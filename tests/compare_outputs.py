#!/usr/bin/env python3
"""
compare_outputs.py - compare two or more archives produced by run_examples.sh.

For every file that is present in all archives we run two passes:

  1. Binary pass    : after stripping volatile lines (banner, command line,
                      random seed, total enumeration time) and normalising
                      CRLF -> LF and the embedded output-path prefix, the
                      file content must be byte identical across archives.
  2. Almost pass    : same as binary, except any numeric token may differ
                      within a relative tolerance (--rtol, default 1e-9) or
                      an absolute tolerance (--atol, default 1e-12). This is
                      wide enough to absorb the "1 unit in the last printed
                      decimal place" drift that legitimately happens between
                      builds, and tight enough that a real numerical bug
                      will still trip the check. For coulomb-energy files
                      the comparison is multiset-based: two lines with the
                      same basename and a tolerance-close energy are
                      considered a tie regardless of order (this lets us
                      tolerate tie-break order swaps).

Exit codes / RESULT line
  0  RESULT: binary_equivalent
  1  RESULT: almost_equivalent
  2  RESULT: not_equivalent
  3  RESULT: error                       (something blew up - bad archive etc)

The final line on stdout is always "RESULT: <label> tags=<a>,<b>,..." so the
script is easy to consume from CI / shell pipelines.

Standard library only (Python 3.8+).
"""

import argparse
import math
import os
import re
import sys
import tarfile
import tempfile
from pathlib import Path

# ---------------------------------------------------------------------------
# constants
# ---------------------------------------------------------------------------

RESULT_BINARY = "binary_equivalent"
RESULT_ALMOST = "almost_equivalent"
RESULT_DIFF = "not_equivalent"
RESULT_ERROR = "error"

EXIT = {
    RESULT_BINARY: 0,
    RESULT_ALMOST: 1,
    RESULT_DIFF: 2,
    RESULT_ERROR: 3,
}

# Volatile log lines that get stripped before comparison.
LOG_DROP_PATTERNS = [
    re.compile(r"^\s*-{5,}\s*$"),                       # banner ruler
    re.compile(r"^\s*-\s.*Supercell program.*-\s*$"),   # version line
    re.compile(r"^\s*-\s*https?://.*-\s*$"),            # url line
    re.compile(r"^\s*-\s*Authors:.*-\s*$"),
    re.compile(r"^\s*-\s+\*\s.*-\s*$"),                 # author bullets
    re.compile(r"^\s*-\s+\(.*@.*\).*-\s*$"),            # author email lines
    re.compile(r"^\s*-\s*please cite:\s*-\s*$"),
    re.compile(r"^\s*-\s+[A-Z]\.\s.*-\s*$"),            # citation author line
    re.compile(r"^\s*-\s+J\..*-\s*$"),                  # citation journal line
    re.compile(r"^Command line:.*$"),
    re.compile(r"^Total enumeration time:.*$"),
    # Diagnostic line printed by cleanup_output_files() when -v 2 is set.
    # Pre-fix v2.1.x emits it inside write_files (after Coulomb calc); post-fix
    # builds emit it at the top of process() (before Coulomb starts). Same
    # information either way, different position in the log — strip it so the
    # run.log compares clean across the cleanup-order fix.
    re.compile(r"^Total \d+ output files? was deleted successfully\s*$"),
]

# Filenames that are pure timing/diagnostic noise - not compared.
IGNORE_FILENAMES = {"time.txt"}

# Regex used both to find numeric tokens and to recognise the energy-line
# inside CIFs.  Captures a signed float / int, optionally with exponent.
NUM_RE = re.compile(r"[-+]?\d+\.\d+(?:[eE][-+]?\d+)?|[-+]?\d+(?:[eE][-+]?\d+)?")

# ---------------------------------------------------------------------------
# helpers - tar extraction, file classification
# ---------------------------------------------------------------------------


def extract_archive(tar_path: Path, dest: Path) -> Path:
    """Extract tar_path into dest. Return the single top-level directory."""
    with tarfile.open(tar_path, "r:*") as tf:
        # safety: refuse absolute paths / parent traversal
        for m in tf.getmembers():
            if m.name.startswith("/") or ".." in Path(m.name).parts:
                raise ValueError(f"unsafe path in archive: {m.name!r}")
        tf.extractall(dest)
    entries = [p for p in dest.iterdir() if p.is_dir()]
    if len(entries) != 1:
        raise ValueError(
            f"{tar_path.name}: expected exactly one top-level dir, got {len(entries)}"
        )
    return entries[0]


def walk_files(root: Path):
    """Yield paths relative to root for every regular file under root."""
    for path in sorted(root.rglob("*")):
        if path.is_file():
            yield path.relative_to(root)


def classify(rel: Path) -> str:
    """Return one of: 'log', 'cif', 'coulomb', 'ignore', 'other'."""
    name = rel.name
    if name in IGNORE_FILENAMES:
        return "ignore"
    if name == "run.log":
        return "log"
    if name.endswith(".cif"):
        return "cif"
    if "coulomb_energy" in name and name.endswith(".txt"):
        return "coulomb"
    return "other"


# ---------------------------------------------------------------------------
# normalisation - turns raw bytes into a canonical comparable form
# ---------------------------------------------------------------------------


def normalise_text(raw: bytes) -> str:
    """Decode bytes (utf-8 with replacement) and normalise CRLF to LF."""
    return raw.decode("utf-8", errors="replace").replace("\r\n", "\n").replace("\r", "\n")


def strip_log(text: str) -> str:
    """Drop volatile log lines (banner, command line, seed, total time)."""
    out = []
    for line in text.splitlines():
        if any(p.match(line) for p in LOG_DROP_PATTERNS):
            continue
        out.append(line)
    return "\n".join(out) + "\n"


def strip_coulomb_paths(text: str) -> str:
    """Replace the path column with just the basename.

    Each line is '<full path>\\t<value> eV'. Different runs embed different
    output directories, so we strip the path down to the leaf filename.
    """
    out = []
    for line in text.splitlines():
        if "\t" in line:
            head, rest = line.split("\t", 1)
            out.append(os.path.basename(head) + "\t" + rest)
        else:
            out.append(line)
    return "\n".join(out) + "\n"


def canonical_text(rel: Path, raw: bytes) -> str:
    """Produce the text we compare for the binary pass."""
    text = normalise_text(raw)
    kind = classify(rel)
    if kind == "log":
        return strip_log(text)
    if kind == "coulomb":
        return strip_coulomb_paths(text)
    return text


# ---------------------------------------------------------------------------
# numeric comparison
# ---------------------------------------------------------------------------


def numbers_close(a: str, b: str, rtol: float, atol: float) -> bool:
    try:
        fa = float(a)
        fb = float(b)
    except ValueError:
        return False
    return math.isclose(fa, fb, rel_tol=rtol, abs_tol=atol)


def line_almost_equal(la: str, lb: str, rtol: float, atol: float) -> bool:
    """True if la and lb match once numeric tokens are compared by tolerance."""
    if la == lb:
        return True
    sa = NUM_RE.split(la)
    sb = NUM_RE.split(lb)
    if sa != sb:
        return False
    na = NUM_RE.findall(la)
    nb = NUM_RE.findall(lb)
    if len(na) != len(nb):
        return False
    return all(numbers_close(x, y, rtol, atol) for x, y in zip(na, nb))


def almost_equal_text(a: str, b: str, rtol: float, atol: float) -> bool:
    """Line-by-line tolerance-aware comparison."""
    la, lb = a.splitlines(), b.splitlines()
    if len(la) != len(lb):
        return False
    return all(line_almost_equal(x, y, rtol, atol) for x, y in zip(la, lb))


def almost_equal_coulomb(a: str, b: str, rtol: float, atol: float) -> bool:
    """Multiset comparison for coulomb energy tables.

    Each line is '<basename>\\t<value> eV'. We pair lines from a and b by
    basename, then check the energies are tolerance-close. Ordering inside
    the file doesn't matter (tie-break order can swap legitimately).
    """
    def parse(text):
        rows = []
        for line in text.splitlines():
            if not line.strip():
                continue
            if "\t" not in line:
                return None
            name, rest = line.split("\t", 1)
            m = NUM_RE.search(rest)
            if not m:
                return None
            rows.append((name, float(m.group(0))))
        return rows

    ra = parse(a)
    rb = parse(b)
    if ra is None or rb is None:
        return almost_equal_text(a, b, rtol, atol)
    if len(ra) != len(rb):
        return False
    from collections import defaultdict
    bucket_b = defaultdict(list)
    for name, val in rb:
        bucket_b[name].append(val)
    for name, val in ra:
        bucket = bucket_b.get(name)
        if not bucket:
            return False
        idx = None
        for i, v in enumerate(bucket):
            if math.isclose(val, v, rel_tol=rtol, abs_tol=atol):
                idx = i
                break
        if idx is None:
            return False
        bucket.pop(idx)
    return all(not v for v in bucket_b.values())


# ---------------------------------------------------------------------------
# per-file comparison
# ---------------------------------------------------------------------------

STATUS_EQUAL = "equal"
STATUS_ALMOST = "almost"
STATUS_DIFF = "differ"
STATUS_MISSING = "missing"


def compare_file(rel: Path, contents: list, rtol: float, atol: float,
                 strict_order: bool = False) -> str:
    """contents is a list of bytes (one per archive). Return a STATUS_*."""
    if any(c is None for c in contents):
        return STATUS_MISSING

    canon = [canonical_text(rel, c) for c in contents]
    ref = canon[0]
    if all(c == ref for c in canon[1:]):
        return STATUS_EQUAL

    kind = classify(rel)
    if kind == "coulomb" and not strict_order:
        ok = all(almost_equal_coulomb(ref, c, rtol, atol) for c in canon[1:])
    else:
        ok = all(almost_equal_text(ref, c, rtol, atol) for c in canon[1:])
    return STATUS_ALMOST if ok else STATUS_DIFF


# ---------------------------------------------------------------------------
# top-level driver
# ---------------------------------------------------------------------------


def main() -> int:
    ap = argparse.ArgumentParser(description="Compare supercell output archives.")
    ap.add_argument("archives", nargs="+", type=Path,
                    help="two or more .tar.gz files produced by run_examples.sh")
    ap.add_argument("--rtol", type=float, default=1e-9,
                    help="relative tolerance for almost-equal numbers (default: 1e-9)")
    ap.add_argument("--atol", type=float, default=1e-12,
                    help="absolute tolerance for almost-equal numbers (default: 1e-12)")
    ap.add_argument("--strict-order", action="store_true",
                    help="require line-by-line ordering even for "
                         "coulomb-energy files (default is multiset, which "
                         "tolerates legitimate tie-break swaps)")
    ap.add_argument("--quiet", action="store_true",
                    help="suppress per-file detail rows")
    ap.add_argument("--show-equal", action="store_true",
                    help="include equal files in detail output (off by default)")
    args = ap.parse_args()

    if len(args.archives) < 2:
        print("error: need at least two archives", file=sys.stderr)
        print(f"RESULT: {RESULT_ERROR} tags=")
        return EXIT[RESULT_ERROR]

    for p in args.archives:
        if not p.is_file():
            print(f"error: archive not found: {p}", file=sys.stderr)
            print(f"RESULT: {RESULT_ERROR} tags=")
            return EXIT[RESULT_ERROR]

    with tempfile.TemporaryDirectory(prefix="sc_compare_") as tmp:
        tmp = Path(tmp)
        roots = []
        tags = []
        try:
            for i, archive in enumerate(args.archives):
                sub = tmp / f"a{i}"
                sub.mkdir()
                root = extract_archive(archive, sub)
                roots.append(root)
                tags.append(root.name)
        except (tarfile.TarError, ValueError, OSError) as exc:
            print(f"error extracting archive: {exc}", file=sys.stderr)
            print(f"RESULT: {RESULT_ERROR} tags=")
            return EXIT[RESULT_ERROR]

        # union of all relative file paths, minus ignored noise.
        all_rels = set()
        for root in roots:
            for rel in walk_files(root):
                if classify(rel) == "ignore":
                    continue
                all_rels.add(rel)

        rows = []
        n_equal = n_almost = n_diff = n_missing = 0
        for rel in sorted(all_rels):
            contents = []
            for root in roots:
                p = root / rel
                contents.append(p.read_bytes() if p.is_file() else None)
            status = compare_file(rel, contents, args.rtol, args.atol,
                                  strict_order=args.strict_order)
            rows.append((rel, status))
            if status == STATUS_EQUAL:
                n_equal += 1
            elif status == STATUS_ALMOST:
                n_almost += 1
            elif status == STATUS_DIFF:
                n_diff += 1
            else:
                n_missing += 1

        # ------------------------------------------------------------------
        # human readable report
        # ------------------------------------------------------------------
        print(f"archives  : {len(args.archives)}")
        for tag, archive in zip(tags, args.archives):
            print(f"  {tag:14s}  {archive}")
        print(f"files     : {len(all_rels)}  "
              f"(equal={n_equal} almost={n_almost} differ={n_diff} missing={n_missing})")
        print(f"tolerance : rtol={args.rtol:g} atol={args.atol:g}"
              + (" strict-order" if args.strict_order else ""))

        if not args.quiet:
            print()
            print(f"{'STATUS':8s} {'KIND':8s} FILE")
            print("-" * 60)
            for rel, status in rows:
                if status == STATUS_EQUAL and not args.show_equal:
                    continue
                print(f"{status:8s} {classify(rel):8s} {rel}")

    # ------------------------------------------------------------------
    # overall verdict
    # ------------------------------------------------------------------
    if n_diff or n_missing:
        verdict = RESULT_DIFF
    elif n_almost:
        verdict = RESULT_ALMOST
    else:
        verdict = RESULT_BINARY

    print()
    print(f"RESULT: {verdict} tags={','.join(tags)}")
    return EXIT[verdict]


if __name__ == "__main__":
    sys.exit(main())
