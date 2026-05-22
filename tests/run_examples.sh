#!/usr/bin/env bash
#
# run_examples.sh - run supercell on the six tutorial cases A..F (see
#                   compare.tar/SUMMARY.md) and pack the outputs into a
#                   single tar.gz archive that the comparator script can
#                   later diff against archives produced by other binaries.
#
# Usage:
#     run_examples.sh <supercell-binary> <prefix> <examples-dir> [output-dir]
#
# Arguments:
#     supercell-binary  Path to the supercell executable to test. If it ends
#                       in .exe it is launched via `wine`.
#     prefix            Tag for this run (e.g. "deps_bump", "windows").
#                       Used both as the top-level directory inside the tar
#                       and as the archive name <prefix>.tar.gz.
#     examples-dir      Path to the directory that holds the supercell
#                       tutorial inputs (typically <repo>/data/examples).
#     output-dir        Optional. Where to place the tar.gz and the working
#                       tree. Defaults to the current directory.
#
# Environment:
#     RUN_WRAP          Optional. A command (and arguments, space-separated)
#                       to prepend to every supercell invocation. Used by the
#                       parallel test to pin core count via
#                         RUN_WRAP="taskset -c 0-2"
#                       which restricts TBB to three cores via CPU affinity.
#
# Result: <output-dir>/<prefix>.tar.gz containing
#     <prefix>/<CASE>/run.log    - combined stdout+stderr of the run
#     <prefix>/<CASE>/time.txt   - elapsed/user/sys/maxrss/rc
#     <prefix>/<CASE>/...        - all output files produced by supercell

set -u

if [[ $# -lt 3 || $# -gt 4 ]]; then
    echo "usage: $0 <supercell-binary> <prefix> <examples-dir> [output-dir]" >&2
    exit 64
fi

BIN_RAW=$1
PREFIX=$2
EX_RAW=$3
OUT_RAW=${4:-.}

# -- resolve absolute paths -------------------------------------------------
# We deliberately do NOT follow symlinks: /snap/bin/<app> is a symlink to
# /usr/bin/snap that re-execs itself based on argv[0], so resolving the
# symlink would turn "/snap/bin/supercell" into "/usr/bin/snap" and the
# snap dispatcher wouldn't know which snap to launch.
abspath() {
    case "$1" in
        /*) printf '%s\n' "$1" ;;
        *)  printf '%s/%s\n' "$PWD" "$1" ;;
    esac
}

BIN=$(abspath "$BIN_RAW")
EX=$(abspath "$EX_RAW")
mkdir -p -- "$OUT_RAW"
OUT=$(abspath "$OUT_RAW")

if [[ ! -x "$BIN" && "${BIN##*.}" != "exe" ]]; then
    echo "error: $BIN is not executable" >&2
    exit 66
fi
if [[ ! -d "$EX" ]]; then
    echo "error: examples dir $EX not found" >&2
    exit 66
fi

# Optional wrapper (e.g. "taskset -c 0-1") is split on whitespace and
# prepended to every invocation. Empty by default. The `${arr[@]+...}`
# idiom is needed because macOS still ships bash 3.2, which treats an
# empty array as "unset" under set -u and errors on a bare ${arr[@]}.
# shellcheck disable=SC2206
WRAP=(${RUN_WRAP:-})

# binaries ending in .exe go through wine
RUN=(${WRAP[@]+"${WRAP[@]}"} "$BIN")
if [[ "${BIN##*.}" == "exe" ]]; then
    if ! command -v wine >/dev/null; then
        echo "error: $BIN is a Windows binary but 'wine' is not on PATH" >&2
        exit 69
    fi
    RUN=(${WRAP[@]+"${WRAP[@]}"} wine "$BIN")
fi

# Locate GNU time (-f support). macOS ships BSD time, so prefer `gtime`
# (brew install gnu-time) when /usr/bin/time doesn't grok -f.
GNU_TIME=""
if /usr/bin/time -f '' true >/dev/null 2>&1; then
    GNU_TIME=/usr/bin/time
elif command -v gtime >/dev/null 2>&1; then
    GNU_TIME=$(command -v gtime)
else
    echo "error: GNU time not found (need /usr/bin/time -f support, or gtime)" >&2
    exit 69
fi

WORK="$OUT/$PREFIX"
rm -rf -- "$WORK"
mkdir -p -- "$WORK"

# -- case definitions -------------------------------------------------------
# Each case is: <id>|<input cif (relative to examples-dir)>|<output prefix>|<extra args>
# --random-seed 42 is pinned on every case so the "Random SEED:" line in
# run.log is byte-deterministic across runs of the same binary, even on the
# enumeration-only cases that would otherwise pick a fresh seed every run.
CASES=(
  "A|Ca2Al2SiO7/Ca2Al2SiO7.cif|struct|-s 1x1x2 -m --random-seed 42"
  "B|Ca2Al2SiO7/Ca2Al2SiO7.cif|struct|-s 2x2x2 -m -n r100 -v 1 --random-seed 42"
  "C|gamma-Fe2O3/Fe2O3-P4332.cif|cell113|-s 1x1x3 -m -q -g -v 2 --random-seed 42"
  "D|gamma-Fe2O3/Fe2O3-P4332.cif|cell123|-s 1x2x3 -m -q -v 2 -n l50 -n r100 -n h20 --random-seed 42"
  "E|alpha-SiGeO2/alpha-SiGeO2.cif|SiGeO2_112|-s 1x1x2 -p Si1:p=2 -p Ge1:p=4 -m --random-seed 42"
  "F|PZT/PZT-PbZr05Ti05O3.cif|PZT421|-s 4x2x1 -m --random-seed 42"
)

run_case() {
    local id=$1 input=$2 oprefix=$3 args=$4
    local case_dir="$WORK/$id"
    mkdir -p -- "$case_dir"

    local input_abs="$EX/$input"
    if [[ ! -f "$input_abs" ]]; then
        echo "[$id] SKIP - input $input_abs not found" >&2
        printf 'elapsed=0 user=0 sys=0 maxrss=0\nrc=127\n' >"$case_dir/time.txt"
        return
    fi

    # shellcheck disable=SC2206
    local extra=($args)
    local cmd=("${RUN[@]}" -i "$input_abs" "${extra[@]}" -o "$case_dir/$oprefix")

    printf '[%s] %s\n' "$id" "${cmd[*]}"

    local rc=0
    "$GNU_TIME" -f 'elapsed=%e user=%U sys=%S maxrss=%M' \
        -o "$case_dir/time.txt" \
        "${cmd[@]}" >"$case_dir/run.log" 2>&1 || rc=$?
    printf 'rc=%d\n' "$rc" >>"$case_dir/time.txt"
}

for spec in "${CASES[@]}"; do
    IFS='|' read -r id input oprefix args <<<"$spec"
    run_case "$id" "$input" "$oprefix" "$args"
done

# -- archive ---------------------------------------------------------------
TARBALL="$OUT/$PREFIX.tar.gz"
rm -f -- "$TARBALL"
tar -czf "$TARBALL" -C "$OUT" "$PREFIX"

echo "wrote $TARBALL"
