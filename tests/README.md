# supercell regression tests

Two scripts that exercise a supercell binary on a fixed set of tutorial
inputs and check that two binaries produce equivalent output. They are
designed to run from CI but work the same way locally.

| script              | role                                                                 |
|---------------------|----------------------------------------------------------------------|
| `run_examples.sh`   | runs six cases (A..F) and packs the outputs into a single `tar.gz`   |
| `compare_outputs.py`| compares two or more such archives and prints a machine-readable verdict |

## Cases

All six cases come from `supercell_tutorial.pdf`. Every case is invoked with
`--random-seed 42` so the `Random SEED:` line in `run.log` is byte
deterministic across runs — even on the enumeration-only cases that do not
otherwise consume the seed.

| ID | Input                                | Args                                                                    |
|----|--------------------------------------|-------------------------------------------------------------------------|
| A  | Ca2Al2SiO7/Ca2Al2SiO7.cif            | `-s 1x1x2 -m`                                                           |
| B  | Ca2Al2SiO7/Ca2Al2SiO7.cif            | `-s 2x2x2 -m -n r100 -v 1`                                              |
| C  | gamma-Fe2O3/Fe2O3-P4332.cif          | `-s 1x1x3 -m -q -g -v 2`                                                |
| D  | gamma-Fe2O3/Fe2O3-P4332.cif          | `-s 1x2x3 -m -q -v 2 -n l50 -n r100 -n h20`                             |
| E  | alpha-SiGeO2/alpha-SiGeO2.cif        | `-s 1x1x2 -p Si1:p=2 -p Ge1:p=4 -m`                                     |
| F  | PZT/PZT-PbZr05Ti05O3.cif             | `-s 4x2x1 -m`                                                           |

## Running

```sh
# 1) Run a binary on every case and pack the result.
tests/run_examples.sh <supercell-bin> <prefix> data/examples [out-dir]

# 2) Compare two archives (or more).
python3 tests/compare_outputs.py out-dir/a.tar.gz out-dir/b.tar.gz
```

`run_examples.sh` writes `<out-dir>/<prefix>.tar.gz` and a working
directory `<out-dir>/<prefix>/<CASE>/` containing one `run.log`,
`time.txt` and the structure / coulomb-energy files supercell emits.

Windows `.exe` binaries are auto-detected and launched through `wine`.

### `RUN_WRAP` environment variable

A space-separated wrapper that is prepended to every invocation. The
parallel test uses it to pin the supercell process to a fixed number of
cores via `taskset`:

```sh
RUN_WRAP="taskset -c 0-2" tests/run_examples.sh ...   # 3 cores
```

TBB picks the affinity mask up via `sched_getaffinity`, so this controls
the worker count without any code change in supercell.

## The comparator

```sh
python3 tests/compare_outputs.py a.tar.gz b.tar.gz [c.tar.gz ...] \
        [--rtol 1e-9] [--atol 1e-12] [--strict-order] [--quiet]
```

For every file present in all archives the comparator does two passes:

1. **Binary pass.** After stripping a known set of volatile lines from
   `run.log` (the banner, the `Command line:` echo, `Total enumeration
   time:`), normalising `CRLF -> LF` and replacing the embedded
   output-directory prefix inside `*_coulomb_energy_*.txt` with the bare
   basename, the file content must be byte identical across all
   archives. `time.txt` is ignored entirely.
2. **Almost pass.** Any numeric token may differ within `--rtol` (default
   `1e-9`) or `--atol` (default `1e-12`) — wide enough to absorb the
   "1 unit in the last printed decimal place" drift that legitimately
   happens between builds, tight enough that a real numerical bug still
   trips the check.

For `*_coulomb_energy_*.txt` the almost pass is **multiset based**: two
lines with the same basename and a tolerance-close energy match
regardless of order. This tolerates legitimate tie-break order swaps in
`_h.txt`. `--strict-order` disables that and reverts to line-by-line
comparison — the parallel test uses it.

### Verdict and exit code

The final line of stdout is always

```
RESULT: <label> tags=<a>,<b>,...
```

with the matching exit code:

| label                | exit | meaning                                                |
|----------------------|------|--------------------------------------------------------|
| `binary_equivalent`  | 0    | all files match byte-for-byte after canonicalisation   |
| `almost_equivalent`  | 1    | only numeric drift within tolerance / tie-break swaps  |
| `not_equivalent`     | 2    | at least one file differs structurally or is missing   |
| `error`              | 3    | could not extract / read an archive                    |

CI gates are typically `exit_code <= 1` ("at least almost"). The
parallel test requires `exit_code == 0` with `--strict-order` — anything
weaker would mean two runs of the same binary disagree on a single bit,
which is exactly what we are checking against.

## How CI uses these scripts

The GitHub Actions workflow has three groups of jobs:

* **`reference-data`** — downloads supercell v2.1.0 from the GitHub
  release page, runs `run_examples.sh` on it, then installs the latest
  supercell from snap and runs it again. The two archives must compare
  at least `almost_equivalent`. The v2.1.0 archive is then uploaded as
  the artifact `reference-test-data` for the downstream jobs.
* **All deploy and develop builds** — download `reference-test-data`, run
  `run_examples.sh` against the freshly built binary, and compare. The
  job passes only if the verdict is `almost_equivalent` or stricter.
* **`parallel-test`** — builds with default flags, then runs the test
  four times under `taskset -c 0..N-1` for N=1..4. All four archives
  must compare `binary_equivalent` with `--strict-order`. This guards
  against thread-count-dependent non-determinism, including last-bit
  floating-point drift.

## Local quick check

```sh
# build supercell first, then:
cd /path/to/supercell
tests/run_examples.sh build/src/sc_cli/supercell run1 data/examples /tmp/sc
tests/run_examples.sh build/src/sc_cli/supercell run2 data/examples /tmp/sc
python3 tests/compare_outputs.py /tmp/sc/run1.tar.gz /tmp/sc/run2.tar.gz
# expected: RESULT: binary_equivalent tags=run1,run2  (exit 0)
```
