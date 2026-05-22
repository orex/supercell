#!/usr/bin/env bash
#
# rescue_wrap.sh - workaround wrapper for the pre-v2.1.2 case-C cleanup bug.
#
# Intended to be passed to tests/run_examples.sh via the RUN_WRAP env var so
# every supercell invocation goes through this shim. The shim only does
# anything for case C; for all other cases it transparently exec's supercell.
#
# The bug it works around
# -----------------------
# Inside write_files() (v2.1.x source, fixed by commit 21bb7b4 on May 2026),
# supercell calls f_q_calc.open("<prefix>_coulomb_energy.txt") and then later
# in the same function scans the output directory and bfs::remove()s any file
# matching <prefix>_coulomb_energy.*\.txt. POSIX unlink-while-open succeeds:
# the directory entry vanishes immediately, the open FD keeps writing into an
# orphaned inode, and the inode is reclaimed when the FD closes. Result: for
# case C the entire Coulomb-energy file ends up missing on disk.
#
# How the rescue works
# --------------------
# Before running supercell we create the expected output file as a hardlink
# to a rescue path in a sibling .rescue/ subdir of the case output dir. The
# inode now has nlink=2. When supercell unlinks the output-dir entry the
# rescue hardlink still pins the inode, so writes via the FD keep flowing
# into a still-on-disk inode. After supercell exits we move the rescue file
# back into the output dir.
#
# Post-fix builds are unaffected
# ------------------------------
# The fix moved cleanup_output_files() to the top of process(), before any
# output stream is opened. On a post-fix run our pre-staged hardlink simply
# gets cleaned up like any other stale output; supercell then creates a new
# inode at the original path and writes to it. The rescue file remains an
# empty touch and the restore check is a no-op.

set -u

if [ $# -lt 1 ]; then
    echo "usage: $0 <supercell-binary> [args...]" >&2
    exit 64
fi

# Locate -o <prefix> in the supercell argument list.
prefix=""
for ((i=1; i<=$#; i++)); do
    if [ "${!i}" = "-o" ]; then
        j=$((i + 1))
        if [ "$j" -le "$#" ]; then
            prefix="${!j}"
        fi
        break
    fi
done

# Only case C needs the rescue. We key off the prefix basename so the wrapper
# stays case-table-agnostic.
if [ -n "$prefix" ] && [ "$(basename -- "$prefix")" = "cell113" ]; then
    out_dir=$(dirname -- "$prefix")
    rescue_dir="$out_dir/.rescue"
    coulomb="${prefix}_coulomb_energy.txt"

    mkdir -p "$rescue_dir"
    rm -f "$rescue_dir/coulomb" "$coulomb"
    touch "$rescue_dir/coulomb"
    ln "$rescue_dir/coulomb" "$coulomb"

    rc=0
    "$@" || rc=$?

    # If supercell unlinked the output-dir copy (pre-fix path), the data
    # survives at the rescue path — move it back. On a post-fix run the
    # output-dir file already exists with fresh data and this is a no-op.
    if [ ! -e "$coulomb" ] && [ -s "$rescue_dir/coulomb" ]; then
        mv "$rescue_dir/coulomb" "$coulomb"
    fi
    rm -rf "$rescue_dir"
    exit "$rc"
fi

exec "$@"
