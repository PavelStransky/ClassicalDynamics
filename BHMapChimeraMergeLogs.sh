#!/bin/bash
# Combines the per-task BHMapChimera.sh logs into two files, once the array
# sweep (or a batch of it) has finished. Safe to re-run; it always rebuilds
# the combined files from whatever per-task logs currently exist.
#
# Usage (from anywhere, e.g. the login node):
#   ./BHMapChimeraMergeLogs.sh

set -euo pipefail

LOGDIR="$HOME/results/bh/lyapunov/5/logs"

shopt -s nullglob
outFiles=("$LOGDIR"/bhmap_*.out)
errFiles=("$LOGDIR"/bhmap_*.err)

if [ ${#outFiles[@]} -eq 0 ] && [ ${#errFiles[@]} -eq 0 ]; then
    echo "No per-task logs found in $LOGDIR"
    exit 0
fi

# Sort numerically by array task id (bhmap_<id>.out/.err) rather than lexically.
sortById() {
    printf '%s\n' "$@" | sed -E 's#.*/bhmap_([0-9]+)\.(out|err)$#\1 &#' | sort -n | cut -d' ' -f2-
}

if [ ${#outFiles[@]} -gt 0 ]; then
    sortById "${outFiles[@]}" | xargs cat > "$LOGDIR/bhmap_all.out"
    echo "Wrote $LOGDIR/bhmap_all.out (${#outFiles[@]} task logs)"
fi

if [ ${#errFiles[@]} -gt 0 ]; then
    sortById "${errFiles[@]}" | xargs cat > "$LOGDIR/bhmap_all.err"
    echo "Wrote $LOGDIR/bhmap_all.err (${#errFiles[@]} task logs)"
fi
