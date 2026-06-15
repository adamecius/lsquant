#!/usr/bin/env bash
# Phase 4 -- regression for the CorrOp state-file bug.
#
# inline_compute-kpm-CorrOp used to do `numStates = atoi(argv[7])` and never load the state
# file (a dead argc==5 branch), so passing "exact" silently became 0 and the run fell back to a
# single default random vector -- which is why a deterministic VAC<->MSD golden was impossible.
# The fix loads argv[7] as a STATE FILE (like the MSD/nonEqOp drivers). This asserts the exact-
# trace state spec is honored: the output is tagged ...stateexact... (the bug tagged it
# ...statedefault...) and is byte-reproducible across two runs (exact trace is deterministic).
#
# Usage: corrop_state_honored.sh <BUILD_DIR>
set -euo pipefail
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BUILD="$(cd "${1:-$REPO/build}" && pwd)"
DRV="$BUILD/inline_compute-kpm-CorrOp"
OPD="$REPO/test/golden/chain1d_ballistic/operators"

[ -x "$DRV" ] || { echo "FAIL: CorrOp driver missing"; exit 1; }
[ -f "$OPD/chain1d.HAM.CSR" ] && [ -f "$OPD/chain1d.VX.CSR" ] || { echo "FAIL: chain operators missing"; exit 1; }

N=$(head -1 "$OPD/chain1d.HAM.CSR" | awk '{print $1}')
run() {  # runs CorrOp in a fresh tmp dir; echoes ONLY that dir
    local T; T="$(mktemp -d)"; mkdir -p "$T/operators"
    cp "$OPD/chain1d.HAM.CSR" "$OPD/chain1d.VX.CSR" "$T/operators/"
    cp "$OPD/chain1d.HAM.desc" "$T/operators/" 2>/dev/null || true
    ( cd "$T"
      { echo local; echo "$N"; seq 0 $((N-1)); } > exact
      "$DRV" chain1d VX VX 32 16 10 exact >/dev/null 2>&1 )
    echo "$T"
}
T1="$(run)"; f1="$(cd "$T1" && ls TimeCorr*chebmomTD)"
T2="$(run)"; f2="$(cd "$T2" && ls TimeCorr*chebmomTD)"
trap 'rm -rf "$T1" "$T2"' EXIT

echo "produced: $f1"
case "$f1" in
  *stateexact*) echo "PASS(honored): exact-trace state file is loaded (output tagged stateexact)";;
  *statedefault*) echo "FAIL: state file ignored -- output is statedefault (atoi bug regressed)"; exit 1;;
  *) echo "FAIL: unexpected output name '$f1'"; exit 1;;
esac

if diff -q "$T1/$f1" "$T2/$f2" >/dev/null; then
    echo "PASS(deterministic): exact-trace VAC moments are byte-identical across runs"
else
    echo "FAIL: exact-trace VAC moments not reproducible"; exit 1
fi
echo "PASS: CorrOp honors the state file (atoi bug fixed)"
