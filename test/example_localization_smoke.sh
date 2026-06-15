#!/usr/bin/env bash
# Smoke test for examples/02_localization_1d_chain (matplotlib-free, numpy-free).
#
# Runs the tutorial's MSD pipeline on a small strongly-disordered chain (fast localization) and a
# clean chain of the same size, and asserts the disordered MSD plateau at late time is finite and
# well below the clean-chain ballistic value. This pins that disorder localizes the spreading
# (plateau) while the clean chain keeps spreading (MSD ~ t^2), through the real binaries.
#
# Usage: example_localization_smoke.sh <BUILD_DIR>
set -euo pipefail
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BUILD="$(cd "${1:-$REPO/build}" && pwd)"
EX="$REPO/examples/02_localization_1d_chain"
PY="${PYTHON:-python3}"
MSD="$BUILD/inline_compute-kpm-MeanSquareDisplacement"
REC="$BUILD/inline_timeCorrelationsFromChebmom"

[ -x "$MSD" ] && [ -x "$REC" ] || { echo "FAIL: MSD/reconstruction drivers missing"; exit 1; }
[ -f "$EX/make_disordered_chain.py" ] || { echo "FAIL: make_disordered_chain.py missing"; exit 1; }
command -v "$PY" >/dev/null 2>&1 || { echo "SKIP: no python3"; exit 0; }

# late-time MSD(E=0) for a chain built with disorder strength W (N, tmax fixed); echoes the value
late_msd() {  # args: W
  local W="$1" T; T="$(mktemp -d)"
  cp "$EX/make_disordered_chain.py" "$T/"
  ( cd "$T"; export KPM_SEED=1
    "$PY" make_disordered_chain.py 256 "$W" 1 >/dev/null
    L="chain1d_dis_N256_W$(printf '%g' "$W")"
    "$MSD" "$L" VX 48 48 100 >/dev/null 2>&1
    M="$(ls Correlation*chebmomTD | head -1)"
    "$REC" "$(basename "$M")" 10 0.0 >/dev/null 2>&1
    awk 'tolower($2)!~/nan|inf/{v=$2} END{print v}' "$(ls mean*EF*.dat | head -1)" )
  rm -rf "$T"
}

dis="$(late_msd 3)"     # strong disorder -> localizes, MSD plateaus low
clean="$(late_msd 0)"   # no disorder -> ballistic, MSD ~ t^2 stays large
echo "late MSD(E=0, t=100):  disordered(W=3) = $dis   clean(W=0) = $clean"

awk -v d="$dis" -v c="$clean" 'BEGIN{
  if (d ~ /nan|inf/ || d+0 <= 0)      { print "FAIL: disordered MSD plateau is not finite/positive"; exit 1 }
  if (!(d+0 < 0.5*(c+0)))             { print "FAIL: disordered MSD not well below clean ballistic"; exit 1 }
  print "PASS: disorder localizes (finite plateau " d ") far below clean ballistic " c
}'
