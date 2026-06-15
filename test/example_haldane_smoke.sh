#!/usr/bin/env bash
# Smoke test for examples/03_haldane_hall (matplotlib-free, numpy-free).
#
# Exercises the tutorial's MODERN pipeline end to end through the single lsquant binary:
#   1. make_haldane.py stages the committed operators and writes operators/haldane.HAM.desc
#      (spectral bounds in the .desc; NO BOUNDS file).
#   2. lsquant compute --config run.json computes the VX-VY (Hall) and VX-VX (longitudinal)
#      exact-trace moments, taking the spectral bounds straight from the descriptor.
#   3. lsquant reconstruct ... bastin gives sigma_xy; the gap plateau gives Chern C = +1.
#   4. lsquant reconstruct ... greenwood gives sigma_xx; it vanishes in the gap.
#
# This complements hall_haldane.sh (which drives the inline binaries on the golden operators
# with a BOUNDS file) by pinning the example wrapper, the .desc-driven bounds, and the single
# binary together. Oracle is the topological invariant, not another run.
#
# Usage: example_haldane_smoke.sh <BUILD_DIR> [M] [CHERN_TOL] [GAP_HALFWIDTH]
set -euo pipefail
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BUILD="$(cd "${1:-$REPO/build}" && pwd)"
M="${2:-128}"
TOL="${3:-0.05}"
GAP_HW="${4:-0.5}"
EX="$REPO/examples/03_haldane_hall"
PY="${PYTHON:-python3}"
LSQ="$BUILD/lsquant"

[ -x "$LSQ" ] || { echo "FAIL: lsquant binary missing ($LSQ)"; exit 1; }
[ -f "$EX/make_haldane.py" ] || { echo "FAIL: make_haldane.py missing"; exit 1; }
command -v "$PY" >/dev/null 2>&1 || { echo "SKIP: no python3"; exit 0; }

T="$(mktemp -d)"; trap 'rm -rf "$T"' EXIT
cp "$EX/make_haldane.py" "$T/"
# make_haldane.py copies operators from ../../test/golden/haldane relative to its own folder,
# so run it from the example dir but stage into the scratch dir.
( cd "$EX" && HALDANE_OPDIR="$T/operators" "$PY" -c "import make_haldane; make_haldane.build('$T')" ) \
  || { echo "FAIL: make_haldane.build crashed"; exit 1; }

[ -f "$T/operators/haldane.HAM.desc" ] || { echo "FAIL: descriptor not written"; exit 1; }
[ ! -f "$T/BOUNDS" ] || { echo "FAIL: a BOUNDS file was written (the .desc should carry bounds)"; exit 1; }

N=$(head -1 "$T/operators/haldane.HAM.CSR" | awk '{print $1}')

( cd "$T"
  # bounds must resolve from the descriptor (no BOUNDS file present)
  "$LSQ" inspect operators/haldane.HAM.desc 2>&1 | grep -q '\[-4.5, 4.5\]' \
    || { echo "FAIL: lsquant inspect did not report the descriptor bounds"; exit 1; }

  cat > run_hall.json <<JSON
{ "label": "haldane", "operator_right": "VX", "operator_left": "VY",
  "num_moments": $M, "kernel": "jackson", "state": "exact" }
JSON
  cat > run_long.json <<JSON
{ "label": "haldane", "operator_right": "VX", "operator_left": "VX",
  "num_moments": $M, "kernel": "jackson", "state": "exact" }
JSON

  echo ">> lsquant compute (VX-VY, VX-VX; exact trace, dim=$N, bounds from .desc)"
  "$LSQ" compute --config run_hall.json >hall.log 2>&1 \
    || { echo "FAIL: Hall moment computation crashed"; tail -15 hall.log; exit 1; }
  "$LSQ" compute --config run_long.json >long.log 2>&1 \
    || { echo "FAIL: longitudinal moment computation crashed"; tail -15 long.log; exit 1; }
  grep -q 'descriptor' hall.log \
    || { echo "FAIL: compute did not take bounds from the descriptor"; tail -5 hall.log; exit 1; }

  XY="$(ls NonEqOpVX-VYhaldane*chebmom2D | head -1)"
  XX="$(ls NonEqOpVX-VXhaldane*chebmom2D | head -1)"
  "$LSQ" reconstruct "$XY" bastin 10 >bastin.log 2>&1 \
    || { echo "FAIL: Kubo-Bastin reconstruction crashed"; tail -15 bastin.log; exit 1; }
  "$LSQ" reconstruct "$XX" greenwood 10 >green.log 2>&1 \
    || { echo "FAIL: Kubo-Greenwood reconstruction crashed"; tail -15 green.log; exit 1; }
  SXY="$(ls KuboBastin*haldane*.dat | head -1)"
  SXX="$(ls KuboGreenWood*haldane*.dat | head -1)"

  for f in "$SXY" "$SXX"; do
    NANS=$(grep -ciE 'nan|inf' "$f" || true)
    [ "$NANS" -eq 0 ] || { echo "FAIL: $NANS NaN/Inf in $f"; exit 1; }
  done

  # Chern from the Hall gap plateau: C = (plateau/A)*2pi, A = (dim/2)*A_cell, A_cell = 3*sqrt(3)/2
  awk -v n="$N" -v tol="$TOL" -v gw="$GAP_HW" '
    function abs(x){return x<0?-x:x}
    abs($1) < gw { s += $2; c++ }
    END {
      if (c == 0) { print "FAIL: no Hall points inside the gap |E|<" gw; exit 1 }
      Acell = 3*sqrt(3)/2; A = (n/2)*Acell;
      C = (s/c/A) * 2 * 3.141592653589793;
      printf "Chern C = (plateau/A)*2pi = %+.4f   (target +1, tol %.3f)\n", C, tol;
      if (abs(C - 1.0) > tol) { print "FAIL: Hall plateau off target beyond tolerance"; exit 1 }
    }' "$SXY"

  # Longitudinal sigma_xx must vanish in the gap relative to its in-band value.
  awk -v gw="$GAP_HW" '
    function abs(x){return x<0?-x:x}
    { v=abs($2); if (v>peak) peak=v; if (abs($1)<gw){ g+=v; gc++ } }
    END {
      if (gc==0) { print "FAIL: no sigma_xx points in the gap"; exit 1 }
      gapmean = g/gc;
      printf "sigma_xx gap mean = %.4g   in-band peak = %.4g\n", gapmean, peak;
      if (peak<=0)               { print "FAIL: sigma_xx is identically zero"; exit 1 }
      if (gapmean > 0.05*peak)   { print "FAIL: sigma_xx does not vanish in the gap"; exit 1 }
    }' "$SXX"

  echo "PASS: .desc-driven compute, quantized Hall plateau (Chern +1), sigma_xx vanishes in the gap"
)
