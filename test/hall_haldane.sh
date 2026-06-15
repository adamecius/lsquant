#!/usr/bin/env bash
# Hall conductivity / Chern number on the Haldane model (Phase 4 #5).
#
# The Haldane model (t1=-1, t2=0.15, phi=pi/2) is a Chern insulator: with the Fermi level
# in the gap the Hall conductivity is quantized at sigma_xy = C*e^2/h with Chern C = +1.
# This is LinQT's finite, quantized "conductance quantum" test and the end-to-end exercise
# of the Kubo-Bastin route (which is NaN-free only with the alpha safeguard, see
# docs/spectral_scales.md).
#
# Method (oracle = the topological invariant, not another run):
#   1. EXACT-TRACE VX-VY Chebyshev moments from the committed operators (local state over the
#      full 512-dim basis -> deterministic, machine precision, no random vectors).
#   2. Kubo-Bastin reconstruction (inline_kuboBastinFromChebmom).
#   3. Assert the reconstruction is finite, and the gap plateau gives
#        C = (plateau / A) * 2pi,   A = N_cells * A_cell,  A_cell = 3*sqrt(3)/2,  N_cells = dim/2
#      equal to +1 within a finite-size/finite-M tolerance.
#
# Units: e = hbar = 1, so e^2/h = 1/(2pi); the (plateau/A)*2pi factor returns the pure Chern.
#
# Usage: hall_haldane.sh <BUILD_DIR> [M] [CHERN_TOL] [GAP_HALFWIDTH]
set -euo pipefail
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BUILD="$(cd "${1:-$REPO/build}" && pwd)"
M="${2:-128}"
TOL="${3:-0.05}"
GAP_HW="${4:-0.5}"
OPDIR="$REPO/test/golden/haldane/operators"
NONEQ="$BUILD/inline_compute-kpm-nonEqOp"
BASTIN="$BUILD/inline_kuboBastinFromChebmom"

[ -x "$NONEQ" ]  || { echo "FAIL: nonEqOp driver missing ($NONEQ)"; exit 1; }
[ -x "$BASTIN" ] || { echo "FAIL: Kubo-Bastin driver missing ($BASTIN)"; exit 1; }
[ -f "$OPDIR/haldane.HAM.CSR" ] || { echo "FAIL: committed Haldane operators missing"; exit 1; }

N=$(head -1 "$OPDIR/haldane.HAM.CSR" | awk '{print $1}')   # SystemSize = dim = 2*N_cells

T="$(mktemp -d)"; trap 'rm -rf "$T"' EXIT
mkdir -p "$T/operators"
cp "$OPDIR"/haldane.{HAM,VX,VY}.CSR "$T/operators/"
( cd "$T"
  printf -- "-4.5 4.5\n" > BOUNDS                       # encloses the Haldane spectrum (|E|max~3)
  { echo local; echo "$N"; seq 0 $((N-1)); } > exact    # exact trace over the full basis
  echo ">> VX-VY moments (M=$M, exact trace, dim=$N)"
  "$NONEQ" haldane VX VY "$M" exact >noneq.log 2>&1 \
    || { echo "FAIL: moment computation crashed"; tail -15 noneq.log; exit 1; }
  MOM="$(ls NonEqOpVX-VYhaldane*chebmom2D | head -1)"
  "$BASTIN" "$MOM" 10 >bastin.log 2>&1 \
    || { echo "FAIL: Kubo-Bastin reconstruction crashed"; tail -15 bastin.log; exit 1; }
  DAT="$(ls KuboBastin*haldane*.dat | head -1)"

  NANS=$(grep -ciE 'nan|inf' "$DAT" || true)
  [ "$NANS" -eq 0 ] || { echo "FAIL: $NANS NaN/Inf in the Hall reconstruction (alpha safeguard inactive?)"; exit 1; }

  awk -v n="$N" -v tol="$TOL" -v gw="$GAP_HW" '
    function abs(x){return x<0?-x:x}
    abs($1) < gw { s += $2; c++ }
    END {
      if (c == 0) { print "FAIL: no points inside the gap |E|<" gw; exit 1 }
      Acell = 3*sqrt(3)/2; Ncells = n/2; A = Ncells*Acell;
      plateau = s/c;
      C = (plateau/A) * 2 * 3.141592653589793;
      printf "gap plateau sigma_xy = %.4f  (mean over %d pts)\n", plateau, c;
      printf "A = N_cells*A_cell = %d*%.4f = %.3f\n", Ncells, Acell, A;
      printf "Chern C = (plateau/A)*2pi = %+.4f   (target +1, tol %.3f)\n", C, tol;
      if (abs(C - 1.0) <= tol) { print "PASS: quantized Hall plateau, Chern = +1"; exit 0 }
      else { print "FAIL: Chern off target beyond tolerance"; exit 1 }
    }' "$DAT"
)
