#!/usr/bin/env bash
# Graphene regression + spectral-bounds-safeguard gate.
#
# Guards three things from the committed graphene golden (VX.VX Kubo-Greenwood moments):
#
#   (1) DETERMINISM: same KPM_SEED => the two fresh runs are byte-identical.
#   (2) REGRESSION: the fresh moments match the committed golden within TOL (numerical,
#       not byte -- toolchain reduction order differs at ~1e-15). A real src/ refactor
#       bug would move moments at O(1e-2)+, so TOL=1e-6 catches regressions while
#       staying immune to compiler/Eigen/OpenMP noise.
#   (3) BOUNDS SAFEGUARD: with no BOUNDS file, the solver must fall back to the
#       Gershgorin + 10% estimate (NOT the old arbitrary [-100,100]). We assert the
#       estimate is announced and its magnitude is sane (1 < |E|_max < 50 for this
#       fixture; the true Gershgorin radius is ~6.83, padded to ~7.52).
#
# Usage: graphene_reproduce.sh <BUILD_DIR> [TOL]
set -euo pipefail
export KPM_SEED=12345

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BUILD="$(cd "${1:-$REPO/build}" && pwd)"
TOL="${2:-1e-6}"
OPDIR="$REPO/test/golden/operators"
GOLD="$REPO/test/golden/NonEqOpVX-VXgraphene_goldenKPM_M64x64_statedefault.chebmom2D"

[ -f "$OPDIR/graphene_golden.HAM.CSR" ] || { echo "FAIL: committed graphene operators missing"; exit 1; }
[ -f "$GOLD" ] || { echo "FAIL: committed graphene golden missing"; exit 1; }

# echoes "<moments-path>|<solver-log-path>"
gen() {
  local T; T="$(mktemp -d)"; mkdir -p "$T/operators"
  cp "$OPDIR/graphene_golden.HAM.CSR" "$OPDIR/graphene_golden.VX.CSR" "$T/operators/"
  # NO BOUNDS file on purpose -> exercises the Gershgorin safeguard
  ( cd "$T"; "$BUILD/inline_compute-kpm-nonEqOp" graphene_golden VX VX 64 >log 2>&1 )
  echo "$T/$(cd "$T" && ls NonEqOp*graphene*chebmom2D | head -1)|$T/log"
}

A="$(gen)"; B="$(gen)"
AMOM="${A%%|*}"; ALOG="${A##*|}"
BMOM="${B%%|*}"

# (1) determinism
diff -q "$AMOM" "$BMOM" >/dev/null \
  || { echo "FAIL(determinism): two runs differ under fixed KPM_SEED"; exit 1; }
echo "PASS(determinism): two fresh runs are byte-identical"

# (3) bounds safeguard: Gershgorin announced + sane magnitude, and NOT [-100,100]
if ! grep -qi "gershgorin" "$ALOG"; then
  echo "FAIL(bounds): Gershgorin safeguard not announced (no BOUNDS file present)"; cat "$ALOG"; exit 1
fi
HI="$(grep -oE '\-?[0-9]+\.[0-9]+' "$ALOG" | awk '{v=$1<0?-$1:$1; if(v>m)m=v} END{print m+0}')"
awk -v hi="$HI" 'BEGIN{ if (hi>1 && hi<50) exit 0; else exit 1 }' \
  || { echo "FAIL(bounds): estimate |E|max=$HI out of sane range (1,50) -- check for [-100,100] fallback"; exit 1; }
echo "PASS(bounds): Gershgorin+10% safeguard active, |E|max~=$HI (not [-100,100])"

# (2) regression vs committed golden, numerical tolerance
"$BUILD/chebmom_compare" "$AMOM" "$GOLD" "$TOL" \
  || { echo "FAIL(regression): fresh graphene moments differ from committed golden beyond TOL=$TOL"; exit 1; }
echo "PASS(regression): fresh graphene moments match committed golden within $TOL"
echo "PASS: graphene golden is deterministic, regression-clean, and uses the bounds safeguard"
