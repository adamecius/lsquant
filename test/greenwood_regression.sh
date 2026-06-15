#!/usr/bin/env bash
# Kubo-Greenwood reconstruction regression -- pins the Bug-1 fix (Phase 0).
#
# Bug 1 (fixed in PR #19): kuboGreenwoodFromChebmom.cpp multiplied the kernel by
# ScaleFactor()*ShiftFactor(); ShiftFactor() = -BandCenter/HalfWidth*CUTOFF is identically 0
# for any band centred at E=0 (clean chain, graphene, every particle-hole-symmetric system),
# so sigma_KG came out IDENTICALLY ZERO. The fix makes the second factor a second ScaleFactor()
# (prefactor -SystemSize*ScaleFactor^2).
#
# Coverage gap this closes: chain1d_analytic guards the KG *moment matrix* (B.15) and
# graphene_regression guards the *moments* -- but NEITHER exercises the reconstruction
# prefactor, so a regression back to sigma_KG==0 would pass CI silently. This test reconstructs
# sigma_KG(E) from the committed graphene VX-VX moments and asserts it is FINITE and NOT
# identically zero in the band interior (the bug gives exactly 0 everywhere; the fix gives
# O(1)). Reference is the physical fact that a centred metal's longitudinal conductivity is not
# identically zero -- not another run.
#
# Usage: greenwood_regression.sh <BUILD_DIR>
set -euo pipefail
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BUILD="$(cd "${1:-$REPO/build}" && pwd)"
DRV="$BUILD/inline_kuboGreenwoodFromChebmom"
GMOM="$REPO/test/golden/NonEqOpVX-VXgraphene_goldenKPM_M64x64_statedefault.chebmom2D"

[ -x "$DRV" ] || { echo "FAIL: Kubo-Greenwood driver missing ($DRV)"; exit 1; }
[ -f "$GMOM" ] || { echo "FAIL: committed graphene VX-VX moments missing ($GMOM)"; exit 1; }

T="$(mktemp -d)"; trap 'rm -rf "$T"' EXIT
cp "$GMOM" "$T/"
( cd "$T"
  MOM="$(ls *.chebmom2D | head -1)"
  "$DRV" "$MOM" 10 >driver.log 2>&1 \
    || { echo "FAIL: Kubo-Greenwood reconstruction crashed"; cat driver.log; exit 1; }
  DAT="$(ls KuboGreenWood_*.dat 2>/dev/null | head -1)"
  [ -n "$DAT" ] || { echo "FAIL: no reconstruction .dat produced"; cat driver.log; exit 1; }

  # Finite everywhere, and a non-negligible interior magnitude (the bug => exactly 0).
  NANS=$(awk 'tolower($2) ~ /nan|inf/ {c++} END{print c+0}' "$DAT")
  MAXABS=$(awk 'function abs(x){return x<0?-x:x}
                tolower($2) !~ /nan|inf/ && abs($1) < 2.0 {if(abs($2)>m)m=abs($2)} END{printf "%.6g", m+0}' "$DAT")
  echo "interior max|sigma_KG| = $MAXABS   (NaN/Inf rows = $NANS)"
  [ "$NANS" -eq 0 ] || { echo "FAIL: $NANS NaN/Inf rows in the Greenwood reconstruction"; exit 1; }
  awk -v m="$MAXABS" 'BEGIN{ if (m > 1e-3) exit 0; else exit 1 }' \
    || { echo "FAIL: sigma_KG is ~0 in the interior (max|sigma_KG|=$MAXABS) -- Bug 1 regressed (ShiftFactor)"; exit 1; }
  echo "PASS: Kubo-Greenwood sigma_KG(E) is finite and non-zero (Bug 1 fix intact)"
)
