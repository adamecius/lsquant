#!/usr/bin/env bash
# Golden-gate for the Kubo-Bastin Sea/Surf reconstruction drivers (plan Step 0.5b).
#
# kuboBastinSeaFromChebmom / kuboBastinSurfFromChebmom carry INLINE reconstruction logic
# (Green-difference cumulative integral for Sea; double-delta POINTWISE accumulation for
# Surf) that is a different algebraic decomposition of Kubo-Bastin than reconstruct_kubo --
# so it cannot be guarded by the unified-kernel goldens. This pins both drivers by
# byte-reproduction of a committed golden BEFORE the A2 migration folds them onto the
# shared A1 grid/integral/write helpers. The kernel loops are parallelized over independent
# energy points (each kernel[i] summed sequentially), so output is byte-deterministic across
# thread counts -- a byte diff is the right gate.
#
# Reconstructs Sea and Surf from the committed graphene VX-VX moments (broadening=10) and
# byte-compares to the committed .dat goldens.
#
# Usage: sea_surf_regression.sh <BUILD_DIR>
set -euo pipefail
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BUILD="$(cd "${1:-$REPO/build}" && pwd)"
GMOM="$REPO/test/golden/NonEqOpVX-VXgraphene_goldenKPM_M64x64_statedefault.chebmom2D"
GDIR="$REPO/test/golden/sea_surf"

[ -f "$GMOM" ] || { echo "FAIL: committed graphene moments missing"; exit 1; }
for d in "$BUILD/inline_kuboBastinSeaFromChebmom" "$BUILD/inline_kuboBastinSurfFromChebmom"; do
  [ -x "$d" ] || { echo "FAIL: driver missing ($d)"; exit 1; }
done

T="$(mktemp -d)"; trap 'rm -rf "$T"' EXIT
cp "$GMOM" "$T/"
( cd "$T"
  M="$(ls *.chebmom2D | head -1)"
  "$BUILD/inline_kuboBastinSeaFromChebmom"  "$M" 10 >/dev/null 2>&1
  "$BUILD/inline_kuboBastinSurfFromChebmom" "$M" 10 >/dev/null 2>&1

  diff -q "$(ls KuboBastinSea_*.dat)"  "$GDIR/graphene_KuboBastinSea.dat"  >/dev/null \
    && echo "PASS: Kubo-Bastin Sea byte-identical to golden" \
    || { echo "FAIL: Kubo-Bastin Sea differs from golden"; exit 1; }
  diff -q "$(ls KuboBastinSurf_*.dat)" "$GDIR/graphene_KuboBastinSurf.dat" >/dev/null \
    && echo "PASS: Kubo-Bastin Surf byte-identical to golden" \
    || { echo "FAIL: Kubo-Bastin Surf differs from golden"; exit 1; }
)
echo "PASS: Sea/Surf reconstruction drivers reproduce their goldens byte-for-byte"
