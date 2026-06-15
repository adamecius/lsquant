#!/usr/bin/env bash
# Phase 5 (slice 3) -- regression net for the FFT Chebyshev-node reconstruction backend.
#
# The FFT-node path (cos(pi(2i+1/2)/M) grid, M points, rearranged for the Bastin integral) is
# coarse and structurally different from the uniform grid -- the two do NOT agree numerically
# (grid resolution), so the only faithful check is byte-reproduction of a committed golden. This
# path previously had NO test; this pins it, so the unification (and any future change) is safe.
#
# Reconstructs the FFT-grid Kubo-Bastin and Kubo-Greenwood from the committed graphene VX-VX
# moments and byte-compares to the committed .dat goldens.
#
# Usage: fft_grid_regression.sh <BUILD_DIR>
set -euo pipefail
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BUILD="$(cd "${1:-$REPO/build}" && pwd)"
GMOM="$REPO/test/golden/NonEqOpVX-VXgraphene_goldenKPM_M64x64_statedefault.chebmom2D"
GDIR="$REPO/test/golden/fft_grid"

[ -f "$GMOM" ] || { echo "FAIL: committed graphene moments missing"; exit 1; }
for d in "$BUILD/inline_kuboBastinFromChebmom_FFTgrid" "$BUILD/inline_kuboGreenwoodFromChebmom_FFTgrid"; do
  [ -x "$d" ] || { echo "FAIL: driver missing ($d)"; exit 1; }
done

T="$(mktemp -d)"; trap 'rm -rf "$T"' EXIT
cp "$GMOM" "$T/"
( cd "$T"
  M="$(ls *.chebmom2D | head -1)"
  "$BUILD/inline_kuboBastinFromChebmom_FFTgrid"    "$M" 10 >/dev/null 2>&1
  "$BUILD/inline_kuboGreenwoodFromChebmom_FFTgrid" "$M" 10 >/dev/null 2>&1

  diff -q "$(ls KuboBastin_*FFTgrid.dat)"    "$GDIR/graphene_KuboBastin_FFTgrid.dat"    >/dev/null \
    && echo "PASS: FFT-grid Kubo-Bastin byte-identical to golden" \
    || { echo "FAIL: FFT-grid Kubo-Bastin differs from golden"; exit 1; }
  diff -q "$(ls KuboGreenWood_*FFTgrid.dat)" "$GDIR/graphene_KuboGreenwood_FFTgrid.dat" >/dev/null \
    && echo "PASS: FFT-grid Kubo-Greenwood byte-identical to golden" \
    || { echo "FAIL: FFT-grid Kubo-Greenwood differs from golden"; exit 1; }
)
echo "PASS: FFT-node reconstruction backend reproduces its goldens byte-for-byte"
