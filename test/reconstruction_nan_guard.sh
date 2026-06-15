#!/usr/bin/env bash
# Reconstruction NaN/Inf guard for the alpha safeguard (Phase 1).
#
# The Kubo-Bastin reconstruction kernels carry 1/sqrt(1 - x^2) (delta_chebF /
# DgreenR_chebF), which is NaN at the band edges x = +-1. The cumulative Fermi-sea
# integral starts at that edge, so a single leading NaN poisons the ENTIRE curve --
# this is the historical "all-NaN graphene Kubo .dat" artifact. The fix is the
# KPM_ALPHA = 0.95 reconstruction-grid cutoff (chebyshev_moments.hpp): the grid runs
# only on [-alpha, alpha], never touching the singularity.
#
# This test reconstructs the committed graphene VX-VX moments through the Kubo-Bastin
# driver and asserts the output is finite everywhere. Before the alpha fix this .dat
# was 100% NaN; with alpha it must be 0% NaN/Inf.
#
# Usage: reconstruction_nan_guard.sh <BUILD_DIR>
set -euo pipefail
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BUILD="$(cd "${1:-$REPO/build}" && pwd)"
DRIVER="$BUILD/inline_kuboBastinFromChebmom"
GMOM="$REPO/test/golden/NonEqOpVX-VXgraphene_goldenKPM_M64x64_statedefault.chebmom2D"

[ -x "$DRIVER" ] || { echo "FAIL: Kubo-Bastin driver missing ($DRIVER)"; exit 1; }
[ -f "$GMOM" ]   || { echo "FAIL: committed graphene Kubo moments missing ($GMOM)"; exit 1; }

T="$(mktemp -d)"; trap 'rm -rf "$T"' EXIT
cp "$GMOM" "$T/"
( cd "$T"
  MOM="$(ls *.chebmom2D | head -1)"
  "$DRIVER" "$MOM" 10 >driver.log 2>&1 \
    || { echo "FAIL: Kubo-Bastin reconstruction crashed"; cat driver.log; exit 1; }
  DAT="$(ls KuboBastin*.dat 2>/dev/null | head -1)"
  [ -n "$DAT" ] || { echo "FAIL: no reconstruction .dat produced"; cat driver.log; exit 1; }

  TOTAL=$(wc -l < "$DAT")
  NANS=$(grep -ciE 'nan|inf' "$DAT" || true)
  RANGE=$(awk 'NR==1{f=$1} END{printf "[%g, %g]", f, $1}' "$DAT")
  echo "reconstructed $TOTAL points over E in $RANGE ; NaN/Inf lines = $NANS"
  [ "$NANS" -eq 0 ] \
    || { echo "FAIL: $NANS/$TOTAL reconstruction points are NaN/Inf (alpha safeguard not active?)"; exit 1; }
  echo "PASS: Kubo-Bastin reconstruction is finite everywhere (alpha=0.95 grid keeps off the 1/sqrt(1-x^2) edge)"
)
