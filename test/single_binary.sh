#!/usr/bin/env bash
# Phase 7 -- the single `lsquant` dispatcher binary.
#
# Asserts the dispatcher subcommands match the legacy path:
#   (1) `lsquant reconstruct ... bastin|greenwood` is BYTE-IDENTICAL to the legacy
#       inline_kuboBastin/Greenwood drivers (same unified routine, one binary);
#   (2) `lsquant inspect run.json` reports the run; and the standalone lsquant_inspect (thin
#       wrapper sharing lsquant::inspect_path) gives the SAME output.
#
# Usage: single_binary.sh <BUILD_DIR>
set -euo pipefail
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BUILD="$(cd "${1:-$REPO/build}" && pwd)"
LSQ="$BUILD/lsquant"
GMOM="$REPO/test/golden/NonEqOpVX-VXgraphene_goldenKPM_M64x64_statedefault.chebmom2D"

[ -x "$LSQ" ] || { echo "FAIL: lsquant dispatcher missing"; exit 1; }
[ -f "$GMOM" ] || { echo "FAIL: graphene moments missing"; exit 1; }

T="$(mktemp -d)"; trap 'rm -rf "$T"' EXIT
cp "$GMOM" "$T/"
( cd "$T"
  M="$(ls *.chebmom2D | head -1)"

  # (1) reconstruct == legacy, byte-for-byte
  "$BUILD/inline_kuboBastinFromChebmom" "$M" 10 >/dev/null 2>&1;  mv KuboBastin_*JACKSON.dat legacy_b.dat
  "$LSQ" reconstruct "$M" bastin 10 >/dev/null 2>&1;             mv KuboBastin_*JACKSON.dat disp_b.dat
  diff -q legacy_b.dat disp_b.dat >/dev/null \
    && echo "PASS(reconstruct/bastin): lsquant == legacy driver byte-for-byte" \
    || { echo "FAIL: lsquant reconstruct bastin differs from legacy"; exit 1; }

  "$BUILD/inline_kuboGreenwoodFromChebmom" "$M" 10 >/dev/null 2>&1; mv KuboGreenWood_*JACKSON.dat legacy_g.dat
  "$LSQ" reconstruct "$M" greenwood 10 >/dev/null 2>&1;            mv KuboGreenWood_*JACKSON.dat disp_g.dat
  diff -q legacy_g.dat disp_g.dat >/dev/null \
    && echo "PASS(reconstruct/greenwood): lsquant == legacy driver byte-for-byte" \
    || { echo "FAIL: lsquant reconstruct greenwood differs from legacy"; exit 1; }

  # (2) inspect via dispatcher == standalone lsquant_inspect
  printf '{ "label":"graphene_golden", "operator_left":"VX", "operator_right":"VX", "num_moments":64 }\n' > run.json
  "$LSQ" inspect run.json > insp_dispatch.txt 2>&1
  grep -q "M = 64" insp_dispatch.txt || { echo "FAIL: lsquant inspect missing run info"; exit 1; }
  echo "PASS(inspect): lsquant inspect reports the run"
  if [ -x "$BUILD/lsquant_inspect" ]; then
    "$BUILD/lsquant_inspect" run.json > insp_standalone.txt 2>&1
    diff -q insp_dispatch.txt insp_standalone.txt >/dev/null \
      && echo "PASS(inspect-share): dispatcher and standalone share lsquant::inspect_path" \
      || { echo "FAIL: dispatcher inspect != standalone inspect"; exit 1; }
  fi

  # (3) compute == legacy driver, byte-for-byte; then an end-to-end compute -> reconstruct.
  OPD="$REPO/test/golden/operators"
  if [ -f "$OPD/graphene_golden.HAM.CSR" ] && [ -f "$OPD/graphene_golden.VX.CSR" ]; then
    export KPM_SEED=12345
    mkdir -p operators; cp "$OPD/graphene_golden.HAM.CSR" "$OPD/graphene_golden.VX.CSR" operators/
    printf '{ "label":"graphene_golden","operator_left":"VX","operator_right":"VX","num_moments":64 }\n' > comp.json
    "$BUILD/inline_compute-kpm-nonEqOp" graphene_golden VX VX 64 >/dev/null 2>&1; mv NonEqOp*chebmom2D legacy.m
    "$LSQ" compute --config comp.json >/dev/null 2>&1;                            mv NonEqOp*chebmom2D disp.m
    diff -q legacy.m disp.m >/dev/null \
      && echo "PASS(compute/noneq): lsquant compute == legacy driver byte-for-byte" \
      || { echo "FAIL: lsquant compute differs from legacy driver"; exit 1; }

    # all three compute modes go through the one dispatcher, each == its legacy inline driver
    "$BUILD/inline_compute-kpm-spectralOp" graphene_golden VX 64 >/dev/null 2>&1; mv SpectralOp*chebmom1D legacy_s.m
    printf '{ "mode":"spectral","label":"graphene_golden","operator":"VX","num_moments":64 }\n' > spec.json
    "$LSQ" compute --config spec.json >/dev/null 2>&1;                            mv SpectralOp*chebmom1D disp_s.m
    diff -q legacy_s.m disp_s.m >/dev/null \
      && echo "PASS(compute/spectral): lsquant compute == legacy spectral driver byte-for-byte" \
      || { echo "FAIL: lsquant compute spectral differs from legacy driver"; exit 1; }

    "$BUILD/inline_compute-kpm-MeanSquareDisplacement" graphene_golden VX 64 8 10 >/dev/null 2>&1; mv Correlation*chebmomTD legacy_m.m
    printf '{ "mode":"msd","label":"graphene_golden","operator":"VX","num_moments":64,"num_times":8,"tmax":10 }\n' > msd.json
    "$LSQ" compute --config msd.json >/dev/null 2>&1;                             mv Correlation*chebmomTD disp_m.m
    diff -q legacy_m.m disp_m.m >/dev/null \
      && echo "PASS(compute/msd): lsquant compute == legacy msd driver byte-for-byte" \
      || { echo "FAIL: lsquant compute msd differs from legacy driver"; exit 1; }
    # end-to-end: the moments lsquant just produced reconstruct cleanly through the same binary
    "$LSQ" reconstruct disp.m bastin 10 >/dev/null 2>&1
    e2e="$(ls KuboBastin_*JACKSON.dat | head -1)"
    nans=$(grep -ciE 'nan|inf' "$e2e" || true)
    [ "$nans" -eq 0 ] \
      && echo "PASS(e2e): lsquant compute -> reconstruct produces a finite curve ($e2e)" \
      || { echo "FAIL: e2e compute->reconstruct has $nans NaN/Inf rows"; exit 1; }
  fi
)
echo "PASS: single lsquant binary subsumes compute + reconstruct + inspect (byte-identical to legacy)"
