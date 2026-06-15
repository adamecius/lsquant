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
)
echo "PASS: single lsquant binary subsumes reconstruct + inspect (byte-identical to legacy)"
