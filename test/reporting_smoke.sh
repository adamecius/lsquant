#!/usr/bin/env bash
# Reporting layer (reporting_system.md) -- two guarantees:
#
#   (1) NUMERICS UNPERTURBED: `lsquant compute` produces byte-identical moments with
#       the reporting layer silent (default), verbose (-v), and quiet (-q). Logging,
#       timing, and resource sampling must never change the computed output.
#   (2) FACADE BEHAVES: the run summary (peak RSS / timing) appears at the default
#       level; --log-level error silences it; -v adds the resource backend line.
#
# "Goldens first, always": this is the smoke gate for the reporting layer.
#
# Usage: reporting_smoke.sh <BUILD_DIR>
set -euo pipefail
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BUILD="$(cd "${1:-$REPO/build}" && pwd)"
LSQ="$BUILD/lsquant"
OPD="$REPO/test/golden/operators"

[ -x "$LSQ" ] || { echo "FAIL: lsquant dispatcher missing"; exit 1; }
[ -f "$OPD/graphene_golden.HAM.CSR" ] || { echo "FAIL: graphene operators missing"; exit 1; }

export KPM_SEED=12345
T="$(mktemp -d)"; trap 'rm -rf "$T"' EXIT
mkdir -p "$T/operators"
cp "$OPD/graphene_golden.HAM.CSR" "$OPD/graphene_golden.VX.CSR" "$T/operators/"

( cd "$T"
  printf '{ "mode":"spectral","label":"graphene_golden","operator":"VX","num_moments":64 }\n' > spec.json

  # (1) byte-identical moments across log levels
  "$LSQ" compute --config spec.json            >/dev/null 2>&1; mv SpectralOp*chebmom1D quiet_default.m
  "$LSQ" compute --config spec.json -v          >/dev/null 2>&1; mv SpectralOp*chebmom1D verbose.m
  "$LSQ" compute --config spec.json -q          >/dev/null 2>&1; mv SpectralOp*chebmom1D quiet.m
  "$LSQ" compute --config spec.json --log-level trace >/dev/null 2>&1; mv SpectralOp*chebmom1D trace.m
  diff -q quiet_default.m verbose.m >/dev/null && diff -q quiet_default.m quiet.m >/dev/null \
     && diff -q quiet_default.m trace.m >/dev/null \
    && echo "PASS(numerics): moments byte-identical across -v / -q / --log-level trace" \
    || { echo "FAIL: reporting flags changed the computed moments"; exit 1; }

  # (2a) default level: a resource summary line is emitted (stderr)
  "$LSQ" compute --config spec.json 2>summary.txt >/dev/null
  grep -q "peak RSS" summary.txt \
    && echo "PASS(summary): default run emits a resource summary line" \
    || { echo "FAIL: no resource summary at default level"; cat summary.txt; exit 1; }

  # (2b) --log-level error silences info-level summary
  "$LSQ" compute --config spec.json --log-level error 2>silent.txt >/dev/null
  if grep -q "peak RSS" silent.txt; then
    echo "FAIL: --log-level error did not silence the summary"; exit 1
  fi
  echo "PASS(level): --log-level error silences the info summary"

  # (2c) verbose exposes the resource backend
  "$LSQ" compute --config spec.json -v 2>verbose.txt >/dev/null
  grep -q "resource backend=" verbose.txt \
    && echo "PASS(verbose): -v reports the resource backend" \
    || { echo "FAIL: -v did not report the resource backend"; exit 1; }
)
echo "PASS: reporting layer is numerics-transparent and level-controllable"
