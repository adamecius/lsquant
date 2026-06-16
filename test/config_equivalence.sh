#!/usr/bin/env bash
# Phase 4 -- JSON run-config equivalence + inspect.
#
# (1) A run.json drives inline_compute-kpm-nonEqOp to produce moments BYTE-IDENTICAL to the
#     legacy positional-argv invocation (same KPM_SEED) -- the config path changes how a run is
#     described, not what it computes.
# (2) `lsquant inspect run.json` prints the run and its resolved numerical provenance
#     (recon alpha, bounds source, hbar/units) -- greppable, no need to decode argv/filenames.
#
# Usage: config_equivalence.sh <BUILD_DIR>
set -euo pipefail
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BUILD="$(cd "${1:-$REPO/build}" && pwd)"
DRV="$BUILD/inline_compute-kpm-nonEqOp"
INSPECT="$BUILD/lsquant_inspect"
OPD="$REPO/test/golden/operators"

[ -x "$DRV" ]     || { echo "FAIL: nonEqOp driver missing"; exit 1; }
[ -x "$INSPECT" ] || { echo "FAIL: lsquant_inspect missing"; exit 1; }
[ -f "$OPD/graphene_golden.HAM.CSR" ] || { echo "FAIL: graphene operators missing"; exit 1; }

export KPM_SEED=12345
T="$(mktemp -d)"; trap 'rm -rf "$T"' EXIT
mkdir -p "$T/operators"
cp "$OPD/graphene_golden.HAM.CSR" "$OPD/graphene_golden.VX.CSR" "$T/operators/"
cp "$OPD/graphene_golden.HAM.desc" "$T/operators/" 2>/dev/null || true   # for the inspect bounds line
cat > "$T/run.json" <<'JSON'
{
  "label":          "graphene_golden",
  "operator_left":  "VX",
  "operator_right": "VX",
  "num_moments":    64,
  "kernel":         "jackson",
  "state":          "default"
}
JSON

( cd "$T"
  # (1) equivalence: argv vs --config, byte-identical
  "$DRV" graphene_golden VX VX 64 >/dev/null 2>&1
  mv NonEqOpVX-VXgraphene_golden*.chebmom2D argv.chebmom2D
  "$DRV" --config run.json >/dev/null 2>&1
  mv NonEqOpVX-VXgraphene_golden*.chebmom2D config.chebmom2D
  if diff -q argv.chebmom2D config.chebmom2D >/dev/null; then
      echo "PASS(equivalence): run.json moments are byte-identical to the argv path"
  else
      echo "FAIL(equivalence): run.json moments differ from the argv path"; exit 1
  fi

  # (1b) DOS (spectral): inline argv vs `lsquant compute --config` (mode spectral), byte-identical
  LSQ="$BUILD/lsquant"
  "$BUILD/inline_compute-kpm-spectralOp" graphene_golden VX 64 >/dev/null 2>&1
  mv SpectralOpVXgraphene_golden*.chebmom1D dos_argv.chebmom1D
  printf '{ "mode":"spectral", "label":"graphene_golden", "operator":"VX", "num_moments":64 }\n' > dos.json
  "$LSQ" compute --config dos.json >/dev/null 2>&1
  mv SpectralOpVXgraphene_golden*.chebmom1D dos_config.chebmom1D
  if diff -q dos_argv.chebmom1D dos_config.chebmom1D >/dev/null; then
      echo "PASS(equivalence/DOS): spectral run.json moments byte-identical to the argv path"
  else
      echo "FAIL(equivalence/DOS): spectral run.json moments differ from the argv path"; exit 1
  fi

  # (1c) MSD: inline argv vs `lsquant compute --config` (mode msd), byte-identical
  "$BUILD/inline_compute-kpm-MeanSquareDisplacement" graphene_golden VX 64 8 10 >/dev/null 2>&1
  mv Correlation*graphene_golden*.chebmomTD msd_argv.chebmomTD
  printf '{ "mode":"msd", "label":"graphene_golden", "operator":"VX", "num_moments":64, "num_times":8, "tmax":10 }\n' > msd.json
  "$LSQ" compute --config msd.json >/dev/null 2>&1
  mv Correlation*graphene_golden*.chebmomTD msd_config.chebmomTD
  if diff -q msd_argv.chebmomTD msd_config.chebmomTD >/dev/null; then
      echo "PASS(equivalence/MSD): msd run.json moments byte-identical to the argv path"
  else
      echo "FAIL(equivalence/MSD): msd run.json moments differ from the argv path"; exit 1
  fi

  # (2) inspect prints the run + numerical provenance
  rep="$("$INSPECT" run.json)"
  echo "$rep" | sed 's/^/    /'
  echo "$rep" | grep -q "graphene_golden"                  || { echo "FAIL(inspect): label missing"; exit 1; }
  echo "$rep" | grep -q "M = 64"                           || { echo "FAIL(inspect): M missing"; exit 1; }
  echo "$rep" | grep -qiE "recon alpha"                    || { echo "FAIL(inspect): alpha missing"; exit 1; }
  echo "$rep" | grep -qiE "hbar"                           || { echo "FAIL(inspect): hbar missing"; exit 1; }
  echo "$rep" | grep -qiE "bounds from"                    || { echo "FAIL(inspect): bounds source missing"; exit 1; }
  echo "PASS(inspect): report includes run parameters and numerical provenance"
)
echo "PASS: JSON config is equivalent to argv and inspect reports provenance"
