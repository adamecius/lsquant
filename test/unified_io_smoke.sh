#!/usr/bin/env bash
# Unified I/O (unified_io.md) -- the entry/output story end to end:
#   (1) system FINGERPRINT is stable + content-sensitive;
#   (2) `lsquant run` validates with intuitive hints (Hall needs VY) and aborts;
#   (3) `lsquant run` computes several observables for one system into ONE results.json;
#   (4) `lsquant merge` combines same-system results and REFUSES different systems.
#
# "Goldens first, always": the gate for the unified-I/O layer.
# Usage: unified_io_smoke.sh <BUILD_DIR>
set -uo pipefail
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BUILD="$(cd "${1:-$REPO/build}" && pwd)"
LSQ="$BUILD/lsquant"
OPD="$REPO/test/golden/operators"

[ -x "$LSQ" ] || { echo "FAIL: lsquant missing"; exit 1; }
[ -f "$OPD/graphene_golden.HAM.CSR" ] || { echo "FAIL: operators missing"; exit 1; }

export KPM_SEED=12345
T="$(mktemp -d)"; trap 'rm -rf "$T"' EXIT
mkdir -p "$T/operators"
cp "$OPD/graphene_golden.HAM.CSR" "$OPD/graphene_golden.VX.CSR" "$T/operators/"

( cd "$T"
  # (1) fingerprint: stable across calls, changes when an operator changes
  f1="$("$LSQ" fingerprint graphene_golden operators)"
  f2="$("$LSQ" fingerprint graphene_golden operators)"
  [ "$f1" = "$f2" ] && [ ${#f1} -eq 64 ] \
    && echo "PASS(fingerprint): stable 64-hex digest" \
    || { echo "FAIL: fingerprint unstable/wrong length ($f1)"; exit 1; }
  cp operators/graphene_golden.HAM.CSR operators/graphene_golden.VX.CSR  # perturb VX
  f3="$("$LSQ" fingerprint graphene_golden operators)"
  [ "$f3" != "$f1" ] \
    && echo "PASS(fingerprint): content-sensitive (perturbing an operator changes it)" \
    || { echo "FAIL: fingerprint not content-sensitive"; exit 1; }
  cp "$OPD/graphene_golden.VX.CSR" operators/   # restore

  # (2) validation hint + abort: Hall asked, VY absent
  cat > bad.json <<'JSON'
{ "system": {"label":"graphene_golden","operators_dir":"operators"},
  "numerics": {"num_moments":64,"broadening_meV":10,"grid_points":20},
  "observables": [ {"kind":"conductivity","component":"xy"} ] }
JSON
  if "$LSQ" run --config bad.json 2>err.txt; then echo "FAIL: run with missing VY should abort"; exit 1; fi
  grep -q "VY is missing" err.txt \
    && echo "PASS(validate): missing-VY gives an intuitive hint and aborts" \
    || { echo "FAIL: no intuitive hint for missing VY"; cat err.txt; exit 1; }

  # (3) unified run: DOS + sigma_xx -> one results.json
  cat > good.json <<'JSON'
{ "system": {"label":"graphene_golden","operators_dir":"operators"},
  "numerics": {"num_moments":64,"broadening_meV":10,"grid_points":25},
  "observables": [ {"kind":"dos"}, {"kind":"conductivity","component":"xx"} ],
  "output": {"results":"g.results.json"} }
JSON
  "$LSQ" run --config good.json >/dev/null 2>&1 || { echo "FAIL: unified run failed"; exit 1; }
  grep -q '"dos"' g.results.json && grep -q '"sigma_xx"' g.results.json \
    && grep -q '"fingerprint"' g.results.json && grep -q '"timing_s"' g.results.json \
    && echo "PASS(run): one results.json holds DOS + sigma_xx + timing + fingerprint" \
    || { echo "FAIL: unified results.json missing expected content"; exit 1; }

  # (4a) merge same system: a DOS-only run + a sigma-only run -> combined
  cat > dos.json <<'JSON'
{ "system": {"label":"graphene_golden","operators_dir":"operators"},
  "numerics": {"num_moments":64,"broadening_meV":10,"grid_points":25},
  "observables": [ {"kind":"dos"} ], "output": {"results":"a.results.json"} }
JSON
  cat > sig.json <<'JSON'
{ "system": {"label":"graphene_golden","operators_dir":"operators"},
  "numerics": {"num_moments":64,"broadening_meV":10,"grid_points":25},
  "observables": [ {"kind":"conductivity","component":"xx"} ], "output": {"results":"b.results.json"} }
JSON
  "$LSQ" run --config dos.json >/dev/null 2>&1
  "$LSQ" run --config sig.json >/dev/null 2>&1
  "$LSQ" merge a.results.json b.results.json -o merged.json >/dev/null 2>&1 || { echo "FAIL: merge of same system failed"; exit 1; }
  grep -q '"dos"' merged.json && grep -q '"sigma_xx"' merged.json \
    && echo "PASS(merge): two same-system runs merge into one file" \
    || { echo "FAIL: merged file missing an observable"; exit 1; }

  # (4b) merge refuses different systems (forge a different fingerprint)
  sed 's/"fingerprint": "[0-9a-f]*"/"fingerprint": "deadbeef"/' b.results.json > other.json
  if "$LSQ" merge a.results.json other.json -o no.json 2>merr.txt; then
    echo "FAIL: merge should refuse different systems"; exit 1
  fi
  grep -q "different systems" merr.txt \
    && echo "PASS(merge): refuses to merge different systems (fingerprint mismatch)" \
    || { echo "FAIL: merge mismatch not reported clearly"; cat merr.txt; exit 1; }
)
echo "PASS: unified I/O -- fingerprint, validated multi-observable run, fingerprint-checked merge"
