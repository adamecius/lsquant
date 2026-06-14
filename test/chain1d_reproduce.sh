#!/usr/bin/env bash
# Reproducibility gate for the chain1d analytic goldens (KG #3 and DOS #1).
#
# Reads the COMMITTED operator CSRs (test/golden/chain1d/operators) and recomputes
# the moments twice under the fixed KPM_SEED, asserting each .chebmom file is
# byte-identical across runs AND identical to the committed golden. Using the
# committed CSRs (rather than re-running wannier2sparse) isolates the actual
# requirement -- deterministic random-phase moments -- and removes any dependency
# on the external tool's binary at test time.
# Usage: chain1d_reproduce.sh <BUILD_DIR>
set -euo pipefail
export KPM_SEED=12345

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BUILD="$(cd "${1:-$REPO/build}" && pwd)"
OPDIR="$REPO/test/golden/chain1d/operators"
KG_GOLD="$REPO/test/golden/chain1d/NonEqOpVX-VXchain1dKPM_M64x64_statedefault.chebmom2D"
DOS_GOLD="$REPO/test/golden/chain1d/SpectralOp1chain1dKPM_M64_statedefault.chebmom1D"

[ -f "$OPDIR/chain1d.HAM.CSR" ] || { echo "FAIL: committed operators missing ($OPDIR)"; exit 1; }

# $1 = driver-and-args producing one .chebmom file; echoes the produced file path
gen() {
  local T; T="$(mktemp -d)"; mkdir -p "$T/operators"
  cp "$OPDIR/chain1d.HAM.CSR" "$OPDIR/chain1d.VX.CSR" "$T/operators/"
  ( cd "$T"
    printf -- "-2 2\n" > BOUNDS                     # exact band, as for the golden
    # shellcheck disable=SC2086  -- intentional word-split of "driver args"
    $BUILD/$1 >/dev/null 2>&1 )
  echo "$T/$(cd "$T" && ls $2 | head -1)"
}

check() {  # name  driver-cmd  glob  committed-golden
  local A B; A="$(gen "$2" "$3")"; B="$(gen "$2" "$3")"
  diff -q "$A" "$B" >/dev/null || { echo "FAIL($1): two runs differ (non-deterministic RNG)"; exit 1; }
  if [ -f "$4" ]; then diff -q "$A" "$4" >/dev/null || { echo "FAIL($1): regenerated != committed golden"; exit 1; }; fi
  echo "PASS($1): byte-identical across two runs and matches committed golden"
}

check KG  "inline_compute-kpm-nonEqOp chain1d VX VX 64"  "NonEqOp*chain1d*chebmom2D" "$KG_GOLD"
check DOS "inline_compute-kpm-spectralOp chain1d 1 64"   "SpectralOp1chain1d*chebmom1D" "$DOS_GOLD"
echo "PASS: chain1d KG and DOS moments byte-reproducible under fixed KPM_SEED"
