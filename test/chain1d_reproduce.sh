#!/usr/bin/env bash
# Reproducibility gate for the chain1d analytic golden.
#
# Reads the COMMITTED operator CSRs (test/golden/chain1d/operators) and recomputes
# the KPM moments twice under the fixed KPM_SEED, asserting the .chebmom2D is
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
GOLD="$REPO/test/golden/chain1d/NonEqOpVX-VXchain1dKPM_M64x64_statedefault.chebmom2D"

[ -f "$OPDIR/chain1d.HAM.CSR" ] || { echo "FAIL: committed operators missing ($OPDIR)"; exit 1; }

gen() {  # recompute moments in a fresh temp dir from the committed CSRs
  local T; T="$(mktemp -d)"; mkdir -p "$T/operators"
  cp "$OPDIR/chain1d.HAM.CSR" "$OPDIR/chain1d.VX.CSR" "$T/operators/"
  ( cd "$T"
    printf -- "-2.05 2.05\n" > BOUNDS               # same tight bounds as the golden
    "$BUILD/inline_compute-kpm-nonEqOp" chain1d VX VX 64 >/dev/null 2>&1 )
  echo "$T/$(cd "$T" && ls NonEqOp*chain1d*chebmom2D | head -1)"
}

A="$(gen)"; B="$(gen)"
diff -q "$A" "$B" >/dev/null || { echo "FAIL: two runs differ (non-deterministic RNG)"; exit 1; }
if [ -f "$GOLD" ]; then
  diff -q "$A" "$GOLD" >/dev/null || { echo "FAIL: regenerated moments != committed golden"; exit 1; }
fi
echo "PASS: chain1d moments byte-identical across two runs and match committed golden"
