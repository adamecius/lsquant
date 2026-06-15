#!/usr/bin/env bash
# Phase 2 -- bounds ingested from the .desc sidecar (authority rule).
#
# Proves lsquant takes the spectral window (a,b) from the operator descriptor instead of a
# hack. The clean-chain DOS is reconstructed with NO BOUNDS file present, so without the
# descriptor the solver would fall back to the loose Gershgorin estimate. With the committed
# chain1d.HAM.desc (Lanczos bounds [-2,2]) the rescaling is exact and the reconstructed DOS
# matches the closed form rho(E) = 1/(pi sqrt(4 - E^2)) at the interior energies.
#
# Exact trace (local state over the full basis) -> deterministic. Reference is the oracle
# closed form, not another run.
#
# Usage: descriptor_bounds.sh <BUILD_DIR> [TOL]
set -euo pipefail
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BUILD="$(cd "${1:-$REPO/build}" && pwd)"
TOL="${2:-0.01}"
DOS="$BUILD/inline_kpm-DOS-standalone"
OPD="$REPO/test/golden/chain1d_ballistic/operators"
DESC="$OPD/chain1d.HAM.desc"

[ -x "$DOS" ]  || { echo "FAIL: DOS driver missing ($DOS)"; exit 1; }
[ -f "$DESC" ] || { echo "FAIL: committed chain1d.HAM.desc missing ($DESC)"; exit 1; }
[ -f "$OPD/chain1d.HAM.CSR" ] || { echo "FAIL: committed chain operators missing"; exit 1; }
# sanity: the descriptor must actually carry bounds (this is what we are testing)
grep -q 'spectral_min' "$DESC" || { echo "FAIL: $DESC has no spectral bounds"; exit 1; }

N=$(head -1 "$OPD/chain1d.HAM.CSR" | awk '{print $1}')
T="$(mktemp -d)"; trap 'rm -rf "$T"' EXIT
mkdir -p "$T/operators"
cp "$OPD/chain1d.HAM.CSR" "$T/operators/"
cp "$DESC" "$T/operators/"                       # the descriptor travels with the operator
( cd "$T"
  # *** intentionally NO BOUNDS file *** -> bounds must come from chain1d.HAM.desc
  { echo local; echo "$N"; seq 0 $((N-1)); } > exact
  "$DOS" chain1d 1 64 exact >dos.log 2>&1 || { echo "FAIL: DOS run crashed"; tail -8 dos.log; exit 1; }
  grep -q 'from descriptor' dos.log || { echo "FAIL: solver did not take bounds from the descriptor"; grep -i 'SpectralBounds' dos.log; exit 1; }
  echo "  $(grep -i 'from descriptor' dos.log | head -1 | sed 's/^[[:space:]]*//')"
  D=$(ls mean*JACKSON.dat | head -1)
  awk -v n="$N" 'tolower($2) !~ /nan|inf/ {printf "%.10g %.10g\n", $1, $2/n}' "$D" > dos.dat
  "$BUILD/dos_reconstructed_test" dos.dat "$TOL"
)
