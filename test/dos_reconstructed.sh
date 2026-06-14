#!/usr/bin/env bash
# Regenerate the reconstructed-DOS golden (#2) for the clean 1D chain by EXACT TRACE.
#
# The DOS Chebyshev moments are delta_{m,0} (chain1d_dos), so the KPM-reconstructed DOS is
# the Jackson-smoothed closed form rho(E)=1/(pi sqrt(4-E^2)) (oracle dos()). A local-state
# set over the full basis makes inline_kpm-DOS-standalone evaluate the EXACT trace in one
# run -- no random vector, no N_r averaging, deterministic. EXACT bounds [-2,2] are
# mandatory (the rescaled band must fill [-1,1]); the driver applies the Jackson kernel.
#
# Writes test/golden/chain1d/DOS_reconstructed.dat (columns: E   rho_per_site).
# Usage: dos_reconstructed.sh [BUILD_DIR]
set -euo pipefail
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BUILD="$(cd "${1:-$REPO/build}" && pwd)"
DOS="$BUILD/inline_kpm-DOS-standalone"
OPDIR="$REPO/test/golden/chain1d/operators"
M=64
OUT="$REPO/test/golden/chain1d/DOS_reconstructed.dat"

[ -f "$OPDIR/chain1d.HAM.CSR" ] || { echo "ERROR: committed chain1d HAM missing"; exit 1; }
N=$(head -1 "$OPDIR/chain1d.HAM.CSR" | awk '{print $1}')   # SystemSize (chain length)

T="$(mktemp -d)"; mkdir -p "$T/operators"
cp "$OPDIR/chain1d.HAM.CSR" "$T/operators/"
( cd "$T"
  printf -- "-2 2\n" > BOUNDS
  { echo local; echo "$N"; seq 0 $((N-1)); } > exact     # exact trace over the full basis
  "$DOS" chain1d 1 "$M" exact >/dev/null 2>&1
  D=$(ls mean*JACKSON.dat | head -1)
  awk -v n="$N" 'tolower($2) !~ /nan|inf/ {printf "%.10g %.10g\n", $1, $2/n}' "$D" > "$OUT"   # per-site, drop band-edge nan/inf
)
rm -rf "$T"
echo "wrote $OUT  (M=$M, EXACT trace, bounds [-2,2], per-site DOS)"
