#!/usr/bin/env bash
# Regenerate the reconstructed-DOS golden (#2) for the clean 1D chain.
#
# The DOS Chebyshev moments are delta_{m,0} (chain1d_dos test), so the KPM-reconstructed
# DOS converges to the closed form rho(E) = 1/(pi sqrt(4-E^2)) (oracle: lsquant_chain_reference.dos).
# inline_kpm-DOS-standalone uses ONE random vector per run (N_r=1); to push the stochastic
# floor below 1% at the interior we average N_r independent runs (distinct KPM_SEED) --
# legitimate N_r averaging with the existing driver, no source change. EXACT bounds [-2,2]
# are mandatory (the rescaled band must fill [-1,1]); Jackson kernel is applied by the driver.
#
# Writes test/golden/chain1d/DOS_reconstructed.dat  (columns: E   rho_per_site).
# Usage: dos_reconstructed.sh [BUILD_DIR]
set -euo pipefail
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BUILD="$(cd "${1:-$REPO/build}" && pwd)"
DOS="$BUILD/inline_kpm-DOS-standalone"
OPDIR="$REPO/test/golden/chain1d/operators"
M=64; NR=128
OUT="$REPO/test/golden/chain1d/DOS_reconstructed.dat"

[ -f "$OPDIR/chain1d.HAM.CSR" ] || { echo "ERROR: committed chain1d HAM missing"; exit 1; }
N=$(head -1 "$OPDIR/chain1d.HAM.CSR" | awk '{print $1}')   # SystemSize (chain length)

T="$(mktemp -d)"; mkdir -p "$T/operators"
cp "$OPDIR/chain1d.HAM.CSR" "$T/operators/"
( cd "$T"
  printf -- "-2 2\n" > BOUNDS                              # EXACT band
  for s in $(seq 1 "$NR"); do
    KPM_SEED="$s" "$DOS" chain1d 1 "$M" >/dev/null 2>&1
    mv mean*JACKSON.dat "seed_$s.dat"
  done
  # average column 2 across the NR runs, normalize to per-site (divide by N)
  awk -v n="$N" '{e[FNR]=$1; s[FNR]+=$2; c[FNR]++} END{for(i=1;i<=FNR;i++) printf "%.10g %.10g\n", e[i], s[i]/(c[i]*n)}' seed_*.dat > "$OUT"
)
rm -rf "$T"
echo "wrote $OUT  (M=$M, N_r=$NR, exact bounds [-2,2], per-site DOS)"
