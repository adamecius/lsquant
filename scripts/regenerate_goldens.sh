#!/usr/bin/env bash
#
# Regenerate LSQUANT golden masters, deterministically.
#
# Pipeline (mirrors the reference flow, corrected to the real tool/driver CLIs):
#   external wannier2sparse  ->  inline_compute-kpm-nonEqOp  ->  inline_kuboBastinFromChebmom
#
#   L-A (inputs):  test/golden/operators/graphene_golden.<OP>.CSR   (from wannier2sparse)
#   L-B (outputs): test/golden/NonEqOp*...chebmom2D                 (KPM moments)
#                  test/golden/KuboBastin_*...dat                   (Kubo-Bastin conductivity)
#
# Determinism: random-phase vectors are pinned by KPM_SEED (honoured by
# qstates::generator). Two runs with the same KPM_SEED are bitwise identical.
#
# Usage:  scripts/regenerate_goldens.sh [BUILD_DIR]      (BUILD_DIR default: build)
set -euo pipefail

export KPM_SEED="${KPM_SEED:-12345}"            # MANDATORY: pins random-phase vectors

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO"

W2S="$(bash scripts/fetch_wannier2sparse.sh)"   # cache-aware; offline-safe on cache hit
BUILD="$(cd "${1:-build}" && pwd)"              # absolute path to the LSQUANT build dir

SEED=graphene_golden
OUT="$REPO/test/golden"
NX=8; NY=8; NZ=1; M=64

[ -f "$OUT/${SEED}.xyz" ] || { echo "ERROR: committed seed $OUT/${SEED}.xyz missing"; exit 1; }

cd "$OUT"
rm -rf operators; mkdir -p operators
# L-A — Hamiltonian + operators (deterministic; pure model expansion)
"$W2S" "$SEED" "$NX" "$NY" "$NZ" all -o operators
# L-B — KPM moments (deterministic under KPM_SEED) then Kubo-Bastin conductivity
"$BUILD/inline_compute-kpm-nonEqOp" "$SEED" VX VX "$M"
MOM="$(ls NonEqOp*"${SEED}"*KPM_M"${M}"x"${M}"*.chebmom2D | head -1)"
"$BUILD/inline_kuboBastinFromChebmom" "$MOM" "$M"

# --- chain1d ANALYTIC goldens: DOS moments (#1) + Kubo-Greenwood matrix (#3) ---
# EXACT spectral bounds = the true band [-2,2] are MANDATORY. For the DOS the rescaled
# band must fill [-1,1] exactly so the moments collapse to delta_{m,0} (a pad leaves the
# divergent band-edge weight out and ruins it); the KG (v^2-weighted, edge weight -> 0)
# is happy at [-2,2] too. The driver reads bounds from a BOUNDS file (else [-100,100]).
# NX=2048 keeps the single random vector's stochastic residual (~1/sqrt(N)) small while
# staying tiny to commit.
CHAIN_NX=2048
COUT="$OUT/chain1d"
[ -f "$COUT/chain1d.xyz" ] || { echo "ERROR: committed chain1d seed $COUT/chain1d.xyz missing"; exit 1; }
( cd "$COUT"
  rm -rf operators; mkdir -p operators
  "$W2S" chain1d "$CHAIN_NX" 1 1 VX -o operators   # emits HAM + VX only (minimal fixture)
  printf -- "-2 2\n" > BOUNDS                       # EXACT band
  # EXACT TRACE: a local-state set over the full basis e_0..e_{D-1} makes
  # SpectralMoments compute Sum_i <e_i|T_m|e_i>/D = Tr[T_m]/D exactly (no random
  # vector, deterministic, machine precision). Uses LSQUANT's own kernel; no source change.
  DIM=$(head -1 operators/chain1d.HAM.CSR | awk '{print $1}')
  { echo local; echo "$DIM"; seq 0 $((DIM-1)); } > exact
  "$BUILD/inline_compute-kpm-nonEqOp"  chain1d VX VX "$M" exact   # #3 KG moments (exact) -> .chebmom2D
  "$BUILD/inline_compute-kpm-spectralOp" chain1d 1 "$M" exact     # #1 DOS moments (exact) -> .chebmom1D
  rm -f exact )
CMOM="$(cd "$COUT" && ls NonEqOp*chain1d*KPM_M"${M}"x"${M}"_stateexact.chebmom2D | head -1)"
DMOM="$(cd "$COUT" && ls SpectralOp1chain1d*KPM_M"${M}"_stateexact.chebmom1D | head -1)"
# (The analytic reference is the closed form embedded in the C++ comparators -- the
#  proven oracle formula -- so there is nothing to (re)generate here and no Python dep.)

# #2 reconstructed-DOS golden (averaged over N_r runs; vs closed form in the comparator).
bash "$REPO/test/dos_reconstructed.sh" "$BUILD"

cd "$REPO"
{
  echo "graphene: KPM_SEED=$KPM_SEED w2s_sha=26cab4a lsquant_head=$(git rev-parse HEAD) Nx=$NX Ny=$NY Nz=$NZ M=$M ops=VX,VX"
  echo "chain1d:  KPM_SEED=$KPM_SEED w2s_sha=26cab4a lsquant_head=$(git rev-parse HEAD) Nx=$CHAIN_NX Ny=1 Nz=1 M=$M bounds=[-2,2] KG=$CMOM(ops=VX,VX) DOS=$DMOM(op=1) ref=oracle:lsquant_chain_reference.py"
} > GOLDEN_PROVENANCE.txt
echo "Goldens regenerated: graphene ($MOM, KuboBastin), chain1d KG ($CMOM) + DOS ($DMOM) vs oracle."
