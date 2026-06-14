#!/usr/bin/env bash
#
# Regenerate LinQT golden masters, deterministically.
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
BUILD="$(cd "${1:-build}" && pwd)"              # absolute path to the LinQT build dir

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

# --- chain1d ANALYTIC golden: Kubo-Greenwood moments vs thesis Eq. (4.9) -----
# Tight spectral bounds around the true band [-2,2] are MANDATORY: loose bounds
# distort the moment matrix (thesis Fig. 4.2). The driver reads them from a
# BOUNDS file (else defaults to [-100,100]). NX=2048 keeps the single random
# vector's stochastic-trace residual small (~1/sqrt(N)) while staying tiny to commit.
CHAIN_NX=2048
COUT="$OUT/chain1d"
[ -f "$COUT/chain1d.xyz" ] || { echo "ERROR: committed chain1d seed $COUT/chain1d.xyz missing"; exit 1; }
( cd "$COUT"
  rm -rf operators; mkdir -p operators
  "$W2S" chain1d "$CHAIN_NX" 1 1 VX -o operators   # emits HAM + VX only (minimal fixture)
  printf -- "-2.05 2.05\n" > BOUNDS
  "$BUILD/inline_compute-kpm-nonEqOp" chain1d VX VX "$M" )
CMOM="$(cd "$COUT" && ls NonEqOp*chain1d*KPM_M"${M}"x"${M}"*.chebmom2D | head -1)"

cd "$REPO"
{
  echo "graphene: KPM_SEED=$KPM_SEED w2s_sha=511f4db linqt_head=$(git rev-parse HEAD) Nx=$NX Ny=$NY Nz=$NZ M=$M ops=VX,VX"
  echo "chain1d:  KPM_SEED=$KPM_SEED w2s_sha=511f4db linqt_head=$(git rev-parse HEAD) Nx=$CHAIN_NX Ny=1 Nz=1 M=$M ops=VX,VX bounds=[-2.05,2.05] ref=test/golden/chain1d_analytic_M.txt(thesis Eq.4.9)"
} > GOLDEN_PROVENANCE.txt
echo "Goldens regenerated: graphene ($MOM, KuboBastin), chain1d ($CMOM vs Eq. 4.9)."
