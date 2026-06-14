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

cd "$REPO"
echo "KPM_SEED=$KPM_SEED  w2s_sha=511f4db  linqt_head=$(git rev-parse HEAD)  Nx=$NX Ny=$NY Nz=$NZ M=$M  ops=VX,VX" \
    > GOLDEN_PROVENANCE.txt
echo "Goldens regenerated under test/golden (operators/*.CSR, $MOM, KuboBastin_*.dat)."
