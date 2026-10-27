#!/usr/bin/env bash
cd "$HOME/lsquant-wrapup-20261024"
export PATH="$HOME/lsq-env/bin:$PATH" LD_LIBRARY_PATH="$HOME/lsq-env/lib" PYTHONPATH="$PWD/candidate-v5/python:$PWD/build-v1-golden/lsquant_core"
export OMP_NUM_THREADS=8 OMP_PROC_BIND=false OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONUNBUFFERED=1 MPLBACKEND=Agg
for material in bi pdse2; do
 python review-figures-20261025/final/stage_9_material_review.py --root "$HOME/lsquant-verify-20261023-1100" --repo "$PWD/candidate-v5" --out "$PWD/review-figures-20261025/final" --material "$material" > "review-figures-20261025/final/$material.log" 2>&1 &
 child=$!; echo "$child" > "review-figures-20261025/final/$material.pid"
 wait "$child"; code=$?; echo "$code" > "review-figures-20261025/final/$material.exit"
 if [ "$code" != 0 ]; then exit "$code"; fi
done
