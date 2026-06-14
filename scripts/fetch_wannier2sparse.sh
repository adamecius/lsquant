#!/usr/bin/env bash
set -euo pipefail
SHA=a7dac9290c3f4ba5b68719561cd4f2520ef02fa5
CACHE="${W2S_CACHE:-$PWD/.cache/wannier2sparse}"
BIN="$CACHE/build/wannier2sparse"
[ -x "$BIN" ] && { echo "$BIN"; exit 0; }                      # cache hit, offline-safe
git clone https://github.com/adamecius/wannier2sparse.git "$CACHE/src" 2>/dev/null || true
git -C "$CACHE/src" fetch --depth 1 origin "$SHA"
git -C "$CACHE/src" checkout -q "$SHA"
cmake -S "$CACHE/src" -B "$CACHE/build" -DW2SP_BUILD_TESTS=OFF ${EIGEN_DIR:+-DEIGEN_DIR=$EIGEN_DIR}
cmake --build "$CACHE/build" -j
echo "$BIN"
