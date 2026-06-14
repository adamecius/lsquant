#!/usr/bin/env bash
set -euo pipefail
SHA=511f4dbf83dfb031ac7d91dd7fe39c55b2026f5d
CACHE="${W2S_CACHE:-$PWD/.cache/wannier2sparse}"
BIN="$CACHE/build/wannier2sparse"
[ -x "$BIN" ] && { echo "$BIN"; exit 0; }                      # cache hit, offline-safe
git clone https://github.com/adamecius/wannier2sparse.git "$CACHE/src" 2>/dev/null || true
git -C "$CACHE/src" fetch --depth 1 origin "$SHA"
git -C "$CACHE/src" checkout -q "$SHA"
cmake -S "$CACHE/src" -B "$CACHE/build" -DW2SP_BUILD_TESTS=OFF ${EIGEN_DIR:+-DEIGEN_DIR=$EIGEN_DIR}
cmake --build "$CACHE/build" -j
echo "$BIN"
