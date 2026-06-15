#!/usr/bin/env bash
# Smoke test for examples/01_dos_1d_chain (matplotlib-free).
#
# Builds a small clean chain with the tutorial's make_chain.py, computes the DOS moments with
# the real binary, and asserts mu_m = delta_{m,0} via the existing chain1d_dos comparator. Uses
# M < N (M=32, N=64) so there is no aliasing revival and the full moment range is delta_{m,0}.
# This pins that the tutorial pipeline runs end-to-end and produces the analytic chain moments.
#
# Usage: example_dos_smoke.sh <BUILD_DIR>
set -euo pipefail
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BUILD="$(cd "${1:-$REPO/build}" && pwd)"
EX="$REPO/examples/01_dos_1d_chain"
PY="${PYTHON:-python3}"

[ -x "$BUILD/inline_compute-kpm-spectralOp" ] || { echo "FAIL: spectralOp driver missing"; exit 1; }
[ -x "$BUILD/chain1d_dos_test" ]             || { echo "FAIL: chain1d_dos_test missing"; exit 1; }
[ -f "$EX/make_chain.py" ]                   || { echo "FAIL: example make_chain.py missing"; exit 1; }
command -v "$PY" >/dev/null 2>&1             || { echo "SKIP: no python3"; exit 0; }

T="$(mktemp -d)"; trap 'rm -rf "$T"' EXIT
cp "$EX/make_chain.py" "$T/"
( cd "$T"
  "$PY" make_chain.py 64 >/dev/null            # writes operators/, BOUNDS, exact_chain1d_N64
  "$BUILD/inline_compute-kpm-spectralOp" chain1d_N64 1 32 exact_chain1d_N64 >/dev/null 2>&1
  MOM="$(ls SpectralOp1*chain1d_N64*.chebmom1D | head -1)"
  echo "  produced $MOM"
  # M=32 < N=64 -> no aliasing -> moments must equal delta_{m,0} to ~machine precision
  "$BUILD/chain1d_dos_test" "$MOM" 1e-6
)
echo "PASS: tutorial DOS pipeline runs and yields mu_m = delta_{m,0}"
