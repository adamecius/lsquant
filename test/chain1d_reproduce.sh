#!/usr/bin/env bash
# Determinism + fresh-run physics gate for the chain1d analytic goldens (KG #3, DOS #1).
#
# This test guards TWO distinct properties, kept separate on purpose:
#
#   (1) DETERMINISM (strict, byte-exact): the same binary under a fixed KPM_SEED must
#       produce byte-identical moment files across two independent runs. This is the
#       only thing for which a byte comparison is legitimate -- it is the same binary
#       on the same machine, so any difference means a non-deterministic RNG.
#
#   (2) PHYSICS OF A FRESH RUN (tolerance-based): the freshly recomputed moments must
#       satisfy the analytic reference (B.15 KG matrix / delta_{m,0} DOS) within the
#       stochastic-trace tolerance, checked by the SAME C++ comparators used on the
#       committed goldens. No Python, no data file.
#
# It deliberately does NOT byte-compare the fresh output against the *committed* golden:
# that golden was produced on a different toolchain, and the moments agree only to ~1e-15
# (floating-point reduction order across compiler/Eigen/OpenMP versions). A byte diff
# there fails for a non-physics reason. The committed golden's physics is already
# guarded with tolerance by the chain1d_analytic / chain1d_dos tests.
#
# Usage: chain1d_reproduce.sh <BUILD_DIR> [TOL]
set -euo pipefail
export KPM_SEED=12345

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BUILD="$(cd "${1:-$REPO/build}" && pwd)"
TOL="${2:-0.10}"
OPDIR="$REPO/test/golden/chain1d/operators"

[ -f "$OPDIR/chain1d.HAM.CSR" ] || { echo "FAIL: committed operators missing ($OPDIR)"; exit 1; }

# $1 = driver-and-args producing one .chebmom file; $2 = glob; echoes the produced path
gen() {
  local T; T="$(mktemp -d)"; mkdir -p "$T/operators"
  cp "$OPDIR/chain1d.HAM.CSR" "$OPDIR/chain1d.VX.CSR" "$T/operators/"
  ( cd "$T"
    printf -- "-2 2\n" > BOUNDS                     # exact band, as for the golden
    # shellcheck disable=SC2086  -- intentional word-split of "driver args"
    $BUILD/$1 >/dev/null 2>&1 )
  echo "$T/$(cd "$T" && ls $2 | head -1)"
}

# name  driver-cmd  glob  physics-comparator(reads <mom> <TOL>)
check() {
  local A B; A="$(gen "$2" "$3")"; B="$(gen "$2" "$3")"
  diff -q "$A" "$B" >/dev/null \
    || { echo "FAIL($1): two runs differ (non-deterministic RNG)"; exit 1; }
  "$BUILD/$4" "$A" "$TOL" >/dev/null \
    || { echo "FAIL($1): fresh run fails the analytic reference (TOL=$TOL)"; exit 1; }
  echo "PASS($1): deterministic across two runs AND fresh run matches the analytic reference"
}

check KG  "inline_compute-kpm-nonEqOp chain1d VX VX 64"  "NonEqOp*chain1d*chebmom2D"    chain1d_analytic_test
check DOS "inline_compute-kpm-spectralOp chain1d 1 64"   "SpectralOp1chain1d*chebmom1D" chain1d_dos_test
echo "PASS: chain1d KG and DOS moments are deterministic under KPM_SEED and physically correct on a fresh run"
