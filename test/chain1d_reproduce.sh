#!/usr/bin/env bash
# Exact-trace determinism + physics gate for the chain1d analytic goldens (KG #3, DOS #1),
# plus a minimal random-phase determinism check so that code path stays exercised.
#
# The goldens are now computed by EXACT TRACE: a local-state set over the full basis
# e_0..e_{D-1} makes SpectralMoments evaluate Sum_i <e_i|T_m|e_i>/D = Tr[T_m]/D exactly,
# through LSQUANT's own Chebyshev kernel, with no random vector. So:
#   (1) two runs are byte-identical (deterministic by construction), and
#   (2) the fresh moments satisfy the analytic reference at MACHINE precision (TOL ~ 1e-9),
#       checked by the C++ comparators (no Python). No byte-diff against the committed
#       golden (toolchain FP differs at ~1e-15; physics is guarded with tolerance).
# Separately, (3) the random-phase path is still deterministic under a fixed KPM_SEED.
#
# Usage: chain1d_reproduce.sh <BUILD_DIR> [TOL]
set -euo pipefail
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BUILD="$(cd "${1:-$REPO/build}" && pwd)"
TOL="${2:-1e-9}"
OPDIR="$REPO/test/golden/chain1d/operators"
[ -f "$OPDIR/chain1d.HAM.CSR" ] || { echo "FAIL: committed operators missing ($OPDIR)"; exit 1; }
DIM=$(head -1 "$OPDIR/chain1d.HAM.CSR" | awk '{print $1}')

# exact-trace run: $1 driver-and-args (no state file), $2 glob -> echoes produced path
gen_exact() {
  local T; T="$(mktemp -d)"; mkdir -p "$T/operators"
  cp "$OPDIR/chain1d.HAM.CSR" "$OPDIR/chain1d.VX.CSR" "$T/operators/"
  ( cd "$T"
    printf -- "-2 2\n" > BOUNDS
    { echo local; echo "$DIM"; seq 0 $((DIM-1)); } > exact
    # shellcheck disable=SC2086
    $BUILD/$1 exact >/dev/null 2>&1 )
  echo "$T/$(cd "$T" && ls $2 | head -1)"
}

check_exact() {  # name  driver-and-args  glob  comparator
  local A B; A="$(gen_exact "$2" "$3")"; B="$(gen_exact "$2" "$3")"
  diff -q "$A" "$B" >/dev/null \
    || { echo "FAIL($1): exact-trace runs differ (should be deterministic)"; exit 1; }
  "$BUILD/$4" "$A" "$TOL" >/dev/null \
    || { echo "FAIL($1): exact-trace moments fail the analytic reference at TOL=$TOL"; exit 1; }
  echo "PASS($1): exact trace, deterministic AND matches the analytic reference to $TOL"
}

check_exact KG  "inline_compute-kpm-nonEqOp chain1d VX VX 32"  "NonEqOp*chain1d*chebmom2D"    chain1d_analytic_test
check_exact DOS "inline_compute-kpm-spectralOp chain1d 1 32"   "SpectralOp1chain1d*chebmom1D" chain1d_dos_test

# (3) random-phase path still deterministic under a fixed seed (keeps the RNG tested)
gen_rand() {
  local T; T="$(mktemp -d)"; mkdir -p "$T/operators"
  cp "$OPDIR/chain1d.HAM.CSR" "$OPDIR/chain1d.VX.CSR" "$T/operators/"
  ( cd "$T"; printf -- "-2 2\n" > BOUNDS
    KPM_SEED=12345 "$BUILD/inline_compute-kpm-nonEqOp" chain1d VX VX 32 >/dev/null 2>&1 )
  echo "$T/$(cd "$T" && ls NonEqOp*chain1d*chebmom2D | head -1)"
}
RA="$(gen_rand)"; RB="$(gen_rand)"
diff -q "$RA" "$RB" >/dev/null \
  || { echo "FAIL(rng): random-phase path non-deterministic under fixed KPM_SEED"; exit 1; }
echo "PASS(rng): random-phase path is deterministic under fixed KPM_SEED"
echo "PASS: chain1d moments are exact-trace deterministic, physically correct, and the RNG path is seed-stable"
