#!/usr/bin/env bash
# Phase 1 invariant guard: the reconstruction cutoff (alpha) must stay OUT of the moment
# rescale, HBAR must have a single home, and alpha must be a runtime parameter.
#
# Three checks:
#   (1) SACRED MOMENT RESCALE. Every ScaleFactor()/ShiftFactor() definition uses a CUTOFF
#       constant and NEVER the reconstruction cutoff (safety_factors / recon_cutoff / KPM_ALPHA /
#       LSQUANT_RECON_ALPHA). Coupling alpha into the moment rescale under-fills [-1,1] and
#       breaks the analytic moment identities (chain1d_analytic / chain1d_dos).
#   (2) SINGLE HBAR. The literal 0.6582119624 appears only in include/units.hpp.
#   (3) RUNTIME ALPHA. The reconstruction grid responds to LSQUANT_RECON_ALPHA at runtime:
#       a larger alpha yields a wider reconstruction energy range (alpha reaches the grid).
#
# Usage: safety_factors_guard.sh <BUILD_DIR>
set -euo pipefail
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BUILD="$(cd "${1:-$REPO/build}" && pwd)"
fail=0

echo "[1] moment rescale (ScaleFactor/ShiftFactor) must use CUTOFF only, never alpha"
defs=$(grep -rnE 'double (ScaleFactor|ShiftFactor)\(\) *const' "$REPO/include" || true)
echo "$defs" | sed 's/^/    /'
while IFS= read -r line; do
  [ -z "$line" ] && continue
  if echo "$line" | grep -qiE 'recon_cutoff|safety_factors|KPM_ALPHA|RECON_ALPHA'; then
    echo "    FAIL: alpha leaked into a moment-rescale factor: $line"; fail=1
  fi
  if ! echo "$line" | grep -qE 'CUTOFF'; then
    echo "    FAIL: moment-rescale factor does not reference a CUTOFF constant: $line"; fail=1
  fi
done <<< "$defs"
[ "$fail" -eq 0 ] && echo "    OK: moment rescale is CUTOFF-only (alpha decoupled)"

echo "[2] HBAR single-sourced in include/units.hpp"
hits=$(grep -rln '0\.6582119624' "$REPO/include" "$REPO/src" || true)
echo "$hits" | sed 's/^/    /'
nonunits=$(echo "$hits" | grep -v 'include/units.hpp' || true)
if [ -n "$nonunits" ]; then echo "    FAIL: HBAR literal outside units.hpp:"; echo "$nonunits" | sed 's/^/      /'; fail=1
else echo "    OK: only include/units.hpp carries the HBAR literal"; fi

echo "[3] alpha is runtime (reconstruction energy range scales with LSQUANT_RECON_ALPHA)"
GMOM="$REPO/test/golden/NonEqOpVX-VXgraphene_goldenKPM_M64x64_statedefault.chebmom2D"
DRV="$BUILD/inline_kuboBastinFromChebmom"
if [ -x "$DRV" ] && [ -f "$GMOM" ]; then
  rng() { # max |E| of the reconstruction ENERGY axis (col 1, always finite) for a given alpha.
          # Note col 2 (sigma) is all-NaN at alpha=1.0 -- that IS Bug 2, the reason alpha<1 ships.
    local a="$1"; local T; T="$(mktemp -d)"; cp "$GMOM" "$T/"
    ( cd "$T"; LSQUANT_RECON_ALPHA="$a" "$DRV" "$(ls *.chebmom2D)" 10 >/dev/null 2>&1
      awk 'function abs(x){return x<0?-x:x} {if(abs($1)>m)m=abs($1)} END{printf "%.4f", m}' KuboBastin*.dat )
    rm -rf "$T"
  }
  r95=$(rng 0.95); r100=$(rng 1.00)
  echo "    max|E|: alpha=0.95 -> $r95 ;  alpha=1.00 -> $r100"
  awk -v a="$r95" -v b="$r100" 'BEGIN{ if (b > a*1.02) exit 0; else exit 1 }' \
    || { echo "    FAIL: alpha did not widen the grid at runtime (env not honored?)"; fail=1; }
  [ "$fail" -eq 0 ] && echo "    OK: alpha=1.00 grid wider than alpha=0.95 -> runtime-tunable"
else
  echo "    SKIP: driver or graphene moments unavailable"
fi

[ "$fail" -eq 0 ] && { echo "PASS: alpha decoupled from moments, HBAR single-sourced, alpha runtime"; exit 0; }
echo "FAIL: Phase-1 safety-factor invariants violated"; exit 1
