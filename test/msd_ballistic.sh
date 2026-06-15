#!/usr/bin/env bash
# Ballistic transport on the clean 1D chain via the mean-square displacement (Phase 4 #4).
#
# For the clean chain the wavepacket spreads ballistically: MSD(E,t) ~ t^2 (no scattering),
# and the energy-resolved transport coefficient tracks the conductivity slope
#     dsigma/dt(E) = (1/pi)*sqrt(4 - E^2)              (oracle: lsquant_derive_transport.py)
# i.e. its energy dependence is sqrt(4 - E^2) (0.636620 at E=0, 0.551329 at E=1).
#
# The MSD reconstruction (inline_timeCorrelationsFromChebmom) returns a(E,t) whose ballistic
# coefficient a(E) = MSD/t^2 is proportional to that conductivity slope. We assert, by EXACT
# TRACE (local state over the full basis -> deterministic, machine precision, NO random
# vectors), two oracle-grounded facts that are independent of the absolute normalization:
#   (1) ballistic scaling: a(E,t) = MSD/t^2 is constant in t (rel. spread < TOL_T);
#   (2) energy law: a(E)/a(0) = sqrt(4 - E^2)/2 at E in {0.5,1.0,1.5} (rel. err < TOL_E).
#
# (The velocity-autocorrelation route gives the SAME energy profile -- representation
# equivalence -- but inline_compute-kpm-CorrOp's state-file/numStates interface is unreliable
# for a deterministic golden, so VAC is not committed here. See the PR notes.)
#
# Usage: msd_ballistic.sh <BUILD_DIR> [M] [TOL_E] [TOL_T]
set -euo pipefail
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BUILD="$(cd "${1:-$REPO/build}" && pwd)"
M="${2:-256}"
TOL_E="${3:-0.02}"
TOL_T="${4:-0.02}"
OPDIR="$REPO/test/golden/chain1d_ballistic/operators"
MSD="$BUILD/inline_compute-kpm-MeanSquareDisplacement"
RECON="$BUILD/inline_timeCorrelationsFromChebmom"

[ -x "$MSD" ]   || { echo "FAIL: MSD driver missing ($MSD)"; exit 1; }
[ -x "$RECON" ] || { echo "FAIL: time-correlation reconstruction missing ($RECON)"; exit 1; }
[ -f "$OPDIR/chain1d.HAM.CSR" ] || { echo "FAIL: committed chain operators missing"; exit 1; }

N=$(head -1 "$OPDIR/chain1d.HAM.CSR" | awk '{print $1}')   # SystemSize

T="$(mktemp -d)"; trap 'rm -rf "$T"' EXIT
mkdir -p "$T/operators"
cp "$OPDIR"/chain1d.{HAM,VX}.CSR "$T/operators/"
( cd "$T"
  printf -- "-2 2\n" > BOUNDS                            # exact band of the clean chain
  { echo local; echo "$N"; seq 0 $((N-1)); } > exact     # exact trace over the full basis
  echo ">> MSD moments (M=$M, numTimes=64, tmax=20, exact trace, dim=$N)"
  "$MSD" chain1d VX "$M" 64 20 exact >msd.log 2>&1 \
    || { echo "FAIL: MSD computation crashed"; tail -15 msd.log; exit 1; }
  MOM="$(ls Correlation*chain1d*chebmomTD | head -1)"

  # reconstruct a(E,t) = MSD/t^2 at each Fermi energy; check ballistic constancy in t,
  # then the energy profile against the analytic sqrt(4-E^2)/2.
  a0=""; fail=0
  for E in 0 0.5 1.0 1.5; do
    "$RECON" "$MOM" 10 "$E" >/dev/null 2>&1
    DAT="$(ls mean*EF*.dat | head -1)"
    # ballistic: a=MSD/t^2 over the window 4<t<12; report mean and relative spread
    aval=$(awk '$1>4 && $1<12 { a=$2/($1*$1); s+=a; n++ } END{ printf "%.8f", s/n }' "$DAT")
    spread=$(awk '$1>4 && $1<12 { a=$2/($1*$1); s+=a; ss+=a*a; n++ }
      END{ m=s/n; v=ss/n-m*m; sd=(v>0?sqrt(v):0); printf "%.5f", (m!=0? sd/m : 1) }' "$DAT")
    rm -f "$DAT"
    awk -v e="$E" -v a="$aval" -v sp="$spread" -v tt="$TOL_T" 'BEGIN{
      printf "E=%-4s a=MSD/t^2 = %-12s  (t-spread %.4f, tol %.3f)%s\n", e, a, sp, tt, (sp>tt?"  <-- NOT ballistic":"") }'
    awk -v sp="$spread" -v tt="$TOL_T" 'BEGIN{exit (sp>tt?1:0)}' || fail=1
    [ -z "$a0" ] && a0="$aval"
    ana=$(awk -v e="$E" 'BEGIN{printf "%.6f", sqrt(4-e*e)/2}')
    rel=$(awk -v a="$aval" -v a0="$a0" -v an="$ana" 'BEGIN{r=(a/a0); printf "%.5f", (r-an>0?r-an:an-r)/an}')
    awk -v e="$E" -v r="$(awk -v a=$aval -v a0=$a0 'BEGIN{printf "%.4f", a/a0}')" -v an="$ana" -v rel="$rel" -v te="$TOL_E" 'BEGIN{
      printf "        a(E)/a(0) = %-8s  analytic sqrt(4-E^2)/2 = %-8s  rel.err %.4f (tol %.3f)%s\n", r, an, rel, te, (rel>te?"  <-- OFF":"") }'
    awk -v rel="$rel" -v te="$TOL_E" 'BEGIN{exit (rel>te?1:0)}' || fail=1
  done
  [ "$fail" -eq 0 ] && echo "PASS: ballistic MSD ~ t^2 and the sigma-slope energy law sqrt(4-E^2) hold (exact trace)" \
                    || { echo "FAIL: ballistic scaling or energy law outside tolerance"; exit 1; }
)
