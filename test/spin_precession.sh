#!/usr/bin/env bash
# Spin precession on the magnetic chain (Phase 4 #6).
#
# The magnetic chain is a constant Zeeman field H = H_hop - J_ex*sigma_z (J_ex=0.1). A spin
# prepared along +x precesses in the x-y plane at the Larmor frequency set by the Zeeman
# splitting 2*J_ex. In LinQT's PHYSICAL units (hbar = 0.6582119624 eV*fs, time in fs):
#
#     <S_x(t)> =  cos(omega t),  <S_y(t)> = -sin(omega t),  <S_z(t)> = 0,
#     omega = 2*J_ex/hbar = 0.2/0.6582119624 = 0.30385 rad/fs   (period 20.68 fs),
#
# independent of E_F. This follows the thesis spin-dynamics formalism (Eqs. 175/176):
# |phi(0)> = 1/2 (1 + x.s)|phi_r> is the +x-projected state and
# S(E,t) = <phi(t)|s delta(E-H)|phi(t)> / <phi(t)|delta(E-H)|phi(t)>. Since the precession is
# E-independent, the ENERGY-INTEGRATED expectation is exactly the m=0 Chebyshev moment:
# integral dE of the delta-reconstruction picks T_0, so <phi(t)|s|phi(t)> = mu_s(0,n), and the
# (constant) denominator <phi(t)|phi(t)> = mu_Sx(0,0). No energy reconstruction is needed.
#
# Method: the DENSITY projection path (chebyshev::TimeEvolvedProjectedOperator), which realises
# |phi(t)><phi(t)| with |phi(0)> = P_x|phi_r>, P_x = 1/2(I+S_x). EXACT TRACE (local state over
# the full basis) -> deterministic, machine precision, no random vectors. We read mu_OP(0,n) for
# OP in {S_x,S_y,S_z}, normalise by mu_Sx(0,0), and compare to the analytic precession.
#
# Usage: spin_precession.sh <BUILD_DIR> [M] [NTIME] [TMAX_fs] [TOL]
set -euo pipefail
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BUILD="$(cd "${1:-$REPO/build}" && pwd)"
M="${2:-32}"; NTIME="${3:-96}"; TMAX="${4:-95}"; TOL="${5:-0.01}"
DRV="$BUILD/inline_compute-kpm-TimeEvProjetedOp"
OPD="$REPO/test/models/chain1d_mag/operators"

[ -x "$DRV" ] || { echo "FAIL: TimeEvProjetedOp driver missing ($DRV)"; exit 1; }
for f in HAM SX SY SZ PX; do
  [ -f "$OPD/chain1d_mag.$f.CSR" ] || { echo "FAIL: committed operator chain1d_mag.$f.CSR missing"; exit 1; }
done
N=$(head -1 "$OPD/chain1d_mag.HAM.CSR" | awk '{print $1}')

T="$(mktemp -d)"; trap 'rm -rf "$T"' EXIT
mkdir -p "$T/operators"
cp "$OPD"/chain1d_mag.{HAM,SX,SY,SZ,PX}.CSR "$T/operators/"
( cd "$T"
  { echo local; echo "$N"; seq 0 $((N-1)); } > exact     # exact trace over the full basis
  for OP in SX SY SZ; do
    "$DRV" chain1d_mag "$OP" PX "$M" "$NTIME" "$TMAX" exact >"log_$OP" 2>&1 \
      || { echo "FAIL: TimeEvProjetedOp $OP crashed"; tail -8 "log_$OP"; exit 1; }
    mv "$(ls TimeProj${OP}-PX*stateexact*chebmomTD | head -1)" "mom_$OP.chebmomTD"
  done

  # read mu_OP(0,n): m=0 row of the m-major file (line1 header, line2 "M Nt", then re/im pairs)
  # one moment (re im) per line, m-major (index m*Nt+n) -> first Nt lines are the m=0 row.
  read_m0() { awk 'NR==1{next} NR==2{Nt=$2; next} {v[c++]=$1} END{for(n=0;n<Nt;n++) print v[n]}' "$1"; }
  DT=$(awk 'NR==1{print $5; exit}' mom_SX.chebmomTD)
  paste <(read_m0 mom_SX.chebmomTD) <(read_m0 mom_SY.chebmomTD) <(read_m0 mom_SZ.chebmomTD) \
  | awk -v dt="$DT" -v tol="$TOL" '
      function abs(x){return x<0?-x:x}
      BEGIN{ hbar=0.6582119624; w=2*0.1/hbar }
      NR==1{ D=$1 }                                   # normalisation = S_x(0)
      { n=NR-1; t=n*dt;
        sx=$1/D; sy=$2/D; sz=$3/D;
        ex=abs(sx-cos(w*t)); ey=abs(sy-(-sin(w*t))); ez=abs(sz);
        if(ex>mex)mex=ex; if(ey>mey)mey=ey; if(ez>mez)mez=ez;
        if(n%16==0) printf "  t=%6.2f fs  Sx=%+.4f(cos=%+.4f)  Sy=%+.4f(-sin=%+.4f)  Sz=%+.4f\n", t, sx, cos(w*t), sy, -sin(w*t), sz;
      }
      END{
        printf "omega=2Jex/hbar=%.5f rad/fs (period %.2f fs); max|dSx|=%.4f max|dSy|=%.4f max|Sz|=%.4f  tol=%.3f\n", w, 6.283185/w, mex, mey, mez, tol;
        if(mex<=tol && mey<=tol && mez<=tol){ print "PASS: spin precession <S(t)>=cos(wt)x - sin(wt)y, S_z=0 (density flavor, exact trace)"; }
        else { print "FAIL: precession outside tolerance"; exit 1; }
      }'
)
