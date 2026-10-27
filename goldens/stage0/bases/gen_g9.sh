#!/bin/bash
# Generate the G9 goldens (stage 5, handoff 0005 codex.md section 2.2, revised after run 1)
# from the independent oracle wannier2sparse/tools/hr_exactdiag.py.
# usage: goldens/stage0/bases/gen_g9.sh <repo>   (writes <repo>/goldens/stage0/bases/<model>/)
#
#   G9 (gated):        hr_exactdiag dos and sigma --comp xx on the Gamma-centred 48 x 48 mesh
#                      k = n/48 (through <repo>/scripts/hr_exactdiag_gamma.py), Gaussian eta = 40 meV
#   finite size (report only): the same on hr_exactdiag's own shifted 240 x 240 mesh (n + 1/2)/240
set -euo pipefail
REPO="$(cd "$1" && pwd)"; W2S=~/wannier2sparse; TOOL="$W2S/tools/hr_exactdiag.py"
NK=48; NKFS=240; ETA=0.040; EMIN=-3.6; EMAX=3.6; NG=2881      # grid spacing 2.5 meV on [-3.6, 3.6]
OUT="$REPO/goldens/stage0/bases"
WORK=$(mktemp -d /tmp/g9.XXXXXX); cd "$WORK"
python3 "$W2S/examples/gen_models.py" >/dev/null
for m in graphene haldane; do
  cd "$WORK/$m"
  "$W2S/build/wannier2sparse" "$m" 1 1 1 VX VY --mode bundle -o bund --quiet >/dev/null
  cp bund/$m.w2sp/operators/VX.hr.dat ${m}_vx_hr.dat
  cp bund/$m.w2sp/operators/VY.hr.dat ${m}_vy_hr.dat
  D="$OUT/$m"; rm -rf "$D"; mkdir -p "$D/model"
  cp ${m}_hr.dat $m.uc $m.xyz ${m}_vx_hr.dat ${m}_vy_hr.dat "$D/model/"
  A="--eta $ETA --ngrid $NG --emin $EMIN --emax $EMAX"
  python3 "$REPO/scripts/hr_exactdiag_gamma.py" "$TOOL" dos   $m --nk $NK $A --out "$D/G9_dos"
  python3 "$REPO/scripts/hr_exactdiag_gamma.py" "$TOOL" sigma $m --nk $NK $A --comp xx --out "$D/G9_sigma_xx"
  python3 "$TOOL" dos   $m --nk $NKFS $A --out "$D/finite_size_shifted240_dos"
  python3 "$TOOL" sigma $m --nk $NKFS $A --comp xx --out "$D/finite_size_shifted240_sigma_xx"
  printf "mesh gamma_centred_k=n/N\nnk %s\neta_eV %s\nemin %s\nemax %s\nngrid %s\nfinite_size_mesh shifted_k=(n+1/2)/N\nfinite_size_nk %s\nsigma_area_hr_A2 |a1_x*a2_y|\n" \
    $NK $ETA $EMIN $EMAX $NG $NKFS > "$D/G9_params.txt"
  echo "$m: G9 Gamma 48x48 and shifted 240x240, eta $ETA eV, grid [$EMIN, $EMAX] x $NG"
done
rm -rf "$WORK"
