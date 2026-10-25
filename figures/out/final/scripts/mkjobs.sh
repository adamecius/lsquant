#!/bin/bash
# writes ~/fv2/jobs.txt
R=$HOME/lsquant-final2; B=$HOME/lsquant-final2-build; F=$HOME/fv2
ENVX="export PATH=$HOME/lsq-env/bin:\$PATH LC_ALL=C KPM_SEED=12345 W2S_BIN=$HOME/from_fatnode/wannier2sparse/build/wannier2sparse TMPDIR=$F/tmp; cd $R"
GATE="GOLDEN_TOL=1 bash scripts/check_goldens.sh"
LB="g1 g2 g3 g4 g4c g5 g6 g7 g7_sx g7_sy g7_sz tutorial01 tutorial02 tutorial03 tutorial04 tutorial05 $(cd $R/runs/goldens && ls g9_*.json | sed 's/\.json//' | tr '\n' ' ')"
LC="$(cd $R/runs/goldens && ls g10_*.json | sed 's/\.json//' | tr '\n' ' ')"
CSR="grep:python -m lsquant run runs/goldens@$F/sched/logs/gate7_b.log"
LONG="tutorial06 g10_chain1d_complex_lowdin g10_chain1d_real_lowdin g10_chain1d_real_lowdin_chebyshev g10_chain1d_real_generalised_chebyshev"
REST=$(cd $R/runs/goldens && ls *.json | sed 's/\.json$//' | grep -vxF -f <(echo $LONG | tr ' ' '\n'))
one() { echo "$ENVX; export PYTHONPATH=$B/lsquant_core:$R/python; mkdir -p $F/thr/T$1; rm -rf $F/thr/T$1/$2; s=\$(date +%s); OMP_NUM_THREADS=$1 python -m lsquant run runs/goldens/$2.json --out $F/thr/T$1/$2 > $F/thr/T$1/$2.out 2> $F/thr/T$1/$2.err; rc=\$?; echo \"T=$1 $2 exit \$rc wall \$((\$(date +%s) - s)) s\" >> $F/thr/runs.log"; }
{
echo "gate7_b|48|-||$ENVX; S7_THREADS=48 S7_SET=list S7_LIST='$LB' $GATE 7runs $B"
echo "gate7_a|48|-||$ENVX; S7_THREADS=48 S7_SET=list S7_LIST='tutorial06' $GATE 7runs $B"
echo "gate7_c|48|-||$ENVX; S7_THREADS=48 S7_SET=list S7_LIST='$LC' $GATE 7runs $B"
echo "gate7_tail|64|done:gate7_a gate7_b gate7_c|R|$ENVX; S7_LEGACY_WT=$F/legacy_wt $GATE 7tail $B"
echo "gate6only|16|-||$ENVX; S6_THREADS=16 S6_KEEP=$F/s6keep $GATE 6only $B"
echo "selfcheck|16|$CSR||bash $F/selfcheck.sh"
for n in $LONG; do echo "t1_$n|1|$CSR||$(one 1 $n)"; done
for n in $LONG; do echo "t16_$n|16|$CSR||$(one 16 $n)"; done
for n in $REST; do echo "t1_$n|1|$CSR||$(one 1 $n)"; done
for n in $REST; do echo "t16_$n|16|$CSR||$(one 16 $n)"; done
} > $F/jobs.txt
wc -l $F/jobs.txt; cut -d'|' -f1-4 $F/jobs.txt | head -8
