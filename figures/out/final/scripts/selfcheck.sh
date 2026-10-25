#!/bin/bash
export PATH=$HOME/lsq-env/bin:$PATH LC_ALL=C KPM_SEED=12345 OMP_NUM_THREADS=16
R=$HOME/lsquant-final2; B=$HOME/lsquant-final2-build; L=$HOME/fv2/self; mkdir -p $L
export LSQUANT_BUILD_DIR=$B PYTHONPATH=$B/lsquant_core:$R/python
cd $R
s=$(date +%s); $B/lsquant_core/stage_1_check > $L/stage_1.log 2>&1; echo "stage_1_check rc=$? wall $(( $(date +%s) - s )) s" >> $L/summary.txt
for k in 2 3 4 5 6 7; do
  s=$(date +%s); python tests/dev/stage_${k}_check.py > $L/stage_$k.log 2>&1; echo "stage_${k}_check.py rc=$? wall $(( $(date +%s) - s )) s" >> $L/summary.txt
done
echo SELF_DONE >> $L/summary.txt
