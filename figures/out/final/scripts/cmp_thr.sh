#!/bin/bash
# item 3: cmp every file written at T=1 against T=16 (Aron, same build), per run file; run.resolved.json excluded
export PATH=$HOME/lsq-env/bin:$PATH LC_ALL=C
O=$HOME/fv2/thr; R=$HOME/lsquant-final2; ti=0; td=0
for d in $O/T16/*/; do
  n=$(basename $d); i=0; x=0; miss=0
  while IFS= read -r f; do
    [ "$f" = run.resolved.json ] && continue
    a=$O/T16/$n/$f; b=$O/T1/$n/$f
    if [ ! -f "$b" ]; then miss=$((miss + 1)); echo "DIFFERS   $n/$f (not written at T=1)"
    elif cmp -s "$a" "$b"; then i=$((i + 1))
    else x=$((x + 1)); echo "DIFFERS   $n/$f: $(python $R/scripts/golden_diff.py "$b" "$a" 2>&1 | tr '\n' ' ' | cut -c1-400)"; fi
  done < <(cd $O/T16/$n && find . -type f | sed 's|^\./||' | sort)
  m1=$(cd $O/T1/$n 2>/dev/null && find . -type f | sed 's|^\./||' | grep -vxF -f <(cd $O/T16/$n && find . -type f | sed 's|^\./||') | wc -l)
  printf "%-42s T1 vs T16: %4d identical %4d differ %3d missing %3d extra-at-T1\n" "$n" $i $x $miss $m1
  ti=$((ti + i)); td=$((td + x + miss + m1))
done
echo "TOTAL T1 vs T16: $ti identical, $td differ, missing or extra; run dirs T16 $(ls -d $O/T16/*/ | wc -l), T1 $(ls -d $O/T1/*/ | wc -l)"
