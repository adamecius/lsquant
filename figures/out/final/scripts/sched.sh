#!/bin/bash
# Weighted job scheduler for the final validation on Aron: sum of OMP_NUM_THREADS over my running
# jobs <= CAP. Jobs file lines: name|weight|precondition|flags|command  (precondition: - , file:<path>, done:<name>;
# flags: R = reserve, i.e. when eligible and not fitting, launch nothing below it). Only PIDs started here are tracked.
CAP=${CAP:-90}; J=$1; D=$HOME/fv2/sched; mkdir -p $D/logs
declare -A PID ST W
mapfile -t LINES < <(grep -v '^#' $J | grep -v '^$')
log() { echo "$(date -u +%FT%TZ) $*" >> $D/sched.log; }
log "START cap=$CAP jobs=${#LINES[@]} scheduler pid $$"
while :; do
  cur=0; pend=0
  for l in "${LINES[@]}"; do
    IFS='|' read -r n w pre fl cmd <<<"$l"
    if [ "${ST[$n]:-}" = run ]; then
      if kill -0 ${PID[$n]} 2>/dev/null; then cur=$((cur + w)); else wait ${PID[$n]} 2>/dev/null; rc=$?; ST[$n]=done; log "END $n pid ${PID[$n]} rc=$rc"; fi
    fi
  done
  for l in "${LINES[@]}"; do
    IFS='|' read -r n w pre fl cmd <<<"$l"
    [ -z "${ST[$n]:-}" ] || continue
    pend=$((pend + 1))
    case "$pre" in
      file:*) [ -e "${pre#file:}" ] || continue ;;
      grep:*) g="${pre#grep:}"; grep -qF -- "${g%@*}" "${g##*@}" 2>/dev/null || continue ;;
      done:*) ok=1; for x in ${pre#done:}; do [ "${ST[$x]:-}" = done ] || ok=0; done; [ $ok = 1 ] || continue ;;
    esac
    if [ $((cur + w)) -le $CAP ]; then
      nohup bash -c "$cmd" > $D/logs/$n.log 2>&1 &
      PID[$n]=$!; ST[$n]=run; cur=$((cur + w)); echo "${PID[$n]} $n" >> $D/pids; log "START $n pid ${PID[$n]} weight $w (load $cur)"
    elif [ "$fl" = R ]; then break; fi
  done
  nrun=0; for l in "${LINES[@]}"; do IFS='|' read -r n _ <<<"$l"; [ "${ST[$n]:-}" = run ] && nrun=$((nrun + 1)); done
  [ $pend -eq 0 ] && [ $nrun -eq 0 ] && break
  sleep 15
done
log "ALL_DONE"
