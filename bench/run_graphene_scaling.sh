#!/usr/bin/env bash
# Phase-F scaling benchmark of the THREADED lsquant SpMV on a graphene (complex)
# Hamiltonian. Two sweeps, both parsing the reporting-layer region times that Phase A
# added ("spmv" and "compute_spectral/SpectralMoments"):
#   (1) thread scaling : fixed N, threads {1..256}, DOS recursion (M moments)
#   (2) size scaling   : fixed threads, N swept, confirms linear-in-N (KPM)
# GB/s for the SpMV stream is reported as % of the measured STREAM read roofline.
set -u
cd "$(dirname "$0")"
LSQ=/home/jgarcia/lsquant/build/lsquant
GEN="python3 make_graphene.py"
WORK=gr; mkdir -p "$WORK"
OUT_T=graphene_threads.csv
OUT_N=graphene_size.csv
# Measured STREAM read roofline (perf/STEP2_ceilings.txt), GB/s
ROOF_SOCKET=348.2 ; ROOF_DUAL=700.8

run_one() { # run_one <label> <M> <threads> <numa>  -> echoes "spmv_ms spmv_calls recur_ms"
  local label="$1" M="$2" T="$3" numa="$4"
  local nc=""
  case "$numa" in sock0) nc="numactl --cpunodebind=0 --membind=0";; interleave) nc="numactl --interleave=all";; esac
  printf '{ "mode":"spectral", "label":"%s", "operator":"1", "num_moments":%d }\n' "$label" "$M" > "$WORK/dos.json"
  local log
  log=$(cd "$WORK" && OMP_NUM_THREADS=$T OMP_PROC_BIND=close OMP_PLACES=cores $nc \
         "$LSQ" compute --config dos.json 2>&1)
  rm -f "$WORK"/SpectralOp1*.chebmom1D
  local spmv_ms spmv_calls recur_ms
  spmv_ms=$(echo "$log"   | grep -oE 'spmv[[:space:]]+[0-9.]+ ms' | grep -oE '[0-9.]+' | head -1)
  spmv_calls=$(echo "$log"| grep -E 'spmv ' | grep -oE '[0-9]+ calls' | grep -oE '[0-9]+' | head -1)
  recur_ms=$(echo "$log"  | grep -oE 'SpectralMoments[[:space:]]+[0-9.]+ ms' | grep -oE '[0-9.]+' | head -1)
  echo "${spmv_ms:-NA} ${spmv_calls:-NA} ${recur_ms:-NA}"
}

gbs() { # gbs <nnz> <N> <spmv_ms> <calls> -> matrix-stream GB/s
  awk -v nnz="$1" -v N="$2" -v ms="$3" -v c="$4" 'BEGIN{
    if(ms=="NA"||c=="NA"){print "NA"; exit}
    bytes=(nnz*(16.0+4.0)+(N+1)*8.0)*c; printf "%.1f", bytes/(ms/1000.0)/1e9 }'
}

echo "=== (1) THREAD SCALING : graphene fixed N, DOS M=512 ==="
TL=1000   # N = 2*L^2 = 2,000,000 sites
read LBL N NNZ < <($GEN $TL 0.10 0.5 "$WORK" | sed -E 's/built ([^:]*): N=([0-9]+) sites, nnz=([0-9]+).*/\1 \2 \3/')
echo "  fixed system: $LBL N=$N nnz=$NNZ"
echo "system,N,nnz,M,threads,numa,spmv_ms,spmv_calls,spmv_per_call_ms,recursion_ms,gbs,pct_roofline" > "$OUT_T"
for T in 1 2 4 8 16 32 64 128 256; do
  numa=sock0; roof=$ROOF_SOCKET
  if [ $T -gt 128 ]; then numa=interleave; roof=$ROOF_DUAL; fi
  read sp ca rc < <(run_one "$LBL" 512 "$T" "$numa")
  g=$(gbs "$NNZ" "$N" "$sp" "$ca")
  pc=$(awk -v g="$g" -v r="$roof" 'BEGIN{ if(g=="NA"){print "NA"}else printf "%.1f", 100*g/r }')
  pcm=$(awk -v s="$sp" -v c="$ca" 'BEGIN{ if(s=="NA"||c=="NA"){print "NA"}else printf "%.4f", s/c }')
  echo "$LBL,$N,$NNZ,512,$T,$numa,$sp,$ca,$pcm,$rc,$g,$pc" | tee -a "$OUT_T"
done

echo "=== (2) SIZE SCALING : graphene N swept, DOS M=256, fixed 128 threads sock0 ==="
echo "system,N,nnz,M,threads,numa,spmv_ms,spmv_calls,spmv_per_call_ms,recursion_ms,gbs,pct_roofline" > "$OUT_N"
for L in 80 160 320 480 640 800 1000; do
  read LBL N NNZ < <($GEN $L 0.10 0.5 "$WORK" | sed -E 's/built ([^:]*): N=([0-9]+) sites, nnz=([0-9]+).*/\1 \2 \3/')
  read sp ca rc < <(run_one "$LBL" 256 128 sock0)
  g=$(gbs "$NNZ" "$N" "$sp" "$ca")
  pc=$(awk -v g="$g" -v r="$ROOF_SOCKET" 'BEGIN{ if(g=="NA"){print "NA"}else printf "%.1f", 100*g/r }')
  pcm=$(awk -v s="$sp" -v c="$ca" 'BEGIN{ if(s=="NA"||c=="NA"){print "NA"}else printf "%.4f", s/c }')
  echo "$LBL,$N,$NNZ,256,128,sock0,$sp,$ca,$pcm,$rc,$g,$pc" | tee -a "$OUT_N"
  rm -f "$WORK"/operators/${LBL}.HAM.CSR
done
echo "=== DONE: $OUT_T , $OUT_N ==="