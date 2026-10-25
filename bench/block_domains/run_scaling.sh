#!/usr/bin/env bash
# Thread-scaling campaign of block_kpm: block sizes 5, 10, 30 at dimension 1e6,
# global vs local, without and with the neighbour-block coupling.
# Writes results.jsonl and per-run moment files in out/.
set -u
cd "$(dirname "$0")"
export OMP_PROC_BIND=spread OMP_PLACES=cores OMP_DYNAMIC=false
DIM=1000000
mkdir -p out
: > results.jsonl
for coupled in 0 1; do
  for w in 5 10 30; do
    for mode in dos cond; do
      if [ "$mode" = dos ]; then M=256; R=4; else M=128; R=1; fi
      for method in global local; do
        # With coupling the local method is not exact: run it once for its error only.
        if [ "$coupled" = 1 ] && [ "$method" = local ]; then threads="64"; else threads="1 2 4 8 16 32 64"; fi
        for t in $threads; do
          OMP_NUM_THREADS=$t ./block_kpm "$mode" "$method" "$w" "$DIM" "$M" "$R" "$coupled" \
            "out/${mode}-${method}-w${w}-c${coupled}-t${t}.bin" >> results.jsonl
        done
      done
    done
  done
done
echo done > campaign.done
