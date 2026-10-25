# Stage 8 measurements on Aron

These drivers are measurement tools, not gates. Run from the checkout root in
`~/lsq-env`, with the selected build's `lsquant_core` directory and `python/` on
`PYTHONPATH`. Keep outputs outside the checkout. Never overwrite a baseline.

Item 1, before any kernel change (repeat after each numbered item):

```bash
export PYTHONPATH="$PWD/python:$PWD/build-golden/lsquant_core"
export KPM_SEED=12345 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export OMP_NUM_THREADS=48 OMP_PROC_BIND=close OMP_PLACES=cores
python bench/profile_run.py runs/goldens/g5.json --out /tmp/stage8/item1-g5 --perf
python bench/profile_run.py "$BI_RUN" --moments 128 --states 2 --out /tmp/stage8/item1-bi --perf
```

`BI_RUN` must name the tester's 6b Bi block spin Hall run file, including its
local SIESTA input paths. Keep its trace mode, word, metric and reconstruction
unchanged; record its path, contents and input checksums with the artifacts.
The short profile changes only M and vector count in memory.

Record `profile.txt` (top 20 cumulative entries), `profile.json` (C++ attributed
self time and fraction), and `occupancy.jsonl` (thread counts and ps CPU percent,
which is the process lifetime average, not an instantaneous sample). Inspect the
attributed entry names: CPython labels pybind11 calls as `PyCapsule` builtins;
unlabelled extension calls are excluded, so this fraction is an attribution,
not a sampling profiler's measurement of all native work. Imports are excluded.
`perf.csv` is an additional unprofiled replay, with its own exit status; missing
perf or denied counters must be recorded, never treated as zero measurements.
Record compiler, dependency versions, affinity, host load and governor as well.

Item 3: configure separate clean directories with the same dependencies and build
type. These commands require only the packages already used by LSQUANT:

```bash
cmake -S . -B build-golden -DCMAKE_BUILD_TYPE=Release -DCMAKE_PREFIX_PATH="$CONDA_PREFIX" -DLSQUANT_BUILD_PROFILE=golden -DCMAKE_EXPORT_COMPILE_COMMANDS=ON
cmake --build build-golden -j 24
cmake -S . -B build-native -DCMAKE_BUILD_TYPE=Release -DCMAKE_PREFIX_PATH="$CONDA_PREFIX" -DLSQUANT_BUILD_PROFILE=native -DCMAKE_EXPORT_COMPILE_COMMANDS=ON
cmake --build build-native -j 24
export PYTHONPATH="$PWD/python:$PWD/build-native/lsquant_core"
```

Golden adds no flags/defines to the historical command lines. Native appends
`-O3 -march=native -mtune=native -fno-math-errno` and `LSQUANT_NATIVE_BUILD=1`.
Neither profile enables fast-math or unsafe-math optimizations. Record the actual
compile commands (including inherited Release and pybind11 LTO flags). Native
headers carry `# build_profile=native` at the end of the first line; the existing
golden headers are untagged and mean golden. No numeric layout or precision changes.
The reader retains source provenance when a moment file is rewritten by another
profile. Receipts always name the loaded engine profile and reject replay under a
different one. Compare moment data only within one profile and one machine.

Item 4: import LSQUANT before another library initializes OpenMP. Defaults fill
only absent variables: close binding, core places, physical cores in the current
CPU affinity (so `numactl`/scheduler cpusets constrain the count). Set BLAS threads
to one to avoid nested teams. A static contiguous block group stays with the same
thread for initialization and multiplication; keep the team size fixed per process.
First touch applies to owned native CSR, block and recursion/batch buffers; on
Linux complete interior pages of newly allocated buffers are discarded before
the parallel fill. It never discards borrowed data. Allocation reuse and policy
still matter: record `numastat -p PID` and `/proc/PID/numa_maps`, not just affinity.

```bash
for threads in 24 48; do
  OMP_NUM_THREADS=$threads numactl --cpunodebind=0 --membind=0 python bench/profile_run.py "$BI_RUN" --moments 128 --states 2 --out /tmp/stage8/item4-bi-socket0-t$threads
done
OMP_NUM_THREADS=96 numactl --cpunodebind=0,1 --membind=0,1 python bench/profile_run.py "$BI_RUN" --moments 128 --states 2 --out /tmp/stage8/item4-bi-dual-t96
OMP_NUM_THREADS=48 numactl --cpunodebind=0 --membind=0 python bench/profile_run.py "$BI_RUN" --moments 128 --states 2 --out /tmp/stage8/item4-bi-p0 > /tmp/stage8/item4-p0.log 2>&1 &
echo $! > /tmp/stage8/item4-p0.pid
OMP_NUM_THREADS=48 numactl --cpunodebind=1 --membind=1 python bench/profile_run.py "$BI_RUN" --moments 128 --states 2 --out /tmp/stage8/item4-bi-p1 > /tmp/stage8/item4-p1.log 2>&1 &
echo $! > /tmp/stage8/item4-p1.pid
wait
```

The last recipe measures two independent jobs and aggregate throughput, not the
latency of one distributed trace. Keep their outputs separate; do not average
already normalized moment files to claim one exact trace. Compare the same total
work (two serial jobs versus these two concurrent jobs). Use the recorded child
PIDs in each output directory for numastat or stopping work, never a process name.

Item 5: enable the optional measurement executable in both build directories:

```bash
cmake -S . -B build-golden -DLSQUANT_BUILD_BENCHMARKS=ON
cmake --build build-golden -j 24
cmake -S . -B build-native -DLSQUANT_BUILD_BENCHMARKS=ON
cmake --build build-native -j 24
mkdir -p /tmp/stage8/kernels
for profile in golden native; do
  for threads in 1 24 48 96; do
    OMP_NUM_THREADS=$threads build-$profile/lsquant_core/stage8_kernels spmv 1000000 1 100 > /tmp/stage8/kernels/$profile-spmv-t$threads.json
    OMP_NUM_THREADS=$threads build-$profile/lsquant_core/stage8_kernels spmm 1000000 16 30 > /tmp/stage8/kernels/$profile-spmm-t$threads.json
    OMP_NUM_THREADS=$threads build-$profile/lsquant_core/stage8_kernels block 2304 32 30 > /tmp/stage8/kernels/$profile-block-t$threads.json
    OMP_NUM_THREADS=$threads build-$profile/lsquant_core/stage8_kernels dot 165888 128 3 > /tmp/stage8/kernels/$profile-dot-t$threads.json
    OMP_NUM_THREADS=$threads build-$profile/lsquant_core/stage8_kernels contract 256 7680 3 > /tmp/stage8/kernels/$profile-contract-t$threads.json
  done
done
for cutoff in 0 4096 16384 65536 262144; do
  for rows in 64 256 1024 4096 16384 65536; do
    OMP_NUM_THREADS=48 LSQUANT_SPMV_MIN_NNZ=$cutoff numactl --cpunodebind=0 --membind=0 build-native/lsquant_core/stage8_kernels spmv $rows 1 1000 > /tmp/stage8/kernels/cutoff-$cutoff-rows-$rows.json
  done
done
for batch in 1 2 4 8 16 32 64 128; do
  OMP_NUM_THREADS=48 LSQUANT_KERNEL_BATCH=$batch numactl --cpunodebind=0 --membind=0 build-native/lsquant_core/stage8_kernels block 2304 128 30 > /tmp/stage8/kernels/block-batch-$batch.json
  OMP_NUM_THREADS=48 LSQUANT_KERNEL_BATCH=$batch numactl --cpunodebind=0 --membind=0 build-native/lsquant_core/stage8_kernels spmm 1000000 128 10 > /tmp/stage8/kernels/csr-batch-$batch.json
done
```

Repeat at least three times with distinct output paths; report medians and range.
The cutoff default is an **unmeasured candidate** until this sweep establishes a
crossover on Aron. Hashes must agree between thread counts within a profile and
configuration; report any discrepancy, never bless a tolerance automatically.
Also compare the full production moment and curve files between thread counts.

Kernel microbatches and storage windows serve different purposes. For a CSR row
shard `q=ceil(N/T)`, the candidate microbatch is
`max(1,floor((512 KiB/(2*16)-3*q)/(2*q)))`, capped at the live vector count. If a
single stream exceeds L2, this returns one and the kernel streams from memory.
For block size W, the microbatch is
`max(1,floor((512 KiB/(2*16)-W*W)/(2*W)))`, likewise capped. The half-L2 budget
leaves space for Eigen packing, stack and other streams; it is a tuning starting
point. `LSQUANT_KERNEL_BATCH` overrides either for measurement. Chebyshev storage
windows keep the requested `expansion.batch` (default M), preventing repeated
right-side recursions merely to satisfy a cache estimate. The moment table packs
the stored window and evaluates fixed 16x16 output tiles by GEMM; contraction
uses fixed groups of 32 energy points. Geometry and reduction trees of each dense
product do not depend on the team size.

FFT is optional via `kubo_fft_grid(..., backend="fft")` or
`LSQUANT_FFT_BACKEND=fft` for `bastin_fft`/`greenwood_fft` reconstruction. The default
`tables` uses the historical scalar contraction in golden, GEMM in native. Compare
FFT and tables on **identical FFT nodes**; uniform-grid curves have a different
grid and cannot be directly subtracted. Frequency folding uses FFTW already linked
by the project, with single-thread transforms and deterministic coefficient sums.
It changes floating-point association, so equality to scalar goldens needs tester
evidence and reviewer acceptance under AGENTS.md, never a silent tolerance change.

Kernel `modeled_gflops` and `modeled_gbs` use the FLOP/byte counts printed by the
driver: complex multiply-add = 8 FLOP, CSR value/index = 24 bytes (64-bit index),
row pointer = 8 bytes, vector scalar = 16 bytes. They exclude setup and include
packing in timed GEMMs. Bytes are an algorithmic stream model, not measured DRAM
traffic (cache reuse/write allocation can change hardware traffic). For contraction
and dots they are lower-bound unique input/output streams; label them accordingly.

Item 6 production workloads (run on Aron; no input generation is timed):

```bash
export CONDA_PREFIX="$HOME/lsq-env"
export PATH="$CONDA_PREFIX/bin:$PATH"
export OMP_DYNAMIC=FALSE OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 KPM_SEED=12345
export OMP_PROC_BIND=close OMP_PLACES=cores
mkdir -p /tmp/stage8/inputs /tmp/stage8/results
export W2S_BIN="$HOME/wannier2sparse/build/wannier2sparse"
git -C "$HOME/wannier2sparse" rev-parse HEAD > /tmp/stage8/inputs/w2s-revision.txt
"$W2S_BIN" graphene 2000 2000 1 VX -p "$PWD/goldens/stage0/bases/graphene/model" -o /tmp/stage8/inputs/graphene
python bench/prepare_runs.py --repo "$PWD" --graphene-csr /tmp/stage8/inputs/graphene --bi-run "$BI_RUN" --out /tmp/stage8/inputs/runs
sha256sum /tmp/stage8/inputs/graphene/graphene.HAM.CSR /tmp/stage8/inputs/graphene/graphene.out "$BI_RUN" > /tmp/stage8/inputs/checksums.txt
cp -a bench /tmp/stage8/drivers
```

Require wannier2sparse revision `26cab4affff8266f8a76ec73f3201e6c063a1ba6` from
the brief; use the tester's existing binary and inputs, no fetching in the driver.
The graphene CSR has 8,000,000 orbitals and is generated from the committed two
orbital nearest-neighbor model. It runs through the CSR front end to avoid building
a clean eigen reference at that size. M=2048, R=1; the moment computation is timed
separately from input parsing and output. Haldane uses 48x48 Bloch blocks, M=256,
R=16 random vectors, explicit window, Greenwood on both uniform and FFT nodes.
Bi keeps the supplied 6b word/metric/reference/oracle settings and uses 2304 blocks,
M=1024, R=16 in the orbital Bloch representation. Supply absolute paths for any
6b model fields other than `dir`/`csr_dir`; those two are resolved by the preparer.
The 6b implementation is absent from this base checkout; Bi measurements require
the owning branch to be integrated by Fable first. Profile the short Bi run with
M=128, R=2 as well. Report reference/oracle costs separately and keep their settings
the same before/after; do not silently remove them to claim a kernel improvement.

For each built checkout, set `PYTHONPATH` to **that checkout's** package and build.
These commands write new artifacts per workload, profile, thread count and repeat:

```bash
revision=$(git rev-parse HEAD)
for profile in golden native; do
  export PYTHONPATH="$PWD/python:$PWD/build-$profile/lsquant_core"
  for threads in 24 48 96; do
    if [ "$threads" -eq 96 ]; then nodes=0,1; else nodes=0; fi
    for repeat in 1 2 3; do
      for workload in graphene haldane bi; do
        dest=/tmp/stage8/results/$revision-$profile-t$threads-r$repeat-$workload
        OMP_NUM_THREADS=$threads numactl --cpunodebind=$nodes --membind=$nodes python /tmp/stage8/drivers/benchmark_run.py /tmp/stage8/inputs/runs/$workload.json --revision "$revision" --out "$dest" > "$dest.log" 2>&1 &
        job=$!
        echo "$job" > "$dest.launcher.pid"
        wait "$job" || exit 1
      done
    done
  done
done
```

Before item 3, build and measure golden only. Record the active engine profile
from the benchmark, not just the build-directory name. To reconstruct before/after
evidence, the tester uses separate detached worktrees at `a6db913`, item 1, item 3,
item 4, item 5, and item 6 commits listed in the delivery. Copy the final `bench/`
drivers outside those worktrees as above; the production timing driver supports
the pre-stage-8 API. Keep the Bi integration revision identical in each comparison
and record it; any conflicts when combining branches belong to Fable's integration.
The optional kernel target and FFT comparison are available only after item 5.
Repeat item-1 profiling at every milestone, using separate output paths.

For one process per socket, run two independent copies of the same workload
concurrently, with disjoint output paths:

```bash
export PYTHONPATH="$PWD/python:$PWD/build-native/lsquant_core"
for socket in 0 1; do
  dest=/tmp/stage8/results/$revision-native-socket$socket-bi
  OMP_NUM_THREADS=48 numactl --cpunodebind=$socket --membind=$socket python /tmp/stage8/drivers/benchmark_run.py /tmp/stage8/inputs/runs/bi.json --revision "$revision" --out "$dest" > "$dest.log" 2>&1 &
  echo $! > "$dest.launcher.pid"
done
wait
```

Repeat the socket recipe for graphene and Haldane, three times under distinct
destinations. Record pair makespan and jobs/hour alongside individual latency.
For a single socket comparison, run these same two jobs sequentially first.
Inspect memory placement during a run with the PID in its `run.pid`:

```bash
numastat -p PID
cat /proc/PID/numa_maps
```

Record total wall time, moment-engine wall time, reconstruction wall time, CPU
seconds, peak RSS, thread occupancy, NUMA pages, per-profile output hashes, compiler
and dependency versions, flags, process affinity, host load and governor. The
driver calculates modeled GFLOP/s/GB/s from measured phase time and explicit
leading-work formulas. Rank 1 uses input CSR actions and dots (omits any added
rescaling diagonal and scalar updates); rank 2 uses dot-table work only (omits
Hamiltonian, metric and static-operator work), with a unique-stream byte lower
bound. These are comparable workload rates, **not measured hardware FLOP or DRAM
traffic**. Read the `model` field of each phase. Actual rates need architecture
counters with their event definitions archived. FFT actual FLOPs stay null: do
not divide a cubic table-contraction work estimate by FFT time and call it GFLOP/s.

Reference one-socket roofline uses the brief's approximate 48x35 GFLOP/s and
190 GB/s: `min(1680,190*modeled_intensity)` GFLOP/s. It is a hardware specification
reference, not a measured ceiling. Report 96-thread results as two-socket work,
even when quoting their ratio to this one-socket reference. If STREAM/compute
ceilings are measured on Aron, archive those separately instead of substituting
numbers from the older `run_graphene_scaling.sh` (which targets a different node).

Compare both contraction backends on an output rank-2 moment file (the tester
sets `MOMENT_FILE` to the exact produced pathname; do not pick an arbitrary glob):

```bash
OMP_NUM_THREADS=48 python bench/compare_contraction.py "$MOMENT_FILE" --formula greenwood --out /tmp/stage8/fft-haldane-t48
OMP_NUM_THREADS=24 python bench/compare_contraction.py "$MOMENT_FILE" --formula greenwood --out /tmp/stage8/fft-haldane-t24
OMP_NUM_THREADS=96 python bench/compare_contraction.py "$MOMENT_FILE" --formula greenwood --out /tmp/stage8/fft-haldane-t96
OMP_NUM_THREADS=48 python bench/compare_contraction.py "$HALL_MOMENT_FILE" --formula bastin --out /tmp/stage8/fft-hall-t48
```

For the FFT/scalar golden comparison, use committed G4/G5 moment files under the
golden engine and record first differing element, max absolute/relative difference
and whether bytes match; no tolerance is set by these benchmark drivers. The
default golden gate runs on the fat node (not Aron), by the tester:

```bash
KPM_SEED=12345 scripts/check_goldens.sh 7 build-golden > /tmp/stage8/gate7.log 2>&1
```

The tester owns the stage-8 gate extension and figure; these drivers do not edit
the harness, figures or goldens. Skip optional item 7 until the requested gcc
profiles have measurements; no oneAPI profile or dependency is added here.
