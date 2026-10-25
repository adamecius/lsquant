# Handoff 0008 (developer) — stage 8: performance on the KPM nodes (Aron class)

DRAFT v0 (2026-10-22). Roles: Codex (effort high) implements while Agi's quota is out, an
Opus tester measures on Aron, Fable reviews. Read `AGENTS.md` (sections 2, 3, 8, 10, 11,
12 binding), `handoff/0006b-siesta-spin-hall/review-1.md` item 5 (the first symptom), and
this brief completely. Branch `handoff/0008-dev` from `origin/cpp-unified` after the 6b
merge (Fable gives the commit). Nothing here changes a formula; every change is measured.

## 0. The machine this is for (user, 2026-10-22: "these are the nodes where most KPM runs")

`tcn-aronnode`: 2 x AMD EPYC 7552 (Zen 2 "Rome"), 48 cores per socket, no SMT, 96 cores,
2 NUMA nodes (cores 0-47 and 48-95), 2.2 GHz base (boost 3.2), AVX2 + FMA (no AVX-512),
L1d 32 KB and L2 512 KB per core, L3 16 MB per 4-core CCX (384 MB total), 2 TB RAM,
transparent huge pages `madvise`, governor `ondemand`, NUMA balancing on. Toolchain: gcc
9.4 (system), the conda env (`~/lsq-env`: Eigen 3.4, FFTW 3.3, GSL, OpenBLAS numpy), Intel
oneAPI 2025 (icpx, MKL) under `/opt/intel/oneapi`. Byte identity with the fat-node goldens
does not hold here (other gcc, Eigen, glibc), so this node is the place for a
performance build profile; the golden build profile stays as it is.

Symptom that opened the stage: the Bi spin Hall run of 6b (block basis, 2304 blocks of
72 x 72, `M = 1024`, 16 vectors) takes 71 minutes with `OMP_NUM_THREADS = 48` while the
process uses 2.6 cores (264 % CPU, 95 threads alive). The build carries no `-march` flag
(SSE2 code on an AVX2 machine).

## 1. Deliverables, in order, one commit each, each with a number

1. **Profile first.** `bench/profile_run.py`: runs one run file under `cProfile` (Python)
   and reports the top 20 cumulative entries and the fraction inside `_core`; plus the
   OpenMP occupancy (`/proc/<pid>/status` threads, `ps` CPU %, `perf stat` if available).
   Run on the 6b Bi spin Hall file at `M = 128`, 2 vectors, and on `runs/goldens/g5.json`.
   This is the "before" of every later item.
2. **The per-block Python loop is gone** (6b review-1 item 5): one C++ call per moment
   step on the block basis (`BlockDiagonalOperator` batched GEMM over the live vectors,
   the per-block exact `S^-1` / `S^-1/2` as block-diagonal operators applied in C++); the
   eigen reference of a run is optional (`"reference": false`), not always computed.
3. **Build profiles.** CMake option `LSQUANT_BUILD_PROFILE = golden | native`: `golden` is
   the present flag set (byte identity with the stage-0 goldens on the fat node); `native`
   adds `-O3 -march=native -mtune=native -fno-math-errno` (no `-ffast-math`, no
   `-funsafe-math-optimizations`: reproducibility within a machine and between thread
   counts stays a rule), reported in `run.resolved.json` and in every moment-file header
   so a golden can never be compared with the wrong profile by accident.
4. **Threads and NUMA.** Defaults set by the Python driver when the environment does not:
   `OMP_PROC_BIND=close`, `OMP_PLACES=cores`, `OMP_NUM_THREADS` = physical cores of the
   NUMA nodes the process is allowed on; first-touch initialisation of every large array
   (CSR arrays, Chebyshev vector batches, block arrays) inside the parallel region that
   will use it; the block loop partitioned so that each thread's blocks stay on its NUMA
   node; a documented recipe for "one process per socket" with `numactl --cpunodebind
   --membind`. Measure 24 / 48 / 96 threads and both recipes.
5. **Kernels.** (a) `SparseMatrix::Multiply`: the parallel threshold (removed in stage 1)
   is back as a measured constant, static schedule, row partition balanced by nnz, and the
   batched form (several live vectors per sweep: the CSR is read once per batch, not once
   per vector); (b) the batched dot products of the moment accumulator as one GEMM
   (`Eigen` on the native profile uses AVX2 FMA); (c) the rank-2 contraction as GEMM with
   the FFT backend of the original plan (`Kubo_solver_FFT` contraction) as an option,
   measured against the GEMM; (d) `BlockDiagonalOperator`: blocks of `W x W` complex
   (72 x 72 = 83 KB) times the batch of live vectors as one GEMM per block, blocks
   grouped per thread to stay in L2; (e) `ChebyshevVectors` batch size chosen so that the
   working set of one batch fits the L2 of the cores that sweep it (document the formula).
6. **Benchmarks** (`bench/`, run by the tester on Aron, numbers in the report): graphene
   real-space 2000 x 2000 cells DOS `M = 2048` (SpMV bound), Haldane rank-2 `M = 256` 48 x 48
   (contraction bound), Bi block basis spin Hall as in 6b (block GEMM bound); each reported
   as wall time, GFLOP/s, GB/s, and against the roofline of one EPYC 7552 socket
   (about 35 GFLOP/s per core in double with FMA, about 190 GB/s per socket). Before and
   after each item.
7. **Optional, report only:** the same binaries with icpx + MKL from oneAPI
   (`LSQUANT_BUILD_PROFILE = oneapi`), if it builds; no default change.

## 2. Restrictions

Every number in the report comes from the tester's runs on Aron, not from estimates. No
change to any formula, any file format or any golden; the golden profile build on the fat
node must still reproduce every golden byte for byte (the tester reruns `check_goldens.sh
7` there at the closure). Thread-count byte identity on one machine is kept (static
schedules, deterministic reductions). No new dependency. The engine stays basis-agnostic.

## 3. Report

`handoff/0008-performance/report.md`: the profiles before and after each item, the
benchmark table, the NUMA recipe, the flags of each build profile, deviations, open
questions. No "done", no gate claims, no push.
