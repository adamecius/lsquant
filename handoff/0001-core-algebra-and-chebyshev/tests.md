# Handoff 0001 tests (Codex role, Opus subagent) — stage 1 gate on TCN, final

**Result: GATE PASS (stage 1).**

The state tested is `handoff/0001-codex`: `cpp-unified` (90b0bae) + `origin/handoff/0001-agi`,
merged twice:
- first in merge commit a3a3d8e (review 1);
- then in a097a58, which brings in the developer's review-2 fixes, 34173e2.

Both merges had no conflicts. The run was on TCN (`fatnode-tacn`, checkout `~/lsquant-unified`)
with g++ 11.4.0, OpenMP 4.5, Python 3.10, numpy 2.2.6 and matplotlib 3.10.9. Nothing under
`lsquant/`, `tests/dev/` or the legacy tree was modified by the tester.

## Summary

| Item | Result |
|---|---|
| Fresh build | exit 0. `stage1_moments` and `stage_1_check` are regular files in `build/lsquant_core/` (no symlinks). 0 warnings under `lsquant/` and `tests/dev/` |
| `stage_1_check` at 1 and 8 threads | dos 2.47e-15, kubo 1.14e-13, 1 vs 8 bit-identical YES (the parallel kernel is now always active) |
| G1, G2 at T = 1, 8, 64 | IDENTICAL |
| G2 with `--batch 16` at T = 8 | IDENTICAL to the unbatched G2 |
| Random-phase determinism | IDENTICAL (4 files) |
| `grep -rn "basis\|std::cout\|getenv" lsquant/core` | empty |
| Legacy `ctest` | 31/31 |
| Figure | `figures/out/stage_1/fig_stage_1.png`; all residuals exactly 0 |

## 1. Build

```
cd ~/lsquant-unified && rm -rf build
cmake -S . -B build -DLSQUANT_BUILD_NEW_CORE=ON -DLSQUANT_BUILD_INPUT_TOOLS=ON > cmake.log 2>&1
cmake --build build -j 32 > build.log 2>&1        # BUILD_EXIT=0
```

`build.log` has 101 lines containing `warning`, all in the legacy tree:

| File | Warnings |
|---|---|
| `include/chebyshev_moments.hpp` | 62 |
| `include/units.hpp` | 32 |
| `src/chebyshev_momentsLocal.cpp` | 7 |

**Under `lsquant/` and `tests/dev/` there are none** (`-Wall -Wextra`).

`ls -la build/lsquant_core`:

```
drwxr-xr-x 5 jgarcia Theoretical   4096 Oct 20 13:57 CMakeFiles
-rw-r--r-- 1 jgarcia Theoretical   1263 Oct 20 13:57 cmake_install.cmake
-rw-r--r-- 1 jgarcia Theoretical  25942 Oct 20 13:57 Makefile
-rwxr-xr-x 1 jgarcia Theoretical 279888 Oct 20 13:58 stage_1_check
-rwxr-xr-x 1 jgarcia Theoretical 175144 Oct 20 13:58 stage1_moments
```

`build/stage1_moments` and `build/stage_1_check` no longer exist.

## 2. `stage_1_check` (`figures/out/stage_1/stage_1_check.log`)

```
### OMP_NUM_THREADS=1 build/lsquant_core/stage_1_check
dos max_abs_error: 2.47025e-15
kubo max_abs_error: 1.13687e-13
OMP_NUM_THREADS=1 vs 8 bit-identical: YES
exit=0
### OMP_NUM_THREADS=8 build/lsquant_core/stage_1_check
dos max_abs_error: 2.47025e-15
kubo max_abs_error: 1.13687e-13
OMP_NUM_THREADS=1 vs 8 bit-identical: YES
exit=0
```

Both errors are below 1e-11. With the threshold removed (F3), the 1 vs 8 comparison now
goes through the row-parallel kernel.

## 3. The gate (`figures/out/stage_1/check_goldens_stage1.log`)

```
$ OMP_NUM_THREADS=16 bash scripts/check_goldens.sh 1 build 2>&1 | tee figures/out/stage_1/check_goldens_stage1.log
== stage 1: /home/ICN2/jgarcia/lsquant-unified/build/lsquant_core/stage1_moments vs goldens/stage0/chain1d (KPM_SEED=12345)
IDENTICAL [T=1] chain1d/SpectralOp1chain1dKPM_M64_stateexact.chebmom1D (G1)
IDENTICAL [T=1] chain1d/NonEqOpVX-VXchain1dKPM_M64x64_stateexact.chebmom2D (G2)
IDENTICAL [T=8] chain1d/SpectralOp1chain1dKPM_M64_stateexact.chebmom1D (G1)
IDENTICAL [T=8] chain1d/NonEqOpVX-VXchain1dKPM_M64x64_stateexact.chebmom2D (G2)
IDENTICAL [T=64] chain1d/SpectralOp1chain1dKPM_M64_stateexact.chebmom1D (G1)
IDENTICAL [T=64] chain1d/NonEqOpVX-VXchain1dKPM_M64x64_stateexact.chebmom2D (G2)
IDENTICAL [T=8, batch 16] chain1d/NonEqOpVX-VXchain1dKPM_M64x64_stateexact.chebmom2D (G2, batched)
IDENTICAL [random, run2 vs run1] NonEqOpVX-VXchain1dKPM_M64x64_statedefault.chebmom
IDENTICAL [random, run2 vs run1] NonEqOpVX-VXchain1dKPM_M64x64_statedefault.chebmom2D
IDENTICAL [random, run2 vs run1] SpectralOp1chain1dKPM_M64_statedefault.chebmom
IDENTICAL [random, run2 vs run1] SpectralOp1chain1dKPM_M64_statedefault.chebmom1D
GREP      grep -rn "basis\|std::cout\|getenv" lsquant/core: empty
== 11 files compared; failures: 0
GATE PASS (stage 1)
```

Exit code 0. Every `cmp` was silent, so no section-3 localisation was needed.

**Batch 16.** The new core with `--batch 16` (M = 64, four blocks) at T = 8 is byte-identical
to the unbatched G2. That is the correct reference: each element gets one dot product per
state, in state order, whatever the batch size. The legacy code with `BATCH_SIZE=16` is wrong:
block k is shifted by k orders, because `Vectors::IterateAll` does not carry the recursion
across batches. That output is kept as evidence only, with its analysis, under
`goldens/stage0/chain1d/legacy_defects/batch16/` (its README.md has the details), and no
harness reads it. The new core does not port the defect.

The harness gained this batch-16 comparison in commit 151caba. The executable lookup is
unchanged (`build_dir/lsquant_core`, then `build_dir/lsquant`), with no fallback to `build/`.

## 4. Legacy `ctest` with the new core on (`figures/out/stage_1/ctest.log`)

```
OMP_NUM_THREADS=16 KPM_SEED=12345 ctest --test-dir build --output-on-failure
100% tests passed, 0 tests failed out of 31
Total Test time (real) =  91.56 sec
```

Known flaw, still present and not fixed (legacy tree): `test/graphene_reproduce.sh` line 51
takes the largest decimal number anywhere in the solver log. The log includes the spdlog
timestamp `[hh:mm:SS.mmm]`. A run that ends at second 50 or later of a minute therefore fails
`graphene_regression` with `|E|max` equal to that second. The real Gershgorin bound is
7.51834. This happened in two of three runs during review 1 (58.682, 53.694). It does not
depend on the new core.

## 5. Figure

`figures/stage_1.py`, run from the fresh build:

```
MPLBACKEND=Agg python3 figures/stage_1.py --threads 1 8 64
```

It exited 0 after 594 s; the log is `figures/out/stage_1/stage_1_figure.log`. The script:
- runs `build/lsquant_core/stage1_moments` itself at each thread count;
- byte-compares the outputs with `goldens/stage0/chain1d`;
- reads the new tensors from the `.chebmom` files.

The residuals are still exactly zero: G1 0, the G2 block 0, and |mu/plateau − Mref| 0.

- Figure: `figures/out/stage_1/fig_stage_1.png`
- Caption: `figures/out/stage_1/fig_stage_1.caption.txt`

> FIG. 1. Stage 1 of the unified C++ core (handoff 0001): Chebyshev moments of the periodic tight-binding chain chain1d, N = 2048 sites, spectral window [-2, 2], exact trace over all 2048 site states, M = 64 moments, KPM_SEED = 12345, computed by build/stage1_moments on goldens/stage0/chain1d/operators at OMP_NUM_THREADS = 1, 8, 64. (a) Rank-1 moments mu_m (G1, word dos): legacy golden (open circles), new core (dashed line) and the analytic reference delta_m0 (triangles); right axis: |new - golden| (dash-dot), max 0.0e+00. (b) |new - golden| over the m, n <= 10 block of the rank-2 moments mu_mn (G2, word kubo, KPM weights and Hermitian symmetrisation applied), max 0.0e+00. (c) The closed-form block B.15, Mref(m, n) of test/chain1d_analytic_test.cpp times the diagonal plateau mean(mu_mm, m = 2..9) = 4 of the new tensor; max |Re mu_mn / plateau - Mref| over m, n < 10 is 0.0e+00, max |Im mu_mn| over the block is 0.0e+00. Byte comparison of the legacy-layout files with goldens/stage0/chain1d (cmp criterion): T=1 SpectralOp1chain1dKPM_M64_stateexact.chebmom1D: IDENTICAL; T=1 NonEqOpVX-VXchain1dKPM_M64x64_stateexact.chebmom2D: IDENTICAL; T=8 SpectralOp1chain1dKPM_M64_stateexact.chebmom1D: IDENTICAL; T=8 NonEqOpVX-VXchain1dKPM_M64x64_stateexact.chebmom2D: IDENTICAL; T=64 SpectralOp1chain1dKPM_M64_stateexact.chebmom1D: IDENTICAL; T=64 NonEqOpVX-VXchain1dKPM_M64x64_stateexact.chebmom2D: IDENTICAL. The new-format .chebmom tensors are bit-identical across the thread counts and agree with the 15-digit legacy-layout files to 0.0e+00. All files IDENTICAL.

## 6. Smells and findings

| # | Finding | Status |
|---|---|---|
| F1 | Executables landed in `build/`, not `build/lsquant_core/`; the developer's earlier runs relied on hand-made symlinks | **Resolved** (review 2, 34173e2): per-target `RUNTIME_OUTPUT_DIRECTORY`. A fresh build puts both in `build/lsquant_core/`, and the official gate ran from there |
| F2 | `apps/stage1_moments.cpp` read `KPM_SEED` and `BATCH_SIZE` from the environment | **Resolved**: both `getenv` calls are removed; `--seed` and `--batch` are the only sources |
| F3 | The parallel SpMV was never exercised (`kMinParallelNnz = 30000`) | **Resolved**: the threshold is removed and the row loop is always `omp parallel for schedule(static)`. G1, G2 at 1/8/64 threads and `stage_1_check` 1 vs 8 now test real parallel runs and are bit-identical |
| F4 | Batching was not compared with anything | **Resolved by a new gate check**: `--batch 16` against G2, IDENTICAL. The legacy batched path is defective (evidence in `goldens/stage0/chain1d/legacy_defects/batch16/`) |
| F5 | `ReadCsrText` Hermitian-folds every operator file | Accepted by Fable (legacy reader convention, documented) |
| F6 | Four hand-written rank paths in `MomentEngine::Run` | Accepted by Fable (a generic loop comes with `ApplySeries` in stage 4) |
| F7 | Real-weight special case in `Trace`, needed to reproduce the legacy signed zeros | Accepted by Fable |
| F8 | `GershgorinBounds` is symmetric about zero | Accepted by Fable |
| F9 (new, for information) | Cost of the always-parallel SpMV on small matrices; numbers in section 7 | Open, for Fable's information; no correctness impact |

Also checked and clean:
- no `TODO` or `FIXME` anywhere under `lsquant/` or `tests/dev/`;
- no `getenv` anywhere under `lsquant/` now (`apps/` included);
- no `std::cout` or `printf` under `core/` or `io/`;
- file streams only in `io/`;
- the out-of-scope classes are absent;
- Eigen and OpenMP are the only link dependencies of `lsquant_core`.

## 7. F9: wall time per thread count

Wall time of `build/lsquant_core/stage1_moments` with the gate arguments (exact trace,
N = 2048, M = 64, `--seed 12345`), one run each. The machine is shared: the load average
was 22 to 46 during these runs.

| Word | T = 1 | T = 8 | T = 64 |
|---|---|---|---|
| dos | 3.3 s | 2.6 s | 124.0 s |
| kubo | 29.5 s | 31.3 s | 400.4 s |
| legacy `inline_compute-kpm-nonEqOp` (kubo, threshold kept) | | | 29.6 s |

The official gate took more than 6 minutes for its T = 64 kubo run, and the whole figure run
took 594 s. On a 2048-row matrix,
every matvec now opens a parallel region, and at 64 threads the fork/join cost exceeds the
work. The legacy kept the 30000-nnz threshold for this reason. The threshold was removed so
that thread-count determinism is actually tested, which was the right call for the gate. For
stage 8 (performance), a threshold that a test can override, or a gate case large enough to
cross it, would keep both properties. No action requested in this stage.

## 8. Tester commits on `handoff/0001-codex`

- ed0e47c: figures, the stage-1 figure script, the figure, the caption and the run logs (review 1)
- d9becf6: tests report (review 1)
- 151caba: harness, batch-16 kubo run compared with G2
- 738da80: goldens, legacy batched defect recorded as evidence only
- a097a58: merge of `origin/handoff/0001-agi` at 34173e2 (review 2)
- this commit: final logs, figure and this report

No golden reference, seed, tolerance or harness criterion was changed.
