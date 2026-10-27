# Handoff 0003 — tester report (stage 3 gate: G4 and G5 from Python; tutorial 03)

Tester: Codex role (Opus subagent driven by Fable). Machine: TCN (`fatnode-tacn`),
checkout `~/lsquant-unified`, branch `handoff/0003-codex` cut from `cpp-unified`
(`e1e6b28`). Nothing under `lsquant/`, `python/lsquant/`, `tests/dev/`, the legacy tree,
`goldens/`, a seed or a tolerance was edited. Every run used an explicit
`OMP_NUM_THREADS` (16, and 1 where the gate says so; the stage-1 part of the harness
keeps its own 1, 8 and 64).

The developer branch was tested three times, each after a merge, a fresh build and a
rerun of the self-check, the gate and the figure script:

| Run | `origin/handoff/0003-agi` | merge into `handoff/0003-codex` | gate | tutorial 03 | logs in `figures/out/stage_3/` |
|---|---|---|---|---|---|
| 1 | `c01b026` | fast-forward | PASS, 41/41 | 1 of 4 IDENTICAL | `*_run1_c01b026.log` |
| 2 | `0774db7` (review 1: keyword-only `kubo_plan`, `word` in `dos_plan`) | fast-forward | PASS, 41/41 | 1 of 4 IDENTICAL | `*_run2_0774db7.log`, `fig_stage_3_run2_0774db7.caption.txt`, `merge_run2.log` |
| **3 (final)** | **`cac5da6`** (review 2: `tutorial03` and `stage3_kubo` form 2 `--formula` round-trip the moments through `WriteLegacyChebmom2D`/`ReadChebmom`; `tutorial03 --operators` defaults to `goldens/stage0/haldane/operators`) | merge commit `3dfa53c`, no conflict (`merge_run3.log`) | **PASS, 41/41** | **4 of 4 IDENTICAL** | `build.log`, `stage_3_check.log`, `check_goldens_stage3.log`, `figure_run.log`, `fig_stage_3.caption.txt` |

**Result (final run, `cac5da6`): GATE PASS (stage 3), 41 of 41 files IDENTICAL at
OMP_NUM_THREADS 16 and 1; both greps empty. Tutorial 03 from Python: 4 of 4 files
IDENTICAL. C = +1.017016743089, equal to the golden's C. `figures/stage_3.py` exits 0.**
In runs 1 and 2, tutorial 03 was 1 of 4 IDENTICAL. The localisation that led to the
review-2 change is kept in section 6.

## 1. Build (fresh, every run)

```
$ rm -rf build
$ cmake -S . -B build -DLSQUANT_BUILD_NEW_CORE=ON -DLSQUANT_BUILD_PYTHON=ON -Dpybind11_DIR="$(python3 -c 'import pybind11; print(pybind11.get_cmake_dir())')"
-- Found pybind11: /home/ICN2/jgarcia/.local/lib/python3.10/site-packages/pybind11/include (found version "3.1.0")
-- Build files have been written to: /home/ICN2/jgarcia/lsquant-unified/build
cmake exit: 0
$ cmake --build build -j 32
[100%] Built target stage_1_check
[100%] Built target _core
build exit: 0
```

Full logs: `figures/out/stage_3/build.log` (final run), `build_run2_0774db7.log` and
`build_run1_c01b026.log`. `compat/legacy_kernels.cpp` compiles into `_core` with the
same `-O3 -Wall -Wextra` flags (`-flto=auto` on `_core`). Warnings: **0 under
`lsquant/`** (core, io, api, compat, apps) in all three runs. The 101 warnings in each
log are all in the legacy tree, as in stages 1 and 2: 62 in
`include/chebyshev_moments.hpp`, 32 in `include/units.hpp`, 7 in
`src/chebyshev_momentsLocal.cpp`.

## 2. Developer self-check (final run, `cac5da6`)

```
$ PYTHONPATH=build/lsquant_core:python OMP_NUM_THREADS=16 python3 tests/dev/stage_3_check.py
(0) SparseMatrix.from_csr (scipy) dos array_equal vs ReadCsrText: True
(a) Haldane VX-VY array_equal vs G5: True
(a) Haldane VX-VX array_equal vs G5: True
(b) Graphene Gershgorin BandWidth = 15.0366847652085, array_equal vs G4: True
(c) GreenRetarded.table vs compat max_abs_diff           = 1.148853e-14
(c) GreenRetardedDerivative.table vs compat max_abs_diff = 2.702165e-13 (max_rel_diff = 7.173405e-15)
(d) G5 Bastin (10 meV)           differing lines: 0 (total lines: 3839)
(d) G5 Sea (10 meV)              differing lines: 0 (total lines: 3839)
(d) G5 Surf (10 meV)             differing lines: 0 (total lines: 3839)
(d) G5 Greenwood (10 meV)        differing lines: 0 (total lines: 127)
(d) G4 Bastin (64 meV)           differing lines: 0 (total lines: 1919)
(d) G4 Sea (10 meV)              differing lines: 0 (total lines: 1919)
(d) G4 Surf (10 meV)             differing lines: 0 (total lines: 1919)
(d) G4 Bastin FFTgrid (10 meV)   differing lines: 0 (total lines: 63)
(d) G4 Greenwood FFTgrid (10 meV) differing lines: 0 (total lines: 63)
wall 26.66 s
exit: 0
```

Exit 0. Runs 1 and 2 printed the same lines (only the wall time differs: 27.42 s and
27.12 s). It matches the developer's report.

## 3. Item 0 of the brief (`SparseMatrix.from_csr`), tester check

Done with a scratch script that is not committed (the brief asks for this check by name;
it is not a test suite). It parses `goldens/stage0/chain1d/operators/chain1d.HAM.CSR`
into a `scipy.sparse.csr_matrix` (scipy gives `int32` indices), passes the arrays
converted to `int64`/`complex128` to `from_csr`, and then tries the wrong buffers.
Output at `0774db7` (review 2 does not touch the binding):

```
scipy csr: shape (2048, 2048) nnz 4096 dtypes int32 int32 complex128
int64/complex128 from_csr + hermitian_fold: dos moments array_equal vs ReadCsrText: True
indptr int32                              : TypeError: SparseMatrix.from_csr: argument 'row_pointers' must be a 1D C-contiguous int64 numpy.ndarray
indices int32                             : TypeError: SparseMatrix.from_csr: argument 'columns' must be a 1D C-contiguous int64 numpy.ndarray
scipy native int32 indptr/indices         : TypeError: SparseMatrix.from_csr: argument 'row_pointers' must be a 1D C-contiguous int64 numpy.ndarray
indices non-contiguous int64 (stride 2)   : TypeError: SparseMatrix.from_csr: argument 'columns' must be a 1D C-contiguous int64 numpy.ndarray
values complex64                          : TypeError: SparseMatrix.from_csr: argument 'values' must be a 1D C-contiguous complex128 numpy.ndarray
values float64                            : TypeError: SparseMatrix.from_csr: argument 'values' must be a 1D C-contiguous complex128 numpy.ndarray
indptr python list                        : TypeError: SparseMatrix.from_csr: argument 'row_pointers' must be a 1D C-contiguous int64 numpy.ndarray
```

Source check: the binding validates each buffer (`ValidateCsrBuffer1D`, forcecast
off) and passes the raw pointers to the new pointer overload
`SparseMatrix::FromCsr(dim, rpSize, rpPtr, colSize, colPtr, valSize, valPtr)`. That
overload assigns an `Eigen::Map<const SparseMatrix>` into the owned matrix, which is
the single accepted copy. There is no `std::vector` staging on the Python path. The
by-value `std::vector` overload remains for C++ callers and delegates to the pointer
overload.

## 4. Gate (`scripts/check_goldens.sh 3 build`)

Harness extension (this handoff, `scripts/check_goldens.sh`): stage 3 runs stage 1 and
stage 2 in full (their comparisons and greps). It then runs
`python3 -m lsquant.apps.stage3_kubo` at `OMP_NUM_THREADS=16` and again at `1`, with
`PYTHONPATH=build/lsquant_core:python`:

1. G5 moments: Haldane `--left VX --right VY` and `--right VX`, `--window -4.5 4.5
   --trace exact --states 512 --moments 128`; `cmp` both `.chebmom2D` against
   `goldens/stage0/haldane/`.
2. G5 curves from the golden moment files at 10 meV: `bastin`, `sea` and `surf` on
   `VX-VY`, `greenwood` on `VX-VX`.
3. G4 moments: graphene `VX VX`, `--window gershgorin --trace random --seed 12345
   --states 1 --moments 64`.
4. G4 curves from the golden moment file: `bastin` at 64 meV, `sea` and `surf` at
   10 meV (`sea_surf/`), `greenwood_fft` and `bastin_fft` at 10 meV (`fft_grid/`).
5. `grep -rn "compat/" lsquant/core lsquant/io` must be empty. The stage-2 code-line
   grep for `delta|DOS|Kubo|basis|std::cout|getenv` on `lsquant/core` and `lsquant/io`
   is part of the stage-2 rerun.

Command: `OMP_NUM_THREADS=16 bash scripts/check_goldens.sh 3 build 2>&1 | tee
figures/out/stage_3/check_goldens_stage3.log`. (In the final run, the wrapper appended
the exit status as the last line of the log.)

### 4.1 Final run (`cac5da6`), verbatim — exit 0

```
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
== stage 2: python3 -m lsquant.apps.stage2_dos, PYTHONPATH=/home/ICN2/jgarcia/lsquant-unified/build/lsquant_core:python, vs goldens/stage0/chain1d (KPM_SEED=12345)
          [T=1] (a) printed: meanSpectralOp1chain1dKPM_M64_stateexactJACKSON.dat
IDENTICAL [T=1] (a) --moments G1 -> chain1d/meanSpectralOp1chain1dKPM_M64_stateexactJACKSON.dat (G3)
          [T=1] (b) printed: SpectralOp1chain1dKPM_M64_stateexact.chebmom1D
          [T=1] (b) printed: meanSpectralOp1chain1dKPM_M64_stateexactJACKSON.dat
IDENTICAL [T=1] (b) --operators -> chain1d/SpectralOp1chain1dKPM_M64_stateexact.chebmom1D (G1)
IDENTICAL [T=1] (b) --operators -> chain1d/meanSpectralOp1chain1dKPM_M64_stateexactJACKSON.dat (G3)
          [T=16] (a) printed: meanSpectralOp1chain1dKPM_M64_stateexactJACKSON.dat
IDENTICAL [T=16] (a) --moments G1 -> chain1d/meanSpectralOp1chain1dKPM_M64_stateexactJACKSON.dat (G3)
          [T=16] (b) printed: SpectralOp1chain1dKPM_M64_stateexact.chebmom1D
          [T=16] (b) printed: meanSpectralOp1chain1dKPM_M64_stateexactJACKSON.dat
IDENTICAL [T=16] (b) --operators -> chain1d/SpectralOp1chain1dKPM_M64_stateexact.chebmom1D (G1)
IDENTICAL [T=16] (b) --operators -> chain1d/meanSpectralOp1chain1dKPM_M64_stateexactJACKSON.dat (G3)
GREP      grep -rn "delta\|DOS\|Kubo\|basis\|std::cout\|getenv" lsquant/core lsquant/io, code lines: empty
          excluded comment line: lsquant/core/chebyshev/chebyshev_moments.hpp:27: * recording whether the `(2 - delta_{m0})` KPM weights (`kpmWeighted`) and Hermitian
          excluded comment line: lsquant/core/chebyshev/chebyshev_moments.hpp:32: *      Never apply KPM `(2 - delta_{m0})` weights or symmetrisation twice; check
          excluded comment line: lsquant/core/chebyshev/chebyshev_moments.hpp:158: * @brief Apply the legacy rank-2 KPM `(2 - delta_{m0})` weights and optional Hermitian symmetrisation.
== stage 3: python3 -m lsquant.apps.stage3_kubo, PYTHONPATH=/home/ICN2/jgarcia/lsquant-unified/build/lsquant_core:python, vs goldens/stage0/{haldane,graphene}
          [T=16] G5 VX-VY moments printed: NonEqOpVX-VYhaldaneKPM_M128x128_stateexact.chebmom2D
          [T=16] G5 VX-VY moments printed: NonEqOpVX-VYhaldaneKPM_M128x128_stateexact.chebmom
          [T=16] G5 VX-VX moments printed: NonEqOpVX-VXhaldaneKPM_M128x128_stateexact.chebmom2D
          [T=16] G5 VX-VX moments printed: NonEqOpVX-VXhaldaneKPM_M128x128_stateexact.chebmom
IDENTICAL [T=16] haldane/NonEqOpVX-VYhaldaneKPM_M128x128_stateexact.chebmom2D (G5 moments)
IDENTICAL [T=16] haldane/NonEqOpVX-VXhaldaneKPM_M128x128_stateexact.chebmom2D (G5 moments)
          [T=16] G5 bastin printed: KuboBastin_NonEqOpVX-VYhaldaneKPM_M128x128_stateexactJACKSON.dat
          [T=16] G5 sea printed: KuboBastinSea_NonEqOpVX-VYhaldaneKPM_M128x128_stateexactJACKSON.dat
          [T=16] G5 surf printed: KuboBastinSurf_NonEqOpVX-VYhaldaneKPM_M128x128_stateexactJACKSON.dat
          [T=16] G5 greenwood printed: KuboGreenWood_NonEqOpVX-VXhaldaneKPM_M128x128_stateexactJACKSON.dat
IDENTICAL [T=16] haldane/KuboBastin_NonEqOpVX-VYhaldaneKPM_M128x128_stateexactJACKSON.dat (G5 curve)
IDENTICAL [T=16] haldane/KuboBastinSea_NonEqOpVX-VYhaldaneKPM_M128x128_stateexactJACKSON.dat (G5 curve)
IDENTICAL [T=16] haldane/KuboBastinSurf_NonEqOpVX-VYhaldaneKPM_M128x128_stateexactJACKSON.dat (G5 curve)
IDENTICAL [T=16] haldane/KuboGreenWood_NonEqOpVX-VXhaldaneKPM_M128x128_stateexactJACKSON.dat (G5 curve)
          [T=16] G4 VX-VX moments printed: NonEqOpVX-VXgraphene_goldenKPM_M64x64_statedefault.chebmom2D
          [T=16] G4 VX-VX moments printed: NonEqOpVX-VXgraphene_goldenKPM_M64x64_statedefault.chebmom
IDENTICAL [T=16] graphene/NonEqOpVX-VXgraphene_goldenKPM_M64x64_statedefault.chebmom2D (G4 moments)
          [T=16] G4 bastin 64 meV printed: KuboBastin_NonEqOpVX-VXgraphene_goldenKPM_M64x64_statedefaultJACKSON.dat
          [T=16] G4 sea 10 meV printed: KuboBastinSea_NonEqOpVX-VXgraphene_goldenKPM_M64x64_statedefaultJACKSON.dat
          [T=16] G4 surf 10 meV printed: KuboBastinSurf_NonEqOpVX-VXgraphene_goldenKPM_M64x64_statedefaultJACKSON.dat
          [T=16] G4 greenwood_fft 10 meV printed: KuboGreenWood_NonEqOpVX-VXgraphene_goldenKPM_M64x64_statedefaultJACKSON_FFTgrid.dat
          [T=16] G4 bastin_fft 10 meV printed: KuboBastin_NonEqOpVX-VXgraphene_goldenKPM_M64x64_statedefaultJACKSON_FFTgrid.dat
IDENTICAL [T=16] graphene/KuboBastin_NonEqOpVX-VXgraphene_goldenKPM_M64x64_statedefaultJACKSON.dat (G4 curve, 64 meV)
IDENTICAL [T=16] graphene/sea_surf/KuboBastinSea_NonEqOpVX-VXgraphene_goldenKPM_M64x64_statedefaultJACKSON.dat (G4 curve, 10 meV)
IDENTICAL [T=16] graphene/sea_surf/KuboBastinSurf_NonEqOpVX-VXgraphene_goldenKPM_M64x64_statedefaultJACKSON.dat (G4 curve, 10 meV)
IDENTICAL [T=16] graphene/fft_grid/KuboGreenWood_NonEqOpVX-VXgraphene_goldenKPM_M64x64_statedefaultJACKSON_FFTgrid.dat (G4 curve, 10 meV)
IDENTICAL [T=16] graphene/fft_grid/KuboBastin_NonEqOpVX-VXgraphene_goldenKPM_M64x64_statedefaultJACKSON_FFTgrid.dat (G4 curve, 10 meV)
          [T=1] G5 VX-VY moments printed: NonEqOpVX-VYhaldaneKPM_M128x128_stateexact.chebmom2D
          [T=1] G5 VX-VY moments printed: NonEqOpVX-VYhaldaneKPM_M128x128_stateexact.chebmom
          [T=1] G5 VX-VX moments printed: NonEqOpVX-VXhaldaneKPM_M128x128_stateexact.chebmom2D
          [T=1] G5 VX-VX moments printed: NonEqOpVX-VXhaldaneKPM_M128x128_stateexact.chebmom
IDENTICAL [T=1] haldane/NonEqOpVX-VYhaldaneKPM_M128x128_stateexact.chebmom2D (G5 moments)
IDENTICAL [T=1] haldane/NonEqOpVX-VXhaldaneKPM_M128x128_stateexact.chebmom2D (G5 moments)
          [T=1] G5 bastin printed: KuboBastin_NonEqOpVX-VYhaldaneKPM_M128x128_stateexactJACKSON.dat
          [T=1] G5 sea printed: KuboBastinSea_NonEqOpVX-VYhaldaneKPM_M128x128_stateexactJACKSON.dat
          [T=1] G5 surf printed: KuboBastinSurf_NonEqOpVX-VYhaldaneKPM_M128x128_stateexactJACKSON.dat
          [T=1] G5 greenwood printed: KuboGreenWood_NonEqOpVX-VXhaldaneKPM_M128x128_stateexactJACKSON.dat
IDENTICAL [T=1] haldane/KuboBastin_NonEqOpVX-VYhaldaneKPM_M128x128_stateexactJACKSON.dat (G5 curve)
IDENTICAL [T=1] haldane/KuboBastinSea_NonEqOpVX-VYhaldaneKPM_M128x128_stateexactJACKSON.dat (G5 curve)
IDENTICAL [T=1] haldane/KuboBastinSurf_NonEqOpVX-VYhaldaneKPM_M128x128_stateexactJACKSON.dat (G5 curve)
IDENTICAL [T=1] haldane/KuboGreenWood_NonEqOpVX-VXhaldaneKPM_M128x128_stateexactJACKSON.dat (G5 curve)
          [T=1] G4 VX-VX moments printed: NonEqOpVX-VXgraphene_goldenKPM_M64x64_statedefault.chebmom2D
          [T=1] G4 VX-VX moments printed: NonEqOpVX-VXgraphene_goldenKPM_M64x64_statedefault.chebmom
IDENTICAL [T=1] graphene/NonEqOpVX-VXgraphene_goldenKPM_M64x64_statedefault.chebmom2D (G4 moments)
          [T=1] G4 bastin 64 meV printed: KuboBastin_NonEqOpVX-VXgraphene_goldenKPM_M64x64_statedefaultJACKSON.dat
          [T=1] G4 sea 10 meV printed: KuboBastinSea_NonEqOpVX-VXgraphene_goldenKPM_M64x64_statedefaultJACKSON.dat
          [T=1] G4 surf 10 meV printed: KuboBastinSurf_NonEqOpVX-VXgraphene_goldenKPM_M64x64_statedefaultJACKSON.dat
          [T=1] G4 greenwood_fft 10 meV printed: KuboGreenWood_NonEqOpVX-VXgraphene_goldenKPM_M64x64_statedefaultJACKSON_FFTgrid.dat
          [T=1] G4 bastin_fft 10 meV printed: KuboBastin_NonEqOpVX-VXgraphene_goldenKPM_M64x64_statedefaultJACKSON_FFTgrid.dat
IDENTICAL [T=1] graphene/KuboBastin_NonEqOpVX-VXgraphene_goldenKPM_M64x64_statedefaultJACKSON.dat (G4 curve, 64 meV)
IDENTICAL [T=1] graphene/sea_surf/KuboBastinSea_NonEqOpVX-VXgraphene_goldenKPM_M64x64_statedefaultJACKSON.dat (G4 curve, 10 meV)
IDENTICAL [T=1] graphene/sea_surf/KuboBastinSurf_NonEqOpVX-VXgraphene_goldenKPM_M64x64_statedefaultJACKSON.dat (G4 curve, 10 meV)
IDENTICAL [T=1] graphene/fft_grid/KuboGreenWood_NonEqOpVX-VXgraphene_goldenKPM_M64x64_statedefaultJACKSON_FFTgrid.dat (G4 curve, 10 meV)
IDENTICAL [T=1] graphene/fft_grid/KuboBastin_NonEqOpVX-VXgraphene_goldenKPM_M64x64_statedefaultJACKSON_FFTgrid.dat (G4 curve, 10 meV)
GREP      grep -rn "compat/" lsquant/core lsquant/io: empty
== 41 files compared; failures: 0
GATE PASS (stage 3)
gate exit: 0
```

### 4.2 Earlier runs

- Run 2 (`0774db7`): `figures/out/stage_3/check_goldens_stage3_run2_0774db7.log`,
  exit 0, GATE PASS, 41/41.
- Run 1 (`c01b026`): `figures/out/stage_3/check_goldens_stage3_run1_c01b026.log`,
  exit 0, GATE PASS, 41/41.

The verdict lines (IDENTICAL/DIFFERS/GREP) of the three runs are identical line by
line.

### 4.3 Extra check of the review-2 path that the gate does not exercise

Form 2 with `--formula` now reconstructs from the round-tripped moments. Running
`stage3_kubo --operators goldens/stage0/haldane/operators --label haldane --left VX
--right VY|VX --window -4.5 4.5 --trace exact --states 512 --moments 128 --formula
<f> --broadening-meV 10` at `cac5da6`, OMP_NUM_THREADS=16, and running `cmp` against
`goldens/stage0/haldane/` gives:

```
IDENTICAL KuboBastin_NonEqOpVX-VYhaldaneKPM_M128x128_stateexactJACKSON.dat
IDENTICAL KuboBastinSea_NonEqOpVX-VYhaldaneKPM_M128x128_stateexactJACKSON.dat
IDENTICAL KuboBastinSurf_NonEqOpVX-VYhaldaneKPM_M128x128_stateexactJACKSON.dat
IDENTICAL KuboGreenWood_NonEqOpVX-VXhaldaneKPM_M128x128_stateexactJACKSON.dat
```

## 5. Figure (`figures/stage_3.py`)

`MPLBACKEND=Agg python3 figures/stage_3.py --build build --threads 16` runs
`python3 -m lsquant.apps.tutorial03 --out <tmp>`. It compares the four `.dat` files
by name with `goldens/stage0/tutorials/03/`, which are byte-identical to the G5
`.dat` files in `goldens/stage0/haldane/` (md5 checked). It then draws the two panels.
Output: **`figures/out/stage_3/fig_stage_3.png`** and `fig_stage_3.caption.txt`; run
log in `figure_run.log`. The run-2 caption is kept as
`fig_stage_3_run2_0774db7.caption.txt`.

tutorial03 output and `cmp` list (final run, `cac5da6`):

```
tutorial03 printed (5 lines):
  <tmp>/KuboBastin_NonEqOpVX-VYhaldaneKPM_M128x128_stateexactJACKSON.dat
  <tmp>/KuboBastinSea_NonEqOpVX-VYhaldaneKPM_M128x128_stateexactJACKSON.dat
  <tmp>/KuboBastinSurf_NonEqOpVX-VYhaldaneKPM_M128x128_stateexactJACKSON.dat
  <tmp>/KuboGreenWood_NonEqOpVX-VXhaldaneKPM_M128x128_stateexactJACKSON.dat
  Chern from gap plateau: C = +1.017 (+1.017016743089)
cmp list:
  KuboBastin_NonEqOpVX-VYhaldaneKPM_M128x128_stateexactJACKSON.dat: IDENTICAL
  KuboBastinSea_NonEqOpVX-VYhaldaneKPM_M128x128_stateexactJACKSON.dat: IDENTICAL
  KuboBastinSurf_NonEqOpVX-VYhaldaneKPM_M128x128_stateexactJACKSON.dat: IDENTICAL
  KuboGreenWood_NonEqOpVX-VXhaldaneKPM_M128x128_stateexactJACKSON.dat: IDENTICAL
C new = +1.017016743089, C golden = +1.017016743089
wrote /home/ICN2/jgarcia/lsquant-unified/figures/out/stage_3/fig_stage_3.png
exit: 0
```

**C = +1.017016743089** from the new Bastin curve (tutorial03 prints `C = +1.017`),
and **C = +1.017016743089** from the golden curve.

Caption (verbatim, `fig_stage_3.caption.txt`):

> FIG. 3. Stage 3 (handoff 0003): tutorial 03, the Hall and longitudinal conductivity of the Haldane model (t1 = -1, t2 = 0.15, phi = pi/2), reproduced from Python by python3 -m lsquant.apps.tutorial03 (lsquant package over the pybind11 module lsquant._core, OMP_NUM_THREADS = 16): N = 512 orbitals, window [-4.5, 4.5] eV, exact trace over all N site states, M = 128 x 128 Chebyshev moments (VX-VY for sigma_xy, VX-VX for sigma_xx), Jackson kernel at eta = 10 meV; conductivities in units of e^2/h as lsqhall.py normalises them (divided by (N/2) A_cell e^2/h, A_cell = 3 sqrt(3)/2, e = hbar = 1). (a) sigma_xy: Kubo-Bastin total, Fermi-sea and Fermi-surface parts; new as lines, the stage-0 golden (goldens/stage0/tutorials/03) as open markers (48 evenly spaced points per curve), dotted line at 1; residual new - golden below. (b) sigma_xx (Kubo-Greenwood), same encoding. Chern number from the gap plateau, C = (mean over |E| < 0.5 eV of sigma_xy / A) 2 pi with A = (N/2) A_cell: C = +1.017016743089 from the new curve (tutorial03 prints +1.017016743089), C = +1.017016743089 from the golden curve. Max |new - golden| over the four curves: 0.0e+00 e^2/h, max energy-axis difference 0.0e+00 eV. Byte comparison with goldens/stage0/tutorials/03 by name (cmp criterion), 4 of 4 IDENTICAL: KuboBastin_NonEqOpVX-VYhaldaneKPM_M128x128_stateexactJACKSON.dat: IDENTICAL; KuboBastinSea_NonEqOpVX-VYhaldaneKPM_M128x128_stateexactJACKSON.dat: IDENTICAL; KuboBastinSurf_NonEqOpVX-VYhaldaneKPM_M128x128_stateexactJACKSON.dat: IDENTICAL; KuboGreenWood_NonEqOpVX-VXhaldaneKPM_M128x128_stateexactJACKSON.dat: IDENTICAL. All files IDENTICAL.

## 6. Runs 1 and 2: localisation of the tutorial-03 DIFFERS (history; fixed by review 2)

At `c01b026` and `0774db7`, tutorial 03 gave 1 of 4 IDENTICAL: Surf was IDENTICAL;
Bastin, Sea and Greenwood DIFFERED, with max |new - golden| 9.4e-13 e^2/h and an
identical energy axis. `golden_diff.py` at `0774db7`:

```
KuboBastin:    first-diff line 1 col 2, new=2.18661e-08 golden=2.1866e-08; max abs diff 1.000e-10 at line 186 col 2; max rel diff 9.711e-06
KuboBastinSea: first-diff line 2 col 2, new=4.60745e-08 golden=4.60743e-08; max abs diff 1.000e-10 at line 411 col 2; max rel diff 3.567e-01 at line 3839 (near-zero value)
KuboGreenWood: first-diff line 1 col 2, new=5.06063e-08 golden=5.06082e-08; max abs diff 1.900e-12 at line 1 col 2; max rel diff 3.754e-05
```

Diagnosis at that time:
- The staged operators were byte-identical to the goldens, and the `.desc` window was
  exactly [-4.5, 4.5].
- The gate showed the moments and the formulas byte-identical.
- The in-memory moments from `run_plan` differed from the golden file in the last bits
  (`VX-VY` max |dIm| 4.9e-15; `VX-VX` max |dRe| 5.3e-15; about 15,400 of 16,384
  elements). The golden file holds them at legacy text precision.
- The legacy `lsqhall.py` reconstructed from the written `.chebmom2D`, while
  `tutorial03.py` reconstructed from the unrounded in-memory moments.
- Surf was IDENTICAL because it uses only `real(mu)` of `VX-VY` (about 1e-31).
- Check: in-memory moments gave 1 of 4 IDENTICAL; the same moments passed through
  `WriteLegacyChebmom2D` + `ReadChebmom` gave 4 of 4 IDENTICAL.

Review 2 adopted that round trip, and run 3 confirms it.

## 7. Smells (final state)

- **Formula or grid in `core/`: none found.** `contraction.cpp` gains only the three
  rank-2 kernels (`m0` outer, `m1` inner, OpenMP over `i`, tables passed in). A grep of
  `lsquant/core` and `lsquant/io` for `pow|acos|0.95|M_PI|cos(|sqrt(` on code lines
  finds only the stage-1 random-phase generator (`state_generator.cpp`) and the
  Lanczos norms (`spectral_bounds.cpp`). The stage-2 code-line grep is empty (three
  comment lines excluded, listed in the gate output).
- **`compat` included from `core/`: no.** `grep -rn "compat/" lsquant/core lsquant/io`
  is empty. `compat/legacy_kernels.cpp` is compiled only into `_core`.
- **numpy where `math` was required: none found.** The prefactors, `ceil(pi/eta)`,
  the uniform and FFT-node grids, the 2D Jackson weights and the cumulative integrals
  are scalar `math` loops in `reconstruct.py` and `kernels.py`. numpy is used only for
  containers and for the elementwise table algebra of Sea (`np.conj(Gr)`, `Gr - Ga`,
  `DGr + DGa`). Those operations are exact IEEE operations, the same ones the legacy
  performs, and the G4/G5 Sea curves are byte-identical.
- **Table built with Python's complex power: no.** `greenR`/`DgreenR` come from
  `_core.compat.legacy_green_tables` (libstdc++ `std::pow(complex<double>, double)`).
  The sympy `GreenRetarded`/`GreenRetardedDerivative` are only compared (1.1e-14 abs;
  7.2e-15 rel).
- **Copy in `from_csr`:** only the accepted Eigen `Map` copy remains; strict dtypes and
  contiguity are enforced with a `TypeError` naming the argument (section 3).
- Self-check (a) still compares moments only after `_roundtrip_legacy_2d`. That is
  correct now that every reconstruction path round-trips, but it is the reason the
  run-1/2 tutorial-03 difference did not show in the self-check.
- Resolved by review 2: tutorial 03 now reconstructs from round-tripped moments, and
  `--operators` defaults to `goldens/stage0/haldane/operators`, so it no longer
  depends on the legacy `test/` tree.
- No golden, seed or tolerance was touched; `goldens/` is unchanged.
