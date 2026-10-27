# Handoff 0002 — tester report (stage 2 gate, tutorial 01 from Python)

Tester: Codex role (Opus subagent driven by Fable). Machine: TCN (`fatnode-tacn`),
checkout `~/lsquant-unified`, branch `handoff/0002-codex` cut from `cpp-unified`
(`3fecd12`) with `origin/handoff/0002-agi` (`656c6f2`) merged without conflict (merge
commit `fa17a2d`). Nothing under `lsquant/`, `python/lsquant/`, the legacy tree,
`goldens/`, a seed or a tolerance was edited.

**Result: GATE PASS (stage 2), 17 of 17 files IDENTICAL; tutorial 01 from Python 22 of
22 files IDENTICAL; code-line grep empty.**

## 1. Versions

```
pybind11: 3.1.0 /home/ICN2/jgarcia/.local/lib/python3.10/site-packages/pybind11/share/cmake/pybind11
sympy: 1.14.0
numpy: 2.2.6
scipy: 1.15.3
matplotlib: 3.10.9
Python 3.10.12
g++ (Ubuntu 11.4.0-1ubuntu1~22.04.3) 11.4.0
cmake version 3.22.1
```

pybind11 and sympy were already installed for the user (`pip install --user`), no
install was needed.

## 2. Build (fresh)

```
$ rm -rf build
$ cmake -S . -B build -DLSQUANT_BUILD_NEW_CORE=ON -DLSQUANT_BUILD_PYTHON=ON -Dpybind11_DIR="$(python3 -c 'import pybind11; print(pybind11.get_cmake_dir())')"
-- Found pybind11: /home/ICN2/jgarcia/.local/lib/python3.10/site-packages/pybind11/include (found version "3.1.0")
-- Build files have been written to: /home/ICN2/jgarcia/lsquant-unified/build
cmake exit: 0
$ cmake --build build -j 32
[ 57%] Linking CXX static library ../lib/liblsquant_core.a
[ 57%] Built target lsquant_core
[ 60%] Linking CXX executable stage1_moments
[ 98%] Linking CXX executable stage_1_check
[100%] Linking CXX shared module _core.cpython-310-x86_64-linux-gnu.so
[100%] Built target _core
build exit: 0
```

Full log: `figures/out/stage_2/build.log`. Flags: `lsquant_core` `-O3 -DNDEBUG -fPIC
-Wall -Wextra -fopenmp -std=c++17`; `_core` adds `-fvisibility=hidden -flto=auto`.
Warnings: **0 under `lsquant/`** (core, io, api, apps, tests/dev). The 101 warnings in
the log are all in the legacy tree (62 `include/chebyshev_moments.hpp` narrowing, 32
`include/units.hpp` inline variables, 7 `src/chebyshev_momentsLocal.cpp` unused
`fread` result), as in stage 1.

## 3. Developer self-check

```
$ PYTHONPATH=build/lsquant_core:python python3 tests/dev/stage_2_check.py
(a) dos_plan array_equal vs G1:  True
(a) kubo_plan array_equal vs G2: True
(b) DeltaFunction.table vs legacy_delta_table max_abs_diff = 2.220446e-16
(c) dos_curve vs G3 differing lines: 0 (total lines: 1590)
```

Exit 0. Identical to the developer's report. (Wall time about 10 min: TCN was at load
average 150 from other jobs and the script runs at the default 127 OpenMP threads.)

## 4. Gate (`scripts/check_goldens.sh 2 build`), verbatim

Harness extension (this handoff): stage 2 runs stage 1 in full (its comparisons and its
grep), then `python3 -m lsquant.apps.stage2_dos` in both forms at `OMP_NUM_THREADS=1`
and `16` with `PYTHONPATH=build/lsquant_core:python`, compares form (a) against G3 and
form (b) against G1 and G3, and runs `grep -rn "delta\|DOS\|Kubo\|basis\|std::cout\|getenv"
lsquant/core lsquant/io` with comment lines (first non-blank characters `//`, `*` or
`/*`) excluded from the verdict and listed below it. Log:
`figures/out/stage_2/check_goldens_stage2.log`.

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
== 17 files compared; failures: 0
GATE PASS (stage 2)
```

Exit code 0. No DIFFERS, so no localisation was needed. The only excluded lines are
the three stage-1 Doxygen comments with the notation `delta_{m0}` that review 1 rules
out as a finding.

## 5. Tutorial 01 from Python: `cmp` list

`figures/stage_2.py` runs `python3 -m lsquant.apps.tutorial01 --out <tmp>`
(`OMP_NUM_THREADS=16`) and matches the outputs to `goldens/stage0/tutorials/01/captured/`
as codex.md section 3 says: moment files by name after the 4-digit prefix (against every
captured copy of that name), the curves `eta20/N256, N1024, N4096` to `0013_, 0015_,
0017_` and `eta80/, eta24/, eta10/` `N256` to `0019_, 0020_, 0021_`. The app printed 12
paths (the 6 moment files at the top of `<tmp>`, the 6 curves under `<tmp>/eta<meV>/`),
exactly the list in the developer's report section 5; no prefixed files, no `captured/`
sequence.

```
SpectralOp1chain1d_N256KPM_M512_stateexact_chain1d_N256.chebmom1D   vs captured/0000_...: IDENTICAL
moments_chain1d_N256_M512.chebmom1D                                 vs captured/0001_...: IDENTICAL
SpectralOp1chain1d_N1024KPM_M512_stateexact_chain1d_N1024.chebmom1D vs captured/0002_...: IDENTICAL
moments_chain1d_N1024_M512.chebmom1D                                vs captured/0003_...: IDENTICAL
SpectralOp1chain1d_N4096KPM_M512_stateexact_chain1d_N4096.chebmom1D vs captured/0004_...: IDENTICAL
moments_chain1d_N4096_M512.chebmom1D                                vs captured/0005_...: IDENTICAL
SpectralOp1chain1d_N256KPM_M512_stateexact_chain1d_N256.chebmom1D   vs captured/0006_...: IDENTICAL
moments_chain1d_N256_M512.chebmom1D                                 vs captured/0007_...: IDENTICAL
SpectralOp1chain1d_N1024KPM_M512_stateexact_chain1d_N1024.chebmom1D vs captured/0008_...: IDENTICAL
moments_chain1d_N1024_M512.chebmom1D                                vs captured/0009_...: IDENTICAL
SpectralOp1chain1d_N4096KPM_M512_stateexact_chain1d_N4096.chebmom1D vs captured/0010_...: IDENTICAL
moments_chain1d_N4096_M512.chebmom1D                                vs captured/0011_...: IDENTICAL
SpectralOp1chain1d_N256KPM_M512_stateexact_chain1d_N256.chebmom1D   vs captured/0012_...: IDENTICAL
SpectralOp1chain1d_N1024KPM_M512_stateexact_chain1d_N1024.chebmom1D vs captured/0014_...: IDENTICAL
SpectralOp1chain1d_N4096KPM_M512_stateexact_chain1d_N4096.chebmom1D vs captured/0016_...: IDENTICAL
SpectralOp1chain1d_N256KPM_M512_stateexact_chain1d_N256.chebmom1D   vs captured/0018_...: IDENTICAL
eta20/meanmoments_chain1d_N256_M512JACKSON.dat   vs captured/0013_meanmoments_chain1d_N256_M512JACKSON.dat: IDENTICAL
eta20/meanmoments_chain1d_N1024_M512JACKSON.dat  vs captured/0015_meanmoments_chain1d_N1024_M512JACKSON.dat: IDENTICAL
eta20/meanmoments_chain1d_N4096_M512JACKSON.dat  vs captured/0017_meanmoments_chain1d_N4096_M512JACKSON.dat: IDENTICAL
eta80/meanmoments_chain1d_N256_M512JACKSON.dat   vs captured/0019_meanmoments_chain1d_N256_M512JACKSON.dat: IDENTICAL
eta24/meanmoments_chain1d_N256_M512JACKSON.dat   vs captured/0020_meanmoments_chain1d_N256_M512JACKSON.dat: IDENTICAL
eta10/meanmoments_chain1d_N256_M512JACKSON.dat   vs captured/0021_meanmoments_chain1d_N256_M512JACKSON.dat: IDENTICAL
```

22 of 22 IDENTICAL (full names in `figures/out/stage_2/stage_2_figure.log`). The
`N = 256` curve at 10 meV uses `M_eff = 512` (capped at `M`), so the whole moment set
and the longest grid (15360 points) are covered byte for byte.

## 6. Figure

`figures/out/stage_2/fig_stage_2.png`, caption `figures/out/stage_2/fig_stage_2.caption.txt`
(produced by `MPLBACKEND=Agg python3 figures/stage_2.py --build build --threads 16`,
exit 0; log `figures/out/stage_2/stage_2_figure.log`). Caption, with the per-file `cmp`
list that closes it shortened to the count (the list is section 5):

> FIG. 2. Stage 2 (handoff 0002): tutorial 01, the density of states of the periodic
> nearest-neighbour chain (t = -1 eV), reproduced from Python by python3 -m
> lsquant.apps.tutorial01 (lsquant package over the pybind11 module lsquant._core,
> OMP_NUM_THREADS = 16): N = 256, 1024, 4096 sites, window [-2, 2] eV, exact trace over
> all N site states, M = 512 Chebyshev moments, Jackson kernel, reconstruction on
> 30*M_eff points over |x| <= 0.95, per-site DOS rho(E) = value / N. (a) N = 256, 1024,
> 4096 at eta = 20 meV (M_eff = 315): new as lines, the stage-0 golden
> (goldens/stage0/tutorials/01/captured 0013, 0015, 0017) as open markers (48 evenly
> spaced points per curve), the analytic bulk DOS 1/(pi sqrt(4 - E^2)) as the thin grey
> line; residual new - golden below. (b) N = 256 at eta = 80, 24, 10 meV (M_eff = 79,
> 262, 512; captured 0019, 0020, 0021), same encoding. Max |new - golden| over all six
> curves: 0.0e+00 per site, max energy-axis difference 0.0e+00 eV. Byte comparison with
> the captured goldens (cmp criterion), 22 of 22 IDENTICAL: [...]. All files IDENTICAL.

The residual panels are identically zero. The 10 meV curve in (b) resolves the 256
discrete levels of the small chain (M_eff = 512 > N), as in the stage-0 tutorial figure.

## 7. Smells

1. **`SparseMatrix.from_csr` copies the three CSR buffers.** `lsquant/api/bindings.cpp`
   (the `from_csr` lambda) builds `std::vector<Index> rp(rowPointers.data(), ...)`,
   `cols(...)` and `std::vector<Scalar> vals(...)` and moves them into `FromCsr`. Brief
   section 2.5 says the bindings copy nothing and `FromCsr` uses the `int64` and
   `complex128` buffers as given; review 1 records `from_csr` as zero-copy, which the
   code does not do. In addition the arguments are `py::array_t<..., c_style>` with the
   default `forcecast`, so a scipy CSR with `int32` indices or a non-contiguous array is
   silently converted (a second copy) instead of rejected. Not exercised by the gate
   (`stage2_dos` and `tutorial01` read the CSR with `ReadCsrText`). `LocalStates` and
   `UserStates` also copy their inputs (small; the brief does not address them).
2. **numpy where `math` was required: none found.** `jackson_damping` and
   `legacy_delta_table` use `math.cos/sin/tan/acos/sqrt/pi` element by element in the
   legacy order; the grid and the energy axis are Python float expressions. The only
   numpy arithmetic on the byte path is the damping `raw[m] * g[m]` (a complex128 scalar
   times a float64 scalar, one IEEE multiply per component, the same as
   `std::complex<double> *= double`). `DeltaFunction.table` lambdifies to numpy by design
   (the general path; agrees to 2.2e-16).
3. **Physics name in C++: none.** Code-line grep empty (section 4); `contraction` and
   `curve_text` speak of moments, table, x and y.
4. **Table computed in C++: none.** `ContractRank1` receives the table from Python;
   `WriteCurveText` only formats.
5. **`print` in the package: only in the two apps.** `stage2_dos.py` prints one line per
   output file and the error on stderr (as brief section 6 asks); `tutorial01.py` prints
   its file list and the error on stderr. No `print` in the library modules.
6. Minor: `tutorial01` also writes `operators/`, `BOUNDS` and `exact_chain1d_N*` into
   `<out>` without listing them, and locates `make_chain.py` through a path relative to
   the installed package (`../../../examples/01_dos_1d_chain`), so it works only from a
   source checkout. `dos_curve` takes the band width from `moments.window`, not from a
   file header; equal for every golden here.

## 8. Files committed on `handoff/0002-codex`

`scripts/check_goldens.sh` (stage 2), `figures/stage_2.py`, `figures/out/stage_2/`
(`build.log`, `stage_2_check.log`, `check_goldens_stage2.log`, `stage_2_figure.log`,
`fig_stage_2.png`, `fig_stage_2.caption.txt`), this file.
