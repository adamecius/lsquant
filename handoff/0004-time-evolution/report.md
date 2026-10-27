# Handoff 0004 Report — Stage 4: Apply Mode, Time Evolution, MSD, and Projected Words (G6, G7, G8)

## 1. Files Added and Modified

- `lsquant/core/chebyshev/plan.hpp` (modified: added `Register = {Ket, Bra, Aux}`, `ApplyDerivedSeries`, `EvolveRegister`, `EvolveDerived`, `CopyRegister`, `AxisStep`, `TimeAxis`, and optional `scale` field on `ApplySeries` and `ApplyDerivedSeries` for the post-recursion band-center phase `scal`)
- `lsquant/core/chebyshev/chebyshev_moments.hpp` (modified: added `timeStep` and `numberOfTimeSteps` metadata fields to `ChebyshevMoments`)
- `lsquant/core/chebyshev/moment_engine.hpp`, `lsquant/core/chebyshev/moment_engine.cpp` (modified: implemented `ApplySeries`, `ApplyDerivedSeries`, generalised `MirrorToBra`, and `TimeAxis` with `{Ket, Bra, Aux}` registers in `MomentEngine::Run`)
- `lsquant/io/chebmom_file.hpp`, `lsquant/io/chebmom_file.cpp` (modified: added `WriteLegacyChebmomTD` and `.chebmomTD` auto-detection in `ReadChebmom`)
- `lsquant/compat/legacy_kernels.hpp`, `lsquant/compat/legacy_kernels.cpp` (modified: added `ComputeLegacyBesselCoefficients` using `gsl_sf_bessel_Jn`)
- `lsquant/CMakeLists.txt` (modified: linked `GSL::gsl` to `_core`; `lsquant_core` remains free of GSL)
- `lsquant/api/bindings.cpp` (modified: bound `Register`, `ApplySeries`, `ApplyDerivedSeries`, `EvolveRegister`, `EvolveDerived`, `CopyRegister`, `TimeAxis`, `ChebyshevMoments.time_step`, `ChebyshevMoments.number_of_time_steps`, `WriteLegacyChebmomTD`, and `_core.compat.bessel_series_coefficients`)
- `python/lsquant/_core.pyi` (modified: added type stubs for stage-4 bindings)
- `python/lsquant/spectral_functions.py` (modified: added `TimeEvolutionFunction`, `TimeEvolution`, `bessel_series_coefficients`, and `legacy_bessel_coefficients`)
- `python/lsquant/moments.py` (modified: added `HBAR`, `time_step_from_tmax`, `time_evolution_series`, `msd_plan`, `projected_evolution_plan`, `time_correlation_plan`, and updated `run_plan` for `TimeAxis` tensors)
- `python/lsquant/reconstruct.py` (modified: added `time_correlation_curve`)
- `python/lsquant/__init__.py` (modified: exported stage-4 symbols)
- `python/lsquant/apps/stage4_time.py` (added: stage-4 gate CLI)
- `python/lsquant/apps/tutorial02.py` (added: Tutorial 02 1D Anderson localisation application)
- `python/lsquant/apps/tutorial04.py` (added: Tutorial 04 Larmor spin precession application)
- `python/lsquant/apps/tutorial05.py` (added: Tutorial 05 2D square-lattice charge transport application)
- `tests/dev/stage_4_check.py` (added: stage-4 developer self-check)

---

## 2. CMake Configure and Build on TCN (`~/lsquant-agi-build`)

```
$ cd ~/lsquant-agi-build
$ cmake -S . -B build -DLSQUANT_BUILD_NEW_CORE=ON -DLSQUANT_BUILD_PYTHON=ON -DLSQUANT_BUILD_INPUT_TOOLS=OFF
$ cmake --build build --target lsquant_core stage1_moments stage_1_check _core -j 16
```

Verbatim configure and build output:

```
-- Found pybind11: /home/ICN2/jgarcia/.local/lib/python3.10/site-packages/pybind11/include (found version "3.1.0")
-- Configuring done
-- Generating done
-- Build files have been written to: /home/ICN2/jgarcia/lsquant-agi-build/build
Consolidate compiler generated dependencies of target lsquant_core
[  0%] Building CXX object lsquant_core/CMakeFiles/lsquant_core.dir/core/algebra/sparse_matrix.cpp.o
[  8%] Building CXX object lsquant_core/CMakeFiles/lsquant_core.dir/core/algebra/composite_operators.cpp.o
[ 16%] Building CXX object lsquant_core/CMakeFiles/lsquant_core.dir/core/chebyshev/rescaled_operator.cpp.o
[ 25%] Building CXX object lsquant_core/CMakeFiles/lsquant_core.dir/core/algebra/spectral_bounds.cpp.o
[ 33%] Building CXX object lsquant_core/CMakeFiles/lsquant_core.dir/core/chebyshev/chebyshev_moments.cpp.o
[ 41%] Building CXX object lsquant_core/CMakeFiles/lsquant_core.dir/core/chebyshev/moment_engine.cpp.o
[ 50%] Building CXX object lsquant_core/CMakeFiles/lsquant_core.dir/io/csr_text.cpp.o
[ 66%] Building CXX object lsquant_core/CMakeFiles/lsquant_core.dir/core/chebyshev/contraction.cpp.o
[ 66%] Building CXX object lsquant_core/CMakeFiles/lsquant_core.dir/io/chebmom_file.cpp.o
[ 66%] Linking CXX static library ../lib/liblsquant_core.a
[100%] Built target lsquant_core
Consolidate compiler generated dependencies of target lsquant_core
[ 92%] Built target lsquant_core
Consolidate compiler generated dependencies of target stage1_moments
[ 92%] Building CXX object lsquant_core/CMakeFiles/stage1_moments.dir/apps/stage1_moments.cpp.o
[100%] Linking CXX executable stage1_moments
[100%] Built target stage1_moments
[ 85%] Built target lsquant_core
Consolidate compiler generated dependencies of target stage_1_check
[ 92%] Building CXX object lsquant_core/CMakeFiles/stage_1_check.dir/__/tests/dev/stage_1_check.cpp.o
[100%] Linking CXX executable stage_1_check
[100%] Built target stage_1_check
[ 85%] Built target lsquant_core
Consolidate compiler generated dependencies of target _core
[ 85%] Building CXX object lsquant_core/CMakeFiles/_core.dir/api/bindings.cpp.o
[ 92%] Building CXX object lsquant_core/CMakeFiles/_core.dir/compat/legacy_kernels.cpp.o
[100%] Linking CXX shared module _core.cpython-310-x86_64-linux-gnu.so
[100%] Built target _core
```

Zero compiler warnings were emitted under `lsquant/` with `-Wall -Wextra`.

---

## 3. Bounds Source of G7 (`chain1d_mag`) and Legacy `RescaleCon2ChebyshevDomain` Defect Note

### Bounds Source of G7 (`chain1d_mag`)

In `src/inline_compute-kpm-TimeEvProjetedOp.cpp` (line 64), the legacy driver calls:

```cpp
spectral_bounds = chebyshev::utility::SpectralBounds(OP[0]);
```

without passing a `descriptor_path` argument (which defaults to `""`), and `test/spin_precession.sh` stages only `operators/chain1d_mag.{HAM,SX,SY,SZ,PX}.CSR` without a `BOUNDS` file. Consequently, `chebyshev::utility::SpectralBounds` takes branch 3 (the Gershgorin row-sum enclosure with `PAD = 0.10`): for `chain1d_mag.HAM.CSR` (`t = -1.0`, Zeeman half-splitting `0.1`), the maximum absolute row sum is `rho = 2 * 1.0 + 0.1 = 2.1`, so `[-rho * (1 + 0.10), +rho * (1 + 0.10)] = [-2.31, 2.31]` and `BandWidth = 4.62`, `BandCenter = 0`, matching the G7 golden header `128 4.62 0 96 1 ` (`window_from_gershgorin(ham)` reproduces this window directly).

### Legacy `RescaleCon2ChebyshevDomain` Defect (for `AGENTS.md` section 10)

In `include/chebyshev_moments.hpp` (lines 422-428), `MomentsTD::RescaleCon2ChebyshevDomain()` rescales the commutator / velocity matrix `Conmuter()` via:

```cpp
this->Conmuter().Rescale(this->ScaleFactor(), this->ShiftFactor());
```

which applies `V <- ScaleFactor * V + ShiftFactor * I` (`RescaledOperator(V, window)`), adding `ShiftFactor * I = -(BandCenter / HalfWidth) * I` to the diagonal of `[X, H~]`. Because `[X, a H + b I] = a [X, H]`, only the multiplicative `ScaleFactor` should rescale the commutator; adding `ShiftFactor * I` is only harmless when `BandCenter == 0` (`ShiftFactor == 0`, as in symmetric windows `[-W, W]`). For byte-for-byte compatibility with legacy MSD goldens, `stage4_time.py` and the tutorial apps wrap the velocity CSR matrix in `_core.RescaledOperator(vel_mat, window)`.

---

## 4. Developer Self-Check Output (`stage_1_check` and `tests/dev/stage_4_check.py`)

```
$ ./build/lsquant_core/stage_1_check
dos max_abs_error: 2.47025e-15
kubo max_abs_error: 1.13687e-13
OMP_NUM_THREADS=1 vs 8 bit-identical: YES

$ OMP_NUM_THREADS=1 PYTHONPATH=build/lsquant_core:python python3 -u tests/dev/stage_4_check.py
(a) G7 chain1d_mag Gershgorin window = [-2.31, 2.31], BandWidth = 4.62
(a) G7 TimeProjSX-PX array_equal vs golden: True
(a) G7 TimeProjSY-PX array_equal vs golden: True
(a) G7 TimeProjSZ-PX array_equal vs golden: True
(b) G6 CorrelationVXchain1d array_equal vs golden: True
(b) G6 EF=+0.0 curve differing lines: 0 (total lines: 64)
(b) G6 EF=+0.5 curve differing lines: 0 (total lines: 64)
(b) G6 EF=+1.0 curve differing lines: 0 (total lines: 64)
(b) G6 EF=+1.5 curve differing lines: 0 (total lines: 64)
(c) compat GSL Bessel vs SciPy jv max_abs_diff = 5.551115e-16
(d) time_correlation_plan chain1d VX-VX (M=64, steps=16, tmax=5, exact) array_equal vs inline_compute-kpm-CorrOp: True
```

### Exact Legacy Command Used in Self-Check (d)

In a temporary working directory staging `goldens/stage0/chain1d_ballistic/operators/chain1d.{HAM,VX}.CSR`, `BOUNDS` (`-2 2`), and an `exact` local-state file over all `256` sites (`local\n256\n0\n...\n255\n`):

```
KPM_SEED=12345 OMP_NUM_THREADS=1 /home/ICN2/jgarcia/lsquant-unified/build/inline_compute-kpm-CorrOp chain1d VX VX 64 16 5 exact
```

which writes `TimeCorrVX-VXchain1dKPM_M64_stateexact.chebmomTD`.

### Review 1 Self-Check Rerun (`tests/dev/stage_4_check.py`)

After updating `time_evolution_series` to default to `use_legacy_bessel=False` (with `stage4_time.py`, `tutorial02.py`, `tutorial04.py`, `tutorial05.py`, and the legacy-golden checks in `stage_4_check.py` passing `use_legacy_bessel=True` explicitly) and documenting the empty-`coeffs` zero-map semantics on `ApplySeries`, `ApplyDerivedSeries`, and `msd_plan`:

```
$ OMP_NUM_THREADS=1 PYTHONPATH=build/lsquant_core:python python3 -u tests/dev/stage_4_check.py
(a) G7 chain1d_mag Gershgorin window = [-2.31, 2.31], BandWidth = 4.62
(a) G7 TimeProjSX-PX array_equal vs golden: True
(a) G7 TimeProjSY-PX array_equal vs golden: True
(a) G7 TimeProjSZ-PX array_equal vs golden: True
(b) G6 CorrelationVXchain1d array_equal vs golden: True
(b) G6 EF=+0.0 curve differing lines: 0 (total lines: 64)
(b) G6 EF=+0.5 curve differing lines: 0 (total lines: 64)
(b) G6 EF=+1.0 curve differing lines: 0 (total lines: 64)
(b) G6 EF=+1.5 curve differing lines: 0 (total lines: 64)
(c) compat GSL Bessel vs SciPy jv max_abs_diff = 5.551115e-16
(d) time_correlation_plan chain1d VX-VX (M=64, steps=16, tmax=5, exact) array_equal vs inline_compute-kpm-CorrOp: True
```

### Review 2 Self-Check Rerun (`tests/dev/stage_4_check.py`)

After applying the three review-2 items (`build_chain` / `build_square` calls in `tutorial02.py` and `tutorial05.py`; renaming `TimeEvolutionStep` to `AxisStep` and removing physics terms from `plan.hpp`; defaulting `tutorial04.py --operators` to `goldens/stage0/chain1d_mag/operators` and documenting `--batch` in `stage4_time.py`):

```
$ OMP_NUM_THREADS=1 PYTHONPATH=build/lsquant_core:python python3 -u tests/dev/stage_4_check.py
(a) G7 chain1d_mag Gershgorin window = [-2.31, 2.31], BandWidth = 4.62
(a) G7 TimeProjSX-PX array_equal vs golden: True
(a) G7 TimeProjSY-PX array_equal vs golden: True
(a) G7 TimeProjSZ-PX array_equal vs golden: True
(b) G6 CorrelationVXchain1d array_equal vs golden: True
(b) G6 EF=+0.0 curve differing lines: 0 (total lines: 64)
(b) G6 EF=+0.5 curve differing lines: 0 (total lines: 64)
(b) G6 EF=+1.0 curve differing lines: 0 (total lines: 64)
(b) G6 EF=+1.5 curve differing lines: 0 (total lines: 64)
(c) compat GSL Bessel vs SciPy jv max_abs_diff = 5.551115e-16
(d) time_correlation_plan chain1d VX-VX (M=64, steps=16, tmax=5, exact) array_equal vs inline_compute-kpm-CorrOp: True
```

### Review 3 Self-Check (`tests/dev/stage_3_check.py`) and Off-Centre Window `stage3_kubo` Run

In `python/lsquant/reconstruct.py`, replaced all five occurrences of `energies[i] = x_grid[i] / scale + shift` (in `dos_curve`, `kubo_greenwood`, `kubo_bastin`, `kubo_bastin_sea`, and `kubo_bastin_surf`; `kubo_fft_grid` already mapped `energies[i] = x_grid[i] * hw + bc`) with `energies[i] = (x_grid[i] - shift) / scale` and updated the docstrings at lines 42 and 175.

```
$ OMP_NUM_THREADS=1 PYTHONPATH=build/lsquant_core:python python3 -u tests/dev/stage_3_check.py
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

$ OMP_NUM_THREADS=1 PYTHONPATH=build/lsquant_core:python python3 -u -m lsquant.apps.stage3_kubo \
    --operators goldens/stage0/graphene/operators --label graphene_golden \
    --left VX --right VX --window -3.5 2.5 --trace exact --moments 128 \
    --formula bastin --broadening-meV 10 --out /tmp/agi_r3_g4c
/tmp/agi_r3_g4c/NonEqOpVX-VXgraphene_goldenKPM_M128x128_stateexact.chebmom2D
/tmp/agi_r3_g4c/NonEqOpVX-VXgraphene_goldenKPM_M128x128_stateexact.chebmom
/tmp/agi_r3_g4c/KuboBastin_NonEqOpVX-VXgraphene_goldenKPM_M128x128_stateexactJACKSON.dat

$ head -n 3 /tmp/agi_r3_g4c/KuboBastin_NonEqOpVX-VXgraphene_goldenKPM_M128x128_stateexactJACKSON.dat
-3.35 1.00303e+150
-3.34852 1.23827e+150
-3.34703 6.85813e+149

$ tail -n 3 /tmp/agi_r3_g4c/KuboBastin_NonEqOpVX-VXgraphene_goldenKPM_M128x128_stateexactJACKSON.dat
2.34555 2.44004e+152
2.34703 2.00768e+152
2.34852 1.46351e+152
```

For window `[-3.5, 2.5]` (`HalfWidth = 3.0`, `Center = -0.5`, `Scale = 1/3`, `Shift = 1/6`), the first and last grid points map to `E_0 = (-0.95 - 1/6) / (1/3) = -3.35` and `E_{3838} = 2.34852`, both strictly inside `[-3.5, 2.5]` (whereas the old `x / scale + shift` produced `[-2.68333, 3.01518]`, exceeding `+2.5`). Note for G4c closure: in the legacy C++ tree, `src/spectralFunctionFromChebmom.cpp:54` writes `grid[i].first / mu.ScaleFactor() + mu.ShiftFactor()` and `src/reconstruction.cpp:148,200,224` writes `e / sf - shf` (both missing `/ sf` on `shf`, whereas `src/inline_kpm-DOS-standalone.cpp:144` and `src/reconstruction.cpp:151` use `e * hw + bc`), and `[-3.5, 2.5]` is narrower than the `graphene_golden` spectrum (Gershgorin radius ~7.5).

---

## 5. Gate CLI (`stage4_time.py`) and Tutorials 02, 04, and 05 Runs on TCN

```
$ OMP_NUM_THREADS=1 PYTHONPATH=build/lsquant_core:python python3 -m lsquant.apps.stage4_time \
    --operators goldens/stage0/chain1d_ballistic/operators --label chain1d --velocity VX \
    --window -2 2 --trace exact --states 256 --moments 256 --steps 64 --tmax 20 --word msd \
    --out /tmp/agi_s4_cli
/tmp/agi_s4_cli/CorrelationVXchain1dKPM_M256_state20exact.chebmomTD
/tmp/agi_s4_cli/CorrelationVXchain1dKPM_M256_state20exact.chebmom

$ OMP_NUM_THREADS=1 PYTHONPATH=build/lsquant_core:python python3 -m lsquant.apps.stage4_time \
    --moments /tmp/agi_s4_cli/CorrelationVXchain1dKPM_M256_state20exact.chebmomTD \
    --broadening-meV 10 --fermi-energy 0.0 --out /tmp/agi_s4_cli
/tmp/agi_s4_cli/meanCorrelationVXchain1dKPM_M256_state20exactEF0.000000JACKSON.dat

$ OMP_NUM_THREADS=1 PYTHONPATH=build/lsquant_core:python python3 -m lsquant.apps.tutorial04 \
    --out /tmp/agi_r2_t04
/tmp/agi_r2_t04/TimeProjSX-PXchain1d_magKPM_M32_stateexact.chebmomTD
/tmp/agi_r2_t04/TimeProjSY-PXchain1d_magKPM_M32_stateexact.chebmomTD
/tmp/agi_r2_t04/TimeProjSZ-PXchain1d_magKPM_M32_stateexact.chebmomTD

$ OMP_NUM_THREADS=1 PYTHONPATH=build/lsquant_core:python python3 -u -m lsquant.apps.tutorial02 \
    --out /tmp/agi_r2_t02
/tmp/agi_r2_t02/meanCorrelationVXchain1d_dis_N512_W1KPM_M256_state600defaultEF0.000000JACKSON.dat
/tmp/agi_r2_t02/captured/0000_meanCorrelationVXchain1d_dis_N256_W1KPM_M256_state600defaultEF0.000000JACKSON.dat
/tmp/agi_r2_t02/captured/0001_meanCorrelationVXchain1d_dis_N256_W1KPM_M256_state600defaultEF1.000000JACKSON.dat
/tmp/agi_r2_t02/captured/0002_meanCorrelationVXchain1d_dis_N256_W1KPM_M256_state600defaultEF-1.000000JACKSON.dat
/tmp/agi_r2_t02/captured/0003_meanCorrelationVXchain1d_dis_N256_W1KPM_M256_state600defaultEF1.900000JACKSON.dat
/tmp/agi_r2_t02/captured/0004_meanCorrelationVXchain1d_dis_N256_W1KPM_M256_state600defaultEF0.000000JACKSON.dat
/tmp/agi_r2_t02/captured/0005_meanCorrelationVXchain1d_dis_N256_W1KPM_M256_state600defaultEF1.000000JACKSON.dat
/tmp/agi_r2_t02/captured/0006_meanCorrelationVXchain1d_dis_N256_W1KPM_M256_state600defaultEF-1.000000JACKSON.dat
/tmp/agi_r2_t02/captured/0007_meanCorrelationVXchain1d_dis_N256_W1KPM_M256_state600defaultEF1.900000JACKSON.dat
/tmp/agi_r2_t02/captured/0008_meanCorrelationVXchain1d_dis_N256_W1KPM_M256_state600defaultEF0.000000JACKSON.dat
/tmp/agi_r2_t02/captured/0009_meanCorrelationVXchain1d_dis_N256_W1KPM_M256_state600defaultEF1.000000JACKSON.dat
/tmp/agi_r2_t02/captured/0010_meanCorrelationVXchain1d_dis_N256_W1KPM_M256_state600defaultEF-1.000000JACKSON.dat
/tmp/agi_r2_t02/captured/0011_meanCorrelationVXchain1d_dis_N256_W1KPM_M256_state600defaultEF1.900000JACKSON.dat
/tmp/agi_r2_t02/captured/0012_meanCorrelationVXchain1d_dis_N256_W1KPM_M256_state600defaultEF0.000000JACKSON.dat
/tmp/agi_r2_t02/captured/0013_meanCorrelationVXchain1d_dis_N256_W1KPM_M256_state600defaultEF1.000000JACKSON.dat
/tmp/agi_r2_t02/captured/0014_meanCorrelationVXchain1d_dis_N256_W1KPM_M256_state600defaultEF-1.000000JACKSON.dat
/tmp/agi_r2_t02/captured/0015_meanCorrelationVXchain1d_dis_N256_W1KPM_M256_state600defaultEF1.900000JACKSON.dat
/tmp/agi_r2_t02/captured/0016_meanCorrelationVXchain1d_dis_N256_W1KPM_M256_state600defaultEF0.000000JACKSON.dat
/tmp/agi_r2_t02/captured/0017_meanCorrelationVXchain1d_dis_N256_W1KPM_M256_state600defaultEF1.000000JACKSON.dat
/tmp/agi_r2_t02/captured/0018_meanCorrelationVXchain1d_dis_N256_W1KPM_M256_state600defaultEF-1.000000JACKSON.dat
/tmp/agi_r2_t02/captured/0019_meanCorrelationVXchain1d_dis_N256_W1KPM_M256_state600defaultEF1.900000JACKSON.dat
/tmp/agi_r2_t02/captured/0020_meanCorrelationVXchain1d_dis_N256_W1KPM_M256_state600defaultEF0.000000JACKSON.dat
/tmp/agi_r2_t02/captured/0021_meanCorrelationVXchain1d_dis_N256_W1KPM_M256_state600defaultEF1.000000JACKSON.dat
/tmp/agi_r2_t02/captured/0022_meanCorrelationVXchain1d_dis_N256_W1KPM_M256_state600defaultEF-1.000000JACKSON.dat
/tmp/agi_r2_t02/captured/0023_meanCorrelationVXchain1d_dis_N256_W1KPM_M256_state600defaultEF1.900000JACKSON.dat
/tmp/agi_r2_t02/captured/0024_meanCorrelationVXchain1d_dis_N256_W1KPM_M256_state600defaultEF0.000000JACKSON.dat
/tmp/agi_r2_t02/captured/0025_meanCorrelationVXchain1d_dis_N256_W1KPM_M256_state600defaultEF1.000000JACKSON.dat
/tmp/agi_r2_t02/captured/0026_meanCorrelationVXchain1d_dis_N256_W1KPM_M256_state600defaultEF-1.000000JACKSON.dat
/tmp/agi_r2_t02/captured/0027_meanCorrelationVXchain1d_dis_N256_W1KPM_M256_state600defaultEF1.900000JACKSON.dat
/tmp/agi_r2_t02/captured/0028_meanCorrelationVXchain1d_dis_N256_W1KPM_M256_state600defaultEF0.000000JACKSON.dat
/tmp/agi_r2_t02/captured/0029_meanCorrelationVXchain1d_dis_N256_W1KPM_M256_state600defaultEF1.000000JACKSON.dat
/tmp/agi_r2_t02/captured/0030_meanCorrelationVXchain1d_dis_N256_W1KPM_M256_state600defaultEF-1.000000JACKSON.dat
/tmp/agi_r2_t02/captured/0031_meanCorrelationVXchain1d_dis_N256_W1KPM_M256_state600defaultEF1.900000JACKSON.dat
/tmp/agi_r2_t02/captured/0032_meanCorrelationVXchain1d_dis_N256_W1KPM_M256_state600defaultEF0.000000JACKSON.dat
/tmp/agi_r2_t02/captured/0033_meanCorrelationVXchain1d_dis_N256_W1KPM_M256_state600defaultEF1.000000JACKSON.dat
/tmp/agi_r2_t02/captured/0034_meanCorrelationVXchain1d_dis_N256_W1KPM_M256_state600defaultEF-1.000000JACKSON.dat
/tmp/agi_r2_t02/captured/0035_meanCorrelationVXchain1d_dis_N256_W1KPM_M256_state600defaultEF1.900000JACKSON.dat
/tmp/agi_r2_t02/captured/0036_meanCorrelationVXchain1d_dis_N256_W1KPM_M256_state600defaultEF0.000000JACKSON.dat
/tmp/agi_r2_t02/captured/0037_meanCorrelationVXchain1d_dis_N256_W1KPM_M256_state600defaultEF1.000000JACKSON.dat
/tmp/agi_r2_t02/captured/0038_meanCorrelationVXchain1d_dis_N256_W1KPM_M256_state600defaultEF-1.000000JACKSON.dat
/tmp/agi_r2_t02/captured/0039_meanCorrelationVXchain1d_dis_N256_W1KPM_M256_state600defaultEF1.900000JACKSON.dat
/tmp/agi_r2_t02/captured/0040_meanCorrelationVXchain1d_dis_N256_W1KPM_M256_state600defaultEF0.000000JACKSON.dat
/tmp/agi_r2_t02/captured/0041_meanCorrelationVXchain1d_dis_N256_W1KPM_M256_state600defaultEF1.000000JACKSON.dat
/tmp/agi_r2_t02/captured/0042_meanCorrelationVXchain1d_dis_N256_W1KPM_M256_state600defaultEF-1.000000JACKSON.dat
/tmp/agi_r2_t02/captured/0043_meanCorrelationVXchain1d_dis_N256_W1KPM_M256_state600defaultEF1.900000JACKSON.dat
/tmp/agi_r2_t02/captured/0044_meanCorrelationVXchain1d_dis_N256_W1KPM_M256_state600defaultEF0.000000JACKSON.dat
/tmp/agi_r2_t02/captured/0045_meanCorrelationVXchain1d_dis_N256_W1KPM_M256_state600defaultEF1.000000JACKSON.dat
/tmp/agi_r2_t02/captured/0046_meanCorrelationVXchain1d_dis_N256_W1KPM_M256_state600defaultEF-1.000000JACKSON.dat
/tmp/agi_r2_t02/captured/0047_meanCorrelationVXchain1d_dis_N256_W1KPM_M256_state600defaultEF1.900000JACKSON.dat
/tmp/agi_r2_t02/captured/0048_meanCorrelationVXchain1d_dis_N256_W1KPM_M256_state600defaultEF0.000000JACKSON.dat
/tmp/agi_r2_t02/captured/0049_meanCorrelationVXchain1d_dis_N256_W1KPM_M256_state600defaultEF1.000000JACKSON.dat
/tmp/agi_r2_t02/captured/0050_meanCorrelationVXchain1d_dis_N256_W1KPM_M256_state600defaultEF-1.000000JACKSON.dat
/tmp/agi_r2_t02/captured/0051_meanCorrelationVXchain1d_dis_N256_W1KPM_M256_state600defaultEF1.900000JACKSON.dat
/tmp/agi_r2_t02/captured/0052_meanCorrelationVXchain1d_dis_N256_W1KPM_M256_state600defaultEF0.000000JACKSON.dat
/tmp/agi_r2_t02/captured/0053_meanCorrelationVXchain1d_dis_N256_W1KPM_M256_state600defaultEF1.000000JACKSON.dat
/tmp/agi_r2_t02/captured/0054_meanCorrelationVXchain1d_dis_N256_W1KPM_M256_state600defaultEF-1.000000JACKSON.dat
/tmp/agi_r2_t02/captured/0055_meanCorrelationVXchain1d_dis_N256_W1KPM_M256_state600defaultEF1.900000JACKSON.dat
/tmp/agi_r2_t02/captured/0056_meanCorrelationVXchain1d_dis_N256_W1KPM_M256_state600defaultEF0.000000JACKSON.dat
/tmp/agi_r2_t02/captured/0057_meanCorrelationVXchain1d_dis_N256_W1KPM_M256_state600defaultEF1.000000JACKSON.dat
/tmp/agi_r2_t02/captured/0058_meanCorrelationVXchain1d_dis_N256_W1KPM_M256_state600defaultEF-1.000000JACKSON.dat
/tmp/agi_r2_t02/captured/0059_meanCorrelationVXchain1d_dis_N256_W1KPM_M256_state600defaultEF1.900000JACKSON.dat
/tmp/agi_r2_t02/captured/0060_meanCorrelationVXchain1d_dis_N256_W1KPM_M256_state600defaultEF0.000000JACKSON.dat
/tmp/agi_r2_t02/captured/0061_meanCorrelationVXchain1d_dis_N256_W1KPM_M256_state600defaultEF1.000000JACKSON.dat
/tmp/agi_r2_t02/captured/0062_meanCorrelationVXchain1d_dis_N256_W1KPM_M256_state600defaultEF-1.000000JACKSON.dat
/tmp/agi_r2_t02/captured/0063_meanCorrelationVXchain1d_dis_N256_W1KPM_M256_state600defaultEF1.900000JACKSON.dat
/tmp/agi_r2_t02/captured/0064_meanCorrelationVXchain1d_dis_N1024_W1KPM_M256_state600defaultEF0.000000JACKSON.dat
/tmp/agi_r2_t02/captured/0065_meanCorrelationVXchain1d_dis_N1024_W1KPM_M256_state600defaultEF1.000000JACKSON.dat
/tmp/agi_r2_t02/captured/0066_meanCorrelationVXchain1d_dis_N1024_W1KPM_M256_state600defaultEF-1.000000JACKSON.dat
/tmp/agi_r2_t02/captured/0067_meanCorrelationVXchain1d_dis_N1024_W1KPM_M256_state600defaultEF1.900000JACKSON.dat
/tmp/agi_r2_t02/captured/0068_meanCorrelationVXchain1d_dis_N1024_W1KPM_M256_state600defaultEF0.000000JACKSON.dat
/tmp/agi_r2_t02/captured/0069_meanCorrelationVXchain1d_dis_N1024_W1KPM_M256_state600defaultEF1.000000JACKSON.dat
/tmp/agi_r2_t02/captured/0070_meanCorrelationVXchain1d_dis_N1024_W1KPM_M256_state600defaultEF-1.000000JACKSON.dat
/tmp/agi_r2_t02/captured/0071_meanCorrelationVXchain1d_dis_N1024_W1KPM_M256_state600defaultEF1.900000JACKSON.dat
/tmp/agi_r2_t02/captured/0072_meanCorrelationVXchain1d_dis_N1024_W1KPM_M256_state600defaultEF0.000000JACKSON.dat
/tmp/agi_r2_t02/captured/0073_meanCorrelationVXchain1d_dis_N1024_W1KPM_M256_state600defaultEF1.000000JACKSON.dat
/tmp/agi_r2_t02/captured/0074_meanCorrelationVXchain1d_dis_N1024_W1KPM_M256_state600defaultEF-1.000000JACKSON.dat
/tmp/agi_r2_t02/captured/0075_meanCorrelationVXchain1d_dis_N1024_W1KPM_M256_state600defaultEF1.900000JACKSON.dat
/tmp/agi_r2_t02/captured/0076_meanCorrelationVXchain1d_dis_N1024_W1KPM_M256_state600defaultEF0.000000JACKSON.dat
/tmp/agi_r2_t02/captured/0077_meanCorrelationVXchain1d_dis_N1024_W1KPM_M256_state600defaultEF1.000000JACKSON.dat
/tmp/agi_r2_t02/captured/0078_meanCorrelationVXchain1d_dis_N1024_W1KPM_M256_state600defaultEF-1.000000JACKSON.dat
/tmp/agi_r2_t02/captured/0079_meanCorrelationVXchain1d_dis_N1024_W1KPM_M256_state600defaultEF1.900000JACKSON.dat
/tmp/agi_r2_t02/captured/0080_meanCorrelationVXchain1d_dis_N1024_W1KPM_M256_state600defaultEF0.000000JACKSON.dat
/tmp/agi_r2_t02/captured/0081_meanCorrelationVXchain1d_dis_N1024_W1KPM_M256_state600defaultEF1.000000JACKSON.dat
/tmp/agi_r2_t02/captured/0082_meanCorrelationVXchain1d_dis_N1024_W1KPM_M256_state600defaultEF-1.000000JACKSON.dat
/tmp/agi_r2_t02/captured/0083_meanCorrelationVXchain1d_dis_N1024_W1KPM_M256_state600defaultEF1.900000JACKSON.dat
/tmp/agi_r2_t02/captured/0084_meanCorrelationVXchain1d_dis_N1024_W1KPM_M256_state600defaultEF0.000000JACKSON.dat
/tmp/agi_r2_t02/captured/0085_meanCorrelationVXchain1d_dis_N1024_W1KPM_M256_state600defaultEF1.000000JACKSON.dat
/tmp/agi_r2_t02/captured/0086_meanCorrelationVXchain1d_dis_N1024_W1KPM_M256_state600defaultEF-1.000000JACKSON.dat
/tmp/agi_r2_t02/captured/0087_meanCorrelationVXchain1d_dis_N1024_W1KPM_M256_state600defaultEF1.900000JACKSON.dat
/tmp/agi_r2_t02/captured/0088_meanCorrelationVXchain1d_dis_N1024_W1KPM_M256_state600defaultEF0.000000JACKSON.dat
/tmp/agi_r2_t02/captured/0089_meanCorrelationVXchain1d_dis_N1024_W1KPM_M256_state600defaultEF1.000000JACKSON.dat
/tmp/agi_r2_t02/captured/0090_meanCorrelationVXchain1d_dis_N1024_W1KPM_M256_state600defaultEF-1.000000JACKSON.dat
/tmp/agi_r2_t02/captured/0091_meanCorrelationVXchain1d_dis_N1024_W1KPM_M256_state600defaultEF1.900000JACKSON.dat
/tmp/agi_r2_t02/captured/0092_meanCorrelationVXchain1d_dis_N1024_W1KPM_M256_state600defaultEF0.000000JACKSON.dat
/tmp/agi_r2_t02/captured/0093_meanCorrelationVXchain1d_dis_N1024_W1KPM_M256_state600defaultEF1.000000JACKSON.dat
/tmp/agi_r2_t02/captured/0094_meanCorrelationVXchain1d_dis_N1024_W1KPM_M256_state600defaultEF-1.000000JACKSON.dat
/tmp/agi_r2_t02/captured/0095_meanCorrelationVXchain1d_dis_N1024_W1KPM_M256_state600defaultEF1.900000JACKSON.dat
/tmp/agi_r2_t02/captured/0096_meanCorrelationVXchain1d_dis_N1024_W1KPM_M256_state600defaultEF0.000000JACKSON.dat
/tmp/agi_r2_t02/captured/0097_meanCorrelationVXchain1d_dis_N1024_W1KPM_M256_state600defaultEF1.000000JACKSON.dat
/tmp/agi_r2_t02/captured/0098_meanCorrelationVXchain1d_dis_N1024_W1KPM_M256_state600defaultEF-1.000000JACKSON.dat
/tmp/agi_r2_t02/captured/0099_meanCorrelationVXchain1d_dis_N1024_W1KPM_M256_state600defaultEF1.900000JACKSON.dat
/tmp/agi_r2_t02/captured/0100_meanCorrelationVXchain1d_dis_N1024_W1KPM_M256_state600defaultEF0.000000JACKSON.dat
/tmp/agi_r2_t02/captured/0101_meanCorrelationVXchain1d_dis_N1024_W1KPM_M256_state600defaultEF1.000000JACKSON.dat
/tmp/agi_r2_t02/captured/0102_meanCorrelationVXchain1d_dis_N1024_W1KPM_M256_state600defaultEF-1.000000JACKSON.dat
/tmp/agi_r2_t02/captured/0103_meanCorrelationVXchain1d_dis_N1024_W1KPM_M256_state600defaultEF1.900000JACKSON.dat
/tmp/agi_r2_t02/captured/0104_meanCorrelationVXchain1d_dis_N1024_W1KPM_M256_state600defaultEF0.000000JACKSON.dat
/tmp/agi_r2_t02/captured/0105_meanCorrelationVXchain1d_dis_N1024_W1KPM_M256_state600defaultEF1.000000JACKSON.dat
/tmp/agi_r2_t02/captured/0106_meanCorrelationVXchain1d_dis_N1024_W1KPM_M256_state600defaultEF-1.000000JACKSON.dat
/tmp/agi_r2_t02/captured/0107_meanCorrelationVXchain1d_dis_N1024_W1KPM_M256_state600defaultEF1.900000JACKSON.dat
/tmp/agi_r2_t02/captured/0108_meanCorrelationVXchain1d_dis_N1024_W1KPM_M256_state600defaultEF0.000000JACKSON.dat
/tmp/agi_r2_t02/captured/0109_meanCorrelationVXchain1d_dis_N1024_W1KPM_M256_state600defaultEF1.000000JACKSON.dat
/tmp/agi_r2_t02/captured/0110_meanCorrelationVXchain1d_dis_N1024_W1KPM_M256_state600defaultEF-1.000000JACKSON.dat
/tmp/agi_r2_t02/captured/0111_meanCorrelationVXchain1d_dis_N1024_W1KPM_M256_state600defaultEF1.900000JACKSON.dat
/tmp/agi_r2_t02/captured/0112_meanCorrelationVXchain1d_dis_N1024_W1KPM_M256_state600defaultEF0.000000JACKSON.dat
/tmp/agi_r2_t02/captured/0113_meanCorrelationVXchain1d_dis_N1024_W1KPM_M256_state600defaultEF1.000000JACKSON.dat
/tmp/agi_r2_t02/captured/0114_meanCorrelationVXchain1d_dis_N1024_W1KPM_M256_state600defaultEF-1.000000JACKSON.dat
/tmp/agi_r2_t02/captured/0115_meanCorrelationVXchain1d_dis_N1024_W1KPM_M256_state600defaultEF1.900000JACKSON.dat
/tmp/agi_r2_t02/captured/0116_meanCorrelationVXchain1d_dis_N1024_W1KPM_M256_state600defaultEF0.000000JACKSON.dat
/tmp/agi_r2_t02/captured/0117_meanCorrelationVXchain1d_dis_N1024_W1KPM_M256_state600defaultEF1.000000JACKSON.dat
/tmp/agi_r2_t02/captured/0118_meanCorrelationVXchain1d_dis_N1024_W1KPM_M256_state600defaultEF-1.000000JACKSON.dat
/tmp/agi_r2_t02/captured/0119_meanCorrelationVXchain1d_dis_N1024_W1KPM_M256_state600defaultEF1.900000JACKSON.dat
/tmp/agi_r2_t02/captured/0120_meanCorrelationVXchain1d_dis_N1024_W1KPM_M256_state600defaultEF0.000000JACKSON.dat
/tmp/agi_r2_t02/captured/0121_meanCorrelationVXchain1d_dis_N1024_W1KPM_M256_state600defaultEF1.000000JACKSON.dat
/tmp/agi_r2_t02/captured/0122_meanCorrelationVXchain1d_dis_N1024_W1KPM_M256_state600defaultEF-1.000000JACKSON.dat
/tmp/agi_r2_t02/captured/0123_meanCorrelationVXchain1d_dis_N1024_W1KPM_M256_state600defaultEF1.900000JACKSON.dat
/tmp/agi_r2_t02/captured/0124_meanCorrelationVXchain1d_dis_N1024_W1KPM_M256_state600defaultEF0.000000JACKSON.dat
/tmp/agi_r2_t02/captured/0125_meanCorrelationVXchain1d_dis_N1024_W1KPM_M256_state600defaultEF1.000000JACKSON.dat
/tmp/agi_r2_t02/captured/0126_meanCorrelationVXchain1d_dis_N1024_W1KPM_M256_state600defaultEF-1.000000JACKSON.dat
/tmp/agi_r2_t02/captured/0127_meanCorrelationVXchain1d_dis_N1024_W1KPM_M256_state600defaultEF1.900000JACKSON.dat
/tmp/agi_r2_t02/captured/0128_meanCorrelationVXchain1d_dis_N4096_W1KPM_M256_state600defaultEF0.000000JACKSON.dat
/tmp/agi_r2_t02/captured/0129_meanCorrelationVXchain1d_dis_N4096_W1KPM_M256_state600defaultEF1.000000JACKSON.dat
/tmp/agi_r2_t02/captured/0130_meanCorrelationVXchain1d_dis_N4096_W1KPM_M256_state600defaultEF-1.000000JACKSON.dat
/tmp/agi_r2_t02/captured/0131_meanCorrelationVXchain1d_dis_N4096_W1KPM_M256_state600defaultEF1.900000JACKSON.dat
/tmp/agi_r2_t02/captured/0132_meanCorrelationVXchain1d_dis_N4096_W1KPM_M256_state600defaultEF0.000000JACKSON.dat
/tmp/agi_r2_t02/captured/0133_meanCorrelationVXchain1d_dis_N4096_W1KPM_M256_state600defaultEF1.000000JACKSON.dat
/tmp/agi_r2_t02/captured/0134_meanCorrelationVXchain1d_dis_N4096_W1KPM_M256_state600defaultEF-1.000000JACKSON.dat
/tmp/agi_r2_t02/captured/0135_meanCorrelationVXchain1d_dis_N4096_W1KPM_M256_state600defaultEF1.900000JACKSON.dat
/tmp/agi_r2_t02/captured/0136_meanCorrelationVXchain1d_dis_N4096_W1KPM_M256_state600defaultEF0.000000JACKSON.dat
/tmp/agi_r2_t02/captured/0137_meanCorrelationVXchain1d_dis_N4096_W1KPM_M256_state600defaultEF1.000000JACKSON.dat
/tmp/agi_r2_t02/captured/0138_meanCorrelationVXchain1d_dis_N4096_W1KPM_M256_state600defaultEF-1.000000JACKSON.dat
/tmp/agi_r2_t02/captured/0139_meanCorrelationVXchain1d_dis_N4096_W1KPM_M256_state600defaultEF1.900000JACKSON.dat
/tmp/agi_r2_t02/captured/0140_meanCorrelationVXchain1d_dis_N4096_W1KPM_M256_state600defaultEF0.000000JACKSON.dat
/tmp/agi_r2_t02/captured/0141_meanCorrelationVXchain1d_dis_N4096_W1KPM_M256_state600defaultEF1.000000JACKSON.dat
/tmp/agi_r2_t02/captured/0142_meanCorrelationVXchain1d_dis_N4096_W1KPM_M256_state600defaultEF-1.000000JACKSON.dat
/tmp/agi_r2_t02/captured/0143_meanCorrelationVXchain1d_dis_N4096_W1KPM_M256_state600defaultEF1.900000JACKSON.dat
/tmp/agi_r2_t02/captured/0144_meanCorrelationVXchain1d_dis_N4096_W1KPM_M256_state600defaultEF0.000000JACKSON.dat
/tmp/agi_r2_t02/captured/0145_meanCorrelationVXchain1d_dis_N4096_W1KPM_M256_state600defaultEF1.000000JACKSON.dat
/tmp/agi_r2_t02/captured/0146_meanCorrelationVXchain1d_dis_N4096_W1KPM_M256_state600defaultEF-1.000000JACKSON.dat
/tmp/agi_r2_t02/captured/0147_meanCorrelationVXchain1d_dis_N4096_W1KPM_M256_state600defaultEF1.900000JACKSON.dat
/tmp/agi_r2_t02/captured/0148_meanCorrelationVXchain1d_dis_N4096_W1KPM_M256_state600defaultEF0.000000JACKSON.dat
/tmp/agi_r2_t02/captured/0149_meanCorrelationVXchain1d_dis_N4096_W1KPM_M256_state600defaultEF1.000000JACKSON.dat
/tmp/agi_r2_t02/captured/0150_meanCorrelationVXchain1d_dis_N4096_W1KPM_M256_state600defaultEF-1.000000JACKSON.dat
/tmp/agi_r2_t02/captured/0151_meanCorrelationVXchain1d_dis_N4096_W1KPM_M256_state600defaultEF1.900000JACKSON.dat
/tmp/agi_r2_t02/captured/0152_meanCorrelationVXchain1d_dis_N4096_W1KPM_M256_state600defaultEF0.000000JACKSON.dat
/tmp/agi_r2_t02/captured/0153_meanCorrelationVXchain1d_dis_N4096_W1KPM_M256_state600defaultEF1.000000JACKSON.dat
/tmp/agi_r2_t02/captured/0154_meanCorrelationVXchain1d_dis_N4096_W1KPM_M256_state600defaultEF-1.000000JACKSON.dat
/tmp/agi_r2_t02/captured/0155_meanCorrelationVXchain1d_dis_N4096_W1KPM_M256_state600defaultEF1.900000JACKSON.dat
/tmp/agi_r2_t02/captured/0156_meanCorrelationVXchain1d_dis_N4096_W1KPM_M256_state600defaultEF0.000000JACKSON.dat
/tmp/agi_r2_t02/captured/0157_meanCorrelationVXchain1d_dis_N4096_W1KPM_M256_state600defaultEF1.000000JACKSON.dat
/tmp/agi_r2_t02/captured/0158_meanCorrelationVXchain1d_dis_N4096_W1KPM_M256_state600defaultEF-1.000000JACKSON.dat
/tmp/agi_r2_t02/captured/0159_meanCorrelationVXchain1d_dis_N4096_W1KPM_M256_state600defaultEF1.900000JACKSON.dat
/tmp/agi_r2_t02/captured/0160_meanCorrelationVXchain1d_dis_N4096_W1KPM_M256_state600defaultEF0.000000JACKSON.dat
/tmp/agi_r2_t02/captured/0161_meanCorrelationVXchain1d_dis_N4096_W1KPM_M256_state600defaultEF1.000000JACKSON.dat
/tmp/agi_r2_t02/captured/0162_meanCorrelationVXchain1d_dis_N4096_W1KPM_M256_state600defaultEF-1.000000JACKSON.dat
/tmp/agi_r2_t02/captured/0163_meanCorrelationVXchain1d_dis_N4096_W1KPM_M256_state600defaultEF1.900000JACKSON.dat
/tmp/agi_r2_t02/captured/0164_meanCorrelationVXchain1d_dis_N4096_W1KPM_M256_state600defaultEF0.000000JACKSON.dat
/tmp/agi_r2_t02/captured/0165_meanCorrelationVXchain1d_dis_N4096_W1KPM_M256_state600defaultEF1.000000JACKSON.dat
/tmp/agi_r2_t02/captured/0166_meanCorrelationVXchain1d_dis_N4096_W1KPM_M256_state600defaultEF-1.000000JACKSON.dat
/tmp/agi_r2_t02/captured/0167_meanCorrelationVXchain1d_dis_N4096_W1KPM_M256_state600defaultEF1.900000JACKSON.dat
/tmp/agi_r2_t02/captured/0168_meanCorrelationVXchain1d_dis_N4096_W1KPM_M256_state600defaultEF0.000000JACKSON.dat
/tmp/agi_r2_t02/captured/0169_meanCorrelationVXchain1d_dis_N4096_W1KPM_M256_state600defaultEF1.000000JACKSON.dat
/tmp/agi_r2_t02/captured/0170_meanCorrelationVXchain1d_dis_N4096_W1KPM_M256_state600defaultEF-1.000000JACKSON.dat
/tmp/agi_r2_t02/captured/0171_meanCorrelationVXchain1d_dis_N4096_W1KPM_M256_state600defaultEF1.900000JACKSON.dat
/tmp/agi_r2_t02/captured/0172_meanCorrelationVXchain1d_dis_N4096_W1KPM_M256_state600defaultEF0.000000JACKSON.dat
/tmp/agi_r2_t02/captured/0173_meanCorrelationVXchain1d_dis_N4096_W1KPM_M256_state600defaultEF1.000000JACKSON.dat
/tmp/agi_r2_t02/captured/0174_meanCorrelationVXchain1d_dis_N4096_W1KPM_M256_state600defaultEF-1.000000JACKSON.dat
/tmp/agi_r2_t02/captured/0175_meanCorrelationVXchain1d_dis_N4096_W1KPM_M256_state600defaultEF1.900000JACKSON.dat
/tmp/agi_r2_t02/captured/0176_meanCorrelationVXchain1d_dis_N4096_W1KPM_M256_state600defaultEF0.000000JACKSON.dat
/tmp/agi_r2_t02/captured/0177_meanCorrelationVXchain1d_dis_N4096_W1KPM_M256_state600defaultEF1.000000JACKSON.dat
/tmp/agi_r2_t02/captured/0178_meanCorrelationVXchain1d_dis_N4096_W1KPM_M256_state600defaultEF-1.000000JACKSON.dat
/tmp/agi_r2_t02/captured/0179_meanCorrelationVXchain1d_dis_N4096_W1KPM_M256_state600defaultEF1.900000JACKSON.dat
/tmp/agi_r2_t02/captured/0180_meanCorrelationVXchain1d_dis_N4096_W1KPM_M256_state600defaultEF0.000000JACKSON.dat
/tmp/agi_r2_t02/captured/0181_meanCorrelationVXchain1d_dis_N4096_W1KPM_M256_state600defaultEF1.000000JACKSON.dat
/tmp/agi_r2_t02/captured/0182_meanCorrelationVXchain1d_dis_N4096_W1KPM_M256_state600defaultEF-1.000000JACKSON.dat
/tmp/agi_r2_t02/captured/0183_meanCorrelationVXchain1d_dis_N4096_W1KPM_M256_state600defaultEF1.900000JACKSON.dat
/tmp/agi_r2_t02/captured/0184_meanCorrelationVXchain1d_dis_N4096_W1KPM_M256_state600defaultEF0.000000JACKSON.dat
/tmp/agi_r2_t02/captured/0185_meanCorrelationVXchain1d_dis_N4096_W1KPM_M256_state600defaultEF1.000000JACKSON.dat
/tmp/agi_r2_t02/captured/0186_meanCorrelationVXchain1d_dis_N4096_W1KPM_M256_state600defaultEF-1.000000JACKSON.dat
/tmp/agi_r2_t02/captured/0187_meanCorrelationVXchain1d_dis_N4096_W1KPM_M256_state600defaultEF1.900000JACKSON.dat
/tmp/agi_r2_t02/captured/0188_meanCorrelationVXchain1d_dis_N4096_W1KPM_M256_state600defaultEF0.000000JACKSON.dat
/tmp/agi_r2_t02/captured/0189_meanCorrelationVXchain1d_dis_N4096_W1KPM_M256_state600defaultEF1.000000JACKSON.dat
/tmp/agi_r2_t02/captured/0190_meanCorrelationVXchain1d_dis_N4096_W1KPM_M256_state600defaultEF-1.000000JACKSON.dat
/tmp/agi_r2_t02/captured/0191_meanCorrelationVXchain1d_dis_N4096_W1KPM_M256_state600defaultEF1.900000JACKSON.dat

$ OMP_NUM_THREADS=16 PYTHONPATH=build/lsquant_core:python python3 -u -m lsquant.apps.tutorial05 \
    --out /tmp/agi_r2_t05
/tmp/agi_r2_t05/meanCorrelationVXsquare_L200_W0.5KPM_M256_state150defaultEF2.000000JACKSON.dat
/tmp/agi_r2_t05/meanSpectralOp1square_L200_W0.5KPM_M1024_statedefaultJACKSON.dat
/tmp/agi_r2_t05/meanCorrelationVXsquare_L200_W1KPM_M256_state150defaultEF2.000000JACKSON.dat
/tmp/agi_r2_t05/meanSpectralOp1square_L200_W1KPM_M1024_statedefaultJACKSON.dat
/tmp/agi_r2_t05/meanCorrelationVXsquare_L200_W2KPM_M256_state150defaultEF2.000000JACKSON.dat
/tmp/agi_r2_t05/meanSpectralOp1square_L200_W2KPM_M1024_statedefaultJACKSON.dat
/tmp/agi_r2_t05/meanCorrelationVXsquare_L200_W4KPM_M256_state150defaultEF2.000000JACKSON.dat
/tmp/agi_r2_t05/meanSpectralOp1square_L200_W4KPM_M1024_statedefaultJACKSON.dat
```

---

## 6. Deviations from the Brief

1. **No `ZeroRegister` instruction needed for MSD seed initialisation.**
   In `ApplySeries`, an empty coefficient vector (`ApplySeries{{}}`) scales the active seed vector by zero (`linalg::scal(0.0, phi)`), matching the `linalg::scal(0, Phi)` pre-step of `MomentsTD::Evolve`. Because `MomentEngine::Run` initialises `seed[Aux] = |r>` at the start of each state, `msd_plan` prepares the initial MSD seeds (`seed[Ket] = 0`, `seed[Aux] = |r>`) with `ApplySeries([])` followed by `MirrorToBra()` without adding a separate `ZeroRegister` instruction.
2. **Optional `scale` scalar on `ApplySeries`, `ApplyDerivedSeries`, `EvolveRegister`, and `EvolveDerived`.**
   Legacy `MomentsTD::Evolve` and `MomentsTD::Conv_Evolve` accumulate `sum_n c_n T_n(H~) v` with pure real/imaginary Bessel coefficients `c_n` and apply the global band-center phase `exp0 = exp(-i BandCenter dt / hbar)` once at the end via `linalg::scal(exp0, Phi)`. Carrying `scale` (default `(1.0, 0.0)`) on the series instructions allows Python to pass `exp0` so the final `linalg::scal(scale, phi)` happens after the Chebyshev sum without introducing physical units or time parameters into `lsquant/core/`.
3. **G6 `EF` energies in `tests/dev/stage_4_check.py`.**
   The brief's text in section 2 listed `0, 0.5, -1.5, 1.9` while `goldens/stage0/PROVENANCE.md` (line 77) and the committed files in `goldens/stage0/chain1d_ballistic/` have `EF = 0.0, 0.5, 1.0, 1.5` (`meanCorrelationVXchain1dKPM_M256_state20exactEF{0.000000,0.500000,1.000000,1.500000}JACKSON.dat`). `tests/dev/stage_4_check.py` checks against the four committed G6 golden files `(0.0, 0.5, 1.0, 1.5)`.

---

## 7. Open Questions

None.
