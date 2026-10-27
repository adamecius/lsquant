# Handoff 0001 Report (Agi) — Stage 1: Core Algebra and Chebyshev Engine

## 1. Files Added and Modified

### Added Files
- `lsquant/CMakeLists.txt`
- `lsquant/core/algebra/types.hpp`
- `lsquant/core/algebra/linear_operator.hpp`
- `lsquant/core/algebra/sparse_matrix.hpp`
- `lsquant/core/algebra/sparse_matrix.cpp`
- `lsquant/core/algebra/composite_operators.hpp`
- `lsquant/core/algebra/composite_operators.cpp`
- `lsquant/core/algebra/linalg.hpp`
- `lsquant/core/algebra/linalg.cpp`
- `lsquant/core/algebra/state_generator.hpp`
- `lsquant/core/algebra/state_generator.cpp`
- `lsquant/core/algebra/spectral_bounds.hpp`
- `lsquant/core/algebra/spectral_bounds.cpp`
- `lsquant/core/chebyshev/spectral_window.hpp`
- `lsquant/core/chebyshev/rescaled_operator.hpp`
- `lsquant/core/chebyshev/rescaled_operator.cpp`
- `lsquant/core/chebyshev/chebyshev_recursion.hpp`
- `lsquant/core/chebyshev/chebyshev_recursion.cpp`
- `lsquant/core/chebyshev/chebyshev_vectors.hpp`
- `lsquant/core/chebyshev/chebyshev_vectors.cpp`
- `lsquant/core/chebyshev/chebyshev_moments.hpp`
- `lsquant/core/chebyshev/chebyshev_moments.cpp`
- `lsquant/core/chebyshev/plan.hpp`
- `lsquant/core/chebyshev/moment_engine.hpp`
- `lsquant/core/chebyshev/moment_engine.cpp`
- `lsquant/io/csr_text.hpp`
- `lsquant/io/csr_text.cpp`
- `lsquant/io/chebmom_file.hpp`
- `lsquant/io/chebmom_file.cpp`
- `lsquant/apps/stage1_moments.cpp`
- `tests/dev/stage_1_check.cpp`

### Modified Files
- `CMakeLists.txt` (`option(LSQUANT_BUILD_NEW_CORE "" ON)` and guarded `add_subdirectory(lsquant lsquant_core)`)

---

## 2. Build Commands and Output on TCN (`fatnode-tacn`)

### Commands
```bash
rsync -a --delete --exclude .git --exclude build /root/lsquant-agi/ jh-x:~/lsquant-agi-build/
ssh jh-x "cd ~/lsquant-agi-build && rm -rf build && cmake -S . -B build -DLSQUANT_BUILD_NEW_CORE=ON && cmake --build build -j 16 && ls build/lsquant_core"
```

### Output
```text
-- The CXX compiler identification is GNU 11.4.0
-- Detecting CXX compiler ABI info
-- Detecting CXX compiler ABI info - done
-- Check for working CXX compiler: /usr/bin/c++ - skipped
-- Detecting CXX compile features
-- Detecting CXX compile features - done
-- Found OpenMP_CXX: -fopenmp (found version "4.5")
-- Found OpenMP: TRUE (found version "4.5")
-- Found PkgConfig: /usr/bin/pkg-config (found version "0.29.2")
-- Found GSL: /usr/include (found version "2.7.1")
-- Looking for C++ include pthread.h
-- Looking for C++ include pthread.h - found
-- Performing Test CMAKE_HAVE_LIBC_PTHREAD
-- Performing Test CMAKE_HAVE_LIBC_PTHREAD - Success
-- Found Threads: TRUE
-- Checking for module 'fftw3'
--   Found fftw3, version 3.3.8
-- wannier2sparse: /home/ICN2/jgarcia/lsquant-agi-build/.cache/wannier2sparse/build/wannier2sparse
-- Configuring done
-- Generating done
-- Build files have been written to: /home/ICN2/jgarcia/lsquant-agi-build/build
[ 11%] Building CXX object lsquant_core/CMakeFiles/lsquant_core.dir/core/algebra/sparse_matrix.cpp.o
[ 16%] Building CXX object lsquant_core/CMakeFiles/lsquant_core.dir/core/algebra/composite_operators.cpp.o
[ 21%] Building CXX object lsquant_core/CMakeFiles/lsquant_core.dir/core/algebra/linalg.cpp.o
[ 24%] Building CXX object lsquant_core/CMakeFiles/lsquant_core.dir/core/algebra/state_generator.cpp.o
[ 26%] Building CXX object lsquant_core/CMakeFiles/lsquant_core.dir/core/algebra/spectral_bounds.cpp.o
[ 26%] Building CXX object lsquant_core/CMakeFiles/lsquant_core.dir/core/chebyshev/rescaled_operator.cpp.o
[ 27%] Building CXX object lsquant_core/CMakeFiles/lsquant_core.dir/core/chebyshev/chebyshev_recursion.cpp.o
[ 30%] Building CXX object lsquant_core/CMakeFiles/lsquant_core.dir/core/chebyshev/chebyshev_vectors.cpp.o
[ 32%] Building CXX object lsquant_core/CMakeFiles/lsquant_core.dir/core/chebyshev/chebyshev_moments.cpp.o
[ 32%] Building CXX object lsquant_core/CMakeFiles/lsquant_core.dir/core/chebyshev/moment_engine.cpp.o
[ 34%] Building CXX object lsquant_core/CMakeFiles/lsquant_core.dir/io/csr_text.cpp.o
[ 34%] Building CXX object lsquant_core/CMakeFiles/lsquant_core.dir/io/chebmom_file.cpp.o
[ 37%] Linking CXX static library ../lib/liblsquant_core.a
[ 37%] Built target lsquant_core
[ 76%] Building CXX object lsquant_core/CMakeFiles/stage1_moments.dir/apps/stage1_moments.cpp.o
[ 79%] Building CXX object lsquant_core/CMakeFiles/stage_1_check.dir/__/tests/dev/stage_1_check.cpp.o
[ 88%] Linking CXX executable stage1_moments
[ 88%] Built target stage1_moments
[ 98%] Linking CXX executable stage_1_check
[ 98%] Built target stage_1_check
[100%] Linking CXX executable ../lsquant
[100%] Built target lsquant
```
Compiler warnings under `lsquant/` and `tests/dev/stage_1_check.cpp` with `-Wall -Wextra`: 0.

### `ls build/lsquant_core` (fresh build from scratch)
```text
CMakeFiles
cmake_install.cmake
Makefile
stage_1_check
stage1_moments
```
(`ls -la build/lsquant_core` confirms both executables are regular files produced directly by CMake in `build/lsquant_core/`.)

### Legacy CTest Summary (`ctest --test-dir build --output-on-failure` with `LSQUANT_BUILD_NEW_CORE=ON`)
```text
Test project /home/ICN2/jgarcia/lsquant-agi-build/build
      Start  1: chain1d_analytic
 1/31 Test  #1: chain1d_analytic .................   Passed    0.01 sec
      Start  2: chain1d_dos
 2/31 Test  #2: chain1d_dos ......................   Passed    0.01 sec
      Start  3: dos_reconstructed
 3/31 Test  #3: dos_reconstructed ................   Passed    0.00 sec
      Start  4: chain1d_reproducible
 4/31 Test  #4: chain1d_reproducible .............   Passed   22.77 sec
      Start  5: spin_algebra
 5/31 Test  #5: spin_algebra .....................   Passed    0.02 sec
      Start  6: graphene_regression
 6/31 Test  #6: graphene_regression ..............   Passed    0.13 sec
      Start  7: reconstruction_nan_guard
 7/31 Test  #7: reconstruction_nan_guard .........   Passed    0.20 sec
      Start  8: msd_ballistic
 8/31 Test  #8: msd_ballistic ....................   Passed   19.69 sec
      Start  9: hall_haldane
 9/31 Test  #9: hall_haldane .....................   Passed    8.46 sec
      Start 10: spin_precession
10/31 Test #10: spin_precession ..................   Passed    7.94 sec
      Start 11: recon_kernel_oracle
11/31 Test #11: recon_kernel_oracle ..............   Passed    0.10 sec
      Start 12: spmv_threaded_bitexact
12/31 Test #12: spmv_threaded_bitexact ...........   Passed    0.01 sec
      Start 13: spmm_batchmultiply_equiv
13/31 Test #13: spmm_batchmultiply_equiv .........   Passed    0.01 sec
      Start 14: greenwood_regression
14/31 Test #14: greenwood_regression .............   Passed    0.05 sec
      Start 15: safety_factors_guard
15/31 Test #15: safety_factors_guard .............   Passed    0.37 sec
      Start 16: descriptor_conformance
16/31 Test #16: descriptor_conformance ...........   Passed    0.00 sec
      Start 17: descriptor_bounds
17/31 Test #17: descriptor_bounds ................   Passed    0.15 sec
      Start 18: tb_descriptor
18/31 Test #18: tb_descriptor ....................   Passed    0.00 sec
      Start 19: config_equivalence
19/31 Test #19: config_equivalence ...............   Passed    0.26 sec
      Start 20: corrop_state_honored
20/31 Test #20: corrop_state_honored .............   Passed    2.27 sec
      Start 21: fft_grid_regression
21/31 Test #21: fft_grid_regression ..............   Passed    0.09 sec
      Start 22: sea_surf_regression
22/31 Test #22: sea_surf_regression ..............   Passed    0.33 sec
      Start 23: jackson_kernel
23/31 Test #23: jackson_kernel ...................   Passed    0.00 sec
      Start 24: jackson_g_scalar
24/31 Test #24: jackson_g_scalar .................   Passed    0.00 sec
      Start 25: phase6_wired
25/31 Test #25: phase6_wired .....................   Passed    0.00 sec
      Start 26: single_binary
26/31 Test #26: single_binary ....................   Passed    0.69 sec
      Start 27: reporting_smoke
27/31 Test #27: reporting_smoke ..................   Passed    0.08 sec
      Start 28: unified_io_smoke
28/31 Test #28: unified_io_smoke .................   Passed    0.24 sec
      Start 29: example_dos_smoke
29/31 Test #29: example_dos_smoke ................   Passed    0.08 sec
      Start 30: example_localization_smoke
30/31 Test #30: example_localization_smoke .......   Passed    0.29 sec
      Start 31: example_haldane_smoke
31/31 Test #31: example_haldane_smoke ............   Passed   16.71 sec

100% tests passed, 0 tests failed out of 31

Total Test time (real) =  80.98 sec
```

---

## 3. `stage_1_check` Output at 1 and 8 Threads on TCN

### Commands
```bash
OMP_NUM_THREADS=1 build/lsquant_core/stage_1_check
OMP_NUM_THREADS=8 build/lsquant_core/stage_1_check
```

### Output (`OMP_NUM_THREADS=1`)
```text
dos max_abs_error: 2.47025e-15
kubo max_abs_error: 1.13687e-13
OMP_NUM_THREADS=1 vs 8 bit-identical: YES
```

### Output (`OMP_NUM_THREADS=8`)
```text
dos max_abs_error: 2.47025e-15
kubo max_abs_error: 1.13687e-13
OMP_NUM_THREADS=1 vs 8 bit-identical: YES
```

---

## 4. `stage1_moments` Output on `test/golden/chain1d/operators`

(`goldens/stage0` is not present in this branch yet; executed against `test/golden/chain1d/operators` with `--window -2 2 --trace exact --seed 12345 --states 2048 --moments 64`, plus the `--word kubo --batch 16` run into `/tmp/stage1_out_b16`.)

### Commands
```bash
OMP_NUM_THREADS=8 build/lsquant_core/stage1_moments \
  --operators test/golden/chain1d/operators --label chain1d --window -2 2 \
  --trace exact --seed 12345 --states 2048 --moments 64 --word dos --out /tmp/stage1_out

OMP_NUM_THREADS=8 build/lsquant_core/stage1_moments \
  --operators test/golden/chain1d/operators --label chain1d --window -2 2 \
  --trace exact --seed 12345 --states 2048 --moments 64 --word kubo --out /tmp/stage1_out

OMP_NUM_THREADS=8 build/lsquant_core/stage1_moments \
  --operators test/golden/chain1d/operators --label chain1d --window -2 2 \
  --trace exact --seed 12345 --states 2048 --moments 64 --word kubo --batch 16 --out /tmp/stage1_out_b16
```

### Standard Output
```text
stage1_moments: word=dos dim=2048 states=2048 (exact) M=64 -> /tmp/stage1_out/SpectralOp1chain1dKPM_M64_stateexact.chebmom1D
stage1_moments: word=kubo dim=2048 states=2048 (exact) M=64x64 -> /tmp/stage1_out/NonEqOpVX-VXchain1dKPM_M64x64_stateexact.chebmom2D
stage1_moments: word=kubo dim=2048 states=2048 (exact) M=64x64 -> /tmp/stage1_out_b16/NonEqOpVX-VXchain1dKPM_M64x64_stateexact.chebmom2D
```

### Leading Lines of Produced Files
```text
==> /tmp/stage1_out/SpectralOp1chain1dKPM_M64_stateexact.chebmom1D <==
2048 4 0
64
1 0
0 0
0 0
0 0
0 0
0 0
0 0
0 0

==> /tmp/stage1_out/SpectralOp1chain1dKPM_M64_stateexact.chebmom <==
chebmom 1
rank 1
dims 64
window -2 2
system_size 2048
states exact 2048
word 1 delta
kpm_weighted 1
symmetrised 0
data

==> /tmp/stage1_out/NonEqOpVX-VXchain1dKPM_M64x64_stateexact.chebmom2D <==
2048 4  0
64 64 
2 -0
0 0
-2 0
0 0
0 0
0 0
0 0
0 0

==> /tmp/stage1_out/NonEqOpVX-VXchain1dKPM_M64x64_stateexact.chebmom <==
chebmom 1
rank 2
dims 64 64
window -2 2
system_size 2048
states exact 2048
word VX delta VX delta
kpm_weighted 1
symmetrised 1
data

==> /tmp/stage1_out_b16/NonEqOpVX-VXchain1dKPM_M64x64_stateexact.chebmom2D <==
2048 4  0
64 64 
2 -0
0 0
-2 0
0 0
0 0
0 0
0 0
0 0
0 0
0 0

==> /tmp/stage1_out_b16/NonEqOpVX-VXchain1dKPM_M64x64_stateexact.chebmom <==
chebmom 1
rank 2
dims 64 64
window -2 2
system_size 2048
states exact 2048
word VX delta VX delta
kpm_weighted 1
symmetrised 1
data
2 -0
0 0
```

---

## 5. Legacy Numerical Conventions Reproduced and Deviations

### Legacy Numerical Conventions Reproduced for Byte Identity
1. **Hermiticity Fold in `ReadCsrText` via `SparseMatrix::HermitianFold()` (`0.5 * (transp + matrix_)`):** Legacy `SparseMatrixType::ConvertFromCSR` folds CSR file inputs unconditionally via `matrix_ = 0.5 * (transp + matrix_)`. In `.CSR` files such as `chain1d.VX.CSR`, off-diagonal entries contain `-0` real components (`-0 -1`) whose conjugate sum with `+0 1` converts `-0.0` into `+0.0`. `SparseMatrix::FromCsr` and `SparseMatrix::FromTriplets` preserve their inputs untouched without folding, while `ReadCsrText` calls `SparseMatrix::HermitianFold()` unconditionally as the legacy reader does.
2. **Real Scalar Weight Multiplication in `MomentEngine::Run` (`Trace`):** Legacy `SpectralMoments` evaluates `scal * linalg::vdot(...)` where `scal` is a `double`. Multiplying `double * std::complex<double>` scales the real and imaginary parts independently (`(scal * re, scal * im)`), whereas multiplying `std::complex<double>(scal, 0.0) * std::complex<double>(re, im)` evaluates `scal * im + 0.0 * re`, turning `-0.0` imaginary parts into `+0.0`. `MomentEngine::Run` checks `weight.imag() == 0.0` and uses `weight.real() * d` so IEEE-754 signed zeros match legacy `vdot` scaling and raw accumulation (`2 -0` in `kubo`).
3. **Legacy `.chebmom2D` Header Whitespace:** `WriteLegacyChebmom2D` emits two spaces between `BandWidth` and `BandCenter` (`"2048 4  0"`) and a trailing space after each dimension on line 2 (`"64 64 "`), matching `chebyshev::Moments2D::saveIn`.

### Deviations from the Brief
1. **Binary Subdirectory `add_subdirectory(lsquant lsquant_core)` and Target `RUNTIME_OUTPUT_DIRECTORY` (and Review-1 Symlink Episode):** To avoid a path collision between the legacy `build/lsquant` executable and the `lsquant/` subdirectory build tree, the root `CMakeLists.txt` uses `add_subdirectory(lsquant lsquant_core)` instead of `add_subdirectory(lsquant)`. Because the root `CMakeLists.txt` sets `CMAKE_RUNTIME_OUTPUT_DIRECTORY` to `${CMAKE_BINARY_DIR}`, removing the directory-wide override in review 1 caused CMake to link `stage1_moments` and `stage_1_check` into `build/` rather than `build/lsquant_core/`. During review 1, Agi created manual symlinks in `~/lsquant-agi-build/build/lsquant_core/` pointing to `../stage1_moments` and `../stage_1_check` on TCN instead of setting the target property in CMake, so a clean build did not match the paths stated in the report. This is now fixed in `lsquant/CMakeLists.txt` via `set_target_properties(<target> PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}")` on both `stage1_moments` and `stage_1_check`; all manual symlinks and the old `build/` directory on TCN were deleted, and a fresh configure and build from scratch places both executables directly in `build/lsquant_core/`.
2. **Const Accessors on `ChebyshevMoments`:** In addition to `std::vector<Scalar>& Data()`, `const std::vector<Scalar>& Data() const` was added alongside the `const` overloads of `At(...)` so `WriteChebmom`, `WriteLegacyChebmom1D`, and `WriteLegacyChebmom2D` can take `const ChebyshevMoments&`.

---

## 6. Open Questions
- None for Stage 1.
