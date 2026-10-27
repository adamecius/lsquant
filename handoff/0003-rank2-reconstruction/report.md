# Handoff 0003 Report — Stage 3: Rank-2 Reconstruction (Kubo-Greenwood, Kubo-Bastin, Sea/Surf), Graphene, and Haldane

## 1. Files Added and Modified

- `lsquant/core/algebra/sparse_matrix.hpp`, `lsquant/core/algebra/sparse_matrix.cpp` (modified for Item 0: added pointer-based `SparseMatrix::FromCsr` overload so CSR buffers are mapped into Eigen without intermediate `std::vector` staging while preserving the existing by-value overload)
- `lsquant/core/chebyshev/contraction.hpp`, `lsquant/core/chebyshev/contraction.cpp` (modified: added `ContractRank2Imag`, `ContractRank2Real`, `ContractRank2RealReal`)
- `lsquant/compat/legacy_kernels.hpp`, `lsquant/compat/legacy_kernels.cpp` (added: gate-only `greenR_chebF`, `DgreenR_chebF`, and `ComputeLegacyGreenTables` using `libstdc++` `std::pow(std::complex<double>, double)`)
- `lsquant/api/bindings.cpp` (modified: fixed `SparseMatrix.from_csr` for Item 0, bound `SpectralBounds`, `GershgorinBounds`, `LanczosBounds`, `ContractRank2Imag`, `ContractRank2Real`, `ContractRank2RealReal`, and `_core.compat.legacy_green_tables`)
- `lsquant/CMakeLists.txt` (modified: compiled `compat/legacy_kernels.cpp` into `_core`)
- `python/lsquant/_core.pyi` (modified: added stubs for the stage-3 bindings and `_core.compat`)
- `python/lsquant/spectral_functions.py` (modified: added `GreenRetarded`, `GreenRetardedDerivative`, `legacy_green_tables`, and `complex128` output support in `SpectralFunction.table`)
- `python/lsquant/kernels.py` (modified: added `jackson_damping_2d`)
- `python/lsquant/moments.py` (modified: generalised `kubo_plan` for distinct left/right operator indices and added `window_from_gershgorin`)
- `python/lsquant/reconstruct.py` (modified: added `kubo_greenwood`, `kubo_bastin`, `kubo_bastin_sea`, `kubo_bastin_surf`, `kubo_fft_grid`)
- `python/lsquant/__init__.py` (modified: exported stage-3 symbols)
- `python/lsquant/apps/stage3_kubo.py` (added: stage-3 gate CLI)
- `python/lsquant/apps/tutorial03.py` (added: Tutorial 03 Haldane Hall application)
- `tests/dev/stage_3_check.py` (added: stage-3 developer self-check)

---

## 2. Item 0: `SparseMatrix.from_csr` Buffer Mapping and Single Remaining Copy

`SparseMatrix.from_csr` in `lsquant/api/bindings.cpp` now validates each of `row_pointers`, `columns`, and `values` via `ValidateCsrBuffer1D` with `forcecast` disabled:
- `row_pointers` and `columns` must be 1D C-contiguous `int64` `numpy.ndarray` objects (`dtype.kind == 'i'`, `itemsize == 8`, `ndim == 1`, `c_style` flag set, element stride `8` bytes);
- `values` must be a 1D C-contiguous `complex128` `numpy.ndarray` (`dtype.kind == 'c'`, `itemsize == 16`, `ndim == 1`, `c_style` flag set, element stride `16` bytes);
- any violation raises a Python `TypeError` (`py::type_error`) explicitly naming the offending argument (`'row_pointers'`, `'columns'`, or `'values'`).

The validated buffer pointers are passed directly to `SparseMatrix::FromCsr(dim, rpSize, rpPtr, colSize, colPtr, valSize, valPtr)` without constructing intermediate `std::vector<Index>` or `std::vector<Scalar>` objects. The existing by-value `SparseMatrix::FromCsr(Index, std::vector<Index>, std::vector<Index>, std::vector<Scalar>)` signature is retained for C++ callers and delegates to the pointer overload.

**Single remaining copy:** `Eigen::SparseMatrix<Scalar, Eigen::RowMajor, Index>` owns its internal `CompressedStorage` buffers (`StorageIndex*` and `Scalar*` allocated by Eigen's aligned allocator) and cannot adopt external NumPy or `std::vector` memory via move. Assigning `Eigen::Map<const Eigen::SparseMatrix<Scalar, Eigen::RowMajor, Index>>(dim, dim, nnz, rowPointers, columns, values)` into `result.matrix_` performs that single copy from the mapped NumPy CSR buffers into Eigen's owned compressed storage.

---

## 3. Fresh CMake Configure and Build on TCN (`~/lsquant-agi-build`)

```
$ cd ~/lsquant-agi-build && rm -rf build
$ cmake -S . -B build -DLSQUANT_BUILD_NEW_CORE=ON -DLSQUANT_BUILD_PYTHON=ON -Dpybind11_DIR="$(python3 -c 'import pybind11; print(pybind11.get_cmake_dir())')"
$ cmake --build build --target lsquant_core _core stage1_moments stage_1_check -j 16
```

Verbatim configure and build output:

```
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
-- Found PythonInterp: /usr/bin/python3 (found suitable version "3.10.12", minimum required is "3.9")
-- Found PythonLibs: /usr/lib/x86_64-linux-gnu/libpython3.10.so
-- Performing Test HAS_FLTO_AUTO
-- Performing Test HAS_FLTO_AUTO - Success
-- Found pybind11: /home/ICN2/jgarcia/.local/lib/python3.10/site-packages/pybind11/include (found version "3.1.0")
-- wannier2sparse: /home/ICN2/jgarcia/lsquant-agi-build/.cache/wannier2sparse/build/wannier2sparse
-- Configuring done
-- Generating done
-- Build files have been written to: /home/ICN2/jgarcia/lsquant-agi-build/build
[  8%] Building CXX object lsquant_core/CMakeFiles/lsquant_core.dir/core/chebyshev/rescaled_operator.cpp.o
[ 16%] Building CXX object lsquant_core/CMakeFiles/lsquant_core.dir/core/algebra/spectral_bounds.cpp.o
[ 25%] Building CXX object lsquant_core/CMakeFiles/lsquant_core.dir/core/algebra/composite_operators.cpp.o
[ 25%] Building CXX object lsquant_core/CMakeFiles/lsquant_core.dir/core/algebra/sparse_matrix.cpp.o
[ 25%] Building CXX object lsquant_core/CMakeFiles/lsquant_core.dir/core/chebyshev/chebyshev_vectors.cpp.o
[ 41%] Building CXX object lsquant_core/CMakeFiles/lsquant_core.dir/core/algebra/linalg.cpp.o
[ 41%] Building CXX object lsquant_core/CMakeFiles/lsquant_core.dir/core/algebra/state_generator.cpp.o
[ 50%] Building CXX object lsquant_core/CMakeFiles/lsquant_core.dir/core/chebyshev/chebyshev_recursion.cpp.o
[ 58%] Building CXX object lsquant_core/CMakeFiles/lsquant_core.dir/core/chebyshev/moment_engine.cpp.o
[ 75%] Building CXX object lsquant_core/CMakeFiles/lsquant_core.dir/io/csr_text.cpp.o
[ 91%] Building CXX object lsquant_core/CMakeFiles/lsquant_core.dir/io/curve_text.cpp.o
[ 91%] Building CXX object lsquant_core/CMakeFiles/lsquant_core.dir/core/chebyshev/contraction.cpp.o
[ 91%] Building CXX object lsquant_core/CMakeFiles/lsquant_core.dir/core/chebyshev/chebyshev_moments.cpp.o
[100%] Building CXX object lsquant_core/CMakeFiles/lsquant_core.dir/io/chebmom_file.cpp.o
[100%] Linking CXX static library ../lib/liblsquant_core.a
[100%] Built target lsquant_core
[ 85%] Built target lsquant_core
[ 85%] Building CXX object lsquant_core/CMakeFiles/_core.dir/api/bindings.cpp.o
[ 92%] Building CXX object lsquant_core/CMakeFiles/_core.dir/compat/legacy_kernels.cpp.o
[100%] Linking CXX shared module _core.cpython-310-x86_64-linux-gnu.so
[100%] Built target _core
[ 92%] Built target lsquant_core
[ 92%] Building CXX object lsquant_core/CMakeFiles/stage1_moments.dir/apps/stage1_moments.cpp.o
[100%] Linking CXX executable stage1_moments
[100%] Built target stage1_moments
[ 85%] Built target lsquant_core
[ 92%] Building CXX object lsquant_core/CMakeFiles/stage_1_check.dir/__/tests/dev/stage_1_check.cpp.o
[100%] Linking CXX executable stage_1_check
[100%] Built target stage_1_check
```

Zero compiler warnings were emitted under `lsquant/` with `-Wall -Wextra`.

---

## 4. Operator Order Verification on G5 (`VX-VY`) and Gershgorin Window on G4

In legacy `inline_compute-kpm-nonEqOp <label> <OP_R> <OP_L> <M>`, the first CLI operator (`OP_R`, passed as `OPL` to `CorrelationExpansionMoments`) seeds the bra recursion (`chevVecL.SetInitVectors(OPL, PhiL)` -> `left_m = T_m OPL |r>`), and the second CLI operator (`OP_L`, passed as `OPR` to `CorrelationExpansionMoments`) is applied after the ket recursion (`chevVecR.IterateAll(); chevVecR.Multiply(OPR)` -> `right_n = OPR T_n |r>`), so `mu_{m,n} = <T_m OPL r | OPR T_n r>`.

`moments.kubo_plan(num_moments: int, *, left_index: int, right_index: int, word: str)` takes `left_index`, `right_index`, and `word` as required keyword-only parameters (no default indices, no argument-order heuristics, and no hard-coded operator names; callers pass `word = f"{left} delta {right} delta"`). `moments.dos_plan(num_moments: int, num_states: int, word: str)` likewise takes `word` from the caller (`"1 delta"`). `kubo_plan` emits:
- ket: `Expand(num_moments)`, `ApplyOperator(right_index)`
- `MirrorToBra()`
- bra: `ApplyOperator(left_index)`, `Expand(num_moments)`
- `Trace([1])`

Verified on G5 (`haldane` with `left_index = 0` for `VX`, `right_index = 1` for `VY`, and `word = "VX delta VY delta"`, `M = 128`, exact trace on `[-4.5, 4.5]`): the resulting `.chebmom2D` tensor matches `goldens/stage0/haldane/NonEqOpVX-VYhaldaneKPM_M128x128_stateexact.chebmom2D` with `numpy.array_equal == True` (and identical file bytes).

Verified `moments.window_from_gershgorin(H)` on `graphene_golden.HAM.CSR`: `_core.GershgorinBounds(H, 0.10)` produces `BandWidth = (upper - lower) * 1.0 = 15.0366847652085`, matching the G4 `.chebmom2D` header `256 15.0366847652085  0` to the last digit.

---

## 5. Developer Self-Check Output (`tests/dev/stage_3_check.py`)

```
$ OMP_NUM_THREADS=16 PYTHONPATH=build/lsquant_core:python python3 tests/dev/stage_3_check.py
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
```

---

## 6. Gate Program Runs (`python -m lsquant.apps.stage3_kubo`)

### Form 2: Compute Moments from CSR Operators

```
$ OMP_NUM_THREADS=16 PYTHONPATH=build/lsquant_core:python python3 -m lsquant.apps.stage3_kubo \
    --operators goldens/stage0/haldane/operators --label haldane --left VX --right VY \
    --window -4.5 4.5 --trace exact --states 512 --moments 128 --out /tmp/stage3_gate
/tmp/stage3_gate/NonEqOpVX-VYhaldaneKPM_M128x128_stateexact.chebmom2D
/tmp/stage3_gate/NonEqOpVX-VYhaldaneKPM_M128x128_stateexact.chebmom

$ OMP_NUM_THREADS=16 PYTHONPATH=build/lsquant_core:python python3 -m lsquant.apps.stage3_kubo \
    --operators goldens/stage0/haldane/operators --label haldane --left VX --right VX \
    --window -4.5 4.5 --trace exact --states 512 --moments 128 --out /tmp/stage3_gate
/tmp/stage3_gate/NonEqOpVX-VXhaldaneKPM_M128x128_stateexact.chebmom2D
/tmp/stage3_gate/NonEqOpVX-VXhaldaneKPM_M128x128_stateexact.chebmom

$ OMP_NUM_THREADS=16 PYTHONPATH=build/lsquant_core:python python3 -m lsquant.apps.stage3_kubo \
    --operators goldens/stage0/graphene/operators --label graphene_golden --left VX --right VX \
    --window gershgorin --trace random --seed 12345 --states 1 --moments 64 --out /tmp/stage3_gate
/tmp/stage3_gate/NonEqOpVX-VXgraphene_goldenKPM_M64x64_statedefault.chebmom2D
/tmp/stage3_gate/NonEqOpVX-VXgraphene_goldenKPM_M64x64_statedefault.chebmom
```

### Form 1: Reconstruct Each Formula from Moment File

```
$ OMP_NUM_THREADS=16 PYTHONPATH=build/lsquant_core:python python3 -m lsquant.apps.stage3_kubo \
    --moments goldens/stage0/haldane/NonEqOpVX-VYhaldaneKPM_M128x128_stateexact.chebmom2D \
    --formula bastin --broadening-meV 10 --out /tmp/stage3_gate
/tmp/stage3_gate/KuboBastin_NonEqOpVX-VYhaldaneKPM_M128x128_stateexactJACKSON.dat

$ OMP_NUM_THREADS=16 PYTHONPATH=build/lsquant_core:python python3 -m lsquant.apps.stage3_kubo \
    --moments goldens/stage0/haldane/NonEqOpVX-VYhaldaneKPM_M128x128_stateexact.chebmom2D \
    --formula sea --broadening-meV 10 --out /tmp/stage3_gate
/tmp/stage3_gate/KuboBastinSea_NonEqOpVX-VYhaldaneKPM_M128x128_stateexactJACKSON.dat

$ OMP_NUM_THREADS=16 PYTHONPATH=build/lsquant_core:python python3 -m lsquant.apps.stage3_kubo \
    --moments goldens/stage0/haldane/NonEqOpVX-VYhaldaneKPM_M128x128_stateexact.chebmom2D \
    --formula surf --broadening-meV 10 --out /tmp/stage3_gate
/tmp/stage3_gate/KuboBastinSurf_NonEqOpVX-VYhaldaneKPM_M128x128_stateexactJACKSON.dat

$ OMP_NUM_THREADS=16 PYTHONPATH=build/lsquant_core:python python3 -m lsquant.apps.stage3_kubo \
    --moments goldens/stage0/haldane/NonEqOpVX-VXhaldaneKPM_M128x128_stateexact.chebmom2D \
    --formula greenwood --broadening-meV 10 --out /tmp/stage3_gate
/tmp/stage3_gate/KuboGreenWood_NonEqOpVX-VXhaldaneKPM_M128x128_stateexactJACKSON.dat

$ OMP_NUM_THREADS=16 PYTHONPATH=build/lsquant_core:python python3 -m lsquant.apps.stage3_kubo \
    --moments goldens/stage0/graphene/NonEqOpVX-VXgraphene_goldenKPM_M64x64_statedefault.chebmom2D \
    --formula bastin --broadening-meV 64 --out /tmp/stage3_gate
/tmp/stage3_gate/KuboBastin_NonEqOpVX-VXgraphene_goldenKPM_M64x64_statedefaultJACKSON.dat

$ OMP_NUM_THREADS=16 PYTHONPATH=build/lsquant_core:python python3 -m lsquant.apps.stage3_kubo \
    --moments goldens/stage0/graphene/NonEqOpVX-VXgraphene_goldenKPM_M64x64_statedefault.chebmom2D \
    --formula bastin_fft --broadening-meV 10 --out /tmp/stage3_gate
/tmp/stage3_gate/KuboBastin_NonEqOpVX-VXgraphene_goldenKPM_M64x64_statedefaultJACKSON_FFTgrid.dat

$ OMP_NUM_THREADS=16 PYTHONPATH=build/lsquant_core:python python3 -m lsquant.apps.stage3_kubo \
    --moments goldens/stage0/graphene/NonEqOpVX-VXgraphene_goldenKPM_M64x64_statedefault.chebmom2D \
    --formula greenwood_fft --broadening-meV 10 --out /tmp/stage3_gate
/tmp/stage3_gate/KuboGreenWood_NonEqOpVX-VXgraphene_goldenKPM_M64x64_statedefaultJACKSON_FFTgrid.dat
```

---

## 7. Tutorial 03 Run (`python -m lsquant.apps.tutorial03`)

`tutorial03` accepts `--operators` defaulting to `goldens/stage0/haldane/operators` (resolved from the repository root, using `make_haldane` only for the `.desc` sidecar and `exact` state file, without reading `test/golden`). To reproduce the legacy `lsqhall.py` pipeline artefact (where `.chebmom2D` was written at 15-digit `digits10` ASCII precision and re-read before reconstruction), `tutorial03` (and the `--formula` branch of `stage3_kubo` form 2) writes the `.chebmom2D` via `_core.WriteLegacyChebmom2D`, reads it back via `_core.ReadChebmom`, and reconstructs the curves from the re-read tensor.

```
$ OMP_NUM_THREADS=16 PYTHONPATH=build/lsquant_core:python python3 -m lsquant.apps.tutorial03 --out /tmp/stage3_tut03
/tmp/stage3_tut03/KuboBastin_NonEqOpVX-VYhaldaneKPM_M128x128_stateexactJACKSON.dat
/tmp/stage3_tut03/KuboBastinSea_NonEqOpVX-VYhaldaneKPM_M128x128_stateexactJACKSON.dat
/tmp/stage3_tut03/KuboBastinSurf_NonEqOpVX-VYhaldaneKPM_M128x128_stateexactJACKSON.dat
/tmp/stage3_tut03/KuboGreenWood_NonEqOpVX-VXhaldaneKPM_M128x128_stateexactJACKSON.dat
Chern from gap plateau: C = +1.017 (+1.017016743089)
```

---

## 8. Deviations from the Brief

1. **Both `numpy.ndarray` and `ChebyshevMoments` overloads for `ContractRank2Imag`, `ContractRank2Real`, and `ContractRank2RealReal`.**
   As in stage 2 (`ContractRank1`), `reconstruct.py` truncates and applies 2D Jackson damping to `moments.data` in Python before calling the C++ contraction kernels. Binding both a 2D `complex128` `numpy.ndarray` overload and a `ChebyshevMoments` overload for each of the three rank-2 contraction functions allows callers to pass either the damped 2D NumPy array directly or a `ChebyshevMoments` instance.

---

## 9. Open Questions

None.
