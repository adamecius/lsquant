# Handoff 0002 Report — Stage 2: Python Library, pybind11 Bindings, SymPy Spectral Functions, and Rank-1 Reconstruction

## 1. Files Added and Modified

- `lsquant/core/chebyshev/contraction.hpp`, `lsquant/core/chebyshev/contraction.cpp`
- `lsquant/io/curve_text.hpp`, `lsquant/io/curve_text.cpp`
- `lsquant/api/bindings.cpp`
- `lsquant/CMakeLists.txt` (modified: added `LSQUANT_BUILD_PYTHON` option, `POSITION_INDEPENDENT_CODE ON` on `lsquant_core`, `contraction.cpp`, `curve_text.cpp`, `find_package(pybind11 CONFIG REQUIRED)`, and `pybind11_add_module(_core MODULE api/bindings.cpp)`)
- `python/pyproject.toml`
- `python/lsquant/__init__.py`
- `python/lsquant/_core.pyi`
- `python/lsquant/kernels.py`
- `python/lsquant/spectral_functions.py`
- `python/lsquant/moments.py`
- `python/lsquant/reconstruct.py`
- `python/lsquant/apps/__init__.py`
- `python/lsquant/apps/stage2_dos.py`
- `python/lsquant/apps/tutorial01.py`
- `tests/dev/stage_2_check.py`

---

## 2. Python Dependencies and CMake Configuration on TCN (`~/lsquant-agi-build`)

### `pip` Installation and Package Versions

```
$ pip install --user pybind11 sympy
Successfully installed mpmath-1.3.0 pybind11-3.1.0 sympy-1.14.0

$ python3 -c "import pybind11, sympy, numpy; print('pybind11:', pybind11.__version__, pybind11.get_cmake_dir()); print('sympy:', sympy.__version__); print('numpy:', numpy.__version__)"
pybind11: 3.1.0 /home/ICN2/jgarcia/.local/lib/python3.10/site-packages/pybind11/share/cmake/pybind11
sympy: 1.14.0
numpy: 2.2.6
```

### Fresh CMake Configure and Build on TCN

```
$ cd ~/lsquant-agi-build && rm -rf build
$ cmake -S . -B build -DLSQUANT_BUILD_NEW_CORE=ON -DLSQUANT_BUILD_PYTHON=ON -Dpybind11_DIR="$(python3 -c 'import pybind11; print(pybind11.get_cmake_dir())')"
$ cmake --build build -j 32
```

Verbatim configure and build output for the `lsquant_core` and `_core` targets:

```
-- The CXX compiler identification is GNU 11.4.0
-- Found OpenMP_CXX: -fopenmp (found version "4.5")
-- Found OpenMP: TRUE (found version "4.5")
-- Found PkgConfig: /usr/bin/pkg-config (found version "0.29.2")
-- Found GSL: /usr/include (found version "2.7.1")
-- Found Threads: TRUE
-- Checking for module 'fftw3'
--   Found fftw3, version 3.3.8
-- Found PythonInterp: /usr/bin/python3 (found suitable version "3.10.12", minimum required is "3.9")
-- Found PythonLibs: /usr/lib/x86_64-linux-gnu/libpython3.10.so
-- Performing Test HAS_FLTO_AUTO
-- Performing Test HAS_FLTO_AUTO - Success
-- Found pybind11: /home/ICN2/jgarcia/.local/lib/python3.10/site-packages/pybind11/include (found version "3.1.0")
-- wannier2sparse: /home/ICN2/jgarcia/lsquant-agi-build/.cache/wannier2sparse/build/wannier2sparse
-- Build files have been written to: /home/ICN2/jgarcia/lsquant-agi-build/build
[  2%] Building CXX object lsquant_core/CMakeFiles/lsquant_core.dir/core/algebra/spectral_bounds.cpp.o
[  2%] Building CXX object lsquant_core/CMakeFiles/lsquant_core.dir/core/algebra/composite_operators.cpp.o
[  5%] Building CXX object lsquant_core/CMakeFiles/lsquant_core.dir/core/chebyshev/moment_engine.cpp.o
[  6%] Building CXX object lsquant_core/CMakeFiles/lsquant_core.dir/io/chebmom_file.cpp.o
[  6%] Building CXX object lsquant_core/CMakeFiles/lsquant_core.dir/io/curve_text.cpp.o
[  7%] Building CXX object lsquant_core/CMakeFiles/lsquant_core.dir/core/algebra/state_generator.cpp.o
[  7%] Building CXX object lsquant_core/CMakeFiles/lsquant_core.dir/core/algebra/sparse_matrix.cpp.o
[  9%] Building CXX object lsquant_core/CMakeFiles/lsquant_core.dir/core/algebra/linalg.cpp.o
[  9%] Building CXX object lsquant_core/CMakeFiles/lsquant_core.dir/core/chebyshev/chebyshev_recursion.cpp.o
[ 10%] Building CXX object lsquant_core/CMakeFiles/lsquant_core.dir/core/chebyshev/chebyshev_vectors.cpp.o
[ 11%] Building CXX object lsquant_core/CMakeFiles/lsquant_core.dir/core/chebyshev/rescaled_operator.cpp.o
[ 15%] Building CXX object lsquant_core/CMakeFiles/lsquant_core.dir/core/chebyshev/contraction.cpp.o
[ 17%] Building CXX object lsquant_core/CMakeFiles/lsquant_core.dir/core/chebyshev/chebyshev_moments.cpp.o
[ 20%] Building CXX object lsquant_core/CMakeFiles/lsquant_core.dir/io/csr_text.cpp.o
[ 57%] Linking CXX static library ../lib/liblsquant_core.a
[ 57%] Built target lsquant_core
[ 58%] Building CXX object lsquant_core/CMakeFiles/stage1_moments.dir/apps/stage1_moments.cpp.o
[ 58%] Building CXX object lsquant_core/CMakeFiles/_core.dir/api/bindings.cpp.o
[ 59%] Building CXX object lsquant_core/CMakeFiles/stage_1_check.dir/__/tests/dev/stage_1_check.cpp.o
[ 60%] Linking CXX executable stage1_moments
[ 60%] Built target stage1_moments
[ 97%] Linking CXX shared module _core.cpython-310-x86_64-linux-gnu.so
[100%] Linking CXX executable stage_1_check
[100%] Built target stage_1_check
[100%] Built target _core
```

Zero compiler warnings were emitted under `lsquant/` with `-Wall -Wextra`.

### `ls build/lsquant_core` and `_core*.so` Path

```
$ ls -la build/lsquant_core
total 1048
drwxr-xr-x 3 jgarcia Theoretical   4096 Oct 20 14:54 .
drwxr-xr-x 8 jgarcia Theoretical   4096 Oct 20 14:54 ..
drwxr-xr-x 6 jgarcia Theoretical   4096 Oct 20 14:53 CMakeFiles
-rw-r--r-- 1 jgarcia Theoretical   1265 Oct 20 14:53 cmake_install.cmake
-rwxr-xr-x 1 jgarcia Theoretical 566344 Oct 20 14:54 _core.cpython-310-x86_64-linux-gnu.so
-rw-r--r-- 1 jgarcia Theoretical  30174 Oct 20 14:53 Makefile
-rwxr-xr-x 1 jgarcia Theoretical 275520 Oct 20 14:54 stage_1_check
-rwxr-xr-x 1 jgarcia Theoretical 174880 Oct 20 14:53 stage1_moments
```

Module path: `build/lsquant_core/_core.cpython-310-x86_64-linux-gnu.so`.

---

## 3. Developer Self-Check Output (`tests/dev/stage_2_check.py`)

```
$ PYTHONPATH=build/lsquant_core:python python3 tests/dev/stage_2_check.py
(a) dos_plan array_equal vs G1:  True
(a) kubo_plan array_equal vs G2: True
(b) DeltaFunction.table vs legacy_delta_table max_abs_diff = 2.220446e-16
(c) dos_curve vs G3 differing lines: 0 (total lines: 1590)
```

---

## 4. Gate Program Runs (`python -m lsquant.apps.stage2_dos`)

### Form 1: Reconstruct from Moment File

```
$ PYTHONPATH=build/lsquant_core:python python3 -m lsquant.apps.stage2_dos \
    --moments goldens/stage0/chain1d/SpectralOp1chain1dKPM_M64_stateexact.chebmom1D \
    --broadening-meV 120 --out /tmp/stage2_out_a
/tmp/stage2_out_a/meanSpectralOp1chain1dKPM_M64_stateexactJACKSON.dat
```

### Form 2: Compute Moments from CSR Operators and Reconstruct

```
$ PYTHONPATH=build/lsquant_core:python python3 -m lsquant.apps.stage2_dos \
    --operators goldens/stage0/chain1d/operators --label chain1d --window -2 2 \
    --trace exact --states 2048 --moments 64 --broadening-meV 120 --out /tmp/stage2_out_b
/tmp/stage2_out_b/SpectralOp1chain1dKPM_M64_stateexact.chebmom1D
/tmp/stage2_out_b/meanSpectralOp1chain1dKPM_M64_stateexactJACKSON.dat
```

---

## 5. Tutorial 01 Run (`python -m lsquant.apps.tutorial01`)

```
$ OMP_NUM_THREADS=1 PYTHONPATH=build/lsquant_core:python python3 -m lsquant.apps.tutorial01 --out /tmp/stage2_tut01
/tmp/stage2_tut01/SpectralOp1chain1d_N256KPM_M512_stateexact_chain1d_N256.chebmom1D
/tmp/stage2_tut01/moments_chain1d_N256_M512.chebmom1D
/tmp/stage2_tut01/SpectralOp1chain1d_N1024KPM_M512_stateexact_chain1d_N1024.chebmom1D
/tmp/stage2_tut01/moments_chain1d_N1024_M512.chebmom1D
/tmp/stage2_tut01/SpectralOp1chain1d_N4096KPM_M512_stateexact_chain1d_N4096.chebmom1D
/tmp/stage2_tut01/moments_chain1d_N4096_M512.chebmom1D
/tmp/stage2_tut01/eta20/meanmoments_chain1d_N256_M512JACKSON.dat
/tmp/stage2_tut01/eta20/meanmoments_chain1d_N1024_M512JACKSON.dat
/tmp/stage2_tut01/eta20/meanmoments_chain1d_N4096_M512JACKSON.dat
/tmp/stage2_tut01/eta80/meanmoments_chain1d_N256_M512JACKSON.dat
/tmp/stage2_tut01/eta24/meanmoments_chain1d_N256_M512JACKSON.dat
/tmp/stage2_tut01/eta10/meanmoments_chain1d_N256_M512JACKSON.dat
```

---

## 6. Deviations from the Brief

1. **Both real-vector and `ChebyshevMoments` overloads for `ContractRank1`.**
   Section 3 of the brief writes `ContractRank1(moments, table: float64[M, N]) -> float64[N]`, while section 4 specifies `_core.ContractRank1(moments_real, table)` taking a real `float64[M]` moment vector (specifically the `.real` view of the damped complex moments). To satisfy both sections without copying buffers, `contraction.hpp` and `bindings.cpp` provide both overloads: one accepting any 1D strided `float64` array (so `damped.real` with 16-byte stride is read directly in place without copying) and one accepting a rank-1 `ChebyshevMoments` instance directly.
2. **Disambiguation of the four `N = 256` broadening curves in `tutorial01.py`.**
   In `examples/01_dos_1d_chain/lsqplot.py`, `moments_chain1d_N256_M512.chebmom1D` is reconstructed four times (at `20, 80, 24, 10` meV) inside separate `tempfile.TemporaryDirectory` contexts, so all four produce the legacy basename `meanmoments_chain1d_N256_M512JACKSON.dat`. `tutorial01.py` writes the moment files at the top of `<out>` and places the reconstructed curves under `<out>/eta<meV>/meanmoments_chain1d_N<N>_M512JACKSON.dat` so that all four `N = 256` curves keep their legacy basename without overwriting each other.

---

## 7. Open Questions

1. **Stage-1 comment lines in `lsquant/core/chebyshev/chebyshev_moments.hpp` under extended `grep`.**
   All stage-2 C++ files (`contraction.hpp`, `contraction.cpp`, `curve_text.hpp`, `curve_text.cpp`, `api/bindings.cpp`) have zero matches for `delta|DOS|Kubo|basis|std::cout|getenv`. Note that the unmodified stage-1 header `lsquant/core/chebyshev/chebyshev_moments.hpp` contains the mathematical notation `(2 - delta_{m0})` in three Doxygen comment lines (matching `AGENTS.md` section 3), which would match a case-sensitive `grep` for `delta` across `lsquant/core` if comments are not excluded or if that stage-1 comment is not reworded. Per instructions, no stage-1 header was modified in this branch.
