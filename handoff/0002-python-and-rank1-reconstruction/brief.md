# Handoff 0002 (Agi) — stage 2: Python is the library (pybind11, sympy, rank-1 reconstruction)

Roles: Agi implements, the tester (`codex.md`) runs the gate and the figure on TCN, Fable
reviews. Read `AGENTS.md` first (sections 3, 4, 8, 9, 11 are binding), then this brief
completely. Nothing is optional. Deviations go in `report.md`. Stage 1 is closed; your
starting point is `origin/cpp-unified` after the stage-1 merge.

## 0. Scope

Make Python the library and C++ the engine, as `AGENTS.md` section 3 says, for the
rank-1 case: a pybind11 module `lsquant._core` over the stage-1 engine, a Python package
`lsquant` where spectral functions are sympy expressions turned into coefficient tables,
and the reconstruction of the density of states as a contraction. Deliverable: the chain1d
DOS curve G3 (`meanSpectralOp1chain1dKPM_M64_stateexactJACKSON.dat`, Jackson 120 meV)
reproduced **byte for byte** from the Python API, and tutorial 01 reproduced from Python.

Out of scope, do not start, no hooks: rank-2 contraction, Green's functions, Fermi-sea
integrals, `ApplySeries`, time evolution, Bloch or eigen bases, metric, the word parser,
CLI, MCP.

## 1. Where things go

```
lsquant/
  core/chebyshev/contraction.hpp .cpp     ContractRank1(moments, table) -> values (sequential sum, section 4)
  io/curve_text.hpp .cpp                  WriteCurveText(filename, x, y): legacy ostream formatting (section 4)
  api/bindings.cpp                        pybind11 module lsquant._core (section 3)
  CMakeLists.txt                          + option LSQUANT_BUILD_PYTHON (ON), find_package(pybind11 CONFIG), pybind11_add_module
python/
  pyproject.toml                          package "lsquant" (setuptools, no build of C++ here: the module is found on PYTHONPATH from the build dir)
  lsquant/__init__.py
  lsquant/_core.pyi                       stubs for the bound names (hand-written, exact signatures)
  lsquant/spectral_functions.py           SpectralFunction (sympy), DeltaFunction, legacy_delta_table (section 4)
  lsquant/kernels.py                      jackson_damping (section 4)
  lsquant/moments.py                      run_plan helpers: dos_plan(M, num_states), kubo_plan(M) (the two stage-1 plans, as data)
  lsquant/reconstruct.py                  dos_curve(moments, broadening_meV, alpha=0.95) -> (E, rho) (section 4)
  lsquant/apps/stage2_dos.py              the gate program (section 6)
  tests/dev/stage_2_check.py              the single self-check (section 5)
```

Legacy tree untouched. `pybind11` and `sympy` come from `pip install --user pybind11 sympy`
on TCN (record the versions in the report); `find_package(pybind11 CONFIG)` uses
`pybind11.get_cmake_dir()` via `-Dpybind11_DIR=` passed on the cmake line (document the
exact line). No `FetchContent`, nothing vendored.

## 2. Restrictions

1. Everything of stage 1 (`brief.md` section 2) still holds for `core/` and `io/`.
2. The Python layer never reimplements the recursion, the matvec or the accumulation;
   it only builds operators, plans and tables and calls `_core`.
3. No physics name enters C++: `contraction.hpp` knows "moments" and "table", not
   "delta" or "DOS". `curve_text` knows two columns of doubles.
4. The coefficient tables that must reproduce the legacy are computed with Python's
   `math` module in the legacy expression order (section 4), not with numpy ufuncs
   (numpy's SIMD `cos`/`acos` can differ from libm in the last bit, which breaks byte
   identity). The sympy path is the general one and must agree with the legacy table to
   1e-14; both are kept.
5. The bindings copy nothing: `SparseMatrix` from scipy CSR arrays uses the `int64`
   index and `complex128` value buffers as given (`FromCsr` validates them);
   `ChebyshevMoments.data` is returned as a numpy view with the tensor's `dims`.
6. Names: `lsquant._core` exposes exactly the C++ names (`SparseMatrix`,
   `SparseMatrixOperator`, `IdentityOperator`, `RandomPhaseStates`, `LocalStates`,
   `UserStates`, `SpectralWindow`, `RescaledOperator`, `ChebyshevMoments`,
   `ApplyKpmWeights`, `ApplyOperator`, `Expand`, `ApplySeries`, `MirrorToBra`, `Trace`,
   `Plan`, `MomentEngine`, `ReadCsrText`, `ReadChebmom`, `WriteChebmom`,
   `WriteLegacyChebmom1D`, `WriteLegacyChebmom2D`, `ContractRank1`, `WriteCurveText`).
   Python-side functions are snake_case.
7. No `print` in the package; the gate program prints one line per output file.

## 3. The binding (api/bindings.cpp)

```cpp
PYBIND11_MODULE(_core, m) {
  // SparseMatrix: from_csr(dim, row_pointers: int64[dim+1], columns: int64[nnz], values: complex128[nnz]) -> SparseMatrix
  //               .dimension, .non_zeros, .hermitian_fold(), .is_hermitian(tol)
  // SparseMatrixOperator(SparseMatrix)  (keeps the matrix alive: py::keep_alive<1,2>)
  // IdentityOperator(dim)
  // RandomPhaseStates(dim, num_states, seed=12345), LocalStates(dim, sites: int64[]), UserStates(list of complex128[])
  // SpectralWindow(lower, upper) with .scale, .shift, .half_width, .center
  // RescaledOperator(SparseMatrix, SpectralWindow)   (keep_alive)
  // ApplyOperator(k), Expand(M), ApplySeries(coeffs: complex128[]), MirrorToBra(), Trace(weights: complex128[]), Plan(code: list, word: str)
  // MomentEngine(rescaled_h, operators: list[LinearOperator], states, batch_size).run(plan) -> ChebyshevMoments
  // ChebyshevMoments: .rank, .dims, .data (numpy complex128 view, shape dims), .window, .system_size, .number_of_states, .state_label, .word, .kpm_weighted, .symmetrised
  // ApplyKpmWeights(moments, num_states, symmetrise)
  // ReadCsrText(path), ReadChebmom(path), WriteChebmom(path, moments), WriteLegacyChebmom1D/2D(path, moments)
  // ContractRank1(moments, table: float64[M, N]) -> float64[N]
  // WriteCurveText(path, x: float64[N], y: float64[N])
}
```

`LinearOperator` is exposed as the common base so a Python list of operators can be
passed; a Python-defined operator (trampoline) is **not** in this stage.

## 4. Reconstruction, and what byte identity with G3 requires

Read `src/spectralFunctionFromChebmom.cpp`, `src/reconstruction.cpp`
(`reconstruct_density_grid`, `make_uniform_grid`, `write_pairs`),
`src/chebyshev_moments1D.cpp` (`ApplyJacksonKernel`, `MomentNumber`),
`src/chebyshev_moments.cpp` (`JacksonKernel`), `include/chebyshev_coefficients.hpp`
(`delta_chebF`) and `include/recon_grid.hpp` before coding. G3 is
`inline_spectralFunctionFromChebmom SpectralOp1chain1dKPM_M64_stateexact.chebmom1D 120`.
The legacy does, in this order:

1. `eta = 2.0*broad/1000/BandWidth; M_eff = ceil(pi/eta)`, capped at `M`; moments
   truncated to `M_eff`.
2. `mu[m] *= g(m, M_eff)` with
   `g = ((M-m+1)*cos(phi*m) + sin(phi*m)/tan(phi)) * phi/pi`, `phi = pi/(M+1)`.
   The kernel is applied to the **complex** moment, then `mu_real[m] = mu[m].real()`.
3. `num_div = 30*M_eff`; `x_i = -alpha + i*(2*alpha)/(num_div-1)`, `alpha = 0.95`.
4. `rho_i = sum_{m=0}^{M_eff-1} delta_chebF(x_i, m) * mu_real[m]`, sequential in `m`,
   with `delta_chebF(x, m) = cos(m*acos(x)) / sqrt(1 - x*x) / pi` (that exact order of
   operations: `fm/f0/M_PI`).
5. Written with `std::ofstream` at default precision (6 significant digits, `%g`
   style): `E_i = x_i/ScaleFactor + ShiftFactor`, `value_i = rho_i*(SystemSize*ScaleFactor)`,
   one `E value` line per point, `std::endl`.

The new pipeline reproduces this exactly, split as the architecture wants:

- `kernels.jackson_damping(M_eff)` returns `g[m]` computed with `math` in that order;
  `reconstruct.dos_curve` truncates and damps the moments in Python (complex
  multiply, then `.real`), as the legacy does.
- `spectral_functions.legacy_delta_table(x: list[float], M_eff) -> float64[M_eff, N]`
  evaluates `delta_chebF` with `math.cos`, `math.acos`, `math.sqrt`, `math.pi` in the
  legacy order, element by element (123k evaluations; speed is irrelevant here).
- `spectral_functions.DeltaFunction` is the sympy object: expression
  `cos(m*acos(x))/(pi*sqrt(1-x**2))` with symbols `x, m`; `table(x_grid, M)` lambdifies to
  numpy. The self-check compares it with `legacy_delta_table` to 1e-14.
- `_core.ContractRank1(moments_real, table)` does `values[i] = sum_m table[m,i]*mu[m]`
  as a sequential loop over `m` for each `i` (OpenMP over `i` allowed; the `m` sum is
  sequential), returning `float64[N]`. It takes a real `float64[M]` moment vector; the
  binding accepts the `.real` of the damped moments.
- `_core.WriteCurveText(path, E, values)` writes with `std::ofstream` default precision
  and `std::endl`, exactly as the legacy loop; the energy axis and the prefactor are
  computed in Python with the same expressions (`x/scale + shift`, `rho*(N*scale)`).

`reconstruct.dos_curve(moments: ChebyshevMoments, broadening_meV: float, alpha=0.95)`
returns `(E, values)` as numpy arrays and is the only public function of the stage.

## 5. The one self-check (tests/dev/stage_2_check.py)

A Python script, no framework: (a) builds the stage-1 `dos` and `kubo` plans from
`moments.py`, runs them on `test/golden/chain1d/operators` through `_core` and asserts
the tensors equal `ReadChebmom` of the committed legacy goldens to the last bit
(`numpy.array_equal`); (b) `DeltaFunction.table` vs `legacy_delta_table` max abs
difference printed (expected below 1e-14); (c) `dos_curve` on G1 at 120 meV vs the
golden G3 file read as text: prints the number of differing lines (expected 0). It is
Agi's tool, not a gate, and is not extended.

## 6. The gate program (python/lsquant/apps/stage2_dos.py)

```
python -m lsquant.apps.stage2_dos --moments <file.chebmom1D|.chebmom> --broadening-meV 120 --out <dir>
python -m lsquant.apps.stage2_dos --operators <dir> --label chain1d --window -2 2 --trace exact --states 2048 --moments 64 --broadening-meV 120 --out <dir>
```

The first form reconstructs from a moment file; the second computes the moments through
`MomentEngine` first (the stage-1 `dos` plan) and then reconstructs. Both write
`<out>/mean<label>JACKSON.dat` with the legacy name (`mean` + system label +
`JACKSON.dat`, where the label is the moment file name without extension, as
`Moments1D(momfilename)` sets it) and, for the second form, also the `.chebmom1D`. The
tester compares against `goldens/stage0/chain1d/meanSpectralOp1chain1dKPM_M64_stateexactJACKSON.dat`
and `SpectralOp1chain1dKPM_M64_stateexact.chebmom1D` with `cmp`. Exit non-zero on any
exception with the message on stderr.

## 7. Tutorial 01 through Python (for the stage figure)

`python/lsquant/apps/tutorial01.py` reproduces `examples/01_dos_1d_chain` with the new
API: builds the chains `N = 256, 1024, 4096` (reuse `examples/01_dos_1d_chain/make_chain.py`
by importing it, do not copy it), computes `M = 512` exact-trace moments, writes the
`.chebmom1D` files with the legacy names, and reconstructs at 20, 80, 24 and 10 meV as the
tutorial does. The tester compares its outputs with `goldens/stage0/tutorials/01/captured/`
(same file names) and draws the stage-2 figure from them.

## 8. Report

`handoff/0002-python-and-rank1-reconstruction/report.md`: files added; pip and cmake
lines with versions; build output on TCN (fresh build, `ls build/lsquant_core`, the
`_core*.so` path); `stage_2_check.py` output; the two gate-program runs with their
output lines; tutorial01 run; deviations with reasons; open questions. No "done", no
`cmp` against the goldens, no push.
