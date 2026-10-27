# Handoff 0003 (Agi) — stage 3: rank-2 reconstruction (Kubo-Greenwood, Kubo-Bastin, Sea/Surf), graphene and Haldane

Roles: Agi implements, the tester (`codex.md`) runs the gate and the figure on TCN, Fable
reviews. Read `AGENTS.md` (sections 3, 4, 8, 9, 11) and this brief completely. Nothing
is optional; deviations go in `report.md`. Starting point: `origin/cpp-unified` after the
stage-2 merge; branch `handoff/0003-agi`.

## 0. Scope

The rank-2 half of the Python library: Green's-function spectral functions as sympy
objects, their coefficient tables, the rank-2 contraction, the Fermi-sea integrals, and
the four conductivity reconstructions of the legacy (`KuboGreenWood`, `KuboBastin`,
`KuboBastinSea`, `KuboBastinSurf`) plus the two FFT-node-grid variants. Deliverables:
G4 (graphene, random phase) and G5 (Haldane, exact) reproduced **byte for byte** from
Python, moments and curves, and tutorial 03 from Python.

Out of scope: `ApplySeries`, optical conductivity (`GammaAC`), the FFT *contraction*
backend (stage 8), any basis or metric work, parser, CLI.

**Item 0, carried from stage 2 (required, first commit of this handoff).**
`SparseMatrix.from_csr` in `lsquant/api/bindings.cpp` copies the three CSR buffers into
`std::vector`s and silently converts non-`int64` or non-contiguous arrays
(`py::array_t<..., c_style>` with forcecast). Make it what the stage-2 brief required:
accept only C-contiguous `int64` row pointers and columns and `complex128` values
(`py::array::forcecast` off; raise `TypeError` naming the offending argument otherwise),
and build the matrix from those buffers without an intermediate copy. `FromCsr` may keep
its by-value signature if a move from a mapped buffer is not possible; then document
in the report which single copy remains and why (Eigen's `Map` of CSR into an owned
`SparseMatrix` is itself one copy; that one is accepted, the `std::vector` staging is
not). `tests/dev/stage_3_check.py` adds one line: load `chain1d.HAM.CSR` with scipy,
pass its `indptr`, `indices`, `data` (converted to `int64`/`complex128`) to `from_csr`,
and assert the `dos` moments equal those from `ReadCsrText` after `hermitian_fold()`.

## 1. Where things go

```
lsquant/core/chebyshev/contraction.hpp .cpp     + ContractRank2Imag, ContractRank2Real, ContractRank2RealReal (section 4)
lsquant/compat/legacy_kernels.hpp .cpp          greenR_chebF, DgreenR_chebF tables with the legacy std::pow (section 4) -- gate-only module
lsquant/api/bindings.cpp                         + the above, + GershgorinBounds
python/lsquant/spectral_functions.py            + GreenRetarded, GreenRetardedDerivative (sympy), legacy_green_tables (via _core.compat)
python/lsquant/kernels.py                        + jackson_damping_2d(M0, M1)
python/lsquant/reconstruct.py                    + kubo_greenwood, kubo_bastin, kubo_bastin_sea, kubo_bastin_surf, kubo_fft_grid (section 4)
python/lsquant/moments.py                        + kubo_plan(left_index, right_index) (generalised), + window_from_gershgorin
python/lsquant/apps/stage3_kubo.py               the gate program (section 6)
python/lsquant/apps/tutorial03.py                tutorial 03 from Python (section 7)
tests/dev/stage_3_check.py                       the single self-check (section 5)
```

`lsquant/compat/` is a new directory with its own rule: it may contain legacy formulas
(it exists only so the gate can be byte-identical to files the legacy wrote), it is
used by nothing under `core/`, and it is deleted in stage 7 with the legacy drivers.
The sympy path is the real one and is checked against it to 1e-12.

## 2. Restrictions

All of stage 1 and stage 2 section 2. In addition: `core/` gains only the three
contraction kernels, which know "a rank-2 tensor and two tables" and nothing else; no
energy grid, no prefactor, no kernel cutoff in C++ outside `compat/`; the Python
reconstruction functions take a `ChebyshevMoments` and a broadening and return
`(E, values)` arrays, nothing else writes files but the apps.

## 3. Operator order in the legacy rank-2 moments (G4, G5)

`inline_compute-kpm-nonEqOp <label> <OP_R> <OP_L> M` maps to
`CorrelationExpansionMoments(OPL = OP_R, OPR = OP_L)`: the **first** CLI operator starts
the bra recursion (`left_m = T_m OP_R r`), the **second** is applied to the ket batch
(`right_n = OP_L T_n r`), so `mu_mn = Tr[ OP_R T_m OP_L T_n ]`. G5's `haldane VX VY`
therefore is `bra: ApplyOperator(VX), Expand; ket: Expand, ApplyOperator(VY)`.
`moments.kubo_plan(left_index, right_index)` takes the two operator indices and emits
that plan; verify on G5 (`VX-VY`, byte-identical) and record it in the report.

G4 uses random phase with seed 12345, `M = 64`, and the window from the Gershgorin
safeguard with pad 0.10 (`lower = -rho (1 + pad)`, `upper = rho (1 + pad)`, `rho` the
max absolute row sum of `H`): `moments.window_from_gershgorin(H)` through
`_core.GershgorinBounds`. The legacy header `BandWidth = 15.0366847652085` must come
out exactly.

## 4. Reconstruction, byte for byte

Read `src/reconstruction.cpp` entirely (`reconstruct_kubo`, `make_uniform_grid`,
`make_fft_node_grid`, `rearrange_crescent_order`, `accumulate`, `accumulate_trapezoid`,
`accumulate_midpoint`, `reconstruct_kubo_bastin_sea`, `reconstruct_kubo_bastin_surf`),
`include/observable.hpp`, `include/chebyshev_coefficients.hpp`, `src/chebyshev_moments2D.cpp`
(`ApplyJacksonKernel(b0, b1)`, `MomentNumber`) and the four `*FromChebmom.cpp` drivers
plus the two `*_FFTgrid.cpp`.

Kernel on a 2D tensor (`jackson_damping_2d`): `eta_i = 2 b_i/1000/BandWidth`,
`M_i = ceil(pi/eta_i)` capped at `dims[i]` (and set to `dims[i]` if `<= 0`), truncate the
tensor to `M_0 x M_1`, then `mu[m0,m1] *= (g(m0, M_0) * g(m1, M_1))`: the two damping
factors are multiplied first (double), then the complex moment is scaled. Same
`math`-module rule as stage 2.

Grids: uniform `num_div = mult * M` points on `[-alpha, alpha]`, `alpha = 0.95`,
`x_i = -alpha + i (2 alpha)/(num_div - 1)`, with `mult = 1` (Greenwood) or `30` (Bastin,
Sea, Surf); `M = max(M_0, M_1)`. FFT-node grid: `x_i = cos(pi (2i + 0.5)/M)`, `M` points,
rearranged to ascending by `rearrange_crescent_order` when the formula integrates.

Tables (`compat` for the gate, sympy for the library; both exposed):

- `delta(x, m) = cos(m acos x) / f0 / pi`, `f0 = sqrt(1 - x x)` (stage 2).
- `greenR(x, m) = -i * pow(x - i f0, m) / f0` with **`std::pow(std::complex<double>, double)`**
  of libstdc++: Python's complex power uses repeated multiplication for small integer
  exponents and does not round the same way, so this table is computed in
  `lsquant/compat/legacy_kernels.cpp` and exposed as `_core.compat.legacy_green_tables`.
- `DgreenR(x, m) = -i * pow(x - i f0, m) * (x + i m f0) / f0 / f0 / f0` (three sequential
  divisions), same module.
- sympy: `GreenRetarded = -I exp(-I m acos x) / sqrt(1 - x^2)` and its `E`-derivative,
  lambdified; the self-check shows they agree with the compat tables to 1e-12.

Contractions (C++, sequential `m0` outer, `m1` inner, exactly the legacy loops; OpenMP
over `i` allowed):

- `ContractRank2Imag(mu, A_real[M0,N], B_complex[M1,N]) -> real[N]`:
  `out_i = sum_m0 sum_m1 A[m0,i] * imag( B[m1,i] * mu[m0,m1] )` (Greenwood, Bastin).
- `ContractRank2Real(mu, A_complex[M0,N], B_complex[M1,N]) -> real[N]`:
  `out_i = sum sum real( (A[m0,i] * B[m1,i]) * mu[m0,m1] )` (Sea, with
  `A = Gr - Ga`, `B = DGr + DGa` built in Python as the legacy does: `Ga = conj(Gr)`,
  `DGr = real(DgreenR)`, `DGa = conj(DGr)` as a complex with zero imaginary part).
- `ContractRank2RealReal(mu, A_real, B_real) -> real[N]`:
  `out_i = sum sum (A[m0,i] * B[m1,i]) * real(mu[m0,m1])` (Surf).

Prefactors and integrals (Python, `math`, legacy order):

- Greenwood: `kernel_i = -1.0 * out_i * pref`, `pref = N * scale * scale`; values are
  `kernel_i` for `i < num_div - 1` (the last point is dropped); `E_i = x_i/scale + shift`.
- Bastin: `kernel_i = -2.0 * out_i * pref`; cumulative rectangle
  `acc += kernel_i (x_{i+1} - x_i)`, values `acc` for `i < num_div - 1`; same axis.
- Sea: `kernel_i = -out_i`, then `*= N scale scale / 2 / pi`; cumulative trapezoid
  `acc += 0.5 (k_i + k_{i+1}) (x_{i+1} - x_i)`.
- Surf: `kernel_i = out_i * (pi * N * scale * scale)`; midpoint `(k_{i+1} + k_i)/2`.
- FFT-node grid (Greenwood and Bastin): `pref = N / hw / hw`, `E_i = x_i hw + bc`,
  Bastin rearranged ascending before the integral.

Written by `_core.WriteCurveText` (default precision, `std::endl`), names
`KuboGreenWood_<label>JACKSON.dat`, `KuboBastin_<label>JACKSON.dat`,
`KuboBastinSea_<label>JACKSON.dat`, `KuboBastinSurf_<label>JACKSON.dat`,
`KuboBastin_<label>JACKSON_FFTgrid.dat`, `KuboGreenWood_<label>JACKSON_FFTgrid.dat`, where `<label>` is the moment file name without
extension.

## 5. The one self-check (`tests/dev/stage_3_check.py`)

(a) `kubo_plan` on Haldane `VX, VY` and `VX, VX` by exact trace vs the committed G5
`.chebmom2D` files: `array_equal`; (b) graphene random phase, seed 12345, Gershgorin
window: `array_equal` with G4; (c) sympy Green tables vs compat tables, max abs
difference (expected below 1e-12); (d) the four reconstructions of G5 and G4 vs the
golden `.dat` files read as text: number of differing lines (expected 0). Not a gate,
not extended.

## 6. The gate program (`python -m lsquant.apps.stage3_kubo`)

```
--moments <file.chebmom2D> --formula greenwood|bastin|sea|surf|greenwood_fft|bastin_fft --broadening-meV <b> --out <dir>
--operators <dir> --label <label> --left VX --right VY --window <lo> <hi> | --window gershgorin
  --trace exact|random --seed 12345 --states <N> --moments <M> --out <dir>
```

The first form reconstructs; the second computes the moments (writes the legacy
`NonEqOp<L>-<R><label>KPM_M<M>x<M>_state<label>.chebmom2D` and the `.chebmom`). One
stdout line per file written; non-zero exit with the message on stderr on error.

## 7. Tutorial 03 through Python (`python/lsquant/apps/tutorial03.py`)

Reproduces `examples/03_haldane_hall/lsqhall.py` with the new API: stages the Haldane
operators (import `make_haldane.py`), computes `VX-VY` and `VX-VX` moments (`M = 128`,
exact), reconstructs Bastin, Sea, Surf and Greenwood at 10 meV, writes the four `.dat`
with the legacy names into `<out>`, and prints the Chern number from the gap plateau as
`lsqhall.py` computes it (`C = (plateau/A) 2 pi`, `A = (dim/2) 3 sqrt(3)/2`,
`|E| < 0.5`). No plotting in the app; the tester draws.

## 8. Report

As in stages 1 and 2: files; build; `stage_3_check.py` output; one run of each gate
form and formula with the output lines; the tutorial-03 run with its `C`; the operator
order check of section 3; deviations; open questions. No "done", no `cmp` against
goldens beyond the self-check, no push.
