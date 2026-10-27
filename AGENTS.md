# LSQUANT — project brief

Branch `cpp-unified` (cut from `main`). This file is the single source of truth for
the project: what exists, what we are building, how it is organised, how it is
validated, and the conventions. It will be refined as the work advances. The
`rust-port` branch and its AIWORK material are abandoned and must not be used as a
reference.

Everything in the repository is written in English: code, identifiers, comments,
commit messages, documentation. The maintainer may write in Spanish; the output is
English regardless.

---

## 1. Purpose

LSQUANT computes spectral and transport properties of large tight-binding models
with the Kernel Polynomial Method (KPM): density of states, spectral functions,
Kubo-Greenwood and Kubo-Bastin conductivities (charge and spin, longitudinal and
Hall), optical conductivity, time-dependent correlations, mean-square displacement,
spin dynamics and relaxation. Cost is linear in system size; nothing is diagonalised.

The new code unifies three existing codes into one C++ core with a Python/sympy
front end:

| Source | Location | What it contributes |
|--------|----------|---------------------|
| outer LSQUANT | `src/`, `include/`, `test/`, `examples/` | the refactored real-space KPM core, the `lsquant` dispatcher, goldens, six tutorials |
| inner LSQUANT | `lsquant-official/` (nested repo) | k-space Hamiltonians with real-space disorder via FFT (`SparseMatrixType_kQuant*`), non-orthogonal bases (`*_nonOrth*`, overlap preconditioners), the FFT Kubo contraction (`Kubo_solver_FFT`) |
| k_quant | `k_quant/` (nested repo, branch `python_implementation`) | the concatenator: a correlator written as an ordered word of operators and functions of `H`, evaluated exactly in the eigenbasis of each `H(k)` block |
| wannier2sparse | `wannier2sparse/` (nested repo) | the operator producer (`_hr.dat` -> supercell CSR, `.desc` sidecars, bundles) and the exact-diagonalisation oracle `tools/hr_exactdiag.py` |

Dropped on purpose: `Device/`, `Graphene_NNN*`, `Graphene_KaneMele*`, `*_ChrisVel*`,
`*_Bismuthene*`, `device-postProcess*`, `MeanSquareDisplacement_DeviceSim` (demo and
model-specific code), the 41 `inline_*` executables, the `chebmom1D/2D/TD` triplet of
formats, the C++ reconstruction kernels (`delta_chebF`, `greenR_chebF`,
`DgreenR_chebF`, `GammaAC_chebF`, Bessel code), and the hard-coded block sizes in the
inner non-orthogonal classes.

---

## 2. The unifying idea: every observable is a word

Every quantity in the three codes is the trace of an ordered word of static
operators `A_i` and spectral functions `f_i(H)`:

    C = Tr[ A_0 f_1(H) A_1 f_2(H) A_2 ... f_p(H) A_p ]

- static operators: identity, `VX`, `VY`, `VZ`, `SX`, `SY`, `SZ`, `VXSZ`, projectors
  such as `(1 + SX)/2`, the position commutator `[X, .]`;
- spectral functions: `delta(E - H)`, `G^+(E)`, `G^-(E)`, `dG/dE`, `exp(-iHt/hbar)`,
  the Fermi function, the identity.

The existing entry points in that language (`p` spectral slots, `rank` = slots whose
Chebyshev index is kept):

| Original driver / class | Word | p | rank |
|---|---|---|---|
| `inline_compute-kpm-spectralOp`, `Moments1D` | `Tr[A delta(E-H)]` | 1 | 1 |
| `inline_compute-kpm-nonEqOp`, `Moments2D`, Kubo-Greenwood | `Tr[VX delta(E-H) VX delta(E-H)]` | 2 | 2 |
| same moments, `kuboBastinFromChebmom` | `Tr[VX delta(E-H) VY dG^-(E)/dE]` + Fermi-sea integral | 2 | 2 |
| `MomentsAC`, optical conductivity | `Tr[VX delta(E-H) VX G(E+w)]` | 2 | 2 |
| `inline_compute-kpm-CorrOp`, `MomentsTD` | `Tr[A delta(E-H) U(t)^+ B U(t)]` | 3 | 1 |
| `inline_compute-kpm-MeanSquareDisplacement` | `Tr[[X,U(t)]^+ delta(E-H) [X,U(t)]]` | 3 | 1 |
| `inline_compute-kpm-TimeEvProjetedOp` | `Tr[P U(t)^+ S_a U(t) delta(E-H) P]` | 3 | 1 |
| `MomentsLocal`, `GetDeltaPhi` | `<i| delta(E-H) |i>`, `delta(E-H) |phi>` | 1 | 0 |
| `Kubo_solver_FFT` (inner) | Kubo-Greenwood / Bastin contracted by FFT | 2 | 2 |

A spectral slot is either **expanded** (its Chebyshev index is kept, producing one
dimension of a moment tensor) or **applied** (contracted on the fly with known
coefficients: `U(t)` with Bessel coefficients in `MomentsTD::Evolve`, `delta(E_0 - H)`
in `Deltaiter`). That choice is the only real difference between the 1D, 2D, TD and
Local classes; everything else they duplicate. In real space the word is therefore
evaluated as the moment tensor

    mu[m_1 ... m_p] = Tr[ A_0 T_{m_1}(H~) A_1 T_{m_2}(H~) ... T_{m_p}(H~) A_p ]

and contracted with the Chebyshev coefficient tables `c^{(i)}_m` of each `f_i`. In the
eigenbasis the same word is evaluated exactly (k_quant's path), which makes the
eigenbasis the oracle of the whole code.

---

## 3. Architecture

Four layers with one-way dependencies: parsers -> interface -> chebyshev -> algebra.
A lower layer never knows about files, strings or physics names.

```
lsquant/
  core/algebra/      LinearOperator (y = a*O*x + b*y), SparseMatrix (CSR, Eigen backend),
                     BlockDiagonalMatrix (k-space), DiagonalMatrix, OperatorSum, OperatorProduct,
                     ScaledOperator, AdjointOperator, BloclTransform (B, B^+ via FFTW batched plans),
                     MetricPolicy (S, S^-1, S^-1/2 appliers), StateGenerator (random phase, exact,
                     local, user), linalg (axpy, dot, batched dot = GEMM), SpectralBounds
                     (Lanczos, Gershgorin).
  core/chebyshev/    SpectralWindow (scale, shift), RescaledOperator, ChebyshevRecursion (single),
                     ChebyshevVectors and SlicedChebyshevVectors (batched, sectioned), ChebyshevMoments
                     (rank r tensor + metadata), MomentAccumulator (expand -> index, apply -> contract),
                     CoefficientTable (c_m(E_j) received as arrays), kernels (Jackson, Lorentz),
                     MomentContraction (rank 1 GEMV, rank 2 GEMM, FFT backend).
  core/basis/        Basis = {dim, layout, metric, transforms}: RealSpaceBasis, BlochBasis (k-grid,
                     W orbitals, FFT plans, disorder as a real-space diagonal), EigenBasis (per-block
                     U(k) and eigenvalues; the exact evaluation path).
  core/physics/      Correlator (vector of slots: operator | spectral function), ExpansionPlan (which
                     slots expand/apply, M per slot), presets as data (SpectralExpansion,
                     CorrelationExpansion, TimeCorrelation, MeanSquareDisplacement,
                     ProjectedTimeEvolution), units and prefactors, velocity rules per basis (bare
                     dH/dk, Berry/covariant, Lowdin, one-sided non-orthogonal).
  io/                CSR text, OperatorDescriptor (.desc), Wannier90 _hr.dat, BLOCH_PHASES, w2s bundle,
                     ChebyshevMoments file (.chebmom, rank + dims header), results JSON. Readers return
                     plain structs; no numerics.
  parsers/           word DSL ("VX * delta(E-H) * VY * dG(E)"), run-config JSON, CLI arguments.
                     Pure functions, unit-tested against an AST; they produce physics::Correlator.
  api/               C++ facade (Session: load basis, register operators, run word -> ChebyshevMoments)
                     and the pybind11 module lsquant._core (zero-copy numpy views).
  python/lsquant/    sympy layer (SpectralFunction -> CoefficientTable, Fermi-sea integrals, kernels),
                     drivers, plotting, MCP tool surface.
  apps/lsquant       thin CLI over api/ (compute | reconstruct | inspect | run).
  tests/             unit (recursion vs dense, operator algebra), golden (examples), oracles.
```

### Python is the library, C++ is the engine

The physics (words, spectral functions, bases, presets, units, names such as `VX`
or `delta`) lives in Python. C++ is a small engine that knows five instructions on
vectors and nothing else. Python compiles every word into a `Plan`, a list of those
instructions, and the engine executes it. A new observable is a Python function,
never a new C++ driver.

```cpp
// The whole C++ instruction set. Nothing physical is named here.
struct ApplyOperator { int operatorIndex; };            // v <- O_k v, for every live vector
struct Expand        { int numMoments; };               // v -> {T_0 v, ..., T_{M-1} v}: opens one tensor index
struct ApplySeries   { std::vector<Scalar> coeffs; };   // v <- sum_m c_m T_m v (Clenshaw), no index opened
struct MirrorToBra   {};                                // from here on build the bra side from a fresh copy of |r>
struct Trace         { std::vector<Scalar> weights; };  // mu[...] += w[...] * <bra_{m..}|ket_{n..}> per trace vector
using Instruction = std::variant<ApplyOperator, Expand, ApplySeries, MirrorToBra, Trace>;
struct Plan { std::vector<Instruction> code; std::string word; };

class MomentEngine {
public:
  MomentEngine(const LinearOperator& rescaledH, std::vector<const LinearOperator*> operators,
               StateGenerator& states, int batchSize);
  ChebyshevMoments Run(const Plan&);                     // rank = number of Expand instructions
};
```

Examples of compilation (Python side):

```
Tr[A delta(E-H)]                  ket: ApplyOperator(A), Expand(M); Trace(w)            -> mu[m]
Tr[VX delta VY dG^-]              ket: Expand(M), ApplyOperator(VY); MirrorToBra;
                                  bra: ApplyOperator(VX), Expand(M); Trace(1)          -> mu[m,n]
Tr[P U(t)^+ S_a U(t) delta P]     ket: ApplyOperator(P), Expand(M); MirrorToBra;
                                  bra: ApplyOperator(P), ApplySeries(bessel(t)),
                                       ApplyOperator(S_a), ApplySeries(bessel(-t)); Trace(w) -> mu[m](t)
```

Rank is whatever the plan opens; the cost (`M^rank` dot products, `M^(rank-1)`
live vectors) is estimated in Python before `Run` is called. Reconstruction never
happens in C++: the tensor goes back to Python, where coefficient tables contract it.

### No code is specific to a basis

The engine, `Plan`, `ChebyshevMoments`, the contraction and the I/O contain no
notion of basis; the word "basis" does not appear under `core/chebyshev/` (this is
checked by `grep` in every gate). What a basis is lives only in adapters that build
`LinearOperator`s and `StateGenerator`s:

| Basis | Adapter output | What the engine sees |
|---|---|---|
| real space | CSR `N x N` | a `LinearOperator` |
| Bloch | `BlockDiagonalMatrix(H(k))`, `BlochTransform` `B`, `B^+`, and the composite `H_k + B^+ V B` | a `LinearOperator` |
| eigen | `DiagonalMatrix(eps)`, `U`, `U^+`, `U^+ O U` per block | a `LinearOperator` |
| non-orthogonal | `S^-1/2 H S^-1/2`, `S^-1 (aH + bS)` as composites | a `LinearOperator` |

The exact (eigenbasis) concatenator of k_quant is therefore not a second code path:
with `H = DiagonalMatrix(eps)` the same plan gives exact moments, and with
`f(H) = DiagonalMatrix(f(eps))` the plan reduces to `ApplyOperator` and `Trace`.
Bases are mixed by inserting transform operators (`B^+ A_k B`) into the word; the
parser only checks dimensions. The basis is a per-calculation choice, never a
program-wide one.

`ChebyshevMoments` carries rank, dims, the window `(a, b)`, system size, number of
trace vectors, the `(2 - delta_m0)` weighting convention, the symmetrisation flag,
and the word as text. Rank 1 is the old `Moments1D`, rank 2 the old `Moments2D`,
rank 1 with a time axis the old `MomentsTD`, rank 0 (vector valued) the old
`MomentsLocal`.

Kept from the originals, moved into these layers: the CSR operator and Eigen matvec,
the three-term recursion, `Vectors`/`Vectors_sliced` (batching and memory
sectioning), `ComputeMomTable` (GEMM), `qstates::generator`, the `.desc` reader,
Gershgorin/Lanczos bounds, the Bloch FFT transform with real-space diagonal disorder,
the three `S^-1` appliers (PCG + ISCHO, Chebyshev of `1/x`, exact per block), the
Chebyshev `S^-1/2` applier, the Jackson/Lorentz formulas, the Kubo-Bastin Fermi-sea
integral and the Sea/Surf decompositions (as coefficient-level recipes), the FFT
contraction of `Kubo_solver_FFT`.

### State at the stage-7 closure (2026-10-21)

The legacy tree is gone from the working tree: `src/`, `include/`, `test/`, `legacy/`,
the legacy CMake targets, `scripts/regenerate_goldens.sh` and the legacy shell scripts
under `examples/`. The tag `legacy-final` (commit `05e2f65`) is the provenance of every
stage-0 golden: `git worktree add /tmp/legacy legacy-final` and the `build_legacy` recipe
of `goldens/stage0/PROVENANCE.md` rebuild the six legacy drivers and reproduce G1, G2 and
G5 byte for byte (verified at the closure, `handoff/0007-words-cli-mcp/tests.md`).
`lsquant/compat/` stays: it holds the two legacy-exact kernels (the GSL Bessel series of
the time evolution and the text round trip of the moments) that byte identity with the
stage-0 captures needs, and nothing else may be added to it. The layout as built differs
from the plan above in two places, both deliberate: there is no `core/physics/` and no
`parsers/` in C++, because the words, the catalogue, the compiler to `Plan`s, the run
file, the CLI and the MCP server live in `python/lsquant/` (`words/`, `runfile.py`,
`cli.py`, `mcp_server.py`); and the bases are Python classes (`python/lsquant/bases.py`)
that build `LinearOperator`s and `StateGenerator`s for the engine. Every golden of stage 0
is reproduced from a run file under `runs/goldens/` through `python -m lsquant run`.

---

## 4. Python and sympy

- Binding: pybind11, fetched by CMake. `lsquant._core` exposes `Operator` (from
  scipy CSR arrays, from composition, or a Python-defined matrix-free operator
  through a trampoline for prototyping), `StateGenerator`, `Plan`,
  `MomentEngine.run(plan)` (numpy view of the tensor) and the metric and transform
  operators. No Python object is touched inside the recursion.
- sympy owns every function of `H`. A `SpectralFunction` is an expression in the
  rescaled variable `x` with parameters (`E`, `eta`, `t`, `mu`, `T`).
  `chebyshev_coefficients(f, M, params, kernel)` returns an `(M, N_E)` table: closed
  form when sympy finds it (`delta`: `T_m(E)/(pi sqrt(1-E^2))`; `G^+`:
  `-i e^{-i m arccos E}/sqrt(1-E^2)`; `exp(-ixt)`: `2(-i)^m J_m(t)`), else a DCT on
  Chebyshev nodes via `lambdify`. The kernel `g_m` and the KPM weight are folded into
  the table, so truncation and damping happen in exactly one place.
- Reconstruction is a contraction: rank 1 `A(E_j) = sum_m c_m(E_j) mu_m` (GEMV);
  rank 2 `C(E_j) = sum_mn c1_m(E_j) mu_mn c2_n(E_j)` (diagonal of `C1 mu C2^T`, GEMM).
  The Kubo-Bastin Fermi-sea integral is an integral of the coefficient tables,
  `int^{E_F} c1_m(E) c2_n(E) dE`, done once symbolically (this is the Sea/Surf split).
- Apply mode receives the same tables: time evolution is `apply(U, coeffs for dt)`,
  local DOS at `E_0` is `apply(delta, coeffs at E_0)`.
- Word syntax from sympy noncommutative symbols:

```python
VX, VY, H = sp.symbols("VX VY H", commutative=False)
word  = VX * delta(E - H) * VY * dG(E - H)
corr  = lsq.Correlator.from_expr(word)      # factor is spectral iff free symbols are {H} + parameters
sigma = lsq.evaluate(corr, basis, M=512, trace="random_phase", R=10, energies=E_grid)
```

  The classifier is the parser layer; the MCP server is one tool
  `evaluate(word: str, model: str, numerics: dict)` over this API.

---

## 5. The concatenator in real space

Evaluate the word right to left on each trace vector:

1. `phi = A_p |r>`.
2. Rightmost spectral slot: **Expand** -> run the recursion and keep the `M` vectors
   `T_n(H~) phi`, batched in sections (`SlicedChebyshevVectors`); **Apply** ->
   contract on the fly (Clenshaw) with the supplied coefficients, keep one vector.
3. Apply the next static operator to every kept vector; continue leftwards.
4. The bra side uses the adjoint word, `<r| A_0 T_m = (T_m A_0^+ |r>)^+`, so static
   operators provide `Adjoint` (a Hermitian CSR is its own adjoint).
5. Final step: dot-product table between left and right batches (rank 1 dot, rank 2
   one GEMM per section pair, rank 3 `M` GEMMs; the planner warns on `O(M^3)`).

The planner in Python chooses expand/apply from the cost model (expanding costs `M`
matvecs per kept vector to its right; applying costs `M_apply` matvecs once). Time
presets apply both `U(t)` slots and expand `delta`; Hall expands both slots;
optical conductivity expands both and applies the shifted `G(E+w)` per frequency.
Symmetrisation `mu_mn <- (mu_mn + conj(mu_nm))/2` (valid for Hermitian `A_i`) and the
`(2 - delta_m0)/N_r` weights are recorded in the tensor header, never applied twice.

---

## 6. The three bases

| Basis | H | operators | trace vectors | transforms |
|---|---|---|---|---|
| real space | `SparseMatrix` `N x N` from wannier2sparse or scripts | CSR (`VX` from hopping vectors, `SX` from basis tags) | random phase over sites; exact `e_i` | `B` (real -> k): FFT over the supercell with the atom gauge |
| Bloch (k-space) | `BlockDiagonalMatrix` of `H(k)`, `N_k` blocks `W x W`; disorder `V` diagonal in real space applied as the composite `B^+ V B` (inner `SparseMatrixType_kQuant`, recomposed in the algebra layer, no subclass override) | block diagonal (`dH/dk` from sympy or finite differences; spin blocks) | random phase per `(k, orbital)`; exact per orbital | `B^+` (k -> real); `U(k)` (Bloch -> eigen) |
| eigen | diagonal `eps_n(k)` from `eigh(H(k))` | `U^+ O(k) U` | exact | `U^+` |

Rules: the eigenbasis is the exact evaluation path and the oracle for the other two
and for the coefficient tables (`sum_m c_m(E) T_m(eps)` must reproduce `f(eps; E)`);
a word may mix bases only through explicit transforms (`A_real = B^+ A_k B`), the
parser refuses anything else; one model file (`_hr.dat` + `.uc` + `.xyz`) feeds all
three bases.

---

## 7. Non-orthogonal bases

`H c = eps S c` enters only through the basis `MetricPolicy`, with three appliers:
`S`, `S^-1`, `S^-1/2`. Two strategies, both composites, so the recursion never
changes:

- **Lowdin (default).** `H' = S^-1/2 H S^-1/2`, `O' = S^-1/2 O S^-1/2`, velocity
  `V'_a = S^-1/2 dH_a S^-1/2 - 1/2 {H', dS_a S^-1}` (the formula in
  `utilities/siesta2sparse/siesta2lsquant.py`). Hermitian, so the standard engine,
  kernels and tables apply unchanged. `S^-1/2` by the Chebyshev expansion of
  `x^-1/2` on `spec(S)` (`ChebyshevInverseSqrtOverlap`, generalised from per-block to
  any sparse `S` with Lanczos bounds on `S`).
- **Generalised (optional).** Recursion on `S^-1 H`, rescaled as
  `H~ = S^-1 (a H + b S) = a S^-1 H + b` (not `S^-1 (a H + b I)`, see section 10),
  with `S^-1` from PCG + incomplete Cholesky (`IschoInverseOverlap`), the Chebyshev
  expansion of `1/x` (`ChebyshevInverseOverlap`), or exact per block
  (`ExactBlockInverseOverlap`). One-sided velocity
  `V_a = S^-1 dH_a - S^-1 dS_a S^-1 H`.

| Basis | S | `S^-1`, `S^-1/2` | notes |
|---|---|---|---|
| real space | sparse CSR (Siesta via sisl, or w2s `S_W(R)`) | PCG + ISCHO, or Chebyshev `1/x`, `x^-1/2` with Lanczos bounds | real-space diagonal disorder leaves `S` untouched |
| Bloch | block-diagonal `S(k)` | exact per block (dense `W x W` LLT) or Chebyshev per block | cheapest; `B^+ V B` disorder keeps `S` block diagonal |
| eigen | solved away: `eigh(H(k), S(k))` | not needed | exact oracle for both strategies |

The policy is chosen once per basis (`Orthogonal | Lowdin(applier) |
Generalised(solver)`); words and presets do not change.

---

## 8. Verification discipline: golden masters and one figure per stage

The project has suffered before from tests piling on tests and nothing moving.
The rule here is the opposite, and it is strict:

- **Golden masters are generated once, from the legacy code, and committed**
  (`goldens/`, stage 0). They are the reference for every later stage. Nobody
  regenerates them from the new code; a golden changes only when Fable re-blesses it
  with a written reason.
- **Every closed change is checked against the goldens, byte for byte.** The gate is
  one script, `scripts/check_goldens.sh <stage>`, that runs the new code on the
  golden inputs and `cmp`s the produced files with the committed ones. A byte
  mismatch that is a pure floating-point reordering (max relative difference below
  1e-12) is reported as such and needs Fable's explicit acceptance; anything larger
  is a defect.
- **Every stage ends in one figure** that shows the thing working: the new result
  overlaid on the golden (or on the analytic curve) with a residual panel, produced
  by `figures/stage_N.py` from the committed goldens and the new output. The figure
  and its caption are the deliverable of the stage, not a test count.
- **No unit-test suites.** A stage may carry at most one developer self-check
  executable (`tests/dev/stage_N_check.cpp`, dense linear algebra on a small random
  matrix) used while implementing; it is not a gate and it is not extended. Codex
  writes no tests beyond the golden harness and the figure scripts unless Fable asks
  for one by name.

The goldens, by stage (section 9 says which stage closes which):

| Golden | Inputs (committed) | Reference output | Oracle |
|---|---|---|---|
| G1 chain1d DOS moments | `test/golden/chain1d/operators`, window `[-2, 2]`, exact trace, `M = 64` | `SpectralOp1chain1d...chebmom1D` | `mu_m = delta_m0` |
| G2 chain1d Kubo-Greenwood moments | same, `VX`,`VX` | `NonEqOpVX-VXchain1d...chebmom2D` | closed form B.15 (`python/lsquant/oracle/chain.py`, `python -m lsquant.oracle.chain`) |
| G3 chain1d DOS curve | G1 + Jackson, broadening 120 meV | `mean...JACKSON.dat`, tutorial 01 figure data | `1/(pi sqrt(4-E^2))` |
| G4 graphene `VX-VX` moments and Kubo-Bastin curve | `test/golden/operators`, seed 12345, `M = 64` | `.chebmom2D`, `KuboBastin...dat` | regression only |
| G5 Haldane Hall | `test/golden/haldane`, exact trace, `M = 128` | Bastin, Sea, Surf, Greenwood `.dat`; tutorial 03 figure data | `C = +1` |
| G6 chain ballistic MSD | `test/golden/chain1d_ballistic`, `M = 256` | `.chebmomTD`, `mean...EF*.dat`; tutorial 02 and 05 figure data | `MSD ~ t^2`, `sqrt(4-E^2)` |
| G7 spin precession | `test/models/chain1d_mag`, `M = 32` | `TimeProj...chebmomTD`; tutorial 04 figure data | `cos(wt)`, `-sin(wt)` |
| G8 tutorials 05, 06 | their generators with fixed seeds | their `.dat` curves | qualitative |
| G9 three bases | graphene and Haldane `_hr.dat` (wannier2sparse) | DOS and `sigma_xx` from `tools/hr_exactdiag.py` (k-mesh exact) | eigenbasis exact |
| G10 non-orthogonal | an orthogonal model and a fixed `T`; `lsquant-official/examples/Siesta_to_LSQuant` | legacy `inline_compute-kpm-spectralOp-nonOrth` DOS; `eigh(H(k), S(k))` DOS | spectrum invariance |

wannier2sparse (pinned, `scripts/fetch_wannier2sparse.sh`, sha 26cab4a) is the only
producer of real-space inputs and, through `tools/hr_exactdiag.py`, the eigenbasis
oracle for G9. Its examples 01 to 05 and the PdSe2 spin Hall case are the models for
the later stages.

---

## 9. Development plan: stages, each closed by a figure and a golden gate

Ordered from fastest and safest to hardest. A stage is closed when its gate passes
byte for byte and its figure is produced and reviewed; only then the next stage
starts. Handoff briefs follow the stage numbers (`handoff/000N-...`).

| Stage | Deliverable | Gate (`check_goldens.sh N`) | Figure (`figures/stage_N.py`) |
|---|---|---|---|
| 0 | Legacy build on the verification machine; generate and commit goldens G1 to G8 with `scripts/regenerate_goldens.sh` and the six tutorial scripts; legacy `ctest` green | legacy suite green; goldens committed with provenance | the six tutorial figures regenerated from the legacy build (proves the reference pipeline) |
| 1 | `core/algebra` + `core/chebyshev`: rank 1 and rank 2 moments, exact and random-phase trace, real space, orthogonal; CSR reader; `.chebmom` writer with the legacy layouts | G1, G2 moments byte-identical to the legacy files | `mu_m` (G1) and the `m,n <= 10` block of `mu_mn` (G2), new vs golden, with residual |
| 2 | pybind11 binding + sympy `SpectralFunction` -> `CoefficientTable` (delta, Jackson) + rank-1 contraction | G3 curve byte-identical | tutorial 01: DOS of the chain vs `1/(pi sqrt(4-E^2))` and vs golden, three sizes |
| 3 | rank-2 contraction (GEMM), `G^-`, `dG/dE`, Fermi-sea integral and Sea/Surf on the coefficient tables; graphene regression | G4 and G5 byte-identical | tutorial 03: `sigma_xy`, sea, surface, `sigma_xx` of Haldane vs golden; `C` in the caption |
| 4 | apply mode: time evolution (`exp(-iHt)` tables), local projection, MSD and projected words | G6, G7, G8 byte-identical | tutorials 02, 04, 05: localisation length, spin precession, diffusion curves vs golden |
| 5 | `BlochBasis` (block `H(k)` + `B^+ V B` disorder via FFT) and `EigenBasis` (exact concatenator) from `_hr.dat` | G9: real, Bloch and eigen within finite-size tolerance of the exact k-mesh curves; G1 to G8 still byte-identical | graphene DOS and `sigma_xx` in the three bases on one panel, residual below |
| 6 | `MetricPolicy`: Lowdin and generalised, the three overlap appliers | G10: transformed model equals the orthogonal DOS; Siesta example vs `eigh(H, S)`; G1 to G9 unchanged | DOS of the transformed model under both strategies vs the orthogonal golden; Siesta DOS vs exact |
| 7 | word parser, CLI, MCP tool; delete the legacy drivers | every golden reproduced from word syntax, byte-identical | the six tutorial figures regenerated from words, side by side with stage 0 |
| 8 | FFT contraction backend and performance | G4, G5 byte-identical through the FFT path | wall time vs system size, new vs legacy |

Why this order: stages 1 to 4 only re-express what the legacy code already does,
so each has a byte-level golden and cannot drift unnoticed; the basis work (5) and
the metric work (6) add physics with no legacy golden, so they come after the engine
is frozen and use the exact eigenbasis as their oracle; the parser (7) is cheap but
only meaningful once every path it drives exists.

The legacy tree (`src/`, `include/`, `test/`, `examples/`) stays untouched and
buildable until stage 7: it is the reference the goldens come from.

---

## 10. Known defects in the original code (do not port)

- `src/chebyshev_solver.cpp`, `MomentosDelta`: the loop does `y = 2 H r_{n-1}` without
  `- r_{n-2}` and reads `mu[2n-3] = <psi|r_n>`; not a Chebyshev recursion. The inner
  copy has the correct product form `mu_{2n-1} = 2<r_{n-1}|r_n> - mu_1`,
  `mu_{2n} = 2<r_n|r_n> - mu_0`.
- Non-orthogonal rescale: `Moments1D_nonOrth::SetInitVectors_nonOrthogonal` rescales
  `H` in place to `aH + bI` and adds `b*phi` again; `Iterate_nonOrthogonal` adds
  nothing. Correct: `S^-1 (aH + bS)` (the inner `Rescale_nonOrth`).
- `Iterate_nonOrthogonal_test`, `IterateAllSliced_test`: `error` accumulated inside
  `omp parallel for` without reduction.
- `Vectors_sliced_nonOrthogonal::MultiplySliced` divides by the scale factor to undo
  the in-place rescale inside the velocity; the velocity formula is buried in the
  recursion class.
- `linalg::orthogonalize`: PCG tolerance 1e-7 hard-coded, prints every iteration.
- Inner `SparseMatrixType_kQuant_nonOrth::set_S` hard-codes `duplicate_spin_blocks(S, 46)`;
  `ChrisVel` hard-codes `nk = 25`, `norbs = 46`.
- Outer `Kubo_solver_FFT::compute` runs a single random vector (`while (r < 1)`).
- `Moments2D::ApplyJacksonKernel` truncates rank-2 tensors differently from rank 1.
- `qstates::FillWithRandomPhase` seeds from `time(0)`.
- `recon_grid.hpp`: `alpha = 0.95` and `bounds_pad = 0.10` are environment knobs;
  they become explicit fields of `SpectralWindow` and of the coefficient-table grid.
- **Batched rank-2 moments are wrong whenever `BATCH_SIZE < M`** (found by the tester in
  stage 1): `Vectors::IterateAll` (`src/chebyshev_vectors.cpp`) leaves the recursion at
  `(T_{B-1}, T_B)` after filling a block of `B` vectors, and `CorrelationExpansionMoments`
  (`src/chebyshev_solver.cpp`) calls it again without restarting, so block `k` holds
  `T_{k(B-1)}` onwards instead of `T_{kB}`: the off-diagonal blocks of `mu_mn` are
  shifted by `k` columns (chain1d, `BATCH_SIZE=16`, `M=64`: 36 elements wrong, max
  error 4.0). The default `BATCH_SIZE = M` hides it, so every committed golden is
  unaffected. Evidence: `goldens/stage0/chain1d/legacy_defects/batch16/`. The batched
  reference for the new engine is the unbatched G2 itself; a correct engine gives
  byte-identical output for any batch size.
- `test/graphene_reproduce.sh` parses the spectral bound as "the largest decimal in
  the solver log", which includes the `[hh:mm:ss.mmm]` timestamp, so
  `graphene_regression` fails about one run in six (whenever the run ends at second
  50 or later). Legacy test; retire with the legacy tree.
- The legacy pipelines write moment files at `precision(digits10)` (15 significant digits) and
  reconstruct from the re-read file, so every legacy curve carries a text round trip that
  loses the last two bits of the moments (tutorial 03: 1e-13 in `sigma_xy`). The new library
  keeps full precision in memory; the tutorial apps round-trip on purpose to stay byte-identical
  with the stage-0 captures, and say so in a comment.
- `test/golden/operators/graphene_golden.HAM.desc` carries bounds `[-3, 3]` inside the
  graphene spectrum (Gershgorin radius about 7.5); with the sidecar present the
  moments diverge. The regeneration script deletes it before running.
- **The legacy uniform-grid reconstruction maps the energy axis wrongly for off-centre
  windows** (found at the stage-4 closure, `handoff/0004-time-evolution/tests.md`):
  `src/reconstruction.cpp` writes `E = e/sf - shf` with `shf = -BandCenter/HalfWidth`,
  i.e. `x HalfWidth + BandCenter/HalfWidth` instead of `x HalfWidth + BandCenter`
  (Bastin, Sea, Surf, Greenwood uniform-grid paths; the FFT-grid path is right). Every
  stage-0 golden has a centred window (`BandCenter = 0`), so none is affected. The new
  library had the same class of defect in `python/lsquant/reconstruct.py` (`x/scale +
  shift`), fixed in handoff 0004 review 3 to `(x - shift)/scale`. Golden G4c (graphene,
  `BOUNDS -8.5 7.5`, `goldens/stage0/graphene/offcentre_enclosing/`) keeps the legacy
  moments and the curves blessed from the new code; the legacy curves are kept as
  evidence in `goldens/stage0/legacy_defects/offcentre_axis/`.
- The MSD velocity is rescaled like the Hamiltonian (`MomentsTD::RescaleCon2ChebyshevDomain`:
  `V <- a V + b I`), which is only right for `b = 0`; reproduced by the stage-4 apps for
  byte identity with the (centred) goldens; the correct form for an off-centre window is
  `V <- a V`.
- `~/lsquant/build` on TCN is an older legacy build (`c1d9318`, reconstruction grid at
  0.99) that does not reproduce G5; goldens are generated only from this tree's legacy
  sources built as `goldens/stage0/PROVENANCE.md` describes (`build_legacy`).
- **The legacy operator loader returns silently on a missing file** (found at the stage-7
  closure, `handoff/0007-words-cli-mcp/tests.md`): `BuildOPFromCSRFile`
  (`legacy-final:include/sparse_matrix.hpp:132`) does nothing when `operators/<label>.<op>.CSR`
  does not exist, so a word with an operator that was never written runs on an empty
  matrix. Tutorial 06 asked the legacy driver for the identity operator `1`, which
  `make_ey.py` and `make_rashba.py` never write; its 139 `TimeProj1-P*` curves (the
  normalising denominators), `ey_ratelaw.dat` and `captured/index.tsv` came from that
  empty operator (values about 1e-6, noisy) and the stage-0 figure was right only because
  the tutorial normalises by the same curve. Re-blessed from the new code (which applies
  the identity) at the stage-7 closure; the legacy files are kept as evidence in
  `goldens/stage0/legacy_defects/tutorial06_identity_operator/`. The new code refuses a
  missing operator file.
- **The CSR readers fold every operator to `1/2 (M + M^+)` silently** (`ReadCsrText`,
  `lsquant::SparseMatrix::HermitianFold`, legacy `ConvertFromCSR`; kept in the new code
  for byte identity). A non-Hermitian input, for example a velocity built from a
  non-Hermitian position matrix (0006b notes section 11), enters the engine as its
  Hermitian part without any warning. Clean-up item: report the anti-Hermitian norm and
  refuse above a tolerance.

---

## 11. Conventions

**Naming.** Identifiers are English, expressive, and derived from the original
names so that someone who knows the old code recognises the new one. Classes
PascalCase, methods PascalCase (as in the originals: `SystemSize()`,
`MomentVector()`), free functions and files snake_case, namespaces lower case
(`chebyshev`, `linalg`, `states`, `basis`, `physics`, `io`). No abbreviations that
are not in the originals.

| Original | New |
|---|---|
| `SparseMatrixType` | `SparseMatrix` |
| `HamiltonianOp`, `SparseHamiltonianOp` | `LinearOperator`, `SparseMatrixOperator` |
| `SparseMatrixType_kQuant` (block `H(k)` + `apply_B`, `apply_Bdagger`) | `BlockDiagonalMatrix` + `BlochTransform` |
| `chebyshev::Moments`, `Moments1D`, `Moments2D`, `MomentsTD`, `MomentsLocal`, `MomentsAC` | `ChebyshevMoments` (rank generic) |
| `chebyshev::Vectors`, `Vectors_sliced` | `ChebyshevVectors`, `SlicedChebyshevVectors` |
| `ScaleFactor()`, `ShiftFactor()`, `CUTOFF` | `SpectralWindow` (`Scale()`, `Shift()`, bounds) |
| `utility::SpectralBounds` | `SpectralBounds` |
| `qstates::generator` | `StateGenerator` (`RandomPhaseStates`, `LocalStates`, `UserStates`) |
| `ComputeMomTable` | `MomentTable` |
| `CorrelationExpansionMoments`, `SpectralMoments`, `TimeDependentCorrelations`, `MeanSquareDisplacement`, `TimeEvolvedProjectedOperator` | presets `CorrelationExpansion`, `SpectralExpansion`, `TimeCorrelation`, `MeanSquareDisplacement`, `ProjectedTimeEvolution` |
| `MetricPolicy`, `OrthogonalMetric`, `NonOrthogonalMetric` | kept |
| `IschoPreconditioner`, `ChebyshevSInvPreconditioner`, `ChebyshevSInvSqrtPreconditioner`, `PerfectPreconditioner` | `IschoInverseOverlap`, `ChebyshevInverseOverlap`, `ChebyshevInverseSqrtOverlap`, `ExactBlockInverseOverlap` |
| `JacksonKernel`, `ApplyJacksonKernel`, `ApplyLorentzKernel` | `JacksonKernel`, `LorentzKernel` (Python, folded into `CoefficientTable`) |
| `delta_chebF`, `greenR_chebF`, `DgreenR_chebF`, `GammaAC_chebF`, `calculate_cn_bessel` | `SpectralFunction` (Python/sympy) -> `CoefficientTable` |
| `Kubo_solver_FFT`, `Greenwood_FFTs`, `Bastin_FFTs` | `FftMomentContraction` |
| `OperatorDescriptor`, `.desc` | kept |
| k_quant `Trace.compute(...)`, `SpectralOperator`, `SpectralRepresentations` | `Correlator`, `SpectralFunction`, `ChebyshevRepresentation` |
| `.chebmom1D`, `.chebmom2D`, `.chebmomTD` | `.chebmom` with a rank header |

**Change discipline.** Goldens before refactors; the suite runs with `KPM_SEED=12345`
before every commit; exact-match goldens stay byte-identical unless the change is
meant to move them. No build caches, binaries or data blobs in the repository.
Commit and push only when asked.

**Documentation.** In-code comments describe the problem solved, the approach,
assumptions, units and conventions, and complexity; not a restatement of the lines.
Templates for comments, PR text, tutorials and figures live in `agent_template/`.

---

## 12. Roles and workflow

Three agents, three roles, never mixed:

| Role | Agent | Owns | Never does |
|---|---|---|---|
| Orchestrator and reviewer | Fable (Claude) | this file, the `handoff/NNNN-<slug>/brief.md` and `codex.md` instructions, code review, checking Codex's `tests.md` and relaying the findings to Agi, the verdict on each handoff | write project source, the harness or the figure scripts |
| Developer | Agi (on Argon 4) | all source under `lsquant/`, `python/`, `apps/`, CMake; the one `tests/dev/stage_N_check.cpp` self-check; `handoff/NNNN/report.md` | edit `scripts/check_goldens.sh`, `figures/` or `goldens/`; certify its own work; touch the legacy tree |
| Tester | Codex, or an Opus subagent driven by Fable when Codex cannot be launched | on TCN: the stage-0 golden generation from the legacy build, the golden harness `scripts/check_goldens.sh`, the figure scripts `figures/stage_N.py`, running every gate byte for byte, and when a gate fails, localising the defect (first differing file, element, step and thread count); `handoff/NNNN/tests.md` | edit source under `lsquant/`; write tests Fable did not ask for by name; relax a tolerance, a seed or a golden |

The cycle for every work item:

1. Fable writes `handoff/NNNN-<slug>/brief.md`: scope, exact names and interfaces,
   restrictions, numbered implementation steps, the stage gate and the stage figure.
   Nothing in a brief is optional.
2. Agi implements on `cpp-unified`, step by step, running each step's check before
   the next, and writes `report.md`: files added, build command and output, every
   step's check output verbatim, deviations from the brief, open questions. Agi does
   not declare the item done.
3. Codex, on TCN and following `handoff/NNNN/codex.md`, writes or updates the golden
   harness and the stage figure script from the brief alone, runs the gate
   (`scripts/check_goldens.sh N`, `KPM_SEED=12345`) and the figure script, and writes
   `tests.md`: for every golden, byte-identical or not; when not, the first differing
   file, line and element, the max relative difference, and a diagnosis of where in
   the implementation the difference is introduced (which step, which thread count,
   which accumulation); the figure path and caption; anything that smells
   (hard-coded values, tolerance drift, a golden touched). Codex never fixes source.
4. Fable checks `tests.md` against the diff and the brief, and returns to Agi a
   numbered list of what is wrong and what to change. Back to 2 until the list is
   empty; only then the item is done and may be committed.

A brief is the contract: if the implementation needs to deviate, the deviation goes
in `report.md` and Fable amends the brief; the code never silently diverges.

**How instructions travel.** Agi and Codex each read their instructions from the
repository (`handoff/NNNN-<slug>/brief.md` for Agi, `codex.md` for Codex); nothing
is passed by chat. Agi works in its own clone of `origin/cpp-unified` on the
maintainer's machine (`C:\Users\adame\LSQUANT\lsquant-agi`, branch
`handoff/NNNN-agi`), commits there and never pushes. Fable reviews that clone,
carries the branch to `origin` (through TCN), and Codex pulls it in
`~/lsquant-unified` on TCN, runs the gate and commits the harness, figures and
`tests.md` on `handoff/NNNN-codex`. When the gate passes and the review is clean,
Fable commits the closed handoff into `cpp-unified` and pushes it.

**Verification machine.** Every gate, golden generation and figure runs on TCN
(`fatnode-tacn`, 64 cores, 2 TB RAM, g++ 11.4, cmake 3.22, Python 3.10), reached
through the Ionos jump host: `ssh ionos` then `ssh jh-x` (or `ssh jose-tcn` from the
maintainer's machine). The checkout there is `~/lsquant`. Local machines are for
editing only; a result produced elsewhere is not a gate result.
**Shared machine rule.** Several sessions (the developer's builds, two or three testers,
the maintainer's own jobs) run the same programs on TCN at the same time. Nobody kills a
process by program name or pattern (`pkill -f stage6_nonorth`, `pkill -f check_goldens`):
every launcher records the PIDs it starts (`nohup ... & echo $! > <stage>.pid`) and stops
only those. Each session uses its own worktree and build directory. A session that stops
another session's run has broken the gate of that session and says so in its report.
(Introduced after one tester session killed another session's gate by name on 2026-10-21.)
