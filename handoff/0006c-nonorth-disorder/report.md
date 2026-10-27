# Stage 6c developer report — nonorthogonal Anderson disorder

Implementation for review, based on `5c54542` (`origin/handoff/0006b-dev`).
Numerical validation is pending the tester on Aron. No build, Python execution,
self-check, material run, golden gate, or timing measurement was performed here,
as explicitly instructed. No numerical output in this report is a measurement.
`git diff --check` is the only source check run locally; it reports no whitespace
errors. This is not a stage closure or a golden certification.

The full 6c brief, AGENTS.md, stage-5 brief, stage-6 theory, stage-6b brief and
reviews 1–4 were read. The incoming staged brief/review documents are preserved
in the first commit. The work does not use the abandoned Rust port.

## Files and integration

- `python/lsquant/metric_forms.py`: `anderson_forms`,
  `disorder_velocity_forms`, `periodic_position_velocity_proof`. These return
  SymPy expressions, Hermiticity verdicts and zero-residual identities for the
  prescribed self-check to execute.
- `lsquant/core/algebra/composite_operators.{hpp,cpp}`: one general algebra
  primitive, `HermitianTripleProduct`, plus its binding and `_core.pyi` entry.
  No physics names, engine instructions, recursion, contraction or state-generator
  implementation changed. The added C++ uses C++17 facilities only; gcc 9.4/11.4
  compilation remains to be checked on the verification machines.
- `python/lsquant/models/disorder.py`: direct periodic CSR assembly from the
  SIESTA `H(R)`, `S(R)`, derivative and midpoint-dipole blocks; checked CSR
  Hermiticity; explicit `anderson_potential(S, w)`.
- `python/lsquant/disorder_bases.py`: `DisorderedBlochMetricBasis` and
  `RealSpaceMetricBasis`. Both implement the existing compiler/engine interface.
  All phase vectors cover the full finite system. Both metric routes are carried.
- `python/lsquant/apps/stage6c_disorder.py`: run-file executor and thin CLI;
  phase timings, moments, curves, saved disorder, Fermi-energy values and receipt.
- `python/lsquant/runfile.py`: accepts the boolean
  `model.params.nonorthogonal_disorder`; dispatches opted-in W=0 runs,
  nonzero-disorder SIESTA runs and real-space SIESTA runs to the new executor.
  Existing clean Bloch/eigen SIESTA files retain their 6b executor.
- `tests/dev/stage_6c_check.py`: the single requested self-check, (a)–(c) below.
- `runs/stage6c/bi_covariant_{bloch,real}_{lowdin,generalised}_W{0,0p1,0p3}.json`:
  twelve files, each computing DOS and spin Hall, 24 x 24, M=1024, 16 states,
  phase seed 12345, disorder seed 7, requested broadening 10 meV. Bloch uses
  exact block metrics; real space uses the stage-6 Chebyshev appliers (order
  1024, coefficient tolerance 1e-14, probe tolerance 1e-10). Each file has its
  own output directory. These are unmeasured gate inputs, not new goldens.

`bases.py`, `models/bloch.py`, `oracle/spin_hall.py`, the 6b app and self-check,
the golden harness, figures, committed golden data and compat are untouched.
The adapters live in a separate module to leave the concurrent 6b vectorisation
merge straightforward. The public imports are
`from lsquant.disorder_bases import DisorderedBlochMetricBasis, RealSpaceMetricBasis`.

## SymPy derivation and potential

`anderson_forms()` starts with the noncommuting matrices S and D (both
Hermitian) and the endpoint ansatz `alpha D S + (1-alpha) S D`. Assigning an
orbital's w to the ket endpoint gives SD; assigning it to the bra gives DS.
Hermitian conjugation exchanges them. Matching their real coefficients gives
`alpha=1/2`, so

    V = (D S + S D)/2,          V_ij = (w_i+w_j) S_ij/2.

This is the symmetric endpoint approximation to a local random potential,
not an assertion about the exact integral of every continuum random field.
The same function derives the physical, generalised, Löwdin and anticommutator
forms, and checks Hermiticity using the existing symbolic reducer. Numerically
this stage selects only generalised and Löwdin, as the 6c brief asks.

The constant-potential identity is decisive: `D=w I` gives `V=w S`, so the
physical Hamiltonian is shifted by w in every state. D alone would instead
produce `S^-1 D=w S^-1` in coefficient representation and would not be a constant
energy shift. Substituting `S=I` in the derived expression gives exactly D.
The functions return both residuals for assertion in the self-check.

The values come from the unchanged stage-5 `make_anderson_disorder`: NumPy
`default_rng(seed).uniform(-W/2,W/2,dimension)`, converted to complex128 with zero
imaginary part. Real ordering is cell-major, n1 fastest, then the interleaved
orbital/spinor index. W=0 avoids random draws, as stage 5 does. The overlap is
never changed. Runtime receipts save the values in `disorder.npy`, their SHA256,
the RNG convention and their layout.

## FFT route, metric and window

Let B denote the actual stage-5 `RealToBloch` implementation:

    B[(k,a),(R,b)] = delta_ab exp(-i k.(R+tau_a))/sqrt(N_cells).

It uses FFTW_FORWARD followed by the unnormalised atom phase and the constructor's
`1/sqrt(N_cells)`; its adjoint uses FFTW_BACKWARD. `SandwichOperator(B,D)` in
the existing API applies `D_k = B D B^+`. The 6c brief calls the real-directed
transform B, so its `B^+ D B` is this same map, with the names interchanged.
No FFT sign or normalisation implementation is changed.

The new adapter applies

    H_dis = H_k + (D_k S_k + S_k D_k)/2.

H, S and the exact cached inverse/inverse square root come from the 6b adapter.
No supercell sparse/dense matrix is assembled by this route. Clean arrays and
eigenvalues retained by the parent are used only for the clean operators,
metric and window bounds. `eigenbasis()` refuses the stage-6c call; disordered
Hamiltonian, velocity and current blocks also refuse `operator_blocks()`.
The executor refuses a per-k oracle request.

The applied Hamiltonian is `Y H_dis Y` for Löwdin (`Y=S^-1/2`) or `X H_dis` for
generalised (`X=S^-1`). Rescaling is always `a X H_dis+b I = X(a H_dis+b S)`.
For W=0 the Bloch Hamiltonian and clean physical operators retain the 6b cached
operation order. For exactly identity S, the potential and applied Hamiltonian
use the stage-5 sandwich directly; summing two half-potentials into H could
otherwise change bytes even though the mathematical matrices are identical.

The window uses the clean finite-mesh generalised spectrum, enlarged by

    ||Y V Y|| <= sqrt(s_max/s_min) max_i |w_i|.

Indeed `Y V Y=(Y D S^1/2+S^1/2 D Y)/2` and each product has this norm bound.
W/2 alone is not a safe enclosure for noncommuting D and S. The 10% window pad
then applies to the enlarged width. A resolved file stores the actual numeric
window. For an exact same-sample comparison, replay both layouts with that same
window; their independently computed clean extrema can differ in the last bit.
An explicit window is required to enclose the conservative norm enclosure.

## Direct real-space reference and velocity

`supercell_csr` places each M_ij(R) at `(cell,i),(cell+R,j)` and wraps only the
indices, summing duplicate bonds. It never Fourier-inverts a Bloch matrix.
Before wrapping, a derivative multiplies by
`i (R.lat+tau_j-tau_i)_a`, or by `i (R.lat)_a` for the bare rung. The same
midpoint connection A(R) certified by 6b is assembled as CSR. Thus boundary bonds
retain their physical displacement. A commutator with a wrapped finite-cell
position diagonal would instead introduce a spurious boundary cut.

The symbolic `periodic_position_velocity_proof()` starts with equation (3),
`i(H X R-R X H)`, puts `R={P,S}/2+A`, and uses

    P H = H P+i dH,    P S = S P+i dS,
    P X = X P-i X dS X.

Its normal-form residual against the 6b covariant velocity is returned for
assertion. The resulting numerical ladder, in both layouts, is

    velocity = dH - (dS X H + H X dS)/2 + i(H X A-A X H).

The connection term is present only in covariant mode; bare uses the cell
derivative pair in the atom-dressed representation. Covariant without dipoles
is refused. Since D commutes with the unwrapped orbital coordinates,

    dV = (D dS + dS D)/2.

`disorder_velocity_forms()` substitutes H+V and dH+dV into the symbolic ladder,
checks the residual and Hermiticity, and derives the correction to the clean
velocity used by the FFT route:

    correction = dV - (dS X V+V X dS)/2 + i(V X A-A X V).

This avoids incorrectly freezing the clean velocity for the generally nonlocal
matrix representation of the potential. The real route evaluates the full
ladder from the independently assembled CSR pieces. Spin is `Sigma=S q` and
spin current is `(velocity X Sigma+Sigma X velocity)/2`, in Pauli units, as in
`spin_forms()` and `spin_current_forms()` of 6b. Every numerical physical factor
is Hermitian before the compiler inserts the metric.

`HermitianTripleProduct(A,M,B,c)` applies
`c A M B+conj(c) B M A`, with dimension/Hermiticity validation of all factors.
It represents the anticommutators and commutators above without a Python matvec
or an intermediate matrix. Its local scratch and alias copy preserve the
`y=a O x+b y` contract; existing C++ matvecs retain their static schedules and
fixed summation order. No new reduction or parallel execution policy was added.

The real metric is the actual supercell S, with stage-6 exact inverse/square
root below `max_dim`, or stage-6 Chebyshev `1/x` and `x^-1/2` above it. Finite-mesh
clean overlap eigenvalues supply the interval for the polynomial (S remains
translation invariant). The clean Hamiltonian spectrum supplies only window
bounds, not reference moments. Both polynomial effective orders and deterministic
inverse/inverse-square-root probe residuals are written to the receipt. A probe
above `metric.tol` refuses the run before moments; it is a diagnostic, not a
certified operator-norm bound. The Bi polynomial order must be measured for
convergence rather than accepted because the probe passed.

## Words and stochastic trace

The existing catalogue/compiler, stage-6b prepare/compute functions and Python
reconstruction are reused:

    Tr[delta(E-H)]
    Tr[JXSZ delta(E-H) VY dG-(E)/dE]

Random phases are drawn once over the full orbital Bloch layout, exactly the
6b review-2 production convention. The real-space generator maps the same
stream by B^+. On the generalised route both original copies are prepared by Y,
with the dual S on the bra, before the compiler's asymmetric equation-(2a–c)
plan. Löwdin uses ordinary orthonormal vectors. No independent bra draw occurs.
Rank-2 weights and the 6b Hermitian trace projection are applied once. DOS is
reported per cell; spin Hall uses the unchanged `2*pi/(N_cells*A_cell)` Pauli
normalisation. No disordered per-k eigentrace is used, including at W=0 in these
files. If an exact trace is explicitly selected, it spans the FULL dimension.

Each phase prints its wall time immediately and persists it to
`run.resolved.json`, including failure timing. The receipt and `fermi_values.json`
contain the interpolated values at the requested Fermi energies. There is no
second production/reference/oracle pass hidden in a material run: the real
reference is a separate named run file with the same W, seed and phase stream.

## Self-check requested in section 1

The script fixes the 6b rotated Kane-Mele fixture (same local T, midpoint A,
interleaved spin), a 12 x 12 finite system, M=64, two states, phase seed 12345 and
batch size 17. It uses both metric routes and prints the symbolic residuals.
It deliberately does not add a unit-test suite or modify the 6b self-check.

(a) W=0: Bloch vs 6b moments to round-off (relative maximum 1e-12), and direct
real vs Bloch moments with maximum absolute difference 1e-10, DOS and spin Hall.

(b) W=1 eV, disorder seed 7, 12 x 12: real vs Bloch DOS and spin Hall moments
with relative maximum 1e-9. The relative denominator is max absolute reference
tensor element, avoiding singular ratios at symmetry zeros. Actual generator
vectors are checked through B before moments are compared.

(c) Hermiticity of CSR and FFT V, their Fourier agreement, the symbolic
constant-potential and orthogonal limits, identical stage-5 disorder bytes,
identical orthogonal sandwich actions and byte-identical DOS moment buffers
through the stage-5 Hamiltonian and the new orthogonal fast path. Both engines
receive the same window and actual orbital-Bloch phase stream; stage 5's default
`B*r_real` generator is a different sample from the 6b production convention.

No self-check output is available here. Expected successful output ends with
`stage_6c_check: (a)-(c) passed; Bi gate and thread-count checks remain with the tester`.
This line is a description of the script, not a claimed observed result.

## Deviations, modelling choices and open questions

1. The brief's literal list includes CSR velocity and spin current. H, S, dH,
   dS, A, potential and spin are CSR in the real route; velocity and current
   remain algebraic compositions of those CSR operators and the metric. S^-1
   is generally dense, so materialising the physical products as CSR would
   lose the reference route's linear storage. The formulas are equation (3)
   expanded into its periodic derivative form, derived symbolically above.
2. The two adapters live in `disorder_bases.py`, rather than modifying the
   concurrently edited `bases.py`. They reuse the existing basis interface.
3. The stage-5 RNG draws once per spinor coefficient. This preserves the brief's
   orthogonal stage-5 byte requirement and is recorded explicitly. If the intended
   material sweep instead requires the SAME draw on each up/down pair (a strictly
   nonmagnetic scalar disorder), that is a different ensemble. It should be
   specified and reviewed before blessing the Bi sweep; no hidden spin pairing
   or reseeding is introduced here. Explicit disorder arrays can already be
   supplied to both adapters for an agreed alternative ensemble.
4. The real adapter supports spin-scalar S, as the SIESTA route and prescribed
   rotated fixture do. An arbitrary supplied `physical_spin` array needs real
   translation blocks; the real adapter refuses it rather than guessing them.
5. Only a single deterministic probe checks polynomial metrics in a material
   run. Exact/Chebyshev convergence and the Bi plateau shape are unmeasured.
   No claim is made that M=1024 resolves the requested 10 meV on Bi's full
   spectral width. Effective KPM resolution must accompany the plateau plot.
6. No anticommutator approximation route is added to this side branch. Both
   exact metric representations requested in the 6c brief are implemented.

## What the tester must run and measure on Aron

Use NumPy 2.2, sisl 0.16 and the existing C++17 Aron build recipe. Build once
with gcc 9.4/11.4 as available and report diagnostics. Do not regenerate or
re-bless any old golden. Keep BLAS threading fixed at one when varying engine
OpenMP threads so that LAPACK/BLAS preprocessing is held constant:

```bash
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1
export KPM_SEED=12345
export PYTHONPATH=python:build:build/lsquant_core
OMP_NUM_THREADS=1 python tests/dev/stage_6c_check.py
OMP_NUM_THREADS=16 python tests/dev/stage_6c_check.py
```

Report all (a)–(c) residuals and the orthogonal byte verdict, without relaxing
the prescribed tolerances. Run the existing golden gate as required for the
shared algebra/binding change. Compare the W=0 Bi Bloch moments against 6b with
the same mesh, M, numeric window and actual orbital-phase stream; an exact block
trace is a different finite stochastic sample and is not a same-sample verdict.

Run every file under `runs/stage6c/`, for example:

```bash
OMP_NUM_THREADS=16 python -m lsquant run runs/stage6c/bi_covariant_bloch_lowdin_W0p1.json
OMP_NUM_THREADS=16 python -m lsquant run runs/stage6c/bi_covariant_real_lowdin_W0p1.json
```

For W=0, 0.1, 0.3 eV, report the DOS and spin Hall moments and reconstructed
curves for both layouts and both metric routes, disorder hashes, numeric
windows, metric effective orders and residuals, per-phase and total wall times,
and memory. For the strict layout comparison use the same numeric window from
one receipt in the paired run. Refine the real polynomial if necessary; report
the change rather than loosening the moment tolerance. The 12 x 12 exact check
is mandatory independently of any production polynomial diagnostics.

Plot the plateau against W with the S. Gimenez sweep shape as the qualitative
comparison requested by the brief, and show the paired-layout residuals and
wall times. Do not infer a clean per-k topological integer for a disordered
sample. Repeat representative W=0.1 outputs with OMP=1 and OMP=16 and compare
moment/curve bytes (timing receipts necessarily differ). No thread-count byte
verdict is claimed in this developer report.

## Commits and read-only Git fallback

The workspace `.git` is read-only: the first `git add` failed with
`Unable to create .../.git/index.lock: Read-only file system`. The original
branch/ref/index are unchanged. A copy of its metadata at
`/tmp/lsquant-stage6c-git` records the commits against this SAME workspace and
branch name `handoff/0006c-dev`; no alternate source branch or checkout was
created. Identity is local to that temporary metadata (`Codex <codex@openai.com>`).
The requested deliverable is `/tmp/lsquant-stage6c.bundle`, with prerequisite
`5c54542` and head `handoff/0006c-dev`. No push was attempted.

Coherent steps before this report:

- `bb52ff6`: incoming brief/reviews and SymPy Anderson/velocity derivation.
- `06bbd88`: Hermitian triple-product algebra and binding.
- `630d38a`: periodic-position proof, FFT adapter and direct CSR reference.
- `59c37ae`: run-file dispatch and physical-word executor.
- `6ec1acd`: required self-check, orthogonal byte-preserving operation order,
  twelve Bi sweep files.
- `e47db8f`: reject invalid disordered per-k blocks, finish the binding stub,
  enforce the absolute W=0 reference tolerance.

The final documentation commit records this report. Import the bundle into a
writable repository that already contains `5c54542`, then inspect/apply its
`handoff/0006c-dev` head there. The local source files match that bundled head;
the original workspace's Git status still compares against its unwritable base.
