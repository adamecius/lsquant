# Handoff 0006 (Agi) — stage 6, run in parallel: non-orthogonal bases through `MetricPolicy`

Roles: Agi implements, the tester (`codex.md`) verifies the fixtures and runs the gate and
the figure on TCN, Fable reviews. Read `AGENTS.md` first (sections 3, 6, 7, 8, 10, 11 are
binding), then this brief completely. Nothing is optional. Deviations go in `report.md`.
This stage depends only on stage 1 (closed) and on the Python binding of stage 2 for the
model generator; it runs in parallel with stages 3 to 5. Branch `handoff/0006-agi` from
`origin/cpp-unified`.

## 0. Scope and the idea of the test

A non-orthogonal basis enters the engine only through operators. Given an orthogonal
model `(H, V_a)` and an invertible, local change of basis `T`, the same physics in the
new basis is `H' = T^+ H T`, `S' = T^+ T`, `V'_a = T^+ V_a T`. The generalised problem
`H' c = eps S' c` has exactly the spectrum of `H`, and traces are basis-invariant, so the
DOS moments, the Kubo moments and `sigma_xy` of `(H', S', V')` must equal those of
`(H, V)`: the goldens G1, G2 and G5 are the reference, with no new oracle. Two routes
must both reproduce them and agree with each other:

- **Löwdin**: `H_L = S^-1/2 H' S^-1/2`, `V_L = S^-1/2 V' S^-1/2`; Hermitian; the stage-1
  engine unchanged.
- **Generalised**: recursion on `H_G = S^-1 (a H' + b S') = a S^-1 H' + b`, with
  `S^-1 V'` as the velocity in the word and the adjoint bra of `theory.md` equations
  2a to 2c. (Note the rescale: `S^-1 (aH' + bS')`, never `S^-1 (aH' + bI)`; AGENTS.md
  section 10.)

Read `theory.md` in this folder before coding: it derives both routes, the Hermitian
velocity matrix, and lists the defects of the legacy non-orthogonal code that must not
be reproduced (one-sided velocity, symmetric estimator with `S^-1 V` on both sides,
`S^-1 (aH + bI)` rescale).

Deliverables: the `MetricPolicy` and its four appliers, the two composites, the model
generator in the Python library, the fixtures (chain1d and Haldane, rotated), the gate
program, and the figure data.

Out of scope: real Siesta inputs (`dHk`, `dSk`, the Löwdin velocity formula with the
`dS` terms; that is handoff 0006b after this one), Bloch per-block `S(k)`, any change to
`MomentEngine`, `Plan` or the stage-1 files.

## 1. Where things go

```
lsquant/core/algebra/
  metric_policy.hpp .cpp            MetricPolicy (interface), OrthogonalMetric, NonOrthogonalMetric
  overlap_appliers.hpp .cpp         ExactInverseOverlap, IschoInverseOverlap, ChebyshevInverseOverlap,
                                    ChebyshevInverseSqrtOverlap  (each a LinearOperator: y = S^-1 x or S^-1/2 x)
  composite_operators.hpp .cpp      + LowdinOperator(S^-1/2, O) = S^-1/2 O S^-1/2 ; + MetricInverseProduct(S^-1, O) = S^-1 O
lsquant/core/chebyshev/
  rescaled_operator.hpp .cpp        + GeneralisedRescaledOperator(H', S', S^-1 applier, window) = S^-1 (a H' + b S')
lsquant/api/bindings.cpp            + the above
lsquant/apps/stage6_nonorth.cpp     the gate program (section 6)
python/lsquant/models/nonorthogonal.py   rotate_basis(H, V_list, T) and make_local_transform(...) (section 5)
python/lsquant/apps/make_nonorth_fixtures.py   writes goldens/stage0/nonorth/<model>/<variant>/*.CSR + manifest (section 5)
tests/dev/stage_6_check.py          the single self-check (section 7)
```

## 2. Restrictions

Everything of stage 1 section 2 holds. In addition: no applier reads a file or prints;
tolerances and iteration caps are constructor parameters with no defaults hidden in
`.cpp`; nothing in `core/chebyshev` learns the word "metric" except the one
`GeneralisedRescaledOperator`, which is a `LinearOperator` like any other; the engine
does not change (if it has to, stop and report).

## 3. Interfaces

```cpp
namespace lsquant {

class MetricPolicy {
public:
  virtual ~MetricPolicy() = default;
  virtual bool IsOrthogonal() const = 0;
  virtual const LinearOperator& Overlap() const = 0;          // S
  virtual const LinearOperator& InverseOverlap() const = 0;   // S^-1 applier
  virtual const LinearOperator& InverseSqrtOverlap() const = 0;   // S^-1/2 applier
};
class OrthogonalMetric : public MetricPolicy { /* all three are IdentityOperator */ };
class NonOrthogonalMetric : public MetricPolicy {
public:
  NonOrthogonalMetric(const SparseMatrix& S, const LinearOperator& inverse, const LinearOperator& inverseSqrt);
};

// Appliers (all LinearOperator, all over a SparseMatrix S that is Hermitian positive definite)
class ExactInverseOverlap : public LinearOperator {         // dense LLT of S; Dimension <= a cap given in the ctor; S^-1 x
public:  ExactInverseOverlap(const SparseMatrix& S, Index maxDimension);  };
class IschoInverseOverlap : public LinearOperator {         // PCG with Eigen::IncompleteCholesky (the legacy ISCHO), S^-1 x
public:  IschoInverseOverlap(const SparseMatrix& S, double tolerance, int maxIterations);
         int LastIterations() const;  double LastError() const;  };
class ChebyshevInverseOverlap : public LinearOperator {     // Chebyshev expansion of 1/x on [a,b] ⊃ spec(S), S^-1 x
public:  ChebyshevInverseOverlap(const SparseMatrix& S, SpectralBounds boundsOfS, int order);  };
class ChebyshevInverseSqrtOverlap : public LinearOperator { // same with x^-1/2, S^-1/2 x
public:  ChebyshevInverseSqrtOverlap(const SparseMatrix& S, SpectralBounds boundsOfS, int order);  };
// The Chebyshev coefficients of 1/x and x^-1/2 on [a,b] are computed in C++ by the Chebyshev-Gauss
// quadrature the legacy chebyshev_sinv_preconditioner.hpp uses (N = order+1 nodes), with the
// adaptive truncation |c_m| < coeffTol |c_0| as a ctor parameter. This is the ONE place a function
// of an operator is expanded in C++, because it is an applier, not an observable; say so in the comment.

class LowdinOperator : public LinearOperator {               // S^-1/2 O S^-1/2
public:  LowdinOperator(const LinearOperator& inverseSqrt, const LinearOperator& O);  };
class MetricInverseProduct : public LinearOperator {         // S^-1 O
public:  MetricInverseProduct(const LinearOperator& inverse, const LinearOperator& O);  };

class GeneralisedRescaledOperator : public LinearOperator {  // S^-1 (a H' + b S') ; a = window.Scale(), b = window.Shift()
public:  GeneralisedRescaledOperator(const SparseMatrix& H, const SparseMatrix& S,
                                     const LinearOperator& inverse, const SpectralWindow& window);
         // holds its own CSR copy of (a H' + b S') built with the same Eigen expression style as SparseMatrix::Rescaled
};
}
```

`IsHermitian()`: true for `LowdinOperator` of a Hermitian `O`; false for
`MetricInverseProduct` and `GeneralisedRescaledOperator` (they are similar to Hermitian,
not Hermitian). The engine does not care; the comment must.

## 4. What the Löwdin route changes in the stage-1 plans: nothing

With `rescaledH = RescaledOperator(H_L)`? No: `S^-1/2 H' S^-1/2` is not a `SparseMatrix`.
Therefore `RescaledOperator` gains a second constructor over a `LinearOperator`
(`RescaledOperator(const LinearOperator& H, const SpectralWindow&)`, matrix-free
`y = a (H x) + b x`). Byte identity with G1/G2 is not required on this route (the
applier is iterative or truncated), so matrix-free is fine here; the stage-1 constructor
over a `SparseMatrix` stays as it is. The `dos` and `kubo` plans of stage 1 are reused
unchanged, with `VX` replaced by `LowdinOperator(S^-1/2, V'_X)`.

For the generalised route the Hamiltonian is `rescaledH =
GeneralisedRescaledOperator(H', S', S^-1, window)`, but the plan is **not** the
symmetric stage-1 one: `S^-1 V'` and `S^-1 H'` are not Hermitian matrices, so the bra
side must carry the adjoint explicitly (`theory.md`, equations 2a to 2c):

    ket: Expand(M), ApplyOperator(S^-1 V')
    MirrorToBra
    bra: ApplyOperator(S^-1), ApplyOperator(S^-1 V'), Expand(M), ApplyOperator(S)
    Trace(1)

(`ApplyOperator(S)` after `Expand` acts on the whole batch; the engine supports it.)
For the DOS no `S` enters the bra at all. A third route, `naive`, is the symmetric plan
with `S^-1 V'` on both sides; it is the legacy estimator and it is wrong (`theory.md`
section 4.1): the gate runs it to demonstrate the defect, never to pass.

`MetricPolicy::Overlap()` therefore matters in the plan: `S` itself is a
`SparseMatrixOperator` registered in the engine's operator list like any other.

## 5. The fixtures: rotated chain1d and rotated Haldane (Python, part of the library)

`python/lsquant/models/nonorthogonal.py`:

- `make_local_transform(H_csr, epsilon, kind, seed) -> T_csr`: `T = I + epsilon * A` with
  `A` built on the sparsity pattern of `H` (nearest neighbours), `kind = "real"`: `A`
  real symmetric with entries `±1` on the bonds (deterministic from `seed`);
  `kind = "complex"`: additionally a per-site phase `diag(exp(i theta_j))` multiplied
  on the right. `epsilon` chosen so that `cond(S)` is between 5 and 10 (the function
  reports `cond(S)` from dense `eigvalsh` for `N <= 4096`, and refuses `S` not
  positive definite).
- `rotate_basis(H_csr, operators: dict, T_csr) -> (H', S', {name: O'})` with
  `X' = T^+ X T` as sparse products, `S' = T^+ T`, Hermitian parts enforced exactly
  (`(X + X^+)/2`) before writing.
- `python/lsquant/apps/make_nonorth_fixtures.py` writes, for `chain1d` (from
  `goldens/stage0/chain1d/operators`, `HAM`, `VX`) and `haldane`
  (`goldens/stage0/haldane/operators`, `HAM`, `VX`, `VY`), the variants `real` and
  `complex`, as `goldens/stage0/nonorth/<model>/<variant>/<model>.{HAM,S,VX[,VY]}.CSR`
  in the legacy CSR text layout (`dim nnz` / values as `re im` pairs / columns / row
  pointers; read `lsquant/io/csr_text.cpp` and `examples/01_dos_1d_chain/make_chain.py`
  for the exact layout) written with scipy and Python's `repr` of doubles, so this stage
  does not wait for the stage-2 binding; plus `manifest.json`
  (epsilon, kind, seed, `cond(S)`, `lambda_min(S)`, `lambda_max(S)`, the spectral
  window of `H`, and the sha256 of each file). Deterministic: seed 12345.

The tester verifies independently (dense) that `S'` is positive definite with the
stated condition number and that `eigvalsh(H', S')` equals `eigvalsh(H)` to 1e-10.

## 6. The gate program (`apps/stage6_nonorth`)

```
stage6_nonorth --fixture goldens/stage0/nonorth/<model>/<variant> --label <model> --window <lo> <hi>
               --route lowdin|generalised|naive --applier exact|ischo|chebyshev --order <P> --tol <t>
               --trace exact --states <N> --moments <M> --word dos|kubo[xy] --out <dir>
```

Writes the legacy-layout moment files with the legacy names (as `stage1_moments`
does), plus the `.chebmom`, plus one line `applier=<name> order=<P> tol=<t>
max_iterations=<n> residual=<r>` on stdout. `kuboxy` is the Haldane `VX`,`VY` word of
G5 (`left = VY`? read `compute_noneq` and the G5 provenance: `inline_compute-kpm-nonEqOp
haldane VX VY 128 exact` means `op_right = VX`, `op_left = VY`, mapped to the stage-1
plan as ket `Expand, ApplyOperator(VX)`, bra `ApplyOperator(VY), Expand`; verify this
mapping against the committed G5 `.chebmom2D` on the orthogonal Haldane first, byte for
byte, and record it in the report).

The gate (tester) compares with `scripts/golden_diff.py` against G1, G2 and G5 at these
tolerances, both routes, and prints the max relative difference:

| route / applier | expected max relative difference vs golden |
|---|---|
| `lowdin` or `generalised`, `exact` (dense) | <= 1e-10 |
| `ischo` (tol 1e-12) | <= 1e-9 |
| `naive` (legacy estimator), any applier | reported, expected to FAIL on `kubo`/`kuboxy` and to pass on `dos`; this row is evidence, not a gate |
| `chebyshev` (order chosen by the ctor truncation at coeffTol 1e-12) | <= 1e-8 |

and the two routes against each other at the same levels. Then the Haldane Hall plateau:
the tester reconstructs `sigma_xy` with the stage-3 tools if they exist by then, else with
the legacy `inline_kuboBastinFromChebmom` on the produced `.chebmom2D` (it reads the
legacy layout), and asserts `C = +1` within 5%, both routes.

## 7. The one self-check (`tests/dev/stage_6_check.py`)

Rotated chain1d, `N = 2048`: `dos` moments through `lowdin/exact` and
`generalised/exact` vs the committed G1 (`delta_m0`), printed max abs error (expected
1e-10); the same with `ischo` and `chebyshev`, printing the applier residuals; and
`ExactInverseOverlap` vs `ChebyshevInverseOverlap` on random vectors. Not a gate, not
extended.

## 8. Report

`handoff/0006-nonorthogonal-metric/report.md`: files; the fixture manifests; the G5
operator-mapping check (byte-identical on the orthogonal Haldane); the self-check output;
one run of `stage6_nonorth` per route, applier and model with its stdout line;
deviations; open questions. No "done", no comparison against goldens beyond the
self-check, no push.
