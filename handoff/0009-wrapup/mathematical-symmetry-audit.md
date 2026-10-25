# Mathematical and spatial-symmetry audit of moments

Status: symbolic audit complete for the identities listed below; material symmetry
certification, optimizations and project closure remain pending.

This report complements [profiling-results.md](profiling-results.md). It answers
whether the current calculation exploits the mathematical structure, especially
how assumed spatial symmetries constrain Chebyshev moments. It is not an
exhaustive proof of the entire code or a claim that a particular point group
has been verified for the frozen Bi/PdSe2 models.

## Findings and priorities

The code already uses important identities correctly: metric similarity,
physical operator projection, exact-trace Hermitian symmetrization and the
frequency-folded Bastin reconstruction. It does **not** exploit all available
structure in the profiled clean projected workflow.

| Property | Current use in the inspected path | Opportunity / qualification |
|---|---|---|
| Chebyshev frequency folding | Implemented as optional `backend="folded"`; profiling used `tables` | Avoid the cubic reconstruction work on the same finite series/grid |
| Selected-subspace trace | `shift` narrows the spectrum but retains the full dimension | Remove inactive complement states for clean, projected words |
| Exact eigenbasis pair formula | Available in exact oracle / diagnostics | Compare with production stochastic KPM for small clean selected blocks |
| Exact-trace Hermiticity | Applied after both triangles have been computed | Cannot simply discard one sampled triangle without changing the estimator |
| DOS product identities | Generic rank-1 engine advances through every order | Specialized identity-DOS path can obtain two orders per recursion step |
| Spatial covariance and Brillouin-zone stars | No point-group mesh reduction found in the inspected material path | Add only after certifying symmetry of H, S, all operators and projection |
| Metric-aware transformations | Existing symbolic Lowdin/generalised proofs pass | Preserve them; Euclidean shortcuts in the generalised basis are unsafe |

Highest priority is reducing mathematically unnecessary work, before tuning
thread tiles. None of these proposals is a measured new end-to-end speedup.

## 1. Definitions and exact moment identities

Work first in an orthonormal representation, with Hermitian rescaled Hamiltonian
\(\widetilde H\) and Hermitian static operators A and B. Define raw moments

\[
\mu^{AB}_{mn}=\operatorname{Tr}
 [A T_m(\widetilde H) B T_n(\widetilde H)].
\]

Normalization by dimension, real Chebyshev weights and Jackson damping can be
restored afterward. Exact trace cyclicity and Hermitian conjugation give

\[
\mu^{AB}_{mn}=\mu^{BA}_{nm},\qquad
(\mu^{AB}_{mn})^*=\mu^{AB}_{nm}.
\]

Thus the real part is symmetric in m,n and the imaginary part is antisymmetric.
For A=B the exact moments are real and symmetric. A generic spin Hall word
has different operators, so its imaginary part must not be discarded.

These are exact-trace identities. For a finite stochastic estimate
\(\widehat\mu=R^{-1}\sum_r r^\dagger W r\), replacing the trace by one
quadratic form does not preserve cyclicity.

**Exact SymPy counterexample:** A=I+sigma_x, B=I+sigma_y,
H=sigma_z/2, r=(1,i)/sqrt(2) gives

\[
\widehat\mu^{AB}_{01}
 -(\widehat\mu^{AB}_{10})^*=-i.
\]

The current `ApplyKpmWeights(..., symmetrise=True)` averages the two sampled
entries. Computing only one triangle and mirroring it would produce a
different finite-R estimator. Exact sampling, a genuinely identical
bra/ket Gram construction, or a separately specified symmetry-adapted
estimator is needed to justify skipping work.

## 2. What an ordinary spatial symmetry does to moments

Let U_g be a unitary symmetry of the Hamiltonian and let the operators transform
in a real representation D:

\[
U_g\widetilde H U_g^\dagger=\widetilde H,\qquad
U_g O_a U_g^\dagger=\sum_\alpha D_{a\alpha}(g)O_\alpha.
\]

A polynomial respects similarity:
\(U_g T_m(\widetilde H)U_g^\dagger=T_m(\widetilde H)\).
Conjugating the complete word and using trace invariance gives

\[
\boxed{\mu^{ab}_{mn}
 =\sum_{\alpha\beta}D^{(A)}_{a\alpha}(g)
 D^{(B)}_{b\beta}(g)\mu^{\alpha\beta}_{mn}.}
\]

For p static slots, there is one representation factor per slot; the same
argument applies to the full ordered word. **An ordinary spatial symmetry
does not swap m and n and does not create an odd/even Chebyshev sign.**

If A and B each have a definite symmetry parity, U_g A U_g^dagger=eta_A A
and similarly for B, then

\[
\mu^{AB}_{mn}=\eta_A\eta_B\mu^{AB}_{mn}.
\]

When eta_A eta_B=-1, **every exact moment of that component is zero**,
not only its odd moments.

### Reflection example: charge versus spin Hall

For reflection x -> -x, polar velocity transforms as
(V_x,V_y,V_z) -> (-V_x,V_y,V_z). Spin is an axial vector:

\[
S_a\mapsto \det(R) R_{ab}S_b.
\]

Consequently S_z -> -S_z and
J_x^z={V_x,S_z}/2 -> +J_x^z.

| Operator / word component | Mirror-x sign | Exact consequence if the symmetry holds |
|---|---:|---|
| V_x | - | Polar vector |
| V_y | + | Polar vector |
| S_z | - | Axial vector |
| J_x^z | + | Product of two odd factors |
| Charge moment V_x,T_m,V_y,T_n | - | Zero for every m,n |
| Spin Hall moment J_x^z,T_m,V_y,T_n | + | Allowed; mirror does not force zero |

A polar-vector treatment of spin would wrongly eliminate this allowed spin
Hall component.

Under inversion, velocities are odd, spin is even and spin currents are odd.
Both the charge two-velocity word and the spin-current/velocity word are even.
**Inversion alone does not forbid spin Hall response or Hall response.**

### Rotation example

For a fourfold rotation about z, the charge tensor at each fixed m,n obeys

\[
\mu^{xx}_{mn}=\mu^{yy}_{mn},\qquad
\mu^{xy}_{mn}=-\mu^{yx}_{mn}.
\]

Combining the second relation with trace cyclicity gives
\(\mu^{xy}_{mn}=-\mu^{xy}_{nm}\). Combining also with Hermiticity makes this
exact component purely imaginary. A fourfold rotation alone allows an
antisymmetric Hall component; time reversal is a separate constraint.

SymPy's group average of a generic 3x3 charge tensor is

\[
\frac14\sum_{g\in C_4}R_g M R_g^T
 =\begin{pmatrix}
 (M_{xx}+M_{yy})/2&(M_{xy}-M_{yx})/2&0\\
 (M_{yx}-M_{xy})/2&(M_{xx}+M_{yy})/2&0\\
 0&0&M_{zz}
 \end{pmatrix}.
\]

For spin Hall the full object has three spatial indices:

\[
\mu^{i a,j}_{mn}
 =\det(R)\sum_{p b q}R_{ip}R_{ab}R_{jq}\mu^{p b,q}_{mn}.
\]

Here i is the spin-current flow direction, a is its spin polarization, and j
is the charge-velocity direction. For a z-preserving proper C4 rotation the
in-plane z-spin tensor has analogous spatial relations. Do not infer a
generic interchange of spin-current and velocity labels from that example.

## 3. Time reversal and chiral symmetry are different constraints

### Antiunitary time reversal

For an antiunitary symmetry Theta, trace conjugation adds a complex conjugate.
With Theta H Theta^-1=H and definite operator parities:

\[
\mu^{AB}_{mn}=\eta_A\eta_B(\mu^{AB}_{mn})^*.
\]

Velocity and spin are time-reversal odd, while their anticommutator spin
current is even. Therefore:

- Charge two-velocity moments are real in a time-reversal invariant exact trace.
- Spin-current/velocity moments are purely imaginary in that setting.
- Together with exact moment Hermiticity, the latter have zero diagonal m=n.

These are statements about raw Chebyshev moments and complete invariant
traces. They are not a blanket prescription to conjugate an entire
retarded/advanced response formula without transforming its spectral slots.

### Chiral / spectral-reflection symmetry

If a unitary Gamma instead satisfies
Gamma H_tilde Gamma^-1=-H_tilde, then

\[
T_m(-\widetilde H)=(-1)^m T_m(\widetilde H),\qquad
\mu^{AB}_{mn}=\eta_A\eta_B(-1)^{m+n}\mu^{AB}_{mn}.
\]

Now parity of m+n matters. For identity DOS, odd moments vanish. For a
higher-rank word the factor involves the sum of expanded polynomial orders.

This requires a true spectral-reflection symmetry. Spatial inversion usually
commutes with H and does not provide it. If the physical spectral symmetry
is centered at E_0, the rescaling center must also be E_0. Otherwise the
transformed dimensionless Hamiltonian carries the residual
2(E_0-center)/half_width instead of changing sign exactly.

## 4. Applying spatial symmetry to k-space, random traces and real space

At each k, a space-group operation relates H(k) to H(gk) by a unitary sewing
matrix in the orthonormal representation. With the convention

\[
U_g(k)O_a(k)U_g(k)^\dagger
 =\sum_\alpha D_{a\alpha}(g)O_\alpha(gk),
\]

the corresponding moment relation is
\(\mu^{ab}_{mn}(k)=\sum D^{(A)}D^{(B)}\mu^{\alpha\beta}_{mn}(gk)\).
The global relation follows only if the sampled mesh and its weights are
closed under the operation.

This permits irreducible-zone evaluation after the symmetry is certified.
For a scalar DOS, orbit multiplicities are sufficient. For tensor components
one must rotate and sum the components over the orbit, including stabilizers;
multiplying one representative component by the number of k points is
generally wrong. For example, SymPy obtains for a C4 xy orbit
[M_xy,-M_yx,M_xy,-M_yx], whose sum is 2(M_xy-M_yx), not 4 M_xy.

Before exploiting this for Bi or PdSe2:

1. Establish the actual operations and their spin/orbital sewing matrices,
   including atom gauge and any nonsymmorphic translation phases.
2. Verify H and S covariance, and covariance of every physical operator
   used in the word, including the chosen velocity construction.
3. Verify projection covariance P(gk)=U_g P(k) U_g^dagger in the chosen
   orthonormal convention. Energy-window selection should include complete
   symmetry-related multiplets; finite tolerances must not split them.
4. Verify mesh/weights, boundary conditions and hopping-cutoff closure.
   A real-space cutoff must retain complete symmetry-related hopping sets.
5. For stochastic runs, map the same states or explicitly define group
   averaging. A finite independent random sample is not invariant.

In a nonorthogonal basis, use the existing Lowdin transformation or a
metric-isometric sewing map; an ordinary Euclidean unitary test on physical
covariant matrices is not the correct substitute.

For finite stochastic states, define rho_R=R^-1 sum_r |r><r|. Exact trace
symmetry descends to that estimate only when rho_R is invariant under the
same group action, or through an explicitly symmetrized estimator.

The spatial SymPy counterexample has H=diag(1/3,2/3), U=diag(1,-1),
A=sigma_x and r=(1,1)/sqrt(2). The exact trace of A is zero, but r^dagger A r=1.
Averaging r and U r restores zero. Postprocessing a complete tensor with a
group projector can enforce the identities when the physical symmetry is
valid; it changes the original sample and is not byte-preserving.

A single generic disorder realization usually breaks crystal symmetries.
A symmetry of the disorder distribution constrains the **ensemble average**,
not every realization. Clean irreducible-zone reduction cannot simply be
reused after disorder mixes k points.

## 5. Chebyshev and reconstruction properties already available

SymPy verified the general trigonometric identity

\[
2T_m(x)T_n(x)=T_{m+n}(x)+T_{|m-n|}(x).
\]

For identity DOS with a Hermitian H_tilde and phi_m=T_m(H_tilde)r:

\[
\mu_{2m}(r)=2\langle\phi_m|\phi_m\rangle-\mu_0(r),\qquad
\mu_{2m+1}(r)=2\langle\phi_m|\phi_{m+1}\rangle-\mu_1(r).
\]

A specialized raw-moment path can reach roughly twice the order with the
same recursion depth. Use mapped projected seeds for selected DOS; in
generalised coordinates the inner product must respect the metric.
This is not a universal rule for arbitrary inserted operators.

For example, with H=sigma_z/2 and A=B=sigma_x, SymPy gives

\[
\operatorname{Tr}(A H B H)-\operatorname{Tr}(A B H^2)=-1.
\]

Thus one cannot collapse a generic rank-2 word by moving B through H.

For Bastin reconstruction, the implemented folded backend exploits

\[
\cos(m\theta)e^{-in\theta}
 =\tfrac12[e^{-i(n-m)\theta}+e^{-i(n+m)\theta}].
\]

Summing frequency diagonals and converting sine terms via second-kind
Chebyshev polynomials reduces the same finite-series reconstruction from
O(M^2 N_E) to O(M^2+M N_E). Since the current uniform grid has N_E=30M,
this changes cubic work to quadratic work. It also avoids M-by-N_E tables.

SymPy verified the derivative sign and the folded imaginary numerator for
arbitrary complex coefficients. Direct numerical corroboration of the
actual implementation used complex rectangular and square tensors up to
8x8; maximum absolute residual was 2.85e-13 and maximum residual normalized
by the reference peak was 1.08e-15. These small cases verify implementation
consistency, not high-order stability near the spectral endpoints.

The route already exists. The material profiling run did not use it.
Changing summation order is mathematically valid but not a promise of byte
identity with legacy tables.

## 6. Spectral compression: the main unused clean-case opportunity

The recorded projection receipts show:

| Material | Original orbitals per k | Retained per k |
|---|---:|---:|
| Bi | 72 | 12 |
| PdSe2 | 220 | 36-40, mean 39.5347 |

The current hard/shift adapter preserves the original dimension. It reduces
the spectral span, not the vector length. The distinct optional
`orthogonal` tool does reduce dimension, but requires a constant isolated
rank and a valid smooth gauge for its real-space representation.

In the orthonormal frame write P=Q Q^dagger and Q^dagger Q=I. For clean
spectral P, [P,H]=0. If every relevant static slot is projected, block
multiplication gives

\[
\operatorname{Tr}_N[PAP f(H_s) PBP g(H_s)]
 =\operatorname{Tr}_d[A_d f(H_d) B_d g(H_d)],
\]

where A_d=Q^dagger A Q and similarly for B and H. Complement energies do not
contribute. SymPy verified a symbolic retained 2x2 block inside a 3x3 matrix,
with independent operator entries and arbitrary spectral values, and the
same equality for a mapped arbitrary sample r_d=Q^dagger r. The argument
extends by block multiplication to arbitrary retained dimension and words.

Important restrictions:

- Keep the physical projected observable, including the spin current itself.
- Preserve normalization. If the old moments are Tr/N and a new generator
  returns Tr/d, multiply by d/N before retaining the old N prefactor, or
  consistently change both conventions. Do not silently mix them.
- To reproduce a finite stochastic run, map the actual vectors without
  renormalizing or generating a new reduced random stream.
- Use the retained projector/original eigenframe, not an energy test on
  the shifted Hamiltonian: discarded modes sit at the window edges and can
  share degenerate eigenspaces with retained modes. Resolve H and P
  consistently inside those eigenspaces before defining the reduced states.
- Unprojected identity-DOS intentionally counts complement modes and is not
  the same observable.
- PdSe2 has variable retained rank; the existing constant-rank orthogonal
  downfolder is not a drop-in solution. Per-k variable blocks or explicitly
  masked padding to 40 can be considered for clean k-space computation.
- Smooth real-space interpolation, finite radius truncation and disorder
  require separate validation. A generic disorder term does not commute with P.

The band counts suggest a roughly five-to-sixfold reduction in vector
length for a suitable clean implementation, **not** a measured time gain.

Projection also does not commute with forming composite operators:

\[
P\{V,S\}P/2-\{PVP,PSP\}/2
 = [PV(1-P)SP+PS(1-P)VP]/2.
\]

SymPy confirmed this missing term, and the existing code correctly projects
the original spin current rather than rebuilding it from projected factors.

## 7. Exact eigenbasis factorization and implementation choices

For a known clean eigenbasis, exact moments obey

\[
\mu^{AB}_{mn}=\sum_{ij} A_{ij}B_{ji}
 T_m(e_j)T_n(e_i).
\]

This factorizes into polynomial tables and an elementwise operator-pair
matrix; in matrix notation it is a pair of dense products with the appropriate
transpose convention. It avoids stochastic error for this exact clean route.
The project already uses the idea in diagnostics/oracles; it is not a new
physical approximation.

Benchmark it on retained blocks against the production estimator before
choosing a route. It is not universally cheaper, and eigenbasis availability
is specific to these clean preprocessed models. Nor does the formula imply a
useful low rank: summing many k blocks can make the moment tensor full rank.

For general real-space/disordered KPM, do not introduce diagonalization as a
replacement for the linear-scaling engine. The existing metric trace proofs
passed for the covered words; the anticommutator metric approximation remains
unequal for generic products.

## 8. Evidence, scope and next work

Executed on Aron using SymPy 1.14.0 and candidate-v5 Python with the existing
golden compiled engine:

- [sympy_audit.py](sympy_audit.py) and [sympy-results.json](sympy-results.json):
  48 explicit zero-residual records, including finite-degree corroborations,
  general identities, counterexamples and actual folded implementation checks.
  Existing metric trace proofs and all 12 spectral projection residuals passed.
- [sympy_spatial_audit.py](sympy_spatial_audit.py) and
  [sympy-spatial-results.json](sympy-spatial-results.json):
  C4 and reflection tensor projectors, polar/axial signs, time-reversal
  constraints, stochastic counterexample and chiral rescaling conditions.
- Remote reproducibility directory:
  `~/lsquant-wrapup-20261024/sympy-audit-20261025/`.
- Material dimension evidence:
  `~/lsquant-wrapup-20261024/profiling-20261025/{bi,pdse2}/run/run.json`.
- Inspected implementation:
  [projection](../../python/lsquant/projection/bases.py),
  [physical projection](../../python/lsquant/projection/spectral.py),
  [reconstruction](../../python/lsquant/reconstruct.py),
  [metric identities](../../python/lsquant/metric_forms.py),
  [moment engine](../../lsquant/core/chebyshev/moment_engine.cpp),
  [moment weights](../../lsquant/core/chebyshev/chebyshev_moments.cpp).

The finite polynomial checks are corroboration, not an infinite-order proof;
the general product/induction identities and block/trace arguments state the
corresponding general results. No material space group was inferred from
its name, and no symmetry tolerance or golden was changed.

Next implementation order: establish consistent trace/source comparisons;
compare folded reconstruction on identical moments; benchmark selected clean
blocks with mapped states; certify actual model spatial symmetries before
irreducible-zone reduction; then tune remaining engine costs. Preserve the
independent full-model/disorder requirements.

No production source optimization was made in this audit. No new performance
gain is claimed. Project closure still requires the user's explicit authorization.

