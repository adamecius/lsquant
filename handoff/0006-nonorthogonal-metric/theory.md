# Kubo and KPM in a non-orthogonal basis: what is correct, and what the legacy does

Fable, 2026-10-20. Derivation from first principles, then a line-by-line review of the
non-orthogonal implementations in `lsquant-official/` (Santiago's branch). The
literature cross-check is appended in section 6 when it arrives; the derivation below
does not depend on it.

## 1. Setting

Basis `{|phi_i>}`, overlap `S_ij = <phi_i|phi_j>`, Hamiltonian matrix `H_ij = <phi_i|H|phi_j>`,
any observable `O_ij = <phi_i|O|phi_j>`. All three matrices are Hermitian; `S` is positive
definite. A state `|psi> = sum_i c_i |phi_i>` is a coefficient vector `c`.

The matrix that represents an operator **acting on coefficient vectors** is
`O_c = S^-1 O` (because `O|psi> = sum_j (S^-1 O c)_j |phi_j>`). Consequences:

- Schrödinger: `H c = E S c`, i.e. `H_c c = E c` with `H_c = S^-1 H`.
- Products: `(AB)_c = A_c B_c = S^-1 A S^-1 B`.
- Functions: `f(H)_c = f(S^-1 H)`; Chebyshev polynomials `T_m(H~_c)` with the rescaled
  `H~_c = S^-1 (a H + b S) = a S^-1 H + b`. (`S^-1 (a H + b I)` is wrong: it shifts by
  `b S^-1`, not by `b`.)
- Traces are basis-independent: `Tr_op[F] = Tr[F_c]`. A stochastic trace uses plain
  random vectors `r` with `E[r r^+] = I`: `Tr[F_c] = E[r^+ F_c r]`. No `S`-weighting is
  needed **as long as the whole product is written in coefficient representation**.
- A change of basis from an orthonormal one, `T` invertible: `H' = T^+ H T`, `S' = T^+ T`,
  `O' = T^+ O T`; then `S'^-1 H' = T^-1 H T`, `S'^-1 O' = T^-1 O T`, and every trace of
  a product of coefficient-representation matrices is invariant.

## 2. The Kubo moments

Kubo-Greenwood (and the Bastin pair with `dG/dE`) needs
`mu_mn = Tr[ v T_m(H~) v T_n(H~) ]` as operators, i.e. in matrices

    mu_mn = Tr[ S^-1 V  T_m(S^-1 H~)  S^-1 V  T_n(S^-1 H~) ]                    (1)

with `V` the **Hermitian** matrix of the velocity operator in the basis. Two consistent
ways to evaluate (1):

**Löwdin route.** `H_L = S^-1/2 H S^-1/2`, `V_L = S^-1/2 V S^-1/2` (Hermitian), and (1)
becomes the orthogonal formula with `H_L`, `V_L`: `Tr[V_L T_m(H~_L) V_L T_n(H~_L)]`,
symmetric left/right vectors, plain dot products. The stage-1 engine unchanged.

**Generalised route** (recursion on `S^-1 H`, no square roots). The product in (1) is not
Hermitian-symmetric, so the left/right split must be done with the adjoint:

    r^+ A r,   A = S^-1 V T_m S^-1 V T_n,        T_k == T_k(S^-1 H~)
    (S^-1 r)^+ V T_m (S^-1 V T_n r)  =  [ T_m^+ V S^-1 r ]^+ [ S^-1 V T_n r ]
    T_m^+ = T_m(H~ S^-1) = S T_m(S^-1 H~) S^-1

so

    left_m  = S * T_m(S^-1 H~) * ( S^-1 V S^-1 r )                                (2a)
    right_n = S^-1 V * T_n(S^-1 H~) * r                                           (2b)
    mu_mn   = < left_m | right_n >  (plain dot)                                   (2c)

In words: the ket side is "recursion, then `S^-1 V`"; the bra side is "`S^-1`, then
`S^-1 V`, then recursion, then `S`" (or, equivalently, the `S`-weighted inner product
`<left|S|right>` with `left_m = T_m (S^-1 V S^-1 r)`). The two sides are **not** mirror
images, because `S^-1 V` and `S^-1 H~` are not Hermitian matrices.

Both routes give (1) exactly; they differ only in cost (`S^-1/2` applications in every
matvec vs. `S^-1` solves). In the rotated-basis test `(H', S', V')` both must reproduce
the orthogonal moments to the applier's precision.

In the engine's instruction set the two routes are plans, nothing else:

    Löwdin:       ket: Expand(M), ApplyOperator(V_L);  MirrorToBra;  bra: ApplyOperator(V_L), Expand(M);  Trace
    generalised:  ket: Expand(M), ApplyOperator(S^-1 V);  MirrorToBra;
                  bra: ApplyOperator(S^-1), ApplyOperator(S^-1 V), Expand(M), ApplyOperator(S);  Trace

(`ApplyOperator` after `Expand` acts on the whole batch, which the stage-1 engine already
does.) Rank 1 with an operator `A`: `Tr[S^-1 A T_m]` = ket `Expand`, bra
`ApplyOperator(S^-1), ApplyOperator(A)` (A Hermitian). The DOS needs no `S` at all on
the bra: `Tr[T_m(S^-1 H~)] = E[r^+ T_m r]`.

## 3. The velocity matrix

`V` must be the matrix of the operator `v = (i/hbar)[H, r]` in the basis, which is
Hermitian by construction. Two situations:

**Real space, position matrix known.** With `R_ij = <phi_i| r |phi_j>` (diagonal
approximation `R = diag(r_i)` when intra-atomic dipoles are neglected):

    V = (i/hbar) ( H S^-1 R - R S^-1 H )                                          (3)

Hermitian (`(H S^-1 R)^+ = R S^-1 H`). Note the `S^-1` inside: in a non-orthogonal basis
the velocity is **not** `i H_ij (r_j - r_i)`; that formula is the `S = I` limit of (3).
This is what wannier2sparse must emit when it gets an `S` (handoff 0006b).

**Bloch representation, derivatives known.** With `H(k)`, `S(k)`, `dH = dH/dk_a`,
`dS = dS/dk_a` (sisl `dHk`, `dSk`, gauge `r`), and neglecting the intra-cell position
term (Berry connection `A`), the Hermitian velocity matrix is

    V(k) = dH - 1/2 ( dS S^-1 H + H S^-1 dS )                                     (4)

Check on eigenvectors (`H c_n = E_n S c_n`): `c_n^+ V c_m = c_n^+ dH c_m - 1/2 (E_n + E_m) c_n^+ dS c_m`,
the standard non-orthogonal tight-binding matrix element. The one-sided
`dH - dS S^-1 H` gives `c_n^+ dH c_m - E_m c_n^+ dS c_m`: not Hermitian, and its
Hermitian part is (4), so the two differ by `1/2 (E_n - E_m) c_n^+ dS c_m`, an interband
term that enters `sigma_xy` and the optical response directly. With the intra-cell
position term, (4) gains `i (H S^-1 A - A^+ S^-1 H)`, i.e. (3) with `R -> A` in Bloch form.

Löwdin transform of (4): `V_L = S^-1/2 V S^-1/2 = S^-1/2 dH S^-1/2 - 1/2 ( D H_L + H_L D )`
with `D = S^-1/2 dS S^-1/2`. Note `D`, not `dS S^-1`.

## 4. Review of the legacy implementations (`lsquant-official/`)

### 4.1 Real-space one-sided: `Moments1D_nonOrth`, `Vectors_sliced_nonOrthogonal`, `Kubo_solver_FFT` + `inline_compute-Kubo-kpm-FFT-nonOrth`

- DOS: recursion on `S^-1 H` through CG (`linalg::orthogonalize`). Correct in structure.
  Defect: the rescale. `SetInitVectors_nonOrthogonal` rescales `H` in place to
  `aH + bI` and adds `b phi` again; `Iterate_nonOrthogonal` adds nothing. The correct
  operator is `S^-1 (aH + bS)`. For `b = 0` (centred band) the error vanishes, which is
  why chain-like tests looked fine. (AGENTS.md section 10.)
- Velocity (`Vectors_sliced_nonOrthogonal::SetInitVectors_2(OP, T0)`, `MultiplySliced`):
  `v_c = S^-1 dH - S^-1 dS S^-1 H / a`. One-sided (non-Hermitian `V`), see section 3.
  The `/a` undoes the in-place rescale of `H` by hand; with `b != 0` the `b S^-1`
  leftover is not removed, so the velocity carries a spurious `-S^-1 dS S^-1 (b/a)` term.
- Trace estimator: `left_m = T_m(S^-1 H~) v_c r`, `right_n = v_c T_n(S^-1 H~) r`, plain
  dot. This evaluates `Tr[ v_c^+ T_m(H~ S^-1) v_c T_n(S^-1 H~) ]`. With a Hermitian `V`
  (`v_c = S^-1 V`) it equals `Tr[ V T_m S^-2 V T_n ]`, whereas (1) is
  `Tr[ V T_m S^-1 V T_n S^-1 ]`. Under `T`: the estimator gives
  `Tr[ V T_m (T T^+)^-1 V T_n (T T^+) ]`, which depends on the basis; (1) does not.
  So the legacy Kubo moments are **not basis-covariant**; they coincide with the truth
  only when `S = I`. The bra side is missing one `S^-1` on `r` and one `S` after the
  recursion (equations 2a-2c).
- Driver (`inline_compute-Kubo-kpm-FFT_kQuant_nonOrth.cpp`): loads `dHk_x`, `dSk_x` for
  **both** directions (`OP_dHk_2`, `OP_dSk_2` read the `_x` files), so `sigma_xy` is
  actually `sigma_xx`; prints "tweaked for the bismuthene case only"; `set_S` hard-codes
  `duplicate_spin_blocks(S, 46)`.

### 4.2 k-space one-sided: `SparseMatrixType_kQuant_nonOrth::vel_i_nonOrth`, `Hk_clean_nonOrth`

Same formulas as 4.1 in Bloch form (`S^-1 dH - S^-1 dS S^-1 H / a`), same non-Hermitian
velocity, same non-covariant estimator, same `/a` patch. `Hk_clean_nonOrth` applies
`S^-1 H` with `H` already rescaled in place as `aH + bI` (`SetAndRescaleHamiltonian`),
i.e. the `S^-1 (aH + bI)` defect again.

### 4.3 Löwdin "ChrisVel": `SparseMatrixType_kQuant_nonOrth_ChrisVel::vel_i_nonOrth_Loewdin`

- `H_L` through the Chebyshev `S^-1/2` applier and `Rescale_nonOrth` (`aH + bS`): the
  rescale is the correct one here.
- Velocity: `V_L = S^-1/2 dH S^-1/2 + i ( H~_L A_L - A_L^+ H~_L )`, with `A` read from
  files `Ak_x.CSR` that no script in the repository produces (provenance unknown), and
  **no `dS` term**. If `A` is the full Bloch position operator including the
  `i d/dk` part this could be complete; if it is only the intra-cell dipole, the whole
  `-1/2 (D H_L + H_L D)` term of section 3 is missing. Cannot be verified from the
  repository. The driver loads `dHk_x` and `Ak_x` for both directions, so again
  `sigma_xy = sigma_xx`.
- The estimator is symmetric (Löwdin), so with a Hermitian `V_L` it would be
  covariant. Whether the `V_L` used is Hermitian depends on `A`.

### 4.4 Converters

- `utilities/siesta2sparse/siesta2lsquant.py` (Löwdin, dense, k-space):
  `V = S^-1/2 dH S^-1/2 - 1/2 H_L dS S^-1 - 1/2 dS S^-1 H_L`. Not the Löwdin transform
  of (4) (which has `D = S^-1/2 dS S^-1/2` in place of `dS S^-1`) and **not Hermitian**
  (`(H_L dS S^-1)^+ = S^-1 dS H_L`). Agrees with the correct form only when `S` and `dS`
  commute (`S = I`, or both diagonal).
- `siesta2linqt_kQuant_nonOrth.py::orth_operators_k`: `H_orth = S^-1 H` (mis-named: it is
  the generalised, non-Hermitian operator, not Löwdin) and `V = S^-1 dH - S^-1 dS S^-1 H`
  (one-sided). Consistent with 4.1/4.2 and inherits their defects.
- `siesta2lsquant_realspace.py`: Löwdin with `V = i [H_L, X_L]`, `X_L = S^-1/2 X S^-1/2`,
  `X` diagonal. This one is **correct** (it is (3) transformed), but it builds dense
  `S^-1/2` (`O(N^3)`), so it is a tool for small checks, not production.

## 5. Verdict

The legacy non-orthogonal code is **not correct in general**:

1. the Kubo trace estimator of the one-sided routes (4.1, 4.2) is not basis-covariant
   (missing `S^-1` on the bra start vector and `S` after the bra recursion);
2. the one-sided velocity `dH - dS S^-1 H` is not Hermitian; the correct Bloch form is
   (4), the correct real-space form is (3);
3. the Löwdin converter velocity is not the transform of (4) and is not Hermitian;
4. the DOS path has the `S^-1 (aH + bI)` rescale defect (harmless only for centred
   bands);
5. the drivers compute `sigma_xx` where they claim `sigma_xy`, and the Löwdin driver
   depends on `A` files of unknown provenance.

What survives: the appliers (`ISCHO` PCG, Chebyshev `1/x` and `x^-1/2`, exact per block)
are sound and reusable; the recursion on `S^-1 H` for the DOS is sound once the rescale
is fixed; `siesta2lsquant_realspace.py` is a correct small-scale reference.

Prediction for the stage-6 gate, used as the figure's punch line: on the rotated
Haldane fixture the correct routes (Löwdin, and generalised with (2a-2c)) reproduce G5
and `C = +1`; the legacy estimator (`--route naive`: symmetric plan with `S^-1 V'` on both
sides) does not, and its deviation grows with `cond(S)`.

## 6. Literature cross-check

(To be appended from the OpenAlex / arXiv survey: the generalised-eigenproblem KPM of
Weiße et al. and its `S`-weighted trace; the non-orthogonal tight-binding velocity
`c_n^+ [dH - 1/2 (E_n + E_m) dS] c_m`; Siesta/LCAO linear-scaling conductivity papers;
the ISCHO reference. Any disagreement with sections 1 to 3 will be stated here
explicitly.)

### 6.1 Result of the survey (`literature.md`, 19 DOIs; OpenAlex, arXiv, OSTI, sisl source)

No published work writes Kubo-Greenwood or Kubo-Bastin in Chebyshev/KPM form with
`S != I`. The closest, Röder, Silver, Drabold, Dong (PRB 55, 15382, 1997; 1996 conference
precursor), runs KPM on `S^-1 H` with PCG and the incomplete stabilised Cholesky
("ISCHO") for the **DOS only**, with a plain stochastic trace; that is the source the
legacy appliers quote. Weiße et al. RMP 2006 has no generalised-eigenproblem section;
Fan et al. Phys. Rep. 2021 is orthogonal only (`V_ij = (i/hbar) t_ij (R_i - R_j)`); the
ICN2/TCN output of 2021-2026 does not treat `S`. So sections 1 to 3 above are our own
algebra on top of published ingredients, and the stage-6 gate is the test of it.

### 6.2 Agreement with the sources

- Trace and functions of `H`: `Tr O = Tr[S^-1 O]`, `f(H) -> f(S^-1 H)`, plain random
  vectors: as in section 1 (Röder for the DOS; the rest is algebra).
- Velocity in `k`: Lee, Lee, Fukuda, Ozaki (PRB 98, 115115, 2018, Eq. 7), Jin et al.
  (PYATB, 2021, Eq. 24) and Esteve-Paredes, Palacios (SciPost Phys. Core 6, 002, 2023,
  Eq. 45) give `hbar v_mn = C_m^+ (dH - eps_m dS) C_n + i (eps_m - eps_n) C_m^+ A C_n`
  with `A(k) = sum_R e^{ikR} <0a| r |Rb>`. Its eigenvector-free Hermitian form is
  `V = dH - 1/2 (H S^-1 dS + dS S^-1 H) + i (H S^-1 A - A S^-1 H)`: equation (4) plus the
  connection term, as stated in section 3.
- Real space: `V = (i/hbar)(H S^-1 R - R S^-1 H)`, covariant under `T`: equation (3).
- The one-sided form `C_m^+ (dH - eps_m dS) C_n` is what sisl computes (`state.py`): it
  differs from the symmetric one by `1/2 (eps_m - eps_n) C_m^+ dS C_n`. **Refinement of
  section 3**: that difference is harmless for the on-shell Kubo-Greenwood (`eps_m =
  eps_n`), but it changes the Fermi-sea term and therefore `sigma_xy`. The legacy's
  `sigma_xx` with the one-sided velocity is thus less wrong than its `sigma_xy`; the
  estimator defect of section 4.1 affects both.
- Löwdin: build `V_L = S^-1/2 V S^-1/2` from `V`; `d/dk` of `H_L` is not the velocity
  for interband elements (literature item 17). Section 3 already says this.
- Green's functions with `S`: `G^± = ((E ± i eta) S - H)^-1 S`, i.e. `eta` multiplies
  `S`; this is automatically what the coefficient tables of `f(S^-1 H)` produce, so
  stage 3 needs no special case.

### 6.3 Two warnings for the fixtures (new, from the survey)

- Rebuilding a velocity from `(H', S', positions)` with the tight-binding convention
  `R ~ 1/2 (X S + S X)` does **not** give `T^+ V T` (the transformation induces
  intra-atomic `r` elements, Sandu 2005 Eqs. 10-13). The stage-6 fixture therefore uses
  `V' = T^+ V T` as the velocity matrix, never a velocity rebuilt from positions. A
  real Siesta model (handoff 0006b) needs the `A` term or accepts the diagonal
  approximation knowingly.
- `Tr[V' delta(E - H') V' delta(E - H')]` with no `S` at all is wrong even for the DOS
  (`T^+ H T` is congruent to `H`, not similar). The `naive` route of the brief keeps
  `S^-1` on `H` and on `V` and is wrong only through the estimator; both are recorded.

Not confirmed by the survey (and not needed by the gate): Graf-Vogl 1995 and Boykin 2001
were catalogued, not read; no source compares the variance of plain, Löwdin and
Cholesky-weighted random-vector estimators.
