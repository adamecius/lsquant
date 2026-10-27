# Kubo conductivity and velocity operator in a non-orthogonal basis (S != I): literature synthesis

Date: 2026-10-20. Sources: OpenAlex (16 queries, plus forward citations of Fan et al. 2021 and Roeder et al. 1997), arXiv full texts (pdftotext), the OSTI scan of Roeder et al. 1996, and the sisl source.

Legend: **[SRC]** means stated in the cited paper (equation numbers are given). **[DER]** means my own algebra built on the sourced formulas: it is short and checkable, but no paper states it in this form. **[UNCONFIRMED]** means I could not verify it from a source.

---

## 0. Summary of what the literature does and does not contain

* **KPM with S is old and only covers the DOS.** Roeder, Silver, Kress, Landrum (1996 LANL conference paper) and Röder, Silver, Drabold, Dong (PRB 55, 15382, 1997) multiply H c = E S c by S^-1 and run ordinary KPM on H~ = S^-1 H. Each S^-1 application is a PCG solve of S x = H|m-1>, preconditioned with an **incomplete stabilized Cholesky (ISCHO)**. The trace is the plain stochastic trace <r|T_m(X)|r> with no S-weighting [SRC, Eqs. 1, 5-9 of the 1996 paper]. The legacy `ischo_preconditioner.hpp` quotes exactly this text, so **ISCHO = Roeder et al. 1996**.
* **The major KPM/LSQT reviews do not treat S.** Weiße et al. RMP 78, 275 (2006) has no section on generalized eigenproblems: a full-text search for "overlap" and "generalized eigen" found nothing. (The brief expected such a section; it does not exist.) Fan et al. Phys. Rep. 903 (2021) writes the velocity only for an orthogonal basis, <R_i|V|R_j> = (i/hbar) t_ij (R_i - R_j) (Eq. 22), "the expression used in all the LSQT methodologies discussed in this review".
* **No paper found writes the Kubo-Greenwood or Kubo-Bastin formula with Chebyshev/KPM and S != I.** I searched OpenAlex, arXiv and the web, and checked the 2021-2026 output of J. H. García, A. W. Cummings and S. Giménez de Castro. None of it handles S. Two recent uses of "non-orthogonal KPM" exist, both for the DOS / Fermi level only: Liu et al. arXiv:2609.37058 (2026), and the Roeder papers. The ML-TB code GPUTB (arXiv:2509.06525) avoids S altogether by learning an *orthogonal* Hamiltonian, and its stated reason is that LSQT for non-orthogonal bases is underdeveloped. **[UNCONFIRMED]** that no ICN2/TCN paper on Siesta-Hamiltonian KPM with S exists: I found none.
* **The velocity in a non-orthogonal basis is well established in k-space** (optics and Berry-curvature literature): Lee, Lee, Fukuda, Ozaki PRB 98, 115115 (2018); Wang et al. NJP 21, 093001 (2019); Jin, Zheng, He JPCM 33, 325503 (2021); Esteve-Paredes and Palacios SciPost Phys. Core 6, 002 (2023). For the Löwdin route and intra-atomic position terms: Sandu PRB 72, 125105 (2005). For position "hoppings" in orthogonal TB: Pedersen et al. PRB 63, 201101 (2001); Ibañez-Azpiroz, de Juan, Souza SciPost Phys. 12, 070 (2022).

The real-space, trace-based (KPM-ready) Kubo formulas in §1-§2 are therefore **[DER]**. They follow from the sourced eigenproblem, normalisation and velocity matrix elements by standard dual-basis algebra.

---

## 1. Change of basis from an orthogonal model: H' = T^+ H T, S' = T^+ T, V' = T^+ V T

### 1.1 Algebra [DER]
For a non-orthogonal basis {|a>} with S_ab = <a|b>, an operator O has covariant matrix O_ab = <a|O|b>. Restricted to the span, the operator is

    O^ = sum_ab |a> (S^-1 O S^-1)_ab <b|

so that:

    Tr O^            = Tr[ S^-1 O ]
    Tr (A^ B^)       = Tr[ S^-1 A  S^-1 B ]
    matrix of f(H^)  = S f(S^-1 H) = f(H S^-1) S        (since H c_n = e_n S c_n, c_n^+ S c_m = delta_nm)

This normalisation is the one stated by Wang et al. (Eq. 8: v_n^+ S v_m = delta_nm) and Lee et al. (Eq. 3) [SRC].

The natural objects are the **mixed (S^-1-weighted) matrices** h = S^-1 H and v = S^-1 V. For T invertible (square):

    S'^-1 H' = (T^+T)^-1 T^+ H T = T^-1 H T          (similar to H: same spectrum)
    S'^-1 V' = T^-1 V T

### 1.2 Correct Kubo-Greenwood in the primed basis [DER]

    sigma_xx(E) = (pi hbar e^2 / Omega) Tr[ v'_x delta(E - h') v'_x delta(E - h') ],
    h' = S'^-1 H',  v'_x = S'^-1 V'_x
    = (pi hbar e^2/Omega) Tr[ T^-1 V_x delta(E-H) V_x delta(E-H) T ]  = orthogonal result, exactly.

### 1.3 Correct Kubo-Bastin (Fermi sea + Fermi surface, Hall) in the primed basis [DER]
Take the Bastin formula in the orthogonal basis, as written in Fan et al. Eq. 20 and García et al. PRL 114, 116602 (2015) [SRC]:

    sigma_xy(mu,T) = (i hbar e^2/Omega) ∫ dE f(E) Tr[ V_x delta(E-H) V_y dG+/dE - V_x dG-/dE V_y delta(E-H) ]

Replace H -> h' and V -> v', with

    G±'(E) = (E - h' ± i eta)^-1 = ((E ± i eta) S' - H')^-1 S'

Then every product inside the trace is T^-1 (...) T, and the result equals the orthogonal one for every chemical potential and every eta.

Note [DER]: the finite-eta regulator multiplies S', giving (E+i eta)S' - H'. A Lorentz-kernel KPM on h' produces this automatically; the form (E S' - H' + i eta) does not, and it is not basis-covariant at finite eta.

### 1.4 Chebyshev moments
Rescale h' with the generalized-eigenvalue bounds (the spectrum of h' is real and equals spec(H)). Then

    mu_mn = Tr[ v'_x T_m(h~') v'_y T_n(h~') ],    T_m(h~')|r> by the recursion |m> = 2 h~'|m-1> - |m-2>

where each h' application is one S'^-1 solve (PCG+ISCHO, Roeder 1996 Eq. 8-9 [SRC]) or a Chebyshev/exact S^-1 applier.

### 1.5 Trace and random vectors
* **Generalised route (plain vectors).** Tr[X] = E[r^+ X r] with E[r r^+] = I, applied to the mixed product X = v' T_m(h') v' T_n(h'). **No S-weighting of the random vectors is needed.** This is the convention of Roeder 1996 Eq. 5 for the DOS [SRC]. The DOS Tr delta(E - h') counts states per orbital dimension without S. The diagonal [delta(E - h')]_aa is the **Mulliken PDOS**, sum_n c_an (S c_n)_a delta(E - e_n). That is exactly sisl's `PDOS` definition: psi* [S psi] [SRC sisl] / [DER identity].
* **Löwdin route (symmetric).**

      H_L = S'^-1/2 H' S'^-1/2,   V_L = S'^-1/2 V' S'^-1/2
      Tr[S'^-1 A S'^-1 B ...] = Tr[A_L B_L ...]

  This is fully Hermitian, so KPM can run unchanged with Hermitian moments. It is equivalent to drawing metric-weighted vectors xi = S'^-1/2 r in the generalised formulas.
* **Cholesky variant.** With S' = L L^+, the operators are L^-1 O L^-+ (also Hermitian), and the vectors are xi = L^-+ r. **[DER]**
* The estimators have the same mean and different variance. **[UNCONFIRMED]** which has lower variance: no source compares them.

### 1.6 What goes wrong with the naive formula (d)
**[DER]**, supported by [SRC] statements:
* Tr[V' delta(E - H') V' delta(E - H')] with the orthogonal formula is wrong already at the DOS level. H' = T^+ H T is congruent to H, not similar, so its eigenvalues differ from those of H unless T is unitary (Sylvester's law preserves only the inertia). Roeder 1996 states the correct problem is H|Psi> = lambda S|Psi>, solved through S^-1 H (Eq. 1, 7) [SRC].
* Mixing conventions is also wrong. Using h' = S'^-1 H' with V' (no S'^-1 on V'), or putting S'^-1 on only one of the two velocities, breaks the T^-1(...)T covariance. Each V' must carry its own S'^-1. The equivalent cyclic placements, e.g. Tr[V' S'^-1 ... ], are fine.
* Rebuilding the velocity of the transformed model from H', S' and the atomic positions (the "TB convention", §2.3) does **not** reproduce V' = T^+ V T. The position operator of the transformed basis, R' = T^+ X T, is not of the form ½(X'S' + S'X'). Sandu [SRC, Eqs. 10-13] shows explicitly that a basis transformation induces intra-atomic position matrix elements and longer-range Hamiltonian terms. So V' must be transported with T, or the full position matrix R' must be kept.

---

## 2. Natively non-orthogonal basis (Siesta / NAO / non-orthogonal Wannier)

### 2.1 k-space velocity [SRC]
Generalized eigenproblem: H(k) C_n = e_n S(k) C_n, with C_n^+ S C_m = delta (Lee Eq. 3; Wang Eq. 7-8).

Intraband (Hellmann-Feynman), from Lee Eq. 6, Wang Eq. 13 and Jin Eq. 24 with m = n:

    hbar v_nn = C_n^+ [ dH/dk - e_n dS/dk ] C_n

Full matrix, from Lee Eq. 7, Jin Eq. 24 and Esteve-Paredes & Palacios Eq. 45:

    hbar v_mn(k) = C_m^+ [ dH/dk - e_m dS/dk ] C_n + i (e_m - e_n) C_m^+ A^R(k) C_n,
    A^R_ab(k) = sum_R e^{ik.R} <0a| r |Rb>       (dipole / position matrix of the basis)

Lee Eq. 8 gives an equivalent form that needs only H(R) and <0a|r|Rb>, with no explicit S [SRC].

The basis Berry connection is A_nm = i C_n^+ S d_k C_m + C_n^+ A^(L) C_m, with A^(L) = sum_R e^{ikR}(<0n|r|Rm> - R<0n|Rm>) (Wang Eq. 9-10; Esteve-Paredes Eq. 20, 44).

Hermiticity relation: <0m|r|-R n> = (<0n|r|Rm>)* - R<0m|Rn> (Wang Eq. 11).

**sisl** (the reader used by this project) computes `derivative()` as <psi_i| dH/dk - e_i dS/dk |psi_j>. This is the one-sided form *without* the A^R term [SRC: sisl/physics/state.py docstring].

### 2.2 Symmetric, basis-level (KPM-ready) form [DER]
Use e_m C_m^+ S = C_m^+ H and S C_n e_n = H C_n, split -e_m dS = -½(e_m+e_n) dS - ½(e_m-e_n) dS, and invoke the Hermiticity relation of Wang Eq. 11. The velocity becomes a **k-dependent matrix** that no longer refers to eigenvectors:

    hbar V(k) = dH - ½ ( H S^-1 dS + dS S^-1 H ) + i ( H S^-1 Ã - Ã S^-1 H ),
    Ã(k) = A^R(k) + (i/2) dS = sum_R e^{ikR} <0a| r - R/2 |Rb>      (Hermitian)

V(k) is Hermitian, and C_m^+ V C_n reproduces Lee Eq. 7. This is the "dH/dk - ½{S^-1 dS, H}"-type expression asked about, written with the overlap inverse placed between the factors.

In Kubo, use v(k) = S^-1 V(k) (generalised route) or S^-1/2 V S^-1/2 (Löwdin route). Taking Ã = 0 is the "TB convention", which drops the r-hoppings of Ibañez-Azpiroz et al. [SRC: those are sizeable for shift current and minimal bases].

### 2.3 Real-space velocity [DER]
Let R_ab = <a|r|b> be the full position matrix (two-center integrals; Lee and Jin note these are computed like the overlap [SRC]). Projecting v = (i/hbar)[H, r] onto the span gives

    V = (i/hbar) ( H S^-1 R - R S^-1 H ),          v = S^-1 V = (i/hbar) [ h, x ],  h = S^-1 H, x = S^-1 R

Under a basis change, h -> T^-1 h T and x -> T^-1 x T, so v is covariant.

With the TB convention R_ab ≈ ½ (x_a + x_b) S_ab (diagonal atomic positions X, R = ½(XS+SX)):

    V = (i/2hbar) ( HX - XH + H S^-1 X S - S X S^-1 H )

This is exactly the real-space image of §2.2 with Ã = 0, via d_k H <-> i[H,X] and d_k S <-> i[S,X]. With S = I it reduces to Fan et al. Eq. 22.

This route needs one extra S^-1 application per velocity application.

Caveat [SRC, Jin et al. §2.4; Esteve-Paredes and Palacios §II; Lee §II.C]: in an incomplete basis, P[H,r]P != [PHP, PrP]. Jin et al. find explicit "correction terms to the naive Kubo formula" for the Berry curvature, which are negligible for full NAO bases and non-negligible for reduced bases. Non-local pseudopotentials also break p = i m [H, r] (Lee §II.C).

### 2.4 One-sided versus symmetric velocity [DER]
The one-sided matrix dH - H S^-1 dS (Lee/sisl without A^R) differs from the symmetric one by ½(e_m - e_n) C_m^+ dS C_n. This term vanishes on the energy shell, so Kubo-Greenwood longitudinal (Fermi surface, eta -> 0) is unaffected. It does **not** vanish in interband / Fermi-sea terms (Kubo-Bastin sigma_xy, optical), where the full Hermitian V with Ã is required.

### 2.5 Two consistent routes, side by side [DER]

| | Generalised (S^-1 H) | Löwdin (S^-1/2) |
|---|---|---|
| Hamiltonian in recursion | h = S^-1 H (non-Hermitian in the Euclidean metric, Hermitian in the S metric) | H_L = S^-1/2 H S^-1/2 |
| velocity | v = S^-1 V | V_L = S^-1/2 V S^-1/2 |
| trace / vectors | plain trace, E[r r^+] = I | plain trace, Hermitian moments |
| per step cost | 1 S^-1 solve (PCG+ISCHO, Roeder) | 2 S^-1/2 applications (Chebyshev x^-1/2) |
| k-space note | d_k acts on H, S | d_k H_L ≠ S^-1/2 (dH - ½{...}) S^-1/2. They differ by i[H_L, X_L]-type terms that cancel in intraband elements but not interband ones, so build V_L from V, not from d_k H_L. |

Sandu (Eq. 10-12) is the source for "transform every operator, including r, with (1+S)^-1/2" [SRC]. The warning in the last row is [DER].

---

## 3. Published implementations

* **Roeder, Silver, Kress, Landrum (1996); Röder, Silver, Drabold, Dong (PRB 1997).** KPM DOS on S^-1 H, PCG with ISCHO preconditioner, plain Gaussian random vectors. Extended Hückel PbTiO3/Pt (N = 868) and amorphous diamond. About 30x MVM overhead from PCG. No conductivity.
* **Liu, Zhong, Guo, Gong, Xiang, arXiv:2609.37058 (2026), HamGNN-NEGF.** "Non-orthogonal KPM" with an implicit S^-1 (linear solves) for the Fermi level only. Transport by RGF/Landauer. Cites Roeder 1997 and Silver 1996.
* **OpenMX (Lee, Lee, Fukuda, Ozaki 2018).** k-space Kubo-Greenwood optical / anomalous Hall conductivity with the non-orthogonal momentum matrix (Eq. 7), by diagonalisation.
* **ABACUS / PYATB (Wang et al. 2019; Jin, Zheng, He 2021; Jin et al. arXiv:2303.18004).** Berry connection, shift current and Berry curvature/AHC in NAO bases with dS/dk and the dipole matrix, by diagonalisation.
* **Esteve-Paredes and Palacios (2023).** Unified derivation; CRYSTAL/Gaussian-type local orbitals.
* **sisl.** Velocity as <psi_i| dH - e_i dS |psi_j> (no dipole term); Mulliken PDOS.
* **LSQT / KITE / GPUQT / LSQUANT orthogonal codes** (Fan et al. 2021 review; García, Covaci, Rappoport 2015). Orthogonal only: V from t_ij(R_i - R_j). GPUTB orthogonalises by construction.
* **Lorentzen et al. 2D Mater. 11, 035014 (2024)** (arXiv:2401.06556). Siesta S_ky and H_ky; Kubo-Greenwood via Green's functions with v = d_k H. **[UNCONFIRMED]** how S enters their trace and velocity: the main text does not say.

---

## 4. Bibliography

1. P.-O. Löwdin, J. Chem. Phys. 18, 365 (1950). doi:10.1063/1.1747632. Symmetric orthogonalisation S^-1/2.
2. H. Roeder, R. N. Silver, J. D. Kress, G. A. Landrum (+ D. Drabold, J. J. Dong on the scan), "The Kernel Polynomial Method for non-orthogonal electronic structure calculations", High Performance Computing 1996 (SCS), pp. 125-129; LA-UR-96-3302, OSTI 380322 (https://www.osti.gov/biblio/380322). No DOI. Source of ISCHO.
3. H. Röder, R. N. Silver, D. A. Drabold, J. J. Dong, PRB 55, 15382 (1997). doi:10.1103/PhysRevB.55.15382. Journal version: amorphous diamond. Full text not read; formulas taken from item 2.
4. R. N. Silver, H. Roeder, A. F. Voter, J. D. Kress, J. Comput. Phys. 124, 115 (1996). doi:10.1006/jcph.1996.0048.
5. A. Weiße, G. Wellein, A. Alvermann, H. Fehske, RMP 78, 275 (2006). doi:10.1103/RevModPhys.78.275; arXiv:cond-mat/0504627. No generalized-eigenproblem section (checked).
6. Z. Fan, J. H. García, A. W. Cummings, J. E. Barrios-Vargas, M. Panhans, A. Harju, F. Ortmann, S. Roche, Phys. Rep. 903, 1 (2021). doi:10.1016/j.physrep.2020.12.001; arXiv:1811.07387. Orthogonal velocity, Eq. 22; Bastin, Eq. 20.
7. J. H. García, L. Covaci, T. G. Rappoport, PRL 114, 116602 (2015). doi:10.1103/PhysRevLett.114.116602; arXiv:1410.8140. KPM Kubo-Bastin, orthogonal.
8. C.-C. Lee, Y.-T. Lee, M. Fukuda, T. Ozaki, PRB 98, 115115 (2018). doi:10.1103/PhysRevB.98.115115; arXiv:1805.08918. Eqs. 3, 6, 7, 8.
9. C. Wang, S. Zhao, X. Guo, X. Ren, B.-L. Gu, Y. Xu, W. Duan, New J. Phys. 21, 093001 (2019). doi:10.1088/1367-2630/ab3c9c; arXiv:1806.07018. Eqs. 7-13.
10. G. Jin, D. Zheng, L. He, J. Phys.: Condens. Matter 33, 325503 (2021). doi:10.1088/1361-648X/ac05e5; arXiv:2105.14662. Eqs. 6-10, 22-25.
11. J. J. Esteve-Paredes, J. J. Palacios, SciPost Phys. Core 6, 002 (2023). doi:10.21468/SciPostPhysCore.6.1.002; arXiv:2201.12290. Eqs. 19-21, 43-45.
12. T. Sandu, PRB 72, 125105 (2005). doi:10.1103/PhysRevB.72.125105; arXiv:cond-mat/0507204. Eqs. 8-13, Löwdin route.
13. T. G. Pedersen, K. Pedersen, T. B. Kristensen, PRB 63, 201101(R) (2001). doi:10.1103/PhysRevB.63.201101. Intra-atomic momentum term. Read via Sandu's Eq. 7, not directly.
14. J. Ibañez-Azpiroz, F. de Juan, I. Souza, SciPost Phys. 12, 070 (2022). doi:10.21468/SciPostPhys.12.2.070; arXiv:1910.06172. Position r-hoppings.
15. M. Graf, P. Vogl, PRB 51, 4940 (1995). doi:10.1103/PhysRevB.51.4940. Peierls coupling in TB. **Not read; content not verified for S != I.**
16. T. B. Boykin, R. C. Bowen, G. Klimeck, PRB 63, 245314 (2001). doi:10.1103/PhysRevB.63.245314. Gauge invariance in ETB. **Not read.**
17. B. Liu, Y. Zhong, Z.-X. Guo, X.-G. Gong, H. Xiang, arXiv:2609.37058 (2026). doi:10.48550/arXiv.2609.37058.
18. Y. Wang et al., GPUTB, arXiv:2509.06525 (2025); Computational Materials Today (2025), doi:10.1016/j.commt.2025.100039.
19. A. B. Lorentzen, F. Gao, P. Bøggild, A.-P. Jauho, M. Brandbyge, 2D Mater. (2024), doi:10.1088/2053-1583/ad349d; arXiv:2401.06556.
20. C. F. A. Negre et al., JCTC 12, 3063 (2016), doi:10.1021/acs.jctc.6b00154. Recursive factorization of S^-1 for linear scaling. Abstract only; volume/page from memory: **[UNCONFIRMED]**.
21. sisl, src/sisl/physics/state.py (`derivative`) and electron.py (`PDOS`), github.com/zerothi/sisl (main, read 2026-10-20).
22. C.-H. Lin, J. J. Moré, SIAM J. Sci. Comput. 21, 24 (1999). Eigen's IncompleteCholesky, cited in the legacy header. **Not verified here.**
