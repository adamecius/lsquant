# HANDOFF — Golden-master references for the LinQT (lsquant) refactor: 1D chain analytics

## Mission context
We are refactoring `lsquant`/LinQT (real-space O(N) quantum transport, KPM + stochastic
trace). The refactor uses a strangler-fig strategy with **golden-master tests** as the
safety net. Your job is the golden-master suite. This package gives you the **analytical
ground truth** for the clean 1D monatomic tight-binding chain — every value below is
closed-form and symbolically verified (sympy), so the tests compare LinQT output against
exact numbers, not against another numerical run.

## Files in this package (all prefixed `lsquant_`)
- **`lsquant_new_physics.md`** — READ FIRST. Full physics reference: conventions, every
  closed form, the two-flavor analysis, and the test spec (its Section 5). Written in
  Spanish; formulas are language-neutral.
- **`lsquant_chain_reference.py`** — The oracle. Closed-form functions
  (`dos`, `v2`, `dos_moment`, `kg_moment`/`kg_moment_matrix`, `vac`, `msd`, `sigma_t`,
  `spin_precession`, `precession_frequency`) + a `__main__` self-test. Import this and
  use it directly as the golden master; do not re-derive.
- **`lsquant_derive_core.py`** — sympy proof of the DOS moments and the Kubo-Greenwood
  moment matrix. Run it to regenerate/verify.
- **`lsquant_derive_spin.py`** — sympy proof of the spin precession, the flavor
  equivalence, and the flavor-difference toy. Run it to regenerate/verify.

## Conventions (must match these exactly)
hopping `γ>0` (the thesis' `t`); lattice `a`; **time is `t`** (no clash);
`ε(k) = -2γ cos(ka)`; `ΔE = 2γ`, `Ē = 0`; `ε̃(k) = -cos(ka)`;
`v(k) = (2γa/ħ) sin(ka)`. All quantities **per site**.

## Closed-form oracle values (also in `lsquant_chain_reference.py`)

CHARGE SECTOR (spinless):
- DOS:            ρ(E) = 1/(π√(4γ²−E²)),  |E|<2γ
- velocity²:      v²(E) = (a²/ħ²)(4γ²−E²)
- DOS moments:    C̄ᴰᴼˢ_m = δ_{m,0}   (exact; only m=0 survives)
- KG moments:     C^KG_mn = (γ²a²/ħ²)[ δ_mn/(1+δ_m0+δ_m1) − ½(1−δ_mn)δ_{|m−n|,2} ]
                  → sparse: only diagonal and 2nd sub/super-diagonal nonzero,
                    values {1/2, 1, −1/2}·(γ²a²/ħ²).
- ballistic:      C_vv(E,t)=v²(E) (const);  ΔX²(E,t)=v²(E)t²;
                  σ(E,t) = (e²a²/πħ²)√(4γ²−E²)·t   (linear in t, NO DC limit)
- conductance:    flat plateau in energy (ρv = a/πħ is E-independent), value =
                  conductance quantum (overshoot only at van Hove edges E=±2γ).

SPIN SECTOR (magnetic chain H = H₀ − J_ex σ_z, initial spin +x̂):
- precession:     ⟨S_x⟩ = (ħ/2)cos(2J_ex t/ħ), ⟨S_y⟩ = −(ħ/2)sin(2J_ex t/ħ), ⟨S_z⟩=0
- frequency:      |f| = J_ex/(πħ)
- **E_F-independent and γ-independent.**

## The two flavors (CRITICAL for the spin tests)
LinQT computes the energy/time-resolved spin polarization in two flavors:
- `density-projection` (José): ⟨S⟩ = Tr[U†SU · P†δ(H−E_F)P] / Tr[P†δ(H−E_F)P].
  Polarization lives in the spin **projector** P.
- `state-projection` (Aron): ⟨S⟩ = [½⟨ψ(t)|s δ(E−H)|ψ(t)⟩ + h.c.] / ⟨ψ(t)|δ(E−H)|ψ(t)⟩,
  |ψ(0)⟩ = ½(𝟙+ĵ·ŝ)|φ_r⟩. Polarization lives in the **state**.

They share the same hot kernel ⟨v| s_H(t) T_m(H̃) |v⟩; only the source vector,
normalization, and (single-state) h.c. symmetrization differ. Two test consequences:

- **Equivalence test (exchange chain):** the spin block is k-independent, so BOTH flavors
  must give exactly cos(2J_ex t/ħ)·x̂ − sin(2J_ex t/ħ)·ŷ, E_F-independent, equal to each
  other within stochastic noise. Disagreement ⇒ the refactor introduced a spurious
  asymmetry.
- **Discrimination test (k-dependent field):** toy model H_k = ε(k)𝟙 + B(k)σ_z,
  B(k)=J_ex+δcos(ka). Then they MUST diverge:
    - density (sharp E_F): stays coherent, ω_dens(E_F) = 2(J_ex − δE_F/2γ)/ħ (1D shell is
      just ±k₀, same |B|).
    - state (energy window η): dephases (averages cos(2B(k)t/ħ) over the window).
  Assert (a) they coincide at δ=0 and (b) differ in a controlled way at δ≠0.

## Golden-master test list (mirrors `lsquant_new_physics.md` §5)
Compare LinQT output vs `lsquant_chain_reference.py`. Tolerances scale like the
stochastic error ~1/√(N·N_r).

CHARGE:
1. DOS-moments: C_m ≈ δ_{m,0}; |C_m|<tol for m≥1.
2. DOS reconstruction vs ρ(E) on a grid away from edges; rel. err < 1%.
3. KG matrix vs closed form (sparse pattern + values). ← fix the normalization here.
4. VAC↔MSD equivalence: σ(E,t) from VAC and from MSD must agree, slope
   (e²a²/πħ²)√(4γ²−E²).
5. Ballistic conductance: flat plateau in E; value = conductance quantum.

SPIN:
6. Precession (equivalence/consolidation): both flavors → cos/−sin, S_z=0,
   |f|=J_ex/(πħ), E_F-independent; flavor-to-flavor equal within noise.
7. Flavor discrimination (B(k) toy): coincide at δ=0, diverge controllably at δ≠0
   with density coherent at ω_dens(E_F) and state dephasing.

## OPEN DECISION — must resolve before freezing test #3
The KG moment normalization. The raw integral Tr[T_m V T_n V]/N gives C^KG₀₀ = 2;
the thesis B.15 closed form gives 1/2. They reconcile via the (1+δ_m0)(1+δ_n0) factor in
the DENOMINATOR of the thesis definition (eq. 4.7/B.15). Before freezing test #3, confirm
which normalization LinQT actually uses in its reconstruction path (thesis eq. 4.6).
`lsquant_chain_reference.kg_moment` currently uses the B.15 (normalized) convention. If
LinQT's reconstruction expects the raw convention, switch the oracle, not the test
tolerance.

## Suggested deliverable from you
A pytest module that (i) imports `lsquant_chain_reference` as oracle, (ii) builds the
small toy systems (clean chain; magnetic chain; B(k) toy) in LinQT, (iii) runs each
representation/flavor, (iv) asserts the 7 checks above. Start with #1, #3, #6 (cheapest
and most diagnostic).
