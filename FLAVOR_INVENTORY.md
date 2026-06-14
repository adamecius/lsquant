# Flavor inventory — spin projection / time-evolution / magnetic-chain paths in `main`

Task D (discovery only, read-only). Input to the G2 spin tests (#6/#7) and the
consolidation decision. **Nothing was modified.** Evidence is `file:line` on
`main` at the time of writing.

## 1. The two projection flavors

| Flavor | Where | Evidence |
|---|---|---|
| **State / wavefunction projection** (project the source vector, asymmetric — "WF") | `chebyshev::TimeEvolvedProjectedOperatorWF` | [`src/chebyshev_solver.cpp:965`](src/chebyshev_solver.cpp#L965): `PhiR=State`, `PhiL=State`; `OPRJ.Multiply(PhiR,PhiL)` then `OP.Multiply(PhiL,PhiL)` → projector applied to build the **bra** only (ket evolves raw). Driver: [`src/inline_compute-kpm-TimeEvProjetedOp-Projecting_WF.cpp:79`](src/inline_compute-kpm-TimeEvProjetedOp-Projecting_WF.cpp#L79) (target `inline_compute-kpm-TimeEvProjetedOpWF`). |
| **Density projection** (projector on both sides, P\|φ⟩⟨φ\|P) | `chebyshev::TimeEvolvedProjectedOperator` | [`src/chebyshev_solver.cpp:1012`](src/chebyshev_solver.cpp#L1012): `PhiR=State`, `PhiL=PhiR`; `OPPRJ.Multiply(PhiR,PhiL)` then `copy(PhiL,PhiR)` → **both** bra and ket carry the projector (symmetric). Driver: [`src/inline_compute-kpm-TimeEvProjetedOp.cpp:79`](src/inline_compute-kpm-TimeEvProjetedOp.cpp#L79) (target `inline_compute-kpm-TimeEvProjetedOp`). |
| **Density projection in the device boundary** (spin projector applied inside the evolution) | `Device::spin_project(...)` | [`src/chebyshev_solver_Device.cpp:37`](src/chebyshev_solver_Device.cpp#L37): `chebMoms.device().spin_project(PhiR.data(), 2)`. Implementations: [`src/Device/Graphene_NNN.cpp:231`](src/Device/Graphene_NNN.cpp#L231), [`src/Device/Graphene_KaneMele.cpp:673`](src/Device/Graphene_KaneMele.cpp#L673); generic `projector(...)` [`src/Device/Graphene_KaneMele.cpp:894`](src/Device/Graphene_KaneMele.cpp#L894) (selected by `projector_option_`). **Device-only** (Graphene_NNN/KaneMele), not the 1D chain. |

**Caveat / uncertainty.** Neither solver function contains an *explicit* `+h.c.`
symmetrization line; the distinction here is structural — `…WF` projects one side
(state/wavefunction), the non-WF projects both sides (density). The exact mapping to
the thesis' `½(𝟙+ĵ·ŝ)|φ_r⟩ +h.c.` (state) vs `P|φ_r⟩` boundary (density) definitions
should be confirmed against the thesis before the G2 spin tests assert on it. The
device boundary path (`spin_project`) is the clearest "density projector in the
boundary" instance, but it lives in the **out-of-scope** Device/Graphene code.

## 2. Time-evolution driver — PRESENT ✅

Time-dependent machinery exists (`chebyshev::MomentsTD`, `Evolve()` in
[`src/chebyshev_momentsTD.cpp`](src/chebyshev_momentsTD.cpp)). Drivers:
`inline_compute-kpm-TimeEvOp`, `…TimeEvOpwithWF`, `…TimeEvProjetedOp`,
`…TimeEvProjetedOp-Projecting_WF`, `…TimeEvOp_Device`. So G2's spin-precession
(#6) has an evolution path to build on.

## 3. Magnetized / exchange 1D chain — ABSENT ❌ (G2 blocker)

`grep -rniE 'J_ex|Jex|exchange|magneti|sigma_z|zeeman'` over `src/` + `include/`
returns **nothing**. There is no `H = H₀ − J_ex σ_z` chain nor the `B(k)` toy model
in `main`. The clean `chain1d` (w2s) is spinless.

**Consequence for G2:** tests #6 (spin precession + flavor equivalence) and #7
(flavor discrimination on the `B(k)` toy) are **not reachable** as-is — they require
a new magnetized-chain model. Per scope that must arrive as **new files** in the
lsquant test area (or a new w2s generator), never as edits to existing physics. The
flavor *solver* paths (§1) and the time-evolution machinery (§2) exist; only the
**magnetic model** is missing.

## Summary for the consolidation decision
- Two distinct projection code paths exist (state/WF vs density), both reachable via
  their own executables — good (G2 #6/#7 can exercise both flavors **once a magnetic
  model exists**).
- Time evolution: present.
- Magnetic chain model: **missing** — the single prerequisite to build before G2.
- The clearest density-in-boundary projector is Device-only (out of scope this pass).
