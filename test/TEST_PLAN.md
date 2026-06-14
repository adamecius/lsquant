# Test plan & handoff — evaluating LinQT against the analytic oracle

This document is the **entry gate** for the agent continuing the test suite. It records
(a) what is now green, (b) a brittle test that was fixed and why, and (c) the precise,
oracle-grounded specification for every remaining test (#2, #4, #5, #6, #7), with the
**verified target numbers** so no reference is hand-waved.

Ground truth is `lsquant_chain_reference.py` (the oracle, symbolically proven in
`lsquant_derive_core.py` / `lsquant_derive_spin.py`). Units throughout: `γ=a=ħ=e=1`.
The clean 1D chain band is `[-2,2]`; `J_ex=0.1`.

## Current status (ctest)

lsquant **5/5**, wannier2sparse **20/20** on a clean toolchain (g++ 13, Eigen apt, no MKL).

| ctest | spec | what it guards |
|---|---|---|
| `chain1d_analytic` | #3 KG-matrix | committed KG moments vs closed form B.15, tolerance |
| `chain1d_dos` | #1 DOS-moments | committed DOS moments vs `δ_{m,0}`, tolerance |
| `chain1d_reproducible` | #1/#3 | **determinism** (two runs byte-identical) + **fresh-run physics** (tolerance) |
| `spin_algebra` | #6 prereq | magnetic-chain `S_x,S_y,S_z` form su(2); `S_x≠0` and `[H,S_x]≠0` |
| `graphene_regression` | regression | fresh graphene moments vs committed golden (tol) + **bounds safeguard** active |

New test sources: `test/spin_algebra_test.cpp`, `test/chebmom_compare.cpp` (generic
numeric moment diff), `test/graphene_reproduce.sh`.

## Fixed: `chain1d_reproducible` was failing for a non-physics reason

It byte-compared a fresh run against the **committed** golden. Across toolchains the
moments agree only to ~`1e-15` (floating-point reduction order in compiler/Eigen/OpenMP),
so the byte diff failed though every element matched to `<1e-12` and the KG structure was
identical. This is the project's "test fails for a non-physics reason" trap.

Fix (principled, not a tolerance bump on physics): the reproducibility test now keeps the
**byte-identical** check only where it is legitimate — two runs of the **same** binary
(RNG determinism) — and validates the fresh run's **physics with the existing tolerance
comparators** instead of byte-diffing it against a foreign-toolchain golden. The committed
golden's physics is still guarded (with tolerance) by `chain1d_analytic`/`chain1d_dos`.

## Remaining tests — spec, verified targets, driver path, caveats

These need **reconstruction / time-evolution** drivers (not just moments), which carry
kernel/broadening/energy-grid choices and are more fragile than the moment-level checks.
The oracle reference numbers below are all verified.

### #2 — DOS reconstructed (charge)
Assert the KPM-reconstructed DOS equals `ρ(E)=1/(π√(4-E²))` with **relative error < 1%**
on a grid **away from the band edges** (edges diverge; KPM smears them).
Verified targets: `ρ(0)=0.159155`, `ρ(0.5)=0.164375`, `ρ(1.0)=0.183776`, `ρ(1.5)=0.240620`.
Driver: `inline_kpm-DOS-standalone` (reconstructs and writes the `.dat` curve). Use exact
bounds `[-2,2]`, a Jackson kernel, enough moments (e.g. M≥512) and `N_r` so the stochastic
floor is below 1%. Compare the `.dat` columns to `ρ(E)` on `|E|≤1.6`.

### #4 — VAC ↔ MSD equivalence (charge)
Assert `σ(E,t)` from the velocity autocorrelation and from the mean-square displacement
**agree** (representation equivalence) and grow **linearly in t** with slope
`(e²a²/πħ²)√(4γ²-E²)`. Verified targets (oracle, exact): slope at `E=0` is
`(1/π)√4 = 0.636620`; at `E=1`, `0.551329`. The oracle's `sigma_t(0,t)` is exactly
`0.636620·t` (e.g. `t=10→6.3662`, `t=40→25.4648`), confirming VAC↔MSD consistency.
Drivers: the MSD and correlation paths (`inline_compute-kpm-MeanSquareDisplacement*`,
`inline_compute-kpm-CorrOp`). The assertion is two-sided: VAC-σ vs MSD-σ within the
stochastic tolerance, and each vs the analytic slope.

### #5 — Ballistic conductance (charge)
Assert a **flat plateau in energy** whose value is the conductance quantum (per channel,
in units `2e²/h`; the 1-channel clean chain gives one quantum). Check flatness (relative
variation across an interior energy window below a small tolerance) and the plateau value.
Driver: the Kubo/conductance reconstruction path. Caveat: define the unit convention
explicitly (the spec fixes normalization in #3 — reuse it).

### #6 — Spin precession (spin) — prerequisites now satisfied
Assert `⟨S(t)⟩ = cos(ω t) x̂ − sin(ω t) ŷ`, `S_z=0`, with **angular** frequency
`ω = 2J_ex/ħ = 0.200`; equivalently the oracle's ordinary frequency
`precession_frequency(0.1) = J_ex/(πħ) = 0.031831` (`f=ω/2π` — state both to avoid the
ω-vs-f trap). The result must be **independent of E_F**, and the two projection flavors
must agree flavor-to-flavor within the stochastic noise.
Oracle samples: `⟨S_x⟩(0)=1`, at `t=π/(2ω)=7.854 → S_x=0, S_y=−1`.
Prereq DONE: `spin_algebra` proves `S_x,S_y,S_z` are su(2) and `[H,S_x]≠0` (`=0.2`, the
precession source). Remaining work: prepare the x-polarized initial state, drive
`MomentsTD`, and **resolve the state-vs-density mapping** (`FLAVOR_INVENTORY.md`) against
the thesis *before asserting* — neither solver path has an explicit `+h.c.` line, so which
of `TimeEvolvedProjectedOperator(WF)` is "state" vs "density" must be pinned first. **This
is a physics decision for the human, not the agent.**

### #7 — Flavor discrimination (spin) — model still missing
Toy `B(k)=J_ex+δ cos(ka)`: **density** stays coherent at
`ω(E_F)=2(J_ex − δE_F/2γ)/ħ`; **state** dephases. Assert (a) the two flavors **agree at
δ=0** and (b) **differ controllably at δ≠0**. Blocker: this `B(k)` model does **not exist**
in w2s or lsquant yet — it must be added as a generator
(`wannier2sparse/examples/gen_models.py`), analogous to `magnetic_chain()`, before #7 is
reachable.

## Recommended order
1. **#2** (DOS reconstructed) — cheapest reconstruction test; validates the KPM
   reconstruction path against a closed form. Low risk.
2. **#4** (VAC↔MSD) — exercises the transport reconstruction and the representation
   equivalence; the analytic slope is a strong, exact target.
3. **#5** (ballistic conductance) — needs the conductance path + a unit-convention call.
4. **#6** (precession) — only after the human pins the state-vs-density mapping.
5. **#7** (flavor discrimination) — only after the `B(k)` model is added upstream.

## Principles carried from prior work
- Reference numbers come from the oracle, never from another run or a fitted tolerance.
- Byte comparison is legitimate **only** for determinism (same binary, twice); cross-run /
  cross-toolchain comparisons use a numerical tolerance above ~`1e-12` noise and well below
  the stochastic-trace error.
- A test that "fails for a non-physics reason" is a bug in the test, not a reason to relax
  the physics.
