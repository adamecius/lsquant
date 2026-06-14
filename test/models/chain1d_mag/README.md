# `chain1d_mag` — magnetized 1D tight-binding chain (G2 prerequisite fixture)

Hand-authored Wannier90-style model for the **magnetic chain**

    H = -2 gamma cos(k a) ⊗ I_spin  -  J_ex sigma_z ,   gamma = 1, J_ex = 0.1

A 1-site, 2-orbital chain. **Spin convention (sup/sdw):** the spin channel is set by the
orbital LABEL in `chain1d_mag.xyz` — `_s+_` = spin↑ (s_z=+1), `_s-_` = spin↓ (s_z=-1)
(w2s `map_id2spin`). Here orbital 1 = `A_s+_` (↑), orbital 2 = `A_s-_` (↓); nearest-
neighbour hop `-gamma` (diagonal in spin), on-site exchange `diag(-J_ex, +J_ex)`.
This is the model the G2 spin tests (#6 precession, #7 flavor discrimination) need —
[`FLAVOR_INVENTORY.md`](../../../FLAVOR_INVENTORY.md) flagged it as absent from both
`main` and wannier2sparse.

## Files
- `chain1d_mag_hr.dat` / `.uc` / `.xyz` — the model (static text, no Python, no generator).
- `operators/chain1d_mag.{HAM,VX,SZ}.CSR` — built via wannier2sparse (NX=64 demo).

## Build & verification
```
wannier2sparse chain1d_mag <NX> 1 1 SZ -o operators   # HAM, VX, SZ
```
Verified (NX=64 ⇒ 128×128):
- **HAM**: Hermitian, exchange-split spectrum `[-2.1, 2.1]` (↑ `[-2.1,1.9]`, ↓ `[-1.9,2.1]`,
  split `2·J_ex=0.2`) — matches analytics.
- **SZ**: a proper sigma_z — eigenvalues `±1`, `S_z² = I`, `[H, S_z] = 0`, traceless.

## IMPORTANT — `SZ` requires the wannier2sparse spin-density fix
The committed `SZ.CSR` was produced with **wannier2sparse PR #12**
(`fix/onsite-spin-density`). The currently *pinned* w2s SHA (`511f4db`,
`scripts/fetch_wannier2sparse.sh`) has a bug where `S_z` re-fills the chain's ±γ bonds
with ±1 (not involutory; eigenvalues escape ±1). **HAM/VX are unaffected** by that bug.
→ When PR #12 merges, **bump the pinned w2s SHA**; until then, do not regenerate `SZ`
with the pinned binary and trust it.

## Open gap for G2 transverse spin (precession #6) — upstream into wannier2sparse
`sigma_x / sigma_y` come out **identically zero** for this collinear model: w2s derives
spin-operator support from H's on-site block, which here is diagonal (no up↔down term),
so there is nothing to host `σ_x/σ_y`. Building `σ_x/σ_y` support from the spin doubling
itself (independent of H) is needed before #6/#7 can measure `⟨σ_x(t)⟩`. This — and the
`magnetic_chain()` generator (so w2s's `validate.py` band-edge–checks it) — belong
upstream in `wannier2sparse/examples/gen_models.py`; tracked in w2s PR #12's notes.
