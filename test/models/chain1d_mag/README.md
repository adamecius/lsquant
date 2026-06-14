# `chain1d_mag` — magnetized 1D tight-binding chain (G2 prerequisite fixture)

Hand-authored Wannier90-style model for the **magnetic chain**

    H = -2 gamma cos(k a) ⊗ I_spin  -  J_ex sigma_z ,   gamma = 1, J_ex = 0.1

i.e. a 1-site, 2-orbital chain where orbital 1 = spin↑, orbital 2 = spin↓:
nearest-neighbour hop `-gamma` (diagonal in spin), on-site exchange `diag(-J_ex, +J_ex)`.
This is the model the G2 spin tests (#6 precession, #7 flavor discrimination) need —
[`FLAVOR_INVENTORY.md`](../../../FLAVOR_INVENTORY.md) flagged it as **absent from both
`main` and wannier2sparse**, so it is provided here.

## Files
- `chain1d_mag_hr.dat` / `.uc` / `.xyz` — the model (no generator; static text, no Python).

## Build the Hamiltonian (via wannier2sparse)
```
wannier2sparse chain1d_mag <NX> 1 1 VX -o operators     # -> operators/chain1d_mag.{HAM,VX}.CSR
```
**Verified** (NX=64 ⇒ 128×128): the HAM is Hermitian and its spectrum is the
exchange-split band `[-2.1, 2.1]` (↑ band `[-2.1, 1.9]`, ↓ band `[-1.9, 2.1]`,
splitting `2·J_ex = 0.2`) — exactly the analytic prediction.

## NOTE — upstream into wannier2sparse
This model should be added as a generator (`magnetic_chain()` / `chain1d_mag()`) in
**`wannier2sparse/examples/gen_models.py`**, next to `chain1d`, so it is produced and
band-edge–validated by w2s's own `validate.py` rather than carried as a hand-authored
fixture here.

**Spin operators are the open piece.** The HAM/VX build correctly, but w2s's spin
operators (`SX/SY/SZ`, `createHoppingSpinDensity_list`, `exact_spin_operator`) follow a
spinor / `.spn` convention that must be set up for this 2-orbital spin basis before the
G2 spin tests can assert on `⟨sigma_z(t)⟩`. That spin-operator support belongs in
wannier2sparse as well. Until it is confirmed, do **not** emit `SZ` for this model and
silently trust it (the project's "test fails for a non-physics reason" trap).
