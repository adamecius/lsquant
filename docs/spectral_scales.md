# Spectral scales, the α safeguard, and the conductance convention   {#spectral_scales}

Single reference for the Chebyshev rescaling, the reconstruction-grid `α` safeguard, and the
conductance / Hall unit convention. Natural units **e = ħ = γ = a = 1** unless stated.

---

## 1. The spectral rescaling — how the scales are defined

KPM needs the Hamiltonian rescaled so its spectrum sits inside the Chebyshev domain
`[-1, 1]`. With band bounds `[E_lo, E_hi]`:

    a = HalfWidth  = (E_hi - E_lo)/2          # BandWidth/2
    b = BandCenter = (E_hi + E_lo)/2
    x = (E - b)/a                              # rescaled energy, x in [-1, 1] for the band
    H̃ = (H - b·I)/a                            # rescaled Hamiltonian, eigenvalues = x

Exposed in `include/chebyshev_moments.hpp`:

    HalfWidth()   = BandWidth()/2
    ScaleFactor() = CUTOFF / HalfWidth()
    ShiftFactor() = -BandCenter()/HalfWidth() · CUTOFF

`CUTOFF` is a single constant, **`1.00`**. The moments must be rescaled so the band fills
`[-1, 1]` *exactly*: the analytic moment identities the test oracle checks — DOS moments
`= δ_{m,0}`, the Kubo-Greenwood matrix `= B.15` — hold only at that filling. **Do not** lower
`CUTOFF` to add an in-domain cushion; that shifts the moments off the oracle (it makes the
exact-trace golden tests fail for a real, physical reason). The cushion belongs on the
*reconstruction* grid instead — see §2.

The spectral bounds `a, b` come from a `BOUNDS` file if present, otherwise from a
Gershgorin-disk estimate padded by 10% (never the old arbitrary `[-100,100]`).

---

## 2. The α safeguard — on the reconstruction grids only

The reconstruction kernels carry `1/√(1 - x²)` (`delta_chebF`, `greenR_chebF`,
`DgreenR_chebF` in `chebyshev_coefficients.hpp`), which is **NaN at x = ±1**. Two failure
modes follow if a reconstruction grid reaches the band edge:

1. **Recursion margin.** Eigenvalues at exactly ±1 sit on the Chebyshev-recursion stability
   edge; rounding can push them past ±1.
2. **Band-edge poisoning.** In the Kubo-Bastin route the output is a *cumulative* integral
   starting at the edge, so one leading NaN poisons the whole curve (the historical all-NaN
   graphene / Haldane Kubo `.dat`).

**Fix.** Reconstruct only on `[-α, α]` with **`α = KPM_ALPHA = 0.95`**, defined alongside
`CUTOFF` in `include/chebyshev_moments.hpp`. This is **decoupled from `CUTOFF` on purpose**:

| quantity            | scale used      | why                                              |
|---------------------|-----------------|--------------------------------------------------|
| Hamiltonian moments | `CUTOFF = 1.00` | band must fill `[-1,1]` for the oracle identities |
| reconstruction grids| `KPM_ALPHA = 0.95` | stay off the `1/√(1-x²)` edge singularity     |

Two reconstruction paths consume `α`:
- the **uniform grids** of the `*FromChebmom` drivers — `xbound = KPM_ALPHA`, so the energy
  grid runs over `[-α, α]` instead of `[-1, 1]`;
- the **FFT route** (`Kubo_solver_FFT_postProcess`) — the integration edge cushion is
  `safety_epsilon = 1 - KPM_ALPHA = 0.05`, so the Fermi-sea integral starts just inside the
  edge.

Validated (Haldane t1=-1, t2=0.15, φ=π/2; 16×16; M=256; exact trace): the Kubo-Bastin driver
goes from **all-NaN (7679/7679)** to **0 NaN**, and the Chern plateau is preserved (`C ≈ +1`).
The `reconstruction_nan_guard` ctest pins this: the Kubo-Bastin reconstruction of the
committed graphene moments must be finite everywhere.

---

## 3. Conductance / Hall unit convention

Code is in natural units (`e = ħ = 1`). Since `ħ = h/2π`, the quantum is

    e²/h = 1/(2π) ≈ 0.159155        # this is where the 2π enters

The `*FromChebmom` reconstructions output **area-/length-extensive** quantities; divide by the
total area (2D) or length (1D) to get the intensive conductivity.

### Hall (Kubo-Bastin) on Haldane — quantized, the reference test
    σ_xy = C · e²/h,    C = (plateau / A) · 2π = Chern number,   A = N_cells · A_cell

Verified (Haldane 16×16, M=256, exact trace): gap plateau ≈ 106, `A_cell = 3√3/2 = 2.598`,
`A = 256·2.598 = 665.1` → **C = +1.004 ≈ +1** (the residual is finite-size / finite-M). This
pins the Kubo-Bastin prefactor `-2·SystemSize·ScaleFactor²` to `e²/h`. `lsquant_hall_reconstruct.py`
reproduces it end-to-end.

### Longitudinal (Kubo-Greenwood) — needs a DC limit
The perfectly clean chain has **no DC limit** — σ_DC diverges (ballistic); its bulk "plateau"
grows ∝ M. So the longitudinal prefactor is calibrated on a **disordered** chain where σ_DC
saturates (W=1, N=512, exact trace): the Greenwood σ_xx(E) matches the Chern-validated Bastin
route to ~5%, and the intensive `σ_xx(0) ≈ 40 e²/h` agrees with the Drude–Born estimate
~48/W². The Greenwood prefactor is `-SystemSize·ScaleFactor²` — half of Bastin's, since
Bastin's Fermi-sea term doubles the Fermi-surface contribution.

---

## 4. Two reconstruction bugs (both fixed)

**Bug 1 — Kubo-Greenwood ≡ 0 for centred bands.** `kuboGreenwoodFromChebmom.cpp` multiplied
the kernel by `ScaleFactor()·ShiftFactor()`; `ShiftFactor() = -BandCenter/HalfWidth·CUTOFF = 0`
for any band centred at E=0 (clean chain, graphene, every particle-hole-symmetric system), so
σ_KG was identically 0. **Fixed:** the second factor is a second `ScaleFactor()` (the prefactor
is `-SystemSize·ScaleFactor²`); `ShiftFactor()` remains, correctly, on the output energy axis.

**Bug 2 — Kubo-Bastin all-NaN (band-edge poisoning).** The reconstruction grid reached the
`1/√(1-x²)` singularity at x=±1 and the cumulative integral propagated the NaN. **Fixed** by
the `α` safeguard (§2): the grid never touches the edge.
