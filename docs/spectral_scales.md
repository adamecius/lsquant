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

**Fix.** Reconstruct only on `[-α, α]` with **`α = recon_cutoff = 0.95`**, defined alongside
`CUTOFF` in `include/chebyshev_moments.hpp`. This is **decoupled from `CUTOFF` on purpose**:

| quantity            | scale used      | why                                              |
|---------------------|-----------------|--------------------------------------------------|
| Hamiltonian moments | `CUTOFF = 1.00` | band must fill `[-1,1]` for the oracle identities |
| reconstruction grids| `safety_factors().recon_cutoff` (alpha = 0.95) | stay off the `1/√(1-x²)` edge singularity     |

Two reconstruction paths consume `α`:
- the **uniform grids** of the `*FromChebmom` drivers — `xbound = recon_cutoff`, so the energy
  grid runs over `[-α, α]` instead of `[-1, 1]`;
- the **FFT route** (`Kubo_solver_FFT_postProcess`) — the integration edge cushion is
  `safety_epsilon = 1 - recon_cutoff = 0.05`, so the Fermi-sea integral starts just inside the
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

---

## 5. Time evolution and spin dynamics — the physical-ℏ convention

Unlike the static spectral quantities above (natural units, e = ħ = 1), LSQUANT's **time
evolution runs in physical units**:

    HBAR = 0.6582119624 eV·fs        # Planck constant, defined in chebyshev_moments.hpp
    time t is in femtoseconds (fs)

The time-evolution operator `Û(Δt) = e^{-iĤΔt/ħ}` is built by a Bessel–Chebyshev expansion
(`MomentsTD::Evolve`), `Û = Σ_m (2-δ_{m0})(-i)^m J_m(ω₀Δt) T_m(H̃)`, with
`ω₀ = ChebyshevFreq() = HalfWidth/CUTOFF/ħ`. This is exact and unitary (‖φ(t)‖ is conserved to
machine precision). **Consequence:** any oscillation/relaxation frequency the code produces is
in **rad/fs with the physical ħ** — *not* the ħ=1 convention of the static oracle. Forgetting
this is the classic "wrong period" trap.

### Spin precession (the magnetic chain)
The magnetic chain is a constant Zeeman field `Ĥ = Ĥ_hop − J_ex σ_z` (J_ex = 0.1 eV). A spin
prepared along +x precesses in the x–y plane (thesis spin-dynamics formalism):

    |φ(0)⟩ = ½(𝟙 + x̂·ŝ)|φ_r⟩                                   # +x-projected state
    S(E,t) = ½⟨φ(t)|ŝ δ(E−Ĥ)|φ(t)⟩ / ⟨φ(t)|δ(E−Ĥ)|φ(t)⟩ + h.c.   # DOS-normalised ratio

    ⟨S_x(t)⟩ =  cos(ω t),  ⟨S_y(t)⟩ = −sin(ω t),  ⟨S_z(t)⟩ = 0,
    ω = 2·J_ex/ħ = 0.2/0.6582119624 = 0.30385 rad/fs   (period 20.68 fs),  E_F-independent.

Two practical points the `spin_precession` test relies on:
- **Energy-integrated shortcut.** When the observable is E-independent (as here), no energy
  reconstruction is needed: `∫dE` of the δ-reconstruction picks out `T₀`, so
  `⟨φ(t)|ŝ|φ(t)⟩` is exactly the **m=0 Chebyshev moment** `μ_ŝ(0,n)`, and the (constant)
  denominator `⟨φ(t)|φ(t)⟩ = μ_Sx(0,0)`.
- **Projection flavor.** The **density** path (`TimeEvolvedProjectedOperator`, both bra and
  ket carry the projector → `|φ(t)⟩⟨φ(t)|` with `|φ(0)⟩=P_x|φ_r⟩`) reproduces ⟨S(t)⟩ exactly.
  The **WF** path (`…ProjectedOperatorWF`) projects one side only — a structurally different,
  asymmetric correlator, *not* a drop-in for ⟨S(t)⟩.

### Ballistic transport (VAC / MSD)
The same expansion underlies the velocity-autocorrelation and mean-square-displacement routes.
For the clean chain, MSD(E,t) ∝ t² and the energy-resolved coefficient tracks the conductivity
slope `dσ/dt = (1/π)√(4−E²)` (the `msd_ballistic` test checks the √(4−E²) energy law by exact
trace). Note `inline_compute-kpm-CorrOp`'s state-file/`numStates` argument handling is currently
unreliable (`atoi("exact")=0`), so the VAC route is not yet a deterministic golden.

---

## 6. Where this is exercised (ctest)

| test | guards |
|---|---|
| `chain1d_analytic`, `chain1d_dos` | moment identities (§1): KG matrix B.15, DOS `δ_{m,0}` |
| `dos_reconstructed` | reconstructed DOS vs `1/(π√(4−E²))` on the α grid (§2) |
| `reconstruction_nan_guard` | Kubo-Bastin reconstruction finite everywhere (§2, Bug 2) |
| `hall_haldane` | quantized Hall plateau, Chern C=+1 (§3) |
| `msd_ballistic` | ballistic MSD, σ-slope energy law √(4−E²) (§5) |
| `spin_precession` | Larmor precession ω=2J_ex/ħ, S_z=0 (§5) |
| `graphene_regression`, `chain1d_reproducible`, `spin_algebra` | regression / determinism / su(2) algebra |
