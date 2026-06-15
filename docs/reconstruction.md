# Reconstruction architecture (Phase 5)   {#reconstruction}

How KPM moment tables become `.dat` spectra/conductivities, after the Phase-5 unification.
"A formula is a piece of data, not an executable" — the reconstruction logic lives in one
place; each driver is a thin wrapper that picks a formula (and grid backend) and writes a file.

---

## Unified core — `include/observable.hpp`, `src/reconstruction.cpp`

- **`KuboObservable`** — the formula as data: `{green-kernel, prefactor, num_div_mult, integrate}`.
  `KUBO_GREENWOOD = {greenR, −1, 1, false}` is the Fermi-surface (N×1) degenerate case of
  `KUBO_BASTIN = {DgreenR, −2, 30, true}`.
- **`reconstruct_kubo(mu, obs, grid)`** — the shared double-Chebyshev integrand
  `Σ_{m,n} delta_chebF(E,m) Im[green(E,n) μ_{mn}]`, with two grid backends:
  - `GRID_UNIFORM` — 30·M points on `[-α, α]` (α = `safety_factors().recon_cutoff`),
  - `GRID_FFT_NODES` — M Chebyshev nodes (rearranged ascending for the Bastin integral).
- **`reconstruct_density_grid(μ_real, α, num_div)`** — the shared 1-D DOS/spectral inner sum
  `Σ_m delta_chebF(E,m) μ_m`. Callers add their scalar prefactor and physical-energy axis.
- The **`recon_kernel_oracle`** test (Phase R) *drives* `reconstruct_density_grid` with a
  synthetic δ of known moments, so the production reconstruction is guarded operator-free.

## Drivers — what is unified, and the gate

| driver | formula / route | gated by |
|---|---|---|
| `inline_kpm-DOS-standalone` | DOS, `reconstruct_density_grid` | `dos_reconstructed`, `recon_kernel_oracle` |
| `spectralFunctionFromChebmom` | spectral, density grid | byte-identity (1-D moments) |
| `spectralFunctionFromChebmom_Lorentz` | spectral + Lorentz kernel | byte-identity |
| `kuboGreenwoodFromChebmom` | `reconstruct_kubo(KUBO_GREENWOOD, UNIFORM)` | `greenwood_regression` |
| `kuboBastinFromChebmom` | `reconstruct_kubo(KUBO_BASTIN, UNIFORM)` | `hall_haldane`, `reconstruction_nan_guard` |
| `kuboGreenwoodFromChebmom_FFTgrid` | `KUBO_GREENWOOD, FFT_NODES` | `fft_grid_regression` |
| `kuboBastinFromChebmom_FFTgrid` | `KUBO_BASTIN, FFT_NODES` | `fft_grid_regression` |

All seven are byte-identical to their pre-refactor binaries.

---

## Secondary reconstructions — NOT yet unified, NO regression golden

These are alternative / experimental formulations with **distinct integrands** (not just
`green` vs `Dgreen`), and none has a golden test. They are **out of scope for the structural
unification** and are flagged for a post-refactor evaluation: for each, decide *unify* (needs a
generalized `Integrand` and a golden), *rewrite*, or *remove*. Do not rely on any of them until
it has a regression test.

| driver | distinct integrand | note |
|---|---|---|
| `kuboBastinIFromChebmom` | Greenwood integrand `delta·Im(greenR·μ)` + trapezoid output | near-duplicate of Greenwood |
| `kuboBastinIIFromChebmom` | `−Re[(GrL·DGrR − DGrL·GrR)·μ]`, prefactor `/2π` | products of green & d-green on both indices |
| `kuboBastinSeaFromChebmom` | `−Re[(Gr−Ga)(DGr+DGa)·μ]` | retarded/advanced Fermi-sea form |
| `kuboBastinSurfFromChebmom` | `delta·delta·μ`, prefactor `π·N·SF²` | double-delta Fermi-surface form |
| `opticalConductivityFromChebmom` | frequency-shifted `greenR(x ± w)` | finite-frequency; needs a frequency arg |
| `localDensitiesFromChebmom` | local (site-resolved) density | different observable shape |
| `timeCorrelationsFromChebmom` | time-domain `MomentsTD` projection | belongs with the TD family |
| `spectralFunctionFromChebmom_FFTgrid` | 1-D density on `cos(π(i+½)/M)` nodes, writes `abs(ρ)` | quirky node formula + abs |

Compute-side: `Kubo_solver_FFT` still carries its own `enum formula {KUBO_GREENWOOD, KUBO_BASTIN}`
(the FFT solver that computes *and* reconstructs); folding it onto the shared formula identity is
part of the same evaluation, and is gated only once that path has a test.

### A latent finding for the evaluation
The **uniform-grid** Kubo drivers write the energy axis as `E/ScaleFactor − ShiftFactor`
= `E·HalfWidth + BandCenter/HalfWidth`, whereas the **FFT-grid** and DOS drivers use
`E·HalfWidth + BandCenter`. These agree only for a band centred at E=0 (all current goldens), so
the `+BandCenter/HalfWidth` term looks like a latent bug for off-centre bands. Confirm and fix
under a non-centred golden during the evaluation.
