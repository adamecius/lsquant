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

Exposed on `lsquant::chebyshev::SpectralWindow` (`lsquant/core/chebyshev/spectral_window.hpp`):

    HalfWidth()   = BandWidth()/2
    ScaleFactor() = 1.0 / HalfWidth()
    ShiftFactor() = -BandCenter()/HalfWidth()

The moments are rescaled so `[E_lo, E_hi]` maps onto `[-1, 1]` directly: the analytic moment
identities checked by the chain oracle (`lsquant.oracle.chain`) — DOS moments `= δ_{m,0}`,
Kubo-Greenwood matrix `= B.15` — hold at that scaling. Spectral bounds come from the model
descriptor (`spectral_min`, `spectral_max`), a `BOUNDS` file, or Lanczos / Gershgorin
estimation (`lsquant::algebra::SpectralBounds`).

---

## 2. The α safeguard — on the reconstruction grids only

The reconstruction kernels carry `1/√(1 - x²)` (`delta_cheb_value`, `green_r_cheb`,
`dgreen_r_cheb` in `python/lsquant/spectral_functions.py`), which diverges at `x = ±1`. Two
failure modes follow if a reconstruction grid reaches the band edge:

1. **Recursion margin.** Eigenvalues at `±1` sit on the Chebyshev-recursion stability edge.
2. **Band-edge poisoning.** In the Kubo-Bastin route the output is a cumulative integral starting
   near the edge, so an edge NaN poisons the entire curve.

**Convention:** Reconstruct on `[-α, α]` with **`α = 0.95`** (`RECON_CUTOFF` in
`python/lsquant/reconstruct.py`), keeping the energy grid strictly inside `(-1, 1)` while the
Hamiltonian moments use the exact `[E_lo, E_hi] -> [-1, 1]` window.

---

## 3. Conductance / Hall unit convention

In natural units (`e = ħ = 1`), since `ħ = h/(2π)`, the conductance quantum is

    e²/h = 1/(2π) ≈ 0.159155

The extensive Kubo reconstructions divide by the total area (2D) or length (1D) to yield
intensive conductivities:

    σ_xy = C · e²/h,    C = (plateau / A) · 2π = Chern number,   A = N_cells · A_cell

On the Haldane honeycomb model (`examples/03_haldane_hall/`), the gap plateau gives `C = +1`.

---

## 4. Time evolution and spin dynamics — the physical-ℏ convention

Unlike the static spectral quantities above, time evolution runs in physical units:

    HBAR = 0.6582119624 eV·fs
    time t is in femtoseconds (fs)

The unitary evolution operator `Û(Δt) = e^{-iĤΔt/ħ}` is expanded in Chebyshev polynomials with
Bessel coefficients `(-i)^m (2 - δ_{m0}) J_m(ω₀ Δt) e^{-i ω_c Δt}`, where `ω₀ = HalfWidth / ħ`
and `ω_c = BandCenter / ħ`, evaluated in the engine via `ApplySeries`.
