# Reconstruction architecture   {#reconstruction}

How KPM moment tensors become `.dat` spectra and conductivities in the unified C++/Python pipeline.
Every observable is evaluated as a word of static operators and spectral functions, compiled to a `Plan` executed by `lsquant::chebyshev::MomentEngine`, and reconstructed in Python (`python/lsquant/reconstruct.py`, `python/lsquant/kernels.py`, `python/lsquant/spectral_functions.py`).

---

## Unified reconstruction core — `python/lsquant/reconstruct.py`

- **Rank-1 density / spectral curves (`reconstruct_dos_curve`)** — contracts the 1D Chebyshev moments $\mu_m$ with the damped delta coefficients (`delta_cheb_value` with Jackson or Lorentz kernel) on an energy grid inside $[-\alpha, \alpha]$ ($\alpha = 0.95$).
- **Rank-2 Kubo conductivities (`reconstruct_kubo_greenwood`, `reconstruct_kubo_bastin`, `reconstruct_kubo_bastin_sea`, `reconstruct_kubo_bastin_surf`)** — contracts the 2D Chebyshev moment matrix $\mu_{mn}$ with the Green's function and derivative-Green's function Chebyshev series (`green_r_cheb`, `dgreen_r_cheb`) on either:
  - the uniform rescaled energy grid inside $[-\alpha, \alpha]$, or
  - Chebyshev nodes for fast Fourier / cosine-transform evaluation.
- **Spin Hall conductivity (`reconstruct_spin_hall_curve`, `reconstruct_spin_hall_bastin_grid`)** — evaluates the zero-temperature Bastin Fermi-sea integral of $\mathrm{Tr}[J_a^s\,\delta(E-H)\,V_b\,dG^-(E)/dE]$ using exact analytical coefficient antiderivatives or energy-grid quadrature.
- **Time-resolved observables (`python/lsquant/apps/stage4_time.py`)** — contracts rank-1 time-Indexed moments $\mu_m(t)$ for velocity autocorrelation, mean-square displacement, and projected spin dynamics.
