# Tutorial 1 - Density of states of a 1D chain

This tutorial builds a one-dimensional tight-binding chain, computes its
Chebyshev moments with LinQT, and reconstructs the density of states (DOS). It
shows two things that hold for every later calculation: the moments are an
intensive, size-independent fingerprint of the spectrum, and the resolution of
the reconstructed DOS follows the number of moments.

This tutorial assumes you are inside the repository, in this example folder,
with LinQT already built and installed (the executables are on your `PATH`, or
the directory that holds them is in `LSQUANT_BIN`). The Python helpers need
`numpy` and `matplotlib` (`pip install -r requirements.txt`).

## The physics

The chain has `N` sites, one orbital per site, periodic boundaries, and a
single hopping `t = -1 eV`:

```
H = t * sum_i ( |i><i+1| + |i+1><i| ),   index i+1 taken modulo N.
```

The eigenvalues are `E_k = 2 t cos(2 pi k / N)`, so the spectrum fills the band
`[-2, 2] eV`. In the large-`N` limit the density of states takes the closed
form

```
rho(E) = 1 / ( pi * sqrt(4 - E^2) ),
```

with van Hove peaks at the band edges `E = ±2 eV`.

LinQT reaches large systems by expanding spectral functions in Chebyshev
polynomials of the rescaled Hamiltonian, the kernel polynomial method (KPM)
[1]. The DOS is the spectral function of the identity operator, so its
expansion coefficients are the moments

```
mu_m = (1/N) Tr T_m(H~),     H~ = (H - b) / a,
```

where `T_m` is the Chebyshev polynomial of order `m` and `(a, b)` rescale the
band into `[-1, 1]`. For this chain the trace gives `mu_m = delta_{m,0}`: the
moments collapse onto a single nonzero value at `m = 0`, which is exactly the
Chebyshev weight that reconstructs `rho(E)` [1].

## Step 1 - build the chain

Build three sizes. The generator writes the sparse Hamiltonian, its physical
sidecar, the spectral bounds, and an exact-trace state set.

```bash
python make_chain.py 256
python make_chain.py 1024
python make_chain.py 4096
```

Each call writes, for label `chain1d_N<N>`:

```
operators/chain1d_N<N>.HAM.CSR    the Hamiltonian, complex-CSR text
operators/chain1d_N<N>.HAM.desc   the physical sidecar
BOUNDS                            the (a b) = (-2 2) line read at compute time
exact_chain1d_N<N>                the exact-trace state set
```

The `.desc` sidecar holds what makes the matrix a physical object: the
observable identity (`hamiltonian`), units (`eV`), basis tags (one `s` orbital
per site), the spectral bounds, and the provenance. The numerical CSR holds the
matrix. They travel together, so a result stays traceable to the operator that
produced it. The sidecar stays beside the kernel and keeps clear of the
Chebyshev recursion.

## Step 2 - compute the moments

Run the recursion for the smallest chain with 512 moments. Operator label `1`
selects the identity, so the spectral function is the DOS.

```bash
inline_compute-kpm-spectralOp chain1d_N256 1 512 exact_chain1d_N256
```

This reads `operators/chain1d_N256.HAM.CSR` and `BOUNDS`, runs the threaded
Chebyshev recursion, and writes a moments file:

```
SpectralOp1chain1d_N256KPM_M512_stateexact_chain1d_N256.chebmom1D
```

The state label in the name is the state-file name you passed (`exact_chain1d_N256`).
The file starts with a header line `SystemSize BandWidth BandCenter`, then the
moment count, then the moments as `real imag` pairs. The exact-trace state set
makes LinQT evaluate `sum_i <e_i|T_m|e_i>/N = Tr[T_m]/N` over the full basis, so
the moments are deterministic and equal the analytic trace.

## Step 3 - the moments are intensive

The compute step is the expensive part, so it stays in C++. Reading the few
hundred moments back is light, so Python reads them directly and plots. Produce
the moments figure:

```bash
python -c "import lsqplot; lsqplot.figure_moments([256,1024,4096], M=256)"
```

This writes `fig_moments.png`. The three curves land on top of one another at
`mu_m = delta_{m,0}`. The moments measure a per-site (intensive) property, so
they keep their value as the system grows. A finite chain enters through
aliasing: the trace of `T_m` revives at `m = N`, so a size `N` shows extra
weight once `m` reaches `N`. For `m < N` the three sizes agree to machine
precision.

## Step 4 - reconstruct the DOS

The reconstruction evaluates

```
rho(E) = sum_m g_m mu_m T_m(E)
```

over a fine energy grid and every moment, with `g_m` the Jackson kernel
coefficients [1]. This sum is heavy, so it stays in C++. Run it once at a
120 meV broadening (the binary is `inline_spectralFunctionFromChebmom` in the
build directory; an install may expose it as `spectralFunctionFromChebmom`):

```bash
inline_spectralFunctionFromChebmom SpectralOp1chain1d_N256KPM_M512_state*.chebmom1D 120
```

This writes `mean...JACKSON.dat` with two columns, energy in eV and the total
DOS. The driver returns the extensive DOS `N * rho(E)`; dividing by the system
size `N` (from the moment-file header) gives the per-site density that matches
`rho(E) = 1/(pi sqrt(4 - E^2))` and overlaps across sizes — `lsqplot.py` does
this division. The kernel keeps `M_eff = ceil( pi * BandWidth / (2 * broadening) )`
moments, which sets a broadening

```
eta = pi * HalfWidth / M_eff.
```

For this chain `HalfWidth = 2 eV`, so `eta = 2 pi / M_eff eV`. The broadening
follows directly from how many moments the kernel keeps [1].

## Step 5 - size, spacing, and resolution

Produce the two-panel figure:

```bash
python lsqplot.py
```

This writes `fig_dos.png`.

A finite sample of `N` sites carries `N` energy levels, a finite spectrum with
mean spacing `Delta ~ BandWidth / N`. Exact-trace moments encode that finite
spectrum through one feature: the trace of `T_m` revives at `m = N` (the
aliasing term). So the discrete levels of a chain of size `N` become visible
once the kernel keeps `M_eff >= N` moments; for fewer moments the chain shows
the smooth bulk form. With `N_small = 256` that threshold is a broadening near
24 meV.

**Top panel - DOS at fixed broadening (20 meV, `M_eff ~ 315`), three sizes.**
At this broadening the kernel resolves below the level spacing of the `N = 256`
chain (`M_eff > 256`), so its curve carries the discrete structure of its 256
levels. The `N = 1024` and `N = 4096` chains keep `M_eff < N`, so they follow
the smooth bulk form `rho(E) = 1/(pi sqrt(4 - E^2))` with its van Hove peaks at
`±2 eV`. The curve sharpens toward the bulk form as `N` grows.

**Bottom panel - DOS at fixed size (`N = 256`), three broadenings
(80, 24, 10 meV).** Here the chain stays fixed and the broadening sweeps across
the `M_eff = N` threshold. At 80 meV (`M_eff ~ 79`) the curve is the smooth
bulk band; at 24 meV (`M_eff ~ 262`) the discrete levels begin to emerge; at
10 meV (`M_eff ~ 512`) they stand out one by one. The legend prints the
`M_eff` behind each curve, which makes the link explicit: in KPM the broadening
is the number of moments [1].

## What to take away

- The DOS moments are intensive, so the same expansion describes every size.
- A finite sample has a finite, discrete spectrum; its level spacing sets the
  maximum meaningful resolution.
- The KPM broadening equals `pi * HalfWidth / M`, so the moment count is the
  resolution knob [1].
- The operator and its sidecar stay together, so each DOS stays traceable to
  the chain that produced it.

The next tutorial reuses this chain and its velocity operator to compute a Kubo
conductivity, where the same moments feed a two-dimensional expansion.

## References

[1] A. Weiße, G. Wellein, A. Alvermann, H. Fehske,
*The kernel polynomial method*, Rev. Mod. Phys. **78**, 275 (2006).
https://journals.aps.org/rmp/abstract/10.1103/RevModPhys.78.275

[2] J. H. García et al., *LinQT / linear-scaling quantum transport*,
arXiv:1811.07387. https://arxiv.org/abs/1811.07387
