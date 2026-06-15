# Tutorial 3 - Quantum Hall plateau of the Haldane model

This tutorial computes the Hall and longitudinal conductivities of the Haldane
model with LinQT. It shows three things: the Hall conductivity locks onto a
quantized plateau in the gap, that plateau comes entirely from the Fermi sea,
and the longitudinal conductivity vanishes in the gap while it stays finite in
the bands.

This tutorial assumes you are inside the repository, in this example folder,
with LinQT built and installed (executables on `PATH`, or `LSQUANT_BIN` set to
the build directory). Read Tutorials 1 and 2 first.

## The physics

The Haldane model lives on a honeycomb lattice with a real nearest-neighbour
hopping `t1`, a complex next-nearest-neighbour hopping `t2 e^{i phi}` that
carries a staggered flux, and an optional sublattice mass `M`:

```
H = t1 * sum_<ij> c_i^dag c_j  +  t2 * sum_<<ij>> e^{i nu_ij phi} c_i^dag c_j
    + M * sum_i xi_i c_i^dag c_i,      nu_ij, xi_i = ±1.
```

At `t1 = -1`, `t2 = 0.15`, `phi = pi/2`, `M = 0` the lower band carries Chern
number `C = +1`, so the model is a Chern insulator [5]. With the Fermi level in
the gap the Hall conductivity is quantized,

```
sigma_xy = C * e^2 / h,
```

a bulk topological response. LinQT measures it from the velocity-velocity
Chebyshev moments through the Kubo-Bastin route [2].

The Kubo-Bastin conductivity splits into two pieces, a Fermi-surface part built
from states at the Fermi level and a Fermi-sea part built from every filled
state below it:

```
sigma_xy = sigma^surface + sigma^sea.
```

The Kubo-Greenwood conductivity is a Fermi-surface quantity, so it gives the
longitudinal `sigma_xx` and stays zero where the Fermi level sits in a gap [2].
Inside the gap the Fermi surface is empty, so the surface part of `sigma_xy`
vanishes there and the whole quantized plateau comes from the Fermi sea. This is
why a topological insulator carries a Hall current while it carries no
longitudinal current.

## Step 1 - get the operators

The committed Haldane operators are the verified model (`C = +1`). Stage them
into this folder with the helper, which copies the operators and writes the
descriptor and an exact-trace state set:

```bash
python make_haldane.py
```

This writes `operators/haldane.{HAM,VX,VY}.CSR` (copied from
`test/golden/haldane`), `operators/haldane.HAM.desc` (the physical sidecar for
the Hamiltonian), and `exact`. The descriptor carries the spectral bounds
(`spectral_min = -4.5`, `spectral_max = 4.5`), which enclose the spectrum
(`|E|max ~ 3`). `lsquant` reads these bounds straight from the descriptor, so
this run needs no BOUNDS file. The exact-trace state set evaluates `Tr[...]/dim`
over the full basis, so the moments are deterministic.

Confirm the descriptor and the run before computing:

```bash
lsquant inspect operators/haldane.HAM.desc
```

This prints the observable, the units, the provenance, and the spectral bounds
`[-4.5, 4.5]` that `lsquant compute` will use.

## Step 2 - the velocity-velocity moments

Compute the off-diagonal moments `VX-VY` for the Hall response and the diagonal
`VX-VX` for the longitudinal response through the single `lsquant` binary. Each
run is a small `run.json`; this is the heavy step, so it runs in C++.

```bash
cat > run_hall.json <<'JSON'
{ "label": "haldane", "operator_right": "VX", "operator_left": "VY",
  "num_moments": 128, "kernel": "jackson", "state": "exact" }
JSON
cat > run_long.json <<'JSON'
{ "label": "haldane", "operator_right": "VX", "operator_left": "VX",
  "num_moments": 128, "kernel": "jackson", "state": "exact" }
JSON
lsquant compute --config run_hall.json
lsquant compute --config run_long.json
```

Each call reads `operators/haldane.{HAM,VX,VY}.CSR`, takes the spectral bounds
from `operators/haldane.HAM.desc`, runs the two-operator Chebyshev recursion
with 128 moments per direction, and writes a `NonEqOp...chebmom2D` moment file.

## Step 3 - the conductivities

The reconstruction sums run in C++ through the single `lsquant` binary for the
totals and through the dedicated drivers for the Fermi-sea and Fermi-surface
parts:

```bash
lsquant reconstruct NonEqOpVX-VYhaldane*.chebmom2D bastin 10        # sigma_xy total
inline_kuboBastinSeaFromChebmom  NonEqOpVX-VYhaldane*.chebmom2D 10  # Fermi sea
inline_kuboBastinSurfFromChebmom NonEqOpVX-VYhaldane*.chebmom2D 10  # Fermi surface
lsquant reconstruct NonEqOpVX-VXhaldane*.chebmom2D greenwood 10     # sigma_xx
```

The Kubo-Bastin route stays finite at the band edges through the alpha safeguard
on the reconstruction grid (Tutorial 1, and `docs/spectral_scales.md`). The
wrapper runs all of this, reads the Chern number off the gap plateau, and draws
the figure:

```bash
python lsqhall.py
```

This writes `fig_haldane_hall.png` with two panels.

**Top panel - Hall conductivity in units of `e^2/h`.** The total Kubo-Bastin
curve sits on a flat plateau at `1` through the gap, so the Chern number reads
`C = (plateau / A) * 2 pi = +1` with `A = N_cells * A_cell` and
`A_cell = 3 sqrt(3) / 2`. The Fermi-sea curve carries the entire plateau. The
Fermi-surface curve stays at zero in the gap and turns on once the Fermi level
enters the bands. The plateau is a Fermi-sea property.

**Bottom panel - longitudinal conductivity.** The Kubo-Greenwood `sigma_xx`
sits at zero through the gap, the insulating interior, and rises inside the
bands where states are available at the Fermi level. The gap region carries a
quantized Hall current and a vanishing longitudinal current at the same time.

## Step 4 - why the plateau is quantized

The plateau value is the Chern number of the filled band, a bulk k-space
invariant. The companion script builds the two-band Bloch Hamiltonian, reads the
Chern number off a k-grid, and sweeps the sublattice mass to show the
topological transition:

```bash
python haldane_kspace.py
```

This reports `C = +1` at the model parameters and `C = 0` once the flux is
removed, and writes `fig_phase_diagram.png`. Raising `M` closes the gap at
`M = 3 sqrt(3) t2 sin(phi) ~ 0.78 eV`, where the Chern number drops from `+1` to
`0`. Past that point the real-space `sigma_xy` plateau also disappears: the
bulk invariant and the transport plateau move together.

## What to take away

- The Hall conductivity of the Haldane model locks onto `C e^2 / h` in the gap,
  the quantized topological response [2, 5].
- The plateau comes entirely from the Fermi sea; the Fermi-surface part stays at
  zero through the gap.
- The longitudinal `sigma_xx` from Kubo-Greenwood vanishes in the gap and turns
  on in the bands, the insulator-then-metal behaviour.
- The plateau value equals the k-space Chern number, and both follow the same
  topological transition as the gap closes.

## References

[1] A. Weiße, G. Wellein, A. Alvermann, H. Fehske,
*The kernel polynomial method*, Rev. Mod. Phys. **78**, 275 (2006).
https://journals.aps.org/rmp/abstract/10.1103/RevModPhys.78.275

[2] J. H. García et al., *LinQT / linear-scaling quantum transport*,
arXiv:1811.07387. https://arxiv.org/abs/1811.07387

[5] F. D. M. Haldane, *Model for a quantum Hall effect without Landau levels*,
Phys. Rev. Lett. **61**, 2015 (1988).
