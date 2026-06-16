# LinQT — Linear-scaling Quantum Transport

**LinQT** (repository `lsquant`) is a C++ library for **large-scale quantum
transport** in tight-binding models. It computes spectral and transport
properties — density of states, conductivities, time evolution, spin dynamics —
by expanding them in **Chebyshev moments of the Hamiltonian** (the Kernel
Polynomial Method, KPM). The cost grows *linearly* with system size and nothing
is ever diagonalized. The numerics run on a single Eigen backend, matrix-free and
OpenMP-parallel.

If you have a tight-binding model as sparse operators in CSR format, LinQT turns
it into a density of states or a Kubo conductivity in two steps: **compute the
moments**, then **reconstruct the spectrum**.

## What it computes

LinQT expands a spectral or transport quantity in Chebyshev polynomials of the
rescaled Hamiltonian, so its cost scales linearly with system size (KPM [1]):

- **Density of states and spectral functions**, from the 1D Chebyshev moments
  $\mu_m = \mathrm{Tr}\,T_m(\tilde H)$ [1].
- **Kubo–Greenwood and Kubo–Bastin conductivities** — longitudinal $\sigma_{xx}$
  and Hall $\sigma_{xy}$, charge and spin, with the Fermi-sea / Fermi-surface
  split [2].
- **Time evolution and time-dependent correlations** — velocity autocorrelation
  and mean-square displacement, for ballistic transport and localization [2].
- **Spin dynamics** — precession and relaxation; non-orthogonal bases [2].

The methods and their derivations are reviewed in [2]; the KPM kernel itself in
[1].

#### References

1. A. Weisse, G. Wellein, A. Alvermann, H. Fehske, *The kernel polynomial
   method*, Rev. Mod. Phys. **78**, 275 (2006).
2. Z. Fan, J. H. García, A. W. Cummings, J. E. Barrios-Vargas, M. Panhans,
   A. Harju, F. Ortmann, S. Roche, *Linear Scaling Quantum Transport
   Methodologies*, arXiv:1811.07387.

## Installation

Requirements: a C++11 compiler (e.g. g++ 13) and CMake ≥ 3.16, plus **Eigen3**
(≥ 3.3), **FFTW3**, **GSL**, **Boost.Math** (header-only), **OpenMP**, and
**pthreads**. Intel MKL is not used; Eigen is the single numeric backend.

On Debian/Ubuntu:

```bash
sudo apt-get install g++ cmake pkg-config libeigen3-dev libfftw3-dev libgsl-dev libboost-dev
```

Dependencies resolve through `find_package` / `pkg-config` as imported targets,
so no hand-passed toolchain file is needed:

```bash
git clone https://github.com/adamecius/lsquant.git
cd lsquant
cmake -S . -B build
cmake --build build -j        # executables land in build/
ctest --test-dir build        # optional: 24 checks, all should pass
```

Build options:

- `-DLINQT_BUILD_INPUT_TOOLS=ON|OFF` (default **ON**) — fetch and build the
  pinned `wannier2sparse` used to generate the test fixtures; it clones at
  configure time and degrades gracefully offline.
- `-DLINQT_BUILD_UTILITIES=ON|OFF` (default **OFF**) — build the Python helper
  utilities.

## The `lsquant` interface

One binary drives the pipeline. A run is described by a small JSON file rather
than positional arguments:

```bash
# 1. compute Chebyshev moments  (mode: spectral | noneq | msd)
lsquant compute --config run.json

# 2. reconstruct a spectrum  (dos from a .chebmom1D; bastin/greenwood from a .chebmom2D)
lsquant reconstruct SpectralOp*.chebmom1D     dos       120      # broadening in meV
lsquant reconstruct NonEqOpVX-VY*.chebmom2D   bastin    10

# inspect a run or an operator and its resolved numerical provenance
lsquant inspect run.json
```

A density-of-states run is the spectral function of the identity operator:

```json
{ "mode": "spectral", "label": "chain1d_N256", "operator": "1",
  "num_moments": 512, "state": "exact_chain1d_N256" }
```

A Hall-conductivity run uses two velocity operators:

```json
{ "label": "haldane", "operator_right": "VX", "operator_left": "VY",
  "num_moments": 128, "state": "exact" }
```

`compute` writes a `.chebmom*` moment table; `reconstruct` turns it into a curve.
Every operator travels with a `.desc` sidecar (its identity, units, and spectral
bounds), so a result stays traceable to the operator that produced it, and
`inspect` reports that provenance. The legacy per-task drivers
(`inline_compute-kpm-*`, `*FromChebmom`) remain in `build/` and give
byte-identical results.

## Tutorials

Three self-contained physics stories under [`examples/`](examples/), each with a
script that reproduces its figures, all driven through the `lsquant` interface:

| | Tutorial | What it teaches |
|---|---|---|
| 1 | [Density of states of a 1D chain](examples/01_dos_1d_chain) | KPM moments, the resolution dial, the finite-size revival |
| 2 | [Localization on a disordered chain](examples/02_localization_1d_chain) | the localization length from mean-square displacement |
| 3 | [Quantized Hall plateau of the Haldane model](examples/03_haldane_hall) | Kubo–Bastin $\sigma_{xy}$ and the Fermi-sea topological plateau |

Start with Tutorial 1; each reuses the previous chain.

## How it works

KPM expands a spectral quantity in Chebyshev polynomials of the rescaled
Hamiltonian $\tilde H = (H - b)/a$, which maps the band onto $[-1, 1]$. The
coefficients are the moments $\mu_m = \mathrm{Tr}\,[\,O\,T_m(\tilde H)\,]$,
estimated by a stochastic (or exact) trace and a three-term recursion whose only
contact with $H$ is a sparse matrix–vector product. That recursion is the
expensive, linear-scaling step; reconstruction is a cheap sum over moments with a
Jackson kernel that sets the broadening $\eta = \pi W_{1/2}/M$, so the moment
count is the resolution dial.

Internally the recursion speaks a single operator interface (`HamiltonianOp`),
with the sparse matrix as the live implementation and the basis metric as a
separate policy, so the hot loop pays one virtual call per matrix–vector product
and the backend stays swappable.

## Inputs

LinQT consumes tight-binding models as sparse operators in **CSR text format**,
one file per operator (`operators/<label>.<OP>.CSR`), produced by external tools
such as **wannier2sparse**, **Wannier90**, or **KWANT**. Each operator carries a
`.desc` sidecar describing what it physically is (observable, units, spectral
bounds, provenance). The `examples/` generators (`make_chain.py`, …) show the
format end to end.

## Tests

The suite validates LinQT against a symbolically-proven analytic oracle for the
clean 1D chain, plus regression goldens for graphene, the Haldane model,
localization, and the FFT reconstruction path:

```bash
ctest --test-dir build --output-on-failure
```

Golden masters live in `test/golden/`; regenerate them deliberately, and only
when behavior changes on purpose, with `scripts/regenerate_goldens.sh` (seed
pinned). The specification is in [`test/TEST_PLAN.md`](test/TEST_PLAN.md).

## Documentation

- Physics notes under [`docs/`](docs/) (`main.md`, `time_evolution.md`,
  `spinrelaxation.md`).
- API reference via Doxygen: `cd docs && doxygen linqt.dox` writes `docs/html/`.

## Repository layout

| Path | Contents |
|---|---|
| `include/`, `src/` | KPM core library, the `lsquant` binary, and the legacy drivers |
| `examples/` | the three worked tutorials and their figure scripts |
| `test/` | analytic-oracle and regression tests, golden masters, `TEST_PLAN.md` |
| `docs/` | physics pages and the Doxygen config |
| `scripts/`, `cmake/` | golden regeneration; the `wannier2sparse` fetch module |

## License and citation

See [LICENSE](LICENSE). If LinQT contributes to your work, please cite the
methodology, reference [2] above (arXiv:1811.07387).
