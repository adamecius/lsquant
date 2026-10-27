# LSQUANT — Linear-Scaling Quantum Transport

**LSQUANT** is a library for large-scale spectral and quantum transport calculations in tight-binding and first-principles models. It uses the Kernel Polynomial Method (KPM) to compute densities of states, longitudinal and Hall conductivities, spin Hall conductivity, real-time dynamics, and mean-square displacement at a cost that grows linearly with system size. It does not diagonalize the Hamiltonian.

It works with real-space sparse models, Bloch $k$-space Hamiltonians, non-orthogonal bases (including SIESTA), and can also perform exact reference calculations on small systems for validation.

## Quick start

```bash
cmake --preset default
cmake --build --preset default
export PYTHONPATH="$PWD/build/default/lsquant_core:$PWD/python"

python3 -m lsquant run examples/01_dos_1d_chain/run_dos.json --out /tmp/chain_dos
python3 -m lsquant inspect /tmp/chain_dos
```

This computes the density of states of a simple 1D chain. You should see the characteristic van Hove singularities.

## What you can calculate

- Density of states (total and projected)
- Longitudinal conductivity (Kubo–Greenwood)
- Hall and spin Hall conductivity (Kubo–Bastin)
- Real-time spin dynamics and mean-square displacement
- Calculations on SIESTA and Wannier models

See the full list in the [Observable Catalogue](docs/reference/observable-words.md).

## Documentation

- **[Getting Started](docs/tutorials/getting-started.md)** — install and run your first calculation
- **[How it works](docs/explanation/unifying-idea.md)** — the core idea behind the calculations
- **[Observable Catalogue](docs/reference/observable-words.md)** — supported physical quantities
- **[Documentation index](docs/index.md)** — complete guide

## Examples

Ready-to-run examples are in [`examples/`](examples/), including Haldane model Hall conductivity, Anderson localization, spin relaxation, and first-principles SIESTA spin Hall calculations.

## Installation

See the [Getting Started](docs/tutorials/getting-started.md) page for dependencies and build instructions.

## License

See [LICENSE](LICENSE).
