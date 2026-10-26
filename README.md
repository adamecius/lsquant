# LSQUANT — Linear-Scaling Quantum Transport

**LSQUANT** is a C++17 / Python library for large-scale spectral and quantum transport calculations using the Kernel Polynomial Method (KPM). It computes densities of states, longitudinal and Hall conductivities, spin Hall conductivity, real-time dynamics, and mean-square displacement at linear cost in system size, without diagonalizing the Hamiltonian.

It supports real-space sparse models, Bloch $k$-space Hamiltonians, non-orthogonal bases (including SIESTA), and an exact eigenbasis oracle.

## Quick start

```bash
cmake --preset golden
cmake --build --preset golden
export PYTHONPATH="$PWD/build/golden/lsquant_core:$PWD/python"

python3 -m lsquant run examples/01_dos_1d_chain/run_dos.json --out /tmp/chain_dos
python3 -m lsquant inspect /tmp/chain_dos
```

## Documentation

The full documentation is organized following the Diátaxis framework:

- **[Getting Started](docs/tutorials/getting-started.md)** — build and run your first calculation
- **[Unifying Idea](docs/explanation/unifying-idea.md)** — the observable-word abstraction
- **[Observable Word Catalogue](docs/reference/observable-words.md)** — supported quantities
- **[Documentation index](docs/index.md)** — complete structure

## Examples

Eight self-contained examples are available under [`examples/`](examples/), covering DOS, localization, Haldane Hall, spin dynamics, and first-principles SIESTA spin Hall.

## Installation

See the [Getting Started](docs/tutorials/getting-started.md) tutorial for prerequisites and build instructions.

## License

See [LICENSE](LICENSE).
