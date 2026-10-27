# Getting Started

This tutorial shows how to build LSQUANT and compute the density of states of a simple 1D tight-binding chain.

## Prerequisites

You need:
- A C++17 compiler (g++ ≥ 11 or clang++ ≥ 14)
- CMake ≥ 3.16
- Eigen3, FFTW3, OpenMP
- Python ≥ 3.10 with numpy, scipy, sympy, and pybind11

On Ubuntu/Debian:

```bash
sudo apt-get install g++ cmake pkg-config libeigen3-dev libfftw3-dev libgsl-dev python3-dev
pip install numpy scipy sympy pybind11 matplotlib
```

## Build

```bash
git clone https://github.com/adamecius/lsquant.git
cd lsquant
cmake --preset default
cmake --build --preset default
export PYTHONPATH="$PWD/build/default/lsquant_core:$PWD/python"
```

## Run a density of states calculation

```bash
python3 -m lsquant run examples/01_dos_1d_chain/run_dos.json --out /tmp/chain_dos
python3 -m lsquant inspect /tmp/chain_dos
```

You should obtain the density of states of a 1D tight-binding chain, showing the expected van Hove singularities.

## Next steps

- Try the Haldane model Hall conductivity example: `examples/03_haldane_hall/`
- See how to write your own calculation: [JSON run file guide](../how-to/run-json-file.md) (coming soon)
- Load a SIESTA Hamiltonian for spin Hall calculations
