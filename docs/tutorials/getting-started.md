# Getting Started

This tutorial shows how to build LSQUANT and compute the density of states of a simple 1D chain.

## Prerequisites

You need a C++17 compiler, CMake ≥ 3.16, Eigen3, FFTW3, OpenMP, and Python ≥ 3.10 with numpy, scipy, sympy, and pybind11.

On Ubuntu/Debian:

```bash
sudo apt-get install g++ cmake pkg-config libeigen3-dev libfftw3-dev libgsl-dev python3-dev
pip install numpy scipy sympy pybind11 matplotlib
```

## Build

```bash
git clone https://github.com/adamecius/lsquant.git
cd lsquant
cmake --preset golden
cmake --build --preset golden
export PYTHONPATH="$PWD/build/golden/lsquant_core:$PWD/python"
```

## Run a DOS calculation

```bash
python3 -m lsquant run examples/01_dos_1d_chain/run_dos.json --out /tmp/chain_dos
python3 -m lsquant inspect /tmp/chain_dos
```

You should see a density of states with the expected van Hove singularities of a 1D tight-binding chain.

Next: see the [Haldane Hall tutorial](haldane-hall.md) or the [JSON run file how-to](../how-to/run-json-file.md).
