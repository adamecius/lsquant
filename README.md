# LSQUANT — Linear-Scaling Quantum Transport

**LSQUANT** is a C++17 and Python/SymPy library for **large-scale spectral and
quantum transport calculations** in tight-binding and first-principles local-orbital
models using the **Kernel Polynomial Method (KPM)** [1, 2]. It evaluates densities
of states, longitudinal and Hall conductivities (Kubo–Greenwood, Kubo–Bastin, and
intrinsic/extrinsic spin Hall), real-time wave-packet and spin dynamics, and
mean-square displacement at a cost **linear in system size** without diagonalising
the Hamiltonian. One observable word grammar runs across real-space sparse
supercells, Bloch $k$-space Hamiltonians with FFT-sandwiched real-space disorder,
non-orthogonal SIESTA and Wannier bases (Löwdin and generalised metric routes),
spectral window projections, and an exact per-$k$ eigenbasis oracle.

## The unifying idea

Every spectral or transport observable in LSQUANT is written as the trace of an
ordered word of static operators $A_i$ and spectral functions $f_i(H)$:

$$
C = \mathrm{Tr}\!\left[\,A_0\,f_1(H)\,A_1\,f_2(H)\,A_2\,\cdots\,f_p(H)\,A_p\,\right]
$$

1. **Python/SymPy parses the word** and classifies each slot as a static operator
   (`1`, `VX`, `VY`, `VZ`, `SX`, `SY`, `SZ`, `JX_SZ`, projectors `(1+SX)/2`, or
   commutators `[X, U(t)]`) or a spectral function of $H$ (`delta(E-H)`, `G+(E)`,
   `G-(E)`, `dG+(E)/dE`, `dG-(E)/dE`, `U(t)`).
2. **The compiler lowers the word into a five-instruction bytecode `Plan`**:
   `ApplyOperator`, `Expand` (opens a Chebyshev moment axis via the three-term
   recursion), `ApplySeries` (contracts a Clenshaw Chebyshev series on the fly,
   such as the Bessel series for $U(t)$), `MirrorToBra`, and `Trace`.
3. **The C++ `MomentEngine` executes the `Plan`** on any `LinearOperator` and
   `StateGenerator`, returning the rank-$r$ `ChebyshevMoments` tensor
   $\mu_{m_1 \dots m_r}$. The engine has no notion of basis or observable names.
4. **Python contracts $\mu$ with damped Chebyshev coefficient tables** (GEMV for
   rank 1, GEMM or 2D FFT for rank 2, symbolic Fermi-sea integrals for
   Kubo–Bastin and spin Hall) to reconstruct the physical curve.

## Quickstart

```bash
# 1. Configure and build the host-tuned C++ core and Python extension
cmake --preset native
cmake --build --preset native
export PYTHONPATH="$PWD/build/native/lsquant_core:$PWD/python"

# 2. Run Tutorial 03: quantized Hall plateau of the Haldane model
python3 -m lsquant run examples/03_haldane_hall/run_haldane.json --out /tmp/haldane

# 3. Inspect the output and generate the publication figure
python3 -m lsquant inspect /tmp/haldane
python3 examples/03_haldane_hall/plot.py /tmp/haldane
```

## Installation

### Prerequisites

- **C++17 compiler** (`g++ >= 11` or `clang++ >= 14`) and **CMake** (`>= 3.21` for
  `CMakePresets.json`, or `>= 3.16` when invoking `cmake -S . -B build` directly)
- **Eigen3** (`>= 3.3`), **FFTW3** (double precision, plus `libfftw3_omp` when
  available), **OpenMP**, and **GSL**
- **Python** (`>= 3.10`) with `numpy`, `scipy`, `sympy`, `pybind11`, and
  `matplotlib` (plus `sisl` for SIESTA `.fdf`/`.TSHS`/`.HSX` models and `mcp` for
  the optional MCP server)

On Debian/Ubuntu:

```bash
sudo apt-get install g++ cmake pkg-config libeigen3-dev libfftw3-dev libgsl-dev python3-dev
pip install numpy scipy sympy pybind11 matplotlib sisl
```

### Build presets and install

[`CMakePresets.json`](CMakePresets.json) provides two build profiles:

- **`native`** (`cmake --preset native && cmake --build --preset native`): enables
  `-O3 -march=native -mtune=native -fno-math-errno` and `LSQUANT_NATIVE_BUILD=1`
  (SIMD-aligned scratch buffers, NUMA-aware first-touch allocation, and
  cache-blocked `BlockDiagonalOperator` evaluation). Recommended for all runs.
- **`golden`** (`cmake --preset golden && cmake --build --preset golden`): portable
  IEEE-754 build without machine-specific flags, used for byte-for-byte golden
  verification.

To install a self-contained tree (`lib/`, `include/lsquant/`, `bin/`, and
`python/lsquant/`):

```bash
cmake --install build/native --prefix ~/.local/lsquant
export PYTHONPATH="$HOME/.local/lsquant/python"
python3 -c "import lsquant, lsquant._core"
```

## What LSQUANT computes

### Observable word catalogue

| Shape ID | Word syntax | Rank | Physical observable |
|---|---|---:|---|
| `dos` | `Tr[delta(E-H)]` | 1 | Density of states $\rho(E)$ |
| `spectral_1d` | `Tr[A delta(E-H)]` | 1 | Operator spectral density $A(E)$ |
| `local_dos` | `<i\| delta(E-H) \|i>` | 0 / 1 | Local density of states on site/orbital $i$ |
| `kubo_greenwood` | `Tr[VX delta(E-H) VX delta(E-H)]` | 2 | Longitudinal Kubo–Greenwood conductivity $\sigma_{\alpha\alpha}(E)$ |
| `kubo_bastin` | `Tr[VX delta(E-H) VY dG-(E)/dE]` | 2 | Kubo–Bastin Hall $\sigma_{\alpha\beta}(E)$ (total, Fermi-sea, Fermi-surface) and spin Hall $\sigma^{\gamma}_{\alpha\beta}(E)$ (`JX_SZ` velocity) |
| `time_correlation` | `Tr[A delta(E-H) U(t)+ B U(t)]` | 1 | Energy-resolved real-time correlator $C_{AB}(E, t)$ |
| `msd` | `Tr[[X, U(t)]+ delta(E-H) [X, U(t)]]` | 1 | Mean-square displacement $\Delta X^2(E, t)$ and diffusivity $D(E, t)$ |
| `projected_evolution` | `Tr[P U(t)+ S_a U(t) delta(E-H) P]` | 1 | Spin precession and relaxation $\langle S_\alpha(E, t) \rangle$ from initial state projector $P$ |

### Representations, metrics, and spectral projection

| Capability | Options | Description |
|---|---|---|
| **Bases** | `real`, `bloch`, `eigen` | Real-space CSR supercells (`RealBasis`); Bloch $k$-space `BlockDiagonalOperator` with optional FFT-sandwiched real-space disorder $H_k + B^\dagger V B$ (`BlochBasis`); exact per-$k$ eigenbasis oracle (`EigenBasis`). |
| **Non-orthogonal metrics** | `orthogonal`, `nonorthogonal` (`lowdin` or `generalised`) | Overlap matrix $S$ handled via `BlochMetricBasis` / `EigenMetricBasis` or real-space appliers (`exact` block LLT/eigh, `chebyshev` polynomial expansion of $S^{-1}$ and $S^{-1/2}$, or `ischo` incomplete-Cholesky PCG), with `covariant` (dipole-corrected), `reduced`, or `truncated` velocity gauges. |
| **Spectral projection** | `shift`, `truncate`, `orthogonal` | Energy-window projection ($P_W H P_W + \sigma Q_W$, hard rank truncation, or isolated-manifold SCDM polar downfolding to an orthogonal effective model). |

## Interfaces

### 1. Command-line interface (`python -m lsquant`)

```bash
# List all catalogued observable word shapes
python3 -m lsquant words

# Parse a word and inspect its AST, bytecode Plan, and allowed reconstruction formulas
python3 -m lsquant words describe "Tr[JX_SZ delta(E-H) VY dG-(E)/dE]"

# Validate and execute a JSON run file
python3 -m lsquant run examples/01_dos_1d_chain/run_dos.json --out /tmp/chain_dos

# Reconstruct a curve from saved moments with a different broadening or formula
python3 -m lsquant reconstruct /tmp/chain_dos --formula dos --broadening-meV 15

# Inspect a run file, moment file (.chebmom*), or output directory
python3 -m lsquant inspect /tmp/chain_dos
```

### 2. Python API

```python
import numpy as np
import lsquant as lsq

basis = lsq.RealBasis("examples/03_haldane_hall/operators", "haldane",
                      operators=["VX", "VY"])
compiled = lsq.words.compile("Tr[VX delta(E-H) VY dG-(E)/dE]",
                             num_moments=128, known_operators=["VX", "VY"])
moments = lsq.run_plan(compiled.plan, basis.hamiltonian(),
                       [basis.operator(op) for op in compiled.operator_names],
                       basis.states("exact"), basis.window())
energies, sigma_xy = lsq.kubo_bastin(moments, np.linspace(-3.0, 3.0, 601), eta=0.1)
```

### 3. JSON run file

Every calculation is reproducible from a single strictly validated JSON run file
(unknown keys are rejected before any operator is loaded). Example
(`examples/03_haldane_hall/run_haldane.json`):

```json
{
  "model": {
    "kind": "csr",
    "dir": "examples/03_haldane_hall/operators",
    "label": "haldane",
    "operators": ["VX", "VY"]
  },
  "basis": { "name": "real" },
  "word": "Tr[VX delta(E-H) VY dG-(E)/dE]",
  "expansion": { "moments": 128, "window": "gershgorin", "trace": "exact" },
  "reconstruction": {
    "kernel": "jackson",
    "broadening_meV": 100.0,
    "formulas": ["bastin", "sea", "surf", "greenwood"]
  }
}
```

See [`python/lsquant/README.md`](python/lsquant/README.md) for the complete
run-file schema and Python API reference.

### 4. MCP server (`python -m lsquant.mcp_server`)

LSQUANT includes a Model Context Protocol (MCP) server over `stdio` exposing
`words_describe`, `run`, `status`, `reconstruct`, and `inspect` for AI coding and
research agents:

```bash
PYTHONPATH="build/native/lsquant_core:python" python3 -m lsquant.mcp_server --root "$PWD"
```

## Tutorials and examples

Eight self-contained examples under [`examples/`](examples/) cover 1D chains, 2D
topological insulators, diffusive/localised lattices, spin relaxation, and
first-principles SIESTA spin-Hall calculations. Every example uses the shared
publication figure stylesheet [`examples/prl_style.py`](examples/prl_style.py).

| # | Example | Physics & observable | Basis & metric | Word | Typical runtime |
|---:|---|---|---|---|---|
| 01 | [`01_dos_1d_chain`](examples/01_dos_1d_chain) | 1D chain DOS, Jackson broadening dial, finite-size revival | `real` (orthogonal) | `Tr[delta(E-H)]` | ~2 s |
| 02 | [`02_localization_1d_chain`](examples/02_localization_1d_chain) | Anderson localisation length $\xi(E)$ from saturated MSD | `real` (orthogonal) | `Tr[[X, U(t)]+ delta(E-H) [X, U(t)]]` | ~15 s |
| 03 | [`03_haldane_hall`](examples/03_haldane_hall) | Quantised Chern plateau $\sigma_{xy} = e^2/h$ and Sea/Surf split | `real` (orthogonal) | `Tr[VX delta(E-H) VY dG-(E)/dE]` | ~3 s |
| 04 | [`04_spin_precession_magnetic_chain`](examples/04_spin_precession_magnetic_chain) | Larmor spin precession in an exchange field | `real` (orthogonal) | `Tr[(1+SX)/2 U(t)+ S_a U(t) delta(E-H) (1+SX)/2]` | ~1 s |
| 05 | [`05_charge_transport_square`](examples/05_charge_transport_square) | Ballistic $\to$ diffusive $\to$ localised transport on a 2D square lattice | `real` (orthogonal) | `Tr[delta(E-H)]` + `Tr[[X, U(t)]+ delta(E-H) [X, U(t)]]` | ~10 s |
| 06 | [`06_spin_relaxation_square`](examples/06_spin_relaxation_square) | Elliott–Yafet vs. Dyakonov–Perel spin relaxation in 2D | `real` (orthogonal) | `Tr[(1+S_a)/2 U(t)+ S_a U(t) delta(E-H) (1+S_a)/2]` | ~25 s |
| 07 | [`07_spin_hall_bi_bilayer`](examples/07_spin_hall_bi_bilayer) | Intrinsic spin Hall conductivity $\sigma^z_{xy}(E)$ of 2D $\mathrm{Bi}(111)$ bilayer from SIESTA SOC | `bloch` + `eigen` (non-orthogonal, Löwdin + generalised, covariant dipoles) | `Tr[JX_SZ delta(E-H) VY dG-(E)/dE]` | ~8 s (`quick`) / ~2 min (`production`) |
| 08 | [`08_spin_hall_pdse2_monolayer`](examples/08_spin_hall_pdse2_monolayer) | Intrinsic spin Hall conductivity $\sigma^z_{xy}(E)$ of pentagonal $\mathrm{PdSe_2}$ monolayer from SIESTA SOC | `bloch` + `eigen` (non-orthogonal, Löwdin + generalised, covariant dipoles) | `Tr[JX_SZ delta(E-H) VY dG-(E)/dE]` | ~15 s (`quick`) / ~4 min (`production`) |

## Project layout

| Path | Contents |
|---|---|
| [`lsquant/core/algebra/`](lsquant/core/algebra) | `LinearOperator`, `SparseMatrix`, `BlockDiagonalOperator`, `BlochTransform` (batched FFTW), `MetricPolicy`, `OverlapAppliers` (`exact`, `chebyshev`, `ischo`), `StateGenerator`, `SpectralBounds`, `FftMomentContraction` |
| [`lsquant/core/chebyshev/`](lsquant/core/chebyshev) | `SpectralWindow`, `RescaledOperator`, `ChebyshevRecursion`, `ChebyshevVectors`, `ChebyshevMoments`, `Plan`, `MomentEngine`, `Contraction` |
| [`lsquant/io/`](lsquant/io) | CSR text reader, `.chebmom` and legacy `.chebmom1D/2D/TD` readers/writers, curve text writer |
| [`lsquant/api/`](lsquant/api), [`lsquant/compat/`](lsquant/compat) | `pybind11` module `lsquant._core` and legacy GSL Bessel / text round-trip helpers |
| [`python/lsquant/`](python/lsquant) | Python library: `words/` parser & compiler, SymPy `spectral_functions`, `kernels`, `reconstruct`, `bases`, `disorder_bases`, `metric_forms`, `models/` (`wannier`, `siesta`, `nonorthogonal`, `disorder`), `projection/`, `oracle/`, `runfile`, `cli`, `mcp_server` |
| [`examples/`](examples) | Worked tutorials `01`–`06` and first-principles SIESTA spin-Hall examples `07`–`08` |
| [`runs/goldens/`](runs/goldens), [`goldens/`](goldens) | Committed golden-master inputs, reference outputs, and JSON run files (`G1`–`G10`, tutorials `01`–`06`, SIESTA/disorder/projection gates) |
| [`figures/`](figures) | Stage verification figure scripts (`stage_1.py` … `stage_9_local_projection.py`) |
| [`docs/`](docs), [`bench/`](bench) | Physics notes, executable SymPy derivations ([`docs/derivations/`](docs/derivations)), performance study ([`docs/performance/`](docs/performance/README.md)), and kernel benchmarks (`bench/kernel_bench.cpp`, `bench/block_domains/`) |

## Verification and performance

Every stage of LSQUANT is gated against committed golden masters and analytic or
exact-diagonalisation oracles by [`scripts/check_goldens.sh`](scripts/check_goldens.sh),
and developer self-checks live under [`tests/dev/`](tests/dev/). OpenMP reductions
use a fixed per-state accumulation order so that multi-threaded runs are
**bit-for-bit identical** between `OMP_NUM_THREADS=1` and `OMP_NUM_THREADS=16`.
In the `native` build profile, SIMD-aligned scratch buffers, NUMA first-touch
allocation, cache-blocked `BlockDiagonalOperator` evaluation, and 2D FFT Kubo–Bastin
contractions accelerate large k-space and real-space runs; see
[`docs/performance/README.md`](docs/performance/README.md) for the full benchmark
report and roofline analysis.

## References, citation, and licence

1. A. Weiße, G. Wellein, A. Alvermann, and H. Fehske, *The kernel polynomial
   method*, Rev. Mod. Phys. **78**, 275 (2006).
2. Z. Fan, J. H. García, A. W. Cummings, J. E. Barrios-Vargas, M. Panhans,
   A. Harju, F. Ortmann, and S. Roche, *Linear scaling quantum transport
   methodologies*, Phys. Rep. **893**, 1–123 (2021) ([arXiv:1811.07387](https://arxiv.org/abs/1811.07387)).

See [`LICENSE`](LICENSE) for licence terms. If LSQUANT contributes to your work,
please cite Reference [2] above.
