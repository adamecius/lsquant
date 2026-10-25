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
   (`1`, `VX`, `VY`, `VZ`, `SX`, `SY`, `SZ`, `JXSZ`, projectors `(1+SX)/2` or `P:PX`,
   or commutators `[X, U(t)]`) or a spectral function of $H$ (`delta(E-H)`,
   `dG-(E)/dE`, `U(t)`).
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

## What is new

With respect to the previous LSQUANT release, this unified architecture provides:

- **Unified five-instruction C++ engine (`MomentEngine`):** replaces the separate
  1D, 2D, time-dependent, and local C++ solver classes and the 41 `inline_*`
  driver executables with a single bytecode interpreter (`ApplyOperator`,
  `Expand`, `ApplySeries`, `MirrorToBra`, `Trace`).
- **Observable word language and compiler:** observables are specified as
  algebraic trace words `Tr[...]` validated against a formal grammar and
  catalogue (`python -m lsquant words`), compiled to `Plan` objects, and driven
  from a single schema-validated JSON run file (`python -m lsquant run`), CLI,
  or MCP server (`python -m lsquant.mcp_server`).
- **Three representations (`real`, `bloch`, `eigen`):** real-space CSR
  supercells (`RealBasis`), $k$-space `BlockDiagonalOperator` Hamiltonians with
  real-space Anderson disorder applied through batched FFTW transforms
  $H_k + B^\dagger V B$ (`BlochBasis`), and an exact per-$k$ diagonalisation
  path (`EigenBasis`) serving as the internal oracle.
- **Non-orthogonal metric policies (`lowdin` and `generalised`):** overlap
  matrices $S$ are handled by composite `LinearOperator` adapters (`exact` dense
  block LLT/eigh, `chebyshev` polynomial expansions of $S^{-1}$ and $S^{-1/2}$,
  or `ischo` incomplete-Cholesky PCG) without modifying the Chebyshev recursion.
- **First-principles SIESTA spin Hall (`sisl` loader) and spectral projection:**
  direct loading of spin–orbit coupled SIESTA `.fdf` + `.TSHS`/`.HSX` and
  `.ion.xml` orbital tables with `bare`, `berry_connection`, and `covariant`
  (certified radial dipole quadrature) velocity gauges, plus spectral window
  projection (`hard`, `shift`, and `orthogonal` SCDM downfolding).
- **Build profiles (`golden` and `native`) and optional FFT contraction:**
  `golden` provides portable IEEE-754 arithmetic for byte-identical golden
  verification; `native` is an opt-in performance profile under validation
  (`-O3 -march=native -mtune=native -fno-math-errno`, NUMA first-touch placement,
  and tiled GEMM moment tables); `FftMomentContraction` provides an opt-in
  FFTW-based rank-2 Chebyshev-node contraction (`LSQUANT_FFT_BACKEND=fft`).
- **Corrected legacy defects:**
  1. *Batched rank-2 moment shift:* in the legacy code, running rank-2 moments
     with `BATCH_SIZE < M` did not restart the Chebyshev recursion between
     sections, shifting off-diagonal blocks by $k$ columns; the new
     `SlicedChebyshevVectors` engine gives byte-identical moments for any batch
     size.
  2. *Off-centre uniform-grid energy mapping:* legacy uniform-grid Kubo
     reconstruction mapped energies as $x\,\mathrm{HalfWidth} +
     \mathrm{BandCenter}/\mathrm{HalfWidth}$ instead of $x\,\mathrm{HalfWidth} +
     \mathrm{BandCenter} = (x - \mathrm{shift})/\mathrm{scale}$; fixed in
     `lsquant.reconstruct`.
  3. *Explicit identity operator in Tutorial 06:* the legacy CSR loader returned
     an empty matrix silently when `<label>.1.CSR` was absent, producing noisy
     denominators in Tutorial 06; the new loader refuses missing operator files
     and maps `1` to `IdentityOperator`.

## Quickstart

```bash
# 1. Configure and build the C++ core and Python extension (golden preset)
cmake --preset golden
cmake --build --preset golden
export PYTHONPATH="$PWD/build/golden/lsquant_core:$PWD/python"

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
- **Eigen3** (`>= 3.3`), **FFTW3** (double precision), **OpenMP**, and **GSL**
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

- **`golden`** (`cmake --preset golden && cmake --build --preset golden`, default):
  portable IEEE-754 build without machine-specific flags, validated byte for byte
  against the committed golden masters.
- **`native`** (`cmake --preset native && cmake --build --preset native`): opt-in
  performance profile under validation (0009 wrap-up item W08), enabling
  `-O3 -march=native -mtune=native -fno-math-errno` and `LSQUANT_NATIVE_BUILD=1`
  (NUMA-aware first-touch allocation, batched operator dispatch, and tiled
  `MomentTable` / contraction kernels).

To install a self-contained tree (`lib/`, `include/lsquant/`, `bin/`, and
`python/lsquant/`):

```bash
cmake --install build/golden --prefix ~/.local/lsquant
export PYTHONPATH="$HOME/.local/lsquant/python"
python3 -c "import lsquant, lsquant._core"
```

## What LSQUANT computes

### Observable word catalogue

| Shape ID | Canonical word | Rank | Plan builder & default reconstruction |
|---|---|---:|---|
| `dos_identity` | `Tr[delta(E-H)]` | 1 | `dos_plan` $\to$ `dos_curve` ($\rho(E)$) |
| `dos_operator` | `Tr[A delta(E-H)]` | 1 | `dos_plan` $\to$ `dos_curve` ($A(E)$) |
| `kubo_greenwood` | `Tr[VX delta(E-H) VX delta(E-H)]` | 2 | `kubo_plan` $\to$ `kubo_greenwood` ($\sigma_{\alpha\alpha}(E)$) |
| `kubo_bastin` | `Tr[VX delta(E-H) VY dG-(E)/dE]` | 2 | `kubo_plan` $\to$ `kubo_bastin` ($\sigma_{\alpha\beta}(E)$, Sea/Surf, spin Hall `shc` with `JXSZ`) |
| `msd` | `Tr[[X,U(t)]^+ delta(E-H) [X,U(t)]]` | 1 | `msd_plan` $\to$ `time_correlation_curve` ($\Delta X^2(E, t)$) |
| `projected_evolution` | `Tr[P U(t)^+ SX U(t) delta(E-H) P]` | 1 | `projected_evolution_plan` $\to$ `row_m0` / `time_correlation_curve` |
| `time_correlation` | `Tr[A delta(E-H) U(t)^+ B U(t)]` | 1 | `time_correlation_plan` $\to$ `time_correlation_curve` ($C_{AB}(E, t)$) |

### Representations, metrics, and spectral projection

| Capability | Options | Description |
|---|---|---|
| **Bases** (`basis.name`) | `real`, `bloch`, `eigen` | Real-space CSR supercells (`RealBasis`); Bloch $k$-space `BlockDiagonalOperator` with optional FFT-sandwiched real-space disorder $H_k + B^\dagger V B$ (`BlochBasis`); exact per-$k$ eigenbasis oracle (`EigenBasis`). |
| **Non-orthogonal metrics** (`metric.policy`, `metric.route`, `metric.applier`) | `orthogonal`, `nonorthogonal` (`lowdin` or `generalised`); appliers `exact`, `chebyshev`, `ischo` | Overlap matrix $S$ handled via `BlochMetricBasis` / `EigenMetricBasis` or real-space appliers (`exact` dense LLT/eigh, `chebyshev` polynomial expansion of $S^{-1}$ and $S^{-1/2}$, or `ischo` incomplete-Cholesky PCG). |
| **SIESTA velocity modes** (`model.params.velocity_mode`) | `bare`, `berry_connection`, `covariant` | Lattice derivative $\partial_k H$ (`bare`), intracell atom-phase + overlap Pulay connection (`berry_connection`), or full two-centre orbital dipole connection integrated from `.ion.xml` (`covariant`). |
| **Spectral projection** (`projection.kind`) | `hard`, `shift`, `orthogonal` | Energy-window projection ($P_W H P_W$ `hard` rank truncation, $P_W H P_W + \sigma Q_W$ `shift`, or isolated-manifold SCDM polar downfolding to an `orthogonal` effective model). |

## Interfaces

### 1. Command-line interface (`python -m lsquant`)

Subcommands and flags (matching `python3 -m lsquant <cmd> --help`):

- `python3 -m lsquant run [-h] [--out OUT] [--json] runfile`
- `python3 -m lsquant reconstruct [-h] [--formula FORMULA] [--broadening-meV BROADENING_MEV] [--fermi-energies FERMI_ENERGIES [FERMI_ENERGIES ...]] [--out OUT] [--json] target [pos_formula] [pos_broadening]`
- `python3 -m lsquant inspect [-h] path`
- `python3 -m lsquant words [-h] [action] [word]`

```bash
# List all catalogued observable word shapes
python3 -m lsquant words

# Parse a word and inspect its AST, bytecode Plan, and allowed reconstruction formulas
python3 -m lsquant words describe "Tr[JXSZ delta(E-H) VY dG-(E)/dE]"

# Validate and execute a JSON run file
python3 -m lsquant run examples/01_dos_1d_chain/run_dos.json --out /tmp/chain_dos

# Reconstruct a curve from saved moments with an explicit formula and broadening (in meV)
python3 -m lsquant reconstruct /tmp/chain_dos --formula dos --broadening-meV 15

# Inspect a run file, moment file (.chebmom*), or output directory
python3 -m lsquant inspect /tmp/chain_dos
```

### 2. Python API

High-level execution via `lsquant.execute_runfile(dict)` and low-level execution
via `_core.ReadCsrText` + `kubo_plan` + `run_plan` + `kubo_bastin` (as in
[`python/lsquant/apps/stage3_kubo.py`](python/lsquant/apps/stage3_kubo.py)):

```python
import lsquant as lsq

# 1. High-level: validate and execute a run-file dictionary
res = lsq.execute_runfile(
    {
        "model": {
            "kind": "csr",
            "dir": "goldens/stage0/haldane/operators",
            "label": "haldane",
            "operators": ["VX", "VY"],
        },
        "basis": {"name": "real"},
        "word": "Tr[VX delta(E-H) VY dG-(E)/dE]",
        "expansion": {"moments": 128, "window": [-4.5, 4.5], "trace": "exact"},
        "reconstruction": {
            "kernel": "jackson",
            "broadening_meV": 10.0,
            "formulas": ["bastin"],
        },
        "output": {"dir": "/tmp/haldane_api"},
    }
)
print("High-level curve_files:", res["curve_files"])

# 2. Low-level: ReadCsrText + kubo_plan + run_plan + kubo_bastin
ham = lsq._core.ReadCsrText("goldens/stage0/haldane/operators/haldane.HAM.CSR")
vx = lsq._core.SparseMatrixOperator(
    lsq._core.ReadCsrText("goldens/stage0/haldane/operators/haldane.VX.CSR")
)
vy = lsq._core.SparseMatrixOperator(
    lsq._core.ReadCsrText("goldens/stage0/haldane/operators/haldane.VY.CSR")
)
window = lsq._core.SpectralWindow(-4.5, 4.5)
states = lsq._core.LocalStates(ham.dimension, list(range(ham.dimension)))
plan = lsq.kubo_plan(128, left_index=0, right_index=1, word="VX delta VY delta")
moments = lsq.run_plan(
    plan, ham, window, states, operators=[vx, vy], batch_size=128, symmetrise=True
)
energies, sigma_xy = lsq.kubo_bastin(moments, 10.0)
print("Low-level energies.shape, sigma_xy.shape:", energies.shape, sigma_xy.shape)
```

### 3. JSON run file

Every calculation is reproducible from a single strictly validated JSON run file
(unknown keys are rejected before any operator is loaded). Example
([`examples/03_haldane_hall/run_haldane.json`](examples/03_haldane_hall/run_haldane.json)):

```json
{
  "model": {
    "kind": "csr",
    "dir": "goldens/stage0/haldane/operators",
    "label": "haldane",
    "operators": ["VX", "VY"]
  },
  "basis": { "name": "real" },
  "metric": { "policy": "orthogonal" },
  "word": "Tr[VX delta(E-H) VY dG-(E)/dE]",
  "words": [
    "Tr[VX delta(E-H) VY dG-(E)/dE]",
    "Tr[VX delta(E-H) VX delta(E-H)]"
  ],
  "expansion": {
    "moments": 128,
    "window": [-4.5, 4.5],
    "trace": "exact",
    "states": 512
  },
  "reconstruction": {
    "kernel": "jackson",
    "broadening_meV": 10.0,
    "formulas": ["bastin", "sea", "surf", "greenwood"],
    "legacy_text_round_trip": true
  },
  "output": {
    "dir": ".",
    "legacy_names": true,
    "legacy_layout": true,
    "capture_tutorial": false
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
PYTHONPATH="build/golden/lsquant_core:python" python3 -m lsquant.mcp_server --root "$PWD"
```

## Tutorials and examples

Eight self-contained examples under [`examples/`](examples/) cover 1D chains, 2D
topological insulators, diffusive/localised lattices, spin relaxation, and
first-principles SIESTA spin-Hall calculations. Every example uses the shared
publication figure stylesheet [`examples/prl_style.py`](examples/prl_style.py).

| # | Example | Physics & observable | Basis & metric | Word |
|---:|---|---|---|---|
| 01 | [`01_dos_1d_chain`](examples/01_dos_1d_chain) | 1D chain DOS, Jackson broadening dial, finite-size revival | `real` (orthogonal) | `Tr[delta(E-H)]` |
| 02 | [`02_localization_1d_chain`](examples/02_localization_1d_chain) | Anderson localisation length $\xi(E)$ from saturated MSD | `real` (orthogonal) | `Tr[[X,U(t)]^+ delta(E-H) [X,U(t)]]` |
| 03 | [`03_haldane_hall`](examples/03_haldane_hall) | Quantised Chern plateau $\sigma_{xy} = e^2/h$ and Sea/Surf split | `real` (orthogonal) | `Tr[VX delta(E-H) VY dG-(E)/dE]` |
| 04 | [`04_spin_precession_magnetic_chain`](examples/04_spin_precession_magnetic_chain) | Larmor spin precession in an exchange field | `real` (orthogonal) | `Tr[P:PX U(t)^+ SX U(t) delta(E-H) P:PX]` |
| 05 | [`05_charge_transport_square`](examples/05_charge_transport_square) | Ballistic $\to$ diffusive $\to$ localised transport on a 2D square lattice | `real` (orthogonal) | `Tr[delta(E-H)]` + `Tr[[X,U(t)]^+ delta(E-H) [X,U(t)]]` |
| 06 | [`06_spin_relaxation_square`](examples/06_spin_relaxation_square) | Elliott–Yafet vs. Dyakonov–Perel spin relaxation in 2D | `real` (orthogonal) | `Tr[P:PZ U(t)^+ SZ U(t) delta(E-H) P:PZ]` |
| 07 | [`07_spin_hall_bi_bilayer`](examples/07_spin_hall_bi_bilayer) | Intrinsic spin Hall conductivity $\sigma^z_{xy}(E)$ of 2D $\mathrm{Bi}(111)$ bilayer from SIESTA SOC | `eigen` (non-orthogonal, Löwdin, `covariant` dipoles) | `Tr[JXSZ delta(E-H) VY dG-(E)/dE]` |
| 08 | [`08_spin_hall_pdse2_monolayer`](examples/08_spin_hall_pdse2_monolayer) | Intrinsic spin Hall conductivity $\sigma^z_{xy}(E)$ of pentagonal $\mathrm{PdSe_2}$ monolayer from SIESTA SOC | `eigen` (non-orthogonal, Löwdin, `covariant` dipoles) | `Tr[JXSZ delta(E-H) VY dG-(E)/dE]` |

## Performance

All performance numbers below are taken directly from committed benchmark and
profiling reports in the repository:

1. **Block-domain vs. global sparse KPM scaling on 64 cores**
   (source: [`bench/block_domains/RESULTS.md`](bench/block_domains/RESULTS.md),
   measured on `jh-tcn`, $2\times$ AMD EPYC 7551, 64 cores across 8 NUMA nodes,
   dimension $N = 10^6$, DOS $M = 256, R = 4$, conductivity $M = 128, R = 1$,
   `OMP_PROC_BIND=spread`):
   - *Without neighbour-block coupling:*
     - Block size $5$ ($\mathrm{nnz}(H) = 5\times 10^6$): global DOS scales from
       $9.81\ \mathrm{s}$ (1 core) to $0.61\ \mathrm{s}$ (64 cores, $16.0\times$
       speedup, $25\%$ efficiency); local-domain DOS scales from
       $11.59\ \mathrm{s}$ to $0.24\ \mathrm{s}$ ($47.5\times$ speedup, $74\%$
       efficiency). Global conductivity scales from $89.42\ \mathrm{s}$ to
       $3.54\ \mathrm{s}$ ($25.2\times$, $39\%$); local-domain conductivity scales
       from $83.97\ \mathrm{s}$ to $1.63\ \mathrm{s}$ ($51.6\times$, $81\%$).
     - Block size $10$ ($\mathrm{nnz}(H) = 10^7$): global DOS $16.44\ \mathrm{s}
       \to 0.74\ \mathrm{s}$ ($22.1\times$, $35\%$); local DOS $17.47\ \mathrm{s}
       \to 0.35\ \mathrm{s}$ ($50.1\times$, $78\%$). Global conductivity
       $92.64\ \mathrm{s} \to 3.71\ \mathrm{s}$ ($25.0\times$, $39\%$); local
       conductivity $86.88\ \mathrm{s} \to 1.63\ \mathrm{s}$ ($53.2\times$,
       $83\%$).
     - Block size $30$ ($\mathrm{nnz}(H) = 29{,}999{,}700$): global DOS
       $44.32\ \mathrm{s} \to 1.28\ \mathrm{s}$ ($34.7\times$, $54\%$); local DOS
       $41.65\ \mathrm{s} \to 0.88\ \mathrm{s}$ ($47.4\times$, $74\%$). Global
       conductivity $109.20\ \mathrm{s} \to 4.50\ \mathrm{s}$ ($24.3\times$,
       $38\%$); local conductivity $100.04\ \mathrm{s} \to 1.87\ \mathrm{s}$
       ($53.6\times$, $84\%$).
     - Relative moment difference $\max|\Delta\mu|/\max|\mu|$ between local and
       global at 64 cores is $4.4\times 10^{-16}$ to $4.2\times 10^{-15}$, and
       between 1 core and 64 cores (global) is $2.4\times 10^{-14}$ to
       $6.9\times 10^{-14}$.
   - *With neighbour-block coupling* (one random inter-block link per orbital):
     global DOS at 64 cores runs in $0.71\ \mathrm{s}$ (block 5, $18.6\times$),
     $0.81\ \mathrm{s}$ (block 10, $24.7\times$), and $1.34\ \mathrm{s}$ (block
     30, $35.6\times$); global conductivity at 64 cores runs in $3.57\ \mathrm{s}$
     ($25.5\times$), $3.86\ \mathrm{s}$ ($24.3\times$), and $4.66\ \mathrm{s}$
     ($23.8\times$).

2. **Native kernel scaling and reduced SIESTA material profiles**
   (source: [`handoff/0009-wrapup/profiling-results.md`](handoff/0009-wrapup/profiling-results.md),
   measured on `Aron`, $2\times$ AMD EPYC 7552, 96 physical cores across 2 sockets):
   - *Unprofiled `native` kernel scaling* (median of 3 processes, close core
     binding, 1 BLAS thread):
     - `linalg::MomentTable` (dimension $165{,}888$, batch $128$):
       $1835.867\ \mathrm{ms}$ (1 core), $290.141\ \mathrm{ms}$ (8 cores),
       $147.750\ \mathrm{ms}$ (24 cores, $12.43\times$), $133.418\ \mathrm{ms}$
       (48 cores, $13.76\times$), and $174.241\ \mathrm{ms}$ (96 cores,
       $10.54\times$; 96 cores is $30.6\%$ slower than 48 cores because a
       $128\times 128$ table tiled at $16\times 16$ exposes 64 tasks).
     - `ContractRank2Imag` ($M = 256$, $7680$ points): $181.227\ \mathrm{ms}$
       (1 core), $23.438\ \mathrm{ms}$ (8 cores), $9.050\ \mathrm{ms}$ (24 cores,
       $20.03\times$), $6.475\ \mathrm{ms}$ (48 cores, $27.99\times$), and
       $4.606\ \mathrm{ms}$ (96 cores, $39.35\times$).
   - *Reduced SIESTA material `cProfile` breakdown* (`golden` build, $24\times 24$
     mesh, `covariant` velocity, `exact` Löwdin metric, `shift` projection
     $[-5, 5]\ \mathrm{eV}$, $M = 256$ followed by $M = 1024$, $R = 16$, 8 OpenMP
     threads):
     - $\mathrm{Bi}$: total $304.753\ \mathrm{s}$ — moments + preparation
       $207.748\ \mathrm{s}$ ($68.2\%$), SHC reconstruction $85.793\ \mathrm{s}$
       ($28.2\%$, of which `_core.ContractRank2Imag` is $66.844\ \mathrm{s}$ and
       `legacy_delta_table` is $15.206\ \mathrm{s}$), other $11.212\ \mathrm{s}$,
       peak RSS $2.79\ \mathrm{GiB}$.
     - $\mathrm{PdSe_2}$: total $1070.146\ \mathrm{s}$ — moments + preparation
       $859.114\ \mathrm{s}$ ($80.3\%$), SHC reconstruction $84.509\ \mathrm{s}$
       ($7.9\%$, of which `_core.ContractRank2Imag` is $65.343\ \mathrm{s}$ and
       `legacy_delta_table` is $15.399\ \mathrm{s}$), other $126.523\ \mathrm{s}$,
       peak RSS $17.00\ \mathrm{GiB}$.

## Verification

Every stage of LSQUANT is gated against committed golden masters and analytic or
exact-diagonalisation oracles by [`scripts/check_goldens.sh`](scripts/check_goldens.sh),
with stage verification figures generated by [`figures/stage_*.py`](figures/) and
developer self-checks under [`tests/dev/`](tests/dev/). OpenMP reductions in the
`golden` profile use a fixed per-state accumulation order so that multi-threaded
runs are **bit-for-bit identical** across thread counts (`OMP_NUM_THREADS=1` vs.
`8`).

## Project layout

| Path | Contents |
|---|---|
| [`lsquant/core/algebra/`](lsquant/core/algebra) | `LinearOperator`, `SparseMatrix`, `BlockDiagonalOperator`, `BlochTransform` (batched FFTW), `MetricPolicy`, `OverlapAppliers` (`exact`, `chebyshev`, `ischo`), `StateGenerator`, `SpectralBounds`, `FftMomentContraction` |
| [`lsquant/core/chebyshev/`](lsquant/core/chebyshev) | `SpectralWindow`, `RescaledOperator`, `ChebyshevRecursion`, `ChebyshevVectors`, `ChebyshevMoments`, `Plan`, `MomentEngine`, `Contraction` |
| [`lsquant/io/`](lsquant/io) | CSR text reader, `.chebmom` and legacy `.chebmom1D/2D/TD` readers/writers, curve text writer |
| [`lsquant/api/`](lsquant/api), [`lsquant/compat/`](lsquant/compat) | `pybind11` module `lsquant._core` and legacy GSL Bessel / text round-trip helpers |
| [`lsquant/apps/`](lsquant/apps) | C++ golden reproduction drivers (`stage1_moments`, `stage6_nonorth`) |
| [`python/lsquant/`](python/lsquant) | Python library: `words/` parser & compiler, SymPy `spectral_functions`, `kernels`, `reconstruct`, `bases`, `disorder_bases`, `metric_forms`, `models/` (`wannier`, `siesta`, `nonorthogonal`, `disorder`), `projection/`, `oracle/`, `apps/` golden reproduction drivers, `runfile`, `cli`, `mcp_server` |
| [`examples/`](examples) | Worked tutorials `01`–`06` and first-principles SIESTA spin-Hall examples `07`–`08` |
| [`runs/goldens/`](runs/goldens), [`goldens/`](goldens) | Committed golden-master inputs, reference outputs, and JSON run files (`G1`–`G10`, tutorials `01`–`06`, SIESTA/disorder/projection gates) |
| [`figures/`](figures) | Stage verification figure scripts (`stage_1.py` … `stage_9_local_projection.py`) |
| [`docs/`](docs), [`bench/`](bench) | Physics notes, executable SymPy derivations ([`docs/derivations/`](docs/derivations)), performance documentation ([`docs/performance/`](docs/performance/README.md)), and kernel benchmarks (`bench/kernel_bench.cpp`, `bench/block_domains/`) |

## References, citation, and licence

1. A. Weiße, G. Wellein, A. Alvermann, and H. Fehske, *The kernel polynomial
   method*, Rev. Mod. Phys. **78**, 275 (2006).
2. Z. Fan, J. H. García, A. W. Cummings, J. E. Barrios-Vargas, M. Panhans,
   A. Harju, F. Ortmann, and S. Roche, *Linear scaling quantum transport
   methodologies*, Phys. Rep. **893**, 1–123 (2021) ([arXiv:1811.07387](https://arxiv.org/abs/1811.07387)).

See [`LICENSE`](LICENSE) for licence terms. If LSQUANT contributes to your work,
please cite Reference [2] above.
