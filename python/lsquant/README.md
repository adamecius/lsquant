# `lsquant` — Python API, Run-File Schema, and MCP Reference

For an overview of LSQUANT, installation instructions, the observable word
catalogue, and the worked examples `01`–`08`, see the main repository
[**README.md**](../../README.md). This page documents the Python package API,
the JSON run-file schema accepted by `python -m lsquant run`, and the MCP server
surface.

---

## Python API Summary

- **High-level execution** (`lsquant.runfile`):
  - `execute_runfile(runfile_path_or_dict, *, out_dir_override=None)` $\to$ `dict` (`output_dir`, `resolved_runfile`, `moment_files`, `curve_files`, `written_files`)
  - `reconstruct_target(target_path, *, formula=None, broadening_meV=None, fermi_energies=None, out_dir=None)`
  - `inspect_target(path)`
- **Bases** (`lsquant.bases`, `lsquant.disorder_bases`):
  - `RealBasis(model, supercell=(1, 1, 1), *, disorder_w=0.0, seed=12345, disorder_seed=7, csr_dir=None)`
  - `BlochBasis(model, supercell=(1, 1, 1), *, disorder_w=0.0, seed=12345, disorder_seed=7)`
  - `EigenBasis(model, supercell=(1, 1, 1))`
  - `BlochMetricBasis(model, supercell, *, route="lowdin", applier="exact", velocity_mode="covariant", ...)`
  - `EigenMetricBasis(model, supercell, *, velocity_mode="covariant")`
  - `make_basis(basis_name, model, supercell=(1, 1, 1), ...)`
- **Models** (`lsquant.models`):
  - `load_wannier_model(directory, seed)` — reads Wannier90 `<seed>_hr.dat`, `.uc`, and `.xyz`.
  - `load_siesta_model(directory, label=None, *, velocity_mode="covariant", ...)` — reads SIESTA `.fdf` + `.TSHS`/`.HSX` and `.ion.xml` via `sisl`, with certified Gauss–Legendre midpoint quadrature for orbital dipole matrices $A_\alpha(\mathbf{R})$.
- **Words & Compilation** (`lsquant.words`):
  - `lsquant.words.parse(word_str)` $\to$ `Word` AST
  - `lsquant.words.compile(word, basis=None, params=None)` $\to$ `CompiledWord` (`plan`, `operator_names`, `operators`, `reconstruction`, `dt`)
  - `lsquant.words.describe(word, basis=None, params=None)` $\to$ `dict`
  - `lsquant.words.list_catalogue()` $\to$ `list[dict]`
- **Engine Execution & Reconstruction** (`lsquant.moments`, `lsquant.reconstruct`):
  - `dos_plan(num_moments, num_states, word, *, operator_index=None)` $\to$ `_core.Plan`
  - `kubo_plan(num_moments, *, left_index, right_index, word)` $\to$ `_core.Plan`
  - `run_plan(plan, hamiltonian, window, states, operators=(), batch_size=None, symmetrise=True, dt=None)` $\to$ `_core.ChebyshevMoments`
  - `dos_curve(moments, broadening_meV, alpha=0.95)`
  - `kubo_greenwood(moments, broadening_meV, alpha=0.95)`
  - `kubo_bastin(moments, broadening_meV, alpha=0.95, *, backend="tables")` (plus `kubo_bastin_sea`, `kubo_bastin_surf`, `kubo_fft_grid`)
  - `shc(moments, broadening_meV, num_cells, cell_area, *, backend="tables")`
  - `time_correlation_curve(moments, broadening_meV, fermi_energy)`
- **Spectral Projection & Downfolding** (`lsquant.projection`):
  - `spectral_projection(hamiltonian, overlap, window, *, kind="shift")` (`kind` in `"hard"` | `"shift"`)
  - `downfold_orthogonal(hamiltonian, overlap, window)` and `LocalMatrix`
  - `compare_windows(source, words, windows, energies, *, cell_area, ...)`

---

## Command-Line Interface (`python -m lsquant`)

```bash
python3 -m lsquant words
python3 -m lsquant words describe "Tr[VX delta(E-H) VY dG-(E)/dE]"
python3 -m lsquant run <run.json> [--out <output_dir>] [--json]
python3 -m lsquant reconstruct <moment_file_or_dir> [--formula <name>] [--broadening-meV <val>] [--fermi-energies <E...>] [--out <dir>] [--json]
python3 -m lsquant inspect <run_json_or_moment_file_or_dir>
```

---

## Run-File Schema (`python/lsquant/runfile.py`)

Every calculation is configured by a strictly validated JSON run file. Unknown keys
at the top level or inside any sub-object are rejected before any operator file is
loaded.

### Top-level sections

| Section / Field | Allowed Values / Type | Default | Description |
|---|---|---|---|
| `model.kind` | `"csr"` \| `"wannier"` \| `"siesta"` \| `"generator"` | `"csr"` | Model input source. |
| `model.dir` | path | `null` | Directory containing `<label>.<OP>.CSR` (`csr`), `<seed>_hr.dat` (`wannier`), or `<label>.fdf` + `.TSHS`/`.HSX` (`siesta`). |
| `model.label` | string | `null` | File prefix `<label>` for `csr` or `siesta` models. |
| `model.seed` | string | `null` | Wannier90 seed name when `model.kind == "wannier"`. |
| `model.operators` | `list[str]` | `[]` | Static operator names available in `model.dir`. |
| `model.csr_dir` | path | `null` | Pre-built supercell CSR directory for `model.kind == "wannier"`. |
| `model.name` | string | `null` | Generator name from `GENERATOR_CATALOGUE` when `model.kind == "generator"`. |
| `model.tutorial` | `1..6` | `null` | Tutorial driver selector when `model.kind == "generator"`. |
| `model.params` | `dict` | `{}` | Generator or SIESTA model parameters (e.g. `velocity_mode`: `"bare"` \| `"berry_connection"` \| `"covariant"`, `oracle`: `bool`). |
| `basis.name` | `"real"` \| `"bloch"` \| `"eigen"` | `"real"` | Representation basis. |
| `basis.supercell` | `[N1, N2, N3]` | `[1, 1, 1]` | Supercell or $k$-mesh dimensions. |
| `basis.disorder` | `{"W": float, "seed": int}` | `{"W": 0.0, "seed": 7}` | Real-space Anderson box disorder parameters. |
| `metric.policy` | `"orthogonal"` \| `"nonorthogonal"` | `"orthogonal"` | Basis overlap policy. |
| `metric.route` | `"lowdin"` \| `"generalised"` | `"lowdin"` | Non-orthogonal evaluation route ($S^{-1/2} H S^{-1/2}$ vs. $S^{-1} H$). |
| `metric.applier` | `"exact"` \| `"chebyshev"` \| `"ischo"` | `"exact"` | Overlap inverse / inverse-square-root applier. |
| `metric.order` | `int` | `60` | Chebyshev polynomial order for overlap appliers. |
| `metric.coeff_tol` | `float` | `1e-12` | Chebyshev coefficient truncation threshold for overlap appliers. |
| `metric.tol` | `float` | `1e-13` | PCG solver tolerance for `"ischo"`. |
| `metric.max_iter` | `int` | `200` | Maximum PCG iterations for `"ischo"`. |
| `metric.max_dim` | `int` | `4096` | Maximum dense matrix dimension for `"exact"` overlap appliers. |
| `projection` | `dict` \| `null` | `null` | Spectral window projection (`{"kind": "hard" \| "shift" \| "orthogonal", "window": [emin, emax], ...}`). |
| `word` | string | *(required)* | Primary observable word `Tr[...]`. |
| `words` | `list[str]` | `[word]` | Multiple observable words evaluated on the same basis in one run. |
| `expansion.moments` | `int` | `64` | Number of Chebyshev moments $M$ per expanded slot. |
| `expansion.window` | `"gershgorin"` \| `"eigen"` \| `[lo, hi]` | `"gershgorin"` (`real`), `"eigen"` (`bloch`, `eigen`) | Spectral window mode or explicit `[lo, hi]` bounds in eV. |
| `expansion.trace` | `"exact"` \| `"random"` | `"exact"` | Trace evaluation mode. |
| `expansion.states` | `int` \| `null` | `null` | Number of trace states (`dim` for `exact`, `1` for `random` when omitted). |
| `expansion.seed` | `int` | `12345` | Random-phase state generator seed (`KPM_SEED`). |
| `expansion.batch` | `int` \| `null` | `moments` | Engine vector batch size. |
| `time.steps` | `int` | — | Number of time steps (`> 1`) for time-dependent words. |
| `time.tmax` | `float` \| `null` | `20.0` | Maximum evolution time in fs. |
| `time.dt` | `float` \| `null` | `null` | Explicit time step in fs as an alternative to `tmax = dt * (steps - 1)`. |
| `time.bessel` | `"scipy"` \| `"legacy"` | `"scipy"` | Bessel coefficient backend (`"legacy"` uses GSL values for legacy byte identity). |
| `reconstruction.kernel` | `"jackson"` \| `"lorentz"` \| `"none"` | `"jackson"` | Chebyshev damping kernel. |
| `reconstruction.broadening_meV` | `float` \| `null` | `null` | Target broadening in meV. |
| `reconstruction.fermi_energies` | `list[float]` | `[]` | Fermi energies in eV for time-correlation reconstruction. |
| `reconstruction.formulas` | `list[str]` \| `null` | shape default | Reconstruction formulas (`dos`, `greenwood`, `bastin`, `sea`, `surf`, `greenwood_fft`, `bastin_fft`, `shc`, `time_correlation`, `row_m0`). |
| `reconstruction.num_points` | `int` \| `null` | `null` | Explicit energy grid point count. |
| `reconstruction.fft_grid` | `bool` | `false` | Reconstruct on the Chebyshev-node FFT grid. |
| `reconstruction.legacy_text_round_trip` | `bool` | `false` | Round-trip moments through the 15-digit legacy ASCII format before reconstruction. |
| `reconstruction.per_site_dos` | `bool` | `false` | Also write standalone per-site `DOS_reconstructed.dat`. |
| `reconstruction.targets` | `list[dict]` | `[]` | Additional `(formulas, broadening_meV, subdir, also_top_dir)` reconstruction targets from the same moments. |
| `output.dir` | path | `"."` | Output directory (always receives `run.resolved.json`). |
| `output.legacy_names` | `bool` | `true` | Use legacy file naming conventions when writing `.chebmom*` and `.dat` files. |
| `output.legacy_layout` | `bool` | `true` | Write `.chebmom1D`/`.chebmom2D`/`.chebmomTD` text files. |
| `output.write_chebmom` | `bool` | `false` | Also write the unified `.chebmom` file. |
| `output.write_basis_json` | `bool` | `false` | Write `<stem>.basis.json` provenance metadata. |
| `output.capture_tutorial` | `bool` | `false` | Copy outputs into `captured/NNNN_*` for tutorial golden verification. |

### Generator catalogue (`model.kind = "generator"`)

| `model.name` | Module path |
|---|---|
| `make_chain` | `examples/01_dos_1d_chain/make_chain.py` |
| `make_disordered_chain` | `examples/02_localization_1d_chain/make_disordered_chain.py` |
| `make_haldane` | `examples/03_haldane_hall/make_haldane.py` |
| `make_chain_mag` | `examples/04_spin_precession_magnetic_chain/make_chain_mag.py` |
| `make_square` | `examples/05_charge_transport_square/make_square.py` |
| `make_ey` | `examples/06_spin_relaxation_square/make_ey.py` |
| `make_rashba` | `examples/06_spin_relaxation_square/make_rashba.py` |
| `gen_models` | `python/lsquant/models/wannier.py` |

---

## MCP Server (`python -m lsquant.mcp_server`)

`lsquant.mcp_server` exposes five tools over `stdio` using `FastMCP` (`pip install mcp`):

- `words_describe(word)` — parses and compiles `word`, returning the same JSON structure as `python -m lsquant words describe "<word>"`.
- `run(runfile_json)` — validates `runfile_json` and dispatches `python -m lsquant run` in a background subprocess, returning `{"run_id", "log", "output_dir"}`.
- `status(run_id)` — polls a background calculation started by `run`.
- `reconstruct(output_dir, formula, broadening_meV, fermi_energies)` — reconstructs curves from moment files in `output_dir`.
- `inspect(path)` — inspects a moment file, run file, or output directory.

All file operations are restricted to `LSQUANT_MCP_ROOT` (defaults to the repository root or `--root <dir>`).

```json
{
  "mcpServers": {
    "lsquant": {
      "command": "python3",
      "args": ["-m", "lsquant.mcp_server"],
      "env": {
        "PYTHONPATH": "/path/to/lsquant/build/golden/lsquant_core:/path/to/lsquant/python",
        "LSQUANT_MCP_ROOT": "/path/to/lsquant",
        "OMP_NUM_THREADS": "16",
        "KPM_SEED": "12345"
      }
    }
  }
}
```
