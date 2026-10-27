# Stage 7 Report — Observable Word Language, Run File, CLI, and MCP Server

Branch: `handoff/0007-agi` (merged `origin/handoff/0006-agi` and `origin/cpp-unified`).

---

## 1. Files Created and Modified

| File | Purpose |
|------|---------|
| `python/lsquant/words/grammar.py` | Hand-written tokenizer (`tokenize`) and recursive-descent parser (`parse`) producing `Word` AST (`OperatorSlot`, `SpectralSlot`, `CommutatorSlot`) with basis-transform sandwich folding (`B^+ O B`, `U^+ O U`). |
| `python/lsquant/words/catalogue.py` | Catalogue table (`CATALOGUE`) of the 7 recognised observable word shapes, basis-consistency validator (`validate_word_bases`), and nearest-shape diagnostic matcher (`match_catalogue`, `nearest_catalogue_entry`). |
| `python/lsquant/words/compiler.py` | Word compiler (`compile`, `describe`, `list_catalogue`) mapping `(Word, basis, params)` to `(Plan, operators, ReconstructionSpec)`, including orthogonal and non-orthogonal (`lowdin`, `generalised`) operator wrapping. |
| `python/lsquant/words/__init__.py` | Public exports of `lsquant.words` (`parse`, `compile`, `describe`, `CATALOGUE`, AST classes). |
| `python/lsquant/moments.py` | Extended `dos_plan` with optional `operator_index` keyword argument for `Tr[A delta(E-H)]` (`A != 1`). |
| `python/lsquant/runfile.py` | `RunFile` dataclass, strict schema validator (`validate_runfile_dict`) refusing unknown keys before any file load, defaults table (`DEFAULT_RUNFILE_DEFAULTS`), execution engine (`execute_runfile`), curve reconstructor (`reconstruct_target`), and inspector (`inspect_target`). |
| `python/lsquant/cli.py` | `argparse` CLI (`run`, `reconstruct`, `inspect`, `words`). |
| `python/lsquant/__main__.py` | Entry point for `python -m lsquant`. |
| `python/lsquant/mcp_server.py` | Stdio MCP server (`words_describe`, `run`, `status`, `reconstruct`, `inspect`) sharing validation and execution logic with the CLI and enforcing root-directory confinement. |
| `python/lsquant/README.md` | CLI usage and Claude Code `mcp.json` configuration entry. |
| `python/pyproject.toml` | Added `[project.scripts]` console entry `lsquant = "lsquant.cli:main"`. |
| `python/lsquant/__init__.py` | Re-exported `words`, `runfile`, and `cli`. |
| `python/lsquant/apps/tutorial01.py`, `python/lsquant/apps/tutorial03.py` | Added top-level 120 meV curve and `captured/NNNN_*` outputs for Tutorial 01 and Tutorial 03. |
| `runs/goldens/*.json` | 36 run files covering G1, G2, G3, G4, G4c, G5, G6, G7, G9 (both models and all bases/disorders), G10 (all 4 fixtures on both routes), and Tutorials 01–05. |
| `tests/dev/stage_7_check.py` | Single developer self-check covering (a) catalogue canonical round-trip, (b) 5 malformed/unsupported refusals, (c) CLI `g1.json` and `g5.json` byte identity vs G1 and G5, and (d) MCP server subprocess equivalence vs CLI. |

Zero C++ files under `lsquant/core/`, `lsquant/io/`, or `lsquant/api/` were modified. The legacy tree (`src/`, `include/`, `test/`) was untouched.

---

## 2. Word Language Grammar (EBNF)

```ebnf
word          ::= "Tr[" body "]" | "sea[" "Tr[" product "]" "]" ;
body          ::= "sea[" product "]" | product ;
product       ::= factor ( ["*"] factor )* ;
factor        ::= spectral | commutator | transform_seq | operator ;
spectral      ::= "delta(E-H)"
                | "G+(E)" | "G-(E)"
                | "dG-(E)/dE" | "dG+(E)/dE"
                | "U(t)^+" | "U(t)" ;
commutator    ::= "[" axis "," "U(t)" "]^+"
                | "[" axis "," "U(t)" "]" ;
axis          ::= "X" | "Y" | "Z" | "VX" | "VY" | "VZ" ;
transform_seq ::= "B^+" operator "B"
                | "B" operator "B^+"
                | "U^+" operator "U"
                | "U" operator "U^+" ;
operator      ::= ("1" | "P:" ident | ident) ["^+"] ["@" basis_tag] ;
ident         ::= LETTER (LETTER | DIGIT | "_")* ;
basis_tag     ::= "real" | "bloch" | "eigen" ;
```

---

## 3. Catalogue of Recognised Word Shapes

| `shape_id` | Canonical Word | Plan Builder | Default Reconstruction | Legacy Moment Prefix | Legacy Curve Prefix | Rank | Time-Dependent |
|---|---|---|---|---|---|---:|---|
| `dos_identity` | `Tr[delta(E-H)]` | `dos_plan` | `dos_curve` | `SpectralOp1` | `mean` | 1 | No |
| `dos_operator` | `Tr[A delta(E-H)]` | `dos_plan` (`operator_index=0`) | `dos_curve` | `SpectralOp<A>` | `mean` | 1 | No |
| `kubo_greenwood` | `Tr[VX delta(E-H) VX delta(E-H)]` | `kubo_plan` | `kubo_greenwood` | `NonEqOp<A>-<B>` | `KuboGreenWood_` | 2 | No |
| `kubo_bastin` | `Tr[VX delta(E-H) VY dG-(E)/dE]` (or `Tr[sea[...]]`) | `kubo_plan` | `kubo_bastin` (`sea`, `surf`) | `NonEqOp<A>-<B>` | `KuboBastin*_` | 2 | No |
| `msd` | `Tr[[X,U(t)]^+ delta(E-H) [X,U(t)]]` | `msd_plan` | `time_correlation_curve` | `Correlation<V>` | `mean...EF` | 1 | Yes |
| `projected_evolution` | `Tr[P U(t)^+ SX U(t) delta(E-H) P]` | `projected_evolution_plan` | `row_m0` | `TimeProj<A>-<P>` | `TimeProj` | 1 | Yes |
| `time_correlation` | `Tr[A delta(E-H) U(t)^+ B U(t)]` | `time_correlation_plan` | `time_correlation_curve` | `TimeCorr<A>-<B>` | `mean...EF` | 1 | Yes |

---

## 4. Run-File Schema and Defaults

Unknown keys at the top level or inside any sub-object raise `ValueError` during `validate_runfile_dict` before any file is loaded.

```json
{
  "model": {
    "kind": "csr | wannier | generator",
    "dir": "<path>",
    "label": "<csr_label>",
    "seed": "<wannier_seed>",
    "operators": ["VX", "VY"],
    "csr_dir": null,
    "name": "make_chain | make_disordered_chain | make_haldane | make_square | make_chain_mag | gen_models",
    "tutorial": null,
    "params": {}
  },
  "basis": {
    "name": "real | bloch | eigen",
    "supercell": [1, 1, 1],
    "disorder": {"W": 0.0, "seed": 7}
  },
  "metric": {
    "policy": "orthogonal | nonorthogonal",
    "route": "lowdin | generalised",
    "applier": "exact | chebyshev | ischo",
    "order": 60,
    "coeff_tol": 1e-12,
    "tol": 1e-13,
    "max_iter": 200,
    "max_dim": 4096
  },
  "word": "Tr[...]",
  "words": ["Tr[...]"],
  "expansion": {
    "moments": 64,
    "window": "gershgorin (on real) | eigen (on bloch, eigen) | [lo, hi]",
    "trace": "exact | random",
    "states": null,
    "seed": 12345,
    "batch": null,
    "state_label": null
  },
  "time": {
    "steps": 64,
    "tmax": 20.0,
    "dt": null,
    "bessel": "scipy | legacy"
  },
  "reconstruction": {
    "kernel": "jackson | lorentz | none",
    "broadening_meV": null,
    "fermi_energies": [],
    "formulas": null,
    "num_points": null,
    "fft_grid": false,
    "legacy_text_round_trip": false,
    "per_site_dos": false,
    "targets": []
  },
  "output": {
    "dir": ".",
    "legacy_names": true,
    "legacy_layout": true,
    "write_chebmom": false,
    "write_basis_json": false,
    "capture_tutorial": false
  }
}
```

Every run writes `run.resolved.json` in its output directory (including the resolved `expansion.window` value from `default_window_for_basis(basis.name)` when not explicitly supplied).

---

## 5. Run Files under `runs/goldens/` and Golden Mapping

| Run File | Golden / Tutorial Reproduced |
|---|---|
| `runs/goldens/g1.json` | **G1**: `goldens/stage0/chain1d/SpectralOp1chain1dKPM_M64_stateexact.chebmom1D` |
| `runs/goldens/g2.json` | **G2**: `goldens/stage0/chain1d/NonEqOpVX-VXchain1dKPM_M64x64_stateexact.chebmom2D` |
| `runs/goldens/g3.json` | **G3**: `goldens/stage0/chain1d/meanSpectralOp1chain1dKPM_M64_stateexactJACKSON.dat`, `DOS_reconstructed.dat`, and G1 `.chebmom1D` |
| `runs/goldens/g4.json` | **G4**: `goldens/stage0/graphene/NonEqOpVX-VXgraphene_goldenKPM_M64x64_statedefault.chebmom2D`, `KuboBastin_...JACKSON.dat` (64 meV), `sea_surf/KuboBastin{Sea,Surf}_...JACKSON.dat` (10 meV), `fft_grid/Kubo{GreenWood,Bastin}_...JACKSON_FFTgrid.dat` (10 meV) |
| `runs/goldens/g4c.json` | **G4c**: `goldens/stage0/graphene/offcentre_enclosing/NonEqOpVX-VYgraphene_goldenKPM_M128x128_stateexact.chebmom2D` and `KuboBastin_...JACKSON.dat` (10 meV) |
| `runs/goldens/g5.json` | **G5**: `goldens/stage0/haldane/NonEqOpVX-{VY,VX}haldaneKPM_M128x128_stateexact.chebmom2D` and `Kubo{Bastin,BastinSea,BastinSurf,GreenWood}_...JACKSON.dat` (10 meV) |
| `runs/goldens/g6.json` | **G6**: `goldens/stage0/chain1d_ballistic/CorrelationVXchain1dKPM_M256_state20exact.chebmomTD` and the 4 `mean...EF{0.000000,0.500000,1.000000,1.500000}JACKSON.dat` curves |
| `runs/goldens/g7.json` (`g7_sx.json`, `g7_sy.json`, `g7_sz.json`) | **G7**: `goldens/stage0/chain1d_mag/TimeProj{SX,SY,SZ}-PXchain1d_magKPM_M32_stateexact.chebmomTD` |
| `runs/goldens/g9_graphene_{real,bloch,bloch_exact,eigen,real_w1,bloch_w1}.json` | **G9 (graphene)**: `goldens/stage0/bases/graphene/` DOS and `sigma_xx` across `real`, `bloch`, and `eigen` bases at `W = 0` and `W = 1` eV |
| `runs/goldens/g9_haldane_{real,bloch,bloch_exact,eigen,real_w1,bloch_w1}.json` | **G9 (haldane)**: `goldens/stage0/bases/haldane/` DOS and `sigma_xx` across `real`, `bloch`, and `eigen` bases at `W = 0` and `W = 1` eV |
| `runs/goldens/g10_chain1d_{real,complex}_{lowdin,generalised}.json` | **G10 (chain1d)**: `goldens/stage0/nonorth/chain1d/{real,complex}` DOS and Kubo moments on `lowdin` and `generalised` routes vs G1/G2 |
| `runs/goldens/g10_haldane_{real,complex}_{lowdin,generalised}.json` | **G10 (haldane)**: `goldens/stage0/nonorth/haldane/{real,complex}` `VX-VY` and `VX-VX` moments and Hall curves on `lowdin` and `generalised` routes vs G5 |
| `runs/goldens/tutorial01.json` | **Tutorial 01**: `goldens/stage0/tutorials/01/` (`captured/0000_*`..`0021_*`, `eta{10,20,24,80}/`, top-level `.chebmom1D` and 120 meV `.dat`) |
| `runs/goldens/tutorial02.json` | **Tutorial 02**: `goldens/stage0/tutorials/02/` (`captured/0000_*`..`0191_*` and top-level `N = 512` `.dat`) |
| `runs/goldens/tutorial03.json` | **Tutorial 03**: `goldens/stage0/tutorials/03/` (`captured/0000_*`..`0003_*` and top-level 4 Haldane `.dat` curves) |
| `runs/goldens/tutorial04.json` | **Tutorial 04**: `goldens/stage0/tutorials/04/` (`captured/0000_*`..`0002_*` and top-level 3 `TimeProj` `.chebmomTD` files) |
| `runs/goldens/tutorial05.json` | **Tutorial 05**: `goldens/stage0/tutorials/05/` (`captured/0000_*`..`0007_*` and top-level 8 MSD and DOS `.dat` curves) |

---

## 6. Self-Check Output and MCP Transcript (`tests/dev/stage_7_check.py` on TCN)

Executed on TCN (`~/lsquant-agi-build`):

```text
$ OMP_NUM_THREADS=16 KPM_SEED=12345 PYTHONPATH=build/lsquant_core:python python3 tests/dev/stage_7_check.py
== (a) Catalogue canonical parse/print round-trip ==
  PASS  Tr[delta(E-H)]                             -> shape=dos_identity         plan_builder=dos_plan
  PASS  Tr[A delta(E-H)]                           -> shape=dos_operator         plan_builder=dos_plan
  PASS  Tr[VX delta(E-H) VX delta(E-H)]            -> shape=kubo_greenwood       plan_builder=kubo_plan
  PASS  Tr[VX delta(E-H) VY dG-(E)/dE]             -> shape=kubo_bastin          plan_builder=kubo_plan
  PASS  Tr[[X,U(t)]^+ delta(E-H) [X,U(t)]]         -> shape=msd                  plan_builder=msd_plan
  PASS  Tr[P U(t)^+ SX U(t) delta(E-H) P]          -> shape=projected_evolution  plan_builder=projected_evolution_plan
  PASS  Tr[A delta(E-H) U(t)^+ B U(t)]             -> shape=time_correlation     plan_builder=time_correlation_plan
  PASS  Tr[sea[VX delta(E-H) VY dG-(E)/dE]]        -> shape=kubo_bastin          plan_builder=kubo_plan
  PASS  Tr[P:PX U(t)^+ SX U(t) delta(E-H) P:PX]    -> shape=projected_evolution  plan_builder=projected_evolution_plan
  PASS  Tr[B^+ VX@bloch B delta(E-H)]              -> shape=dos_operator         plan_builder=dos_plan

== (b) Refusal of 5 malformed / unsupported inputs ==
  PASS  [unknown_operator] refused with: Unknown operator 'FOO' in word 'Tr[FOO delta(E-H)]' (registered operators: ['1', 'VX', 'VY'])
  PASS  [rank_3] refused with: Unsupported word 'Tr[VX delta(E-H) VY delta(E-H) VZ delta(E-H)]' of rank 3 (3 expanded spectral slots): rank >= 3 moment tensors are out of scope; nearest recognised shape is 'Tr[VX delta(E-H) VX delta(E-H)]' (kubo_greenwood)
  PASS  [mixed_basis_without_transform] refused with: Mixed basis without explicit transform: operator 'VX@bloch' is in 'bloch' basis while calculation/active basis is 'real' (insert explicit B^+ ... B or U^+ ... U transform)
  PASS  [ut_without_time] refused with: Word 'Tr[A delta(E-H) U(t)^+ B U(t)]' contains time evolution U(t) without 'time' parameters (steps and tmax required)
  PASS  [unknown_key_in_runfile] refused with: Unknown key 'unknown_extra_key' in run file (top-level); allowed keys are ['basis', 'expansion', 'metric', 'model', 'output', 'reconstruction', 'time', 'word', 'words']
  PASS  [unknown_generator] refused with: Unknown generator name 'evil_script' in run file (model.name); allowed generators are ['gen_models', 'make_chain', 'make_chain_mag', 'make_disordered_chain', 'make_haldane', 'make_square']
  PASS  [window_default_by_basis] real->gershgorin, bloch->eigen, eigen->eigen
  PASS  [runfile_schema] all 36 files in runs/goldens/*.json validated

== (c) CLI reproduction of G1 (runs/goldens/g1.json) and G5 (runs/goldens/g5.json) ==
  PASS  cmp IDENTICAL: G1 (SpectralOp1chain1dKPM_M64_stateexact.chebmom1D)
  PASS  cmp IDENTICAL: G5 (NonEqOpVX-VYhaldaneKPM_M128x128_stateexact.chebmom2D)
  PASS  cmp IDENTICAL: G5 (NonEqOpVX-VXhaldaneKPM_M128x128_stateexact.chebmom2D)
  PASS  cmp IDENTICAL: G5 (KuboBastin_NonEqOpVX-VYhaldaneKPM_M128x128_stateexactJACKSON.dat)
  PASS  cmp IDENTICAL: G5 (KuboBastinSea_NonEqOpVX-VYhaldaneKPM_M128x128_stateexactJACKSON.dat)
  PASS  cmp IDENTICAL: G5 (KuboBastinSurf_NonEqOpVX-VYhaldaneKPM_M128x128_stateexactJACKSON.dat)
  PASS  cmp IDENTICAL: G5 (KuboGreenWood_NonEqOpVX-VXhaldaneKPM_M128x128_stateexactJACKSON.dat)

== (d) MCP server subprocess vs CLI JSON equivalence ==
  MCP tools registered: ['words_describe', 'run', 'status', 'reconstruct', 'inspect']
  MCP -> call_tool: {"tool": "words_describe", "arguments": {"word": "Tr[VX delta(E-H) VY dG-(E)/dE]"}}
  MCP <- words_describe response:
    {
      "word": {
        "canonical": "Tr[VX delta(E-H) VY dG-(E)/dE]",
        "sea": false,
        "rank": 2,
        "time_dependent": false,
        "slots": [
          {
            "type": "operator",
            "name": "VX",
            "operator_key": "VX"
          },
          {
            "type": "spectral",
            "kind": "delta",
            "token": "delta(E-H)",
            "mode": "expand"
          },
          {
            "type": "operator",
            "name": "VY",
            "operator_key": "VY"
          },
          {
            "type": "spectral",
            "kind": "dG-/dE",
            "token": "dG-(E)/dE",
            "mode": "expand"
          }
        ]
      },
      "catalogue": {
        "shape_id": "kubo_bastin",
        "canonical_example": "Tr[VX delta(E-H) VY dG-(E)/dE]",
        "plan_builder": "kubo_plan",
        "default_reconstruction": "kubo_bastin",
        "default_formulas": [
          "bastin"
        ],
        "engine_word": "VX delta VY delta",
        "legacy_moment_prefix": "NonEqOpVX-VY",
        "legacy_moment_ext": ".chebmom2D",
        "operators": [
          {
            "type": "operator",
            "name": "VX",
            "operator_key": "VX"
          },
          {
            "type": "operator",
            "name": "VY",
            "operator_key": "VY"
          }
        ]
      },
      "plan": {
        "word": "VX delta VY delta",
        "instructions": [
          {
            "op": "Expand",
            "num_moments": 64
          },
          {
            "op": "ApplyOperator",
            "operator_index": 1,
            "operator_name": "VY"
          },
          {
            "op": "MirrorToBra"
          },
          {
            "op": "ApplyOperator",
            "operator_index": 0,
            "operator_name": "VX"
          },
          {
            "op": "Expand",
            "num_moments": 64
          },
          {
            "op": "Trace",
            "num_weights": 1
          }
        ]
      },
      "operators": [
        "VX",
        "VY"
      ],
      "reconstruction": {
        "shape_id": "kubo_bastin",
        "default_reconstruction": "kubo_bastin",
        "formulas": [
          "bastin"
        ],
        "kernel": "jackson",
        "broadening_meV": null,
        "fermi_energies": [],
        "num_points": null,
        "use_fft_grid": false,
        "legacy_moment_prefix": "NonEqOpVX-VY",
        "legacy_moment_ext": ".chebmom2D",
        "legacy_text_round_trip": false
      },
      "dt": null
    }
  PASS  MCP words_describe JSON == CLI words describe JSON
  MCP -> call_tool: {"tool": "inspect", "arguments": {"path": "goldens/stage0/chain1d/SpectralOp1chain1dKPM_M64_stateexact.chebmom1D"}}
  MCP <- inspect response:
    {
      "kind": "moments",
      "path": "/home/ICN2/jgarcia/lsquant-agi-build/goldens/stage0/chain1d/SpectralOp1chain1dKPM_M64_stateexact.chebmom1D",
      "rank": 1,
      "dims": [
        64
      ],
      "system_size": 2048,
      "number_of_states": 0,
      "state_label": "exact",
      "word": "",
      "window": [
        -2.0,
        2.0
      ],
      "band_width": 4.0,
      "half_width": 2.0,
      "band_center": 0.0,
      "kpm_weighted": true,
      "symmetrised": false,
      "number_of_time_steps": 0,
      "time_step": 0.0
    }
  PASS  MCP inspect JSON == CLI inspect JSON

STAGE 7 SELF-CHECK PASS
```

---

## 7. Deviations and Open Questions

- **Deviations (schema extensions beyond `brief.md` section 4 and tutorial app additions)**:
  - **`words` (top-level list of observable words)**: Added so a single run file can evaluate multiple observable words on the same loaded model/basis in one output directory. Required by **G7** (`runs/goldens/g7.json`, evaluating `SX`, `SY`, `SZ`), **G9** (`runs/goldens/g9_*.json`, evaluating DOS and `sigma_xx`), **G10 chain1d** (`runs/goldens/g10_chain1d_*.json`, evaluating DOS and Kubo `VX-VX`), **Tutorial 04** (`runs/goldens/tutorial04.json`), and **Tutorial 05** (`runs/goldens/tutorial05.json`, evaluating MSD and DOS).
  - **`model.kind = "generator"` with `name`, `tutorial`, and `params`**: Added to invoke the repository's model generators from a fixed catalogue (`GENERATOR_CATALOGUE` in `runfile.py`: `make_chain`, `make_disordered_chain`, `make_haldane`, `make_square`, `make_chain_mag`, `gen_models`, resolved from the repository root) for multi-size / multi-disorder tutorial sweeps. Required by **Tutorial 01** (`runs/goldens/tutorial01.json`, `N = 256, 1024, 4096`), **Tutorial 02** (`runs/goldens/tutorial02.json`, `N = 256, 1024, 4096` and `16` disorder realisations), and **Tutorial 05** (`runs/goldens/tutorial05.json`, `L = 200` and `W = 0.5, 1.0, 2.0, 4.0`).
  - **`model.csr_dir`**: Added to pass the pre-built `wannier2sparse` supercell CSR directory when `model.kind == "wannier"` and `basis.name == "real"`. Required by **G9** (`runs/goldens/g9_*_real*.json`) to match the Stage-5 `wannier2sparse` real-space CSR operators byte-for-byte.
  - **`metric.coeff_tol`, `metric.max_iter`, and `metric.max_dim`**: Added to expose the Chebyshev truncation threshold, PCG iteration cap, and dense solver dimension cap established in Stage 6 review 2. Used by **G10** (`runs/goldens/g10_*.json`).
  - **`expansion.state_label`**: Added to allow overriding the `_state<label>` filename token when a golden's filename label differs from the generator's default.
  - **`time.dt`**: Added to allow specifying the time step `dt` directly as an alternative to `tmax = dt * (steps - 1)`.
  - **`reconstruction.kernel = "lorentz" | "none"`**: Added alongside `"jackson"` to expose the Lorentz and bare (undamped) kernels implemented in `python/lsquant/kernels.py`.
  - **`reconstruction.num_points` and `reconstruction.fft_grid`**: Added to expose explicit energy grid point counts and the Chebyshev-node FFT grid reconstruction path (`Kubo*_JACKSON_FFTgrid.dat`).
  - **`reconstruction.per_site_dos`**: Added to emit the standalone per-site `DOS_reconstructed.dat` curve alongside the legacy `meanSpectralOp1...JACKSON.dat` file. Required by **G3** (`runs/goldens/g3.json`).
  - **`reconstruction.targets`**: Added to run multiple reconstruction passes (with distinct `formulas`, `broadening_meV`, `subdir`, and `also_top_dir`) from a single computed moment tensor. Required by **G4** (`runs/goldens/g4.json`), which reconstructs the same `VX-VX` moment tensor at `64` meV in the top directory and at `10` meV in `sea_surf/` (`sea`, `surf`) and `fft_grid/` (`greenwood_fft`, `bastin_fft`).
  - **`output.write_chebmom`**: Added so a run with `legacy_layout: true` can also write the unified `.chebmom` file in the same directory. Required by **G9** (`runs/goldens/g9_*.json`) and **G10** (`runs/goldens/g10_*.json`).
  - **`output.write_basis_json`**: Added to write `<stem>.basis.json` provenance metadata (`window_source`, `window`, `supercell`, `operators`, `disorder_w`). Required by **G9** (`runs/goldens/g9_*.json`).
  - **`output.capture_tutorial`**: Added to copy generated tutorial outputs into `captured/NNNN_*` matching the Stage-0 tutorial capture layout. Required by **Tutorials 01–05** (`runs/goldens/tutorial01.json`..`tutorial05.json`).
  - **Tutorial apps `python/lsquant/apps/tutorial01.py` and `python/lsquant/apps/tutorial03.py` output additions**: Updated `tutorial01.py` to also write the README step-3 top-level `120` meV curve (`meanSpectralOp1chain1d_N256KPM_M512_stateexact_chain1d_N256JACKSON.dat`), the `moments_chain1d_N<N>_M512.chebmom1D` copies, and `captured/0000_*`..`0021_*`, and updated `tutorial03.py` to also populate `captured/0000_*`..`0003_*`, so that running `tutorial01.json` and `tutorial03.json` reproduces the full `goldens/stage0/tutorials/01/` and `goldens/stage0/tutorials/03/` golden trees.
  - **`mcp` package class name**: `pip install --user mcp` on TCN installed `mcp 2.3.0`, where `FastMCP` was renamed to `MCPServer` (`from mcp.server.mcpserver import MCPServer`). `python/lsquant/mcp_server.py` imports `MCPServer` with a fallback to `FastMCP` so both `mcp 2.x` and `mcp 1.x` work transparently.
  - **Legacy tree deletion (`brief.md` section 8)**: Executed in commit `7687487` after merging `origin/cpp-unified` (`3265265`) and verifying tag `legacy-final` (`05e2f65`). See section **Legacy tree deleted** below for the fresh CMake build and Stage 1–7 self-check outputs on the deleted tree.
- **Open Questions**: None.

---

## Legacy tree deleted

Following Fable's authorisation of `brief.md` section 8 (`legacy-final` tag at `05e2f65` on `origin`), `origin/cpp-unified` was merged (`3265265`) and the legacy tree was removed in a single commit (`7687487`):

- Deleted `src/`, `include/`, `test/`, `scripts/regenerate_goldens.sh`, and `examples/utilities/kwant2wannier/run_example.sh` (keeping all Python model generators used by the tutorials and `GENERATOR_CATALOGUE` in `python/lsquant/runfile.py`, as well as `lsquant/compat/`).
- Removed the legacy targets (`add_subdirectory(src)`, `add_subdirectory(test)`, `add_subdirectory(utilities)`, `add_subdirectory(docs)`, and `include_directories(${PROJECT_SOURCE_DIR}/include)`) from the top-level `CMakeLists.txt`, retaining the `lsquant/` targets, Eigen3, OpenMP, GSL, Threads, FFTW3, pybind11, and the stage apps.
- Replaced `scripts/regenerate_goldens.sh` with a section in `goldens/stage0/PROVENANCE.md` documenting how to regenerate Stage-0 goldens from a detached worktree of the tag `legacy-final`.
- Updated `tests/dev/stage_2_check.py` and `tests/dev/stage_3_check.py` to load `chain1d` operators from `goldens/stage0/chain1d/operators` instead of the deleted `test/golden/chain1d/operators`, and updated `tests/dev/stage_4_check.py` check (d) to verify against `inline_compute-kpm-CorrOp` if present and against the exact SHA-256 of the `legacy-final` `.chebmomTD` output (`84047151f40de0933e725fbfbf56627bc1e5b33f1bba7ec9671da1908e661445`) when the legacy binary is absent.

After `rsync --delete` to `~/lsquant-agi-build` on TCN (`ssh jh-x`), a fresh build (`rm -rf build`) with FFTW and all seven stage self-checks (`stage_1_check`, `stage_2_check.py`, `stage_3_check.py`, `stage_4_check.py`, `stage_5_check.py`, `stage_6_check.py`, `stage_7_check.py`) passed:

```text
=== 0. Verify deleted paths on TCN ===
ABSENT: src
ABSENT: include
ABSENT: test
ABSENT: scripts/regenerate_goldens.sh
ABSENT: examples/utilities/kwant2wannier/run_example.sh

=== 1. Fresh CMake configure and build from scratch (rm -rf build) ===
$ cmake -S . -B build -DLSQUANT_BUILD_NEW_CORE=ON -DLSQUANT_BUILD_PYTHON=ON -Dpybind11_DIR=/home/ICN2/jgarcia/.local/lib/python3.10/site-packages/pybind11/share/cmake/pybind11
-- The CXX compiler identification is GNU 11.4.0
-- Detecting CXX compiler ABI info
-- Detecting CXX compiler ABI info - done
-- Check for working CXX compiler: /usr/bin/c++ - skipped
-- Detecting CXX compile features
-- Detecting CXX compile features - done
-- Found OpenMP_CXX: -fopenmp (found version "4.5") 
-- Found OpenMP: TRUE (found version "4.5")  
-- Found PkgConfig: /usr/bin/pkg-config (found version "0.29.2") 
-- Found GSL: /usr/include (found version "2.7.1") 
-- Looking for C++ include pthread.h
-- Looking for C++ include pthread.h - found
-- Performing Test CMAKE_HAVE_LIBC_PTHREAD
-- Performing Test CMAKE_HAVE_LIBC_PTHREAD - Success
-- Found Threads: TRUE  
-- Checking for module 'fftw3'
--   Found fftw3, version 3.3.8
-- Found PythonInterp: /usr/bin/python3 (found suitable version "3.10.12", minimum required is "3.9") 
-- Found PythonLibs: /usr/lib/x86_64-linux-gnu/libpython3.10.so
-- Performing Test HAS_FLTO_AUTO
-- Performing Test HAS_FLTO_AUTO - Success
-- Found pybind11: /home/ICN2/jgarcia/.local/lib/python3.10/site-packages/pybind11/include (found version "3.1.0")
-- Configuring done
-- Generating done
-- Build files have been written to: /home/ICN2/jgarcia/lsquant-agi-build/build
$ cmake --build build -j 32
[  3%] Building CXX object lsquant_core/CMakeFiles/lsquant_core.dir/core/algebra/block_diagonal_operator.cpp.o
[ 10%] Building CXX object lsquant_core/CMakeFiles/lsquant_core.dir/core/algebra/overlap_appliers.cpp.o
[  7%] Building CXX object lsquant_core/CMakeFiles/lsquant_core.dir/core/algebra/sparse_matrix.cpp.o
[ 17%] Building CXX object lsquant_core/CMakeFiles/lsquant_core.dir/core/algebra/bloch_transform.cpp.o
[ 21%] Building CXX object lsquant_core/CMakeFiles/lsquant_core.dir/core/algebra/spectral_bounds.cpp.o
[ 25%] Building CXX object lsquant_core/CMakeFiles/lsquant_core.dir/core/algebra/composite_operators.cpp.o
[ 25%] Building CXX object lsquant_core/CMakeFiles/lsquant_core.dir/core/algebra/state_generator.cpp.o
[ 28%] Building CXX object lsquant_core/CMakeFiles/lsquant_core.dir/core/chebyshev/moment_engine.cpp.o
[ 35%] Building CXX object lsquant_core/CMakeFiles/lsquant_core.dir/io/csr_text.cpp.o
[ 32%] Building CXX object lsquant_core/CMakeFiles/lsquant_core.dir/core/algebra/linalg.cpp.o
[ 39%] Building CXX object lsquant_core/CMakeFiles/lsquant_core.dir/io/chebmom_file.cpp.o
[ 46%] Building CXX object lsquant_core/CMakeFiles/lsquant_core.dir/core/chebyshev/chebyshev_recursion.cpp.o
[ 46%] Building CXX object lsquant_core/CMakeFiles/lsquant_core.dir/core/algebra/metric_policy.cpp.o
[ 53%] Building CXX object lsquant_core/CMakeFiles/lsquant_core.dir/core/chebyshev/rescaled_operator.cpp.o
[ 53%] Building CXX object lsquant_core/CMakeFiles/lsquant_core.dir/core/chebyshev/chebyshev_moments.cpp.o
[ 57%] Building CXX object lsquant_core/CMakeFiles/lsquant_core.dir/core/chebyshev/contraction.cpp.o
[ 60%] Building CXX object lsquant_core/CMakeFiles/lsquant_core.dir/core/chebyshev/chebyshev_vectors.cpp.o
[ 64%] Building CXX object lsquant_core/CMakeFiles/lsquant_core.dir/io/curve_text.cpp.o
[ 67%] Linking CXX static library ../lib/liblsquant_core.a
[ 67%] Built target lsquant_core
[ 71%] Building CXX object lsquant_core/CMakeFiles/stage1_moments.dir/apps/stage1_moments.cpp.o
[ 75%] Building CXX object lsquant_core/CMakeFiles/stage6_nonorth.dir/apps/stage6_nonorth.cpp.o
[ 78%] Building CXX object lsquant_core/CMakeFiles/stage_1_check.dir/__/tests/dev/stage_1_check.cpp.o
[ 82%] Building CXX object lsquant_core/CMakeFiles/_core.dir/api/bindings.cpp.o
[ 85%] Building CXX object lsquant_core/CMakeFiles/_core.dir/compat/legacy_kernels.cpp.o
[ 89%] Linking CXX executable stage1_moments
[ 89%] Built target stage1_moments
[ 92%] Linking CXX executable stage6_nonorth
[ 92%] Built target stage6_nonorth
[ 96%] Linking CXX executable stage_1_check
[ 96%] Built target stage_1_check
[100%] Linking CXX shared module _core.cpython-310-x86_64-linux-gnu.so
[100%] Built target _core

=== 2. Stage 1 check (build/lsquant_core/stage_1_check) ===
dos max_abs_error: 2.47025e-15
kubo max_abs_error: 1.13687e-13
OMP_NUM_THREADS=1 vs 8 bit-identical: YES

=== 3. Stage 2 check (tests/dev/stage_2_check.py) ===
(a) dos_plan array_equal vs G1:  True
(a) kubo_plan array_equal vs G2: True
(b) DeltaFunction.table vs legacy_delta_table max_abs_diff = 2.220446e-16
(c) dos_curve vs G3 differing lines: 0 (total lines: 1590)

=== 4. Stage 3 check (tests/dev/stage_3_check.py) ===
(0) SparseMatrix.from_csr (scipy) dos array_equal vs ReadCsrText: True
(a) Haldane VX-VY array_equal vs G5: True
(a) Haldane VX-VX array_equal vs G5: True
(b) Graphene Gershgorin BandWidth = 15.0366847652085, array_equal vs G4: True
(c) GreenRetarded.table vs compat max_abs_diff           = 1.148853e-14
(c) GreenRetardedDerivative.table vs compat max_abs_diff = 2.702165e-13 (max_rel_diff = 7.173405e-15)
(d) G5 Bastin (10 meV)           differing lines: 0 (total lines: 3839)
(d) G5 Sea (10 meV)              differing lines: 0 (total lines: 3839)
(d) G5 Surf (10 meV)             differing lines: 0 (total lines: 3839)
(d) G5 Greenwood (10 meV)        differing lines: 0 (total lines: 127)
(d) G4 Bastin (64 meV)           differing lines: 0 (total lines: 1919)
(d) G4 Sea (10 meV)              differing lines: 0 (total lines: 1919)
(d) G4 Surf (10 meV)             differing lines: 0 (total lines: 1919)
(d) G4 Bastin FFTgrid (10 meV)   differing lines: 0 (total lines: 63)
(d) G4 Greenwood FFTgrid (10 meV) differing lines: 0 (total lines: 63)

=== 5. Stage 4 check (tests/dev/stage_4_check.py) ===
(a) G7 chain1d_mag Gershgorin window = [-2.31, 2.31], BandWidth = 4.62
(a) G7 TimeProjSX-PX array_equal vs golden: True
(a) G7 TimeProjSY-PX array_equal vs golden: True
(a) G7 TimeProjSZ-PX array_equal vs golden: True
(b) G6 CorrelationVXchain1d array_equal vs golden: True
(b) G6 EF=+0.0 curve differing lines: 0 (total lines: 64)
(b) G6 EF=+0.5 curve differing lines: 0 (total lines: 64)
(b) G6 EF=+1.0 curve differing lines: 0 (total lines: 64)
(b) G6 EF=+1.5 curve differing lines: 0 (total lines: 64)
(c) compat GSL Bessel vs SciPy jv max_abs_diff = 5.551115e-16
(d) time_correlation_plan chain1d VX-VX (M=64, steps=16, tmax=5, exact) array_equal vs inline_compute-kpm-CorrOp: True

=== 6. Stage 5 check (tests/dev/stage_5_check.py) ===
[check a] BlochTransform vs dense B (6,4,1) graphene: rel(B x) = 3.563e-15, rel(B^+ x) = 3.221e-15
[check a] BlochTransform unitarity (8,8,1) graphene: max|B B^+ x - x| = 5.648e-17, max|B^+ B x - x| = 5.722e-17
[check b] (graphene 6x4x1) max|H_real_basis - H_w2s_csr| = 0.000e+00
[check b] (graphene 6x4x1) right-size CSR with wrong receipt raises ValueError: PASS
[check b] (graphene 6x4x1) max_x ||H_k x - B H_csr B^+ x||_inf = 4.780e-16
[check c] (graphene 6x4x1) max_x ||v_x x - B VX_csr B^+ x||_inf = 2.483e-16
[check b] (haldane 6x4x1) max|H_real_basis - H_w2s_csr| = 0.000e+00
[check b] (haldane 6x4x1) max_x ||H_k x - B H_csr B^+ x||_inf = 4.100e-16
[check c] (haldane 6x4x1) max_x ||v_x x - B VX_csr B^+ x||_inf = 2.582e-16
[check c] (haldane 48x48x1) rel ||v_x x - B VX_w2s B^+ x||_inf = 6.990e-16
[check d] Anderson disorder W=1.0 (seed=12345): max|mu_real - mu_bloch| DOS(64) = 4.719e-16, VX-VX(64x64) = 1.556e-15
[check e] (16,16,1) real vs bloch: DOS max|d_mu|=7.439e-15, max_rel|d_curve|=2.943e-15; KuboXX max|d_mu|=4.352e-14, max|d_curve/N|=1.399e-14, max_rel|d_curve|=9.851e-15
[check e] (16,16,1) real vs eigen: DOS max|d_mu|=7.438e-15, max_rel|d_curve|=2.968e-15; KuboXX max|d_mu|=4.308e-14, max|d_curve/N|=1.343e-14, max_rel|d_curve|=9.457e-15
[check e] (16,16,1) bloch vs eigen: DOS max|d_mu|=2.109e-15, max_rel|d_curve|=9.730e-16; KuboXX max|d_mu|=7.994e-15, max|d_curve/N|=1.278e-15, max_rel|d_curve|=9.001e-16
[check e] Jackson KPM (M=128) vs Gaussian exact-diag: rel_L2(DOS)=1.150e-02, rel_L2(KuboXX)=4.129e-02
[check e] (graphene 12x12x1) EigenBasis.exact_dos/exact_sigma_xx vs hr_exactdiag (Gamma mesh): rel_dos=1.822e-15, rel_sigma_xx=2.758e-15
[check f] G1 chain1d DOS byte-identical to golden: PASS
[check f] G4 graphene Kubo byte-identical to golden: PASS
[check f] G7 chain1d_mag spin precession byte-identical to golden: PASS
stage_5_check: ALL CHECKS PASSED

=== 7. Stage 6 check (tests/dev/stage_6_check.py) ===
=== Stage 6 Self-Check: Rotated chain1d (N = 2048) ===
[1] real: max |eigvalsh(H', S') - eigvalsh(H)| = 1.154632e-14
[4] real: max ||S^-1_cheb v - S^-1_exact v|| / ||S^-1_exact v|| (M_S=60) = 1.452447e-12
[2] real route=lowdin      applier=exact     applier_res=2.687185e-14 max|mu - G1|=7.959603e-13 max|mu - delta_m0|=7.959603e-13
[2] real route=generalised applier=exact     applier_res=3.918055e-16 max|mu - G1|=4.664319e-16 max|mu - delta_m0|=4.664319e-16
[3] real route=lowdin      applier=chebyshev (M_S=12, eff_order=12) applier_res=7.094662e-05 max|mu - delta_m0|=5.308428e-04
[3] real route=generalised applier=chebyshev (M_S=12, eff_order=12) applier_res=2.580411e-04 max|mu - delta_m0|=1.859680e-03
[3] real route=lowdin      applier=chebyshev (M_S=60, eff_order=38) applier_res=1.781045e-12 max|mu - delta_m0|=9.506879e-12
[3] real route=generalised applier=chebyshev (M_S=60, eff_order=41) applier_res=1.687868e-12 max|mu - delta_m0|=1.517385e-12
[3] real route=generalised applier=ischo     (tol=1e-06) applier_res=5.071275e-07 max|mu - delta_m0|=1.675476e-08
[3] real route=generalised applier=ischo     (tol=1e-13) applier_res=1.345871e-16 max|mu - delta_m0|=2.737152e-15
[1] complex: max |eigvalsh(H', S') - eigvalsh(H)| = 7.993606e-15
[4] complex: max ||S^-1_cheb v - S^-1_exact v|| / ||S^-1_exact v|| (M_S=60) = 9.664000e-13
[2] complex route=lowdin      applier=exact     applier_res=3.004363e-14 max|mu - G1|=6.202940e-13 max|mu - delta_m0|=6.202940e-13
[2] complex route=generalised applier=exact     applier_res=4.162287e-16 max|mu - G1|=1.022115e-15 max|mu - delta_m0|=1.022115e-15
[3] complex route=lowdin      applier=chebyshev (M_S=60, eff_order=33) applier_res=1.782907e-12 max|mu - delta_m0|=9.336730e-13
[3] complex route=generalised applier=chebyshev (M_S=60, eff_order=36) applier_res=1.138146e-12 max|mu - delta_m0|=2.915651e-12
[3] complex route=generalised applier=ischo     (tol=1e-13) applier_res=1.491737e-16 max|mu - delta_m0|=2.860823e-15
=== Stage 6 Self-Check PASSED ===

=== 8. Stage 7 check (tests/dev/stage_7_check.py) ===
== (a) Catalogue canonical parse/print round-trip ==
  PASS  Tr[delta(E-H)]                             -> shape=dos_identity         plan_builder=dos_plan
  PASS  Tr[A delta(E-H)]                           -> shape=dos_operator         plan_builder=dos_plan
  PASS  Tr[VX delta(E-H) VX delta(E-H)]            -> shape=kubo_greenwood       plan_builder=kubo_plan
  PASS  Tr[VX delta(E-H) VY dG-(E)/dE]             -> shape=kubo_bastin          plan_builder=kubo_plan
  PASS  Tr[[X,U(t)]^+ delta(E-H) [X,U(t)]]         -> shape=msd                  plan_builder=msd_plan
  PASS  Tr[P U(t)^+ SX U(t) delta(E-H) P]          -> shape=projected_evolution  plan_builder=projected_evolution_plan
  PASS  Tr[A delta(E-H) U(t)^+ B U(t)]             -> shape=time_correlation     plan_builder=time_correlation_plan
  PASS  Tr[sea[VX delta(E-H) VY dG-(E)/dE]]        -> shape=kubo_bastin          plan_builder=kubo_plan
  PASS  Tr[P:PX U(t)^+ SX U(t) delta(E-H) P:PX]    -> shape=projected_evolution  plan_builder=projected_evolution_plan
  PASS  Tr[B^+ VX@bloch B delta(E-H)]              -> shape=dos_operator         plan_builder=dos_plan

== (b) Refusal of 5 malformed / unsupported inputs ==
  PASS  [unknown_operator] refused with: Unknown operator 'FOO' in word 'Tr[FOO delta(E-H)]' (registered operators: ['1', 'VX', 'VY'])
  PASS  [rank_3] refused with: Unsupported word 'Tr[VX delta(E-H) VY delta(E-H) VZ delta(E-H)]' of rank 3 (3 expanded spectral slots): rank >= 3 moment tensors are out of scope; nearest recognised shape is 'Tr[VX delta(E-H) VX delta(E-H)]' (kubo_greenwood)
  PASS  [mixed_basis_without_transform] refused with: Mixed basis without explicit transform: operator 'VX@bloch' is in 'bloch' basis while calculation/active basis is 'real' (insert explicit B^+ ... B or U^+ ... U transform)
  PASS  [ut_without_time] refused with: Word 'Tr[A delta(E-H) U(t)^+ B U(t)]' contains time evolution U(t) without 'time' parameters (steps and tmax required)
  PASS  [unknown_key_in_runfile] refused with: Unknown key 'unknown_extra_key' in run file (top-level); allowed keys are ['basis', 'expansion', 'metric', 'model', 'output', 'reconstruction', 'time', 'word', 'words']
  PASS  [unknown_generator] refused with: Unknown generator name 'evil_script' in run file (model.name); allowed generators are ['gen_models', 'make_chain', 'make_chain_mag', 'make_disordered_chain', 'make_haldane', 'make_square']
  PASS  [window_default_by_basis] real->gershgorin, bloch->eigen, eigen->eigen
  PASS  [runfile_schema] all 36 files in runs/goldens/*.json validated

== (c) CLI reproduction of G1 (runs/goldens/g1.json) and G5 (runs/goldens/g5.json) ==
  PASS  cmp IDENTICAL: G1 (SpectralOp1chain1dKPM_M64_stateexact.chebmom1D)
  PASS  cmp IDENTICAL: G5 (NonEqOpVX-VYhaldaneKPM_M128x128_stateexact.chebmom2D)
  PASS  cmp IDENTICAL: G5 (NonEqOpVX-VXhaldaneKPM_M128x128_stateexact.chebmom2D)
  PASS  cmp IDENTICAL: G5 (KuboBastin_NonEqOpVX-VYhaldaneKPM_M128x128_stateexactJACKSON.dat)
  PASS  cmp IDENTICAL: G5 (KuboBastinSea_NonEqOpVX-VYhaldaneKPM_M128x128_stateexactJACKSON.dat)
  PASS  cmp IDENTICAL: G5 (KuboBastinSurf_NonEqOpVX-VYhaldaneKPM_M128x128_stateexactJACKSON.dat)
  PASS  cmp IDENTICAL: G5 (KuboGreenWood_NonEqOpVX-VXhaldaneKPM_M128x128_stateexactJACKSON.dat)

== (d) MCP server subprocess vs CLI JSON equivalence ==
  MCP tools registered: ['words_describe', 'run', 'status', 'reconstruct', 'inspect']
  MCP -> call_tool: {"tool": "words_describe", "arguments": {"word": "Tr[VX delta(E-H) VY dG-(E)/dE]"}}
  MCP <- words_describe response:
    {
      "word": {
        "canonical": "Tr[VX delta(E-H) VY dG-(E)/dE]",
        "sea": false,
        "rank": 2,
        "time_dependent": false,
        "slots": [
          {
            "type": "operator",
            "name": "VX",
            "operator_key": "VX"
          },
          {
            "type": "spectral",
            "kind": "delta",
            "token": "delta(E-H)",
            "mode": "expand"
          },
          {
            "type": "operator",
            "name": "VY",
            "operator_key": "VY"
          },
          {
            "type": "spectral",
            "kind": "dG-/dE",
            "token": "dG-(E)/dE",
            "mode": "expand"
          }
        ]
      },
      "catalogue": {
        "shape_id": "kubo_bastin",
        "canonical_example": "Tr[VX delta(E-H) VY dG-(E)/dE]",
        "plan_builder": "kubo_plan",
        "default_reconstruction": "kubo_bastin",
        "default_formulas": [
          "bastin"
        ],
        "engine_word": "VX delta VY delta",
        "legacy_moment_prefix": "NonEqOpVX-VY",
        "legacy_moment_ext": ".chebmom2D",
        "operators": [
          {
            "type": "operator",
            "name": "VX",
            "operator_key": "VX"
          },
          {
            "type": "operator",
            "name": "VY",
            "operator_key": "VY"
          }
        ]
      },
      "plan": {
        "word": "VX delta VY delta",
        "instructions": [
          {
            "op": "Expand",
            "num_moments": 64
          },
          {
            "op": "ApplyOperator",
            "operator_index": 1,
            "operator_name": "VY"
          },
          {
            "op": "MirrorToBra"
          },
          {
            "op": "ApplyOperator",
            "operator_index": 0,
            "operator_name": "VX"
          },
          {
            "op": "Expand",
            "num_moments": 64
          },
          {
            "op": "Trace",
            "num_weights": 1
          }
        ]
      },
      "operators": [
        "VX",
        "VY"
      ],
      "reconstruction": {
        "shape_id": "kubo_bastin",
        "default_reconstruction": "kubo_bastin",
        "formulas": [
          "bastin"
        ],
        "kernel": "jackson",
        "broadening_meV": null,
        "fermi_energies": [],
        "num_points": null,
        "use_fft_grid": false,
        "legacy_moment_prefix": "NonEqOpVX-VY",
        "legacy_moment_ext": ".chebmom2D",
        "legacy_text_round_trip": false
      },
      "dt": null
    }
  PASS  MCP words_describe JSON == CLI words describe JSON
  MCP -> call_tool: {"tool": "inspect", "arguments": {"path": "goldens/stage0/chain1d/SpectralOp1chain1dKPM_M64_stateexact.chebmom1D"}}
  MCP <- inspect response:
    {
      "kind": "moments",
      "path": "/home/ICN2/jgarcia/lsquant-agi-build/goldens/stage0/chain1d/SpectralOp1chain1dKPM_M64_stateexact.chebmom1D",
      "rank": 1,
      "dims": [
        64
      ],
      "system_size": 2048,
      "number_of_states": 0,
      "state_label": "exact",
      "word": "",
      "window": [
        -2.0,
        2.0
      ],
      "band_width": 4.0,
      "half_width": 2.0,
      "band_center": 0.0,
      "kpm_weighted": true,
      "symmetrised": false,
      "number_of_time_steps": 0,
      "time_step": 0.0
    }
  PASS  MCP inspect JSON == CLI inspect JSON

STAGE 7 SELF-CHECK PASS
```


---

## Review 2 applied

1. **Item 1 (`2055207`)**: Replaced all eleven `MetricPolicy` method calls (`overlap()`, `inverse_overlap()`, `inverse_sqrt_overlap()`) with property accesses (`overlap`, `inverse_overlap`, `inverse_sqrt_overlap`) and annotated them with `@property` in the type stub — `python/lsquant/runfile.py`, `python/lsquant/words/compiler.py`, `python/lsquant/_core.pyi`.
2. **Item 2 (`44ee28e`)**: Added `_get_or(mapping, key, default)` so keys present with value `None` (including `expansion.states`) fall back to their defaults across compilation and run-file resolution — `python/lsquant/words/compiler.py`, `python/lsquant/runfile.py`.
3. **Item 3 (`65a40c6`)**: Restricted the `kubo_bastin` catalogue entry to `dG-/dE` only (refusing `G+(E)`, `G-(E)`, and `dG+(E)/dE` with the nearest-shape message) and refused `expansion.window='eigen'` on `model.kind='csr'` models during validation and execution — `python/lsquant/words/catalogue.py`, `python/lsquant/runfile.py`.
4. **Item 4 (`71ddce5`)**: Parameterised `run_tutorial01`, `run_tutorial02`, and `run_tutorial05` to accept and honour all run-file `word`, `expansion`, `time`, `reconstruction`, and `model.params` settings while validating allowed keys per generator in `validate_runfile_dict` — `python/lsquant/apps/tutorial01.py`, `python/lsquant/apps/tutorial02.py`, `python/lsquant/apps/tutorial05.py`, `python/lsquant/runfile.py`, `runs/goldens/tutorial02.json`, `runs/goldens/tutorial05.json`.
5. **Item 5 (`3d69bab`)**: Removed the implicit Greenwood word injection from `execute_runfile`, refused any `reconstruction.formulas` entry whose required word is not declared in `word`/`words`, and listed both words explicitly in multi-formula run files — `python/lsquant/words/catalogue.py`, `python/lsquant/words/compiler.py`, `python/lsquant/runfile.py`, `runs/goldens/g5.json`, `runs/goldens/tutorial03.json`, `runs/goldens/tutorial05.json`, `runs/goldens/g10_haldane_real_lowdin.json`, `runs/goldens/g10_haldane_real_generalised.json`, `runs/goldens/g10_haldane_complex_lowdin.json`, `runs/goldens/g10_haldane_complex_generalised.json`.
6. **Item 6 (`bf49444`)**: Set `"csr_dir": "build/w2s_csr/<seed>"` in all four G9 real-basis run files and updated `RealBasis` to invoke `wannier2sparse` automatically when `csr_dir` is not yet populated — `python/lsquant/bases.py`, `runs/goldens/g9_graphene_real.json`, `runs/goldens/g9_graphene_real_w1.json`, `runs/goldens/g9_haldane_real.json`, `runs/goldens/g9_haldane_real_w1.json`.
7. **Item 7 (`7ae3510`)**: Wrapped MCP tool handler exceptions in `mcp.server.fastmcp.exceptions.ToolError(str(exc))` so error reasons are returned to MCP clients in `CallToolResult.content[0].text` (`isError=True`) — `python/lsquant/mcp_server.py`.
8. **Item 8 (`12255e9`)**: Added `python/lsquant/apps/tutorial06.py` and `runs/goldens/tutorial06.json` (reproducing `goldens/stage0/tutorials/06/`), added G10 run files for the `chebyshev` and `ischo` overlap appliers (`g10_chain1d_real_lowdin_chebyshev.json`, `g10_chain1d_real_generalised_chebyshev.json`, `g10_chain1d_real_generalised_ischo.json`), and extended `tests/dev/stage_7_check.py` to test the new refusals, MCP error reasons, `g9_graphene_eigen.json`, and `g10_chain1d_real_generalised.json` — `python/lsquant/apps/tutorial06.py`, `python/lsquant/runfile.py`, `runs/goldens/tutorial06.json`, `runs/goldens/g10_chain1d_real_lowdin_chebyshev.json`, `runs/goldens/g10_chain1d_real_generalised_chebyshev.json`, `runs/goldens/g10_chain1d_real_generalised_ischo.json`, `tests/dev/stage_7_check.py`.

### `tests/dev/stage_7_check.py` Output on TCN (`~/lsquant-agi-build`)

```text
== (a) Catalogue canonical parse/print round-trip ==
  PASS  Tr[delta(E-H)]                             -> shape=dos_identity         plan_builder=dos_plan
  PASS  Tr[A delta(E-H)]                           -> shape=dos_operator         plan_builder=dos_plan
  PASS  Tr[VX delta(E-H) VX delta(E-H)]            -> shape=kubo_greenwood       plan_builder=kubo_plan
  PASS  Tr[VX delta(E-H) VY dG-(E)/dE]             -> shape=kubo_bastin          plan_builder=kubo_plan
  PASS  Tr[[X,U(t)]^+ delta(E-H) [X,U(t)]]         -> shape=msd                  plan_builder=msd_plan
  PASS  Tr[P U(t)^+ SX U(t) delta(E-H) P]          -> shape=projected_evolution  plan_builder=projected_evolution_plan
  PASS  Tr[A delta(E-H) U(t)^+ B U(t)]             -> shape=time_correlation     plan_builder=time_correlation_plan
  PASS  Tr[sea[VX delta(E-H) VY dG-(E)/dE]]        -> shape=kubo_bastin          plan_builder=kubo_plan
  PASS  Tr[P:PX U(t)^+ SX U(t) delta(E-H) P:PX]    -> shape=projected_evolution  plan_builder=projected_evolution_plan
  PASS  Tr[B^+ VX@bloch B delta(E-H)]              -> shape=dos_operator         plan_builder=dos_plan

== (b) Refusal of 5 malformed / unsupported inputs ==
  PASS  [unknown_operator] refused with: Unknown operator 'FOO' in word 'Tr[FOO delta(E-H)]' (registered operators: ['1', 'VX', 'VY'])
  PASS  [rank_3] refused with: Unsupported word 'Tr[VX delta(E-H) VY delta(E-H) VZ delta(E-H)]' of rank 3 (3 expanded spectral slots): rank >= 3 moment tensors are out of scope; nearest recognised shape is 'Tr[VX delta(E-H) VX delta(E-H)]' (kubo_greenwood)
  PASS  [mixed_basis_without_transform] refused with: Mixed basis without explicit transform: operator 'VX@bloch' is in 'bloch' basis while calculation/active basis is 'real' (insert explicit B^+ ... B or U^+ ... U transform)
  PASS  [ut_without_time] refused with: Word 'Tr[A delta(E-H) U(t)^+ B U(t)]' contains time evolution U(t) without 'time' parameters (steps and tmax required)
  PASS  [unknown_key_in_runfile] refused with: Unknown key 'unknown_extra_key' in run file (top-level); allowed keys are ['basis', 'expansion', 'metric', 'model', 'output', 'reconstruction', 'time', 'word', 'words']
  PASS  [unsupported_green_word] refused with: Unrecognised word shape 'Tr[VX delta(E-H) VY G+(E)]'; nearest recognised shape is 'Tr[VX delta(E-H) VY dG-(E)/dE]' (kubo_bastin)
  PASS  [eigen_window_on_csr] refused with: expansion.window='eigen' is not supported for model.kind='csr' (no k blocks; use gershgorin or an explicit window)
  PASS  [unmatched_formula_for_word] refused with: Reconstruction formula 'greenwood' in reconstruction.formulas does not match any word in ['Tr[VX delta(E-H) VY dG-(E)/dE]'] (allowed formulas for given word(s): ['bastin', 'bastin_fft', 'kubo_bastin', 'kubo_bastin_fft', 'kubo_bastin_sea', 'kubo_bastin_surf', 'sea', 'surf'])
  PASS  [unknown_generator_param] refused with: Unknown key 'bogus_param' in run file (model.params (make_chain)); allowed keys are ['hopping_t', 'sizes']
  PASS  [unknown_generator] refused with: Unknown generator name 'evil_script' in run file (model.name); allowed generators are ['gen_models', 'make_chain', 'make_chain_mag', 'make_disordered_chain', 'make_ey', 'make_haldane', 'make_rashba', 'make_square']
  PASS  [window_default_by_basis] real->gershgorin, bloch->eigen, eigen->eigen
  PASS  [runfile_schema] all 40 files in runs/goldens/*.json validated

== (c) CLI reproduction of G1, G5, G9 (g9_graphene_eigen.json), and G10 (g10_chain1d_real_generalised.json) ==
  PASS  cmp IDENTICAL: G1 (SpectralOp1chain1dKPM_M64_stateexact.chebmom1D)
  PASS  cmp IDENTICAL: G5 (NonEqOpVX-VYhaldaneKPM_M128x128_stateexact.chebmom2D)
  PASS  cmp IDENTICAL: G5 (NonEqOpVX-VXhaldaneKPM_M128x128_stateexact.chebmom2D)
  PASS  cmp IDENTICAL: G5 (KuboBastin_NonEqOpVX-VYhaldaneKPM_M128x128_stateexactJACKSON.dat)
  PASS  cmp IDENTICAL: G5 (KuboBastinSea_NonEqOpVX-VYhaldaneKPM_M128x128_stateexactJACKSON.dat)
  PASS  cmp IDENTICAL: G5 (KuboBastinSurf_NonEqOpVX-VYhaldaneKPM_M128x128_stateexactJACKSON.dat)
  PASS  cmp IDENTICAL: G5 (KuboGreenWood_NonEqOpVX-VXhaldaneKPM_M128x128_stateexactJACKSON.dat)
  PASS  cmp IDENTICAL: G9 g9_graphene_eigen (KuboGreenWood_NonEqOpVX-VXgraphene_eigen_KPM_M256x256_stateexactJACKSON.dat)
  PASS  cmp IDENTICAL: G9 g9_graphene_eigen (NonEqOpVX-VXgraphene_eigen_KPM_M256x256_stateexact.basis.json)
  PASS  cmp IDENTICAL: G9 g9_graphene_eigen (NonEqOpVX-VXgraphene_eigen_KPM_M256x256_stateexact.chebmom)
  PASS  cmp IDENTICAL: G9 g9_graphene_eigen (NonEqOpVX-VXgraphene_eigen_KPM_M256x256_stateexact.chebmom2D)
  PASS  cmp IDENTICAL: G9 g9_graphene_eigen (SpectralOp1graphene_eigen_KPM_M256_stateexact.basis.json)
  PASS  cmp IDENTICAL: G9 g9_graphene_eigen (SpectralOp1graphene_eigen_KPM_M256_stateexact.chebmom)
  PASS  cmp IDENTICAL: G9 g9_graphene_eigen (SpectralOp1graphene_eigen_KPM_M256_stateexact.chebmom1D)
  PASS  cmp IDENTICAL: G9 g9_graphene_eigen (meanSpectralOp1graphene_eigen_KPM_M256_stateexactJACKSON.dat)
  PASS  cmp IDENTICAL: G10 g10_chain1d_real_generalised (KuboGreenWood_NonEqOpVX-VXchain1dKPM_M64x64_stateexactJACKSON.dat)
  PASS  cmp IDENTICAL: G10 g10_chain1d_real_generalised (NonEqOpVX-VXchain1dKPM_M64x64_stateexact.chebmom)
  PASS  cmp IDENTICAL: G10 g10_chain1d_real_generalised (NonEqOpVX-VXchain1dKPM_M64x64_stateexact.chebmom2D)
  PASS  cmp IDENTICAL: G10 g10_chain1d_real_generalised (SpectralOp1chain1dKPM_M64_stateexact.chebmom)
  PASS  cmp IDENTICAL: G10 g10_chain1d_real_generalised (SpectralOp1chain1dKPM_M64_stateexact.chebmom1D)
  PASS  cmp IDENTICAL: G10 g10_chain1d_real_generalised (meanSpectralOp1chain1dKPM_M64_stateexactJACKSON.dat)

== (d) MCP server subprocess vs CLI JSON equivalence ==
  MCP tools registered: ['words_describe', 'run', 'status', 'reconstruct', 'inspect']
  MCP -> call_tool: {"tool": "words_describe", "arguments": {"word": "Tr[VX delta(E-H) VY dG-(E)/dE]"}}
  MCP <- words_describe response:
    {
      "word": {
        "canonical": "Tr[VX delta(E-H) VY dG-(E)/dE]",
        "sea": false,
        "rank": 2,
        "time_dependent": false,
        "slots": [
          {
            "type": "operator",
            "name": "VX",
            "operator_key": "VX"
          },
          {
            "type": "spectral",
            "kind": "delta",
            "token": "delta(E-H)",
            "mode": "expand"
          },
          {
            "type": "operator",
            "name": "VY",
            "operator_key": "VY"
          },
          {
            "type": "spectral",
            "kind": "dG-/dE",
            "token": "dG-(E)/dE",
            "mode": "expand"
          }
        ]
      },
      "catalogue": {
        "shape_id": "kubo_bastin",
        "canonical_example": "Tr[VX delta(E-H) VY dG-(E)/dE]",
        "plan_builder": "kubo_plan",
        "default_reconstruction": "kubo_bastin",
        "default_formulas": [
          "bastin"
        ],
        "engine_word": "VX delta VY delta",
        "legacy_moment_prefix": "NonEqOpVX-VY",
        "legacy_moment_ext": ".chebmom2D",
        "operators": [
          {
            "type": "operator",
            "name": "VX",
            "operator_key": "VX"
          },
          {
            "type": "operator",
            "name": "VY",
            "operator_key": "VY"
          }
        ]
      },
      "plan": {
        "word": "VX delta VY delta",
        "instructions": [
          {
            "op": "Expand",
            "num_moments": 64
          },
          {
            "op": "ApplyOperator",
            "operator_index": 1,
            "operator_name": "VY"
          },
          {
            "op": "MirrorToBra"
          },
          {
            "op": "ApplyOperator",
            "operator_index": 0,
            "operator_name": "VX"
          },
          {
            "op": "Expand",
            "num_moments": 64
          },
          {
            "op": "Trace",
            "num_weights": 1
          }
        ]
      },
      "operators": [
        "VX",
        "VY"
      ],
      "reconstruction": {
        "shape_id": "kubo_bastin",
        "default_reconstruction": "kubo_bastin",
        "formulas": [
          "bastin"
        ],
        "kernel": "jackson",
        "broadening_meV": null,
        "fermi_energies": [],
        "num_points": null,
        "use_fft_grid": false,
        "legacy_moment_prefix": "NonEqOpVX-VY",
        "legacy_moment_ext": ".chebmom2D",
        "legacy_text_round_trip": false
      },
      "dt": null
    }
  PASS  MCP words_describe JSON == CLI words describe JSON
  MCP -> call_tool: {"tool": "inspect", "arguments": {"path": "goldens/stage0/chain1d/SpectralOp1chain1dKPM_M64_stateexact.chebmom1D"}}
  MCP <- inspect response:
    {
      "kind": "moments",
      "path": "/home/ICN2/jgarcia/lsquant-agi-build/goldens/stage0/chain1d/SpectralOp1chain1dKPM_M64_stateexact.chebmom1D",
      "rank": 1,
      "dims": [
        64
      ],
      "system_size": 2048,
      "number_of_states": 0,
      "state_label": "exact",
      "word": "",
      "window": [
        -2.0,
        2.0
      ],
      "band_width": 4.0,
      "half_width": 2.0,
      "band_center": 0.0,
      "kpm_weighted": true,
      "symmetrised": false,
      "number_of_time_steps": 0,
      "time_step": 0.0
    }
  PASS  MCP inspect JSON == CLI inspect JSON
  PASS  MCP error responses preserve reason: 'Error executing tool words_describe: Unrecognised word shape 'Tr[VX delta(E-H) VY G+(E)]'; nearest recognised shape is 'Tr[VX delta(E-H) VY dG-(E)/dE]' (kubo_bastin)'

STAGE 7 SELF-CHECK PASS
```

### G9 and G10 Run Files (`python -m lsquant run` + `cmp` on TCN)

```text
== G9 run files (12 run files under runs/goldens/g9_*.json) ==
IDENTICAL  runs/goldens/g9_graphene_bloch.json -> KuboGreenWood_NonEqOpVX-VXgraphene_bloch_KPM_M256x256_statedefaultJACKSON.dat
IDENTICAL  runs/goldens/g9_graphene_bloch.json -> meanSpectralOp1graphene_bloch_KPM_M256_statedefaultJACKSON.dat
IDENTICAL  runs/goldens/g9_graphene_bloch.json -> NonEqOpVX-VXgraphene_bloch_KPM_M256x256_statedefault.basis.json
IDENTICAL  runs/goldens/g9_graphene_bloch.json -> NonEqOpVX-VXgraphene_bloch_KPM_M256x256_statedefault.chebmom
IDENTICAL  runs/goldens/g9_graphene_bloch.json -> NonEqOpVX-VXgraphene_bloch_KPM_M256x256_statedefault.chebmom2D
IDENTICAL  runs/goldens/g9_graphene_bloch.json -> SpectralOp1graphene_bloch_KPM_M256_statedefault.basis.json
IDENTICAL  runs/goldens/g9_graphene_bloch.json -> SpectralOp1graphene_bloch_KPM_M256_statedefault.chebmom
IDENTICAL  runs/goldens/g9_graphene_bloch.json -> SpectralOp1graphene_bloch_KPM_M256_statedefault.chebmom1D
IDENTICAL  runs/goldens/g9_graphene_bloch_exact.json -> KuboGreenWood_NonEqOpVX-VXgraphene_bloch_KPM_M256x256_stateexactJACKSON.dat
IDENTICAL  runs/goldens/g9_graphene_bloch_exact.json -> meanSpectralOp1graphene_bloch_KPM_M256_stateexactJACKSON.dat
IDENTICAL  runs/goldens/g9_graphene_bloch_exact.json -> NonEqOpVX-VXgraphene_bloch_KPM_M256x256_stateexact.basis.json
IDENTICAL  runs/goldens/g9_graphene_bloch_exact.json -> NonEqOpVX-VXgraphene_bloch_KPM_M256x256_stateexact.chebmom
IDENTICAL  runs/goldens/g9_graphene_bloch_exact.json -> NonEqOpVX-VXgraphene_bloch_KPM_M256x256_stateexact.chebmom2D
IDENTICAL  runs/goldens/g9_graphene_bloch_exact.json -> SpectralOp1graphene_bloch_KPM_M256_stateexact.basis.json
IDENTICAL  runs/goldens/g9_graphene_bloch_exact.json -> SpectralOp1graphene_bloch_KPM_M256_stateexact.chebmom
IDENTICAL  runs/goldens/g9_graphene_bloch_exact.json -> SpectralOp1graphene_bloch_KPM_M256_stateexact.chebmom1D
IDENTICAL  runs/goldens/g9_graphene_bloch_w1.json -> KuboGreenWood_NonEqOpVX-VXgraphene_bloch_KPM_M256x256_statedefaultJACKSON.dat
IDENTICAL  runs/goldens/g9_graphene_bloch_w1.json -> meanSpectralOp1graphene_bloch_KPM_M256_statedefaultJACKSON.dat
IDENTICAL  runs/goldens/g9_graphene_bloch_w1.json -> NonEqOpVX-VXgraphene_bloch_KPM_M256x256_statedefault.basis.json
IDENTICAL  runs/goldens/g9_graphene_bloch_w1.json -> NonEqOpVX-VXgraphene_bloch_KPM_M256x256_statedefault.chebmom
IDENTICAL  runs/goldens/g9_graphene_bloch_w1.json -> NonEqOpVX-VXgraphene_bloch_KPM_M256x256_statedefault.chebmom2D
IDENTICAL  runs/goldens/g9_graphene_bloch_w1.json -> SpectralOp1graphene_bloch_KPM_M256_statedefault.basis.json
IDENTICAL  runs/goldens/g9_graphene_bloch_w1.json -> SpectralOp1graphene_bloch_KPM_M256_statedefault.chebmom
IDENTICAL  runs/goldens/g9_graphene_bloch_w1.json -> SpectralOp1graphene_bloch_KPM_M256_statedefault.chebmom1D
IDENTICAL  runs/goldens/g9_graphene_eigen.json -> KuboGreenWood_NonEqOpVX-VXgraphene_eigen_KPM_M256x256_stateexactJACKSON.dat
IDENTICAL  runs/goldens/g9_graphene_eigen.json -> meanSpectralOp1graphene_eigen_KPM_M256_stateexactJACKSON.dat
IDENTICAL  runs/goldens/g9_graphene_eigen.json -> NonEqOpVX-VXgraphene_eigen_KPM_M256x256_stateexact.basis.json
IDENTICAL  runs/goldens/g9_graphene_eigen.json -> NonEqOpVX-VXgraphene_eigen_KPM_M256x256_stateexact.chebmom
IDENTICAL  runs/goldens/g9_graphene_eigen.json -> NonEqOpVX-VXgraphene_eigen_KPM_M256x256_stateexact.chebmom2D
IDENTICAL  runs/goldens/g9_graphene_eigen.json -> SpectralOp1graphene_eigen_KPM_M256_stateexact.basis.json
IDENTICAL  runs/goldens/g9_graphene_eigen.json -> SpectralOp1graphene_eigen_KPM_M256_stateexact.chebmom
IDENTICAL  runs/goldens/g9_graphene_eigen.json -> SpectralOp1graphene_eigen_KPM_M256_stateexact.chebmom1D
IDENTICAL  runs/goldens/g9_graphene_real.json -> KuboGreenWood_NonEqOpVX-VXgraphene_real_KPM_M256x256_statedefaultJACKSON.dat
IDENTICAL  runs/goldens/g9_graphene_real.json -> meanSpectralOp1graphene_real_KPM_M256_statedefaultJACKSON.dat
IDENTICAL  runs/goldens/g9_graphene_real.json -> NonEqOpVX-VXgraphene_real_KPM_M256x256_statedefault.basis.json
IDENTICAL  runs/goldens/g9_graphene_real.json -> NonEqOpVX-VXgraphene_real_KPM_M256x256_statedefault.chebmom
IDENTICAL  runs/goldens/g9_graphene_real.json -> NonEqOpVX-VXgraphene_real_KPM_M256x256_statedefault.chebmom2D
IDENTICAL  runs/goldens/g9_graphene_real.json -> SpectralOp1graphene_real_KPM_M256_statedefault.basis.json
IDENTICAL  runs/goldens/g9_graphene_real.json -> SpectralOp1graphene_real_KPM_M256_statedefault.chebmom
IDENTICAL  runs/goldens/g9_graphene_real.json -> SpectralOp1graphene_real_KPM_M256_statedefault.chebmom1D
IDENTICAL  runs/goldens/g9_graphene_real_w1.json -> KuboGreenWood_NonEqOpVX-VXgraphene_real_KPM_M256x256_statedefaultJACKSON.dat
IDENTICAL  runs/goldens/g9_graphene_real_w1.json -> meanSpectralOp1graphene_real_KPM_M256_statedefaultJACKSON.dat
IDENTICAL  runs/goldens/g9_graphene_real_w1.json -> NonEqOpVX-VXgraphene_real_KPM_M256x256_statedefault.basis.json
IDENTICAL  runs/goldens/g9_graphene_real_w1.json -> NonEqOpVX-VXgraphene_real_KPM_M256x256_statedefault.chebmom
IDENTICAL  runs/goldens/g9_graphene_real_w1.json -> NonEqOpVX-VXgraphene_real_KPM_M256x256_statedefault.chebmom2D
IDENTICAL  runs/goldens/g9_graphene_real_w1.json -> SpectralOp1graphene_real_KPM_M256_statedefault.basis.json
IDENTICAL  runs/goldens/g9_graphene_real_w1.json -> SpectralOp1graphene_real_KPM_M256_statedefault.chebmom
IDENTICAL  runs/goldens/g9_graphene_real_w1.json -> SpectralOp1graphene_real_KPM_M256_statedefault.chebmom1D
IDENTICAL  runs/goldens/g9_haldane_bloch.json -> KuboGreenWood_NonEqOpVX-VXhaldane_bloch_KPM_M256x256_statedefaultJACKSON.dat
IDENTICAL  runs/goldens/g9_haldane_bloch.json -> meanSpectralOp1haldane_bloch_KPM_M256_statedefaultJACKSON.dat
IDENTICAL  runs/goldens/g9_haldane_bloch.json -> NonEqOpVX-VXhaldane_bloch_KPM_M256x256_statedefault.basis.json
IDENTICAL  runs/goldens/g9_haldane_bloch.json -> NonEqOpVX-VXhaldane_bloch_KPM_M256x256_statedefault.chebmom
IDENTICAL  runs/goldens/g9_haldane_bloch.json -> NonEqOpVX-VXhaldane_bloch_KPM_M256x256_statedefault.chebmom2D
IDENTICAL  runs/goldens/g9_haldane_bloch.json -> SpectralOp1haldane_bloch_KPM_M256_statedefault.basis.json
IDENTICAL  runs/goldens/g9_haldane_bloch.json -> SpectralOp1haldane_bloch_KPM_M256_statedefault.chebmom
IDENTICAL  runs/goldens/g9_haldane_bloch.json -> SpectralOp1haldane_bloch_KPM_M256_statedefault.chebmom1D
IDENTICAL  runs/goldens/g9_haldane_bloch_exact.json -> KuboGreenWood_NonEqOpVX-VXhaldane_bloch_KPM_M256x256_stateexactJACKSON.dat
IDENTICAL  runs/goldens/g9_haldane_bloch_exact.json -> meanSpectralOp1haldane_bloch_KPM_M256_stateexactJACKSON.dat
IDENTICAL  runs/goldens/g9_haldane_bloch_exact.json -> NonEqOpVX-VXhaldane_bloch_KPM_M256x256_stateexact.basis.json
IDENTICAL  runs/goldens/g9_haldane_bloch_exact.json -> NonEqOpVX-VXhaldane_bloch_KPM_M256x256_stateexact.chebmom
IDENTICAL  runs/goldens/g9_haldane_bloch_exact.json -> NonEqOpVX-VXhaldane_bloch_KPM_M256x256_stateexact.chebmom2D
IDENTICAL  runs/goldens/g9_haldane_bloch_exact.json -> SpectralOp1haldane_bloch_KPM_M256_stateexact.basis.json
IDENTICAL  runs/goldens/g9_haldane_bloch_exact.json -> SpectralOp1haldane_bloch_KPM_M256_stateexact.chebmom
IDENTICAL  runs/goldens/g9_haldane_bloch_exact.json -> SpectralOp1haldane_bloch_KPM_M256_stateexact.chebmom1D
IDENTICAL  runs/goldens/g9_haldane_bloch_w1.json -> KuboGreenWood_NonEqOpVX-VXhaldane_bloch_KPM_M256x256_statedefaultJACKSON.dat
IDENTICAL  runs/goldens/g9_haldane_bloch_w1.json -> meanSpectralOp1haldane_bloch_KPM_M256_statedefaultJACKSON.dat
IDENTICAL  runs/goldens/g9_haldane_bloch_w1.json -> NonEqOpVX-VXhaldane_bloch_KPM_M256x256_statedefault.basis.json
IDENTICAL  runs/goldens/g9_haldane_bloch_w1.json -> NonEqOpVX-VXhaldane_bloch_KPM_M256x256_statedefault.chebmom
IDENTICAL  runs/goldens/g9_haldane_bloch_w1.json -> NonEqOpVX-VXhaldane_bloch_KPM_M256x256_statedefault.chebmom2D
IDENTICAL  runs/goldens/g9_haldane_bloch_w1.json -> SpectralOp1haldane_bloch_KPM_M256_statedefault.basis.json
IDENTICAL  runs/goldens/g9_haldane_bloch_w1.json -> SpectralOp1haldane_bloch_KPM_M256_statedefault.chebmom
IDENTICAL  runs/goldens/g9_haldane_bloch_w1.json -> SpectralOp1haldane_bloch_KPM_M256_statedefault.chebmom1D
IDENTICAL  runs/goldens/g9_haldane_eigen.json -> KuboGreenWood_NonEqOpVX-VXhaldane_eigen_KPM_M256x256_stateexactJACKSON.dat
IDENTICAL  runs/goldens/g9_haldane_eigen.json -> meanSpectralOp1haldane_eigen_KPM_M256_stateexactJACKSON.dat
IDENTICAL  runs/goldens/g9_haldane_eigen.json -> NonEqOpVX-VXhaldane_eigen_KPM_M256x256_stateexact.basis.json
IDENTICAL  runs/goldens/g9_haldane_eigen.json -> NonEqOpVX-VXhaldane_eigen_KPM_M256x256_stateexact.chebmom
IDENTICAL  runs/goldens/g9_haldane_eigen.json -> NonEqOpVX-VXhaldane_eigen_KPM_M256x256_stateexact.chebmom2D
IDENTICAL  runs/goldens/g9_haldane_eigen.json -> SpectralOp1haldane_eigen_KPM_M256_stateexact.basis.json
IDENTICAL  runs/goldens/g9_haldane_eigen.json -> SpectralOp1haldane_eigen_KPM_M256_stateexact.chebmom
IDENTICAL  runs/goldens/g9_haldane_eigen.json -> SpectralOp1haldane_eigen_KPM_M256_stateexact.chebmom1D
IDENTICAL  runs/goldens/g9_haldane_real.json -> KuboGreenWood_NonEqOpVX-VXhaldane_real_KPM_M256x256_statedefaultJACKSON.dat
IDENTICAL  runs/goldens/g9_haldane_real.json -> meanSpectralOp1haldane_real_KPM_M256_statedefaultJACKSON.dat
IDENTICAL  runs/goldens/g9_haldane_real.json -> NonEqOpVX-VXhaldane_real_KPM_M256x256_statedefault.basis.json
IDENTICAL  runs/goldens/g9_haldane_real.json -> NonEqOpVX-VXhaldane_real_KPM_M256x256_statedefault.chebmom
IDENTICAL  runs/goldens/g9_haldane_real.json -> NonEqOpVX-VXhaldane_real_KPM_M256x256_statedefault.chebmom2D
IDENTICAL  runs/goldens/g9_haldane_real.json -> SpectralOp1haldane_real_KPM_M256_statedefault.basis.json
IDENTICAL  runs/goldens/g9_haldane_real.json -> SpectralOp1haldane_real_KPM_M256_statedefault.chebmom
IDENTICAL  runs/goldens/g9_haldane_real.json -> SpectralOp1haldane_real_KPM_M256_statedefault.chebmom1D
IDENTICAL  runs/goldens/g9_haldane_real_w1.json -> KuboGreenWood_NonEqOpVX-VXhaldane_real_KPM_M256x256_statedefaultJACKSON.dat
IDENTICAL  runs/goldens/g9_haldane_real_w1.json -> meanSpectralOp1haldane_real_KPM_M256_statedefaultJACKSON.dat
IDENTICAL  runs/goldens/g9_haldane_real_w1.json -> NonEqOpVX-VXhaldane_real_KPM_M256x256_statedefault.basis.json
IDENTICAL  runs/goldens/g9_haldane_real_w1.json -> NonEqOpVX-VXhaldane_real_KPM_M256x256_statedefault.chebmom
IDENTICAL  runs/goldens/g9_haldane_real_w1.json -> NonEqOpVX-VXhaldane_real_KPM_M256x256_statedefault.chebmom2D
IDENTICAL  runs/goldens/g9_haldane_real_w1.json -> SpectralOp1haldane_real_KPM_M256_statedefault.basis.json
IDENTICAL  runs/goldens/g9_haldane_real_w1.json -> SpectralOp1haldane_real_KPM_M256_statedefault.chebmom
IDENTICAL  runs/goldens/g9_haldane_real_w1.json -> SpectralOp1haldane_real_KPM_M256_statedefault.chebmom1D

== G10 run files (11 run files under runs/goldens/g10_*.json) ==
IDENTICAL  runs/goldens/g10_chain1d_complex_generalised.json -> KuboGreenWood_NonEqOpVX-VXchain1dKPM_M64x64_stateexactJACKSON.dat
IDENTICAL  runs/goldens/g10_chain1d_complex_generalised.json -> meanSpectralOp1chain1dKPM_M64_stateexactJACKSON.dat
IDENTICAL  runs/goldens/g10_chain1d_complex_generalised.json -> NonEqOpVX-VXchain1dKPM_M64x64_stateexact.chebmom
IDENTICAL  runs/goldens/g10_chain1d_complex_generalised.json -> NonEqOpVX-VXchain1dKPM_M64x64_stateexact.chebmom2D
IDENTICAL  runs/goldens/g10_chain1d_complex_generalised.json -> SpectralOp1chain1dKPM_M64_stateexact.chebmom
IDENTICAL  runs/goldens/g10_chain1d_complex_generalised.json -> SpectralOp1chain1dKPM_M64_stateexact.chebmom1D
IDENTICAL  runs/goldens/g10_chain1d_complex_lowdin.json -> KuboGreenWood_NonEqOpVX-VXchain1dKPM_M64x64_stateexactJACKSON.dat
IDENTICAL  runs/goldens/g10_chain1d_complex_lowdin.json -> meanSpectralOp1chain1dKPM_M64_stateexactJACKSON.dat
IDENTICAL  runs/goldens/g10_chain1d_complex_lowdin.json -> NonEqOpVX-VXchain1dKPM_M64x64_stateexact.chebmom
IDENTICAL  runs/goldens/g10_chain1d_complex_lowdin.json -> NonEqOpVX-VXchain1dKPM_M64x64_stateexact.chebmom2D
IDENTICAL  runs/goldens/g10_chain1d_complex_lowdin.json -> SpectralOp1chain1dKPM_M64_stateexact.chebmom
IDENTICAL  runs/goldens/g10_chain1d_complex_lowdin.json -> SpectralOp1chain1dKPM_M64_stateexact.chebmom1D
IDENTICAL  runs/goldens/g10_chain1d_real_generalised.json -> KuboGreenWood_NonEqOpVX-VXchain1dKPM_M64x64_stateexactJACKSON.dat
IDENTICAL  runs/goldens/g10_chain1d_real_generalised.json -> meanSpectralOp1chain1dKPM_M64_stateexactJACKSON.dat
IDENTICAL  runs/goldens/g10_chain1d_real_generalised.json -> NonEqOpVX-VXchain1dKPM_M64x64_stateexact.chebmom
IDENTICAL  runs/goldens/g10_chain1d_real_generalised.json -> NonEqOpVX-VXchain1dKPM_M64x64_stateexact.chebmom2D
IDENTICAL  runs/goldens/g10_chain1d_real_generalised.json -> SpectralOp1chain1dKPM_M64_stateexact.chebmom
IDENTICAL  runs/goldens/g10_chain1d_real_generalised.json -> SpectralOp1chain1dKPM_M64_stateexact.chebmom1D
IDENTICAL  runs/goldens/g10_chain1d_real_generalised_chebyshev.json -> KuboGreenWood_NonEqOpVX-VXchain1dKPM_M64x64_stateexactJACKSON.dat
IDENTICAL  runs/goldens/g10_chain1d_real_generalised_chebyshev.json -> meanSpectralOp1chain1dKPM_M64_stateexactJACKSON.dat
IDENTICAL  runs/goldens/g10_chain1d_real_generalised_chebyshev.json -> NonEqOpVX-VXchain1dKPM_M64x64_stateexact.chebmom
IDENTICAL  runs/goldens/g10_chain1d_real_generalised_chebyshev.json -> NonEqOpVX-VXchain1dKPM_M64x64_stateexact.chebmom2D
IDENTICAL  runs/goldens/g10_chain1d_real_generalised_chebyshev.json -> SpectralOp1chain1dKPM_M64_stateexact.chebmom
IDENTICAL  runs/goldens/g10_chain1d_real_generalised_chebyshev.json -> SpectralOp1chain1dKPM_M64_stateexact.chebmom1D
IDENTICAL  runs/goldens/g10_chain1d_real_generalised_ischo.json -> KuboGreenWood_NonEqOpVX-VXchain1dKPM_M64x64_stateexactJACKSON.dat
IDENTICAL  runs/goldens/g10_chain1d_real_generalised_ischo.json -> meanSpectralOp1chain1dKPM_M64_stateexactJACKSON.dat
IDENTICAL  runs/goldens/g10_chain1d_real_generalised_ischo.json -> NonEqOpVX-VXchain1dKPM_M64x64_stateexact.chebmom
IDENTICAL  runs/goldens/g10_chain1d_real_generalised_ischo.json -> NonEqOpVX-VXchain1dKPM_M64x64_stateexact.chebmom2D
IDENTICAL  runs/goldens/g10_chain1d_real_generalised_ischo.json -> SpectralOp1chain1dKPM_M64_stateexact.chebmom
IDENTICAL  runs/goldens/g10_chain1d_real_generalised_ischo.json -> SpectralOp1chain1dKPM_M64_stateexact.chebmom1D
IDENTICAL  runs/goldens/g10_chain1d_real_lowdin.json -> KuboGreenWood_NonEqOpVX-VXchain1dKPM_M64x64_stateexactJACKSON.dat
IDENTICAL  runs/goldens/g10_chain1d_real_lowdin.json -> meanSpectralOp1chain1dKPM_M64_stateexactJACKSON.dat
IDENTICAL  runs/goldens/g10_chain1d_real_lowdin.json -> NonEqOpVX-VXchain1dKPM_M64x64_stateexact.chebmom
IDENTICAL  runs/goldens/g10_chain1d_real_lowdin.json -> NonEqOpVX-VXchain1dKPM_M64x64_stateexact.chebmom2D
IDENTICAL  runs/goldens/g10_chain1d_real_lowdin.json -> SpectralOp1chain1dKPM_M64_stateexact.chebmom
IDENTICAL  runs/goldens/g10_chain1d_real_lowdin.json -> SpectralOp1chain1dKPM_M64_stateexact.chebmom1D
IDENTICAL  runs/goldens/g10_chain1d_real_lowdin_chebyshev.json -> KuboGreenWood_NonEqOpVX-VXchain1dKPM_M64x64_stateexactJACKSON.dat
IDENTICAL  runs/goldens/g10_chain1d_real_lowdin_chebyshev.json -> meanSpectralOp1chain1dKPM_M64_stateexactJACKSON.dat
IDENTICAL  runs/goldens/g10_chain1d_real_lowdin_chebyshev.json -> NonEqOpVX-VXchain1dKPM_M64x64_stateexact.chebmom
IDENTICAL  runs/goldens/g10_chain1d_real_lowdin_chebyshev.json -> NonEqOpVX-VXchain1dKPM_M64x64_stateexact.chebmom2D
IDENTICAL  runs/goldens/g10_chain1d_real_lowdin_chebyshev.json -> SpectralOp1chain1dKPM_M64_stateexact.chebmom
IDENTICAL  runs/goldens/g10_chain1d_real_lowdin_chebyshev.json -> SpectralOp1chain1dKPM_M64_stateexact.chebmom1D
IDENTICAL  runs/goldens/g10_haldane_complex_generalised.json -> KuboBastin_NonEqOpVX-VYhaldaneKPM_M128x128_stateexactJACKSON.dat
IDENTICAL  runs/goldens/g10_haldane_complex_generalised.json -> KuboBastinSea_NonEqOpVX-VYhaldaneKPM_M128x128_stateexactJACKSON.dat
IDENTICAL  runs/goldens/g10_haldane_complex_generalised.json -> KuboBastinSurf_NonEqOpVX-VYhaldaneKPM_M128x128_stateexactJACKSON.dat
IDENTICAL  runs/goldens/g10_haldane_complex_generalised.json -> KuboGreenWood_NonEqOpVX-VXhaldaneKPM_M128x128_stateexactJACKSON.dat
IDENTICAL  runs/goldens/g10_haldane_complex_generalised.json -> NonEqOpVX-VXhaldaneKPM_M128x128_stateexact.chebmom
IDENTICAL  runs/goldens/g10_haldane_complex_generalised.json -> NonEqOpVX-VXhaldaneKPM_M128x128_stateexact.chebmom2D
IDENTICAL  runs/goldens/g10_haldane_complex_generalised.json -> NonEqOpVX-VYhaldaneKPM_M128x128_stateexact.chebmom
IDENTICAL  runs/goldens/g10_haldane_complex_generalised.json -> NonEqOpVX-VYhaldaneKPM_M128x128_stateexact.chebmom2D
IDENTICAL  runs/goldens/g10_haldane_complex_lowdin.json -> KuboBastin_NonEqOpVX-VYhaldaneKPM_M128x128_stateexactJACKSON.dat
IDENTICAL  runs/goldens/g10_haldane_complex_lowdin.json -> KuboBastinSea_NonEqOpVX-VYhaldaneKPM_M128x128_stateexactJACKSON.dat
IDENTICAL  runs/goldens/g10_haldane_complex_lowdin.json -> KuboBastinSurf_NonEqOpVX-VYhaldaneKPM_M128x128_stateexactJACKSON.dat
IDENTICAL  runs/goldens/g10_haldane_complex_lowdin.json -> KuboGreenWood_NonEqOpVX-VXhaldaneKPM_M128x128_stateexactJACKSON.dat
IDENTICAL  runs/goldens/g10_haldane_complex_lowdin.json -> NonEqOpVX-VXhaldaneKPM_M128x128_stateexact.chebmom
IDENTICAL  runs/goldens/g10_haldane_complex_lowdin.json -> NonEqOpVX-VXhaldaneKPM_M128x128_stateexact.chebmom2D
IDENTICAL  runs/goldens/g10_haldane_complex_lowdin.json -> NonEqOpVX-VYhaldaneKPM_M128x128_stateexact.chebmom
IDENTICAL  runs/goldens/g10_haldane_complex_lowdin.json -> NonEqOpVX-VYhaldaneKPM_M128x128_stateexact.chebmom2D
IDENTICAL  runs/goldens/g10_haldane_real_generalised.json -> KuboBastin_NonEqOpVX-VYhaldaneKPM_M128x128_stateexactJACKSON.dat
IDENTICAL  runs/goldens/g10_haldane_real_generalised.json -> KuboBastinSea_NonEqOpVX-VYhaldaneKPM_M128x128_stateexactJACKSON.dat
IDENTICAL  runs/goldens/g10_haldane_real_generalised.json -> KuboBastinSurf_NonEqOpVX-VYhaldaneKPM_M128x128_stateexactJACKSON.dat
IDENTICAL  runs/goldens/g10_haldane_real_generalised.json -> KuboGreenWood_NonEqOpVX-VXhaldaneKPM_M128x128_stateexactJACKSON.dat
IDENTICAL  runs/goldens/g10_haldane_real_generalised.json -> NonEqOpVX-VXhaldaneKPM_M128x128_stateexact.chebmom
IDENTICAL  runs/goldens/g10_haldane_real_generalised.json -> NonEqOpVX-VXhaldaneKPM_M128x128_stateexact.chebmom2D
IDENTICAL  runs/goldens/g10_haldane_real_generalised.json -> NonEqOpVX-VYhaldaneKPM_M128x128_stateexact.chebmom
IDENTICAL  runs/goldens/g10_haldane_real_generalised.json -> NonEqOpVX-VYhaldaneKPM_M128x128_stateexact.chebmom2D
IDENTICAL  runs/goldens/g10_haldane_real_lowdin.json -> KuboBastin_NonEqOpVX-VYhaldaneKPM_M128x128_stateexactJACKSON.dat
IDENTICAL  runs/goldens/g10_haldane_real_lowdin.json -> KuboBastinSea_NonEqOpVX-VYhaldaneKPM_M128x128_stateexactJACKSON.dat
IDENTICAL  runs/goldens/g10_haldane_real_lowdin.json -> KuboBastinSurf_NonEqOpVX-VYhaldaneKPM_M128x128_stateexactJACKSON.dat
IDENTICAL  runs/goldens/g10_haldane_real_lowdin.json -> KuboGreenWood_NonEqOpVX-VXhaldaneKPM_M128x128_stateexactJACKSON.dat
IDENTICAL  runs/goldens/g10_haldane_real_lowdin.json -> NonEqOpVX-VXhaldaneKPM_M128x128_stateexact.chebmom
IDENTICAL  runs/goldens/g10_haldane_real_lowdin.json -> NonEqOpVX-VXhaldaneKPM_M128x128_stateexact.chebmom2D
IDENTICAL  runs/goldens/g10_haldane_real_lowdin.json -> NonEqOpVX-VYhaldaneKPM_M128x128_stateexact.chebmom
IDENTICAL  runs/goldens/g10_haldane_real_lowdin.json -> NonEqOpVX-VYhaldaneKPM_M128x128_stateexact.chebmom2D

== TOTAL: 12 G9 + 11 G10 run files, 170 files compared with cmp: 0 failures ==
```

---

## Review 3 applied

The four items of `handoff/0007-words-cli-mcp/review-3.md` were applied in order on `handoff/0007-agi` after `5e8186d`, one commit per item:

1. **`1b91cd7` — Tutorial 06 honours the run file or refuses it (`python/lsquant/apps/tutorial06.py`, `python/lsquant/runfile.py`)**: `_eval_and_reconstruct_pair` now compiles the run file's `words` directly (`pz_words` and `px_words` partitioned from `rf.words`) and derives moment/curve filename stems from `compiled.reconstruction.legacy_moment_prefix`, while forwarding `expansion` (`moments`, `window`, `trace`, `states`, `seed`), `time` (`steps`, `tmax`/`dt`, `bessel`), `reconstruction` (`kernel`, `broadening_meV`, `fermi_energies`, `formulas`, `legacy_text_round_trip`), and `output.capture_tutorial` from the run file (and refusing `window: "eigen"` on `model.kind == "generator"`).
2. **`bafbb84` — Remove `ensure_w2s_csr_dir` and refuse missing or mismatched `csr_dir` before any load (`python/lsquant/bases.py`, `python/lsquant/runfile.py`)**: deleted `find_w2s_binary` and `ensure_w2s_csr_dir` and replaced them with `validate_w2s_csr_dir` (called in `execute_runfile` before `os.makedirs(out_dir)` and before `load_wannier_model`, as well as in `RealBasis.__init__`), which refuses when `<seed>.HAM.CSR` or `<seed>.out` is missing or disagrees with the requested supercell/seed/operators with a message naming `csr_dir`, the mismatch, and the exact `wannier2sparse <seed> <N1> <N2> <N3> VX -p <model_dir> -o <csr_dir>` command.
3. **`e281749` — Resolve `LSQUANT_BUILD_DIR` (default `build/`) in `tests/dev/stage_7_check.py`**: replaced all hard-coded `os.path.join(REPO, "build", ...)` paths (including `stage6_nonorth` lookup at lines 370–372 and `sys.path`/`PYTHONPATH` setup) with `BUILD_DIR` resolved from `LSQUANT_BUILD_DIR` (default `build/`).
4. **`b84d763` — Remove `legacy/` and all `inline_*` driver references in `docs/` and `examples/`**: deleted `legacy/` (19 tracked legacy driver sources, preserved by the `legacy-final` tag), rewrote `examples/06_spin_relaxation_square/lsqrelax.py` and `examples/06_spin_relaxation_square/rp_spin.py` (along with the other tutorial scripts) to execute via `lsquant.runfile.execute_runfile` modelled on `runs/goldens/tutorial06.json`, and removed all `inline_*` driver mentions across `docs/` and `examples/` (`grep -rn "inline_" docs/ examples/` is empty).

### `tests/dev/stage_7_check.py` output on TCN (`OMP_NUM_THREADS=16 LSQUANT_BUILD_DIR=build`)

```text
== (a) Catalogue canonical parse/print round-trip ==
  PASS  Tr[delta(E-H)]                             -> shape=dos_identity         plan_builder=dos_plan
  PASS  Tr[A delta(E-H)]                           -> shape=dos_operator         plan_builder=dos_plan
  PASS  Tr[VX delta(E-H) VX delta(E-H)]            -> shape=kubo_greenwood       plan_builder=kubo_plan
  PASS  Tr[VX delta(E-H) VY dG-(E)/dE]             -> shape=kubo_bastin          plan_builder=kubo_plan
  PASS  Tr[[X,U(t)]^+ delta(E-H) [X,U(t)]]         -> shape=msd                  plan_builder=msd_plan
  PASS  Tr[P U(t)^+ SX U(t) delta(E-H) P]          -> shape=projected_evolution  plan_builder=projected_evolution_plan
  PASS  Tr[A delta(E-H) U(t)^+ B U(t)]             -> shape=time_correlation     plan_builder=time_correlation_plan
  PASS  Tr[sea[VX delta(E-H) VY dG-(E)/dE]]        -> shape=kubo_bastin          plan_builder=kubo_plan
  PASS  Tr[P:PX U(t)^+ SX U(t) delta(E-H) P:PX]    -> shape=projected_evolution  plan_builder=projected_evolution_plan
  PASS  Tr[B^+ VX@bloch B delta(E-H)]              -> shape=dos_operator         plan_builder=dos_plan

== (b) Refusal of 5 malformed / unsupported inputs ==
  PASS  [unknown_operator] refused with: Unknown operator 'FOO' in word 'Tr[FOO delta(E-H)]' (registered operators: ['1', 'VX', 'VY'])
  PASS  [rank_3] refused with: Unsupported word 'Tr[VX delta(E-H) VY delta(E-H) VZ delta(E-H)]' of rank 3 (3 expanded spectral slots): rank >= 3 moment tensors are out of scope; nearest recognised shape is 'Tr[VX delta(E-H) VX delta(E-H)]' (kubo_greenwood)
  PASS  [mixed_basis_without_transform] refused with: Mixed basis without explicit transform: operator 'VX@bloch' is in 'bloch' basis while calculation/active basis is 'real' (insert explicit B^+ ... B or U^+ ... U transform)
  PASS  [ut_without_time] refused with: Word 'Tr[A delta(E-H) U(t)^+ B U(t)]' contains time evolution U(t) without 'time' parameters (steps and tmax required)
  PASS  [unknown_key_in_runfile] refused with: Unknown key 'unknown_extra_key' in run file (top-level); allowed keys are ['basis', 'expansion', 'metric', 'model', 'output', 'reconstruction', 'time', 'word', 'words']
  PASS  [unsupported_green_word] refused with: Unrecognised word shape 'Tr[VX delta(E-H) VY G+(E)]'; nearest recognised shape is 'Tr[VX delta(E-H) VY dG-(E)/dE]' (kubo_bastin)
  PASS  [eigen_window_on_csr] refused with: expansion.window='eigen' is not supported for model.kind='csr' (no k blocks; use gershgorin or an explicit window)
  PASS  [unmatched_formula_for_word] refused with: Reconstruction formula 'greenwood' in reconstruction.formulas does not match any word in ['Tr[VX delta(E-H) VY dG-(E)/dE]'] (allowed formulas for given word(s): ['bastin', 'bastin_fft', 'kubo_bastin', 'kubo_bastin_fft', 'kubo_bastin_sea', 'kubo_bastin_surf', 'sea', 'surf'])
  PASS  [unknown_generator_param] refused with: Unknown key 'bogus_param' in run file (model.params (make_chain)); allowed keys are ['hopping_t', 'sizes']
  PASS  [unknown_generator] refused with: Unknown generator name 'evil_script' in run file (model.name); allowed generators are ['gen_models', 'make_chain', 'make_chain_mag', 'make_disordered_chain', 'make_ey', 'make_haldane', 'make_rashba', 'make_square']
  PASS  [window_default_by_basis] real->gershgorin, bloch->eigen, eigen->eigen
  PASS  [runfile_schema] all 40 files in runs/goldens/*.json validated

== (c) CLI reproduction of G1, G5, G9 (g9_graphene_eigen.json), and G10 (g10_chain1d_real_generalised.json) ==
  PASS  cmp IDENTICAL: G1 (SpectralOp1chain1dKPM_M64_stateexact.chebmom1D)
  PASS  cmp IDENTICAL: G5 (NonEqOpVX-VYhaldaneKPM_M128x128_stateexact.chebmom2D)
  PASS  cmp IDENTICAL: G5 (NonEqOpVX-VXhaldaneKPM_M128x128_stateexact.chebmom2D)
  PASS  cmp IDENTICAL: G5 (KuboBastin_NonEqOpVX-VYhaldaneKPM_M128x128_stateexactJACKSON.dat)
  PASS  cmp IDENTICAL: G5 (KuboBastinSea_NonEqOpVX-VYhaldaneKPM_M128x128_stateexactJACKSON.dat)
  PASS  cmp IDENTICAL: G5 (KuboBastinSurf_NonEqOpVX-VYhaldaneKPM_M128x128_stateexactJACKSON.dat)
  PASS  cmp IDENTICAL: G5 (KuboGreenWood_NonEqOpVX-VXhaldaneKPM_M128x128_stateexactJACKSON.dat)
  PASS  cmp IDENTICAL: G9 g9_graphene_eigen (KuboGreenWood_NonEqOpVX-VXgraphene_eigen_KPM_M256x256_stateexactJACKSON.dat)
  PASS  cmp IDENTICAL: G9 g9_graphene_eigen (NonEqOpVX-VXgraphene_eigen_KPM_M256x256_stateexact.basis.json)
  PASS  cmp IDENTICAL: G9 g9_graphene_eigen (NonEqOpVX-VXgraphene_eigen_KPM_M256x256_stateexact.chebmom)
  PASS  cmp IDENTICAL: G9 g9_graphene_eigen (NonEqOpVX-VXgraphene_eigen_KPM_M256x256_stateexact.chebmom2D)
  PASS  cmp IDENTICAL: G9 g9_graphene_eigen (SpectralOp1graphene_eigen_KPM_M256_stateexact.basis.json)
  PASS  cmp IDENTICAL: G9 g9_graphene_eigen (SpectralOp1graphene_eigen_KPM_M256_stateexact.chebmom)
  PASS  cmp IDENTICAL: G9 g9_graphene_eigen (SpectralOp1graphene_eigen_KPM_M256_stateexact.chebmom1D)
  PASS  cmp IDENTICAL: G9 g9_graphene_eigen (meanSpectralOp1graphene_eigen_KPM_M256_stateexactJACKSON.dat)
  PASS  cmp IDENTICAL: G10 g10_chain1d_real_generalised (KuboGreenWood_NonEqOpVX-VXchain1dKPM_M64x64_stateexactJACKSON.dat)
  PASS  cmp IDENTICAL: G10 g10_chain1d_real_generalised (NonEqOpVX-VXchain1dKPM_M64x64_stateexact.chebmom)
  PASS  cmp IDENTICAL: G10 g10_chain1d_real_generalised (NonEqOpVX-VXchain1dKPM_M64x64_stateexact.chebmom2D)
  PASS  cmp IDENTICAL: G10 g10_chain1d_real_generalised (SpectralOp1chain1dKPM_M64_stateexact.chebmom)
  PASS  cmp IDENTICAL: G10 g10_chain1d_real_generalised (SpectralOp1chain1dKPM_M64_stateexact.chebmom1D)
  PASS  cmp IDENTICAL: G10 g10_chain1d_real_generalised (meanSpectralOp1chain1dKPM_M64_stateexactJACKSON.dat)

== (d) MCP server subprocess vs CLI JSON equivalence ==
  MCP tools registered: ['words_describe', 'run', 'status', 'reconstruct', 'inspect']
  PASS  MCP words_describe JSON == CLI words describe JSON
  PASS  MCP inspect JSON == CLI inspect JSON
  PASS  MCP error responses preserve reason: 'Error executing tool words_describe: Unrecognised word shape 'Tr[VX delta(E-H) VY G+(E)]'; nearest recognised shape is 'Tr[VX delta(E-H) VY dG-(E)/dE]' (kubo_bastin)'

STAGE 7 SELF-CHECK PASS
```

### Tutorial 06 (`runs/goldens/tutorial06.json`) verification and `cmp` status

- `git fetch origin handoff/0007-codex` was run: `origin/handoff/0007-codex` is at `22f92a8` (`handoff 0007 tester: closure run on 5e8186d (deleted legacy tree), FIG. 7 with tutorial 06`), and the tester has not yet committed/pushed the re-blessed `goldens/stage0/tutorials/06` commit to `origin/handoff/0007-codex`.
- File counts produced by `runs/goldens/tutorial06.json`: **281** golden-compared files (**13** top-level `.dat` curves + **255** `captured/0000_*.dat`..`0254_*.dat` files + **1** `captured/index.tsv` + **12** `.chebmomTD` moment files, plus `run.resolved.json`).
- After Item 1 (`1b91cd7`), executing `tutorial06` via `execute_runfile` produces `meanTimeProjSZ-PZey_L8_D0.1KPM_M96_statedefaultEF2.000000JACKSON.dat` and `meanTimeProj1-PZey_L8_D0.1KPM_M96_statedefaultEF2.000000JACKSON.dat` byte-identical (`cmp` PASS) to the tester's `~/s7keepc/T16/tutorial06` output on TCN; the full 281-file `cmp` against `goldens/stage0/tutorials/06` on `origin/handoff/0007-codex` is pending the tester's re-blessing commit.

