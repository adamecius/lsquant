# Unified I/O & results system — plan

**Status:** proposal / research (rev 2, reconciled with the sanitized codebase).
No code changed by this document. Companion to
[`reporting_system.md`](reporting_system.md) — timing, error, and resource data
that land in the unified output are *produced* by the reporting layer; this plan
defines the **input schema**, **system identity / merge**, **unified output
container**, and **presentation (Python + HTML dashboard)**. Everything is
**additive** (new schema, new files, opt-in), so it can land beside in-flight work.

## 0. Reconciliation with the sanitize effort (rev 2)

Sanitize Phases A–D changed the ground this plan stands on:

- **Data ban (Phase C — `.gitignore` + `.githooks/pre-commit`).** Generated
  transport artifacts (`*.CSR`, `*.chebmom*`) are a **regenerable cache that must
  never be committed**; a pre-commit hook rejects >200 KB files outside
  `test/golden/`|`test/models/`. **This reclassifies the `<system>.results.json`
  container: it is *output*, not a source artifact.** Rev-1 sold "git-diffable"
  as a virtue — that is now *wrong*. Corrected in §4.4: `results.json` is
  **gitignored** alongside the other outputs; "human-readable" means *easy to
  read/transmit/inspect*, not *commit it to the repo*. Small thin-grid spectra
  also keep it under the 200 KB hook limit if one ever is staged deliberately.
- **Supported surface = `lsquant <command>` (Phase D).** `compute`,
  `reconstruct`, `inspect` are *the* interface; 18 drivers were quarantined to
  `legacy/` (read-only). This **reinforces** the plan's direction: the unified
  entry point *is* `lsquant compute` (extended schema), and merge becomes a new
  first-class subcommand **`lsquant merge`** — not a side script.
- **Pruned docs.** `ROADMAP.md` / `TEST_PLAN.md` are gone; references re-anchored
  to `README.md` and [`legacy/README.md`](legacy/README.md).
- **"Goldens first, always."** Each phase below ships with golden/test coverage
  before it is wired into the build (see §6).
- **Spin relaxation is now real** (Tutorial 6, `examples/06_spin_relaxation_square`,
  driven through `lsquant compute msd`). The `spin_relaxation` observable in the
  schema (§4.1) maps onto that existing `msd`/time-correlation path, not a new
  engine.
- **The embryonic merge scripts are partly legacy now.** The C++ `*FromChebmom`
  drivers split: `kuboBastin`/`greenwood`/`spectral` stay active in `src/`; the
  I/II/optical/localDensities/FFTgrid-spectral variants moved to `legacy/`. The
  Python `add/average/substract_chebmom2D.py` remain under `utilities/python` and
  are still the header-string-match merge to be superseded.

## 1. Goal (the vision, restated)

1. **One unified entry point.** For a *single system*, submit *one* config that
   requests several observables at once — e.g. DOS **and** conductivity **and**
   spin relaxation — and get *one* unified output.
2. **Merge of separate runs.** If two runs target the *same, unchanged* system
   (e.g. a DOS run, then later a conductivity run), I can **merge** them into a
   single results file. The merger proves "same system" by comparing a
   **fingerprint hash**, not by trusting filenames.
3. **Validated input with intuitive hints.** A malformed/incomplete config tells
   me *exactly what is missing* to define a simulation, in plain language.
4. **Unified output that is both human- and machine-readable.** Final spectra
   (σ, DOS, …) stored on a **thin presentation grid**, plus **timing**, **errors/
   warnings**, and **provenance** — all in one machine-friendly file.
5. **Use mature OSS** for the two steps (config-in, results-out) where it exists;
   build only the domain-specific glue (the physics merge) ourselves.
6. **HTML dashboard** generated from the output (self-contained file now; online
   site later).
7. **A Python presenter** in `utilities/python` that ingests the output and shows
   the results nicely (tables + plots + the HTML).

## 2. What exists today (audit)

| Step | Current state | Gap vs. vision |
|---|---|---|
| Input | Flat `run.json`, **one mode per run** (`noneq`/`spectral`/`msd`), parsed by a hand-rolled dependency-free reader (`run_config.cpp`). Validation returns terse strings (`"missing 'label'"`). | No multi-observable run; no schema; hints are minimal. |
| System identity | A system = `operators/<label>.<ID>.CSR` + `<label>.<ID>.desc` descriptors (observable/units/provenance/bounds/basis). **No hash anywhere.** | Nothing proves two runs share the *same* system. |
| Output (moments) | `.chebmom1D/2D`: 2-line plaintext header (`size bandwidth bandcenter` / `num_moments`) + data (`chebyshev_moments1D.cpp:112`, `..2D.cpp:53`). | Per-step, ad-hoc; no timing/errors/provenance inside. |
| Output (spectra) | `reconstruct` writes `mean*JACKSON.dat` / `Kubo*JACKSON.dat`, two-column plaintext on the *full* recon grid (`lsquant.cpp`). | Not unified; full grid (heavy); no metadata. |
| Merge | **Already in embryo**: `utilities/python/bin/{add,average,substract}_chebmom2D.py` merge by requiring the **header string to be byte-identical** (`add_chebmom2D.py`). | This is exactly the "same system?" check — but brittle (string match, not content hash) and only for one file type. |
| Inspect | `lsquant inspect` resolves numerical *provenance* (alpha, bounds source, hbar, descriptor) — `inspect.cpp`. | Doesn't hash, time, or emit machine-readable output. |
| Python utils | `utilities/python/` has `sciplots.py` + the merge scripts, but carries **committed `build/` dirs** and py3.6 `__pycache__` cruft. | Stale; needs a clean presenter entry point. |

Conventions to honour (same as the reporting plan): **C++11**; dependencies via
`find_package` first, then a pinned `FetchContent` (offline-graceful like
`cmake/FetchWannier2Sparse.cmake`), then a vendored single header; keep the
**offline-safe** posture the hand-rolled JSON reader was written for.

## 3. Architecture overview

```
            unified config (run.json / run.toml)
                       │
              ┌────────▼─────────┐   schema validation → intuitive hints
              │  parse + validate │   (mature OSS: nlohmann/json + JSON-Schema)
              └────────┬─────────┘
                       │  expands to N observable jobs sharing ONE system
              ┌────────▼─────────┐
              │  system fingerprint│  hash(CSR ops + descriptors + numerics)
              └────────┬─────────┘   → stamped into every output
                       │
        ┌──────────────┼───────────────┐
     compute DOS   compute σ        compute spin-relax     (each timed; errors captured)
        └──────────────┼───────────────┘
                       │  thin-grid final spectra + timing + errors + provenance
              ┌────────▼─────────┐
              │ unified results   │  <system>.results.json   (merge-by-fingerprint)
              │   container       │
              └────────┬─────────┘
              ┌────────┴───────────────┐
        utilities/python              static HTML dashboard
        lsquant_report.py  ──────────▶  <system>.report.html   (→ online site later)
```

## 4. Component plans

### 4.1 Input: unified config + schema validation — recommend **nlohmann/json + JSON-Schema**

The hand-rolled flat parser was a deliberate offline-safe choice, but it cannot
express "a system + a list of observables," and its validation can't produce rich
hints. Adopt the de-facto mature stack, kept offline-safe by vendoring:

- **`nlohmann/json`** — the standard C++ JSON library (header-only, MIT, single
  amalgamated header, vendorable → zero network, zero link). Replaces
  `parse_flat_json` *behind the existing `read_run_config` API* so callers don't
  change — additive.
- **JSON Schema** (validator: `pboettch/json-schema-validator`, header-ish, MIT)
  — declare the run schema once; the validator yields messages like
  *"`/observables/1`: required property `operator_left` is missing"*. Wrap these
  into **intuitive hints** keyed to physics (see 4.2).
- **Schema v2 (the unified entry point)** — extend the flat v1 into:
  ```jsonc
  {
    "system": {                      // defines identity (→ fingerprint)
      "label": "graphene_50nm",
      "operators_dir": "operators",  // <label>.HAM.CSR, <label>.VX.CSR, ...
      "bounds": "auto"               // auto | descriptor | [a,b]
    },
    "numerics": { "num_moments": 512, "kernel": "jackson", "broadening_meV": 10 },
    "observables": [                 // request several at once, same system
      { "kind": "dos" },
      { "kind": "conductivity", "component": "xx" },
      { "kind": "conductivity", "component": "xy" },
      { "kind": "spin_relaxation", "spin_axis": "z" }
    ],
    "output": { "results": "graphene_50nm.results.json", "grid_points": 600 }
  }
  ```
  Each `observables[]` entry expands to one existing `compute_*` call (`noneq` /
  `spectral` / `msd`), reusing today's dispatch in `lsquant.cpp` — the unified
  config is a **thin front-end** over the existing engine, not a rewrite.
- **Backward compatible:** a v1 flat `run.json` (single mode) is accepted as "a
  system with one observable." Detect by presence of `observables` vs `mode`.
- **Format choice:** keep **JSON** as canonical (continuity with `run.json`, and
  it doubles as the *output* format and the web payload). Optionally also accept
  **TOML** (via `toml++`/`toml11`) since physicists often prefer it for hand
  editing — both deserialize to the same internal struct. *Decision point E1.*

### 4.2 Validation with intuitive, physics-aware hints

JSON-Schema gives the *structural* "what's missing"; layer a thin **semantic**
check that speaks the domain, e.g.:

- `conductivity` with `component: "xy"` but only one velocity operator present →
  *"Hall conductivity σ_xy needs two velocity operators; found VX but not VY in
  `operators/`. Provide graphene_50nm.VY.CSR or set component to xx."*
- `spin_relaxation` but no spin descriptor on the Hamiltonian → *"spin relaxation
  needs a spin-resolved basis; `graphene_50nm.HAM.desc` has no `basis`/spin
  block."*
- `num_moments` missing → *"num_moments not set: it controls energy resolution
  (≈ bandwidth / num_moments). Typical 256–2048."*

This is a small rules table (`include/io/validate.hpp`), additive, returning a
list of `Hint{severity, where, message, suggestion}`. Severity `error` blocks the
run; `warning` is logged via the reporting layer and embedded in the output.

### 4.3 System identity / fingerprint — recommend **SHA-256 (picosha2)**

A run can be merged with another **iff they share the same physical system**.
Define a **system fingerprint** = a stable hash over exactly the inputs that
change the physics:

```
fingerprint = SHA256(
    content(HAM.CSR) ⊕ content(each operator .CSR used) ⊕
    content(relevant .desc files) ⊕
    canonical(numerics that affect moments: bounds, scaling, rescaling) )
```

- **Hash lib:** **`picosha2`** — single-header, public-domain, deterministic,
  hex output that's human-quotable (`fp: 9f3a…`). Identity, not speed, matters
  here. If hashing multi-GB CSR files becomes a bottleneck, swap to **xxHash**
  (header-only, BSD, ~GB/s) — note as *decision point E2*; the interface
  (`io::fingerprint(system)`) hides the choice.
- The fingerprint is **stamped into every output** (`.results.json` and, as a new
  optional header comment, into `.chebmom` files — additive, ignored by old
  readers). It **supersedes the brittle header-string match** the python merge
  scripts use today.
- `lsquant inspect` gains a line: `system fingerprint: 9f3a…` so a human can eyeball
  whether two runs match before merging.
- **Note:** the fingerprint deliberately **excludes** `num_moments`,
  `broadening`, and which *observable* — those differ between mergeable runs (a
  64-moment DOS and a 512-moment σ_xy on the same system still merge into one
  results file; each observable keeps its own moments/resolution metadata).

### 4.4 Unified output container — recommend **JSON for results, keep `.chebmom` for raw moments**

Two tiers, by data size:

- **Raw moments** (large): stay in the existing `.chebmom1D/2D` files (optionally
  gain an HDF5 backend later for big systems — *decision point E3*). Add the
  fingerprint + a provenance comment to their header (additive).
- **Unified results** (`<system>.results.json`) — the machine+human-readable
  container the user described. One file per system, **merge-by-fingerprint**:
  ```jsonc
  {
    "schema": "lsquant.results/1",
    "system": { "label": "...", "fingerprint": "9f3a…", "size": 2048,
                "bandwidth": 4.0, "bandcenter": 0.0 },
    "results": {
      "dos":            { "energy_eV": [...], "dos": [...],        // THIN grid
                          "num_moments": 512, "kernel": "jackson",
                          "timing_s": 31.2, "warnings": [] },
      "sigma_xx":       { "energy_eV": [...], "sigma": [...], "units": "e^2/h",
                          "num_moments": 512, "timing_s": 88.0, "warnings": [] },
      "sigma_xy":       { ... },
      "spin_relax_z":   { "time_fs": [...], "Sz": [...], "tau_fs": 412.0,
                          "timing_s": 54.7, "warnings": [] }
    },
    "run": { "host": "...", "started": "...", "wall_s": 174,
             "peak_rss_bytes": 4509715660,        // from reporting_system.md
             "errors": [] }
  }
  ```
- **Thin grid:** `output.grid_points` (default ~600) downsamples the final spectra
  for storage/plotting; the dense recon grid stays available on demand. Keeps the
  unified file small (web-friendly, well under the 200 KB pre-commit limit) — but
  note it is still a **regenerable output**: add `*.results.json` to `.gitignore`
  (Phase C data-ban philosophy); it is produced, transmitted, and rendered, never
  committed (except, like goldens, a deliberately curated fixture under
  `test/golden/`).
- **Machine-friendly + human-readable:** pretty-printed JSON satisfies both; it
  carries **timing** and **errors/warnings** inline (sourced from the reporting
  layer), and it **merges trivially** — merge = load both, verify
  `system.fingerprint` equal, deep-merge `results`, concatenate `run.errors`. If
  fingerprints differ, refuse with a clear diff ("HAM.CSR hash differs").
- **`lsquant merge a.results.json b.results.json -o out.results.json`** — a new,
  small C++ subcommand (and/or python), the matured successor to the
  header-string python scripts. This is the "merge function" with the hash check.

> "Use mature OSS if it exists." For the *container + validation* steps, the
> mature pieces are **nlohmann/json** (parse/emit) and **JSON-Schema** (validate).
> No off-the-shelf tool knows how to merge *KPM results by physical system* — that
> domain glue (fingerprint + deep-merge) is the only part we build, and it's thin.

### 4.5 Presentation: Python tool in `utilities/python` + HTML dashboard

- **`utilities/python/lsquant_report.py`** (new) — ingests one or many
  `*.results.json`, and:
  1. prints a summary table (system, fingerprint, observables present, timing,
     peak RSS, any errors);
  2. plots DOS / σ_xx / σ_xy / spin-relaxation (reuse `lib/sciplots.py`);
  3. emits a **self-contained HTML dashboard** `<system>.report.html`.
- **HTML dashboard** — recommend a **single self-contained HTML file** via
  **Plotly** (offline mode embeds the JS) or a Jinja2 template — interactive plots,
  the metadata/timing/error panel, opens in any browser, nothing to host. This
  directly satisfies "report everything in HTML." *Decision point E4: Plotly
  (interactive, ~3 MB JS embedded) vs. Matplotlib-PNG-in-HTML (tiny, static).*
- **Online site (later).** The same `results.json` is the API payload, so a future
  dashboard (static site reading a directory of results, or a small FastAPI/Flask
  service) is an additive consumer — no format change. Out of scope for v1.
- **Cleanup (flagged, not bundled):** `utilities/python` has committed `build/`
  and `__pycache__/` directories and py3.6 artifacts. Recommend a separate
  housekeeping pass (`.gitignore` them, target Python ≥3.8); do **not** mix it
  into this feature to keep the change auditable.

## 5. Build & dependency integration (CMake)

Same ladder as the reporting plan — system package → pinned `FetchContent` →
vendored header — so the project stays offline-safe.

- `nlohmann/json`: `find_package(nlohmann_json)` (Ubuntu `nlohmann-json3-dev`) →
  fetch → vendored single header. Define `LSQ_WITH_JSON`.
- `json-schema-validator`: find → fetch → vendor; gated `LSQUANT_WITH_JSON_SCHEMA`
  (default ON; if absent, fall back to the semantic-only validator in 4.2).
- `picosha2`: vendored single header (no option needed).
- New sources under `src/io/` (`config.cpp`, `validate.cpp`, `fingerprint.cpp`,
  `results.cpp`, `merge.cpp`) folded into `KPM_LIB_SOURCES` so every driver gets
  them; new public headers under `include/io/`.
- Python presenter: optional `LSQUANT_BUILD_UTILITIES` already exists — hook
  `lsquant_report.py` there; add a `requirements.txt` (numpy, plotly).

## 6. Rollout (additive, auditable)

Per *"Goldens first, always"*: each phase lands with test coverage before it is
wired in — schema/validation get unit tests on good/bad configs; the fingerprint
gets a stability test (same inputs → same hash, perturbed CSR → different hash);
`results.json` and `merge` get a golden fixture under `test/golden/`.

- **Phase 0:** vendor nlohmann/json; swap `parse_flat_json` internals behind the
  unchanged `read_run_config` API. No behavior change (v1 configs still work).
- **Phase 1:** JSON-Schema + semantic hints; `lsquant inspect` prints hints.
- **Phase 2:** schema v2 (`system`/`numerics`/`observables`); `lsquant compute`
  expands a multi-observable config into the existing `compute_*` calls.
- **Phase 3:** `io::fingerprint`; stamp it into outputs; show in `inspect`.
- **Phase 4:** `<system>.results.json` writer; populate timing/errors from the
  reporting layer; thin-grid downsampling.
- **Phase 5:** `lsquant merge` (fingerprint-checked) — supersede the python merge
  scripts (keep them until parity is proven).
- **Phase 6:** `utilities/python/lsquant_report.py` + static HTML dashboard.
- **Phase 7 (later):** online site; optional HDF5 moment backend.

## 7. Open decisions (need owner confirmation)

| # | Decision | Recommendation |
|---|---|---|
| E1 | Config format: JSON only vs also accept TOML | **JSON canonical**, TOML optional later |
| E2 | Fingerprint hash: SHA-256 (picosha2) vs xxHash | **SHA-256** (identity); xxHash if too slow |
| E3 | Raw moments: keep `.chebmom` vs add HDF5 | **Keep `.chebmom`**, HDF5 later if needed |
| E4 | HTML dashboard: Plotly (interactive) vs Matplotlib PNG | **Plotly** self-contained |
| E5 | Merge tool: C++ `lsquant merge` vs Python | **C++ subcommand** (reuses fingerprint), python wrapper for convenience |
| E6 | utilities/python cleanup (build/ cruft, py3.6) | separate housekeeping pass, not in this feature |

## 8. Relationship to `reporting_system.md`

This plan **consumes** the reporting layer: every `timing_s`, `wall_s`,
`peak_rss_bytes`, and `errors[]`/`warnings[]` field in `results.json` is produced
by `lsquant::util::Timer`, `lsquant::diag::ResourceProbe`, and the
`lsquant::Error` hierarchy defined there. Land the reporting layer's Phases 1–4
first (or in parallel); this plan's Phase 4 is where they surface to the user.

## 9. Out of scope

Hosted multi-user dashboard/auth, a job queue/scheduler, distributed (MPI) result
aggregation, and GPU outputs — infra deferred for now, consistent with the
current `lsquant`-only supported surface (`legacy/README.md`). The `results.json`
schema is versioned (`"schema": "lsquant.results/1"`) so these can extend it
without breaking consumers.
