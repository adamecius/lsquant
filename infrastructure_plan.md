# lsquant infrastructure plan — master

**Status:** proposal / research (rev 2, reconciled with the sanitized codebase).
No code changed by this document. This is the umbrella that ties together two
detailed plans; read them for depth:

- [`reporting_system.md`](reporting_system.md) — logging, timing/profiling, error
  handling, and memory/CPU resource accounting.
- [`unified_io.md`](unified_io.md) — unified config + validation, system
  fingerprint/merge, unified results container, and the Python/HTML presentation.

Everything is **additive** (new `util/` and `io/` modules, opt-in CMake options,
backward-compatible schema) so it can land beside other in-flight work without
disturbing the existing numerics. Targets the **final architecture** (the single
`lsquant` library + driver), not the legacy per-executable layout.

### What changed in rev 2 (Sanitize Phases A–D)

- **Supported surface narrowed** to `lsquant compute|reconstruct|inspect` + ~13
  oracle drivers; 18 drivers quarantined to read-only [`legacy/`](legacy/README.md)
  (`# [LEGACY]` in `src/CMakeLists.txt`). The two new modules wire into
  `kpm_lib` + the `lsquant` driver only; **`legacy/` is never touched.**
- **Data ban** (`.gitignore` + `.githooks/pre-commit`, 200 KB gate): generated
  artifacts are a regenerable cache. **Log files and `results.json` are outputs →
  gitignored, never committed** (this corrects rev-1's "git-diffable results"
  framing).
- **"Goldens first, always"**: every phase below now lands with test coverage
  before it is wired into the build.
- **Pruned docs**: `ROADMAP.md`/`TEST_PLAN.md` removed; references re-anchored to
  `README.md` and `legacy/README.md`.
- **Non-orth throwing placeholder** (preserved at `5cc0924`) is the first client
  of the `lsquant::Error` hierarchy.

## 1. The vision in one picture

```
   unified config ─▶ validate (hints) ─▶ fingerprint ─▶ compute N observables ─▶ unified results.json ─▶ HTML / Python
        │                 │                  │            │  (each timed,            │  (timing, peak RSS,      dashboard
        │                 │                  │            │   errors captured)       │   errors, thin spectra)
        └─────────────────┴──────────────────┴───────────┴── reporting layer: log · timer · error · resource ──┘
```

One config defines a **system** and asks for several **observables** (DOS,
conductivity, spin relaxation). The run is validated with physics-aware hints,
fingerprinted, and executed; the reporting layer instruments every phase; results
land in one machine- and human-readable file that **merges with prior runs on the
same system** and renders to an HTML dashboard.

## 2. Two layers, one stack

| Layer | Module | Concerns | Detail doc |
|---|---|---|---|
| **Reporting** (produces) | `lsquant::util` / `diag` — `include/util/`, `src/util/` | log levels+sinks, scoped timers, `Error` hierarchy + `Result`, Linux-first resource probe behind a backend interface | `reporting_system.md` |
| **Unified I/O** (consumes) | `lsquant::io` — `include/io/`, `src/io/` | schema v2 (system + observables), JSON-Schema validation + hints, SHA-256 system fingerprint, `results.json` container, `lsquant merge`, Python presenter + HTML | `unified_io.md` |

The I/O layer surfaces what the reporting layer produces: every `timing_s`,
`wall_s`, `peak_rss_bytes`, and `errors[]`/`warnings[]` field in `results.json`
comes from `util::Timer`, `diag::ResourceProbe`, and the `lsquant::Error`
hierarchy. **Build reporting first (or in parallel); I/O depends on it.**

## 3. Where it goes

`include/` is currently flat; `src/` is flat with a `Device/` subdir. Introduce
two new modules (additive, the start of proper modularization for the target
architecture):

```
include/util/   log.hpp timer.hpp error.hpp resource.hpp report.hpp
include/io/     config.hpp validate.hpp fingerprint.hpp results.hpp merge.hpp
src/util/       log.cpp timer.cpp error.cpp resource.cpp resource_linux.cpp
src/io/         config.cpp validate.cpp fingerprint.cpp results.cpp merge.cpp
utilities/python/  lsquant_report.py   (+ reuse lib/sciplots.py)
```

Both `src/util/*.cpp` and `src/io/*.cpp` fold into `KPM_LIB_SOURCES`
(`src/CMakeLists.txt`) so every driver linking `kpm_lib` inherits them with no
per-target wiring.

## 4. Dependencies (mature OSS, offline-safe)

All resolved with the project's existing ladder: **`find_package` (system pkg) →
pinned `FetchContent` (offline-graceful, like `cmake/FetchWannier2Sparse.cmake`)
→ vendored single header**. C++11 throughout.

| Need | Pick | Mode | CMake gate |
|---|---|---|---|
| Logging | **spdlog** (`libspdlog-dev`) | find → fetch → in-house fallback | `LSQUANT_WITH_SPDLOG` (ON) |
| Error vocab | **tl::expected** (C++11 backport of `std::expected`) | vendored header | — |
| JSON parse/emit | **nlohmann/json** (`nlohmann-json3-dev`) | find → fetch → vendored header | `LSQ_WITH_JSON` |
| Schema validation | **pboettch/json-schema-validator** | find → fetch → semantic-only fallback | `LSQUANT_WITH_JSON_SCHEMA` (ON) |
| Fingerprint hash | **picosha2** (SHA-256) | vendored header | — |
| Deep profiling (later) | **Tracy** | opt-in | `LSQUANT_WITH_TRACY` (OFF) |
| HW counters (later) | **PAPI** / `perf_event_open` | opt-in | `LSQUANT_WITH_PAPI` (OFF) |
| Plots/HTML (Python) | **Plotly** + numpy | `requirements.txt` | `LSQUANT_BUILD_UTILITIES` |

In-house (no good tiny cross-platform dep exists): the resource probe interface +
Linux backend, the scoped timers, and the domain-specific **physics merge**
(fingerprint + deep-merge of results) — all thin.

## 5. Unified phased roadmap

Each phase is independently mergeable, touches **new files only** until the
explicit migration steps, and — per *"Goldens first, always"* — **lands with test
coverage before it is wired into the build** (smoke test that goldens stay
byte-identical for R-phases; unit/golden fixtures for I-phases).

| Phase | Deliverable | From |
|---|---|---|
| **R0** | `util/` scaffolding, CMake options, `init_reporting()` no-op, linked into `kpm_lib` | reporting |
| **R1** | Logging: spdlog wiring + `LSQ_LOG_*` facade; `-v/--log-file` in `lsquant.cpp` | reporting |
| **R2** | Timer: `Timer`/`ScopedTimer`/`RegionRegistry`; `time_station*` forward to it; instrument `compute.cpp` phases | reporting |
| **R3** | Errors: `lsquant::Error` hierarchy + `tl::expected`; top-level `try/catch` in driver; convert existing bare throws/`cerr` | reporting |
| **R4** | Resources: `ResourceProbe` + Linux backend; end-of-run summary line | reporting |
| **I0** | Vendor nlohmann/json; swap `parse_flat_json` internals behind unchanged `read_run_config` API (no behavior change) | I/O |
| **I1** | JSON-Schema + physics-aware hints; `lsquant inspect` prints hints | I/O |
| **I2** | Schema v2 (`system`/`numerics`/`observables`); `compute` expands multi-observable configs into existing `compute_*` | I/O |
| **I3** | `io::fingerprint` (SHA-256); stamp into outputs; show in `inspect` | I/O |
| **I4** | `<system>.results.json` writer; populate timing/errors/peak-RSS from the reporting layer; thin-grid spectra | I/O (needs R2–R4) |
| **I5** | `lsquant merge` (fingerprint-checked); supersede python merge scripts after parity | I/O |
| **I6** | `utilities/python/lsquant_report.py` + self-contained Plotly HTML dashboard | I/O |
| **later** | online dashboard site; optional HDF5 moment backend; Tracy/PAPI backends | both |

Suggested order: **R0→R4 then I0→I6**, or run the R and I tracks in parallel and
gate **I4** on R2–R4 being in.

## 6. Consolidated decision points

From the two detail docs — recommendations given, owner confirmation wanted.

| # | Decision | Recommendation |
|---|---|---|
| D1 | `util` vs `utils` directory/namespace | **`util`** (matches `lsquant::util`) |
| D2 | Logging: spdlog vs in-house | **spdlog** (facade keeps escape hatch) |
| D3 | Error vocab: `tl::expected` vs in-house `Result` | **tl::expected** (vendored, C++11) |
| D4 | Resource probe: in-house vs third-party | **in-house** interface + Linux backend |
| D5 | Migration of existing `cout`/`time_station` | **opportunistic**, with back-compat shims |
| E1 | Config format: JSON only vs also TOML | **JSON canonical**, TOML optional later |
| E2 | Fingerprint: SHA-256 vs xxHash | **SHA-256** (identity); xxHash if too slow |
| E3 | Raw moments: keep `.chebmom` vs HDF5 | **keep `.chebmom`**, HDF5 later |
| E4 | HTML dashboard: Plotly vs Matplotlib PNG | **Plotly** self-contained |
| E5 | Merge tool: C++ subcommand vs Python | **C++ `lsquant merge`** + python wrapper |
| E6 | `utilities/python` cleanup (build/ cruft, py3.6) | separate housekeeping pass |

## 7. Out of scope (both plans)

GPU/CUDA accounting and outputs, MPI/distributed result aggregation, a hosted
multi-user dashboard with auth, and a job scheduler — GPU/parallel infra stays
out of scope, consistent with the current `lsquant`-only supported surface and
the `legacy/` quarantine (`legacy/README.md`). Interfaces are designed so these
extend later without reworking call sites: `ResourceProbe` admits a per-rank/
per-device backend, and `results.json` is versioned
(`"schema": "lsquant.results/1"`).
