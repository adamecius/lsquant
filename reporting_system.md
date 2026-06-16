# Reporting system — plan

**Status:** proposal / research (rev 2, reconciled with the sanitized codebase).
No code changed by this document. Everything below is **additive** (new files
under a new `util/` namespace, opt-in CMake options) so it can land alongside
other in-flight work without touching existing sources until each integration
step is explicitly taken.

## 0. Reconciliation with the sanitize effort (rev 2)

Since rev 1, the other agent landed Sanitize Phases A–D. What changed and how this
plan adapts:

- **Pruned planning docs.** `ROADMAP.md`, `FLAVOR_INVENTORY.md`, `HANDOFF_*`, and
  `test/TEST_PLAN.md` were removed (Phase A). Earlier rev-1 references to
  `ROADMAP.md` are now re-anchored to the live policy in
  [`legacy/README.md`](legacy/README.md) and `README.md`.
- **`legacy/` quarantine (Phase D).** Eighteen single-purpose / non-orthogonal
  drivers moved to `legacy/` (commented `# [LEGACY]` in `src/CMakeLists.txt`,
  never built). The **supported surface shrank** to `lsquant <command>`
  (`compute`, `reconstruct`, `inspect`) plus ~13 oracle drivers still in `src/`.
  This *simplifies* integration: the reporting layer now needs to reach **one
  library (`kpm_lib`) + one primary driver (`lsquant`) + a handful of oracle
  drivers**, not 40+ `main()`s. **`legacy/` is read-only — this plan never
  touches it.**
- **"Goldens first, always" (legacy/README.md).** Any capability wired into the
  build must arrive with golden/test coverage first. Adopted here as an explicit
  **gate on every phase** (see §6) — even logging/timing land with a smoke test
  proving they don't perturb golden numerics.
- **Non-orthogonal path is now an intentional throwing placeholder** in the
  library (working impl preserved at `5cc0924`). That throw is the **first
  natural client** of the `lsquant::Error` hierarchy in §4.3 (a `NotImplemented`/
  `NumericError` carrying the "restore from 5cc0924" hint).
- **Data-ban philosophy (Phase C).** Generated artifacts are a regenerable cache,
  never committed. The reporting layer's log files (`lsquant.log`) are such
  artifacts → they must be **gitignored**, never written into `test/golden/`.

## 1. Goal

Give `lsquant` a single, coherent **reporting layer** covering four concerns it
currently lacks:

1. **Logging** — levelled, timestamped, routable (console + file), thread-aware.
2. **Timing / profiling** — one canonical scoped timer + named-region
   accumulators, replacing the two ad-hoc `time_station` classes.
3. **Error handling & reporting** — a consistent vocabulary for failures
   (exception hierarchy at API boundaries + a lightweight `Result`/`Status` for
   hot paths) that funnels into the logger.
4. **Resource accounting** — peak/current **memory (RSS)** and **CPU time**
   (user/sys, wall), sampled cheaply, **Ubuntu/Linux first** but behind an
   interface so macOS/Windows/HPC backends can slot in later.

Design constraint from the owner: **build on existing, proven products** rather
than writing a framework from scratch; design for the **target architecture**
(the single `lsquant` library + `lsquant` driver, the strangler-fig endpoint),
not the legacy per-executable layout.

## 2. What exists today (audit)

| Concern | Current state | Verdict |
|---|---|---|
| Timing | Two header-only classes: `include/time_station.hpp` (`time_station`, `int` µs, prints to `std::cout`) and `src/time_station_2.hpp` (`time_station_2`, `long int` µs, `+=` accumulation, depends on `static_vars.hpp` for `r_type`). Overlapping, divergent, both hard-print to stdout. | **Consolidate** into one util; keep a thin back-compat shim. |
| Logging | None. ~40 files write directly to `std::cout`/`std::cerr`/`printf` (see `grep -c std::cout`). No levels, no timestamps, no file sink, no thread id. | **New.** |
| Errors | Mixed: CLI layer returns codes + `std::cerr` (`src/lsquant.cpp`, `compute.cpp`); library throws bare `std::runtime_error` in a few spots (`eigen_sparse_matrix.cpp:137`, `chebyshev_solver.cpp:137,286`); much of the numeric core has **no** error checking. | **Standardize**, additively. |
| Resources | None. No RSS/CPU/peak-memory reporting anywhere. | **New.** |

Relevant conventions to honour:

- C++ standard is **C++11** (`CMakeLists.txt:13`). This rules out `std::optional`
  (C++17), `std::expected` (C++23), `std::filesystem`, and `[[nodiscard]]`. Any
  vocabulary type must be C++11-clean.
- Dependencies are resolved with `find_package` / `pkg_config` as **imported
  targets**; no hand-passed toolchain (`CMakeLists.txt:30-60`).
- There is already a precedent for pulling an **external pinned product** with
  graceful offline degradation: `cmake/FetchWannier2Sparse.cmake` +
  `scripts/fetch_wannier2sparse.sh`. The reporting deps should follow the **same
  pattern** (system package first → pinned fetch → vendored header fallback).
- Namespaces in the target library are `lsquant::` (driver/config/compute) over
  `chebyshev::` / `qstates::` (numerics). New utilities should live under
  `lsquant::`.

## 3. Where it goes (directory & namespace)

`include/` is currently **flat**; `src/` is flat with a `Device/` subdir. Neither
`src/utils` nor `include/utils` exists yet. For the target architecture, introduce
a dedicated module — **`util/`** (singular, matching the `lsquant::util` namespace
we'll use; "utils" plural is fine too, pick one — this doc uses `util`):

```
include/util/            # public headers (consumed by library + drivers)
  log.hpp                #   lsquant::log  — logging facade
  timer.hpp              #   lsquant::util::Timer / ScopedTimer / RegionTimer
  error.hpp              #   lsquant::Error hierarchy + Result<T>/Status
  resource.hpp           #   lsquant::diag::ResourceProbe (interface) + snapshot
  report.hpp             #   umbrella include + one-line init (init_reporting())
src/util/
  log.cpp
  timer.cpp
  error.cpp
  resource.cpp
  resource_linux.cpp     #   Linux/proc + getrusage backend (default on Ubuntu)
  # resource_macos.cpp   #   (future) mach task_info
  # resource_win.cpp     #   (future) PSAPI GetProcessMemoryInfo
```

Rationale for `include/util/` (public) **+** `src/util/` (impl), not header-only:

- The logger and resource probe carry **state** (sinks, config, cached `/proc`
  fds) and a backend choice — a `.cpp` keeps that out of every translation unit
  and avoids ODR/duplicate-init hazards across the built `main()`s (the `lsquant`
  driver + the ~13 active oracle drivers).
- It mirrors how `kpm_lib` is already built (`src/CMakeLists.txt`): add the new
  `src/util/*.cpp` to `KPM_LIB_SOURCES` so every *built* driver linking `kpm_lib`
  gets it for free, no per-target wiring. Post-quarantine that is the `lsquant`
  driver + the ~13 active oracle drivers; the commented `# [LEGACY]` entries are
  untouched.
- Pure-template helpers (the scoped-timer RAII, `Result<T>`) stay header-only in
  `include/util/`; only the stateful parts get a `.cpp`.

> Decision point D1: `util` vs `utils` directory/namespace name. Recommendation:
> **`util`** (singular) to match `lsquant::util`. Cheap to change now, costly later.

## 4. Component plans

### 4.1 Logging — recommend **spdlog**

The de-facto C++ logging library: fast, MIT-licensed, header-only *or* compiled,
levels (`trace…critical`), multiple sinks (stdout-color, basic/rotating file),
`fmt`-style formatting, thread-safe, async option. Packaged on Ubuntu as
`libspdlog-dev`.

- **Facade, not raw spdlog everywhere.** Expose `lsquant::log` (thin wrapper in
  `include/util/log.hpp`) — `LSQ_LOG_INFO(...)`, `_WARN`, `_ERROR`, `_DEBUG`,
  `_TRACE`. Internals call spdlog. This means call sites never include
  `<spdlog/...>` directly, so spdlog can be swapped or compiled out behind a
  macro (`LINQT_WITH_SPDLOG`) without touching numerics.
- **Default sinks:** colored stderr (human) at `info`; optional rotating file
  `lsquant.log` when `--log-file` / `LSQUANT_LOG_FILE` is set. Level controlled by
  `LSQUANT_LOG_LEVEL` env + a `--verbose/-v` CLI flag wired in `lsquant.cpp`.
- **Format:** `[time] [level] [thread] message`. KPM loops are OpenMP-parallel,
  so thread id matters; keep per-iteration logging at `trace` to avoid hot-loop
  cost (level check is a cheap atomic load when disabled).
- **C++11:** spdlog supports C++11. Pin a C++11-compatible release in the fetch
  (e.g. the 1.x line) to stay within the project standard.

Fallback if a third-party dep is unwanted: a ~150-line in-house logger with the
same `LSQ_LOG_*` macro surface. The facade makes this a drop-in. (Recommend
spdlog; note this as **decision point D2**.)

Alternatives considered: **glog** (heavier, gflags coupling), **Boost.Log**
(compiled, heavy — and the project only uses Boost.Math header-only today),
**loguru**/**Quill** (fine, smaller ecosystem). spdlog wins on ubiquity + Ubuntu
packaging + the header-only escape hatch.

### 4.2 Timing — in-house, consolidate onto spdlog output

No external dep needed here; the existing `time_station` classes already do 90%.
Modernize and unify into `include/util/timer.hpp`:

- `lsquant::util::Timer` — start/stop/accumulate (`long long` µs via
  `std::chrono::steady_clock`), `elapsed_ms()/elapsed_us()`; supersedes both
  `time_station` and `time_station_2`.
- `lsquant::util::ScopedTimer` — RAII; logs `region X: 12.3 ms` at `debug` on
  scope exit. One-liner instrumentation: `LSQ_SCOPED("CorrelationExpansion");`.
- `lsquant::util::RegionRegistry` — named accumulators so repeated regions (per
  Chebyshev iteration) sum and print a summary table at shutdown (count, total,
  mean). This is the natural home for "where did the time go" in a KPM run.
- **Output through the logger**, not `std::cout` — fixes the current hard-coded
  printing and lets timing be silenced by log level.
- **Back-compat shim:** keep `time_station.hpp`/`time_station_2.hpp` as thin
  aliases that forward to `util::Timer` (deprecation comment), so the ~dozen
  current users keep compiling. Migrate call sites opportunistically, not in a
  big bang — keeps this change additive and auditable.

Optional deep-profiling backend (future, behind `LINQT_WITH_TRACY`): **Tracy**
frame profiler for flame-graph analysis of the hot loops. Out of scope for v1;
the `Timer` API is designed so a Tracy zone can be emitted alongside.

### 4.3 Error handling & reporting

Two layers, matching where failures actually occur:

1. **Exception hierarchy at API boundaries** (`include/util/error.hpp`):
   ```
   lsquant::Error : std::runtime_error      // base; carries context string
     ├─ ConfigError      (bad run.json / missing field)   // replaces ad-hoc cerr in lsquant.cpp
     ├─ IOError          (file open/read/write failures)  // replaces bare runtime_error in chebyshev_solver.cpp
     ├─ InputDataError   (non-Hermitian CSR, bad operator) // eigen_sparse_matrix.cpp:137
     ├─ NumericError     (NaN/Inf, non-convergence, OOB bounds)
     └─ ResourceError    (alloc failure, probe failure)
   ```
   Each carries a message + optional `where` (file:line via a `LSQ_THROW(Type,
   msg)` macro that captures `__FILE__/__LINE__/__func__`). The `lsquant` driver
   gets **one** top-level `try/catch(const lsquant::Error&)` in `main` that logs
   the error and returns a stable exit code — replacing scattered
   `cerr + return 1`.

2. **Lightweight `Result`/`Status` for hot paths** where throwing is undesirable
   (tight numeric loops, OpenMP regions — exceptions can't cross `#pragma omp`
   boundaries cleanly). Two options:
   - **`tl::expected`** — a single-header, **C++11-compatible** backport of
     `std::expected`. Vendor the one header under `include/util/third_party/`
     (zero build cost, no link). Gives `expected<T, Error>`.
   - A minimal in-house `Result<T>` / `Status` (≈80 lines) if vendoring is
     unwanted.
   Recommendation: **`tl::expected`** (decision point D3) — standard-shaped,
   migrates to `std::expected` for free when the project moves past C++11.

3. **Reporting:** errors don't just propagate — the throw-site macro and the
   top-level catch both emit through `LSQ_LOG_ERROR`, so every failure is
   timestamped and (optionally) file-logged. A `NumericError` for NaN/Inf in
   moments is the kind of thing that's currently silent and should be loud.

Policy (document in `error.hpp` header): **throw across module boundaries,
`Result` within a hot loop, assert only for invariants that indicate a bug.**

### 4.4 Resource accounting (memory / CPU) — Linux first, backend interface

Define the interface once; ship the Linux backend; leave stubs for others.

`include/util/resource.hpp`:
```cpp
namespace lsquant { namespace diag {
struct ResourceSnapshot {
  long long  rss_bytes;        // current resident set
  long long  peak_rss_bytes;   // high-water mark
  double     cpu_user_sec;     // user CPU
  double     cpu_sys_sec;      // system CPU
  double     wall_sec;         // since process/probe start
  long long  available_bytes;  // system free mem (optional, -1 if unknown)
};
struct ResourceProbe {                 // backend interface
  virtual ResourceSnapshot sample() = 0;
  virtual ~ResourceProbe() {}
};
ResourceProbe& probe();                // process-wide, backend chosen at build
ResourceSnapshot sample();             // convenience
}}
```

- **Linux/Ubuntu backend** (`src/util/resource_linux.cpp`, default):
  - Memory: `/proc/self/statm` or `/proc/self/status` (`VmRSS`, `VmHWM` for peak),
    *and* `getrusage(RUSAGE_SELF).ru_maxrss` (peak RSS; KiB on Linux). Cross-check
    cheap; reading `/proc` is a couple syscalls.
  - CPU: `getrusage` `ru_utime`/`ru_stime`, or `/proc/self/stat` utime/stime.
  - System free memory: `/proc/meminfo` `MemAvailable` (optional).
  - All POSIX/Linux; no extra dependency.
- **Backend selection** in CMake: `if(UNIX AND NOT APPLE)` → compile
  `resource_linux.cpp`; `APPLE`/`WIN32` compile a stub that returns `-1` fields
  and logs a one-time `warn` ("resource probe: backend not implemented on this
  platform"). The interface guarantees callers degrade gracefully, never break.
- **Future backends** (interface-ready, not built in v1): macOS `task_info`
  (`mach/task.h`), Windows `GetProcessMemoryInfo` (PSAPI). HPC extensions —
  **PAPI** or Linux `perf_event_open` for hardware counters, jemalloc/tcmalloc
  `mallctl` stats for allocator-level detail — slot in as additional `ResourceProbe`
  implementations behind `LINQT_WITH_PAPI` etc.
- **Usage pattern:** a `ScopedResourceReport` RAII (sibling of `ScopedTimer`)
  logs `ΔRSS` and CPU at the end of a major phase (moment computation,
  reconstruction). The `lsquant` driver prints a final one-line summary:
  `peak RSS 4.2 GiB · 312 s user · 18 s sys · 95 s wall`. Cheap, high-value for
  the large-system runs this library targets.

> No mature, tiny, truly-cross-platform RSS/CPU library exists worth a dependency
> here — every option is platform-specific glue under the hood. A ~150-line
> interface + Linux backend is the right call and keeps us dependency-light, which
> matches the project's posture. (decision point D4: confirm in-house vs pulling
> something like a process-stats helper — recommend **in-house**.)

## 5. Build integration (CMake)

Follow the existing imported-target + pinned-fetch + graceful-offline pattern.

- New options (default the heavy/optional bits **OFF** or auto):
  - `LINQT_WITH_SPDLOG` (default **ON**, auto-off if neither system pkg nor fetch
    available → falls back to in-house logger).
  - `LINQT_WITH_TRACY` (default **OFF**).
  - `LINQT_WITH_PAPI` (default **OFF**).
- **spdlog resolution order**, mirroring Eigen's fallback ladder
  (`CMakeLists.txt:32-46`):
  1. `find_package(spdlog QUIET)` → `spdlog::spdlog` (Ubuntu `libspdlog-dev`).
  2. else `cmake/FetchSpdlog.cmake` (pinned commit, `FetchContent`, offline-graceful
     like `FetchWannier2Sparse.cmake`).
  3. else header-only vendored / in-house logger; define `LSQ_NO_SPDLOG`.
- `tl::expected`: vendor the single header (no fetch, no link).
- Wire `src/util/*.cpp` into `KPM_LIB_SOURCES` (`src/CMakeLists.txt`) so all
  drivers inherit it; link `spdlog::spdlog` PUBLIC on `kpm_lib`.
- Apt one-liner addition to README: `libspdlog-dev` (optional).

## 6. Rollout (additive, auditable, non-blocking to other work)

Each phase is independently mergeable and touches **new files only** until the
explicit "migrate call sites" steps, which are mechanical and reviewable. Per
*"Goldens first, always"*, **every phase ships with a test** — at minimum a smoke
check that `ctest` golden output is byte-identical with the reporting layer
active (logging/timing must not perturb numerics; resource sampling must not
reorder OpenMP work).

- **Phase 0 — scaffolding (no behavior change):** create `include/util/` +
  `src/util/`, the headers/stubs, CMake options, and an `init_reporting()` no-op.
  Link into `kpm_lib`. Nothing calls it yet.
- **Phase 1 — logging:** land `log.hpp` + spdlog wiring + `init_reporting()` real
  body. Add `-v/--log-file` to `lsquant.cpp`. Do **not** mass-convert existing
  `std::cout` yet.
- **Phase 2 — timer:** land `timer.hpp/.cpp`; make `time_station*` forward to it.
  Instrument the 2–3 top-level phases in `compute.cpp` with `ScopedTimer`.
- **Phase 3 — errors:** land `error.hpp` + `tl::expected`; add the top-level
  `try/catch` in `lsquant.cpp`; convert the existing bare `runtime_error`/`cerr`
  sites (small, enumerated above) to the hierarchy.
- **Phase 4 — resources:** land `resource.hpp` + Linux backend; add the
  end-of-run summary line and `ScopedResourceReport` around moment computation.
- **Phase 5 — opportunistic migration:** convert remaining `std::cout`/`printf`
  call sites to `LSQ_LOG_*` as files are touched for other reasons. Never a
  big-bang rewrite.

## 7. Open decisions (need owner confirmation)

| # | Decision | Recommendation |
|---|---|---|
| D1 | Directory/namespace name `util` vs `utils` | `util` (matches `lsquant::util`) |
| D2 | Logging: spdlog dependency vs in-house | **spdlog** (facade keeps escape hatch) |
| D3 | Error vocabulary: `tl::expected` vs in-house `Result` | **tl::expected** (vendored header, C++11) |
| D4 | Resource probe: in-house interface vs third-party | **in-house** interface + Linux backend |
| D5 | How aggressively to migrate existing `cout`/`time_station` | **lazy/opportunistic**, with back-compat shims |

## 8. Out of scope (explicitly)

GPU/CUDA resource accounting, domain-decomposition/MPI-aware aggregation, and
distributed tracing — GPU/parallel infra remains out of scope for now (the
`legacy/` quarantine and the `lsquant`-only supported surface set the current
boundary). The `ResourceProbe` interface is designed so a future per-rank/
per-device probe can be added without reworking call sites.
