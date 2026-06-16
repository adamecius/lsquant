# CLAUDE.md — working agreement for agents on `lsquant`

LinQT/`lsquant` is a C++11 KPM quantum-transport library. Build:
`cmake -S . -B build && cmake --build build -j`. Test: `ctest --test-dir build`.

## Reporting requirement (MANDATORY before opening a PR)

Every **new implementation** — a new `compute_*` path, a new `lsquant` subcommand,
a new library routine that does I/O or non-trivial work — MUST be wired into the
reporting layer (`include/util/`, see [`reporting_system.md`](reporting_system.md))
**before** it goes into a PR. Concretely:

1. **Log through the facade, never raw streams.** Use `LSQ_LOG_INFO/_WARN/_ERROR/
   _DEBUG/_TRACE` (`util/log.hpp`). Do **not** add `std::cout`/`std::cerr`/
   `printf` to new code. User-facing progress is `info`; per-iteration detail is
   `trace`/`debug`.
2. **Time the significant phases.** Wrap heavy kernels with `LSQ_SCOPED("name")`
   (`util/timer.hpp`). Long-running work surfaces in the end-of-run summary.
3. **Fail through the error hierarchy.** Throw `lsquant::ConfigError` /
   `IOError` / `InputDataError` / `NumericError` / `ResourceError` /
   `NotImplemented` via `LSQ_THROW(Type, msg)` (`util/error.hpp`) — not bare
   `std::runtime_error`, `exit()`, or silent `return 1`. In OpenMP hot loops use
   `Result<T>`/`Status` instead of throwing across the parallel region.
4. **Surface resources for long runs.** A new long-running subcommand calls
   `lsquant::report_run_summary()` on success (peak RSS / CPU / wall + timing).
5. **Numerics must stay transparent.** Reporting must never change computed
   output. `ctest` (incl. `reporting_smoke` and `single_binary`) must be green;
   moment/spectrum files stay byte-identical with `-v`, `-q`, default.

Self-check before a PR: run **`scripts/check_reporting.sh`** (flags raw I/O and
bare throws in changed C++ outside `legacy/`) and ensure `ctest` passes. The
build is the authority; the script is an advisory pre-flight.

## Other standing rules (from the sanitize effort)

- **`legacy/` is read-only.** Never compile, modify, or "fix" it
  ([`legacy/README.md`](legacy/README.md)). The supported surface is
  `lsquant compute|reconstruct|inspect|run|merge|fingerprint` plus the active
  oracle drivers in `src/`.
- **Goldens first, always.** Any capability wired into the build arrives with
  golden/test coverage first.
- **Data ban.** Generated artifacts (`*.CSR`, `*.chebmom*`, log files,
  `*.results.json`) are a regenerable cache — never commit them. A pre-commit
  hook enforces a 200 KB size gate (`git config core.hooksPath .githooks` once per
  clone). Curated fixtures live only under `test/golden/` and `test/models/`.

## Figures

Tutorial figures follow the shared publication (APS/PRL) style in
[`examples/prl_style.py`](examples/prl_style.py): fixed column widths, math
labels, no titles, colour + linestyle + marker redundancy (greyscale- and
colour-blind-safe). Present `results.json` with
[`utilities/python/lsquant_report.py`](utilities/python/lsquant_report.py).

## Planning docs

[`infrastructure_plan.md`](infrastructure_plan.md) is the master plan; it links the
reporting and unified-I/O detail plans. [`TODO.md`](TODO.md) tracks open items.
