# TODO — leftovers before the PR

Tracks the open items from the reporting-layer + unified-I/O + PRL-figure work.
Everything below is **additive** and the build was green at the last C++ change;
these are finish-and-verify items, not known breakage.

## Threaded-SpMV effort (landed, branch `feat/threaded-spmv-phaseA`)

Stacked on base HEAD `86d6946`. Full suite **31/31 green at KPM_SEED=12345**, all 29
goldens **byte-identical** (the recursion uses only exact scalars a∈{1,2}, b∈{0,−1},
so the row-parallel SpMV is bit-for-bit the serial result).

- [x] **Phase A** (`79e1a7e`) — thread `SparseMatrixType::Multiply` (row-parallel fused
      OpenMP). New gate `spmv_threaded_bitexact`. Measured ~16–20× SpMV speedup/socket.
- [x] **Phase B** (`43055f5`) — implement `BatchMultiply` (block SpMM); new gate
      `spmm_batchmultiply_equiv` (R=1 bit-identical to `Multiply`).
- [x] **Phase C** (NUMA) — `schedule(static)` + autonuma converge to NUMA-local; numactl
      guidance in `docs/performance/README.md` (no first-touch refactor: Eigen/std::vector
      own storage; autonuma covers steady state).
- [x] **Phase F** (`b7f5945`) — graphene scaling benchmark (`bench/`), `docs/performance/README.md`,
      README front-page "Performance & scaling" section. SpMV reaches ~65% of the 1-socket
      STREAM read roofline; cost linear in N.
- [ ] **Optional W5** — threaded moment reductions (`vdot`) behind a default-OFF flag
      (`LSQ_THREADED_REDUCTION`); MUST stay off the byte-exact golden path.

## Must do before opening the PR

- [ ] **Run the full suite green.** `cmake --build build -j && ctest --test-dir build`.
      Targeted tests already pass (`reporting_smoke`, `unified_io_smoke`,
      `single_binary`, `config_equivalence`); the *full* 25+2-test run was not
      re-executed after the `io/` layer + `lsquant run|merge|fingerprint` +
      `compute_*` signature change. Confirm the example smoke tests still pass.
- [ ] **Regenerate tutorial 06 figures** (`fig_ey.png`, `fig_dp.png`). The script
      `examples/06_spin_relaxation_square/lsqrelax.py` is already restyled to
      `prl_style`, but the PNGs are still the old style — the run is slow because
      it uses the `exact` trace (loops over all basis states). Regenerate when
      convenient:
      `cd examples/06_spin_relaxation_square && LSQUANT_BIN=../../build python3 lsqrelax.py`
      Then confirm both PNGs are < 200 KB (the pre-commit hook limit).
- [ ] **Clean generated cruft from `examples/`** before staging: tutorial runs
      leave gitignored data (`operators/`, `*.chebmom*`, `BOUNDS`,
      `exact_chain1d_*`, `moments_*.chebmom1D`, `*.results.json`, `*.report.html`,
      staged `*.CSR`). They are gitignored, but tidy the tree.
- [ ] **Remove (or gitignore) `examples/_style_preview/`** — throwaway 600-dpi
      demo PNGs/PDFs from the style evaluation (some > 200 KB; regenerable via
      `python3 examples/style_demo.py`). Keep `examples/prl_style.py` and
      `examples/style_demo.py`.
- [ ] **Run `scripts/check_reporting.sh`** (per `CLAUDE.md`) and confirm the only
      advisory hits are the `lsquant.cpp` argv-usage `std::cerr` handlers (those
      are the program's argument-error channel, acceptable).
- [ ] **Verify all regenerated figures < 200 KB** (all current ones are; `fig_dos`
      is capped at 110 dpi for this reason).

## Follow-ups (can be separate PRs)

- [ ] **`lsquant run` spin_relaxation: time-domain reconstruction.** The
      `spin_relaxation` observable currently records the `.chebmomTD` moments file
      + parameters in `results.json` but does not yet reconstruct `S(t)` / the
      relaxation time `tau` into the doc (see the `note` field). Wire the
      `timeCorrelations` route so the thin `S(t)` curve and `tau` land in the
      unified output.
- [ ] **Unified-flow notes in the other tutorial READMEs.** Only
      `03_haldane_hall/README.md` shows the `lsquant run` + `results.json` +
      `lsquant_report.py` path. Add the equivalent short note to tutorials
      01/02/04/05/06 where it fits.
- [ ] **Optional dependency accelerators.** The JSON reader and logger are
      in-house (offline-safe, hook-compliant). If desired, add
      `find_package(nlohmann_json)` / `find_package(spdlog)` as *optional*
      accelerators behind the existing `LSQ_WITH_JSON` / `LSQUANT_WITH_SPDLOG`
      gates, keeping the in-house fallback (see `reporting_system.md` §5,
      `unified_io.md` §5).
- [ ] **Plan-doc decision points.** Confirm the open D1–D5 / E1–E6 choices in
      `infrastructure_plan.md` §6 now that the in-house route is implemented
      (e.g. E2 SHA-256 vs xxHash, E3 HDF5 for large moments).
- [ ] **macOS/Windows resource backends.** `include/util/resource.hpp` defines the
      `ResourceProbe` interface; only the Linux backend is implemented
      (`src/util/resource_linux.cpp`); others fall back to the stub. Add
      `task_info` (macOS) / PSAPI (Windows) backends if needed.

## Done this pass (for the PR description)

- Reporting layer (`include/util/`, `src/util/`): levelled logging, scoped
  timers, `lsquant::Error` hierarchy + `Result`, Linux resource probe; wired into
  `kpm_lib`; integrated into the `lsquant` driver (`-v`/`--log-file`, top-level
  catch, end-of-run summary) and `compute.cpp`. Test: `reporting_smoke`.
- Unified I/O (`include/io/`, `src/io/`): in-house JSON, SHA-256 system
  fingerprint, schema-v2 config + physics-aware validation hints, `results.json`
  container + fingerprint-checked merge; new subcommands `lsquant run|merge|
  fingerprint`; fingerprint shown in `inspect`. Test: `unified_io_smoke`.
- Presenter `utilities/python/lsquant_report.py` (summary + plots +
  self-contained HTML).
- PRL publication figure style `examples/prl_style.py` (+ `style_demo.py`); all
  six tutorial plot scripts restyled; figures regenerated for tutorials 01–05
  (06 pending, above).
- Agent directive: `CLAUDE.md` reporting requirement + `scripts/check_reporting.sh`.
- Plans reconciled: `infrastructure_plan.md`, `reporting_system.md`,
  `unified_io.md`.
