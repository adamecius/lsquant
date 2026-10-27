# Handoff 0004 — review 2 (Fable to Agi)

Tester result on `handoff/0004-codex` at `3dac161` (fresh build of `be54b79`): **gate PASS,
65 of 65 byte-identical** (stages 1 to 3: 41; stage 4: G6 moments, the four G6 curves from
golden and from new moments, G7 `SX`, `SY`, `SZ`, at 16 and 1 threads), the four greps
empty, the self-check all `array_equal`. Tutorial 04: 3 of 3 IDENTICAL. Review-1 items are
in (`6c457c9`).

Three items before the handoff closes; nothing else is to be touched.

1. **Tutorials 02 and 05 exit 1 before computing anything.** `tutorial02.py` calls
   `make_mod.build(stage_dir, n, W, seed)` and `tutorial05.py` calls
   `make_square.build(stage_dir, L, W, seed)`; the generators define
   `build_chain(N, W, seed, outdir=".")` and `build_square(L, W, seed, outdir=".")`. Call
   the real functions with the real argument order. Then run both apps on TCN (02 at 1
   thread, about 95 s; 05 at 16 threads, about 10 min; nohup and poll) and paste their
   stdout into the report: the tester showed with a scratch alias that the outputs are
   then 193 of 193 and 16 of 16 IDENTICAL, so nothing else should move. The report
   section 5 must show every app the brief names actually run, not one smoke run.

2. **A physics word in the engine.** `TimeEvolutionStep` (`plan.hpp`, used in
   `moment_engine.cpp`) is outside the allowed names (`TimeAxis`, `EvolveRegister`,
   `EvolveDerived`, `ApplyDerivedSeries`). Rename it `AxisStep`, and rewrite the `plan.hpp`
   comments that say "time-evolution" and "mean-square" in the engine's own terms (a
   trailing tensor axis, register updates between slices); the physics stays in the
   docstrings of `msd_plan`, `projected_evolution_plan` and `time_correlation_plan`.
   Update the binding and `_core.pyi` names accordingly.

3. **`tutorial04.py` reads the legacy `test/models/chain1d_mag/operators`.** As in stage 3
   (`tutorial03 --operators`): add `--operators` with default
   `goldens/stage0/chain1d_mag/operators` resolved from the repository root, and never read
   `test/`. Also say in the `stage4_time` help text that `--batch` does not apply to the
   time-axis path (the tester noted it is silently ignored).

Carried to closure by Fable: the MSD velocity rescaled like the Hamiltonian (legacy
convention, AGENTS.md section 10); the thread-count slowdown of the time-evolution path
(unconditional parallel SpMV on small systems) goes to the stage-8 performance item with
F9.

Procedure: one commit per item on `handoff/0004-agi`; rsync to `~/lsquant-agi-build` on TCN
(`ssh jh-x`), rebuild, rerun `stage_4_check.py` and the three tutorials, paste the outputs
into `report.md`, commit. Do not push.
