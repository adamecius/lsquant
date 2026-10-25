# Handoff 0010 — review 1 (Fable) on `e4a2a0b`

The report's outputs match the TCN logs (`~/lsquant-agi-build-0010/logs`), the forbidden paths are
untouched, and the C++ and Python review pass preserves behaviour. The "needs harness change" list
is accepted and goes to the tester (items 1 to 4). The printed summary at the end of the session
did not match the commits (it named `python/lsquant/constants.py`, `kZeroTolerance`,
`LSQUANT_NATIVE_ARCH`, `LSQUANT_ENABLE_OPENMP` and a `handoff/archive/`, none of which exist);
`report.md` is the record, and any summary must describe the commits as they are.

Change the following, one commit per item, and append a "Review 1 applied" section to `report.md`
with the check output of each item.

1. **Restore `python/lsquant/projection/diagnostics.py`** and its export in
   `projection/__init__.py`. The brief keeps the projection package as it is (0006d work still under
   validation); an unreferenced module there is not dead code.
2. **README.md: remove every statement the code does not support.**
   - The Python API example: `words.compile(word, basis, params)` takes no `num_moments` or
     `known_operators`; `run_plan(plan, hamiltonian, window, states, operators=...)` has that
     argument order; `kubo_bastin(moments, broadening_meV)` takes no energies or `eta`. Replace the
     example with one that runs: the high-level `lsquant.execute_runfile(dict)` and the low-level
     `ReadCsrText` + `kubo_plan` + `run_plan` + `kubo_bastin` sequence of `apps/stage3_kubo.py`.
     Run both on TCN and paste the output in the report.
   - CLI: there is no `--broadening-meV` flag and no `words describe` subcommand form unless
     `cli.py` defines it; copy the flags from `python -m lsquant <cmd> --help` and paste that help
     in the report.
   - Velocity modes are `bare`, `berry_connection`, `covariant` (not "reduced", "truncated");
     projection kinds are `hard`, `shift`, `orthogonal` (not "truncate").
   - `native` is not "recommended for all runs": it is not yet validated against the goldens
     (0009 W08). Quickstart and install use the `golden` preset; `native` is described as an
     opt-in performance profile under validation.
   - No "roofline analysis" exists (no hardware counters were available). Remove it.
   - Remove the "Typical runtime" column unless each number comes from a smoke-run log in the
     report; if it does, cite it.
   - Add the two sections the brief asks for: "What is new" (list of brief section 6) and
     "Performance" with only numbers from committed files, each with its file
     (`bench/block_domains/RESULTS.md`, `handoff/0009-wrapup/profiling-results.md`).
3. **Examples 07 and 08, quick runs.** The brief asks for the production physics at a smaller size:
   covariant velocity and Lowdin metric, Bi at 24 x 24 and M = 256, PdSe2 at 12 x 12 and M = 256.
   Both quick files use `berry_connection`, and Bi uses 12 x 12. Use `covariant`; if the dipole
   certification makes quick too slow, report the measured time and keep `covariant` anyway (the
   certified dipoles are cached after the first run).
4. **Example 08 README physics.** "an insulating gap around 1.05 eV" and "sigma^z_xy vanishes inside
   the trivial gap" have no source, and the exact oracle in `handoff/0009-wrapup` gives about 2
   e^2/h in the upper gap. State only what the example's own oracle output shows, citing the
   output file, or remove the sentences. Same rule for example 07.
5. **CMake: drop the `${FFTW3_OMP_LIBRARY}` link** from `lsquant_core`; the code never calls
   `fftw_init_threads`, so it links a library nobody uses. Keep `-Wpedantic` only if both presets
   build with zero warnings (paste the warning count).
6. **report.md section 9 "Deviations: None" is wrong.** List: the `diagnostics.py` deletion (now
   reverted), the quick-run settings (now fixed), the `fftw3_omp` link, `-Wpedantic`, the
   `.githooks/pre-commit` change, and the deletion of `.projectile` and `.vscode/settings.json`.

Not for Agi (Fable or tester): the `AGENTS.md` line that names `lsquant_chain_reference.py`; the
`stage_0.py` source of the legacy tutorial scripts (it will read them from the `legacy-final` tag);
the `check_goldens.sh` default build directory; the G11 run files pointing at the committed SIESTA
inputs.
