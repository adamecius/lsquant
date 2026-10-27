# Handoff 0001 — review 1 (Fable to Agi)

Gate result on TCN (`~/lsquant-agi-build`, commit `a00d51c`): G1 and G2 **byte-identical**
at `OMP_NUM_THREADS = 1, 8, 64`; random-phase run twice identical; `grep` for
`basis|std::cout|getenv` under `lsquant/core` and `lsquant/io` empty; zero warnings.
The engine, the recursion, the generators and the writers are correct and the report
is complete. The following must change before the handoff closes. Numbers 1 to 3 are
defects; 4 to 6 are required cleanups. Nothing else is to be touched.

1. **`lsquant/CMakeLists.txt` moves the legacy `lsquant` executable to
   `build/legacy_bin`.** This breaks the legacy tree, which must stay buildable and
   runnable unchanged: `test/single_binary.sh`, `test/config_equivalence.sh`, the
   `example_*_smoke.sh` tests and every tutorial (`LSQUANT_BIN=build`,
   `_bin("lsquant")`) look for `build/lsquant`. Fix: in the root `CMakeLists.txt`
   use `add_subdirectory(lsquant lsquant_core)` so the new tree's binary directory is
   `build/lsquant_core/` and there is no name clash; delete the
   `if(TARGET lsquant) set_target_properties(...)` block and the
   `CMAKE_RUNTIME_OUTPUT_DIRECTORY` override from `lsquant/CMakeLists.txt`. The
   new executables then live in `build/lsquant_core/stage1_moments` and
   `build/lsquant_core/stage_1_check`; update `report.md` paths accordingly. Verify
   on TCN that the legacy `ctest --test-dir build` is fully green with
   `LSQUANT_BUILD_NEW_CORE=ON` and paste the summary in the report.

2. **Implicit Hermitian fold with a magic threshold in `SparseMatrix::FromCsr` and
   `FromTriplets`.** A generic matrix container must never alter the values it is
   given on a heuristic (`||M - M^+||/||M|| <= 1e-8`). The fold is a legacy
   convention of the CSR *file reader* (`ConvertFromCSR` folds `0.5 (M + M^+)`), so
   it belongs there, explicitly. Fix: remove the fold from both factories; add
   `SparseMatrix SparseMatrix::HermitianFold() const` (returns `0.5 (M + M^+)` with
   the same Eigen expression as the legacy, documented as the legacy convention
   that also normalises `-0` to `+0`); make `ReadCsrText` apply `HermitianFold()`
   unconditionally, as the legacy reader does, and say so in its comment and in the
   report. `FromTriplets` never folds. Re-run the gate: it must stay byte-identical.

3. **`MomentEngine::Run` sniffs the operator type with
   `dynamic_cast<const RescaledOperator*>` to fill `mu.window`.** The engine must
   know only `LinearOperator`; this is exactly the kind of coupling AGENTS.md
   section 3 forbids. Fix: remove the `dynamic_cast` and the include of
   `rescaled_operator.hpp` from `moment_engine.cpp`; the engine leaves `window`,
   `systemSize` (set it from `Dimension()`, that is fine), `numberOfStates`,
   `stateLabel` and `word` as it does now, and the caller sets `window`
   (`stage1_moments.cpp` already does).

4. **`ReadChebmom` legacy layouts leave the metadata incomplete.** The brief
   (section 7) requires `kpmWeighted = true` and `symmetrised = (rank == 2)` for
   `.chebmom1D` / `.chebmom2D`; also set `stateLabel` from the file name token after
   `_state` when the name has it, else leave it empty. Add this.

5. **`ChebyshevMoments` refuses rank > 2 in its constructor.** The container is
   rank-generic by design (`FlatIndex` already is); the stage-1 limit belongs to the
   engine only (it already throws on a second `Expand` per side). Remove the
   constructor check.

6. **`plan.hpp` re-exports the instruction types in a `lsquant::chebyshev`
   namespace.** Not in the brief; remove the aliases. One name per thing.

Notes, no action required: the `kMinParallelNnz = 30000` threshold in
`SparseMatrix::Multiply` is acceptable for now because determinism is per row; it
will become a constructor parameter when the Bloch operators arrive. The report's
"deviation 2" (const accessors) is accepted.

Procedure: fix 1 to 6 in that order, one commit each, on `handoff/0001-agi`; rsync
to `~/lsquant-agi-build` on TCN; rebuild; run `stage_1_check` at 1 and 8 threads;
run `stage1_moments` for `dos` and `kubo` on `test/golden/chain1d/operators` with
`--window -2 2 --trace exact --states 2048 --moments 64` into `/tmp/stage1_out`;
run the legacy `ctest --test-dir build --output-on-failure` and paste its summary;
update `report.md` (paths, the fold note, the ctest summary) and commit. Do not run
`cmp` against the goldens yourself; Fable does. Do not push.
