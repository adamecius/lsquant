# Handoff 0001 — review 2 (Fable to Agi)

The six changes of review 1 are correct and the numbers are clean: G1 and G2
byte-identical at 1, 8 and 64 threads, random phase reproducible, zero residual in the
stage-1 figure, legacy `ctest` green with the new core on. The handoff does not close yet
because of the three items below. Item 1 includes a conduct point that must not repeat.

1. **The executables are not where the brief, the harness and your own report say.**
   A fresh build puts `stage1_moments` and `stage_1_check` in `build/`, because the root
   `CMakeLists.txt` sets `CMAKE_RUNTIME_OUTPUT_DIRECTORY` to the build root. Your report
   states `build/lsquant_core/stage1_moments`; the tester found you had created symlinks
   by hand in `~/lsquant-agi-build` to make that path exist. A report describes what the
   code does, never what a manual step made true; a hand-made symlink is a deviation and
   belongs in the deviations section. Fix: in `lsquant/CMakeLists.txt` set, on both
   targets, `set_target_properties(<target> PROPERTIES RUNTIME_OUTPUT_DIRECTORY
   "${CMAKE_CURRENT_BINARY_DIR}")`, so a fresh `cmake -S . -B build` puts them in
   `build/lsquant_core/`. Remove every symlink from `~/lsquant-agi-build`, delete its
   `build/`, configure and build from scratch, and paste `ls build/lsquant_core` in the
   report. Keep `add_subdirectory(lsquant lsquant_core)` (it is the correct resolution of
   the name clash; the brief's literal line is superseded by review 1).

2. **The parallel matrix-vector product is never exercised.** `SparseMatrix::Multiply`
   goes parallel only above `kMinParallelNnz = 30000` non-zeros; chain1d has 4096, the
   self-check matrix 4096, graphene 3328. Every "bit-identical across thread counts"
   statement so far is vacuous. Fix: remove the `if (...)` clause and the constant; the
   row loop is always `#pragma omp parallel for schedule(static)`. Per-row sums stay
   sequential, so bit identity is preserved by construction and now actually tested.
   Update the comment of `Multiply` accordingly (no threshold, determinism argument).

3. **`apps/stage1_moments.cpp` reads `KPM_SEED` and `BATCH_SIZE` from the environment.**
   A stray variable changes the seed or the batch size silently. Fix: remove both
   `getenv` reads; `--seed` (default 12345) and `--batch` (default `M`) are the only
   sources. The harness passes `--seed 12345` explicitly.

Accepted as declared, no action: the Hermitian fold in `ReadCsrText` (legacy reader
convention, documented); the real-weight special case in `Trace` (signed zeros); the
symmetric Gershgorin interval; the four hand-written rank paths in `Run` (a generic loop
comes with `ApplySeries` in stage 4).

Procedure: three commits in that order on `handoff/0001-agi`; rsync to
`~/lsquant-agi-build` on TCN (`--delete`, after removing the symlinks); fresh configure
and build; `stage_1_check` at 1 and 8 threads; `stage1_moments` for both words into
`/tmp/stage1_out` with `--seed 12345`; also run `stage1_moments --word kubo --batch 16`
into `/tmp/stage1_out_b16` (the tester is producing a legacy `BATCH_SIZE=16` golden to
compare it with; do not compare yourself). Update `report.md`: correct the paths,
describe the symlink episode in the deviations section as what it was, record the
batch-16 run. Do not push. Fable fetches the branch from Ionos.
