# Handoff 0007 (tester) — stage 7 gate: every golden from words, byte for byte; legacy gone

Read `AGENTS.md` sections 2, 3, 8, 9, 12 and `brief.md`. TCN, your own worktree of
`~/lsquant-unified`, branch `handoff/0007-codex` from `cpp-unified` after merging
`origin/handoff/0007-agi`. Edit only `scripts/check_goldens.sh`, `figures/`, `tests.md`.

## 1. Build and self-check
Stage-2 cmake line (plus FFTW); fresh build; `pip install --user mcp` if missing;
`tests/dev/stage_7_check.py`, output recorded.

## 2. The gate (`scripts/check_goldens.sh 7`)
1. For every run file under `runs/goldens/`, `python -m lsquant run <file>` into a
   temporary directory and `cmp` each produced file against its golden (G1, G2, G3, G4,
   G4c, G5, G6, G7, G9, G10 and the tutorial captures 01 to 05): one line per file,
   IDENTICAL or DIFFERS, totals per golden. At 16 and 1 threads.
2. The same through the MCP server: start it, call `run` for `g1.json` and `g5.json`,
   poll `status`, `cmp` the outputs (2 files). Record the JSON exchanged.
3. Refusals: the five malformed words of the self-check through the CLI, exit code
   non-zero, message on stderr; a run file with an unknown key refused before any load.
4. Stages 1 to 6 as before (now driven by the same run files; the old stage sections of
   the harness are replaced by section 1, keep the counts comparable: list the mapping).
5. After the developer's deletion commit: the tree builds without `src/`, `include/`,
   `test/`; `git ls-files src include test` is empty; the `legacy-final` tag exists and
   `git worktree add /tmp/legacy legacy-final && cmake` there still builds the legacy
   binaries (provenance stays reproducible); the gate passes on the deleted tree.
Localisation on `DIFFERS`: a name or layout difference points at `output` options; a
value difference points at the compiler's choice of operators or window (compare the
resolved run file with the stage app's arguments); a whole-golden failure points at a
missing `legacy_text_round_trip` or `bessel` option.

## 3. The figure (`figures/stage_7.py`)
The six tutorial figures regenerated from the run files, side by side with the stage-0
captures, residuals zero; caption `FIG. 7.` with the run file of each panel.
`figures/out/stage_7/`.

## 4. Report
`handoff/0007-words-cli-mcp/tests.md`: the per-golden table, the MCP transcript, the
refusals, the deleted-tree checks, smells (a word approximated instead of refused, a
basis word under `core`, an `eval`, a legacy file left behind, a golden re-blessed).
Commit on `handoff/0007-codex`, push, report to Fable.
