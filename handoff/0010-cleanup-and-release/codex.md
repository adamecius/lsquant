# Handoff 0010 (tester) — final tests on `experimental` before the pull request to `main`

Read `AGENTS.md` sections 8, 9, 12, `brief.md`, `report.md` (section 8 is the harness list) and
`review-1.md`. Run everything on TCN in your own worktree and build directory (shared machine
rule: record the PIDs you start, stop only those). Branch `handoff/0010-codex` from
`origin/experimental` at the commit Fable names. Edit only `scripts/check_goldens.sh`, `figures/`,
`runs/goldens/`, the `fig_*` files under `examples/`, and `tests.md` in this directory.
`KPM_SEED=12345` everywhere.

## 1. Harness changes approved by Fable (the only ones)

1. `stage7_deleted`: add `AGENTS.md` and `CLAUDE.md` to the record exclusion of the `inline_*` text
   check, next to `handoff/`.
2. `figures/stage_0.py`: copy the legacy tutorial scripts from the `legacy-final` tag
   (`git worktree add <tmp> legacy-final`, or `git show legacy-final:examples/...`) instead of the
   current `examples/`; stage 0 is legacy provenance and must not depend on the current tree.
3. `scripts/check_goldens.sh`: when `$ROOT/build/lsquant_core` is absent, fall back to
   `$ROOT/build/golden/lsquant_core` (the `golden` preset).
4. `runs/goldens/g11_*.json` and the `6b`/`6c` sections of the harness: point the SIESTA model
   directories at `examples/07_spin_hall_bi_bilayer/siesta` and
   `examples/08_spin_hall_pdse2_monolayer/siesta`. First check that these files are byte-identical
   to `~/shc_bi/run` and `~/shc_pdse2/run` (`SHA256SUMS` in each directory) and record it.

No tolerance, seed or golden changes.

## 2. Build

Fresh `cmake --preset golden` and `cmake --preset native`, build both, record the warning count
(expected 0). `cmake --install` into a temporary prefix and `python -c "import lsquant,
lsquant._core"` from `/tmp`.

## 3. Gates — "the same as before"

The reference is the 2026-10-25 state: candidate v5 (`bbf81cf`) and the gate logs in
`~/lsquant-wrapup-20261024/logs` on TCN.

1. `scripts/check_goldens.sh 7` with the golden build, at 16 and at 1 thread: the same numeric files
   IDENTICAL as on 2026-10-25 (70), every text check PASS (the `inline_*` check passed nowhere
   before; it must pass now).
2. `scripts/check_goldens.sh 6only` with the golden build: same verdicts as on `bbf81cf`.
3. Developer checks `tests/dev/stage_6b_check.py`, `stage_6c_check.py`,
   `stage_6d_check.py --moments 8192 --backend folded`: same output as on `bbf81cf` (build `bbf81cf`
   once in a separate directory for the comparison).
4. The 24 G11 run files under `runs/goldens/g11_*` with the golden build: compare every output
   byte for byte with the same run files on `bbf81cf`.
5. The `native` build is not a gate: run gate 7 with it once and report identical/differing counts.

## 4. Examples

For each of `examples/01` to `08`, from a clean copy: run every run file as its README says, run
`plot.py`, commit the `fig_*` files. For 01 to 06 overlay the new figure data on the stage-0
captures and report the max abs difference per curve. For 07 and 08 run `quick` and `production`;
report `sigma^z_xy` at the README's reference energies from KPM and from the oracle, and the wall
time.

## 5. tests.md

Per gate: IDENTICAL or the first differing file, element and relative difference, and where in the
code the difference comes from. Build logs, install check, example numbers and timings, figure
paths with a one-line caption each, anything that smells. A verdict line at the top:
"same as before: YES/NO".
