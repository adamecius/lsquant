# Handoff 0010 (tester) — release gate: goldens unchanged, examples regenerated

Read `AGENTS.md` sections 8, 9, 12 and `brief.md`. TCN, your own worktree and build directory
(shared machine rule: record your PIDs, stop only those). Branch `handoff/0010-codex` from
`origin/handoff/0010-agi` at the commit Fable names. Edit only `scripts/check_goldens.sh`,
`figures/`, the `fig_*` files under `examples/`, and `tests.md`.

## 1. Build
Fresh configure and build from each preset (`cmake --preset golden`, `cmake --preset native`);
record warnings (expected: none). `cmake --install` into a temporary prefix, `pip install ./python`
into a fresh venv, `python -c "import lsquant, lsquant._core"` from `/tmp`.

## 2. Harness changes (the only ones)
1. `stage7_deleted`: add `AGENTS.md` and `CLAUDE.md` to the record exclusion of the `inline_*`
   text check, next to `handoff/`.
2. Every path listed under "needs harness change" in `report.md` that Fable approved.
Nothing else in the harness moves; no tolerance, seed or golden changes.

## 3. Gate
`KPM_SEED=12345 scripts/check_goldens.sh 7` with the golden build, at 16 and at 1 thread. Expected:
the same 70 numeric files identical as on 2026-10-25, every text check PASS. Then the 24 G11 run
files under `runs/goldens/g11_*` with the golden build: compare with the outputs of the same run files on the
unchanged starting point `bbf81cf` (candidate v5), built the same way: byte for byte, report the
first difference if any. The developer checks `tests/dev/stage_6b_check.py`, `stage_6c_check.py`
and `stage_6d_check.py --moments 8192 --backend folded` must give the same output as on `bbf81cf`.
The native build is not a gate; run gate 7 with it once and report identical/differing counts.

## 4. Examples
For each of `examples/01` to `08`: follow its README from a clean copy (no files from earlier
runs), run every run file, run `plot.py`, and commit the `fig_*` files. For 01 to 06, overlay the
new figure data on the stage-0 tutorial captures and report the max abs difference per curve. For
07 and 08, run `quick` and `production`; report `sigma^z_xy(E_F = 0)` from KPM and from the
oracle for both, and the wall time.

## 5. tests.md
Per golden: IDENTICAL or the first differing file, element and relative difference. Build logs,
install check, example timings and numbers, figure paths with a one-line caption each, and
anything that smells (a number in a README with no source, a hard-coded path, a figure that
changed shape).
