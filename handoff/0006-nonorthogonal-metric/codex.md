# Handoff 0006 (tester) — stage 6 gate: rotated models must reproduce G1, G2, G5

Read `AGENTS.md` sections 7, 8, 9, 12 and `brief.md` of this handoff. Everything on TCN,
`~/lsquant-unified`, branch `handoff/0006-codex` from `cpp-unified` after merging
`origin/handoff/0006-agi`. You edit only `scripts/check_goldens.sh`, `figures/`,
`goldens/stage0/PROVENANCE.md` and `tests.md`.

## 1. Fixture verification (independent of the developer's generator)

For each `goldens/stage0/nonorth/<model>/<variant>`: load `HAM`, `S`, `VX` (`VY`) and
the orthogonal originals; with dense numpy (`N = 2048` and `512`, fine) check:
`S` Hermitian and positive definite (`lambda_min > 0`), `cond(S)` equals the manifest
to 1e-6 and lies in `[5, 10]`; `eigvalsh(H', S')` equals `eigvalsh(H)` to 1e-10;
`H'`, `V'` Hermitian to 1e-14; sha256 of each file equals the manifest. Record the
numbers. If a check fails, the fixture is rejected and reported; nothing is regenerated
by you.

## 2. The gate (`scripts/check_goldens.sh 6`)

For `model in chain1d haldane`, `variant in real complex`, `route in lowdin generalised`,
`applier in exact ischo chebyshev`: run `stage6_nonorth` (window `[-2, 2]` for chain1d,
`[-4.5, 4.5]` for Haldane as G5, exact trace, `M = 64` chain1d / `128` Haldane, words
`dos` for chain1d, `kubo` for chain1d, `kuboxy` and `kubo` for Haldane) and compare
each produced legacy-layout file with the orthogonal golden (G1, G2, G5 `VX-VY` and
`VX-VX`) using `scripts/golden_diff.py`, printing the max relative difference and
PASS/FAIL against the brief's table (section 6). Also compare the two routes with each
other. Then the Hall plateau: reconstruct `sigma_xy` from each produced Haldane
`VX-VY` file (legacy `inline_kuboBastinFromChebmom <file> 10` unless stage 3 has
landed) and compute `C = (plateau/A) 2 pi` as `test/hall_haldane.sh` does; PASS if
`|C - 1| <= 0.05`. Stages 1 and 2 gates still run as part of this.

`DIFFERS` localisation (do not fix): the DOS wrong but the spectrum check of section 1
right points at the applier (residual line on stdout) or the rescale
(`S^-1 (aH' + bS')` vs `S^-1 (aH' + bI)`: a constant shift of the whole DOS); Kubo wrong
with DOS right points at the velocity (`V'` used without `S^-1` on the generalised
route, or without the two `S^-1/2` on Löwdin); only the `complex` variant wrong points
at a missing conjugation in `T^+`; both routes equal to each other but not to the golden
points at the fixture.

## 3. The figure (`figures/stage_6.py`)

Two rows. Top: chain1d DOS (reconstructed with stage-2 `lsquant.reconstruct.dos_curve`)
for orthogonal golden, Löwdin and generalised (exact applier), residual below. Bottom:
Haldane `sigma_xy(E)` for the three cases with the plateau and `C` in the caption, plus
the max relative moment differences per applier as a small table in the caption.
`FIG. 6.` Save under `figures/out/stage_6/`.

## 4. Report

`handoff/0006-nonorthogonal-metric/tests.md`: fixture checks with numbers; the full
gate table (model, variant, route, applier, max relative difference, PASS/FAIL,
applier residual and iterations); the `C` values; figure path and caption; smells (a
tolerance hidden in a `.cpp`, an applier printing, a change under `core/chebyshev`
other than the one allowed class, the engine touched). Commit on `handoff/0006-codex`,
push, report to Fable.
