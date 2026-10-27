# Handoff 0003 (tester) — stage 3 gate: G4 and G5 byte for byte from Python

Read `AGENTS.md` sections 8, 9, 12 and `brief.md` of this handoff. Everything on TCN,
`~/lsquant-unified`, branch `handoff/0003-codex` from `cpp-unified` after merging
`origin/handoff/0003-agi`. You edit only `scripts/check_goldens.sh`, `figures/`,
`tests.md`.

## 1. Build and self-check

Fresh build with the new core and the Python module (same cmake line as stage 2);
`PYTHONPATH=build/lsquant_core:python python3 tests/dev/stage_3_check.py`, output
recorded.

## 2. The gate (`scripts/check_goldens.sh 3`)

All with `PYTHONPATH=build/lsquant_core:python`, `KPM_SEED` irrelevant (the apps take
`--seed`), `OMP_NUM_THREADS=16` and once at `1`:

1. G5 moments: `stage3_kubo --operators goldens/stage0/haldane/operators --label haldane --left VX --right VY --window -4.5 4.5 --trace exact --states 512 --moments 128 --out $T` and the same with `--right VX`; `cmp` both `.chebmom2D` against `goldens/stage0/haldane/`.
2. G5 curves from the golden moment files: `--formula bastin|sea|surf` on the `VX-VY` file and `greenwood` on the `VX-VX` file, 10 meV; `cmp` the four `.dat`.
3. G4 moments: `--operators goldens/stage0/graphene/operators --label graphene_golden --left VX --right VX --window gershgorin --trace random --seed 12345 --states 1 --moments 64`; `cmp` against `goldens/stage0/graphene/NonEqOpVX-VX...chebmom2D`.
4. G4 curves: `bastin` at 64 meV (the regeneration script used `64`) against
   `KuboBastin_...JACKSON.dat`; `sea` and `surf` at 10 meV against `sea_surf/`;
   `greenwood_fft` and `bastin_fft` at 10 meV against `fft_grid/`.
5. Stages 1 and 2 gates rerun; grep on code lines of `lsquant/core` and `lsquant/io`
   (`lsquant/compat` is exempt by design, but must not be included from `core/`:
   `grep -rn "compat/" lsquant/core lsquant/io` empty).

Localisation on `DIFFERS` (do not fix): only the last significant digit on scattered
lines points at a table (`pow`, sympy instead of compat) or a summation order; a
constant factor points at the prefactor; a shift points at the axis; a different line
count points at `num_div`, `M_eff` or the dropped last point; the Bastin curve wrong
with Greenwood right points at the cumulative integral or `DgreenR`; `VX-VY` wrong with
`VX-VX` right points at the operator order of brief section 3.

## 3. The figure (`figures/stage_3.py`)

Tutorial 03 from Python: run `python -m lsquant.apps.tutorial03 --out <tmp>`, `cmp` its
four `.dat` with `goldens/stage0/tutorials/03/` (same names), and draw the tutorial's
two panels: (a) `sigma_xy` total, sea and surface in units of `e^2/h` (normalise as
`lsqhall.py` does), golden as markers, new as lines, residual below, `C` in the
legend; (b) `sigma_xx` (Greenwood). Caption `FIG. 3.` with `N = 512`, window
`[-4.5, 4.5]`, exact trace, `M = 128`, 10 meV, the `cmp` results and `C`. Save under
`figures/out/stage_3/`.

## 4. Report

`handoff/0003-rank2-reconstruction/tests.md`: build, self-check, gate output verbatim
(stages 1 to 3), tutorial-03 `cmp` list and `C`, figure path and caption, smells (a
formula or grid in `core/`, `compat` included from `core/`, numpy where `math` was
required, a table built with Python's complex power). Commit on `handoff/0003-codex`,
push, report to Fable.
