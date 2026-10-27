# Handoff 0004 (tester) — stage 4 gate: G6 and G7 byte for byte, tutorials 02, 04, 05

Read `AGENTS.md` sections 8, 9, 12 and `brief.md`. TCN, `~/lsquant-unified`, branch
`handoff/0004-codex` from `cpp-unified` after merging `origin/handoff/0004-agi`. Edit only
`scripts/check_goldens.sh`, `figures/`, `tests.md`.

## 1. Build and self-check

Stage-2 cmake line; fresh build; `tests/dev/stage_4_check.py` with
`OMP_NUM_THREADS=16`, output recorded (it runs one legacy binary from `build/`).

## 2. The gate (`scripts/check_goldens.sh 4`)

1. G6 moments: `stage4_time --word msd` on `goldens/stage0/chain1d_ballistic/operators`
   (window `-2 2`, exact, 256 states, `M = 256`, 64 steps, `tmax 20`); `cmp` the
   `.chebmomTD`. Then the four curves (`--fermi-energy 0 0.5 1.0 1.5`, 10 meV) `cmp`
   against `goldens/stage0/chain1d_ballistic/mean...EF*.dat`.
2. G7 moments: `--word projected` for `SX`, `SY`, `SZ` with `PX` on
   `goldens/stage0/chain1d_mag/operators` (`M = 32`, 96 steps, `tmax 95`, exact,
   window as the brief's section 6 resolves it); `cmp` the three `.chebmomTD`.
3. At `OMP_NUM_THREADS = 1` and `16`.
4. Stages 1 to 3 rerun; code-line grep on `core` and `io` extended with
   `time|velocity|spin|projector|bessel|gsl` (the identifiers `TimeAxis`,
   `EvolveRegister`, `EvolveDerived`, `ApplyDerivedSeries` are allowed: match whole
   words, lowercase, in code lines).

Localisation on `DIFFERS`: time slice 0 wrong only (seed initialisation), all slices
scaled (the `(2-delta)/N` weights or the global phase), drift growing with `n`
(series length `L`, the `|J_n| <= 1e-15` cut, or the commutator operator not rescaled
with `b`), `SY` wrong with `SX` right (a conjugation in `MirrorToBra` copying), curves
wrong with moments right (reconstruction axis `n*dt` or `std::to_string` of `E_F`).

## 3. The figure (`figures/stage_4.py`)

Three rows from the tutorial apps run into temporary directories and `cmp`'d against
`goldens/stage0/tutorials/02|04|05/captured/` (match by name; several captures of the
same name are compared each): (a) tutorial 04, `<S_x>`, `<S_y>`, `<S_z>` vs `cos`,
`-sin`, 0, residual; (b) tutorial 02, the localisation panel (`MSD` plateau) new vs
golden; (c) tutorial 05, `D(E,t)` for the disorder sweep new vs golden. Caption
`FIG. 4.` with every parameter and the `cmp` counts. `figures/out/stage_4/`.

## 4. Report

`handoff/0004-time-evolution/tests.md`: as before, plus the wall time of tutorial 05
and 02 (they are the slow ones) and smells (a physics word in the engine, GSL included
from `core`, a Clenshaw sum where forward accumulation was required, scipy Bessel used
on the gate path). Commit on `handoff/0004-codex`, push, report to Fable.
