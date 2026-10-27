# Handoff 0004 — tester report (stage 4 gate: G6 and G7 from Python; tutorials 02, 04, 05; G4c)

Tester: Codex role (Opus subagent driven by Fable). Machine: TCN (`fatnode-tacn`).

Nothing under `lsquant/`, `python/lsquant/`, `tests/dev/` or the legacy tree was edited, and
no seed or tolerance was changed. Files edited:

- `scripts/check_goldens.sh`
- `figures/stage_4.py`
- `figures/out/stage_4/`
- this file
- for the closure run only, the new golden G4c (`goldens/stage0/graphene/offcentre*/`) and
  its `PROVENANCE.md` entry, as Fable asked. No existing golden file was touched.

The developer branch was tested twice:

| Run | `origin/handoff/0004-agi` | Checkout | Gate | Tutorials 04 / 02 / 05 | Logs in `figures/out/stage_4/` |
|---|---|---|---|---|---|
| 1 | `be54b79` | `~/lsquant-unified`, fast-forward | PASS, 65/65 | 3/3, 0/193, 0/16 (02 and 05 exit 1) | `*_run1_be54b79.*`, `fig_stage_4_diag_alias.*` |
| **2 (closure)** | **`aba179a`** (review 2: tutorials 02/05 call fix, `TimeEvolutionStep` renamed `AxisStep`, `tutorial04 --operators` default, `--batch` help; review 3: `reconstruct.py` axis `(x - shift)/scale`) | worktree `~/lsquant-unified-0004`, merge commit `83be4f7` (`merge_run2.log`) | **FAIL, 69 of 77 IDENTICAL**. The 8 failures are all G4c Bastin curves: energy column only, from a legacy defect (section 4) | **3/3, 193/193, 16/16** | `build.log`, `stage_4_check.log`, `check_goldens_stage4.log`, `figure_run.log`, `fig_stage_4.*`, `g4c_generation.log` |

| **2b (final)** | `aba179a`, same build; G4c re-blessed by Fable's decision (section 0) | worktree `~/lsquant-unified-0004` | **PASS, 71 of 71 IDENTICAL** | (unchanged from run 2) | `check_goldens_stage4.log`; the unblessed run-2 gate log is `check_goldens_stage4_run2_g4c_unblessed.log` |

## 0. Final run (2b): G4c re-blessed, gate PASS 71/71

Fable's decision on G4c: the legacy uniform-grid axis is the defect (D2), and the new code is
the reference for the energy column. It was applied as follows:

1. **The non-enclosing window is dropped.** `goldens/stage0/graphene/offcentre/` (`-3.5 2.5`)
   is removed, together with its gate lines. A window that excludes the spectrum is not a
   golden; S9 is accepted.
2. **`goldens/stage0/graphene/offcentre_enclosing/` (`BOUNDS` `-8.5 7.5`) is the G4c golden.**
   - Its `BOUNDS` and `.chebmom2D` stay the legacy files.
   - The legacy `KuboBastin_...JACKSON.dat` moved to
     `goldens/stage0/legacy_defects/offcentre_axis/`, together with a copy of its `BOUNDS` and
     a `README.md` stating the D2 mapping, the FFT-grid path being correct, the identical value
     column and the +0.4375 eV shift.
   - The curve golden is now the new-code output, blessed:
     `OMP_NUM_THREADS=16 python3 -m lsquant.apps.stage3_kubo --operators goldens/stage0/graphene/operators --label graphene_golden --left VX --right VY --window -8.5 7.5 --trace exact --states 256 --moments 128`
     (wall 32 s), then `--moments <that .chebmom2D> --formula bastin --broadening-meV 10`
     (wall 1 s), at `aba179a`. sha1 of the blessed curve: `318fad25b2279b3d518041074d65070e12332841`.
   - Checks before blessing:
     - the new moments are byte-identical to the legacy moments;
     - the curve from the legacy moments file is byte-identical to the curve from the new
       moments;
     - the value column is byte-identical to the legacy curve.
   - `PROVENANCE.md` G4c says exactly this, and its "Legacy defects" list gains the
     `offcentre_axis` entry.
3. **Gate.** The stage-3 section runs the enclosing window only, at 16 and 1 threads:
   moments are compared with the legacy moments, and the Bastin curve (from the golden moments
   and from the produced moments) with the blessed file.

Command: `OMP_NUM_THREADS=16 bash scripts/check_goldens.sh 4 build`, fresh log
`figures/out/stage_4/check_goldens_stage4.log`, wall 4207 s:

| Section | Files | IDENTICAL | DIFFERS |
|---|---|---|---|
| stage 1 | 11 | 11 | 0 |
| stage 2 | 6 | 6 | 0 |
| stage 3, G5 and G4 | 24 | 24 | 0 |
| stage 3, G4c (moments, Bastin from golden moments, Bastin from new moments; at 16 and 1) | 6 | 6 | 0 |
| stage 4 (G6 moments and 8 curves, G7 × 3; at 16 and 1) | 24 | 24 | 0 |
| **total** | **71** | **71** | **0** |

`== 71 files compared; failures: 0`, `GATE PASS (stage 4)`, `gate exit: 0`. All four greps
are empty.

Stage-4 walls in this run:

| Run | 16 threads | 1 thread |
|---|---|---|
| G6 moments | 744 s | 23 s |
| G7 SX | 388 s | 4 s |
| G7 SY | 202 s | 4 s |
| G7 SZ | 83 s | 4 s |

The load average was 60 to 95.

Sections 3 and 4 below describe run 2 before the re-blessing: G4c with both windows and the
legacy curves, 69 of 77. They are kept as the record of the localisation.

About the checkout for run 2: the shared checkout `~/lsquant-unified` had been switched by
another session to `handoff/0006-codex`, with uncommitted stage-6 work in it. I did not touch
it. `git worktree add ~/lsquant-unified-0004 handoff/0004-codex` gave this handoff its own
tree and its own fresh `build/`. My first command in run 2 ran before I noticed the branch
switch. Its `cd` into the missing `figures/out/stage_4` failed, so its `git mv`, `commit`,
`merge` and build steps did nothing. I checked afterwards that the stage-6 checkout was
unchanged: HEAD `a5108e3`, same `git status`.

**Result (run 2, `aba179a`):**

- **Stages 1, 2 and 4 are all IDENTICAL**, at the thread counts of each stage; all greps
  are empty.
- **Stage 3 has every pre-existing golden IDENTICAL**, including the centred G4 and G5
  curves after the `reconstruct.py` axis change. The new G4c moments are IDENTICAL.
- **The new G4c Bastin curves DIFFER, in the energy column only.** The conductivity column is
  byte-identical as text, at 16 and 1 threads, whether reconstructed from the golden moments
  or from the new ones.
- **The defect is in the legacy, not in the new code.** The legacy rank-2 reconstruction
  writes `E = x HalfWidth + BandCenter/HalfWidth`. The new code writes the correct
  `x HalfWidth + BandCenter`.
- **The tutorials from Python are all IDENTICAL:** 04 is 3 of 3, 02 is 193 of 193, and 05 is
  16 of 16.

## 1. Build (run 2, fresh)

```
$ rm -rf build
$ cmake -S . -B build -DLSQUANT_BUILD_NEW_CORE=ON -DLSQUANT_BUILD_PYTHON=ON -Dpybind11_DIR=...
cmake exit: 0
$ cmake --build build -j 32
build exit: 0
```

`figures/out/stage_4/build.log` shows 0 warnings under `lsquant/`. As in every run, all 101
warnings are in the legacy tree.

For G4c a second, legacy-only build was made exactly as `PROVENANCE.md` records the stage-0
build: `cmake -S . -B build_legacy -DLSQUANT_BUILD_NEW_CORE=OFF -DLSQUANT_BUILD_INPUT_TOOLS=ON`,
then `cmake --build build_legacy -j 32`. It is not committed. The legacy tree is unchanged
since `8ecbf83` except for CMake options (`git diff --stat 8ecbf83 HEAD -- src include` is
empty).

## 2. Developer self-check (run 2)

```
$ PYTHONPATH=build/lsquant_core:python OMP_NUM_THREADS=1 python3 -u tests/dev/stage_4_check.py
(a) G7 chain1d_mag Gershgorin window = [-2.31, 2.31], BandWidth = 4.62
(a) G7 TimeProjSX-PX array_equal vs golden: True
(a) G7 TimeProjSY-PX array_equal vs golden: True
(a) G7 TimeProjSZ-PX array_equal vs golden: True
(b) G6 CorrelationVXchain1d array_equal vs golden: True
(b) G6 EF=+0.0 curve differing lines: 0 (total lines: 64)
(b) G6 EF=+0.5 curve differing lines: 0 (total lines: 64)
(b) G6 EF=+1.0 curve differing lines: 0 (total lines: 64)
(b) G6 EF=+1.5 curve differing lines: 0 (total lines: 64)
(c) compat GSL Bessel vs SciPy jv max_abs_diff = 5.551115e-16
(d) time_correlation_plan chain1d VX-VX (M=64, steps=16, tmax=5, exact) array_equal vs inline_compute-kpm-CorrOp: True
wall 34.50 s
exit: 0
```

I first started it at OMP_NUM_THREADS = 16, as in run 1. Load average reached 140, and after
49 minutes it was still in part (b), so I stopped it. Its partial output, (a) all `True`, is
kept in `stage_4_check_T16_aborted.log`. The full check was then run at 1 thread, the
developer's own setting; it took 34.5 s. The gate covers the 16-thread case.

## 3. Gate (run 2): `OMP_NUM_THREADS=16 bash scripts/check_goldens.sh 4 build`

Log: `figures/out/stage_4/check_goldens_stage4.log`. Wall time: 6810 s. The harness now
includes G4c in the stage-3 section (section 4).

| Section | Files | IDENTICAL | DIFFERS |
|---|---|---|---|
| stage 1 (G1, G2 at 1, 8, 64 threads; batched G2; random-phase rerun) | 11 | 11 | 0 |
| stage 2 (G1, G3 from Python at 1 and 16 threads) | 6 | 6 | 0 |
| stage 3, old goldens (G5 moments and 4 curves, G4 moments and 5 curves, at 16 and 1) | 24 | 24 | 0 |
| stage 3, **G4c** (2 windows × (moments + Bastin from golden moments + Bastin from new moments), at 16 and 1) | 12 | 4 (moments) | **8 (Bastin curves)** |
| stage 4 (G6 moments, 8 G6 curves, 3 G7 moments, at 16 and 1) | 24 | 24 | 0 |
| **total** | **77** | **69** | **8** |

Greps, all empty on code lines:
- `basis|std::cout|getenv` in core;
- `delta|DOS|Kubo|basis|std::cout|getenv` in core and io;
- `compat/` in core and io;
- `time|velocity|spin|projector|bessel|gsl` as whole lowercase words in core and io.

The `plan.hpp` comment lines that said "time-evolution" and "mean-square" are gone from
the excluded-comment list.

Stage-4 walls at 16 threads: G6 moments 3258 s at 16 threads vs 23 s at 1; G7 53, 52, 177 s at 16 vs 4, 3, 4 s at 1.

The verdict line is `GATE FAIL (stage 4)`. Only the 8 G4c Bastin lines are DIFFERS; each is
followed by the harness's numeric localisation and an informational column note:

```
DIFFERS   [T=16] graphene/offcentre/KuboBastin_NonEqOpVX-VYgraphene_goldenKPM_M128x128_stateexactJACKSON.dat (G4c curve, 10 meV, from golden moments)
          first-diff line 1 col 1, new=-3.35 golden=-3.01667; ...
          note: value column identical as text; energy column golden - new: min 0.333328 max 0.33334 mean 0.333333 over 3839 lines
DIFFERS   [T=16] graphene/offcentre_enclosing/KuboBastin_...JACKSON.dat (G4c curve, 10 meV, from golden moments)
          first-diff line 1 col 1, new=-8.1 golden=-7.6625; ...
          note: value column identical as text; energy column golden - new: min 0.437495 max 0.437505 mean 0.4375 over 3839 lines
```

The same lines appear for the curves built from the new moments, and again at T=1.

## 4. G4c: the off-centre golden

### 4.1 What was generated

`goldens/stage0/graphene/offcentre/` was asked for: `BOUNDS` `-3.5 2.5`. I also added
`offcentre_enclosing/` (`BOUNDS` `-8.5 7.5`, the same centre of `-0.5`).

Each directory holds `BOUNDS`, `NonEqOpVX-VYgraphene_goldenKPM_M128x128_stateexact.chebmom2D`
and `KuboBastin_NonEqOpVX-VYgraphene_goldenKPM_M128x128_stateexactJACKSON.dat`. The
`PROVENANCE.md` entry is "G4c". The generation log is `g4c_generation.log`.

In a temporary directory I staged `operators/graphene_golden.{HAM,VX,VY}.CSR`, the `BOUNDS`
file and an exact local-state file over the 256 sites, then ran, with `KPM_SEED=12345
OMP_NUM_THREADS=16`:

```
inline_compute-kpm-nonEqOp graphene_golden VX VY 128 exact
inline_kuboBastinFromChebmom NonEqOpVX-VYgraphene_goldenKPM_M128x128_stateexact.chebmom2D 10
```

The legacy honoured `BOUNDS`: the headers are `256 6  -0.5` and `256 16  -0.5`.

Three decisions, recorded here:

1. **Legacy build: not `~/lsquant/build`.** That directory, suggested in the instructions,
   holds an older and unrelated legacy build: checkout `c1d9318`, binaries dated June 2025.
   It reconstructs on the grid `(-0.99, 0.99)` and **does not reproduce G5**: reconstructing
   the committed G5 VX-VY moments gives `-4.455...`, where the golden first line is
   `-4.275`. The repository's own legacy build does reproduce G5 Bastin byte for byte, so
   G4c comes from that build, made as in `PROVENANCE.md` (`build_legacy`, section 1). The
   moments from both legacy builds are identical (same sha1).
2. **The requested window does not enclose the spectrum.** By dense diagonalisation the
   graphene spectrum is `[-6.834, 6.834]` eV, and `[-3.5, 2.5]` covers less than half of it.
   The rescaled eigenvalues reach -2.11 and +2.44, so the Chebyshev moments grow
   geometrically: the curve values are up to 1e139. This golden still tests the arithmetic
   and the energy axis byte for byte, but it is not a physical result. I therefore added
   `offcentre_enclosing/` (`-8.5 7.5`), which encloses the spectrum. Its graphene
   `sigma_xy` is about 1e-16, i.e. 0, as expected.
3. **Gate coverage.** The gate compares the moments, and the Bastin curve both from the
   golden moments and from the new moments, at 16 and 1 threads.

### 4.2 Result and localisation

| Item | offcentre (`-3.5 2.5`) | offcentre_enclosing (`-8.5 7.5`) |
|---|---|---|
| moments (`.chebmom2D`) | IDENTICAL | IDENTICAL |
| Bastin, value column (text) | identical, 3839 lines | identical, 3839 lines |
| Bastin, energy column | DIFFERS: golden - new = +0.333333 on every line | DIFFERS: golden - new = +0.4375 on every line |
| new first / last E | -3.35 / 2.34852 (= -0.95·3 - 0.5 …) | -8.1 / 7.09604 |
| legacy first / last E | -3.01667 / 2.68185 | -7.6625 / 7.53354 |

The localisation is in the energy axis, not in the moments or the values. The offset is
exactly `BandCenter/HalfWidth - BandCenter`: -0.5/3 + 0.5 = 0.333333, and -0.5/8 + 0.5 =
0.4375.

- **New code (review 3):** `E = (x - shift)/scale = x·HalfWidth + BandCenter`. This is
  correct, with a max deviation of 4e-16.
- **Legacy, this repository's tree:** `src/reconstruction.cpp` writes the uniform-grid Kubo
  curves with `write_pairs(values, energies, [sf, shf](double e){ return e / sf - shf; })`,
  where `sf = ScaleFactor() = 1/HalfWidth` and `shf = ShiftFactor() = -BandCenter/HalfWidth`
  (`include/chebyshev_moments.hpp:85,88`). That gives `E = x·HalfWidth + BandCenter/HalfWidth`,
  verified with a max deviation of 7e-6, the 6-digit text precision. **This is a legacy
  defect.** Review 3's statement that "the legacy maps `x = (E - BandCenter)/HalfWidth`" holds
  for the FFT-grid path (`e * hw + bc`, `reconstruction.cpp:150`) and for the old
  `~/lsquant` build, but not for this tree's uniform-grid Bastin, Sea, Surf and Greenwood
  (`reconstruction.cpp:147-148`, `199`, `223`).
- The same pattern appears in the unified `lsquant` binary's DOS output
  (`src/lsquant.cpp:108,240`: `x / ScaleFactor + ShiftFactor`, which is the sign the new code
  had before review 3). I did not test it.
- In the centred goldens `BandCenter = 0`, so both forms agree; that is why nothing else
  moved.

Consequence: G4c's curves faithfully encode a legacy defect. Two ways to close this, for
Fable to decide:

- (a) accept the energy-column difference with a written reason and re-bless G4c's
  `KuboBastin` files with the corrected axis; or
- (b) keep G4c's moments as the golden and compare only the value column of the curves.

No golden was edited and no code was changed. The gate therefore reports FAIL on these 8
files.

## 5. Figure (run 2): `MPLBACKEND=Agg OPENBLAS_NUM_THREADS=1 python3 figures/stage_4.py --build build --threads 1`

Thread count: run 1 used `--threads 16`. In run 2 the figure ran at 1 thread, because 16
threads on this load is about 20 times slower (S1): the run-1 attempt of tutorial 02 at 16
threads had managed only 28 of 192 realisations in 10 minutes. Thread invariance of the
engine is what the gate checks.

```
tutorial04: wall 14.0 s, 3 lines printed, 3 of 3 IDENTICAL (3 new captures)
tutorial02: wall 132.6 s, 193 lines printed, 193 of 193 IDENTICAL (192 new captures)
tutorial05: wall 223.1 s, 8 lines printed, 16 of 16 IDENTICAL (8 new captures)
max |new - golden|: (a) 0.0e+00  (b) 0.0e+00  (c) 0.0e+00; (a) max |golden - analytic| 3.997e-15
exit: 0
```

- `tutorial04` now reads `goldens/stage0/chain1d_mag/operators` by default.
- Wall times of the slow tutorials at 1 thread, with load average 90 to 140: tutorial 02
  took 132.6 s and tutorial 05 took 223.1 s. For comparison, tutorial 05 took 624 s at 16
  threads in run 1's diagnostic run.
- The figure is `figures/out/stage_4/fig_stage_4.png`. Its caption is in
  `fig_stage_4.caption.txt` and the cmp list of all 212 files in `fig_stage_4.cmp.txt`.

Caption (verbatim):

> FIG. 4. Stage 4 (handoff 0004): tutorials 04, 02 and 05 reproduced from Python by python3 -m lsquant.apps.tutorial04|tutorial02|tutorial05 (lsquant package over the pybind11 module lsquant._core, OMP_NUM_THREADS = 1), each run in a temporary directory; time evolution U(dt) = exp(-i H dt/hbar) by the Bessel-Chebyshev series (GSL J_n via lsquant._core.compat, |J_n| <= 1e-15 cut, hbar = 0.6582119624 eV fs), dt = tmax/(steps - 1). (a) Tutorial 04, Larmor precession on chain1d_mag (N = 128, t = -1, J_ex = 0.1 eV): m = 0 row of TimeProj{SX,SY,SZ}-PX (projector PX = (1 + S_x)/2, exact trace over the 128 site states, M = 32, 96 time steps, tmax = 95 fs, Gershgorin window [-2.31, 2.31] eV), normalised by <S_x(0)> as lsqspin.py does; new as filled markers, golden as open grey markers (every 4th step), cos(wt), -sin(wt), 0 with w = 2 J_ex/hbar as thin lines; residual new - analytic below (golden - analytic max 4.00e-15, round-off level); max |new - golden| = 0.0e+00. (b) Tutorial 02, localisation in the 1D Anderson chain (box disorder W = 1 eV, disorder seed and KPM_SEED = r for r = 1..16, one random-phase vector, window from BOUNDS, M = 256, 128 steps, tmax = 600 fs, MSD(E,t) at 10 meV): xi(E) = cal sqrt(<MSD>) over the last 30% of the time window, averaged over the 16 realisations, N = 256, 1024, 4096, E = 0, +-1, 1.9 eV, cal fixed once against the Born length 24(4 - E^2)/W^2 at N = 4096, E = 1 eV (lsqloc.py); new as lines, golden as open markers, Born dashed; residual new - golden below (max 0.0e+00 sites). (c) Tutorial 05, 2D square lattice L = 200 (N = 40000), disorder W = 0.5, 1, 2, 4 eV (seed 1, KPM_SEED = 1, one random-phase vector): D(E,t) = DeltaX^2/(2t) at E = 2 eV from the MSD (M = 256, 160 steps, tmax = 150 fs, 80 meV; the app also writes the DOS, M = 1024, 80 meV); new as lines, golden as open markers (every 8th step); residual new - golden below (max 0.0e+00). Byte comparison with goldens/stage0/tutorials/0N (captured/NNNN_<name> against the capture with the same index and name, top-level .dat by name; cmp criterion): tutorial 04 3 of 3 IDENTICAL; tutorial 02 193 of 193 IDENTICAL; tutorial 05 16 of 16 IDENTICAL. Wall times (TCN, load 60 to 100 on 64 cores): tutorial 04 14.0 s, tutorial 02 132.6 s, tutorial 05 223.1 s. All files IDENTICAL.

## 6. Defects

- **D1 (run 1): fixed.** `tutorial02` and `tutorial05` now call `build_chain` and
  `build_square` with the generators' argument order. They produce 193 of 193 and 16 of 16
  IDENTICAL.
- **D2 (new, legacy, not the new code; resolved for the gate in run 2b by blessing the
  new-code curve, with the legacy curve kept in `legacy_defects/offcentre_axis/`):** the energy axis of this tree's legacy uniform-grid
  rank-2 reconstruction is `x·HalfWidth + BandCenter/HalfWidth` (section 4.2). It needs to be
  recorded in AGENTS.md section 10, and G4c's curve comparison needs a decision.

## 7. Smells (run 2)

- **S1 (open, carried by Fable to stage 8): thread count and speed.** The time-evolution path
  is 10 to 20 times slower at 16 threads than at 1 on loaded TCN. The cause is the
  unconditional `#pragma omp parallel for` per SpMV in `sparse_matrix.cpp:126`. Run 2
  figures:
  - self-check: 34.5 s at 1 thread; still unfinished after 49 min at 16 threads;
  - tutorial04: 14 s at 1 thread;
  - stage-4 gate at 16 threads: G6 moments 3258 s at 16 threads vs 23 s at 1; G7 53, 52, 177 s at 16 vs 4, 3, 4 s at 1.

  Results are unaffected.
- **S2 (gone).** `TimeEvolutionStep` is renamed `AxisStep` in `plan.hpp`, `moment_engine.cpp`,
  `bindings.cpp` and `_core.pyi`. The `plan.hpp` comments no longer say "time-evolution" or
  "mean-square". The stage-4 grep is empty, and no lowercase physics word remains in the
  excluded comment lines of `plan.hpp`.
- **S3 (unchanged, fine):** GSL appears only in `compat/`.
- **S4 (unchanged, minor):** forward accumulation, not Clenshaw. One unused recursion step
  remains after the last coefficient of each series.
- **S5 (gone; it was a non-finding):** scipy's `jv` is not on any gate path; the default
  coefficients come from the GSL compat routine.
- **S6 (gone).** `tutorial04 --operators` defaults to `goldens/stage0/chain1d_mag/operators`,
  and no app reads `test/`.
- **S7 (carried by Fable to AGENTS.md section 10):** the MSD velocity is rescaled like H.
  This is harmless because every stage-4 window is centred.
- **S8 (resolved as documentation):** the `stage4_time --batch` help now says it does not
  apply to the time-axis path. The extra `.chebmom` output remains and is ignored by the gate.
- **S9 (accepted by Fable; resolved in run 2b by dropping `offcentre/`): the legacy G4c window.** As asked for, `[-3.5, 2.5]` excludes the graphene
  spectrum, so the moments diverge to about 1e138 (section 4.1). As a golden it is
  arithmetic-only, and it amplifies any future reordering by a huge factor: a single ulp
  change becomes a large absolute difference. `offcentre_enclosing/` is the physically
  meaningful off-centre case.
- **S10 (new): two legacy builds on TCN.** `~/lsquant/build` (`c1d9318`, June 2025) is not
  the build behind the stage-0 goldens: it does not reproduce G5. Any instruction that points
  at it for golden generation should point at a build of this repository's legacy tree
  instead.
- **No existing golden, seed or tolerance was touched.** The goldens added are G4c (two
  directories) and a `PROVENANCE.md` section.

## 8. Run 1 (`be54b79`) record

The run-1 sections follow unchanged, apart from heading levels. Their logs and figure were
renamed with the suffix `_run1_be54b79`: `build_run1_be54b79.log`, `stage_4_check_run1_be54b79.log`,
`check_goldens_stage4_run1_be54b79.log`, `figure_run_run1_be54b79.log`,
`fig_stage_4_run1_be54b79.{png,caption.txt,cmp.txt}`. The diagnostic files are
`fig_stage_4_diag_alias_run1_be54b79.*`. Run 1's smells S1 to S8 are restated, with their
current status, in section 7.

### A.1 Build (fresh)

```
$ rm -rf build
$ cmake -S . -B build -DLSQUANT_BUILD_NEW_CORE=ON -DLSQUANT_BUILD_PYTHON=ON -Dpybind11_DIR="$(python3 -c 'import pybind11; print(pybind11.get_cmake_dir())')"
cmake exit: 0
$ cmake --build build -j 32
[100%] Built target _core
build exit: 0
```

Full log: `figures/out/stage_4/build.log`. Warnings: **0 under `lsquant/`**. All 101
warnings are in the legacy tree, as in stages 1 to 3: 62 in `include/chebyshev_moments.hpp`,
32 in `include/units.hpp` and 7 in `src/chebyshev_momentsLocal.cpp`. `_core` links
`GSL::gsl` (`lsquant/CMakeLists.txt:51`). `lsquant_core` links only Eigen and OpenMP.

### A.2 Developer self-check (`tests/dev/stage_4_check.py`, OMP_NUM_THREADS = 16)

```
$ PYTHONPATH=build/lsquant_core:python OMP_NUM_THREADS=16 python3 -u tests/dev/stage_4_check.py
(a) G7 chain1d_mag Gershgorin window = [-2.31, 2.31], BandWidth = 4.62
(a) G7 TimeProjSX-PX array_equal vs golden: True
(a) G7 TimeProjSY-PX array_equal vs golden: True
(a) G7 TimeProjSZ-PX array_equal vs golden: True
(b) G6 CorrelationVXchain1d array_equal vs golden: True
(b) G6 EF=+0.0 curve differing lines: 0 (total lines: 64)
(b) G6 EF=+0.5 curve differing lines: 0 (total lines: 64)
(b) G6 EF=+1.0 curve differing lines: 0 (total lines: 64)
(b) G6 EF=+1.5 curve differing lines: 0 (total lines: 64)
(c) compat GSL Bessel vs SciPy jv max_abs_diff = 5.551115e-16
(d) time_correlation_plan chain1d VX-VX (M=64, steps=16, tmax=5, exact) array_equal vs inline_compute-kpm-CorrOp: True
wall 655.12 s
exit: 0
```

Log: `figures/out/stage_4/stage_4_check.log`. The lines match the developer's report, which
ran at OMP_NUM_THREADS = 1. Item (d) ran the legacy binary from `build/` in this session.

### A.3 Gate (`scripts/check_goldens.sh 4 build`)

Harness extension (this handoff): stage 4 runs stages 1, 2 and 3 in full, with their
comparisons and greps. It then runs `python3 -m lsquant.apps.stage4_time` at
`OMP_NUM_THREADS=16` and again at `1`, with `PYTHONPATH=build/lsquant_core:python`:

1. **G6 moments.** Command: `--operators goldens/stage0/chain1d_ballistic/operators --label
   chain1d --velocity VX --window -2 2 --trace exact --states 256 --moments 256 --steps 64
   --tmax 20 --word msd`. The produced `.chebmomTD` is compared with `cmp`.
2. **G6 curves.** Command: `--moments <file> --broadening-meV 10 --fermi-energy 0 0.5 1.0
   1.5`. The curves are reconstructed twice: from the committed golden `.chebmomTD`, and from
   the `.chebmomTD` produced in step 1. Both sets are compared with
   `goldens/stage0/chain1d_ballistic/mean...EF{0.000000,0.500000,1.000000,1.500000}JACKSON.dat`.
   This also checks that `--fermi-energy` formats `E_F` the way `std::to_string` does.
3. **G7 moments.** Command: `--operators goldens/stage0/chain1d_mag/operators --label
   chain1d_mag --operator SX|SY|SZ --projector PX --window gershgorin --trace exact --states
   128 --moments 32 --steps 96 --tmax 95 --word projected`. The three `.chebmomTD` files are
   compared with `cmp`.

   About the window: the brief's section 6 resolves it as the legacy `SpectralBounds` without
   a `BOUNDS` file. The developer identifies that as branch 3, the Gershgorin enclosure with
   pad 0.10, which gives `[-2.31, 2.31]`, BandWidth 4.62, matching the golden header
   `128 4.62 0 96 1`. The gate therefore passes `--window gershgorin`, not the literal
   numbers.
4. **Grep.** `grep -rnwE "time|velocity|spin|projector|bessel|gsl" lsquant/core lsquant/io`
   matches whole lowercase words. Comment lines (first non-blank characters `//`, `*` or `/*`)
   are excluded from the verdict and listed. The identifiers `TimeAxis`, `EvolveRegister`,
   `EvolveDerived` and `ApplyDerivedSeries` cannot match: the pattern is case-sensitive and
   whole-word.

Command: `OMP_NUM_THREADS=16 bash scripts/check_goldens.sh 4 build`. Wall time: 2458 s. Log:
`figures/out/stage_4/check_goldens_stage4.log`. The stage 1 to 3 part (41 files and three
greps) is identical, line for line in its verdict lines, to the stage-3 final run.

#### A.3.1 Stage-4 part, verbatim (the "printed:" lines are shortened here; the full text is in the log)

```
== stage 4: python3 -m lsquant.apps.stage4_time, PYTHONPATH=/home/ICN2/jgarcia/lsquant-unified/build/lsquant_core:python, vs goldens/stage0/{chain1d_ballistic,chain1d_mag}
          [T=16] G6 msd moments wall 140 s
IDENTICAL [T=16] chain1d_ballistic/CorrelationVXchain1dKPM_M256_state20exact.chebmomTD (G6 moments)
IDENTICAL [T=16] chain1d_ballistic/meanCorrelationVXchain1dKPM_M256_state20exactEF0.000000JACKSON.dat (G6 curve from golden moments)
IDENTICAL [T=16] chain1d_ballistic/meanCorrelationVXchain1dKPM_M256_state20exactEF0.500000JACKSON.dat (G6 curve from golden moments)
IDENTICAL [T=16] chain1d_ballistic/meanCorrelationVXchain1dKPM_M256_state20exactEF1.000000JACKSON.dat (G6 curve from golden moments)
IDENTICAL [T=16] chain1d_ballistic/meanCorrelationVXchain1dKPM_M256_state20exactEF1.500000JACKSON.dat (G6 curve from golden moments)
IDENTICAL [T=16] chain1d_ballistic/meanCorrelationVXchain1dKPM_M256_state20exactEF0.000000JACKSON.dat (G6 curve from new moments)
IDENTICAL [T=16] chain1d_ballistic/meanCorrelationVXchain1dKPM_M256_state20exactEF0.500000JACKSON.dat (G6 curve from new moments)
IDENTICAL [T=16] chain1d_ballistic/meanCorrelationVXchain1dKPM_M256_state20exactEF1.000000JACKSON.dat (G6 curve from new moments)
IDENTICAL [T=16] chain1d_ballistic/meanCorrelationVXchain1dKPM_M256_state20exactEF1.500000JACKSON.dat (G6 curve from new moments)
          [T=16] G7 SX-PX moments wall 45 s
          [T=16] G7 SY-PX moments wall 52 s
          [T=16] G7 SZ-PX moments wall 36 s
IDENTICAL [T=16] chain1d_mag/TimeProjSX-PXchain1d_magKPM_M32_stateexact.chebmomTD (G7 moments)
IDENTICAL [T=16] chain1d_mag/TimeProjSY-PXchain1d_magKPM_M32_stateexact.chebmomTD (G7 moments)
IDENTICAL [T=16] chain1d_mag/TimeProjSZ-PXchain1d_magKPM_M32_stateexact.chebmomTD (G7 moments)
          [T=1] G6 msd moments wall 23 s
IDENTICAL [T=1] chain1d_ballistic/CorrelationVXchain1dKPM_M256_state20exact.chebmomTD (G6 moments)
IDENTICAL [T=1] chain1d_ballistic/meanCorrelationVXchain1dKPM_M256_state20exactEF0.000000JACKSON.dat (G6 curve from golden moments)
IDENTICAL [T=1] chain1d_ballistic/meanCorrelationVXchain1dKPM_M256_state20exactEF0.500000JACKSON.dat (G6 curve from golden moments)
IDENTICAL [T=1] chain1d_ballistic/meanCorrelationVXchain1dKPM_M256_state20exactEF1.000000JACKSON.dat (G6 curve from golden moments)
IDENTICAL [T=1] chain1d_ballistic/meanCorrelationVXchain1dKPM_M256_state20exactEF1.500000JACKSON.dat (G6 curve from golden moments)
IDENTICAL [T=1] chain1d_ballistic/meanCorrelationVXchain1dKPM_M256_state20exactEF0.000000JACKSON.dat (G6 curve from new moments)
IDENTICAL [T=1] chain1d_ballistic/meanCorrelationVXchain1dKPM_M256_state20exactEF0.500000JACKSON.dat (G6 curve from new moments)
IDENTICAL [T=1] chain1d_ballistic/meanCorrelationVXchain1dKPM_M256_state20exactEF1.000000JACKSON.dat (G6 curve from new moments)
IDENTICAL [T=1] chain1d_ballistic/meanCorrelationVXchain1dKPM_M256_state20exactEF1.500000JACKSON.dat (G6 curve from new moments)
          [T=1] G7 SX-PX moments wall 4 s
          [T=1] G7 SY-PX moments wall 4 s
          [T=1] G7 SZ-PX moments wall 4 s
IDENTICAL [T=1] chain1d_mag/TimeProjSX-PXchain1d_magKPM_M32_stateexact.chebmomTD (G7 moments)
IDENTICAL [T=1] chain1d_mag/TimeProjSY-PXchain1d_magKPM_M32_stateexact.chebmomTD (G7 moments)
IDENTICAL [T=1] chain1d_mag/TimeProjSZ-PXchain1d_magKPM_M32_stateexact.chebmomTD (G7 moments)
GREP      grep -rnwE "time|velocity|spin|projector|bessel|gsl" lsquant/core lsquant/io, code lines: empty
          excluded comment line: (27 comment lines, listed in the log; e.g. plan.hpp:7 "time-evolution step descriptors", plan.hpp:26 "... recursions such as mean-square", sparse_matrix.hpp:208 "purely imaginary velocity operators", and "O(N) time" complexity notes)
== 65 files compared; failures: 0
GATE PASS (stage 4)
gate exit: 0 (wall 2458 s)
```

Per-file count: stage 1 has 11, stage 2 has 6, stage 3 has 24, and stage 4 has 24
(12 at T=16 and 12 at T=1). Total 65, all IDENTICAL. Greps: `basis|std::cout|getenv` (core)
empty; `delta|DOS|Kubo|basis|std::cout|getenv` (core, io, code lines) empty; `compat/`
(core, io) empty; stage-4 physics words (core, io, code lines) empty. No `DIFFERS`, so the
localisation checklist of codex.md section 2 was not needed for the gate.

### A.4 Figure (`figures/stage_4.py`)

Command: `MPLBACKEND=Agg python3 figures/stage_4.py --build build --threads 16`.

The script runs `python3 -m lsquant.apps.tutorial04|tutorial02|tutorial05 --out <tmp>`, each
in its own temporary directory, timed, with `KPM_SEED=12345` in the environment. Each app
sets its own trace seeds, as the legacy tutorial scripts do.

Matching rule for the byte comparison: every golden `captured/NNNN_<name>` is compared
(`filecmp`, `shallow=False`) with the new `captured/NNNN_<name>`, that is, the same index and
the same name. Tutorial 02 has 16 captures per name (the 16 disorder realisations), and each
is compared with its own counterpart. If a capture differs, it is also checked against every
new capture of the same name, to tell an ordering difference from a value difference. The
golden top-level `*.dat` files are compared by name.

A failing app does not stop the script: its files count as `DIFFERS (not produced)` and its
panel shows the golden only, with the error printed in the residual axis. The script exits
1 if any file differs.

Outputs:
- `figures/out/stage_4/fig_stage_4.png`
- `fig_stage_4.caption.txt`
- `fig_stage_4.cmp.txt` (the per-file cmp list, 212 entries)
- run log `figure_run.log` (exit 1)

#### A.4.1 Production run (the deliverable)

```
tutorial04: wall 217.9 s, 3 lines printed, 3 of 3 IDENTICAL (3 new captures)
tutorial02: wall 1.0 s, 0 lines printed, 0 of 193 IDENTICAL (0 new captures); tutorial02 exited 1: tutorial02 error: module 'lsquant_tutorial02_make_disordered_chain' has no attribute 'build'
tutorial05: wall 1.1 s, 0 lines printed, 0 of 16 IDENTICAL (0 new captures); tutorial05 exited 1: tutorial05 error: module 'lsquant_tutorial05_make_square' has no attribute 'build'
max |new - golden|: (a) 0.0e+00  (b) n/a (no new data)  (c) n/a (no new data); (a) max |golden - analytic| 3.997e-15
exit: 1
```

Tutorial 04 cmp list: `captured/0000_TimeProjSX-PX...chebmomTD`, `0001_...SY...` and
`0002_...SZ...` are all IDENTICAL. Tutorial 02: all 192 captures and the top-level
`meanCorrelationVXchain1d_dis_N512_W1KPM_M256_state600defaultEF0.000000JACKSON.dat` are
DIFFERS (not produced). Tutorial 05: the 8 captures and the 8 top-level `.dat` files are
DIFFERS (not produced).

Caption (verbatim, `fig_stage_4.caption.txt`):

> FIG. 4. Stage 4 (handoff 0004): tutorials 04, 02 and 05 reproduced from Python by python3 -m lsquant.apps.tutorial04|tutorial02|tutorial05 (lsquant package over the pybind11 module lsquant._core, OMP_NUM_THREADS = 16), each run in a temporary directory; time evolution U(dt) = exp(-i H dt/hbar) by the Bessel-Chebyshev series (GSL J_n via lsquant._core.compat, |J_n| <= 1e-15 cut, hbar = 0.6582119624 eV fs), dt = tmax/(steps - 1). (a) Tutorial 04, Larmor precession on chain1d_mag (N = 128, t = -1, J_ex = 0.1 eV): m = 0 row of TimeProj{SX,SY,SZ}-PX (projector PX = (1 + S_x)/2, exact trace over the 128 site states, M = 32, 96 time steps, tmax = 95 fs, Gershgorin window [-2.31, 2.31] eV), normalised by <S_x(0)> as lsqspin.py does; new as filled markers, golden as open grey markers (every 4th step), cos(wt), -sin(wt), 0 with w = 2 J_ex/hbar as thin lines; residual new - analytic below (golden - analytic max 4.00e-15, round-off level); max |new - golden| = 0.0e+00. (b) Tutorial 02, localisation in the 1D Anderson chain (box disorder W = 1 eV, disorder seed and KPM_SEED = r for r = 1..16, one random-phase vector, window from BOUNDS, M = 256, 128 steps, tmax = 600 fs, MSD(E,t) at 10 meV): xi(E) = cal sqrt(<MSD>) over the last 30% of the time window, averaged over the 16 realisations, N = 256, 1024, 4096, E = 0, +-1, 1.9 eV, cal fixed once against the Born length 24(4 - E^2)/W^2 at N = 4096, E = 1 eV (lsqloc.py); new as lines, golden as open markers, Born dashed; residual new - golden below (max n/a (no new data) sites). (c) Tutorial 05, 2D square lattice L = 200 (N = 40000), disorder W = 0.5, 1, 2, 4 eV (seed 1, KPM_SEED = 1, one random-phase vector): D(E,t) = DeltaX^2/(2t) at E = 2 eV from the MSD (M = 256, 160 steps, tmax = 150 fs, 80 meV; the app also writes the DOS, M = 1024, 80 meV); new as lines, golden as open markers (every 8th step); residual new - golden below (max n/a (no new data)). Byte comparison with goldens/stage0/tutorials/0N (captured/NNNN_<name> against the capture with the same index and name, top-level .dat by name; cmp criterion): tutorial 04 3 of 3 IDENTICAL; tutorial 02 0 of 193 IDENTICAL (app failed: tutorial02 exited 1: tutorial02 error: module 'lsquant_tutorial02_make_disordered_chain' has no attribute 'build'); tutorial 05 0 of 16 IDENTICAL (app failed: tutorial05 exited 1: tutorial05 error: module 'lsquant_tutorial05_make_square' has no attribute 'build'). Wall times (TCN, load 60 to 100 on 64 cores): tutorial 04 217.9 s, tutorial 02 1.0 s (failed), tutorial 05 1.1 s (failed). AT LEAST ONE FILE DIFFERS OR IS NOT PRODUCED.

#### A.4.2 Diagnostic run with the generator entry point aliased (localisation of D1; not a gate result)

Purpose: find out whether D1 is the only thing wrong with `tutorial02` and `tutorial05`.

Method: a scratch wrapper, `/tmp/s4_shim/shim_run.py` on TCN, not committed and not part of
the harness. It imports the app module unchanged and wraps its `_load_make_*()` so that the
loaded generator module gains one alias:

- `build(dir, n, W, seed) = build_chain(n, W, seed, outdir=dir)`
  (`make_disordered_chain.py`)
- `build(dir, L, W, seed) = build_square(L, W, seed, outdir=dir)` (`make_square.py`)

It then calls `main(["--out", <dir>])`. Nothing else is changed. The outputs were fed back
through the same figure script with `--reuse 04=... --reuse 02=... --reuse 05=... --name
fig_stage_4_diag_alias`.

Runs:
- tutorial05 at OMP_NUM_THREADS = 16: wall 624.4 s, concurrent with the gate.
- tutorial02 at OMP_NUM_THREADS = 1: wall 94.9 s. A first attempt at 16 threads had done
  only 28 of 192 realisations after about 10 minutes and was stopped (see S1).
- tutorial04 at OMP_NUM_THREADS = 1, unaliased, for comparison: wall 9 s.

```
tutorial04: reused /tmp/s4_shim/out04_T1, 3 of 3 IDENTICAL
tutorial02: reused /tmp/s4_shim/out02, 193 of 193 IDENTICAL (192 new captures)
tutorial05: reused /tmp/s4_shim/out05, 16 of 16 IDENTICAL (8 new captures)
max |new - golden|: (a) 0.0e+00  (b) 0.0e+00  (c) 0.0e+00
exit: 0
```

Files: `figures/out/stage_4/fig_stage_4_diag_alias.png`, `.caption.txt`, `.cmp.txt` and
`fig_stage_4_diag_alias_run.log`. That caption says OMP_NUM_THREADS = 16 because the figure
script was run with `--threads 16`, but the reused tutorial 02 and 04 outputs came from
1-thread runs, as listed above. The physics is visible in this figure:
- (a) `<S_x>`, `<S_y>` and `<S_z>` follow cos, -sin and 0, with golden - analytic at about
  4e-15.
- (b) xi(E) follows the Born parabola for N = 1024 and 4096, while N = 256 is finite-size
  limited at E = 0.
- (c) D(E,t) goes ballistic, then diffusive, then localised as W grows.

All of it is byte-identical to the goldens.

### A.5 Defects found (localised; not fixed)

- **D1 (blocking for the tutorial deliverable): `tutorial02.py` and `tutorial05.py` call a
  generator function that does not exist.**
  - `python/lsquant/apps/tutorial02.py`, `_run_one_chain`, calls
    `make_mod.build(stage_dir, n_sites, w_val, disorder_seed)`.
  - `python/lsquant/apps/tutorial05.py`, `run_tutorial05`, calls
    `make_square.build(stage_dir, L_SIDE, w_val, DISORDER_SEED)`.
  - The example generators on `cpp-unified` define only
    `build_chain(N, W, seed, outdir=".")` (`examples/02_localization_1d_chain/make_disordered_chain.py:108`)
    and `build_square(L, W, seed, outdir=".")` (`examples/05_charge_transport_square/make_square.py:100`).
    Both take the directory as the keyword argument `outdir`, last, and return the label.
  - Both apps exit 1 with `AttributeError ... has no attribute 'build'` before any
    computation.
  - Section 5.2 shows that with the call alone corrected, both apps are byte-identical to
    the goldens: 193 of 193 and 16 of 16.
  - The developer's report shows only a `tutorial04` smoke run, so tutorials 02 and 05 were
    apparently never executed.

