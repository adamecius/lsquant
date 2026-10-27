# Handoff 0006 — tester report (stage 6 gate: rotated models reproduce G1, G2, G5)

Tester: Codex role (Opus subagent driven by Fable). Machine: TCN (`fatnode-tacn`),
checkout `~/lsquant-unified`, branch `handoff/0006-codex` cut from `cpp-unified`
(`552de6e`) and merged with `origin/handoff/0006-agi` at `a5108e3` (fast-forward: the
developer branch already contains `552de6e`). Nothing under `lsquant/`, `python/lsquant/`,
`tests/dev/`, the legacy tree, `goldens/` (data), a seed or a tolerance was edited. Files
edited or added: `scripts/check_goldens.sh` (stage 6), `figures/stage_6.py`,
`figures/stage_6_fixtures.py` (section 1 verification), `figures/out/stage_6/*`,
`goldens/stage0/PROVENANCE.md` (a G10 paragraph), this file.

**Closure run (section 6, at the end; supersedes sections 2 to 4 where they differ):**
`origin/handoff/0006-agi` at `8305f11` (review 2, plus `cpp-unified` with stages 4 and 5),
merged into `handoff/0006-codex`. Stages 1 to 5: 103 files IDENTICAL, 52 stage-5 rows PASS,
5 greps empty, 0 failures. Stage 6: 74 of 74 rows PASS (40 vs goldens, 24 route vs route,
10 Hall plateaus), C = +1.017017 on every route and applier. Naive evidence: DOS PASS
(4.7e-16, 1.0e-15), Kubo FAIL by O(1) (0.17 to 0.20; Hall C = 1.19, 1.20). 1-thread
repeat of the 24 chain1d runs: 20 IDENTICAL, **4 DIFFERS, all `lowdin/exact`** (rounding
level, 3e-16 to 1.1e-14; confirmed on a second machine). Harness exit: `GATE FAIL
(stage 6only)`, 98 compared, 4 failures, these 4 thread-count rows.

**First run, result (sections 1 to 5, at `a5108e3`).** Fixtures: 4 of 4 ACCEPTED. Stage-6 rows: 74 of 74 PASS (40 route/applier rows
vs the goldens, 24 route-vs-route rows, 10 Hall plateaus); 8 `lowdin/ischo` rows N/A; 10
`naive` rows reported as evidence. Stages 1 to 3 inside the same run: 41 of 41 files
IDENTICAL, the stage-2 comment-aware grep and the stage-3 `compat/` grep empty. **The
harness exit is nevertheless `GATE FAIL (stage 6)`, 115 compared, 1 failure: the stage-1
grep `grep -rn "basis\|std::cout\|getenv" lsquant/core` (no comment exclusion) finds the
word `basis` in 9 Doxygen comment lines of the new `metric_policy.hpp` and
`composite_operators.hpp` (section 2.4).** No physics defect was found. C = +1.017017 for
every route and applier, equal to the golden G5 value.

## 0. Build

```
$ cmake -S . -B build6 -DLSQUANT_BUILD_NEW_CORE=ON -DLSQUANT_BUILD_PYTHON=ON -Dpybind11_DIR="$(python3 -c 'import pybind11; print(pybind11.get_cmake_dir())')"
$ cmake --build build6 -j 24
BUILD_EXIT 0
```

Separate build directory `build6` so that other tester sessions using `build/` are not
disturbed. Full log: `figures/out/stage_6/build.log`.

## 1. Fixture verification (`figures/stage_6_fixtures.py`, independent of the generator)

Own reader for the legacy CSR text layout, dense numpy/scipy (`eigvalsh`, `scipy.linalg.eigh(H', S')`),
originals from `goldens/stage0/{chain1d,haldane}/operators`. Log:
`figures/out/stage_6/fixture_verification.log`.

| fixture | N | sha256 (all files vs manifest, and file list = manifest list) | max \|S - S^+\| | max \|H' - H'^+\| | max \|V' - V'^+\| | lambda_min(S) | lambda_max(S) | cond(S) | rel. diff vs manifest cond | max \|eig(H',S') - eig(H)\| | extra: max \|eig(H'+0.37VX'+0.21VY', S') - eig(same, orth.)\| | verdict |
|---|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|
| chain1d/real | 2048 | all match | 0 | 0 | 0 | 0.31377193275553 | 2.10355392820201 | 6.704085702404 | 2.0e-14 | 1.3e-14 | 5.2e-15 | ACCEPTED |
| chain1d/complex | 2048 | all match | 0 | 0 | 0 | 0.36043374446001 | 2.02063482826285 | 5.606120013236 | 1.7e-15 | 6.7e-15 | 6.2e-15 | ACCEPTED |
| haldane/real | 512 | all match | 0 | 0 | 0 | 0.35741617850047 | 2.06916590808509 | 5.789234042975 | 1.5e-15 | 1.0e-14 | 5.8e-15 | ACCEPTED |
| haldane/complex | 512 | all match | 0 | 0 | 0 | 0.38717798606128 | 2.02066370188938 | 5.218952974174 | 8.8e-15 | 1.5e-14 | 4.9e-15 | ACCEPTED |

All checks of codex.md section 1 hold: S Hermitian and positive definite, cond(S) equal
to the manifest to 1e-6 (actually 2e-14) and inside [5, 10], the generalised spectrum
equal to the orthogonal one to 1e-10 (actually 1.5e-14), H' and V' Hermitian to 1e-14
(exactly), sha256 of each file equal to the manifest. The extra row (tester's addition,
cheap) checks that every operator was rotated with the same T: the generalised spectrum
of a fixed combination of H', VX', VY' equals that of H, VX, VY. The `complex` variants
really are complex (max |Im S| = 0.50 chain1d, 0.25 Haldane; the `real` S has Im = 0).
Note: the manifest names the strength `alpha` (0.55), the brief `epsilon`.

## 2. The gate (`scripts/check_goldens.sh 6 build6`)

Command (TCN, `KPM_SEED=12345`, stage-6 runs at `OMP_NUM_THREADS=16`):

```
S6_KEEP=$HOME/lsquant-unified/figures/out/stage_6/data scripts/check_goldens.sh 6 build6
```

Stage 6 runs stages 1, 2 and 3 in full first (decision: stage 3 has landed and G4, G5
must stay unchanged), then the stage-6 block. Full log:
`figures/out/stage_6/check_goldens_stage6.log`; tables: `figures/out/stage_6/data/gate_table.tsv`,
`figures/out/stage_6/data/chern.tsv` (the run directories under `data/runs` are 37 MB and are not committed).

### 2.1 Decisions recorded

1. **Measure.** These runs cannot be byte-identical, so the criterion is the brief's
   tolerance (section 6). The measure is `rel = max_k |mu_new[k] - mu_gold[k]| /
   max_k |mu_gold[k]|` over all complex moments of the file (the moment vector's max-norm
   relative difference), plus header lines (N, BandWidth, BandCenter; M) numerically
   identical. `scripts/golden_diff.py` is still printed under every row; its
   per-element relative difference is ~1 on moments that vanish by symmetry (odd DOS
   moments, 1e-17 vs 0) and cannot be the criterion.
2. **Applier parameters.** `exact`: defaults. `ischo`: `--tol 1e-12` (brief table),
   `max_iterations` 200 (program default). `chebyshev`: `--order 60` (program default) with
   the constructor truncation `coeffTol = 1e-12` (header default; `stage6_nonorth` does not
   expose it); the effective order is printed.
3. **`lowdin/ischo` is N/A**: ISCHO provides `S^-1` only, the Löwdin route needs `S^-1/2`;
   `stage6_nonorth` rejects the combination with `std::invalid_argument` (developer report,
   section 5). Not counted.
4. **Route vs route**: `generalised/exact vs lowdin/exact` (1e-10),
   `generalised/chebyshev vs lowdin/chebyshev` (1e-8), `generalised/ischo vs lowdin/exact`
   (1e-9, as there is no `lowdin/ischo`).
5. **`naive`** run with `exact` only (the applier is not used on this route), evidence, not counted.
6. **Hall plateau**: `python3 -m lsquant.apps.stage3_kubo --moments <VX-VY file> --formula
   bastin --broadening-meV 10`, then `C = (mean_{|E|<0.5} sigma_xy / A) 2 pi`,
   `A = (512/2) 3 sqrt(3)/2` (identical to `tutorial03.chern_from_plateau` and
   `test/hall_haldane.sh`), PASS if `|C - 1| <= 0.05`.
7. **Environment workaround.** On TCN another session runs a loop
   `while true; do pkill -f "stage6_nonorth.*--route lowdi[n]"; sleep 20; done` (seen in
   `ps`, not mine, left alone). It would kill every Löwdin run of the gate. The harness
   therefore calls the same binary through a symlink `$TMP/s6_nonorth_exe`, so the command
   line does not match. The first gate launch was restarted for this before it reached stage 6.
8. A stage-6-only copy of the harness (`case 6) stage6`) was run in parallel in `~/s6`
   while the stage-1 part (64 threads, machine load ~70) was slow; its table is identical to
   the official run's (`diff` of the first eight columns empty). Only the official run is reported.

### 2.2 Stages 1 to 3 (inside the stage-6 run)

41 of 41 files IDENTICAL: G1, G2 at 1, 8, 64 threads, batched G2, random-phase
determinism, G3 forms (a) and (b) at 1 and 16 threads, G5 moments and curves, G4 moments
and curves at 16 and 1 threads. Stage-2 grep (code lines): empty. Stage-3 grep
`compat/`: empty. **Stage-1 grep: FOUND** (section 2.4).

### 2.3 Stage-6 table

`rel` as defined in 2.1; residual = `residual=` of the stdout line (`||S S^-1 v - v||/||v||`
on one random vector, or for ISCHO the `LastError()` of the last solve); it. = ISCHO
`LastIterations()` of the last solve; eff. = Chebyshev effective order.

| model | variant | word | route/applier | rel vs golden | tol | verdict | residual | it. | eff. |
|---|---|---|---|---:|---:|---|---:|---:|---:|
| chain1d | real | dos | lowdin/exact | 7.960e-13 | 1e-10 | PASS | 2.64e-14 | | |
| chain1d | real | dos | lowdin/chebyshev | 9.507e-12 | 1e-8 | PASS | 1.77e-12 | | 38 |
| chain1d | real | dos | generalised/exact | 4.664e-16 | 1e-10 | PASS | 3.74e-16 | | |
| chain1d | real | dos | generalised/ischo | 3.269e-14 | 1e-9 | PASS | 3.52e-14 | 0 | |
| chain1d | real | dos | generalised/chebyshev | 1.517e-12 | 1e-8 | PASS | 1.65e-12 | | 41 |
| chain1d | real | kubo | lowdin/exact | 1.510e-14 | 1e-10 | PASS | 2.64e-14 | | |
| chain1d | real | kubo | lowdin/chebyshev | 1.675e-13 | 1e-8 | PASS | 1.77e-12 | | 38 |
| chain1d | real | kubo | generalised/exact | 2.442e-15 | 1e-10 | PASS | 3.74e-16 | | |
| chain1d | real | kubo | generalised/ischo | 2.838e-15 | 1e-9 | PASS | 2.01e-16 | 0 | |
| chain1d | real | kubo | generalised/chebyshev | 1.026e-13 | 1e-8 | PASS | 1.65e-12 | | 41 |
| chain1d | complex | dos | lowdin/exact | 6.203e-13 | 1e-10 | PASS | 2.79e-14 | | |
| chain1d | complex | dos | lowdin/chebyshev | 9.337e-13 | 1e-8 | PASS | 1.76e-12 | | 33 |
| chain1d | complex | dos | generalised/exact | 1.022e-15 | 1e-10 | PASS | 4.01e-16 | | |
| chain1d | complex | dos | generalised/ischo | 2.735e-14 | 1e-9 | PASS | 1.46e-14 | 0 | |
| chain1d | complex | dos | generalised/chebyshev | 2.916e-12 | 1e-8 | PASS | 1.13e-12 | | 36 |
| chain1d | complex | kubo | lowdin/exact | 1.243e-14 | 1e-10 | PASS | 2.79e-14 | | |
| chain1d | complex | kubo | lowdin/chebyshev | 6.750e-14 | 1e-8 | PASS | 1.76e-12 | | 33 |
| chain1d | complex | kubo | generalised/exact | 2.554e-15 | 1e-10 | PASS | 4.01e-16 | | |
| chain1d | complex | kubo | generalised/ischo | 7.466e-15 | 1e-9 | PASS | 1.91e-16 | 0 | |
| chain1d | complex | kubo | generalised/chebyshev | 8.260e-14 | 1e-8 | PASS | 1.13e-12 | | 36 |
| haldane | real | kuboxy | lowdin/exact | 6.041e-14 | 1e-10 | PASS | 1.98e-14 | | |
| haldane | real | kuboxy | lowdin/chebyshev | 7.662e-13 | 1e-8 | PASS | 1.55e-12 | | 34 |
| haldane | real | kuboxy | generalised/exact | 1.340e-14 | 1e-10 | PASS | 1.04e-15 | | |
| haldane | real | kuboxy | generalised/ischo | 3.116e-13 | 1e-9 | PASS | 9.82e-13 | 6 | |
| haldane | real | kuboxy | generalised/chebyshev | 5.878e-13 | 1e-8 | PASS | 1.05e-12 | | 37 |
| haldane | real | kubo | lowdin/exact | 8.285e-14 | 1e-10 | PASS | 1.98e-14 | | |
| haldane | real | kubo | lowdin/chebyshev | 7.667e-13 | 1e-8 | PASS | 1.55e-12 | | 34 |
| haldane | real | kubo | generalised/exact | 1.148e-14 | 1e-10 | PASS | 1.04e-15 | | |
| haldane | real | kubo | generalised/ischo | 3.854e-13 | 1e-9 | PASS | 9.67e-13 | 6 | |
| haldane | real | kubo | generalised/chebyshev | 7.574e-13 | 1e-8 | PASS | 1.05e-12 | | 37 |
| haldane | complex | kuboxy | lowdin/exact | 4.036e-14 | 1e-10 | PASS | 1.52e-14 | | |
| haldane | complex | kuboxy | lowdin/chebyshev | 6.433e-13 | 1e-8 | PASS | 2.21e-12 | | 31 |
| haldane | complex | kuboxy | generalised/exact | 1.340e-14 | 1e-10 | PASS | 1.03e-15 | | |
| haldane | complex | kuboxy | generalised/ischo | 5.101e-14 | 1e-9 | PASS | 8.56e-14 | 6 | |
| haldane | complex | kuboxy | generalised/chebyshev | 3.179e-13 | 1e-8 | PASS | 1.14e-12 | | 34 |
| haldane | complex | kubo | lowdin/exact | 5.196e-14 | 1e-10 | PASS | 1.52e-14 | | |
| haldane | complex | kubo | lowdin/chebyshev | 6.069e-13 | 1e-8 | PASS | 2.21e-12 | | 31 |
| haldane | complex | kubo | generalised/exact | 1.148e-14 | 1e-10 | PASS | 1.03e-15 | | |
| haldane | complex | kubo | generalised/ischo | 6.636e-14 | 1e-9 | PASS | 9.98e-14 | 6 | |
| haldane | complex | kubo | generalised/chebyshev | 4.133e-13 | 1e-8 | PASS | 1.14e-12 | | 34 |

`lowdin/ischo`: N/A in all 8 model/variant/word cells (2.1, item 3).

Worst case per route/applier over all cells: lowdin/exact 7.96e-13, lowdin/chebyshev
9.51e-12, generalised/exact 1.34e-14, generalised/ischo 3.85e-13, generalised/chebyshev
2.92e-12: two to four orders of magnitude inside the tolerances. On the generalised/exact
route the agreement is at round-off (the dense `LLT` inverse); the Löwdin/exact route sits
at ~1e-12 on the chain DOS (the dense eigen-decomposition `S^-1/2`, error growing with m,
max at m = 62, 63).

Route vs route (24 rows, all PASS): generalised/exact vs lowdin/exact max 7.96e-13
(chain1d/real dos); generalised/chebyshev vs lowdin/chebyshev max 7.99e-12 (chain1d/real
dos); generalised/ischo vs lowdin/exact max 8.29e-13 (chain1d/real dos). Every row is in
the log and in `gate_table.tsv`.

Naive route (evidence, not counted): chain1d dos 1.1e+27 (real) / 2.3e+23 (complex), kubo
4.7e+49 / 2.9e+47 (the recursion diverges); Haldane kuboxy 0.85 / 0.86, kubo 1.00 / 1.00.
All would FAIL. The brief expected the naive route to pass on `dos`; it does not, because
the implemented `naive` is not the brief's naive (smell 5).

### 2.4 The one failure: stage-1 grep on comment lines

```
GREP      grep -rn "basis\|std::cout\|getenv" lsquant/core: FOUND
          lsquant/core/algebra/metric_policy.hpp:23,35,45,47,99,157,187  (Doxygen comments: "basis non-orthogonality", "orthonormal basis", ...)
          lsquant/core/algebra/composite_operators.hpp:321,382            (Doxygen comments: "non-orthogonal basis")
```

On `cpp-unified` the same grep is empty (0 lines). All 9 hits are comment lines; the
stage-2 grep, which excludes comment lines, is empty. Not relaxed here: the stage-1 grep
is as Fable wrote it. Either the comments say "metric" / "overlap" instead of "basis", or
Fable decides the stage-1 grep should exclude comments as the stage-2 one does.

### 2.5 Hall plateau (stage-3 Bastin, 10 meV)

| source | C | verdict |
|---|---:|---|
| golden G5 Bastin (legacy, reference) | +1.017017 | REF |
| haldane/real lowdin/exact, lowdin/chebyshev, generalised/exact, generalised/ischo, generalised/chebyshev | +1.017017 (each) | PASS (5) |
| haldane/complex lowdin/exact, lowdin/chebyshev, generalised/exact, generalised/ischo, generalised/chebyshev | +1.017017 (each) | PASS (5) |
| haldane/real naive | +1.015324 | evidence |
| haldane/complex naive | +1.017298 | evidence |

plateau = 107.656458 (450 points with |E| < 0.5) for all ten routes. The naive C is close
to 1 only by accident (its moments differ by O(1) from G5).

## 3. Figure

`figures/out/stage_6/fig_stage_6.png`, produced by
`MPLBACKEND=Agg OMP_NUM_THREADS=16 python3 figures/stage_6.py --build build6` from the kept gate runs (exit 0).
Caption (`figures/out/stage_6/fig_stage_6.caption.txt`):

> FIG. 6. Stage 6 (handoff 0006): a non-orthogonal basis through MetricPolicy reproduces the orthogonal goldens. Fixtures: the committed chain1d and Haldane operators rotated by a local transform T (H' = T^+ H T, S' = T^+ T, V' = T^+ V T; cond(S') 5.2 to 6.7; variants real and complex T), goldens/stage0/nonorth. Moments from build/lsquant_core/stage6_nonorth (exact trace, OMP_NUM_THREADS = 16). (a) chain1d DOS per site, window [-2, 2] eV, M = 64, reconstructed with the stage-2 lsquant.reconstruct.dos_curve (Jackson, 120 meV): orthogonal golden G1 (open markers, 40 points), Lowdin route (S^-1/2 H' S^-1/2) and generalised route (S^-1 (a H' + b S')) with the exact applier on the complex fixture; residual route - G1 below for both fixtures, max |residual| = 3.9e-14 1/eV. (b) Haldane sigma_xy(E), window [-4.5, 4.5] eV, M = 128, Kubo-Bastin at 10 meV by python3 -m lsquant.apps.stage3_kubo, in units of e^2/h per (N/2) A_cell (A_cell = 3 sqrt(3)/2) so the gap plateau (shaded, |E| < 0.5 eV) is C; golden G5 as markers; residual below, max |residual| = 9.4e-11. Chern numbers C = (mean_{|E|<0.5} sigma_xy / A) 2 pi: golden_G5 C = 1.017017 (REF); haldane_complex_kuboxy_generalised_chebyshev C = 1.017017 (PASS); haldane_complex_kuboxy_generalised_exact C = 1.017017 (PASS); haldane_complex_kuboxy_generalised_ischo C = 1.017017 (PASS); haldane_complex_kuboxy_lowdin_chebyshev C = 1.017017 (PASS); haldane_complex_kuboxy_lowdin_exact C = 1.017017 (PASS); haldane_complex_kuboxy_naive_exact C = 1.017298 (EVIDENCE(PASS)); haldane_real_kuboxy_generalised_chebyshev C = 1.017017 (PASS); haldane_real_kuboxy_generalised_exact C = 1.017017 (PASS); haldane_real_kuboxy_generalised_ischo C = 1.017017 (PASS); haldane_real_kuboxy_lowdin_chebyshev C = 1.017017 (PASS); haldane_real_kuboxy_lowdin_exact C = 1.017017 (PASS); haldane_real_kuboxy_naive_exact C = 1.015324 (EVIDENCE(PASS)). Max relative moment difference vs G1, G2, G5 (max_k |mu - mu_gold| / max_k |mu_gold|, over models, variants and words) per route/applier: generalised/chebyshev 2.92e-12; generalised/exact 1.34e-14; generalised/ischo 3.85e-13; lowdin/chebyshev 9.51e-12; lowdin/exact 7.96e-13 (tolerances: exact 1e-10, ischo 1e-9, chebyshev 1e-8; Lowdin/ischo not defined). Generalised vs Lowdin, same applier: chebyshev 7.99e-12; exact 7.96e-13; ischo 8.29e-13. Naive route (legacy estimator, evidence only): chain1d/real dos 1.149e+27; chain1d/real kubo 4.664e+49; chain1d/complex dos 2.306e+23; chain1d/complex kubo 2.914e+47; haldane/real kuboxy 8.497e-01; haldane/real kubo 9.956e-01; haldane/complex kuboxy 8.595e-01; haldane/complex kubo 1.002e+00. Gate rows: 64 PASS, 0 FAIL.

("Gate rows" in the caption counts the 40 + 24 moment rows of `gate_table.tsv`; the 10
Hall rows are listed separately in the same caption.) The sigma_xy residual of 9.4e-11
(in e^2/h per cell area) is concentrated in two spikes at |E| ~ 2.9 eV, where the Bastin
integrand is large; it does not touch the plateau.

## 4. Smells

1. **Hidden tolerance in a `.cpp`**: `overlap_appliers.cpp:134`
   `SparseMatrix CompressDenseMatrix(const Eigen::MatrixXcd& dense, double dropThreshold = 1e-18)`:
   the dense `S^-1` and `S^-1/2` are pruned at 1e-18 inside the `.cpp`, not a constructor
   parameter (reported by the developer as deviation 2; brief section 2 forbids it).
2. **Defaults on tolerances and caps in the headers**, not in the brief's interfaces:
   `ExactInverseOverlap(S, maxDimension = 4096)`, `ExactInverseSqrtOverlap(S, maxDimension = 4096)`,
   `IschoInverseOverlap(S, tolerance = 1e-13, maxIterations = 200)`,
   `ChebyshevInverseOverlap(..., coeffTol = 1e-12)`, `ChebyshevInverseSqrtOverlap(..., coeffTol = 1e-12)`.
   `stage6_nonorth` never passes `maxDimension` or `coeffTol` (no flag), so the gate's
   Chebyshev rows rely on the header default 1e-12.
3. **A fifth applier class** not in the brief: `ExactInverseSqrtOverlap` (dense
   eigen-decomposition), used by `lowdin/exact`. Consequence: there is no `S^-1/2` for
   ISCHO and `lowdin/ischo` is undefined (N/A rows).
4. **ISCHO diagnostics are last-call state written from a `const Multiply`**
   (`mutable lastIterations_, lastError_`): the `iterations=` and `residual=` printed by the
   gate program are those of the last solve, not a maximum; if the engine calls the
   operator from several OpenMP threads these members are written concurrently (a data
   race on diagnostics; results unaffected here). `iterations=0` on chain1d means the
   incomplete Cholesky of the banded chain `S` is already the exact Cholesky. The ISCHO
   non-convergence `throw` (commit `28d4def`) is raised inside `Multiply`, possibly inside an
   OpenMP region, where an exception cannot propagate.
5. **`--route naive` is not the brief's naive.** Brief section 4: naive = the symmetric
   plan with `S^-1 V'` on both sides (the legacy estimator, on the generalised H). The
   program's naive uses `OrthogonalMetric` and `RescaledOperator(H', window)`: it ignores
   `S` entirely, so the recursion runs on `H'`, whose spectrum leaves the window, and the
   chain1d moments diverge (1e+23 to 1e+49) instead of "passing on dos". The row does not
   demonstrate the legacy defect it is meant to show.
6. **Dead option**: `stage6_nonorth --drop-tol` (default 1e-3) is parsed and never used.
7. **Window padding constants in the app**: `ComputePaddedOverlapBounds` (Lanczos 80
   steps, seed 12345, pad 5% of the span, floor at half the lower Ritz value): in the app,
   not in core; acceptable but not in the brief.
8. **No applier prints** (`grep cout|cerr|printf|getenv|fopen|ifstream` on the four new
   or changed core `.cpp` files: empty); **no file reading in core**; the word `metric`
   does not appear under `core/chebyshev` (only `GeneralisedRescaledOperator` was added
   there, plus the brief-mandated `RescaledOperator(const LinearOperator&, window)`
   constructor; the stage-1 `SparseMatrix` member became `std::optional<SparseMatrix>`,
   and G1/G2 stay byte-identical). **The engine is untouched**: `git diff
   cpp-unified...a5108e3` changes no file of `MomentEngine`, `Plan`, or `io/`.
9. Stage-1 grep hits on comment lines (section 2.4).
10. Environment, not code: the foreign `pkill -f "stage6_nonorth.*--route lowdin"` loop
    on TCN (2.1 item 7). Whoever owns it should know it also kills gate runs.

## 5. Files

- `scripts/check_goldens.sh`: stage 6 (`6) stage1; core_grep; stage2; stage2_grep; stage3; stage3_grep; stage6`).
- `figures/stage_6.py`, `figures/stage_6_fixtures.py`.
- `figures/out/stage_6/`: `build.log`, `check_goldens_stage6.log`, `fixture_verification.log`,
  `fig_stage_6.png`, `fig_stage_6.caption.txt`, `stage_6_figure.log`,
  `data/gate_table.tsv`, `data/chern.tsv`.
- `goldens/stage0/PROVENANCE.md`: G10 paragraph (fixtures verified, not regenerated).

## 6. Closure run (`8305f11`)

### 6.1 Branch, build, self-check

- `handoff/0006-codex`: merged `c05fd7b` (review 2, as `bace579`), then `8305f11` (the
  developer's merge `092da48` of `origin/cpp-unified` with stages 4 and 5, plus the report)
  as `4fb389c`. Two conflicts, resolved by keeping both sides: `scripts/check_goldens.sh`
  (the stage-4/5 sections from `cpp-unified`, plus my stage-6 section; the stage-6 case now
  runs stages 1 to 5 in full, then stage 6) and `goldens/stage0/PROVENANCE.md` (the
  `offcentre_axis` bullet of `cpp-unified`, then the G10 paragraph). `85046a7` adds a
  `6only` case (stage-6 part alone, as stage 5 has `5only`) and keeps the first run's
  outputs under `*_run1_a5108e3` names.
- Fresh `build6` (TCN, g++ 11.4): `BUILD_EXIT 0` (`figures/out/stage_6/build.log`).
- Self-check `tests/dev/stage_6_check.py --exe build6/lsquant_core/stage6_nonorth`, 16
  threads: `Stage 6 Self-Check PASSED`, wall 47 min (`figures/out/stage_6/stage_6_check.log`).
  Its first launch failed: the script defaults to `build/lsquant_core/stage6_nonorth`
  (smell 11); I passed `--exe`.
- Harness changes for review 2:
  - `--max-dim 4096` (exact), `--tol 1e-12 --max-iter 200` (ischo) and
    `--order 60 --coeff-tol 1e-12` (chebyshev) are now explicit on every command line.
    The program echoes them as `coeff_tol=`, `max_dim=`, `max_iterations=`; the fields the
    harness reads (`residual=`, `iterations=`, `effective_order=`) are unchanged.
  - The symlink workaround is gone.
  - New: the chain1d runs are repeated at `OMP_NUM_THREADS=1` and `cmp`'d with the
    16-thread runs. Haldane is not repeated, to save time.

### 6.2 Runs and wall times

| run | what | result | wall |
|---|---|---|---|
| run 1 of the closure | `check_goldens.sh 6 build6`, stages 1 to 5 complete, stage 6 up to chain1d/complex dos | killed by SIGTERM (`GATE_EXIT 143`); stages 1 to 5 complete and all clean | 8078 s until killed |
| run 2 of the closure | `check_goldens.sh 6only build6` (stage 6 only, plus the 1-thread chain1d repeat) | 98 compared, 4 failures (section 6.7) | 37702 s (10.5 h) |
| Aron cross-check | 24 chain1d runs at T = 16 (one after another) and T = 1 (24 processes, one thread each) | 20 IDENTICAL, 4 DIFFERS | 7001 s |

- **The SIGTERM in run 1** was not mine. Another session's (the stage-7 tester's) `pkill`
  matched my gate by program name; the coordinator confirmed it. Run 1's log is
  `figures/out/stage_6/check_goldens_stage6_stages1to5.log`.
- **Run 2 is on the same build.** It reruns only the stage-6 part, because stages 1 to 5
  had already finished on this build.
- **The orphaned subshell.** After the kill, one stage-6 subshell of run 1 was still
  running. I stopped it by its PID (4000895, 4000896), after checking from its command
  line that it was my own process.
- **Slow 1-thread repeat.** With the fat node under load (load average 35 to 70), each
  1-thread chain1d run took about 10 minutes. The lowdin/exact runs, with a dense
  2048 x 2048 `S^-1/2`, were the slowest.

### 6.3 Stages 1 to 5 (closure run 1)

| stage | files IDENTICAL | other rows | grep |
|---|---:|---|---|
| 1 | 11 | | `basis\|std::cout\|getenv` on `lsquant/core`: **empty** (the 9 comment lines are reworded) |
| 2 | 6 | | code lines: empty |
| 3 | 30 | | `compat/`: empty |
| 4 | 24 | | `time\|velocity\|spin\|projector\|bessel\|gsl`, code lines: empty |
| 5 | 32 | 52 PASS, 0 FAIL (198 INFO, report only) | `fftw`: empty |

Total: 103 IDENTICAL, 0 DIFFERS, 52 PASS, 0 FAIL.

### 6.4 Stage 6: rotated fixtures against G1, G2, G5, both routes

Fixtures unchanged since the first run, so the sha256 values of section 1 still hold.
`rel` = max|mu - mu_gold| / max|mu_gold|. Each cell gives `rel`, then the stdout `residual`
in brackets, then the ISCHO iteration count (`it`) or the Chebyshev effective order
(`eff`).

| fixture, word | lowdin/exact | lowdin/chebyshev | generalised/exact | generalised/ischo | generalised/chebyshev |
|---|---|---|---|---|---|
| chain1d/real dos | 7.96e-13 (2.6e-14) | 9.51e-12 (1.8e-12; eff 38) | 4.66e-16 (3.7e-16) | 3.27e-14 (3.5e-14; it 0) | 1.52e-12 (1.7e-12; eff 41) |
| chain1d/real kubo | 1.51e-14 | 1.68e-13 | 2.44e-15 | 2.84e-15 (it 0) | 1.03e-13 |
| chain1d/complex dos | 6.20e-13 (2.8e-14) | 9.34e-13 (1.8e-12; eff 33) | 1.02e-15 (4.0e-16) | 2.74e-14 (1.5e-14; it 0) | 2.92e-12 (1.1e-12; eff 36) |
| chain1d/complex kubo | 1.24e-14 | 6.75e-14 | 2.55e-15 | 7.47e-15 (it 0) | 8.26e-14 |
| haldane/real kuboxy | 6.04e-14 (2.0e-14) | 7.66e-13 (1.6e-12; eff 34) | 1.34e-14 (1.0e-15) | 3.12e-13 (9.8e-13; it 6) | 5.88e-13 (1.0e-12; eff 37) |
| haldane/real kubo | 8.29e-14 | 7.67e-13 | 1.15e-14 | 3.85e-13 (it 6) | 7.57e-13 |
| haldane/complex kuboxy | 4.04e-14 (1.5e-14) | 6.43e-13 (2.2e-12; eff 31) | 1.34e-14 (1.0e-15) | 5.10e-14 (8.6e-14; it 6) | 3.18e-13 (1.1e-12; eff 34) |
| haldane/complex kubo | 5.20e-14 | 6.07e-13 | 1.15e-14 | 6.64e-14 (it 6) | 4.13e-13 |

- **All 40 PASS** (tolerances: exact 1e-10, ischo 1e-9, chebyshev 1e-8). The values equal
  the first run's to the last printed digit.
- **lowdin/ischo is N/A** in all 8 rows.
- **Route against route: 24 of 24 PASS.** The worst cases are exact 7.96e-13,
  chebyshev 7.99e-12, and generalised/ischo against lowdin/exact 8.29e-13.

### 6.5 The generalised-route asymmetric bra, against the naive (symmetric) bra

The generalised plan carries the adjoint explicitly on the bra:
`ApplyOperator(S^-1), ApplyOperator(S^-1 V'), Expand, ApplyOperator(S)`.
The review-2 `naive` route is the brief's legacy estimator: the same recursion on
`S^-1 H~` and the same `S^-1 V'` on the ket, but the symmetric bra (`S^-1 V'` and `Expand`,
with no `S^-1` and no `S`). The two routes differ only in the bra, so this pair is the
asymmetric-bra check:

| fixture, word | generalised/exact (asymmetric bra) | naive/exact (symmetric bra), evidence |
|---|---:|---:|
| chain1d/real dos | 4.664e-16 | 4.664e-16 (same plan: DOS PASS) |
| chain1d/complex dos | 1.022e-15 | 1.022e-15 (DOS PASS) |
| chain1d/real kubo | 2.44e-15 | **1.989e-01** (FAIL, O(1)) |
| chain1d/complex kubo | 2.55e-15 | **1.809e-01** |
| haldane/real kuboxy | 1.34e-14 | **1.772e-01** |
| haldane/real kubo | 1.15e-14 | **1.847e-01** |
| haldane/complex kuboxy | 1.34e-14 | **1.745e-01** |
| haldane/complex kubo | 1.15e-14 | **1.765e-01** |
| Hall C (real, complex) | +1.017017, +1.017017 | **+1.197943, +1.193531** (would FAIL) |

The evidence rows now read "DOS PASS, Kubo FAIL by O(1)", as the brief expects. The
dropped `S^-1 ... S` on the bra costs about 18 to 20 % of the largest moment and 19 to 20 %
of the Hall plateau. The first run's naive route had no metric at all and diverged
(smell 5); it is fixed.

### 6.6 Hall plateau (stage-3 Bastin, 10 meV)

Golden G5: C = +1.017017. All 10 lowdin and generalised routes, real and complex:
C = +1.017017, plateau 107.656458 over 450 points. **10 of 10 PASS.**

### 6.7 Thread-count repeat (chain1d, T = 1 against T = 16): the finding of this closure

**Fat node (TCN, g++ 11.4):** 20 of 24 IDENTICAL. The 4 DIFFERS rows are exactly the
`lowdin/exact` runs:

| run | first differing value (T=1 vs T=16) | max abs diff | golden_diff class |
|---|---|---:|---|
| chain1d_real_dos_lowdin_exact | line 4, m=1 re: 3.54669635674321e-17 vs 3.37457926186113e-17 | 2.97e-16 (m = 62) | rounding (DEFECT label is relative on tiny values) |
| chain1d_complex_dos_lowdin_exact | line 4, m=1 re: -2.82529533622566e-16 vs -2.83410447887711e-16 | 8.88e-16 (m = 62) | rounding |
| chain1d_real_kubo_lowdin_exact | line 4, (0,1) re: -3.87751354462285e-16 vs -3.86789125034204e-16 | 1.07e-14 at (11,11) | rounding |
| chain1d_complex_kubo_lowdin_exact | line 4, (0,1) re: -6.50934655569563e-16 vs -6.2560498231487e-16 | 9.99e-15 at (54,56) | rounding |

**Aron (tcn-aronnode, g++ 9.4, Eigen from `~/lsq-env`; `85046a7` built in
`~/lsquant-codex-6`, target `stage6_nonorth` only).** The legacy `kpm_lib` does not
compile with that Eigen; it is not needed here. Result: the same 4 rows DIFFER and the
other 20 are IDENTICAL. Max abs diff: 3.1e-16 and 6.9e-16 (dos), 1.02e-14 (both kubo
runs). Log: `figures/out/stage_6/aron_T1_vs_T16_chain1d.log`. The comparison is
machine-local; Aron's files are not compared with the fat node's.

**What these numbers rule out.** The differences are at the level of floating-point
reordering: at most 1.1e-14 on moments of order 1 to 4. They are far inside the 1e-10
tolerance against G1 and G2, and every lowdin/exact row passes against its golden at
both thread counts.

They are still a thread-count dependence, and only the exact Löwdin applier shows it:
- lowdin/chebyshev uses the same `LowdinOperator`, the same `RescaledOperator` and the same
  engine, and is byte-identical.
- generalised/exact uses the same dense-factorisation idea (`ExactInverseOverlap`, `LLT`
  plus `solve(Identity)`) and is byte-identical.

So the dependence is introduced in the `ExactInverseSqrtOverlap` constructor:
`overlap_appliers.cpp:206-217`, `SelfAdjointEigenSolver` followed by the dense product
`V * diag(lambda^-1/2) * V^+`.

**Likely cause (not proved):** Eigen parallelises dense matrix-matrix products with
OpenMP, and the blocking depends on the thread count. The LLT path uses triangular solves
and stays thread-invariant. A rebuild with `EIGEN_DONT_PARALLELIZE`, or computing that
product single-threaded, would decide it. I did not change source.

The engine itself is thread-invariant on every other row: generalised and Löwdin, dos and
kubo, chebyshev and ischo.

### 6.8 Figure

`figures/out/stage_6/fig_stage_6.png`, produced by `figures/stage_6.py` from closure
run 2 (exit 0). Panels and caption are as in section 3; the caption now carries
naive C = 1.197943 (real) and 1.193531 (complex), and naive Kubo 0.17 to 0.20.

### 6.9 Smells after review 2

Gone:
- **Smell 1**, the hidden 1e-18 drop threshold. It is now a named constant
  `kExactInverseDropThreshold` in the header.
- **Smell 2**, applier parameters never set by the program. `coeffTol` and `maxDimension`
  are now command-line flags, printed on stdout, and set explicitly by the harness. The
  header defaults remain.
- **Smell 5**, the naive route. It is now the brief's naive.
- **Smell 6**, the dead `--drop-tol` option. It is removed.
- **Smell 9**, the stage-1 grep hitting comment lines. The grep is now empty.
- **Smell 10**, the `pkill` loop. Nobody runs it any more. Closure run 1 was still killed
  by another session's `pkill` (6.2).

Remaining:
- **Smell 2 (partly):** default arguments on tolerances and caps remain in the appliers'
  constructors in the header (`maxDimension = 4096`, `tolerance = 1e-13`,
  `maxIterations = 200`, `coeffTol = 1e-12`). They are visible, but the brief's interfaces
  have none.
- **Smell 3:** the fifth applier `ExactInverseSqrtOverlap`; `lowdin/ischo` is undefined.
- **Smell 4:** the ISCHO diagnostics are last-call state, written from a `const Multiply`.
  The developer documented thread safety in review 2.
- **Smell 7:** the window padding constants are in the app.
- **New, smell 11:** `tests/dev/stage_6_check.py` defaults to `build/` for both the
  executable and `sys.path`. It needs `--exe` when the build directory is not `build/`.
- **New, smell 12 (the finding of 6.7):** `ExactInverseSqrtOverlap` is not thread-count
  invariant at the last bits. Byte identity across thread counts does not hold for
  lowdin/exact, while it holds for every other route and applier.
- **Still clean:**
  - No applier prints.
  - The engine, `Plan` and `io/` are not touched by the stage-6 commits.
  - `core/chebyshev` has only `GeneralisedRescaledOperator` and the brief-mandated
    `RescaledOperator` constructor.

### 6.10 Files (closure)

- Harness and figure scripts: `scripts/check_goldens.sh` (stages 4 and 5 plus 6, `6only`)
  and `figures/stage_6.py` (unchanged).
- Main logs in `figures/out/stage_6/`: `build.log`, `stage_6_check.log`,
  `check_goldens_stage6_stages1to5.log`, `check_goldens_stage6only.log`,
  `aron_T1_vs_T16_chain1d.log`.
- Figure and tables in `figures/out/stage_6/`: `fig_stage_6.png`,
  `fig_stage_6.caption.txt`, `stage_6_figure.log`, `data/gate_table.tsv`, `data/chern.tsv`.
- First-run outputs are kept under `*_run1_a5108e3` names.
