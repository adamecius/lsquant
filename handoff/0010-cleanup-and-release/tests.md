# Handoff 0010 (tester): final tests, `tests.md`

**same as before: YES.** No numeric output differs from v5 (`bbf81cf`), and the only failing text
check is fixed by `904c145` (details below).

Candidate: `handoff/0010-codex` at `f0792f8`, plus the harness fix `904c145`. Reference: v5 = `bbf81cf`.
`KPM_SEED=12345` everywhere. Sources: on TCN (`jose-fatnode`), the logs in
`~/lsquant-tester-0010-logs/`; on `jose-sov-01`, the results in `~/lsq-sov01/results/`. I re-read
both read-only for this report.

## 1. Harness changes

| Commit | Change |
|---|---|
| on `handoff/0010-codex` up to `f0792f8` | the four approved changes in `codex.md` section 1 (AGENTS.md/CLAUDE.md record exclusion, `stage_0.py` from `legacy-final`, `build/golden` fallback, G11 model dirs at `examples/0{7,8}_*/siesta`) |
| `904c145` | `stage7_deleted`: the `inline_*` text check also excludes `figures/out/final/legacy_vs_new_on_aron.txt`. That file is the recorded legacy-vs-new comparison, and its line 11 names `inline_compute-kpm-nonEqOp` as the provenance of G4c. The file is added to both `ex` and `bad`, and nothing else is excluded. |

## 2. Build and install (TCN)

| Build | Tree | Exit | Warnings (`grep -ci warning`) |
|---|---|---|---|
| `cmake --preset golden` | `f0792f8` | 0 | 0 |
| `cmake --preset native` (LSQUANT build profile: native) | `f0792f8` | 0 | 0 |
| v5, golden profile, separate build dir | `bbf81cf` | 0 | 0 |

`cmake --install` (golden) into `/tmp/lsq-tester-0010-install.AijM` completed: configuration
"Release", 91 install lines and no errors. The log has no output from the
`python -c "import lsquant, lsquant._core"` check, so that check is not recorded here.

## 3. Gates

### 3.1 Gate 7, golden build (64 run files; the 24 G11 run files are deferred to 3.4)

| Threads | Wall (s) | Exit | IDENTICAL lines | PASS lines | DIFFERS lines | FAIL lines | Summary line |
|---|---|---|---|---|---|---|---|
| 16 (`g7-T16-r2`) | 45551 | 1 | 722 | 55 | 0 | 1 | `808 files compared; failures: 1` |
| 1 (`g7-T1-r2`) | 77788 | 1 | 722 | 55 | 0 | 1 | `808 files compared; failures: 1` |

These are counts of lines that begin with each verdict. The per-golden totals table is the same at
16 and at 1 thread (only the `T16`/`T1` label differs), and every golden has 0 differing:

| G1 | G2 | G3 | G4 | G4c | G5 | G6 | G7 | Tut01 | Tut02 | Tut03 | Tut04 | Tut05 | Tut06 | G9 csr_dir | G9 | G9 table | G10 | G10 vs G1/G2/G5 | MCP |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | 1 | 3 | 6 | 2 | 6 | 5 | 6 | 24 | 194 | 9 | 4 | 17 | 281 | 4 | 96 | 16 | 44 | 22 | 7 |

The only FAIL, at both thread counts, is the text check
`no inline_* mention outside goldens/stage0/**/legacy_defects and PROVENANCE.md: 1 files: figures/out/final/legacy_vs_new_on_aron.txt`.
`904c145` fixes it. I verified the fix by re-running only that `git grep` pipeline on `904c145`
(`bad` is empty, and `ex` lists 36 files including `figures/out/final/legacy_vs_new_on_aron.txt`),
together with `bash -n scripts/check_goldens.sh`. I did not re-run the full gate. Every other text
check is PASS.

### 3.2 Gate 6only, golden build (`g6only`, 16 threads)

Exit 0, wall 48859 s, `98 files compared; failures: 0`, `GATE PASS (stage 6only)`. The log has
24 IDENTICAL lines, 74 PASS lines and 0 FAIL or DIFFERS lines. It also has 10 EVIDENCE lines, which
the gate does not count: these are the naive/exact rows, for example
`chain1d/real kubo naive/exact rel=1.989e-01` and `Hall haldane_real_kuboxy_naive_exact C=1.197943`,
as expected. There is no `bbf81cf` 6only log in the result set, so this verdict is compared with the
gate's own PASS criteria only.

### 3.3 Developer checks, HEAD (`f0792f8`) vs v5 (`bbf81cf`), 16 threads

| Check | HEAD rc / wall (s) | v5 rc / wall (s) | Logs (`cmp`) |
|---|---|---|---|
| `stage_6b_check.py` | 0 / 1859 | 0 / 1841 | identical (167 lines) |
| `stage_6c_check.py` | 0 / 11 | 0 / 12 | identical (35 lines) |
| `stage_6d_check.py --moments 8192 --backend folded` | 0 / 4858 | 0 / 4813 | identical (65 lines) |

### 3.4 G11 run files, HEAD vs v5, byte comparison with `cmp` over `out/`

| Set | Machine, threads | Pairs | Files | Identical | Differing files |
|---|---|---|---|---|---|
| Bi + PdSe2 DOS (3 + 3) | TCN, 16 | 6 | 54 | 42 | `basis.json`, `run.resolved.json` in each pair |
| Bi SHC: bare x3, berry_connection_anticommutator | TCN, 16 | 4 | 36 | 28 | `basis.json`, `run.resolved.json` in each pair |
| Bi SHC: berry_connection_{generalised,lowdin}, covariant x3 (`sov-bi-head` vs `sov-bi-v5`) | sov-01, 8 | 5 | 45 | 37 | `run.resolved.json` x5, `basis.json` x3 (covariant) |
| PdSe2 SHC reduced, 24x24, M=256, all 9 (`sov-red-head` vs `sov-red-v5`) | sov-01, 8 | 9 | 81 | 69 | `run.resolved.json` x9, `basis.json` x3 (covariant) |

Every moment file (`.chebmom`, `.chebmom1D`/`2D`), every curve (`dos.dat`, `shc.dat`) and every
oracle and comparison JSON (`exact_shc.json`, `spin_chern.json`, `eigen_*.json`) is byte-identical.
I compared the differing JSON files key by key:

* sov-01, all 14 pairs: the only differences are `output/dir`, `receipt/phases[*]/wall_seconds`,
  `dipole_cache/wall_seconds`, `dipole_quadrature/wall_seconds` and `dipole_cache/path`.
* TCN, all 10 pairs: the only differences are `output/dir`, `wall_seconds`, and the input paths
  `model/dir`, `fdf` and `input_file`. HEAD reads `examples/0N_*/siesta`; v5 reads the copy
  `~/lsquant-tester-bbf81cf-siesta/{bi,pdse2}`.

**Behaviour difference (not numeric): dipole cache location.** HEAD writes the dipole cache to
`examples/0N_*/.lsquant-cache/`, and v5 writes it to `examples/0N_*/siesta/.lsquant-cache/`. The
cache file name, and with it the content hash, is the same in both:
`dipoles-58b06a67…986a.npz` for Bi and `dipoles-6433b5be…f643.npz` for PdSe2. On sov-01 both v5
paths are under `~/lsq-sov01/head/examples/...`, so both runs read the same model inputs.

### 3.5 Native build

I did not run gate 7 with the `native` build: the result set has no native gate log. The native
build compiles with 0 warnings (section 2).

## 4. Examples (TCN, golden build, 16 threads, from a clean copy `~/lsquant-tester-0010-ex`)

| Example | Run files | rc | Wall (s) | plot.py rc |
|---|---|---|---|---|
| 01 DOS chain | `run_dos.json` | 0 | 654 | 0 |
| 02 localisation | `run_localization.json` | 0 | 565 | 0 |
| 03 Haldane Hall | `run_haldane.json` | 0 | 76 | 0 |
| 04 spin precession | `run_spin.json` | 0 | 339 | 0 |
| 05 charge transport | `run_transport.json` | 0 | 446 | 0 |
| 06 spin relaxation | `run_relaxation.json` | 0 | 13759 | 0 |
| 07 Bi quick | `run_dos_quick` / `run_shc_quick` | 0 / 0 | 3085 / 165 | 0 |
| 07 Bi production | `run_dos_production` / `run_shc_production` | 0 / 0 | 132 / 15114 | 0 |
| 08 PdSe2 quick | `run_dos_quick` / `run_shc_quick` | 0 / 0 | 2576 / 553 | 0 |
| 08 PdSe2 production | `run_dos_production` | 0 (run twice: 1580, 1344) | | not run |
| 08 PdSe2 production | `run_shc_production` | **not run to completion** | ~25-30 h per run | deferred |

Example 05 (from its log): at `W = 0.5, 1, 2, 4`, `sigma = 2.489, 0.770, 0.185, 0.044` and
`ell = 19.66, 6.16, 1.54, 0.35`. Not in this report: the overlay of examples 01-06 on the stage-0
captures with the max abs difference per curve, and `sigma^z_xy` at the README reference energies
for 07 and 08 (the result set does not contain them). The figure files are
`fig_dos_shc.png` in `bi_quick`, `bi_prod` and `pdse2_quick` under
`~/lsquant-tester-0010-ex-out/`.

## 5. Deviations from `codex.md`

1. **Full-size PdSe2 G11 SHC (48x48, M=1024) replaced, with the owner's approval,** by a reduced
   comparison (24x24, M=256). All 9 PdSe2 SHC variants were run with HEAD and with v5 on sov-01
   (section 3.4).
2. **5 Bi G11 SHC pairs rerun on sov-01** (HEAD and v5 on the same machine and environment,
   OMP 8) to get a same-environment comparison. The other 4 Bi SHC pairs and the 6 DOS pairs are
   from TCN.
3. Gate 7 was not re-run after `904c145`; only the changed pipeline was checked (section 3.1).
4. Gate 7 with the native build was not run, the install import check is not logged, and the
   example overlays and `sigma^z_xy` numbers are not reported (sections 2, 3.5, 4).
5. Example 08 production SHC was deferred (about 25-30 h per run).

## 6. Smells

* The default dipole cache location moved from `siesta/.lsquant-cache` to the example directory
  (section 3.4). The bytes are the same, but a cache that v5 left behind will not be found by HEAD.
* In `~/lsquant-tester-0010-logs/examples/summary.txt`, `08p` lists `run_dos_production` twice and
  has no `run_shc_production` line. `08p.log` shows that the SHC moments started twice and never
  completed.

Same as before: YES. Every numeric output (gate 7 at 16 and 1 thread, 6only, developer checks
6b/6c/6d, G11 DOS and SHC) is byte-identical to v5 or to the goldens. The only differences are in
paths and timings. The only FAIL was the `inline_*` text check, fixed by `904c145`.
