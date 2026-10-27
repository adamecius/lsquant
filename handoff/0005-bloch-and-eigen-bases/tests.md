# Handoff 0005 — tester report (stage 5 gate: Bloch and eigen bases, G9; G1 to G8 unchanged)

Tester: Codex role (Opus subagent driven by Fable). Machine: TCN (`fatnode-tacn`), worktree
`~/lsquant-unified-0005` of the clone `~/lsquant-unified`, branch `handoff/0005-codex`, own
build directory `~/lsquant-unified-build5` (the other testers' branches and builds were not
touched). Nothing under `lsquant/`, `python/lsquant/`, `tests/dev/`, the legacy tree, a seed, a
tolerance or a stage-0 to stage-4 golden was edited. Edited or added by the tester:
`scripts/check_goldens.sh` (stage 5), `scripts/hr_exactdiag_gamma.py` (codex.md 2.2 wrapper),
`figures/stage_5.py`, `figures/stage_5_table.py`, `figures/out/stage_5/`,
`goldens/stage0/PROVENANCE.md` (G9 section), `goldens/stage0/bases/` (G9), this file.

| Run | developer head | tester head before the run | G9 | gate | logs (`figures/out/stage_5/`) |
|---|---|---|---|---|---|
| 1 | `7109c20` | `81ad680` (fast-forward of `7109c20`) | shifted 48 x 48 mesh, Gaussian pi a / M (first codex.md) | FAIL 32/137, then FAIL 37/149 with the convention rows counted | `check_goldens_stage5_run1_7109c20.log`, `check_goldens_stage5_run1b_7109c20.log`, `build_run1.log`, `stage_5_check_run1.log`, `fig_stage_5_run1.*` |
| 2 | `dd989b0` (review 2, items 1 to 6) | merge `001d41a` of `origin/handoff/0005-agi` (contains `origin/cpp-unified` `5cb02fe`, so G4c and the stage-4 harness of `209dfd6`) | Gamma-centred 48 x 48, Gaussian 40 meV (revised codex.md 2.2) | PASS, 155/155 | `check_goldens_stage5_run2_dd989b0.log`, `stage_5_check_run2_dd989b0.log`, `table_T{16,1}_run2_dd989b0.txt`, `figure_run_run2_dd989b0.log`, `fig_stage_5_run2_dd989b0.*`, `build.log`, `gen_g9.log` |
| **3 (final)** | **`642955b`** (review 3: strict CSR receipt, operator sources in `basis.json`, `<seed>_<basis>` in the stems; Python only, no rebuild needed: the build of run 2 is the build of the same C++) | merge `2e7e5e7` | same G9 as run 2 (unchanged) | **PASS, 155/155** | `check_goldens_stage5.log`, `stage_5_check.log`, `table_T16.txt`, `table_T1.txt`, `figure_run.log`, `fig_stage_5.*` |

**Result (final run 3, `642955b`): GATE PASS (stage 5), 155 of 155 rows, 0 failures, exit 0
(28 min 21 s).** Every gated row and every convention row of run 3 is identical in text to
run 2 (diff of the PASS/FAIL lines empty), and the figure caption is byte-identical. The
numbers below are those of run 3 (= run 2). The self-check adds one line (the strict receipt,
PASS) and is otherwise unchanged.

- Stages 1 to 4 (real-space path): **71 of 71 IDENTICAL**, with G4c. All greps are empty.
- Real vs Bloch moments, tolerance 1e-9: **2.4e-14 or better at W = 0, 1.9e-14 or better at
  W = 1** (disorder seed 7 through `--disorder-seed`).
- Bloch (exact trace) vs eigen curves, tolerance 1e-8: **4.6e-15 or better**.
- `EigenBasis.exact_dos/exact_sigma_xx` vs G9 (hr_exactdiag on the same Gamma-centred mesh with
  the same Gaussian), tolerance 1e-8: **4.7e-15 or better**.
- Convention rows (20): the real basis is now the wannier2sparse CSR. That holds whether
  RealBasis loads the CSR or builds it, at 6 x 4 x 1 and at 48 x 48 x 1, for both models, to
  1.8e-16 or exactly 0. `B VX_w2s B^+ = v_x` holds to 7.6e-16 for graphene and Haldane at
  48 x 48.
- Thread counts: every stage-5 file is byte-identical between OMP_NUM_THREADS = 1 and 16. That
  is 96 files: 32 real-space files counted, 64 Bloch and eigen files reported.
- Report-only rows: KPM vs G9 (kernel shape) is 7.6e-3 to 1.8e-1. Finite size (G9 vs the
  shifted 240 x 240 mesh) is 1.1e-2 to 4.4e-1.
- All five run-1 defects (D1 to D5) are fixed, as section 9 shows.

## 1. Build (fresh, closure run)

```
$ rm -rf ~/lsquant-unified-build5        # HEAD 001d41a
$ cmake -S . -B ~/lsquant-unified-build5 -DLSQUANT_BUILD_NEW_CORE=ON -DLSQUANT_BUILD_PYTHON=ON -Dpybind11_DIR="$(python3 -c 'import pybind11; print(pybind11.get_cmake_dir())')"
--   Found fftw3, version 3.3.8
cmake exit: 0
$ cmake --build ~/lsquant-unified-build5 -j 32
[100%] Built target _core
build exit: 0
```

FFTW reaches `lsquant_core` through `pkg_check_modules(FFTW3 REQUIRED IMPORTED_TARGET fftw3)`
(top-level `CMakeLists.txt`) and `target_link_libraries(lsquant_core PUBLIC ... PkgConfig::FFTW3)`;
no extra flag. `fftw3.h` is included only in `lsquant/core/algebra/bloch_transform.cpp`.
Warnings: 101, all in the legacy tree; **0 under `lsquant/`**.

## 2. Developer self-check, `tests/dev/stage_5_check.py`, OMP_NUM_THREADS = 16 (run 3, exit 0)

```
[check a] BlochTransform vs dense B (6,4,1) graphene: rel(B x) = 3.563e-15, rel(B^+ x) = 3.221e-15
[check a] BlochTransform unitarity (8,8,1) graphene: max|B B^+ x - x| = 5.648e-17, max|B^+ B x - x| = 5.722e-17
[check b] (graphene 6x4x1) max|H_real_basis - H_w2s_csr| = 0.000e+00
[check b] (graphene 6x4x1) right-size CSR with wrong receipt raises ValueError: PASS
[check b] (graphene 6x4x1) max_x ||H_k x - B H_csr B^+ x||_inf = 4.780e-16
[check c] (graphene 6x4x1) max_x ||v_x x - B VX_csr B^+ x||_inf = 2.483e-16
[check b] (haldane 6x4x1) max|H_real_basis - H_w2s_csr| = 0.000e+00
[check b] (haldane 6x4x1) max_x ||H_k x - B H_csr B^+ x||_inf = 4.100e-16
[check c] (haldane 6x4x1) max_x ||v_x x - B VX_csr B^+ x||_inf = 2.582e-16
[check c] (haldane 48x48x1) rel ||v_x x - B VX_w2s B^+ x||_inf = 6.990e-16
[check d] Anderson disorder W=1.0 (seed=12345): max|mu_real - mu_bloch| DOS(64) = 4.719e-16, VX-VX(64x64) = 1.556e-15
[check e] (16,16,1) real vs bloch: DOS max|d_mu|=7.439e-15, max_rel|d_curve|=2.943e-15; KuboXX max|d_mu|=4.352e-14, max|d_curve/N|=1.399e-14, max_rel|d_curve|=9.851e-15
[check e] (16,16,1) real vs eigen: DOS max|d_mu|=7.438e-15, max_rel|d_curve|=2.968e-15; KuboXX max|d_mu|=4.308e-14, max|d_curve/N|=1.343e-14, max_rel|d_curve|=9.457e-15
[check e] (16,16,1) bloch vs eigen: DOS max|d_mu|=2.109e-15, max_rel|d_curve|=9.730e-16; KuboXX max|d_mu|=7.994e-15, max|d_curve/N|=1.278e-15, max_rel|d_curve|=9.001e-16
[check e] Jackson KPM (M=128) vs Gaussian exact-diag: rel_L2(DOS)=1.150e-02, rel_L2(KuboXX)=4.129e-02
[check e] (graphene 12x12x1) EigenBasis.exact_dos/exact_sigma_xx vs hr_exactdiag (Gamma mesh): rel_dos=1.822e-15, rel_sigma_xx=2.758e-15
[check f] G1 chain1d DOS byte-identical to golden: PASS
[check f] G4 graphene Kubo byte-identical to golden: PASS
[check f] G7 chain1d_mag spin precession byte-identical to golden: PASS
stage_5_check: ALL CHECKS PASSED

real	0m35.836s
EXIT:0
```

The review-2 additions are present: the dense-`B` check on 6 x 4 x 1, the non-square
6 x 4 x 1 equality with the wannier2sparse CSR, Haldane 48 x 48 `B VX_w2s B^+`, and the
12 x 12 Gamma-centred hr_exactdiag check.

## 3. G9 (regenerated per revised codex.md 2.2; PROVENANCE section G9)

`goldens/stage0/bases/gen_g9.sh` (wannier2sparse `3f553db`). For graphene and Haldane from
`examples/gen_models.py` (files committed unchanged under `<model>/model/`, with the bundle-mode
velocities `<model>_vx_hr.dat`, `_vy_hr.dat`), it runs:

- **G9 (the gated oracle):** `hr_exactdiag.py dos` and `sigma --comp xx` with `--nk 48`, on
  the **Gamma-centred mesh k = n/48**. The run goes through `scripts/hr_exactdiag_gamma.py`,
  a thin wrapper that loads the tool's source and replaces exactly its four
  `(np.arange(nk) + .5) / nk` mesh expressions by `np.arange(nk) / nk`. It aborts if the
  expected occurrences are not found, then runs the tool's own `main`. Parameters:
  **Gaussian eta = 40 meV**, grid [-3.6, 3.6] eV with 2881 points (2.5 meV spacing).
- **Finite size (report only):** the same commands on the tool's own shifted 240 x 240 mesh.

**sigma_xx factor.** `hr_exactdiag sigma` divides by A_hr = |a1_x a2_y| = 1.299038 A^2. That is
half the honeycomb cell |a1 x a2| = 2.598076 A^2, so the stored `sigma_xx_e2h` of these models
is twice the e^2/h value. The gate multiplies `EigenBasis.exact_sigma_xx`, which uses the
stage-3 `kubo_greenwood` normalisation pi Tr[v delta v delta], by
**2 pi / (A_hr N_k) = 2.0993e-3** (N_k = 2304). It then compares against the stored values
unchanged. The figure shows both in e^2/h with the true cell.

A second generation into a scratch tree reproduced all 20 files byte for byte; the sha256 sums
are in PROVENANCE.

## 4. Gate: `OMP_NUM_THREADS=16 S5_KEEP=~/s5work/s5out bash scripts/check_goldens.sh 5 ~/lsquant-unified-build5`

Log: `figures/out/stage_5/check_goldens_stage5.log` (`== 155 files compared; failures: 0`,
`GATE PASS (stage 5)`).

### 4.1 Stages 1 to 4 (the real-space path did not move)

| Stage | Files | Result |
|---|---|---|
| 1 (G1, G2 at T = 1, 8, 64; batched G2; random-phase determinism) | 11 | 11 IDENTICAL |
| 2 (G1, G3 from Python, T = 1, 16) | 6 | 6 IDENTICAL |
| 3 (G4, G5, **G4c**, T = 16, 1) | 30 | 30 IDENTICAL |
| 4 (G6 moments and 8 curves, G7 moments, T = 16, 1) | 24 | 24 IDENTICAL |
| **total** | **71** | **71 IDENTICAL** |

The stage-1 to stage-4 harness is the one in `cpp-unified` (`209dfd6`, with G4c). Greps:
`basis|std::cout|getenv` (core), `delta|DOS|Kubo|basis|std::cout|getenv` (core, io, code
lines), `compat/`, `time|velocity|spin|projector|bessel|gsl` (whole words, code lines), and
stage 5 `fftw` in `lsquant/core/chebyshev` and `lsquant/io`: all empty.

### 4.2 Stage-5 inputs as run

`stage5_bases` uses the brief's CLI:

```
--model goldens/stage0/bases/<m>/model --seed <m> --seed 12345 --supercell 48 48 1
--word dos|kuboxx --moments 256 --broadening-meV 40 --csr-dir <wannier2sparse CSR dir>
--basis real|bloch --trace random --states 32       (W = 0)
--basis bloch|eigen --trace exact                   (W = 0)
--basis real|bloch --trace random --states 32 --disorder 1 --disorder-seed 7   (W = 1)
```

- **Model files:** the unchanged gen_models.py files; the run-1 format adapter is gone.
- **Window:** `--window eigen` (default), recorded in `basis.json` as [-3.6, 3.6] eV at
  W = 0 and [-4.2, 4.2] eV at W = 1.
- **Real-space CSR:** the real basis reads the wannier2sparse CSR
  (`wannier2sparse <m> 48 48 1 VX`, written by the gate into its temporary directory).

### 4.3 Convention rows (counted; all PASS; `rel = max|a - b| / max|b|` over 4 random vectors)

| Row | graphene | Haldane |
|---|---|---|
| FFTW `B` vs dense `exp(-i k.(R + tau))/sqrt(N_k)`, iR = n1 + N1 n2, 6 x 4 x 1 (`B`, `B^+`) | 2.9e-15, 2.5e-15 | 2.9e-15, 2.5e-15 |
| (INFO) the opposite sign | 1.64 | 1.64 |
| 6 x 4 x 1: RealBasis loaded from `--csr-dir` vs wannier2sparse CSR (H, VX) | 0, 0 | 0, 0 |
| 6 x 4 x 1: RealBasis built by lsquant vs wannier2sparse CSR | 0, 0 | 0, 1.4e-16 |
| 6 x 4 x 1: `WannierModel.hamiltonian_csr / velocity_csr` vs wannier2sparse CSR | 0, 0 | 0, 1.4e-16 |
| 48 x 48: `B^+ B x = x` | 3.9e-16 | 3.9e-16 |
| 48 x 48: the three CSR rows above | 0 (all) | 0 / 1.8e-16 (VX built) |
| 48 x 48: `B H_w2s B^+` vs `H_k` | 7.6e-16 | 7.7e-16 |
| 48 x 48: `B VX_w2s B^+` vs `v_x(k)` (brief section 3, tol 1e-10) | 7.6e-16 | 7.4e-16 |

**FFTW convention, as verified:**

- `B` (RealToBloch) applies `FFTW_FORWARD` (exponent -1) on the (N1, N2, N3) grid, cells
  numbered n1 fastest, `howmany = W`, stride W. It then multiplies by `exp(-i k.tau_alpha)`
  and 1/sqrt(N_k).
- `B^+` (BlochToReal) multiplies by the conjugate phase, then applies `FFTW_BACKWARD` (+1).
- So `B_(k,a),(R,b) = delta_ab exp(-i k.(R + tau_a)) / sqrt(N_k)`, which is unitary. This is
  the corrected `brief.md` section 3.

### 4.4 The codex.md section 2 table (T = 16; the T = 1 table is identical in every gated row)

| Model, word | real vs Bloch moments W = 0 (1e-9) | Bloch(exact) vs eigen curve (1e-8) | eigen exact vs G9 (1e-8) | real vs Bloch moments W = 1, seed 7 (1e-9) |
|---|---|---|---|---|
| graphene DOS | PASS 2.28e-14 | PASS 2.20e-15 | PASS 3.90e-15 | PASS 1.89e-14 |
| graphene sigma_xx | PASS 2.15e-15 | PASS 3.27e-15 | PASS 4.74e-15 | PASS 2.96e-16 |
| Haldane DOS | PASS 2.38e-14 | PASS 1.18e-15 | PASS 3.68e-15 | PASS 1.94e-14 |
| Haldane sigma_xx | PASS 1.74e-15 | PASS 4.59e-15 | PASS 4.12e-15 | PASS 2.77e-16 |

Report-only rows. KPM = Jackson, M = 256, 40 meV. Masks: graphene |E| < 2.5, Haldane
0.5 < |E| < 2.5.

| Model, word | KPM real / Bloch (random 32) vs G9 | KPM Bloch(exact) / eigen vs G9 | G9 vs shifted 240 x 240 (finite size) | independent exact moments vs eigen / Bloch(exact) | random 32 vs eigen curve (stochastic) |
|---|---|---|---|---|---|
| graphene DOS | 4.74e-2 | 2.76e-2 | 7.18e-2 | 2.27e-14 / 2.27e-14 | 2.95e-2 |
| graphene sigma_xx | 1.43e-1 | 1.33e-1 | 4.43e-1 | 1.62e-14 / 1.46e-14 | 4.49e-2 |
| Haldane DOS | 1.62e-2 | 7.64e-3 | 1.08e-2 | 2.37e-14 / 2.37e-14 | 1.37e-2 |
| Haldane sigma_xx | 1.82e-1 | 1.71e-1 | 1.12e-1 | 9.32e-15 / 9.32e-15 | 4.82e-2 |

The "independent exact moments" are the Chebyshev moments of the Gamma-mesh system computed in
numpy. They use hr_exactdiag's `read_hr`/`opk` H(k) and the wannier2sparse velocity, with no
lsquant model code. Real vs Bloch curves, W = 0 and 1: 1.7e-15 to 2.0e-14. Pointwise `prel`
of the moment rows: 9.6e-13 to 3.8e-10, all at near-zero moments.

Decisions for the table:

- The 1e-8 Bloch-vs-eigen row uses Bloch with exact trace (BlockLocalStates). The random
  Bloch row differs from eigen by the stochastic error, so it is reported only.
- Curves are recomputed at full precision from the `.chebmom` files at the run's broadening
  (`basis.json`). The written `.dat` files have 6 significant digits (the legacy text);
  `written .dat vs recomputed` is 5.5e-7 to 3.6e-6 nrel.
- The gated quantity is `nrel = max|a - b| / max|b|`.

### 4.5 Thread counts

All 96 stage-5 output files are byte-identical between OMP_NUM_THREADS = 1 and 16:
`.chebmom`, `.chebmom1D`/`.chebmom2D`, `.dat` and `.basis.json`, for 2 models x 2 words x 6
runs. The 32 real-space files are counted; the 64 Bloch and eigen files are reported. FFTW
runs single-threaded, and the OpenMP loops (phases, blocks) write each element from one
thread.

### 4.6 Wall times (run 2; process wall time, seconds, T = 16 / T = 1; node shared; run 3 equal or up to 15 % faster)

| Run | graphene DOS | graphene kuboxx | Haldane DOS | Haldane kuboxx |
|---|---|---|---|---|
| W = 0 real (random 32) | 2.0 / 2.1 | 18.0 / 19.4 | 2.0 / 2.6 | 18.3 / 20.4 |
| W = 0 Bloch (random 32) | 1.9 / 2.0 | 17.7 / 19.1 | 1.9 / 2.2 | 17.1 / 17.9 |
| W = 0 Bloch (exact, 2 states) | 1.8 / 1.7 | 2.3 / 2.4 | 1.7 / 1.7 | 3.1 / 2.7 |
| W = 0 eigen (exact, 2 states) | 1.7 / 1.6 | 2.8 / 2.4 | 1.6 / 1.5 | 2.6 / 2.4 |
| W = 1 real (random 32) | 2.1 / 2.4 | 19.1 / 18.7 | 2.2 / 2.9 | 17.8 / 21.0 |
| W = 1 Bloch (random 32) | 4.0 / 4.0 | 23.8 / 21.8 | 4.0 / 3.7 | 20.5 / 23.1 |

With disorder, each Bloch product adds two FFTs. W = 1 Bloch Kubo costs 1.1 to 1.25 times the
real-space time. On this loaded node 16 threads give no gain at N = 4608.

## 5. Figure: `figures/stage_5.py`

Command:
`MPLBACKEND=Agg PYTHONPATH=~/lsquant-unified-build5/lsquant_core:python python3 figures/stage_5.py --runs ~/s5work/s5out/T16 --table figures/out/stage_5/table_T16.txt`.
It exits 0 and writes `figures/out/stage_5/fig_stage_5.png` and `fig_stage_5.caption.txt`.

What the figure shows:

- The three KPM bases lie on one another.
- The EigenBasis exact curve sits on the G9 markers, with a residual of 1e-14.
- The Jackson curves differ from the Gaussian ones only in kernel shape.
- At W = 1 the real - Bloch residual is about 1e-14, and disorder cuts the clean sigma_xx to
  about a third.

Caption:

> FIG. 5. Graphene (wannier2sparse examples/gen_models.py, t = -1 eV, nearest-neighbour distance 1 A) on the 48 x 48 x 1 supercell (4608 orbitals; Bloch mesh k = n/48), evaluated by lsquant.apps.stage5_bases in three bases from one model file, M = 256 Chebyshev moments, Jackson kernel, broadening 40.0 meV. Top: DOS per cell; bottom: sigma_xx (Kubo-Greenwood) in e^2/h with the true cell area 2.5981 A^2. Left, W = 0, window [-3.6, 3.6] eV (--window eigen): real space (the wannier2sparse CSR, random phase, 32 states, seed 12345, solid), Bloch (block-diagonal H(k), the same 32 vectors mapped by B, dashed), eigen (exact trace in the eigenbasis of each H(k), dotted); markers: G9, tools/hr_exactdiag.py (wannier2sparse 3f553db) dos and sigma --comp xx on the same Gamma-centred 48 x 48 mesh (scripts/hr_exactdiag_gamma.py), Gaussian eta = 40 meV; thin black line: EigenBasis.exact_dos/exact_sigma_xx with the same Gaussian. Residual below: exact - G9 (units 1e-14; max relative 3.9e-15 DOS, 4.7e-15 sigma_xx; the gated row). KPM eigen - G9 (Jackson vs Gaussian kernel shape, report only, not drawn): max relative 4.14e-02 DOS, 1.33e-01 sigma_xx. Right, W = 1 eV Anderson disorder (uniform in [-0.5, 0.5] eV, --disorder-seed 7), window [-4.2, 4.2] eV: real space (CSR + diag V) and Bloch (H(k) + B^+ V B through FFTW); clean G9 markers for reference; residual real - Bloch (units 1e-14; max relative on the curve 2.0e-14 DOS, 2.1e-15 sigma_xx); eigen skipped (no per-k eigenbasis with disorder). Gate rows (figures/stage_5_table.py, T = 16; identical at T = 1), graphene: DOS: real vs Bloch moments 2.284e-14 (W = 0), 1.891e-14 (W = 1), tol 1e-9; Bloch (exact trace) vs eigen curve 2.199e-15, tol 1e-8; EigenBasis exact vs G9 3.902e-15, tol 1e-8; sigma_xx: real vs Bloch moments 2.146e-15 (W = 0), 2.964e-16 (W = 1), tol 1e-9; Bloch (exact trace) vs eigen curve 3.274e-15, tol 1e-8; EigenBasis exact vs G9 4.735e-15, tol 1e-8.

## 6. Defects of run 1 and their state after review 2

| Run 1 | Now |
|---|---|
| D1 `read_hr` dense layout only | fixed: the unchanged gen_models.py `_hr.dat` loads (gate runs on the committed files) |
| D2 `read_xyz` needs a comment line | fixed: the unchanged single-header `.xyz` loads |
| D3 cell index n1 slowest vs wannier2sparse n1 fastest | fixed: all 20 convention rows PASS; RealBasis reads the wannier2sparse CSR (`--csr-dir`, or `<model_dir>`, `<model_dir>/operators`) and builds the identical one otherwise |
| D4 CLI differs from the brief | fixed: `--model`, `--disorder`, `--disorder-seed` (independent of the state seed), `--window eigen` default, legacy `.chebmom1D/2D` and legacy curve names plus `.chebmom` and `<stem>.basis.json`; the old spellings (`--model-dir`, `--disorder-w`, `--label`, `--rng-seed`) remain as aliases |
| D5 no exact evaluation in EigenBasis | fixed: `exact_dos`, `exact_sigma_xx` (GaussianDelta `SpectralFunction`), the G9 row at 4.7e-15 |

The brief corrections (cell numbering n1 fastest, sign of `B`) are verified by section 4.3.

## 7. Smells (closure run)

Gone since run 1:

- Gershgorin on a Bloch run: `--window eigen` is the default for every basis, and
  `BlochBasis`/`EigenBasis.window()` use the block eigenvalues. `--window gershgorin` remains
  an explicit option. On `--basis bloch|eigen` it builds a RealBasis only to compute the
  bound.
- All the CLI items.
- The format adapter.
- The disorder-seed wrapper.
- `report.md`'s "no open questions": section 5 now lists the known limitations.
- Silent fallbacks in the wannier2sparse CSR loading (run 2), gone in run 3 (review 3 item 1).
  A CSR of the right size with a wrong receipt now raises; the self-check line is
  `[check b] (graphene 6x4x1) right-size CSR with wrong receipt raises ValueError: PASS`.
  `basis.json` records the source and sha256 of every operator. For example, run 3's graphene
  real kuboxx has `"HAM": {"source": "wannier2sparse .../graphene.HAM.CSR 2f9b8614..."}` and
  the same for `VX`.
- Output names without the basis (run 2), gone in run 3 (review 3 item 2). The legacy stem
  carries `<seed>_<basis>` again, e.g. `SpectralOp1graphene_bloch_KPM_M256_statedefault.chebmom1D`
  and `NonEqOpVX-VXgraphene_real_KPM_M256x256_statedefault.chebmom2D`. The harness needed no
  change: it globs per run directory.

Remaining:

- **A basis name in the engine:** none. `grep basis lsquant/core` is empty. BlochTransform,
  BlockDiagonalOperator, SandwichOperator, BlockLocalStates and TransformedStates live in
  `core/algebra`.
- **A block format file:** none. Blocks travel as one `(N_k, W, W)` array.
- **FFTW included from `core/chebyshev`:** no.
- **The engine changed:** `lsquant/core/chebyshev/rescaled_operator.{hpp,cpp}` keeps the
  constructor from a general `LinearOperator` (on-the-fly `scale H + shift`), against "the
  engine does not change" (brief section 2). `MomentEngine`, `Plan` and `ChebyshevRecursion`
  are untouched, and G1 to G8 plus G4c are byte-identical.
- **Thread safety:** `SandwichOperator` (`mutable t1_, t2_`) and `BlochTransform` (one FFTW
  workspace per instance) are not safe for concurrent `Multiply` on one object. This is now
  documented in the headers and the report, and is fine today: the engine parallelises inside
  one `Multiply`.
- **Output precision:** curves are written with 6 significant digits (legacy), so 1e-8
  comparisons go through the moments.
- **Oracle side (wannier2sparse, not LSQUANT):**
  - `hr_exactdiag.py sigma` uses half the honeycomb cell area (factor stated in section 3).
  - Its mesh is hard-coded (hence the wrapper).
  - The pinned wannier2sparse `26cab4a` has no `hr_exactdiag.py`; G9 uses `3f553db`.

## 8. Run 1 (`7109c20`), summary

Run 1 used the first codex.md: G9 on hr_exactdiag's shifted 48 x 48 mesh with the Gaussian
pi a / M = 40.497 meV, and tolerances 2e-2 / 5e-2.

- Stages 1 to 4 were 65 of 65 IDENTICAL (before G4c).
- All three-basis rows passed at 1e-14.
- The 16 G9 rows failed at each thread count (2.9e-2 to 0.81). The mismatch was localised to
  the shifted vs Gamma mesh plus Jackson vs Gaussian; independent exact moments agreed with
  lsquant to 1e-14.
- The convention rows found D3: 5 failures, cell index n1 N2 + n2 against wannier2sparse's
  n1 + N1 n2.
- D1 and D2 forced a format adapter, and D4 a seed wrapper.
- Logs: `check_goldens_stage5_run1_7109c20.log` (137 rows, 32 failures) and
  `check_goldens_stage5_run1b_7109c20.log` (149 rows, 37 failures, with the convention rows),
  `build_run1.log`, `stage_5_check_run1.log`, `fig_stage_5_run1.png`, `.caption.txt`.
- The run-1 G9 files are superseded by the regenerated set in section 3 and are no longer in
  the tree.
