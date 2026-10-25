# Handoff 0006b — tester report (SIESTA non-orthogonal spin Hall, first execution of b51a8c6)

Tester: Opus subagent in the Codex role, driven by Fable. Machine: Aron (`tcn-aronnode`,
96 cores), toolchain `~/lsq-env` (python 3.10, numpy 2.2, scipy 1.15, sympy, sisl 0.16.4,
pybind11 3.1). Worktree `~/lsquant-codex-6b` at `origin/handoff/0006b-dev` = `b51a8c6`,
build directory `~/lsquant-codex-6b-build`, all logs under `~/t6b/` on Aron. Nothing under
`lsquant/`, `python/`, `tests/`, `runs/` or `goldens/` was edited; no golden was touched;
only processes started by this session were ever stopped (none were). Byte verdicts are not
issued (Aron toolchain, and stage 6b has no goldens yet). Diagnostics that needed a changed
behaviour were run from wrapper scripts outside the repository and are labelled as such.

## Summary

| Item | Verdict |
|---|---|
| Build | PASS (cmake rc 0, build rc 0, no warnings from the compiler) |
| Self-check `tests/dev/stage_6b_check.py` | **FAIL** at the dense-mesh step: tolerance defect in `checked_hermitian` at Gamma (section 2). With that one criterion patched in a wrapper, every assertion passes: covariance 8.9e-16, Hermiticity 0, HF 2.0e-10, plateau 2.000000000000, spin Chern 2, KPM 2.0055 vs oracle 1.9999 |
| SIESTA reader on Bi | PASS: sisl layout/gauge 2e-15, generalised eigenvalues vs sisl to 1e-6, eq.-4 velocity Hermitian to 3e-7 raw, HF 3e-5 |
| Dipole quadrature (Bi) | PASS, certified: refinement 9.87e-6 (< 1e-5 at order 96, the maximum), overlap vs stored S 9.2e-7, pair Hermiticity 8e-17, diag A(0) 2e-16; 256 s serial |
| Covariant velocity on Bi | Hermitian to 8e-13 raw, HF 2.7e-5; dipole term up to 20 eV A |
| Bi oracle | flat plateaus in the whole gap: **1.5527** (bare = berry_connection), **1.1747** (covariant), e^2/h Pauli; spin Chern 2, Z2 1 in every mode |
| Bi KPM (16 vectors, the run files) | +7.07 (berry_connection), -5.21 (covariant) at E=0: **noise-dominated**; production = eigen reference to 1e-11 |
| berry_connection shc run files | **FAIL** (rc 1) in the eigen-reference pass: `Operator 'JXSZ' ... is not Hermitian` (absolute 1e-12 test on C^+JC, section 5) |
| Random-phase trace of the SHC word in the 6b Bloch path | **variance defect**: std 11 per 16-vector sample on 60x60 Kane-Mele, against std 0.0041 for the stage-3 real-space path on the same model (section 6) |

## 1. Build

```
cmake -S . -B ~/lsquant-codex-6b-build -DCMAKE_BUILD_TYPE=Release -DLSQUANT_BUILD_NEW_CORE=ON -DLSQUANT_BUILD_PYTHON=ON \
  -Dpybind11_DIR=$(python -c 'import pybind11; print(pybind11.get_cmake_dir())') \
  -DCMAKE_DISABLE_FIND_PACKAGE_Eigen3=ON -DEIGEN3_INCLUDE_DIR=$HOME/lsq-env/include/eigen3 \
  -DFFTW3_DIR=$HOME/lsq-env/lib/cmake/fftw3 -DGSL_ROOT_DIR=$HOME/lsq-env      -> rc 0
cmake --build ~/lsquant-codex-6b-build -j 48                                     -> rc 0, 31.9 s
```

cmake notes only: the pybind11 rpath warning for `libgomp` (as in every Aron build) and
"FFTW3_DIR unused". `grep -iE "warning|error"` on the build log is empty. Products:
`_core.cpython-310-x86_64-linux-gnu.so`, `stage_1_check`, `stage1_moments`, `stage6_nonorth`.

## 2. Self-check (`LSQUANT_BUILD_DIR=... OMP_NUM_THREADS=16 KPM_SEED=12345`)

### 2.1 As committed: FAIL, rc 1, 2.1 s

Printed before the failure (symbolic part abbreviated; the full output is `~/t6b/selfcheck.log`):

- Metric forms and Hermiticity verdicts: V, J physical Hermitian; generalised S-Hermitian,
  not Hermitian; Löwdin Hermitian; anticommutator Hermitian. Sigma_z: all Hermitian.
- Symbolic trace proofs: Hamiltonian similarity True; DOS, spectral operator, Greenwood,
  Bastin, time correlation, projected evolution, MSD term: trace residual 0, equal True;
  anticommutator: single operator equal True, generic products equal False (pair residual
  `A*B*X**2/4 - A*X*B*X/2 + A*X**2*B/4`).
- Dipole term norms (axis y, 5 random k): 2.1355, 1.5501, 2.1024, 1.0633, 4.9778 (axis x: 0).
- `(a)-(d) covariance=8.88e-16 Hermiticity=0 HF=1.986e-10`

Then:

```
File "python/lsquant/bases.py", line 1281, in operator_blocks
File "python/lsquant/models/nonorthogonal.py", line 327, in bloch_velocity
File "python/lsquant/models/nonorthogonal.py", line 284, in checked_hermitian
ValueError: covariant velocity: anti-Hermitian relative norm 0.13285 exceeds 1e-10
```

**Localisation.** `checked_hermitian` (`python/lsquant/models/nonorthogonal.py`, lines
271-285) divides each block's anti-Hermitian norm by **that block's own norm** and exempts
only exact zeros. On the 96x96 mesh the covariant V_y at Gamma is physically zero: its
norm is 1.41e-15 and its anti-Hermitian part 1.9e-16 (round-off), relative 0.133. Every
other block is below 2.6e-13 relative; the absolute maximum over the mesh is 1.4e-15. The
formula is right; the criterion is not. A block-local relative test needs an absolute floor
or the largest block norm of the call as its scale. The same function guards SIESTA Bloch
sums, the overlap, the connection and the spin current, so any physically vanishing block
(high-symmetry k) can trip it.

### 2.2 Rerun with the criterion patched in a wrapper: PASS, rc 0, 11.8 s

Wrapper `/tmp/t6b_patched_check.py` (not in the repo) replaces `checked_hermitian` by the
same test with the largest block norm of the call as the scale, then runs the unchanged
script:

```
(a)-(d) covariance=8.8817841970012523e-16 Hermiticity=0 HF=1.9860313393849083e-10
(e) sigma_xy=2.000000000000 spin_chern=2.000000000000 Z2=1
(e) explicit Lowdin oracle plateau=2.000000000000
oracle eigenphase-minus-exact=-1.0004441719502211e-11
(f) generalised exact KPM M=256 Nk=12x12 plateau=2.005512511421 oracle=1.999923700485
SHC eigenphase Bloch-minus-eigen moments=1.8803272174855169e-14; random-minus-exact=0.10399082741868143
random phase plateaus generalised=-10.772790144273 Lowdin=-10.772790144273
anticommutator plateau=-4.557517873631; equality is NOT required by the derivation
DOS single eigenphase-minus-exact=1.0011542351459481e-15
stage_6b_check: developer self-check passed (not a golden gate)
```

Note the unasserted line: the 16-vector random-phase plateau is -10.77 against the exact
2.0055 on the same 12x12 mesh. Section 6 shows this is not bad luck.

## 3. SIESTA reader on Bi(111) (`~/from_fatnode/shc_bi/run`, `bi.fdf`, `bi.TSHS`)

Script `/tmp/t6b_siesta_check.py` (diagnostic, outside the repo), log `~/t6b/siesta_bi.log`.

- Model: 36 orbitals (W = 72 with spin), N_R = 25, lattice a = 4.34 A, printed E_F
  -3.0585 eV (not subtracted; stored H already relative to E_F). Load without dipoles 0.1 s.
- Built-in public-API checks (`basis.json`): H and S, gauges `R` and `r`, max abs
  difference 2.2e-16 to 2.0e-15; anti-Hermitian relative part of the Bloch sums 1.2e-10
  (H), 2.0e-10 (S), i.e. TSHS round-off. This pins the SOC component order
  `(Re Huu, Re Hdd, Re Hud, Im Hud, Im Huu, Im Hdd, Re Hdu, Im Hdu)`, the interleaved
  spinor order (orbital up, orbital down), the phase sign and the spin duplication of S.
- Generalised eigenvalues `eigh(H(k), S(k))` vs sisl's own `eigh(k)` at 6 k (Gamma, M,
  K, 3 random): 5.4e-8 to 1.1e-6 eV; min eig S 9.6e-4 to 4.3e-3 (cond up to ~1e3). Gap at
  Gamma [-0.186, +0.176] eV.
- Velocity (3 random k, x and y), Hellmann-Feynman by central differences of the
  generalised eigenvalues (delta 1e-4 1/A, degenerate multiplets diagonalised):

| velocity | max abs (V-V^+), before the fold | max abs V (eV A) | HF max error |
|---|---|---|---|
| eq. (4), berry_connection | 5.3e-8 to 2.9e-7 | 10.4-15.9 | 3.5e-6 to 2.9e-5 |
| bare (cell derivatives dressed in the atom gauge) | (built from folded inputs) | 12.9-16.3 | 3.5e-6 to 2.9e-5 |
| covariant (with the certified A) | 5.7e-14 to 8.3e-13 | 10.3-10.7 | 9e-9 to 2.7e-5 |
| legacy one-sided `dH - dS S^-1 H` | 19.3 to 33.0 | — | — |

  These reproduce notes section 11 (eq.-4 Hermitian to ~1e-7 from TSHS round-off, HF
  ~1e-5, legacy non-Hermitian by 20-30 eV A).
- Dipole quadrature (`load_siesta_model(..., dipoles=True)`, log `~/t6b/dip_bi.log`):
  256 s serial, certificate `largest_order 96, tolerance 1e-5, max_refinement_error
  9.870e-6, max_overlap_error 9.218e-7, pair_hermiticity_error 8.3e-17`; post-fold
  max|A(R)-A(-R)^+| = 0; max|diag A(0)| = 8e-17, 1.8e-16, 2.7e-17 (x, y, z); max|A| 0.776 A.
  Spot pairs: S converges to the stored overlap by order 48 (e.g. 0.053200, 0.032040,
  0.014646, -0.099019); same-centre s-s and identical-radial s-s pairs give A = 0 by
  parity; the s-p pair (i=1, j=20) gives A = (-0.023916, -0.031486, -0.018736) A, stable
  to 1e-6 from order 24. The refinement error sits just under the tolerance at the maximum
  order: the certificate passes with no margin.
- Size of the dipole term on Bi: up to 14.7-19.9 eV A per element; projected on the
  eigenstates within |E| < 3 eV it changes the velocity matrix by 6 % to 20 % (Frobenius).

## 4. Bi oracle and the plateau value

Oracle units and prefactor, as the code defines them (`oracle/spin_hall.py`):
`Omega_n = -2 Im sum_(m != n) J_nm V_mn / (E_n - E_m)^2`, `J = 1/2 {V_x, sigma_z}` with Pauli
sigma_z (eigenvalues +-1), `sigma = 2 pi / (N_k A_cell) sum_kn f(E_kn) Omega_kn`, e = hbar = 1,
reported as `e^2/h, sigma_z Pauli` (Kane-Mele gives 2.000, section 2.2). `reconstruct.shc`
uses the same `2 pi / (N_k A_cell)`.

Fine-grid oracle on 48x48 (`/tmp/t6b_bior.py`, Löwdin eigenvectors, eta = 0, step 0.02 eV):

| velocity | sigma^z_xy in the gap (-0.11 to +0.17 eV) | at E=0, eta 10 meV | outside: -0.13 / +0.19 eV |
|---|---|---|---|
| bare | 1.5527 (constant) | 1.55268 | 1.2849 / 1.5170 |
| berry_connection | 1.5527 (constant) | 1.55268 | 1.2849 / 1.5170 |
| covariant | 1.1747 (constant) | 1.17473 | 0.9718 / 1.1385 |

Spin Chern (FHS on P sigma_z P): C+ = +1.000, C- = -1.000, C_spin = 2.000, Z2 = 1, 30
occupied bands, energy gap 0.3032 eV. `projected_spin_gap` = 0.513 is the min |eigenvalue|
of P sigma_z P, **dimensionless (Pauli units), not eV**.

Spin mixing (developer's open question 5): eigenvalues s of `C_occ^+ Sigma_z C_occ` over the
mesh: min |s| 0.513, mean |s| 0.991, mean (1-|s|) 0.0092. At the gap edge the mixing is
strong: the highest occupied Kramers pair has |<sigma_z>| mean 0.614, min 0.099.

Reading: the plateau is flat (no grid artefact) and non-integer while the spin Chern number
is 2; the same prefactor gives 2.000 on Kane-Mele. The non-quantised value is the expected
consequence of strong spin mixing near the gap, not a prefactor. Against notes section 1:
the 6.09 there (and 0.97 = 6.09/2 pi vs the legacy 1.05) is not reproduced by this code for
any mode; the units or the operators of that day-20 exact calculation must be re-checked
before the legacy comparison is made. The dipole term lowers the Bi plateau by 24 %
(1.5527 -> 1.1747). The `exact_shc.json` the app writes uses 1001 points over the whole
97 eV window (0.097 eV step, 2 points in the gap): adequate for the value at E_F, useless for
plateau flatness; a gap-restricted grid would be better.

## 5. Bi run files (`python -m lsquant run`, 48 threads, one at a time)

All: 48x48x1, M = 1024, `trace random`, **16 eigenbasis phase vectors**, seed 12345,
10 meV, window from the generalised spectrum [-32.80, 64.49] eV. KPM plateau = `shc.dat` at
E = 0 and mean +- std over |E| < 0.1.

| run file | rc | wall (s) | oracle plateau (gap) | spin Chern | KPM at E=0 | KPM mean +- std, abs(E)<0.1 | production vs eigen reference |
|---|---|---|---|---|---|---|---|
| g11_bi_berry_connection_lowdin_dos | 0 | 187 | 1.5527 | 2 (Z2 1) | DOS | — | moments 5.6e-13, curve 8.4e-13 |
| g11_bi_berry_connection_lowdin_shc | **1** | 4273 | (1.5527, bior) | 2 | +7.070 | +7.083 +- 0.341 | not reached |
| g11_bi_berry_connection_generalised_shc | **1** | 4200 | (1.5527, bior) | 2 | +7.070 | +7.083 +- 0.341 | not reached |
| g11_bi_covariant_lowdin_shc | 0 | 8783 | 1.1747 | 2 (Z2 1) | -5.210 | -5.264 +- 0.100 | moments 3.9e-12, curve 5.9e-11 |
| g11_bi_covariant_generalised_shc | 0 | 10853 | 1.1747 | 2 (Z2 1) | -5.210 | -5.264 +- 0.100 | moments 4.2e-12, curve 7.3e-11 |

- Löwdin and generalised routes agree with each other to round-off, and the production
  Bloch path equals the independent eigenbasis engine run to 1e-11: the metric routes and
  the eigen binding are consistent. The KPM numbers themselves are noise (+7.07 vs 1.55,
  -5.21 vs 1.17): section 6.
- **berry_connection shc failure.** Both runs wrote the production moments and `shc.dat`,
  then stopped in the eigen-reference pass: `lsquant error: Operator 'JXSZ' in word
  'Tr[JXSZ delta(E-H) VY dG-(E)/dE]' is not Hermitian` (`python/lsquant/words/compiler.py:413`).
  `EigenMetricBasis.operator` (`python/lsquant/bases.py`) forms `C^+ J C` without a
  Hermitian fold, and `BlockDiagonalOperator` (`lsquant/core/algebra/block_diagonal_operator.cpp:37`)
  tests Hermiticity with an **absolute** 1e-12; on Bi the high-energy elements of C^+ J C
  reach hundreds of eV A, so round-off alone exceeds it. The covariant runs passed this
  test by chance. Because the oracle is written after the reference, these two runs have no
  `exact_shc.json`/`spin_chern.json`; their oracle values come from section 4.
- Cost: the covariant runs include the 256 s quadrature; the remainder is the two rank-2
  passes (production and reference) at M = 1024 x 1024, 16 vectors, dimension 165 888.
  PdSe2 was not run.

## 6. The random-phase trace of the spin Hall word

Brief section 5 makes the eigenbasis random phase (16 vectors) the trace of every gate
comparison. Measured spread of the 16-vector estimator of sigma^z_xy (value at E = 0,
one value per seed):

| model, mesh, M, broadening | trace | seeds | mean | s.e.m. | std per 16-vector sample |
|---|---|---|---|---|---|
| rotated KM, 12x12, 256, 10 meV (exact KPM 2.0055) | eigenbasis phase | 100 | 2.513 | 2.505 | 25.05 |
| orthogonal KM, 12x12, 256, 10 meV | eigenbasis phase | 100 | 3.166 | 2.608 | 26.08 |
| orthogonal KM, 12x12, 256, 10 meV | orbital-Bloch phase | 100 | -0.041 | 1.737 | 17.37 |
| orthogonal KM, 48x48, 256, 10 meV | eigenbasis phase | 20 | 3.822 | 1.491 | 6.67 |
| orthogonal KM, 48x48, 256, 10 meV | orbital-Bloch phase | 20 | 1.339 | 0.885 | 3.96 |
| orthogonal KM, 60x60, 512, 30 meV | eigenbasis phase | 10 | 1.726 | 3.531 | 11.17 |
| orthogonal KM, 60x60, 512, 30 meV | orbital-Bloch phase | 10 | 5.675 | 3.452 | 10.92 |
| orthogonal KM, 60x60, 512, 30 meV, **stage-3 real-space CSR path** | real-space phase | 10 | 2.0027 | 0.0013 | **0.0041** |

Earlier two-seed series (rotated KM, Löwdin, 12x12 / 48x48, M = 256): 16 vectors -10.8 / 63.5
and 0.07 / -2.26; 64 vectors 6.3 / 25.4 and 1.30 / -2.02; 256 vectors 11.5 / 7.4 and
1.48 / -0.75.

Stage-3 rows (`python -m lsquant.apps.stage3_kubo --operators ~/from_fatnode/shc_km/model/out
--label kanemele --left VXSZ --right VY --window gershgorin --trace random --states 16
--moments 512`, then `--formula bastin --broadening-meV 30`, built from this same worktree;
the wannier2sparse 60x60 Kane-Mele CSR of notes section 1, berry_connection mode, Hermitian):
seeds 1-10: 2.0014, 2.0047, 2.0041, 2.0040, 2.0036, 1.9963, 2.0019, 2.0061, 1.9951, 2.0098
(value at E = 0; the mean over |E| < 0.25 is the same to 4 digits); seed 12345: 1.9957, which
reproduces notes section 1 (1.9957). 113-141 s per seed at 32 threads.

"orbital-Bloch phase" = the same 6b `BlochMetricBasis` and compiled plan without
`bind_eigen_trace` (wrapper `/tmp/t6b_diag4.py`, outside the repo).

**Conclusion.** No estimator shows a significant bias (every 6b mean is within 1.3 s.e.m. of
the exact value). But on the same 60x60 Kane-Mele, with the same 16 vectors, M and
broadening, the stage-3 real-space path gives a spread of 0.0041 and reproduces the day-20
1.996, while both 6b Bloch-path estimators give a spread of 11. For a translation-invariant
(k-block-diagonal) operator M the variance of a random-phase sample is
`sum_(i != j) |M_ij|^2 = ||M||_F^2 - sum_i |M_ii|^2`; `||M||_F` is basis independent and the
real-space diagonal is the k-average of the k-space diagonal, so a phase vector in the
orbital-Bloch or eigen basis must have **at most** the real-space variance. A ratio of
~2700 in the standard deviation therefore cannot be a property of the trace basis: it is a
**defect of the 6b sampling path** (variance, not bias), common to both phase choices and
present with S = I, hence not in the metric. Each 6b sample is evidently not the quadratic
form `r^+ M r` of one random vector; the place to look is the rank-2 plan on block
operators: `RandomPhaseStates` over the block layout, the `MirrorToBra` bra seed in rank-2
Bastin plans, `bind_eigen_trace` (`python/lsquant/bases.py`), and `ApplyKpmWeights(mu,
nstates, False)` in `apps/stage6b_siesta.py`. Not localised to a line. Cheap discriminating
test for the developer: one sample, M = 8, 4x4 mesh, compare the engine's mu_nm with
`r^+ (T_n-side ... ) r` computed in numpy for the same r.

Consequence for brief section 5: as implemented, the eigenbasis random phase cannot be the
trace of the spin Hall gate (the self-check's own unasserted line, -10.77 vs 2.0055, and the
Bi KPM values +7.07 / -5.21 against 1.55 / 1.17 are this defect). Exactness for diagonal
words holds (DOS single eigenphase vs exact 1e-15; Bi DOS production vs reference 6e-13).
The section can stand only after the sampling defect is fixed and re-measured with these
tables; until then no Bi or PdSe2 KPM plateau can be judged.

## 7. Profile of a short Bi run (review-1 item 5, "before")

Copy of `g11_bi_berry_connection_lowdin_shc.json` with `moments 128`, `states 2`
(`~/t6b/prof/prof_run.json`), `OMP_NUM_THREADS=48`, `python -m cProfile -o prof.out -m lsquant
run`. There is no switch to disable the eigen reference, so the run ends in the reference
pass with the JXSZ Hermiticity error of section 5 (rc 1); the profile covers loading, basis
construction, operator construction and the production moments. Wall 86 s, profiled 81.2 s.

Top cumulative entries (abridged to the non-trivial ones of the top 20):

| cumtime (s) | tottime (s) | ncalls | function |
|---|---|---|---|
| 80.50 | 0.00 | 1 | `runfile.execute_runfile` -> `apps/stage6b_siesta.execute_siesta_runfile` |
| 53.83 | 16.75 | 2 | `apps/stage6b_siesta.compute_moments` (tottime = the C++ engine call, pybind methods are not profiled separately) |
| 50.67 | 0.07 | 6 | `models/siesta.SiestaModel._sum` |
| 48.92 | 0.02 | 10 | `models/bloch.bloch_sum` |
| 48.67 | **47.82** | 10 | `models/bloch.py:36 <listcomp>` (the per-k Python loop of the Bloch sum) |
| 37.08 | 0.00 | 2 | `compute_moments <dictcomp>` (operator construction) |
| 36.35 | 0.24 | 2 | `bases.BlochMetricBasis.operator` |
| 36.10 | 0.00 | 5/4 | `bases.BlochMetricBasis.operator_blocks` |
| 33.74 | 0.00 | 2 | `models/siesta.derivatives` |
| 25.60 | 2.10 | 1 | `bases.BlochMetricBasis.__init__` |
| 8.48 / 8.45 | 0.00 | 1 / 1 | `SiestaModel.Hk` / `SiestaModel.Sk` |
| 5.80 | 5.80 | 2 | `numpy.linalg.eigh` |
| 3.48 | 0.35 | 1 | `nonorthogonal.exact_block_metric` |
| 2.90 | 1.13 | 10 | `nonorthogonal.checked_hermitian` |
| 1.76 | 1.18 | 2 | `nonorthogonal.bloch_velocity` |

Fraction of time: the Python Bloch sums alone take 48.7 s of 81.2 s (60 %); the engine call
(C++ `_core`, inside the 16.75 s tottime of `compute_moments`) at most 21 %; eigh/metric
~7 %. A filename-based count of `_core` frames gives 0.2 % only because pybind calls are
attributed to their Python caller. The `bloch_sum` loop is a serial Python list comprehension
over 2304 k points with an (N_R, W, W, 3) phase product per k.

`ps` every 5 s (elapsed s, %CPU, threads, RSS kB); note that `ps` %CPU is the
lifetime average (CPU time / elapsed), not instantaneous:

```
 0     0   1    10672
 5   216  48   279492
10   158  48   385616
15   139  48   502244
21   385  48   865504
26  1005  48  2034944
31  1004  48  2328568
36   875  48  2565580
41   778  48  2552440
46   723  48  2705852
51   661  48  3035496
61   567  48  3151128
66   655  95  3903660
71   929  95  3903660
76   948  95  3903660
81   939  95  3429472
```

Reading: during the Bloch-sum and operator-construction phase (5-61 s) the lifetime average
falls from ~1000 % to 567 %, i.e. that phase runs on few cores (the loop is serial Python;
the numpy calls inside it are partly threaded); the rise after 61 s (thread count 95, two
OpenMP pools) is the engine, which then runs at roughly 40 cores of 48. For the full M = 1024,
16-vector runs (section 5), the 70-180 min wall times are dominated by the two rank-2 engine
passes; at small M the serial Python Bloch sums dominate.

## 8. Smells

- Wannier import in the SIESTA route: none in `models/siesta.py`; `runfile.py` imports
  `models.wannier` at module level (pre-existing, not used on the SIESTA path).
- Silent fallback: none found. Covariant without dipoles is refused (self-check and
  `BlochMetricBasis`); missing radial functions or an uncertified quadrature raise.
- Zero-dipole shortcut: none. The zero connection in the self-check is the true connection
  of the unrotated orthogonal model (Q = diag(tau), A = Q - tau S = 0).
- Golden touched: none (`git diff --stat` on `goldens/`, `figures/`, `scripts/`,
  `lsquant/core` vs ab3adeb is empty). The 24 `runs/goldens/g11_*.json` are new run files,
  not golden data.
- `eval`/`exec`: none; `sympy.sympify` only on internal expressions in `metric_forms.py`.
- Tolerance policy: SIESTA Bloch sums accept and fold relative round-off up to 1e-7 (declared
  in report.md, measured here at 1e-10 to 3.5e-10). Two Hermiticity tests have the wrong
  scale: block-relative with no floor (`nonorthogonal.py` 271-285, section 2) and absolute
  1e-12 on large blocks (`block_diagonal_operator.cpp:37` via `EigenMetricBasis`, section 5).
- The quadrature certificate passes only at the maximum order with the refinement error
  just under the tolerance (9.87e-6 vs 1e-5): no margin; 256 s serial per covariant load,
  repeated in every covariant run.
- Spin-Chern convention: `C_spin = C+ - C-` = 2 for KM, `Z2 = (C+ - C-)/2 mod 2` = 1. This
  deviates from the brief's `Z2 = C_s mod 2` (which would give 0 for KM), as the report says.

## 9. Files

On Aron (`~/t6b/`): `cmake.log`, `build.log`, `selfcheck.log`, `selfcheck_patched.log`,
`siesta_bi.log`, `dip_bi.log`, `dipoles_bi.npz`, `bior.log`, `diag3.log`, `diag4.log`,
`job2.log`, `km3/s*/`, `prof/`, `runs_bi.log`, `runs/*.log`, `runs/out/g11_bi_*`.
Diagnostic scripts (outside the repo): `/tmp/t6b_patched_check.py`, `/tmp/t6b_siesta_check.py`,
`/tmp/t6b_dip.py`, `/tmp/t6b_cov.py`, `/tmp/t6b_bior.py`, `/tmp/t6b_diag2.py`,
`/tmp/t6b_diag3.py`, `/tmp/t6b_diag4.py`, `/tmp/t6b_plat.py`, `~/t6b/job2.sh`.
