# Stage 5 Report — Bloch and Eigen Bases (`handoff/0005-agi`)

## 1. Files Created / Modified

`git diff --stat origin/handoff/0004-agi -- lsquant python tests` (stage-5 code delta on top of stage 4):

```text
 lsquant/CMakeLists.txt                           |    4 +-
 lsquant/api/bindings.cpp                         |  145 ++-
 lsquant/core/algebra/bloch_transform.cpp         |  217 +++++
 lsquant/core/algebra/bloch_transform.hpp         |  228 +++++
 lsquant/core/algebra/block_diagonal_operator.cpp |  121 +++
 lsquant/core/algebra/block_diagonal_operator.hpp |  126 +++
 lsquant/core/algebra/composite_operators.cpp     |   91 +-
 lsquant/core/algebra/composite_operators.hpp     |  176 +++-
 lsquant/core/algebra/state_generator.cpp         |   70 +-
 lsquant/core/algebra/state_generator.hpp         |  147 +++
 lsquant/core/chebyshev/rescaled_operator.cpp     |   62 +-
 lsquant/core/chebyshev/rescaled_operator.hpp     |   42 +-
 python/lsquant/__init__.py                       |   10 +-
 python/lsquant/_core.pyi                         |   95 ++
 python/lsquant/apps/stage5_bases.py              |  402 ++++++++
 python/lsquant/bases.py                          | 1092 ++++++++++++++++++++++
 python/lsquant/models/__init__.py                |   33 +
 python/lsquant/models/wannier.py                 |  633 +++++++++++++
 python/lsquant/moments.py                        |    2 +-
 python/lsquant/spectral_functions.py             |   19 +
 tests/dev/stage_5_check.py                       |  741 +++++++++++++++
 21 files changed, 4413 insertions(+), 43 deletions(-)
```

Full `git diff --stat main...handoff/0005-agi`:

```text
 .gitignore                                         |    4 +
 AGENTS.md                                          |  766 ++++++++++++++
 CMakeLists.txt                                     |    5 +
 agent_template/AGENTS.md                           |   62 ++
 agent_template/README.md                           |   44 +
 agent_template/install.sh                          |   80 ++
 agent_template/language_rules/c_cpp.md             |  153 +++
 agent_template/language_rules/python.md            |  176 ++++
 figures/stage_0.py                                 |  366 +++++++
 figures/stage_1.py                                 |  157 +++
 figures/stage_2.py                                 |  175 ++++
 figures/stage_3.py                                 |  226 +++++
 figures/stage_4.py                                 |  191 ++++
 handoff/0001-algebra-and-chebyshev-core/brief.md   |  152 +++
 handoff/0001-algebra-and-chebyshev-core/codex.md   |   89 ++
 handoff/0001-algebra-and-chebyshev-core/report.md  |  117 +++
 .../0002-pybind-sympy-rank1-contraction/brief.md   |  177 ++++
 .../0002-pybind-sympy-rank1-contraction/codex.md   |  116 +++
 .../0002-pybind-sympy-rank1-contraction/report.md  |  209 ++++
 handoff/0003-rank2-kubo-bastin-sea-surf/brief.md   |  215 ++++
 handoff/0003-rank2-kubo-bastin-sea-surf/codex.md   |  134 +++
 handoff/0003-rank2-kubo-bastin-sea-surf/report.md  |  246 +++++
 handoff/0004-time-evolution/brief.md               |  218 ++++
 handoff/0004-time-evolution/codex.md               |  111 ++
 handoff/0004-time-evolution/report.md              |  217 ++++
 handoff/0005-bloch-and-eigen-bases/brief.md        |  154 +++
 handoff/0005-bloch-and-eigen-bases/codex.md        |  103 ++
 handoff/0006-non-orthogonal-bases/brief.md         |  200 ++++
 handoff/0006-non-orthogonal-bases/codex.md         |  133 +++
 lsquant/CMakeLists.txt                             |   54 +
 lsquant/api/bindings.cpp                           |  731 ++++++++++++++
 lsquant/apps/stage1_moments.cpp                    |  302 ++++++
 lsquant/compat/legacy_kernels.cpp                  |  129 +++
 lsquant/compat/legacy_kernels.hpp                  |   67 ++
 lsquant/core/algebra/bloch_transform.cpp           |  202 ++++
 lsquant/core/algebra/bloch_transform.hpp           |  202 ++++
 lsquant/core/algebra/block_diagonal_operator.cpp   |  121 +++
 lsquant/core/algebra/block_diagonal_operator.hpp   |  132 +++
 lsquant/core/algebra/composite_operators.cpp       |  239 +++++
 lsquant/core/algebra/composite_operators.hpp       |  472 +++++++++
 lsquant/core/algebra/linalg.cpp                    |   98 ++
 lsquant/core/algebra/linalg.hpp                    |  120 +++
 lsquant/core/algebra/linear_operator.hpp           |  125 +++
 lsquant/core/algebra/sparse_matrix.cpp             |  466 +++++++++
 lsquant/core/algebra/sparse_matrix.hpp             |  304 ++++++
 lsquant/core/algebra/spectral_bounds.cpp           |  124 +++
 lsquant/core/algebra/spectral_bounds.hpp           |   84 ++
 lsquant/core/algebra/state_generator.cpp           |  186 ++++
 lsquant/core/algebra/state_generator.hpp           |  399 ++++++++
 lsquant/core/algebra/types.hpp                     |   34 +
 lsquant/core/chebyshev/chebyshev_moments.cpp       |  148 +++
 lsquant/core/chebyshev/chebyshev_moments.hpp       |  171 ++++
 lsquant/core/chebyshev/chebyshev_recursion.cpp     |   90 ++
 lsquant/core/chebyshev/chebyshev_recursion.hpp     |  145 +++
 lsquant/core/chebyshev/chebyshev_vectors.cpp       |  186 ++++
 lsquant/core/chebyshev/chebyshev_vectors.hpp       |  279 ++++++
 lsquant/core/chebyshev/contraction.cpp             |  153 +++
 lsquant/core/chebyshev/contraction.hpp             |  144 +++
 lsquant/core/chebyshev/moment_engine.cpp           |  439 ++++++++
 lsquant/core/chebyshev/moment_engine.hpp           |  108 ++
 lsquant/core/chebyshev/plan.hpp                    |  226 +++++
 lsquant/core/chebyshev/rescaled_operator.cpp       |   78 ++
 lsquant/core/chebyshev/rescaled_operator.hpp       |  120 +++
 lsquant/core/chebyshev/spectral_window.hpp         |  105 ++
 lsquant/io/chebmom_file.cpp                        |  275 +++++
 lsquant/io/chebmom_file.hpp                        |   96 ++
 lsquant/io/csr_text.cpp                            |  134 +++
 lsquant/io/csr_text.hpp                            |   43 +
 lsquant/io/curve_text.cpp                          |   53 +
 lsquant/io/curve_text.hpp                          |   45 +
 python/lsquant/__init__.py                         |   97 ++
 python/lsquant/_core.pyi                           |  382 +++++++
 python/lsquant/apps/__init__.py                    |    1 +
 python/lsquant/apps/stage2_dos.py                  |  280 ++++++
 python/lsquant/apps/stage3_kubo.py                 |  355 +++++++
 python/lsquant/apps/stage4_time.py                 |  412 ++++++++
 python/lsquant/apps/stage5_bases.py                |  272 +++++
 python/lsquant/bases.py                            |  408 ++++++++
 python/lsquant/kernels.py                          |  141 +++
 python/lsquant/models/__init__.py                  |   36 +
 python/lsquant/models/wannier.py                   |  467 +++++++++
 python/lsquant/moments.py                          |  470 +++++++++
 python/lsquant/reconstruct.py                      |  644 ++++++++++++
 python/lsquant/spectral_functions.py               |  586 +++++++++++
 scripts/check_goldens.sh                           |  395 ++++++++
 scripts/fetch_wannier2sparse.sh                    |   40 +
 scripts/golden_diff.py                             |  135 +++
 scripts/regenerate_goldens.sh                      |  285 ++++++
 tests/dev/stage_1_check.cpp                        |  145 +++
 tests/dev/stage_2_check.py                         |  168 ++++
 tests/dev/stage_3_check.py                         |  337 ++++++
 tests/dev/stage_4_check.py                         |  308 ++++++
 tests/dev/stage_5_check.py                         |  517 ++++++++++
```

---

## 2. Verified FFTW Sign, Normalisation, and Atom-Gauge Velocity Conventions

1. **Real-space hopping and cell-numbering convention (`wannier2sparse`)**:
   In `wannier2sparse` (`src/hopping_list.cpp`), a primitive hopping entry $(dR, H_{\alpha\beta}(dR))$ in `<seed>_hr.dat` connects source cell $R_1$ orbital $\alpha$ (row index $(i_{R_1} W) + \alpha$) to target cell $R_2 = (R_1 + dR) \bmod N$ orbital $\beta$ (column index $(i_{R_2} W) + \beta$):
   $$(H_{\text{CSR}})_{(R_1, \alpha), (R_1 + dR, \beta)} = H_{\alpha\beta}(dR),$$
   with $n_1$-fastest cell index $i_R = n_1 + N_1 n_2 + N_1 N_2 n_3$ (`index_aliasing` in `wannier2sparse`).

2. **Atom-gauge Bloch Hamiltonian (`brief.md` Section 3)**:
   $$H(k)_{\alpha\beta} = \sum_{dR} \exp\!\big(+i \mathbf{k} \cdot (\mathbf{dR} + \boldsymbol{\tau}_\beta - \boldsymbol{\tau}_\alpha)\big) H_{\alpha\beta}(dR).$$
   For $(B H_{\text{CSR}} B^\dagger)_{(k,\alpha),(k,\beta)} = H(k)_{\alpha\beta}$ to hold for all hoppings $dR$ and sublattice offsets $\boldsymbol{\tau}_\alpha$:
   $$\sum_{R_1, dR} B_{(k,\alpha),(R_1,\alpha)} H_{\alpha\beta}(dR) (B^\dagger)_{(R_1+dR,\beta),(k,\beta)} = \sum_{dR} e^{+i \mathbf{k} \cdot (\mathbf{dR} + \boldsymbol{\tau}_\beta - \boldsymbol{\tau}_\alpha)} H_{\alpha\beta}(dR)$$
   uniquely fixes:
   - **`RealToBloch` ($B$: real $\to$ Bloch)**:
     $$B_{(k,\alpha),(R,\beta)} = \delta_{\alpha\beta} \frac{1}{\sqrt{N_k}} \exp\!\big(-i \mathbf{k} \cdot (\mathbf{R} + \boldsymbol{\tau}_\alpha)\big),$$
     evaluated via `fftw_plan_many_dft` on `{N_3, N_2, N_1}` ($n_1$ fastest) with `howmany = W`, `istride = ostride = W`, `idist = odist = 1`, and `FFTW_FORWARD` ($\sum_R e^{-i \mathbf{k} \cdot \mathbf{R}}$), followed by pointwise multiplication by the precomputed phase table $\exp(-i \mathbf{k} \cdot \boldsymbol{\tau}_\alpha) / \sqrt{N_k}$.
   - **`BlochToReal` ($B^\dagger$: Bloch $\to$ real)**:
     $$(B^\dagger)_{(R,\beta),(k,\alpha)} = \delta_{\alpha\beta} \frac{1}{\sqrt{N_k}} \exp\!\big(+i \mathbf{k} \cdot (\mathbf{R} + \boldsymbol{\tau}_\alpha)\big),$$
     evaluated by pointwise multiplication by the conjugate phase $\exp(+i \mathbf{k} \cdot \boldsymbol{\tau}_\alpha) / \sqrt{N_k}$ followed by `FFTW_BACKWARD` ($\sum_k e^{+i \mathbf{k} \cdot \mathbf{R}}$).
   - This matches `lsquant-official/src/Hamiltonian/Graphene_KaneMele.cpp` (`to_kSpace_pruned` using `FFTW_FORWARD` and `to_rSpace_pruned` using `FFTW_BACKWARD`).

3. **Atom-gauge velocity sign convention (`wannier2sparse/include/tbmodel.hpp:236`)**:
   In `wannier2sparse`, `berry_connection` multiplies each hopping $H_{\alpha\beta}(dR)$ by $(-i)(\mathbf{dR} + \boldsymbol{\tau}_\beta - \boldsymbol{\tau}_\alpha)_a$. Consequently, its Bloch transform is
   $$v_a(k)_{\alpha\beta} = (B \, \text{V}_{a,\text{CSR}} \, B^\dagger)_{(k,\alpha),(k,\beta)} = -i \sum_{dR} (\mathbf{dR} + \boldsymbol{\tau}_\beta - \boldsymbol{\tau}_\alpha)_a \exp\!\big(+i \mathbf{k} \cdot (\mathbf{dR} + \boldsymbol{\tau}_\beta - \boldsymbol{\tau}_\alpha)\big) H_{\alpha\beta}(dR),$$
   giving $\max_x \| v_x x - B \, \text{VX}_{\text{CSR}} \, B^\dagger x \|_\infty \approx 7 \times 10^{-16}$ on both graphene and Haldane.

---

## 3. Verbatim Outputs on TCN

### 3.1 `tests/dev/stage_5_check.py` (`OMP_NUM_THREADS=16`)

```text
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
EXIT=0
```

### 3.2 `python -m lsquant.apps.stage5_bases` CLI Suite (`OMP_NUM_THREADS=16`, default `--window eigen`)

```text
Opening File: /tmp/s5_cli_r3_5tx9od8x/graphene/graphene.xyz
Velocity mode: berry_connection
Creating the Current matrix J0
$ python -m lsquant.apps.stage5_bases --model <graphene> --seed graphene --basis real --supercell 16 16 1 --word dos --moments 128 --trace exact --broadening-meV 200 --csr-dir <w2s_csr> --out <out/w0_dos>
/tmp/s5_cli_r3_5tx9od8x/w0_dos/SpectralOp1graphene_real_KPM_M128_stateexact.chebmom1D
/tmp/s5_cli_r3_5tx9od8x/w0_dos/SpectralOp1graphene_real_KPM_M128_stateexact.chebmom
/tmp/s5_cli_r3_5tx9od8x/w0_dos/meanSpectralOp1graphene_real_KPM_M128_stateexactJACKSON.dat
/tmp/s5_cli_r3_5tx9od8x/w0_dos/SpectralOp1graphene_real_KPM_M128_stateexact.basis.json
  operators in SpectralOp1graphene_real_KPM_M128_stateexact.basis.json: {"HAM": {"source": "wannier2sparse /tmp/s5_cli_r3_5tx9od8x/w2s_graphene_16x16/graphene.HAM.CSR beb91b9064cd91450a2ae4d0593f22213ee00ede582f53a20409b543acb9e692"}}
$ python -m lsquant.apps.stage5_bases --model <graphene> --seed graphene --basis bloch --supercell 16 16 1 --word dos --moments 128 --trace exact --broadening-meV 200 --out <out/w0_dos>
/tmp/s5_cli_r3_5tx9od8x/w0_dos/SpectralOp1graphene_bloch_KPM_M128_stateexact.chebmom1D
/tmp/s5_cli_r3_5tx9od8x/w0_dos/SpectralOp1graphene_bloch_KPM_M128_stateexact.chebmom
/tmp/s5_cli_r3_5tx9od8x/w0_dos/meanSpectralOp1graphene_bloch_KPM_M128_stateexactJACKSON.dat
/tmp/s5_cli_r3_5tx9od8x/w0_dos/SpectralOp1graphene_bloch_KPM_M128_stateexact.basis.json
  operators in SpectralOp1graphene_bloch_KPM_M128_stateexact.basis.json: {"HAM": {"source": "built by lsquant"}}
$ python -m lsquant.apps.stage5_bases --model <graphene> --seed graphene --basis eigen --supercell 16 16 1 --word dos --moments 128 --trace exact --broadening-meV 200 --out <out/w0_dos>
/tmp/s5_cli_r3_5tx9od8x/w0_dos/SpectralOp1graphene_eigen_KPM_M128_stateexact.chebmom1D
/tmp/s5_cli_r3_5tx9od8x/w0_dos/SpectralOp1graphene_eigen_KPM_M128_stateexact.chebmom
/tmp/s5_cli_r3_5tx9od8x/w0_dos/meanSpectralOp1graphene_eigen_KPM_M128_stateexactJACKSON.dat
/tmp/s5_cli_r3_5tx9od8x/w0_dos/SpectralOp1graphene_eigen_KPM_M128_stateexact.basis.json
  operators in SpectralOp1graphene_eigen_KPM_M128_stateexact.basis.json: {"HAM": {"source": "built by lsquant"}}
  [dos W=0 real vs bloch] max|d_mu|=7.439e-15, max|d_dat|=0.000e+00
  [dos W=0 real vs eigen] max|d_mu|=7.438e-15, max|d_dat|=0.000e+00
  [dos W=0 bloch vs eigen] max|d_mu|=2.109e-15, max|d_dat|=0.000e+00
$ python -m lsquant.apps.stage5_bases --model <graphene> --seed graphene --basis real --supercell 16 16 1 --word kuboxx --moments 128 --trace exact --broadening-meV 200 --csr-dir <w2s_csr> --out <out/w0_kuboxx>
/tmp/s5_cli_r3_5tx9od8x/w0_kuboxx/NonEqOpVX-VXgraphene_real_KPM_M128x128_stateexact.chebmom2D
/tmp/s5_cli_r3_5tx9od8x/w0_kuboxx/NonEqOpVX-VXgraphene_real_KPM_M128x128_stateexact.chebmom
/tmp/s5_cli_r3_5tx9od8x/w0_kuboxx/KuboGreenWood_NonEqOpVX-VXgraphene_real_KPM_M128x128_stateexactJACKSON.dat
/tmp/s5_cli_r3_5tx9od8x/w0_kuboxx/NonEqOpVX-VXgraphene_real_KPM_M128x128_stateexact.basis.json
  operators in NonEqOpVX-VXgraphene_real_KPM_M128x128_stateexact.basis.json: {"HAM": {"source": "wannier2sparse /tmp/s5_cli_r3_5tx9od8x/w2s_graphene_16x16/graphene.HAM.CSR beb91b9064cd91450a2ae4d0593f22213ee00ede582f53a20409b543acb9e692"}, "VX": {"source": "wannier2sparse /tmp/s5_cli_r3_5tx9od8x/w2s_graphene_16x16/graphene.VX.CSR 8eeffd788d5afe0a4c78473c7e0cb7da8d37fac5c85daf8c5b1ba5d2a2f13368"}}
$ python -m lsquant.apps.stage5_bases --model <graphene> --seed graphene --basis bloch --supercell 16 16 1 --word kuboxx --moments 128 --trace exact --broadening-meV 200 --out <out/w0_kuboxx>
/tmp/s5_cli_r3_5tx9od8x/w0_kuboxx/NonEqOpVX-VXgraphene_bloch_KPM_M128x128_stateexact.chebmom2D
/tmp/s5_cli_r3_5tx9od8x/w0_kuboxx/NonEqOpVX-VXgraphene_bloch_KPM_M128x128_stateexact.chebmom
/tmp/s5_cli_r3_5tx9od8x/w0_kuboxx/KuboGreenWood_NonEqOpVX-VXgraphene_bloch_KPM_M128x128_stateexactJACKSON.dat
/tmp/s5_cli_r3_5tx9od8x/w0_kuboxx/NonEqOpVX-VXgraphene_bloch_KPM_M128x128_stateexact.basis.json
  operators in NonEqOpVX-VXgraphene_bloch_KPM_M128x128_stateexact.basis.json: {"HAM": {"source": "built by lsquant"}, "VX": {"source": "built by lsquant"}}
$ python -m lsquant.apps.stage5_bases --model <graphene> --seed graphene --basis eigen --supercell 16 16 1 --word kuboxx --moments 128 --trace exact --broadening-meV 200 --out <out/w0_kuboxx>
/tmp/s5_cli_r3_5tx9od8x/w0_kuboxx/NonEqOpVX-VXgraphene_eigen_KPM_M128x128_stateexact.chebmom2D
/tmp/s5_cli_r3_5tx9od8x/w0_kuboxx/NonEqOpVX-VXgraphene_eigen_KPM_M128x128_stateexact.chebmom
/tmp/s5_cli_r3_5tx9od8x/w0_kuboxx/KuboGreenWood_NonEqOpVX-VXgraphene_eigen_KPM_M128x128_stateexactJACKSON.dat
/tmp/s5_cli_r3_5tx9od8x/w0_kuboxx/NonEqOpVX-VXgraphene_eigen_KPM_M128x128_stateexact.basis.json
  operators in NonEqOpVX-VXgraphene_eigen_KPM_M128x128_stateexact.basis.json: {"HAM": {"source": "built by lsquant"}, "VX": {"source": "built by lsquant"}}
  [kuboxx W=0 real vs bloch] max|d_mu|=4.352e-14, max|d_dat|=0.000e+00
  [kuboxx W=0 real vs eigen] max|d_mu|=4.308e-14, max|d_dat|=0.000e+00
  [kuboxx W=0 bloch vs eigen] max|d_mu|=7.994e-15, max|d_dat|=0.000e+00
$ python -m lsquant.apps.stage5_bases --model <graphene> --seed graphene --basis real --supercell 16 16 1 --word dos --moments 128 --trace random --states 8 --seed 12345 --disorder 1.0 --disorder-seed 7 --broadening-meV 200 --out <out/w1_dos>
/tmp/s5_cli_r3_5tx9od8x/w1_dos/SpectralOp1graphene_real_KPM_M128_statedefault.chebmom1D
/tmp/s5_cli_r3_5tx9od8x/w1_dos/SpectralOp1graphene_real_KPM_M128_statedefault.chebmom
/tmp/s5_cli_r3_5tx9od8x/w1_dos/meanSpectralOp1graphene_real_KPM_M128_statedefaultJACKSON.dat
/tmp/s5_cli_r3_5tx9od8x/w1_dos/SpectralOp1graphene_real_KPM_M128_statedefault.basis.json
  operators in SpectralOp1graphene_real_KPM_M128_statedefault.basis.json: {"HAM": {"source": "built by lsquant"}}
$ python -m lsquant.apps.stage5_bases --model <graphene> --seed graphene --basis bloch --supercell 16 16 1 --word dos --moments 128 --trace random --states 8 --seed 12345 --disorder 1.0 --disorder-seed 7 --broadening-meV 200 --out <out/w1_dos>
/tmp/s5_cli_r3_5tx9od8x/w1_dos/SpectralOp1graphene_bloch_KPM_M128_statedefault.chebmom1D
/tmp/s5_cli_r3_5tx9od8x/w1_dos/SpectralOp1graphene_bloch_KPM_M128_statedefault.chebmom
/tmp/s5_cli_r3_5tx9od8x/w1_dos/meanSpectralOp1graphene_bloch_KPM_M128_statedefaultJACKSON.dat
/tmp/s5_cli_r3_5tx9od8x/w1_dos/SpectralOp1graphene_bloch_KPM_M128_statedefault.basis.json
  operators in SpectralOp1graphene_bloch_KPM_M128_statedefault.basis.json: {"HAM": {"source": "built by lsquant"}}
  [dos W=1 (disorder-seed=7, state-seed=12345) real vs bloch] max|d_mu|=7.080e-16, max|d_dat|=0.000e+00
$ python -m lsquant.apps.stage5_bases --model <graphene> --seed graphene --basis real --supercell 16 16 1 --word kuboxx --moments 128 --trace random --states 8 --seed 12345 --disorder 1.0 --disorder-seed 7 --broadening-meV 200 --out <out/w1_kuboxx>
/tmp/s5_cli_r3_5tx9od8x/w1_kuboxx/NonEqOpVX-VXgraphene_real_KPM_M128x128_statedefault.chebmom2D
/tmp/s5_cli_r3_5tx9od8x/w1_kuboxx/NonEqOpVX-VXgraphene_real_KPM_M128x128_statedefault.chebmom
/tmp/s5_cli_r3_5tx9od8x/w1_kuboxx/KuboGreenWood_NonEqOpVX-VXgraphene_real_KPM_M128x128_statedefaultJACKSON.dat
/tmp/s5_cli_r3_5tx9od8x/w1_kuboxx/NonEqOpVX-VXgraphene_real_KPM_M128x128_statedefault.basis.json
  operators in NonEqOpVX-VXgraphene_real_KPM_M128x128_statedefault.basis.json: {"HAM": {"source": "built by lsquant"}, "VX": {"source": "built by lsquant"}}
$ python -m lsquant.apps.stage5_bases --model <graphene> --seed graphene --basis bloch --supercell 16 16 1 --word kuboxx --moments 128 --trace random --states 8 --seed 12345 --disorder 1.0 --disorder-seed 7 --broadening-meV 200 --out <out/w1_kuboxx>
/tmp/s5_cli_r3_5tx9od8x/w1_kuboxx/NonEqOpVX-VXgraphene_bloch_KPM_M128x128_statedefault.chebmom2D
/tmp/s5_cli_r3_5tx9od8x/w1_kuboxx/NonEqOpVX-VXgraphene_bloch_KPM_M128x128_statedefault.chebmom
/tmp/s5_cli_r3_5tx9od8x/w1_kuboxx/KuboGreenWood_NonEqOpVX-VXgraphene_bloch_KPM_M128x128_statedefaultJACKSON.dat
/tmp/s5_cli_r3_5tx9od8x/w1_kuboxx/NonEqOpVX-VXgraphene_bloch_KPM_M128x128_statedefault.basis.json
  operators in NonEqOpVX-VXgraphene_bloch_KPM_M128x128_statedefault.basis.json: {"HAM": {"source": "built by lsquant"}, "VX": {"source": "built by lsquant"}}
  [kuboxx W=1 (disorder-seed=7, state-seed=12345) real vs bloch] max|d_mu|=1.465e-15, max|d_dat|=0.000e+00
EXIT=0
```

---

## 4. Deviations from the Brief

1. **`TransformedStates(transform, inner)` in `lsquant/core/algebra/state_generator.hpp`**:
   Added `TransformedStates` alongside `BlockLocalStates` so that when `BlochBasis.states("random", num_states, seed)` (or `EigenBasis.states("random", ...)`) is requested, it applies $B$ (`RealToBloch`) to `RandomPhaseStates(dim, num_states, seed)`. Since $B$ is unitary, $|r_k\rangle = B |r_R\rangle$ is a properly normalised isotropic random state in Bloch space and represents the **exact same physical state vector** as $|r_R\rangle$ in real space, allowing finite-sample stochastic traces in `RealBasis` and `BlochBasis` (at both $W = 0$ and $W > 0$) to match to $\sim 10^{-15}$ for any number of states $R$.
2. **`--seed` argument handling in `python/lsquant/apps/stage5_bases.py`**:
   The CLI example in Section 6 of the brief passes both `--seed graphene` (model seed prefix) and `--seed 12345` (PRNG integer seed) on the same command line. `stage5_bases.py` accepts `--seed` with `action="append"` (plus optional `--label` and `--rng-seed`) and disambiguates integer state seeds from model prefix strings (unless a `<tok>_hr.dat` file literally exists for that token) so both single-`--seed` and dual-`--seed` invocations work seamlessly.

---

## 5. Known Limitations and Resolved Defect Notes

1. **Conflict with the first `brief.md` on the sign of $B$ and cell numbering (D3, resolved in `brief.md` Section 3 and Review 2 Item 1)**:
   - The first version of `brief.md` wrote $B_{(k,\alpha),(R,\beta)} = \delta_{\alpha\beta} N_k^{-1/2} \exp(+i \mathbf{k} \cdot (\mathbf{R} + \boldsymbol{\tau}_\alpha))$ (`FFTW_BACKWARD` for `RealToBloch`) and C row-major cell index $i_R = n_1 N_2 N_3 + n_2 N_3 + n_3$.
   - Because `wannier2sparse` builds $(H_{\text{CSR}})_{(R,\alpha),(R+dR,\beta)} = H_{\alpha\beta}(dR)$ with $n_1$-fastest cell index $i_R = n_1 + N_1 n_2 + N_1 N_2 n_3$, satisfying $B H_{\text{CSR}} B^\dagger = H(k)$ requires $B_{(k,\alpha),(R,\beta)} = \delta_{\alpha\beta} N_k^{-1/2} \exp(-i \mathbf{k} \cdot (\mathbf{R} + \boldsymbol{\tau}_\alpha))$ (`FFTW_FORWARD` on FFTW dimension array `{N_3, N_2, N_1}`) and $i_R = n_1 + N_1 n_2 + N_1 N_2 n_3$ (`n_1` fastest) everywhere. Both are now aligned with `wannier2sparse` and the corrected `brief.md` Section 3, and `RealBasis` loads `<seed>.HAM.CSR` and `<seed>.VX.CSR` written by `wannier2sparse` when present (or builds an identical CSR in Python otherwise).
2. **`read_hr` sparse vs. dense `_hr.dat` files (D1, resolved in Review 2 Item 2)**:
   - `read_hr` now parses `R_1 R_2 R_3 i j re im` hopping lines by content, assigning `ndegen[ir]` in the order of first appearance of each `R` (matching `wannier2sparse` and `hr_exactdiag.py`), so both sparse hand-written/`gen_models.py` files and dense Wannier90 files are supported. If a file lists fewer distinct `R` vectors than `nrpts` in the header (e.g. a vector whose entire $W \times W$ block is omitted), only the appearing `R` vectors are stored; duplicate `(R, i, j)` lines overwrite earlier ones in file order.
3. **`read_xyz` single-header vs. two-header `.xyz` files (D2, resolved in Review 2 Item 3)**:
   - `read_xyz` now accepts both `wannier2sparse`'s single-header format (orbital count line followed immediately by `<label> <x> <y> <z>` lines, no comment line) and the standard two-header XYZ format (orbital count line, comment line, orbital lines). Note that `wannier2sparse`'s own C++ reader (`read_xyz_file`) requires the single-header format without a comment line.
4. **Thread-safety of `BlochTransform` and `SandwichOperator` (Review 2 Item 6)**:
   - Each `BlochTransform` instance owns a single FFTW workspace (`inBuffer_`, `outBuffer_`, `plan_`), and each `SandwichOperator` instance owns mutable scratch vectors (`t1_`, `t2_`). Concurrent calls to `Multiply` on the *same* `BlochTransform` or `SandwichOperator` instance from multiple threads are therefore not thread-safe. This is safe in `MomentEngine`, which calls `Multiply` sequentially per state and parallelises with OpenMP only *inside* a single `Multiply` call.
5. **`EigenBasis` scope (`W = 0` only)**:
   - `EigenBasis` diagonalises the $N_k$ decoupled $W \times W$ blocks $H(k)$ of a translationally invariant crystal (`disorder_w == 0`). For non-zero real-space disorder (`disorder_w > 0`), $k$-blocks are coupled across the full $(N_k W) \times (N_k W)$ supercell and `EigenBasis` raises `ValueError`.
6. **Comparison with `wannier2sparse/tools/hr_exactdiag.py` (`EigenBasis.exact_dos` and `exact_sigma_xx`, Review 2 Item 4)**:
   - `hr_exactdiag.py` uses a half-step shifted mesh $k_i = (n_i + 1/2)/N_i$ by default (whereas the periodic $N_1 \times N_2 \times N_3$ supercell has the $\Gamma$-centred mesh $k_i = n_i / N_i$), and its `cmd_sigma` uses $A_{\text{hr}} = |a_{1x} a_{2y}|$ (half the honeycomb unit-cell area for `gen_models.py` graphene/Haldane) with prefactor $(2\pi^2 / A_{\text{hr}}) / N_k$. `EigenBasis.exact_sigma_xx` uses the stage-3 `kubo_greenwood` normalisation $\pi \sum_{k,n,m} |v_{x,nm}(k)|^2 g(E-\varepsilon_n) g(E-\varepsilon_m)$, so multiplying `exact_sigma_xx` by $2\pi / (A_{\text{hr}} N_k)$ matches `hr_exactdiag.py` evaluated on the same $\Gamma$-centred mesh to $\sim 2.8 \times 10^{-15}$.
