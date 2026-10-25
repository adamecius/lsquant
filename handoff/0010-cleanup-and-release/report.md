# Handoff 0010 — Clean-up, Build, README, and Examples Report

## 0. Starting Point and Merge Resolution

- **Branch:** `handoff/0010-dev` cut from `origin/handoff/0010-agi` (`77df68db4b03594cb250920c5c751bc416a468ca`).
- **Merge resolution:** None required. `origin/handoff/0010-agi` already contains the `handoff/0009-wrapup` integration commit (`bbf81cf`: stages 6b, 6c, 6d, and 8, tree `bf70c9b8`), the wrap-up reports, this handoff's brief, and the committed SIESTA inputs (`SHA256SUMS`, `.fdf`, `.TSHS`, `.ion.xml`, `.psml`) for `examples/07_spin_hall_bi_bilayer/siesta/` and `examples/08_spin_hall_pdse2_monolayer/siesta/`.
- **Verification host:** All builds, developer self-checks, and example smoke runs were executed on TCN (`jh-x`, `Ubuntu 22.04`, `g++ 11.4.0`, 64 logical cores) in `~/lsquant-agi-build-0010`. No compilation or execution was performed on Ionos.

---

## 1. Step 1 — Baseline Build and Self-Check on TCN

### Commands executed on TCN (`~/lsquant-agi-build-0010`, commit `77df68d`)

```bash
cmake -S . -B build-step1-golden -DCMAKE_BUILD_TYPE=Release -DLSQUANT_BUILD_PROFILE=golden \
  -DLSQUANT_BUILD_BENCHMARKS=ON -Dpybind11_DIR=$(python3 -c "import pybind11; print(pybind11.get_cmake_dir())")
cmake --build build-step1-golden -j 8

cmake -S . -B build-step1-native -DCMAKE_BUILD_TYPE=Release -DLSQUANT_BUILD_PROFILE=native \
  -DLSQUANT_BUILD_BENCHMARKS=ON -Dpybind11_DIR=$(python3 -c "import pybind11; print(pybind11.get_cmake_dir())")
cmake --build build-step1-native -j 8

./build-step1-golden/lsquant_core/stage_1_check
./build-step1-native/lsquant_core/stage_1_check
```

### Baseline build log (`LSQUANT_BUILD_PROFILE=golden`, 0 compiler warnings)

```text
-- The CXX compiler identification is GNU 11.4.0
-- Detecting CXX compiler ABI info
-- Detecting CXX compiler ABI info - done
-- Check for working CXX compiler: /usr/bin/c++ - skipped
-- Detecting CXX compile features
-- Detecting CXX compile features - done
-- LSQUANT build profile: golden
-- Found OpenMP_CXX: -fopenmp (found version "4.5") 
-- Found OpenMP: TRUE (found version "4.5")  
-- Found PkgConfig: /usr/bin/pkg-config (found version "0.29.2") 
-- Found GSL: /usr/include (found version "2.7.1") 
-- Looking for C++ include pthread.h
-- Looking for C++ include pthread.h - found
-- Performing Test CMAKE_HAVE_LIBC_PTHREAD
-- Performing Test CMAKE_HAVE_LIBC_PTHREAD - Success
-- Found Threads: TRUE  
-- Checking for module 'fftw3'
--   Found fftw3, version 3.3.8
-- Found PythonInterp: /usr/bin/python3 (found suitable version "3.10.12", minimum required is "3.9") 
-- Found PythonLibs: /usr/lib/x86_64-linux-gnu/libpython3.10.so
-- Performing Test HAS_FLTO_AUTO
-- Performing Test HAS_FLTO_AUTO - Success
-- Found pybind11: /home/ICN2/jgarcia/.local/lib/python3.10/site-packages/pybind11/include (found version "3.1.0")
-- Configuring done
-- Generating done
-- Build files have been written to: /home/ICN2/jgarcia/lsquant-agi-build-0010/build-step1-golden
[  3%] Building CXX object lsquant_core/CMakeFiles/lsquant_core.dir/core/algebra/sparse_matrix.cpp.o
[  6%] Building CXX object lsquant_core/CMakeFiles/lsquant_core.dir/core/algebra/bloch_transform.cpp.o
[  9%] Building CXX object lsquant_core/CMakeFiles/lsquant_core.dir/core/algebra/block_diagonal_operator.cpp.o
[ 12%] Building CXX object lsquant_core/CMakeFiles/lsquant_core.dir/core/algebra/linalg.cpp.o
[ 19%] Building CXX object lsquant_core/CMakeFiles/lsquant_core.dir/core/algebra/state_generator.cpp.o
[ 19%] Building CXX object lsquant_core/CMakeFiles/lsquant_core.dir/core/algebra/composite_operators.cpp.o
[ 22%] Building CXX object lsquant_core/CMakeFiles/lsquant_core.dir/core/algebra/spectral_bounds.cpp.o
[ 25%] Building CXX object lsquant_core/CMakeFiles/lsquant_core.dir/core/algebra/metric_policy.cpp.o
[ 29%] Building CXX object lsquant_core/CMakeFiles/lsquant_core.dir/core/algebra/overlap_appliers.cpp.o
[ 32%] Building CXX object lsquant_core/CMakeFiles/lsquant_core.dir/core/chebyshev/rescaled_operator.cpp.o
[ 41%] Building CXX object lsquant_core/CMakeFiles/lsquant_core.dir/core/chebyshev/chebyshev_recursion.cpp.o
[ 41%] Building CXX object lsquant_core/CMakeFiles/lsquant_core.dir/core/chebyshev/contraction.cpp.o
[ 41%] Building CXX object lsquant_core/CMakeFiles/lsquant_core.dir/core/chebyshev/chebyshev_vectors.cpp.o
[ 45%] Building CXX object lsquant_core/CMakeFiles/lsquant_core.dir/core/chebyshev/moment_engine.cpp.o
[ 51%] Building CXX object lsquant_core/CMakeFiles/lsquant_core.dir/core/algebra/fft_moment_contraction.cpp.o
[ 51%] Building CXX object lsquant_core/CMakeFiles/lsquant_core.dir/core/chebyshev/chebyshev_moments.cpp.o
[ 54%] Building CXX object lsquant_core/CMakeFiles/lsquant_core.dir/io/csr_text.cpp.o
[ 58%] Building CXX object lsquant_core/CMakeFiles/lsquant_core.dir/io/chebmom_file.cpp.o
[ 61%] Building CXX object lsquant_core/CMakeFiles/lsquant_core.dir/io/curve_text.cpp.o
[ 64%] Linking CXX static library ../lib/liblsquant_core.a
[ 64%] Built target lsquant_core
[ 67%] Building CXX object lsquant_core/CMakeFiles/stage1_moments.dir/apps/stage1_moments.cpp.o
[ 70%] Building CXX object lsquant_core/CMakeFiles/stage_1_check.dir/__/tests/dev/stage_1_check.cpp.o
[ 77%] Building CXX object lsquant_core/CMakeFiles/_core.dir/api/bindings.cpp.o
[ 77%] Building CXX object lsquant_core/CMakeFiles/stage6_nonorth.dir/apps/stage6_nonorth.cpp.o
[ 80%] Building CXX object lsquant_core/CMakeFiles/_core.dir/compat/legacy_kernels.cpp.o
[ 83%] Building CXX object lsquant_core/CMakeFiles/stage8_kernels.dir/__/bench/kernel_bench.cpp.o
[ 87%] Linking CXX executable stage8_kernels
[ 90%] Linking CXX executable stage1_moments
[ 90%] Built target stage8_kernels
[ 90%] Built target stage1_moments
[ 93%] Linking CXX executable stage6_nonorth
[ 93%] Built target stage6_nonorth
[ 96%] Linking CXX executable stage_1_check
[ 96%] Built target stage_1_check
[100%] Linking CXX shared module _core.cpython-310-x86_64-linux-gnu.so
[100%] Built target _core
```

### Baseline build log (`LSQUANT_BUILD_PROFILE=native`, 0 compiler warnings)

```text
-- The CXX compiler identification is GNU 11.4.0
-- Detecting CXX compiler ABI info
-- Detecting CXX compiler ABI info - done
-- Check for working CXX compiler: /usr/bin/c++ - skipped
-- Detecting CXX compile features
-- Detecting CXX compile features - done
-- LSQUANT build profile: native
-- Found OpenMP_CXX: -fopenmp (found version "4.5") 
-- Found OpenMP: TRUE (found version "4.5")  
-- Found PkgConfig: /usr/bin/pkg-config (found version "0.29.2") 
-- Found GSL: /usr/include (found version "2.7.1") 
-- Looking for C++ include pthread.h
-- Looking for C++ include pthread.h - found
-- Performing Test CMAKE_HAVE_LIBC_PTHREAD
-- Performing Test CMAKE_HAVE_LIBC_PTHREAD - Success
-- Found Threads: TRUE  
-- Checking for module 'fftw3'
--   Found fftw3, version 3.3.8
-- Found PythonInterp: /usr/bin/python3 (found suitable version "3.10.12", minimum required is "3.9") 
-- Found PythonLibs: /usr/lib/x86_64-linux-gnu/libpython3.10.so
-- Performing Test HAS_FLTO_AUTO
-- Performing Test HAS_FLTO_AUTO - Success
-- Found pybind11: /home/ICN2/jgarcia/.local/lib/python3.10/site-packages/pybind11/include (found version "3.1.0")
-- Configuring done
-- Generating done
-- Build files have been written to: /home/ICN2/jgarcia/lsquant-agi-build-0010/build-step1-native
[  3%] Building CXX object lsquant_core/CMakeFiles/lsquant_core.dir/core/algebra/composite_operators.cpp.o
[  9%] Building CXX object lsquant_core/CMakeFiles/lsquant_core.dir/core/algebra/spectral_bounds.cpp.o
[  9%] Building CXX object lsquant_core/CMakeFiles/lsquant_core.dir/core/algebra/sparse_matrix.cpp.o
[ 12%] Building CXX object lsquant_core/CMakeFiles/lsquant_core.dir/core/algebra/block_diagonal_operator.cpp.o
[ 16%] Building CXX object lsquant_core/CMakeFiles/lsquant_core.dir/core/chebyshev/chebyshev_vectors.cpp.o
[ 22%] Building CXX object lsquant_core/CMakeFiles/lsquant_core.dir/core/algebra/metric_policy.cpp.o
[ 22%] Building CXX object lsquant_core/CMakeFiles/lsquant_core.dir/core/algebra/bloch_transform.cpp.o
[ 25%] Building CXX object lsquant_core/CMakeFiles/lsquant_core.dir/core/algebra/state_generator.cpp.o
[ 29%] Building CXX object lsquant_core/CMakeFiles/lsquant_core.dir/core/algebra/linalg.cpp.o
[ 32%] Building CXX object lsquant_core/CMakeFiles/lsquant_core.dir/core/algebra/overlap_appliers.cpp.o
[ 35%] Building CXX object lsquant_core/CMakeFiles/lsquant_core.dir/core/chebyshev/chebyshev_recursion.cpp.o
[ 41%] Building CXX object lsquant_core/CMakeFiles/lsquant_core.dir/core/chebyshev/rescaled_operator.cpp.o
[ 45%] Building CXX object lsquant_core/CMakeFiles/lsquant_core.dir/core/chebyshev/chebyshev_moments.cpp.o
[ 48%] Building CXX object lsquant_core/CMakeFiles/lsquant_core.dir/core/chebyshev/moment_engine.cpp.o
[ 48%] Building CXX object lsquant_core/CMakeFiles/lsquant_core.dir/core/chebyshev/contraction.cpp.o
[ 51%] Building CXX object lsquant_core/CMakeFiles/lsquant_core.dir/core/algebra/fft_moment_contraction.cpp.o
[ 54%] Building CXX object lsquant_core/CMakeFiles/lsquant_core.dir/io/csr_text.cpp.o
[ 58%] Building CXX object lsquant_core/CMakeFiles/lsquant_core.dir/io/chebmom_file.cpp.o
[ 61%] Building CXX object lsquant_core/CMakeFiles/lsquant_core.dir/io/curve_text.cpp.o
[ 64%] Linking CXX static library ../lib/liblsquant_core.a
[ 64%] Built target lsquant_core
[ 67%] Building CXX object lsquant_core/CMakeFiles/stage6_nonorth.dir/apps/stage6_nonorth.cpp.o
[ 70%] Building CXX object lsquant_core/CMakeFiles/stage1_moments.dir/apps/stage1_moments.cpp.o
[ 74%] Building CXX object lsquant_core/CMakeFiles/_core.dir/api/bindings.cpp.o
[ 77%] Building CXX object lsquant_core/CMakeFiles/stage_1_check.dir/__/tests/dev/stage_1_check.cpp.o
[ 80%] Building CXX object lsquant_core/CMakeFiles/_core.dir/compat/legacy_kernels.cpp.o
[ 83%] Building CXX object lsquant_core/CMakeFiles/stage8_kernels.dir/__/bench/kernel_bench.cpp.o
[ 87%] Linking CXX executable stage8_kernels
[ 87%] Built target stage8_kernels
[ 90%] Linking CXX executable stage1_moments
[ 90%] Built target stage1_moments
[ 93%] Linking CXX executable stage6_nonorth
[ 93%] Built target stage6_nonorth
[ 96%] Linking CXX executable stage_1_check
[ 96%] Built target stage_1_check
[100%] Linking CXX shared module _core.cpython-310-x86_64-linux-gnu.so
[100%] Built target _core
```

### Baseline `stage_1_check` output

- `golden` profile:
```text
dos max_abs_error: 2.47025e-15
kubo max_abs_error: 1.13687e-13
OMP_NUM_THREADS=1 vs 8 bit-identical: YES
```

- `native` profile:
```text
dos max_abs_error: 1.77636e-15
kubo max_abs_error: 9.94766e-14
OMP_NUM_THREADS=1 vs 8 bit-identical: YES
```

---

## 2. Step 2 — Remove Dead Files (Commit `da6239f`)

### Files changed in `da6239f` (`git show --name-status da6239f`)

```text
M	.githooks/pre-commit
M	.gitignore
D	.projectile
D	.vscode/settings.json
D	GOLDEN_PROVENANCE.txt
M	README.md
D	TODO.md
M	agent_template/tutorial_style.md
M	agent_template/tutorial_template.md
D	bench/CAPTIONS.md
A	bench/block_domains/RESULTS.md
A	bench/block_domains/analyze.py
A	bench/block_domains/block_kpm.cpp
A	bench/block_domains/compare_moments.py
A	bench/block_domains/results.jsonl
A	bench/block_domains/run_scaling.sh
D	bench/make_graphene.py
D	bench/plot_graphene_scaling.py
D	bench/run_graphene_scaling.sh
M	docs/CMakeLists.txt
A	docs/derivations/README.md
R100	lsquant_derive_core.py	docs/derivations/lsquant_derive_core.py
R100	lsquant_derive_spin.py	docs/derivations/lsquant_derive_spin.py
R100	lsquant_derive_transport.py	docs/derivations/lsquant_derive_transport.py
D	docs/main.cpp
M	docs/performance/README.md
M	docs/reconstruction.md
M	docs/spectral_scales.md
M	examples/.gitignore
D	examples/style_demo.py
D	examples/utilities/kwant2wannier/graphene.py
D	examples/utilities/kwant2wannier/graphene_example.ipynb
D	infrastructure_plan.md
D	install.py
M	lsquant/core/algebra/bloch_transform.hpp
M	lsquant/core/algebra/composite_operators.hpp
M	python/lsquant/moments.py
R097	lsquant_chain_reference.py	python/lsquant/oracle/chain.py
D	reporting_system.md
D	scripts/check_reporting.sh
M	tests/dev/stage_4_check.py
D	unified_io.md
D	utilities/CMakeLists.txt
D	utilities/postprocessing/average_moments1D.cpp
D	utilities/python/.ipynb_checkpoints/qstate_generator-checkpoint.ipynb
D	utilities/python/CMakeLists.txt
D	utilities/python/bin/add_chebmom2D.py
D	utilities/python/bin/average_chebmom2D.py
D	utilities/python/bin/substract_chebmom2D.py
D	utilities/python/lib/__pycache__/sciplots.cpython-36.pyc
D	utilities/python/lib/sciplots.py
D	utilities/python/lsquant_report.py
D	utilities/python/test/bin/add_chebmom2D.py
D	utilities/python/test/bin/average_chebmom2D.py
D	utilities/python/test/bin/substract_chebmom2D.py
D	utilities/python/test/lib/python3.6/__pycache__/sciplots.cpython-36.pyc
D	utilities/python/test/lib/python3.6/sciplots.py
D	utilities/siesta2sparse/siesta2lsquant.py
D	utilities/siesta2sparse/siesta2lsquant_realspace.py
```

### `git grep -n <basename>` evidence at baseline (`77df68d`) outside `handoff/` and `goldens/`

```text
=== git grep -n TODO.md ===
(no matches)

=== git grep -n infrastructure_plan.md ===
TODO.md:73:      `infrastructure_plan.md` §6 now that the in-house route is implemented
TODO.md:96:- Plans reconciled: `infrastructure_plan.md`, `reporting_system.md`,

=== git grep -n reporting_system.md ===
README.md:139:[`reporting_system.md`](reporting_system.md)): `-v`/`--log-file` and
TODO.md:70:      gates, keeping the in-house fallback (see `reporting_system.md` §5,
TODO.md:96:- Plans reconciled: `infrastructure_plan.md`, `reporting_system.md`,
infrastructure_plan.md:7:- [`reporting_system.md`](reporting_system.md) — logging, timing/profiling, error
infrastructure_plan.md:53:| **Reporting** (produces) | `lsquant::util` / `diag` — `include/util/`, `src/util/` | log levels+sinks, scoped timers, `Error` hierarchy + `Result`, Linux-first resource probe behind a backend interface | `reporting_system.md` |
unified_io.md:5:[`reporting_system.md`](reporting_system.md) — timing, error, and resource data
unified_io.md:223:             "peak_rss_bytes": 4509715660,        // from reporting_system.md
unified_io.md:316:## 8. Relationship to `reporting_system.md`

=== git grep -n unified_io.md ===
TODO.md:71:      `unified_io.md` §5).
TODO.md:97:  `unified_io.md`.
infrastructure_plan.md:9:- [`unified_io.md`](unified_io.md) — unified config + validation, system
infrastructure_plan.md:54:| **Unified I/O** (consumes) | `lsquant::io` — `include/io/`, `src/io/` | schema v2 (system + observables), JSON-Schema validation + hints, SHA-256 system fingerprint, `results.json` container, `lsquant merge`, Python presenter + HTML | `unified_io.md` |

=== git grep -n GOLDEN_PROVENANCE.txt ===
(no matches)

=== git grep -n install.py ===
docs/main.md:18:installer helper is `install.py` at the repository root, and the full
install.py:4:python3 -m pip install pybind11

=== git grep -n lsquant_derive_core.py ===
(no matches)

=== git grep -n lsquant_derive_spin.py ===
(no matches)

=== git grep -n lsquant_derive_transport.py ===
(no matches)

=== git grep -n lsquant_chain_reference.py ===
AGENTS.md:353:| G2 chain1d Kubo-Greenwood moments | same, `VX`,`VX` | `NonEqOpVX-VXchain1d...chebmom2D` | closed form B.15 (`lsquant_chain_reference.py`) |
GOLDEN_PROVENANCE.txt:2:chain1d:  KPM_SEED=12345 w2s_sha=26cab4a lsquant_head=5dcce0f200aef6dad39bde257060cbb1c60d1d37 Nx=2048 Ny=1 Nz=1 M=64 bounds=[-2,2] KG=NonEqOpVX-VXchain1dKPM_M64x64_stateexact.chebmom2D(ops=VX,VX) DOS=SpectralOp1chain1dKPM_M64_stateexact.chebmom1D(op=1) ref=oracle:lsquant_chain_reference.py

=== git grep -n utilities/python / chebmom2D.py / sciplots.py / lsquant_report.py ===
README.md:136:[`utilities/python/lsquant_report.py`](utilities/python/lsquant_report.py).
TODO.md:64:      `lsquant_report.py` path. Add the equivalent short note to tutorials
TODO.md:90:- Presenter `utilities/python/lsquant_report.py` (summary + plots +
examples/03_haldane_hall/README.md:110:> python ../../utilities/python/lsquant_report.py haldane.results.json   # table + HTML
examples/utilities/kwant2wannier/graphene_example.ipynb:20:    "sys.path.append('/data/jgarcia/codes/lsquant/utilities/python/tight-binding/')\n",
infrastructure_plan.md:72:utilities/python/  lsquant_report.py   (+ reuse lib/sciplots.py)
infrastructure_plan.md:120:| **I6** | `utilities/python/lsquant_report.py` + self-contained Plotly HTML dashboard | I/O |
infrastructure_plan.md:142:| E6 | `utilities/python` cleanup (build/ cruft, py3.6) | separate housekeeping pass |
unified_io.md:40:  Python `add/average/substract_chebmom2D.py` remain under `utilities/python` and
unified_io.md:61:7. **A Python presenter** in `utilities/python` that ingests the output and shows
unified_io.md:72:| Merge | **Already in embryo**: `utilities/python/bin/{add,average,substract}_chebmom2D.py` merge by requiring the **header string to be byte-identical** (`add_chebmom2D.py`). |
unified_io.md:74:| Python utils | `utilities/python/` has `sciplots.py` + the merge scripts, but carries **committed `build/` dirs** and py3.6 `__pycache__` cruft. |
unified_io.md:103:        utilities/python              static HTML dashboard
unified_io.md:104:        lsquant_report.py  ──────────▶  <system>.report.html   (→ online site later)
unified_io.md:248:### 4.5 Presentation: Python tool in `utilities/python` + HTML dashboard
unified_io.md:250:- **`utilities/python/lsquant_report.py`** (new) — ingests one or many
unified_io.md:254:  2. plots DOS / σ_xx / σ_xy / spin-relaxation (reuse `lib/sciplots.py`);
unified_io.md:264:- **Cleanup (flagged, not bundled):** `utilities/python` has committed `build/`
unified_io.md:283:  `lsquant_report.py` there; add a `requirements.txt` (numpy, plotly).
unified_io.md:302:- **Phase 6:** `utilities/python/lsquant_report.py` + static HTML dashboard.
unified_io.md:314:| E6 | utilities/python cleanup (build/ cruft, py3.6) | separate housekeeping pass, not in this feature |
utilities/python/lsquant_report.py:15:    lsquant_report.py graphene.results.json                # summary + HTML
utilities/python/lsquant_report.py:16:    lsquant_report.py *.results.json -o dashboard.html     # several systems
utilities/python/lsquant_report.py:17:    lsquant_report.py run.results.json --no-html           # summary only
Binary file utilities/python/lib/__pycache__/sciplots.cpython-36.pyc matches
Binary file utilities/python/test/lib/python3.6/__pycache__/sciplots.cpython-36.pyc matches

=== git grep -n average_moments1D.cpp ===
(no matches)

=== git grep -n utilities/CMakeLists.txt ===
(no matches)

=== git grep -n siesta2sparse / siesta2lsquant.py / siesta2lsquant_realspace.py ===
AGENTS.md:301:  `utilities/siesta2sparse/siesta2lsquant.py`). Hermitian, so the standard engine,
utilities/CMakeLists.txt:3:# github.com/adamecius/wannier2sparse (see LSQUANT_BUILD_INPUT_TOOLS). siesta2sparse
utilities/siesta2sparse/siesta2lsquant.py:3:siesta2lsquant.py
utilities/siesta2sparse/siesta2lsquant.py:29:    python siesta2lsquant.py <fdf_file> <kx> <ky> [--prefix NAME]

=== git grep -n lsqplot.py / lsqloc.py / lsqhall.py / lsqspin.py / lsqtransport.py / lsqrelax.py / plot_ratelaw.py ===
agent_template/tutorial_style.md:85:by calling the compiled kernels through `lsqplot.py` or, where the physics has a
agent_template/tutorial_template.md:94:[[python -c "import lsqplot; lsqplot.figure_*( ... )"   or   python lsqplot.py]]
examples/01_dos_1d_chain/README.md:150:python lsqplot.py
examples/02_localization_1d_chain/README.md:104:python lsqloc.py
examples/03_haldane_hall/README.md:126:python lsqhall.py
examples/03_haldane_hall/haldane_kspace.py:11:The real-space LSQUANT calculation (see lsqhall.py / README) reproduces the same
examples/04_spin_precession_magnetic_chain/README.md:92:python lsqspin.py
examples/05_charge_transport_square/README.md:118:python lsqtransport.py
examples/05_charge_transport_square/make_square.py:20:functions are importable: lsqtransport.py calls build_square with a fresh seed per
examples/06_spin_relaxation_square/README.md:70:python lsqrelax.py        # runs Part A (this figure) and Part B (Step 4)
examples/06_spin_relaxation_square/README.md:97:`lsqrelax.py` also sweeps the disorder and compares the spin directions:
examples/06_spin_relaxation_square/README.md:145:python3 plot_ratelaw.py                     # writes fig_ey_ratelaw
examples/06_spin_relaxation_square/README.md:167:> engine (`python -m lsquant run ../../runs/goldens/tutorial06.json`, or `lsqrelax.py`
figures/stage_0.py:71:            ("py", ["lsqplot.py"]),
figures/stage_0.py:79:        "size": "N = 512 (README step 2); N = 256, 1024, 4096 x 16 realisations (lsqloc.py)",
figures/stage_0.py:81:        "seed": "KPM_SEED=12345 for the README CLI step; lsqloc.py sets KPM_SEED = disorder seed 1..16",
figures/stage_0.py:88:            ("py", ["lsqloc.py"]),
figures/stage_0.py:110:            ("py", ["lsqhall.py"]),
figures/stage_0.py:126:            ("py", ["lsqspin.py"]),
figures/stage_0.py:138:        "seed": "KPM_SEED=12345 for the README CLI step; lsqtransport.py sets KPM_SEED=1",
figures/stage_0.py:145:            ("py", ["lsqtransport.py"]),
figures/stage_0.py:149:                 "defines; it is not run standalone. lsqtransport.py writes and runs that "
figures/stage_0.py:160:            ("py", ["lsqrelax.py"]),
figures/stage_0.py:163:            ("py", ["plot_ratelaw.py"]),
figures/stage_2.py:9,16,43: comments referencing lsqplot.py normalisation and capture order
figures/stage_3.py:11,59,179: comments referencing lsqhall.py normalisation and chern_from_plateau
figures/stage_4.py:20,22,55,56,128,137,142,164,175,373,380: comments referencing lsqspin.py, lsqloc.py, lsqtransport.py
python/lsquant/README.md:33:| `make_chain_mag` | `examples/04_spin_precession_magnetic_chain/lsqspin.py` |
python/lsquant/runfile.py:76:    "make_chain_mag": "examples/04_spin_precession_magnetic_chain/lsqspin.py",
python/lsquant/apps/tutorial01.py:180: comment on capture order
python/lsquant/apps/tutorial02.py:6,210: comments on lsqloc.py ensemble
python/lsquant/apps/tutorial03.py:103,177,238: comments on lsqhall.py normalisation

=== git grep -n style_demo.py ===
TODO.md:46:      `python3 examples/style_demo.py`). Keep `examples/prl_style.py` and
TODO.md:47:      `examples/style_demo.py`.
TODO.md:92:- PRL publication figure style `examples/prl_style.py` (+ `style_demo.py`); all

=== git grep -n prl_style.py ===
README.md:170:[`examples/prl_style.py`](examples/prl_style.py): fixed column widths, math
agent_template/plotting_style.md:8:The reference implementation is `prl_style.py` (presets `single` / `double` /
examples/01..06/lsq*.py: # Shared publication (APS/PRL) plot style: examples/prl_style.py.
figures/stage_0.py:11,12,239: copies examples/prl_style.py

=== git grep -n kwant2wannier ===
bench/make_graphene.py:4:benchmark. Self-contained (no kwant): the official examples/utilities/kwant2wannier
docs/performance/README.md:30:`examples/utilities/kwant2wannier` path needs kwant, absent here). Thread sweep at
examples/utilities/kwant2wannier/graphene_example.ipynb:23:    "from kwant2wannier import syst2txt, geo2txt\n",

=== git grep -n run_graphene_scaling.sh / make_graphene.py / plot_graphene_scaling.py / CAPTIONS.md ===
bench/README.md:248:numbers from the older `run_graphene_scaling.sh` (which targets a different node).
bench/make_graphene.py:19:Usage:  make_graphene.py <L> [t2] [phi_over_pi] [outdir]
bench/plot_graphene_scaling.py:3,79,102: Reads graphene_threads.csv + graphene_size.csv (from run_graphene_scaling.sh), writes CAPTIONS.md
bench/run_graphene_scaling.sh:11:GEN="python3 make_graphene.py"
docs/performance/README.md:7,8,29,72: references bench/run_graphene_scaling.sh, make_graphene.py, plot_graphene_scaling.py

=== git grep -n Doxyfile.in / docs/main.cpp / docs/main.md / docs/CMakeLists.txt / docs/requirements.txt ===
docs/CMakeLists.txt:1,5,19,29,36,51: configures docs/Doxyfile.in and docs/main.md
docs/Doxyfile.in:1,6: configured by docs/CMakeLists.txt
(no matches for docs/main.cpp or docs/requirements.txt)

=== git grep -n downfold_origin / lsquant/projection / runs/stage6c / runs/projection / stage_9_ ===
python/lsquant/projection/origin.py:19,21,29: loads tools/downfold_origin
tools/downfold_origin/README.md:20,21,25: references tools/downfold_origin and python/lsquant/projection/
tools/downfold_origin/SHA256SUMS:1..32: checksums of tools/downfold_origin/*.py
figures/out/review_20261025/run.sh:6: stage_9_material_review.py
```

### Items Required by Step 2 Table

1. **`GOLDEN_PROVENANCE.txt` contents (deleted; preserved here for Fable):**
   ```text
   graphene: KPM_SEED=12345 w2s_sha=26cab4a lsquant_head=5dcce0f200aef6dad39bde257060cbb1c60d1d37 Nx=8 Ny=8 Nz=1 M=64 ops=VX,VX
   chain1d:  KPM_SEED=12345 w2s_sha=26cab4a lsquant_head=5dcce0f200aef6dad39bde257060cbb1c60d1d37 Nx=2048 Ny=1 Nz=1 M=64 bounds=[-2,2] KG=NonEqOpVX-VXchain1dKPM_M64x64_stateexact.chebmom2D(ops=VX,VX) DOS=SpectralOp1chain1dKPM_M64_stateexact.chebmom1D(op=1) ref=oracle:lsquant_chain_reference.py
   ```
   (`goldens/stage0/PROVENANCE.md` already documents `graphene` and `chain1d` at `legacy-final` `05e2f65`; the only additional datum in `GOLDEN_PROVENANCE.txt` was the earlier commit SHA `5dcce0f200aef6dad39bde257060cbb1c60d1d37`.)

2. **`lsquant_chain_reference.py` -> `python/lsquant/oracle/chain.py` (`AGENTS.md` update for Fable):**
   Moved to `python/lsquant/oracle/chain.py` with CLI entry point `python -m lsquant.oracle.chain`. In `AGENTS.md` section 8 (line 353), replace:
   ```markdown
   | G2 chain1d Kubo-Greenwood moments | same, `VX`,`VX` | `NonEqOpVX-VXchain1d...chebmom2D` | closed form B.15 (`lsquant_chain_reference.py`) |
   ```
   with:
   ```markdown
   | G2 chain1d Kubo-Greenwood moments | same, `VX`,`VX` | `NonEqOpVX-VXchain1d...chebmom2D` | closed form B.15 (`python/lsquant/oracle/chain.py`, `python -m lsquant.oracle.chain`) |
   ```

3. **`utilities/siesta2sparse/` replacement mapping:**
   Deleted both `utilities/siesta2sparse/siesta2lsquant.py` and `utilities/siesta2sparse/siesta2lsquant_realspace.py` because they are completely superseded by the library:
   - `siesta2lsquant.py` (k-space Bloch matrices $H(k)$, $S(k)$, $\partial_a H(k)$, $\partial_a S(k)$, spin matrices, and Löwdin / covariant velocities) is replaced by `python/lsquant/models/siesta.py`: `load_siesta_model`, `_evaluate_bloch_matrices`, `_compute_dipole_correction_orbital_blocks`, `_build_spin_operators`, and `python/lsquant/bases.py::BlochBasis` / `python/lsquant/metric_forms.py`.
   - `siesta2lsquant_realspace.py` (real-space supercell tiling of SIESTA $H$, $S$, and velocity operators) is replaced by `python/lsquant/models/nonorthogonal.py::build_siesta_rspace_supercell` and `python/lsquant/disorder_bases.py::build_real_supercell_session`.

4. **`docs/` disposition:**
   - Kept and updated: `docs/Doxyfile.in` (already points `INPUT` to `lsquant/` and `python/lsquant/`), `docs/CMakeLists.txt` (removed deleted `docs/main.cpp`), `docs/main.md` (removed `install.py` reference and updated to CMake presets and `python -m lsquant`), `docs/requirements.txt`, `docs/performance/README.md`, `docs/reconstruction.md`, and `docs/spectral_scales.md` (removed all mentions of `.chebmom1D/2D/TD` as user formats, `inline_*`, `lsquant compute`, and `lsquant reconstruct <file> dos`).
   - Deleted: `docs/main.cpp` (empty legacy stub).

5. **Whole-tree search for untracked/dead artefacts:**
   - `git ls-files '*.pyc' '*.orig' '*checkpoint*'`: found `utilities/python/lib/__pycache__/sciplots.cpython-36.pyc`, `utilities/python/test/lib/python3.6/__pycache__/sciplots.cpython-36.pyc`, and `utilities/python/.ipynb_checkpoints/qstate_generator-checkpoint.ipynb` — all deleted with `utilities/python/`.
   - Editor files in `git ls-files`: found `.projectile` and `.vscode/settings.json` (which contained CMake/legacy header paths from the old tree) — both deleted.
   - `scripts/check_reporting.sh` and `.githooks/pre-commit`: `scripts/check_reporting.sh` checked the deleted `src/` and `include/` directories for the abandoned reporting plan; deleted `scripts/check_reporting.sh` and updated `.githooks/pre-commit` so it only invokes `check_reporting.sh` if present.
   - Tracked files larger than 1 MB outside `goldens/` and `examples/*/siesta/`: only `figures/out/review_20261025/graphene_disorder_compare.npz` (1.6 MB, part of the 0009 wrap-up review figures in `figures/`; left untouched per Section 1 restriction on `figures/`).
   - Unreferenced Python modules under `python/lsquant/`: found `python/lsquant/projection/diagnostics.py` (never imported anywhere in `python/lsquant/`, `tests/`, `runs/`, or `figures/`; only listed in `projection/__init__.py` `__all__`). Deleted in Step 3 (`a8252cc`) and removed from `projection/__init__.py`.
   - Uncalled C++ functions under `lsquant/`: none; every public method in `lsquant/core/` and `lsquant/io/` is bound in `lsquant/api/bindings.cpp` or called by `MomentEngine`, `stage1_moments`, `stage6_nonorth`, or `stage_1_check`.

6. **`inline_*` mentions removed from the 5 non-`AGENTS.md` files:**
   Removed all mentions of `inline_compute-kpm-*`, `device-postProcess`, `MeanSquareDisplacement_DeviceSim`, and `kuboBastinFromChebmom` from `README.md`, `agent_template/tutorial_style.md`, `agent_template/tutorial_template.md`, `python/lsquant/moments.py`, and `tests/dev/stage_4_check.py`. Running:
   ```bash
   git grep -n -E '(inline_compute-kpm-|device-postProcess|MeanSquareDisplacement_DeviceSim|kuboBastinFromChebmom)' -- ':!handoff/' ':!goldens/' ':!AGENTS.md'
   ```
   now returns 0 matches.

---

## 3. Step 3 — Code Review Pass (Commit `a8252cc`)

### Files changed in `a8252cc` (`git show --name-status a8252cc`)

```text
M	lsquant/api/bindings.cpp
M	lsquant/core/algebra/block_diagonal_operator.cpp
M	lsquant/core/algebra/build_profile.hpp
M	lsquant/core/algebra/fft_moment_contraction.cpp
M	lsquant/core/algebra/fft_moment_contraction.hpp
M	lsquant/core/algebra/fftw_planner.hpp
M	lsquant/core/algebra/linalg.hpp
M	lsquant/core/algebra/linear_operator.hpp
M	lsquant/core/algebra/memory_placement.hpp
M	lsquant/core/algebra/sparse_matrix.hpp
M	lsquant/core/algebra/work_partition.hpp
M	lsquant/core/chebyshev/chebyshev_moments.hpp
M	lsquant/core/chebyshev/contraction.cpp
M	lsquant/core/chebyshev/plan.hpp
M	lsquant/io/chebmom_file.cpp
M	lsquant/io/chebmom_file.hpp
M	lsquant/io/csr_text.hpp
M	python/lsquant/apps/stage6b_siesta.py
M	python/lsquant/apps/stage6c_disorder.py
M	python/lsquant/apps/tutorial01.py
M	python/lsquant/bases.py
M	python/lsquant/disorder_bases.py
M	python/lsquant/kernels.py
M	python/lsquant/metric_forms.py
M	python/lsquant/models/disorder.py
M	python/lsquant/models/nonorthogonal.py
M	python/lsquant/models/siesta.py
M	python/lsquant/models/wannier.py
M	python/lsquant/moments.py
M	python/lsquant/oracle/spin_hall.py
M	python/lsquant/projection/__init__.py
M	python/lsquant/projection/bases.py
D	python/lsquant/projection/diagnostics.py
M	python/lsquant/projection/origin.py
M	python/lsquant/projection/orthogonal.py
M	python/lsquant/projection/real_space.py
M	python/lsquant/projection/runfile.py
M	python/lsquant/projection/words.py
M	python/lsquant/reconstruct.py
M	python/lsquant/runfile.py
M	python/lsquant/runtime.py
M	python/lsquant/words/catalogue.py
M	python/lsquant/words/compiler.py
M	scripts/golden_diff.py
```

### 1. Identifiers and Named Constants (Before -> After)

- **Library module / symbol cleanup:**
  - `python/lsquant/projection/diagnostics.py`: deleted unused review-only diagnostic module and removed `"diagnostics"` from `python/lsquant/projection/__init__.py`.
  - Gate entry points (`stage1_moments`, `stage6_nonorth`, `stage_1_check`, `python -m lsquant.apps.stage*`, `python -m lsquant.apps.tutorial0*`, `tests/dev/stage_6b_check.py`, `stage_6c_check.py`, `stage_6d_check.py`) were preserved unchanged per Section 1 of the brief.
- **Named constants extracted from magic numbers:**

| File | Before (hard-coded literal) | After (named constant) | Value |
|---|---|---|---|
| `python/lsquant/models/siesta.py` | `1e-12` (quadrature grid radius cutoff) | `DEFAULT_DIPOLE_MIN_RADIUS_BOHR` | `1.0e-12` |
| `python/lsquant/models/siesta.py` | `48` (radial Gauss-Legendre points) | `DEFAULT_DIPOLE_NRADIAL` | `48` |
| `python/lsquant/models/siesta.py` | `16` (polar Gauss-Legendre points) | `DEFAULT_DIPOLE_NTHETA` | `16` |
| `python/lsquant/models/siesta.py` | `32` (azimuthal trapezoidal points) | `DEFAULT_DIPOLE_NPHI` | `32` |
| `python/lsquant/models/siesta.py` | `root / ".lsquant_dipole_cache"` (inside `siesta/`) | `_default_dipole_cache_dir(root)` (`root.parent / ".lsquant-cache"` when `root.name == "siesta"`) | outside `siesta/` |
| `python/lsquant/oracle/spin_hall.py` | `1e-12` (degenerate eigenvalue cutoff in Kubo-Bastin sum) | `DEFAULT_DEGENERACY_TOL_EV` | `1.0e-12` |
| `python/lsquant/metric_forms.py` | `1e-13` (default PCG tolerance) | `DEFAULT_PCG_TOLERANCE` | `1.0e-13` |
| `python/lsquant/metric_forms.py` | `200` (default PCG max iterations) | `DEFAULT_PCG_MAX_ITER` | `200` |
| `python/lsquant/metric_forms.py` | `60` (default Chebyshev overlap polynomial order) | `DEFAULT_OVERLAP_CHEBYSHEV_ORDER` | `60` |
| `python/lsquant/metric_forms.py` | `1e-12` (default Chebyshev coefficient cutoff) | `DEFAULT_OVERLAP_COEFF_TOL` | `1.0e-12` |
| `python/lsquant/metric_forms.py` | `4096` (default maximum dimension for dense exact overlap inverse) | `DEFAULT_EXACT_OVERLAP_MAX_DIM` | `4096` |

### 2. Comments and Error Messages

- Removed all stage/handoff/review/tester scaffolding comments (`Stage-1`, `Stage-6b`, `Stage-6c`, `Stage-6d`, `Stage-8`, `handoff 000*`) across 43 C++ and Python files while preserving all mathematical, physical, unit, and legacy-quirk explanations (AGENTS.md section 10).
- Improved `throw` / `raise` error messages in `lsquant/io/chebmom_file.cpp`, `python/lsquant/projection/origin.py`, and `python/lsquant/models/siesta.py` so each states what was expected and what was received.
- Removed unused imports flagged by `pyflakes` across `python/lsquant/apps/stage6b_siesta.py`, `python/lsquant/apps/stage6c_disorder.py`, `python/lsquant/disorder_bases.py`, `python/lsquant/projection/bases.py`, and `python/lsquant/projection/real_space.py`. `python3 -m pyflakes python/lsquant` exits 0 with zero output.

---

## 4. Step 4 — CMake, Presets, and Install Rules (Commit `26313dc`)

### Files changed in `26313dc` (`git show --name-status 26313dc`)

```text
M	CMakeLists.txt
A	CMakePresets.json
M	lsquant/CMakeLists.txt
M	tests/dev/stage_6_check.py
M	tests/dev/stage_7_check.py
```

### Summary of CMake changes

- Single top-level `CMakeLists.txt` defining options `LSQUANT_BUILD_PROFILE` (`golden` | `native`, default `golden`), `LSQUANT_BUILD_PYTHON` (`ON`), `LSQUANT_BUILD_BENCHMARKS` (`OFF`), and `LSQUANT_BUILD_DEV_CHECKS` (`ON`), with `CMAKE_EXPORT_COMPILE_COMMANDS ON`.
- Removed `LSQUANT_BUILD_NEW_CORE` and `LSQUANT_BUILD_INPUT_TOOLS`.
- Dependencies (`OpenMP`, `FFTW3`, `Eigen3`, and when `LSQUANT_BUILD_PYTHON=ON`, `GSL` and `pybind11`) are found once at the top level.
- Added `lsquant_warnings` interface target (`-Wall -Wextra -Wpedantic`) shared across all targets without per-target flag repetition.
- Added `install()` rules for `lsquant_core`, public headers under `include/lsquant/`, `stage1_moments` and `stage6_nonorth` under `bin/`, and `_core` plus the `lsquant` Python package under `python/lsquant`.
- Created `CMakePresets.json` with `golden` and `native` configure and build presets (`build/golden` and `build/native`).
- Updated `tests/dev/stage_6_check.py` and `tests/dev/stage_7_check.py` to detect `build/golden` or `build/native` automatically when `build/lsquant_core` is not present.

### Verbatim preset build, `cmake --install`, and out-of-tree import output on TCN

Commands:
```bash
cmake --preset golden
cmake --build --preset golden
cmake --preset native
cmake --build --preset native
rm -rf /tmp/lsquant-install-test
cmake --install build/native --prefix /tmp/lsquant-install-test
cd /tmp && PYTHONPATH=/tmp/lsquant-install-test/python python3 -c \
  "import lsquant, lsquant._core; print('Out-of-tree import OK:', lsquant.__file__, lsquant._core.__file__)"
```

Output (`~/lsquant-agi-build-0010/logs/verify_part1.log` lines 1–208):
```text
=== 1. Presets build & install ===
Preset CMake variables:

  CMAKE_BUILD_TYPE="Release"
  LSQUANT_BUILD_BENCHMARKS="OFF"
  LSQUANT_BUILD_DEV_CHECKS="ON"
  LSQUANT_BUILD_PROFILE="golden"
  LSQUANT_BUILD_PYTHON="ON"

-- The CXX compiler identification is GNU 11.4.0
-- Detecting CXX compiler ABI info
-- Detecting CXX compiler ABI info - done
-- Check for working CXX compiler: /usr/bin/c++ - skipped
-- Detecting CXX compile features
-- Detecting CXX compile features - done
-- LSQUANT build profile: golden
-- Found FFTW3 OpenMP library: /usr/lib/x86_64-linux-gnu/libfftw3_omp.so
-- Found PkgConfig: /usr/bin/pkg-config (found version "0.29.2") 
-- Found PythonLibs: /usr/lib/x86_64-linux-gnu/libpython3.10.so
-- Performing Test HAS_FLTO_AUTO
-- Performing Test HAS_FLTO_AUTO - Success
-- Configuring done
-- Generating done
-- Build files have been written to: /home/ICN2/jgarcia/lsquant-agi-build-0010/build/golden
[  3%] Building CXX object lsquant_core/CMakeFiles/lsquant_core.dir/core/algebra/block_diagonal_operator.cpp.o
[  6%] Building CXX object lsquant_core/CMakeFiles/lsquant_core.dir/core/algebra/bloch_transform.cpp.o
[ 17%] Building CXX object lsquant_core/CMakeFiles/lsquant_core.dir/core/algebra/linalg.cpp.o
[ 17%] Building CXX object lsquant_core/CMakeFiles/lsquant_core.dir/core/algebra/sparse_matrix.cpp.o
[ 17%] Building CXX object lsquant_core/CMakeFiles/lsquant_core.dir/core/algebra/state_generator.cpp.o
[ 20%] Building CXX object lsquant_core/CMakeFiles/lsquant_core.dir/core/algebra/spectral_bounds.cpp.o
[ 24%] Building CXX object lsquant_core/CMakeFiles/lsquant_core.dir/core/algebra/composite_operators.cpp.o
[ 27%] Building CXX object lsquant_core/CMakeFiles/lsquant_core.dir/core/algebra/metric_policy.cpp.o
[ 31%] Building CXX object lsquant_core/CMakeFiles/lsquant_core.dir/core/algebra/overlap_appliers.cpp.o
[ 34%] Building CXX object lsquant_core/CMakeFiles/lsquant_core.dir/core/chebyshev/rescaled_operator.cpp.o
[ 37%] Building CXX object lsquant_core/CMakeFiles/lsquant_core.dir/core/chebyshev/chebyshev_recursion.cpp.o
[ 41%] Building CXX object lsquant_core/CMakeFiles/lsquant_core.dir/core/chebyshev/chebyshev_vectors.cpp.o
[ 44%] Building CXX object lsquant_core/CMakeFiles/lsquant_core.dir/core/chebyshev/chebyshev_moments.cpp.o
[ 48%] Building CXX object lsquant_core/CMakeFiles/lsquant_core.dir/core/chebyshev/moment_engine.cpp.o
[ 51%] Building CXX object lsquant_core/CMakeFiles/lsquant_core.dir/core/chebyshev/contraction.cpp.o
[ 55%] Building CXX object lsquant_core/CMakeFiles/lsquant_core.dir/core/algebra/fft_moment_contraction.cpp.o
[ 58%] Building CXX object lsquant_core/CMakeFiles/lsquant_core.dir/io/csr_text.cpp.o
[ 62%] Building CXX object lsquant_core/CMakeFiles/lsquant_core.dir/io/chebmom_file.cpp.o
[ 65%] Building CXX object lsquant_core/CMakeFiles/lsquant_core.dir/io/curve_text.cpp.o
[ 68%] Linking CXX static library ../lib/liblsquant_core.a
[ 68%] Built target lsquant_core
[ 72%] Building CXX object lsquant_core/CMakeFiles/_core.dir/api/bindings.cpp.o
[ 75%] Building CXX object lsquant_core/CMakeFiles/stage1_moments.dir/apps/stage1_moments.cpp.o
[ 79%] Building CXX object lsquant_core/CMakeFiles/_core.dir/compat/legacy_kernels.cpp.o
[ 86%] Building CXX object lsquant_core/CMakeFiles/stage_1_check.dir/__/tests/dev/stage_1_check.cpp.o
[ 86%] Building CXX object lsquant_core/CMakeFiles/stage6_nonorth.dir/apps/stage6_nonorth.cpp.o
[ 89%] Linking CXX executable stage1_moments
[ 89%] Built target stage1_moments
[ 93%] Linking CXX executable stage6_nonorth
[ 93%] Built target stage6_nonorth
[ 96%] Linking CXX executable stage_1_check
[ 96%] Built target stage_1_check
[100%] Linking CXX shared module _core.cpython-310-x86_64-linux-gnu.so
[100%] Built target _core
Preset CMake variables:

  CMAKE_BUILD_TYPE="Release"
  LSQUANT_BUILD_BENCHMARKS="OFF"
  LSQUANT_BUILD_DEV_CHECKS="ON"
  LSQUANT_BUILD_PROFILE="native"
  LSQUANT_BUILD_PYTHON="ON"

-- The CXX compiler identification is GNU 11.4.0
-- Detecting CXX compiler ABI info
-- Detecting CXX compiler ABI info - done
-- Check for working CXX compiler: /usr/bin/c++ - skipped
-- Detecting CXX compile features
-- Detecting CXX compile features - done
-- LSQUANT build profile: native
-- Found FFTW3 OpenMP library: /usr/lib/x86_64-linux-gnu/libfftw3_omp.so
-- Found PkgConfig: /usr/bin/pkg-config (found version "0.29.2") 
-- Found PythonLibs: /usr/lib/x86_64-linux-gnu/libpython3.10.so
-- Performing Test HAS_FLTO_AUTO
-- Performing Test HAS_FLTO_AUTO - Success
-- Configuring done
-- Generating done
-- Build files have been written to: /home/ICN2/jgarcia/lsquant-agi-build-0010/build/native
[  3%] Building CXX object lsquant_core/CMakeFiles/lsquant_core.dir/core/algebra/sparse_matrix.cpp.o
[  6%] Building CXX object lsquant_core/CMakeFiles/lsquant_core.dir/core/algebra/block_diagonal_operator.cpp.o
[ 10%] Building CXX object lsquant_core/CMakeFiles/lsquant_core.dir/core/algebra/bloch_transform.cpp.o
[ 17%] Building CXX object lsquant_core/CMakeFiles/lsquant_core.dir/core/algebra/linalg.cpp.o
[ 17%] Building CXX object lsquant_core/CMakeFiles/lsquant_core.dir/core/algebra/composite_operators.cpp.o
[ 20%] Building CXX object lsquant_core/CMakeFiles/lsquant_core.dir/core/algebra/metric_policy.cpp.o
[ 24%] Building CXX object lsquant_core/CMakeFiles/lsquant_core.dir/core/algebra/spectral_bounds.cpp.o
[ 27%] Building CXX object lsquant_core/CMakeFiles/lsquant_core.dir/core/algebra/state_generator.cpp.o
[ 31%] Building CXX object lsquant_core/CMakeFiles/lsquant_core.dir/core/algebra/overlap_appliers.cpp.o
[ 34%] Building CXX object lsquant_core/CMakeFiles/lsquant_core.dir/core/chebyshev/rescaled_operator.cpp.o
[ 37%] Building CXX object lsquant_core/CMakeFiles/lsquant_core.dir/core/chebyshev/chebyshev_recursion.cpp.o
[ 41%] Building CXX object lsquant_core/CMakeFiles/lsquant_core.dir/core/chebyshev/chebyshev_vectors.cpp.o
[ 44%] Building CXX object lsquant_core/CMakeFiles/lsquant_core.dir/core/chebyshev/chebyshev_moments.cpp.o
[ 48%] Building CXX object lsquant_core/CMakeFiles/lsquant_core.dir/core/chebyshev/moment_engine.cpp.o
[ 51%] Building CXX object lsquant_core/CMakeFiles/lsquant_core.dir/core/chebyshev/contraction.cpp.o
[ 55%] Building CXX object lsquant_core/CMakeFiles/lsquant_core.dir/core/algebra/fft_moment_contraction.cpp.o
[ 58%] Building CXX object lsquant_core/CMakeFiles/lsquant_core.dir/io/csr_text.cpp.o
[ 62%] Building CXX object lsquant_core/CMakeFiles/lsquant_core.dir/io/chebmom_file.cpp.o
[ 65%] Building CXX object lsquant_core/CMakeFiles/lsquant_core.dir/io/curve_text.cpp.o
[ 68%] Linking CXX static library ../lib/liblsquant_core.a
[ 68%] Built target lsquant_core
[ 72%] Building CXX object lsquant_core/CMakeFiles/stage6_nonorth.dir/apps/stage6_nonorth.cpp.o
[ 75%] Building CXX object lsquant_core/CMakeFiles/stage_1_check.dir/__/tests/dev/stage_1_check.cpp.o
[ 79%] Building CXX object lsquant_core/CMakeFiles/stage1_moments.dir/apps/stage1_moments.cpp.o
[ 82%] Building CXX object lsquant_core/CMakeFiles/_core.dir/api/bindings.cpp.o
[ 86%] Building CXX object lsquant_core/CMakeFiles/_core.dir/compat/legacy_kernels.cpp.o
[ 89%] Linking CXX executable stage1_moments
[ 89%] Built target stage1_moments
[ 93%] Linking CXX executable stage6_nonorth
[ 93%] Built target stage6_nonorth
[ 96%] Linking CXX executable stage_1_check
[ 96%] Built target stage_1_check
[100%] Linking CXX shared module _core.cpython-310-x86_64-linux-gnu.so
[100%] Built target _core
-- Install configuration: "Release"
-- Installing: /tmp/lsquant-install-test/lib/liblsquant_core.a
-- Installing: /tmp/lsquant-install-test/include/lsquant/core
...
-- Installing: /tmp/lsquant-install-test/bin/stage1_moments
-- Installing: /tmp/lsquant-install-test/bin/stage6_nonorth
-- Installing: /tmp/lsquant-install-test/python/lsquant
...
-- Installing: /tmp/lsquant-install-test/python/lsquant/_core.cpython-310-x86_64-linux-gnu.so
Out-of-tree import OK: /tmp/lsquant-install-test/python/lsquant/__init__.py /tmp/lsquant-install-test/python/lsquant/_core.cpython-310-x86_64-linux-gnu.so
```

---

## 5. Step 5 — README Rewrite (Commit `0297e79`)

### Files changed in `0297e79` (`git show --name-status 0297e79`)

```text
M	README.md
M	python/lsquant/README.md
```

- Rewrote `README.md` preserving the section structure of the original README (what it computes, references, what is new, performance, installation, interface, tutorials and examples 01–08, how it works, inputs, verification, documentation, repository layout, licence and citation).
- Documented the word language (`Tr[A_0 f_1(H) A_1 ... f_p(H) A_p]`), the catalogue (`python -m lsquant words`), the five-instruction C++ engine (`ApplyOperator`, `Expand`, `ApplySeries`, `MirrorToBra`, `Trace`), the three bases (real, Bloch with FFT disorder, eigenbasis oracle), the non-orthogonal metric (`orthogonal`, `lowdin`, `generalised`), SIESTA and Wannier90 model loaders, build profiles (`golden` and `native`), corrected legacy defects (batched rank-2 shift, off-centre energy axis, tutorial 06 explicit identity operator), and performance numbers cited strictly from `bench/block_domains/RESULTS.md` and `handoff/0008-performance/report.md`.
- Rewrote `python/lsquant/README.md` as a concise API reference pointing to the root `README.md`.

---

## 6. Step 6 — Examples 01 to 08 (Commit `360339a`)

### Files changed in `360339a` (`git show --name-status 360339a`)

```text
M	agent_template/tutorial_style.md
M	agent_template/tutorial_template.md
M	examples/01_dos_1d_chain/README.md
D	examples/01_dos_1d_chain/fig_dos.png
D	examples/01_dos_1d_chain/fig_levels.png
D	examples/01_dos_1d_chain/fig_moments.png
D	examples/01_dos_1d_chain/lsqplot.py
A	examples/01_dos_1d_chain/plot.py
A	examples/01_dos_1d_chain/run_dos.json
M	examples/02_localization_1d_chain/README.md
D	examples/02_localization_1d_chain/fig_anderson.png
D	examples/02_localization_1d_chain/fig_localization.png
D	examples/02_localization_1d_chain/lsqloc.py
M	examples/02_localization_1d_chain/make_disordered_chain.py
A	examples/02_localization_1d_chain/plot.py
A	examples/02_localization_1d_chain/run_localization.json
M	examples/03_haldane_hall/README.md
D	examples/03_haldane_hall/fig_bands.png
D	examples/03_haldane_hall/fig_haldane_hall.png
D	examples/03_haldane_hall/fig_phase_diagram.png
D	examples/03_haldane_hall/haldane_kspace.py
D	examples/03_haldane_hall/lsqhall.py
M	examples/03_haldane_hall/make_haldane.py
A	examples/03_haldane_hall/plot.py
A	examples/03_haldane_hall/run_haldane.json
M	examples/04_spin_precession_magnetic_chain/README.md
D	examples/04_spin_precession_magnetic_chain/fig_precession.png
D	examples/04_spin_precession_magnetic_chain/fig_zeeman.png
D	examples/04_spin_precession_magnetic_chain/lsqspin.py
A	examples/04_spin_precession_magnetic_chain/make_chain_mag.py
A	examples/04_spin_precession_magnetic_chain/plot.py
A	examples/04_spin_precession_magnetic_chain/run_spin.json
M	examples/05_charge_transport_square/README.md
D	examples/05_charge_transport_square/fig_diffusion.png
D	examples/05_charge_transport_square/fig_fermi.png
D	examples/05_charge_transport_square/fig_transport.png
D	examples/05_charge_transport_square/lsqtransport.py
M	examples/05_charge_transport_square/make_square.py
A	examples/05_charge_transport_square/plot.py
A	examples/05_charge_transport_square/run_transport.json
M	examples/06_spin_relaxation_square/README.md
D	examples/06_spin_relaxation_square/ey_ratelaw.dat
D	examples/06_spin_relaxation_square/fig_dp.png
D	examples/06_spin_relaxation_square/fig_ey.png
D	examples/06_spin_relaxation_square/fig_ey_ratelaw.pdf
D	examples/06_spin_relaxation_square/fig_ey_ratelaw.png
D	examples/06_spin_relaxation_square/lsqrelax.py
A	examples/06_spin_relaxation_square/plot.py
D	examples/06_spin_relaxation_square/plot_ratelaw.py
D	examples/06_spin_relaxation_square/rp_spin.py
D	examples/06_spin_relaxation_square/rp_validate.py
A	examples/06_spin_relaxation_square/run_relaxation.json
A	examples/07_spin_hall_bi_bilayer/README.md
A	examples/07_spin_hall_bi_bilayer/plot.py
A	examples/07_spin_hall_bi_bilayer/run_dos_production.json
A	examples/07_spin_hall_bi_bilayer/run_dos_quick.json
A	examples/07_spin_hall_bi_bilayer/run_shc_production.json
A	examples/07_spin_hall_bi_bilayer/run_shc_quick.json
A	examples/08_spin_hall_pdse2_monolayer/README.md
A	examples/08_spin_hall_pdse2_monolayer/plot.py
A	examples/08_spin_hall_pdse2_monolayer/run_dos_production.json
A	examples/08_spin_hall_pdse2_monolayer/run_dos_quick.json
A	examples/08_spin_hall_pdse2_monolayer/run_shc_production.json
A	examples/08_spin_hall_pdse2_monolayer/run_shc_quick.json
M	python/lsquant/runfile.py
```

- **Tutorials 01 to 06:** Replaced `lsq*.py`, `rp_*.py`, and `plot_ratelaw.py` with `run_*.json` and `plot.py` using `examples/prl_style.py`; removed committed legacy `fig_*.png` files (the tester regenerates and commits the figures on TCN). Extracted `make_chain_mag.py` in `examples/04_spin_precession_magnetic_chain/` and updated `GENERATOR_CATALOGUE["make_chain_mag"]` in `python/lsquant/runfile.py`.
- **Examples 07 (`Bi(111)` bilayer) and 08 (`PdSe2` monolayer):** Created `README.md` (with production $\sigma^z_{xy}(E_F = 0)$ marked "value pending" pending Fable's diagnosis of the Aron difference), `run_dos_quick.json`, `run_shc_quick.json`, `run_dos_production.json`, `run_shc_production.json`, and 3-panel `plot.py` (DOS $M=512$ vs $M=2048$; KPM $\sigma^z_{xy}(E)$ vs exact eigenbasis oracle; residual $\sigma^z_{xy,\mathrm{KPM}} - \sigma^z_{xy,\mathrm{oracle}}$).

---

## 7. Verification Outputs on TCN (`~/lsquant-agi-build-0010`)

### 7.1 Basis-isolation check (`! grep -ri basis lsquant/core/chebyshev/`)

```text
=== 2. Basis-isolation grep check ===
PASS: ! grep -ri basis lsquant/core/chebyshev/ (0 matches)
```

### 7.2 Legacy-executable text check (`git grep` outside `handoff/`, `goldens/`, `AGENTS.md`)

```text
$ git grep -n -E 'inline_compute-kpm-|device-postProcess|MeanSquareDisplacement_DeviceSim|kuboBastinFromChebmom' -- ':!handoff' ':!goldens' ':!AGENTS.md'
(0 matches — exit code 1)
```

### 7.3 Developer self-checks (`stage_1_check` through `stage_7_check.py`) verbatim

```text
=== 4. Developer self-checks ===
dos max_abs_error: 2.47025e-15
kubo max_abs_error: 1.13687e-13
OMP_NUM_THREADS=1 vs 8 bit-identical: YES
--- stage_2_check.py ---
ReadCsrText: /home/ICN2/jgarcia/lsquant-agi-build-0010/goldens/stage0/chain1d/operators/chain1d.HAM.CSR anti_hermitian_norm=0 relative=0 legacy_fold=0
ReadCsrText: /home/ICN2/jgarcia/lsquant-agi-build-0010/goldens/stage0/chain1d/operators/chain1d.VX.CSR anti_hermitian_norm=0 relative=0 legacy_fold=0
(a) dos_plan array_equal vs G1:  True
(a) kubo_plan array_equal vs G2: True
(b) DeltaFunction.table vs legacy_delta_table max_abs_diff = 2.220446e-16
(c) dos_curve vs G3 differing lines: 0 (total lines: 1590)
--- stage_3_check.py ---
ReadCsrText: /home/ICN2/jgarcia/lsquant-agi-build-0010/goldens/stage0/chain1d/operators/chain1d.HAM.CSR anti_hermitian_norm=0 relative=0 legacy_fold=0
ReadCsrText: /home/ICN2/jgarcia/lsquant-agi-build-0010/goldens/stage0/haldane/operators/haldane.HAM.CSR anti_hermitian_norm=0 relative=0 legacy_fold=0
ReadCsrText: /home/ICN2/jgarcia/lsquant-agi-build-0010/goldens/stage0/haldane/operators/haldane.VX.CSR anti_hermitian_norm=0 relative=0 legacy_fold=0
ReadCsrText: /home/ICN2/jgarcia/lsquant-agi-build-0010/goldens/stage0/haldane/operators/haldane.VY.CSR anti_hermitian_norm=0 relative=0 legacy_fold=0
ReadCsrText: /home/ICN2/jgarcia/lsquant-agi-build-0010/goldens/stage0/graphene/operators/graphene_golden.HAM.CSR anti_hermitian_norm=0 relative=0 legacy_fold=0
ReadCsrText: /home/ICN2/jgarcia/lsquant-agi-build-0010/goldens/stage0/graphene/operators/graphene_golden.VX.CSR anti_hermitian_norm=0 relative=0 legacy_fold=0
(0) SparseMatrix.from_csr (scipy) dos array_equal vs ReadCsrText: True
(a) Haldane VX-VY array_equal vs G5: True
(a) Haldane VX-VX array_equal vs G5: True
(b) Graphene Gershgorin BandWidth = 15.0366847652085, array_equal vs G4: True
(c) GreenRetarded.table vs compat max_abs_diff           = 1.148853e-14
(c) GreenRetardedDerivative.table vs compat max_abs_diff = 2.702165e-13 (max_rel_diff = 7.173405e-15)
(d) G5 Bastin (10 meV)           differing lines: 0 (total lines: 3839)
(d) G5 Sea (10 meV)              differing lines: 0 (total lines: 3839)
(d) G5 Surf (10 meV)             differing lines: 0 (total lines: 3839)
(d) G5 Greenwood (10 meV)        differing lines: 0 (total lines: 127)
(d) G4 Bastin (64 meV)           differing lines: 0 (total lines: 1919)
(d) G4 Sea (10 meV)              differing lines: 0 (total lines: 1919)
(d) G4 Surf (10 meV)             differing lines: 0 (total lines: 1919)
(d) G4 Bastin FFTgrid (10 meV)   differing lines: 0 (total lines: 63)
(d) G4 Greenwood FFTgrid (10 meV) differing lines: 0 (total lines: 63)
--- stage_4_check.py ---
ReadCsrText: /home/ICN2/jgarcia/lsquant-agi-build-0010/goldens/stage0/chain1d_mag/operators/chain1d_mag.HAM.CSR anti_hermitian_norm=0 relative=0 legacy_fold=0
ReadCsrText: /home/ICN2/jgarcia/lsquant-agi-build-0010/goldens/stage0/chain1d_mag/operators/chain1d_mag.PX.CSR anti_hermitian_norm=0 relative=0 legacy_fold=0
ReadCsrText: /home/ICN2/jgarcia/lsquant-agi-build-0010/goldens/stage0/chain1d_mag/operators/chain1d_mag.SX.CSR anti_hermitian_norm=0 relative=0 legacy_fold=0
ReadCsrText: /home/ICN2/jgarcia/lsquant-agi-build-0010/goldens/stage0/chain1d_mag/operators/chain1d_mag.SY.CSR anti_hermitian_norm=0 relative=0 legacy_fold=0
ReadCsrText: /home/ICN2/jgarcia/lsquant-agi-build-0010/goldens/stage0/chain1d_mag/operators/chain1d_mag.SZ.CSR anti_hermitian_norm=0 relative=0 legacy_fold=0
ReadCsrText: /home/ICN2/jgarcia/lsquant-agi-build-0010/goldens/stage0/chain1d_ballistic/operators/chain1d.HAM.CSR anti_hermitian_norm=0 relative=0 legacy_fold=0
ReadCsrText: /home/ICN2/jgarcia/lsquant-agi-build-0010/goldens/stage0/chain1d_ballistic/operators/chain1d.VX.CSR anti_hermitian_norm=0 relative=0 legacy_fold=0
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
(d) time_correlation_plan chain1d VX-VX (M=64, steps=16, tmax=5, exact) array_equal vs reference SHA-256: True
--- stage_5_check.py ---
ReadCsrText: /tmp/lsquant_stage5_wbrcup6k/w2s_graphene_6x4/graphene.HAM.CSR anti_hermitian_norm=0 relative=0 legacy_fold=0
ReadCsrText: /tmp/lsquant_stage5_wbrcup6k/w2s_haldane_6x4/haldane.HAM.CSR anti_hermitian_norm=0 relative=0 legacy_fold=0
ReadCsrText: /tmp/lsquant_stage5_wbrcup6k/w2s_haldane_48x48/haldane.HAM.CSR anti_hermitian_norm=0 relative=0 legacy_fold=0
ReadCsrText: /tmp/lsquant_stage5_wbrcup6k/w2s_haldane_48x48/haldane.VX.CSR anti_hermitian_norm=0 relative=0 legacy_fold=0
ReadCsrText: /home/ICN2/jgarcia/lsquant-agi-build-0010/goldens/stage0/chain1d/operators/chain1d.HAM.CSR anti_hermitian_norm=0 relative=0 legacy_fold=0
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
ReadCsrText: /home/ICN2/jgarcia/lsquant-agi-build-0010/goldens/stage0/graphene/operators/graphene_golden.HAM.CSR anti_hermitian_norm=0 relative=0 legacy_fold=0
ReadCsrText: /home/ICN2/jgarcia/lsquant-agi-build-0010/goldens/stage0/graphene/operators/graphene_golden.VX.CSR anti_hermitian_norm=0 relative=0 legacy_fold=0
[check f] G4 graphene Kubo byte-identical to golden: PASS
ReadCsrText: /home/ICN2/jgarcia/lsquant-agi-build-0010/goldens/stage0/chain1d_mag/operators/chain1d_mag.HAM.CSR anti_hermitian_norm=0 relative=0 legacy_fold=0
ReadCsrText: /home/ICN2/jgarcia/lsquant-agi-build-0010/goldens/stage0/chain1d_mag/operators/chain1d_mag.PX.CSR anti_hermitian_norm=0 relative=0 legacy_fold=0
ReadCsrText: /home/ICN2/jgarcia/lsquant-agi-build-0010/goldens/stage0/chain1d_mag/operators/chain1d_mag.SX.CSR anti_hermitian_norm=0 relative=0 legacy_fold=0
[check f] G7 chain1d_mag spin precession byte-identical to golden: PASS
stage_5_check: ALL CHECKS PASSED
--- stage_6_check.py ---
read_csr_text: /home/ICN2/jgarcia/lsquant-agi-build-0010/goldens/stage0/chain1d/operators/chain1d.HAM.CSR anti_hermitian_norm=0 relative=0 legacy_fold=False
=== Stage 6 Self-Check: Rotated chain1d (N = 2048) ===
read_csr_text: /home/ICN2/jgarcia/lsquant-agi-build-0010/goldens/stage0/nonorth/chain1d/real/chain1d.HAM.CSR anti_hermitian_norm=0 relative=0 legacy_fold=False
read_csr_text: /home/ICN2/jgarcia/lsquant-agi-build-0010/goldens/stage0/nonorth/chain1d/real/chain1d.S.CSR anti_hermitian_norm=0 relative=0 legacy_fold=False
[1] real: max |eigvalsh(H', S') - eigvalsh(H)| = 1.154632e-14
[4] real: max ||S^-1_cheb v - S^-1_exact v|| / ||S^-1_exact v|| (M_S=60) = 1.452447e-12
[2] real route=lowdin      applier=exact     applier_res=2.686186e-14 max|mu - G1|=7.956630e-13 max|mu - delta_m0|=7.956630e-13
[2] real route=generalised applier=exact     applier_res=3.918055e-16 max|mu - G1|=4.664319e-16 max|mu - delta_m0|=4.664319e-16
[3] real route=lowdin      applier=chebyshev (M_S=12, eff_order=12) applier_res=7.094662e-05 max|mu - delta_m0|=5.308428e-04
[3] real route=generalised applier=chebyshev (M_S=12, eff_order=12) applier_res=2.580411e-04 max|mu - delta_m0|=1.859680e-03
[3] real route=lowdin      applier=chebyshev (M_S=60, eff_order=38) applier_res=1.781045e-12 max|mu - delta_m0|=9.506879e-12
[3] real route=generalised applier=chebyshev (M_S=60, eff_order=41) applier_res=1.687868e-12 max|mu - delta_m0|=1.517385e-12
[3] real route=generalised applier=ischo     (tol=1e-06) applier_res=5.071275e-07 max|mu - delta_m0|=1.675476e-08
[3] real route=generalised applier=ischo     (tol=1e-13) applier_res=1.345871e-16 max|mu - delta_m0|=2.737152e-15
read_csr_text: /home/ICN2/jgarcia/lsquant-agi-build-0010/goldens/stage0/nonorth/chain1d/complex/chain1d.HAM.CSR anti_hermitian_norm=0 relative=0 legacy_fold=False
read_csr_text: /home/ICN2/jgarcia/lsquant-agi-build-0010/goldens/stage0/nonorth/chain1d/complex/chain1d.S.CSR anti_hermitian_norm=0 relative=0 legacy_fold=False
[1] complex: max |eigvalsh(H', S') - eigvalsh(H)| = 7.993606e-15
[4] complex: max ||S^-1_cheb v - S^-1_exact v|| / ||S^-1_exact v|| (M_S=60) = 9.664000e-13
[2] complex route=lowdin      applier=exact     applier_res=3.004640e-14 max|mu - G1|=6.211824e-13 max|mu - delta_m0|=6.211824e-13
[2] complex route=generalised applier=exact     applier_res=4.162287e-16 max|mu - G1|=1.022115e-15 max|mu - delta_m0|=1.022115e-15
[3] complex route=lowdin      applier=chebyshev (M_S=60, eff_order=33) applier_res=1.782907e-12 max|mu - delta_m0|=9.336730e-13
[3] complex route=generalised applier=chebyshev (M_S=60, eff_order=36) applier_res=1.138146e-12 max|mu - delta_m0|=2.915651e-12
[3] complex route=generalised applier=ischo     (tol=1e-13) applier_res=1.491737e-16 max|mu - delta_m0|=2.860823e-15
=== Stage 6 Self-Check PASSED ===
--- stage_7_check.py ---
ReadCsrText: /home/ICN2/jgarcia/lsquant-agi-build-0010/goldens/stage0/nonorth/chain1d/real/chain1d.HAM.CSR anti_hermitian_norm=0 relative=0 legacy_fold=0
ReadCsrText: /home/ICN2/jgarcia/lsquant-agi-build-0010/goldens/stage0/nonorth/chain1d/real/chain1d.S.CSR anti_hermitian_norm=0 relative=0 legacy_fold=0
ReadCsrText: /home/ICN2/jgarcia/lsquant-agi-build-0010/goldens/stage0/nonorth/chain1d/real/chain1d.HAM.CSR anti_hermitian_norm=0 relative=0 legacy_fold=0
ReadCsrText: /home/ICN2/jgarcia/lsquant-agi-build-0010/goldens/stage0/nonorth/chain1d/real/chain1d.S.CSR anti_hermitian_norm=0 relative=0 legacy_fold=0
ReadCsrText: /home/ICN2/jgarcia/lsquant-agi-build-0010/goldens/stage0/nonorth/chain1d/real/chain1d.VX.CSR anti_hermitian_norm=0 relative=0 legacy_fold=0
== (a) Catalogue canonical parse/print round-trip ==
  PASS  Tr[delta(E-H)]                             -> shape=dos_identity         plan_builder=dos_plan
  PASS  Tr[A delta(E-H)]                           -> shape=dos_operator         plan_builder=dos_plan
  PASS  Tr[VX delta(E-H) VX delta(E-H)]            -> shape=kubo_greenwood       plan_builder=kubo_plan
  PASS  Tr[VX delta(E-H) VY dG-(E)/dE]             -> shape=kubo_bastin          plan_builder=kubo_plan
  PASS  Tr[[X,U(t)]^+ delta(E-H) [X,U(t)]]         -> shape=msd                  plan_builder=msd_plan
  PASS  Tr[P U(t)^+ SX U(t) delta(E-H) P]          -> shape=projected_evolution  plan_builder=projected_evolution_plan
  PASS  Tr[A delta(E-H) U(t)^+ B U(t)]             -> shape=time_correlation     plan_builder=time_correlation_plan
  PASS  Tr[sea[VX delta(E-H) VY dG-(E)/dE]]        -> shape=kubo_bastin          plan_builder=kubo_plan
  PASS  Tr[P:PX U(t)^+ SX U(t) delta(E-H) P:PX]    -> shape=projected_evolution  plan_builder=projected_evolution_plan
  PASS  Tr[B^+ VX@bloch B delta(E-H)]              -> shape=dos_operator         plan_builder=dos_plan

== (b) Refusal of 5 malformed / unsupported inputs ==
  PASS  [unknown_operator] refused with: Unknown operator 'FOO' in word 'Tr[FOO delta(E-H)]' (registered operators: ['1', 'VX', 'VY'])
  PASS  [rank_3] refused with: Unsupported word 'Tr[VX delta(E-H) VY delta(E-H) VZ delta(E-H)]' of rank 3 (3 expanded spectral slots): rank >= 3 moment tensors are out of scope; nearest recognised shape is 'Tr[VX delta(E-H) VX delta(E-H)]' (kubo_greenwood)
  PASS  [mixed_basis_without_transform] refused with: Mixed basis without explicit transform: operator 'VX@bloch' is in 'bloch' basis while calculation/active basis is 'real' (insert explicit B^+ ... B or U^+ ... U transform)
  PASS  [ut_without_time] refused with: Word 'Tr[A delta(E-H) U(t)^+ B U(t)]' contains time evolution U(t) without 'time' parameters (steps and tmax required)
  PASS  [unknown_key_in_runfile] refused with: Unknown key 'unknown_extra_key' in run file (top-level); allowed keys are ['basis', 'build_profile', 'expansion', 'metric', 'model', 'output', 'projection', 'receipt', 'reconstruction', 'reference', 'runtime', 'time', 'word', 'words']
  PASS  [unsupported_green_word] refused with: Unrecognised word shape 'Tr[VX delta(E-H) VY G+(E)]'; nearest recognised shape is 'Tr[VX delta(E-H) VY dG-(E)/dE]' (kubo_bastin)
  PASS  [eigen_window_on_csr] refused with: expansion.window='eigen' is not supported for model.kind='csr' (no k blocks; use gershgorin or an explicit window)
  PASS  [unmatched_formula_for_word] refused with: Reconstruction formula 'greenwood' in reconstruction.formulas does not match any word in ['Tr[VX delta(E-H) VY dG-(E)/dE]'] (allowed formulas for given word(s): ['bastin', 'bastin_fft', 'kubo_bastin', 'kubo_bastin_fft', 'kubo_bastin_sea', 'kubo_bastin_surf', 'sea', 'shc', 'surf'])
  PASS  [unknown_generator_param] refused with: Unknown key 'bogus_param' in run file (model.params (make_chain)); allowed keys are ['hopping_t', 'sizes']
  PASS  [unknown_generator] refused with: Unknown generator name 'evil_script' in run file (model.name); allowed generators are ['gen_models', 'make_chain', 'make_chain_mag', 'make_disordered_chain', 'make_ey', 'make_haldane', 'make_rashba', 'make_square']
  PASS  [window_default_by_basis] real->gershgorin, bloch->eigen, eigen->eigen
  PASS  [runfile_schema] all 64 files in runs/goldens/*.json validated

== (c) CLI reproduction of G1, G5, G9 (g9_graphene_eigen.json), and G10 (g10_chain1d_real_generalised.json) ==
  PASS  cmp IDENTICAL: G1 (SpectralOp1chain1dKPM_M64_stateexact.chebmom1D)
  PASS  cmp IDENTICAL: G5 (NonEqOpVX-VYhaldaneKPM_M128x128_stateexact.chebmom2D)
  PASS  cmp IDENTICAL: G5 (NonEqOpVX-VXhaldaneKPM_M128x128_stateexact.chebmom2D)
  PASS  cmp IDENTICAL: G5 (KuboBastin_NonEqOpVX-VYhaldaneKPM_M128x128_stateexactJACKSON.dat)
  PASS  cmp IDENTICAL: G5 (KuboBastinSea_NonEqOpVX-VYhaldaneKPM_M128x128_stateexactJACKSON.dat)
  PASS  cmp IDENTICAL: G5 (KuboBastinSurf_NonEqOpVX-VYhaldaneKPM_M128x128_stateexactJACKSON.dat)
  PASS  cmp IDENTICAL: G5 (KuboGreenWood_NonEqOpVX-VXhaldaneKPM_M128x128_stateexactJACKSON.dat)
  PASS  cmp IDENTICAL: G9 g9_graphene_eigen (KuboGreenWood_NonEqOpVX-VXgraphene_eigen_KPM_M256x256_stateexactJACKSON.dat)
  PASS  cmp IDENTICAL: G9 g9_graphene_eigen (NonEqOpVX-VXgraphene_eigen_KPM_M256x256_stateexact.basis.json)
  PASS  cmp IDENTICAL: G9 g9_graphene_eigen (NonEqOpVX-VXgraphene_eigen_KPM_M256x256_stateexact.chebmom)
  PASS  cmp IDENTICAL: G9 g9_graphene_eigen (NonEqOpVX-VXgraphene_eigen_KPM_M256x256_stateexact.chebmom2D)
  PASS  cmp IDENTICAL: G9 g9_graphene_eigen (SpectralOp1graphene_eigen_KPM_M256_stateexact.basis.json)
  PASS  cmp IDENTICAL: G9 g9_graphene_eigen (SpectralOp1graphene_eigen_KPM_M256_stateexact.chebmom)
  PASS  cmp IDENTICAL: G9 g9_graphene_eigen (SpectralOp1graphene_eigen_KPM_M256_stateexact.chebmom1D)
  PASS  cmp IDENTICAL: G9 g9_graphene_eigen (meanSpectralOp1graphene_eigen_KPM_M256_stateexactJACKSON.dat)
  PASS  cmp IDENTICAL: G10 g10_chain1d_real_generalised (KuboGreenWood_NonEqOpVX-VXchain1dKPM_M64x64_stateexactJACKSON.dat)
  PASS  cmp IDENTICAL: G10 g10_chain1d_real_generalised (NonEqOpVX-VXchain1dKPM_M64x64_stateexact.chebmom)
  PASS  cmp IDENTICAL: G10 g10_chain1d_real_generalised (NonEqOpVX-VXchain1dKPM_M64x64_stateexact.chebmom2D)
  PASS  cmp IDENTICAL: G10 g10_chain1d_real_generalised (SpectralOp1chain1dKPM_M64_stateexact.chebmom)
  PASS  cmp IDENTICAL: G10 g10_chain1d_real_generalised (SpectralOp1chain1dKPM_M64_stateexact.chebmom1D)
  PASS  cmp IDENTICAL: G10 g10_chain1d_real_generalised (meanSpectralOp1chain1dKPM_M64_stateexactJACKSON.dat)

== (d) MCP server subprocess vs CLI JSON equivalence ==
  MCP tools registered: ['words_describe', 'run', 'status', 'reconstruct', 'inspect']
  PASS  MCP words_describe JSON == CLI words describe JSON
  PASS  MCP inspect JSON == CLI inspect JSON
  PASS  MCP error responses preserve reason: 'Error executing tool words_describe: Unrecognised word shape 'Tr[VX delta(E-H) VY G+(E)]'; nearest recognised shape is 'Tr[VX delta(E-H) VY dG-(E)/dE]' (kubo_bastin)'

STAGE 7 SELF-CHECK PASS
```

### 7.4 Example 01 to 08 smoke runs on TCN (`python3 -m lsquant run` + `plot.py`)

All eight examples were executed on TCN into `/tmp/lsquant_examples_smoke/` (figures were not committed to the git repository).

```text
=== 01_dos_1d_chain ===
wrote /tmp/lsquant_examples_smoke/01_dos_1d_chain/run_dos/fig_levels.png
wrote /tmp/lsquant_examples_smoke/01_dos_1d_chain/run_dos/fig_moments.png
wrote /tmp/lsquant_examples_smoke/01_dos_1d_chain/run_dos/fig_dos.png
total 252
drwxr-xr-x  3 jgarcia Theoretical   4096 Oct 25 21:05 .
drwxr-xr-x 10 jgarcia Theoretical   4096 Oct 25 20:57 ..
-rw-r--r--  1 jgarcia Theoretical 159204 Oct 25 21:05 fig_dos.png
-rw-r--r--  1 jgarcia Theoretical  49459 Oct 25 21:05 fig_levels.png
-rw-r--r--  1 jgarcia Theoretical  30642 Oct 25 21:05 fig_moments.png
drwxr-xr-x  8 jgarcia Theoretical   4096 Oct 25 21:05 run_dos
DONE 01_dos_1d_chain

=== 02_localization_1d_chain ===
wrote /tmp/lsquant_examples_smoke/02_localization_1d_chain/run_localization/fig_anderson.png
wrote /tmp/lsquant_examples_smoke/02_localization_1d_chain/run_localization/fig_localization.png
total 192
drwxr-xr-x  3 jgarcia Theoretical   4096 Oct 25 21:06 .
drwxr-xr-x 10 jgarcia Theoretical   4096 Oct 25 20:57 ..
-rw-r--r--  1 jgarcia Theoretical  53051 Oct 25 21:06 fig_anderson.png
-rw-r--r--  1 jgarcia Theoretical 130033 Oct 25 21:06 fig_localization.png
drwxr-xr-x 19 jgarcia Theoretical   4096 Oct 25 21:06 run_localization
DONE 02_localization_1d_chain

=== 03_haldane_hall ===
ReadCsrText: /home/ICN2/jgarcia/lsquant-agi-build-0010/goldens/stage0/haldane/operators/haldane.HAM.CSR anti_hermitian_norm=0 relative=0 legacy_fold=0
ReadCsrText: /home/ICN2/jgarcia/lsquant-agi-build-0010/goldens/stage0/haldane/operators/haldane.VX.CSR anti_hermitian_norm=0 relative=0 legacy_fold=0
ReadCsrText: /home/ICN2/jgarcia/lsquant-agi-build-0010/goldens/stage0/haldane/operators/haldane.VY.CSR anti_hermitian_norm=0 relative=0 legacy_fold=0
/tmp/lsquant_examples_smoke/03_haldane_hall/run_haldane/NonEqOpVX-VYhaldaneKPM_M128x128_stateexact.chebmom2D
/tmp/lsquant_examples_smoke/03_haldane_hall/run_haldane/KuboBastin_NonEqOpVX-VYhaldaneKPM_M128x128_stateexactJACKSON.dat
/tmp/lsquant_examples_smoke/03_haldane_hall/run_haldane/KuboBastinSea_NonEqOpVX-VYhaldaneKPM_M128x128_stateexactJACKSON.dat
/tmp/lsquant_examples_smoke/03_haldane_hall/run_haldane/KuboBastinSurf_NonEqOpVX-VYhaldaneKPM_M128x128_stateexactJACKSON.dat
/tmp/lsquant_examples_smoke/03_haldane_hall/run_haldane/NonEqOpVX-VXhaldaneKPM_M128x128_stateexact.chebmom2D
/tmp/lsquant_examples_smoke/03_haldane_hall/run_haldane/KuboGreenWood_NonEqOpVX-VXhaldaneKPM_M128x128_stateexactJACKSON.dat
wrote /tmp/lsquant_examples_smoke/03_haldane_hall/run_haldane/fig_bands.png
wrote /tmp/lsquant_examples_smoke/03_haldane_hall/run_haldane/fig_haldane_hall.png
Chern from gap plateau: C = +1.017
wrote /tmp/lsquant_examples_smoke/03_haldane_hall/run_haldane/fig_phase_diagram.png
total 192
drwxr-xr-x  3 jgarcia Theoretical   4096 Oct 25 21:06 .
drwxr-xr-x 10 jgarcia Theoretical   4096 Oct 25 20:57 ..
-rw-r--r--  1 jgarcia Theoretical  37845 Oct 25 21:06 fig_bands.png
-rw-r--r--  1 jgarcia Theoretical 101122 Oct 25 21:06 fig_haldane_hall.png
-rw-r--r--  1 jgarcia Theoretical  40283 Oct 25 21:06 fig_phase_diagram.png
drwxr-xr-x  2 jgarcia Theoretical   4096 Oct 25 21:06 run_haldane
DONE 03_haldane_hall

=== 04_spin_precession_magnetic_chain ===
ReadCsrText: /home/ICN2/jgarcia/lsquant-agi-build-0010/goldens/stage0/chain1d_mag/operators/chain1d_mag.HAM.CSR anti_hermitian_norm=0 relative=0 legacy_fold=0
ReadCsrText: /home/ICN2/jgarcia/lsquant-agi-build-0010/goldens/stage0/chain1d_mag/operators/chain1d_mag.SX.CSR anti_hermitian_norm=0 relative=0 legacy_fold=0
ReadCsrText: /home/ICN2/jgarcia/lsquant-agi-build-0010/goldens/stage0/chain1d_mag/operators/chain1d_mag.SY.CSR anti_hermitian_norm=0 relative=0 legacy_fold=0
ReadCsrText: /home/ICN2/jgarcia/lsquant-agi-build-0010/goldens/stage0/chain1d_mag/operators/chain1d_mag.SZ.CSR anti_hermitian_norm=0 relative=0 legacy_fold=0
ReadCsrText: /home/ICN2/jgarcia/lsquant-agi-build-0010/goldens/stage0/chain1d_mag/operators/chain1d_mag.PX.CSR anti_hermitian_norm=0 relative=0 legacy_fold=0
/tmp/lsquant_examples_smoke/04_spin_precession_magnetic_chain/run_spin/TimeProjSX-PXchain1d_magKPM_M32_stateexact.chebmomTD
/tmp/lsquant_examples_smoke/04_spin_precession_magnetic_chain/run_spin/TimeProjSY-PXchain1d_magKPM_M32_stateexact.chebmomTD
/tmp/lsquant_examples_smoke/04_spin_precession_magnetic_chain/run_spin/TimeProjSZ-PXchain1d_magKPM_M32_stateexact.chebmomTD
wrote /tmp/lsquant_examples_smoke/04_spin_precession_magnetic_chain/run_spin/fig_zeeman.png
wrote /tmp/lsquant_examples_smoke/04_spin_precession_magnetic_chain/run_spin/fig_precession.png
total 164
drwxr-xr-x  3 jgarcia Theoretical  4096 Oct 25 21:05 .
drwxr-xr-x 10 jgarcia Theoretical  4096 Oct 25 20:57 ..
-rw-r--r--  1 jgarcia Theoretical 90901 Oct 25 21:05 fig_precession.png
-rw-r--r--  1 jgarcia Theoretical 58507 Oct 25 21:05 fig_zeeman.png
drwxr-xr-x  2 jgarcia Theoretical  4096 Oct 25 21:05 run_spin
DONE 04_spin_precession_magnetic_chain

=== 05_charge_transport_square ===
/tmp/lsquant_examples_smoke/05_charge_transport_square/run_transport/meanCorrelationVXsquare_L200_W0.5KPM_M256_state150defaultEF2.000000JACKSON.dat
/tmp/lsquant_examples_smoke/05_charge_transport_square/run_transport/meanSpectralOp1square_L200_W0.5KPM_M1024_statedefaultJACKSON.dat
/tmp/lsquant_examples_smoke/05_charge_transport_square/run_transport/meanCorrelationVXsquare_L200_W1KPM_M256_state150defaultEF2.000000JACKSON.dat
/tmp/lsquant_examples_smoke/05_charge_transport_square/run_transport/meanSpectralOp1square_L200_W1KPM_M1024_statedefaultJACKSON.dat
/tmp/lsquant_examples_smoke/05_charge_transport_square/run_transport/meanCorrelationVXsquare_L200_W2KPM_M256_state150defaultEF2.000000JACKSON.dat
/tmp/lsquant_examples_smoke/05_charge_transport_square/run_transport/meanSpectralOp1square_L200_W2KPM_M1024_statedefaultJACKSON.dat
/tmp/lsquant_examples_smoke/05_charge_transport_square/run_transport/meanCorrelationVXsquare_L200_W4KPM_M256_state150defaultEF2.000000JACKSON.dat
/tmp/lsquant_examples_smoke/05_charge_transport_square/run_transport/meanSpectralOp1square_L200_W4KPM_M1024_statedefaultJACKSON.dat
wrote /tmp/lsquant_examples_smoke/05_charge_transport_square/run_transport/fig_fermi.png
W=0.5  D_max=21.874  rho(E)=0.1138  sigma= 2.489  ell= 19.66
W=1.0  D_max= 6.855  rho(E)=0.1123  sigma= 0.770  ell=  6.16
W=2.0  D_max= 1.711  rho(E)=0.1079  sigma= 0.185  ell=  1.54
W=4.0  D_max= 0.392  rho(E)=0.1115  sigma= 0.044  ell=  0.35
wrote /tmp/lsquant_examples_smoke/05_charge_transport_square/run_transport/fig_diffusion.png
wrote /tmp/lsquant_examples_smoke/05_charge_transport_square/run_transport/fig_transport.png
total 260
drwxr-xr-x  3 jgarcia Theoretical   4096 Oct 25 21:08 .
drwxr-xr-x 10 jgarcia Theoretical   4096 Oct 25 20:57 ..
-rw-r--r--  1 jgarcia Theoretical 107700 Oct 25 21:08 fig_diffusion.png
-rw-r--r--  1 jgarcia Theoretical  85871 Oct 25 21:07 fig_fermi.png
-rw-r--r--  1 jgarcia Theoretical  53378 Oct 25 21:08 fig_transport.png
drwxr-xr-x  3 jgarcia Theoretical   4096 Oct 25 21:08 run_transport
DONE 05_charge_transport_square

=== 06_spin_relaxation_square ===
EY  Dsf=0.05  T1=51.3 fs  1/T1=0.0195
EY  Dsf=0.07  T1=38.8 fs  1/T1=0.0258
EY  Dsf=0.10  T1=32.5 fs  1/T1=0.0308
wrote /tmp/lsquant_examples_smoke/06_fast/run_relaxation/fig_ey.png
wrote /tmp/lsquant_examples_smoke/06_fast/run_relaxation/fig_dp.png
wrote /tmp/lsquant_examples_smoke/06_fast/run_relaxation/fig_ey_ratelaw.png
wrote fig_ey_ratelaw.png (fitted p = 1.26; theory p = 2)
total 252
drwxr-xr-x  3 jgarcia Theoretical  4096 Oct 25 22:07 .
drwxr-xr-x 13 jgarcia Theoretical  4096 Oct 25 21:33 ..
-rw-r--r--  1 jgarcia Theoretical   273 Oct 25 22:07 ey_ratelaw.dat
-rw-r--r--  1 jgarcia Theoretical 90266 Oct 25 22:07 fig_dp.png
-rw-r--r--  1 jgarcia Theoretical 77342 Oct 25 22:07 fig_ey.png
-rw-r--r--  1 jgarcia Theoretical 66253 Oct 25 22:07 fig_ey_ratelaw.png
drwxr-xr-x  2 jgarcia Theoretical  4096 Oct 25 22:07 run_relaxation
DONE 06_spin_relaxation_square

=== 07_spin_hall_bi_bilayer (run_dos_quick.json + run_shc_quick.json + plot.py + sha256sum) ===
stage6b: dipole quadrature started {'enabled': False}
stage6b: dipole quadrature completed wall_seconds=0.000433 {'enabled': False}
stage6b: metric build started {}
stage6b: metric build completed wall_seconds=14.269572 {}
stage6b: operators started {}
stage6b: operators completed wall_seconds=0.000800 {}
stage6b: moments started {'word': 'Tr[delta(E-H)]'}
stage6b: moments completed wall_seconds=17.069579 {'word': 'Tr[delta(E-H)]'}
stage6b: reconstruction started {'word': 'Tr[delta(E-H)]'}
stage6b: reconstruction completed wall_seconds=9.681301 {'word': 'Tr[delta(E-H)]'}
stage6b: eigen reference skipped: production already computes the exact eigenbasis trace {'word': 'Tr[delta(E-H)]'}
stage6b: oracle started {}
stage6b: oracle completed wall_seconds=18.737633 {}
/tmp/lsquant_examples_smoke/07_spin_hall_bi_bilayer/run_dos_quick/SpectralOp1biKPM_M256_stateexact.chebmom
/tmp/lsquant_examples_smoke/07_spin_hall_bi_bilayer/run_dos_quick/SpectralOp1biKPM_M256_stateexact.chebmom1D
/tmp/lsquant_examples_smoke/07_spin_hall_bi_bilayer/run_dos_quick/basis.json
/tmp/lsquant_examples_smoke/07_spin_hall_bi_bilayer/run_dos_quick/dos.dat
/tmp/lsquant_examples_smoke/07_spin_hall_bi_bilayer/run_dos_quick/eigen_comparison.json
/tmp/lsquant_examples_smoke/07_spin_hall_bi_bilayer/run_dos_quick/eigen_hermiticity.json
/tmp/lsquant_examples_smoke/07_spin_hall_bi_bilayer/run_dos_quick/exact_shc.json
/tmp/lsquant_examples_smoke/07_spin_hall_bi_bilayer/run_dos_quick/run.resolved.json
/tmp/lsquant_examples_smoke/07_spin_hall_bi_bilayer/run_dos_quick/spin_chern.json
stage6b: dipole quadrature started {'enabled': False}
stage6b: dipole quadrature completed wall_seconds=0.000280 {'enabled': False}
stage6b: metric build started {}
stage6b: metric build completed wall_seconds=13.765673 {}
stage6b: operators started {}
eigen JXSZ: Hermitian fold removed anti-Hermitian norm=3.0160545172943591e-12 operator_scale=1874.5901651207578 relative_norm=1.6089140834151715e-15
eigen VY: Hermitian fold removed anti-Hermitian norm=5.301513538420041e-12 operator_scale=1641.0756141631316 relative_norm=3.2305114357107513e-15
stage6b: operators completed wall_seconds=5.863416 {}
stage6b: moments started {'word': 'Tr[JXSZ delta(E-H) VY dG-(E)/dE]'}
stage6b: moments completed wall_seconds=170.727192 {'word': 'Tr[JXSZ delta(E-H) VY dG-(E)/dE]'}
stage6b: reconstruction started {'word': 'Tr[JXSZ delta(E-H) VY dG-(E)/dE]'}
stage6b: reconstruction completed wall_seconds=15.774334 {'word': 'Tr[JXSZ delta(E-H) VY dG-(E)/dE]'}
stage6b: eigen reference skipped: production already computes the exact eigenbasis trace {'word': 'Tr[JXSZ delta(E-H) VY dG-(E)/dE]'}
stage6b: oracle started {}
stage6b: oracle completed wall_seconds=9.644210 {}
/tmp/lsquant_examples_smoke/07_spin_hall_bi_bilayer/run_shc_quick/NonEqOpJXSZ-VYbiKPM_M256x256_stateexact.chebmom
/tmp/lsquant_examples_smoke/07_spin_hall_bi_bilayer/run_shc_quick/NonEqOpJXSZ-VYbiKPM_M256x256_stateexact.chebmom2D
/tmp/lsquant_examples_smoke/07_spin_hall_bi_bilayer/run_shc_quick/basis.json
/tmp/lsquant_examples_smoke/07_spin_hall_bi_bilayer/run_shc_quick/eigen_comparison.json
/tmp/lsquant_examples_smoke/07_spin_hall_bi_bilayer/run_shc_quick/eigen_hermiticity.json
/tmp/lsquant_examples_smoke/07_spin_hall_bi_bilayer/run_shc_quick/exact_shc.json
/tmp/lsquant_examples_smoke/07_spin_hall_bi_bilayer/run_shc_quick/run.resolved.json
/tmp/lsquant_examples_smoke/07_spin_hall_bi_bilayer/run_shc_quick/shc.dat
/tmp/lsquant_examples_smoke/07_spin_hall_bi_bilayer/run_shc_quick/spin_chern.json
wrote /tmp/lsquant_examples_smoke/07_spin_hall_bi_bilayer/fig_dos_shc.png
total 136
drwxr-xr-x  4 jgarcia Theoretical   4096 Oct 25 21:00 .
drwxr-xr-x 10 jgarcia Theoretical   4096 Oct 25 20:57 ..
-rw-r--r--  1 jgarcia Theoretical 121563 Oct 25 21:05 fig_dos_shc.png
drwxr-xr-x  2 jgarcia Theoretical   4096 Oct 25 20:54 run_dos_quick
drwxr-xr-x  2 jgarcia Theoretical   4096 Oct 25 20:59 run_shc_quick
examples/07_spin_hall_bi_bilayer/siesta/Bi.ion.xml: OK
examples/07_spin_hall_bi_bilayer/siesta/Bi.psml: OK
examples/07_spin_hall_bi_bilayer/siesta/bi.TSHS: OK
examples/07_spin_hall_bi_bilayer/siesta/bi.fdf: OK
DONE 07_spin_hall_bi_bilayer

=== 08_spin_hall_pdse2_monolayer (run_dos_quick.json + run_shc_quick.json + plot.py + sha256sum) ===
stage6b: dipole quadrature started {'enabled': False}
stage6b: dipole quadrature completed wall_seconds=0.002004 {'enabled': False}
stage6b: metric build started {}
stage6b: metric build completed wall_seconds=161.687203 {}
stage6b: operators started {}
stage6b: operators completed wall_seconds=0.000819 {}
stage6b: moments started {'word': 'Tr[delta(E-H)]'}
stage6b: moments completed wall_seconds=63.331709 {'word': 'Tr[delta(E-H)]'}
stage6b: reconstruction started {'word': 'Tr[delta(E-H)]'}
stage6b: reconstruction completed wall_seconds=8.161716 {'word': 'Tr[delta(E-H)]'}
stage6b: eigen reference skipped: production already computes the exact eigenbasis trace {'word': 'Tr[delta(E-H)]'}
stage6b: oracle started {}
stage6b: oracle completed wall_seconds=257.026443 {}
/tmp/lsquant_examples_smoke/08_spin_hall_pdse2_monolayer/run_dos_quick/SpectralOp1pdse2KPM_M256_stateexact.chebmom
/tmp/lsquant_examples_smoke/08_spin_hall_pdse2_monolayer/run_dos_quick/SpectralOp1pdse2KPM_M256_stateexact.chebmom1D
/tmp/lsquant_examples_smoke/08_spin_hall_pdse2_monolayer/run_dos_quick/basis.json
/tmp/lsquant_examples_smoke/08_spin_hall_pdse2_monolayer/run_dos_quick/dos.dat
/tmp/lsquant_examples_smoke/08_spin_hall_pdse2_monolayer/run_dos_quick/eigen_comparison.json
/tmp/lsquant_examples_smoke/08_spin_hall_pdse2_monolayer/run_dos_quick/eigen_hermiticity.json
/tmp/lsquant_examples_smoke/08_spin_hall_pdse2_monolayer/run_dos_quick/exact_shc.json
/tmp/lsquant_examples_smoke/08_spin_hall_pdse2_monolayer/run_dos_quick/run.resolved.json
/tmp/lsquant_examples_smoke/08_spin_hall_pdse2_monolayer/run_dos_quick/spin_chern.json
stage6b: dipole quadrature started {'enabled': False}
stage6b: dipole quadrature completed wall_seconds=0.000706 {'enabled': False}
stage6b: metric build started {}
stage6b: metric build completed wall_seconds=123.152524 {}
stage6b: operators started {}
eigen JXSZ: Hermitian fold removed anti-Hermitian norm=4.8478798952421658e-12 operator_scale=1965.050829809747 relative_norm=2.4670506338563923e-15
eigen VY: Hermitian fold removed anti-Hermitian norm=6.1228446661570681e-12 operator_scale=2097.3250013136594 relative_norm=2.9193590227179977e-15
stage6b: operators completed wall_seconds=123.652831 {}
stage6b: moments started {'word': 'Tr[JXSZ delta(E-H) VY dG-(E)/dE]'}
stage6b: moments completed wall_seconds=1008.983832 {'word': 'Tr[JXSZ delta(E-H) VY dG-(E)/dE]'}
stage6b: reconstruction started {'word': 'Tr[JXSZ delta(E-H) VY dG-(E)/dE]'}
stage6b: reconstruction completed wall_seconds=1.085918 {'word': 'Tr[JXSZ delta(E-H) VY dG-(E)/dE]'}
stage6b: eigen reference skipped: production already computes the exact eigenbasis trace {'word': 'Tr[JXSZ delta(E-H) VY dG-(E)/dE]'}
stage6b: oracle started {}
stage6b: oracle completed wall_seconds=12.506428 {}
/tmp/lsquant_examples_smoke/08_spin_hall_pdse2_monolayer/run_shc_quick/NonEqOpJXSZ-VYpdse2KPM_M256x256_stateexact.chebmom
/tmp/lsquant_examples_smoke/08_spin_hall_pdse2_monolayer/run_shc_quick/NonEqOpJXSZ-VYpdse2KPM_M256x256_stateexact.chebmom2D
/tmp/lsquant_examples_smoke/08_spin_hall_pdse2_monolayer/run_shc_quick/basis.json
/tmp/lsquant_examples_smoke/08_spin_hall_pdse2_monolayer/run_shc_quick/eigen_comparison.json
/tmp/lsquant_examples_smoke/08_spin_hall_pdse2_monolayer/run_shc_quick/eigen_hermiticity.json
/tmp/lsquant_examples_smoke/08_spin_hall_pdse2_monolayer/run_shc_quick/exact_shc.json
/tmp/lsquant_examples_smoke/08_spin_hall_pdse2_monolayer/run_shc_quick/run.resolved.json
/tmp/lsquant_examples_smoke/08_spin_hall_pdse2_monolayer/run_shc_quick/shc.dat
/tmp/lsquant_examples_smoke/08_spin_hall_pdse2_monolayer/run_shc_quick/spin_chern.json
wrote /tmp/lsquant_examples_smoke/08_spin_hall_pdse2_monolayer/fig_dos_shc.png
total 136
drwxr-xr-x  4 jgarcia Theoretical   4096 Oct 25 21:23 .
drwxr-xr-x 11 jgarcia Theoretical   4096 Oct 25 21:10 ..
-rw-r--r--  1 jgarcia Theoretical 119474 Oct 25 21:23 fig_dos_shc.png
drwxr-xr-x  2 jgarcia Theoretical   4096 Oct 25 21:02 run_dos_quick
drwxr-xr-x  2 jgarcia Theoretical   4096 Oct 25 21:23 run_shc_quick
examples/08_spin_hall_pdse2_monolayer/siesta/Pd.ion.xml: OK
examples/08_spin_hall_pdse2_monolayer/siesta/Pd.psml: OK
examples/08_spin_hall_pdse2_monolayer/siesta/Se.ion.xml: OK
examples/08_spin_hall_pdse2_monolayer/siesta/Se.psml: OK
examples/08_spin_hall_pdse2_monolayer/siesta/pdse2.TSHS: OK
examples/08_spin_hall_pdse2_monolayer/siesta/pdse2.fdf: OK
DONE 08_spin_hall_pdse2_monolayer
```

---

## 8. "Needs Harness Change" List (For Tester / Fable)

1. **`AGENTS.md` (line 353 and `scripts/check_goldens.sh` stage 7 text check):**
   - `AGENTS.md` line 353 references `lsquant_chain_reference.py`, which was moved to `python/lsquant/oracle/chain.py` (`python -m lsquant.oracle.chain`).
   - `AGENTS.md` is the only remaining file outside `handoff/` and `goldens/` that mentions `inline_compute-kpm-*`, `device-postProcess`, `MeanSquareDisplacement_DeviceSim`, and `kuboBastinFromChebmom`. Per Section 3 of the brief, the tester adds `:!AGENTS.md` to the `inline_*` record exclusion in `scripts/check_goldens.sh`.
2. **`figures/stage_0.py`:**
   - References the deleted legacy tutorial driver filenames (`lsqplot.py`, `lsqloc.py`, `lsqhall.py`, `lsqspin.py`, `lsqtransport.py`, `lsqrelax.py`, `plot_ratelaw.py`). Per Section 1 of the brief (`figures/` is not edited in this handoff), updating or retiring `figures/stage_0.py` requires a harness change.
3. **`scripts/check_goldens.sh` default `BUILD` directory:**
   - `scripts/check_goldens.sh` defaults `BUILD` to `$ROOT/build` (expecting `$ROOT/build/lsquant_core`), whereas `CMakePresets.json` places preset builds in `build/golden` and `build/native`. Updating `scripts/check_goldens.sh` to fall back to `$ROOT/build/golden` when `$ROOT/build/lsquant_core` is absent is a 3-line harness change.
4. **`runs/goldens/g11_bi_*.json`, `runs/goldens/g11_pdse2_*.json`, and `scripts/check_goldens.sh` stages `6b`/`6c`:**
   - The `g11_*` golden run files hard-code `"dir": "~/from_fatnode/shc_bi/run"` and `"~/from_fatnode/shc_pdse2/run"`, and `scripts/check_goldens.sh` looks for `~/shc_bi/run` or `~/from_fatnode/shc_bi/run`. They can now point to the committed `examples/07_spin_hall_bi_bilayer/siesta` and `examples/08_spin_hall_pdse2_monolayer/siesta` directories.
5. **Golden reproduction entry-point names (`stage1_moments`, `stage6_nonorth`, `lsquant.apps.stage*`, `lsquant.apps.tutorial0*`):**
   - Kept in this handoff as required by Section 1 because `scripts/check_goldens.sh` and `figures/stage_*.py` invoke them by name; renaming them is deferred to a follow-up harness change.

---

## 9. Deviations and Open Questions

- **Deviations from the brief:**
  1. **`python/lsquant/projection/diagnostics.py` deletion (now reverted):** Deleted in `a8252cc` as an unimported module; reverted in Review 1 Item 1 (`fae21b2`) and exported in `python/lsquant/projection/__init__.py` because `python/lsquant/projection/` is 0006d work still under validation.
  2. **Quick-run settings in Examples 07 and 08 (now fixed):** In `360339a`, `examples/07_spin_hall_bi_bilayer/run_*_quick.json` and `examples/08_spin_hall_pdse2_monolayer/run_*_quick.json` used `"velocity_mode": "berry_connection"` (and `12x12x1` for Bi instead of `24x24x1`). Fixed in Review 1 Item 3 (`cbd7dc4`) to use `"velocity_mode": "covariant"` with `24x24x1` ($M = 256$) for Bi and `12x12x1` ($M = 256$) for PdSe2.
  3. **`fftw3_omp` link (now removed):** In `26313dc`, `CMakeLists.txt` looked up `FFTW3_OMP_LIBRARY` and linked it into `lsquant_core` even though `fftw_init_threads` is never called. Removed in Review 1 Item 5 (`79cb251`).
  4. **`-Wpedantic` warning flag:** Added `-Wpedantic` alongside `-Wall -Wextra` on the `lsquant_warnings` interface target in `CMakeLists.txt` (`26313dc`); retained per Review 1 Item 5 because both `golden` and `native` presets build with 0 compiler warnings.
  5. **`.githooks/pre-commit` and `scripts/check_reporting.sh`:** Deleted `scripts/check_reporting.sh` (which checked the removed legacy `src/` and `include/` directories for the abandoned reporting plan) and updated `.githooks/pre-commit` in `da6239f` to guard its invocation with `if [ -x scripts/check_reporting.sh ]`.
  6. **Deletion of `.projectile` and `.vscode/settings.json`:** Deleted both tracked editor configuration files (which referenced legacy header and CMake paths) in `da6239f` as part of the Step 2 whole-tree search for editor files in `git ls-files`.
- **Open questions:** None for Handoff 0010. As stated in Section 0 and Section 8 of the brief, the 0009 wrap-up items W02 to W12 (including the numerical difference on Aron for Bi and PdSe2) remain open for Fable's separate diagnosis and maintainer decision.
