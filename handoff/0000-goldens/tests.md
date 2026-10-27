# Handoff 0000 — tests (Codex): stage-0 golden masters from the legacy build

Everything below ran on TCN in `~/lsquant-unified`, branch `handoff/0000-codex`, cut from
`cpp-unified` at `8ecbf836677f8d7a13f72dcdd503225a6ffac820`. Date 2026-10-20 (UTC).
No file under `src/`, `include/`, `test/`, `examples/` or `lsquant/` was modified; the
`test/golden` regeneration was reverted (section 3).

## 1. Machine and versions

| | |
|---|---|
| Host | `fatnode-tacn`, AMD EPYC 7551, 64 logical CPUs, 2 TB RAM |
| OS | Ubuntu 22.04.4 LTS, kernel 5.15.0-161-generic |
| Compiler | g++ (Ubuntu 11.4.0-1ubuntu1~22.04.3) 11.4.0 |
| CMake | 3.22.1 (CMAKE_BUILD_TYPE=Release, default flags) |
| Python | 3.10.12; numpy 2.2.6, scipy 1.15.3, matplotlib 3.10.9 (MPLBACKEND=Agg) |
| Libraries | Eigen 3.4.0, FFTW 3.3.8, GSL 2.7.1, OpenMP 4.5 |
| wannier2sparse | `26cab4affff8266f8a76ec73f3201e6c063a1ba6` (`26cab4a`, confirmed) |

## 2. Legacy build and ctest

```
cmake -S . -B build -DLSQUANT_BUILD_NEW_CORE=OFF -DLSQUANT_BUILD_INPUT_TOOLS=ON
cmake --build build -j 32
OMP_NUM_THREADS=16 KPM_SEED=12345 ctest --test-dir build --output-on-failure
```

Configure exit 0, build exit 0, ctest exit 0. The full ctest summary:

```
 1/31 Test  #1: chain1d_analytic .................   Passed    0.01 sec
 2/31 Test  #2: chain1d_dos ......................   Passed    0.00 sec
 3/31 Test  #3: dos_reconstructed ................   Passed    0.01 sec
 4/31 Test  #4: chain1d_reproducible .............   Passed   22.83 sec
 5/31 Test  #5: spin_algebra .....................   Passed    0.02 sec
 6/31 Test  #6: graphene_regression ..............   Passed    0.16 sec
 7/31 Test  #7: reconstruction_nan_guard .........   Passed    0.27 sec
 8/31 Test  #8: msd_ballistic ....................   Passed   19.66 sec
 9/31 Test  #9: hall_haldane .....................   Passed    9.42 sec
10/31 Test #10: spin_precession ..................   Passed    8.02 sec
11/31 Test #11: recon_kernel_oracle ..............   Passed    0.21 sec
12/31 Test #12: spmv_threaded_bitexact ...........   Passed    0.01 sec
13/31 Test #13: spmm_batchmultiply_equiv .........   Passed    0.01 sec
14/31 Test #14: greenwood_regression .............   Passed    0.06 sec
15/31 Test #15: safety_factors_guard .............   Passed    0.53 sec
16/31 Test #16: descriptor_conformance ...........   Passed    0.00 sec
17/31 Test #17: descriptor_bounds ................   Passed    0.13 sec
18/31 Test #18: tb_descriptor ....................   Passed    0.00 sec
19/31 Test #19: config_equivalence ...............   Passed    0.27 sec
20/31 Test #20: corrop_state_honored .............   Passed    2.26 sec
21/31 Test #21: fft_grid_regression ..............   Passed    0.07 sec
22/31 Test #22: sea_surf_regression ..............   Passed    0.62 sec
23/31 Test #23: jackson_kernel ...................   Passed    0.00 sec
24/31 Test #24: jackson_g_scalar .................   Passed    0.00 sec
25/31 Test #25: phase6_wired .....................   Passed    0.00 sec
26/31 Test #26: single_binary ....................   Passed    1.06 sec
27/31 Test #27: reporting_smoke ..................   Passed    0.16 sec
28/31 Test #28: unified_io_smoke .................   Passed    0.25 sec
29/31 Test #29: example_dos_smoke ................   Passed    0.09 sec
30/31 Test #30: example_localization_smoke .......   Passed    0.31 sec
31/31 Test #31: example_haldane_smoke ............   Passed   17.53 sec

100% tests passed, 0 tests failed out of 31

Total Test time (real) =  84.00 sec
```

Configure-time remarks (not failures, recorded as found):

- `CMake Warning: Manually-specified variables were not used by the project:
  LSQUANT_BUILD_NEW_CORE`. The option does not exist yet at `8ecbf83`, so `=OFF` is a no-op;
  the build is the pure legacy tree.
- `CMake Warning at cmake/FetchWannier2Sparse.cmake:9: wannier2sparse unavailable
  (offline?)`. Not offline: the fetch ran (`* branch 26cab4a... -> FETCH_HEAD` in the log)
  and built `.cache/wannier2sparse/build/wannier2sparse`. The `execute_process` captures
  the stdout of the nested `cmake` configure/build into `W2S_BIN`, so `EXISTS "${W2S_BIN}"`
  is false on a cold cache. Running `bash scripts/fetch_wannier2sparse.sh` afterwards
  (cache hit) prints only the path,
  `/home/ICN2/jgarcia/lsquant-unified/.cache/wannier2sparse/build/wannier2sparse`, at sha
  `26cab4affff8266f8a76ec73f3201e6c063a1ba6`. No test depends on `LSQUANT_W2S_BIN`.
- The build log has 101 compiler warnings in the legacy tree (32 "inline variables are only
  available with -std=c++17", 62 `-Wnarrowing` `size_t -> int` of `m0`/`m1`, 7
  unused `fread` results). Legacy code, not touched.

## 3. Golden regeneration (`test/golden`)

```
KPM_SEED=12345 OMP_NUM_THREADS=16 bash scripts/regenerate_goldens.sh build
git status --short test/golden GOLDEN_PROVENANCE.txt
```

Exit 0 (37 s). **`test/golden` changed**:

```
 M GOLDEN_PROVENANCE.txt
 M test/golden/KuboBastin_NonEqOpVX-VXgraphene_goldenKPM_M64x64_statedefaultJACKSON.dat
 M test/golden/NonEqOpVX-VXgraphene_goldenKPM_M64x64_statedefault.chebmom2D
 D test/golden/operators/graphene_golden.HAM.desc
 D test/golden/operators/graphene_golden.VX.desc
```

The whole chain1d part (operators, `BOUNDS`, G1 `.chebmom1D`, G2 `.chebmom2D`,
`DOS_reconstructed.dat`) and every graphene `.CSR` came out byte-identical. The changes, in
order of the files:

1. `GOLDEN_PROVENANCE.txt`: only `lsquant_head=5dcce0f2...` -> `8ecbf836...` (expected, the
   script writes HEAD).
2. **`NonEqOpVX-VXgraphene_goldenKPM_M64x64_statedefault.chebmom2D` (first differing
   file)**. Header identical (`256 15.0366847652085  0`, `64 64`). `cmp`: `differ: byte 85,
   line 4`. First differing element: line 4, `(m,n) = (0,1)` imaginary part,
   new `-9.39092282799606e-09`, committed `-9.39092282756238e-09`. Over the file:
   max absolute difference `9.992e-16` at `(62,62)` real; **max relative difference
   `1.917e-06`** at `(39,45)` imaginary (an element of order 1e-9). The difference is at
   the level of the last bit of O(1) moments (abs <= 1e-15) but exceeds the 1e-12
   relative threshold of AGENTS.md section 8 on the small elements. It is not thread
   dependent here: the same moments rerun at `OMP_NUM_THREADS = 1, 16, 64` are
   byte-identical to each other and to the regenerated file. So the committed file was
   produced by a build that rounds differently (different machine/compiler, or code
   before `5dcce0f`..`8ecbf83`), not by this build. ctest `graphene_regression` passes
   because it compares with `chebmom_compare` at `TOL=1e-6`, not byte for byte.
3. `KuboBastin_NonEqOpVX-VXgraphene_goldenKPM_M64x64_statedefaultJACKSON.dat`: same
   number of lines (1919) but a **different energy grid**: line 1 new `-7.14243 <finite>`,
   committed `-7.51834 nan`. The committed curve has `nan` rows the current
   reconstruction does not produce (max relative difference infinite: nan vs finite).
   The committed file comes from an older reconstruction (grid/pad and alpha safeguard
   changed since); `inline_kuboBastinFromChebmom` at `8ecbf83` cannot reproduce it.
4. `operators/graphene_golden.{HAM,VX}.desc` deleted: `regenerate_goldens.sh` does
   `rm -rf operators` and wannier2sparse does not write these committed sidecars.
   Probe (not a golden): with the committed `.desc` present the driver takes the
   Lanczos bounds `[-3, 3]` from it (header `256 6 -2.22e-16`) and the moments diverge
   (max `|mu| ~ 2e73` at `(63,63)`), i.e. the descriptor bounds do not enclose the
   spectrum; without it the Gershgorin + 10% pad gives `BandWidth = 15.0367` and finite
   moments. Both the committed graphene golden and the regenerated one were made
   without the `.desc`.

Per instructions nothing of this was committed: `git checkout -- test/golden
GOLDEN_PROVENANCE.txt` restored the tree. The graphene files in `goldens/stage0/graphene/`
are this session's regenerated ones (rule: every stage-0 golden produced by the legacy
build in this session); they therefore **differ from `test/golden`** as described above.
For the derived graphene files:

| this session (`goldens/stage0/graphene/`) | vs committed `test/golden/` |
|---|---|
| `fft_grid/KuboBastin_..._FFTgrid.dat` | IDENTICAL to `fft_grid/graphene_KuboBastin_FFTgrid.dat` |
| `fft_grid/KuboGreenWood_..._FFTgrid.dat` | IDENTICAL to `fft_grid/graphene_KuboGreenwood_FFTgrid.dat` |
| `sea_surf/KuboBastinSurf_...dat` | IDENTICAL to `sea_surf/graphene_KuboBastinSurf.dat` |
| `sea_surf/KuboBastinSea_...dat` | DIFFERS: line 1 col 2, new `-2.49859e-13`, committed `-2.49851e-13`; max abs `1.0e-15`, max rel `1.669e-04` (line 1915), from the moment difference above |
| `NonEqOp...chebmom2D`, `KuboBastin_...dat` | DIFFER as in items 2 and 3 |

Fable decides which graphene reference is blessed for G4.

## 4. Stability of the goldens

Second run of the G1, G2 and G5 pipelines in fresh temporary directories, against
`goldens/stage0/`, then G1 and G2 at 1 and 64 threads. Explicit transcript
(`KPM_SEED=12345`; empty `cmp` output and exit 0 = byte-identical):

```
## G1, G2 at OMP_NUM_THREADS=16 (workdir /tmp/tmp.rmhjC8XPO1)
$ cmp /tmp/tmp.rmhjC8XPO1/SpectralOp1chain1dKPM_M64_stateexact.chebmom1D goldens/stage0/chain1d/SpectralOp1chain1dKPM_M64_stateexact.chebmom1D
  exit 0
$ cmp /tmp/tmp.rmhjC8XPO1/NonEqOpVX-VXchain1dKPM_M64x64_stateexact.chebmom2D goldens/stage0/chain1d/NonEqOpVX-VXchain1dKPM_M64x64_stateexact.chebmom2D
  exit 0
## G1, G2 at OMP_NUM_THREADS=1 (workdir /tmp/tmp.ZztksUm4gR)
$ cmp /tmp/tmp.ZztksUm4gR/SpectralOp1chain1dKPM_M64_stateexact.chebmom1D goldens/stage0/chain1d/SpectralOp1chain1dKPM_M64_stateexact.chebmom1D
  exit 0
$ cmp /tmp/tmp.ZztksUm4gR/NonEqOpVX-VXchain1dKPM_M64x64_stateexact.chebmom2D goldens/stage0/chain1d/NonEqOpVX-VXchain1dKPM_M64x64_stateexact.chebmom2D
  exit 0
## G1, G2 at OMP_NUM_THREADS=64 (workdir /tmp/tmp.dfw8HNOwi6)
$ cmp /tmp/tmp.dfw8HNOwi6/SpectralOp1chain1dKPM_M64_stateexact.chebmom1D goldens/stage0/chain1d/SpectralOp1chain1dKPM_M64_stateexact.chebmom1D
  exit 0
$ cmp /tmp/tmp.dfw8HNOwi6/NonEqOpVX-VXchain1dKPM_M64x64_stateexact.chebmom2D goldens/stage0/chain1d/NonEqOpVX-VXchain1dKPM_M64x64_stateexact.chebmom2D
  exit 0
## G5 at OMP_NUM_THREADS=16 (workdir /tmp/tmp.nnaXKOkjkw)
$ cmp /tmp/tmp.nnaXKOkjkw/NonEqOpVX-VXhaldaneKPM_M128x128_stateexact.chebmom2D goldens/stage0/haldane/NonEqOpVX-VXhaldaneKPM_M128x128_stateexact.chebmom2D
  exit 0
$ cmp /tmp/tmp.nnaXKOkjkw/NonEqOpVX-VYhaldaneKPM_M128x128_stateexact.chebmom2D goldens/stage0/haldane/NonEqOpVX-VYhaldaneKPM_M128x128_stateexact.chebmom2D
  exit 0
$ cmp /tmp/tmp.nnaXKOkjkw/KuboBastin_NonEqOpVX-VYhaldaneKPM_M128x128_stateexactJACKSON.dat goldens/stage0/haldane/KuboBastin_NonEqOpVX-VYhaldaneKPM_M128x128_stateexactJACKSON.dat
  exit 0
$ cmp /tmp/tmp.nnaXKOkjkw/KuboBastinSea_NonEqOpVX-VYhaldaneKPM_M128x128_stateexactJACKSON.dat goldens/stage0/haldane/KuboBastinSea_NonEqOpVX-VYhaldaneKPM_M128x128_stateexactJACKSON.dat
  exit 0
$ cmp /tmp/tmp.nnaXKOkjkw/KuboBastinSurf_NonEqOpVX-VYhaldaneKPM_M128x128_stateexactJACKSON.dat goldens/stage0/haldane/KuboBastinSurf_NonEqOpVX-VYhaldaneKPM_M128x128_stateexactJACKSON.dat
  exit 0
$ cmp /tmp/tmp.nnaXKOkjkw/KuboGreenWood_NonEqOpVX-VXhaldaneKPM_M128x128_stateexactJACKSON.dat goldens/stage0/haldane/KuboGreenWood_NonEqOpVX-VXhaldaneKPM_M128x128_stateexactJACKSON.dat
  exit 0
```

Result: **byte-identical, every file, at 1, 16 and 64 threads.** The harness
(section 5) repeated the same comparison three times during the session (12/12
IDENTICAL each time).

Cross-checks found on the way (all byte-identical):

- tutorial 03 run through the modern `lsquant compute` (bounds from the `.desc`) +
  `lsquant reconstruct bastin|greenwood` + the Sea/Surf drivers reproduces the four G5
  `.dat` files and the VX-VY `.chebmom2D` of the inline drivers;
- tutorial 04's `TimeProjSX-PX...chebmomTD` equals G7;
- the second full run of `figures/stage_0.py` reproduced every captured tutorial data file
  and every PNG of the first run (section 6).

## 5. The harness

```
usage: scripts/check_goldens.sh <stage> [build_dir]
```

Files: `scripts/check_goldens.sh` (gate) and `scripts/golden_diff.py` (numeric
localisation: first differing line, element index mapped to `m`, `(m,n)` or `(m,t)` from
the legacy layouts, both values, max absolute and relative difference over the file, and
the class `rounding level (< 1e-12)` / `DEFECT level`). Stage 0 reruns section 4; stage 1
runs `build/lsquant_core/stage1_moments` (falling back to `build/lsquant/stage1_moments`)
with the arguments of `handoff/0001-core-algebra-and-chebyshev/codex.md` section 2 at
`T = 1, 8, 64`, the random-phase determinism pair (`--trace random --seed 12345 --states 4`,
two runs, `cmp`), and the `grep -rn "basis\|std::cout\|getenv" lsquant/core`. The script
works in `mktemp -d` and never writes into `goldens/`. Exit 0 only if every file is
IDENTICAL (and the grep is empty for stage >= 1); usage errors and missing executables
exit 2. Note for stage 1: the gate passes `--seed 12345` (the brief's CLI option) in
addition to `KPM_SEED=12345` in the environment, also on the exact-trace runs.

Output on stage 0 (`bash scripts/check_goldens.sh 0 build`, 2 min wall):

```
== stage 0: legacy pipelines vs goldens/stage0 (KPM_SEED=12345, build /home/ICN2/jgarcia/lsquant-unified/build)
IDENTICAL [T=16] chain1d/SpectralOp1chain1dKPM_M64_stateexact.chebmom1D (G1)
IDENTICAL [T=16] chain1d/NonEqOpVX-VXchain1dKPM_M64x64_stateexact.chebmom2D (G2)
IDENTICAL [T=16] haldane/NonEqOpVX-VYhaldaneKPM_M128x128_stateexact.chebmom2D (G5)
IDENTICAL [T=16] haldane/NonEqOpVX-VXhaldaneKPM_M128x128_stateexact.chebmom2D (G5)
IDENTICAL [T=16] haldane/KuboBastin_NonEqOpVX-VYhaldaneKPM_M128x128_stateexactJACKSON.dat (G5)
IDENTICAL [T=16] haldane/KuboBastinSea_NonEqOpVX-VYhaldaneKPM_M128x128_stateexactJACKSON.dat (G5)
IDENTICAL [T=16] haldane/KuboBastinSurf_NonEqOpVX-VYhaldaneKPM_M128x128_stateexactJACKSON.dat (G5)
IDENTICAL [T=16] haldane/KuboGreenWood_NonEqOpVX-VXhaldaneKPM_M128x128_stateexactJACKSON.dat (G5)
IDENTICAL [T=1] chain1d/SpectralOp1chain1dKPM_M64_stateexact.chebmom1D (G1)
IDENTICAL [T=1] chain1d/NonEqOpVX-VXchain1dKPM_M64x64_stateexact.chebmom2D (G2)
IDENTICAL [T=64] chain1d/SpectralOp1chain1dKPM_M64_stateexact.chebmom1D (G1)
IDENTICAL [T=64] chain1d/NonEqOpVX-VXchain1dKPM_M64x64_stateexact.chebmom2D (G2)
== 12 files compared; failures: 0
GATE PASS (stage 0)
```

The `DIFFERS` path was exercised once on a throw-away copy of `scripts/` and
`goldens/stage0` in `/tmp` with one digit of G1 changed (the committed goldens were not
touched):

```
DIFFERS   [T=16] chain1d/SpectralOp1chain1dKPM_M64_stateexact.chebmom1D (G1)
          cmp: .../g12_T16/SpectralOp1chain1dKPM_M64_stateexact.chebmom1D /tmp/hx/goldens/stage0/chain1d/SpectralOp1chain1dKPM_M64_stateexact.chebmom1D differ: byte 41, line 10
          cmp -l:  41  60  71
          first-diff line 10 m=7 re, new=0 golden=9; max abs diff 9.000e+00 at line 10 m=7 re; max rel diff 1.000e+00 at line 10 m=7 re; DEFECT level (>= 1e-12)
```

Stage 1 cannot run here yet (`ERROR: stage1_moments not found under .../build/lsquant_core
or .../build/lsquant`, exit 2): this build has no new core.

## 6. The stage-0 figure

`figures/stage_0.py` runs `scripts/check_goldens.sh 0`, then every README command
block of the six tutorials (in a work tree outside the repository; `examples/` is never
written) and writes into `figures/out/stage_0/0N/` the PNGs, a caption file per PNG
(`<name>.caption.txt`, beginning `FIG. 0.N`) and `compare.txt` (figure data vs
`goldens/stage0/tutorials/0N`). It was run twice: once to produce the tutorial goldens
(output outside the repository, then copied into `goldens/stage0/tutorials/`), and once
more, from scratch, to produce `figures/out/stage_0/`:

```
scripts/check_goldens.sh 0: 12 IDENTICAL, 0 DIFFERS (G1, G2, G5 at 16 threads; G1, G2 at 1 and 64 threads) -> GATE PASS (stage 0)
tutorial 01: 3 figures, 24 data files, 445 s; figure data vs goldens/stage0/tutorials/01: 24 IDENTICAL, 0 DIFFERS
tutorial 02: 2 figures, 194 data files, 84 s; figure data vs goldens/stage0/tutorials/02: 194 IDENTICAL, 0 DIFFERS
tutorial 03: 3 figures, 9 data files, 52 s; figure data vs goldens/stage0/tutorials/03: 9 IDENTICAL, 0 DIFFERS
tutorial 04: 2 figures, 4 data files, 17 s; figure data vs goldens/stage0/tutorials/04: 4 IDENTICAL, 0 DIFFERS
tutorial 05: 3 figures, 17 data files, 52 s; figure data vs goldens/stage0/tutorials/05: 17 IDENTICAL, 0 DIFFERS
tutorial 06: 3 figures, 281 data files, 4493 s; figure data vs goldens/stage0/tutorials/06: 281 IDENTICAL, 0 DIFFERS
```

All 16 PNGs of the second run are also byte-identical to those of the first (the ones in
`goldens/stage0/tutorials/`). Every README step exited 0.

Figure paths (each with `<name>.caption.txt` next to it, and `compare.txt` per tutorial):

```
figures/out/stage_0/01/fig_levels.png  fig_moments.png  fig_dos.png             FIG. 0.1a-c
figures/out/stage_0/02/fig_anderson.png  fig_localization.png                   FIG. 0.2a-b
figures/out/stage_0/03/fig_bands.png  fig_haldane_hall.png  fig_phase_diagram.png  FIG. 0.3a-c
figures/out/stage_0/04/fig_precession.png  fig_zeeman.png                       FIG. 0.4a-b
figures/out/stage_0/05/fig_diffusion.png  fig_fermi.png  fig_transport.png       FIG. 0.5a-c
figures/out/stage_0/06/fig_dp.png  fig_ey.png  fig_ey_ratelaw.png               FIG. 0.6a-c
figures/out/stage_0/check_goldens_stage0.log                                     the section-4 run
```

(Letters follow the alphabetical order of the PNG names.) Example caption,
`figures/out/stage_0/03/fig_haldane_hall.caption.txt`:

```
FIG. 0.3b (fig_haldane_hall.png). Tutorial 03 (03_haldane_hall), regenerated by the legacy build (HEAD 8ecbf83) through the README pipeline.
Model: Haldane model t1 = -1, t2 = 0.15, phi = pi/2 (test/golden/haldane operators).
Size: dim = 512 (256 cells), bounds [-4.5, 4.5] eV from the .desc.
Moments: M = 128 x 128 (exact trace), Jackson 10 meV.
Seed: KPM_SEED=12345 (exact trace).
Data: figure data vs goldens/stage0/tutorials/03: 9 IDENTICAL, 0 DIFFERS.
Stability (section 4): scripts/check_goldens.sh 0: 12 IDENTICAL, 0 DIFFERS (G1, G2, G5 at 16 threads; G1, G2 at 1 and 64 threads) -> GATE PASS (stage 0).
```

Physics read off the run log: tutorial 03 `Chern from gap plateau: C = +1.017`, k-space
`C = +1.000` (topological) and `+0.000` (trivial); tutorial 05 `sigma` falls from 2.489 to
0.044 and `ell` from 19.66 to 0.35 as `W` goes 0.5 -> 4.

The stage-0 figure is "the six tutorial figures regenerated from the legacy build": the
tutorial scripts draw them from the legacy drivers' output, and `figures/stage_0.py`
proves that output equals the committed goldens byte for byte.

## 7. Goldens committed (`goldens/stage0/`, 7.1 MiB, 604 files)

Sizes in bytes (`find -printf %s`); the tutorial directories are summarised (their full
content is listed in `captured/index.tsv` and by `git ls-tree`). Total 7 488 556 bytes, 604 files.

```
         5  goldens/stage0/chain1d_ballistic/BOUNDS
    363427  goldens/stage0/chain1d_ballistic/CorrelationVXchain1dKPM_M256_state20exact.chebmomTD
      1002  goldens/stage0/chain1d_ballistic/meanCorrelationVXchain1dKPM_M256_state20exactEF0.000000JACKSON.dat
      1001  goldens/stage0/chain1d_ballistic/meanCorrelationVXchain1dKPM_M256_state20exactEF0.500000JACKSON.dat
       998  goldens/stage0/chain1d_ballistic/meanCorrelationVXchain1dKPM_M256_state20exactEF1.000000JACKSON.dat
      1002  goldens/stage0/chain1d_ballistic/meanCorrelationVXchain1dKPM_M256_state20exactEF1.500000JACKSON.dat
      5372  goldens/stage0/chain1d_ballistic/operators/chain1d.HAM.CSR
      5372  goldens/stage0/chain1d_ballistic/operators/chain1d.VX.CSR
         5  goldens/stage0/chain1d/BOUNDS
     41874  goldens/stage0/chain1d/DOS_reconstructed.dat
      3807  goldens/stage0/chain1d_mag/operators/chain1d_mag.HAM.CSR
      3321  goldens/stage0/chain1d_mag/operators/chain1d_mag.PX.CSR
      1331  goldens/stage0/chain1d_mag/operators/chain1d_mag.SX.CSR
      1395  goldens/stage0/chain1d_mag/operators/chain1d_mag.SY.CSR
      1395  goldens/stage0/chain1d_mag/operators/chain1d_mag.SZ.CSR
     65675  goldens/stage0/chain1d_mag/TimeProjSX-PXchain1d_magKPM_M32_stateexact.chebmomTD
     65713  goldens/stage0/chain1d_mag/TimeProjSY-PXchain1d_magKPM_M32_stateexact.chebmomTD
     66128  goldens/stage0/chain1d_mag/TimeProjSZ-PXchain1d_magKPM_M32_stateexact.chebmomTD
     26785  goldens/stage0/chain1d/meanSpectralOp1chain1dKPM_M64_stateexactJACKSON.dat
     18605  goldens/stage0/chain1d/NonEqOpVX-VXchain1dKPM_M64x64_stateexact.chebmom2D
     48443  goldens/stage0/chain1d/operators/chain1d.HAM.CSR
     48443  goldens/stage0/chain1d/operators/chain1d.VX.CSR
       268  goldens/stage0/chain1d/SpectralOp1chain1dKPM_M64_stateexact.chebmom1D
      1095  goldens/stage0/graphene/fft_grid/KuboBastin_NonEqOpVX-VXgraphene_goldenKPM_M64x64_statedefaultJACKSON_FFTgrid.dat
      1107  goldens/stage0/graphene/fft_grid/KuboGreenWood_NonEqOpVX-VXgraphene_goldenKPM_M64x64_statedefaultJACKSON_FFTgrid.dat
      1623  goldens/stage0/graphene/graphene_golden_hr.dat
       147  goldens/stage0/graphene/graphene_golden.uc
       228  goldens/stage0/graphene/graphene_golden.xyz
     32440  goldens/stage0/graphene/KuboBastin_NonEqOpVX-VXgraphene_goldenKPM_M64x64_statedefaultJACKSON.dat
    167007  goldens/stage0/graphene/NonEqOpVX-VXgraphene_goldenKPM_M64x64_statedefault.chebmom2D
     50725  goldens/stage0/graphene/operators/graphene_golden.HAM.CSR
      2867  goldens/stage0/graphene/operators/graphene_golden.SX.CSR
      2995  goldens/stage0/graphene/operators/graphene_golden.SY.CSR
      2995  goldens/stage0/graphene/operators/graphene_golden.SZ.CSR
     55973  goldens/stage0/graphene/operators/graphene_golden.VX.CSR
     55973  goldens/stage0/graphene/operators/graphene_golden.VXSX.CSR
     56293  goldens/stage0/graphene/operators/graphene_golden.VXSY.CSR
     57125  goldens/stage0/graphene/operators/graphene_golden.VXSZ.CSR
     52389  goldens/stage0/graphene/operators/graphene_golden.VY.CSR
     52389  goldens/stage0/graphene/operators/graphene_golden.VYSX.CSR
     52133  goldens/stage0/graphene/operators/graphene_golden.VYSY.CSR
     53029  goldens/stage0/graphene/operators/graphene_golden.VYSZ.CSR
     27941  goldens/stage0/graphene/operators/graphene_golden.VZ.CSR
     27941  goldens/stage0/graphene/operators/graphene_golden.VZSX.CSR
     27813  goldens/stage0/graphene/operators/graphene_golden.VZSY.CSR
     29349  goldens/stage0/graphene/operators/graphene_golden.VZSZ.CSR
     40158  goldens/stage0/graphene/sea_surf/KuboBastinSea_NonEqOpVX-VXgraphene_goldenKPM_M64x64_statedefaultJACKSON.dat
     32438  goldens/stage0/graphene/sea_surf/KuboBastinSurf_NonEqOpVX-VXgraphene_goldenKPM_M64x64_statedefaultJACKSON.dat
         9  goldens/stage0/haldane/BOUNDS
     69689  goldens/stage0/haldane/KuboBastin_NonEqOpVX-VYhaldaneKPM_M128x128_stateexactJACKSON.dat
     69744  goldens/stage0/haldane/KuboBastinSea_NonEqOpVX-VYhaldaneKPM_M128x128_stateexactJACKSON.dat
     81102  goldens/stage0/haldane/KuboBastinSurf_NonEqOpVX-VYhaldaneKPM_M128x128_stateexactJACKSON.dat
      2292  goldens/stage0/haldane/KuboGreenWood_NonEqOpVX-VXhaldaneKPM_M128x128_stateexactJACKSON.dat
    362347  goldens/stage0/haldane/NonEqOpVX-VXhaldaneKPM_M128x128_stateexact.chebmom2D
    362922  goldens/stage0/haldane/NonEqOpVX-VYhaldaneKPM_M128x128_stateexact.chebmom2D
     50613  goldens/stage0/haldane/operators/haldane.HAM.CSR
     52149  goldens/stage0/haldane/operators/haldane.VX.CSR
     70325  goldens/stage0/haldane/operators/haldane.VY.CSR
      7300  goldens/stage0/PROVENANCE.md
   1220206  goldens/stage0/tutorials/01/  (27 files: 3 PNG, 1 top-level .dat, 23 captured/ incl. index.tsv)
    645128  goldens/stage0/tutorials/02/  (196 files: 2 PNG, 1 top-level .dat, 193 captured/ incl. index.tsv)
    633856  goldens/stage0/tutorials/03/  (12 files: 3 PNG, 4 top-level .dat, 5 captured/ incl. index.tsv)
    352352  goldens/stage0/tutorials/04/  (6 files: 2 PNG, 0 top-level .dat, 4 captured/ incl. index.tsv)
   1049593  goldens/stage0/tutorials/05/  (20 files: 3 PNG, 8 top-level .dat, 9 captured/ incl. index.tsv)
    746045  goldens/stage0/tutorials/06/  (284 files: 3 PNG, 25 top-level .dat, 256 captured/ incl. index.tsv)
```

## 8. Failures and things that smell

1. `test/golden` does not regenerate byte-identically (section 3): graphene moments
   differ at max rel `1.9e-6` (abs `1e-15`), the committed Kubo-Bastin curve has a
   different energy grid and `nan` rows, and the regeneration deletes the two committed
   `.desc` sidecars. Not committed; needs Fable's decision for G4.
2. The committed `graphene_golden.HAM.desc` gives bounds `[-3, 3]` that make the
   graphene moments diverge (`|mu| ~ 2e73`) when the sidecar is present (section 3, item 4).
3. `cmake/FetchWannier2Sparse.cmake` reports "wannier2sparse unavailable (offline?)" on a
   cold cache although the fetch and build succeed (stdout of the nested cmake captured
   into `W2S_BIN`).
4. `-DLSQUANT_BUILD_NEW_CORE=OFF` is ignored at `8ecbf83` (option not defined yet).
5. Tutorial 04's README says the operators are shipped under `examples/04/operators`; the
   directory does not exist (staged from `test/models/chain1d_mag`). Tutorial 05's README
   runs `lsquant compute --config dos.json` without defining `dos.json`.
6. Tutorial 06 (qualitative, G8): the rate-law fit gives `p = 1.257` against the
   Elliott-Yafet `p = 2` printed by `rp_validate.py`, and the exact vs random-phase
   cross-check differs by 20.0 % (`T1 = 1004.1` vs `802.9` fs). Recorded as produced;
   nothing changed. Tutorial 06 is also slow with the legacy drivers (`lsqrelax.py` 68 min
   for twelve exact-trace runs on 800 orbitals, single-threaded at under 50 % CPU).
7. `.gitignore` excludes `*.CSR` and `*.chebmom*` except under `test/golden` and
   `test/models`; the stage-0 goldens were added with `git add -f`. Fable may want
   `!goldens/**` in `.gitignore` (not edited here).
8. The legacy `TimeEvProjetedOp` and MSD drivers print the Bessel coefficients of every
   step on stdout (`n:22 Jn:... x:...`), which is why `figures/out/stage_0/run.log` is
   several MB; it is not committed.
