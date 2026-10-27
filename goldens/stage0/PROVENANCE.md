# Stage-0 golden masters — provenance

Every file under `goldens/stage0/` was produced by the **legacy** LSQUANT build in one
session on TCN, following `handoff/0000-goldens/codex.md`. No file was edited by hand.
Inputs (operators, `BOUNDS`, model seeds) are copied from the committed tree and listed
as such.

| Item | Value |
|---|---|
| Machine | `fatnode-tacn` (TCN), AMD EPYC 7551, 64 logical CPUs, 2 TB RAM, Ubuntu 22.04.4 LTS, kernel 5.15.0-161 |
| Compiler / tools | g++ 11.4.0 (Ubuntu 11.4.0-1ubuntu1~22.04.3), cmake 3.22.1, Python 3.10.12, numpy 2.2.6, scipy 1.15.3, matplotlib 3.10.9; Eigen 3.4.0, FFTW 3.3.8, GSL 2.7.1 |
| Repository | `adamecius/lsquant`, branch `cpp-unified`, HEAD `8ecbf836677f8d7a13f72dcdd503225a6ffac820` (clone `~/lsquant-unified`) |
| wannier2sparse | `26cab4affff8266f8a76ec73f3201e6c063a1ba6` (scripts/fetch_wannier2sparse.sh) |
| Build | `cmake -S . -B build -DLSQUANT_BUILD_NEW_CORE=OFF -DLSQUANT_BUILD_INPUT_TOOLS=ON` (CMAKE_BUILD_TYPE=Release, default flags); `cmake --build build -j 32` |
| Seed | `KPM_SEED=12345` for every command, except where a tutorial script sets its own `KPM_SEED` (noted below) |
| Threads | `OMP_NUM_THREADS=16` for every generation command (tutorial 06 `rp_validate.py`: 64, as its README says; `rp_spin.py` sets 16 itself) |
| Date | 2026-10-20 (UTC) |

**Regenerating from `legacy-final`.** The legacy C++ tree (`src/`, `include/`, `test/`),
`scripts/regenerate_goldens.sh`, and the legacy shell scripts were removed in Stage 7 once
every golden was reproduced from observable words. The git tag `legacy-final` (`05e2f65`)
marks the last commit containing the legacy tree and `scripts/regenerate_goldens.sh`. Any
Stage-0 golden is regenerated from a worktree of the tag `legacy-final`:

```bash
git worktree add /tmp/lsquant-legacy legacy-final
cmake -S /tmp/lsquant-legacy -B /tmp/lsquant-legacy/build -DLSQUANT_BUILD_NEW_CORE=OFF -DLSQUANT_BUILD_INPUT_TOOLS=ON
cmake --build /tmp/lsquant-legacy/build -j 32
cd /tmp/lsquant-legacy && KPM_SEED=12345 OMP_NUM_THREADS=16 bash scripts/regenerate_goldens.sh build
```

## G1, G2, G3 — `chain1d/` (exact trace, `M = 64`, window `[-2, 2]`)

Produced by `KPM_SEED=12345 OMP_NUM_THREADS=16 bash scripts/regenerate_goldens.sh build`
(wannier2sparse `chain1d 2048 1 1 VX`, then `inline_compute-kpm-nonEqOp chain1d VX VX 64 exact`,
`inline_compute-kpm-spectralOp chain1d 1 64 exact`, and `test/dos_reconstructed.sh`).
All files of this regeneration under `test/golden/chain1d/` came out byte-identical to the
committed ones; they were copied here:

- `operators/chain1d.HAM.CSR`, `operators/chain1d.VX.CSR` (N = 2048), `BOUNDS` (`-2 2`)
- `SpectralOp1chain1dKPM_M64_stateexact.chebmom1D` (G1)
- `NonEqOpVX-VXchain1dKPM_M64x64_stateexact.chebmom2D` (G2)
- `DOS_reconstructed.dat` (G3: `inline_kpm-DOS-standalone chain1d 1 64 exact`, per-site, nan/inf rows dropped)
- `meanSpectralOp1chain1dKPM_M64_stateexactJACKSON.dat` (G3 raw curve:
  `inline_spectralFunctionFromChebmom SpectralOp1chain1dKPM_M64_stateexact.chebmom1D 120`, Jackson 120 meV)

## G4 — `graphene/` (random phase, `KPM_SEED=12345`, `M = 64`)

Produced by the same `scripts/regenerate_goldens.sh` run (wannier2sparse
`graphene_golden 8 8 1 all`, `inline_compute-kpm-nonEqOp graphene_golden VX VX 64`,
`inline_kuboBastinFromChebmom <mom> 64`). **These files differ from the committed
`test/golden` ones** (see `handoff/0000-goldens/tests.md`, section 3); the committed
`test/golden` was restored and not changed. The files here are this session's output:

- `operators/graphene_golden.*.CSR` (17 operators, dim 256; byte-identical to the committed CSRs),
  `graphene_golden_hr.dat`, `.uc`, `.xyz` (committed wannier2sparse seed inputs)
- `NonEqOpVX-VXgraphene_goldenKPM_M64x64_statedefault.chebmom2D`
- `KuboBastin_NonEqOpVX-VXgraphene_goldenKPM_M64x64_statedefaultJACKSON.dat`
- `fft_grid/`: `inline_kuboBastinFromChebmom_FFTgrid <mom> 10`, `inline_kuboGreenwoodFromChebmom_FFTgrid <mom> 10`
- `sea_surf/`: `inline_kuboBastinSeaFromChebmom <mom> 10`, `inline_kuboBastinSurfFromChebmom <mom> 10`

Note: `scripts/regenerate_goldens.sh` deletes `test/golden/operators/` and wannier2sparse
does not write the two committed sidecars `graphene_golden.HAM.desc` and `.VX.desc`, so the
moments are computed without them (bounds from the Gershgorin + pad safeguard, header
`BandWidth = 15.0366847652085`).

## G4c — `graphene/offcentre_enclosing/` (off-centre window, exact trace, `M = 128`, 10 meV)

Added at the stage-4 closure (handoff 0004, review 3) so that the gate covers a window with
`BandCenter != 0` (`Shift != 0`); every other golden window is centred.

Inputs: `operators/graphene_golden.{HAM,VX,VY}.CSR` copied from `graphene/operators/`
(dim 256, spectrum `[-6.834, 6.834]` eV by dense diagonalisation); `BOUNDS` = `-8.5 7.5` in the
run directory (window source = `BOUNDS` file, branch 1 of the legacy `SpectralBounds`; centre
`-0.5`, encloses the spectrum; header `256 16  -0.5`); an exact local-state file over the 256
sites.

- `BOUNDS`, `NonEqOpVX-VYgraphene_goldenKPM_M128x128_stateexact.chebmom2D`: **legacy**. Produced on
  TCN, 2026-10-20, by this repository's legacy tree (unchanged since `8ecbf83` apart from CMake
  options), built as in the table above in a separate directory: `cmake -S . -B build_legacy
  -DLSQUANT_BUILD_NEW_CORE=OFF -DLSQUANT_BUILD_INPUT_TOOLS=ON`, `cmake --build build_legacy -j 32`
  (it reproduces G5 Bastin byte for byte), with `KPM_SEED=12345 OMP_NUM_THREADS=16`:
  `inline_compute-kpm-nonEqOp graphene_golden VX VY 128 exact`. (`~/lsquant/build` on TCN,
  `c1d9318`, June 2025, is an older unrelated legacy build that does not reproduce G5; not used.)
- `KuboBastin_NonEqOpVX-VYgraphene_goldenKPM_M128x128_stateexactJACKSON.dat`: **blessed from the
  new code because of legacy defect D2** (Fable's decision, handoff 0004 closure):
  `OMP_NUM_THREADS=16 python3 -m lsquant.apps.stage3_kubo --moments
  NonEqOpVX-VYgraphene_goldenKPM_M128x128_stateexact.chebmom2D --formula bastin --broadening-meV 10`
  at `aba179a` (moments from `stage3_kubo --operators ... --left VX --right VY --window -8.5 7.5
  --trace exact --states 256 --moments 128`, byte-identical to the legacy moments; the curve from
  the legacy moments file is byte-identical too). Its value column is byte-identical to the
  legacy `inline_kuboBastinFromChebmom <mom> 10` output; its energy column is the correct
  `x*HalfWidth + BandCenter`. The legacy curve, whose energy column is shifted by
  `BandCenter/HalfWidth - BandCenter = +0.4375` eV, is kept as evidence in
  `legacy_defects/offcentre_axis/` (see its README).

A first G4c window `-3.5 2.5` was dropped: it excludes the spectrum (rescaled eigenvalues up
to 2.44, moments ~1e138), so it is not a golden. Generation log:
`figures/out/stage_4/g4c_generation.log`.

## G5 — `haldane/` (exact trace, `M = 128`, bounds `[-4.5, 4.5]`, 10 meV)

Inputs: `operators/haldane.{HAM,VX,VY}.CSR` copied from `test/golden/haldane/operators`
(dim 512); `BOUNDS` written as `-4.5 4.5` (test/hall_haldane.sh). In a temporary directory
with an exact local-state file over the full basis:

```
inline_compute-kpm-nonEqOp haldane VX VY 128 exact
inline_compute-kpm-nonEqOp haldane VX VX 128 exact
inline_kuboBastinFromChebmom     NonEqOpVX-VYhaldaneKPM_M128x128_stateexact.chebmom2D 10
inline_kuboBastinSeaFromChebmom  NonEqOpVX-VYhaldaneKPM_M128x128_stateexact.chebmom2D 10
inline_kuboBastinSurfFromChebmom NonEqOpVX-VYhaldaneKPM_M128x128_stateexact.chebmom2D 10
inline_kuboGreenwoodFromChebmom  NonEqOpVX-VXhaldaneKPM_M128x128_stateexact.chebmom2D 10
```

## G6 — `chain1d_ballistic/` (exact trace, `M = 256`, 64 times, `tmax = 20`, 10 meV)

Inputs: `operators/chain1d.{HAM,VX}.CSR` from `test/golden/chain1d_ballistic/operators`
(the committed `chain1d.HAM.desc` is deliberately not copied: test/msd_ballistic.sh stages only
the two CSRs and a `BOUNDS` file); `BOUNDS` written as `-2 2`.

```
inline_compute-kpm-MeanSquareDisplacement chain1d VX 256 64 20 exact
inline_timeCorrelationsFromChebmom CorrelationVXchain1dKPM_M256_state20exact.chebmomTD 10 <E>   for E = 0 0.5 1.0 1.5
```

## G7 — `chain1d_mag/` (exact trace, `M = 32`, 96 times, `tmax = 95`)

Inputs: `operators/chain1d_mag.{HAM,SX,SY,SZ,PX}.CSR` from `test/models/chain1d_mag/operators`.

```
inline_compute-kpm-TimeEvProjetedOp chain1d_mag <OP> PX 32 96 95 exact    for OP = SX SY SZ
```

## G8 and the tutorial figure data — `tutorials/01` ... `tutorials/06`

Produced by `figures/stage_0.py` (run with `--out`/`--work` outside the repository), which runs
every README command block of `examples/0N` in a fresh copy of the tutorial (inputs only:
`*.py`, `*.json`, `README.md`, `requirements.txt`) with `LSQUANT_BIN=build`,
`PATH=build:$PATH`, `KPM_SEED=12345`, `OMP_NUM_THREADS=16`, `MPLBACKEND=Agg`. The exact
command list per tutorial is the `TUTORIALS` table in `figures/stage_0.py`.

Per tutorial directory:

- `fig_*.png`: the figures the tutorial scripts write;
- `*.dat`: curve files left in the tutorial directory by the README commands and scripts;
- `captured/NNNN_<name>`: every `.dat` / `.chebmom1D` / `.chebmom2D` / `.chebmomTD` file the
  tutorial Python scripts read (through `numpy.loadtxt` or `open`), in read order, copied
  by the capture shim of `figures/stage_0.py` before the scripts delete them;
  `captured/index.tsv` maps each to its SHA-1 and its source (`<tmpdir>/` = a
  `tempfile.TemporaryDirectory` of the script).

Parameters a README does not fix, and deviations:

- 01: README step 2 "Repeat for chain1d_N1024 and chain1d_N4096" run as a loop over the
  three labels with the same JSON fields.
- 02: `lsqloc.py` sets `KPM_SEED` to the disorder seed (1..16) for each realisation.
- 04: the README says the operators are shipped under `examples/04/operators`; that
  directory does not exist. Operators staged from `test/models/chain1d_mag/operators` and an
  exact local-state file written (as test/spin_precession.sh does).
- 05: README step 3 references a `dos.json` the README never defines; not run standalone
  (`lsqtransport.py` writes and runs it). `lsqtransport.py` sets `KPM_SEED=1`.
- 06: `rp_validate.py` with `OMP_NUM_THREADS=64` (README); `rp_spin.py` sets `KPM_SEED=1..20`
  and `OMP_NUM_THREADS=16` per random vector itself.

**Tutorial 06 denominators: re-blessed from the new code because of legacy defect D5**
(Fable's decision, handoff 0007 closure, as for G4c / D2). 142 files of `tutorials/06` were
replaced; the legacy versions are kept as evidence in
`legacy_defects/tutorial06_identity_operator/` with the same relative layout (see its
`README.md`):

- the 139 curves `meanTimeProj1-P*` (top level and `captured/`);
- `ey_ratelaw.dat` and `captured/0254_ey_ratelaw.dat`;
- `captured/index.tsv`.

What the files hold:

- The replaced curves are the denominator word `Tr[P U(t)^+ 1 U(t) delta(E-H) P]` of
  `lsqrelax.py` and `rp_spin.py`.
- The legacy pipeline ran `inline_compute-kpm-TimeEvProjetedOp <label> 1 <proj> ...`, which reads
  `operators/<label>.1.CSR` (`legacy-final:src/inline_compute-kpm-TimeEvProjetedOp.cpp:53-62`).
  `make_ey.py` and `make_rashba.py` never write that file.
- `SparseMatrixBuilder::BuildOPFromCSRFile` returns silently when the file cannot be read
  (`legacy-final:include/sparse_matrix.hpp:132-133`). The operator `1` was therefore never built,
  and the legacy denominators are about 1e-6 noise.
- The new code applies the identity: the curve is constant in t and equal to `Tr[P delta P]`.
- The tutorials normalise by `mean(den)` and then by `s[0]`, so the scale cancels. The stage-0
  figures were right, and the new `ey_ratelaw.dat` agrees with the legacy one to 1.8e-16 relative.
- The 139 numerator curves (`TimeProjSZ-*`, `TimeProjSX-*`) were byte-identical between legacy and
  new code and are unchanged.

How the new files were produced:

- Command: `KPM_SEED=12345 OMP_NUM_THREADS=16 PYTHONPATH=<build>/lsquant_core:python python3 -m
  lsquant run runs/goldens/tutorial06.json --out <dir>`, at `5e8186d` (tester merge `9cdbe0f`),
  build `~/lsquant-unified-build7c` on TCN, 2026-10-21. The log is the stage-7 closure gate
  `~/s7logs/gate7c.log`.
- `captured/index.tsv` keeps the stage-0 names and source column; its SHA-1 column is recomputed
  from the new captures (128 of 255 entries change).

## G9 — `bases/graphene/`, `bases/haldane/` (stage 5; independent oracle, not the legacy build)

G9 is the only golden not produced by the legacy LSQUANT build: AGENTS.md section 8 names
`wannier2sparse/tools/hr_exactdiag.py` (exact diagonalisation on a k mesh, independent code)
as its oracle. Regenerated on TCN on 2026-10-21 by the stage-5 tester for the closure run
(handoff `0005-bloch-and-eigen-bases`, codex.md section 2.2 as revised in review 2:
Gamma-centred mesh, Gaussian 40 meV), Python 3.10.12, numpy 2.2.6. The run-1 G9 (shifted
48 x 48 mesh, Gaussian pi a / M) is superseded.

| Item | Value |
|---|---|
| wannier2sparse | `~/wannier2sparse` at `3f553db8dc2f226f44592b6ee2a717d3cd71d8c7` (local edits only in `tools/block_pilot.cpp`, `docs/FINDINGS.md`, `roadmap_plans.md`; `tools/hr_exactdiag.py`, `examples/gen_models.py` and `src/` clean). The pinned `26cab4a` of `scripts/fetch_wannier2sparse.sh` has no `tools/hr_exactdiag.py`. |
| Models | `python3 examples/gen_models.py` (graphene: t = -1 eV, a1 = (1.5, 0.8660254, 0), a2 = (1.5, -0.8660254, 0), tau_B = (1, 0, 0); Haldane: t1 = -1, t2 = 0.15, phi = pi/2). Committed unchanged as `<model>/model/<model>_hr.dat`, `.uc`, `.xyz`. |
| Velocities | `wannier2sparse <model> 1 1 1 VX VY --mode bundle -o bund` (velocity mode `berry_connection`, `-i (R.lat + tau_b - tau_a) H`, eV Angstrom); `bund/<model>.w2sp/operators/VX.hr.dat`, `VY.hr.dat` copied to `<model>/model/<model>_vx_hr.dat`, `_vy_hr.dat` (the names `hr_exactdiag.py sigma` reads). |
| Mesh (G9) | **Gamma-centred 48 x 48, k = n/48** (the k mesh of the 48 x 48 x 1 supercell). `hr_exactdiag.py` hard-codes the shifted mesh (n + 1/2)/N, so it is run through the thin wrapper **`scripts/hr_exactdiag_gamma.py`**, which reads the tool's source, replaces exactly the four expressions `(np.arange(nk) + .5) / nk` / `(np.arange(a.nk) + .5) / a.nk` by `np.arange(nk) / nk` / `np.arange(a.nk) / a.nk` (it aborts if the expected occurrences are not found) and runs the tool's own `main`; nothing else of the tool changes. |
| Broadening | Gaussian **eta = 0.040 eV** (40 meV). |
| Grid | `--emin -3.6 --emax 3.6 --ngrid 2881` (2.5 meV spacing; [-3.6, 3.6] eV is the `--window eigen` of the lsquant runs). |
| sigma_xx units | `hr_exactdiag sigma` divides by A_hr = `abs(a1_x a2_y)` = 1.299038 A^2, half the honeycomb cell `abs(a1 x a2)` = 2.598076 A^2, so `sigma_xx_e2h` of these two models is **twice** the e^2/h value. The gate compares `EigenBasis.exact_sigma_xx` (the stage-3 `kubo_greenwood` normalisation, pi Tr[v delta v delta]) multiplied by **2 pi / (A_hr N_k)** = 2.0993e-3 (N_k = 2304) with the G9 values as stored; the figure shows both in e^2/h with the true cell (G9 x A_hr / A). |
| Finite size (report only) | the same `dos` and `sigma --comp xx` with the tool unwrapped on its own shifted 240 x 240 mesh, same eta and grid: `finite_size_shifted240_dos.json`, `finite_size_shifted240_sigma_xx.json`. No gate row uses them as a reference. |

Commands (in a scratch copy of the generated models, `<m>` = graphene, haldane; `T` =
`~/wannier2sparse/tools/hr_exactdiag.py`; `A` = `--eta 0.040 --ngrid 2881 --emin -3.6 --emax 3.6`):

```
python3 scripts/hr_exactdiag_gamma.py T dos   <m> --nk 48 A            --out goldens/stage0/bases/<m>/G9_dos
python3 scripts/hr_exactdiag_gamma.py T sigma <m> --nk 48 A --comp xx  --out goldens/stage0/bases/<m>/G9_sigma_xx
python3 T dos   <m> --nk 240 A           --out goldens/stage0/bases/<m>/finite_size_shifted240_dos
python3 T sigma <m> --nk 240 A --comp xx --out goldens/stage0/bases/<m>/finite_size_shifted240_sigma_xx
```

The whole generation is `goldens/stage0/bases/gen_g9.sh <repo>` (committed). Files per model:
`G9_dos.json` (`energy_eV`, `dos` per cell, integral 2.000), `G9_sigma_xx.json` (`energy_eV`,
`sigma_xx_e2h`), `G9_params.txt` (mesh, nk, eta, grid; read by `figures/stage_5_table.py`),
the two finite-size files, `model/`. A second generation into a scratch tree reproduced all
20 files byte for byte.

sha256:

```
8c8a857d35291856832f54c6b14ab2614d65247260ff910b3d50c285acb12625  bases/graphene/G9_dos.json
a3ca1a5a32c55717024ce5d8e97ebf21688f1d8babfeb81d1e5eded896afa090  bases/graphene/G9_sigma_xx.json
53711db6669b96996381d19f2ab09e3e07cee9837c49b1d49b787e0a13e0d2af  bases/graphene/G9_params.txt
9ef13d31a06ce8387485107592807fd612dcd65a3e56f73d31261836c8fa289e  bases/graphene/finite_size_shifted240_dos.json
2d25877f358c2cf3493c4d5bcc5dad99ac8301368137b34d1f7380f159dda2a9  bases/graphene/finite_size_shifted240_sigma_xx.json
1dc70fea90f60fc8340ee0ced64639a6d9390f7f636a5eea10c06dc333ab2fde  bases/graphene/model/graphene_hr.dat
0bc2d2218e7805356c364afce524055907c1fccc287ef88cff54e6e4ba022560  bases/graphene/model/graphene.uc
2779e3f2fecf3dcff54fad0817088b4b9ab536293cb7e30f44b525a13563cdf9  bases/graphene/model/graphene.xyz
b3014714d39c0354052b6d0308b6b3a92bbb0da6d0cda0ee6d032a7b443cc2b6  bases/graphene/model/graphene_vx_hr.dat
204b566f8a05dc3405219a49630b204beaed79a2d4cbba06324783951e53d8f5  bases/graphene/model/graphene_vy_hr.dat
18f61df5b4bfa8ccad06831488558d627e9e24a86c720a36e1d4bc5221516c58  bases/haldane/G9_dos.json
f7a8f11c7717eedb60259141495f7fc7e207a431c8488ad49e63fc9295a7764a  bases/haldane/G9_sigma_xx.json
53711db6669b96996381d19f2ab09e3e07cee9837c49b1d49b787e0a13e0d2af  bases/haldane/G9_params.txt
e213a220eebd2abe951b4f7dcac3ed4c644630772d17048d5582bc7dca095be6  bases/haldane/finite_size_shifted240_dos.json
34d3455330a7b5173b4ff9e2374f12ec1f132a0f07d81958aede2236e7fbcdfb  bases/haldane/finite_size_shifted240_sigma_xx.json
d4b53b22c4e33e067ad0ca77a6ac66c0aa76d1379d2d265e0278700dcd83d706  bases/haldane/model/haldane_hr.dat
0bc2d2218e7805356c364afce524055907c1fccc287ef88cff54e6e4ba022560  bases/haldane/model/haldane.uc
2779e3f2fecf3dcff54fad0817088b4b9ab536293cb7e30f44b525a13563cdf9  bases/haldane/model/haldane.xyz
b7bdaaad6a393fc5da52d9d0109f0872dc5c86f5771c2e8a688f3523627d05b8  bases/haldane/model/haldane_vx_hr.dat
ea86f1903319c27467003d92f6e332248ae9a7dfcb991a38571c552054f18c10  bases/haldane/model/haldane_vy_hr.dat
```

## Stability

`scripts/check_goldens.sh 0 build` reran G1, G2, G5 (16 threads) and G1, G2 (1 and 64
threads) in fresh temporary directories: every file byte-identical (see
`handoff/0000-goldens/tests.md`, section 4).

## Legacy defects (evidence only)

Files under `chain1d/legacy_defects/` are legacy outputs known to be wrong. They are kept as
evidence and **are not references**: no harness reads them.

- `chain1d/legacy_defects/batch16/NonEqOpVX-VXchain1dKPM_M64x64_stateexact.chebmom2D`:
  G2 recomputed with `BATCH_SIZE=16`
  (`BATCH_SIZE=16 KPM_SEED=12345 OMP_NUM_THREADS=16 inline_compute-kpm-nonEqOp chain1d VX VX 64 exact`,
  legacy build at d9becf6). It is stable at 16 and 1 threads but wrong: legacy
  `Vectors::IterateAll` does not carry the recursion across batches, so block k is shifted
  by k orders. The command, the evidence and the cause are in its `README.md`. The batched
  run of the new core is gated against G2 itself (`scripts/check_goldens.sh 1`).
- `legacy_defects/offcentre_axis/KuboBastin_NonEqOpVX-VYgraphene_goldenKPM_M128x128_stateexactJACKSON.dat`:
  the legacy Bastin curve of G4c (`BOUNDS` `-8.5 7.5`), whose energy column carries legacy defect
  D2 (uniform-grid axis `x*HalfWidth + BandCenter/HalfWidth`, shifted by +0.4375 eV here); see
  its `README.md`. The G4c curve golden was blessed from the new code instead.
- `legacy_defects/tutorial06_identity_operator/` (142 files): the legacy tutorial-06 denominator
  curves `meanTimeProj1-P*`, `ey_ratelaw.dat` and `captured/index.tsv`. They were computed by the
  legacy `inline_compute-kpm-TimeEvProjetedOp` with an operator `1` that was never built, because
  `operators/<label>.1.CSR` was missing and `BuildOPFromCSRFile` returned silently (legacy defect
  D5). See its `README.md`; the tutorial-06 goldens were re-blessed from the new code instead.

## G10 (part) — `nonorth/` (rotated chain1d and Haldane, handoff 0006)

Not legacy output: generated by the developer with
`python3 python/lsquant/apps/make_nonorth_fixtures.py` (seed 12345, strength `alpha = 0.55`)
from `chain1d/operators` and `haldane/operators`: `H' = T^+ H T`, `S' = T^+ T`,
`V' = T^+ V T`, variants `real` and `complex`; each directory carries `manifest.json` with
the sha256 of every `.CSR`. The reference for these inputs is the orthogonal G1, G2, G5
(traces are basis-invariant), never a file computed from them. Verified independently by
the tester (`figures/stage_6_fixtures.py`, dense numpy/scipy, log
`figures/out/stage_6/fixture_verification.log`, handoff `0006` `tests.md` section 1): sha256
equal to the manifests; S Hermitian, positive definite, cond(S) 5.22 to 6.70 equal to the
manifests to 2e-14; `eigh(H', S')` equal to `eigvalsh(H)` to 1.5e-14; H', V' exactly
Hermitian. Accepted as committed at `a5108e3`; not regenerated.
