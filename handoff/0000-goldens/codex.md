# Handoff 0000 (Codex) — stage 0: golden masters from the legacy build, on TCN

You are the tester. You do not write project source. Read `AGENTS.md` sections 8,
9 and 12 first. Everything below runs on TCN; a result produced anywhere else is
not a result.

## 1. The machine

TCN is `fatnode-tacn` (Ubuntu 22.04, 64 cores, 2 TB RAM, g++ 11.4, cmake 3.22,
Python 3.10 with numpy 2.2, scipy 1.15, matplotlib 3.10; Eigen3, FFTW3, GSL and
Boost dev packages installed; `sympy` not installed, not needed in this stage).
It is reached through the Ionos jump host:

```
ssh ionos            # root@31.70.157.48
ssh jh-x             # from ionos: jgarcia@fatnode-tacn
```

The existing `~/lsquant` there is a different repository (`icn2lsquant/lsquant`,
branch `3.0.0_beta`): do not touch it. Work in `~/lsquant-unified`:

```
git clone -b cpp-unified https://github.com/adamecius/lsquant.git ~/lsquant-unified
cd ~/lsquant-unified
git submodule status || true        # none expected; wannier2sparse is fetched by scripts/fetch_wannier2sparse.sh
```

Use `OMP_NUM_THREADS=16` unless a step says otherwise, and `KPM_SEED=12345` always.

## 2. Build the legacy code

```
cmake -S . -B build -DLSQUANT_BUILD_NEW_CORE=OFF -DLSQUANT_BUILD_INPUT_TOOLS=ON
cmake --build build -j 32
ctest --test-dir build --output-on-failure
```

Record in `tests.md`: compiler and cmake versions, the exact commands, the full
`ctest` summary (expected: every test passes). If a test fails, record it and stop;
do not patch anything.

`-DLSQUANT_BUILD_INPUT_TOOLS=ON` fetches wannier2sparse at the pinned sha
(`scripts/fetch_wannier2sparse.sh`); confirm the sha printed is `26cab4a`.

## 3. Generate the goldens

Run the pinned regeneration, which produces G1, G2, G4 and the chain1d
reconstructed DOS:

```
KPM_SEED=12345 bash scripts/regenerate_goldens.sh build
git status --short test/golden GOLDEN_PROVENANCE.txt
```

Expected: no change under `test/golden` (the committed goldens were produced by the
same pipeline). If `git diff` shows a change, report the first differing file and
the max relative difference; do not commit it.

Then produce the remaining goldens by running the tutorial pipelines exactly as
their READMEs say (`examples/01` to `06`), with `LSQUANT_BIN=$PWD/build`,
`KPM_SEED=12345`, and the sizes stated in each README. Collect, under
`goldens/stage0/`, every `.chebmom1D`, `.chebmom2D`, `.chebmomTD`, every
reconstructed `.dat`, and the figure PNGs, with this layout:

```
goldens/stage0/
  chain1d/operators/*.CSR            copied from test/golden/chain1d/operators
  chain1d/BOUNDS
  chain1d/*.chebmom1D *.chebmom2D    G1, G2 (copied from test/golden/chain1d)
  chain1d/DOS_reconstructed.dat      G3
  graphene/...                       G4: operators, .chebmom2D, KuboBastin .dat, fft_grid, sea_surf
  haldane/...                        G5: operators, NonEqOpVX-VY and VX-VX .chebmom2D, Bastin/Sea/Surf/Greenwood .dat
  chain1d_ballistic/...              G6: operators, .chebmomTD, mean...EF*.dat for E = 0, 0.5, 1.0, 1.5
  chain1d_mag/...                    G7: operators, TimeProjSX/SY/SZ...chebmomTD
  tutorials/01 ... 06/               every .dat the tutorial scripts read for their figures, and the PNGs
  PROVENANCE.md                      commands, seeds, sizes, HEAD sha, machine, thread count, date
```

Every file in `goldens/stage0/` is produced by the legacy build in this session and
nothing else. Do not hand-edit a golden. Do not regenerate one twice with different
settings; if a tutorial needs a parameter the README does not fix, write it in
`PROVENANCE.md`.

## 4. Verify the goldens are stable

Run the G1, G2 and G5 pipelines a second time in a fresh temporary directory and
`cmp` the outputs against `goldens/stage0/`. Expected: byte-identical. Then rerun
G1 and G2 with `OMP_NUM_THREADS=1` and `64`: byte-identical. Report any deviation
with the first differing byte (`cmp -l | head`) and the element it corresponds to.

## 5. The harness

Write `scripts/check_goldens.sh`:

```
usage: scripts/check_goldens.sh <stage> [build_dir]
```

For stage `0` it reruns section 4. For stage `1` and later it runs the stage's gate
executable on the golden inputs (the brief of each stage names the executable and
its arguments) into a temporary directory and `cmp`s each produced file against
`goldens/stage0/...`. Output: one line per file, `IDENTICAL` or `DIFFERS`, and for
`DIFFERS` the first differing line number, the element index, the two values, and
the max relative difference over the file (a small Python helper under
`scripts/` may do the numeric part). Exit code 0 only if every file is
`IDENTICAL`. It also runs `grep -rn "basis\|std::cout\|getenv" lsquant/core` and
fails if anything is found (stage 1 onward).

The script never writes into `goldens/` and never modifies a golden.

## 6. The stage-0 figure

`figures/stage_0.py`: regenerate the six tutorial figures from the goldens (the
tutorial scripts already do this; wrap them) into `figures/out/stage_0/`. Caption file
next to each PNG, beginning `FIG. 0.N`, listing model, size, `M`, seed, and the
`cmp` result of section 4.

## 7. Report

`handoff/0000-goldens/tests.md`: machine, versions, commands, `ctest` summary,
regeneration result (`git status` of `test/golden`), the list of files under
`goldens/stage0/` with sizes, the stability results of section 4 with the exact
`cmp` output, the harness usage line and its output on stage 0, the figure paths.
Then stop; Fable checks the report before the goldens are committed.
