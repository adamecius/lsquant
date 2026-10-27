# Handoff 0001 (Codex) — stage 1 gate: byte-for-byte against G1 and G2, on TCN

You are the tester. You do not fix source; you find and localise. Read `AGENTS.md`
sections 8, 9 and 12 and `brief.md` of this handoff first. Everything runs on TCN
(`ssh ionos`, then `ssh jh-x`; checkout `~/lsquant-unified`, branch `cpp-unified`,
see `handoff/0000-goldens/codex.md` section 1). Stage 0 must be closed and
`goldens/stage0/` committed before you start.

## 1. Build

```
cd ~/lsquant-unified && git pull
cmake -S . -B build -DLSQUANT_BUILD_NEW_CORE=ON
cmake --build build -j 32 2>&1 | tee build.log
```

Record compiler warnings under `lsquant/` (there must be none with `-Wall -Wextra`).
Then run Agi's self-check and record its output verbatim:

```
OMP_NUM_THREADS=1 build/lsquant/stage_1_check
OMP_NUM_THREADS=8 build/lsquant/stage_1_check
```

## 2. The gate

Extend `scripts/check_goldens.sh` for stage 1. It runs, in a temporary directory:

```
OMP_NUM_THREADS=$T KPM_SEED=12345 build/lsquant/stage1_moments \
   --operators goldens/stage0/chain1d/operators --label chain1d --window -2 2 \
   --trace exact --states 2048 --moments 64 --word dos  --out $TMP
OMP_NUM_THREADS=$T KPM_SEED=12345 build/lsquant/stage1_moments \
   --operators goldens/stage0/chain1d/operators --label chain1d --window -2 2 \
   --trace exact --states 2048 --moments 64 --word kubo --out $TMP
```

for `T` in `1 8 64`, and compares:

```
cmp $TMP/SpectralOp1chain1dKPM_M64_stateexact.chebmom1D        goldens/stage0/chain1d/SpectralOp1chain1dKPM_M64_stateexact.chebmom1D
cmp $TMP/NonEqOpVX-VXchain1dKPM_M64x64_stateexact.chebmom2D    goldens/stage0/chain1d/NonEqOpVX-VXchain1dKPM_M64x64_stateexact.chebmom2D
```

plus the random-phase determinism check (`--trace random --states 4`, run twice,
`cmp` of the two outputs; there is no golden for it, identity of the two runs is
the criterion), plus `grep -rn "basis\|std::cout\|getenv" lsquant/core` (must be
empty). The gate passes only if every `cmp` is silent at every `T`.

## 3. When a `cmp` differs: localise, do not fix

Do this in order and write every result in `tests.md`:

1. `cmp -l new golden | head -5`: first differing byte; map it to line number and
   element index `m` (rank 1) or `(m, n)` (rank 2) of the file layout.
2. Numeric diff of the whole file (`scripts/` Python helper): max absolute and max
   relative difference, the index where it occurs, and whether the differences are
   at the 1e-16 level (rounding or summation order) or larger (logic).
3. Classify the first differing element:
   - header line differs: formatting of `WriteLegacyChebmom1D/2D` (spaces,
     precision, `-0`), or the window values (`BandWidth`, `BandCenter` computed
     differently);
   - `m = 0` differs, `m > 0` fine: the `(2 - delta_m0)` weight or the trace
     normalisation (`N_states`, state order);
   - all `m` differ at rounding level: the matvec summation order
     (`SparseMatrix::Multiply` vs `src/eigen_sparse_matrix.cpp`), or
     `RescaledOperator` not using the rescaled-copy expression, or `dot` not being
     `Eigen::Map(x).dot(Map(y))`;
   - all `m` differ at logic level: the recursion start (`T_1 = H~ phi`),
     the `Step` (`2 H~ T_m - T_{m-1}`), or the operator applied at the wrong end
     (bra: `A` before `Expand`; ket: `Expand` then `A`);
   - rank 2 differs only off the diagonal or only for `m != n`: the
     symmetrisation loop (`ApplyKpmWeights`) or the bra/ket assignment of `left`
     and `right`;
   - rank 2 differs at block boundaries of `batchSize`: the batching order or the
     rounding-up of `M`;
   - differs only at `T > 1`: a parallel reduction that is not per-row, or a
     race.
4. Confirm the classification with one targeted run: for example `--moments 2` or
   `--states 1` isolates weights from recursion; comparing the `dos` output of
   `--trace exact --states 1` against `<e_0| T_m |e_0>` computed with the legacy
   `inline_compute-kpm-spectralOp` and a one-state `local` file isolates the
   recursion from the trace.
5. Name the file and function in `lsquant/` where you believe the difference is
   introduced, and the legacy line it should mirror (`src/...:line`). This is a
   hypothesis for Fable and Agi, stated as such.

Never change a golden, a seed, a tolerance, or the harness criterion to make the
gate pass. Never edit `lsquant/`.

## 4. The stage-1 figure

`figures/stage_1.py` reads `goldens/stage0/chain1d/*.chebmom1D|2D` and the new
outputs (via the `.chebmom` files `stage1_moments` writes alongside) and draws one
figure, two panels: (a) `mu_m`, `m = 0..63`: golden as markers, new as a line,
`delta_m0` as the analytic reference, residual `new - golden` on a secondary axis;
(b) the `m, n <= 10` block of `mu_mn`: new minus golden as a heat map, with the
closed-form B.15 block (`test/chain1d_analytic_test.cpp`, `Mref`, times the
`m = 2..9` diagonal plateau) beside it. Save `figures/out/stage_1/fig_stage_1.png`
and a caption file beginning `FIG. 1.` with `N = 2048`, window `[-2, 2]`, exact
trace, `M = 64`, the thread counts tried and the `cmp` result.

## 5. Report

`handoff/0001-core-algebra-and-chebyshev/tests.md`: build log summary and
warnings; `stage_1_check` output at 1 and 8 threads; for each golden and each `T`,
`IDENTICAL` or the full section-3 localisation; the random-phase determinism
result; the `grep` result; the figure path and caption; anything that smells in the
implementation (hard-coded sizes, a physics name under `core/`, an environment
variable read in the core, a `TODO` for an out-of-scope feature). Then stop; Fable
checks and relays to Agi.
