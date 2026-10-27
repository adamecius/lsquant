# Legacy defect: batched `CorrelationExpansionMoments` (evidence only, not a reference)

`NonEqOpVX-VXchain1dKPM_M64x64_stateexact.chebmom2D` in this directory is the legacy
chain1d Kubo-Greenwood moment matrix computed with `BATCH_SIZE=16` (M = 64). **It is wrong.**
It is kept as evidence of the defect. **No harness reads it.** `scripts/check_goldens.sh`
compares the new core's batched run (`--batch 16`) with the unbatched G2,
`goldens/stage0/chain1d/NonEqOpVX-VXchain1dKPM_M64x64_stateexact.chebmom2D`. A correct
batched engine is byte-identical to G2: each element gets one dot product per state,
summed in state order, whatever the batch size.

## How it was produced (TCN, legacy build of `~/lsquant-unified` at d9becf6)

The run used a fresh temporary directory, set up as the chain1d part of
`scripts/regenerate_goldens.sh`:

```
mkdir -p operators && cp goldens/stage0/chain1d/operators/*.CSR operators/
printf -- "-2 2\n" > BOUNDS
{ echo local; echo 2048; seq 0 2047; } > exact
BATCH_SIZE=16 KPM_SEED=12345 OMP_NUM_THREADS=16 build/inline_compute-kpm-nonEqOp chain1d VX VX 64 exact
```

The driver logs `Using user's BATCH_SIZE = 16`.

## Stability

The result is stable:
- a second run at `OMP_NUM_THREADS=16` is byte-identical (`cmp`);
- a run at `OMP_NUM_THREADS=1` is byte-identical.

sha256 `fbf3e69745ae9296d73d5fb7922a88de08373dfa354a25a431107704913c3abf`, 18611 bytes.

## What is wrong

- **Size of the damage:** compared with G2, 36 elements differ, all in the off-diagonal
  16 x 16 blocks (0,1), (1,0), (1,2), (2,1), (2,3), (3,2). The four diagonal blocks equal
  G2 exactly.
- **First differing element:** (m, n) = (13, 16), real part: -2 here, 0 in G2 (byte 3541,
  line 851).
- **Largest difference:** max abs 4.0 at (15, 16).
- **Pattern:** block k is shifted by k orders. Row 15 has 4 at n = 15 and n = 16 and -2 at
  n = 18; the closed form B.15 (and G2) has 4 at n = 15 and -2 at n = 17.
- **Why the diagonal blocks look right:** both sides shift by the same amount, and the
  B.15 matrix depends only on m - n once m >= 2.

## Cause

`chebyshev::Vectors::IterateAll`, `src/chebyshev_vectors.cpp` lines 353-364:
- It writes `Vector(0) = Chebyshev0`, then for m = 1..B-1 writes `Vector(m) = Chebyshev1`
  and steps the recursion.
- When the loop ends, the recursion state is (T_{B-1}, T_B).
- `chebyshev::CorrelationExpansionMoments`, `src/chebyshev_solver.cpp` lines 550-561, calls
  `IterateAll` once per block without restarting the recursion.
- So the next block's `Vector(0)` is T_{B-1} again: block k starts at T_{k(B-1)}, not at
  T_{kB}.
- With the default `BATCH_SIZE` (= M, one block) the defect does not show, so G2 is correct.

Listed in AGENTS.md section 10 (known defects, do not port).
