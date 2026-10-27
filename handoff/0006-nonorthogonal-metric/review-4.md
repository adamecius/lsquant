# Handoff 0006 — review 4 (Fable to Agi): the closure run on `8305f11`

Tester (`handoff/0006-codex` at `06384ce`, fat node `build6`, g++ 11.4, plus Aron for the
thread comparison): stages 1 to 5 byte-identical (103 files, 52 PASS stage-5 rows, greps
empty); stage 6 74 of 74 PASS against G1, G2, G5 and route against route (worst
lowdin/chebyshev 9.5e-12, generalised/exact 1.3e-14); the ten Hall plateaus give
C = +1.017017 exactly as G5 on both routes, real and complex; the naive estimator (evidence)
gives +1.19 and fails the Kubo rows by O(1) while passing the DOS rows, which is the
asymmetric-bra argument of `theory.md` made visible. The closure gate is accepted on those
numbers. One finding stays open and is the only item of this review.

## Finding: `ExactInverseSqrtOverlap` depends on the thread count

The four `lowdin/exact` chain1d runs (dos and kubo, real and complex) are not byte-identical
between `OMP_NUM_THREADS=1` and `16`; the other 20 rows are. Differences are rounding level
(max 1e-14 on moments of order 1), the same on the fat node and on Aron (gcc 9.4, other
Eigen), so it is the algorithm, not the machine. `lowdin/chebyshev` (same `LowdinOperator`,
same engine) and `generalised/exact` (LLT and solve) are byte-identical, which localises it
to the constructor of `ExactInverseSqrtOverlap` (`overlap_appliers.cpp:206-217`): the dense
product `V diag(lambda^-1/2) V^+` after `SelfAdjointEigenSolver`, whose Eigen GEMM blocking
depends on the number of threads. The engine rule (AGENTS.md section 8: moments byte-identical
between thread counts on the real path) applies to the appliers too.

## Items for Agi, one commit each, on `handoff/0006-agi` after `8305f11`

1. **Make the exact `S^-1/2` deterministic.** The product that forms `S^-1/2` (and, for
   symmetry, the one that forms `S^-1` in `ExactInverseOverlap` if it is a GEMM) is
   computed once per basis; compute it single-threaded (a scoped `Eigen::setNbThreads(1)`
   around the product, or an explicit non-parallel loop) so that the matrix is independent
   of `OMP_NUM_THREADS`. Do not touch the apply path. Then the four rows are byte-identical
   at 1 and 16 threads; the values against G1, G2, G5 may move at the 1e-14 level (report
   the worst).
2. **`tests/dev/stage_6_check.py` takes `LSQUANT_BUILD_DIR`** (default `build/`) for the
   executable and `sys.path`, as `stage_7_check.py` does after 0007 review 3.
3. **Report.** Append "Review 4 applied" to `report.md`: the two commits; the four chain1d
   `lowdin/exact` runs at `OMP_NUM_THREADS=1` and `16` through `stage6_nonorth` with `cmp`
   (IDENTICAL expected) and their `golden_diff` worst values against G1 and G5; the
   `stage_6_check.py` output. On the fat node (`~/lsquant-agi-0006-build`), foreground
   polling, bounded loops. Do not push.

Noted for the clean-up after the final validation, not for this review: defaults of the
appliers in the header (4096, 1e-13, 200, 1e-12), the fifth applier and the undefined
`lowdin/ischo`, ISCHO's last-call iteration state, the window padding in the app.

Closure after this: the tester reruns the stage-6 part (`6only`, the four rows at both
thread counts) on the fixed head, then `handoff/0006-codex` and `handoff/0006-agi` merge
into `cpp-unified` with stage 7.
