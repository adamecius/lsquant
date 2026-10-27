# Handoff 0001 — review 3 (Fable): stage 1 closed

Verdict: **closed**, merged into `cpp-unified`.

Gate on TCN from a fresh build of `handoff/0001-agi` at `34173e2` (tester commit
`010555a`): G1 and G2 byte-identical at 1, 8 and 64 threads; the batched run
(`--batch 16`) byte-identical to the unbatched G2, which the legacy could not do;
random-phase runs reproducible; `grep` clean; legacy `ctest` 31/31 with the new core
on; `stage_1_check` 2.5e-15 (dos) and 1.1e-13 (kubo) against dense algebra, bit-identical
across threads through the parallel matvec; stage-1 figure with residual exactly zero.

What this stage established: the five-instruction engine (`Plan`) reproduces the legacy
rank-1 and rank-2 moments to the byte, including the legacy's signed zeros and its
accumulation order, and does so with a basis-free `LinearOperator`, a rank-generic
tensor, and an explicit Hermitian fold confined to the CSR reader.

Carried forward:

- F9 (tester): with the parallel threshold removed, small systems pay the fork/join cost
  at 64 threads (`kubo` 30 s at 8 threads, 400 s at 64 on a loaded machine). Correctness
  is unaffected. Stage 8 gives `SparseMatrix::Multiply` an explicit parallel threshold
  as a constructor parameter (default 0 = always parallel, so the gate keeps exercising
  the parallel path) and the gate runs at 1, 8 and 16 threads.
- The four hand-written rank paths in `MomentEngine::Run` become one generic loop when
  `ApplySeries` arrives (stage 4).
- The legacy batched-moments defect, the flaky `graphene_regression` and the bad graphene
  `.desc` are recorded in `AGENTS.md` section 10.

Conduct note, closed: the review-1 symlink episode is declared in `report.md`; the
rule stands that a report describes what the committed code does from a fresh build.
