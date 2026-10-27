# Handoff 0006 — review 3 (Fable): stage 6 closed on the developer branch; final gate deferred

Developer branch `handoff/0006-agi` at `8305f11` (reviews 1 and 2 applied; merged with
stages 4 and 5, all six self-checks and the orthogonal G5 byte check passing after the
merge). Tester branch `handoff/0006-codex` at `98f7807` (run 1 on `a5108e3`): fixtures
4 of 4 accepted independently (cond(S) 5.2 to 6.7, spectra invariant to 1.5e-14,
sha256 as in the manifests); 40 of 40 rows against G1, G2, G5 PASS (Löwdin exact 8e-13,
Chebyshev 1e-11; generalised exact 1e-14, ISCHO 4e-13, Chebyshev 3e-12); routes against
each other 24 of 24; Hall plateau `C = +1.017017` on all ten route/applier pairs, equal to
the golden; the naive legacy estimator off by O(1) on Haldane (evidence row, as
`theory.md` section 4.1 predicts); stages 1 to 3 in the same run 41 of 41 byte-identical.

Decision (user, 2026-10-21 13:00): the project closes stage 7 first and runs the long
gates once at the end. Stage 6 is therefore merged now on the strength of run 1 and the
post-merge self-checks; the tester's closure run on `8305f11` (in progress: the comment
grep, the corrected naive route, the explicit applier parameters) and its harness,
figure and `tests.md` are merged when it reports, and repeated in the final validation
run together with the deferred stage-7 goldens.

Review-2 items (verified by the developer on TCN, to be confirmed by the closure run):
comments reworded (no `basis` under `core`); `--route naive` is the legacy estimator
(DOS equals G1 to 1e-15, Kubo off by 0.17 to 0.20); `--drop-tol` removed, the dense
compression threshold is a named constant, `--coeff-tol` and `--max-dim` explicit.

Carried forward: `ExactInverseSqrtOverlap` as the fifth applier; the Lanczos padding of
the overlap bounds in the app; the design decision of the notes in
`handoff/0006b-siesta-spin-hall/` (exact `S^-1` per block by default, the certified
single-shot `X` for real-space `S`, study 0006c).
