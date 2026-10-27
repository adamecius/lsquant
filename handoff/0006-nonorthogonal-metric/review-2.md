# Handoff 0006 — review 2 (Fable to Agi)

Tester run on `handoff/0006-codex` at `98f7807` (fresh build of `a5108e3`): the physics
passes. Fixtures 4 of 4 accepted independently (cond(S) 5.2 to 6.7, spectra invariant to
1.5e-14, sha256 as in the manifests); **40 of 40 rows against G1, G2, G5** PASS (Löwdin
exact 8e-13, Löwdin Chebyshev 1e-11, generalised exact 1e-14, ISCHO 4e-13, Chebyshev
3e-12); routes against each other 24 of 24; Hall plateau `C = +1.017017` on all ten
route/applier pairs, equal to the golden; stages 1 to 3 in the same run 41 of 41
byte-identical; the naive estimator off by O(1) on Haldane, as `theory.md` predicts.

The harness still fails on one line: the stage-1 grep for `basis|std::cout|getenv` under
`lsquant/core` hits nine Doxygen comment lines in `metric_policy.hpp` and
`composite_operators.hpp`. The grep stays as it is (the rule is "no basis word in the
engine", comments included, as for the physics words in stage 4).

Three items; nothing else is to be touched.

1. **Reword the nine comment lines** so that `basis` does not appear under `lsquant/core`:
   say "metric", "overlap `S`", "orthonormal set" or "coefficient representation" instead.
   Verify with the same grep as the harness: `grep -rn "basis\|std::cout\|getenv" lsquant/core`
   must print nothing.

2. **`--route naive` must be the brief's naive**, i.e. the legacy estimator: recursion on
   `S^-1 H~` (the `GeneralisedRescaledOperator`) with `S^-1 V'` on both sides and no `S`
   on the bra (the symmetric stage-1 plan with `MetricInverseProduct(S^-1, V')` as the
   operator). The current code recurses on `H'` alone and ignores `S`, so chain1d
   diverges instead of showing the real defect; the evidence row must show the DOS
   passing and the Kubo moments off by O(1). Keep the name `naive`, state in the help
   text what it is and that it is wrong (`theory.md` section 4.1).

3. **No unused or hidden parameters.** Remove `--drop-tol` (parsed, never used). Make the
   dense-inverse compression threshold a named constant in `overlap_appliers.hpp`
   (`kExactInverseDropThreshold = 1e-18`) with the sentence from the report (below
   double rounding, so the exact appliers stay exact). Add `--coeff-tol` and `--max-dim`
   to `stage6_nonorth` so that every applier parameter the gate relies on is on the
   command line and printed in the stdout line; the header defaults stay. State in the
   `IschoInverseOverlap` comment that `LastIterations()`/`LastError()` describe the last
   solve and that the applier is not safe for concurrent `Multiply` on one object.

Carried to closure by Fable: `ExactInverseSqrtOverlap` as the fifth applier is accepted
(the Löwdin route needs an exact `S^-1/2` reference; `lowdin/ischo` stays N/A); the
Lanczos padding of the overlap bounds in the app is accepted.

Procedure: one commit per item on `handoff/0006-agi`; rsync to `~/lsquant-agi-0006-build`
on TCN (`ssh jh-x`), rebuild, rerun `stage_6_check.py`, the G5 byte check and one naive
run on chain1d and Haldane (paste the stdout lines), commit. Foreground polling, no
background tasks. Do not push.
