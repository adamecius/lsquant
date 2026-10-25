# Handoff 0006b — review 2 (Fable to Codex): the sampling path of the block basis has a variance defect

Evidence (tester, Aron, `tests.md` of this handoff): on the same orthogonal Kane-Mele
model, 60 x 60 mesh, `M = 512`, 30 meV, 16 random-phase vectors, exact value 2.0055:

| estimator | spread over seeds |
|---|---|
| stage-3 real-space random phase on the wannier2sparse CSR (`apps/stage3_kubo`) | 2.0014, 2.0047, 2.0041, ... (about 0.003) |
| 6b block-basis sampler, eigenbasis phases (`bind_eigen_trace`) | std 11.2 (10 seeds) |
| 6b block-basis sampler, orbital-Bloch phases | std 10.9 (10 seeds) |

Both 6b samplers are unbiased (100-seed means within the standard error) and show the same
spread on the rotated and on the orthogonal model, so neither the metric nor the phase
basis is the cause. And the spread is impossible for a correct quadratic form: for a
random-phase vector `r` the variance of `r^+ M r` is `||M||_F^2 - sum_i |M_ii|^2`; `||M||_F`
is basis independent, `M` is block-diagonal in `k`, and each real-space diagonal element is
the `k`-average of the `k`-basis diagonal, so by Jensen the `k`-basis variance is at most
the real-space one, never 3000 times larger. Hence each 6b sample is not a single
quadratic form `r^+ M r`: most likely the bra and the ket of the rank-2 plan do not carry
the same random vector (the `MirrorToBra` adapter that "prepends C on the ket and S on the
copied bra", or `RandomPhaseStates` over the block layout regenerating or correlating
phases across blocks, or the normalisation). Consequence seen on Bi: 16-vector KPM
plateaus of +7.07 and -5.21 against oracle values 1.55 and 1.17, while production and
eigen reference agree to 1e-11 (the engine does what the sampler asks).

This also voids the measurement behind review-1 item 6: the noise of a correct eigenbasis
or orbital-basis estimator is unknown until the sampler is fixed. Item 6(a) stands (exact
block trace for tests and goldens); 6(b) is re-measured after this review.

## Items, one commit each, on `handoff/0006b-dev` after the review-1 commits

1. **Find it with the cheap test first.** On a 4 x 4 mesh at small `M`, for one random
   vector, compute the sampled word `r^+ A f(H) B g(H) r` directly in numpy in the eigen
   (or Bloch) basis from the same `r` the engine receives, and compare with the engine's
   sample (`EigenMetricBasis` / `BlochMetricBasis` rank-2 plan, one vector). They must
   agree to round-off. Put this comparison into `tests/dev/stage_6b_check.py` as assertion
   (g): one-vector sample equals the numpy quadratic form (1e-10), on the rotated and the
   orthogonal fixture, generalised and Löwdin routes.
2. **Fix the sampler** so that bra and ket are the same random vector in every plan
   (`MirrorToBra` semantics of AGENTS.md section 2: the bra side is built from a fresh copy
   of the same `|r>`), with the metric inserted as theory.md equations 2a to 2c prescribe
   (`S` on the bra copy is right only if the ket copy is `C r` and the bra copy is the same
   `C r`), and phases drawn once per sample over the whole block layout.
3. **Re-measure** in the self-check: the std over 20 seeds of the 16-vector estimate on
   the orthogonal Kane-Mele, 60 x 60, `M = 512`, 30 meV, for eigenbasis and orbital-Bloch
   phases, printed next to the stage-3 real-space value (expected: at most the real-space
   spread, about 0.003, by the Jensen argument). The number decides review-1 item 6(b).
4. **Hermiticity of `C^+ O C`.** `EigenMetricBasis.operator` builds `C^+ J C` without a
   Hermitian fold and `BlockDiagonalOperator` tests Hermiticity with an absolute 1e-12
   (`block_diagonal_operator.cpp:37`); on Bi the elements reach hundreds of eV A and
   round-off alone fails it, killing the eigen-reference pass after the production pass.
   Make the block test relative to the operator's largest block norm (same rule as
   review-1 item 1), and fold `C^+ O C` Hermitian explicitly when `O` is Hermitian
   (reporting the anti-Hermitian norm it removed).
5. **Report.** Append "Review 2 applied" to `report.md`: the localisation (file, line, what
   the bra carried before), the fix, the self-check assertion (g), and the re-measured
   std that the tester must confirm.

## Addendum from the tester's final report (`handoff/0006b-codex` at `78e0431`, `tests.md`)

- Stage-3 real-space estimator, 10 seeds on the 60 x 60 Kane-Mele: mean 2.0027, spread
  0.0041; seed 12345 gives 1.9957, which reproduces the day-20 number. The 6b samplers
  spread 11 on the same problem: 2700 times more. Item 1 and 2 stand as written; a further
  place to look is `ApplyKpmWeights(mu, nstates, False)` in `apps/stage6b_siesta.py`.
- Profile of the short Bi run (`M = 128`, 2 vectors, 48 threads, wall 86 s): 60 % of the
  time (48.7 s) in the serial Python per-k loop `bloch_sum` (`python/lsquant/models/bloch.py:36`);
  the engine call at most 21 %; `eigh` and the metric 7 %. CPU 567 % of 4800 % during the
  Python phase and 940 % during the engine phase. So review-1 item 5 has its target:
  `bloch_sum` vectorised over k (one `einsum`/FFT over the R list, or the C++ Bloch
  transform of stage 5), the eigen reference optional, and the block engine's own
  occupancy (940 % of 4800 %) goes to stage 8 item 5(d).
- The dipole quadrature certificate has no margin: refinement error 9.87e-6 against the
  1e-5 tolerance reached only at the maximum order 96. Item 6 of this review: raise the
  maximum order (or use a finer radial grid) so the certificate passes with a factor 10 of
  margin, and report the cost; and cache the certified `A(R)` per SIESTA run on disk so the
  256 s are paid once, not once per run file.
- The Z2 convention deviates from the brief as the developer's report declares: state the
  convention used (`Z2 = C_s mod 2` with the Prodan sectors) in the oracle docstring and
  keep it; the brief is amended.
