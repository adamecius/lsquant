# Handoff 0006 — review 1 (Fable to Agi)

Branch `handoff/0006-agi` at `dd9e908` (five commits), published to origin. Read against
`brief.md`, `theory.md` and the report.

What is right and stays: `MetricPolicy` with the two policies; the five appliers (dense
`LLT` and dense eigensolver references, ISCHO-preconditioned CG with the explicit `L^+`
solve, the two Chebyshev-Gauss series with adaptive truncation, coefficients computed in
C++ as the one allowed exception); `LowdinOperator` and `MetricInverseProduct`;
`RescaledOperator` over a `LinearOperator`; `GeneralisedRescaledOperator` holding its
own CSR copy of `a H' + b S'`; the generalised plan with the adjoint bra
(`S^-1`, `S^-1 V'`, `Expand`, `S`) exactly as `theory.md` 2c; `MomentEngine` and `Plan`
untouched; the orthogonal G5 `kuboxy` byte-identity recorded; the self-check output
(generalised/exact at 4e-16 from G1, Löwdin/exact at 6e-13, the applier tables). The
dense inverses stored as compressed CSR below 1e-18 are accepted as the reference
appliers (the dropped entries are below double rounding; say so in the class comment,
which already does). The branch merges cleanly into `cpp-unified`.

Five items before the tester runs the gate. Nothing else is to be touched.

1. **Merge `origin/cpp-unified` (now `552de6e`: stages 2 and 3 closed) into
   `handoff/0006-agi` and add the binding.** The brief lists `lsquant/api/bindings.cpp`
   as a deliverable; the branch was cut from a stage-1 base, which explains the omission
   but does not remove it (Python is the library). Bind `MetricPolicy`,
   `OrthogonalMetric`, `NonOrthogonalMetric`, `ExactInverseOverlap`,
   `ExactInverseSqrtOverlap`, `IschoInverseOverlap` (with `last_iterations`,
   `last_error`), `ChebyshevInverseOverlap`, `ChebyshevInverseSqrtOverlap`
   (`effective_order`), `LowdinOperator`, `MetricInverseProduct`,
   `GeneralisedRescaledOperator`, and the `LinearOperator` overload of
   `RescaledOperator` (a second `py::init`), with `keep_alive` on every reference the
   C++ side stores, following the existing style of `bindings.cpp`; update
   `python/lsquant/_core.pyi`; add `python/lsquant/models/__init__.py` so
   `lsquant.models.nonorthogonal` is a regular package module like `lsquant.apps`.
   Rebuild on TCN with the binding and run the stage-1 to stage-3 self-checks once to
   show nothing moved.

2. **Fixtures: the condition number is too easy.** The brief fixes `cond(S')` between 5
   and 10; the fixtures have 1.8 to 1.9 (`||alpha K||_2 <= 0.35` gives at most
   `cond(S') ~ 4`). Raise the norm target so that `cond(S')` lands in [5, 10] for all
   four variants (around `||alpha K||_2 = 0.5` to `0.55`; keep `R = I + alpha K` on the
   graph of `H`, seed 12345), refuse outside that band, regenerate the four fixtures,
   and add to each `manifest.json` the spectral window of `H` used by the goldens
   (chain1d `[-2, 2]`, Haldane `[-4.5, 4.5]`) and the sha256 of every `.CSR` written, as
   the brief asks. Rerun `stage_6_check.py` on the new fixtures; the Chebyshev order
   needed at `coeffTol 1e-12` will rise, and that is the point of the test.

3. **`stage6_nonorth` stdout line.** The brief asks for `applier=<name> order=<P>
   tol=<t> max_iterations=<n> residual=<r>` on every run; the current line omits all
   but the name. Print, after the existing fields: `order` (the ctor order for
   chebyshev, `EffectiveOrder()` as `effective_order`), `tol`, `max_iterations`, and
   `residual` (ischo: `LastError()` of the last apply with `iterations=LastIterations()`;
   chebyshev and exact: the residual `||S (S^-1 v) - v|| / ||v||` on one random vector
   at the end of the run). The tester reads these for the figure.

4. **`IschoInverseOverlap::Multiply` must not return an unconverged solution
   silently.** When the iteration cap is reached with `relErr > tolerance`, or when
   `p^+ S p <= 0`, throw `std::runtime_error` naming the residual and the cap. The
   caller raises `maxIterations` or loosens `tolerance`; the recursion must never be
   fed a vector whose accuracy is unknown. Keep `LastIterations()` and `LastError()`.

5. **Report.** Update sections 2, 4 and 5 for the new fixtures and the new stdout
   line; add the binding to section 1; record in section 6 that the branch now
   contains the stage-2 and stage-3 merge. Keep the G5 byte-identity check (section 3)
   and rerun it after the merge, since `RescaledOperator` changed internally
   (`std::optional`): the orthogonal path must still be byte-identical.

Procedure: one commit per item on `handoff/0006-agi` (the merge commit counts as part
of item 1); rsync to `~/lsquant-agi-0006-build` on TCN via `ssh jh-x`, rebuild with the
binding, rerun `stage_6_check.py` and the G5 check, paste the outputs into the report,
commit. Do not push.
