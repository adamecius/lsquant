# Handoff 0006b — review 3 (Fable to Codex): round 2 on `99cbda7`

Tester round 2 (Aron, interim): review-2 fix confirmed. Assertion (g) passes (one-vector
sample = NumPy quadratic form, residual 5.6e-16); the 20-seed eigenbasis random-phase
estimate on the orthogonal Kane-Mele 60 x 60 gives 2.002592650725 with std 1e-12, equal
to the exact KPM (the user's statement that eigenbasis random phase equals the exact trace
holds once the sampler is right); Gamma Hermiticity failure gone. PdSe2 through SIESTA:
trivial gap at E_F (spin Chern 0, sigma 0.006 / 0.011), topological gap at +1.01 to
+1.09 eV (spin Chern 2, Z2 1), plateau 2.08 (berry_connection) / 1.98 (covariant),
converged at 24 x 24 to 0.003. Bi unchanged: 1.5527 / 1.1747.

## Items, one commit each, on `handoff/0006b-dev` after `99cbda7`

1. **Oracle degeneracy tolerance.** `oracle/spin_hall.py:38,52` drops a pair only if
   `|E_n - E_m| <= 1e-10` eV absolute. With inversion symmetry (PdSe2) Kramers pairs are
   split by round-off (1e-10 to 1e-6 eV), their `1/dE^2` reaches 1e20 and the result is
   thread-dependent garbage (-25.6, -14.3, +1.39 at 12, 24, 48; reruns differ); the old
   "1.40 -> 1.70" was this. Treat near-degenerate pairs (within a tolerance tied to the
   eigensolver residual and the bandwidth, default 1e-6 eV, a run-file field) as one
   degenerate subspace: their intra-subspace terms are excluded, as for exact degeneracy.
   Report in the oracle output the tolerance used and the number of pairs it merged. The
   tester found the values stable to 5 digits from 1e-6 to 1e-3 eV.
2. **Quadrature certificate target.** The factor-10 margin of review 2 cannot be met at
   the default tolerance 1e-5: the overlap residual against the stored SIESTA `S` levels
   off near 1e-6 (the stored `S` itself carries that precision) while the refinement is
   converged to 5e-8. Separate the two criteria: refinement convergence with margin 10
   (certifies the integral), and agreement with the stored `S` against a tolerance that
   is the stored precision (report it, refuse only above 1e-5). The committed covariant
   run files must run as committed.
3. **The rank-2 contraction is serial.** In the Bi exact-trace run the moments phase runs
   on one core: the `M x M` loop of `linalg::dot` at `lsquant/core/chebyshev/moment_engine.cpp:541-558`.
   Parallelise it with a static schedule over the outer moment index (deterministic:
   each `(m, n)` element is computed by one thread in a fixed order), or call the
   batched dot as one GEMM where the vectors are stored contiguously. Thread-count byte
   identity must hold (the tester checks 1 vs 16 threads on G5).
4. **Report.** Append "Review 3 applied" to `report.md`. Do not push.
