# Handoff 0005 — review 4 (Fable): stage 5 closed

Developer branch `handoff/0005-agi` at `642955b` (reviews 1 to 3 applied); tester branch
`handoff/0005-codex` at `274302a`.

Gate on TCN: **155 of 155 rows PASS** (real vs Bloch moments to 2e-14 with and without
Anderson disorder; Bloch exact trace vs eigen curves to 5e-15; eigen exact Gaussian vs the
independent k-mesh oracle `hr_exactdiag` on the Gamma-centred mesh to 5e-15 for graphene
and Haldane, DOS and `sigma_xx`; the FFTW `B` against a dense `B` to 3e-15; the real basis
equal to the wannier2sparse CSR; every stage-5 output byte-identical between 1 and 16
threads), stages 1 to 4 **71 of 71** byte-identical, greps empty. G9 committed under
`goldens/stage0/bases/` with provenance (`scripts/hr_exactdiag_gamma.py`, Gaussian
40 meV, area factor stated). Figure `figures/out/stage_5/fig_stage_5.png`.

Found and resolved in this handoff: the brief's cell order and sign of `B` (corrected in
the brief from the tester's run 1); the Gershgorin window on Wannier models (eigenvalue
window, app default); the wannier2sparse CSR as the real basis with strict receipt
validation and provenance in `basis.json`; sparse `_hr.dat` and single-header `.xyz`
readers; the exact Gaussian evaluation in the eigenbasis as the oracle of the tables.

Carried forward: `RescaledOperator` over a `LinearOperator` is the agreed stage-5/6
extension (also on `handoff/0006-agi`, reconciled at the stage-6 merge);
`SandwichOperator` and `BlochTransform` are documented as not safe for concurrent
`Multiply` on one object (the engine applies operators to one vector at a time); the KPM
vs Gaussian kernel shape (1e-2 to 1e-1) and the finite-size rows are reported, not gated.

Stage 5 closed.
