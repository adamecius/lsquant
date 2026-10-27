# Handoff 0005 — review 1 (Fable to Agi)

Branch `handoff/0005-agi` at `7109c20` (four commits), published to origin; the tester
starts the gate on it in parallel.

What is right and stays: `BlockDiagonalOperator` (batched block GEMV, Hermiticity cached),
`BlochTransform` with `fftw_plan_many_dft` and the atom-gauge phase table, with the sign
and normalisation conventions derived in the report from `B H_csr B^+ = H(k)` and checked
to 7e-16 against the wannier2sparse `VX`; `SandwichOperator` and `DiagonalOperator`;
`BlockLocalStates` and `TransformedStates` (the same random vectors in both bases, which is
why real and Bloch moments agree to 1e-15: good choice); the Wannier reader with `ndegen`
and `wsvec` folded once; the three bases and the app; the self-check numbers; G1, G4, G7
byte-identical. The eigen path run as KPM on the diagonalised blocks with the exact trace
is accepted as the "same kernel" comparison, and the Jackson-vs-Gaussian 2 percent number
is the right cross-check against the independent oracle.

Two items; nothing else is to be touched.

1. **The window of a Bloch or eigen run must come from the eigenvalues, not from
   Gershgorin** (brief section 4). All three `window()` methods call
   `window_from_gershgorin`. Implement `EigenBasis.window(pad)` from
   `min/max` of the block eigenvalues padded by `pad` of the width, let
   `BlochBasis.window` return the same (it owns the blocks: diagonalise once, cache), and
   keep `RealBasis.window` on Gershgorin only as the stage-1-compatible default. In
   `stage5_bases` add `--window eigen|gershgorin|<lo> <hi>` with default `eigen` for every
   basis, so that the three bases of one comparison share the same window (the app
   computes it once from the blocks and passes it to all three); say in the docstring
   why Gershgorin is not used on Wannier models (several times wider than the spectrum).
   Rerun the self-check: check (e) numbers change (tighter window, same agreement), (f)
   unchanged.

2. **No heuristic scaling in `BlochTransform`.** The constructor decides whether the
   phases "already carry `1/sqrt(N_k)`" by testing `|phase_0|^2 N_k ~ 1`. Remove that:
   the constructor takes unscaled phases `exp(-i k . tau)` (as `wannier.bloch_phases`
   returns them) and always multiplies by `1/sqrt(N_k)`; the `shared_ptr` constructor
   takes scaled phases and says so in its name (`FromScaledPhases` factory or a tag
   type). State in the class comment that `Multiply` uses a per-object FFTW workspace and
   is not safe to call concurrently on the same object (the engine never does; copies get
   their own workspace).

Note for the merge, no action: `RescaledOperator(const LinearOperator&, ...)` and the
`composite_operators` additions also exist on `handoff/0006-agi`; whichever branch merges
second will be reconciled by Fable.

Procedure: one commit per item on `handoff/0005-agi`; rsync to `~/lsquant-agi-build` on
TCN (`ssh jh-x`), rebuild, rerun `stage_5_check.py` and the CLI suite with the new default
window (foreground polling, no background tasks), paste into `report.md`, commit. Do not
push.
