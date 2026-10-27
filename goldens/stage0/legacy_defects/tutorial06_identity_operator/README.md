# Legacy defect D5: tutorial 06 denominators from an operator that was never built (evidence only)

These 142 files are the legacy `goldens/stage0/tutorials/06` files that the new code does not
reproduce, kept with their original relative layout. They are **not references**; no harness
compares against them. The tutorial-06 goldens were re-blessed from the new code (Fable's
decision, handoff 0007 closure; see `PROVENANCE.md`, G8).

## The files

- 139 curves `meanTimeProj1-P*EF2.000000JACKSON.dat` (top level and `captured/`): the
  denominator word `Tr[P U(t)^+ 1 U(t) delta(E-H) P]` of `examples/06_spin_relaxation_square`
  (`lsqrelax.py:57`, `rp_spin.py:51,64`);
- `ey_ratelaw.dat` and `captured/0254_ey_ratelaw.dat`, derived from them;
- `captured/index.tsv`, the stage-0 capture index with the legacy SHA-1 values.

## Cause

`inline_compute-kpm-TimeEvProjetedOp <label> 1 <proj> M Nt tmax exact` builds its operators from
`operators/<label>.<ID>.CSR` (`legacy-final:src/inline_compute-kpm-TimeEvProjetedOp.cpp:53-62`).
The tutorial generators `make_ey.py` and `make_rashba.py` write `HAM`, `SX`, `SY`, `SZ`, `PX` and
`PZ`, but no `<label>.1.CSR`. `SparseMatrixBuilder::BuildOPFromCSRFile` returns without a message
when `Sparse::OPERATOR_FromCSRFile` fails (`legacy-final:include/sparse_matrix.hpp:132-133`), so
`OP[1]` is used without ever being built.

## Effect

The legacy denominators are noise of order 1e-6. An example is
`captured/0001_meanTimeProj1-PZey_L20_D0.05...`, whose first line is `0 9.98496e-07`.

The new code (`python -m lsquant run runs/goldens/tutorial06.json`) applies the identity. Its
curve is constant in t and equal to the numerator at t = 0, which is `Tr[P delta(E-H) P]`
(first line `0 0.0532345`).

The tutorial normalises `s = num / mean(den)` and then `s / s[0]`, so the denominator's scale
cancels. The stage-0 figures were correct, and `ey_ratelaw.dat` agrees to 1.8e-16 relative
(rounding).

All 139 numerator curves are byte-identical between the legacy build and the new code.
