# Handoff 0005 — review 3 (Fable to Agi)

Tester closure run on `handoff/0005-codex` at `d2f739c` (fresh build of `dd989b0`): **155 of
155 rows PASS**: real vs Bloch moments to 2e-14 with and without disorder, Bloch vs eigen
curves to 5e-15, eigen exact (Gaussian) vs the independent oracle on the Gamma-centred
mesh to 5e-15, the FFTW `B` against a dense `B` to 3e-15, the real basis equal to the
wannier2sparse CSR, every stage-5 output byte-identical between 1 and 16 threads; stages 1
to 4: 71 of 71 byte-identical. Review-2 items verified.

Two items before the handoff closes; nothing else is to be touched.

1. **No silent fallback when loading the wannier2sparse CSR.** `RealBasis._try_load_w2s_csr`
   swallows every exception (`except Exception: pass`) and accepts any CSR whose dimension
   matches even when the `.out` receipt is missing or unreadable, so a stale file of the
   right size would be used without a trace. Make it strict: if `--csr-dir` (or the model
   directory) holds `<seed>.HAM.CSR`, it must also hold a readable receipt whose label,
   supercell and operator list match, or the load raises with the reason; no exception is
   caught and discarded. Record in `<stem>.basis.json` where every operator came from
   (`"source": "wannier2sparse <path> <sha256 of the CSR>"` or `"built by lsquant"`).
   Add to the self-check one negative case: a CSR of the right size with a wrong receipt
   must raise.

2. **The basis in the output names.** The stem lost the basis name, so the three bases
   written to one directory overwrite each other. Put it back as in the first version
   (`<seed>_<basis>` as the label part of the legacy stem, e.g.
   `SpectralOp1graphene_bloch_KPM_M256_stateexact.chebmom1D`), for the moment files, the
   curves and `basis.json`; say in the app docstring that the legacy stem is kept
   otherwise.

Carried to closure by Fable: the `RescaledOperator` constructor from a `LinearOperator`
(also on `handoff/0006-agi`) is accepted as the stage-5/6 extension of the engine that the
Bloch and Löwdin paths need; the thread-safety limitation is documented and matches the
engine's single-vector use; curve files keep the legacy 6-digit format, checks go through
the moments.

Procedure: one commit per item on `handoff/0005-agi`; rsync to `~/lsquant-agi-build` on
TCN (`ssh jh-x`), rerun the self-check and the app suite (foreground polling), paste into
`report.md`, commit. Do not push.
