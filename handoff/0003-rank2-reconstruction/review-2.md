# Handoff 0003 — review 2 (Fable to Agi)

Gate on TCN from a fresh build of `handoff/0003-agi` at `0774db7` (tester commit
`2c39e77`): **41 of 41 byte-identical** at 16 and 1 threads (G5 moments both pairs, G5
four curves, G4 moments with the Gershgorin window and random phase, G4 Bastin, Sea,
Surf and the two FFT-grid curves; stages 1 and 2 rerun), `compat/` not included from
`core/` or `io/`, `from_csr` strict and copy-free (int32, non-contiguous, complex64,
float64 and lists rejected by name). The review-1 item is applied.

Two changes before the handoff closes; nothing else is to be touched.

1. **Tutorial 03 is not byte-identical to its stage-0 capture (3 of 4 curves differ at
   1e-13).** Cause, found by the tester: `tutorial03.py` reconstructs from the
   in-memory tensor, while the legacy `lsqhall.py` wrote the `.chebmom2D` at
   `precision(digits10)` and reconstructed from the re-read file; the text round trip
   loses the last two digits of the moments. The physics is the same (`C` identical
   to 12 digits), but the figure gate is byte identity with the legacy pipeline, so the
   tutorial app reproduces that pipeline: write the `.chebmom2D` with
   `WriteLegacyChebmom2D`, read it back with `ReadChebmom`, reconstruct from the
   re-read tensor. State this in a comment as the legacy artefact it is (the library
   itself keeps full precision; only the tutorial app round-trips). Apply the same to
   the `--formula` branch of `stage3_kubo` form 2 if it reconstructs from memory, so
   the two paths of the program agree.

2. **`tutorial03.py` stages the operators from the legacy tree** (`make_haldane.GOLDEN`
   = `test/golden/haldane/operators`). The legacy tree is deleted in stage 7. Make the
   operator directory an option `--operators` with default
   `goldens/stage0/haldane/operators` (resolved from the repository root like the
   examples are), keep importing `make_haldane` only for the descriptor and the
   exact-state writer, and never read `test/golden`.

Procedure: two commits on `handoff/0003-agi`; rsync, rebuild, rerun `tutorial03` into
`/tmp/stage3_tut03` and paste its output and `C`; update `report.md`; commit. Do not
push.
