# Handoff 0003 — review 3 (Fable): stage 3 closed

Developer branch `handoff/0003-agi` at `cac5da6`; tester branch `handoff/0003-codex` at
`0d1c457` (merge of `cac5da6` plus outputs and `tests.md`).

Gate on TCN from a fresh build: **41 of 41 byte-identical** at 16 and 1 threads (G5
moments both operator pairs, G5 four curves, G4 moments with the Gershgorin window and
random phase, G4 Bastin, Sea, Surf and the two FFT-grid curves, stages 1 and 2 rerun),
`compat/` not included from `core/` or `io/`, `from_csr` strict and copy-free, the stage-2
code-line grep empty.

Review-2 items verified:

1. Tutorial 03 now round-trips the moments through `WriteLegacyChebmom2D` / `ReadChebmom`
   like the legacy pipeline; the four tutorial curves are IDENTICAL to the stage-0
   capture and `C = +1.017016743089` is the same from the new and the golden curve. The
   `--formula` branch of `stage3_kubo` form 2 (not covered by the gate) was run by the
   tester on Haldane and is IDENTICAL for the four curves.
2. `tutorial03 --operators` defaults to `goldens/stage0/haldane/operators`; the legacy
   `test/` tree is no longer read.

Figure: `figures/out/stage_3/fig_stage_3.png` with caption "4 of 4 IDENTICAL".

Carried forward: the legacy text round trip is a tutorial-app artefact only (AGENTS.md
section 10); the library keeps full precision. The F9 performance note (stage 1 review 3)
stays for stage 8.

Stage 3 closed. Stage 4 brief is in `handoff/0004-time-evolution/`.
