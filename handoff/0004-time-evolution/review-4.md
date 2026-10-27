# Handoff 0004 — review 4 (Fable): stage 4 closed

Developer branch `handoff/0004-agi` at `aba179a` (reviews 1 to 3 applied); tester branch
`handoff/0004-codex` at `745f5a3`.

Gate on TCN from a fresh build: **71 of 71 byte-identical** at 16 and 1 threads
(stages 1 to 3: 41 including G4 and G5; G4c off-centre window 6; stage 4: G6 moments and
the eight curves, G7 `SX`, `SY`, `SZ`: 24), the four greps empty, the self-check all
`array_equal`. Tutorials 02, 04, 05 from Python: 193 of 193, 3 of 3, 16 of 16 IDENTICAL,
figure `figures/out/stage_4/fig_stage_4.png` with every residual zero.

Found and resolved in this handoff:
- tutorials 02 and 05 called a non-existent generator function (review 2);
- `TimeEvolutionStep` renamed `AxisStep`, physics words out of the engine (review 2);
- `tutorial04` reads the golden operators, not the legacy tree (review 2);
- the stage-3 energy axis for non-centred windows, `(x - shift)/scale` (review 3), now
  gated by G4c (graphene, `BOUNDS -8.5 7.5`); the legacy uniform-grid reconstruction has
  the same defect (D2, AGENTS.md section 10), so the G4c curves are blessed from the new
  code with the legacy moments and the legacy curve kept as evidence;
- the library default no longer depends on `compat/` (review 1).

Carried forward: the thread-count slowdown of the time-evolution path (unconditional
parallel SpMV on small systems, up to 100x at 16 threads under load) joins F9 in stage 8;
the extra `.chebmom` written by `stage4_time` and the unused final recursion step are
cosmetic; `~/lsquant/build` on TCN is not the golden build (AGENTS.md section 10).

Stage 4 closed. Stage 5 (`handoff/0005-agi`, branched from `aba179a`) is in review and test.
