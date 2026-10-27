# Handoff 0006b — review 4 (Fable to Codex): the cost of a material run

Tester round 2 on `99cbda7`, Aron: the fixed sampler works on Bi (16 eigenbasis vectors:
1.66 +- 0.03 against the oracle 1.5527; round 1 gave +7.07), but a Bi spin Hall run
(48 x 48, W = 72, M = 1024, 16 vectors) takes 72 minutes and PdSe2 (W = 220) more than
40 h per exact pass. Profile and timings (short Bi run, M = 128):

| hot spot | time | where |
|---|---|---|
| rank-2 contraction on one core | ~270 s per vector at M = 1024 | `moment_engine.cpp:541-558` (review 3 item 3) |
| oracle `sigma_xy`: erfc evaluated by a Python generator | 82.6 s | `oracle/spin_hall.py:91` |
| `bloch_sum`: serial Python loop over k | 48.9 s (60 % of the short run) | `models/bloch.py:36` |
| reference pass repeated in the eigen basis with exact trace | doubles every run | the app |

Rule for all of it (user): no ninja programming. Plain vectorisation, plain OpenMP with
static schedules, determinism across thread counts kept.

## Items, one commit each, after the review-3 commits

1. **`bloch_sum` vectorised.** One NumPy operation over all k and all R at once (the phase
   matrix `exp(i k.R)` of shape `(N_k, N_R)` times the stacked `(N_R, W, W)` blocks, as
   one `einsum`/`tensordot`, chunked over k so that memory stays below a run-file limit),
   or the stage-5 C++ Bloch transform. Same numbers to round-off (the tester compares).
2. **Oracle `sigma_xy` vectorised.** The broadened occupation is `0.5 erfc(...)` evaluated
   for all bands and energies as one `scipy.special.erfc` call on an array; no Python
   generator in any per-k or per-energy loop of the oracle.
3. **No redundant reference.** When the trace is exact in the eigen basis (or the random
   phase is drawn in the eigen basis), production and reference compute the same thing:
   the reference pass is skipped and the run's receipt says so; `"reference": true` in the
   run file forces it.
4. **Report.** Append "Review 4 applied" to `report.md` with what the tester must time.
