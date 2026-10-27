# Handoff 0004 — review 1 (Fable to Agi)

Branch `handoff/0004-agi` at `be54b79` (five commits), published to origin; the tester is
running the gate on it in parallel.

What is right and stays: the seed registers and the generalised `MirrorToBra`; `ApplySeries`
accumulated forward with the recursion advancing (not Clenshaw) and the global phase as a
final `scal` carried by the `scale` field (accepted deviation 2: the coefficients still come
from outside, the engine applies one scalar); `ApplyDerivedSeries` in the `Conv_Evolve`
order; `TimeAxis` as a trailing tensor dimension with the body run from copies of the
seeds; the `(2-delta)/N` weights replicated over the time slices; GSL confined to
`compat/` and linked only by `_core`; the `.chebmomTD` writer and the reader detection;
the three plans in Python; self-check (a) to (d) all `array_equal`/zero differing lines;
the G7 bounds source (Gershgorin branch, `BandWidth = 4.62`) and the
`RescaleCon2ChebyshevDomain` note, which Fable adds to AGENTS.md section 10 at closure.
Deviation 1 (no `ZeroRegister`) is accepted: the empty sum is the zero map, which is
mathematically consistent, but it must be stated.

Two items; nothing else is to be touched.

1. **The library must not depend on `compat/` by default.** `time_evolution_series`
   defaults to `use_legacy_bessel=True`, so every library caller goes through the GSL
   path that stage 7 deletes. Make the default `False` (the sympy/scipy
   `TimeEvolution` coefficients), and have the gate app `stage4_time.py` and the three
   tutorial apps pass `use_legacy_bessel=True` explicitly with a one-line comment that
   byte identity with the legacy captures requires the GSL Bessel values (the same
   pattern as the text round trip in tutorial 03). Rerun the self-check: (a), (b), (d)
   unchanged because the apps keep the legacy path; (c) unchanged.

2. **State the empty-series semantics where the instruction is declared.** In
   `plan.hpp`, the `ApplySeries` and `ApplyDerivedSeries` comments must say that an empty
   `coeffs` is the zero map (`v <- 0`), that this is how a seed register is cleared, and
   that `msd_plan` relies on it; the same sentence in the `msd_plan` docstring. No code
   change.

Procedure: one commit per item on `handoff/0004-agi`; rsync to `~/lsquant-agi-build` on TCN
(`ssh jh-x`; TCN is loaded, run the self-check with nohup and poll), paste the self-check
output into `report.md`, commit. Do not push.
