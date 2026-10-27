# Handoff 0004 (Agi) — stage 4: applied series, time axis, mean-square displacement, spin dynamics

Roles: Agi implements, the tester (`codex.md`) runs the gate and the figure on TCN, Fable
reviews. Read `AGENTS.md` (sections 2, 3, 8, 9, 11) and this brief completely. Nothing is
optional; deviations go in `report.md`. Starting point: `origin/cpp-unified` after the
stage-3 merge; branch `handoff/0004-agi`.

## 0. Scope

The "apply" half of the engine: a spectral slot contracted on the fly with a coefficient
series instead of opening an index, a derived series for commutators, and a time axis
on the tensor. With them the three legacy time-dependent programs become plans:
`MeanSquareDisplacement` (G6), `TimeEvolvedProjectedOperator` (G7) and
`TimeDependentCorrelations` (no golden; plan only, tested by the self-check against the
legacy driver output on the chain). Deliverables: G6 and G7 moment files and the G6
curves reproduced **byte for byte** from Python; tutorials 02, 04 and 05 from Python.

Out of scope: optical conductivity, bases, metric, parser, CLI, any change to the
rank-1/rank-2 paths beyond what section 3 names.

## 1. Where things go

```
lsquant/core/chebyshev/plan.hpp                 + ApplyDerivedSeries, TimeAxis, Register, EvolveRegister, EvolveDerived, CopyRegister (section 3)
lsquant/core/chebyshev/moment_engine.cpp .hpp   + the semantics of section 3 (ApplySeries implemented; registers; time axis)
lsquant/core/chebyshev/chebyshev_moments.hpp    + the time axis is an ordinary trailing dimension; metadata: timeStep, numberOfTimeSteps
lsquant/io/chebmom_file.cpp .hpp                + WriteLegacyChebmomTD, ReadChebmom for the .chebmomTD layout (section 4)
lsquant/compat/legacy_kernels.cpp .hpp          + bessel_series_coefficients(x, tol) with GSL gsl_sf_bessel_Jn (section 4) -- gate-only
lsquant/api/bindings.cpp                         + all of the above
python/lsquant/spectral_functions.py            + TimeEvolution (sympy: exp(-I x t)), its Bessel coefficients via scipy (checked vs compat)
python/lsquant/moments.py                        + msd_plan, projected_evolution_plan, time_correlation_plan, time-step helpers
python/lsquant/reconstruct.py                    + time_correlation_curve (section 4)
python/lsquant/apps/stage4_time.py               the gate program (section 6)
python/lsquant/apps/tutorial02.py, tutorial04.py, tutorial05.py   (section 7)
tests/dev/stage_4_check.py                       the single self-check (section 5)
```

## 2. Restrictions

All of stages 1 to 3. In addition: the engine gains instructions and registers, never a
notion of "time", "velocity", "spin" or "projector" (the word `TimeAxis` names a tensor
dimension, nothing physical); GSL is linked only by `compat/` (it is already a legacy
dependency on TCN) and `core/` never includes it; the series coefficients always come
from outside the engine as arrays.

## 3. The instruction set, completed

Registers: every side (ket, bra) has a **seed** vector; a third register `Aux` is
available. `Register = { Ket, Bra, Aux }`.

- `MirrorToBra`: the bra seed becomes a copy of the **current ket seed** (which is `r`
  if nothing has touched it). This generalises the stage-1 semantics (there the ket seed
  was always `r` at that point); stage-1 plans are unchanged.
- `ApplyOperator(k)` before `Expand` acts on the active seed; after `Expand` on the batch
  (as today).
- `ApplySeries{coeffs}` (declared in stage 1, implemented now): active seed
  `v <- sum_{n=0}^{L-1} c_n T_n(H~) v`, accumulated **forward** term by term
  (`y = 0; for n: y += c_n T_n v`, with the recursion advancing), exactly as the legacy
  `Evolve`: not Clenshaw. `L = coeffs.size()`.
- `ApplyDerivedSeries{coeffs, operatorIndex k}`: active seed `v <- sum_n c_n D_n v` where
  `D_0 v = 0`, `D_1 v = O_k T_0 v`, `D_{n+1} v = 2 H~ D_n v - D_{n-1} v + 2 O_k T_n v`,
  i.e. the recursion of the legacy `Conv_Evolve`: this is `[X, f(H~)] v` when
  `O_k = [X, H~]`. Operation order per step exactly as `Conv_Evolve` (recursion of `T`,
  then of `D`, then the `+2 O_k T_n` term, then swap; the `n = 0` term contributes
  nothing).
- `TimeAxis{numSteps, evolve}`: the tensor gains a trailing index of length `numSteps`.
  For `n = 0 .. numSteps-1`: the instructions after `TimeAxis` (the body: `Expand`,
  post-`ApplyOperator`, `Trace`) run from copies of the current seeds and accumulate
  into `mu[..., n]`; then `evolve` (a list of `EvolveRegister{target, coeffs}`,
  `EvolveDerived{target, source, coeffs, operatorIndex}` meaning
  `target += sum_n c_n D_n source`, and `CopyRegister{from, to}`) is applied to the
  seeds in order. `TimeAxis` appears at most once, after all seed preparation.

Legacy plans in this language (`U` = the Bessel series of section 4, `V` = the
rescaled commutator operator, `P` the projector; `w` the stage-1 `(2-delta_m0)/N` weights):

```
MeanSquareDisplacement (G6):           seeds: ket = 0 (MSDWF), aux = r, bra = ket
  ApplyOperator? no. Plan: CopyRegister{Ket -> Aux}; ScaleSeed... -> instead: TimeAxis{N, evolve:[EvolveRegister{Ket, U}, EvolveDerived{Ket += D(U) Aux, V}, EvolveRegister{Aux, U}]}
  body: Expand(M); Trace(w)   with bra = ket seed at each step (MirrorToBra inside the body: bra seed = copy of the ket seed of that step)
  Because the ket seed starts as the zero vector, the first time slice is zero, as in the legacy.
TimeEvolvedProjectedOperator (G7):     ket: ApplyOperator(P); MirrorToBra;   TimeAxis{N, evolve:[EvolveRegister{Ket,U}, EvolveRegister{Bra,U}]}
  body: Expand(M); ApplyOperator(S_a) [post]; Trace(w)
TimeDependentCorrelations (VAC):       ket: r;  MirrorToBra; bra: ApplyOperator(OPL);   TimeAxis{N, evolve:[EvolveRegister{Ket,U}, EvolveRegister{Bra,U}]}
  body: ApplyOperator(OPR) [pre, on a copy of the ket seed]; Expand(M); Trace(w)
```

Read `src/chebyshev_solver.cpp` (`MeanSquareDisplacement`, `TimeEvolvedProjectedOperator`,
`TimeDependentCorrelations`) and `src/chebyshev_momentsTD.cpp` (`Evolve`, `Conv_Evolve`,
`MSD_Evolve`, `ResetGaussian`, `RescaleCon2ChebyshevDomain`) before writing the semantics;
the exact order of `Evolve(MSDWF)`, `Conv_Evolve(PhiT)`, `axpy`, `Evolve(PhiR)` in
`MSD_Evolve` is what `evolve` must reproduce. Design the `MSD` seed initialisation so it
needs no new instruction if `MirrorToBra` inside the body and a zero ket seed are
enough; if a `ZeroRegister` instruction is needed, add it and say so.

## 4. What byte identity with G6 and G7 requires

- Series coefficients of `U(dt) = exp(-i H dt/hbar)`: `x = dt * HalfWidth/hbar`
  (`hbar = 0.6582119624` eV fs, `units.hpp`), `c_0 = J_0(x)`, `c_n = (-i)^n 2 J_n(x)`,
  length `L` = the first `n >= 1` with `|J_n(x)| <= 1e-15` (that term excluded), then
  the global phase `exp(-i BandCenter dt/hbar)` applied to the result (legacy
  multiplies the accumulated vector by `exp0` at the end: do the same, as a final
  `scal`, not folded into the coefficients). `J_n` from GSL `gsl_sf_bessel_Jn` in
  `compat/` (scipy's `jv` differs in the last bits; the sympy `TimeEvolution` path uses
  scipy and is checked against compat to 1e-13 in the self-check).
- The derived series uses the same coefficients, with the `n = 0` term zero.
- The commutator operator is the velocity **rescaled like the Hamiltonian**
  (`RescaleCon2ChebyshevDomain`: `V <- a V + b I`), a legacy convention that is only
  right for `b = 0`; reproduce it (`RescaledOperator(V, window)`) and record it in
  `AGENTS.md` section 10 via the report.
- Accumulation: `mu(m, n) += scal * vdot(seed_bra, T_m seed_ket)` per state, `scal =
  2/N_states` halved at `m = 0`, per time slice; states in generator order.
- `.chebmomTD` layout: header `size " " bandwidth " " center " " numTimes " " dt " "` then
  `numMoms " " numTimes " "` (trailing spaces, `precision(digits10)`), then `re im`
  lines in `m`-major order (`m` outer, `n` inner). File names
  `Correlation<V><label>KPM_M<M>_state<tmax><statelabel>.chebmomTD` (MSD; note the
  `tmax` token formatted with the default stream format of the double, e.g. `20`) and
  `TimeProj<OP>-<P><label>KPM_M<M>_state<statelabel>.chebmomTD`.
- Reconstruction (`timeCorrelationsFromChebmom <file> <meV> <E_F>`): Jackson on `m`
  with the 1D rule, `x_F = (E_F - center)/HalfWidth`, `value_n = (sum_m delta(x_F, m)
  mu(m, n).real) * scale`, written as `n*dt value` with default precision, name
  `mean<label>EF<E_F as std::to_string>JACKSON.dat` (`std::to_string(0.5)` is
  `0.500000`).

## 5. The one self-check (`tests/dev/stage_4_check.py`)

(a) G7: the three `TimeProj` tensors by `projected_evolution_plan` vs the goldens,
`array_equal`; (b) G6: the `Correlation` tensor by `msd_plan` vs the golden,
`array_equal`, and the four `EF` curves, zero differing lines; (c) compat Bessel
coefficients vs scipy, max abs difference; (d) `time_correlation_plan` on chain1d
`VX, VX` (`M = 64`, 16 steps, `tmax = 5`, exact) vs the legacy
`inline_compute-kpm-CorrOp chain1d VX VX 64 16 5 exact` run on TCN in the same session
(the legacy binary is in `~/lsquant-unified/build`): `array_equal`, and the command
recorded in the report. Not a gate, not extended.

## 6. The gate program (`python -m lsquant.apps.stage4_time`)

```
--operators <dir> --label chain1d --velocity VX --window -2 2 --trace exact --states 256 --moments 256 --steps 64 --tmax 20 --word msd --out <dir>
--operators <dir> --label chain1d_mag --operator SX --projector PX --window <lo> <hi> --trace exact --states 128 --moments 32 --steps 96 --tmax 95 --word projected --out <dir>
--moments <file.chebmomTD> --broadening-meV 10 --fermi-energy <E> --out <dir>
```

(G7's window is the legacy `SpectralBounds` of `chain1d_mag.HAM.CSR` without `BOUNDS`
file: read `test/spin_precession.sh` and the golden header `BandWidth = 4.62` and
reproduce the bounds source exactly; say which it was.) Legacy names and layouts; one
stdout line per file.

## 7. Tutorials 02, 04, 05 through Python

`tutorial02.py` (localisation: disordered chain, `make_disordered_chain.py`, msd at
`E = 0`), `tutorial04.py` (precession: the three `TimeProj` runs and the `m = 0` row),
`tutorial05.py` (square lattice, `make_square.py`, msd at `E = 2` plus the DOS of stage 2).
Each writes the files the tutorial reads for its figures with the legacy names into
`<out>` (the tester matches them to `goldens/stage0/tutorials/0N/captured/`). No
plotting in the apps.

## 8. Report

As before; plus the exact legacy command used in self-check (d) and the bounds source
of G7. No "done", no `cmp` beyond the self-check, no push.
