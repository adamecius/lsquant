# Handoff 0011 (Agi) — a fixed, engine-independent stochastic test state

Roles as in AGENTS.md section 12 (developer Agi; Codex at reasoning effort `medium` if Agi has no
quota; tester on TCN; Fable reviews). Starts after 0010 is reviewed, from `experimental` with
0010 merged. Read AGENTS.md sections 3, 8, 10, 11 and this brief completely.

## 0. Why

The Aron numbers that disagreed (Bi 1.135 against 1.155, PdSe2 2.52 against 2.34; see
`handoff/0009-wrapup/profiling-results.md`) came from two runs with the same seed and different
trace vectors. The reduced-material study draws random phases in the eigenbasis of each projected
`H(k)` (`ProjectedEigenBasis.states`). Bi and PdSe2 have inversion and time reversal, so every band
is a Kramers pair, and `eigh` returns an arbitrary vector inside each pair that depends on the
round-off of `H(k)`. The Bloch sum changed between 6d (`22965bf`) and the integration (`9e9a40e`
and the 0009 correction 2) at round-off level; with the same seed, the phases landed on other
physical vectors. Evidence: identical resolved configuration and exact oracle (1.0797595 against
1.0797585), first moment `mu[0,0]` -0.00106 against -0.00013.

`RandomPhaseStates` has two further dependences: `std::uniform_real_distribution` is not specified
by the C++ standard (libstdc++, libc++ and MSVC differ), and `cos`/`sin` come from the platform
libm and the build profile. A before/after comparison therefore cannot rely on "same seed".

## 1. What to build

1. `TestPhaseStates` in `lsquant/core/algebra/state_generator.{hpp,cpp}`: component `i` of state
   `s` is `z^q / sqrt(N)` with `z = i` (the imaginary unit) and `q = h(seed, s, i) mod 4`, where
   `h` is SplitMix64 applied to the three integers. No floating-point function, no
   distribution object, no dependence on the order in which components are filled; parallel fill
   allowed. Same `StateGenerator` interface (`Next`, `Reset`, `Dimension`, `NumberOfStates`).
2. Binding in `lsquant/api/bindings.cpp` and `_core.pyi`.
3. The index `i` is always the physical orbital index of the real-space supercell (cell, orbital,
   spin) as the 6c mapping defines it. Bloch and eigen bases receive the state through the
   existing transforms (`B`, then `U(k)^+`); no basis draws its own phases. Write the mapping once
   and use it from `bases.py`, `disorder_bases.py` and `projection/bases.py`.
4. Run file: `expansion.trace = "test_phase"`, valid in every basis. `random` keeps
   `RandomPhaseStates` unchanged, because golden G4 depends on its bits.
5. Receipt: every run records `trace_states_sha256`, the SHA-256 of the first two states as
   little-endian complex128, so two runs prove they used the same vectors.

## 2. Checks (developer)

- The same `trace_states_sha256` from: golden and native profiles, 1 and 16 threads, and the
  real, Bloch and eigen bases of the same small model (graphene 4 x 4).
- The stochastic estimate with `test_phase` agrees with the exact trace within the usual
  `1/sqrt(R N)` error on G1 and G2 inputs (report the numbers; no new threshold).
- Every golden unchanged (`scripts/check_goldens.sh 7`).

## 3. Then (tester, final tests)

Repeat the reduced Bi and PdSe2 workload of `profiling-20261025` with `trace = "test_phase"` in the
eigen and Bloch bases of the final candidate: the moments must agree to round-off between the two
bases and between two builds (golden, native), so any remaining difference is a code difference. Also run Bi at R = 16 with five seeds of
`RandomPhaseStates` in the eigenbasis to show that the 1.135 to 1.155 spread is sampling.

## 4. Report

`handoff/0011-fixed-test-states/report.md` as usual: files, check output verbatim, deviations.
