# Handoff 0006c (developer) — disorder in a non-orthogonal basis, both routes

Non-blocking side branch (user, 2026-10-23): `handoff/0006c-dev` from `origin/handoff/0006b-dev`
(it needs the 6b SIESTA reader and the per-block metric); merged after 6b closes. Codex
(effort high) implements, an Opus tester measures on Aron, Fable reviews. Read AGENTS.md
(sections 2, 3, 6, 8, 10, 11, 12), `../0005-bloch-and-eigen-bases/brief.md` (section 3: `B`,
the atom gauge, disorder as `B^+ V B`), `../0006-nonorthogonal-metric/theory.md`,
`../0006b-siesta-spin-hall/brief.md` (sections 3, 4, 11, 12) and its reviews 1 to 4.
Rules: no ninja programming; every formula derived with sympy in `metric_forms.py` first.

## 0. Scope

1. **Anderson disorder in a non-orthogonal basis.** With on-site energies `w_i` (uniform in
   `[-W/2, W/2]`, numpy `default_rng(seed)` as stage 5), the matrix of the potential in the
   atomic basis is `V = 1/2 (D S + S D)` with `D = diag(w_i)` (the symmetric, Hermitian
   form; state its derivation and why `D` alone is wrong when `S != I`). The disordered
   Hamiltonian is `H + V` with the same `S`.
2. **Bloch route (the FFT trick, no supercell matrix):** `V` acts on a Bloch vector as
   `1/2 (B^+ D B S(k) + S(k) B^+ D B)`, with `B^+ D B` the stage-5 `SandwichOperator`
   (FFTW) and `S(k)` the block-diagonal overlap; the metric (Löwdin `S^-1/2 (.) S^-1/2`
   per block, or generalised `S^-1 (.)`) is applied exactly as for every other operator
   of 6b. Note that with disorder the eigenbasis is no longer per-k: the trace is random
   phase over the Bloch layout (the 6b sampler after review 2) and the reference is the
   real-space route below, not the per-k oracle.
3. **Real-space route (the reference):** the SIESTA model on an `N1 x N2` supercell as CSR
   (`H`, `S`, `dH`/velocity, spin, spin current) built from the same `H(R)`, `S(R)` and the
   6b velocity ladder in real space (theory.md equation 3 with the dipole matrix), plus
   the same `V`; the metric by the stage-6 appliers (exact for small supercells,
   Chebyshev for large).
4. **Words:** DOS and the spin Hall word, both routes, both metric routes.

## 1. Gate and self-check

`tests/dev/stage_6c_check.py`: (a) `W = 0`: both routes reproduce the 6b numbers (Bloch
moments equal to 6b's to round-off; real-space moments equal to Bloch to 1e-10 on the same
finite system); (b) rotated Kane-Mele with disorder `W = 1` eV, seed 7, 12 x 12: Bloch and
real-space DOS and spin Hall moments agree to 1e-9 relative (the stage-5 tolerance), same
random phases mapped by `B`; (c) `V` Hermitian, and `V -> D` when `S -> I` (the orthogonal
limit reproduces stage 5's disorder byte for byte on the orthogonal fixture). Gate on
Aron: Bi at `W = 0, 0.1, 0.3` eV, 24 x 24, DOS and spin Hall, both routes, the plateau
against `W` (the shape of S. Gimenez's W sweep), with wall times.

## 2. Report

`handoff/0006c-nonorth-disorder/report.md`: files, the sympy derivation of `V`, the
formulas, deviations, open questions, what the tester must measure. No push.
