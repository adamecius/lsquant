# Handoff 0006d (developer) — spectral projection toolset for non-orthogonal bases (no orthogonalisation)

Non-blocking side branch (user, 2026-10-23): `handoff/0006d-dev` from `origin/handoff/0006b-dev`.
Codex (effort high) implements, an Opus tester measures on Aron, Fable reviews. Read AGENTS.md
(sections 2, 3, 4, 6, 8, 10, 11, 12), `../0006-nonorthogonal-metric/theory.md`,
`../0006b-siesta-spin-hall/brief.md` (sections 3, 5, 11, 12) and its reviews, and the
existing downfolding code brought in under `tools/downfold_origin/` (from the user's other
project, `/root/scratch/gb` on Ionos: `downfold.py`, `downfold_scdm.py`,
`siesta2w90_sisl.py`, `study.py` and the notes they cite). Rules: no ninja programming;
every formula derived with sympy in `metric_forms.py` first.

## 0. Why

SIESTA multi-zeta bases put spurious states 30 to 60 eV away from the Fermi level: the
Bi spectral window is 97 eV for a 0.3 eV gap, so the KPM resolution at M = 1024 is of the
order of the gap. The user's `downfold.py` solves this by reducing to a small **orthogonal**
tight-binding model (Löwdin/polar gauge, Souza-Marzari-Vanderbilt disentanglement, SCDM).
The user wants the opposite route kept available: **project out the unwanted states while
staying in the non-orthogonal basis**, so that `S`, the dipole matrix, the velocity ladder
and every 6b operator keep their meaning, and only the spectral range shrinks.

## 1. The toolset (`python/lsquant/projection/`)

1. **Bring the origin code in** under `tools/downfold_origin/` unchanged (provenance, with
   a README saying where it came from), and keep its orthogonal path callable as one tool
   (`downfold_orthogonal`) with its own validation (band comparison on off-grid k).
2. **S-orthogonal spectral projector, per k.** With the generalised eigenpairs
   `H c_n = E_n S c_n`, `c^+ S c = 1`, the window `W = {n : E_n in [E_lo, E_hi]}`:
   `P = sum_{n in W} c_n c_n^+ S` (idempotent, `P^2 = P`, S-self-adjoint). Derive with
   sympy and implement the three forms:
   - **hard projection**: `H_P = S P H_tilde P` in the matrix sense, i.e.
     `H_P = S C_W diag(E_W) C_W^+ S`, same `S`; the out-of-window states become zero modes
     of the generalised problem; state where they land in the spectrum and how the metric
     routes see them;
   - **spectral shift** (push, not remove): `H_s = H + S C_O diag(E_shift - E_O) C_O^+ S`
     for the out-of-window set `O`, which moves those states to energies at the window
     edge, keeping every in-window eigenpair exact and the matrix size; the spectral range
     becomes `[E_lo, E_hi]`;
   - **operators**: every operator `O` of a word is projected consistently,
     `O_P = S P O P` in matrix form (velocity, spin, spin current), and the error of each
     word against the unprojected oracle is reported as a function of the window.
3. **Local (real-space) version** for the supercell route of 6c: the same projector built
   from the per-k eigenpairs and Fourier-transformed to `R` with Wigner-Seitz images
   (reuse the origin code's FFT/WS machinery), truncated at a hopping radius, with the
   truncation error reported.
4. **Run-file interface:** `"projection": {"kind": "none" | "hard" | "shift" | "orthogonal",
   "window": [E_lo, E_hi], "radius": ...}`; the resolved window, the number of kept bands
   per k, the new spectral bounds and the oracle error go into the receipt.

## 2. Self-check (`tests/dev/stage_6d_check.py`)

(a) sympy: `P^2 = P`, S-self-adjointness, in-window eigenpairs of `H_s` and `H_P` exact;
(b) rotated Kane-Mele with two extra high-energy orbitals per site added (a toy of the
SIESTA spurious states): with `shift` the spectral range drops to the window and the spin
Hall plateau (oracle and exact-trace KPM at the same M) equals the unprojected one to 1e-3,
while the KPM resolution at fixed M improves by the ratio of the ranges (print it);
(c) the orthogonal tool on the same toy reproduces the in-window bands to 1e-6.

## 3. Gate on Aron (tester)

Bi and PdSe2 (SIESTA): window +-3 eV and +-5 eV around E_F, `shift` and `hard`, 24 x 24
and 48 x 48: kept bands per k, new spectral range, oracle plateau error against the
unprojected 1.5527 / 1.1747 (Bi) and 2.08 / 1.98 (PdSe2), and the KPM plateau at M = 256
and 1024 with 16 eigenbasis vectors, with wall times. The point is the same plateau with
far fewer moments.

## 4. Report

`handoff/0006d-spectral-projection/report.md`: files, the sympy derivations, the formulas,
deviations, open questions, what to measure. No push.
