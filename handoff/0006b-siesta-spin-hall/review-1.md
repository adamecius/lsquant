# Handoff 0006b — review 1 (Fable to Codex): the tester's first run on `b51a8c6`

DRAFT, to be completed with the tester's `tests.md` (branch `handoff/0006b-codex`). Tester:
Opus on Aron, worktree at `b51a8c6`, build clean (0 warnings), logs in `~/t6b/` on Aron.

## What is right

- Rotated Kane-Mele (the exact fixture of brief section 8), with the Gamma tolerance
  patched outside the repo: sympy forms and proofs as the brief asks (generalised forms
  S-Hermitian; Löwdin and anticommutator Hermitian; generalised and Löwdin traces equal
  for every word of the catalogue; the anticommutator equal for one operator, not for
  products); covariant velocity equal to `T^+ V T` to 8.9e-16; Hermiticity 0;
  Hellmann-Feynman 2e-10; oracle `sigma^z_xy = 2.000000000000`, spin Chern 2, Z2 = 1 on
  both exact routes; exact-trace KPM (`M = 256`, 12 x 12) 2.0055 against its oracle 1.9999;
  Bloch and eigen moments of the same sample equal to 1.9e-14; DOS with one eigenphase
  vector exact to 1e-15.
- SIESTA reader on Bi(111) (`bi.TSHS`, 36 orbitals, W = 72): sisl conventions pinned
  (2e-15), generalised eigenvalues against sisl `eigh` to 1e-6, eq. (4) velocity Hermitian
  to 3e-7 before folding and Hellmann-Feynman to 3e-5, the dipole quadrature certified
  (order 96, refinement 9.9e-6 against 1e-5, overlap against stored `S` 9.2e-7, `A` Hermitian
  to 8e-17, diagonal zero), covariant velocity Hermitian to 8e-13. No Wannier import, no
  silent fallback, no zero-dipole shortcut, no eval.
- The midpoint dipole `A = Q - (tau_i + tau_j + R) S / 2` corrects the row-centred
  definition of brief section 3; accepted, the brief is amended (section 13).
- Bi physics from the oracle (48 x 48, Löwdin eigenvectors): spin Chern `C+ = +1`,
  `C- = -1`, Z2 = 1, energy gap 0.303 eV; flat plateau over the whole gap:
  bare 1.5527, berry_connection 1.5527, covariant 1.1747 e^2/h (Pauli); spin mixing
  `min|s| = 0.513` over the occupied `P sigma_z P`, `|<sigma_z>|` of the highest occupied
  Kramers pair down to 0.099. The dipole term moves the plateau by 24 %, the user's
  suspicion about the velocity confirmed; the non-integer value is spin mixing, not a
  prefactor (Kane-Mele gives 2.000 with the same code).

## Items, one commit each, on `handoff/0006b-dev` after `b51a8c6`

1. **Gamma tolerance.** `checked_hermitian` (`python/lsquant/models/nonorthogonal.py`
   about lines 271-285, reached from `bases.py:1281`) divides each block's anti-Hermitian
   norm by that block's own norm, so a block that is physically zero at Gamma (`V_y`,
   norm 1.4e-15) is refused with "relative norm 0.13". Scale by the largest block norm of
   the operator (or an absolute floor tied to it); keep refusing real defects.
2. **The generalised route refuses its own operator.** `g11_bi_berry_connection_generalised_shc`
   died after 70 minutes with "Operator 'JXSZ' in word ... is not Hermitian": on the
   generalised route the applied matrix is `S^-1 J`, which is S-Hermitian, not Hermitian
   (your own report says so). The Hermiticity check must test the physical matrix before
   the metric is inserted (or test S-Hermiticity after), and it must run **before** any
   heavy work: validate every operator of the word at compile time, not after the oracle
   and the eigen reference have been computed.
3. **Exact curve grid.** `exact_shc.json` samples 1001 points over the whole window
   (97 eV for Bi, 0.097 eV spacing, two points in the gap). The oracle grid is a run-file
   choice (`reconstruction.energies` or `[E_min, E_max, n]`), default the gap plus a margin.
4. **Progress and timing.** Each run writes one line per phase (metric build, dipole
   quadrature, operators, moments, reconstruction, oracle, eigen reference) with its wall
   time to the log and to `run.resolved.json`'s receipt; 70-minute runs without a line in
   the log are not acceptable.
5. **Cost.** The Bi spin Hall run (block basis, 2304 blocks of 72, `M = 1024`, 16
   vectors) takes 71 minutes with `OMP_NUM_THREADS = 48` while the process uses 2.6 cores
   (264 % CPU, 95 threads). The arithmetic is small (rank-2 contraction
   16 x 1024^2 dot products of length 166k is minutes on one socket), so the time is in
   serial Python or in a Python-level per-block loop around the engine. Profile one short
   run (`M = 128`, 2 vectors; `cProfile` or `py-spy` if available) and move whatever is in
   Python per block into one C++ call (the stage-5 `BlockDiagonalOperator` and the
   stage-6 appliers already exist for this). Report the profile before and after.
6. **Brief section 5 is revised, the code follows.** The eigenbasis random-phase trace is
   exact for the DOS (confirmed, 1e-15) but not usable for the rank-2 spin Hall word: on
   Kane-Mele, 100 seeds x 16 vectors give an unbiased mean with a standard deviation of
   25 per sample at 12 x 12 and 6.7 at 48 x 48 (noise about 1/sqrt(N_k)); phases in the
   orbital Bloch basis are 0.6 times less noisy at both meshes; the metric is not the
   cause (identical on the orthogonal model). Therefore: (a) every test and gate
   comparison of a rank-2 word uses the **exact block trace in the eigenbasis** (H
   diagonal there; this is what the oracle does and what `EigenMetricBasis` must expose
   as `trace: "exact"`); (b) random phase stays as the production estimator, with the
   phases set in the orbital Bloch basis, and the self-check reports its standard
   deviation at the gate mesh so the vector count is chosen from a number, not a belief;
   (c) the G11 run files for the goldens use the exact trace.
7. **Report.** Append "Review 1 applied" to `report.md`: one line per item, the profile
   numbers of item 5, the self-check output expected from the tester (you cannot run it;
   say so). Do not push.

## Open physics questions (not items)

- Bi covariant plateau 1.17 against the legacy bismuthene 1.05 (one-sided, non-Hermitian
  velocity, other geometry?): which part of the difference is physics; needs S. Gimenez's
  geometry and units.
- PdSe2 needs a denser NSCF grid before it is a test (notes section 10).
