# Handoff 0005 — review 2 (Fable to Agi): the tester's run 1

Tester run on `handoff/0005-codex` at `81ad680` (fresh build of `7109c20`, before review 1):
the three bases agree with each other to 1e-14 (real vs Bloch moments, Bloch vs eigen
curves, with and without disorder), stages 1 to 4 are 65 of 65 byte-identical, all 48
stage-5 outputs are byte-identical between 1 and 16 threads, and the FFTW convention is
confirmed against a dense `B`. Two things in the brief were wrong and are corrected in
`brief.md` (section 3): the cell numbering (wannier2sparse has `n_1` fastest) and the sign
of `B` (your derivation was right). G9 is redefined in `codex.md` (Gamma-centred mesh,
Gaussian, compared with an exact Gaussian evaluation) because the oracle's shifted mesh
and the Jackson-vs-Gaussian shape made the old row untestable; nothing of that is yours.

Six items; nothing else is to be touched. Items 1 to 5 are the tester's D3, D1, D2, D5, D4.

1. **Cell numbering = wannier2sparse's.** `RealBasis` rebuilds the supercell with `n_1`
   slowest, so its CSR is a permutation of the wannier2sparse CSR and `B VX_w2s B^+ = v_x`
   fails on Haldane 48x48 by 0.9. Switch the flat index to `iR = n_1 + N_1 n_2 + N_1 N_2 n_3`
   everywhere (`RealBasis`, `BlochTransform`'s FFTW dimension array `{N_3, N_2, N_1}`,
   `BlockLocalStates`, `wannier.k_mesh`), so that the wannier2sparse supercell CSR *is*
   the real basis: `RealBasis` must read `<seed>.HAM.CSR` and `<seed>.VX.CSR` written by
   wannier2sparse when they exist (and build them itself, identically, otherwise). The
   self-check gains: `max |H_real_basis - H_w2s_csr| = 0` on graphene and Haldane 6x4x1
   (non-square on purpose) and the `B VX_w2s B^+ = v_x` check to 1e-10 on Haldane 48x48.
   Fix the report sentence that says wannier2sparse is row-major.

2. **`read_hr` must read every Wannier90-layout file**, including the sparse ones that
   `gen_models.py` writes (`ndegen` line, then only the non-zero hoppings): parse
   `R i j re im` records by content, assign `ndegen` by the order of first appearance of
   each `R` as Wannier90 does, and read the dense layout as a special case.

3. **`read_xyz` must accept wannier2sparse's `.xyz`** (count line followed directly by
   the orbital lines, no comment line) as well as the two-header form.

4. **`EigenBasis.exact_dos(energies, eta)` and `exact_sigma_xx(energies, eta)`**: the
   exact evaluation on the eigenvalues with a Gaussian of width `eta` (DOS: sum over `k`
   and `n` of the Gaussian, per cell; `sigma_xx`: Kubo-Greenwood `sum |v_nm|^2 g(E - eps_n)
   g(E - eps_m)` with the same normalisation as the stage-3 `kubo_greenwood` so that the
   two curves are comparable), using the stage-2 `SpectralFunction` machinery where a
   Gaussian is available and numpy otherwise. These are the oracle of the three bases
   and the G9 row of the gate (1e-8 against `hr_exactdiag` on the same mesh with the same
   Gaussian). Add them to the self-check: eigen exact vs `hr_exactdiag` wrapped on the
   Gamma-centred mesh of a `12x12` graphene supercell, 1e-10.

5. **The app as the brief says**: `--model`, `--disorder W`, `--disorder-seed S`
   (independent of the state seed), legacy-layout moment files and curves with the legacy
   names plus the new `.chebmom`, and `<stem>.basis.json`; `--window` default `eigen`
   (review 1) so it is not mandatory.

6. **Thread-safety statement for `SandwichOperator`** as for `BlochTransform` (shared
   mutable buffers, not safe for concurrent `Multiply` on one object); and `report.md`
   must mention every limitation you know of (D1, D2, the sign conflict with the first
   brief) instead of "no open questions".

Procedure: one commit per item on `handoff/0005-agi`; rsync to `~/lsquant-agi-build` on
TCN (`ssh jh-x`), rebuild, rerun the self-check and the app suite (foreground polling),
paste into `report.md`, commit. Do not push.
