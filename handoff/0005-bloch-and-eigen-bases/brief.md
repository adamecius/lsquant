# Handoff 0005 (Agi) — stage 5: the Bloch basis and the eigenbasis from one model file

Roles: Agi implements, the tester (`codex.md`) runs the gate and the figure on TCN, Fable
reviews. Read `AGENTS.md` (sections 3, 6, 8, 9, 11 are binding) and this brief completely.
Nothing is optional; deviations go in `report.md`. Starting point: branch `handoff/0005-agi`
from `origin/handoff/0004-agi` at `aba179a` (it carries the stage-4 code and the stage-3
energy-axis fix of `reconstruct.py`), then `git merge origin/cpp-unified` (documents only).

## 0. Scope and the idea of the test

One model file (`<seed>_hr.dat` + `<seed>.uc` + `<seed>.xyz`, the wannier2sparse layout)
feeds three representations of the same Hamiltonian, and the engine computes the same
words in all three:

- **real space**: the supercell CSR that wannier2sparse already writes (stage 1 to 4 path);
- **Bloch**: `N_k` dense blocks `H(k)` of size `W x W` on the `k` mesh commensurate with the
  supercell (`k = n/N`), as one `LinearOperator`; a real-space diagonal disorder `V` enters
  as the composite `B^+ V B` where `B` is the Bloch transform (FFT over the supercell with
  the atom gauge);
- **eigen**: the exact path: `eps_n(k)`, `U(k)` from `eigh(H(k))`, every operator
  `U^+ O(k) U`, every word evaluated exactly in that basis (the "exact concatenator"). It
  is the oracle of the other two and of the coefficient tables.

Traces are basis independent, so the DOS and `sigma_xx` of graphene and of Haldane must
agree in the three bases up to finite size (the real-space supercell and the Bloch mesh
are the same finite system; the eigen path on the same mesh is exact for that system).
G1 to G8 stay byte-identical: nothing in the stage-1 to stage-4 paths changes.

Out of scope: non-orthogonal metric (stage 6), parser and CLI (stage 7), FFT contraction
(stage 8), `_r.dat` covariant velocity (0006b).

## 1. Where things go

```
lsquant/core/algebra/
  block_diagonal_operator.hpp .cpp   BlockDiagonalOperator: N_k dense W x W blocks, Multiply = batched block GEMV (OpenMP over blocks)
  bloch_transform.hpp .cpp           BlochTransform: B and B^+ with FFTW (fftw_plan_many_dft, howmany = W) and the atom-gauge phases
  composite_operators.hpp .cpp       + SandwichOperator(B, V): y = B^+ V B x   (V any LinearOperator; here DiagonalOperator)
                                     + DiagonalOperator(values)
lsquant/core/algebra/state_generator.hpp .cpp   RandomPhaseStates already covers (k, orbital); + BlockLocalStates(block index) for the exact trace per block
python/lsquant/models/wannier.py     read_hr, read_uc, read_xyz (ndegen folded, wsvec optional), Hk(k), velocity dHk(k) in the atom gauge (berry_connection), the k mesh of a supercell, the Bloch phases
python/lsquant/bases.py              RealBasis, BlochBasis, EigenBasis: build the engine operators for a word from the model and a basis name; EigenBasis evaluates words exactly
python/lsquant/apps/stage5_bases.py  the gate program (section 6)
tests/dev/stage_5_check.py           the single self-check (section 5)
lsquant/api/bindings.cpp             + the above classes
```

## 2. Restrictions

All of stages 1 to 4. In addition: the engine (`MomentEngine`, `Plan`, `ChebyshevRecursion`)
does not change; a basis is a way to build `LinearOperator`s and `StateGenerator`s, nothing
else; no file format for blocks (the blocks travel from Python as one contiguous complex
array `(N_k, W, W)`, copied once into the operator); FFTW is linked by `lsquant_core`
(it is a legacy dependency on TCN; `find_package` or `pkg-config fftw3`); the FFT
normalisation and the atom-gauge phases are fixed by section 3 and checked by the
self-check against the dense `B`; the eigen path uses `Eigen::SelfAdjointEigenSolver`
per block in C++ or `numpy.linalg.eigh` in Python, your choice, but the oracle numbers in
the report come from Python.

## 3. Definitions that must be exact

Supercell `N_1 x N_2 x N_3` cells, `W` orbitals per cell, orbital positions `tau_alpha`
(from `.xyz`, Cartesian), lattice vectors `a_i` (from `.uc`). Real-space index
`(R, alpha)` with `R = (n_1, n_2, n_3)`, `n_i in [0, N_i)`, flat `iR = n_1 + N_1 n_2 + N_1 N_2 n_3`
(`n_1` fastest: this is the wannier2sparse cell order, so that its supercell CSR is the real
basis without any permutation; corrected in review 2, the first version of this brief had
`n_1` slowest); Bloch index `(k, alpha)` with `k = sum_i (m_i / N_i) b_i`, `m_i in [0, N_i)`,
same flat order (`m_1` fastest).

    H(k)_{alpha beta} = sum_R exp(+i k . (R + tau_beta - tau_alpha)) H_{alpha beta}(R)       (atom gauge)
    v_a(k)_{alpha beta} = sum_R i (R + tau_beta - tau_alpha)_a exp(+i k . (R + tau_beta - tau_alpha)) H_{alpha beta}(R)
    B_{(k,alpha),(R,beta)} = delta_{alpha beta} exp(-i k . (R + tau_alpha)) / sqrt(N_k)        (unitary; sign fixed by B H_R B^+ = H(k) above, corrected in review 2)

`B` maps real-space coefficient vectors to Bloch coefficient vectors: `psi_k = B psi_R`,
`psi_R = B^+ psi_k`. With FFTW: `B^+` is `fftw FORWARD` on the `(N_1, N_2, N_3)` grid with
`howmany = W`, stride `W`, after multiplying by `conj(atom_phase)`; `B` is `fftw BACKWARD`
followed by `atom_phase` (the `1/sqrt(N_k)` lives in `atom_phase`). Record in the report
which sign of the FFTW exponent corresponds to which of `B`, `B^+`, and verify against the
dense `B` in the self-check (section 5). The velocity `v_a(k)` above is exactly the
`berry_connection` velocity of wannier2sparse (`-i (R.lat + dr) H` in its sign convention):
the block-diagonal `v_a` and the supercell CSR `VX` of wannier2sparse must agree through
`B` to 1e-10 relative (self-check).

Disorder: Anderson `V = diag(w_i)` on the real-space sites with `w_i` uniform in
`[-W/2, W/2]`, `std::mt19937` seed, generated in Python (numpy `default_rng`), passed as
the values of `DiagonalOperator`. In the Bloch basis the Hamiltonian is the composite
`OperatorSum(H_k, SandwichOperator(B, V))`; in real space it is the CSR plus the diagonal.
The two must give the same DOS moments to 1e-10 relative (not byte: the FFT rounding).

Eigenbasis: `H(k) = U diag(eps) U^+` per block; an operator `O(k)` becomes `U^+ O(k) U`;
a word `Tr[A f(H) B g(H)]` on the mesh is `sum_k sum_{nm} A_nm(k) f(eps_m) B_mn(k) g(eps_n)`,
and the Chebyshev form of the same word is `sum_k sum_{nm} A_nm B_mn sum_{ab} c_a(E) c_b(E')
T_a(eps~_m) T_b(eps~_n)`: `EigenBasis` evaluates both, the first with the sympy
`SpectralFunction` evaluated on the eigenvalues, the second with the coefficient tables of
stage 2 and 3, and reports the difference (this is the check of section 6 of AGENTS.md
that the tables reproduce `f(eps; E)`).

## 4. What the engine sees

Nothing new. `BlochBasis.hamiltonian()` returns a `LinearOperator` (`BlockDiagonalOperator`
or the composite with disorder), `BlochBasis.operator("VX")` returns the block-diagonal
velocity, `BlochBasis.states(kind)` returns `RandomPhaseStates` over `N_k W` or, for the
exact trace, the local states; the stage-2/3 plans (`dos_plan`, `kubo_plan`) run unchanged
on them. The spectral window for a Bloch run comes from the eigenvalues of the blocks
(`EigenBasis.window(pad)`), never from Gershgorin (on Wannier models Gershgorin is several
times wider than the spectrum and wastes the moments: say so in the docstring).

## 5. The one self-check (`tests/dev/stage_5_check.py`)

Graphene from `wannier2sparse/examples` models (`gen_models.py` graphene, `N = 24`):
(a) `B` from FFTW vs the dense `B` built in numpy on random vectors: max relative difference
(expect 1e-13); `B^+ B = 1` on random vectors; (b) `H_k = B H_R B^+` on random vectors vs
`BlockDiagonalOperator`: 1e-12; (c) `v_x` block-diagonal vs `B VX_csr B^+`: 1e-10; (d) DOS
moments `M = 64`, random phase seed 12345, real vs Bloch: max relative difference (1e-10)
and with Anderson `W = 1`, same seed for the disorder: 1e-10; (e) the eigen path: DOS
curve from `dos_plan` on the Bloch basis vs `EigenBasis` exact (Jackson at the same `M`
applied to the exact spectrum, i.e. the same kernel on both sides): max abs difference
(expect 1e-12, this is the table check); (f) `sigma_xx` of graphene, `kubo_plan` on the
Bloch basis vs eigen exact with the same kernel: 1e-10 relative at `E = 0.5`. Not a gate,
not extended.

## 6. The gate program (`python -m lsquant.apps.stage5_bases`)

```
--model <dir with seed_hr.dat seed.uc seed.xyz> --seed graphene --supercell 48 48 1 --basis real|bloch|eigen
--word dos|kuboxx --moments 256 --trace random --states 32 --seed 12345 [--disorder W --disorder-seed 7]
--broadening-meV 20 --out <dir>
```
Writes the legacy-layout moment files and curves with the legacy names plus
`<stem>.basis.json` (basis, mesh, window source, disorder); for `eigen` the curve is the
exact evaluation with the same kernel. G9: the tester compares the three bases with
each other and with `wannier2sparse/tools/hr_exactdiag.py dos|sigma` (k-mesh exact,
independent code) at the brief's tolerances (`codex.md`), for graphene and Haldane,
with and without disorder (`W = 0` and `W = 1` eV), and reruns stages 1 to 4 byte for byte.

## 7. Report

`handoff/0005-bloch-and-eigen-bases/report.md`: files; the FFTW sign and normalisation
conventions with the self-check numbers; the self-check output; one `stage5_bases` run per
basis and word with the stdout line; deviations; open questions. No "done", no gate
comparison beyond the self-check, no push.
