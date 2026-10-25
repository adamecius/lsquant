# Handoff 0006b (Agi) — stage 6b: SIESTA inputs, the non-orthogonal velocity ladder, and the spin Hall of a topological insulator

DRAFT v0 (2026-10-21). Becomes binding when stage 7 is closed and the PdSe2 velocity-mode
run on Aron (`notes-2026-10-21.md` section 8) is in. Roles: Agi implements, the tester
(`codex.md`, to be written) reproduces the oracle and the goldens on Aron, Fable reviews.
Read `AGENTS.md` (sections 2, 3, 6, 8, 10, 11, 12 binding), `../0006-nonorthogonal-metric/theory.md`
(section 3 is the velocity; its equation numbers are used below) and `notes-2026-10-21.md`
completely. Nothing is optional; deviations go in `report.md`. Branch `handoff/0006b-agi`
from `origin/cpp-unified` after the stage-7 merge.

Machine: Aron (`tcn-aronnode`, 96 cores, no internet), reached from the fat node:
`ssh jh-x "ssh -i ~/.ssh/socks_aron jgarcia@10.100.200.20 '<cmd>'"`; rsync with
`-e "ssh jh-x ssh -i ~/.ssh/socks_aron"`. Toolchain `~/lsq-env` (`export PATH=$HOME/lsq-env/bin:$PATH`),
recipe `~/lsquant/aron_build.sh`, `sisl` 0.16.4 installed there. Inputs under
`~/from_fatnode/shc_bi/run` (`bi.fdf`, `bi.TSHS`), `~/from_fatnode/shc_pdse2/run` (`pdse2.fdf`,
`pdse2.TSHS`), the Wannier PdSe2 under `~/from_fatnode/pdse2/ex`, wannier2sparse built under
`~/from_fatnode/wannier2sparse/build`. Byte verdicts of goldens are still issued on the fat
node (different toolchain, last-bit differences on Aron: `notes` section 8 and `tests.md` of 0007).

## 0. Scope and the idea of the test

Stage 6 made the engine metric-aware for a basis given as `(H', S', V')`. Stage 6b feeds it a
real non-orthogonal basis from SIESTA and builds the operators it needs from `H(R)`, `S(R)`,
the orbital centres and the intra-cell dipoles. The physics target is the spin Hall
conductivity `sigma^z_xy` of a 2D topological insulator (Bi(111) bilayer, PdSe2 monolayer),
whose gap plateau is an integer in `e^2/h` with `sigma_z` in Pauli units (`C_up - C_down`;
Kane-Mele gives 2 by this convention: `notes` section 1). The integer is the oracle: it does
not depend on the DFT code or the basis, so the same material through the orthogonal
Wannier route (wannier2sparse, `--velocity-mode covariant`, user: "this one has the velocity")
and through the SIESTA route must give the same plateau and the same spin Chern number.

What the user asked for (2026-10-21): implement in the SIESTA route what wannier2sparse
already has for Wannier models, the velocity ladder with its covariant term, and `S^-1`
exact, per block, computed once ("single shot"); the expansion of `S^-1` is a separate
study (0006c), not this stage.

Out of scope: disorder in a non-orthogonal basis (the matrix of an on-site potential is not
diagonal when `S != I`; 0006c), the rank-3 and optical words, performance (stage 8),
Chebyshev or ISCHO `S^-1` on these materials (they were the reason the first runs were
aborted, `notes` section 6).

## 1. Where things go

```
python/lsquant/models/siesta.py        read a SIESTA run with sisl: H(R), S(R) (spin-orbit, 2x2 spinor blocks), orbital centres tau_i,
                                       lattice, E_F, species/orbital metadata; the intra-cell dipole matrix A_a(R) (section 3);
                                       H(k), S(k), dH_a(k), dS_a(k) in the atom gauge on the supercell mesh k = n/N
python/lsquant/models/nonorthogonal.py + bloch_velocity(H, S, dH, dS, A, mode): the ladder of section 3, one function, used by every caller
                                       + spin_operators(S, ordering): Sigma_a = S (x) sigma_a and the applied form 1 (x) sigma_a
                                       + spin_current(V_a, sigma_b): J^b_a = 1/2 (V_a (1 (x) sigma_b) + (1 (x) sigma_b) V_a)
python/lsquant/bases.py                + BlochMetricBasis: the stage-5 BlochBasis with a per-block S(k): hamiltonian(), overlap(),
                                       operator(name) for VX VY VZ SX SY SZ JXSZ ..., metric policy per block (section 4)
python/lsquant/oracle/spin_hall.py     exact evaluation on the generalised eigenproblem per k: spin Berry curvature, sigma^z_xy(E),
                                       spin Chern number (section 5); numpy only
python/lsquant/apps/stage6b_siesta.py  the gate program (section 6)
runs/goldens/g11_*.json                run files (stage-7 words) for every golden of section 7
tests/dev/stage_6b_check.py            the single self-check (section 8)
lsquant/core/...                       no change expected. If the per-block exact S^-1 needs a C++ applier for speed
                                       (dense W x W per k, OpenMP over k), it is BlockMetricApplier in core/algebra, and nothing else.
```

## 2. Restrictions

All of stages 1 to 7. In addition: `sisl` is imported only in `models/siesta.py`; nothing
under `lsquant/core` or `words/` learns about SIESTA; the velocity formula exists once
(`bloch_velocity`) and the three modes carry the wannier2sparse names `bare`,
`berry_connection`, `covariant`; `covariant` without a dipole matrix is refused, never
approximated; the metric applier on these bases is the exact per-block inverse (and inverse
square root for Löwdin), computed once per basis and cached; the word the user writes is
the physical word (`Tr[JXSZ delta(E-H) VY dG-(E)/dE]`), the metric and the route are run-file
choices (`"metric": {"policy": "nonorthogonal", "route": "generalised", "applier": "exact"}`,
`"basis": {"name": "bloch", "supercell": [N, N, 1]}`); the oracle uses only numpy on the
same matrices; units are written into the curve header (`e^2/h, sigma_z Pauli`).

## 3. Definitions that must be exact

Orbital `i` at `tau_i` (Cartesian, from sisl), cell `R`, spinor index inside the orbital
block as sisl stores it (state it in the report). Atom gauge:

    H(k)_ij  = sum_R exp(+i k.(R + tau_j - tau_i)) H_ij(R),   S(k) likewise
    dH_a(k)  = sum_R i (R + tau_j - tau_i)_a exp(...) H_ij(R),   dS_a(k) likewise          (atom gauge)
    dH_a^cell(k) = sum_R i R_a exp(+i k.R) H_ij(R)                                        (cell gauge, for "bare" only)

Metric: `S^-1(k)` per block, dense, exact (`numpy.linalg.inv` on the Hermitian block, or
Cholesky); `S^-1/2(k)` by `eigh`. Both computed once per basis and cached in the basis object.

The velocity ladder (theory.md equations 3 and 4; `[X, Y]` matrices, `A^+` the adjoint):

    bare:             V_a = dH^cell - 1/2 ( dS^cell S^-1 H + H S^-1 dS^cell )
    berry_connection: V_a = dH       - 1/2 ( dS       S^-1 H + H S^-1 dS )                 (4)
    covariant:        V_a = (4) + i ( H S^-1 A_a - A_a^+ S^-1 H )

with `A_a(k) = sum_R exp(+i k.(R + tau_j - tau_i)) A_a,ij(R)` and the intra-cell dipole
`A_a,ij(R) = <phi_i,0 | (r - tau_i)_a | phi_j,R>`, computed by numerical quadrature of the
SIESTA radial basis functions that sisl reads from the `.ion.xml` files (same sparsity as
`S(R)`, Hermitian as an operator: `A_a,ij(R) = conj(A_a,ji(-R)) + (tau_j - tau_i)_a S_ij(R)`
up to the gauge term; derive and state the exact relation in the report; the diagonal
`A_a,ii(0)` vanishes for orbitals of definite parity). All three `V_a` are Hermitian
matrices; the operator the engine applies is `S^-1 V_a` (theory.md section 1), inserted by
the metric policy, never by hand.

Spin: `Sigma_b = S (x) sigma_b` is the matrix of the spin operator; its applied form is
`S^-1 Sigma_b = 1 (x) sigma_b`. Spin current (operator product `AB -> A S^-1 B`):
`J^b_a = 1/2 ( V_a (1 (x) sigma_b) + (1 (x) sigma_b) V_a )`, a Hermitian matrix.

Units: `sigma^z_xy(E) = (1 / (N_k A_cell)) * 2 pi * [Kubo-Bastin sum in e^2/h units of
stage 3]`, with `sigma_z` Pauli (eigenvalues +-1), so that Kane-Mele gives 2.000. Fix this
once in `reconstruct.py` as the formula `shc` and name the convention in the header.

## 4. What the engine sees

Nothing new: a `BlochMetricBasis` returns `LinearOperator`s (`BlockDiagonalOperator` for
`H(k)`, `S(k)`, `V_a(k)`, `J(k)`; the exact per-block `S^-1` as another block-diagonal
operator), `StateGenerator`s as stage 5, and the stage-6 `MetricPolicy` with the exact
applier bound to those blocks. Window from the generalised eigenvalues of the blocks
(`eigh(H(k), S(k))`), never Gershgorin. The Löwdin route uses `S^-1/2` per block.

## 5. The oracle (`oracle/spin_hall.py`)

Per `k`: `eigh(H(k), S(k)) -> E_n, c_n` (`c^+ S c = 1`); matrix elements
`J_nm = c_n^+ J c_m`, `V_nm = c_n^+ V_y c_m`; spin Berry curvature
`Omega_n(k) = -2 Im sum_{m != n} J_nm V_mn / (E_n - E_m)^2`; `sigma^z_xy(E_F)` from the
occupied sum with a Gaussian broadening on `E_n` (as G9); the spin Chern number in the gap
as the plateau value and, independently, Prodan's `P sigma_z P` spin Chern number by the
Fukui-Hatsugai-Suzuki lattice method on the sectors of `P sigma_z P` (integer; `Z2 = C_s mod 2`).
This is the number the KPM plateau is compared with, and the number that must agree
between the SIESTA route and the Wannier route of the same material.

**Trace for every test of this stage (user, 2026-10-21): random phase in the eigenbasis.**
In the basis of the generalised eigenvectors `(k, n)` a random-phase vector has unit
modulus on every `|k n>`, so `<r| A f(H) |r>` equals the exact trace for every word whose
operators are diagonal in that basis (DOS, `Tr[f(H)]`), and for the rank-2 words the only
fluctuating terms are the interband cross terms inside one `k` block (`J`, `V` are
block-diagonal in `k`): `W` terms per block, not `N_k W`. So the eigenbasis random-phase
trace (stage-5 `EigenBasis.states("random")`, seed 12345) is the trace of the oracle and of
every self-check and gate comparison; the exact trace is used only on the rotated
Kane-Mele fixture to show the two agree (report the difference). The production path,
the Bloch basis with the block engine and the per-block `S^-1`, is compared against that
eigenbasis result at the gate sizes (section 6), not against an exact trace.

## 6. The gate program (`python -m lsquant.apps.stage6b_siesta`)

```
--siesta <dir with .fdf .TSHS .ion.xml> --label bi --supercell 48 48 1 --velocity-mode bare|berry_connection|covariant
--route lowdin|generalised --word shc|dos --moments 1024 --trace random --states 16 --seed 12345
--broadening-meV 10 --out <dir>   [--oracle: also write exact_shc.json, spin_chern.json]
```
Legacy-layout moment files and curves, `run.resolved.json`, and the oracle files. Every
run is also expressible as a stage-7 run file (`runs/goldens/g11_*.json`): the app is the
run file's executor, not a second code path.

## 7. The goldens (section 5 of `notes`)

Bi(111) bilayer (`bi.fdf`, DZP+SOC, pseudo-dojo nc-fr-04) and PdSe2 monolayer
(`pdse2.fdf`): the oracle curves and spin Chern numbers on the `48 x 48` mesh for the three
velocity modes; the KPM moments and curves of the validated mode (seed 12345, `M = 1024`,
16 random vectors), byte for byte, as regression goldens; PROVENANCE with SIESTA version,
fdf, pseudopotentials, sisl version, commit of wannier2sparse. Cross-checks the tester
runs: PdSe2 SIESTA vs PdSe2 Wannier `covariant` (same integer), Bi vs the legacy bismuthene
curve of S. Gimenez (`figures/out/shc_bi/santiago_*`; 0.97 vs 1.05, the `2 pi` of
`notes` section 4 must be settled here).

## 8. The one self-check (`tests/dev/stage_6b_check.py`)

Rotated Kane-Mele as the exact fixture (stage-6 `rotate_basis` with a local, intra-cell
`T`): (a) `H'(k) = T^+ H(k) T`, `S'(k) = T^+ T` from the rotated `H'(R)`, `S'(R)` through
`models/siesta.py`-style Bloch sums (1e-12); (b) with `A' = T^+ R T` (`R = diag(tau)`
intra-cell), the `covariant` velocity equals `T^+ V T` of the orthogonal model to 1e-10,
and `berry_connection` differs from it by exactly the dipole term (print the norm);
(c) Hermiticity of all three `V_a` and of `J` (1e-13); (d) Hellmann-Feynman:
`c_n^+ V_a c_n = dE_n/dk_a` by central finite differences (1e-6); (e) the oracle on rotated
Kane-Mele gives `sigma^z_xy = 2.000` in the gap and spin Chern 2 (1e-6); (f) a short KPM
(`M = 256`, exact trace on a `12 x 12` mesh) on the generalised route reproduces the oracle
plateau to 1e-2. Not a gate, not extended.

## 9. Report

`handoff/0006b-siesta-spin-hall/report.md`: files; the sisl conventions you verified
(spinor ordering, gauge, units of `H`, `S`); the dipole quadrature (grid, cutoff, error
estimate); the self-check output; one `stage6b_siesta` run per material and mode with the
oracle numbers; deviations; open questions. No "done", no gate comparison beyond the
self-check, no push.

## 10. Goal restated (user, 2026-10-22) and references

The user's instruction: verify and correct the velocity operator, on both routes, until the
spin Hall plateau of a topological insulator comes out of the KPM; the rotated Kane-Mele
fixture is the exact check, Bi(111) and PdSe2 the materials. Evidence already in hand
(`notes-2026-10-21.md`, sections 8 to 11): the SIESTA eq. (4) velocity is Hermitian and
passes Hellmann-Feynman; the wannier2sparse `covariant` velocity is non-Hermitian because
`_r.dat` is used raw (not hermitised, not divided by `ndegen`), and becomes Hermitian to
1e-16 once `A(R) <- 1/2 [A(R) + A(-R)^+] / ndegen(R)`; the PdSe2 4 x 4 Wannier model shows
no plateau in any mode and the exact reference agrees, so PdSe2 needs a denser NSCF grid
before it can be a test; the legacy one-sided velocity is non-Hermitian by 20 to 30 eV A on
Bi while its diagonal is right. Developer for this stage: Codex (effort high) while Agi's
quota is out; tester: an Opus agent on Aron; Fable reviews.

References (read before coding; `../0006-nonorthogonal-metric/literature.md` has the
stage-6 list):
- X. Wang, J. R. Yates, I. Souza, D. Vanderbilt, Phys. Rev. B 74, 195118 (2006): the
  Wannier-interpolated velocity and Berry curvature; the covariant velocity
  `v_a = dH/dk_a - i [H, A_a]` with the Wannier position matrix `A_a(R) = <0n| r_a |R m>`,
  and the sum rules it must satisfy (Hermiticity, `A_a(R) = A_a(-R)^+`).
- D. Vanderbilt, Berry Phases in Electronic Structure Theory (CUP 2018), chapter 5 on
  Wannier functions and the position operator in a non-orthogonal basis.
- E. Prodan, Phys. Rev. B 80, 125327 (2009): the spin Chern number from `P sigma_z P` when
  `S_z` is not conserved (the oracle of section 5).
- T. Fukui, Y. Hatsugai, H. Suzuki, J. Phys. Soc. Jpn. 74, 1674 (2005): the lattice
  Chern number (integer by construction) used for the oracle.
- C. L. Kane, E. J. Mele, Phys. Rev. Lett. 95, 226801 (2005): the model of the self-check
  (`sigma^z_xy = 2 e^2/h` with Pauli `sigma_z`, section 3 units).
- A. Weisse, G. Wellein, A. Alvermann, H. Fehske, Rev. Mod. Phys. 78, 275 (2006): KPM,
  Jackson kernel, the Kubo-Bastin expansion; and `theory.md` of stage 6 for the
  non-orthogonal trace (`S`-weighted inner product, operators as `S^-1 O`).
- N. Papior et al., sisl (doi 10.5281/zenodo.597181): `Hamiltonian.Hk/Sk/dHk/dSk`, gauge
  `atomic` versus `lattice`, spin-orbit block layout; SIESTA basis functions from
  `.ion.xml` for the dipole quadrature.
- wannier2sparse (https://github.com/adamecius/wannier2sparse): `include/tbmodel.hpp`
  `createCovariantVelocity_list`, `readPositionMatrix` (the `ndegen` pitfall is in its own
  comment), `hopping_product_R`.

## 11. Binding clarification (user, 2026-10-22): non-orthogonal is the target, not a cross-check

The systems this stage is for cannot be Wannierised. The SIESTA route must therefore stand
on its own: every operator (`H`, `S`, `dH`, `dS`, the intra-cell dipole matrix `A`, the
spin and spin-current matrices) comes from the SIESTA output and the SIESTA basis
functions, never from a Wannier90 file; the oracle (section 5) runs on those same SIESTA
operators with the generalised eigenproblem; the KPM runs on them with the per-block exact
`S^-1`. The Wannier route (wannier2sparse, `_r.dat`) is only an optional cross-check on a
material where both exist, and nothing in `models/siesta.py`, `bases.py` or the oracle may
import or require it. A `covariant` velocity on the SIESTA route means the dipole matrix
computed by quadrature of the atomic orbitals (section 3), and the self-check of section 8
(rotated Kane-Mele with `A' = T^+ R T`) is what proves that formula; it needs no Wannier
input either.

## 12. The metric forms of every operator, derived with sympy (user, 2026-10-22)

Every physical operator of this stage (velocity in its three rungs, spin, spin current, and
the words built from them) is to be written in the three metric forms and the code must
carry all three, selected by the run file's `route`:

- **generalised**: operators as `S^-1 O` on the ket side with the asymmetric bra of
  `theory.md` equations 2a to 2c (the matrix the engine applies is not Hermitian, the
  trace is);
- **Löwdin**: `O_L = S^-1/2 O S^-1/2` with `H_L = S^-1/2 H S^-1/2`; Hermitian matrices, the
  stage-1 engine unchanged; the user's preferred default for this reason;
- **symmetrised anticommutator**: `1/2 {S^-1, O}` and its products, Hermitian when `O` is,
  acceptable where the trace identity holds; state exactly for which words it is equal to
  the other two and for which it is only an approximation (it is not a similarity
  transform, so this must be shown, not assumed).

Derive them with sympy (noncommutative symbols for `H`, `S`, `S^-1`, `S^-1/2`, `V`,
`sigma`, with the rules `S^-1/2 S^-1/2 = S^-1`, `S S^-1 = 1`, Hermitian conjugation as a
transformation) in `python/lsquant/metric_forms.py`: a function per operator that returns
the three forms as sympy expressions, a check that each form is Hermitian or not, and the
proof (symbolic expansion to the same normal form) that the traces `Tr[A f(H) B g(H)]` of
the generalised and Löwdin forms coincide for the words of the catalogue; the sympy output
is then what the numerical code implements (the formulas in docstrings cite the sympy
function), as the spectral functions already go from sympy to coefficient tables
(AGENTS.md section 4). The self-check prints the three forms of `V`, `Sigma_z` and `J` and
the Hermiticity verdicts, and runs the rotated Kane-Mele plateau on both the generalised
and the Löwdin route (2.000 on both) and on the anticommutator form (report the number;
PASS only if the derivation says it must be equal).
