# Stage 6d developer report

Branch: `handoff/0006d-dev`, from `5c54542`. This is a non-blocking side
branch. No build, Python execution, numerical check, golden gate or figure
was run here, as instructed. All numerical assertions and material results
below are tasks for the tester on Aron, not developer measurements.

The checkout's `.git` is read-only (`index.lock`: read-only file system).
Commits were therefore made using a copy of that metadata at
`/tmp/lsquant-stage6d.git`, with the original checkout as work tree, retaining
the existing branch name. The deliverable is `/tmp/lsquant-stage6d.bundle`;
the checkout's original HEAD cannot be advanced here. Nothing was pushed.

## Source and implementation

The first commit preserves all supplied Python files under
`tools/downfold_origin/` unchanged, with README provenance and SHA256SUMS.
The source was supplied from `/root/scratch/gb` on Ionos; its upstream commit
was not supplied. No separate downfolder notes or PDFs were present. Read
in full: both downfolders, `siesta2w90_sisl.py`, `study.py`, the 6d brief,
AGENTS.md, the stage-6 theory/literature notes, the 6b brief, its dated notes
and reviews 1–4. The supplied handoff documents are a separate commit.

Files added or extended:

- `python/lsquant/metric_forms.py`: symbolic projection expressions,
  arbitrary-dimension noncommutative identities, polar and Fourier forms.
- `python/lsquant/projection/spectral.py`: checked physical matrices,
  generalised eigensolve, S-projector, hard/shift pencil, projected
  observables, generic ordered-word trace and error diagnostics.
- `projection/origin.py`: lazy loading of the unchanged origin helper and
  subprocess entry point for either complete original reducer.
- `projection/orthogonal.py`: isolated-manifold SCDM trial selection,
  original polar helper, S-isometry, physical observable reduction and
  deterministic off-grid band comparison.
- `projection/real_space.py`: origin FFT/WS interpolation, gauge/order
  conversion, radius errors and periodic cell-gauge CSR export.
- `projection/bases.py`: projection adapters with independent caches,
  original overlap appliers, selected-identity handling and CSR interface.
- `projection/words.py`: preserve an explicit projected identity through
  the existing compiler, which normally simplifies literal `1` away.
- `projection/runfile.py`: clean SIESTA executor, physical word oracles,
  per-word errors, timings, bounds/counts and interpolation receipts.
- `projection/diagnostics.py`: exact word errors over a sequence of windows.
- `python/lsquant/runfile.py`: additive configuration field, validation,
  resolved configuration and dispatch; the `none` executor is unchanged.
- `runs/projection/{bi,pdse2}_shift.json`: editable material-study examples,
  outside `runs/goldens`; no reference values are embedded in calculations.
- `tests/dev/stage_6d_check.py`: the one explicitly requested self-check.

`bases.py`, `models/bloch.py` and `oracle/spin_hall.py` were not edited. The
C++ engine, compatibility code, golden harness, goldens and figures were
not edited. No dependencies beyond the requested NumPy 2.2, SciPy 1.15,
sisl 0.16 environment and the existing SymPy/core extension were added.
No sub-agents or remote commands were used.

## Derivations and conventions

`metric_forms.spectral_projection_forms()` is the source for the numerical
forms. Let `C` be the complete generalised eigenframe, `D=diag(E)`,
`C† S C=I`, `C C†=S^-1`, and `W=diag(0 or 1)` the inclusive window mask.
Then

```
P = C W C† S
P† = S C W C†
P² = P,                 P† S = S P
H = S C D C† S
H_P = S C W D C† S
H_s = S C D_s C† S
    = H + S C (D_s-D) C† S
D_s = diag(clip(E, E_lo, E_hi))
O_P = P† O P = S P (S^-1 O) P
```

Here `O` is a stored PHYSICAL matrix, including a velocity, spin or spin
current. The brief's shorthand `S P O P` is valid when its `O` means the
coefficient action `O_c=S^-1 O_physical`. Applying that shorthand directly
to the stored physical matrix would insert an extra S and is incorrect.
The implementation uses `P† O_physical P`. Physical identity is `S`, so
its projection is `S P`. Every current is projected as an existing physical
matrix, rather than reconstructed as an anticommutator of projected
velocity and projected spin (which inserts an extra internal projector).

`spectral_projection_proof()` expands and reduces these identities using
`C† S C=I`, completeness, `W²=W`, `WD=DW`, and `D_s W=D W`. It returns
residual expressions for idempotency, S-self-adjointness, both retained
pencils, hard complement zero, equivalent operator representations,
projected identity, the hard coefficient form, the shift update and
physical Hermiticity. The tester must verify every residual is zero.
The original metric forms `S^-1 O_P` and `S^-1/2 O_P S^-1/2` then feed the
unchanged exact generalised and Löwdin engine plans, respectively. The
approximate anticommutator route is refused for this spectral tool.

Hard projection keeps the original matrix dimension and S. All discarded
states become EXACT ZERO modes of the generalised pencil. Löwdin sees
these as zero eigenvalues of `S^-1/2 H_P S^-1/2`; the generalised route sees
the same zero eigenvalues of `S^-1 H_P`. In particular a window `[20,40]`
still has complement states at zero, outside that window. Its actual bounds
are reported, and the expansion interval must enclose them. Neither
metric route removes modes or inserts a second shift.

Shift preserves every in-window eigenpair and moves lower/upper discarded
states to the respective window edge. Its spectrum is contained in the
requested interval before local truncation. A window boundary cutting a
near-degenerate multiplet remains a physical selection ambiguity: the
selection is inclusive and literal, never silently enlarged. The oracle
retains the 6b degeneracy policy (default 1e-6 eV and residual floor).

`Tr[delta(E-H)]` counts all zero/edge complement modes. To count only kept
states, use `Tr[P:WINDOW delta(E-H)]`, or use the projection-aware
`compute_projection_moments` with `Tr[1 delta(E-H)]`. That wrapper renames
an explicit identity slot to `P:WINDOW` before calling the existing
compiler. The reported word retains the user's original text. Direct use
of the generic 6b compiler simplifies literal `1` and is therefore not the
selected-DOS API. This distinction is included in receipts.

`orthogonal_projection_forms()` states `Q=C_W U`, `Q† S Q=I`,
`H_W=U† D_W U` and `O_W=Q† O Q`, with the origin polar `U=u v†`. QRCP of
the k-averaged Löwdin subspace density selects deterministic SCDM columns,
as in the original isolated SCDM path. This optional tool requires constant
window rank and full-rank polar trials. It refuses variable rank or trial
singular values below 1e-8 rather than fabricating missing directions.
The original shell, frozen-window and SCDM scans remain available through
`run_origin_downfolder(..., method='polar'|'scdm')`; both original programs
perform their own off-grid band validation. The automatic run-file
orthogonal path uses the narrower isolated-window array tool, not the
original broad candidate search. Off-grid error is recorded, not blessed.

## Local representation and 6c boundary

Origin `Model` is reused without source rewriting: its negative-sign FFT,
WS nearest images, boundary-image degeneracies and positive-sign forward
sum are unchanged. Its reshape assumes the last mesh coordinate varies
fastest; LSQUANT's first coordinate varies fastest. The adapter validates
integer Gamma mesh coordinates and sorts into the origin order.

With `D_ii(k)=exp(i k.tau_i)`, the symbolic gauge identity is
`M_cell=D M_atom D†`. Convert to cell gauge before the origin FFT and
back afterward. Radius cuts use the origin orbital-pair hopping distance
in Angstrom. Full and radius-cut reconstruction errors are measured
separately on the training mesh. When a run requests a radius, fresh
per-k projectors and physical operators at 16 deterministic off-grid points
also measure interpolation and total truncation errors. Orthogonal runs
instead compare reduced bands at 40 independent points, seed 1.

`ProjectedBlochBasis.local_matrix(name)` builds local matrices lazily.
`local_matrix('1')` stores the physical selected identity `S P`; its
coefficient projector is the composite `S^-1 (S P)`. Truncated local
projectors need not remain idempotent. `supercell_operators(names)` exports
periodic CSR H, original S, and projected observables, summing duplicate
WS images and sorting CSR indices. S has no radius cut. With a finite mesh
this is the periodic representation of the sampled original S, not a
claim to reconstruct off-grid S exactly.

With `radius` absent, the engine uses exact per-k projected matrices and
skips local FFT preprocessing on hard/shift. A supplied finite radius
uses local truncated matrices; the receipt then reports ACTUAL eigenvalue
bounds, which can extend beyond the selection window. Truncation is a
measured approximation, not silently certified spectral clipping.
The unchanged origin WS algorithm stores orbital-pair/image distances;
its memory grows with `N_k W²` times the WS image count. Large material
local exports should be measured separately from the per-k gate.

The 6c caller can compose exported physical CSR operators with its own
metric/disorder operators. This branch does not add disorder or alter 6c
files. Adding disorder after clean projection can couple the subspaces;
its error needs a separate measurement. Run-file `basis.name='real'` and
nonzero disorder are refused here; local CSR is the real-space API.

## Run interface and word comparisons

```
"projection": {"kind": "none"|"hard"|"shift"|"orthogonal",
               "window": [E_lo,E_hi], "radius": 12.0}
```

`radius` is optional; `none` takes neither a window nor a radius. Selection
energies refer to the stored SIESTA E_F=0 convention. They are distinct
from `expansion.window`, which is either `eigen` (with the existing 10%
engine padding) or an explicit enclosing interval. Replay preserves the
projection configuration and resolved expansion interval.

Projected run files support clean SIESTA Bloch/eigen bases, exact per-block
metric appliers, the physical SHC word and one-slot spectral words.
They write `.chebmom`, optional legacy layouts, curve data, original and
projected word oracles, `projection.json`, and `run.resolved.json`.
The receipt contains counts per k, resolved selection window, raw and
padded bounds, metric route, matrix dimension, projector residuals,
truncation errors, oracle word errors and phase wall times.
`reconstruction.energies` or `energy_grid` controls a common oracle grid;
otherwise the existing 6b gap-aware default is used.

For the material study, `basis.name='eigen'`, `trace='random'`, `states=16`
uses independent seeded eigenbasis phases. Exact trace uses W block-local
states. Bloch random traces retain the existing orbital Bloch estimator.
No reference KPM pass is added: the reference is the exact unprojected
word oracle. Both original and projected word oracles are always computed
because their difference is the scientific output of this tool.

`projection.diagnostics.compare_windows` reports word errors as a function
of a supplied sequence of windows, holding the source, energy grid and
word parameters fixed. `spectral.word_trace` supplies a generic exact
ordered-word diagnostic with caller-provided spectral arrays. The run
executor covers the existing stage-6b spectral/SHC scope, not rank-three
or time words. KPM/Jackson versus Gaussian oracle errors are labelled as
comparisons of different kernels and do not imply byte equivalence.

Projecting observables removes kept/discarded virtual transitions. Even
though retained eigenpairs are exact, a Hall plateau need NOT be exact for
a finite window on Bi/PdSe2. The code reports the deviation; it neither
adds compensation terms nor adjusts prefactors to recover a target. The
decoupled high-energy toy must pass, but cannot establish this material
invariance. This is the principal open physics question of the gate.

## What the Aron tester must measure

Use the requested NumPy 2.2, SciPy 1.15, sisl 0.16 environment and the built
extension. The developer has not executed even the symbolic self-check.

```
PYTHONPATH=python:build:build/lsquant_core python tests/dev/stage_6d_check.py
python -m lsquant run runs/projection/bi_shift.json
python -m lsquant run runs/projection/pdse2_shift.json
```

Self-check expectations (assertions, NOT observed output):

1. All SymPy residuals zero. Numeric P idempotent and S-self-adjoint,
   retained hard/shift eigenpairs exact to 1e-10; hard complement zero.
2. Rotated spinful Kane-Mele plus two +30/+60 eV orbitals per site:
   hard/shift oracle plateau errors below 1e-3, and both exact metric
   compiler plans agree with their diagonal-H engine moments at 1e-10.
   At equal M (default 1024), wide and shifted exact-trace KPM plateaus
   must agree with the finite-mesh unprojected oracle and each other to
   1e-3. The script prints all plateaus, ranges, and nominal fixed-M
   resolution gain `range_before/range_after`. It never changes the
   tolerance or tunes a window to force a pass. Failure of the wide-band
   KPM at that M must be reported; convergence has not been assumed from
   unavailable measurements. The nominal resolution is `pi*range/(2M)`;
   actual Chebyshev local resolution also depends on the rescaled energy.
3. Optional orthogonal reduction: training bands and S-isometry to 1e-10,
   independent off-grid bands to 1e-6. Local FFT/WS training-grid error,
   selected-identity handling, Hermitian CSR export and run-file replay
   are checked in the same script. No separate test suite was added.

For each material, run windows ±3 and ±5 eV, hard and shift, 24x24 and
48x48, M=256 and 1024, 16 eigenbasis vectors, seed 12345. Repeat for
berry_connection/covariant velocities and use the matching unprojected
oracle: Bi 1.5527/1.1747, PdSe2 2.08/1.98 (from the supplied reviews).
Those numbers are comparison targets, never calculation inputs.
PdSe2's topological gap is around +1.05 eV; its E_F gap is trivial.
The committed examples use three points inside the corresponding gap;
expand the grids if gap edges or window-error profiles are needed.
Give each variant its own output directory. Record retained-band count
range and per-k distribution, actual spectral bounds, all oracle plateau
errors, KPM plateau values, moments/reconstruction/oracle timings, and
wall time. Repeat thread counts to assess determinism in the target BLAS
environment. No parallel reductions were added by this branch; each k
has an independent Hermitian solve, QR/SVD and deterministic FFT order.

The existing golden gate and the stage figure remain tester deliverables.
No goldens were regenerated or accepted in this developer checkout.
