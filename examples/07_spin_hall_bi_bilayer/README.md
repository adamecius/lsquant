# Example 7: Spin Hall conductivity of a 2D $\mathrm{Bi}(111)$ bilayer from SIESTA

This example computes the density of states and the intrinsic Kubo–Bastin spin Hall
conductivity $\sigma^z_{xy}(E)$ of a two-dimensional $\mathrm{Bi}(111)$ bilayer with
spin–orbit coupling (SOC) directly from the SIESTA `.fdf`, `.TSHS`, and `.ion.xml`
files shipped under [`siesta/`](siesta/). In SIESTA's finite-range numeric
atomic-orbital (NAO) basis ($W = 72$ spinor orbitals per unit cell as recorded in
`basis.json`), neighbouring orbitals overlap, so the secular equation is the
generalised eigenproblem

$$ H(\mathbf{k})\,\mathbf{c}_n(\mathbf{k}) = \varepsilon_n(\mathbf{k})\,S(\mathbf{k})\,\mathbf{c}_n(\mathbf{k}). $$

The KPM Chebyshev expansion is compared against the exact per-$k$ eigenbasis
spin-Berry-curvature oracle written to `exact_shc.json`.

## The physics and non-orthogonal formulation

In a non-orthogonal basis with overlap matrix $S$, both the Hamiltonian recursion
and the velocity / spin-current operators depend on the metric policy:

1. **Metric routes (`metric.route`):**
   - **`lowdin` (default):** symmetric orthogonalisation
     $H' = S^{-1/2} H S^{-1/2}$ and $O' = S^{-1/2} O S^{-1/2}$, yielding an
     ordinary Hermitian problem on which the standard Chebyshev recursion runs.
   - **`generalised`:** one-sided contragredient recursion on $S^{-1} H$, rescaled
     as $\tilde{H} = S^{-1}(a H + b S) = a S^{-1} H + b I$ so that the identity
     shift preserves the generalised eigenvectors of $(H, S)$.
2. **Velocity gauges (`model.params.velocity_mode`):**
   - **`bare`:** lattice derivative $\partial_{k_\alpha} H(\mathbf{k})$ without
     intracell position or overlap connections.
   - **`berry_connection`:** includes intracell atom-position phase factors and the
     overlap Pulay connection $-\tfrac12 \{H, S^{-1} \partial_{k_\alpha} S\}$.
   - **`covariant`:** adds the intra-atomic and inter-atomic orbital dipole
     matrix elements $\mathbf{A}(\mathbf{R}) = \langle \phi_{i\mathbf{0}} | \mathbf{r} - \mathbf{R}/2 | \phi_{j\mathbf{R}} \rangle$
     integrated by certified Gauss–Legendre midpoint quadrature from SIESTA's
     radial `.ion.xml` tables (cached in `.lsquant-cache/` outside `siesta/`).
3. **Observable words:**
   - Density of states: `Tr[delta(E-H)]`
   - Intrinsic spin Hall conductivity: `Tr[JXSZ delta(E-H) VY dG-(E)/dE]`, where
     $J_x^{s_z} = \tfrac12 \{v_x, \sigma_z\}_S$ is the metric-consistent spin
     current operator, reconstructed via `shc` in units of $e^2/h$ (Pauli
     $\sigma_z$ convention).

## Step 1: verify the SIESTA inputs

The committed SIESTA run files live under [`siesta/`](siesta/) (`bi.fdf`,
`bi.TSHS`, `Bi.ion.xml`, `Bi.psml`) and are never modified at run time:

```bash
(cd examples/07_spin_hall_bi_bilayer/siesta && sha256sum -c SHA256SUMS)
```

## Step 2: run the quick or production calculation

For interactive testing, run the `24x24x1` ($M = 256$, $\eta = 20\ \mathrm{meV}$,
`covariant` velocity + `lowdin` metric) quick configurations into a single output
directory:

```bash
python3 -m lsquant run examples/07_spin_hall_bi_bilayer/run_dos_quick.json --out /tmp/bi_quick
python3 -m lsquant run examples/07_spin_hall_bi_bilayer/run_shc_quick.json --out /tmp/bi_quick
```

For the production setting (`48x48x1`, $M = 1024$, $\eta = 10\ \mathrm{meV}$,
`covariant` velocity + `lowdin` metric, matching
`runs/goldens/g11_bi_covariant_lowdin_shc.json`):

```bash
python3 -m lsquant run examples/07_spin_hall_bi_bilayer/run_dos_production.json --out /tmp/bi_prod
python3 -m lsquant run examples/07_spin_hall_bi_bilayer/run_shc_production.json --out /tmp/bi_prod
```

## Step 3: plot the DOS and spin Hall conductivity

Generate the three-panel figure comparing the KPM reconstruction against the exact
eigenbasis oracle:

```bash
python3 examples/07_spin_hall_bi_bilayer/plot.py /tmp/bi_quick
```

![Density of states and intrinsic spin Hall conductivity of the 2D Bi(111) bilayer](fig_dos_shc.png)

Panel (a) plots `dos.dat` around $E_F = 0$. Panel (b) plots the KPM Kubo–Bastin
curve `shc.dat` alongside the exact eigenbasis oracle from `exact_shc.json`, and
panel (c) shows the residual $\sigma^z_{xy,\mathrm{KPM}} - \sigma^z_{xy,\mathrm{oracle}}$.
On the `24x24x1` `covariant` + `lowdin` mesh, the exact oracle in `exact_shc.json`
(and [`handoff/0009-wrapup/tests.md`](../../handoff/0009-wrapup/tests.md)) gives
$\sigma^z_{xy}(E_F = 0) \approx 1.080\ e^2/h$; the converged production KPM value
at $E_F = 0$ is **value pending** (pending diagnosis of the reference difference
recorded in
[`handoff/0009-wrapup/profiling-results.md`](../../handoff/0009-wrapup/profiling-results.md)).

## Output files and figures

| File | Description |
|---|---|
| `SpectralOp1biKPM_M*` (`.chebmom`, `.chebmom1D`) | Rank-1 Chebyshev moments for `Tr[delta(E-H)]` |
| `NonEqOpJXSZ-VYbiKPM_M*` (`.chebmom`, `.chebmom2D`) | Rank-2 Chebyshev moments for `Tr[JXSZ delta(E-H) VY dG-(E)/dE]` |
| `dos.dat` | Reconstructed density of states per cell (`states / eV / cell`) |
| `shc.dat` | Reconstructed Kubo–Bastin spin Hall conductivity $\sigma^z_{xy}(E)$ ($e^2/h$, Pauli $\sigma_z$) |
| `exact_shc.json`, `spin_chern.json` | Exact per-$k$ eigenbasis spin Berry curvature oracle and topological diagnostics |
| `basis.json`, `eigen_hermiticity.json`, `run.resolved.json` | Metric residuals, dipole quadrature certificate, and resolved run receipt |
| `fig_dos_shc.png` | Three-panel DOS, $\sigma^z_{xy}(E)$ vs. exact oracle, and residual figure |

## References

- Z. Fan, J. H. García, A. W. Cummings et al., *Linear scaling quantum transport
  methodologies*, Phys. Rep. **893**, 1–123 (2021) ([arXiv:1811.07387](https://arxiv.org/abs/1811.07387)).
- J. H. García, M. Vila, C.-H. Hsu, X. Waintal, V. M. Pereira, and S. Roche,
  *Canted persistent spin texture and quantum spin Hall effect in $\mathrm{WTe_2}$*,
  Phys. Rev. Lett. **125**, 256603 (2020).
