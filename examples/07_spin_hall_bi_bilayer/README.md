# Example 7: Spin Hall conductivity of a 2D $\mathrm{Bi}(111)$ bilayer from SIESTA

Bismuth is the heaviest stable group-V element, and its strong atomic spin–orbit
coupling (SOC) turns a buckled honeycomb $\mathrm{Bi}(111)$ bilayer into a prototypical
two-dimensional quantum spin Hall system. Unlike a toy orthogonal tight-binding
model, a first-principles SIESTA calculation expands the wavefunctions in a finite-range
numeric atomic-orbital (NAO) basis whose basis functions on neighbouring atoms
overlap, so the secular equation is the generalised eigenproblem

$$ H(\mathbf{k})\,\mathbf{c}_n(\mathbf{k}) = \varepsilon_n(\mathbf{k})\,S(\mathbf{k})\,\mathbf{c}_n(\mathbf{k}). $$

This example computes the density of states and the intrinsic Kubo–Bastin spin Hall
conductivity $\sigma^z_{xy}(E)$ directly from the SIESTA `.fdf`, `.TSHS`, and `.ion.xml`
files shipped under [`siesta/`](siesta/), comparing the KPM Chebyshev expansion
against the exact per-$k$ eigenbasis Berry-curvature oracle.

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
     $J_x^{s_z} = \tfrac12 \{v_x, s_z\}_S$ is the metric-consistent spin current
     operator, reconstructed in units of $e / (4\pi)$.

## Step 1: verify the SIESTA inputs

The committed SIESTA run files live under [`siesta/`](siesta/) (`bi.fdf`,
`bi.TSHS`, `Bi.ion.xml`, `Bi.psml`) and are never modified at run time:

```bash
(cd examples/07_spin_hall_bi_bilayer/siesta && sha256sum -c SHA256SUMS)
```

## Step 2: run the quick or production calculation

For interactive testing and CI, run the `12x12x1` ($M = 256$, `berry_connection` +
`lowdin`) quick configurations into a single output directory:

```bash
python3 -m lsquant run examples/07_spin_hall_bi_bilayer/run_dos_quick.json --out /tmp/bi_quick
python3 -m lsquant run examples/07_spin_hall_bi_bilayer/run_shc_quick.json --out /tmp/bi_quick
```

For converged research curves with certified radial dipole quadrature (`48x48x1`,
$M = 1024$, $\eta = 10\ \mathrm{meV}$, `covariant` + `lowdin`):

```bash
python3 -m lsquant run examples/07_spin_hall_bi_bilayer/run_dos_production.json --out /tmp/bi_prod
python3 -m lsquant run examples/07_spin_hall_bi_bilayer/run_shc_production.json --out /tmp/bi_prod
```

## Step 3: plot the DOS and spin Hall conductivity

Generate the two-panel figure comparing the KPM reconstruction against the exact
eigenbasis oracle:

```bash
python3 examples/07_spin_hall_bi_bilayer/plot.py /tmp/bi_quick
```

![Density of states and intrinsic spin Hall conductivity of the 2D Bi(111) bilayer](fig_dos_shc.png)

Panel (a) shows the SOC-split density of states per unit cell around the Fermi
level $E_F = 0$. Panel (b) shows the intrinsic spin Hall conductivity
$\sigma^z_{xy}(E)$ in units of $e/(4\pi)$ together with the exact occupied-band
Berry-curvature oracle (`exact_shc.json`).

## Output files and figures

| File | Description |
|---|---|
| `SpectralOp1biKPM_M*` (`.chebmom`, `.chebmom1D`) | Rank-1 Chebyshev moments for `Tr[delta(E-H)]` |
| `NonEqOpJXSZ-VYbiKPM_M*` (`.chebmom`, `.chebmom2D`) | Rank-2 Chebyshev moments for `Tr[JXSZ delta(E-H) VY dG-(E)/dE]` |
| `dos.dat` | Reconstructed density of states per cell (`states / eV / cell`) |
| `shc.dat` | Reconstructed Kubo–Bastin spin Hall conductivity $\sigma^z_{xy}(E)$ (`e / 4pi`) |
| `exact_shc.json`, `spin_chern.json` | Exact per-$k$ eigenbasis spin Berry curvature oracle and topological diagnostics |
| `basis.json`, `eigen_hermiticity.json`, `run.resolved.json` | Metric residuals, dipole quadrature certificate, and resolved run receipt |
| `fig_dos_shc.png` | Two-panel DOS and $\sigma^z_{xy}(E)$ figure with exact oracle overlay |

## References

- Z. Fan, J. H. García, A. W. Cummings et al., *Linear scaling quantum transport
  methodologies*, Phys. Rep. **893**, 1–123 (2021) ([arXiv:1811.07387](https://arxiv.org/abs/1811.07387)).
- J. H. García, M. Vila, C.-H. Hsu, X. Waintal, V. M. Pereira, and S. Roche,
  *Canted persistent spin texture and quantum spin Hall effect in $\mathrm{WTe_2}$*,
  Phys. Rev. Lett. **125**, 256603 (2020).
