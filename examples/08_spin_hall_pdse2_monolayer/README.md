# Example 8: Spin Hall conductivity of a pentagonal $\mathrm{PdSe_2}$ monolayer from SIESTA

Monolayer $\mathrm{PdSe_2}$ crystallises in a puckered pentagonal Cairo-tiling
lattice (space group $P2_1/c$ in bulk, nonsymmorphic $Pb2_1a$ / $P2_1/c$ symmetry
in monolayer form) with six atoms per unit cell (two $\mathrm{Pd}$ and four
$\mathrm{Se}$) and $156$ spinor atomic orbitals in the SIESTA double-$\zeta$
polarised basis. Its low-symmetry crystal field and $4d$/$4p$ spin–orbit coupling
produce strongly anisotropic valence and conduction bands separated by an
insulating gap around $E \approx 1.05\ \mathrm{eV}$ in the stored SIESTA energy
reference, where $\sigma^z_{xy}(E)$ vanishes inside the trivial gap and develops
sharp spin-Berry-curvature peaks at the valence and conduction band edges.

## The physics and non-orthogonal formulation

As in Example 7, the SIESTA Hamiltonian $H(\mathbf{R})$ and overlap $S(\mathbf{R})$
define a generalised eigenproblem $H(\mathbf{k})\,\mathbf{c}_n(\mathbf{k}) = \varepsilon_n(\mathbf{k})\,S(\mathbf{k})\,\mathbf{c}_n(\mathbf{k})$:

1. **Metric routes (`metric.route`):**
   - **`lowdin` (default):** symmetric orthogonalisation
     $H' = S^{-1/2} H S^{-1/2}$ and $O' = S^{-1/2} O S^{-1/2}$.
   - **`generalised`:** one-sided contragredient recursion on
     $\tilde{H} = S^{-1}(a H + b S) = a S^{-1} H + b I$.
2. **Velocity gauges (`model.params.velocity_mode`):**
   - **`bare`**, **`berry_connection`**, and **`covariant`** (which integrates the
     two-centre orbital dipole matrix elements $\mathbf{A}(\mathbf{R})$ from
     `Pd.ion.xml` and `Se.ion.xml` via certified Gauss–Legendre quadrature, cached
     in `.lsquant-cache/` outside `siesta/`).
3. **Observable words:**
   - Density of states: `Tr[delta(E-H)]`
   - Intrinsic spin Hall conductivity: `Tr[JXSZ delta(E-H) VY dG-(E)/dE]` in units
     of $e / (4\pi)$.

## Step 1: verify the SIESTA inputs

The committed SIESTA run files live under [`siesta/`](siesta/) (`pdse2.fdf`,
`pdse2.TSHS`, `Pd.ion.xml`, `Pd.psml`, `Se.ion.xml`, `Se.psml`) and are never
modified at run time:

```bash
(cd examples/08_spin_hall_pdse2_monolayer/siesta && sha256sum -c SHA256SUMS)
```

## Step 2: run the quick or production calculation

For interactive testing and CI, run the `12x12x1` ($M = 256$, `berry_connection` +
`lowdin`) quick configurations into a single output directory:

```bash
python3 -m lsquant run examples/08_spin_hall_pdse2_monolayer/run_dos_quick.json --out /tmp/pdse2_quick
python3 -m lsquant run examples/08_spin_hall_pdse2_monolayer/run_shc_quick.json --out /tmp/pdse2_quick
```

For converged research curves with certified radial dipole quadrature (`48x48x1`,
$M = 1024$, $\eta = 10\ \mathrm{meV}$, `covariant` + `lowdin`):

```bash
python3 -m lsquant run examples/08_spin_hall_pdse2_monolayer/run_dos_production.json --out /tmp/pdse2_prod
python3 -m lsquant run examples/08_spin_hall_pdse2_monolayer/run_shc_production.json --out /tmp/pdse2_prod
```

## Step 3: plot the DOS and spin Hall conductivity

Generate the two-panel figure comparing the KPM reconstruction against the exact
eigenbasis oracle:

```bash
python3 examples/08_spin_hall_pdse2_monolayer/plot.py /tmp/pdse2_quick
```

![Density of states and intrinsic spin Hall conductivity of the pentagonal PdSe2 monolayer](fig_dos_shc.png)

Panel (a) shows the density of states per unit cell across the valence and
conduction bands and the semiconducting gap around $E \approx 1.05\ \mathrm{eV}$.
Panel (b) shows the intrinsic spin Hall conductivity $\sigma^z_{xy}(E)$ in units
of $e/(4\pi)$ alongside the exact occupied-band Berry-curvature oracle
(`exact_shc.json`).

## Output files and figures

| File | Description |
|---|---|
| `SpectralOp1pdse2KPM_M*` (`.chebmom`, `.chebmom1D`) | Rank-1 Chebyshev moments for `Tr[delta(E-H)]` |
| `NonEqOpJXSZ-VYpdse2KPM_M*` (`.chebmom`, `.chebmom2D`) | Rank-2 Chebyshev moments for `Tr[JXSZ delta(E-H) VY dG-(E)/dE]` |
| `dos.dat` | Reconstructed density of states per cell (`states / eV / cell`) |
| `shc.dat` | Reconstructed Kubo–Bastin spin Hall conductivity $\sigma^z_{xy}(E)$ (`e / 4pi`) |
| `exact_shc.json`, `spin_chern.json` | Exact per-$k$ eigenbasis spin Berry curvature oracle and topological diagnostics |
| `basis.json`, `eigen_hermiticity.json`, `run.resolved.json` | Metric residuals, dipole quadrature certificate, and resolved run receipt |
| `fig_dos_shc.png` | Two-panel DOS and $\sigma^z_{xy}(E)$ figure with exact oracle overlay |

## References

- Z. Fan, J. H. García, A. W. Cummings et al., *Linear scaling quantum transport
  methodologies*, Phys. Rep. **893**, 1–123 (2021) ([arXiv:1811.07387](https://arxiv.org/abs/1811.07387)).
- A. Oyedele et al., *$\mathrm{PdSe_2}$: Pentagonal two-dimensional layers with high
  air stability for electronics*, J. Am. Chem. Soc. **139**, 14090 (2017).
