# Observable Catalogue

LSQUANT can compute the following physical quantities. Each one is specified by a short expression (an “observable word”) that describes the operators and spectral functions involved.

| Quantity | Expression | Typical use |
|----------|------------|-------------|
| Density of states | `Tr[delta(E-H)]` | Total DOS |
| Projected DOS | `Tr[A delta(E-H)]` | Orbital- or spin-projected DOS |
| Longitudinal conductivity | `Tr[VX delta(E-H) VX delta(E-H)]` | $\sigma_{xx}(E)$ (Kubo–Greenwood) |
| Hall / spin Hall conductivity | `Tr[VX delta(E-H) VY dG-(E)/dE]` | $\sigma_{xy}(E)$, spin Hall (with spin current) |
| Mean-square displacement | `Tr[[X,U(t)]^+ delta(E-H) [X,U(t)]]` | Diffusion, localization length |
| Spin dynamics | `Tr[P U(t)^+ S U(t) delta(E-H) P]` | Spin precession, relaxation |
| Time correlation | `Tr[A delta(E-H) U(t)^+ B U(t)]` | General dynamical correlations |

All of these are available for real-space models, Bloch Hamiltonians, and non-orthogonal bases (including SIESTA).

For the complete list of supported operators and functions see the examples in the `examples/` directory.
