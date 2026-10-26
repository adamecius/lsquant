# Observable Word Catalogue

Every spectral or transport quantity in LSQUANT is expressed as a trace of an ordered word of operators and spectral functions.

| Shape ID | Canonical word | Rank | Typical use |
|----------|----------------|------|-------------|
| `dos_identity` | `Tr[delta(E-H)]` | 1 | Density of states |
| `dos_operator` | `Tr[A delta(E-H)]` | 1 | Projected or operator-weighted DOS |
| `kubo_greenwood` | `Tr[VX delta(E-H) VX delta(E-H)]` | 2 | Longitudinal conductivity |
| `kubo_bastin` | `Tr[VX delta(E-H) VY dG-(E)/dE]` | 2 | Hall conductivity, spin Hall |
| `msd` | `Tr[[X,U(t)]^+ delta(E-H) [X,U(t)]]` | 1 | Mean-square displacement |
| `projected_evolution` | `Tr[P U(t)^+ SX U(t) delta(E-H) P]` | 1 | Spin precession / relaxation |
| `time_correlation` | `Tr[A delta(E-H) U(t)^+ B U(t)]` | 1 | General time correlations |

See the [unifying idea](../explanation/unifying-idea.md) for the conceptual background.
