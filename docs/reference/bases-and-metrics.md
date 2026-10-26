# Bases and Metrics

LSQUANT can represent the Hamiltonian in three ways and handle both orthogonal and non-orthogonal bases.

## Representations

| Name | Description | When to use |
|------|-------------|-------------|
| `real` | Real-space sparse supercell | Large disordered systems, localization, real-time dynamics |
| `bloch` | Bloch $k$-space Hamiltonian | Clean crystals, with optional real-space disorder |
| `eigen` | Exact diagonalization of the unit cell (or small system) | Reference calculations, validation, small unit cells |

## Orthogonal vs non-orthogonal

- **Orthogonal**: standard tight-binding models (overlap matrix = identity).
- **Non-orthogonal**: first-principles bases (SIESTA, Wannier). Two common treatments:
  - **Löwdin**: orthogonalize the basis with $S^{-1/2}$.
  - **Generalised**: keep the non-orthogonal metric and use the appropriate velocity and current operators.

In the JSON file:

```json
"basis": { "name": "real" },
"metric": { "policy": "orthogonal" }
```

or for SIESTA:

```json
"basis": { "name": "eigen" },
"metric": { "policy": "nonorthogonal", "route": "lowdin" }
```

## Velocity operators (SIESTA)

When loading SIESTA Hamiltonians you can choose how the velocity operators are constructed:

- `bare`: lattice derivative of the Hamiltonian,
- `berry_connection`: includes the intracell Berry connection,
- `covariant`: full orbital dipole contribution from the `.ion.xml` files.

The choice affects the accuracy of the Hall and spin Hall conductivities.
