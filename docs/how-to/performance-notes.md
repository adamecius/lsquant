# Performance Notes

LSQUANT is designed for systems with $10^5$–$10^7$ orbitals. The cost of a calculation scales linearly with the number of orbitals and with the number of Chebyshev moments.

## Practical guidelines

- **Number of moments**: 128–256 is usually enough for densities of states and conductivities with moderate broadening (10–20 meV). Increase to 512 for sharper features.
- **Broadening**: start with 10–20 meV. Smaller values require more moments.
- **Real-space vs Bloch**: use real-space for disordered systems and dynamics; use Bloch (or eigen) for clean crystals and unit-cell reference calculations.
- **Build type**: the default build is portable and reproducible. A performance-oriented build (compiler flags `-O3 -march=native`) is available if you need maximum speed on a specific machine.

## Scaling

On multi-core machines the sparse matrix–vector products that dominate the calculation parallelize well. For very large systems, NUMA-aware placement (`numactl --interleave=all` or binding to one socket) can improve performance.

Detailed benchmark numbers are available in the repository under `bench/` and in `docs/performance/`.
