# How LSQUANT Works

LSQUANT calculates spectral and transport properties of large systems without diagonalizing the Hamiltonian. Instead of finding all eigenvalues and eigenvectors, it expands the desired physical quantity in Chebyshev polynomials and evaluates the expansion using matrix–vector multiplications. The computational cost therefore grows linearly with the number of atoms (or orbitals).

Any spectral or transport observable can be written as a trace of the form

$$
C = \mathrm{Tr}\bigl[ A_0\, f_1(H)\, A_1\, f_2(H)\, \cdots\, f_p(H)\, A_p \bigr]
$$

where the $A_i$ are static operators (velocity, spin, projectors, etc.) and the $f_i(H)$ are functions of the Hamiltonian (delta functions, energy derivatives of the Green function, time evolution, etc.).

The library:
1. Parses the expression (the “observable word”),
2. Executes the necessary Chebyshev expansions,
3. Reconstructs the physical curve (density of states, conductivity, etc.).

This approach is used for density of states, longitudinal and Hall conductivities, spin Hall conductivity, real-time spin dynamics, and mean-square displacement, on both tight-binding and first-principles (SIESTA) models.

For small systems the same quantities can also be computed by exact diagonalization to validate the results.
