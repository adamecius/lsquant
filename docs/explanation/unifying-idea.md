# The Unifying Idea

LSQUANT computes spectral and transport observables without diagonalizing the Hamiltonian. Every quantity is written as the trace of an ordered product of static operators and spectral functions of the Hamiltonian:

$$
C = \mathrm{Tr}\bigl[ A_0\, f_1(H)\, A_1\, f_2(H)\, \cdots\, f_p(H)\, A_p \bigr]
$$

The library parses this “observable word”, compiles it into a short bytecode plan, executes the plan with the Kernel Polynomial Method on any linear operator, and reconstructs the physical curve from the resulting Chebyshev moments.

This single abstraction covers:
- Density of states
- Longitudinal and Hall conductivities
- Spin Hall conductivity
- Real-time dynamics and mean-square displacement
- Calculations on real-space, Bloch, and eigen bases
- Orthogonal and non-orthogonal metrics
- First-principles SIESTA models

The same word language and engine work for all of these cases.
