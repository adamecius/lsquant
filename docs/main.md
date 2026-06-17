# LSQUANT: A tool for linear quantum transport

LSQUANT is a tools which allows performing quantum transport calculations in tight-binding models. LSQUANT possess interfaces with other tight-binding modelling packages such as **Wannier90,** or **KWANT** and Tbmodel. Its core is written in C++ and uses openMP and MPI for optimized calculations.

LSQUANT allows for the following kind of calculations:

1. Nonequilibrium electrical response of different operators within the linear regime. Some particular examples are:
  1. The conductivity.
  2. The spin Hall conductivity.
  3. The spin susceptibility.
2. Time-correlation which are the basis for computing:
  1. The Mean-square spreading.
  2. \subpage spinrelaxation.

## Getting started

Build LSQUANT with CMake (`cmake -B build && cmake --build build -j`); the
installer helper is `install.py` at the repository root, and the full
prerequisites and build options are described in the top-level `README.md`. This
page does not repeat those steps so they cannot drift out of sync.

Worked, runnable calculations live in the `examples/` directory, one
self-contained subdirectory per topic (1D density of states, localization,
Haldane Hall, spin precession, charge transport, and spin relaxation). Each
example carries its own `README.md` explaining the physics it demonstrates and
how to run it.

