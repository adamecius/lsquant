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

