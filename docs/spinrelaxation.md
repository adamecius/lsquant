\page spinrelaxation The spin relaxation times

**Introduction**                                                                                                                                                                                                                 
Systems described by a spin-dependent Hamiltonians, endows spins with dynamics.  Typically, crystalline systems will posses spin-dependent fields which will make the spins to precess coherently. This is the case of external magnetic fields and spin-orbit coupling fields. However, in the presence of randomness, the spin-coherence is lost due to the irreversibility, making it relax over some relevant time scale \f$\tau_{\rm s}\f$ know as spin relaxation time. 

To obtain the spin relaxation time it is sufficient to solve the Schrodinger equation:
\f[
i\hbar \frac{d}{dt}|\Psi(t)\rangle = H |\Psi\rangle
\f]
where \f$H\f$ the hamiltonian of the system, \f$\hbar\f$ the Planck's constant, and \f$|\Psi(t)\rangle\f$ the state of the system. Then, spin density is computed in the Heisenberg picture as
\f[
S( t) = {\rm Tr}\left[ \rho\frac{s(t)}{\Omega} \right],
\f]
where $\f\rho\f$ the density matrix, \f$s(t)=U^\dagger(t) s U(t)\f$ the time evolved spin operator, with \f$s\f$ its static representation in the Schrodinger picture and \f$U(t)\f$ the time-evolution operator. 

In LSQUANT we deal with tight-binding models of solids described by a time-independent Hamiltonian. Therefore, the evolution operator can be written simply as \f$U(t) = {\rm e}^{-i H t/\hbar}\f$, and the density matrix depends on the Fermi energy \f$\varepsilon_{\rm F}\f$. To compute the evolution of the spin density, we follow the approach done by [Cummings et. al][1], in which the density matrix takes the following form

\f[
\rho(\varepsilon_{\rm F}) = \frac{ P \delta(H-\varepsilon_{\rm F}) P }{ {\rm Tr}[\delta(H-\varepsilon_{\rm F})] }
\f]

with \f$ P \f$ a projector operator.  In [1], the projector operator is choosen such that spin density acts on spin polarizes in a given spatial dimension. If such spatial direction is described by an altitude angle \f$\theta\f$ and a azimutal angle \f$\phi\f$, then the spinorial component of the projector operator takes the form of:
\f[
P_{\pm}(\theta,\phi) =\frac{1}{2} \left(\begin{array}{cc}
1\pm \cos(\theta) & {\rm e}^{-i\phi}\sin\theta \\
{\rm e}^{i\phi}\sin\theta &1\mp \cos(\theta)
\end{array} \right)
\f]
but an arbitrary operator can be choose.  

**The Chebyshev approach**

We rewrite the expression by using the permutation property of the trace, \f[
S( t) = \frac{ {\rm Tr}\left[  PU^\dagger(t) sU(t) P\delta(H-\varepsilon_{\rm F})\right]}{Z},
\f]
where we had defined \f$Z={\rm Tr}\left[  P\delta(H-\varepsilon_{\rm F})\right]\f$. The trace is then approximate by a mean over a random phase vector
\f[
\mu_{m,n} = \frac{1}{Z}\langle \chi_{\rm L}(n)|  s|\chi_{\rm R}(n,m)\rangle,
\f]
where \f$ |\chi_{\rm L}(n) \rangle =  U(t_n)P|\chi\rangle \f$ and \f$ |\chi_{\rm R}(m,n) \rangle =  U(t_n)PT_m(H)|\chi\rangle \f$

Finally, the \subpage time_evolution  is computed in terms of the expansion moments as:
\f[
S(t_n)= \sum_{m}\mu_{m,n} T_m(\varepsilon_{\rm F})
\f]



[1]: https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.119.206601




