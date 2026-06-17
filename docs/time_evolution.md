\page time_evolution The time evolution 

**Introduction**                                                                                                                                                                                                                 
Given an arbitrary state \f$|\Psi\rangle\f$, its evolution in time is given by the Schrodinger equation:
\f[
i\hbar \frac{d}{dt}|\Psi(t)\rangle = H |\Psi\rangle
\f]
where \f$H\f$ the hamiltonian of the system, \f$\hbar\f$ the Planck's constant, and \f$|\Psi(t)\rangle\f$ the state of the system at a given time \f$t\f$.  For time-independent Hamiltonians, the time-evolve state has a simple expression:

\f[
|\Psi(t)\rangle = U(t)|\Psi\rangle
\f]

where \f$U(t)\equiv{\rm e}^{-i H t}\f$ an unitary operator known as time-evolution operator. Typically, the calculation of the exponential of an operator is not a trivial task. In LSQUANT we circumvent this problem by expanding it in terms of Chebyshev polynomials. Following the approach in [1], we first express the evolution operator in terms of the normalized hamiltonian \f$\tilde{H}\f$
\f[
U(t)={\rm e}^{-i \omega_c t} {\rm e}^{-i \omega \tilde{H}}\
\f]
where  \f$\omega_c\equiv E_c/\hbar\f$ with \f$E_c\f$ the band center, and \f$\omega\equiv\frac{W}{2\hbar}\f$ with  \f$W\f$ the bandwidth. 

\f[
U(t) = \sum_{m=0}^{\infty} U_m(t) T_m(\tilde{H})
\f]
where 
\f[
U_m(t)  = (-i)^m(2-\delta_{m0})J_m(\omega t){\rm e}^{-i \omega_c t}
\f]





[1]: LSQT
[1]: https://journals.aps.org/rmp/abstract/10.1103/RevModPhys.78.275


