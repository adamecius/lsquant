#ifndef LSQUANT_UNITS
#define LSQUANT_UNITS

namespace chebyshev
{
	// Reduced Planck constant, in eV*fs -- the SOLE home for the time-evolution unit.
	//
	// The Chebyshev/Bessel propagator advances PHYSICAL time t (in fs): the per-step
	// frequency is ChebyshevFreq() = HalfWidth/CUTOFF/HBAR, so e.g. spin precession runs at
	// omega = 2*J_ex/HBAR rad/fs (period 2*pi*HBAR/(2 J_ex)), NOT the hbar=1 convention of the
	// static moment oracle. See docs/spectral_scales.md section 5.
	//
	// This value was hard-coded x5 (twice in chebyshev_moments.hpp, plus the _Device and
	// _Graphene_NNN forks and chebyshev_solver.cpp). Consolidated here; the forks include this
	// header now and the duplicate copies disappear when the forks merge in Phase 6.
	inline constexpr double HBAR = 0.6582119624; // eV*fs
}

#endif // LSQUANT_UNITS
