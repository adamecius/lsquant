#ifndef LSQUANT_OBSERVABLE
#define LSQUANT_OBSERVABLE

#include <vector>
#include <utility>
#include <complex>
#include "chebyshev_moments.hpp"

// Phase 5 -- "a formula is a piece of data, not an executable."
//
// The Kubo-Greenwood and Kubo-Bastin reconstructions share one inner double-Chebyshev sum;
// they differ ONLY in four data fields. Greenwood is the degenerate (Fermi-surface-only) case
// of Bastin: same sum, half the prefactor, no Fermi-sea integral, coarser grid. Capturing those
// four fields as a KuboObservable collapses the duplicated *FromChebmom drivers into one
// reconstruction routine; each driver becomes a thin wrapper that picks the observable and
// writes the file.

namespace lsquant
{
	struct KuboObservable
	{
		const char* name;
		// Green's-function kernel on the n-index: greenR_chebF (Greenwood) or DgreenR_chebF (Bastin)
		std::complex<double> (*green)(double, const double);
		double prefactor;      // -1.0 (Greenwood, Fermi-surface) | -2.0 (Bastin, + Fermi-sea term)
		int    num_div_mult;   // energy grid = num_div_mult * HighestMomentNumber : 1 | 30
		bool   integrate;      // false (Greenwood) | true (Bastin: cumulative Fermi-sea integral)
	};

	extern const KuboObservable KUBO_GREENWOOD;
	extern const KuboObservable KUBO_BASTIN;

	// Two grid backends under the SAME observable (the v2 plan's "FFT = its own integrand backend"):
	//   GRID_UNIFORM   : 30*M points uniform on [-alpha, alpha]; energy axis E/ScaleFactor - ShiftFactor.
	//   GRID_FFT_NODES : M Chebyshev nodes cos(pi(2i+1/2)/M) (rearranged to ascending for the Bastin
	//                    integral); energy axis E*HalfWidth + BandCenter. Coarse (M, not 30*M) -- the
	//                    legacy *_FFTgrid behaviour, kept byte-identical.
	enum ReconGrid { GRID_UNIFORM = 0, GRID_FFT_NODES = 1 };

	// Reconstruct sigma(E) for the given observable on the chosen grid backend. Returns (physical
	// energy, value) pairs (num_div-1 of them, matching the legacy drivers). Moments must already
	// be kernel-damped (ApplyJacksonKernel).
	std::vector<std::pair<double,double> >
	reconstruct_kubo(chebyshev::Moments2D& mu, const KuboObservable& obs, ReconGrid grid = GRID_UNIFORM);

	// 1-D KPM density reconstruction (DOS / spectral function): the shared inner sum
	//   rho(x) = sum_m delta_chebF(x, m) * mu_real[m]
	// on a uniform grid of `num_div` points over [-alpha, alpha]. Returns (x_rescaled, rho); the
	// caller applies its own scalar prefactor and physical-energy axis (the only per-formula
	// differences between DOS / spectral). `mu_real` must already be kernel-damped and carry the
	// KPM (2 - delta_{m0}) weighting, as the stored moments do. This is the exact primitive the
	// Phase-R oracle drives, so that oracle now guards the production reconstruction.
	std::vector<std::pair<double,double> >
	reconstruct_density_grid(const std::vector<double>& mu_real, double alpha, int num_div);
}

#endif // LSQUANT_OBSERVABLE
