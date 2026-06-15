#include "observable.hpp"
#include "chebyshev_coefficients.hpp"   // delta_chebF, greenR_chebF, DgreenR_chebF
#include "recon_grid.hpp"               // chebyshev::safety_factors().recon_cutoff (alpha)
#include <cmath>                        // std::cos, M_PI (FFT node grid)

namespace lsquant
{
	// The two formulas, as data. Greenwood = N x 1 (Fermi-surface) degenerate case of Bastin.
	const KuboObservable KUBO_GREENWOOD = { "Kubo-Greenwood", &greenR_chebF,  -1.0,  1, false };
	const KuboObservable KUBO_BASTIN    = { "Kubo-Bastin",    &DgreenR_chebF, -2.0, 30, true  };

	namespace {
		// FFT node ordering -> ascending energy (needed for the cumulative Bastin integral).
		// Verbatim from the legacy *_FFTgrid drivers so the output is byte-identical.
		void rearrange_crescent_order(std::vector<double>& a)
		{
			const int nump = (int)a.size();
			std::vector<double> orig(a.begin(), a.end());
			for (int k = 0; k < nump / 2; ++k) { a[2*k] = orig[k]; a[2*k+1] = orig[nump-k-1]; }
			for (int k = 0; k < nump / 2; ++k) { double t = a[k]; a[k] = a[nump-k-1]; a[nump-k-1] = t; }
		}
	}

	std::vector<std::pair<double,double> >
	reconstruct_kubo(chebyshev::Moments2D& mu, const KuboObservable& obs, ReconGrid grid)
	{
		// --- grid backend: uniform [-alpha,alpha] (30*M pts) or FFT Chebyshev nodes (M pts) ---
		std::vector<double> energies;
		if (grid == GRID_UNIFORM)
		{
			const int    num_div = obs.num_div_mult * mu.HighestMomentNumber();
			const double xbound  = chebyshev::safety_factors().recon_cutoff;   // alpha
			energies.resize(num_div);
			for (int i = 0; i < num_div; ++i)
				energies[i] = -xbound + i * (2 * xbound) / (num_div - 1);
		}
		else // GRID_FFT_NODES
		{
			const int num_div = mu.HighestMomentNumber();
			const double disp = 0.5;
			energies.resize(num_div);
			for (int i = 0; i < num_div; ++i)
				energies[i] = std::cos( M_PI * (2 * i + disp) / num_div );
			if (obs.integrate) rearrange_crescent_order(energies);  // ascending for the integral
		}
		const int num_div = (int)energies.size();

		// --- shared integrand: prefactor * N * (1/HalfWidth)^2 * sum_{m,n} delta(E,m) Im[green(E,n) mu] ---
		// The extensive prefactor is expressed exactly as each legacy backend wrote it (the two
		// forms are mathematically equal at CUTOFF=1, but kept distinct for byte-identity).
		const double pref_ext = (grid == GRID_UNIFORM)
			? ( mu.SystemSize() * mu.ScaleFactor() * mu.ScaleFactor() )
			: ( mu.SystemSize() / mu.HalfWidth() / mu.HalfWidth() );
		std::vector<double> kernel(num_div, 0.0);
		#pragma omp parallel for
		for (int i = 0; i < num_div; ++i)
		{
			double out = 0.0;
			const double energ = energies[i];
			for (int m0 = 0; m0 < mu.HighestMomentNumber(0); ++m0)
			for (int m1 = 0; m1 < mu.HighestMomentNumber(1); ++m1)
				out += delta_chebF(energ, m0) * ( obs.green(energ, m1) * mu(m0, m1) ).imag();
			kernel[i] = obs.prefactor * out * pref_ext;
		}

		// --- output: Bastin integrates cumulatively (Fermi sea), Greenwood reports the kernel ---
		std::vector<std::pair<double,double> > result;
		result.reserve(num_div > 0 ? num_div - 1 : 0);
		double acc = 0.0;
		for (int i = 0; i < num_div - 1; ++i)
		{
			const double e_phys = (grid == GRID_UNIFORM)
				? ( energies[i] / mu.ScaleFactor() - mu.ShiftFactor() )
				: ( energies[i] * mu.HalfWidth()   + mu.BandCenter() );
			if (obs.integrate) { acc += kernel[i] * (energies[i+1] - energies[i]); result.push_back(std::make_pair(e_phys, acc)); }
			else               { result.push_back(std::make_pair(e_phys, kernel[i])); }
		}
		return result;
	}

	std::vector<std::pair<double,double> >
	reconstruct_density_grid(const std::vector<double>& mu_real, double alpha, int num_div)
	{
		std::vector<double> energies(num_div, 0.0);
		for (int i = 0; i < num_div; ++i)
			energies[i] = -alpha + i * (2 * alpha) / (num_div - 1);

		std::vector<std::pair<double,double> > out(num_div);
		const int M = static_cast<int>(mu_real.size());
		#pragma omp parallel for
		for (int i = 0; i < num_div; ++i)
		{
			double s = 0.0;
			const double energ = energies[i];
			for (int m = 0; m < M; ++m)
				s += delta_chebF(energ, m) * mu_real[m];
			out[i] = std::make_pair(energ, s);
		}
		return out;
	}
}
