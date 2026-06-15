#include "observable.hpp"
#include "chebyshev_coefficients.hpp"   // delta_chebF, greenR_chebF, DgreenR_chebF
#include "recon_grid.hpp"               // chebyshev::safety_factors().recon_cutoff (alpha)

namespace lsquant
{
	// The two formulas, as data. Greenwood = N x 1 (Fermi-surface) degenerate case of Bastin.
	const KuboObservable KUBO_GREENWOOD = { "Kubo-Greenwood", &greenR_chebF,  -1.0,  1, false };
	const KuboObservable KUBO_BASTIN    = { "Kubo-Bastin",    &DgreenR_chebF, -2.0, 30, true  };

	std::vector<std::pair<double,double> >
	reconstruct_kubo(chebyshev::Moments2D& mu, const KuboObservable& obs)
	{
		const int    num_div = obs.num_div_mult * mu.HighestMomentNumber();
		const double xbound  = chebyshev::safety_factors().recon_cutoff;   // alpha

		std::vector<double> energies(num_div, 0.0);
		for (int i = 0; i < num_div; ++i)
			energies[i] = -xbound + i * (2 * xbound) / (num_div - 1);

		// sigma kernel per energy: prefactor * N * ScaleFactor^2 * sum_{m,n} delta(E,m) Im[green(E,n) mu_{mn}]
		std::vector<double> kernel(num_div, 0.0);
		#pragma omp parallel for
		for (int i = 0; i < num_div; ++i)
		{
			double out = 0.0;
			const double energ = energies[i];
			for (int m0 = 0; m0 < mu.HighestMomentNumber(0); ++m0)
			for (int m1 = 0; m1 < mu.HighestMomentNumber(1); ++m1)
				out += delta_chebF(energ, m0) * ( obs.green(energ, m1) * mu(m0, m1) ).imag();
			kernel[i] = obs.prefactor * out * ( mu.SystemSize() * mu.ScaleFactor() * mu.ScaleFactor() );
		}

		// Output on the physical-energy axis. Bastin integrates the kernel cumulatively
		// (Fermi sea); Greenwood reports the kernel directly (Fermi surface). Both emit num_div-1.
		std::vector<std::pair<double,double> > result;
		result.reserve(num_div > 0 ? num_div - 1 : 0);
		double acc = 0.0;
		for (int i = 0; i < num_div - 1; ++i)
		{
			const double e_phys = energies[i] / mu.ScaleFactor() - mu.ShiftFactor();
			if (obs.integrate) { acc += kernel[i] * (energies[i+1] - energies[i]); result.push_back(std::make_pair(e_phys, acc)); }
			else               { result.push_back(std::make_pair(e_phys, kernel[i])); }
		}
		return result;
	}
}
