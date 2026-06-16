#include "observable.hpp"
#include "chebyshev_coefficients.hpp"   // delta_chebF, greenR_chebF, DgreenR_chebF
#include "recon_grid.hpp"               // chebyshev::safety_factors().recon_cutoff (alpha)
#include <cmath>                        // std::cos, M_PI (FFT node grid)

namespace lsquant
{
	// The two formulas, as data. Greenwood = N x 1 (Fermi-surface) degenerate case of Bastin.
	const KuboObservable KUBO_GREENWOOD = { "Kubo-Greenwood", &greenR_chebF,  -1.0,  1, false };
	const KuboObservable KUBO_BASTIN    = { "Kubo-Bastin",    &DgreenR_chebF, -2.0, 30, true  };

	// --- shared reconstruction scaffolding (internal "detail" helpers) ------------------------
	// Single home for the grid construction, the energy accumulation rule, and the (energy,value)
	// assembly. reconstruct_kubo and reconstruct_density_grid route through these so a change to a
	// grid or the integral is made once; the upcoming Sea/Surf migration (A2) reuses them too.
	// Kept in the anonymous namespace -> internal linkage. NO behaviour change vs the inline code.
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

		// Uniform grid of num_div points on [-alpha, alpha].
		std::vector<double> make_uniform_grid(int num_div, double alpha)
		{
			std::vector<double> energies(num_div);
			for (int i = 0; i < num_div; ++i)
				energies[i] = -alpha + i * (2 * alpha) / (num_div - 1);
			return energies;
		}

		// Chebyshev-node grid cos(pi(2i+1/2)/M) (M points); rearranged to ascending when the
		// caller integrates (Bastin needs a monotone axis for the cumulative sum).
		std::vector<double> make_fft_node_grid(int num_div, bool integrate)
		{
			std::vector<double> energies(num_div);
			const double disp = 0.5;
			for (int i = 0; i < num_div; ++i)
				energies[i] = std::cos( M_PI * (2 * i + disp) / num_div );
			if (integrate) rearrange_crescent_order(energies);  // ascending for the integral
			return energies;
		}

		// Per-point output values aligned to energies[0..N-2]. Bastin: forward-cumulative
		// rectangle integral (Fermi sea). Greenwood: the kernel reported pointwise. This is
		// exactly what the legacy output loop did, factored out of the energy-axis mapping.
		std::vector<double> accumulate(const std::vector<double>& kernel,
		                               const std::vector<double>& energies, bool integrate)
		{
			const int num_div = (int)kernel.size();
			std::vector<double> values;
			values.reserve(num_div > 0 ? num_div - 1 : 0);
			double acc = 0.0;
			for (int i = 0; i < num_div - 1; ++i)
			{
				if (integrate) { acc += kernel[i] * (energies[i+1] - energies[i]); values.push_back(acc); }
				else           { values.push_back(kernel[i]); }
			}
			return values;
		}

		// Cumulative trapezoid integral (Sea Fermi-sea): values aligned to energies[0..N-2].
		std::vector<double> accumulate_trapezoid(const std::vector<double>& kernel,
		                                         const std::vector<double>& energies)
		{
			const int num_div = (int)kernel.size();
			std::vector<double> values;
			values.reserve(num_div > 0 ? num_div - 1 : 0);
			double acc = 0.0;
			for (int i = 0; i < num_div - 1; ++i)
			{
				const double denerg = energies[i+1] - energies[i];
				acc += 0.5 * (kernel[i] + kernel[i+1]) * denerg;
				values.push_back(acc);
			}
			return values;
		}

		// Pointwise midpoint average (Surf): value[i] = (k[i+1]+k[i])/2.
		std::vector<double> accumulate_midpoint(const std::vector<double>& kernel)
		{
			const int num_div = (int)kernel.size();
			std::vector<double> values;
			values.reserve(num_div > 0 ? num_div - 1 : 0);
			for (int i = 0; i < num_div - 1; ++i)
				values.push_back( (kernel[i+1] + kernel[i]) / 2.0 );
			return values;
		}

		// Zip each value with its physical energy (the per-grid axis transform is the caller's).
		template <class ToPhys>
		std::vector<std::pair<double,double> >
		write_pairs(const std::vector<double>& values, const std::vector<double>& energies, ToPhys to_phys)
		{
			std::vector<std::pair<double,double> > result;
			result.reserve(values.size());
			for (size_t i = 0; i < values.size(); ++i)
				result.push_back(std::make_pair(to_phys(energies[i]), values[i]));
			return result;
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
			energies = make_uniform_grid(num_div, xbound);
		}
		else // GRID_FFT_NODES
		{
			energies = make_fft_node_grid(mu.HighestMomentNumber(), obs.integrate);
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
		const std::vector<double> values = accumulate(kernel, energies, obs.integrate);
		if (grid == GRID_UNIFORM)
		{
			const double sf = mu.ScaleFactor(), shf = mu.ShiftFactor();
			return write_pairs(values, energies, [sf, shf](double e){ return e / sf - shf; });
		}
		const double hw = mu.HalfWidth(), bc = mu.BandCenter();
		return write_pairs(values, energies, [hw, bc](double e){ return e * hw + bc; });
	}

	std::vector<std::pair<double,double> >
	reconstruct_density_grid(const std::vector<double>& mu_real, double alpha, int num_div)
	{
		const std::vector<double> energies = make_uniform_grid(num_div, alpha);

		std::vector<double> values(num_div, 0.0);
		const int M = static_cast<int>(mu_real.size());
		#pragma omp parallel for
		for (int i = 0; i < num_div; ++i)
		{
			double s = 0.0;
			const double energ = energies[i];
			for (int m = 0; m < M; ++m)
				s += delta_chebF(energ, m) * mu_real[m];
			values[i] = s;
		}
		return write_pairs(values, energies, [](double e){ return e; });
	}

	// --- Kubo-Bastin Sea: Green-difference form, cumulative Fermi-sea integral ----------------
	std::vector<std::pair<double,double> >
	reconstruct_kubo_bastin_sea(chebyshev::Moments2D& mu)
	{
		const int    num_div = 30 * mu.HighestMomentNumber();
		const double xbound  = chebyshev::safety_factors().recon_cutoff;
		const std::vector<double> energies = make_uniform_grid(num_div, xbound);

		std::vector<double> kernel(num_div, 0.0);
		#pragma omp parallel for
		for (int i = 0; i < num_div; ++i)
		{
			const double energ = energies[i];
			for (int m0 = 0; m0 < mu.HighestMomentNumber(0); ++m0)
			for (int m1 = 0; m1 < mu.HighestMomentNumber(1); ++m1)
			{
				const auto Gr  = greenR_chebF(energ, m0);
				const auto Ga  = std::conj( greenR_chebF(energ, m0) );
				const auto DGr = DgreenR_chebF(energ, m1).real();
				const auto DGa = std::conj( DGr );
				kernel[i] += -( (Gr - Ga) * ( DGr + DGa ) * mu(m0, m1) ).real();
			}
			kernel[i] *= mu.SystemSize() * mu.ScaleFactor() * mu.ScaleFactor() / 2 / M_PI;
		}

		const std::vector<double> values = accumulate_trapezoid(kernel, energies);
		const double sf = mu.ScaleFactor(), shf = mu.ShiftFactor();
		return write_pairs(values, energies, [sf, shf](double e){ return e / sf - shf; });
	}

	// --- Kubo-Bastin Surf: double-delta form, pointwise (midpoint) ----------------------------
	std::vector<std::pair<double,double> >
	reconstruct_kubo_bastin_surf(chebyshev::Moments2D& mu)
	{
		const int    num_div = 30 * mu.HighestMomentNumber();
		const double xbound  = chebyshev::safety_factors().recon_cutoff;
		const std::vector<double> energies = make_uniform_grid(num_div, xbound);

		std::vector<double> kernel(num_div, 0.0);
		#pragma omp parallel for
		for (int i = 0; i < num_div; ++i)
		{
			const double energ = energies[i];
			for (int m0 = 0; m0 < mu.HighestMomentNumber(0); ++m0)
			for (int m1 = 0; m1 < mu.HighestMomentNumber(1); ++m1)
				kernel[i] += delta_chebF(energ, m0) * delta_chebF(energ, m1) * mu(m0, m1).real();
			kernel[i] *= M_PI * mu.SystemSize() * mu.ScaleFactor() * mu.ScaleFactor();
		}

		const std::vector<double> values = accumulate_midpoint(kernel);
		const double sf = mu.ScaleFactor(), shf = mu.ShiftFactor();
		return write_pairs(values, energies, [sf, shf](double e){ return e / sf - shf; });
	}
}
