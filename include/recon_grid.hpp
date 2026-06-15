#ifndef LSQUANT_RECON_GRID
#define LSQUANT_RECON_GRID

#include <cstdlib>   // getenv, atof

namespace chebyshev
{
	// Two DISTINCT, named safety factors. Deliberately NOT a single "alpha": they do different
	// jobs, live at different points, and have different correct values.
	struct SafetyFactors
	{
		// recon_cutoff (alpha): the reconstruction grid runs on [-alpha, alpha] (or, for the
		// FFT route, with an edge cushion 1 - alpha). It keeps the grid off the 1/sqrt(1-x^2)
		// reconstruction-kernel singularity at x = +-1 -- the band-edge NaN that poisoned the
		// cumulative Kubo-Bastin integral. Default 0.95.
		//
		//   *** NOT a moment knob. *** The moment rescale stays at CUTOFF = 1.00 (the band must
		//   fill [-1,1] EXACTLY for the analytic moment identities -- DOS moments delta_{m,0},
		//   Kubo-Greenwood matrix B.15). Under-filling to [-alpha,alpha] in moments moves every
		//   moment and breaks the oracle (residual 0.05-0.375 vs tol 1e-9). alpha = 1.0 selects
		//   the legacy moment-identity reconstruction grid.
		double recon_cutoff;

		// bounds_pad: fractional margin added to the Gershgorin spectral-bounds estimate so the
		// rescaled spectrum stays inside [-1,1] (Chebyshev recursion stability). This is the
		// estimator's safety margin, independent of recon_cutoff. Default 0.10.
		double bounds_pad;
	};

	// Process-wide safety factors, resolved ONCE. Defaults reproduce the shipped goldens
	// bit-for-bit. Override at runtime via environment -- the Phase-1 "alpha as a first-class
	// runtime parameter" step; Phase 4 will route these through run.json:
	//   LSQUANT_RECON_ALPHA  -> recon_cutoff   (e.g. 1.0 for the legacy/edge-unsafe grid)
	//   LSQUANT_BOUNDS_PAD   -> bounds_pad
	inline const SafetyFactors& safety_factors()
	{
		static const SafetyFactors sf = []
		{
			SafetyFactors s{ 0.95, 0.10 };
			if (const char* a = std::getenv("LSQUANT_RECON_ALPHA")) s.recon_cutoff = std::atof(a);
			if (const char* p = std::getenv("LSQUANT_BOUNDS_PAD"))  s.bounds_pad   = std::atof(p);
			return s;
		}();
		return sf;
	}
}

#endif // LSQUANT_RECON_GRID
