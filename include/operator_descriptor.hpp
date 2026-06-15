#ifndef LSQUANT_OPERATOR_DESCRIPTOR
#define LSQUANT_OPERATOR_DESCRIPTOR

#include <string>

// Reader for the operator ".desc" sidecar emitted by wannier2sparse (--bounds).
//
// SCHEMA COUPLING DECISION (refactor plan Phase 2): this is a DUPLICATED-BUT-VERSIONED mirror
// of the producer schema in wannier2sparse/include/descriptor.hpp -- NOT a shared submodule.
// Rationale: a format submodule would re-couple lsquant's build to w2s's layout and undo the
// offline-safe decoupling already banked (FetchWannier2Sparse). The risk of duplication --
// drift between writer and reader -- is closed by the `descriptor_conformance` test, which
// parses a committed w2s-produced .desc and checks every field round-trips.
//
// Schema v1 (plain "key: value" text, one per line):
//   observable:   hamiltonian | velocity | spin | spin_current | external | ...
//   component:    "" | X | Y | Z | XSZ | ...
//   units:        eV | eV*Angstrom | hbar/2 | dimensionless | ...
//   provenance:   free text (how it was produced)
//   spectral_min: <double>   } present ONLY when the operator carries bounds
//   spectral_max: <double>   }  (typically just the Hamiltonian; Lanczos estimate)
//
// The descriptor is the AUTHORITATIVE source of an operator's identity and of the spectral
// window (a,b) used to rescale H into [-1,1]. It never enters the Chebyshev recursion.

namespace lsquant
{
	constexpr int DESCRIPTOR_SCHEMA_VERSION = 1;

	struct OperatorDescriptor
	{
		std::string observable;
		std::string component;
		std::string units;
		std::string provenance;
		bool        has_bounds = false;   // (a,b) meaningful (Hamiltonian)
		double      a = 0.0, b = 0.0;     // spectral_min, spectral_max
		bool        valid = false;        // false if the file was absent / unparseable
	};

	// Parse a .desc file. Returns a descriptor with valid=false if the path cannot be opened.
	OperatorDescriptor read_descriptor(const std::string& path);

	// Convenience: if `desc_path` exists and carries bounds, fill (a,b) and return true.
	bool descriptor_bounds(const std::string& desc_path, double& a, double& b);
}

#endif // LSQUANT_OPERATOR_DESCRIPTOR
