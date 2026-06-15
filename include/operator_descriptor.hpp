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

#include <vector>
#include <complex>

namespace lsquant
{
	constexpr int DESCRIPTOR_SCHEMA_VERSION = 2;   // v2: optional basis block (Phase 3)

	struct OperatorDescriptor
	{
		std::string observable;
		std::string component;
		std::string units;
		std::string provenance;
		bool        has_bounds = false;   // (a,b) meaningful (Hamiltonian)
		double      a = 0.0, b = 0.0;     // spectral_min, spectral_max
		bool        valid = false;        // false if the file was absent / unparseable

		// --- basis block (Phase 3, manual tight-binding traceability) -------------------
		// Minimal orbital-tag basis for a spin-doubled TB model: each unit cell holds
		// `orbitals_per_cell` orbitals whose spins are `orbital_spin` (+1 = up, -1 = down);
		// the global index of (cell c, orbital o) is c*orbitals_per_cell + o. This makes the
		// spin/orbital operators (SX/SY/SZ) DERIVABLE from the descriptor rather than opaque
		// CSR blobs (the build_spin_operator gate below reproduces them). Kept IN the .desc
		// (the Phase-3 "base in .desc" choice); a large H may move it to a shared sidecar later.
		bool             has_basis = false;
		int              orbitals_per_cell = 0;
		std::vector<int> orbital_spin;    // size = orbitals_per_cell; +1/-1 per orbital
		int              num_cells = 0;
	};

	// Parse a .desc file. Returns a descriptor with valid=false if the path cannot be opened.
	OperatorDescriptor read_descriptor(const std::string& path);

	// Convenience: if `desc_path` exists and carries bounds, fill (a,b) and return true.
	bool descriptor_bounds(const std::string& desc_path, double& a, double& b);

	// Build a Pauli spin operator (component 'X'|'Y'|'Z') over the descriptor's spin-doubled
	// basis, as COO triplets. Within each cell the up/down orbital pair carries the standard
	// Pauli matrix: sigma_x=[[0,1],[1,0]], sigma_y=[[0,-i],[i,0]], sigma_z=[[1,0],[0,-1]].
	// Returns false if the descriptor has no usable spin-doubled basis.
	bool build_spin_operator(const OperatorDescriptor& d, char component,
	                         std::vector<int>& rows, std::vector<int>& cols,
	                         std::vector<std::complex<double> >& vals);
}

#endif // LSQUANT_OPERATOR_DESCRIPTOR
