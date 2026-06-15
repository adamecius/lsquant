#ifndef LSQUANT_COMPUTE
#define LSQUANT_COMPUTE

#include <string>

namespace lsquant
{
	// Compute the 2D Kubo non-equilibrium Chebyshev moments mu_{mn} = <phi| OPL T_m(H~) OPR
	// T_n(H~) |phi> from operators/<label>.{HAM,OPR,OPL}.CSR, and save the .chebmom2D. Spectral
	// bounds come from operators/<label>.HAM.desc (Phase 2), else a BOUNDS file / Gershgorin.
	// state_file empty -> default random-phase state; otherwise LoadStateFile(state_file).
	// Returns 0 on success. Shared by the standalone inline_compute-kpm-nonEqOp driver and the
	// `lsquant compute` subcommand (Phase 7).
	int compute_noneq(const std::string& label, const std::string& op_right,
	                  const std::string& op_left, int num_moments, const std::string& state_file);
}

#endif // LSQUANT_COMPUTE
