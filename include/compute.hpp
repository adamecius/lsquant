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
	// out_filename (optional): receives the written .chebmom path, so an orchestrator
	// (lsquant run) can reconstruct it in-memory without re-deriving the name.
	int compute_noneq(const std::string& label, const std::string& op_right,
	                  const std::string& op_left, int num_moments, const std::string& state_file,
	                  std::string* out_filename = 0);

	// Compute the 1D spectral (DOS / spectral-function) Chebyshev moments mu_m = <phi| OP T_m(H~)
	// |phi> from operators/<label>.{HAM,OP}.CSR, and save the .chebmom1D. Same orchestration as the
	// standalone inline_compute-kpm-spectralOp driver and the `lsquant compute --config` (mode
	// "spectral"). state_file empty -> default random-phase state. Returns 0 on success.
	int compute_spectral(const std::string& label, const std::string& op,
	                     int num_moments, const std::string& state_file,
	                     std::string* out_filename = 0);

	// Compute the time-dependent mean-square-displacement Chebyshev moments from
	// operators/<label>.{HAM,VOP}.CSR and save the .chebmomTD. num_times/tmax set the time grid
	// (TimeDiff = tmax/(num_times-1)). Shared by the standalone
	// inline_compute-kpm-MeanSquareDisplacement driver and `lsquant compute --config` (mode "msd").
	// state_file empty -> default random-phase state. Returns 0 on success.
	int compute_msd(const std::string& label, const std::string& vop,
	                int num_moments, int num_times, double tmax, const std::string& state_file,
	                std::string* out_filename = 0);
}

#endif // LSQUANT_COMPUTE
