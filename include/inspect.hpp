#ifndef LSQUANT_INSPECT
#define LSQUANT_INSPECT

#include <string>

namespace lsquant
{
	// Print what a run.json (or an operator .desc) will do, including the resolved numerical
	// provenance (recon alpha, bounds source, hbar/units, kernel). Dispatches on the file suffix
	// (.desc vs .json). Returns 0 on success. Shared by the standalone `lsquant_inspect` binary
	// and the `lsquant inspect` subcommand (Phase 7).
	int inspect_path(const std::string& path);
}

#endif // LSQUANT_INSPECT
