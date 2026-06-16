#ifndef LSQUANT_RUN_CONFIG
#define LSQUANT_RUN_CONFIG

#include <string>
#include <map>

// A run is described by a small, human-readable JSON object instead of positional argv and
// filename-encoded metadata. Example run.json (flat object, the v1 schema):
//
//   {
//     "label":          "graphene_golden",
//     "operator_left":  "VX",
//     "operator_right": "VX",
//     "num_moments":    64,
//     "kernel":         "jackson",
//     "lambda":         0.0,
//     "state":          "default",
//     "alpha":          0.95,
//     "output":         "NonEqOpVX-VXgraphene_goldenKPM_M64x64_statedefault.chebmom2D"
//   }
//
// The parser is a tiny, dependency-FREE flat-JSON reader (string/number/bool values, no nested
// objects or arrays). This keeps lsquant offline-safe -- no fetched JSON library -- consistent
// with the w2s decoupling; a full JSON lib can replace parse_flat_json later without touching
// callers. Numerical knobs (alpha, bounds) still resolve through their existing homes
// (recon_grid.hpp / .desc); run.json just records the intent, and `lsquant inspect` reports it.

namespace lsquant
{
	// Parse a flat JSON object into raw (unquoted) string values keyed by name.
	// Returns empty map if the file is missing or has no '{'.
	std::map<std::string, std::string> parse_flat_json(const std::string& path);

	struct RunConfig
	{
		std::string label;
		std::string mode = "noneq";    // noneq | spectral | msd  (which compute_* to dispatch)
		std::string operator_left;     // OPL (noneq)
		std::string operator_right;    // OPR (noneq)
		std::string op;                // single operator: OP (spectral) / velocity (msd)
		int         num_moments = 0;
		int         num_times   = 0;   // msd time grid size
		double      tmax        = 0.0; // msd max time
		std::string kernel = "jackson";
		double      lambda  = 0.0;
		std::string state   = "default";
		double      alpha   = -1.0;    // <0 means "unset -> use the recon_grid default/env"
		std::string output;            // declared output name (optional)
		bool        valid   = false;
		std::string error;
	};

	RunConfig read_run_config(const std::string& path);
}

#endif // LSQUANT_RUN_CONFIG
