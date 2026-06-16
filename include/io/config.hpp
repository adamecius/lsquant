#ifndef LSQUANT_IO_CONFIG_HPP
#define LSQUANT_IO_CONFIG_HPP

// Unified run configuration, schema v2 (unified_io.md, Phase I1/I2).
//
// One config defines a SYSTEM and a list of OBSERVABLES computed on it, so a single
// `lsquant run --config run.json` produces one results.json holding DOS +
// conductivity + spin relaxation together. The flat v1 run.json (single mode) is
// still accepted by `lsquant compute`; v2 adds the "system"/"observables" shape.
//
//   { "system":   { "label": "...", "operators_dir": "operators", "bounds": "auto" },
//     "numerics": { "num_moments": 512, "kernel": "jackson", "broadening_meV": 10,
//                   "grid_points": 600, "state": "default" },
//     "observables": [ {"kind":"dos"},
//                      {"kind":"conductivity","component":"xx"},
//                      {"kind":"conductivity","component":"xy"},
//                      {"kind":"spin_relaxation","axis":"z","num_times":64,"tmax":500} ],
//     "output": { "results": "graphene.results.json" } }

#include <string>
#include <vector>
#include "io/json.hpp"

namespace lsquant { namespace io {

struct Observable {
	std::string kind;        // "dos" | "conductivity" | "spin_relaxation"
	std::string component;   // conductivity: "xx" | "xy"
	std::string op;          // operator-id override (dos defaults to identity "1")
	std::string axis;        // spin_relaxation: "x"|"y"|"z"
	int         num_times = 0;
	double      tmax      = 0.0;

	// Result-block name, e.g. "dos", "sigma_xx", "spin_relax_z".
	std::string result_name() const;
};

struct RunSpec {
	std::string label;
	std::string operators_dir  = "operators";
	std::string bounds         = "auto";
	int         num_moments    = 0;
	std::string kernel         = "jackson";
	double      broadening_meV = 0.0;
	int         grid_points    = 600;
	std::string state          = "default";
	std::string results_path;            // output path (defaults to <label>.results.json)
	std::vector<Observable> observables;
};

struct Hint {
	std::string severity;   // "error" | "warning"
	std::string where;      // e.g. "observables[1]"
	std::string message;    // human, physics-aware
	std::string suggestion; // optional fix
};

// Structural parse of a v2 document into a RunSpec. Throws lsquant::ConfigError on
// a fundamentally malformed document; per-field problems are reported by validate().
RunSpec parse_run_spec(const json::Value& doc);

// Returns true when the spec is runnable (no "error" hints). Fills `hints` with
// physics-aware errors/warnings; errors block the run, warnings are advisory.
bool validate(const RunSpec& spec, std::vector<Hint>& hints);

}} // namespace lsquant::io

#endif // LSQUANT_IO_CONFIG_HPP
