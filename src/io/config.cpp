#include "io/config.hpp"
#include "util/error.hpp"

#include <fstream>
#include <cctype>

namespace lsquant { namespace io {

namespace {
bool file_exists(const std::string& p) { std::ifstream f(p.c_str()); return f.good(); }
std::string lower(const std::string& s) {
	std::string o; o.reserve(s.size());
	for (char c : s) o += static_cast<char>(std::tolower((unsigned char)c));
	return o;
}
} // namespace

std::string Observable::result_name() const
{
	if (kind == "conductivity") return component.empty() ? "sigma" : ("sigma_" + component);
	if (kind == "spin_relaxation") return axis.empty() ? "spin_relax" : ("spin_relax_" + axis);
	return kind;   // "dos"
}

RunSpec parse_run_spec(const json::Value& doc)
{
	if (!doc.is_object())
		throw lsquant::ConfigError("run config: top level must be a JSON object");

	RunSpec s;
	const json::Value& sys = doc.at("system");
	s.label         = sys.at("label").as_string();
	s.operators_dir = sys.contains("operators_dir") ? sys.at("operators_dir").as_string() : "operators";
	s.bounds        = sys.contains("bounds") ? sys.at("bounds").as_string("auto") : "auto";

	const json::Value& num = doc.at("numerics");
	s.num_moments    = num.at("num_moments").as_int(0);
	s.kernel         = num.contains("kernel") ? num.at("kernel").as_string("jackson") : "jackson";
	s.broadening_meV = num.at("broadening_meV").as_number(0.0);
	s.grid_points    = num.contains("grid_points") ? num.at("grid_points").as_int(600) : 600;
	s.state          = num.contains("state") ? num.at("state").as_string("default") : "default";

	if (doc.contains("output") && doc.at("output").contains("results"))
		s.results_path = doc.at("output").at("results").as_string();
	if (s.results_path.empty() && !s.label.empty())
		s.results_path = s.label + ".results.json";

	const json::Value& obs = doc.at("observables");
	for (size_t k = 0; k < obs.size(); ++k) {
		const json::Value& o = obs[k];
		Observable ob;
		ob.kind      = o.at("kind").as_string();
		ob.component = o.contains("component") ? o.at("component").as_string() : "";
		ob.op        = o.contains("operator")  ? o.at("operator").as_string()  : "";
		ob.axis      = o.contains("axis")       ? o.at("axis").as_string()      : "";
		ob.num_times = o.contains("num_times") ? o.at("num_times").as_int(0) : 0;
		ob.tmax      = o.contains("tmax")       ? o.at("tmax").as_number(0.0) : 0.0;
		s.observables.push_back(ob);
	}
	return s;
}

bool validate(const RunSpec& s, std::vector<Hint>& hints)
{
	const std::string& dir = s.operators_dir;
	const std::string base = dir + "/" + s.label + ".";

	if (s.label.empty())
		hints.push_back(Hint{"error", "system.label", "missing 'label': it names the system and its operator files.", "set system.label to your <label> (operators/<label>.HAM.CSR must exist)"});
	else if (!file_exists(base + "HAM.CSR"))
		hints.push_back(Hint{"error", "system.label", "Hamiltonian not found: " + base + "HAM.CSR is missing.", "generate the operators or fix system.operators_dir / system.label"});

	if (s.num_moments <= 0)
		hints.push_back(Hint{"error", "numerics.num_moments", "num_moments not set: it controls energy resolution (~ bandwidth / num_moments).", "typical 256-2048"});

	if (s.observables.empty())
		hints.push_back(Hint{"error", "observables", "no observables requested: nothing to compute.", "add e.g. {\"kind\":\"dos\"} or {\"kind\":\"conductivity\",\"component\":\"xx\"}"});

	for (size_t k = 0; k < s.observables.size(); ++k) {
		const Observable& o = s.observables[k];
		std::string w = "observables[" + std::to_string(k) + "]";
		if (o.kind == "dos") {
			// identity operator -> DOS; an explicit op is allowed (spectral function).
			if (!o.op.empty() && o.op != "1" && !file_exists(base + o.op + ".CSR"))
				hints.push_back(Hint{"error", w, "spectral operator '" + o.op + "' not found: " + base + o.op + ".CSR missing.", "use \"operator\":\"1\" for DOS or provide the .CSR"});
		}
		else if (o.kind == "conductivity") {
			const std::string comp = lower(o.component);
			if (comp != "xx" && comp != "xy")
				hints.push_back(Hint{"error", w, "conductivity needs component 'xx' or 'xy' (got '" + o.component + "').", "set component to xx (longitudinal) or xy (Hall)"});
			else {
				if (!file_exists(base + "VX.CSR"))
					hints.push_back(Hint{"error", w, "conductivity needs the velocity operator VX: " + base + "VX.CSR missing.", "provide " + s.label + ".VX.CSR"});
				if (comp == "xy" && !file_exists(base + "VY.CSR"))
					hints.push_back(Hint{"error", w, "Hall conductivity sigma_xy needs TWO velocity operators; VY is missing (" + base + "VY.CSR).", "provide " + s.label + ".VY.CSR or use component xx"});
			}
		}
		else if (o.kind == "spin_relaxation") {
			if (o.num_times <= 0)
				hints.push_back(Hint{"error", w, "spin_relaxation needs num_times > 0 (time-grid size).", "e.g. num_times: 64"});
			if (o.tmax <= 0.0)
				hints.push_back(Hint{"error", w, "spin_relaxation needs tmax > 0 (max time, fs).", "e.g. tmax: 500"});
			const std::string ax = o.axis.empty() ? "z" : lower(o.axis);
			const std::string sop = "S" + std::string(1, static_cast<char>(std::toupper(ax[0])));
			if (!file_exists(base + sop + ".CSR"))
				hints.push_back(Hint{"warning", w, "spin operator " + sop + " not found (" + base + sop + ".CSR); spin relaxation needs a spin-resolved basis.", "provide " + s.label + "." + sop + ".CSR"});
		}
		else {
			hints.push_back(Hint{"warning", w, "unknown observable kind '" + o.kind + "' (known: dos, conductivity, spin_relaxation); it will be skipped.", ""});
		}
	}

	for (size_t k = 0; k < hints.size(); ++k)
		if (hints[k].severity == "error") return false;
	return true;
}

}} // namespace lsquant::io
