// lsquant -- single-binary dispatcher (Phase 7).
//
// Strangler-fig: this binary grows BESIDE the legacy per-task executables (41 main()s). It
// already subsumes the parts that the refactor has library-ized -- `inspect` (run.json / .desc)
// and `reconstruct` (the unified Kubo reconstruction) -- and will absorb `compute` and the
// model plugins as Phase 6/7 proceed. Subcommand style: `lsquant <command> [args...]`.
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <complex>

#include "inspect.hpp"
#include "observable.hpp"
#include "compute.hpp"
#include "run_config.hpp"
#include "chebyshev_moments.hpp"
#include "util/report.hpp"
#include "io/json.hpp"
#include "io/config.hpp"
#include "io/fingerprint.hpp"
#include "io/results.hpp"

static int usage(const char* prog)
{
	std::cerr << "usage: " << prog << " <command> [args...]\n\n"
	          << "commands:\n"
	          << "  compute --config <run.json>\n"
	          << "        compute Chebyshev moments for a run (mode: noneq|spectral|msd);\n"
	          << "        writes NonEqOp*.chebmom2D / SpectralOp*.chebmom1D / Correlation*.chebmomTD\n"
	          << "  reconstruct <moments> <dos|bastin|greenwood> <broadening_meV>\n"
	          << "        reconstruct a spectrum via the unified routine: dos reads a .chebmom1D\n"
	          << "        and writes mean*JACKSON.dat; bastin/greenwood read a .chebmom2D and write\n"
	          << "        Kubo*JACKSON.dat\n"
	          << "  inspect <run.json | operator.desc>\n"
	          << "        print a run / operator and its resolved numerical provenance\n"
	          << "  run --config <run.json>\n"
	          << "        unified entry point: one system, many observables (dos / conductivity /\n"
	          << "        spin_relaxation) -> one <label>.results.json (machine- + human-readable)\n"
	          << "  merge <a.results.json> <b.results.json> [...] -o <out.json>\n"
	          << "        merge results for the SAME system (fingerprint-checked)\n"
	          << "  fingerprint <label> [operators_dir]\n"
	          << "        print the system fingerprint (SHA-256 of its operator files)\n"
	          << "\nglobal options (any position):\n"
	          << "  -v/--verbose, -vv, -q/--quiet, --log-level <trace|debug|info|warn|error>,\n"
	          << "  --log-file <path>   (also via env LSQUANT_LOG_LEVEL / LSQUANT_LOG_FILE)\n";
	return 2;
}

// Strip the global reporting flags (already consumed by init_reporting) so the
// per-command argument checks see only their own arguments.
static std::vector<char*> strip_global_flags(int argc, char** argv)
{
	std::vector<char*> out;
	for (int i = 0; i < argc; ++i)
	{
		const std::string a = argv[i];
		if (a == "-v" || a == "--verbose" || a == "-vv" || a == "-q" || a == "--quiet")
			continue;
		if ((a == "--log-level" || a == "--log-file") && i + 1 < argc) { ++i; continue; }
		out.push_back(argv[i]);
	}
	return out;
}

static int cmd_compute(int argc, char** argv)
{
	if (argc != 4 || std::string(argv[2]) != "--config")
	{ std::cerr << "usage: " << argv[0] << " compute --config <run.json>\n"; return 2; }
	const lsquant::RunConfig c = lsquant::read_run_config(argv[3]);
	if (!c.valid) { std::cerr << "config error: " << c.error << std::endl; return 1; }
	const std::string state = (c.state != "default" && !c.state.empty()) ? c.state : std::string();
	if (c.mode == "noneq")
		return lsquant::compute_noneq(c.label, c.operator_right, c.operator_left, c.num_moments, state);
	if (c.mode == "spectral")
		return lsquant::compute_spectral(c.label, c.op, c.num_moments, state);
	if (c.mode == "msd")
		return lsquant::compute_msd(c.label, c.op, c.num_moments, c.num_times, c.tmax, state);
	std::cerr << "compute: unknown mode '" << c.mode << "' (noneq|spectral|msd)\n";
	return 2;
}

static int cmd_reconstruct(int argc, char** argv)
{
	if (argc != 5) { std::cerr << "usage: " << argv[0]
		<< " reconstruct <moments> <dos|bastin|greenwood> <broadening_meV>\n"; return 2; }
	const std::string momfile = argv[2];
	const std::string which   = argv[3];
	const double broadening   = std::stod(argv[4]);

	// DOS / spectral function: the 1-D density-grid route (mirrors the
	// inline_spectralFunctionFromChebmom driver, byte-identical: same kernel, grid, prefactor
	// N*ScaleFactor, and energy axis x/ScaleFactor + ShiftFactor).
	if (which == "dos")
	{
		chebyshev::Moments1D mu(momfile.c_str());
		mu.ApplyJacksonKernel(broadening);
		const int num_div = 30 * mu.HighestMomentNumber();
		std::vector<double> mu_real(mu.HighestMomentNumber());
		for (int m = 0; m < mu.HighestMomentNumber(); ++m) mu_real[m] = mu(m).real();
		const std::vector<std::pair<double,double> > grid =
			lsquant::reconstruct_density_grid(mu_real, chebyshev::safety_factors().recon_cutoff, num_div);

		const std::string outputName = "mean" + mu.SystemLabel() + "JACKSON.dat";
		std::cout << "Saving the data in " << outputName << std::endl;
		std::ofstream f(outputName.c_str());
		for (size_t i = 0; i < grid.size(); ++i)
			f << grid[i].first/mu.ScaleFactor() + mu.ShiftFactor()
			  << " " << grid[i].second*( mu.SystemSize()*mu.ScaleFactor() ) << std::endl;
		return 0;
	}

	const lsquant::KuboObservable* obs = 0;
	std::string tag;
	if      (which == "bastin")    { obs = &lsquant::KUBO_BASTIN;    tag = "KuboBastin_"; }
	else if (which == "greenwood") { obs = &lsquant::KUBO_GREENWOOD; tag = "KuboGreenWood_"; }
	else { std::cerr << "reconstruct: formula must be 'dos', 'bastin' or 'greenwood'\n"; return 2; }

	chebyshev::Moments2D mu(momfile.c_str());
	mu.ApplyJacksonKernel(broadening, broadening);
	const std::vector<std::pair<double,double> > data = lsquant::reconstruct_kubo(mu, *obs);

	const std::string outputName = tag + mu.SystemLabel() + "JACKSON.dat";
	std::cout << "Saving the data in " << outputName << std::endl;
	std::ofstream f(outputName.c_str());
	for (size_t i = 0; i < data.size(); ++i) f << data[i].first << " " << data[i].second << std::endl;
	return 0;
}

// --- unified I/O helpers (unified_io.md) -------------------------------------
namespace {

// Downsample a reconstructed (x,y) curve to ~grid_points and pack it as a JSON
// block {xname:[...], yname:[...]}.
lsquant::json::Value thin_curve(const std::vector<std::pair<double,double> >& g,
                                const std::string& xname, const std::string& yname,
                                int grid_points)
{
	lsquant::json::Value x = lsquant::json::Value::array();
	lsquant::json::Value y = lsquant::json::Value::array();
	const size_t n = g.size();
	const size_t target = grid_points > 0 ? static_cast<size_t>(grid_points) : 600;
	const size_t stride = (n > target && target > 0) ? (n / target) : 1;
	for (size_t i = 0; i < n; i += stride) {
		x.push_back(lsquant::json::Value(g[i].first));
		y.push_back(lsquant::json::Value(g[i].second));
	}
	lsquant::json::Value blk = lsquant::json::Value::object();
	blk.set(xname, x);
	blk.set(yname, y);
	return blk;
}

} // namespace

static int cmd_fingerprint(int argc, char** argv)
{
	if (argc < 3) { std::cerr << "usage: " << argv[0]
		<< " fingerprint <label> [operators_dir]\n"; return 2; }
	const std::string label = argv[2];
	const std::string dir   = (argc >= 4) ? argv[3] : "operators";
	std::cout << lsquant::io::fingerprint_system(label, dir) << "\n";
	return 0;
}

static int cmd_merge(int argc, char** argv)
{
	std::vector<std::string> inputs;
	std::string out;
	for (int i = 2; i < argc; ++i) {
		const std::string a = argv[i];
		if (a == "-o" && i + 1 < argc) out = argv[++i];
		else inputs.push_back(a);
	}
	if (inputs.size() < 2 || out.empty()) {
		std::cerr << "usage: " << argv[0]
		          << " merge <a.results.json> <b.results.json> [...] -o <out.results.json>\n"
		          << "  merges results for the SAME system (fingerprint-checked).\n";
		return 2;
	}
	lsquant::json::Value base = lsquant::io::load_results(inputs[0]);
	for (size_t k = 1; k < inputs.size(); ++k)
		lsquant::io::merge_into(base, lsquant::io::load_results(inputs[k]));   // throws on mismatch
	lsquant::io::save_results(base, out);
	LSQ_LOG_INFO("merged " << inputs.size() << " results into " << out);
	return 0;
}

// `lsquant run --config <v2.json>` -- one system, many observables, one results.json.
static int cmd_run(int argc, char** argv)
{
	std::string cfg;
	for (int i = 2; i < argc; ++i) {
		const std::string a = argv[i];
		if (a == "--config" && i + 1 < argc) cfg = argv[++i];
	}
	if (cfg.empty()) { std::cerr << "usage: " << argv[0] << " run --config <run.json>\n"; return 2; }

	const lsquant::json::Value doc = lsquant::json::parse_file(cfg);
	lsquant::io::RunSpec spec = lsquant::io::parse_run_spec(doc);

	std::vector<lsquant::io::Hint> hints;
	const bool ok = lsquant::io::validate(spec, hints);
	for (size_t k = 0; k < hints.size(); ++k) {
		const lsquant::io::Hint& h = hints[k];
		const std::string tail = h.suggestion.empty() ? "" : ("  (hint: " + h.suggestion + ")");
		if (h.severity == "error") LSQ_LOG_ERROR("[" << h.where << "] " << h.message << tail);
		else                       LSQ_LOG_WARN ("[" << h.where << "] " << h.message << tail);
	}
	if (!ok) { LSQ_LOG_ERROR("run aborted: the config is missing what a simulation needs (see above)."); return 1; }

	const std::string state = (spec.state != "default" && !spec.state.empty()) ? spec.state : std::string();
	const std::string fp = lsquant::io::fingerprint_system(spec.label, spec.operators_dir);
	LSQ_LOG_INFO("system '" << spec.label << "'  fingerprint " << fp.substr(0, 16) << "...  ("
	             << spec.observables.size() << " observable(s))");

	lsquant::json::Value results;
	bool header_set = false;

	for (size_t k = 0; k < spec.observables.size(); ++k) {
		const lsquant::io::Observable& o = spec.observables[k];
		const std::string rname = o.result_name();
		lsquant::util::Timer t;
		lsquant::json::Value block = lsquant::json::Value::object();

		try {
			if (o.kind == "dos") {
				const std::string op = o.op.empty() ? std::string("1") : o.op;
				std::string mom;
				lsquant::compute_spectral(spec.label, op, spec.num_moments, state, &mom);
				chebyshev::Moments1D mu(mom.c_str());
				mu.ApplyJacksonKernel(spec.broadening_meV);
				std::vector<double> mr(mu.HighestMomentNumber());
				for (int m = 0; m < mu.HighestMomentNumber(); ++m) mr[m] = mu(m).real();
				const int num_div = 30 * mu.HighestMomentNumber();
				std::vector<std::pair<double,double> > grid =
					lsquant::reconstruct_density_grid(mr, chebyshev::safety_factors().recon_cutoff, num_div);
				std::vector<std::pair<double,double> > curve; curve.reserve(grid.size());
				for (size_t i = 0; i < grid.size(); ++i)
					curve.push_back(std::make_pair(grid[i].first/mu.ScaleFactor() + mu.ShiftFactor(),
					                               grid[i].second*(mu.SystemSize()*mu.ScaleFactor())));
				block = thin_curve(curve, "energy_eV", "dos", spec.grid_points);
				if (!header_set) {
					results = lsquant::io::new_results_doc(spec.label, fp, mu.SystemSize(),
					                                       mu.BandWidth(), mu.BandCenter());
					header_set = true;
				}
			}
			else if (o.kind == "conductivity") {
				const bool hall = (o.component == "xy");
				const std::string opl = hall ? "VY" : "VX";
				std::string mom;
				lsquant::compute_noneq(spec.label, "VX", opl, spec.num_moments, state, &mom);
				chebyshev::Moments2D mu(mom.c_str());
				mu.ApplyJacksonKernel(spec.broadening_meV, spec.broadening_meV);
				const lsquant::KuboObservable& obs = hall ? lsquant::KUBO_BASTIN : lsquant::KUBO_GREENWOOD;
				std::vector<std::pair<double,double> > curve = lsquant::reconstruct_kubo(mu, obs);
				block = thin_curve(curve, "energy_eV", "sigma", spec.grid_points);
				block.set("units", lsquant::json::Value("e^2/h"));
				block.set("formula", lsquant::json::Value(hall ? "kubo-bastin" : "kubo-greenwood"));
				if (!header_set) {
					results = lsquant::io::new_results_doc(spec.label, fp, mu.SystemSize(),
					                                       mu.BandWidth(), mu.BandCenter());
					header_set = true;
				}
			}
			else if (o.kind == "spin_relaxation") {
				const std::string ax = o.axis.empty() ? "z" : o.axis;
				std::string sop = "S"; sop += static_cast<char>(std::toupper(ax[0]));
				std::string mom;
				lsquant::compute_msd(spec.label, sop, spec.num_moments, o.num_times, o.tmax, state, &mom);
				// Time-domain reconstruction (S(t), tau) is the timeCorrelations route; for now the
				// unified doc records the moments file + parameters so the curve can be reconstructed.
				block.set("moments_file", lsquant::json::Value(mom));
				block.set("spin_axis",    lsquant::json::Value(ax));
				block.set("num_times",    lsquant::json::Value(o.num_times));
				block.set("tmax_fs",      lsquant::json::Value(o.tmax));
				block.set("note", lsquant::json::Value("time-domain reconstruction pending; reconstruct from moments_file"));
				if (!header_set) {
					// minimal header without spectral scalars
					results = lsquant::io::new_results_doc(spec.label, fp, 0, 0.0, 0.0);
					header_set = true;
				}
			}
			else {
				LSQ_LOG_WARN("skipping unknown observable '" << o.kind << "'");
				continue;
			}
		}
		catch (const std::exception& e) {
			LSQ_LOG_ERROR("observable '" << rname << "' failed: " << e.what());
			lsquant::json::Value err = lsquant::json::Value::object();
			err.set("observable", lsquant::json::Value(rname));
			err.set("error", lsquant::json::Value(e.what()));
			if (header_set) {
				lsquant::json::Value run = results.at("run");
				lsquant::json::Value errs = run.at("errors");
				errs.push_back(err);
				run.set("errors", errs);
				results.set("run", run);
			}
			continue;
		}

		block.set("num_moments", lsquant::json::Value(spec.num_moments));
		block.set("kernel",      lsquant::json::Value(spec.kernel));
		block.set("timing_s",    lsquant::json::Value(t.elapsed_s()));
		block.set("warnings",    lsquant::json::Value::array());
		lsquant::io::add_result(results, rname, block);
		LSQ_LOG_INFO("  + " << rname << "  (" << t.elapsed_s() << " s)");
	}

	if (!header_set) { LSQ_LOG_ERROR("run produced no results."); return 1; }

	// run-level resources
	const lsquant::diag::ResourceSnapshot rs = lsquant::diag::sample();
	lsquant::json::Value run = results.at("run");
	run.set("wall_s",         lsquant::json::Value(rs.wall_sec));
	run.set("peak_rss_bytes", lsquant::json::Value(static_cast<double>(rs.peak_rss_bytes)));
	results.set("run", run);

	lsquant::io::save_results(results, spec.results_path);
	LSQ_LOG_INFO("wrote unified results to " << spec.results_path);
	lsquant::report_run_summary();
	return 0;
}

int main(int argc, char** argv)
{
	lsquant::init_reporting(argc, argv);
	std::vector<char*> clean = strip_global_flags(argc, argv);
	int    cargc = static_cast<int>(clean.size());
	char** cargv = clean.data();

	if (cargc < 2) return usage(cargv[0]);
	const std::string cmd = cargv[1];

	// One top-level catch: every library failure surfaces as a logged error with a
	// stable non-zero exit code, instead of an unhandled-exception abort.
	try
	{
		if (cmd == "inspect")
		{
			if (cargc < 3) { std::cerr << "usage: " << cargv[0] << " inspect <run.json | operator.desc>\n"; return 2; }
			return lsquant::inspect_path(cargv[2]);
		}
		if (cmd == "compute")
		{
			const int rc = cmd_compute(cargc, cargv);
			if (rc == 0) lsquant::report_run_summary();   // peak RSS / CPU / wall + timing regions
			return rc;
		}
		if (cmd == "reconstruct") return cmd_reconstruct(cargc, cargv);
		if (cmd == "run")         return cmd_run(cargc, cargv);
		if (cmd == "merge")       return cmd_merge(cargc, cargv);
		if (cmd == "fingerprint") return cmd_fingerprint(cargc, cargv);
		if (cmd == "--help" || cmd == "-h" || cmd == "help") { usage(cargv[0]); return 0; }
	}
	catch (const lsquant::Error& e)
	{
		LSQ_LOG_ERROR("lsquant: " << e.what());
		return 1;
	}
	catch (const std::exception& e)
	{
		LSQ_LOG_ERROR("lsquant: unhandled exception: " << e.what());
		return 1;
	}

	std::cerr << "lsquant: unknown command '" << cmd << "'\n\n";
	return usage(cargv[0]);
}
