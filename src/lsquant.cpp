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

static int usage(const char* prog)
{
	std::cerr << "usage: " << prog << " <command> [args...]\n\n"
	          << "commands:\n"
	          << "  compute --config <run.json>\n"
	          << "        compute Chebyshev moments for a run (mode: noneq|spectral|msd);\n"
	          << "        writes NonEqOp*.chebmom2D / SpectralOp*.chebmom1D / Correlation*.chebmomTD\n"
	          << "  reconstruct <moments.chebmom2D> <bastin|greenwood> <broadening_meV>\n"
	          << "        reconstruct sigma(E) via the unified Kubo routine; writes Kubo*JACKSON.dat\n"
	          << "  inspect <run.json | operator.desc>\n"
	          << "        print a run / operator and its resolved numerical provenance\n";
	return 2;
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
		<< " reconstruct <moments.chebmom2D> <bastin|greenwood> <broadening_meV>\n"; return 2; }
	const std::string momfile = argv[2];
	const std::string which   = argv[3];
	const double broadening   = std::stod(argv[4]);

	const lsquant::KuboObservable* obs = 0;
	std::string tag;
	if      (which == "bastin")    { obs = &lsquant::KUBO_BASTIN;    tag = "KuboBastin_"; }
	else if (which == "greenwood") { obs = &lsquant::KUBO_GREENWOOD; tag = "KuboGreenWood_"; }
	else { std::cerr << "reconstruct: formula must be 'bastin' or 'greenwood'\n"; return 2; }

	chebyshev::Moments2D mu(momfile.c_str());
	mu.ApplyJacksonKernel(broadening, broadening);
	const std::vector<std::pair<double,double> > data = lsquant::reconstruct_kubo(mu, *obs);

	const std::string outputName = tag + mu.SystemLabel() + "JACKSON.dat";
	std::cout << "Saving the data in " << outputName << std::endl;
	std::ofstream f(outputName.c_str());
	for (size_t i = 0; i < data.size(); ++i) f << data[i].first << " " << data[i].second << std::endl;
	return 0;
}

int main(int argc, char** argv)
{
	if (argc < 2) return usage(argv[0]);
	const std::string cmd = argv[1];
	if (cmd == "inspect")
	{
		if (argc < 3) { std::cerr << "usage: " << argv[0] << " inspect <run.json | operator.desc>\n"; return 2; }
		return lsquant::inspect_path(argv[2]);
	}
	if (cmd == "compute")     return cmd_compute(argc, argv);
	if (cmd == "reconstruct") return cmd_reconstruct(argc, argv);
	if (cmd == "--help" || cmd == "-h" || cmd == "help") { usage(argv[0]); return 0; }
	std::cerr << "lsquant: unknown command '" << cmd << "'\n\n";
	return usage(argv[0]);
}
