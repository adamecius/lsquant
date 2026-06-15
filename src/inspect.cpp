#include "inspect.hpp"

#include <iostream>
#include <fstream>
#include <string>
#include "run_config.hpp"
#include "operator_descriptor.hpp"
#include "units.hpp"
#include "recon_grid.hpp"

namespace
{
	bool ends_with(const std::string& s, const std::string& suf) {
		return s.size() >= suf.size() && s.compare(s.size()-suf.size(), suf.size(), suf) == 0;
	}
	bool file_exists(const std::string& p) { std::ifstream f(p.c_str()); return f.good(); }

	void inspect_desc(const std::string& path) {
		lsquant::OperatorDescriptor d = lsquant::read_descriptor(path);
		std::cout << "operator descriptor: " << path << "\n";
		if (!d.valid) { std::cout << "  (could not parse)\n"; return; }
		std::cout << "  observable : " << d.observable << "\n"
		          << "  component  : " << (d.component.empty() ? "(none)" : d.component) << "\n"
		          << "  units      : " << d.units << "\n"
		          << "  provenance : " << d.provenance << "\n";
		if (d.has_bounds) std::cout << "  bounds(a,b): [" << d.a << ", " << d.b << "]\n";
		if (d.has_basis)  { std::cout << "  basis      : " << d.orbitals_per_cell << " orbitals/cell x "
		                              << d.num_cells << " cells, spins ";
		                    for (int s : d.orbital_spin) std::cout << (s>0?'+':'-'); std::cout << "\n"; }
	}

	void inspect_run(const std::string& path) {
		lsquant::RunConfig c = lsquant::read_run_config(path);
		std::cout << "run config: " << path << "\n";
		if (!c.valid) { std::cout << "  INVALID: " << c.error << "\n"; return; }

		const std::string ham_desc = "operators/" + c.label + ".HAM.desc";
		const std::string bounds_src =
			file_exists("BOUNDS")  ? "BOUNDS file (explicit override)" :
			file_exists(ham_desc)  ? ("descriptor " + ham_desc + " (producer Lanczos)") :
			                         "Gershgorin enclosure + pad (fallback)";
		const double alpha = (c.alpha >= 0.0) ? c.alpha : chebyshev::safety_factors().recon_cutoff;

		std::cout << "  what runs   : <" << c.operator_left << "|.|" << c.operator_right
		          << "> correlation on '" << c.label << "', M = " << c.num_moments << "\n"
		          << "  state       : " << c.state << "\n"
		          << "  kernel      : " << c.kernel;
		if (c.kernel == "lorentz") std::cout << " (lambda = " << c.lambda << ")";
		std::cout << "\n";
		std::cout << "  -- numerical provenance --\n"
		          << "  recon alpha : " << alpha << (c.alpha>=0.0 ? " (from run.json)" : " (recon_grid default/env)") << "\n"
		          << "  bounds pad  : " << chebyshev::safety_factors().bounds_pad << "\n"
		          << "  bounds from : " << bounds_src << "\n"
		          << "  hbar (units): " << chebyshev::HBAR << " eV*fs\n";
		if (file_exists(ham_desc)) { std::cout << "  -- Hamiltonian descriptor --\n"; inspect_desc(ham_desc); }
		std::cout << "  output      : " << (c.output.empty() ? "(derived from run parameters)" : c.output) << "\n";
	}
}

namespace lsquant
{
	int inspect_path(const std::string& path)
	{
		if (ends_with(path, ".desc")) inspect_desc(path);
		else                          inspect_run(path);
		return 0;
	}
}
