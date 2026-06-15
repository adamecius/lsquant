// lsquant inspect -- print what a run.json (and an operator .desc) will actually do, including
// the resolved NUMERICAL PROVENANCE (alpha, bounds source, hbar/units, kernel). Replaces having
// to read positional argv and decode filename metadata to understand a run.
//
// Usage:
//   lsquant_inspect <run.json>            inspect a run config
//   lsquant_inspect <operator.desc>       inspect an operator descriptor
#include <iostream>
#include <fstream>
#include <string>
#include "run_config.hpp"
#include "operator_descriptor.hpp"
#include "units.hpp"
#include "recon_grid.hpp"

static bool ends_with(const std::string& s, const std::string& suf) {
    return s.size() >= suf.size() && s.compare(s.size()-suf.size(), suf.size(), suf) == 0;
}
static bool file_exists(const std::string& p) { std::ifstream f(p.c_str()); return f.good(); }

static void inspect_desc(const std::string& path) {
    lsquant::OperatorDescriptor d = lsquant::read_descriptor(path);
    std::cout << "operator descriptor: " << path << "\n";
    if (!d.valid) { std::cout << "  (could not parse)\n"; return; }
    std::cout << "  observable : " << d.observable << "\n"
              << "  component  : " << (d.component.empty() ? "(none)" : d.component) << "\n"
              << "  units      : " << d.units << "\n"
              << "  provenance : " << d.provenance << "\n";
    if (d.has_bounds) std::cout << "  bounds(a,b): [" << d.a << ", " << d.b << "]\n";
    if (d.has_basis)  { std::cout << "  basis      : " << d.orbitals_per_cell << " orbitals/cell x "
                                  << d.num_cells << " cells, spins "; for (int s:d.orbital_spin) std::cout<<(s>0?'+':'-'); std::cout<<"\n"; }
}

static void inspect_run(const std::string& path) {
    lsquant::RunConfig c = lsquant::read_run_config(path);
    std::cout << "run config: " << path << "\n";
    if (!c.valid) { std::cout << "  INVALID: " << c.error << "\n"; return; }

    const std::string ham_desc = "operators/" + c.label + ".HAM.desc";
    const std::string bounds_src =
        file_exists("BOUNDS")  ? "BOUNDS file (explicit override)" :
        file_exists(ham_desc)  ? ("descriptor " + ham_desc + " (producer Lanczos)") :
                                 "Gershgorin enclosure + pad (fallback)";
    // alpha: run.json value if set, else the recon_grid default/env
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
    std::cout << "  output      : "
              << (c.output.empty() ? "(derived from run parameters)" : c.output) << "\n";
}

int main(int argc, char** argv) {
    if (argc < 2) {
        std::cerr << "usage: " << argv[0] << " <run.json | operator.desc>\n";
        return 2;
    }
    const std::string path = argv[1];
    if (ends_with(path, ".desc")) inspect_desc(path);
    else                          inspect_run(path);
    return 0;
}
