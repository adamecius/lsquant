#include "operator_descriptor.hpp"

#include <fstream>
#include <sstream>
#include <string>

namespace
{
	// trim leading/trailing spaces, tabs and CR (tolerate CRLF line endings)
	std::string trim(std::string s)
	{
		const char* ws = " \t\r\n";
		const auto b = s.find_first_not_of(ws);
		if (b == std::string::npos) return "";
		const auto e = s.find_last_not_of(ws);
		return s.substr(b, e - b + 1);
	}
}

namespace lsquant
{
	OperatorDescriptor read_descriptor(const std::string& path)
	{
		OperatorDescriptor d;
		std::ifstream f(path.c_str());
		if (!f.is_open()) return d;   // valid stays false

		std::string line;
		while (std::getline(f, line))
		{
			const auto colon = line.find(':');
			if (colon == std::string::npos) continue;
			const std::string key = trim(line.substr(0, colon));
			const std::string val = trim(line.substr(colon + 1));

			if      (key == "observable")   d.observable  = val;
			else if (key == "component")    d.component   = val;
			else if (key == "units")        d.units       = val;
			else if (key == "provenance")   d.provenance  = val;
			else if (key == "spectral_min") { d.a = std::stod(val); d.has_bounds = true; }
			else if (key == "spectral_max") { d.b = std::stod(val); d.has_bounds = true; }
			else if (key == "basis_orbitals_per_cell") { d.orbitals_per_cell = std::stoi(val); d.has_basis = true; }
			else if (key == "basis_num_cells")         { d.num_cells = std::stoi(val);         d.has_basis = true; }
			else if (key == "basis_orbital_spin")
			{
				// space-separated +/- (or +1/-1) signs, one per orbital in a cell
				std::istringstream ss(val); std::string tok;
				d.orbital_spin.clear();
				while (ss >> tok) d.orbital_spin.push_back((!tok.empty() && tok[0] == '-') ? -1 : +1);
				d.has_basis = true;
			}
			// unknown keys are ignored (forward-compatible with later schema versions)
		}
		d.valid = true;
		return d;
	}

	bool descriptor_bounds(const std::string& desc_path, double& a, double& b)
	{
		const OperatorDescriptor d = read_descriptor(desc_path);
		if (d.valid && d.has_bounds) { a = d.a; b = d.b; return true; }
		return false;
	}

	bool build_spin_operator(const OperatorDescriptor& d, char component,
	                         std::vector<int>& rows, std::vector<int>& cols,
	                         std::vector<std::complex<double> >& vals)
	{
		typedef std::complex<double> cd;
		rows.clear(); cols.clear(); vals.clear();
		if (!d.has_basis || d.orbitals_per_cell <= 0 || d.num_cells <= 0) return false;
		if ((int)d.orbital_spin.size() != d.orbitals_per_cell) return false;

		// Identify the up (+1) and down (-1) orbital within a cell (spin-1/2 doubling).
		int up = -1, dn = -1;
		for (int o = 0; o < d.orbitals_per_cell; ++o)
		{
			if (d.orbital_spin[o] > 0 && up < 0) up = o;
			if (d.orbital_spin[o] < 0 && dn < 0) dn = o;
		}
		if (up < 0 || dn < 0) return false;   // need one up and one down orbital

		const int opc = d.orbitals_per_cell;
		for (int c = 0; c < d.num_cells; ++c)
		{
			const int iu = c * opc + up;   // global index of the up orbital
			const int id = c * opc + dn;   // global index of the down orbital
			switch (component)
			{
				case 'X': case 'x':   // sigma_x = [[0,1],[1,0]]
					rows.push_back(iu); cols.push_back(id); vals.push_back(cd(1, 0));
					rows.push_back(id); cols.push_back(iu); vals.push_back(cd(1, 0));
					break;
				case 'Y': case 'y':   // sigma_y = [[0,-i],[i,0]]
					rows.push_back(iu); cols.push_back(id); vals.push_back(cd(0, -1));
					rows.push_back(id); cols.push_back(iu); vals.push_back(cd(0, +1));
					break;
				case 'Z': case 'z':   // sigma_z = [[1,0],[0,-1]]
					rows.push_back(iu); cols.push_back(iu); vals.push_back(cd(+1, 0));
					rows.push_back(id); cols.push_back(id); vals.push_back(cd(-1, 0));
					break;
				default: return false;
			}
		}
		return true;
	}
}
