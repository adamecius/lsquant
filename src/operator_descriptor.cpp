#include "operator_descriptor.hpp"

#include <fstream>
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
}
