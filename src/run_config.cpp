#include "run_config.hpp"

#include <fstream>
#include <sstream>
#include <cctype>

namespace
{
	// Minimal flat-JSON object reader: extracts top-level "key": value pairs.
	// Values may be quoted strings, numbers, or true/false/null. Nested objects/arrays are
	// NOT supported (the v1 run.json schema is flat). Tolerant of whitespace and trailing commas.
	std::map<std::string, std::string> parse(const std::string& text)
	{
		std::map<std::string, std::string> kv;
		size_t i = text.find('{');
		if (i == std::string::npos) return kv;
		++i;
		const size_t n = text.size();
		auto skipws = [&] { while (i < n && std::isspace((unsigned char)text[i])) ++i; };

		while (i < n)
		{
			skipws();
			if (i >= n || text[i] == '}') break;
			if (text[i] == ',') { ++i; continue; }
			if (text[i] != '"') { ++i; continue; }      // expect a quoted key
			// --- key ---
			++i; std::string key;
			while (i < n && text[i] != '"') { if (text[i]=='\\' && i+1<n) ++i; key += text[i++]; }
			if (i < n) ++i;                              // closing quote
			skipws();
			if (i < n && text[i] == ':') ++i;
			skipws();
			// --- value ---
			std::string val;
			if (i < n && text[i] == '"')
			{
				++i;
				while (i < n && text[i] != '"') { if (text[i]=='\\' && i+1<n) ++i; val += text[i++]; }
				if (i < n) ++i;
			}
			else
			{
				while (i < n && text[i] != ',' && text[i] != '}' && !std::isspace((unsigned char)text[i]))
					val += text[i++];
			}
			if (!key.empty()) kv[key] = val;
		}
		return kv;
	}
}

namespace lsquant
{
	std::map<std::string, std::string> parse_flat_json(const std::string& path)
	{
		std::ifstream f(path.c_str());
		if (!f.is_open()) return {};
		std::stringstream ss; ss << f.rdbuf();
		return parse(ss.str());
	}

	RunConfig read_run_config(const std::string& path)
	{
		RunConfig c;
		auto kv = parse_flat_json(path);
		if (kv.empty()) { c.error = "cannot open or parse " + path; return c; }

		auto get = [&](const char* k, const std::string& dflt) -> std::string {
			auto it = kv.find(k); return it == kv.end() ? dflt : it->second;
		};
		c.label          = get("label", "");
		c.operator_left  = get("operator_left", "");
		c.operator_right = get("operator_right", "");
		c.kernel         = get("kernel", "jackson");
		c.state          = get("state", "default");
		c.output         = get("output", "");
		try {
			if (kv.count("num_moments")) c.num_moments = std::stoi(kv["num_moments"]);
			if (kv.count("lambda"))      c.lambda      = std::stod(kv["lambda"]);
			if (kv.count("alpha"))       c.alpha       = std::stod(kv["alpha"]);
		} catch (...) { c.error = "non-numeric value for num_moments/lambda/alpha"; return c; }

		if (c.label.empty())       { c.error = "missing 'label'"; return c; }
		if (c.num_moments <= 0)    { c.error = "missing/invalid 'num_moments'"; return c; }
		c.valid = true;
		return c;
	}
}
