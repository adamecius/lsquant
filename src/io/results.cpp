#include "io/results.hpp"
#include "util/error.hpp"

#include <fstream>

namespace lsquant { namespace io {

json::Value new_results_doc(const std::string& label, const std::string& fingerprint,
                            long long size, double bandwidth, double bandcenter)
{
	json::Value doc = json::Value::object();
	doc.set("schema", json::Value("lsquant.results/1"));

	json::Value sys = json::Value::object();
	sys.set("label",       json::Value(label));
	sys.set("fingerprint", json::Value(fingerprint));
	sys.set("size",        json::Value(static_cast<double>(size)));
	sys.set("bandwidth",   json::Value(bandwidth));
	sys.set("bandcenter",  json::Value(bandcenter));
	doc.set("system", sys);

	doc.set("results", json::Value::object());

	json::Value run = json::Value::object();
	run.set("errors", json::Value::array());
	doc.set("run", run);
	return doc;
}

void add_result(json::Value& doc, const std::string& name, const json::Value& block)
{
	json::Value results = doc.contains("results") ? doc.at("results") : json::Value::object();
	results.set(name, block);
	doc.set("results", results);
}

void merge_into(json::Value& base, const json::Value& incoming)
{
	const std::string fb = base.at("system").at("fingerprint").as_string();
	const std::string fi = incoming.at("system").at("fingerprint").as_string();
	if (fb.empty() || fi.empty())
		throw lsquant::InputDataError("merge: missing system.fingerprint in one of the documents");
	if (fb != fi)
		throw lsquant::InputDataError(
			"merge refused: different systems (fingerprint " + fb.substr(0, 12) +
			"... vs " + fi.substr(0, 12) + "...). These runs are not on the same system.");

	// Merge observable blocks.
	json::Value results = base.contains("results") ? base.at("results") : json::Value::object();
	const json::Value& inc_results = incoming.at("results");
	for (std::map<std::string, json::Value>::const_iterator it = inc_results.items().begin();
	     it != inc_results.items().end(); ++it)
		results.set(it->first, it->second);
	base.set("results", results);

	// Concatenate run.errors.
	json::Value run = base.contains("run") ? base.at("run") : json::Value::object();
	json::Value errors = run.contains("errors") ? run.at("errors") : json::Value::array();
	const json::Value& inc_errors = incoming.at("run").at("errors");
	for (size_t k = 0; k < inc_errors.size(); ++k) errors.push_back(inc_errors[k]);
	run.set("errors", errors);
	base.set("run", run);
}

json::Value load_results(const std::string& path)
{
	json::Value doc = json::parse_file(path);   // throws IOError/ConfigError
	if (!doc.is_object() || !doc.contains("system"))
		throw lsquant::ConfigError("not a results document (no 'system' object): " + path);
	return doc;
}

void save_results(const json::Value& doc, const std::string& path)
{
	std::ofstream f(path.c_str());
	if (!f.is_open()) throw lsquant::IOError("cannot write " + path);
	f << json::dump(doc, 2) << "\n";
}

}} // namespace lsquant::io
