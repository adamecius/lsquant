#ifndef LSQUANT_IO_RESULTS_HPP
#define LSQUANT_IO_RESULTS_HPP

// Unified results container (unified_io.md, Phase I4/I5).
//
// One <system>.results.json per physical system, machine- and human-readable,
// holding thin-grid final spectra plus timing / errors / provenance. Two docs
// merge IFF their system fingerprints match -- the matured successor to the
// header-string-match python scripts. The doc is a regenerable OUTPUT: gitignored,
// never committed (Phase C data-ban).
//
// Schema ("schema":"lsquant.results/1"):
//   { "schema", "system": {label, fingerprint, size, bandwidth, bandcenter},
//     "results": { "<name>": { ...observable-specific..., "timing_s", "warnings" } },
//     "run": { "wall_s", "peak_rss_bytes", "errors":[] } }

#include <string>
#include "io/json.hpp"

namespace lsquant { namespace io {

// A fresh results doc for a system.
json::Value new_results_doc(const std::string& label, const std::string& fingerprint,
                            long long size, double bandwidth, double bandcenter);

// Insert/replace one observable block under results.<name>.
void add_result(json::Value& doc, const std::string& name, const json::Value& block);

// Merge `incoming` into `base` IN PLACE. Requires equal system.fingerprint;
// throws InputDataError with a clear message otherwise. Observable blocks from
// `incoming` are added (overwriting same-named ones); run.errors are concatenated.
void merge_into(json::Value& base, const json::Value& incoming);

json::Value load_results(const std::string& path);          // throws IOError/ConfigError
void        save_results(const json::Value& doc, const std::string& path);

}} // namespace lsquant::io

#endif // LSQUANT_IO_RESULTS_HPP
