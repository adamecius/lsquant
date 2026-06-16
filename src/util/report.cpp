#include "util/report.hpp"

#include <cstdlib>
#include <cstring>
#include <string>

namespace lsquant {

void init_reporting(int argc, char** argv)
{
	// 1) env
	if (const char* lv = std::getenv("LSQUANT_LOG_LEVEL"))
		log::set_level(log::level_from_string(lv));
	if (const char* lf = std::getenv("LSQUANT_LOG_FILE"))
		if (lf[0]) log::set_log_file(lf);

	// 2) CLI flags (scanned, not consumed)
	for (int i = 1; i < argc && argv; ++i) {
		const std::string a = argv[i];
		if (a == "--verbose" || a == "-v") log::set_level(log::Level::debug);
		else if (a == "-vv")               log::set_level(log::Level::trace);
		else if (a == "--quiet" || a == "-q") log::set_level(log::Level::warn);
		else if (a == "--log-level" && i + 1 < argc)
			log::set_level(log::level_from_string(argv[++i]));
		else if (a == "--log-file" && i + 1 < argc)
			log::set_log_file(argv[++i]);
	}

	LSQ_LOG_DEBUG("reporting initialised (level=" << static_cast<int>(log::level())
	              << ", resource backend=" << diag::probe().backend_name() << ")");
}

void report_run_summary()
{
	util::log_region_summary(log::Level::info);
	const diag::ResourceSnapshot s = diag::sample();
	LSQ_LOG_INFO(diag::format_snapshot(s));
}

} // namespace lsquant
