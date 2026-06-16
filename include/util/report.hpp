#ifndef LSQUANT_UTIL_REPORT_HPP
#define LSQUANT_UTIL_REPORT_HPP

// lsquant reporting layer -- umbrella include + one-line initialisation.
//
//   #include "util/report.hpp"
//   int main(int argc, char** argv) {
//       lsquant::init_reporting(argc, argv);   // honours -v / --log-file + env
//       ... LSQ_LOG_INFO("hello"); ...
//   }
//
// See reporting_system.md. Pulls in logging, timing, errors, and resource
// accounting; each header is independently usable.

#include "util/log.hpp"
#include "util/timer.hpp"
#include "util/error.hpp"
#include "util/resource.hpp"

namespace lsquant {

// Initialise the reporting layer. Resolves the log level and file sink from, in
// increasing precedence: defaults -> env (LSQUANT_LOG_LEVEL, LSQUANT_LOG_FILE)
// -> CLI flags scanned from argv (-v/-vv/--verbose, --quiet, --log-level <x>,
// --log-file <path>). Recognised flags are left in argv (callers may ignore
// them); returns nothing. Idempotent.
void init_reporting(int argc = 0, char** argv = nullptr);

// Log a final resource summary line (peak RSS, CPU, wall) at info level, plus the
// timing-region summary. Safe to call once near process exit.
void report_run_summary();

} // namespace lsquant

#endif // LSQUANT_UTIL_REPORT_HPP
