#ifndef LSQUANT_UTIL_LOG_HPP
#define LSQUANT_UTIL_LOG_HPP

// lsquant logging facade (reporting_system.md, Phase R1).
//
// A thin, levelled, thread-safe logger. Call sites use the LSQ_LOG_* macros and
// never include a backend header, so the implementation (currently an in-house,
// zero-dependency, offline-safe sink set) can be swapped for spdlog later behind
// LINQT_WITH_SPDLOG without touching a single call site.
//
//   LSQ_LOG_INFO("saved " << n << " moments to " << name);
//   LSQ_LOG_WARN("non-Hermitian input, symmetrising");
//
// The argument is an ostream expression, so it composes like std::cout but is
// only evaluated when the level is enabled (a cheap atomic load otherwise).

#include <string>
#include <sstream>

namespace lsquant { namespace log {

enum class Level { trace = 0, debug = 1, info = 2, warn = 3, error = 4, off = 5 };

// Process-wide configuration. Safe to call before/after init_reporting().
void  set_level(Level lvl);
Level level();
bool  enabled(Level lvl);            // lvl >= current threshold

// Parse "trace|debug|info|warn|error|off" (case-insensitive); returns info on
// an unrecognised string. Used by the LSQUANT_LOG_LEVEL env and the -v flag.
Level level_from_string(const std::string& s);

// Route a copy of every record to an additional file sink (appended). Empty path
// disables the file sink. The file is a regenerable artifact -- gitignored, never
// committed (Phase C data-ban). Returns false if the file could not be opened.
bool  set_log_file(const std::string& path);

// Emit a single already-formatted record. Prefer the macros below.
void  emit(Level lvl, const std::string& msg);

}} // namespace lsquant::log

// --- call-site macros --------------------------------------------------------
#define LSQ_LOG(lvl, expr)                                                     \
	do {                                                                       \
		if (::lsquant::log::enabled(lvl)) {                                    \
			std::ostringstream _lsq_oss; _lsq_oss << expr;                     \
			::lsquant::log::emit((lvl), _lsq_oss.str());                       \
		}                                                                      \
	} while (0)

#define LSQ_LOG_TRACE(expr) LSQ_LOG(::lsquant::log::Level::trace, expr)
#define LSQ_LOG_DEBUG(expr) LSQ_LOG(::lsquant::log::Level::debug, expr)
#define LSQ_LOG_INFO(expr)  LSQ_LOG(::lsquant::log::Level::info,  expr)
#define LSQ_LOG_WARN(expr)  LSQ_LOG(::lsquant::log::Level::warn,  expr)
#define LSQ_LOG_ERROR(expr) LSQ_LOG(::lsquant::log::Level::error, expr)

#endif // LSQUANT_UTIL_LOG_HPP
