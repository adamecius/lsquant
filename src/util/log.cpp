#include "util/log.hpp"

#include <atomic>
#include <chrono>
#include <ctime>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <mutex>
#include <thread>
#if defined(__unix__) || defined(__APPLE__)
#include <unistd.h>
#endif

// In-house, zero-dependency logging backend (offline-safe). Thread-safe; emits
// to stderr (optionally coloured on a TTY) and, when configured, appends to a
// file sink. Format: "[HH:MM:SS.mmm] [LEVEL] [tid] message".

namespace lsquant { namespace log {

namespace {

std::atomic<int> g_level{ static_cast<int>(Level::info) };
std::mutex       g_mutex;          // serialises both sinks
std::ofstream    g_file;
bool             g_file_open = false;

const char* level_name(Level l)
{
	switch (l) {
		case Level::trace: return "TRACE";
		case Level::debug: return "DEBUG";
		case Level::info:  return "INFO ";
		case Level::warn:  return "WARN ";
		case Level::error: return "ERROR";
		default:           return "OFF  ";
	}
}

// ANSI colour per level for the stderr sink (only when stderr is a terminal).
const char* level_colour(Level l)
{
	switch (l) {
		case Level::trace: return "\033[90m";   // bright black
		case Level::debug: return "\033[36m";   // cyan
		case Level::info:  return "\033[32m";   // green
		case Level::warn:  return "\033[33m";   // yellow
		case Level::error: return "\033[31m";   // red
		default:           return "";
	}
}

bool stderr_is_tty()
{
#if defined(__unix__) || defined(__APPLE__)
	return ::isatty(2) != 0;
#else
	return false;
#endif
}

std::string timestamp()
{
	using namespace std::chrono;
	const auto now  = system_clock::now();
	const auto secs = time_point_cast<seconds>(now);
	const auto ms   = duration_cast<milliseconds>(now - secs).count();
	const std::time_t t = system_clock::to_time_t(now);
	std::tm tmv;
#if defined(_WIN32)
	localtime_s(&tmv, &t);
#else
	localtime_r(&t, &tmv);
#endif
	char buf[16];
	std::snprintf(buf, sizeof(buf), "%02d:%02d:%02d.%03d",
	              tmv.tm_hour, tmv.tm_min, tmv.tm_sec, static_cast<int>(ms));
	return std::string(buf);
}

} // namespace

void  set_level(Level lvl) { g_level.store(static_cast<int>(lvl)); }
Level level()              { return static_cast<Level>(g_level.load()); }
bool  enabled(Level lvl)   { return static_cast<int>(lvl) >= g_level.load(); }

Level level_from_string(const std::string& s)
{
	std::string v; v.reserve(s.size());
	for (char c : s) v += static_cast<char>(std::tolower((unsigned char)c));
	if (v == "trace") return Level::trace;
	if (v == "debug") return Level::debug;
	if (v == "info")  return Level::info;
	if (v == "warn" || v == "warning") return Level::warn;
	if (v == "error") return Level::error;
	if (v == "off" || v == "none" || v == "silent") return Level::off;
	return Level::info;
}

bool set_log_file(const std::string& path)
{
	std::lock_guard<std::mutex> lk(g_mutex);
	if (g_file.is_open()) g_file.close();
	g_file_open = false;
	if (path.empty()) return true;
	g_file.open(path.c_str(), std::ios::out | std::ios::app);
	g_file_open = g_file.is_open();
	return g_file_open;
}

void emit(Level lvl, const std::string& msg)
{
	if (!enabled(lvl)) return;
	const std::string ts  = timestamp();
	const char* lname     = level_name(lvl);

	std::ostringstream tid; tid << std::this_thread::get_id();
	std::lock_guard<std::mutex> lk(g_mutex);

	static const bool colour = stderr_is_tty();
	if (colour)
		std::cerr << "[" << ts << "] " << level_colour(lvl) << "[" << lname << "]"
		          << "\033[0m [" << tid.str() << "] " << msg << "\n";
	else
		std::cerr << "[" << ts << "] [" << lname << "] [" << tid.str() << "] " << msg << "\n";

	if (g_file_open)
		g_file << "[" << ts << "] [" << lname << "] [" << tid.str() << "] " << msg << "\n" << std::flush;
}

}} // namespace lsquant::log
