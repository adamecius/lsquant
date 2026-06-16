#include "util/timer.hpp"

#include <map>
#include <mutex>
#include <sstream>
#include <iomanip>

namespace lsquant { namespace util {

namespace {
struct Region { long long total_us = 0; long long count = 0; };
std::map<std::string, Region> g_regions;
std::mutex                    g_mutex;
} // namespace

void region_add(const std::string& name, long long microseconds)
{
	std::lock_guard<std::mutex> lk(g_mutex);
	Region& r = g_regions[name];
	r.total_us += microseconds;
	r.count    += 1;
}

void clear_regions()
{
	std::lock_guard<std::mutex> lk(g_mutex);
	g_regions.clear();
}

void log_region_summary(log::Level lvl)
{
	if (!log::enabled(lvl)) return;
	std::map<std::string, Region> snapshot;
	{
		std::lock_guard<std::mutex> lk(g_mutex);
		snapshot = g_regions;
	}
	if (snapshot.empty()) return;

	LSQ_LOG(lvl, "--- timing summary (region: total / calls / mean) ---");
	for (std::map<std::string, Region>::const_iterator it = snapshot.begin();
	     it != snapshot.end(); ++it)
	{
		const double total_ms = static_cast<double>(it->second.total_us) / 1e3;
		const double mean_ms   = it->second.count
			? total_ms / static_cast<double>(it->second.count) : 0.0;
		std::ostringstream line;
		line << std::left << std::setw(40) << it->first
		     << std::right << std::setw(12) << total_ms << " ms"
		     << std::setw(8)  << it->second.count << " calls"
		     << std::setw(12) << mean_ms << " ms/call";
		LSQ_LOG(lvl, line.str());
	}
}

}} // namespace lsquant::util
