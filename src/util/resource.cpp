#include "util/resource.hpp"
#include "util/log.hpp"

#include <chrono>
#include <sstream>
#include <iomanip>
#include <cmath>

// Backend-agnostic parts of the resource probe: the wall-clock origin, the
// formatting helpers, and the probe() selector. The actual memory/CPU readings
// live in the per-platform backends (resource_linux.cpp; stub otherwise).

namespace lsquant { namespace diag {

namespace detail {
// Defined by exactly one backend translation unit.
ResourceProbe& make_probe();
} // namespace detail

namespace {
std::chrono::steady_clock::time_point& origin()
{
	static std::chrono::steady_clock::time_point t = std::chrono::steady_clock::now();
	return t;
}
} // namespace

double wall_since_origin_sec()
{
	const auto now = std::chrono::steady_clock::now();
	return std::chrono::duration_cast<std::chrono::milliseconds>(now - origin()).count() / 1e3;
}

ResourceProbe& probe()
{
	origin();                       // fix wall-clock origin on first use
	return detail::make_probe();
}

ResourceSnapshot sample() { return probe().sample(); }

std::string human_bytes(long long bytes)
{
	if (bytes < 0) return "n/a";
	const char* unit[] = {"B", "KiB", "MiB", "GiB", "TiB"};
	double v = static_cast<double>(bytes);
	int u = 0;
	while (v >= 1024.0 && u < 4) { v /= 1024.0; ++u; }
	std::ostringstream os;
	os << std::fixed << std::setprecision(v < 10 ? 2 : 1) << v << " " << unit[u];
	return os.str();
}

std::string format_snapshot(const ResourceSnapshot& s)
{
	std::ostringstream os;
	os << std::fixed << std::setprecision(1);
	os << "peak RSS " << human_bytes(s.peak_rss_bytes);
	if (s.cpu_user_sec >= 0) os << " | " << s.cpu_user_sec << " s user";
	if (s.cpu_sys_sec  >= 0) os << " | " << s.cpu_sys_sec  << " s sys";
	if (s.wall_sec     >= 0) os << " | " << s.wall_sec     << " s wall";
	return os.str();
}

}} // namespace lsquant::diag
