#ifndef LSQUANT_UTIL_RESOURCE_HPP
#define LSQUANT_UTIL_RESOURCE_HPP

// lsquant resource accounting (reporting_system.md, Phase R4).
//
// Cheap sampling of memory (current/peak RSS) and CPU (user/sys) plus wall time.
// Ubuntu/Linux-first, but behind a ResourceProbe interface so macOS (task_info)
// and Windows (PSAPI) backends slot in later without touching call sites. On a
// platform with no backend, fields are -1 and callers degrade gracefully.

#include <string>

namespace lsquant { namespace diag {

struct ResourceSnapshot {
	long long rss_bytes      = -1;  // current resident set size
	long long peak_rss_bytes = -1;  // high-water mark
	double    cpu_user_sec   = -1;  // user CPU time
	double    cpu_sys_sec    = -1;  // system CPU time
	double    wall_sec       = -1;  // since first probe() / process start
	long long avail_bytes    = -1;  // system MemAvailable (optional, -1 if unknown)
};

struct ResourceProbe {
	virtual ResourceSnapshot sample() = 0;
	virtual const char*      backend_name() const = 0;
	virtual ~ResourceProbe() {}
};

// Process-wide probe (Linux backend on Ubuntu; stub elsewhere). The first call
// fixes the wall-clock origin.
ResourceProbe&   probe();
ResourceSnapshot sample();

// Human-readable one-liner, e.g. "peak RSS 4.20 GiB | 312.0 s user | 18.0 s sys
// | 95.0 s wall". Used by the lsquant driver's end-of-run summary.
std::string format_snapshot(const ResourceSnapshot& s);
std::string human_bytes(long long bytes);

}} // namespace lsquant::diag

#endif // LSQUANT_UTIL_RESOURCE_HPP
