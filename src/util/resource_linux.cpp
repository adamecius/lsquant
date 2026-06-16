#include "util/resource.hpp"

// Linux/Ubuntu resource backend: /proc/self/status (VmRSS, VmHWM) and
// getrusage(RUSAGE_SELF) for CPU and peak RSS, plus /proc/meminfo MemAvailable.
// Compiled only when CMake selects the Linux backend (UNIX AND NOT APPLE).

#include <cstdio>
#include <cstring>
#include <string>
#include <sys/resource.h>

namespace lsquant { namespace diag {

double wall_since_origin_sec();   // defined in resource.cpp

namespace {

// Read a "Key:   <num> kB" line from /proc/self/status -> bytes (-1 if absent).
long long read_status_kb(const char* key)
{
	std::FILE* f = std::fopen("/proc/self/status", "r");
	if (!f) return -1;
	char line[256];
	long long bytes = -1;
	const size_t klen = std::strlen(key);
	while (std::fgets(line, sizeof(line), f)) {
		if (std::strncmp(line, key, klen) == 0) {
			long long kb = 0;
			if (std::sscanf(line + klen, " %lld", &kb) == 1) bytes = kb * 1024LL;
			break;
		}
	}
	std::fclose(f);
	return bytes;
}

long long read_mem_available()
{
	std::FILE* f = std::fopen("/proc/meminfo", "r");
	if (!f) return -1;
	char line[256];
	long long bytes = -1;
	while (std::fgets(line, sizeof(line), f)) {
		if (std::strncmp(line, "MemAvailable:", 13) == 0) {
			long long kb = 0;
			if (std::sscanf(line + 13, " %lld", &kb) == 1) bytes = kb * 1024LL;
			break;
		}
	}
	std::fclose(f);
	return bytes;
}

struct LinuxProbe : ResourceProbe {
	ResourceSnapshot sample() override {
		ResourceSnapshot s;
		s.rss_bytes      = read_status_kb("VmRSS:");
		s.peak_rss_bytes = read_status_kb("VmHWM:");
		s.avail_bytes    = read_mem_available();

		struct rusage ru;
		if (getrusage(RUSAGE_SELF, &ru) == 0) {
			s.cpu_user_sec = ru.ru_utime.tv_sec + ru.ru_utime.tv_usec / 1e6;
			s.cpu_sys_sec  = ru.ru_stime.tv_sec + ru.ru_stime.tv_usec / 1e6;
			// ru_maxrss is KiB on Linux; use it as a peak-RSS fallback.
			if (s.peak_rss_bytes < 0 && ru.ru_maxrss > 0)
				s.peak_rss_bytes = static_cast<long long>(ru.ru_maxrss) * 1024LL;
		}
		s.wall_sec = wall_since_origin_sec();
		return s;
	}
	const char* backend_name() const override { return "linux(/proc+getrusage)"; }
};

} // namespace

namespace detail {
ResourceProbe& make_probe() { static LinuxProbe p; return p; }
} // namespace detail

}} // namespace lsquant::diag
