#include "util/resource.hpp"
#include "util/log.hpp"

// Fallback resource backend for platforms without a native probe (non-Linux).
// Returns -1 fields and warns once. A macOS (task_info) or Windows (PSAPI)
// backend can replace this without touching any call site. Compiled only when
// CMake does NOT select the Linux backend.

namespace lsquant { namespace diag {

double wall_since_origin_sec();   // defined in resource.cpp

namespace {
struct StubProbe : ResourceProbe {
	ResourceSnapshot sample() override {
		static bool warned = false;
		if (!warned) {
			warned = true;
			LSQ_LOG_WARN("resource probe: no backend on this platform; memory/CPU unavailable");
		}
		ResourceSnapshot s;
		s.wall_sec = wall_since_origin_sec();
		return s;
	}
	const char* backend_name() const override { return "stub(none)"; }
};
} // namespace

namespace detail {
ResourceProbe& make_probe() { static StubProbe p; return p; }
} // namespace detail

}} // namespace lsquant::diag
