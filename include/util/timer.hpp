#ifndef LSQUANT_UTIL_TIMER_HPP
#define LSQUANT_UTIL_TIMER_HPP

// lsquant timing utilities (reporting_system.md, Phase R2).
//
// One canonical wall-clock timer plus an RAII scoped timer and a named-region
// registry, all emitting through the logging facade rather than std::cout. This
// supersedes the two legacy time_station classes (include/time_station.hpp and
// src/time_station_2.hpp), which now forward here.

#include <string>
#include <chrono>
#include "util/log.hpp"

namespace lsquant { namespace util {

// Monotonic wall-clock stopwatch. Accumulates across stop()/start() cycles.
class Timer {
public:
	Timer() : accum_us_(0), running_(false) { start(); }

	void start() {
		start_ = std::chrono::steady_clock::now();
		running_ = true;
	}
	// Stop and fold the interval into the accumulator; returns total elapsed us.
	long long stop() {
		if (running_) {
			const auto now = std::chrono::steady_clock::now();
			accum_us_ += std::chrono::duration_cast<std::chrono::microseconds>(now - start_).count();
			running_ = false;
		}
		return accum_us_;
	}
	void reset() { accum_us_ = 0; running_ = false; }

	long long elapsed_us() const {
		long long e = accum_us_;
		if (running_) {
			const auto now = std::chrono::steady_clock::now();
			e += std::chrono::duration_cast<std::chrono::microseconds>(now - start_).count();
		}
		return e;
	}
	double elapsed_ms() const { return static_cast<double>(elapsed_us()) / 1e3; }
	double elapsed_s()  const { return static_cast<double>(elapsed_us()) / 1e6; }

private:
	std::chrono::steady_clock::time_point start_;
	long long accum_us_;
	bool      running_;
};

// Record one named region's elapsed time into the process-wide registry
// (count + total), e.g. for per-iteration accumulation. Summarised by
// log_region_summary() at shutdown.
void  region_add(const std::string& name, long long microseconds);
void  log_region_summary(log::Level lvl = log::Level::info);
void  clear_regions();

// RAII: times its scope and, on destruction, logs "<label>: <ms> ms" at `lvl`
// and folds the interval into the named region registry.
class ScopedTimer {
public:
	explicit ScopedTimer(const std::string& label, log::Level lvl = log::Level::debug)
		: label_(label), lvl_(lvl) {}
	~ScopedTimer() {
		const long long us = timer_.stop();
		region_add(label_, us);
		LSQ_LOG(lvl_, label_ << ": " << (static_cast<double>(us) / 1e3) << " ms");
	}
private:
	std::string label_;
	log::Level  lvl_;
	Timer       timer_;
};

}} // namespace lsquant::util

// Convenience: one-line scope instrumentation.  LSQ_SCOPED("CorrelationExpansion");
#define LSQ_SCOPED(label) ::lsquant::util::ScopedTimer _lsq_scoped_timer((label))

#endif // LSQUANT_UTIL_TIMER_HPP
