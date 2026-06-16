#ifndef LSQUANT_UTIL_ERROR_HPP
#define LSQUANT_UTIL_ERROR_HPP

// lsquant error vocabulary (reporting_system.md, Phase R3).
//
// Two layers:
//   1. An exception hierarchy (lsquant::Error and friends) for failures at module
//      boundaries -- the lsquant driver installs ONE top-level catch.
//   2. A lightweight Result<T>/Status for hot paths where throwing is undesirable
//      (tight numeric loops, across OpenMP regions). C++11-clean; migrates to
//      std::expected for free when the project leaves C++11.
//
// Policy: throw across module boundaries; return Result within a hot loop; assert
// only for invariants that indicate a bug. Throw sites use LSQ_THROW so file:line
// is captured; the driver's catch logs through the reporting facade.

#include <stdexcept>
#include <string>
#include <sstream>
#include <utility>

namespace lsquant {

// --- exception hierarchy -----------------------------------------------------
class Error : public std::runtime_error {
public:
	explicit Error(const std::string& msg, const std::string& where = std::string())
		: std::runtime_error(where.empty() ? msg : (msg + "  [" + where + "]"))
		, where_(where) {}
	const std::string& where() const { return where_; }
private:
	std::string where_;
};

// Bad run.json / missing field / unknown mode.
struct ConfigError    : Error { ConfigError   (const std::string& m, const std::string& w = std::string()) : Error(m, w) {} };
// File open/read/write failures.
struct IOError        : Error { IOError       (const std::string& m, const std::string& w = std::string()) : Error(m, w) {} };
// Non-Hermitian CSR, bad operator, malformed input data.
struct InputDataError : Error { InputDataError(const std::string& m, const std::string& w = std::string()) : Error(m, w) {} };
// NaN/Inf, non-convergence, out-of-bounds spectral bounds.
struct NumericError   : Error { NumericError  (const std::string& m, const std::string& w = std::string()) : Error(m, w) {} };
// Allocation failure, resource-probe failure.
struct ResourceError  : Error { ResourceError (const std::string& m, const std::string& w = std::string()) : Error(m, w) {} };
// Deliberately-unimplemented path (e.g. the non-orthogonal placeholder preserved
// at 5cc0924). Carries the restore hint in its message.
struct NotImplemented : Error { NotImplemented(const std::string& m, const std::string& w = std::string()) : Error(m, w) {} };

// --- lightweight Status / Result<T> ------------------------------------------
// A movable success-or-message. Use where exceptions cannot cross (OpenMP) or
// where the caller wants to branch without a try/catch.
class Status {
public:
	Status() : ok_(true) {}
	static Status success() { return Status(); }
	static Status fail(const std::string& msg) { Status s; s.ok_ = false; s.msg_ = msg; return s; }
	bool ok()  const { return ok_; }
	explicit operator bool() const { return ok_; }
	const std::string& message() const { return msg_; }
private:
	bool        ok_;
	std::string msg_;
};

template <class T>
class Result {
public:
	Result(const T& v)  : ok_(true), value_(v) {}
	Result(T&& v)       : ok_(true), value_(std::move(v)) {}
	static Result error(const std::string& msg) { Result r; r.ok_ = false; r.msg_ = msg; return r; }

	bool ok() const { return ok_; }
	explicit operator bool() const { return ok_; }
	const std::string& message() const { return msg_; }

	const T& value() const { return value_; }   // precondition: ok()
	T&       value()       { return value_; }
	const T& value_or(const T& fallback) const { return ok_ ? value_ : fallback; }
private:
	Result() : ok_(false) {}
	bool        ok_;
	T           value_;
	std::string msg_;
};

} // namespace lsquant

// Capture file:line:func at the throw site.  LSQ_THROW(ConfigError, "missing 'label'");
#define LSQ_WHERE_STR_2(x) #x
#define LSQ_WHERE_STR(x)   LSQ_WHERE_STR_2(x)
#define LSQ_THROW(ExcType, expr)                                               \
	do {                                                                       \
		std::ostringstream _lsq_oss; _lsq_oss << expr;                         \
		throw ExcType(_lsq_oss.str(),                                          \
		              std::string(__FILE__ ":" LSQ_WHERE_STR(__LINE__)));      \
	} while (0)

#endif // LSQUANT_UTIL_ERROR_HPP
