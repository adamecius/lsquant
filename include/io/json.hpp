#ifndef LSQUANT_IO_JSON_HPP
#define LSQUANT_IO_JSON_HPP

// Minimal nested-JSON value + parser/serializer (unified_io.md, Phase I0).
//
// Why in-house rather than nlohmann/json: the project is offline-safe (FetchContent
// is not guaranteed) AND the data-ban pre-commit hook rejects files > 200 KB, while
// the nlohmann single header is ~930 KB -- so it can be neither vendored nor fetched
// reliably here. This compact, dependency-free reader (objects, arrays, strings,
// numbers, bools, null) succeeds the deliberately flat parser in run_config.cpp and
// is enough for the v2 run schema and the results container. A system nlohmann
// (find_package) can replace it later behind this same interface.

#include <string>
#include <vector>
#include <map>
#include <stdexcept>

namespace lsquant { namespace json {

class Value {
public:
	enum Type { Null, Bool, Number, String, Array, Object };

	Value() : type_(Null), bool_(false), num_(0.0) {}
	Value(bool b) : type_(Bool), bool_(b), num_(0.0) {}
	Value(double d) : type_(Number), bool_(false), num_(d) {}
	Value(int i) : type_(Number), bool_(false), num_(i) {}
	Value(const char* s) : type_(String), bool_(false), num_(0.0), str_(s) {}
	Value(const std::string& s) : type_(String), bool_(false), num_(0.0), str_(s) {}

	static Value array()  { Value v; v.type_ = Array;  return v; }
	static Value object() { Value v; v.type_ = Object; return v; }

	Type type() const { return type_; }
	bool is_null()   const { return type_ == Null; }
	bool is_bool()   const { return type_ == Bool; }
	bool is_number() const { return type_ == Number; }
	bool is_string() const { return type_ == String; }
	bool is_array()  const { return type_ == Array; }
	bool is_object() const { return type_ == Object; }

	bool        as_bool(bool d = false)            const { return is_bool()   ? bool_ : d; }
	double      as_number(double d = 0.0)          const { return is_number() ? num_  : d; }
	int         as_int(int d = 0)                  const { return is_number() ? static_cast<int>(num_) : d; }
	std::string as_string(const std::string& d="") const { return is_string() ? str_  : d; }

	// Object access.
	bool contains(const std::string& k) const {
		return type_ == Object && obj_.find(k) != obj_.end();
	}
	const Value& at(const std::string& k) const {
		static const Value null_value;
		std::map<std::string, Value>::const_iterator it = obj_.find(k);
		return it == obj_.end() ? null_value : it->second;
	}
	const Value& operator[](const std::string& k) const { return at(k); }
	Value&       operator[](const std::string& k) { type_ = Object; return obj_[k]; }
	const std::map<std::string, Value>& items() const { return obj_; }

	// Array access.
	size_t size() const { return type_ == Array ? arr_.size() : 0; }
	const Value& operator[](size_t i) const { return arr_[i]; }
	const std::vector<Value>& elements() const { return arr_; }
	void push_back(const Value& v) { type_ = Array; arr_.push_back(v); }
	void set(const std::string& k, const Value& v) { type_ = Object; obj_[k] = v; }

private:
	Type type_;
	bool bool_;
	double num_;
	std::string str_;
	std::vector<Value> arr_;
	std::map<std::string, Value> obj_;
};

// Parse text into a Value. Throws lsquant::ConfigError with a position-tagged
// message on malformed input.
Value parse(const std::string& text);

// Parse a file. Throws lsquant::IOError if it cannot be opened.
Value parse_file(const std::string& path);

// Serialise. indent > 0 pretty-prints with that many spaces per level.
std::string dump(const Value& v, int indent = 2);

}} // namespace lsquant::json

#endif // LSQUANT_IO_JSON_HPP
