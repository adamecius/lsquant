#include "io/json.hpp"
#include "util/error.hpp"

#include <fstream>
#include <sstream>
#include <cctype>
#include <cmath>
#include <cstdio>

namespace lsquant { namespace json {

namespace {

struct Parser {
	const std::string& s;
	size_t i;
	explicit Parser(const std::string& text) : s(text), i(0) {}

	[[noreturn]] void fail(const std::string& msg) {
		// 1-based line/col for a friendly message.
		size_t line = 1, col = 1;
		for (size_t k = 0; k < i && k < s.size(); ++k) {
			if (s[k] == '\n') { ++line; col = 1; } else ++col;
		}
		std::ostringstream os;
		os << "JSON parse error at line " << line << ", col " << col << ": " << msg;
		throw lsquant::ConfigError(os.str());
	}

	void skip_ws() {
		while (i < s.size()) {
			const char c = s[i];
			if (c == ' ' || c == '\t' || c == '\n' || c == '\r') { ++i; continue; }
			// // line comments and /* */ block comments (tolerated, JSONC-style)
			if (c == '/' && i + 1 < s.size()) {
				if (s[i+1] == '/') { i += 2; while (i < s.size() && s[i] != '\n') ++i; continue; }
				if (s[i+1] == '*') { i += 2; while (i + 1 < s.size() && !(s[i]=='*'&&s[i+1]=='/')) ++i; i += 2; continue; }
			}
			break;
		}
	}

	char peek() { return i < s.size() ? s[i] : '\0'; }

	Value parse_value() {
		skip_ws();
		const char c = peek();
		switch (c) {
			case '{': return parse_object();
			case '[': return parse_array();
			case '"': return Value(parse_string());
			case 't': case 'f': return parse_bool();
			case 'n': return parse_null();
			case '\0': fail("unexpected end of input");
			default:  return parse_number();
		}
	}

	Value parse_object() {
		Value v = Value::object();
		++i; // {
		skip_ws();
		if (peek() == '}') { ++i; return v; }
		for (;;) {
			skip_ws();
			if (peek() != '"') fail("expected string key");
			const std::string key = parse_string();
			skip_ws();
			if (peek() != ':') fail("expected ':' after key");
			++i;
			v.set(key, parse_value());
			skip_ws();
			const char d = peek();
			if (d == ',') { ++i; skip_ws(); if (peek() == '}') { ++i; return v; } continue; }
			if (d == '}') { ++i; return v; }
			fail("expected ',' or '}' in object");
		}
	}

	Value parse_array() {
		Value v = Value::array();
		++i; // [
		skip_ws();
		if (peek() == ']') { ++i; return v; }
		for (;;) {
			v.push_back(parse_value());
			skip_ws();
			const char d = peek();
			if (d == ',') { ++i; skip_ws(); if (peek() == ']') { ++i; return v; } continue; }
			if (d == ']') { ++i; return v; }
			fail("expected ',' or ']' in array");
		}
	}

	std::string parse_string() {
		++i; // opening quote
		std::string out;
		while (i < s.size()) {
			const char c = s[i++];
			if (c == '"') return out;
			if (c == '\\') {
				if (i >= s.size()) fail("unterminated escape");
				const char e = s[i++];
				switch (e) {
					case '"':  out += '"';  break;
					case '\\': out += '\\'; break;
					case '/':  out += '/';  break;
					case 'n':  out += '\n'; break;
					case 't':  out += '\t'; break;
					case 'r':  out += '\r'; break;
					case 'b':  out += '\b'; break;
					case 'f':  out += '\f'; break;
					case 'u': {
						if (i + 4 > s.size()) fail("bad \\u escape");
						// Minimal: decode BMP code point to UTF-8.
						unsigned cp = 0;
						for (int k = 0; k < 4; ++k) {
							const char h = s[i++];
							cp <<= 4;
							if (h >= '0' && h <= '9') cp |= (h - '0');
							else if (h >= 'a' && h <= 'f') cp |= (h - 'a' + 10);
							else if (h >= 'A' && h <= 'F') cp |= (h - 'A' + 10);
							else fail("bad hex in \\u escape");
						}
						if (cp < 0x80) out += static_cast<char>(cp);
						else if (cp < 0x800) {
							out += static_cast<char>(0xC0 | (cp >> 6));
							out += static_cast<char>(0x80 | (cp & 0x3F));
						} else {
							out += static_cast<char>(0xE0 | (cp >> 12));
							out += static_cast<char>(0x80 | ((cp >> 6) & 0x3F));
							out += static_cast<char>(0x80 | (cp & 0x3F));
						}
						break;
					}
					default: fail("invalid escape character");
				}
			} else {
				out += c;
			}
		}
		fail("unterminated string");
	}

	Value parse_number() {
		const size_t start = i;
		if (peek() == '-' || peek() == '+') ++i;
		bool any = false;
		while (i < s.size() && (std::isdigit((unsigned char)s[i]) || s[i]=='.' ||
		       s[i]=='e' || s[i]=='E' || s[i]=='+' || s[i]=='-')) { ++i; any = true; }
		if (!any) fail("invalid value");
		const std::string tok = s.substr(start, i - start);
		try { return Value(std::stod(tok)); }
		catch (...) { fail("invalid number '" + tok + "'"); }
	}

	Value parse_bool() {
		if (s.compare(i, 4, "true") == 0)  { i += 4; return Value(true); }
		if (s.compare(i, 5, "false") == 0) { i += 5; return Value(false); }
		fail("invalid literal");
	}

	Value parse_null() {
		if (s.compare(i, 4, "null") == 0) { i += 4; return Value(); }
		fail("invalid literal");
	}
};

void dump_string(std::ostringstream& os, const std::string& s) {
	os << '"';
	for (char c : s) {
		switch (c) {
			case '"':  os << "\\\""; break;
			case '\\': os << "\\\\"; break;
			case '\n': os << "\\n";  break;
			case '\t': os << "\\t";  break;
			case '\r': os << "\\r";  break;
			default:   os << c;      break;
		}
	}
	os << '"';
}

void dump_number(std::ostringstream& os, double d) {
	if (std::floor(d) == d && std::fabs(d) < 1e15) {
		os << static_cast<long long>(d);          // integral -> no trailing ".0"
	} else {
		char buf[32];
		std::snprintf(buf, sizeof(buf), "%.10g", d);
		os << buf;
	}
}

void dump_value(std::ostringstream& os, const Value& v, int indent, int depth) {
	const std::string pad(indent > 0 ? (depth + 1) * indent : 0, ' ');
	const std::string pad0(indent > 0 ? depth * indent : 0, ' ');
	const char* nl = indent > 0 ? "\n" : "";
	const char* sp = indent > 0 ? " " : "";

	switch (v.type()) {
		case Value::Null:   os << "null"; break;
		case Value::Bool:   os << (v.as_bool() ? "true" : "false"); break;
		case Value::Number: dump_number(os, v.as_number()); break;
		case Value::String: dump_string(os, v.as_string()); break;
		case Value::Array: {
			if (v.size() == 0) { os << "[]"; break; }
			os << "[" << nl;
			const std::vector<Value>& a = v.elements();
			for (size_t k = 0; k < a.size(); ++k) {
				os << pad; dump_value(os, a[k], indent, depth + 1);
				if (k + 1 < a.size()) os << ",";
				os << nl;
			}
			os << pad0 << "]";
			break;
		}
		case Value::Object: {
			const std::map<std::string, Value>& o = v.items();
			if (o.empty()) { os << "{}"; break; }
			os << "{" << nl;
			size_t k = 0;
			for (std::map<std::string, Value>::const_iterator it = o.begin(); it != o.end(); ++it, ++k) {
				os << pad; dump_string(os, it->first); os << ":" << sp;
				dump_value(os, it->second, indent, depth + 1);
				if (k + 1 < o.size()) os << ",";
				os << nl;
			}
			os << pad0 << "}";
			break;
		}
	}
}

} // namespace

Value parse(const std::string& text) {
	Parser p(text);
	Value v = p.parse_value();
	p.skip_ws();
	if (p.i != text.size()) p.fail("trailing characters after value");
	return v;
}

Value parse_file(const std::string& path) {
	std::ifstream f(path.c_str());
	if (!f.is_open()) throw lsquant::IOError("cannot open " + path);
	std::stringstream ss; ss << f.rdbuf();
	return parse(ss.str());
}

std::string dump(const Value& v, int indent) {
	std::ostringstream os;
	dump_value(os, v, indent, 0);
	return os.str();
}

}} // namespace lsquant::json
