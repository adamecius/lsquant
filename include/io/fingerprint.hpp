#ifndef LSQUANT_IO_FINGERPRINT_HPP
#define LSQUANT_IO_FINGERPRINT_HPP

// System fingerprint (unified_io.md, Phase I3).
//
// A stable SHA-256 over exactly the inputs that define the physical system: the
// Hamiltonian + operator CSR files and their descriptors. Two runs may be merged
// iff their fingerprints match. The hash deliberately EXCLUDES num_moments,
// broadening, and which observable -- those differ between mergeable runs.
//
// In-house SHA-256 (no crypto-strength requirement; identity + human-quotable hex)
// keeps the project dependency-free and offline-safe.

#include <string>
#include <vector>

namespace lsquant { namespace io {

// SHA-256 of an arbitrary byte string, lowercase hex (64 chars).
std::string sha256_hex(const std::string& data);

// Fingerprint a system given its label and operators directory. Hashes, in a
// fixed order, the bytes of every "<dir>/<label>.<ID>.CSR" present plus any
// "<dir>/<label>.<ID>.desc". Missing files are skipped (recorded as absent) so
// the hash is well-defined whatever subset exists. Returns 64-char hex.
std::string fingerprint_system(const std::string& label, const std::string& operators_dir);

// Lower-level: fingerprint an explicit ordered list of files (each contributes
// its name + length + content). Files that cannot be read contribute a marker.
std::string fingerprint_files(const std::vector<std::string>& paths);

}} // namespace lsquant::io

#endif // LSQUANT_IO_FINGERPRINT_HPP
