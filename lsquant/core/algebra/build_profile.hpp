#ifndef LSQUANT_CORE_ALGEBRA_BUILD_PROFILE_HPP_
#define LSQUANT_CORE_ALGEBRA_BUILD_PROFILE_HPP_

namespace lsquant {
// Item 3: provenance without changing the golden compiler command line.
// Tester: native headers/receipts must be tagged; untagged historical files are golden.
inline const char* BuildProfileName() {
#ifdef LSQUANT_NATIVE_BUILD
  return "native";
#else
  return "golden";
#endif
}
}  // namespace lsquant
#endif
