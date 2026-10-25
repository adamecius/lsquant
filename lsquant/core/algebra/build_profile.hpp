/**
 * @file build_profile.hpp
 * @brief Compile-time build profile identifier ("golden" or "native").
 */
#ifndef LSQUANT_CORE_ALGEBRA_BUILD_PROFILE_HPP_
#define LSQUANT_CORE_ALGEBRA_BUILD_PROFILE_HPP_

namespace lsquant {

inline const char* BuildProfileName() {
#ifdef LSQUANT_NATIVE_BUILD
  return "native";
#else
  return "golden";
#endif
}

}  // namespace lsquant

#endif  // LSQUANT_CORE_ALGEBRA_BUILD_PROFILE_HPP_
