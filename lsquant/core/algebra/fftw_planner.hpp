#ifndef LSQUANT_CORE_ALGEBRA_FFTW_PLANNER_HPP_
#define LSQUANT_CORE_ALGEBRA_FFTW_PLANNER_HPP_
#include <mutex>
namespace lsquant {
namespace detail {
// FFTW planning/destruction is process-global, unlike execution of separate plans.
inline std::mutex& FftwPlannerMutex() {
  static std::mutex mutex;
  return mutex;
}
}  // namespace detail
}  // namespace lsquant
#endif
