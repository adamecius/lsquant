/**
 * @file memory_placement.hpp
 * @brief NUMA-aware page placement helpers for native builds.
 */
#ifndef LSQUANT_CORE_ALGEBRA_MEMORY_PLACEMENT_HPP_
#define LSQUANT_CORE_ALGEBRA_MEMORY_PLACEMENT_HPP_

#include "core/algebra/types.hpp"
#ifdef LSQUANT_NATIVE_BUILD
#include <cstddef>
#include <cstdint>
#ifdef __linux__
#include <sys/mman.h>
#include <unistd.h>
#endif

namespace lsquant {

constexpr std::size_t kMinPageDiscardBytes = 1024 * 1024;

// Discard complete interior pages of a newly allocated buffer before parallel first-touch
// initialization so each OpenMP thread faults its slice on its local NUMA node.
inline void DiscardInteriorPages(void* data, std::size_t bytes) {
#ifdef __linux__
  if (!data || bytes < kMinPageDiscardBytes) return;
  const long pageQuery = sysconf(_SC_PAGESIZE);
  if (pageQuery <= 0) return;
  const auto page = static_cast<std::uintptr_t>(pageQuery);
  const auto begin = reinterpret_cast<std::uintptr_t>(data);
  const auto first = ((begin + page - 1) / page) * page;
  const auto last = ((begin + bytes) / page) * page;
  if (last > first) {
    (void)madvise(reinterpret_cast<void*>(first), last - first, MADV_DONTNEED);
  }
#else
  (void)data;
  (void)bytes;
#endif
}

inline Vector MakeZeroVector(Index size) {
  Vector result(static_cast<std::size_t>(size > 0 ? size : 0));
  DiscardInteriorPages(result.data(), result.size() * sizeof(Scalar));
#pragma omp parallel for schedule(static)
  for (Index i = 0; i < size; ++i) result[static_cast<std::size_t>(i)] = Scalar(0.0, 0.0);
  return result;
}

}  // namespace lsquant
#endif
#endif  // LSQUANT_CORE_ALGEBRA_MEMORY_PLACEMENT_HPP_
