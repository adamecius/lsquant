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
// Item 4: std::vector value-initializes serially. Discard only complete interior
// pages of a newly owned buffer BEFORE filling it; the static owner then faults
// them locally. Never call on live data or shared/borrowed buffers. If madvise is
// unavailable, correctness is unchanged, placement is best effort. Tester must
// inspect numa_maps/numastat as well as timings; rewriting live pages is not a fix.
inline void DiscardInteriorPages(void* data, std::size_t bytes) {
#ifdef __linux__
  if (!data || bytes < 1024 * 1024) return;
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
#endif
