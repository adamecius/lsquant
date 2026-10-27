/**
 * @file work_partition.hpp
 * @brief Row partitioning and cache microbatch sizing for OpenMP sparse kernels.
 */
#ifndef LSQUANT_CORE_ALGEBRA_WORK_PARTITION_HPP_
#define LSQUANT_CORE_ALGEBRA_WORK_PARTITION_HPP_

#include <algorithm>
#include <omp.h>
#include "core/algebra/types.hpp"

namespace lsquant {

constexpr Index kDefaultSpmvThreshold = 16384;
constexpr Index kHalfL2CacheScalars = (512 * 1024) / (2 * sizeof(Scalar));

inline int& KernelBatchSetting() {
  static int width = 0;
  return width;
}

inline Index& SpmvThresholdSetting() {
  static Index threshold = kDefaultSpmvThreshold;
  return threshold;
}

inline int KernelBatchOverride(int maximum) {
  return KernelBatchSetting() > 0 ? std::min(maximum, KernelBatchSetting()) : maximum;
}

inline Index SpmvParallelThreshold() {
  return SpmvThresholdSetting();
}

// Fixed contiguous row ranges with approximately equal nnz; no row is split and
// each row's reduction order is unchanged.
inline Index RowBoundary(const Index* rowPointers, Index rows, int part, int parts) {
  if (part == 0) return 0;
  if (part == parts) return rows;
  const Index nnz = rowPointers[rows];
  if (nnz == 0) return rows * part / parts;
  const Index target = (nnz / parts) * part + (nnz % parts) * part / parts;
  return static_cast<Index>(std::lower_bound(rowPointers, rowPointers + rows, target) - rowPointers);
}

// Bound the per-core operator working set (input/output plus three recursion streams)
// within half of a 512 KiB L2 cache without restarting the Chebyshev recursion.
inline int CacheBatchSize(Index rows, int maximum) {
  const Index q = std::max<Index>(1, (rows + omp_get_max_threads() - 1) / omp_get_max_threads());
  const Index available = kHalfL2CacheScalars - 3 * q;
  const int candidate = static_cast<int>(std::max<Index>(1, std::min<Index>(maximum, available / (2 * q))));
  return KernelBatchSetting() > 0 ? KernelBatchOverride(maximum) : candidate;
}

}  // namespace lsquant

#endif  // LSQUANT_CORE_ALGEBRA_WORK_PARTITION_HPP_
