#ifndef LSQUANT_CORE_ALGEBRA_WORK_PARTITION_HPP_
#define LSQUANT_CORE_ALGEBRA_WORK_PARTITION_HPP_

#include <algorithm>
#include <omp.h>
#include "core/algebra/types.hpp"

namespace lsquant {
inline int& KernelBatchSetting() {
  static int width = 0;
  return width;
}
inline Index& SpmvThresholdSetting() {
  static Index threshold = 16384;
  return threshold;
}
inline int KernelBatchOverride(int maximum) {
  return KernelBatchSetting() > 0 ? std::min(maximum, KernelBatchSetting()) : maximum;
}
// Item 5a: candidate crossover, NOT a measured Aron constant. Tester sweeps
// LSQUANT_SPMV_MIN_NNZ=0,4096,16384,65536,262144 before blessing the default.
inline Index SpmvParallelThreshold() {
  return SpmvThresholdSetting();
}

// Fixed contiguous row ranges, approximately equal nnz; no row is split and
// each row's reduction order is unchanged. Empty rows still belong to one owner.
inline Index RowBoundary(const Index* rowPointers, Index rows, int part, int parts) {
  if (part == 0) return 0;
  if (part == parts) return rows;
  const Index nnz = rowPointers[rows];
  if (nnz == 0) return rows * part / parts;
  const Index target = (nnz / parts) * part + (nnz % parts) * part / parts;
  return static_cast<Index>(std::lower_bound(rowPointers, rowPointers + rows, target) - rowPointers);
}

// Item 5e: input/output plus three recursion streams per core; half of 512 KiB
// L2 reserved for streams. A global batch is a storage window; this microbatch
// bounds the hot operator working set without restarting H. For rows/core q,
// B=max(1,floor((L2/(2*sizeof(Scalar))-3q)/(2q))). Tester sweeps B and measures.
inline int CacheBatchSize(Index rows, int maximum) {
  const Index q = std::max<Index>(1, (rows + omp_get_max_threads() - 1) / omp_get_max_threads());
  const Index available = Index(512 * 1024 / (2 * sizeof(Scalar))) - 3 * q;
  const int candidate = static_cast<int>(std::max<Index>(1, std::min<Index>(maximum, available / (2 * q))));
  return KernelBatchSetting() > 0 ? KernelBatchOverride(maximum) : candidate;
}
}  // namespace lsquant
#endif
