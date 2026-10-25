/**
 * @file block_diagonal_operator.cpp
 * @brief Implementation of BlockDiagonalOperator.
 */
#include "core/algebra/block_diagonal_operator.hpp"

#include <algorithm>
#include <cmath>
#include <complex>
#include <limits>
#include <stdexcept>
#include <utility>
#include "core/algebra/memory_placement.hpp"
#ifdef LSQUANT_NATIVE_BUILD
#include <Eigen/Core>
#include "core/algebra/work_partition.hpp"
#endif

namespace lsquant {

BlockDiagonalOperator::BlockDiagonalOperator(Index numBlocks, Index blockSize,
                                             std::vector<Scalar> blocks)
    : numBlocks_(numBlocks), blockSize_(blockSize), blocks_(std::move(blocks)) {
  if (numBlocks < 0 || blockSize < 0) {
    throw std::invalid_argument(
        "BlockDiagonalOperator: numBlocks and blockSize must be non-negative");
  }
  const std::size_t expected = static_cast<std::size_t>(numBlocks) *
                               static_cast<std::size_t>(blockSize) *
                               static_cast<std::size_t>(blockSize);
  if (blocks_.size() != expected) {
    throw std::invalid_argument(
        "BlockDiagonalOperator: blocks size must equal numBlocks * blockSize^2");
  }

  const std::size_t w = static_cast<std::size_t>(blockSize_);
  const std::size_t w2 = w * w;
#ifdef LSQUANT_NATIVE_BUILD
  // First-touch contiguous block groups on the same static OpenMP partition as Multiply.
  std::vector<Scalar> placed(blocks_.size());
  DiscardInteriorPages(placed.data(), placed.size() * sizeof(Scalar));
#pragma omp parallel for schedule(static)
  for (Index ik = 0; ik < numBlocks_; ++ik) {
    const std::size_t begin = static_cast<std::size_t>(ik) * w2;
    std::copy_n(blocks_.data() + begin, w2, placed.data() + begin);
  }
  blocks_.swap(placed);
#endif
  double largestNorm = 0.0;
  double largestAntiNorm = 0.0;
  bool finite = true;
  for (Index ik = 0; ik < numBlocks_; ++ik) {
    const Scalar* blk = blocks_.data() + static_cast<std::size_t>(ik) * w2;
    double norm = 0.0;
    double antiNorm = 0.0;
    for (std::size_t i = 0; i < w; ++i) {
      for (std::size_t j = 0; j < w; ++j) {
        finite = finite && std::isfinite(blk[i * w + j].real()) &&
                 std::isfinite(blk[i * w + j].imag());
        norm = std::hypot(norm, std::abs(blk[i * w + j]));
        const Scalar diff = blk[i * w + j] - std::conj(blk[j * w + i]);
        antiNorm = std::hypot(antiNorm, 0.5 * std::abs(diff));
      }
    }
    largestNorm = std::max(largestNorm, norm);
    largestAntiNorm = std::max(largestAntiNorm, antiNorm);
  }
  // One common operator scale also accepts round-off in symmetry-zero blocks.
  isHermitian_ = finite && largestAntiNorm <=
      1e-12 * std::max(largestNorm, std::numeric_limits<double>::min());
}

Index BlockDiagonalOperator::Dimension() const { return numBlocks_ * blockSize_; }

void BlockDiagonalOperator::Multiply(Scalar a, const Vector& x, Scalar b, Vector& y) const {
  const Index dim = Dimension();
  if (static_cast<Index>(x.size()) != dim || static_cast<Index>(y.size()) != dim) {
    throw std::invalid_argument("BlockDiagonalOperator::Multiply: vector dimension mismatch");
  }
  if (dim == 0) {
    return;
  }

#ifdef LSQUANT_NATIVE_BUILD
  // A single vector is the B=1 case; materialize the W-entry result before
  // overwriting so aliases retain their documented semantics.
  using Dense = Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
  using DenseVector = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>;
#pragma omp parallel for schedule(static)
  for (Index ik = 0; ik < numBlocks_; ++ik) {
    const Index offset = ik * blockSize_;
    Eigen::Map<const Dense> block(blocks_.data() + ik * blockSize_ * blockSize_, blockSize_, blockSize_);
    Eigen::Map<const DenseVector> input(x.data() + offset, blockSize_);
    const DenseVector product = block * input;
    for (Index row = 0; row < blockSize_; ++row) y[offset + row] = a * product(row) + b * y[offset + row];
  }
#else
  const Index nk = numBlocks_;
  const std::size_t w = static_cast<std::size_t>(blockSize_);
  const std::size_t w2 = w * w;
  const Scalar* blkBase = blocks_.data();
  const Scalar* xPtr = x.data();
  Scalar* yPtr = y.data();
  const bool aliased = (xPtr == yPtr);

  if (aliased) {
    #pragma omp parallel for schedule(static)
    for (Index ik = 0; ik < nk; ++ik) {
      const std::size_t offset = static_cast<std::size_t>(ik) * w;
      const Scalar* blk = blkBase + static_cast<std::size_t>(ik) * w2;
      Vector xLocal(xPtr + offset, xPtr + offset + w);
      Scalar* yk = yPtr + offset;
      for (std::size_t alpha = 0; alpha < w; ++alpha) {
        const Scalar* row = blk + alpha * w;
        Scalar acc(0.0, 0.0);
        for (std::size_t beta = 0; beta < w; ++beta) {
          acc += row[beta] * xLocal[beta];
        }
        yk[alpha] = a * acc + b * yk[alpha];
      }
    }
  } else {
    #pragma omp parallel for schedule(static)
    for (Index ik = 0; ik < nk; ++ik) {
      const std::size_t offset = static_cast<std::size_t>(ik) * w;
      const Scalar* blk = blkBase + static_cast<std::size_t>(ik) * w2;
      const Scalar* xk = xPtr + offset;
      Scalar* yk = yPtr + offset;
      for (std::size_t alpha = 0; alpha < w; ++alpha) {
        const Scalar* row = blk + alpha * w;
        Scalar acc(0.0, 0.0);
        for (std::size_t beta = 0; beta < w; ++beta) {
          acc += row[beta] * xk[beta];
        }
        yk[alpha] = a * acc + b * yk[alpha];
      }
    }
  }
#endif
}

BlockDiagonalOperator BlockDiagonalOperator::Adjoint() const {
  const std::size_t nk = static_cast<std::size_t>(numBlocks_);
  const std::size_t w = static_cast<std::size_t>(blockSize_);
  const std::size_t w2 = w * w;
  std::vector<Scalar> adj(blocks_.size(), Scalar(0.0, 0.0));
  for (std::size_t ik = 0; ik < nk; ++ik) {
    const Scalar* src = blocks_.data() + ik * w2;
    Scalar* dst = adj.data() + ik * w2;
    for (std::size_t i = 0; i < w; ++i) {
      for (std::size_t j = 0; j < w; ++j) {
        dst[i * w + j] = std::conj(src[j * w + i]);
      }
    }
  }
  return BlockDiagonalOperator(numBlocks_, blockSize_, std::move(adj));
}

void BlockDiagonalOperator::MultiplyBatch(Scalar a, const std::vector<Vector>& x, Scalar b,
                                         std::vector<Vector>& y) const {
  if (x.size() != y.size() || &x == &y) {
    throw std::invalid_argument("BlockDiagonalOperator::MultiplyBatch: invalid batch");
  }
  for (std::size_t j = 0; j < x.size(); ++j) {
    if (static_cast<Index>(x[j].size()) != Dimension() || x[j].size() != y[j].size()) {
      throw std::invalid_argument(
          "BlockDiagonalOperator::MultiplyBatch: expected vector length " +
          std::to_string(Dimension()) + ", got x=" + std::to_string(x[j].size()) +
          " and y=" + std::to_string(y[j].size()));
    }
  }
#ifdef LSQUANT_NATIVE_BUILD
  if (x.empty() || blockSize_ == 0) return;
  using Dense = Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
  // Fixed microbatch geometry: one GEMM per block and microbatch, with no team-dependent
  // reduction tree. Bounds per-core working set by half of a 512 KiB L2 cache.
  constexpr Index kHalfL2Scalars = (512 * 1024) / (2 * sizeof(Scalar));
  const int candidate = static_cast<int>(std::max<Index>(1, std::min<Index>(
      static_cast<Index>(x.size()), (kHalfL2Scalars - blockSize_ * blockSize_) / (2 * blockSize_))));
  const int width = KernelBatchSetting() > 0 ? KernelBatchOverride(static_cast<int>(x.size())) : candidate;
#pragma omp parallel
  {
    Dense input(blockSize_, width), output(blockSize_, width);
#pragma omp for schedule(static)
    for (Index ik = 0; ik < numBlocks_; ++ik) {
      const Index offset = ik * blockSize_;
      Eigen::Map<const Dense> block(blocks_.data() + ik * blockSize_ * blockSize_, blockSize_, blockSize_);
      for (std::size_t base = 0; base < x.size(); base += static_cast<std::size_t>(width)) {
        const int count = static_cast<int>(std::min<std::size_t>(width, x.size() - base));
        for (Index row = 0; row < blockSize_; ++row) {
          for (int j = 0; j < count; ++j) input(row, j) = x[base + j][offset + row];
        }
        output.leftCols(count).noalias() = block * input.leftCols(count);
        for (Index row = 0; row < blockSize_; ++row) {
          for (int j = 0; j < count; ++j) {
            y[base + j][offset + row] = a * output(row, j) + b * y[base + j][offset + row];
          }
        }
      }
    }
  }
#else
  LinearOperator::MultiplyBatch(a, x, b, y);
#endif
}

}  // namespace lsquant
