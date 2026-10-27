/**
 * @file chebyshev_vectors.cpp
 * @brief Implementation of ChebyshevVectors.
 */
#include "core/chebyshev/chebyshev_vectors.hpp"

#include <stdexcept>

#include "core/algebra/linalg.hpp"
#include "core/algebra/memory_placement.hpp"

namespace lsquant {

ChebyshevVectors::ChebyshevVectors(Index dimension, int batchSize)
    : dimension_(dimension),
      batchSize_(batchSize),
      firstOrder_(0),
#ifdef LSQUANT_NATIVE_BUILD
      vectors_(static_cast<std::size_t>(batchSize > 0 ? batchSize : 0)),
      scratch_(MakeZeroVector(dimension)) {
#else
      vectors_(static_cast<std::size_t>(batchSize > 0 ? batchSize : 0),
               Vector(static_cast<std::size_t>(dimension > 0 ? dimension : 0), Scalar(0.0, 0.0))),
      scratch_(static_cast<std::size_t>(dimension > 0 ? dimension : 0), Scalar(0.0, 0.0)) {
#endif
  if (dimension <= 0) {
    throw std::invalid_argument("ChebyshevVectors: dimension must be positive");
  }
  if (batchSize <= 0) {
    throw std::invalid_argument("ChebyshevVectors: batchSize must be positive");
  }
#ifdef LSQUANT_NATIVE_BUILD
  for (Vector& vector : vectors_) vector = MakeZeroVector(dimension);
#endif
}

void ChebyshevVectors::Fill(ChebyshevRecursion& recursion) {
  if (static_cast<Index>(recursion.Current().size()) != dimension_) {
    throw std::invalid_argument("ChebyshevVectors::Fill: recursion dimension mismatch");
  }
  firstOrder_ = recursion.Order();
  for (int m = 0; m < batchSize_; ++m) {
    linalg::copy(recursion.Current(), vectors_[static_cast<std::size_t>(m)]);
    recursion.Step();
  }
}

void ChebyshevVectors::ApplyToAll(const LinearOperator& A) {
  if (A.Dimension() != dimension_) {
    throw std::invalid_argument("ChebyshevVectors::ApplyToAll: operator dimension mismatch");
  }
  if (A.IsIdentity()) {
    return;
  }
#ifdef LSQUANT_NATIVE_BUILD
  if (batchScratch_.empty()) {
    for (int m = 0; m < batchSize_; ++m) batchScratch_.push_back(MakeZeroVector(dimension_));
  }
  A.MultiplyBatch(Scalar(1.0, 0.0), vectors_, Scalar(0.0, 0.0), batchScratch_);
  vectors_.swap(batchScratch_);
#else
  for (int m = 0; m < batchSize_; ++m) {
    Vector& vec = vectors_[static_cast<std::size_t>(m)];
    linalg::copy(vec, scratch_);
    A.Multiply(Scalar(1.0, 0.0), scratch_, Scalar(0.0, 0.0), vec);
  }
#endif
}

Vector& ChebyshevVectors::At(int m) {
  if (m < 0 || m >= batchSize_) {
    throw std::out_of_range("ChebyshevVectors::At: index out of range");
  }
  return vectors_[static_cast<std::size_t>(m)];
}

const Vector& ChebyshevVectors::At(int m) const {
  if (m < 0 || m >= batchSize_) {
    throw std::out_of_range("ChebyshevVectors::At: index out of range");
  }
  return vectors_[static_cast<std::size_t>(m)];
}

int ChebyshevVectors::Size() const { return batchSize_; }

int ChebyshevVectors::FirstOrder() const { return firstOrder_; }

}  // namespace lsquant
