/**
 * @file metric_policy.cpp
 * @brief Implementation of OrthogonalMetric and NonOrthogonalMetric.
 */
#include "core/algebra/metric_policy.hpp"

#include <stdexcept>

namespace lsquant {

OrthogonalMetric::OrthogonalMetric(Index dimension) : identity_(dimension) {}

const LinearOperator& OrthogonalMetric::Overlap() const { return identity_; }

const LinearOperator& OrthogonalMetric::InverseOverlap() const { return identity_; }

const LinearOperator& OrthogonalMetric::InverseSqrtOverlap() const { return identity_; }

NonOrthogonalMetric::NonOrthogonalMetric(const SparseMatrix& S, const LinearOperator& inverse,
                                         const LinearOperator& inverseSqrt)
    : overlap_(S), inverse_(&inverse), inverseSqrt_(&inverseSqrt) {
  const Index dim = S.Dimension();
  if (dim <= 0) {
    throw std::invalid_argument("NonOrthogonalMetric: overlap dimension must be positive");
  }
  if (!overlap_.IsHermitian()) {
    throw std::invalid_argument("NonOrthogonalMetric: overlap matrix S must be Hermitian");
  }
  if (inverse.Dimension() != dim || inverseSqrt.Dimension() != dim) {
    throw std::invalid_argument("NonOrthogonalMetric: applier dimension mismatch with S");
  }
}

const LinearOperator& NonOrthogonalMetric::Overlap() const { return overlap_; }

const LinearOperator& NonOrthogonalMetric::InverseOverlap() const { return *inverse_; }

const LinearOperator& NonOrthogonalMetric::InverseSqrtOverlap() const { return *inverseSqrt_; }

}  // namespace lsquant
