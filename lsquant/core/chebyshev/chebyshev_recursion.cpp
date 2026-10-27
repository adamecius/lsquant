/**
 * @file chebyshev_recursion.cpp
 * @brief Implementation of ChebyshevRecursion.
 */
#include "core/chebyshev/chebyshev_recursion.hpp"

#include <stdexcept>

#include "core/algebra/linalg.hpp"
#include "core/algebra/memory_placement.hpp"

namespace lsquant {

ChebyshevRecursion::ChebyshevRecursion(const LinearOperator& rescaledH)
    : rescaledH_(&rescaledH),
#ifdef LSQUANT_NATIVE_BUILD
      current_(MakeZeroVector(rescaledH.Dimension())),
      previous_(MakeZeroVector(rescaledH.Dimension())),
      next_(MakeZeroVector(rescaledH.Dimension())),
#else
      current_(static_cast<std::size_t>(rescaledH.Dimension() > 0 ? rescaledH.Dimension() : 0),
               Scalar(0.0, 0.0)),
      previous_(static_cast<std::size_t>(rescaledH.Dimension() > 0 ? rescaledH.Dimension() : 0),
                Scalar(0.0, 0.0)),
      next_(static_cast<std::size_t>(rescaledH.Dimension() > 0 ? rescaledH.Dimension() : 0),
            Scalar(0.0, 0.0)),
#endif
      order_(-1) {
  if (rescaledH.Dimension() <= 0) {
    throw std::invalid_argument("ChebyshevRecursion: operator dimension must be positive");
  }
}

void ChebyshevRecursion::Start(const Vector& phi) {
  const Index dim = rescaledH_->Dimension();
  if (static_cast<Index>(phi.size()) != dim) {
    throw std::invalid_argument("ChebyshevRecursion::Start: phi dimension mismatch");
  }
  linalg::copy(phi, current_);
  linalg::copy(phi, previous_);
  rescaledH_->Multiply(Scalar(1.0, 0.0), current_, Scalar(0.0, 0.0), next_);
  order_ = 0;
}

void ChebyshevRecursion::Step() {
  if (order_ < 0) {
    throw std::runtime_error("ChebyshevRecursion::Step: Start must be called before Step");
  }
  if (order_ == 0) {
    linalg::copy(current_, previous_);
    current_.swap(next_);
    order_ = 1;
    return;
  }
  rescaledH_->Multiply(Scalar(2.0, 0.0), current_, Scalar(-1.0, 0.0), previous_);
  previous_.swap(current_);
  ++order_;
}

const Vector& ChebyshevRecursion::Current() const {
  if (order_ < 0) {
    throw std::runtime_error("ChebyshevRecursion::Current: Start must be called first");
  }
  return current_;
}

const Vector& ChebyshevRecursion::Previous() const {
  if (order_ < 0) {
    throw std::runtime_error("ChebyshevRecursion::Previous: Start must be called first");
  }
  return previous_;
}

int ChebyshevRecursion::Order() const { return order_; }

}  // namespace lsquant
