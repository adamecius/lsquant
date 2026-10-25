/**
 * @file composite_operators.cpp
 * @brief Implementation of algebraic composite operators.
 */
#include "core/algebra/composite_operators.hpp"

#include <cmath>
#include <complex>
#include <stdexcept>
#include <utility>

#include "core/algebra/bloch_transform.hpp"
#include "core/algebra/block_diagonal_operator.hpp"

namespace lsquant {

IdentityOperator::IdentityOperator(Index dimension) : dimension_(dimension) {
  if (dimension < 0) {
    throw std::invalid_argument("IdentityOperator: negative dimension");
  }
}

Index IdentityOperator::Dimension() const { return dimension_; }

void IdentityOperator::Multiply(Scalar a, const Vector& x, Scalar b, Vector& y) const {
  if (static_cast<Index>(x.size()) != dimension_ ||
      static_cast<Index>(y.size()) != dimension_) {
    throw std::invalid_argument("IdentityOperator::Multiply: vector dimension mismatch");
  }
  for (Index i = 0; i < dimension_; ++i) {
    const std::size_t idx = static_cast<std::size_t>(i);
    y[idx] = a * x[idx] + b * y[idx];
  }
}

DiagonalOperator::DiagonalOperator(std::vector<Scalar> values)
    : values_(std::move(values)) {
  bool herm = true;
  for (const Scalar& v : values_) {
    if (std::abs(v.imag()) > 1e-12) {
      herm = false;
      break;
    }
  }
  isHermitian_ = herm;
}

Index DiagonalOperator::Dimension() const {
  return static_cast<Index>(values_.size());
}

void DiagonalOperator::Multiply(Scalar a, const Vector& x, Scalar b, Vector& y) const {
  const Index dim = Dimension();
  if (static_cast<Index>(x.size()) != dim || static_cast<Index>(y.size()) != dim) {
    throw std::invalid_argument("DiagonalOperator::Multiply: vector dimension mismatch");
  }
  const Scalar* dPtr = values_.data();
  const Scalar* xPtr = x.data();
  Scalar* yPtr = y.data();
  #pragma omp parallel for schedule(static)
  for (Index i = 0; i < dim; ++i) {
    const Scalar xi = xPtr[i];
    yPtr[i] = a * dPtr[i] * xi + b * yPtr[i];
  }
}

DiagonalOperator DiagonalOperator::Adjoint() const {
  std::vector<Scalar> conjVals(values_.size());
  for (std::size_t i = 0; i < values_.size(); ++i) {
    conjVals[i] = std::conj(values_[i]);
  }
  return DiagonalOperator(std::move(conjVals));
}

ScaledOperator::ScaledOperator(Scalar alpha, const LinearOperator& op)
    : alpha_(alpha), op_(&op) {}

Index ScaledOperator::Dimension() const { return op_->Dimension(); }

void ScaledOperator::Multiply(Scalar a, const Vector& x, Scalar b, Vector& y) const {
  op_->Multiply(a * alpha_, x, b, y);
}

bool ScaledOperator::IsHermitian() const {
  return op_->IsHermitian() && (alpha_.imag() == 0.0);
}

OperatorSum::OperatorSum(const LinearOperator& first, const LinearOperator& second)
    : first_(&first), second_(&second) {
  if (first.Dimension() != second.Dimension()) {
    throw std::invalid_argument("OperatorSum: operator dimensions must match");
  }
}

Index OperatorSum::Dimension() const { return first_->Dimension(); }

void OperatorSum::Multiply(Scalar a, const Vector& x, Scalar b, Vector& y) const {
  const Index dim = Dimension();
  if (static_cast<Index>(x.size()) != dim || static_cast<Index>(y.size()) != dim) {
    throw std::invalid_argument("OperatorSum::Multiply: vector dimension mismatch");
  }
  if (&x == &y) {
    Vector xCopy = x;
    first_->Multiply(a, xCopy, b, y);
    second_->Multiply(a, xCopy, Scalar(1.0, 0.0), y);
    return;
  }
  first_->Multiply(a, x, b, y);
  second_->Multiply(a, x, Scalar(1.0, 0.0), y);
}

bool OperatorSum::IsHermitian() const {
  return first_->IsHermitian() && second_->IsHermitian();
}

OperatorProduct::OperatorProduct(const LinearOperator& first, const LinearOperator& second)
    : first_(&first), second_(&second) {
  if (first.Dimension() != second.Dimension()) {
    throw std::invalid_argument("OperatorProduct: operator dimensions must match");
  }
}

Index OperatorProduct::Dimension() const { return first_->Dimension(); }

void OperatorProduct::Multiply(Scalar a, const Vector& x, Scalar b, Vector& y) const {
  const Index dim = Dimension();
  if (static_cast<Index>(x.size()) != dim || static_cast<Index>(y.size()) != dim) {
    throw std::invalid_argument("OperatorProduct::Multiply: vector dimension mismatch");
  }
  Vector temp(static_cast<std::size_t>(dim), Scalar(0.0, 0.0));
  second_->Multiply(Scalar(1.0, 0.0), x, Scalar(0.0, 0.0), temp);
  first_->Multiply(a, temp, b, y);
}

bool OperatorProduct::IsHermitian() const {
  if (first_->IsIdentity()) {
    return second_->IsHermitian();
  }
  if (second_->IsIdentity()) {
    return first_->IsHermitian();
  }
  return false;
}

AdjointOperator::AdjointOperator(const LinearOperator& op) : isHermitian_(op.IsHermitian()) {
  if (const auto* sparseOp = dynamic_cast<const SparseMatrixOperator*>(&op)) {
    adjointMatrix_ = sparseOp->Matrix().Adjoint();
  } else if (const auto* blochOp = dynamic_cast<const BlochTransform*>(&op)) {
    ownedAdjoint_ = std::make_shared<BlochTransform>(blochOp->Adjoint());
  } else if (const auto* blockOp = dynamic_cast<const BlockDiagonalOperator*>(&op)) {
    if (blockOp->IsHermitian()) {
      forwardedOp_ = &op;
    } else {
      ownedAdjoint_ = std::make_shared<BlockDiagonalOperator>(blockOp->Adjoint());
    }
  } else if (const auto* diagOp = dynamic_cast<const DiagonalOperator*>(&op)) {
    if (diagOp->IsHermitian()) {
      forwardedOp_ = &op;
    } else {
      ownedAdjoint_ = std::make_shared<DiagonalOperator>(diagOp->Adjoint());
    }
  } else if (op.IsHermitian()) {
    forwardedOp_ = &op;
  } else {
    throw std::invalid_argument(
        "AdjointOperator: unsupported non-Hermitian LinearOperator type");
  }
}

AdjointOperator::AdjointOperator(const SparseMatrix& matrix)
    : adjointMatrix_(matrix.Adjoint()), isHermitian_(matrix.IsHermitian()) {}

Index AdjointOperator::Dimension() const {
  if (adjointMatrix_.has_value()) {
    return adjointMatrix_->Dimension();
  }
  if (ownedAdjoint_ != nullptr) {
    return ownedAdjoint_->Dimension();
  }
  return forwardedOp_->Dimension();
}

void AdjointOperator::Multiply(Scalar a, const Vector& x, Scalar b, Vector& y) const {
  if (adjointMatrix_.has_value()) {
    adjointMatrix_->Multiply(a, x, b, y);
  } else if (ownedAdjoint_ != nullptr) {
    ownedAdjoint_->Multiply(a, x, b, y);
  } else {
    forwardedOp_->Multiply(a, x, b, y);
  }
}

bool AdjointOperator::IsHermitian() const { return isHermitian_; }

SandwichOperator::SandwichOperator(const LinearOperator& B, const LinearOperator& V)
    : b_(&B),
      v_(&V),
      bAdjoint_(B),
      t1_(static_cast<std::size_t>(B.Dimension()), Scalar(0.0, 0.0)),
      t2_(static_cast<std::size_t>(B.Dimension()), Scalar(0.0, 0.0)) {
  if (B.Dimension() != V.Dimension()) {
    throw std::invalid_argument("SandwichOperator: B and V dimensions must match");
  }
}

Index SandwichOperator::Dimension() const { return b_->Dimension(); }

void SandwichOperator::Multiply(Scalar a, const Vector& x, Scalar b, Vector& y) const {
  const Index dim = Dimension();
  if (static_cast<Index>(x.size()) != dim || static_cast<Index>(y.size()) != dim) {
    throw std::invalid_argument("SandwichOperator::Multiply: vector dimension mismatch");
  }
  bAdjoint_.Multiply(Scalar(1.0, 0.0), x, Scalar(0.0, 0.0), t1_);
  v_->Multiply(Scalar(1.0, 0.0), t1_, Scalar(0.0, 0.0), t2_);
  b_->Multiply(a, t2_, b, y);
}

bool SandwichOperator::IsHermitian() const { return v_->IsHermitian(); }

LowdinOperator::LowdinOperator(const LinearOperator& inverseSqrt, const LinearOperator& O)
    : inverseSqrt_(&inverseSqrt), op_(&O) {
  if (inverseSqrt.Dimension() != O.Dimension()) {
    throw std::invalid_argument("LowdinOperator: operator dimensions must match");
  }
}

Index LowdinOperator::Dimension() const { return op_->Dimension(); }

void LowdinOperator::Multiply(Scalar a, const Vector& x, Scalar b, Vector& y) const {
  const Index dim = Dimension();
  if (static_cast<Index>(x.size()) != dim || static_cast<Index>(y.size()) != dim) {
    throw std::invalid_argument("LowdinOperator::Multiply: vector dimension mismatch");
  }
  Vector t1(static_cast<std::size_t>(dim), Scalar(0.0, 0.0));
  Vector t2(static_cast<std::size_t>(dim), Scalar(0.0, 0.0));
  inverseSqrt_->Multiply(Scalar(1.0, 0.0), x, Scalar(0.0, 0.0), t1);
  op_->Multiply(Scalar(1.0, 0.0), t1, Scalar(0.0, 0.0), t2);
  inverseSqrt_->Multiply(a, t2, b, y);
}

bool LowdinOperator::IsHermitian() const {
  return inverseSqrt_->IsHermitian() && op_->IsHermitian();
}

MetricInverseProduct::MetricInverseProduct(const LinearOperator& inverse,
                                           const LinearOperator& O)
    : inverse_(&inverse), op_(&O) {
  if (inverse.Dimension() != O.Dimension()) {
    throw std::invalid_argument("MetricInverseProduct: operator dimensions must match");
  }
}

Index MetricInverseProduct::Dimension() const { return op_->Dimension(); }

void MetricInverseProduct::Multiply(Scalar a, const Vector& x, Scalar b, Vector& y) const {
  const Index dim = Dimension();
  if (static_cast<Index>(x.size()) != dim || static_cast<Index>(y.size()) != dim) {
    throw std::invalid_argument("MetricInverseProduct::Multiply: vector dimension mismatch");
  }
  Vector temp(static_cast<std::size_t>(dim), Scalar(0.0, 0.0));
  op_->Multiply(Scalar(1.0, 0.0), x, Scalar(0.0, 0.0), temp);
  inverse_->Multiply(a, temp, b, y);
}

HermitianTripleProduct::HermitianTripleProduct(
    const LinearOperator& first, const LinearOperator& middle,
    const LinearOperator& last, Scalar coefficient)
    : first_(&first), middle_(&middle), last_(&last), coefficient_(coefficient) {
  if (first.Dimension() != middle.Dimension() || first.Dimension() != last.Dimension()) {
    throw std::invalid_argument("HermitianTripleProduct: operator dimensions must match");
  }
  if (!first.IsHermitian() || !middle.IsHermitian() || !last.IsHermitian()) {
    throw std::invalid_argument("HermitianTripleProduct: every factor must be Hermitian");
  }
  if (!std::isfinite(coefficient.real()) || !std::isfinite(coefficient.imag())) {
    throw std::invalid_argument("HermitianTripleProduct: coefficient must be finite");
  }
}

Index HermitianTripleProduct::Dimension() const { return first_->Dimension(); }

void HermitianTripleProduct::Multiply(Scalar a, const Vector& x, Scalar b, Vector& y) const {
  const Index dim = Dimension();
  if (static_cast<Index>(x.size()) != dim || static_cast<Index>(y.size()) != dim) {
    throw std::invalid_argument("HermitianTripleProduct::Multiply: vector dimension mismatch");
  }
  Vector alias;
  const Vector* input = &x;
  if (&x == &y) {
    alias = x;
    input = &alias;
  }
  Vector t1(static_cast<std::size_t>(dim)), t2(static_cast<std::size_t>(dim));
  last_->Multiply(Scalar(1.0, 0.0), *input, Scalar(0.0, 0.0), t1);
  middle_->Multiply(Scalar(1.0, 0.0), t1, Scalar(0.0, 0.0), t2);
  first_->Multiply(a * coefficient_, t2, b, y);
  first_->Multiply(Scalar(1.0, 0.0), *input, Scalar(0.0, 0.0), t1);
  middle_->Multiply(Scalar(1.0, 0.0), t1, Scalar(0.0, 0.0), t2);
  last_->Multiply(a * std::conj(coefficient_), t2, Scalar(1.0, 0.0), y);
}

}  // namespace lsquant
