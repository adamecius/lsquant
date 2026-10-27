/**
 * @file overlap_appliers.cpp
 * @brief Implementation of ExactInverseOverlap, ExactInverseSqrtOverlap, IschoInverseOverlap,
 *        ChebyshevInverseOverlap, and ChebyshevInverseSqrtOverlap.
 */
#include "core/algebra/overlap_appliers.hpp"

#include <Eigen/Eigenvalues>
#include <algorithm>
#include <cmath>
#include <sstream>
#include <stdexcept>

namespace lsquant {

namespace {

constexpr double kPi = 3.141592653589793238462643383279502884;

template <typename Func>
void ComputeChebyshevGaussCoeffs(SpectralBounds boundsOfS, int order, double coeffTol, Func func,
                                 double& centerOut, double& halfwidthOut,
                                 std::vector<double>& coeffsOut, int& effectiveOrderOut) {
  if (!(boundsOfS.lower > 0.0) || !(boundsOfS.upper > boundsOfS.lower) ||
      !std::isfinite(boundsOfS.lower) || !std::isfinite(boundsOfS.upper)) {
    throw std::invalid_argument(
        "ComputeChebyshevGaussCoeffs: boundsOfS must satisfy 0 < lower < upper");
  }
  if (order < 1) {
    throw std::invalid_argument("ComputeChebyshevGaussCoeffs: order must be >= 1");
  }
  if (coeffTol < 0.0 || !std::isfinite(coeffTol)) {
    throw std::invalid_argument("ComputeChebyshevGaussCoeffs: coeffTol must be non-negative");
  }

  const double a = boundsOfS.lower;
  const double b = boundsOfS.upper;
  const double center = 0.5 * (a + b);
  const double halfwidth = 0.5 * (b - a);

  const int numNodes = order + 1;
  std::vector<double> coeffsFull(static_cast<std::size_t>(numNodes), 0.0);
  std::vector<double> fvals(static_cast<std::size_t>(numNodes), 0.0);
  std::vector<double> thetas(static_cast<std::size_t>(numNodes), 0.0);

  for (int j = 0; j < numNodes; ++j) {
    const double theta = kPi * (2.0 * static_cast<double>(j) + 1.0) /
                         (2.0 * static_cast<double>(numNodes));
    const double xj = center + halfwidth * std::cos(theta);
    thetas[static_cast<std::size_t>(j)] = theta;
    fvals[static_cast<std::size_t>(j)] = func(xj);
  }

  for (int m = 0; m < numNodes; ++m) {
    double sum = 0.0;
    for (int j = 0; j < numNodes; ++j) {
      sum += fvals[static_cast<std::size_t>(j)] *
             std::cos(static_cast<double>(m) * thetas[static_cast<std::size_t>(j)]);
    }
    coeffsFull[static_cast<std::size_t>(m)] = (2.0 / static_cast<double>(numNodes)) * sum;
  }

  const double c0mag = std::abs(coeffsFull[0]);
  int effOrder = order;
  for (int m = order; m >= 1; --m) {
    if (std::abs(coeffsFull[static_cast<std::size_t>(m)]) > coeffTol * c0mag) {
      effOrder = m;
      break;
    }
    effOrder = 0;
  }
  effOrder = std::max(effOrder, 1);

  centerOut = center;
  halfwidthOut = halfwidth;
  coeffsOut.assign(coeffsFull.begin(), coeffsFull.begin() + effOrder + 1);
  effectiveOrderOut = effOrder;
}

void ApplyChebyshevSeries(const SparseMatrix& S, double center, double halfwidth,
                          const std::vector<double>& coeffs, Scalar a, const Vector& x, Scalar b,
                          Vector& y) {
  const Index dim = S.Dimension();
  if (static_cast<Index>(x.size()) != dim || static_cast<Index>(y.size()) != dim) {
    throw std::invalid_argument("ApplyChebyshevSeries: vector dimension mismatch");
  }
  if (&x == &y) {
    throw std::invalid_argument("ApplyChebyshevSeries: input and output vectors must not alias");
  }

  const std::size_t n = static_cast<std::size_t>(dim);
  Vector v0 = x;
  Vector out(n, Scalar(0.0, 0.0));
  const double c0Half = 0.5 * coeffs[0];
  for (std::size_t i = 0; i < n; ++i) {
    out[i] = c0Half * v0[i];
  }

  if (coeffs.size() > 1) {
    const double invHalfwidth = 1.0 / halfwidth;
    Vector v1(n, Scalar(0.0, 0.0));
    S.Multiply(Scalar(1.0, 0.0), v0, Scalar(0.0, 0.0), v1);
    const double c1 = coeffs[1];
    for (std::size_t i = 0; i < n; ++i) {
      v1[i] = (v1[i] - center * v0[i]) * invHalfwidth;
      out[i] += c1 * v1[i];
    }

    Vector v2(n, Scalar(0.0, 0.0));
    const double twoInvHalfwidth = 2.0 * invHalfwidth;
    for (std::size_t m = 2; m < coeffs.size(); ++m) {
      S.Multiply(Scalar(1.0, 0.0), v1, Scalar(0.0, 0.0), v2);
      const double cm = coeffs[m];
      for (std::size_t i = 0; i < n; ++i) {
        v2[i] = twoInvHalfwidth * (v2[i] - center * v1[i]) - v0[i];
        out[i] += cm * v2[i];
      }
      v0.swap(v1);
      v1.swap(v2);
    }
  }

  if (b == Scalar(0.0, 0.0)) {
    for (std::size_t i = 0; i < n; ++i) {
      y[i] = a * out[i];
    }
  } else {
    for (std::size_t i = 0; i < n; ++i) {
      y[i] = a * out[i] + b * y[i];
    }
  }
}

class ScopedEigenSingleThread {
public:
  ScopedEigenSingleThread() : previousThreads_(Eigen::nbThreads()) {
    Eigen::setNbThreads(1);
  }
  ~ScopedEigenSingleThread() { Eigen::setNbThreads(previousThreads_); }
  ScopedEigenSingleThread(const ScopedEigenSingleThread&) = delete;
  ScopedEigenSingleThread& operator=(const ScopedEigenSingleThread&) = delete;

private:
  int previousThreads_;
};

SparseMatrix CompressDenseMatrix(const Eigen::MatrixXcd& dense) {
  const Index dim = static_cast<Index>(dense.rows());
  const double sqThreshold = kExactInverseDropThreshold * kExactInverseDropThreshold;
  std::vector<Index> rowPointers(static_cast<std::size_t>(dim + 1), 0);
  std::vector<Index> columns;
  std::vector<Scalar> values;

  for (Index i = 0; i < dim; ++i) {
    rowPointers[static_cast<std::size_t>(i)] = static_cast<Index>(columns.size());
    for (Index j = 0; j < dim; ++j) {
      const Scalar val = dense(i, j);
      if (std::norm(val) > sqThreshold) {
        columns.push_back(j);
        values.push_back(val);
      }
    }
  }
  rowPointers[static_cast<std::size_t>(dim)] = static_cast<Index>(columns.size());
  return SparseMatrix::FromCsr(dim, std::move(rowPointers), std::move(columns),
                               std::move(values));
}

}  // namespace

ExactInverseOverlap::ExactInverseOverlap(const SparseMatrix& S, Index maxDimension)
    : dimension_(S.Dimension()) {
  if (maxDimension <= 0) {
    throw std::invalid_argument("ExactInverseOverlap: maxDimension must be positive");
  }
  if (dimension_ <= 0) {
    throw std::invalid_argument("ExactInverseOverlap: S dimension must be positive");
  }
  if (dimension_ > maxDimension) {
    throw std::invalid_argument("ExactInverseOverlap: S dimension exceeds maxDimension");
  }
  if (!S.IsHermitian()) {
    throw std::invalid_argument("ExactInverseOverlap: S must be Hermitian");
  }

  const Eigen::MatrixXcd denseS(S.Eigen());
  Eigen::MatrixXcd denseInv;
  {
    ScopedEigenSingleThread singleThread;
    Eigen::LLT<Eigen::MatrixXcd> llt(denseS);
    if (llt.info() != Eigen::Success) {
      throw std::runtime_error(
          "ExactInverseOverlap: exact Cholesky (LLT) failed; S is not positive definite");
    }
    denseInv = llt.solve(Eigen::MatrixXcd::Identity(dimension_, dimension_));
  }
  inverse_ = CompressDenseMatrix(denseInv);
}

Index ExactInverseOverlap::Dimension() const { return dimension_; }

void ExactInverseOverlap::Multiply(Scalar a, const Vector& x, Scalar b, Vector& y) const {
  inverse_.Multiply(a, x, b, y);
}

ExactInverseSqrtOverlap::ExactInverseSqrtOverlap(const SparseMatrix& S, Index maxDimension)
    : dimension_(S.Dimension()) {
  if (maxDimension <= 0) {
    throw std::invalid_argument("ExactInverseSqrtOverlap: maxDimension must be positive");
  }
  if (dimension_ <= 0) {
    throw std::invalid_argument("ExactInverseSqrtOverlap: S dimension must be positive");
  }
  if (dimension_ > maxDimension) {
    throw std::invalid_argument("ExactInverseSqrtOverlap: S dimension exceeds maxDimension");
  }
  if (!S.IsHermitian()) {
    throw std::invalid_argument("ExactInverseSqrtOverlap: S must be Hermitian");
  }

  const Eigen::MatrixXcd denseS(S.Eigen());
  Eigen::MatrixXcd denseInvSqrt;
  {
    ScopedEigenSingleThread singleThread;
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> es(denseS);
    if (es.info() != Eigen::Success) {
      throw std::runtime_error("ExactInverseSqrtOverlap: dense eigensolver failed on S");
    }
    const Eigen::VectorXd evals = es.eigenvalues();
    if (evals.minCoeff() <= 0.0) {
      throw std::runtime_error("ExactInverseSqrtOverlap: S is not positive definite");
    }
    const Eigen::VectorXd invSqrtEvals = evals.cwiseSqrt().cwiseInverse();
    denseInvSqrt =
        es.eigenvectors() * invSqrtEvals.asDiagonal() * es.eigenvectors().adjoint();
  }
  inverseSqrt_ = CompressDenseMatrix(denseInvSqrt);
}

Index ExactInverseSqrtOverlap::Dimension() const { return dimension_; }

void ExactInverseSqrtOverlap::Multiply(Scalar a, const Vector& x, Scalar b, Vector& y) const {
  inverseSqrt_.Multiply(a, x, b, y);
}

IschoInverseOverlap::IschoInverseOverlap(const SparseMatrix& S, double tolerance, int maxIterations)
    : S_(&S), tolerance_(tolerance), maxIterations_(maxIterations) {
  if (S.Dimension() <= 0) {
    throw std::invalid_argument("IschoInverseOverlap: S dimension must be positive");
  }
  if (!S.IsHermitian()) {
    throw std::invalid_argument("IschoInverseOverlap: S must be Hermitian");
  }
  if (!(tolerance > 0.0) || !std::isfinite(tolerance)) {
    throw std::invalid_argument("IschoInverseOverlap: tolerance must be positive and finite");
  }
  if (maxIterations <= 0) {
    throw std::invalid_argument("IschoInverseOverlap: maxIterations must be positive");
  }

  Eigen::SparseMatrix<Scalar, Eigen::ColMajor, Index> colS(S.Eigen());
  ICholType ichol;
  ichol.compute(colS);
  if (ichol.info() != Eigen::Success) {
    throw std::runtime_error(
        "IschoInverseOverlap: incomplete Cholesky (ISCHO) factorisation failed");
  }

  factorL_ = ichol.matrixL();
  scale_ = ichol.scalingS();
  perm_ = ichol.permutationP();
}

Index IschoInverseOverlap::Dimension() const { return S_->Dimension(); }

void IschoInverseOverlap::ApplyPreconditioner(const Eigen::VectorXcd& rhs,
                                              Eigen::VectorXcd& out) const {
  Eigen::VectorXcd x = perm_ * rhs;
  x = scale_.asDiagonal() * x;
  Eigen::VectorXcd y = factorL_.template triangularView<Eigen::Lower>().solve(x);
  Eigen::VectorXcd z = factorL_.adjoint().template triangularView<Eigen::Upper>().solve(y);
  z = scale_.asDiagonal() * z;
  out = perm_.inverse() * z;
}

void IschoInverseOverlap::Multiply(Scalar a, const Vector& x, Scalar b, Vector& y) const {
  const Index dim = Dimension();
  if (static_cast<Index>(x.size()) != dim || static_cast<Index>(y.size()) != dim) {
    throw std::invalid_argument("IschoInverseOverlap::Multiply: vector dimension mismatch");
  }
  if (&x == &y) {
    throw std::invalid_argument(
        "IschoInverseOverlap::Multiply: input and output vectors must not alias");
  }

  Eigen::Map<const Eigen::VectorXcd> rhs(x.data(), dim);
  Eigen::Map<Eigen::VectorXcd> yMap(y.data(), dim);

  const double rhsNorm = rhs.norm();
  if (rhsNorm == 0.0) {
    lastIterations_ = 0;
    lastError_ = 0.0;
    if (b == Scalar(0.0, 0.0)) {
      yMap.setZero();
    } else {
      yMap *= b;
    }
    return;
  }

  const auto& matS = S_->Eigen();

  Eigen::VectorXcd sol(dim);
  ApplyPreconditioner(rhs, sol);

  Eigen::VectorXcd r = rhs - matS * sol;
  double relErr = r.norm() / rhsNorm;
  if (relErr <= tolerance_) {
    lastIterations_ = 0;
    lastError_ = relErr;
    if (b == Scalar(0.0, 0.0)) {
      yMap = a * sol;
    } else {
      yMap = a * sol + b * yMap;
    }
    return;
  }

  Eigen::VectorXcd z(dim);
  ApplyPreconditioner(r, z);
  Eigen::VectorXcd p = z;
  double rzOld = r.dot(z).real();

  Eigen::VectorXcd Ap(dim);
  int iter = 0;
  for (iter = 1; iter <= maxIterations_; ++iter) {
    Ap.noalias() = matS * p;
    const double pAp = p.dot(Ap).real();
    if (!(pAp > 0.0)) {
      lastIterations_ = iter;
      lastError_ = relErr;
      std::ostringstream oss;
      oss << "IschoInverseOverlap::Multiply: non-positive curvature p^+ S p = " << pAp
          << " at iteration " << iter << " (maxIterations=" << maxIterations_
          << ", residual=" << relErr << ", tolerance=" << tolerance_ << ")";
      throw std::runtime_error(oss.str());
    }
    const double alpha = rzOld / pAp;
    sol.noalias() += alpha * p;
    r.noalias() -= alpha * Ap;

    relErr = r.norm() / rhsNorm;
    if (relErr <= tolerance_) {
      break;
    }

    ApplyPreconditioner(r, z);
    const double rzNew = r.dot(z).real();
    const double beta = rzNew / rzOld;
    p = z + beta * p;
    rzOld = rzNew;
  }

  lastIterations_ = std::min(iter, maxIterations_);
  lastError_ = relErr;
  if (!(relErr <= tolerance_)) {
    std::ostringstream oss;
    oss << "IschoInverseOverlap::Multiply: failed to converge within maxIterations="
        << maxIterations_ << " (residual=" << relErr << " > tolerance=" << tolerance_ << ")";
    throw std::runtime_error(oss.str());
  }

  if (b == Scalar(0.0, 0.0)) {
    yMap = a * sol;
  } else {
    yMap = a * sol + b * yMap;
  }
}

ChebyshevInverseOverlap::ChebyshevInverseOverlap(const SparseMatrix& S, SpectralBounds boundsOfS,
                                                 int order, double coeffTol)
    : S_(&S) {
  if (S.Dimension() <= 0) {
    throw std::invalid_argument("ChebyshevInverseOverlap: S dimension must be positive");
  }
  if (!S.IsHermitian()) {
    throw std::invalid_argument("ChebyshevInverseOverlap: S must be Hermitian");
  }
  ComputeChebyshevGaussCoeffs(
      boundsOfS, order, coeffTol, [](double x) { return 1.0 / x; }, center_, halfwidth_, coeffs_,
      effectiveOrder_);
}

Index ChebyshevInverseOverlap::Dimension() const { return S_->Dimension(); }

void ChebyshevInverseOverlap::Multiply(Scalar a, const Vector& x, Scalar b, Vector& y) const {
  ApplyChebyshevSeries(*S_, center_, halfwidth_, coeffs_, a, x, b, y);
}

ChebyshevInverseSqrtOverlap::ChebyshevInverseSqrtOverlap(const SparseMatrix& S,
                                                         SpectralBounds boundsOfS, int order,
                                                         double coeffTol)
    : S_(&S) {
  if (S.Dimension() <= 0) {
    throw std::invalid_argument("ChebyshevInverseSqrtOverlap: S dimension must be positive");
  }
  if (!S.IsHermitian()) {
    throw std::invalid_argument("ChebyshevInverseSqrtOverlap: S must be Hermitian");
  }
  ComputeChebyshevGaussCoeffs(
      boundsOfS, order, coeffTol, [](double x) { return 1.0 / std::sqrt(x); }, center_, halfwidth_,
      coeffs_, effectiveOrder_);
}

Index ChebyshevInverseSqrtOverlap::Dimension() const { return S_->Dimension(); }

void ChebyshevInverseSqrtOverlap::Multiply(Scalar a, const Vector& x, Scalar b, Vector& y) const {
  ApplyChebyshevSeries(*S_, center_, halfwidth_, coeffs_, a, x, b, y);
}

}  // namespace lsquant
