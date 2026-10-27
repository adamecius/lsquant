/**
 * @file spectral_bounds.cpp
 * @brief Implementation of LanczosBounds and GershgorinBounds.
 */
#include "core/algebra/spectral_bounds.hpp"

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <vector>

#include "core/algebra/linalg.hpp"
#include "core/algebra/state_generator.hpp"

namespace lsquant {

SpectralBounds LanczosBounds(const LinearOperator& H, int iterations, unsigned seed) {
  const Index dim = H.Dimension();
  if (dim <= 0) {
    throw std::invalid_argument("LanczosBounds: operator dimension must be positive");
  }
  if (iterations <= 0) {
    throw std::invalid_argument("LanczosBounds: iterations must be positive");
  }

  const int maxSteps = static_cast<int>(std::min<Index>(iterations, dim));
  RandomPhaseStates gen(dim, 1, seed);
  Vector vCurr;
  gen.Next(vCurr);
  const double initNorm = std::sqrt(linalg::norm2(vCurr));
  if (initNorm > 0.0) {
    linalg::scal(Scalar(1.0 / initNorm, 0.0), vCurr);
  }

  Vector vPrev(static_cast<std::size_t>(dim), Scalar(0.0, 0.0));
  Vector w(static_cast<std::size_t>(dim), Scalar(0.0, 0.0));

  std::vector<double> alphas;
  std::vector<double> betas;
  alphas.reserve(static_cast<std::size_t>(maxSteps));
  betas.reserve(static_cast<std::size_t>(maxSteps));

  double betaPrev = 0.0;
  for (int step = 0; step < maxSteps; ++step) {
    H.Multiply(Scalar(1.0, 0.0), vCurr, Scalar(0.0, 0.0), w);
    if (step > 0) {
      linalg::axpy(Scalar(-betaPrev, 0.0), vPrev, w);
    }
    const double alpha = linalg::dot(vCurr, w).real();
    alphas.push_back(alpha);
    linalg::axpy(Scalar(-alpha, 0.0), vCurr, w);

    const double beta = std::sqrt(linalg::norm2(w));
    if (beta < 1e-14 || step + 1 == maxSteps) {
      break;
    }
    betas.push_back(beta);
    vPrev = vCurr;
    vCurr = w;
    linalg::scal(Scalar(1.0 / beta, 0.0), vCurr);
    betaPrev = beta;
  }

  const int m = static_cast<int>(alphas.size());
  Eigen::VectorXd diag(m);
  Eigen::VectorXd subdiag(std::max(0, m - 1));
  for (int i = 0; i < m; ++i) {
    diag(i) = alphas[static_cast<std::size_t>(i)];
  }
  for (int i = 0; i < m - 1; ++i) {
    subdiag(i) = betas[static_cast<std::size_t>(i)];
  }

  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver;
  solver.computeFromTridiagonal(diag, subdiag, Eigen::EigenvaluesOnly);
  if (solver.info() != Eigen::Success) {
    throw std::runtime_error("LanczosBounds: tridiagonal eigenvalue solver failed");
  }
  const auto& evals = solver.eigenvalues();
  double lower = evals.minCoeff();
  double upper = evals.maxCoeff();
  if (!(upper > lower)) {
    lower -= 1.0;
    upper += 1.0;
  }
  return SpectralBounds{lower, upper};
}

SpectralBounds GershgorinBounds(const SparseMatrix& H, double pad) {
  if (pad < 0.0) {
    throw std::invalid_argument("GershgorinBounds: negative pad");
  }
  const auto& mat = H.Eigen();
  double rho = 0.0;
  for (Index r = 0; r < mat.outerSize(); ++r) {
    double rowSum = 0.0;
    for (Eigen::SparseMatrix<Scalar, Eigen::RowMajor, Index>::InnerIterator it(mat, r); it; ++it) {
      rowSum += std::abs(it.value());
    }
    if (rowSum > rho) {
      rho = rowSum;
    }
  }
  if (rho <= 0.0) {
    rho = 1.0;
  }
  const double bound = rho * (1.0 + pad);
  return SpectralBounds{-bound, bound};
}

}  // namespace lsquant
