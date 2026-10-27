/**
 * @file overlap_appliers.hpp
 * @brief LinearOperator implementations for inverse and inverse-square-root overlap application.
 *
 * Ports the four overlap preconditioner/applier strategies from `lsquant-official`:
 * - `ExactInverseOverlap` (and `ExactInverseSqrtOverlap` for dense Löwdin reference),
 * - `IschoInverseOverlap` (Preconditioned Conjugate Gradient with incomplete stabilised Cholesky),
 * - `ChebyshevInverseOverlap` (Chebyshev expansion of `1/x` on `spec(S)`),
 * - `ChebyshevInverseSqrtOverlap` (Chebyshev expansion of `x^-1/2` on `spec(S)`).
 *
 * Note: The Chebyshev coefficients of `1/x` and `x^-1/2` on `[a, b]` in
 * `ChebyshevInverseOverlap` and `ChebyshevInverseSqrtOverlap` are computed in C++ by
 * Chebyshev-Gauss quadrature with adaptive truncation `|c_m| < coeffTol * |c_0|`.
 * This is the ONE place a function of an operator is expanded in C++, because it is an
 * applier, not an observable.
 */
#ifndef LSQUANT_CORE_ALGEBRA_OVERLAP_APPLIERS_HPP_
#define LSQUANT_CORE_ALGEBRA_OVERLAP_APPLIERS_HPP_

#include <Eigen/Cholesky>
#include <Eigen/Dense>
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/OrderingMethods>
#include <Eigen/Sparse>
#include <vector>

#include "core/algebra/linear_operator.hpp"
#include "core/algebra/sparse_matrix.hpp"
#include "core/algebra/spectral_bounds.hpp"
#include "core/algebra/types.hpp"

namespace lsquant {

/**
 * @brief Magnitude threshold below which entries of the dense `S^-1` and `S^-1/2` matrices
 *        are pruned when compressing into internal CSR storage.
 *
 * Set to `1e-18`, below IEEE-754 double-precision rounding noise (`2.22e-16`), so the exact
 * appliers stay exact while dropping exact structural zeros.
 */
inline constexpr double kExactInverseDropThreshold = 1e-18;

/**
 * @brief Exact dense Cholesky (`LLT`) inverse applier `y = a * S^-1 * x + b * y`.
 *
 * Solves the problem of providing a machine-precision reference applier for `S^-1` on small to
 * moderate system sizes (`Dimension() <= maxDimension`). Materialises `S` as a dense matrix,
 * computes its exact Hermitian positive-definite Cholesky factorisation `Eigen::LLT`, and
 * stores entries of `S^-1 = LLT(S).solve(I)` with magnitude above `kExactInverseDropThreshold`
 * (`1e-18`, below double rounding, so the exact appliers stay exact). Construction costs
 * `O(N^3)` time and `O(N^2)` memory; each `Multiply` costs `O(N^2)` time.
 *
 * @warning Throws `std::invalid_argument` if `S.Dimension() > maxDimension` or `!S.IsHermitian()`,
 *          and `std::runtime_error` if `S` is not positive definite (`LLT` fails).
 * @par Pitfall
 *      Intended only for reference and validation runs where `N <= maxDimension`; use
 *      `IschoInverseOverlap` or `ChebyshevInverseOverlap` for large sparse systems.
 */
class ExactInverseOverlap : public LinearOperator {
public:
  /**
   * @brief Construct the exact dense inverse overlap operator `S^-1`.
   *
   * @param S            Hermitian positive-definite sparse overlap matrix.
   * @param maxDimension Upper bound on `S.Dimension()` allowed for dense inversion.
   *
   * @warning Throws `std::invalid_argument` if `maxDimension <= 0`, `S.Dimension() <= 0`,
   *          `S.Dimension() > maxDimension`, or `!S.IsHermitian()`. Throws `std::runtime_error`
   *          if `Eigen::LLT` fails.
   * @par Pitfall
   *      Materialises an `N x N` dense complex matrix of `16 * N^2` bytes.
   */
  explicit ExactInverseOverlap(const SparseMatrix& S, Index maxDimension = 4096);

  /**
   * @brief Return the operator dimension `N`.
   *
   * @return Square dimension `N > 0`.
   */
  Index Dimension() const override;

  /**
   * @brief Compute `y = a * S^-1 * x + b * y` using the dense inverse matrix.
   *
   * @param a Scalar multiplier for `S^-1 * x`.
   * @param x Input vector of length `Dimension()`.
   * @param b Scalar multiplier for `y`.
   * @param y Input/output vector of length `Dimension()`.
   */
  void Multiply(Scalar a, const Vector& x, Scalar b, Vector& y) const override;

  using LinearOperator::Multiply;

  /**
   * @brief Return `true` because the inverse of a Hermitian positive-definite matrix is Hermitian.
   *
   * @return Always `true`.
   */
  bool IsHermitian() const override { return true; }

private:
  Index dimension_ = 0;
  SparseMatrix inverse_;
};

/**
 * @brief Exact dense inverse square-root overlap applier `y = a * S^-1/2 * x + b * y`.
 *
 * Solves the problem of providing a machine-precision dense reference for `S^-1/2` on the
 * Löwdin route (`--route lowdin --applier exact`) when `Dimension() <= maxDimension`.
 * Diagonalises the dense overlap matrix `S = U diag(lambda) U^+` via
 * `Eigen::SelfAdjointEigenSolver` and stores entries of `S^-1/2 = U diag(lambda^-1/2) U^+`
 * with magnitude above `kExactInverseDropThreshold` (`1e-18`, below double rounding, so the
 * exact appliers stay exact). Construction costs `O(N^3)` time and `O(N^2)` memory; each
 * `Multiply` costs `O(N^2)` time.
 *
 * @warning Throws `std::invalid_argument` if `S.Dimension() > maxDimension` or `!S.IsHermitian()`,
 *          and `std::runtime_error` if `S` has a non-positive eigenvalue.
 * @par Pitfall
 *      Intended for validation and reference gates on systems with `N <= maxDimension`; use
 *      `ChebyshevInverseSqrtOverlap` for linear-scaling sparse execution.
 */
class ExactInverseSqrtOverlap : public LinearOperator {
public:
  /**
   * @brief Construct the exact dense inverse square-root overlap operator `S^-1/2`.
   *
   * @param S            Hermitian positive-definite sparse overlap matrix.
   * @param maxDimension Upper bound on `S.Dimension()` allowed for dense eigendecomposition.
   *
   * @warning Throws `std::invalid_argument` if `maxDimension <= 0`, `S.Dimension() <= 0`,
   *          `S.Dimension() > maxDimension`, or `!S.IsHermitian()`. Throws `std::runtime_error`
   *          if `SelfAdjointEigenSolver` fails or `lambda_min(S) <= 0`.
   * @par Pitfall
   *      Materialises an `N x N` dense complex matrix of `16 * N^2` bytes.
   */
  explicit ExactInverseSqrtOverlap(const SparseMatrix& S, Index maxDimension = 4096);

  /**
   * @brief Return the operator dimension `N`.
   *
   * @return Square dimension `N > 0`.
   */
  Index Dimension() const override;

  /**
   * @brief Compute `y = a * S^-1/2 * x + b * y` using the dense `S^-1/2` matrix.
   *
   * @param a Scalar multiplier for `S^-1/2 * x`.
   * @param x Input vector of length `Dimension()`.
   * @param b Scalar multiplier for `y`.
   * @param y Input/output vector of length `Dimension()`.
   */
  void Multiply(Scalar a, const Vector& x, Scalar b, Vector& y) const override;

  using LinearOperator::Multiply;

  /**
   * @brief Return `true` because the principal inverse square root of a Hermitian positive-definite
   *        matrix is Hermitian.
   *
   * @return Always `true`.
   */
  bool IsHermitian() const override { return true; }

private:
  Index dimension_ = 0;
  SparseMatrix inverseSqrt_;
};

/**
 * @brief Preconditioned Conjugate Gradient (`PCG`) inverse overlap applier with Incomplete
 *        Stabilised Cholesky (`ISCHO`) preconditioning.
 *
 * Solves `S z = x` iteratively for each `Multiply` call using Hermitian Conjugate Gradient
 * preconditioned by `Eigen::IncompleteCholesky<Scalar, Eigen::Lower, Eigen::AMDOrdering<Index>>`,
 * reproducing the legacy `IschoPreconditioner` forward/backward substitution with explicit
 * conjugate-transpose upper solve `L^+` and initial guess `z_0 = M^-1 x`. Complexity per
 * application is `O(k * nnz(S))` where `k <= maxIterations` is the number of PCG iterations
 * required to reach relative residual `||x - S z||_2 / ||x||_2 <= tolerance`.
 * `LastIterations()` and `LastError()` describe the last solve, and the applier is not safe
 * for concurrent `Multiply` on one object.
 *
 * @warning Throws `std::runtime_error` at construction if the incomplete Cholesky factorisation
 *          cannot stabilise `S` to positive definiteness. Because `Multiply` updates mutable
 *          last-solve diagnostics (`LastIterations()`/`LastError()`), this applier is not safe
 *          for concurrent `Multiply` calls on a single object.
 * @par Pitfall
 *      Choosing `tolerance` too loose (e.g., `1e-5`) introduces non-Hermitian drift that
 *      accumulates across Chebyshev recursion steps `T_m(S^-1 H~)`.
 */
class IschoInverseOverlap : public LinearOperator {
public:
  /**
   * @brief Construct the ISCHO-preconditioned PCG solver for `S^-1`.
   *
   * @param S             Hermitian positive-definite sparse overlap matrix (must outlive `*this`).
   * @param tolerance     Relative residual stopping criterion `||x - S z||_2 / ||x||_2 <= tolerance`.
   * @param maxIterations Hard cap on PCG iterations per `Multiply` call.
   *
   * @warning Throws `std::invalid_argument` if `S.Dimension() <= 0`, `!S.IsHermitian()`,
   *          `tolerance <= 0`, or `maxIterations <= 0`. Throws `std::runtime_error` if
   *          `Eigen::IncompleteCholesky` fails.
   * @par Pitfall
   *      `S` is referenced without copying its row-major structure; `S` must outlive `*this`.
   */
  explicit IschoInverseOverlap(const SparseMatrix& S, double tolerance = 1e-13,
                               int maxIterations = 200);

  /**
   * @brief Return the operator dimension `N`.
   *
   * @return Square dimension `N > 0`.
   */
  Index Dimension() const override;

  /**
   * @brief Solve `S z = x` via ISCHO-PCG and accumulate `y = a * z + b * y`.
   *
   * @param a Scalar multiplier for `S^-1 * x`.
   * @param x Input vector of length `Dimension()`.
   * @param b Scalar multiplier for `y`.
   * @param y Input/output vector of length `Dimension()`.
   *
   * @warning Throws `std::runtime_error` if `p^+ S p <= 0` or if `maxIterations` is reached
   *          with relative residual `relErr > tolerance`.
   */
  void Multiply(Scalar a, const Vector& x, Scalar b, Vector& y) const override;

  using LinearOperator::Multiply;

  /**
   * @brief Return `true` because `S^-1` is self-adjoint up to the solver tolerance.
   *
   * @return Always `true`.
   */
  bool IsHermitian() const override { return true; }

  /**
   * @brief Return the PCG iteration count from the most recent `Multiply` call.
   *
   * @return Number of iterations executed in `[0, maxIterations]`.
   */
  int LastIterations() const { return lastIterations_; }

  /**
   * @brief Return the relative residual `||x - S z||_2 / ||x||_2` from the most recent `Multiply`.
   *
   * @return Relative residual (dimensionless).
   */
  double LastError() const { return lastError_; }

private:
  using OrderingType = Eigen::AMDOrdering<Index>;
  using ICholType = Eigen::IncompleteCholesky<Scalar, Eigen::Lower, OrderingType>;
  using FactorType = typename ICholType::FactorType;
  using PermutationType = typename ICholType::PermutationType;

  void ApplyPreconditioner(const Eigen::VectorXcd& rhs, Eigen::VectorXcd& out) const;

  const SparseMatrix* S_;
  double tolerance_;
  int maxIterations_;
  FactorType factorL_;
  Eigen::VectorXd scale_;
  PermutationType perm_;
  mutable int lastIterations_ = 0;
  mutable double lastError_ = 0.0;
};

/**
 * @brief Matrix-free Chebyshev polynomial expansion of `S^-1` on `[a, b] supset spec(S)`.
 *
 * Approximates `S^-1 x` by expanding `f(lambda) = 1 / lambda` in Chebyshev polynomials on the
 * positive spectral interval `[boundsOfS.lower, boundsOfS.upper]`:
 * `S^-1 x ~= (c_0 / 2) v_0 + sum_{m=1}^{P_eff} c_m v_m`, with `v_0 = x`, `v_1 = S_hat v_0`,
 * `v_{m+1} = 2 S_hat v_m - v_{m-1}`, and `S_hat = (S - center * I) / halfwidth`.
 * The Chebyshev coefficients `c_m` are computed in C++ at construction by Chebyshev-Gauss
 * quadrature on `N = order + 1` nodes (`legacy chebyshev_sinv_preconditioner.hpp`) and
 * adaptively truncated once `|c_m| <= coeffTol * |c_0|` for all trailing `m`.
 * This is the ONE place a function of an operator is expanded in C++, because it is an
 * applier, not an observable.
 *
 * @warning `boundsOfS.lower` must be strictly positive (`0 < lower < upper`) and `[lower, upper]`
 *          must enclose the entire spectrum of `S`.
 * @par Pitfall
 *      Because `1/x` has a pole at `x = 0`, the required polynomial order grows with
 *      `sqrt(cond(S))`; if `boundsOfS.lower` is too close to `0`, increase `order`.
 */
class ChebyshevInverseOverlap : public LinearOperator {
public:
  /**
   * @brief Construct the Chebyshev `S^-1` applier on `boundsOfS` up to maximum order `order`.
   *
   * @param S         Hermitian positive-definite sparse overlap matrix (must outlive `*this`).
   * @param boundsOfS Positive interval `[lower, upper]` enclosing `spec(S)` (`0 < lower < upper`).
   * @param order     Maximum Chebyshev polynomial degree `P >= 1` (`N = order + 1` quadrature nodes).
   * @param coeffTol  Relative adaptive truncation threshold `|c_m| < coeffTol * |c_0|` (default `1e-12`).
   *
   * @warning Throws `std::invalid_argument` if `S.Dimension() <= 0`, `!S.IsHermitian()`,
   *          `boundsOfS.lower <= 0`, `boundsOfS.upper <= boundsOfS.lower`, `order < 1`, or
   *          `coeffTol < 0`.
   * @par Pitfall
   *      Ensure `boundsOfS` includes a small safety padding outside the extremal Ritz values of `S`
   *      so that `S_hat` maps `spec(S)` strictly into `(-1, 1)`.
   */
  ChebyshevInverseOverlap(const SparseMatrix& S, SpectralBounds boundsOfS, int order,
                          double coeffTol = 1e-12);

  /**
   * @brief Return the operator dimension `N`.
   *
   * @return Square dimension `N > 0`.
   */
  Index Dimension() const override;

  /**
   * @brief Evaluate the truncated Chebyshev series for `y = a * S^-1 * x + b * y`.
   *
   * @param a Scalar multiplier for `S^-1 * x`.
   * @param x Input vector of length `Dimension()`.
   * @param b Scalar multiplier for `y`.
   * @param y Input/output vector of length `Dimension()`.
   */
  void Multiply(Scalar a, const Vector& x, Scalar b, Vector& y) const override;

  using LinearOperator::Multiply;

  /**
   * @brief Return `true` because a real-coefficient polynomial of a Hermitian matrix is Hermitian.
   *
   * @return Always `true`.
   */
  bool IsHermitian() const override { return true; }

  /**
   * @brief Return the effective polynomial order after adaptive truncation.
   *
   * @return Truncated degree `P_eff` in `[1, order]`.
   */
  int EffectiveOrder() const { return effectiveOrder_; }

private:
  const SparseMatrix* S_;
  double center_ = 0.0;
  double halfwidth_ = 0.0;
  std::vector<double> coeffs_;
  int effectiveOrder_ = 0;
};

/**
 * @brief Matrix-free Chebyshev polynomial expansion of `S^-1/2` on `[a, b] supset spec(S)`.
 *
 * Approximates `S^-1/2 x` by expanding `f(lambda) = lambda^-1/2` in Chebyshev polynomials on the
 * positive spectral interval `[boundsOfS.lower, boundsOfS.upper]`:
 * `S^-1/2 x ~= (c_0 / 2) v_0 + sum_{m=1}^{P_eff} c_m v_m`.
 * The Chebyshev coefficients `c_m` are computed in C++ at construction by Chebyshev-Gauss
 * quadrature on `N = order + 1` nodes (`legacy chebyshev_sinvsqrt_preconditioner.hpp`) and
 * adaptively truncated once `|c_m| <= coeffTol * |c_0|` for all trailing `m`.
 * This is the ONE place a function of an operator is expanded in C++, because it is an
 * applier, not an observable.
 *
 * @warning `boundsOfS.lower` must be strictly positive (`0 < lower < upper`) and `[lower, upper]`
 *          must enclose the entire spectrum of `S`.
 * @par Pitfall
 *      Used twice per matrix-vector product inside `LowdinOperator` (`S^-1/2 O S^-1/2`), so any
 *      truncation error enters with factor `2` per step.
 */
class ChebyshevInverseSqrtOverlap : public LinearOperator {
public:
  /**
   * @brief Construct the Chebyshev `S^-1/2` applier on `boundsOfS` up to maximum order `order`.
   *
   * @param S         Hermitian positive-definite sparse overlap matrix (must outlive `*this`).
   * @param boundsOfS Positive interval `[lower, upper]` enclosing `spec(S)` (`0 < lower < upper`).
   * @param order     Maximum Chebyshev polynomial degree `P >= 1` (`N = order + 1` quadrature nodes).
   * @param coeffTol  Relative adaptive truncation threshold `|c_m| < coeffTol * |c_0|` (default `1e-12`).
   *
   * @warning Throws `std::invalid_argument` if `S.Dimension() <= 0`, `!S.IsHermitian()`,
   *          `boundsOfS.lower <= 0`, `boundsOfS.upper <= boundsOfS.lower`, `order < 1`, or
   *          `coeffTol < 0`.
   * @par Pitfall
   *      Supply padded spectral bounds so all eigenvalues of `S` lie strictly inside `[lower, upper]`.
   */
  ChebyshevInverseSqrtOverlap(const SparseMatrix& S, SpectralBounds boundsOfS, int order,
                              double coeffTol = 1e-12);

  /**
   * @brief Return the operator dimension `N`.
   *
   * @return Square dimension `N > 0`.
   */
  Index Dimension() const override;

  /**
   * @brief Evaluate the truncated Chebyshev series for `y = a * S^-1/2 * x + b * y`.
   *
   * @param a Scalar multiplier for `S^-1/2 * x`.
   * @param x Input vector of length `Dimension()`.
   * @param b Scalar multiplier for `y`.
   * @param y Input/output vector of length `Dimension()`.
   */
  void Multiply(Scalar a, const Vector& x, Scalar b, Vector& y) const override;

  using LinearOperator::Multiply;

  /**
   * @brief Return `true` because a real-coefficient polynomial of a Hermitian matrix is Hermitian.
   *
   * @return Always `true`.
   */
  bool IsHermitian() const override { return true; }

  /**
   * @brief Return the effective polynomial order after adaptive truncation.
   *
   * @return Truncated degree `P_eff` in `[1, order]`.
   */
  int EffectiveOrder() const { return effectiveOrder_; }

private:
  const SparseMatrix* S_;
  double center_ = 0.0;
  double halfwidth_ = 0.0;
  std::vector<double> coeffs_;
  int effectiveOrder_ = 0;
};

}  // namespace lsquant

#endif  // LSQUANT_CORE_ALGEBRA_OVERLAP_APPLIERS_HPP_
