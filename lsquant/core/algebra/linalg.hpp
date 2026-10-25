/**
 * @file linalg.hpp
 * @brief Dense complex vector linear algebra primitives backed by Eigen::Map.
 *
 * Provides thin, zero-copy wrappers around `Eigen::Map` that reproduce the exact arithmetic
 * and summation order of the legacy LSQUANT `linalg` namespace (`axpy`, `scal`, `copy`,
 * `dot`, `norm2`).
 */
#ifndef LSQUANT_CORE_ALGEBRA_LINALG_HPP_
#define LSQUANT_CORE_ALGEBRA_LINALG_HPP_

#include "core/algebra/types.hpp"

namespace lsquant {
namespace linalg {

/**
 * @brief Compute the in-place complex vector update `y = a * x + y`.
 *
 * Maps `x` and `y` as Eigen column vectors and evaluates `eig_y = a * eig_x + eig_y` in
 * `O(N)` time.
 *
 * @param a Multiplicative scalar factor for `x`.
 * @param x Input vector of length `N`.
 * @param y Input/output vector of length `N`.
 *
 * @warning Throws `std::invalid_argument` if `x.size() != y.size()`.
 * @par Pitfall
 *      Both vectors must be pre-allocated to the same length before calling `axpy`.
 */
void axpy(Scalar a, const Vector& x, Vector& y);

/**
 * @brief Scale a complex vector in place: `x = a * x`.
 *
 * Maps `x` as an Eigen vector and multiplies every entry by `a` in `O(N)` time.
 *
 * @param a Scalar multiplier.
 * @param x Input/output vector scaled in place.
 *
 * @warning Scaling by `0` overwrites existing values but can propagate NaNs if `x` is
 *          uninitialised.
 * @par Pitfall
 *      Prefer `std::fill` when zeroing a freshly allocated vector that may contain non-finite
 *      values.
 */
void scal(Scalar a, Vector& x);

/**
 * @brief Copy complex vector `x` into pre-allocated vector `y`.
 *
 * Maps `x` and `y` as Eigen vectors and assigns `eig_y = eig_x` in `O(N)` time, matching
 * legacy `linalg::copy`.
 *
 * @param x Source vector of length `N`.
 * @param y Destination vector of length `N`.
 *
 * @warning Throws `std::invalid_argument` if `x.size() != y.size()`.
 * @par Pitfall
 *      Unlike `std::vector::operator=`, `linalg::copy` requires `y` to be pre-sized to match
 *      `x` and throws on mismatch.
 */
void copy(const Vector& x, Vector& y);

/**
 * @brief Compute the sesquilinear inner product `conj(x) . y = sum_i conj(x[i]) * y[i]`.
 *
 * Evaluates `Eigen::Map(x).dot(Eigen::Map(y))` sequentially on a single thread in `O(N)`
 * time, matching legacy `linalg::vdot(x, y)` bit for bit.
 *
 * @param x Left-hand (bra) vector of length `N`, conjugated in the product.
 * @param y Right-hand (ket) vector of length `N`.
 * @return  Complex inner product `<x | y>`.
 *
 * @warning Throws `std::invalid_argument` if `x.size() != y.size()`.
 * @par Pitfall
 *      Do not parallelise this reduction with OpenMP; changing the summation tree alters
 *      floating-point rounding and breaks byte identity with the golden moment tables.
 */
Scalar dot(const Vector& x, const Vector& y);

/**
 * @brief Compute the squared Euclidean norm `||x||_2^2 = sum_i |x[i]|^2`.
 *
 * Evaluates `Eigen::Map(x).squaredNorm()` in `O(N)` time, matching legacy `linalg::nrm2`.
 *
 * @param x Input vector of length `N`.
 * @return  Non-negative real squared norm `||x||_2^2`.
 *
 * @warning Returns the squared norm `||x||_2^2`, not `||x||_2`.
 * @par Pitfall
 *      Take `std::sqrt(linalg::norm2(x))` when the Euclidean length itself is needed (e.g.
 *      for state normalisation in Lanczos).
 */
double norm2(const Vector& x);

// Item 5b: GEMM table, row-major output <left_m|right_n>. Fixed output tiles
// own the full inner reduction regardless of team size. Tester compares data
// bytes at 1/24/48/96 threads, and numerical residuals against scalar dot.
std::vector<Scalar> MomentTable(const std::vector<const Vector*>& left,
                               const std::vector<const Vector*>& right);

}  // namespace linalg
}  // namespace lsquant

#endif  // LSQUANT_CORE_ALGEBRA_LINALG_HPP_
