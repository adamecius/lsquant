/**
 * @file legacy_kernels.hpp
 * @brief Gate-only compatibility evaluators for legacy retarded Green's function tables.
 *
 * Reproduces `greenR_chebF` and `DgreenR_chebF` from `include/chebyshev_coefficients.hpp`
 * using `std::pow(std::complex<double>, double)` from libstdc++ and the exact three-division
 * arithmetic order so that stage-3 golden curves match the legacy C++ reconstruction byte for byte.
 * This header belongs to `lsquant/compat/` and is never included by `lsquant/core/`.
 */
#ifndef LSQUANT_COMPAT_LEGACY_KERNELS_HPP_
#define LSQUANT_COMPAT_LEGACY_KERNELS_HPP_

#include <cmath>
#include <complex>
#include <vector>

#include "core/algebra/types.hpp"

namespace lsquant {
namespace compat {

/**
 * @brief Evaluate the `m`-th Chebyshev coefficient of the retarded Green's function at `x`.
 *
 * Evaluates `-i * std::pow(x - i * f0, m) / f0` with `f0 = std::sqrt(1.0 - x * x)` using
 * `std::pow(std::complex<double>, double)` to match `greenR_chebF(double, const double)` bit
 * for bit.
 *
 * @param x Rescaled energy coordinate in `(-1, 1)`.
 * @param m Chebyshev polynomial index as a `double`.
 * @return  Complex retarded Green's function Chebyshev factor.
 *
 * @warning Singular at `|x| >= 1.0` where `f0 == 0.0`.
 * @par Pitfall
 *      Keep the exponent `m` typed as `double` so overload resolution selects
 *      `std::pow(const std::complex<double>&, const double&)` rather than integer exponentiation.
 */
inline std::complex<double> greenR_chebF(double x, const double m) {
  const std::complex<double> I(0.0, 1.0);
  const double f0 = std::sqrt(1.0 - x * x);
  const std::complex<double> fm = std::pow(x - I * f0, m);
  return -I * fm / f0;
}

/**
 * @brief Evaluate the `m`-th Chebyshev coefficient of the derivative of the retarded Green's
 *        function at `x`.
 *
 * Evaluates `-i * std::pow(x - i * f0, m) * (x + i * m * f0) / f0 / f0 / f0` with
 * `f0 = std::sqrt(1.0 - x * x)` using `std::pow(std::complex<double>, double)` and three
 * sequential divisions by `f0` to match `DgreenR_chebF(double, const double)` bit for bit.
 *
 * @param x Rescaled energy coordinate in `(-1, 1)`.
 * @param m Chebyshev polynomial index as a `double`.
 * @return  Complex retarded Green's function derivative Chebyshev factor.
 *
 * @warning Singular at `|x| >= 1.0` where `f0 == 0.0`.
 * @par Pitfall
 *      Preserve the three sequential divisions `/ f0 / f0 / f0`; replacing them with
 *      `/ (f0 * f0 * f0)` alters IEEE-754 rounding in the last bit.
 */
inline std::complex<double> DgreenR_chebF(double x, const double m) {
  const std::complex<double> I(0.0, 1.0);
  const double f0 = std::sqrt(1.0 - x * x);
  const std::complex<double> fm = std::pow(x - I * f0, m) * (x + I * m * f0);
  return -I * fm / f0 / f0 / f0;
}

/**
 * @brief Populate row-major `(numMoments0, numPoints)` `greenR` and `(numMoments1, numPoints)`
 *        `DgreenR` tables over `xGrid`.
 *
 * Solves the problem of generating the retarded Green's function and its derivative tables
 * with the exact `libstdc++` `std::pow(std::complex<double>, double)` rounding of the legacy
 * C++ drivers. Complexity is `O((numMoments0 + numMoments1) * numPoints)`.
 *
 * @param xGrid       Pointer to rescaled coordinates `x_i in (-1, 1)` of length `numPoints`.
 * @param numPoints   Number of grid points `N >= 0`.
 * @param numMoments0 Number of moments `M0 > 0` for `greenR`.
 * @param numMoments1 Number of moments `M1 > 0` for `DgreenR`.
 * @param greenR      Output buffer of length `numMoments0 * numPoints` (row-major).
 * @param dgreenR     Output buffer of length `numMoments1 * numPoints` (row-major).
 *
 * @warning Throws `std::invalid_argument` if `numMoments0 <= 0`, `numMoments1 <= 0`,
 *          `numPoints < 0`, or any required pointer is null.
 * @par Pitfall
 *      Tables are stored in row-major `(M, N)` layout with moment order `m` as the slow axis
 *      and grid index `i` as the fast axis.
 */
void ComputeLegacyGreenTables(const double* xGrid, Index numPoints, Index numMoments0,
                              Index numMoments1, Scalar* greenR, Scalar* dgreenR);

/**
 * @brief Compute the Chebyshev Bessel series coefficients `c_n` of `exp(-i x y)` using GSL.
 *
 * Reproduces the coefficient sequence of legacy `MomentsTD::Evolve`:
 * `c_0 = J_0(x)`, `c_n = (-i)^n * 2 * J_n(x)` for `n >= 1`, stopping at the first
 * `n >= 1` where `|J_n(x)| <= tol` (that term excluded), using `gsl_sf_bessel_Jn(n, x)`
 * and the exact complex multiplication order `nIp * Scalar(2.0, 0.0) * Jn`.
 *
 * @param x   Dimensionless time-frequency product `dt * HalfWidth / hbar`.
 * @param tol Cutoff tolerance on `|J_n(x)|` (default `1e-15`).
 * @return    Vector of complex Chebyshev series coefficients `[c_0, c_1, ..., c_{L-1}]`.
 *
 * @warning Throws `std::invalid_argument` if `tol <= 0.0`.
 * @par Pitfall
 *      Does not include the band-center phase `exp(-i BandCenter dt / hbar)`, which is
 *      applied after the Chebyshev recursion via `linalg::scal`.
 */
std::vector<Scalar> ComputeLegacyBesselCoefficients(double x, double tol = 1e-15);

}  // namespace compat
}  // namespace lsquant

#endif  // LSQUANT_COMPAT_LEGACY_KERNELS_HPP_
