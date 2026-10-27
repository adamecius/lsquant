/**
 * @file fft_moment_contraction.hpp
 * @brief Fast Fourier/cosine-transform contraction of rank-2 Chebyshev moment matrices on Chebyshev nodes.
 */
#ifndef LSQUANT_CORE_ALGEBRA_FFT_MOMENT_CONTRACTION_HPP_
#define LSQUANT_CORE_ALGEBRA_FFT_MOMENT_CONTRACTION_HPP_

#include "core/algebra/types.hpp"

namespace lsquant {

/**
 * @brief Evaluate `sum_{m,n} mu[m,n] cos(m*theta) exp(-i*n*theta)` and
 *        `sum_{m,n} n * mu[m,n] cos(m*theta) exp(-i*n*theta)` on the `points`
 *        Chebyshev nodes `theta_j = (2*j + 0.5) * pi / points`.
 *
 * Returns two contiguous rows of `points` complex values (`O(M0*M1 + N log N)` work).
 */
std::vector<Scalar> FftMomentContraction(const Scalar* moments, Index m0, Index m1,
                                         Index points);

}  // namespace lsquant

#endif  // LSQUANT_CORE_ALGEBRA_FFT_MOMENT_CONTRACTION_HPP_
