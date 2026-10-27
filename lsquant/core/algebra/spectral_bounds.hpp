/**
 * @file spectral_bounds.hpp
 * @brief Spectral interval estimators (Lanczos and Gershgorin) for Hermitian operators.
 *
 * Computes `[lower, upper]` enclosures of an operator's spectrum so that the Hamiltonian
 * can be rescaled into the Chebyshev domain `[-1, 1]` without divergence.
 */
#ifndef LSQUANT_CORE_ALGEBRA_SPECTRAL_BOUNDS_HPP_
#define LSQUANT_CORE_ALGEBRA_SPECTRAL_BOUNDS_HPP_

#include "core/algebra/linear_operator.hpp"
#include "core/algebra/sparse_matrix.hpp"

namespace lsquant {

/**
 * @brief Closed interval `[lower, upper]` bounding the spectrum of a Hermitian operator.
 *
 * Holds the minimum and maximum energy bounds in the physical energy units of the operator
 * (typically eV).
 *
 * @warning `upper` must be strictly greater than `lower` before constructing a `SpectralWindow`.
 * @par Pitfall
 *      Underestimating the spectral width causes the Chebyshev recursion `T_m(H~)` to grow
 *      exponentially outside `[-1, 1]`.
 */
struct SpectralBounds {
  double lower;
  double upper;
};

/**
 * @brief Estimate extremal eigenvalues of a Hermitian operator `H` using Lanczos iteration.
 *
 * Runs a `k`-step Hermitian Lanczos tridiagonalisation starting from a deterministic
 * random-phase state seeded with `seed`, then diagonalises the resulting `k x k` symmetric
 * tridiagonal matrix with `Eigen::SelfAdjointEigenSolver` and returns the smallest and
 * largest Ritz values. Complexity is `O(iterations * cost(H) + iterations^2)`.
 *
 * @param H          Hermitian linear operator of dimension `N > 0`.
 * @param iterations Maximum number of Lanczos steps (default `40`).
 * @param seed       PRNG seed for the starting state (default `12345`).
 * @return           Estimated `[lower, upper]` Ritz bounds.
 *
 * @warning Throws `std::invalid_argument` if `H.Dimension() <= 0` or `iterations <= 0`.
 * @par Pitfall
 *      Ritz values are interior approximations to the extreme eigenvalues; add a safety
 *      margin if `iterations` is small and no exact analytical bound is known.
 */
SpectralBounds LanczosBounds(const LinearOperator& H, int iterations = 40, unsigned seed = 12345);

/**
 * @brief Compute a guaranteed symmetric Gershgorin spectral enclosure for a sparse matrix.
 *
 * Evaluates the infinity norm `rho = max_i sum_j |H_{ij}|` across all CSR rows and returns
 * `[-rho * (1 + pad), +rho * (1 + pad)]`, using `rho = 1.0` if the matrix is identically
 * zero. Matches the legacy `utility::SpectralBounds` fallback formula. Complexity is
 * `O(Dimension() + NonZeros())`.
 *
 * @param H   Compressed sparse matrix.
 * @param pad Relative safety padding fraction (default `0.10`, dimensionless).
 * @return    Symmetric spectral enclosure `[-rho*(1+pad), +rho*(1+pad)]`.
 *
 * @warning Throws `std::invalid_argument` if `pad < 0`.
 * @par Pitfall
 *      Gershgorin bounds can be significantly wider than the true bandwidth, reducing energy
 *      resolution for a fixed moment count `M`.
 */
SpectralBounds GershgorinBounds(const SparseMatrix& H, double pad = 0.10);

}  // namespace lsquant

#endif  // LSQUANT_CORE_ALGEBRA_SPECTRAL_BOUNDS_HPP_
