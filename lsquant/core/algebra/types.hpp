/**
 * @file types.hpp
 * @brief Fundamental numeric and container type aliases for the LSQUANT algebra layer.
 *
 * Centralises the 64-bit integer index type, double-precision complex scalar type, and
 * dense vector container used across all operators, state generators, and polynomial
 * recursions. Fixing these aliases in a single header guarantees binary compatibility
 * across sparse storage arrays, dense vector kernels, and moment tensors.
 */
#ifndef LSQUANT_CORE_ALGEBRA_TYPES_HPP_
#define LSQUANT_CORE_ALGEBRA_TYPES_HPP_

#include <complex>
#include <cstdint>
#include <vector>

namespace lsquant {

/**
 * @brief Signed 64-bit integer index used for matrix dimensions, non-zero counts, and CSR offsets.
 *
 * Chosen as a signed 64-bit integer to match Eigen sparse storage index conventions while
 * supporting large real-space tight-binding models whose non-zero count exceeds 32-bit limits.
 *
 * @warning Passing negative values where dimensions or offsets are expected causes runtime
 *          validation exceptions.
 * @par Pitfall
 *      Mixing 32-bit loop counters with CSR row pointer arrays on large supercells overflows
 *      silently; always index CSR arrays with `Index`.
 */
using Index = std::int64_t;

/**
 * @brief Double-precision complex scalar type representing wave-function amplitudes and matrix elements.
 *
 * All tight-binding operators, Chebyshev recursion vectors, and moment tensors operate in
 * IEEE-754 double-precision complex arithmetic.
 *
 * @warning IEEE-754 signed zero (`-0.0`) is preserved by real scalar multiplication but can be
 *          flipped by general complex multiplication with `(a, 0.0)`.
 * @par Pitfall
 *      Avoid constructing temporary complex scalars with zero imaginary part when scaling inner
 *      products if signed-zero preservation is required by downstream text formatters.
 */
using Scalar = std::complex<double>;

/**
 * @brief Contiguous dense vector of complex scalars representing a quantum state in a Hilbert space.
 *
 * Backed by `std::vector<Scalar>` so that `Eigen::Map` views can wrap the underlying buffer with
 * zero copy Overhead.
 *
 * @warning Operations assuming matching vector dimensions throw `std::invalid_argument` if sizes
 *          differ.
 * @par Pitfall
 *      Reallocating a `Vector` invalidates any active `Eigen::Map` or raw pointer view into its
 *      storage.
 */
using Vector = std::vector<Scalar>;

}  // namespace lsquant

#endif  // LSQUANT_CORE_ALGEBRA_TYPES_HPP_
