/**
 * @file linear_operator.hpp
 * @brief Abstract matrix-free linear operator interface for the Chebyshev engine.
 *
 * Decouples the polynomial recursion and moment engine from concrete matrix representations.
 * Every operator (CSR sparse matrix, rescaled Hamiltonian, or algebraic composition) implements
 * the fused multiply-accumulate operation `y = a * O * x + b * y` at whole-vector granularity,
 * paying a single virtual dispatch per matrix-vector product.
 */
#ifndef LSQUANT_CORE_ALGEBRA_LINEAR_OPERATOR_HPP_
#define LSQUANT_CORE_ALGEBRA_LINEAR_OPERATOR_HPP_

#include "core/algebra/types.hpp"
#include <stdexcept>

namespace lsquant {

/**
 * @brief Abstract interface for a linear operator acting on dense complex vectors.
 *
 * Solves the problem of running polynomial recursions and operator words identically over
 * stored sparse matrices and composite operator expressions without templated code bloat.
 * Assumes a finite-dimensional complex Hilbert space of dimension `Dimension()`. Matrix
 * elements carry the caller's physical or dimensionless rescaled units. Virtual dispatch
 * occurs once per whole-vector product in `O(1)` dispatch overhead plus the concrete operator
 * cost.
 *
 * @warning Input and output vectors passed to `Multiply` must already have length `Dimension()`.
 * @par Pitfall
 *      Aliasing `x` and `y` (`&x == &y`) in `Multiply` is invalid unless a concrete subclass
 *      explicitly copies the input buffer first.
 */
class LinearOperator {
public:
  /**
   * @brief Virtual destructor for polymorphic cleanup.
   */
  virtual ~LinearOperator() = default;

  /**
   * @brief Return the Hilbert-space dimension (number of rows and columns) of the operator.
   *
   * Queries the square dimension `N` of the linear map `C^N -> C^N` in `O(1)` time.
   *
   * @return Square operator dimension `N >= 0` (dimensionless count).
   *
   * @warning Never returns a negative value.
   * @par Pitfall
   *      All vectors passed to `Multiply` must match this dimension exactly.
   */
  virtual Index Dimension() const = 0;

  /**
   * @brief Compute the fused matrix-vector operation `y = a * O * x + b * y`.
   *
   * Evaluates the linear action of the operator on `x`, scaled by `a`, and accumulates into
   * `y` scaled by `b`. This four-argument signature directly supports the two-term Chebyshev
   * recursion step `y = 2 * H * x - y` without extra vector passes. Complexity is determined
   * by the concrete operator (`O(nnz)` for sparse CSR matrices).
   *
   * @param a Multiplicative scalar factor applied to `O * x`.
   * @param x Input vector of length `Dimension()`.
   * @param b Multiplicative scalar factor applied to the existing entries of `y`.
   * @param y Input/output vector of length `Dimension()` overwritten with `a * O * x + b * y`.
   *
   * @warning Throws `std::invalid_argument` if `x.size()` or `y.size()` does not equal
   *          `Dimension()`.
   * @par Pitfall
   *      When `b == 0`, existing NaN entries in `y` could propagate if multiplied by zero in
   *      IEEE arithmetic; callers should initialise `y` to finite values.
   */
  virtual void Multiply(Scalar a, const Vector& x, Scalar b, Vector& y) const = 0;

  // Item 5: one dispatch for live vectors. CSR/block implementations reuse the
  // matrix; the fallback preserves other adapters. Tester compares batches
  // against repeated Multiply at every measured thread count.
  virtual void MultiplyBatch(Scalar a, const std::vector<Vector>& x, Scalar b,
                             std::vector<Vector>& y) const {
    if (x.size() != y.size()) throw std::invalid_argument("MultiplyBatch: batch size mismatch");
    for (std::size_t j = 0; j < x.size(); ++j) Multiply(a, x[j], b, y[j]);
  }

  /**
   * @brief Compute the direct matrix-vector product `y = O * x`.
   *
   * Convenience overload forwarding to `Multiply(1, x, 0, y)` via a single virtual call.
   * Complexity matches the four-argument `Multiply`.
   *
   * @param x Input vector of length `Dimension()`.
   * @param y Output vector of length `Dimension()` overwritten with `O * x`.
   *
   * @warning Throws `std::invalid_argument` if `x.size()` or `y.size()` does not equal
   *          `Dimension()`.
   * @par Pitfall
   *      Do not pass the same vector reference for both `x` and `y`; use a temporary vector
   *      when transforming a state in place.
   */
  void Multiply(const Vector& x, Vector& y) const {
    Multiply(Scalar(1.0, 0.0), x, Scalar(0.0, 0.0), y);
  }

  /**
   * @brief Report whether this operator is structurally the identity map.
   *
   * Allows callers to elide redundant operator applications in `O(1)` time when a slot in an
   * operator word is the identity.
   *
   * @return `true` if the operator is guaranteed to be the identity, `false` otherwise.
   *
   * @warning Returns `false` for general sparse matrices even if their entries happen to form
   *          an identity matrix.
   * @par Pitfall
   *      Use this only as an algebraic fast-path; do not rely on it to test numerical identity.
   */
  virtual bool IsIdentity() const { return false; }

  /**
   * @brief Report whether the operator is self-adjoint (`O == O^+`).
   *
   * Used by adjoint wrappers and moment symmetrisation checks to determine whether the
   * operator equals its conjugate transpose.
   *
   * @return `true` if the operator is Hermitian within numerical tolerance, `false` otherwise.
   *
   * @warning For composite products of non-commuting Hermitian operators, this returns `false`.
   * @par Pitfall
   *      Calling an uncached Frobenius norm check inside a tight loop would be expensive;
   *      implementations should precompute or cheaply evaluate this property.
   */
  virtual bool IsHermitian() const = 0;
};

}  // namespace lsquant

#endif  // LSQUANT_CORE_ALGEBRA_LINEAR_OPERATOR_HPP_
