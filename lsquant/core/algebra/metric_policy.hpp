/**
 * @file metric_policy.hpp
 * @brief Overlap metric policy interface and orthogonal/non-orthogonal implementations.
 *
 * Encapsulates the overlap matrix `S` and its inverse appliers (`S^-1` and `S^-1/2`) as
 * `LinearOperator` instances so that non-orthogonal tight-binding models (`H c = E S c`)
 * enter the Chebyshev engine solely through operator compositions without modifying the
 * polynomial recursion or trace state machine.
 */
#ifndef LSQUANT_CORE_ALGEBRA_METRIC_POLICY_HPP_
#define LSQUANT_CORE_ALGEBRA_METRIC_POLICY_HPP_

#include "core/algebra/composite_operators.hpp"
#include "core/algebra/linear_operator.hpp"
#include "core/algebra/sparse_matrix.hpp"
#include "core/algebra/types.hpp"

namespace lsquant {

/**
 * @brief Abstract policy providing the overlap operator `S`, `S^-1`, and `S^-1/2`.
 *
 * Solves the problem of decoupling a non-trivial overlap `S` from the Chebyshev moment
 * engine. Callers query `Overlap()`, `InverseOverlap()`, and `InverseSqrtOverlap()` to
 * assemble either Löwdin-orthogonalised operators (`S^-1/2 O S^-1/2`) or generalised
 * coefficient-representation operators (`S^-1 O`, `S^-1 (a H + b S)`, and the adjoint
 * bra factors `S^-1` and `S`). All three returned operators act on a common Hilbert
 * space of dimension `Overlap().Dimension()`.
 *
 * @warning All operators returned by a `MetricPolicy` reference internal or externally
 *          owned objects whose lifetimes are bound to the `MetricPolicy` instance.
 * @par Pitfall
 *      In the generalised (one-sided) route, `InverseOverlap()` and `Overlap()` must
 *      both be applied on the bra side of rank-2 Kubo words (`S T_m(S^-1 H~) S^-1 V S^-1 |r>`)
 *      to preserve metric covariance of the trace estimator.
 */
class MetricPolicy {
public:
  /**
   * @brief Virtual destructor for polymorphic cleanup.
   */
  virtual ~MetricPolicy() = default;

  /**
   * @brief Report whether the overlap metric is trivial (`S == I`).
   *
   * Allows plan builders to elide redundant metric factors when the representation is
   * already an orthonormal set. Runs in `O(1)` time.
   *
   * @return `true` if `S`, `S^-1`, and `S^-1/2` are all identity operators, `false` otherwise.
   *
   * @warning Returns `false` for `NonOrthogonalMetric` even if the supplied sparse matrix `S`
   *          happens to equal the identity matrix numerically.
   * @par Pitfall
   *      Use this only to select simplified orthogonal plans when the metric is structurally
   *      orthogonal.
   */
  virtual bool IsOrthogonal() const = 0;

  /**
   * @brief Return the overlap operator `S`.
   *
   * @return Const reference to the overlap `LinearOperator` (dimensionless matrix elements).
   *
   * @warning The returned reference is valid only while this `MetricPolicy` lives.
   * @par Pitfall
   *      On the generalised route, `Overlap()` must be applied to the expanded bra batch
   *      before taking Euclidean inner products with the ket batch.
   */
  virtual const LinearOperator& Overlap() const = 0;

  /**
   * @brief Return the inverse overlap applier `S^-1`.
   *
   * @return Const reference to the `LinearOperator` approximating or computing `y = S^-1 x`.
   *
   * @warning The accuracy of generalised-route moments is bounded by the accuracy of this
   *          applier.
   * @par Pitfall
   *      Do not pass a `S^-1/2` applier here; `InverseOverlap()` must approximate `S^-1`.
   */
  virtual const LinearOperator& InverseOverlap() const = 0;

  /**
   * @brief Return the inverse square-root overlap applier `S^-1/2`.
   *
   * @return Const reference to the `LinearOperator` approximating or computing `y = S^-1/2 x`.
   *
   * @warning Used by `LowdinOperator` to form `S^-1/2 O S^-1/2`; requires two applications per
   *          matrix-vector product.
   * @par Pitfall
   *      Do not pass an `S^-1` applier here; `InverseSqrtOverlap()` must satisfy
   *      `(S^-1/2)^2 == S^-1`.
   */
  virtual const LinearOperator& InverseSqrtOverlap() const = 0;
};

/**
 * @brief Trivial metric policy `S = S^-1 = S^-1/2 = I` for an orthonormal set.
 *
 * Holds a single `IdentityOperator` of dimension `N` and returns it from `Overlap()`,
 * `InverseOverlap()`, and `InverseSqrtOverlap()`. Because `IdentityOperator::IsIdentity()`
 * returns `true`, the moment engine automatically elides any identity applications in `O(1)`.
 *
 * @warning Throws `std::invalid_argument` at construction if `dimension < 0`.
 * @par Pitfall
 *      Ensure `dimension` matches the Hamiltonian dimension when registering `Overlap()`
 *      in a `MomentEngine`.
 */
class OrthogonalMetric : public MetricPolicy {
public:
  /**
   * @brief Construct an orthogonal metric of dimension `dimension`.
   *
   * @param dimension Hilbert-space dimension `N >= 0`.
   *
   * @warning Throws `std::invalid_argument` if `dimension < 0`.
   * @par Pitfall
   *      Default-constructing with `dimension = 0` is unsuitable for a positive-dimension
   *      `MomentEngine`.
   */
  explicit OrthogonalMetric(Index dimension = 0);

  /**
   * @brief Return `true` because the metric is the identity.
   *
   * @return Always `true`.
   */
  bool IsOrthogonal() const override { return true; }

  /**
   * @brief Return the `IdentityOperator` representing `S = I`.
   *
   * @return Const reference to the internal `IdentityOperator`.
   */
  const LinearOperator& Overlap() const override;

  /**
   * @brief Return the `IdentityOperator` representing `S^-1 = I`.
   *
   * @return Const reference to the internal `IdentityOperator`.
   */
  const LinearOperator& InverseOverlap() const override;

  /**
   * @brief Return the `IdentityOperator` representing `S^-1/2 = I`.
   *
   * @return Const reference to the internal `IdentityOperator`.
   */
  const LinearOperator& InverseSqrtOverlap() const override;

private:
  IdentityOperator identity_;
};

/**
 * @brief Metric policy for a non-orthogonal coefficient representation with Hermitian
 *        positive-definite overlap `S`.
 *
 * Wraps the sparse overlap matrix `S` in a `SparseMatrixOperator` and binds caller-supplied
 * `LinearOperator` appliers for `S^-1` and `S^-1/2`. Construction takes `O(dim + nnz(S))`
 * to verify Hermiticity of `S` via `SparseMatrixOperator`.
 *
 * @warning `S`, `inverse`, and `inverseSqrt` are referenced without copying and must outlive
 *          this `NonOrthogonalMetric` instance.
 * @par Pitfall
 *      Throws `std::invalid_argument` if `S` is not Hermitian, has non-positive dimension, or
 *      if the dimensions of `S`, `inverse`, and `inverseSqrt` do not all match.
 */
class NonOrthogonalMetric : public MetricPolicy {
public:
  /**
   * @brief Construct a non-orthogonal metric policy from `S` and its `S^-1` and `S^-1/2` appliers.
   *
   * @param S           Hermitian positive-definite sparse overlap matrix (must outlive `*this`).
   * @param inverse     Linear operator applying `S^-1` (must outlive `*this`).
   * @param inverseSqrt Linear operator applying `S^-1/2` (must outlive `*this`).
   *
   * @warning Throws `std::invalid_argument` if `S.Dimension() <= 0`, `!S.IsHermitian()`, or if
   *          `inverse.Dimension()` or `inverseSqrt.Dimension()` differs from `S.Dimension()`.
   * @par Pitfall
   *      All three referenced objects must remain alive for the lifetime of this policy.
   */
  NonOrthogonalMetric(const SparseMatrix& S, const LinearOperator& inverse,
                      const LinearOperator& inverseSqrt);

  /**
   * @brief Return `false` because the metric carries a non-trivial overlap `S`.
   *
   * @return Always `false`.
   */
  bool IsOrthogonal() const override { return false; }

  /**
   * @brief Return the `SparseMatrixOperator` wrapping `S`.
   *
   * @return Const reference to the overlap operator `S`.
   */
  const LinearOperator& Overlap() const override;

  /**
   * @brief Return the `S^-1` applier.
   *
   * @return Const reference to the `inverse` operator supplied at construction.
   */
  const LinearOperator& InverseOverlap() const override;

  /**
   * @brief Return the `S^-1/2` applier.
   *
   * @return Const reference to the `inverseSqrt` operator supplied at construction.
   */
  const LinearOperator& InverseSqrtOverlap() const override;

private:
  SparseMatrixOperator overlap_;
  const LinearOperator* inverse_;
  const LinearOperator* inverseSqrt_;
};

}  // namespace lsquant

#endif  // LSQUANT_CORE_ALGEBRA_METRIC_POLICY_HPP_
