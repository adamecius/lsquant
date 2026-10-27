/**
 * @file composite_operators.hpp
 * @brief Algebraic composite operators implementing LinearOperator.
 *
 * Provides matrix-free operator combinators (`IdentityOperator`, `DiagonalOperator`,
 * `ScaledOperator`, `OperatorSum`, `OperatorProduct`, `AdjointOperator`,
 * `SandwichOperator`, `LowdinOperator`, `MetricInverseProduct`) so that diagonal
 * operators, linear combinations, products, conjugate transposes, unitary similarity
 * transforms, symmetric Löwdin sandwiches, and one-sided metric-inverse products can be
 * passed into the Chebyshev engine without materialising intermediate dense matrices.
 */
#ifndef LSQUANT_CORE_ALGEBRA_COMPOSITE_OPERATORS_HPP_
#define LSQUANT_CORE_ALGEBRA_COMPOSITE_OPERATORS_HPP_

#include <memory>
#include <optional>
#include <vector>

#include "core/algebra/linear_operator.hpp"
#include "core/algebra/sparse_matrix.hpp"
#include "core/algebra/types.hpp"

namespace lsquant {

/**
 * @brief Matrix-free identity operator `I` on a complex vector space of dimension `N`.
 *
 * Implements `y = a * x + b * y` in `O(N)` time and reports `IsIdentity() == true` so
 * that callers can bypass redundant identity applications in `O(1)`.
 *
 * @warning Throws `std::invalid_argument` if constructed with a negative dimension or
 *          applied to vectors of mismatched size.
 * @par Pitfall
 *      Even when `IsIdentity()` is true, calling `Multiply(a, x, b, y)` still applies
 *      the scalar weights `a` and `b`.
 */
class IdentityOperator : public LinearOperator {
public:
  /**
   * @brief Construct an identity operator of dimension `dimension`.
   *
   * @param dimension Vector-space dimension `N >= 0`.
   *
   * @warning Throws `std::invalid_argument` if `dimension < 0`.
   * @par Pitfall
   *      Ensure `dimension` matches the Hamiltonian dimension when mixing operators.
   */
  explicit IdentityOperator(Index dimension);

  /**
   * @brief Return the operator dimension `N`.
   *
   * @return Dimension `N >= 0`.
   */
  Index Dimension() const override;

  /**
   * @brief Evaluate `y = a * x + b * y`.
   *
   * @param a Scalar multiplier for `x`.
   * @param x Input vector of length `Dimension()`.
   * @param b Scalar multiplier for `y`.
   * @param y Input/output vector of length `Dimension()`.
   */
  void Multiply(Scalar a, const Vector& x, Scalar b, Vector& y) const override;

  using LinearOperator::Multiply;

  /**
   * @brief Return `true` to indicate structural identity.
   *
   * @return Always `true`.
   */
  bool IsIdentity() const override { return true; }

  /**
   * @brief Return `true` because the identity map is self-adjoint.
   *
   * @return Always `true`.
   */
  bool IsHermitian() const override { return true; }

private:
  Index dimension_;
};

/**
 * @brief Matrix-free diagonal operator `diag(d_0, ..., d_{N-1})`.
 *
 * Solves the problem of applying diagonal potentials (such as real-space Anderson
 * disorder `V`) or diagonal eigenvalue tables `eps_n(k)` in `O(N)` time and memory
 * via OpenMP-parallel elementwise multiplication `y[i] = a * d[i] * x[i] + b * y[i]`.
 *
 * @warning Throws `std::invalid_argument` if applied to vectors whose length does
 *          not equal `Dimension()`.
 * @par Pitfall
 *      Hermiticity holds if and only if every diagonal entry has imaginary part
 *      within `1e-12` of zero.
 */
class DiagonalOperator : public LinearOperator {
public:
  /**
   * @brief Construct a diagonal operator from its diagonal entries.
   *
   * @param values Vector of `N` complex diagonal elements.
   */
  explicit DiagonalOperator(std::vector<Scalar> values);

  /**
   * @brief Return the operator dimension `N`.
   *
   * @return Dimension `N >= 0`.
   */
  Index Dimension() const override;

  /**
   * @brief Evaluate `y[i] = a * d[i] * x[i] + b * y[i]` across all `i in 0..N-1`.
   *
   * @param a Scalar multiplier for `D * x`.
   * @param x Input vector of length `Dimension()`.
   * @param b Scalar multiplier for `y`.
   * @param y Input/output vector of length `Dimension()`.
   */
  void Multiply(Scalar a, const Vector& x, Scalar b, Vector& y) const override;

  using LinearOperator::Multiply;

  /**
   * @brief Return `true` iff all diagonal entries are real within `1e-12`.
   *
   * @return Cached Hermiticity flag.
   */
  bool IsHermitian() const override { return isHermitian_; }

  /**
   * @brief Return the conjugate diagonal operator `diag(conj(d))`.
   *
   * @return Adjoint `DiagonalOperator`.
   */
  DiagonalOperator Adjoint() const;

private:
  std::vector<Scalar> values_;
  bool isHermitian_ = true;
};

/**
 * @brief Matrix-free scalar multiple `alpha * O` of an existing `LinearOperator`.
 *
 * Forwards `Multiply(a, x, b, y)` to `op.Multiply(a * alpha, x, b, y)` with zero extra
 * vector allocations and `O(1)` overhead beyond the wrapped operator.
 *
 * @warning Does not own the wrapped `LinearOperator`; `op` must outlive this instance.
 * @par Pitfall
 *      A Hermitian operator scaled by a complex number with non-zero imaginary part is
 *      no longer Hermitian.
 */
class ScaledOperator : public LinearOperator {
public:
  /**
   * @brief Construct the scaled operator `alpha * op`.
   *
   * @param alpha Complex scalar factor.
   * @param op    Referenced operator (must outlive `*this`).
   */
  ScaledOperator(Scalar alpha, const LinearOperator& op);

  /**
   * @brief Return the dimension of the wrapped operator.
   *
   * @return Dimension `N >= 0`.
   */
  Index Dimension() const override;

  /**
   * @brief Evaluate `y = (a * alpha) * O * x + b * y`.
   *
   * @param a Scalar multiplier for `(alpha * O) * x`.
   * @param x Input vector of length `Dimension()`.
   * @param b Scalar multiplier for `y`.
   * @param y Input/output vector of length `Dimension()`.
   */
  void Multiply(Scalar a, const Vector& x, Scalar b, Vector& y) const override;

  using LinearOperator::Multiply;

  /**
   * @brief Return `true` iff the wrapped operator is Hermitian and `alpha` is real.
   *
   * @return Hermiticity of `alpha * O`.
   */
  bool IsHermitian() const override;

private:
  Scalar alpha_;
  const LinearOperator* op_;
};

/**
 * @brief Matrix-free sum `O1 + O2` of two linear operators of equal dimension.
 *
 * Evaluates `y = a * (O1 + O2) * x + b * y` as `y = a * O1 * x + b * y` followed by
 * `y = a * O2 * x + 1 * y`. Complexity is `O(cost(O1) + cost(O2))` with no temporary
 * vector allocation unless `x` and `y` alias.
 *
 * @warning Both referenced operators must outlive this `OperatorSum` instance.
 * @par Pitfall
 *      Throws `std::invalid_argument` at construction if `first.Dimension() != second.Dimension()`.
 */
class OperatorSum : public LinearOperator {
public:
  /**
   * @brief Construct the sum `first + second`.
   *
   * @param first  First summand operator.
   * @param second Second summand operator (must match `first.Dimension()`).
   *
   * @warning Throws `std::invalid_argument` if dimensions differ.
   * @par Pitfall
   *      Neither summand is copied; keep both alive for the lifetime of this object.
   */
  OperatorSum(const LinearOperator& first, const LinearOperator& second);

  /**
   * @brief Return the common operator dimension.
   *
   * @return Dimension `N >= 0`.
   */
  Index Dimension() const override;

  /**
   * @brief Evaluate `y = a * (O1 + O2) * x + b * y`.
   *
   * @param a Scalar multiplier for `(O1 + O2) * x`.
   * @param x Input vector of length `Dimension()`.
   * @param b Scalar multiplier for `y`.
   * @param y Input/output vector of length `Dimension()`.
   */
  void Multiply(Scalar a, const Vector& x, Scalar b, Vector& y) const override;

  using LinearOperator::Multiply;

  /**
   * @brief Return `true` if both summands report `IsHermitian() == true`.
   *
   * @return Sufficient condition for Hermiticity of `O1 + O2`.
   */
  bool IsHermitian() const override;

private:
  const LinearOperator* first_;
  const LinearOperator* second_;
};

/**
 * @brief Matrix-free composition `O1 * O2` of two linear operators (`O2` applied first).
 *
 * Evaluates `temp = O2 * x` into a temporary vector of length `Dimension()`, then
 * `y = a * O1 * temp + b * y`. Complexity is `O(cost(O1) + cost(O2) + N)`.
 *
 * @warning Both referenced operators must outlive this `OperatorProduct` instance.
 * @par Pitfall
 *      In operator notation `O1 * O2`, `second` (`O2`) acts on the ket `x` first, and
 *      `first` (`O1`) acts on the intermediate result.
 */
class OperatorProduct : public LinearOperator {
public:
  /**
   * @brief Construct the composition `first * second` (`second` applied first to `x`).
   *
   * @param first  Outer operator `O1`.
   * @param second Inner operator `O2` applied first to the input ket.
   *
   * @warning Throws `std::invalid_argument` if `first.Dimension() != second.Dimension()`.
   * @par Pitfall
   *      A product of two Hermitian operators is generally non-Hermitian unless they commute.
   */
  OperatorProduct(const LinearOperator& first, const LinearOperator& second);

  /**
   * @brief Return the common operator dimension.
   *
   * @return Dimension `N >= 0`.
   */
  Index Dimension() const override;

  /**
   * @brief Evaluate `y = a * O1 * (O2 * x) + b * y`.
   *
   * @param a Scalar multiplier for `(O1 * O2) * x`.
   * @param x Input vector of length `Dimension()`.
   * @param b Scalar multiplier for `y`.
   * @param y Input/output vector of length `Dimension()`.
   */
  void Multiply(Scalar a, const Vector& x, Scalar b, Vector& y) const override;

  using LinearOperator::Multiply;

  /**
   * @brief Report whether the product is guaranteed to be Hermitian.
   *
   * @return `true` if one factor is the identity and the other is Hermitian; `false` otherwise.
   */
  bool IsHermitian() const override;

private:
  const LinearOperator* first_;
  const LinearOperator* second_;
};

/**
 * @brief Conjugate transpose operator `O^+`.
 *
 * Supports `SparseMatrixOperator` (via `SparseMatrix::Adjoint()`), `BlochTransform`
 * (via `BlochTransform::Adjoint()`), `BlockDiagonalOperator` (via
 * `BlockDiagonalOperator::Adjoint()`), `DiagonalOperator` (via
 * `DiagonalOperator::Adjoint()`), and any Hermitian `LinearOperator` (forwarded
 * directly to `op.Multiply`).
 *
 * @warning Throws `std::invalid_argument` if constructed from a non-Hermitian
 *          `LinearOperator` whose concrete type does not provide an adjoint.
 * @par Pitfall
 *      Construct `AdjointOperator` outside tight inner loops so any explicit
 *      adjoint buffers are allocated only once.
 */
class AdjointOperator : public LinearOperator {
public:
  /**
   * @brief Construct the adjoint `op^+` of a `LinearOperator`.
   *
   * @param op Referenced operator.
   *
   * @warning Throws `std::invalid_argument` if `op` is non-Hermitian and not one
   *          of the supported transform/matrix operator types.
   * @par Pitfall
   *      If `op` is forwarded without an owned copy, `op` must outlive `*this`.
   */
  explicit AdjointOperator(const LinearOperator& op);

  /**
   * @brief Construct the adjoint `matrix^+` directly from a `SparseMatrix`.
   *
   * @param matrix Source sparse matrix whose conjugate transpose is stored.
   */
  explicit AdjointOperator(const SparseMatrix& matrix);

  /**
   * @brief Return the operator dimension.
   *
   * @return Dimension `N >= 0`.
   */
  Index Dimension() const override;

  /**
   * @brief Evaluate `y = a * O^+ * x + b * y`.
   *
   * @param a Scalar multiplier for `O^+ * x`.
   * @param x Input vector of length `Dimension()`.
   * @param b Scalar multiplier for `y`.
   * @param y Input/output vector of length `Dimension()`.
   */
  void Multiply(Scalar a, const Vector& x, Scalar b, Vector& y) const override;

  using LinearOperator::Multiply;

  /**
   * @brief Report whether `O^+` is Hermitian.
   *
   * @return `true` if `O` is Hermitian, `false` otherwise.
   */
  bool IsHermitian() const override;

private:
  std::optional<SparseMatrix> adjointMatrix_;
  std::shared_ptr<const LinearOperator> ownedAdjoint_;
  const LinearOperator* forwardedOp_ = nullptr;
  bool isHermitian_ = false;
};

/**
 * @brief Matrix-free unitary similarity / congruence operator `B * V * B^+`.
 *
 * Solves the problem of applying an operator `V` defined in one representation
 * (such as diagonal Anderson disorder in real space) to vectors in another
 * representation (such as Bloch k-space) without assembling the dense matrix
 * `B * V * B^+`. Evaluates `y = a * (B * V * B^+) * x + b * y` in three steps
 * using two internal scratch vectors:
 * `t1 = B^+ * x`, `t2 = V * t1`, `y = a * B * t2 + b * y`.
 *
 * `Multiply` uses per-object mutable scratch buffers (`t1_`, `t2_`) and is not
 * safe to call concurrently on the same `SandwichOperator` object (the engine
 * never does).
 *
 * @warning `B` and `V` must have the same `Dimension()` and must outlive this
 *          `SandwichOperator` instance. `Multiply` is not reentrant across
 *          threads sharing one `SandwichOperator` instance.
 * @par Pitfall
 *      In `SandwichOperator(B, V)`, `B^+` acts on the input ket `x` first, then
 *      `V`, and finally `B` maps back to the output representation.
 */
class SandwichOperator : public LinearOperator {
public:
  /**
   * @brief Construct `B * V * B^+` from transform `B` and inner operator `V`.
   *
   * @param B Outer transform operator (applied last; `B^+` is applied first).
   * @param V Inner operator applied between `B` and `B^+`.
   *
   * @warning Throws `std::invalid_argument` if `B.Dimension() != V.Dimension()`
   *          or if `AdjointOperator(B)` cannot be constructed.
   * @par Pitfall
   *      Both `B` and `V` are held by pointer; keep them alive for the lifetime
   *      of this `SandwichOperator`.
   */
  SandwichOperator(const LinearOperator& B, const LinearOperator& V);

  /**
   * @brief Return the operator dimension `N`.
   *
   * @return Dimension `N >= 0`.
   */
  Index Dimension() const override;

  /**
   * @brief Evaluate `y = a * (B * V * B^+) * x + b * y`.
   *
   * @param a Scalar multiplier for `(B * V * B^+) * x`.
   * @param x Input vector of length `Dimension()`.
   * @param b Scalar multiplier for `y`.
   * @param y Input/output vector of length `Dimension()`.
   */
  void Multiply(Scalar a, const Vector& x, Scalar b, Vector& y) const override;

  using LinearOperator::Multiply;

  /**
   * @brief Return `true` iff the inner operator `V` is Hermitian.
   *
   * @return `V.IsHermitian()`.
   */
  bool IsHermitian() const override;

private:
  const LinearOperator* b_;
  const LinearOperator* v_;
  AdjointOperator bAdjoint_;
  mutable Vector t1_;
  mutable Vector t2_;
};

/**
 * @brief Matrix-free symmetric Löwdin sandwich `S^{-1/2} O S^{-1/2}`.
 *
 * Solves the problem of transforming an operator `O'` with overlap `S` into the
 * orthonormal Löwdin coefficient representation without materialising
 * dense `N x N` matrices. Evaluates `t1 = S^{-1/2} x`, `t2 = O t1`, and
 * `y = a * S^{-1/2} t2 + b * y` using two temporary vectors of length `Dimension()` per
 * `Multiply` call. Complexity is `O(2 * cost(S^{-1/2}) + cost(O) + N)`.
 *
 * @warning Neither `inverseSqrt` nor `O` is owned; both must outlive this instance.
 * @par Pitfall
 *      Both factors must have identical `Dimension()` or the constructor throws
 *      `std::invalid_argument`.
 */
class LowdinOperator : public LinearOperator {
public:
  /**
   * @brief Construct the symmetric sandwich `inverseSqrt * O * inverseSqrt`.
   *
   * @param inverseSqrt Inverse square-root overlap applier `S^{-1/2}`.
   * @param O           Inner operator `O'` in the non-orthogonal representation.
   *
   * @warning Throws `std::invalid_argument` if `inverseSqrt.Dimension() != O.Dimension()`.
   * @par Pitfall
   *      For `IsHermitian()` to return `true`, both `inverseSqrt` and `O` must report
   *      `IsHermitian() == true`.
   */
  LowdinOperator(const LinearOperator& inverseSqrt, const LinearOperator& O);

  /**
   * @brief Return the common operator dimension `N`.
   *
   * @return Dimension `N >= 0`.
   */
  Index Dimension() const override;

  /**
   * @brief Evaluate `y = a * (S^{-1/2} * O * S^{-1/2}) * x + b * y`.
   *
   * @param a Scalar multiplier for `(S^{-1/2} O S^{-1/2}) x`.
   * @param x Input vector of length `Dimension()`.
   * @param b Scalar multiplier for `y`.
   * @param y Input/output vector of length `Dimension()`.
   */
  void Multiply(Scalar a, const Vector& x, Scalar b, Vector& y) const override;

  using LinearOperator::Multiply;

  /**
   * @brief Return `true` iff both `inverseSqrt` and `O` are Hermitian.
   *
   * @return Hermiticity of `S^{-1/2} O S^{-1/2}`.
   */
  bool IsHermitian() const override;

private:
  const LinearOperator* inverseSqrt_;
  const LinearOperator* op_;
};

/**
 * @brief Matrix-free one-sided product `S^{-1} O` for the generalised non-orthogonal route.
 *
 * Solves the problem of representing the effective operator `S^{-1} O'` in the
 * coefficient representation where `O'` maps coordinates to covariant dual components and
 * `S^{-1}` raises the index back to contravariant coordinates. Evaluates `t = O x` into a
 * single temporary vector and then `y = a * S^{-1} t + b * y`. Complexity is
 * `O(cost(S^{-1}) + cost(O) + N)`.
 *
 * @warning Reports `IsHermitian() == false` unconditionally because `S^{-1} O` is
 *          `S`-self-adjoint rather than Euclidean-Hermitian whenever `S != I`.
 * @par Pitfall
 *      Never assume `(S^{-1} O)^+ = S^{-1} O`; the Euclidean adjoint is `O^+ S^{-1}`,
 *      which is why the generalised Kubo bra path requires explicit `S^{-1}` and `S`
 *      insertions.
 */
class MetricInverseProduct : public LinearOperator {
public:
  /**
   * @brief Construct the one-sided product `inverse * O`.
   *
   * @param inverse Inverse overlap applier `S^{-1}`.
   * @param O       Right factor `O'`.
   *
   * @warning Throws `std::invalid_argument` if `inverse.Dimension() != O.Dimension()`.
   * @par Pitfall
   *      Neither referenced operator is copied; both must outlive this object.
   */
  MetricInverseProduct(const LinearOperator& inverse, const LinearOperator& O);

  /**
   * @brief Return the common operator dimension `N`.
   *
   * @return Dimension `N >= 0`.
   */
  Index Dimension() const override;

  /**
   * @brief Evaluate `y = a * S^{-1} * (O * x) + b * y`.
   *
   * @param a Scalar multiplier for `(S^{-1} O) x`.
   * @param x Input vector of length `Dimension()`.
   * @param b Scalar multiplier for `y`.
   * @param y Input/output vector of length `Dimension()`.
   */
  void Multiply(Scalar a, const Vector& x, Scalar b, Vector& y) const override;

  using LinearOperator::Multiply;

  /**
   * @brief Return `false` because `S^{-1} O` is not Euclidean-Hermitian in general.
   *
   * @return Always `false`.
   */
  bool IsHermitian() const override { return false; }

private:
  const LinearOperator* inverse_;
  const LinearOperator* op_;
};

/**
 * @brief Apply c*A*M*B + conj(c)*B*M*A for three Hermitian operators.
 *
 * c=1/2 gives a symmetric product; c=i gives a Hermitian commutator.
 * No intermediate matrix is assembled. Cost is two applications of each
 * factor and two length-N scratch vectors. All references must outlive this
 * operator. Scratch is local, and each factor retains its own reentrancy rules.
 */
class HermitianTripleProduct : public LinearOperator {
public:
  HermitianTripleProduct(const LinearOperator& first, const LinearOperator& middle,
                        const LinearOperator& last, Scalar coefficient = Scalar(0.5, 0.0));
  Index Dimension() const override;
  void Multiply(Scalar a, const Vector& x, Scalar b, Vector& y) const override;
  using LinearOperator::Multiply;
  bool IsHermitian() const override { return true; }

private:
  const LinearOperator* first_;
  const LinearOperator* middle_;
  const LinearOperator* last_;
  Scalar coefficient_;
};

}  // namespace lsquant

#endif  // LSQUANT_CORE_ALGEBRA_COMPOSITE_OPERATORS_HPP_
