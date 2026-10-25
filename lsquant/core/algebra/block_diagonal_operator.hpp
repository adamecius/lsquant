/**
 * @file block_diagonal_operator.hpp
 * @brief Block-diagonal dense linear operator for k-space representations.
 *
 * Stores `N_k` dense `W x W` complex blocks contiguously and evaluates the
 * fused matrix-vector operation `y = a * O * x + b * y` as an OpenMP-parallel
 * batched block GEMV across the `N_k` independent blocks.
 */
#ifndef LSQUANT_CORE_ALGEBRA_BLOCK_DIAGONAL_OPERATOR_HPP_
#define LSQUANT_CORE_ALGEBRA_BLOCK_DIAGONAL_OPERATOR_HPP_

#include <vector>

#include "core/algebra/linear_operator.hpp"
#include "core/algebra/types.hpp"

namespace lsquant {

/**
 * @brief Block-diagonal operator holding `N_k` dense `W x W` complex blocks.
 *
 * Solves the problem of applying k-local operators (`H(k)`, velocity `dH/dk`,
 * or per-k unitary rotations `U(k)`) in `O(N_k * W^2)` time with contiguous
 * memory access and OpenMP parallelism over `k`. Assumes vectors are ordered
 * cell/block-major: entry `(ik, alpha)` resides at flat index `ik * W + alpha`,
 * and block `ik` row `alpha`, column `beta` resides at
 * `blocks[(ik * W + alpha) * W + beta]`.
 *
 * @warning Throws `std::invalid_argument` if `numBlocks < 0`, `blockSize < 0`,
 *          or `blocks.size() != numBlocks * blockSize * blockSize`.
 * @par Pitfall
 *      Each `W x W` block is stored in row-major order; passing column-major
 *      blocks without transposing exchanges the operator with its transpose.
 */
class BlockDiagonalOperator : public LinearOperator {
public:
  /**
   * @brief Construct an empty 0-dimensional block-diagonal operator.
   */
  BlockDiagonalOperator() = default;

  /**
   * @brief Construct a block-diagonal operator from contiguous row-major blocks.
   *
   * @param numBlocks Number of diagonal blocks `N_k >= 0`.
   * @param blockSize Dimension `W >= 0` of each square block.
   * @param blocks    Contiguous array of length `N_k * W * W` in row-major order.
   *
   * @warning Throws `std::invalid_argument` if dimensions are negative or the
   *          buffer size does not equal `numBlocks * blockSize * blockSize`.
   * @par Pitfall
   *      Hermiticity is evaluated and cached once at construction time within
   *      tolerance `1e-12`.
   */
  BlockDiagonalOperator(Index numBlocks, Index blockSize, std::vector<Scalar> blocks);

  /**
   * @brief Return the total Hilbert-space dimension `N_k * W`.
   *
   * @return Dimension `N_k * W >= 0`.
   */
  Index Dimension() const override;

  /**
   * @brief Return the number of diagonal blocks `N_k`.
   *
   * @return Block count `N_k >= 0`.
   */
  Index NumBlocks() const { return numBlocks_; }

  /**
   * @brief Return the size `W` of each square block.
   *
   * @return Block dimension `W >= 0`.
   */
  Index BlockSize() const { return blockSize_; }

  /**
   * @brief Return a const reference to the contiguous `N_k * W * W` block buffer.
   *
   * @return Contiguous block array.
   */
  const std::vector<Scalar>& Blocks() const { return blocks_; }

  /**
   * @brief Evaluate `y = a * O * x + b * y` via OpenMP-parallel per-block GEMV.
   *
   * @param a Multiplicative scalar factor for `O * x`.
   * @param x Input vector of length `Dimension()`.
   * @param b Multiplicative scalar factor for `y`.
   * @param y Input/output vector of length `Dimension()`.
   *
   * @warning Throws `std::invalid_argument` if `x.size()` or `y.size()` does not
   *          equal `Dimension()`.
   * @par Pitfall
   *      Supports in-place evaluation (`&x == &y`) by staging each block's
   *      `W`-element input slice before overwriting the output slice.
   */
  void Multiply(Scalar a, const Vector& x, Scalar b, Vector& y) const override;
  void MultiplyBatch(Scalar a, const std::vector<Vector>& x, Scalar b,
                     std::vector<Vector>& y) const override;

  using LinearOperator::Multiply;

  /**
   * @brief Return cached Hermiticity of the block-diagonal matrix.
   *
   * @return `true` if every block is self-adjoint within `1e-12`, `false` otherwise.
   */
  bool IsHermitian() const override { return isHermitian_; }

  /**
   * @brief Return the conjugate transpose `O^+` as a new `BlockDiagonalOperator`.
   *
   * @return Block-diagonal operator whose `ik`-th block is `B_{ik}^+`.
   */
  BlockDiagonalOperator Adjoint() const;

private:
  Index numBlocks_ = 0;
  Index blockSize_ = 0;
  std::vector<Scalar> blocks_;
  bool isHermitian_ = true;
};

}  // namespace lsquant

#endif  // LSQUANT_CORE_ALGEBRA_BLOCK_DIAGONAL_OPERATOR_HPP_
