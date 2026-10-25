/**
 * @file chebyshev_vectors.hpp
 * @brief Batched window of Chebyshev polynomial vectors `T_m(H~) |phi>`.
 *
 * Stores a contiguous block of `batchSize` Chebyshev vectors `T_m |phi>` for
 * `m in [FirstOrder(), FirstOrder() + Size())` and supports in-place static operator
 * application across the batch, matching legacy `chebyshev::Vectors`.
 */
#ifndef LSQUANT_CORE_CHEBYSHEV_CHEBYSHEV_VECTORS_HPP_
#define LSQUANT_CORE_CHEBYSHEV_CHEBYSHEV_VECTORS_HPP_

#include <vector>

#include "core/algebra/linear_operator.hpp"
#include "core/algebra/types.hpp"
#include "core/chebyshev/chebyshev_recursion.hpp"

namespace lsquant {

/**
 * @brief Contiguous batch of `Size()` Chebyshev vectors `T_m(H~) |phi>`.
 *
 * Solves the problem of holding a block of `batchSize` expanded Chebyshev states in memory
 * so that rank-2 moment contractions (`<bra_m | ket_n>`) can amortise recursion work across
 * blocks of size `batchSize x batchSize`. `Fill` populates the batch from a running
 * `ChebyshevRecursion` (matching legacy `Vectors::IterateAll`), and `ApplyToAll` transforms
 * every vector in the batch by a static operator `A` (matching legacy `Vectors::Multiply`).
 * Memory footprint is `O(batchSize * dimension)`.
 *
 * @warning `Fill` advances the passed `ChebyshevRecursion` by `Size()` steps.
 * @par Pitfall
 *      `At(m)` uses batch-local indexing `0 <= m < Size()`; the global Chebyshev order of
 *      `At(m)` is `FirstOrder() + m`.
 */
class ChebyshevVectors {
public:
  /**
   * @brief Allocate a batch of `batchSize` vectors of length `dimension`.
   *
   * @param dimension Hilbert-space dimension `N > 0`.
   * @param batchSize Number of vectors `B > 0` held in the batch.
   *
   * @warning Throws `std::invalid_argument` if `dimension <= 0` or `batchSize <= 0`.
   * @par Pitfall
   *      Large `batchSize * dimension` products allocate `16 * B * N` bytes of RAM.
   */
  ChebyshevVectors(Index dimension, int batchSize);

  /**
   * @brief Fill the batch with `Size()` consecutive vectors from `recursion`.
   *
   * Records `FirstOrder() = recursion.Order()`, then for `m = 0 .. Size() - 1` copies
   * `recursion.Current()` into `At(m)` and calls `recursion.Step()`, matching legacy
   * `Vectors::IterateAll`.
   *
   * @param recursion Active `ChebyshevRecursion` whose `Current()` is copied and stepped.
   *
   * @warning Throws `std::invalid_argument` if `recursion.Current().size() != dimension`.
   * @par Pitfall
   *      Leaves `recursion` positioned at order `FirstOrder() + Size()`, ready for the
   *      next batch.
   */
  void Fill(ChebyshevRecursion& recursion);

  /**
   * @brief Apply linear operator `A` in place to every vector in the batch (`v_m <- A v_m`).
   *
   * Reproduces legacy `Vectors::Multiply`: if `A.IsIdentity()`, returns immediately;
   * otherwise copies `At(m)` into a scratch vector and evaluates `A.Multiply(scratch, At(m))`
   * for `m = 0 .. Size() - 1`. Complexity is `O(Size() * cost(A))`.
   *
   * @param A Linear operator of matching dimension.
   *
   * @warning Throws `std::invalid_argument` if `A.Dimension() != dimension`.
   * @par Pitfall
   *      Uses an internal scratch buffer to avoid aliasing `x` and `y` in `A.Multiply`.
   */
  void ApplyToAll(const LinearOperator& A);

  /**
   * @brief Mutable access to batch vector `m` (`0 <= m < Size()`).
   *
   * @param m Batch-local vector index.
   * @return  Mutable reference to the `m`-th vector in the batch.
   *
   * @warning Throws `std::out_of_range` if `m < 0` or `m >= Size()`.
   * @par Pitfall
   *      Do not resize the returned `Vector`.
   */
  Vector& At(int m);

  /**
   * @brief Const access to batch vector `m` (`0 <= m < Size()`).
   *
   * @param m Batch-local vector index.
   * @return  Const reference to the `m`-th vector in the batch.
   *
   * @warning Throws `std::out_of_range` if `m < 0` or `m >= Size()`.
   * @par Pitfall
   *      Index `m` is relative to `FirstOrder()`, not necessarily global order `0`.
   */
  const Vector& At(int m) const;

  /**
   * @brief Return the number of vectors in the batch.
   *
   * @return Batch size `B > 0`.
   */
  int Size() const;

  /**
   * @brief Return the global Chebyshev order corresponding to `At(0)`.
   *
   * @return Polynomial order `first >= 0` set by the most recent `Fill` call.
   */
  int FirstOrder() const;

private:
  Index dimension_;
  int batchSize_;
  int firstOrder_ = 0;
  std::vector<Vector> vectors_;
  Vector scratch_;
#ifdef LSQUANT_NATIVE_BUILD
  std::vector<Vector> batchScratch_;
#endif
};

namespace chebyshev {
using lsquant::ChebyshevVectors;
}  // namespace chebyshev

}  // namespace lsquant

#endif  // LSQUANT_CORE_CHEBYSHEV_CHEBYSHEV_VECTORS_HPP_
