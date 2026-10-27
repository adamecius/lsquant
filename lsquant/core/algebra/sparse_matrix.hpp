/**
 * @file sparse_matrix.hpp
 * @brief Compressed Sparse Row (CSR) complex matrix and its LinearOperator adapter.
 *
 * Wraps `Eigen::SparseMatrix<Scalar, Eigen::RowMajor, Index>` with deterministic,
 * thread-invariant row-parallel matrix-vector multiplication and affine rescaling
 * that matches the legacy LSQUANT numerical order bit for bit.
 */
#ifndef LSQUANT_CORE_ALGEBRA_SPARSE_MATRIX_HPP_
#define LSQUANT_CORE_ALGEBRA_SPARSE_MATRIX_HPP_

#include <Eigen/Sparse>
#include <vector>

#include "core/algebra/linear_operator.hpp"
#include "core/algebra/types.hpp"

namespace lsquant {

/**
 * @brief Square complex sparse matrix in Compressed Sparse Row (CSR) format.
 *
 * Stores tight-binding Hamiltonians and static observables in compressed row-major
 * layout (`Eigen::SparseMatrix<Scalar, Eigen::RowMajor, Index>`). Solves the problem
 * of high-throughput sparse matrix-vector multiplication with strict OpenMP thread-count
 * bit invariance: each row reduction is performed sequentially in stored column order by
 * a single thread. Matrix entries carry the physical units of the operator (e.g., eV for
 * Hamiltonians, eV*Angstrom for velocity matrices) or dimensionless units after rescaling.
 * Storage complexity is `O(dim + nnz)` and matrix-vector product complexity is `O(nnz)`.
 *
 * @warning All indices must be within `[0, Dimension())` and CSR storage is always kept
 *          compressed.
 * @par Pitfall
 *      Do not parallelise across columns within a single row in `Multiply`; changing the
 *      reduction order alters IEEE-754 rounding and breaks byte-for-byte golden identity.
 */
class SparseMatrix {
public:
  /**
   * @brief Construct an empty `0 x 0` sparse matrix.
   */
  SparseMatrix() = default;

  /**
   * @brief Build a compressed CSR sparse matrix directly from raw CSR buffers.
   *
   * Maps the supplied row pointer, column index, and value buffers directly into an Eigen
   * row-major sparse matrix via `Eigen::Map<const Eigen::SparseMatrix<...>>` without staging
   * into intermediate `std::vector` containers. Performs a single copy from the mapped buffers
   * into the owned `Eigen::SparseMatrix` storage. Complexity is `O(dim + nnz)`.
   *
   * @param dim             Square matrix dimension `N >= 0`.
   * @param rowPointersSize Number of elements in `rowPointers` (must equal `dim + 1`).
   * @param rowPointers     Pointer to CSR row offset buffer of length `dim + 1`.
   * @param columnsSize     Number of elements in `columns` (must equal `nnz`).
   * @param columns         Pointer to CSR column index buffer of length `nnz`.
   * @param valuesSize      Number of elements in `values` (must equal `nnz`).
   * @param values          Pointer to complex non-zero values buffer of length `nnz`.
   * @return                Constructed and compressed `SparseMatrix`.
   *
   * @warning Throws `std::invalid_argument` if buffer sizes, pointers, row pointer monotonicity,
   *          or column bounds are violated.
   * @par Pitfall
   *      Column indices within each row should be sorted in ascending order for canonical
   *      Eigen CSR layout.
   */
  static SparseMatrix FromCsr(Index dim, Index rowPointersSize, const Index* rowPointers,
                              Index columnsSize, const Index* columns, Index valuesSize,
                              const Scalar* values);

  /**
   * @brief Build a compressed CSR sparse matrix from raw CSR arrays.
   *
   * Maps the supplied row pointer, column index, and value arrays directly into an Eigen
   * row-major sparse matrix without altering stored values. Complexity is `O(dim + nnz)`.
   *
   * @param dim         Square matrix dimension `N >= 0`.
   * @param rowPointers CSR row offset array of length `dim + 1`, starting at `0` and ending at `nnz`.
   * @param columns     CSR column index array of length `nnz`, with entries in `[0, dim)`.
   * @param values      Complex non-zero values array of length `nnz`.
   * @return            Constructed and compressed `SparseMatrix`.
   *
   * @warning Throws `std::invalid_argument` if array sizes, row pointer monotonicity, or
   *          column bounds are violated.
   * @par Pitfall
   *      Column indices within each row should be sorted in ascending order for canonical
   *      Eigen CSR layout.
   */
  static SparseMatrix FromCsr(Index dim, std::vector<Index> rowPointers,
                              std::vector<Index> columns, std::vector<Scalar> values);

  /**
   * @brief Build a compressed CSR sparse matrix from coordinate (COO) triplets.
   *
   * Assembles a `dim x dim` sparse matrix from parallel arrays of row indices, column
   * indices, and complex values, summing duplicate coordinates via Eigen's `setFromTriplets`
   * without folding or altering values. Complexity is `O(nnz log(nnz / dim))`.
   *
   * @param dim  Square matrix dimension `N >= 0`.
   * @param rows Row indices of length `nnz` in `[0, dim)`.
   * @param cols Column indices of length `nnz` in `[0, dim)`.
   * @param vals Complex values of length `nnz`.
   * @return     Constructed and compressed `SparseMatrix`.
   *
   * @warning Throws `std::invalid_argument` if `rows`, `cols`, and `vals` differ in size or
   *          contain out-of-range indices.
   * @par Pitfall
   *      Duplicate `(row, col)` entries are summed, not overwritten.
   */
  static SparseMatrix FromTriplets(Index dim, const std::vector<Index>& rows,
                                   const std::vector<Index>& cols, const std::vector<Scalar>& vals);

  /**
   * @brief Return the square matrix dimension `N`.
   *
   * @return Number of rows (equal to number of columns) in `O(1)` time.
   *
   * @warning Returns `0` for a default-constructed matrix.
   * @par Pitfall
   *      Always check `Dimension()` before allocating state vectors.
   */
  Index Dimension() const;

  /**
   * @brief Return the number of stored non-zero entries (`nnz`).
   *
   * @return Count of compressed non-zero elements in `O(1)` time.
   *
   * @warning Structural zeros explicitly present in CSR arrays count towards `NonZeros()`.
   * @par Pitfall
   *      Rescaling by `a * M + b * I` may increase `NonZeros()` if diagonal entries were absent.
   */
  Index NonZeros() const;

  /**
   * @brief Compute the row-parallel sparse matrix-vector product `y = a * M * x + b * y`.
   *
   * Evaluates each row `i` as an independent sequential sum `acc = sum_p vl[p] * x[cl[p]]`
   * in stored column order, then writes `y[i] = a * acc + b * y[i]`. Rows are always
   * distributed across OpenMP threads with static scheduling (`#pragma omp parallel for
   * schedule(static)`). Because each row reduction is owned by a single thread and summed
   * strictly in stored column order without cross-thread reductions, the IEEE-754 rounding
   * sequence for every `y[i]` is independent of the thread count, guaranteeing bit-identical
   * results across all OpenMP thread counts and matching the Eigen sparse product on recursion
   * scalars `(a, b) in {(1,0), (2,-1)}`. Complexity is `O(nnz)` work and `O(Dimension())`
   * parallel tasks.
   *
   * @param a Scalar multiplier for `M * x`.
   * @param x Input vector of length `Dimension()`.
   * @param b Scalar multiplier for `y`.
   * @param y Input/output vector of length `Dimension()`.
   *
   * @warning Throws `std::invalid_argument` if `x.size()` or `y.size()` differs from `Dimension()`,
   *          or if `&x == &y`.
   * @par Pitfall
   *      Never pass the same vector as both `x` and `y`; reading `x` while updating `y` across
   *      rows causes a data race and mathematical corruption.
   */
  void Multiply(Scalar a, const Vector& x, Scalar b, Vector& y) const;
  void MultiplyBatch(Scalar a, const std::vector<Vector>& x, Scalar b,
                     std::vector<Vector>& y) const;

  /**
   * @brief Return a new sparse matrix equal to the affine transformation `a * M + b * I`.
   *
   * Evaluates the exact Eigen sparse expression `a * matrix_ + bID` used by the
   * legacy CSR rescaling path, without mutating `*this`. Complexity is `O(dim + nnz)`.
   *
   * @param a Multiplicative scale factor applied to `M`.
   * @param b Additive shift factor multiplying the identity matrix `I`.
   * @return  New compressed `SparseMatrix` holding `a * M + b * I`.
   *
   * @warning Evaluates in place on a fresh copy so that diagonal shifts undergo the exact
   *          IEEE-754 rounding of the stored rescaled matrix.
   * @par Pitfall
   *      Replacing this precomputed matrix sum with an on-the-fly `a * (M * x) + b * x` vector
   *      combination alters rounding at the `1e-16` level and breaks byte-identical goldens.
   */
  SparseMatrix Rescaled(Scalar a, Scalar b) const;

  /**
   * @brief Return the conjugate transpose `M^+` as a new compressed sparse matrix.
   *
   * Transposes the sparsity pattern and conjugates all non-zero values via Eigen's
   * `.adjoint()`. Complexity is `O(dim + nnz)`.
   *
   * @return New compressed `SparseMatrix` representing `M^+`.
   *
   * @warning Allocates a full copy of the sparse structure.
   * @par Pitfall
   *      For known Hermitian matrices, `Adjoint()` is redundant except for sign of `-0.0` on
   *      diagonal imaginary parts.
   */
  SparseMatrix Adjoint() const;

  /**
   * @brief Return the Hermitian fold `0.5 * (M + M^+)` using the legacy Eigen expression.
   *
   * Evaluates `transp = matrix_.adjoint(); result = 0.5 * (transp + matrix_)` matching
   * the legacy CSR-reader convention
   * that enforces exact self-adjointness and normalises IEEE-754 `-0.0` components to `+0.0`
   * when paired with `+0.0` in the conjugate transpose. Complexity is `O(dim + nnz)`.
   *
   * @return New compressed `SparseMatrix` holding `0.5 * (M^+ + M)`.
   *
   * @warning Only meaningful for matrices intended to represent Hermitian operators; applying
   *          this to a non-Hermitian matrix projects onto its Hermitian part.
   * @par Pitfall
   *      Omitting this fold when loading legacy `.CSR` files leaves `-0.0` real parts in
   *      purely imaginary velocity operators (`VX`), shifting signed zeros in moment outputs.
   */
  SparseMatrix HermitianFold() const;

  /**
   * @brief Test whether the matrix is self-adjoint within relative Frobenius tolerance `tol`.
   *
   * Computes `||M - M^+||_F / ||M||_F <= tol` (or `||M - M^+||_F <= tol` if `||M||_F == 0`)
   * in `O(dim + nnz)` time.
   *
   * @param tol Relative Frobenius-norm tolerance (default `1e-12`, dimensionless).
   * @return    `true` if `M` equals `M^+` within `tol`, `false` otherwise.
   *
   * @warning Throws `std::invalid_argument` if `tol < 0`.
   * @par Pitfall
   *      This check constructs a temporary adjoint matrix; cache the result if queried
   *      repeatedly.
   */
  bool IsHermitian(double tol = 1e-12) const;

  /**
   * @brief Access the underlying Eigen row-major compressed sparse matrix.
   *
   * @return Const reference to `Eigen::SparseMatrix<Scalar, Eigen::RowMajor, Index>`.
   *
   * @warning The returned reference is valid only while this `SparseMatrix` instance lives.
   * @par Pitfall
   *      Do not cast away constness or uncompress the returned matrix.
   */
  const Eigen::SparseMatrix<Scalar, Eigen::RowMajor, Index>& Eigen() const;

private:
#ifdef LSQUANT_NATIVE_BUILD
  void PlaceOnThreads();
#endif
  Eigen::SparseMatrix<Scalar, Eigen::RowMajor, Index> matrix_;
};

/**
 * @brief Non-owning `LinearOperator` adapter over a `SparseMatrix`.
 *
 * Bridges a stored `SparseMatrix` into the abstract `LinearOperator` hierarchy without
 * copying matrix data. Caches Hermiticity at construction time in `O(dim + nnz)` so subsequent
 * `IsHermitian()` queries run in `O(1)`.
 *
 * @warning The referenced `SparseMatrix` must outlive this `SparseMatrixOperator`.
 * @par Pitfall
 *      Binding a `SparseMatrixOperator` to a temporary `SparseMatrix` leaves a dangling pointer.
 */
class SparseMatrixOperator : public LinearOperator {
public:
  /**
   * @brief Construct a non-owning operator view of `matrix`.
   *
   * @param matrix Referenced sparse matrix (must outlive `*this`).
   *
   * @warning Evaluates `matrix.IsHermitian()` once during construction.
   * @par Pitfall
   *      Do not destroy or move-from `matrix` while this operator is registered in an engine.
   */
  explicit SparseMatrixOperator(const SparseMatrix& matrix);

  /**
   * @brief Return the dimension of the underlying sparse matrix.
   *
   * @return Square matrix dimension `N >= 0`.
   */
  Index Dimension() const override;

  /**
   * @brief Forward `y = a * M * x + b * y` to `SparseMatrix::Multiply`.
   *
   * @param a Multiplicative scalar for `M * x`.
   * @param x Input vector of length `Dimension()`.
   * @param b Multiplicative scalar for `y`.
   * @param y Input/output vector of length `Dimension()`.
   */
  void Multiply(Scalar a, const Vector& x, Scalar b, Vector& y) const override;
  void MultiplyBatch(Scalar a, const std::vector<Vector>& x, Scalar b,
                     std::vector<Vector>& y) const override;

  using LinearOperator::Multiply;

  /**
   * @brief Return cached Hermiticity of the underlying sparse matrix.
   *
   * @return `true` if the underlying matrix is Hermitian within `1e-12`, `false` otherwise.
   */
  bool IsHermitian() const override;

  /**
   * @brief Access the referenced `SparseMatrix`.
   *
   * @return Const reference to the underlying `SparseMatrix`.
   */
  const SparseMatrix& Matrix() const;

private:
  const SparseMatrix* matrix_;
  bool isHermitian_;
};

}  // namespace lsquant

#endif  // LSQUANT_CORE_ALGEBRA_SPARSE_MATRIX_HPP_
