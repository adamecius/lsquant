/**
 * @file sparse_matrix.cpp
 * @brief Implementation of SparseMatrix and SparseMatrixOperator.
 */
#include "core/algebra/sparse_matrix.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include "core/algebra/memory_placement.hpp"
#ifdef LSQUANT_NATIVE_BUILD
#include "core/algebra/work_partition.hpp"
#endif

namespace lsquant {

SparseMatrix SparseMatrix::FromCsr(Index dim, Index rowPointersSize, const Index* rowPointers,
                                   Index columnsSize, const Index* columns, Index valuesSize,
                                   const Scalar* values) {
  if (dim < 0) {
    throw std::invalid_argument("SparseMatrix::FromCsr: negative dimension");
  }
  if (rowPointersSize != dim + 1) {
    throw std::invalid_argument("SparseMatrix::FromCsr: rowPointers size must equal dim + 1");
  }
  if (columnsSize != valuesSize || columnsSize < 0) {
    throw std::invalid_argument("SparseMatrix::FromCsr: columns and values size mismatch");
  }
  if (rowPointers == nullptr) {
    throw std::invalid_argument("SparseMatrix::FromCsr: rowPointers pointer is null");
  }
  const Index nnz = valuesSize;
  if (nnz > 0 && (columns == nullptr || values == nullptr)) {
    throw std::invalid_argument("SparseMatrix::FromCsr: columns or values pointer is null");
  }
  if (rowPointers[0] != 0 || rowPointers[dim] != nnz) {
    throw std::invalid_argument("SparseMatrix::FromCsr: invalid rowPointers endpoints");
  }
  for (Index r = 0; r < dim; ++r) {
    const Index start = rowPointers[r];
    const Index end = rowPointers[r + 1];
    if (start > end) {
      throw std::invalid_argument("SparseMatrix::FromCsr: rowPointers must be non-decreasing");
    }
    for (Index p = start; p < end; ++p) {
      const Index col = columns[p];
      if (col < 0 || col >= dim) {
        throw std::invalid_argument("SparseMatrix::FromCsr: column index out of range");
      }
    }
  }

  SparseMatrix result;
  if (dim == 0) {
    return result;
  }
  if (nnz == 0) {
    result.matrix_.resize(dim, dim);
    result.matrix_.makeCompressed();
    return result;
  }

  result.matrix_ = Eigen::Map<const Eigen::SparseMatrix<Scalar, Eigen::RowMajor, Index>>(
      dim, dim, nnz, rowPointers, columns, values);

  if (!result.matrix_.isCompressed()) {
    result.matrix_.makeCompressed();
  }
#ifdef LSQUANT_NATIVE_BUILD
  result.PlaceOnThreads();
#endif
  return result;
}

SparseMatrix SparseMatrix::FromCsr(Index dim, std::vector<Index> rowPointers,
                                   std::vector<Index> columns, std::vector<Scalar> values) {
  return FromCsr(dim, static_cast<Index>(rowPointers.size()), rowPointers.data(),
                 static_cast<Index>(columns.size()), columns.data(),
                 static_cast<Index>(values.size()), values.data());
}

SparseMatrix SparseMatrix::FromTriplets(Index dim, const std::vector<Index>& rows,
                                        const std::vector<Index>& cols,
                                        const std::vector<Scalar>& vals) {
  if (dim < 0) {
    throw std::invalid_argument("SparseMatrix::FromTriplets: negative dimension");
  }
  if (rows.size() != cols.size() || cols.size() != vals.size()) {
    throw std::invalid_argument("SparseMatrix::FromTriplets: triplet array size mismatch");
  }

  std::vector<Eigen::Triplet<Scalar, Index>> triplets;
  triplets.reserve(vals.size());
  for (std::size_t k = 0; k < vals.size(); ++k) {
    if (rows[k] < 0 || rows[k] >= dim || cols[k] < 0 || cols[k] >= dim) {
      throw std::invalid_argument("SparseMatrix::FromTriplets: index out of range");
    }
    triplets.emplace_back(rows[k], cols[k], vals[k]);
  }

  SparseMatrix result;
  result.matrix_.resize(dim, dim);
  result.matrix_.setFromTriplets(triplets.begin(), triplets.end());

  if (!result.matrix_.isCompressed()) {
    result.matrix_.makeCompressed();
  }
#ifdef LSQUANT_NATIVE_BUILD
  result.PlaceOnThreads();
#endif
  return result;
}

Index SparseMatrix::Dimension() const { return matrix_.rows(); }

Index SparseMatrix::NonZeros() const { return matrix_.nonZeros(); }

void SparseMatrix::Multiply(Scalar a, const Vector& x, Scalar b, Vector& y) const {
  const Index n = Dimension();
  if (static_cast<Index>(x.size()) != n || static_cast<Index>(y.size()) != n) {
    throw std::invalid_argument("SparseMatrix::Multiply: vector dimension mismatch");
  }
  if (&x == &y) {
    throw std::invalid_argument("SparseMatrix::Multiply: input and output vectors must not alias");
  }
  if (n == 0) {
    return;
  }

  const Index* const rp = matrix_.outerIndexPtr();
  const Index* const cl = matrix_.innerIndexPtr();
  const Scalar* const vl = matrix_.valuePtr();
  const Scalar* const xp = x.data();
  Scalar* const yp = y.data();

#ifdef LSQUANT_NATIVE_BUILD
#pragma omp parallel if(NonZeros() >= SpmvParallelThreshold())
  {
    const Index first = RowBoundary(rp, n, omp_get_thread_num(), omp_get_num_threads());
    const Index last = RowBoundary(rp, n, omp_get_thread_num() + 1, omp_get_num_threads());
    for (Index i = first; i < last; ++i) {
      Scalar acc(0.0, 0.0);
      for (Index p = rp[i]; p < rp[i + 1]; ++p) acc += vl[p] * xp[cl[p]];
      yp[i] = a * acc + b * yp[i];
    }
  }
#else
#pragma omp parallel for schedule(static)
  for (Index i = 0; i < n; ++i) {
    Scalar acc(0.0, 0.0);
    for (Index p = rp[i]; p < rp[i + 1]; ++p) {
      acc += vl[p] * xp[cl[p]];
    }
    yp[i] = a * acc + b * yp[i];
  }
#endif
}

void SparseMatrix::MultiplyBatch(Scalar a, const std::vector<Vector>& x, Scalar b,
                                std::vector<Vector>& y) const {
  if (x.size() != y.size() || &x == &y) {
    throw std::invalid_argument("SparseMatrix::MultiplyBatch: mismatched or aliased batch");
  }
  for (std::size_t j = 0; j < x.size(); ++j) {
    if (static_cast<Index>(x[j].size()) != Dimension() || x[j].size() != y[j].size()) {
      throw std::invalid_argument("SparseMatrix::MultiplyBatch: dimension mismatch");
    }
  }
#ifdef LSQUANT_NATIVE_BUILD
  const Index n = Dimension();
  if (n == 0 || x.empty()) return;
  const Index* rp = matrix_.outerIndexPtr();
  const Index* cl = matrix_.innerIndexPtr();
  const Scalar* vl = matrix_.valuePtr();
  const int width = CacheBatchSize(n, static_cast<int>(x.size()));
#pragma omp parallel if(NonZeros() >= SpmvParallelThreshold())
  {
    const Index first = RowBoundary(rp, n, omp_get_thread_num(), omp_get_num_threads());
    const Index last = RowBoundary(rp, n, omp_get_thread_num() + 1, omp_get_num_threads());
    Vector sums(static_cast<std::size_t>(width));
    for (std::size_t base = 0; base < x.size(); base += static_cast<std::size_t>(width)) {
      const std::size_t count = std::min<std::size_t>(width, x.size() - base);
      for (Index row = first; row < last; ++row) {
        std::fill(sums.begin(), sums.end(), Scalar(0.0, 0.0));
        for (Index p = rp[row]; p < rp[row + 1]; ++p) {
          const Scalar value = vl[p];
          const Index column = cl[p];
          for (std::size_t j = 0; j < count; ++j) sums[j] += value * x[base + j][column];
        }
        for (std::size_t j = 0; j < count; ++j) {
          y[base + j][row] = a * sums[j] + b * y[base + j][row];
        }
      }
    }
  }
#else
  for (std::size_t j = 0; j < x.size(); ++j) Multiply(a, x[j], b, y[j]);
#endif
}

SparseMatrix SparseMatrix::Rescaled(Scalar a, Scalar b) const {
  const Index n = Dimension();
  Eigen::SparseMatrix<Scalar, Eigen::RowMajor, Index> bID(n, n);
  bID.setIdentity();
  bID *= b;

  SparseMatrix result;
  result.matrix_ = a * matrix_ + bID;
  if (!result.matrix_.isCompressed()) {
    result.matrix_.makeCompressed();
  }
#ifdef LSQUANT_NATIVE_BUILD
  result.PlaceOnThreads();
#endif
  return result;
}

SparseMatrix SparseMatrix::Adjoint() const {
  SparseMatrix result;
  result.matrix_ = Eigen::SparseMatrix<Scalar, Eigen::RowMajor, Index>(matrix_.adjoint());
  if (!result.matrix_.isCompressed()) {
    result.matrix_.makeCompressed();
  }
#ifdef LSQUANT_NATIVE_BUILD
  result.PlaceOnThreads();
#endif
  return result;
}

SparseMatrix SparseMatrix::HermitianFold() const {
  SparseMatrix result;
  if (Dimension() == 0) {
    return result;
  }
  Eigen::SparseMatrix<Scalar, Eigen::RowMajor, Index> transp(matrix_.adjoint());
  result.matrix_ = 0.5 * (transp + matrix_);
  if (!result.matrix_.isCompressed()) {
    result.matrix_.makeCompressed();
  }
#ifdef LSQUANT_NATIVE_BUILD
  result.PlaceOnThreads();
#endif
  return result;
}

bool SparseMatrix::IsHermitian(double tol) const {
  if (tol < 0.0) {
    throw std::invalid_argument("SparseMatrix::IsHermitian: negative tolerance");
  }
  if (Dimension() == 0) {
    return true;
  }
  Eigen::SparseMatrix<Scalar, Eigen::RowMajor, Index> adj(matrix_.adjoint());
  const double diff = (matrix_ - adj).norm();
  const double norm = matrix_.norm();
  if (norm == 0.0) {
    return diff <= tol;
  }
  return (diff / norm) <= tol;
}

const Eigen::SparseMatrix<Scalar, Eigen::RowMajor, Index>& SparseMatrix::Eigen() const {
  return matrix_;
}

#ifdef LSQUANT_NATIVE_BUILD
void SparseMatrix::PlaceOnThreads() {
  const Index n = Dimension();
  if (n == 0) return;
  const Index nnz = NonZeros();
  Eigen::SparseMatrix<Scalar, Eigen::RowMajor, Index> placed(n, n);
  placed.makeCompressed();
  placed.resizeNonZeros(nnz);
  DiscardInteriorPages(placed.valuePtr(), static_cast<std::size_t>(nnz) * sizeof(Scalar));
  DiscardInteriorPages(placed.innerIndexPtr(), static_cast<std::size_t>(nnz) * sizeof(Index));
  DiscardInteriorPages(placed.outerIndexPtr(), static_cast<std::size_t>(n + 1) * sizeof(Index));
#pragma omp parallel if(nnz >= SpmvParallelThreshold())
  {
    const Index first = RowBoundary(matrix_.outerIndexPtr(), n, omp_get_thread_num(), omp_get_num_threads());
    const Index last = RowBoundary(matrix_.outerIndexPtr(), n, omp_get_thread_num() + 1, omp_get_num_threads());
    for (Index row = first; row < last; ++row) {
      placed.outerIndexPtr()[row] = matrix_.outerIndexPtr()[row];
      for (Index p = matrix_.outerIndexPtr()[row]; p < matrix_.outerIndexPtr()[row + 1]; ++p) {
        placed.innerIndexPtr()[p] = matrix_.innerIndexPtr()[p];
        placed.valuePtr()[p] = matrix_.valuePtr()[p];
      }
    }
  }
  placed.outerIndexPtr()[n] = nnz;
  matrix_.swap(placed);
}
#endif

SparseMatrixOperator::SparseMatrixOperator(const SparseMatrix& matrix)
    : matrix_(&matrix), isHermitian_(matrix.IsHermitian()) {}

Index SparseMatrixOperator::Dimension() const { return matrix_->Dimension(); }

void SparseMatrixOperator::Multiply(Scalar a, const Vector& x, Scalar b, Vector& y) const {
  matrix_->Multiply(a, x, b, y);
}

void SparseMatrixOperator::MultiplyBatch(Scalar a, const std::vector<Vector>& x, Scalar b,
                                        std::vector<Vector>& y) const {
  matrix_->MultiplyBatch(a, x, b, y);
}

bool SparseMatrixOperator::IsHermitian() const { return isHermitian_; }

const SparseMatrix& SparseMatrixOperator::Matrix() const { return *matrix_; }

}  // namespace lsquant
