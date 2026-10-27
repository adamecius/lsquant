/**
 * @file csr_text.cpp
 * @brief Implementation of ASCII CSR file reader and writer.
 */
#include "io/csr_text.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>

namespace lsquant {

SparseMatrix ReadCsrText(const std::string& filename, bool legacyFold) {
  std::ifstream file(filename);
  if (!file.is_open()) {
    throw std::runtime_error("ReadCsrText: cannot open file: " + filename);
  }

  Index dim = 0;
  Index nnz = 0;
  if (!(file >> dim >> nnz) || dim < 0 || nnz < 0) {
    throw std::runtime_error("ReadCsrText: invalid header in file: " + filename);
  }

  std::vector<Scalar> values(static_cast<std::size_t>(nnz));
  for (Index i = 0; i < nnz; ++i) {
    double re = 0.0;
    double im = 0.0;
    if (!(file >> re >> im)) {
      throw std::runtime_error("ReadCsrText: failed reading values in file: " + filename);
    }
    values[static_cast<std::size_t>(i)] = Scalar(re, im);
  }

  std::vector<Index> columns(static_cast<std::size_t>(nnz));
  for (Index i = 0; i < nnz; ++i) {
    Index col = 0;
    if (!(file >> col)) {
      throw std::runtime_error("ReadCsrText: failed reading column indices in file: " + filename);
    }
    columns[static_cast<std::size_t>(i)] = col;
  }

  std::vector<Index> rowPointers(static_cast<std::size_t>(dim + 1));
  for (Index i = 0; i <= dim; ++i) {
    Index rp = 0;
    if (!(file >> rp)) {
      throw std::runtime_error("ReadCsrText: failed reading row pointers in file: " + filename);
    }
    rowPointers[static_cast<std::size_t>(i)] = rp;
  }

  auto matrix = SparseMatrix::FromCsr(dim, std::move(rowPointers), std::move(columns),
                                      std::move(values));
  // Report ||(M-M^+)/2||_F / ||M||_F before the legacy round-off fold.
  const auto& stored = matrix.Eigen();
  Eigen::SparseMatrix<Scalar, Eigen::RowMajor, Index> adjoint = stored.adjoint();
  Eigen::SparseMatrix<Scalar, Eigen::RowMajor, Index> difference = stored - adjoint;
  const double norm = stored.norm();
  const double antiNorm = 0.5 * difference.norm();
  const double relative = antiNorm / std::max(norm, std::numeric_limits<double>::min());
  std::clog << "ReadCsrText: " << filename << " anti_hermitian_norm=" << antiNorm
            << " relative=" << relative << " legacy_fold=" << legacyFold << "\n";
  if (!std::isfinite(norm) || !std::isfinite(relative) ||
      (!legacyFold && relative > 1e-10)) {
    throw std::runtime_error("ReadCsrText: non-Hermitian operator: " + filename +
                             "; use legacy_fold only to reproduce legacy inputs");
  }
  return matrix.HermitianFold();
}

void WriteCsrText(const std::string& filename, const SparseMatrix& matrix) {
  std::ofstream file(filename);
  if (!file.is_open()) {
    throw std::runtime_error("WriteCsrText: cannot open file for writing: " + filename);
  }

  const Index dim = matrix.Dimension();
  const Index nnz = matrix.NonZeros();
  const auto& eig = matrix.Eigen();
  const Index* const rp = eig.outerIndexPtr();
  const Index* const cl = eig.innerIndexPtr();
  const Scalar* const vl = eig.valuePtr();

  file.precision(17);
  file << dim << " " << nnz << "\n";

  for (Index i = 0; i < nnz; ++i) {
    file << vl[i].real() << " " << vl[i].imag() << " ";
  }
  file << "\n";

  for (Index i = 0; i < nnz; ++i) {
    file << cl[i] << " ";
  }
  file << "\n";

  for (Index i = 0; i <= dim; ++i) {
    file << rp[i] << " ";
  }
  file << "\n";

  if (!file) {
    throw std::runtime_error("WriteCsrText: error writing file: " + filename);
  }
}

}  // namespace lsquant
