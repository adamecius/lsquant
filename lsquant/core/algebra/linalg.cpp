/**
 * @file linalg.cpp
 * @brief Implementation of dense vector linear algebra wrappers over Eigen::Map.
 */
#include "core/algebra/linalg.hpp"

#include <Eigen/Core>
#include <algorithm>
#include <stdexcept>
#include "core/algebra/memory_placement.hpp"

namespace lsquant {
namespace linalg {

void axpy(Scalar a, const Vector& x, Vector& y) {
  if (x.size() != y.size()) {
    throw std::invalid_argument("linalg::axpy: vector size mismatch");
  }
  Eigen::Map<const Eigen::Vector<Scalar, -1>> eig_x(x.data(), static_cast<Eigen::Index>(x.size()));
  Eigen::Map<Eigen::Vector<Scalar, -1>> eig_y(y.data(), static_cast<Eigen::Index>(y.size()));
  eig_y = a * eig_x + eig_y;
}

void scal(Scalar a, Vector& x) {
  Eigen::Map<Eigen::Vector<Scalar, -1>> eig_x(x.data(), static_cast<Eigen::Index>(x.size()));
  eig_x *= a;
}

void copy(const Vector& x, Vector& y) {
  if (x.size() != y.size()) {
    throw std::invalid_argument("linalg::copy: vector size mismatch");
  }
  Eigen::Map<const Eigen::Vector<Scalar, -1>> eig_x(x.data(), static_cast<Eigen::Index>(x.size()));
  Eigen::Map<Eigen::Vector<Scalar, -1>> eig_y(y.data(), static_cast<Eigen::Index>(y.size()));
  eig_y = eig_x;
}

Scalar dot(const Vector& x, const Vector& y) {
  if (x.size() != y.size()) {
    throw std::invalid_argument("linalg::dot: vector size mismatch");
  }
  Eigen::Map<const Eigen::Vector<Scalar, -1>> eig_x(x.data(), static_cast<Eigen::Index>(x.size()));
  Eigen::Map<const Eigen::Vector<Scalar, -1>> eig_y(y.data(), static_cast<Eigen::Index>(y.size()));
  return eig_x.dot(eig_y);
}

double norm2(const Vector& x) {
  Eigen::Map<const Eigen::Vector<Scalar, -1>> eig_x(x.data(), static_cast<Eigen::Index>(x.size()));
  return eig_x.squaredNorm();
}

std::vector<Scalar> MomentTable(const std::vector<const Vector*>& left,
                               const std::vector<const Vector*>& right) {
  if (left.empty() || right.empty()) return {};
  if (!left.front()) throw std::invalid_argument("MomentTable: null left vector");
  const Index dimension = static_cast<Index>(left.front()->size());
  for (const Vector* v : left) {
    if (!v || static_cast<Index>(v->size()) != dimension) throw std::invalid_argument("MomentTable: invalid left vector");
  }
  for (const Vector* v : right) {
    if (!v || static_cast<Index>(v->size()) != dimension) throw std::invalid_argument("MomentTable: invalid right vector");
  }
  using Dense = Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>;
  Dense a(dimension, left.size()), b(dimension, right.size());
#ifdef LSQUANT_NATIVE_BUILD
  DiscardInteriorPages(a.data(), static_cast<std::size_t>(a.size()) * sizeof(Scalar));
  DiscardInteriorPages(b.data(), static_cast<std::size_t>(b.size()) * sizeof(Scalar));
#endif
#pragma omp parallel for schedule(static)
  for (Index row = 0; row < dimension; ++row) {
    for (std::size_t j = 0; j < left.size(); ++j) a(row, j) = (*left[j])[row];
    for (std::size_t j = 0; j < right.size(); ++j) b(row, j) = (*right[j])[row];
  }
  using RowDense = Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
  std::vector<Scalar> values(left.size() * right.size());
#ifdef LSQUANT_NATIVE_BUILD
  DiscardInteriorPages(values.data(), values.size() * sizeof(Scalar));
#endif
  Eigen::Map<RowDense> result(values.data(), left.size(), right.size());
  constexpr Index tile = 16;
  const Index rows = static_cast<Index>(left.size()), cols = static_cast<Index>(right.size());
  const Index rowTiles = (rows + tile - 1) / tile, colTiles = (cols + tile - 1) / tile;
#pragma omp parallel for schedule(static)
  for (Index task = 0; task < rowTiles * colTiles; ++task) {
    const Index row = (task / colTiles) * tile, col = (task % colTiles) * tile;
    const Index nr = std::min(tile, rows - row), nc = std::min(tile, cols - col);
    result.block(row, col, nr, nc).noalias() = a.middleCols(row, nr).adjoint() * b.middleCols(col, nc);
  }
  return values;
}

}  // namespace linalg
}  // namespace lsquant
