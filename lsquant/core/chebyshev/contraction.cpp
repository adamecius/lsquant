/**
 * @file contraction.cpp
 * @brief Implementation of deterministic rank-1 and rank-2 Chebyshev series contractions.
 */
#include "core/chebyshev/contraction.hpp"

#include <stdexcept>
#ifdef LSQUANT_NATIVE_BUILD
#include <Eigen/Core>
#include <algorithm>
#endif

namespace lsquant {

namespace {

#ifdef LSQUANT_NATIVE_BUILD
// Evaluate (mu * B) by GEMM, then the diagonal of A^T * (mu * B). Fixed groups of
// 32 energy points keep every reduction independent of the OpenMP team size.
template <typename A, typename B>
void GemmContract(const Scalar* moments, Index m0, Index m1, const A* tableA,
                   const B* tableB, Index points, double* values, bool imaginary) {
  using Dense = Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
  using TableB = Eigen::Matrix<B, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
  Eigen::Map<const Dense> mu(moments, m0, m1);
  Eigen::Map<const TableB> b(tableB, m1, points);
  constexpr Index kPointTileSize = 32;
#pragma omp parallel for schedule(static)
  for (Index first = 0; first < points; first += kPointTileSize) {
    const Index count = std::min(kPointTileSize, points - first);
    Dense partial(m0, count);
    partial.noalias() = mu * b.middleCols(first, count).template cast<Scalar>();
    for (Index j = 0; j < count; ++j) {
      Scalar sum(0.0, 0.0);
      for (Index m = 0; m < m0; ++m) sum += Scalar(tableA[m * points + first + j]) * partial(m, j);
      values[first + j] = imaginary ? sum.imag() : sum.real();
    }
  }
}
#endif

void ValidateRank2Inputs(const char* funcName, const Scalar* moments, Index numMoments0,
                         Index numMoments1, const void* tableA, const void* tableB,
                         Index numPoints, const double* values) {
  if (numMoments0 < 0 || numMoments1 < 0) {
    throw std::invalid_argument(std::string(funcName) + ": moment dimensions must be non-negative");
  }
  if (numPoints < 0) {
    throw std::invalid_argument(std::string(funcName) + ": numPoints must be non-negative");
  }
  if (numMoments0 > 0 && numMoments1 > 0 && moments == nullptr) {
    throw std::invalid_argument(std::string(funcName) + ": moments pointer is null");
  }
  if (numMoments0 > 0 && numPoints > 0 && tableA == nullptr) {
    throw std::invalid_argument(std::string(funcName) + ": tableA pointer is null");
  }
  if (numMoments1 > 0 && numPoints > 0 && tableB == nullptr) {
    throw std::invalid_argument(std::string(funcName) + ": tableB pointer is null");
  }
  if (numPoints > 0 && values == nullptr) {
    throw std::invalid_argument(std::string(funcName) + ": values pointer is null");
  }
}

Index ValidateRank2VectorSizes(const char* funcName, const ChebyshevMoments& moments,
                               std::size_t sizeA, std::size_t sizeB) {
  if (moments.Rank() != 2) {
    throw std::invalid_argument(std::string(funcName) + ": ChebyshevMoments must have rank 2");
  }
  const Index m0 = static_cast<Index>(moments.Dims()[0]);
  const Index m1 = static_cast<Index>(moments.Dims()[1]);
  if (m0 <= 0 || m1 <= 0) {
    throw std::invalid_argument(std::string(funcName) + ": moment dimensions must be positive");
  }
  const Index totalA = static_cast<Index>(sizeA);
  const Index totalB = static_cast<Index>(sizeB);
  if (totalA % m0 != 0 || totalB % m1 != 0) {
    throw std::invalid_argument(std::string(funcName) +
                                ": table sizes must be multiples of moment dimensions");
  }
  const Index numPointsA = totalA / m0;
  const Index numPointsB = totalB / m1;
  if (numPointsA != numPointsB) {
    throw std::invalid_argument(std::string(funcName) +
                                ": tableA and tableB must have the same number of grid points");
  }
  return numPointsA;
}

}  // namespace

void ContractRank1(const double* moments, Index momentStride, Index numMoments,
                   const double* table, Index tableStrideM, Index tableStrideN,
                   Index numPoints, double* values) {
  if (numMoments < 0) {
    throw std::invalid_argument("ContractRank1: numMoments must be non-negative");
  }
  if (numPoints < 0) {
    throw std::invalid_argument("ContractRank1: numPoints must be non-negative");
  }
  if (numMoments > 0 && moments == nullptr) {
    throw std::invalid_argument("ContractRank1: moments pointer is null");
  }
  if (numMoments > 0 && numPoints > 0 && table == nullptr) {
    throw std::invalid_argument("ContractRank1: table pointer is null");
  }
  if (numPoints > 0 && values == nullptr) {
    throw std::invalid_argument("ContractRank1: values pointer is null");
  }

#pragma omp parallel for schedule(static)
  for (Index i = 0; i < numPoints; ++i) {
    double s = 0.0;
    const double* colPtr = table + i * tableStrideN;
    for (Index m = 0; m < numMoments; ++m) {
      s += colPtr[m * tableStrideM] * moments[m * momentStride];
    }
    values[i] = s;
  }
}

std::vector<double> ContractRank1(const std::vector<double>& moments,
                                  const std::vector<double>& table) {
  const Index numMoments = static_cast<Index>(moments.size());
  if (numMoments == 0) {
    if (!table.empty()) {
      throw std::invalid_argument(
          "ContractRank1: table must be empty when moments vector is empty");
    }
    return {};
  }
  const Index totalSize = static_cast<Index>(table.size());
  if (totalSize % numMoments != 0) {
    throw std::invalid_argument(
        "ContractRank1: table size must be a multiple of moments size");
  }
  const Index numPoints = totalSize / numMoments;
  std::vector<double> values(static_cast<std::size_t>(numPoints), 0.0);
  ContractRank1(moments.data(), 1, numMoments, table.data(), numPoints, 1, numPoints,
                values.data());
  return values;
}

std::vector<double> ContractRank1(const ChebyshevMoments& moments,
                                  const std::vector<double>& table) {
  if (moments.Rank() != 1) {
    throw std::invalid_argument("ContractRank1: ChebyshevMoments must have rank 1");
  }
  const Index numMoments = static_cast<Index>(moments.Dims()[0]);
  const Index totalSize = static_cast<Index>(table.size());
  if (numMoments <= 0 || totalSize % numMoments != 0) {
    throw std::invalid_argument(
        "ContractRank1: table size must be a multiple of moments dimension");
  }
  const Index numPoints = totalSize / numMoments;
  std::vector<double> values(static_cast<std::size_t>(numPoints), 0.0);
  const auto* realPtr = reinterpret_cast<const double*>(moments.Data().data());
  ContractRank1(realPtr, 2, numMoments, table.data(), numPoints, 1, numPoints, values.data());
  return values;
}

void ContractRank2Imag(const Scalar* moments, Index numMoments0, Index numMoments1,
                       const double* tableA, const Scalar* tableB, Index numPoints,
                       double* values) {
  ValidateRank2Inputs("ContractRank2Imag", moments, numMoments0, numMoments1, tableA, tableB,
                      numPoints, values);

#ifdef LSQUANT_NATIVE_BUILD
  GemmContract(moments, numMoments0, numMoments1, tableA, tableB, numPoints, values, true);
  return;
#endif

#pragma omp parallel for schedule(static)
  for (Index i = 0; i < numPoints; ++i) {
    double out = 0.0;
    for (Index m0 = 0; m0 < numMoments0; ++m0) {
      const double a = tableA[m0 * numPoints + i];
      const Scalar* muRow = moments + m0 * numMoments1;
      for (Index m1 = 0; m1 < numMoments1; ++m1) {
        out += a * (tableB[m1 * numPoints + i] * muRow[m1]).imag();
      }
    }
    values[i] = out;
  }
}

std::vector<double> ContractRank2Imag(const ChebyshevMoments& moments,
                                      const std::vector<double>& tableA,
                                      const std::vector<Scalar>& tableB) {
  const Index numPoints =
      ValidateRank2VectorSizes("ContractRank2Imag", moments, tableA.size(), tableB.size());
  const Index m0 = static_cast<Index>(moments.Dims()[0]);
  const Index m1 = static_cast<Index>(moments.Dims()[1]);
  std::vector<double> values(static_cast<std::size_t>(numPoints), 0.0);
  ContractRank2Imag(moments.Data().data(), m0, m1, tableA.data(), tableB.data(), numPoints,
                    values.data());
  return values;
}

void ContractRank2Real(const Scalar* moments, Index numMoments0, Index numMoments1,
                       const Scalar* tableA, const Scalar* tableB, Index numPoints,
                       double* values) {
  ValidateRank2Inputs("ContractRank2Real", moments, numMoments0, numMoments1, tableA, tableB,
                      numPoints, values);

#ifdef LSQUANT_NATIVE_BUILD
  GemmContract(moments, numMoments0, numMoments1, tableA, tableB, numPoints, values, false);
  return;
#endif

#pragma omp parallel for schedule(static)
  for (Index i = 0; i < numPoints; ++i) {
    double out = 0.0;
    for (Index m0 = 0; m0 < numMoments0; ++m0) {
      const Scalar a = tableA[m0 * numPoints + i];
      const Scalar* muRow = moments + m0 * numMoments1;
      for (Index m1 = 0; m1 < numMoments1; ++m1) {
        out += ((a * tableB[m1 * numPoints + i]) * muRow[m1]).real();
      }
    }
    values[i] = out;
  }
}

std::vector<double> ContractRank2Real(const ChebyshevMoments& moments,
                                      const std::vector<Scalar>& tableA,
                                      const std::vector<Scalar>& tableB) {
  const Index numPoints =
      ValidateRank2VectorSizes("ContractRank2Real", moments, tableA.size(), tableB.size());
  const Index m0 = static_cast<Index>(moments.Dims()[0]);
  const Index m1 = static_cast<Index>(moments.Dims()[1]);
  std::vector<double> values(static_cast<std::size_t>(numPoints), 0.0);
  ContractRank2Real(moments.Data().data(), m0, m1, tableA.data(), tableB.data(), numPoints,
                    values.data());
  return values;
}

void ContractRank2RealReal(const Scalar* moments, Index numMoments0, Index numMoments1,
                           const double* tableA, const double* tableB, Index numPoints,
                           double* values) {
  ValidateRank2Inputs("ContractRank2RealReal", moments, numMoments0, numMoments1, tableA, tableB,
                      numPoints, values);

#ifdef LSQUANT_NATIVE_BUILD
  GemmContract(moments, numMoments0, numMoments1, tableA, tableB, numPoints, values, false);
  return;
#endif

#pragma omp parallel for schedule(static)
  for (Index i = 0; i < numPoints; ++i) {
    double out = 0.0;
    for (Index m0 = 0; m0 < numMoments0; ++m0) {
      const double a = tableA[m0 * numPoints + i];
      const Scalar* muRow = moments + m0 * numMoments1;
      for (Index m1 = 0; m1 < numMoments1; ++m1) {
        out += (a * tableB[m1 * numPoints + i]) * muRow[m1].real();
      }
    }
    values[i] = out;
  }
}

std::vector<double> ContractRank2RealReal(const ChebyshevMoments& moments,
                                          const std::vector<double>& tableA,
                                          const std::vector<double>& tableB) {
  const Index numPoints =
      ValidateRank2VectorSizes("ContractRank2RealReal", moments, tableA.size(), tableB.size());
  const Index m0 = static_cast<Index>(moments.Dims()[0]);
  const Index m1 = static_cast<Index>(moments.Dims()[1]);
  std::vector<double> values(static_cast<std::size_t>(numPoints), 0.0);
  ContractRank2RealReal(moments.Data().data(), m0, m1, tableA.data(), tableB.data(), numPoints,
                        values.data());
  return values;
}

}  // namespace lsquant
