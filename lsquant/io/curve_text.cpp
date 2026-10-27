/**
 * @file curve_text.cpp
 * @brief Implementation of the two-column ASCII curve writer.
 */
#include "io/curve_text.hpp"

#include <fstream>
#include <stdexcept>

namespace lsquant {

void WriteCurveText(const std::string& filename, const double* x, Index xStride,
                    const double* y, Index yStride, Index count) {
  if (count < 0) {
    throw std::invalid_argument("WriteCurveText: count must be non-negative");
  }
  if (count > 0 && (x == nullptr || y == nullptr)) {
    throw std::invalid_argument("WriteCurveText: x and y pointers must not be null");
  }

  std::ofstream outputfile(filename.c_str());
  if (!outputfile.is_open()) {
    throw std::runtime_error("WriteCurveText: cannot open file for writing: " + filename);
  }

  for (Index i = 0; i < count; ++i) {
    outputfile << x[i * xStride] << " " << y[i * yStride] << std::endl;
  }

  if (!outputfile) {
    throw std::runtime_error("WriteCurveText: error writing file: " + filename);
  }
  outputfile.close();
}

void WriteCurveText(const std::string& filename, const std::vector<double>& x,
                    const std::vector<double>& y) {
  if (x.size() != y.size()) {
    throw std::invalid_argument("WriteCurveText: x and y must have the same size");
  }
  WriteCurveText(filename, x.data(), 1, y.data(), 1, static_cast<Index>(x.size()));
}

}  // namespace lsquant
