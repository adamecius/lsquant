/**
 * @file curve_text.hpp
 * @brief Two-column ASCII curve writer reproducing legacy ostream formatting.
 *
 * Writes paired coordinate and value arrays `(x, y)` to a text file using default
 * `std::ofstream` formatting (`%g`-equivalent, 6 significant digits) and `std::endl`
 * line termination.
 */
#ifndef LSQUANT_IO_CURVE_TEXT_HPP_
#define LSQUANT_IO_CURVE_TEXT_HPP_

#include <string>
#include <vector>

#include "core/algebra/types.hpp"

namespace lsquant {

/**
 * @brief Write two strided arrays of `double` values as a two-column ASCII file.
 *
 * Solves the problem of emitting reconstructed 1D curves that match the legacy C++
 * post-processing output streams byte for byte. Opens a default `std::ofstream` (default
 * precision 6, default general floatfield) and writes `count` lines of the form
 * `x[i * xStride] << " " << y[i * yStride] << std::endl`. Complexity is `O(count)`.
 *
 * @param filename Destination file path.
 * @param x        Pointer to the first coordinate element.
 * @param xStride  Element stride between consecutive entries in `x`.
 * @param y        Pointer to the first value element.
 * @param yStride  Element stride between consecutive entries in `y`.
 * @param count    Number of `(x, y)` pairs to write (`count >= 0`).
 *
 * @warning Throws `std::invalid_argument` if `count < 0` or if `x` or `y` is null when
 *          `count > 0`, and `std::runtime_error` if the file cannot be opened or writing fails.
 * @par Pitfall
 *      Do not set custom stream precision on the output stream here; legacy curve files rely
 *      on the default `std::ofstream` precision of 6 significant digits.
 */
void WriteCurveText(const std::string& filename, const double* x, Index xStride,
                    const double* y, Index yStride, Index count);

/**
 * @brief Write two equal-length `std::vector<double>` arrays as a two-column ASCII file.
 *
 * Forwards contiguous vectors `x` and `y` to the strided `WriteCurveText` writer.
 * Complexity is `O(x.size())`.
 *
 * @param filename Destination file path.
 * @param x        First column values.
 * @param y        Second column values (must have the same length as `x`).
 *
 * @warning Throws `std::invalid_argument` if `x.size() != y.size()`, and `std::runtime_error`
 *          if opening or writing `filename` fails.
 * @par Pitfall
 *      Overwrites any existing file at `filename`.
 */
void WriteCurveText(const std::string& filename, const std::vector<double>& x,
                    const std::vector<double>& y);

namespace io {
using lsquant::WriteCurveText;
}  // namespace io

}  // namespace lsquant

#endif  // LSQUANT_IO_CURVE_TEXT_HPP_
