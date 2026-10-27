/**
 * @file chebmom_file.hpp
 * @brief Readers and writers for the unified `.chebmom` format and legacy `.chebmom1D`/`.chebmom2D` files.
 *
 * Implements the rank-generic `.chebmom` v1 format (with full metadata and 17-digit precision)
 * alongside byte-exact writers for legacy 1D and 2D moment outputs,
 * plus an auto-detecting reader that loads all formats into `ChebyshevMoments`.
 */
#ifndef LSQUANT_IO_CHEBMOM_FILE_HPP_
#define LSQUANT_IO_CHEBMOM_FILE_HPP_

#include <string>

#include "core/chebyshev/chebyshev_moments.hpp"

namespace lsquant {

/**
 * @brief Write a `ChebyshevMoments` tensor in the unified `.chebmom` version 1 text format.
 *
 * Writes the structured header (`chebmom 1`, `rank`, `dims`, `window`, `system_size`,
 * `states`, `word`, `kpm_weighted`, `symmetrised`, `data`) followed by row-major `re im`
 * pairs formatted at `precision(17)`. Complexity is `O(Data().size())`.
 *
 * @param filename Destination file path.
 * @param moments  Moment tensor and metadata to serialise.
 *
 * @warning Throws `std::runtime_error` if `filename` cannot be opened or writing fails.
 * @par Pitfall
 *      Uses 17-digit precision for lossless IEEE-754 double round-tripping; use
 *      `WriteLegacyChebmom1D` or `WriteLegacyChebmom2D` when comparing with `digits10`
 *      legacy goldens.
 */
void WriteChebmom(const std::string& filename, const ChebyshevMoments& moments);

/**
 * @brief Read a moment tensor from a `.chebmom`, `.chebmom1D`, or `.chebmom2D` file.
 *
 * Auto-detects the file format from the first token: if the file begins with `"chebmom"`,
 * parses the version-1 `.chebmom` header and data; otherwise parses the legacy header
 * `systemSize bandwidth bandcenter` followed by one dimension (`.chebmom1D`) or two
 * dimensions (`.chebmom2D`) and row-major `re im` pairs, setting `kpmWeighted = true`,
 * `symmetrised = (rank == 2)`, and `stateLabel` from the filename token after `_state`
 * when present. Complexity is `O(Data().size())`.
 *
 * @param filename Path to the moment file.
 * @return         Populated `ChebyshevMoments` instance.
 *
 * @warning Throws `std::runtime_error` if the file cannot be opened or has an invalid header
 *          or truncated data section.
 * @par Pitfall
 *      Legacy `.chebmom1D` and `.chebmom2D` files do not store `numberOfStates` or `word`;
 *      those fields remain at their defaults when reading legacy files.
 */
ChebyshevMoments ReadChebmom(const std::string& filename);

/**
 * @brief Write a rank-1 `ChebyshevMoments` tensor in the byte-exact legacy `.chebmom1D` layout.
 *
 * Reproduces the legacy 1D moment writer verbatim: sets stream precision to
 * `std::numeric_limits<double>::digits10`, writes `systemSize " " bandwidth " " center "\n"`,
 * writes `Dims()[0] "\n"`, and writes each moment as `re " " im "\n"`.
 *
 * @param filename Destination `.chebmom1D` file path.
 * @param moments  Rank-1 `ChebyshevMoments` tensor.
 *
 * @warning Throws `std::invalid_argument` if `moments.Rank() != 1` and `std::runtime_error`
 *          if the file cannot be opened.
 * @par Pitfall
 *      `bandwidth` and `center` are computed as `(upper - lower) * 1.0` and
 *      `(upper + lower) * 0.5` to match the legacy header convention.
 */
void WriteLegacyChebmom1D(const std::string& filename, const ChebyshevMoments& moments);

/**
 * @brief Write a rank-2 `ChebyshevMoments` tensor in the byte-exact legacy `.chebmom2D` layout.
 *
 * Reproduces the legacy 2D moment writer verbatim: sets stream precision to
 * `std::numeric_limits<double>::digits10`, writes `systemSize " " bandwidth "  " center "\n"`
 * (with two spaces before `center`), writes `Dims()[0] " " Dims()[1] " \n"` (with a trailing
 * space before the newline), and writes each row-major moment as `re " " im "\n"`.
 *
 * @param filename Destination `.chebmom2D` file path.
 * @param moments  Rank-2 `ChebyshevMoments` tensor.
 *
 * @warning Throws `std::invalid_argument` if `moments.Rank() != 2` and `std::runtime_error`
 *          if the file cannot be opened.
 * @par Pitfall
 *      Removing the double space before `center` or the trailing space on the dimensions line
 *      breaks byte-for-byte `cmp` against the legacy `.chebmom2D` goldens.
 */
void WriteLegacyChebmom2D(const std::string& filename, const ChebyshevMoments& moments);

/**
 * @brief Write a rank-2 time-dependent `ChebyshevMoments` tensor in the byte-exact legacy `.chebmomTD` layout.
 *
 * Reproduces the legacy time-dependent moment writer verbatim: sets stream precision to
 * `std::numeric_limits<double>::digits10`, writes
 * `systemSize " " bandwidth " " center " " maxTimeStep " " timeDiff " \n"` (with a trailing
 * space before the newline), writes `numMoms " " maxTimeStep " \n"` (with a trailing space
 * before the newline), and writes each moment in `m`-major (`m` outer, `n` inner) order as
 * `re " " im "\n"`.
 *
 * @param filename Destination `.chebmomTD` file path.
 * @param moments  Rank-2 `(numMoms, numTimes)` `ChebyshevMoments` tensor.
 *
 * @warning Throws `std::invalid_argument` if `moments.Rank() != 2` and `std::runtime_error`
 *          if the file cannot be opened.
 * @par Pitfall
 *      Both header lines carry a single trailing space before the newline in the legacy
 *      writer; omitting either breaks byte-for-byte comparison with G6/G7 goldens.
 */
void WriteLegacyChebmomTD(const std::string& filename, const ChebyshevMoments& moments);

namespace io {
using lsquant::ReadChebmom;
using lsquant::WriteChebmom;
using lsquant::WriteLegacyChebmom1D;
using lsquant::WriteLegacyChebmom2D;
using lsquant::WriteLegacyChebmomTD;
}  // namespace io

}  // namespace lsquant

#endif  // LSQUANT_IO_CHEBMOM_FILE_HPP_
