/**
 * @file csr_text.hpp
 * @brief Reader and writer for LSQUANT ASCII CSR operator files (`.CSR`).
 *
 * Parses and serialises the four-section text CSR format produced by wannier2sparse and
 * legacy LSQUANT input generators: header `dim nnz`, `nnz` pairs of `(real, imag)` values,
 * `nnz` zero-based column indices, and `dim + 1` CSR row pointer offsets.
 */
#ifndef LSQUANT_IO_CSR_TEXT_HPP_
#define LSQUANT_IO_CSR_TEXT_HPP_

#include <string>

#include "core/algebra/sparse_matrix.hpp"

namespace lsquant {

/**
 * @brief Read an ASCII `.CSR` operator file and construct a Hermitian-folded `SparseMatrix`.
 *
 * Parses `dim` and `nnz`, followed by `nnz` complex value pairs `(re, im)`, `nnz` column
 * indices, and `dim + 1` row offsets, assembles the matrix via `SparseMatrix::FromCsr`, and
 * checks the anti-Hermitian norm and applies `SparseMatrix::HermitianFold()` (`0.5 * (M + M^+)`) to match the
 * legacy CSR reader convention (which also normalises `-0.0`
 * to `+0.0`). Complexity is `O(dim + nnz)` in time and memory.
 *
 * Reports `||(M-M^+)/2||_F / ||M||_F` for every operator before folding. Refuses
 * relative values above `1e-10` unless `legacyFold` is explicitly true. The
 * final fold preserves the established rounding of accepted Hermitian inputs.
 *
 * @param filename Path to the `.CSR` text file.
 * @param legacyFold Permit a non-Hermitian input for legacy reproduction only.
 * @return         Compressed `SparseMatrix` built from the file contents and folded with
 *                 `HermitianFold()`.
 *
 * @warning Throws `std::runtime_error` if the file cannot be opened or is truncated/malformed,
 *          and `std::invalid_argument` if CSR invariants fail.
 * @par Pitfall
 *      In the LSQUANT `.CSR` text layout, the values array comes first, followed by column
 *      indices, and the `dim + 1` row pointer array comes last.
 */
SparseMatrix ReadCsrText(const std::string& filename, bool legacyFold = false);

/**
 * @brief Write a `SparseMatrix` to an ASCII `.CSR` operator file.
 *
 * Serialises `Dimension()` and `NonZeros()`, followed by the compressed value pairs at
 * 17-digit precision, the column indices, and the `dim + 1` row pointer offsets. Complexity
 * is `O(dim + nnz)`.
 *
 * @param filename Destination file path.
 * @param matrix   Compressed `SparseMatrix` to write.
 *
 * @warning Throws `std::runtime_error` if `filename` cannot be opened or writing fails.
 * @par Pitfall
 *      Overwrites any existing file at `filename`.
 */
void WriteCsrText(const std::string& filename, const SparseMatrix& matrix);

namespace io {
using lsquant::ReadCsrText;
using lsquant::WriteCsrText;
}  // namespace io

}  // namespace lsquant

#endif  // LSQUANT_IO_CSR_TEXT_HPP_
