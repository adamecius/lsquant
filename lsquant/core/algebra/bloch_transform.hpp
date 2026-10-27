/**
 * @file bloch_transform.hpp
 * @brief Unitary real-space <-> Bloch-space transform via batched 3D FFTW.
 *
 * Implements the unitary lattice Fourier transform `B` (`RealToBloch`) and its
 * conjugate transpose `B^+` (`BlochToReal`) for an `(N_1, N_2, N_3)` supercell
 * with `W` orbitals per unit cell and arbitrary intra-cell atomic positions
 * `tau_alpha` encoded in a precomputed phase table.
 */
#ifndef LSQUANT_CORE_ALGEBRA_BLOCH_TRANSFORM_HPP_
#define LSQUANT_CORE_ALGEBRA_BLOCH_TRANSFORM_HPP_

#include <memory>
#include <vector>

#include "core/algebra/linear_operator.hpp"
#include "core/algebra/types.hpp"

namespace lsquant {

/**
 * @brief Direction of the unitary lattice Fourier transform.
 */
enum class BlochDirection {
  RealToBloch,  ///< Forward transform `B`: real space `(R, alpha)` -> Bloch `(k, alpha)`.
  BlochToReal   ///< Adjoint transform `B^+`: Bloch `(k, alpha)` -> real space `(R, alpha)`.
};

/**
 * @brief Unitary transform `B` (`RealToBloch`) or `B^+` (`BlochToReal`) using FFTW.
 *
 * Solves the problem of mapping state vectors between real-space supercell
 * coordinates and Bloch k-space coordinates in `O(N_k * W * log N_k)` time
 * without assembling dense `N x N` Fourier matrices. Assumes cell-major vector
 * layout `x[(iR * W) + alpha]` with `n_1`-fastest flat cell index
 * `iR = n_1 + N_1 * n_2 + N_1 * N_2 * n_3` (matching `wannier2sparse`) on the
 * real side and `ik = m_1 + N_1 * m_2 + N_1 * N_2 * m_3` on the Bloch side.
 *
 * For `RealToBloch`, executes `FFTW_FORWARD` (`sum_R exp(-i k . R)`) on the
 * `{N_3, N_2, N_1}` FFTW grid (`n_1` fastest) with `howmany = W`, stride `W`,
 * and multiplies each `(ik, alpha)` component by `exp(-i k . tau_alpha) / sqrt(N_k)`.
 * For `BlochToReal`, multiplies each `(ik, alpha)` component by the conjugate
 * factor `exp(+i k . tau_alpha) / sqrt(N_k)` and executes `FFTW_BACKWARD`
 * (`sum_k exp(+i k . R)`).
 *
 * `Multiply` uses a per-object FFTW workspace (`inBuffer`, `outBuffer`, `plan`)
 * and is not safe to call concurrently on the same `BlochTransform` object
 * (the engine never does; copies constructed via the copy constructor or
 * `Adjoint()` get their own workspace).
 *
 * @warning Throws `std::invalid_argument` if grid dimensions or orbital count
 *          are not positive, or if `atomPhases` does not have length `0` or
 *          `N_1 * N_2 * N_3 * W`.
 * @par Pitfall
 *      Both the real-space CSR indexing and the k-mesh indexing must follow
 *      the `n_1`-fastest cell order `(n_1 + N_1 * n_2 + N_1 * N_2 * n_3) * W + alpha`
 *      matching `wannier2sparse` and FFTW dimension array `{N_3, N_2, N_1}`.
 */
class BlochTransform : public LinearOperator {
public:
  using Direction = BlochDirection;

  /**
   * @brief Tag type indicating that a shared phase table already carries `1 / sqrt(N_k)`.
   */
  struct ScaledPhasesTag {};

  /**
   * @brief Construct a `BlochTransform` for an `(n1, n2, n3)` supercell with `numOrbitals` orbitals.
   *
   * @param n1          Supercell grid count along lattice vector 1 (`n1 > 0`).
   * @param n2          Supercell grid count along lattice vector 2 (`n2 > 0`).
   * @param n3          Supercell grid count along lattice vector 3 (`n3 > 0`).
   * @param numOrbitals Number of orbitals `W > 0` per unit cell.
   * @param atomPhases  Unscaled phase table of length `N_k * W` with entries
   *                    `exp(-i k . tau_alpha)` (as returned by `wannier.bloch_phases`),
   *                    which the constructor always multiplies by `1 / sqrt(N_k)`.
   *                    If empty, defaults to zero intra-cell offsets (`tau_alpha = 0`).
   * @param direction   Transform direction (`RealToBloch` for `B`, `BlochToReal` for `B^+`).
   *
   * @warning Throws `std::invalid_argument` on invalid dimensions or phase table length.
   * @par Pitfall
   *      Always pass unscaled `exp(-i k . tau_alpha)` (with the minus sign) for both
   *      directions; the constructor always scales by `1 / sqrt(N_k)`, and `BlochToReal`
   *      conjugates the shared phase table internally.
   */
  BlochTransform(Index n1, Index n2, Index n3, Index numOrbitals,
                 std::vector<Scalar> atomPhases = {},
                 Direction direction = Direction::RealToBloch);

  /**
   * @brief Construct a `BlochTransform` from a pre-scaled phase table `exp(-i k . tau_alpha) / sqrt(N_k)`.
   *
   * @param n1           Supercell grid count along lattice vector 1 (`n1 > 0`).
   * @param n2           Supercell grid count along lattice vector 2 (`n2 > 0`).
   * @param n3           Supercell grid count along lattice vector 3 (`n3 > 0`).
   * @param numOrbitals  Number of orbitals `W > 0` per unit cell.
   * @param scaledPhases Shared phase table of length `N_k * W` already multiplied by `1 / sqrt(N_k)`.
   * @param direction    Transform direction (`RealToBloch` for `B`, `BlochToReal` for `B^+`).
   * @return             `BlochTransform` instance sharing `scaledPhases` with its own FFTW workspace.
   */
  static BlochTransform FromScaledPhases(
      Index n1, Index n2, Index n3, Index numOrbitals,
      std::shared_ptr<const std::vector<Scalar>> scaledPhases,
      Direction direction = Direction::RealToBloch);

  /**
   * @brief Copy constructor creating dedicated FFTW scratch buffers and sharing the phase table.
   *
   * @param other Source transform.
   */
  BlochTransform(const BlochTransform& other);

  /**
   * @brief Copy assignment operator.
   *
   * @param other Source transform.
   * @return Reference to `*this`.
   */
  BlochTransform& operator=(const BlochTransform& other);

  /**
   * @brief Move constructor.
   */
  BlochTransform(BlochTransform&&) noexcept = default;

  /**
   * @brief Move assignment operator.
   */
  BlochTransform& operator=(BlochTransform&&) noexcept = default;

  /**
   * @brief Destructor.
   */
  ~BlochTransform() override = default;

  /**
   * @brief Return the total Hilbert-space dimension `N_1 * N_2 * N_3 * W`.
   *
   * @return Dimension `N_k * W > 0`.
   */
  Index Dimension() const override { return dimension_; }

  /**
   * @brief Return grid count `N_1`.
   */
  Index N1() const { return n1_; }

  /**
   * @brief Return grid count `N_2`.
   */
  Index N2() const { return n2_; }

  /**
   * @brief Return grid count `N_3`.
   */
  Index N3() const { return n3_; }

  /**
   * @brief Return the total number of cells/k-points `N_k = N_1 * N_2 * N_3`.
   */
  Index NumCells() const { return numCells_; }

  /**
   * @brief Return the number of orbitals per cell `W`.
   */
  Index NumOrbitals() const { return numOrbitals_; }

  /**
   * @brief Return the transform direction (`RealToBloch` or `BlochToReal`).
   */
  Direction GetDirection() const { return direction_; }

  /**
   * @brief Evaluate `y = a * B * x + b * y` (or `a * B^+ * x + b * y`).
   *
   * @param a Multiplicative scalar factor for the transformed vector.
   * @param x Input vector of length `Dimension()`.
   * @param b Multiplicative scalar factor for `y`.
   * @param y Input/output vector of length `Dimension()`.
   */
  void Multiply(Scalar a, const Vector& x, Scalar b, Vector& y) const override;

  using LinearOperator::Multiply;

  /**
   * @brief Report whether the operator is self-adjoint.
   *
   * @return `false` for general Fourier transforms.
   */
  bool IsHermitian() const override { return false; }

  /**
   * @brief Return the conjugate transpose transform (`RealToBloch <-> BlochToReal`)
   *        sharing the precomputed phase table.
   *
   * @return Adjoint `BlochTransform` instance.
   */
  BlochTransform Adjoint() const;

private:
  struct FftwWorkspace;

  BlochTransform(ScaledPhasesTag, Index n1, Index n2, Index n3, Index numOrbitals,
                 std::shared_ptr<const std::vector<Scalar>> scaledPhases,
                 Direction direction);

  void InitWorkspace();

  Index n1_ = 0;
  Index n2_ = 0;
  Index n3_ = 0;
  Index numCells_ = 0;
  Index numOrbitals_ = 0;
  Index dimension_ = 0;
  Direction direction_ = Direction::RealToBloch;
  std::shared_ptr<const std::vector<Scalar>> scaledPhases_;
  std::shared_ptr<FftwWorkspace> workspace_;
};

}  // namespace lsquant

#endif  // LSQUANT_CORE_ALGEBRA_BLOCH_TRANSFORM_HPP_
