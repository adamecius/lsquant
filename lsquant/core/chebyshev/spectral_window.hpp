/**
 * @file spectral_window.hpp
 * @brief Spectral energy window and affine map into the Chebyshev interval [-1, 1].
 *
 * Encapsulates the lower and upper energy bounds `[lower, upper]` and derives the
 * half-width, band centre, multiplicative scale, and additive shift that map `[lower, upper]`
 * onto `[-1, 1]` with `CUTOFF = 1.0`.
 */
#ifndef LSQUANT_CORE_CHEBYSHEV_SPECTRAL_WINDOW_HPP_
#define LSQUANT_CORE_CHEBYSHEV_SPECTRAL_WINDOW_HPP_

namespace lsquant {

/**
 * @brief Energy window `[lower, upper]` defining the affine rescale into `[-1, 1]`.
 *
 * Maps an energy or Hamiltonian `H` with spectrum inside `[lower, upper]` to the
 * dimensionless Chebyshev variable `H~ = Scale() * H + Shift() * I` in `[-1, 1]`.
 * Preserves the exact legacy `ScaleFactor()` and `ShiftFactor()` arithmetic with unit
 * cutoff so analytic moment identities hold without distortion. All operations are `O(1)`.
 *
 * @warning `upper` must be strictly greater than `lower`; otherwise `HalfWidth() <= 0` and
 *          `Scale()` divides by zero or flips orientation.
 * @par Pitfall
 *      Never fold reconstruction-grid safety factors into `Scale()` or `Shift()`; moment
 *      rescaling must map `[lower, upper]` onto `[-1, 1]` exactly.
 */
struct SpectralWindow {
  double lower = 0.0;
  double upper = 0.0;

  /**
   * @brief Return the half-bandwidth `(upper - lower) / 2` in energy units.
   *
   * @return Half-width `a = 0.5 * (upper - lower)`.
   */
  double HalfWidth() const { return 0.5 * (upper - lower); }

  /**
   * @brief Return the band centre `(upper + lower) / 2` in energy units.
   *
   * @return Centre `b = 0.5 * (upper + lower)`.
   */
  double Center() const { return 0.5 * (upper + lower); }

  /**
   * @brief Return the multiplicative scale factor `1 / HalfWidth()` (inverse energy units).
   *
   * @return Scale factor mapping `[lower, upper]` to width `2`.
   */
  double Scale() const { return 1.0 / HalfWidth(); }

  /**
   * @brief Return the dimensionless additive shift `-Center() / HalfWidth()`.
   *
   * @return Shift factor centring `[lower, upper]` at `0`.
   */
  double Shift() const { return -Center() / HalfWidth(); }
};

namespace chebyshev {
using lsquant::SpectralWindow;
}  // namespace chebyshev

}  // namespace lsquant

#endif  // LSQUANT_CORE_CHEBYSHEV_SPECTRAL_WINDOW_HPP_
