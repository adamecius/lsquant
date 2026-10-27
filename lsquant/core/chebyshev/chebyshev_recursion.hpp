/**
 * @file chebyshev_recursion.hpp
 * @brief Single-vector three-term Chebyshev polynomial recursion engine.
 *
 * Generates the sequence `T_m(H~) |phi>` via `T_0 = phi`, `T_1 = H~ phi`, and
 * `T_{m+1} = 2 H~ T_m - T_{m-1}` using two rolling state buffers, reproducing the exact
 * call sequence of legacy `Moments::SetInitVectors` and `Moments::Iterate`.
 */
#ifndef LSQUANT_CORE_CHEBYSHEV_CHEBYSHEV_RECURSION_HPP_
#define LSQUANT_CORE_CHEBYSHEV_CHEBYSHEV_RECURSION_HPP_

#include "core/algebra/linear_operator.hpp"
#include "core/algebra/types.hpp"

namespace lsquant {

/**
 * @brief Stateful three-term Chebyshev polynomial recursion `T_{m+1}(H~) = 2 H~ T_m(H~) - T_{m-1}(H~)`.
 *
 * Solves the problem of stepping through Chebyshev polynomials of a rescaled operator `H~`
 * acting on an initial state `|phi>` in `O(Dimension())` memory without storing the full
 * history of vectors. `Start(phi)` initialises `T_0 = phi` and precomputes `T_1 = H~ phi`
 * (matching legacy `SetInitVectors`), setting `Order() == 0` and `Current() == T_0`.
 * Each subsequent `Step()` advances `Order()` by 1 so `Current()` yields `T_{Order()}` and
 * `Previous()` yields `T_{Order()-1}`. Each step costs one `rescaledH.Multiply` call.
 *
 * @warning `Start(phi)` must be called before `Step()`, `Current()`, or `Previous()`.
 * @par Pitfall
 *      `rescaledH` must already be mapped into `[-1, 1]`; applying the recursion to an
 *      unrescaled operator overflows within a few dozen steps.
 */
class ChebyshevRecursion {
public:
  /**
   * @brief Bind a Chebyshev recursion to a rescaled operator `rescaledH`.
   *
   * @param rescaledH Operator with spectrum inside `[-1, 1]` (must outlive `*this`).
   *
   * @warning Throws `std::invalid_argument` if `rescaledH.Dimension() <= 0`.
   * @par Pitfall
   *      Does not copy `rescaledH`; the caller owns the operator lifetime.
   */
  explicit ChebyshevRecursion(const LinearOperator& rescaledH);

  /**
   * @brief Initialise the recursion from state `phi` (`T_0 = phi`, `T_1 = H~ phi`).
   *
   * Copies `phi` into the `T_0` buffer and immediately evaluates `T_1 = H~ T_0` into the
   * lookahead buffer, matching legacy `Moments::SetInitVectors`. Sets `Order() == 0`.
   *
   * @param phi Initial vector of length `rescaledH.Dimension()`.
   *
   * @warning Throws `std::invalid_argument` if `phi.size() != rescaledH.Dimension()`.
   * @par Pitfall
   *      After `Start(phi)`, `Order()` is `0` and `Current()` is `T_0`; `T_1` becomes
   *      `Current()` on the first `Step()`.
   */
  void Start(const Vector& phi);

  /**
   * @brief Advance the recursion by one polynomial order (`m -> m + 1`).
   *
   * Uses the precomputed `T_1` on the transition `0 -> 1`, and for `m >= 1` evaluates
   * `rescaledH.Multiply(2.0, T_m, -1.0, T_{m-1})` followed by a buffer swap, matching
   * legacy `Moments::Iterate` bit for bit.
   *
   * @warning Throws `std::runtime_error` if called before `Start(phi)`.
   * @par Pitfall
   *      `Step()` invalidates references previously obtained from `Current()` and
   *      `Previous()` by swapping internal vector buffers.
   */
  void Step();

  /**
   * @brief Return the current Chebyshev vector `T_{Order()}(H~) |phi>`.
   *
   * @return Const reference to `T_m |phi>`.
   *
   * @warning Throws `std::runtime_error` if called before `Start(phi)`.
   * @par Pitfall
   *      Copy the contents if you need to retain `T_m |phi>` across a subsequent `Step()`.
   */
  const Vector& Current() const;

  /**
   * @brief Return the previous Chebyshev vector `T_{Order()-1}(H~) |phi>`.
   *
   * @return Const reference to `T_{m-1} |phi>` (at `Order() == 0`, returns `T_0`).
   *
   * @warning Throws `std::runtime_error` if called before `Start(phi)`.
   * @par Pitfall
   *      Only meaningful as `T_{m-1}` once `Order() >= 1`.
   */
  const Vector& Previous() const;

  /**
   * @brief Return the polynomial order `m` of `Current()`.
   *
   * @return Non-negative integer order `m` (`-1` before `Start` is called).
   */
  int Order() const;

private:
  const LinearOperator* rescaledH_;
  Vector current_;
  Vector previous_;
  Vector next_;
  int order_ = -1;
};

namespace chebyshev {
using lsquant::ChebyshevRecursion;
}  // namespace chebyshev

}  // namespace lsquant

#endif  // LSQUANT_CORE_CHEBYSHEV_CHEBYSHEV_RECURSION_HPP_
