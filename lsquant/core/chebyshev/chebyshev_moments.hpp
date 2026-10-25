/**
 * @file chebyshev_moments.hpp
 * @brief Rank-generic Chebyshev moment tensor container and KPM weighting routines.
 *
 * Unifies 1D, 2D, time-dependent, and local Chebyshev moment containers into
 * a single rank-generic tensor `ChebyshevMoments` storing row-major
 * complex moments alongside spectral window, trace, and word metadata.
 */
#ifndef LSQUANT_CORE_CHEBYSHEV_CHEBYSHEV_MOMENTS_HPP_
#define LSQUANT_CORE_CHEBYSHEV_CHEBYSHEV_MOMENTS_HPP_

#include <string>
#include <vector>

#include "core/algebra/types.hpp"
#include "core/algebra/build_profile.hpp"
#include "core/chebyshev/spectral_window.hpp"

namespace lsquant {

/**
 * @brief Rank-generic tensor of Chebyshev expansion moments with provenance metadata.
 *
 * Stores `prod_k dims[k]` complex moments in contiguous row-major layout (1 element when
 * `Rank() == 0`). Carries the `SpectralWindow` used to rescale the Hamiltonian, the
 * Hilbert-space dimension `systemSize`, the trace vector count `numberOfStates`, the
 * state generator label `stateLabel`, the evaluated word string `word`, and flags
 * recording whether the `(2 - delta_{m0})` KPM weights (`kpmWeighted`) and Hermitian
 * transpose symmetrisation (`symmetrised`) have been applied.
 *
 * @warning All dimensions in `dims` must be strictly positive.
 * @par Pitfall
 *      Never apply KPM `(2 - delta_{m0})` weights or symmetrisation twice; check
 *      `kpmWeighted` and `symmetrised` before post-processing.
 */
class ChebyshevMoments {
public:
  /**
   * @brief Construct a zero-initialised moment tensor of shape `dims`.
   *
   * @param dims Tensor dimensions per expanded Chebyshev slot (empty for rank 0).
   *
   * @warning Throws `std::invalid_argument` if any dimension `<= 0`.
   * @par Pitfall
   *      For rank 0 (`dims = {}`), `Data()` has size `1`, not `0`.
   */
  explicit ChebyshevMoments(std::vector<int> dims = {});

  /**
   * @brief Return the tensor rank (`dims.size()`).
   *
   * @return Non-negative integer rank.
   */
  int Rank() const;

  /**
   * @brief Return the tensor dimensions along each expanded Chebyshev axis.
   *
   * @return Const reference to the dimension vector.
   */
  const std::vector<int>& Dims() const;

  /**
   * @brief Mutable element access by general multi-index.
   *
   * @param index Vector of `Rank()` zero-based indices.
   * @return      Mutable reference to the corresponding tensor element.
   *
   * @warning Throws `std::invalid_argument` if `index.size() != Rank()` or `std::out_of_range`
   *          if any coordinate is out of bounds.
   * @par Pitfall
   *      Prefer `At(m)` or `At(m, n)` in tight loops to avoid `std::vector` allocation.
   */
  Scalar& At(const std::vector<int>& index);

  /**
   * @brief Const element access by general multi-index.
   *
   * @param index Vector of `Rank()` zero-based indices.
   * @return      Const reference to the corresponding tensor element.
   */
  const Scalar& At(const std::vector<int>& index) const;

  /**
   * @brief Mutable element access for rank-1 moment vectors.
   *
   * @param m Zero-based moment index in `[0, Dims()[0])`.
   * @return  Mutable reference to `mu[m]`.
   *
   * @warning Throws `std::invalid_argument` if `Rank() != 1` or `std::out_of_range` if
   *          `m` is out of bounds.
   * @par Pitfall
   *      Do not call `At(m)` on a rank-2 tensor to access flat storage; use `Data()[idx]`.
   */
  Scalar& At(int m);

  /**
   * @brief Const element access for rank-1 moment vectors.
   *
   * @param m Zero-based moment index in `[0, Dims()[0])`.
   * @return  Const reference to `mu[m]`.
   */
  const Scalar& At(int m) const;

  /**
   * @brief Mutable element access for rank-2 moment matrices (row-major: `m * Dims()[1] + n`).
   *
   * @param m Zero-based first (bra) moment index in `[0, Dims()[0])`.
   * @param n Zero-based second (ket) moment index in `[0, Dims()[1])`.
   * @return  Mutable reference to `mu[m, n]`.
   *
   * @warning Throws `std::invalid_argument` if `Rank() != 2` or `std::out_of_range` if
   *          `m` or `n` is out of bounds.
   * @par Pitfall
   *      Row index `m` corresponds to the bra expansion and column index `n` to the ket
   *      expansion.
   */
  Scalar& At(int m, int n);

  /**
   * @brief Const element access for rank-2 moment matrices.
   *
   * @param m Zero-based first (bra) moment index in `[0, Dims()[0])`.
   * @param n Zero-based second (ket) moment index in `[0, Dims()[1])`.
   * @return  Const reference to `mu[m, n]`.
   */
  const Scalar& At(int m, int n) const;

  /**
   * @brief Mutable access to the contiguous row-major storage buffer.
   *
   * @return Reference to the underlying `std::vector<Scalar>`.
   */
  std::vector<Scalar>& Data();

  /**
   * @brief Const access to the contiguous row-major storage buffer.
   *
   * @return Const reference to the underlying `std::vector<Scalar>`.
   */
  const std::vector<Scalar>& Data() const;

  SpectralWindow window;
  Index systemSize = 0;
  int numberOfStates = 0;
  std::string stateLabel;
  std::string word;
  std::string buildProfile = BuildProfileName();
  bool kpmWeighted = false;
  bool symmetrised = false;
  double timeStep = 0.0;
  int numberOfTimeSteps = 0;

private:
  std::size_t FlatIndex(const std::vector<int>& index) const;

  std::vector<int> dims_;
  std::vector<Scalar> data_;
};

/**
 * @brief Apply the legacy rank-2 KPM `(2 - delta_{m0})` weights and optional Hermitian symmetrisation.
 *
 * Reproduces the post-processing loop of legacy `CorrelationExpansionMoments` verbatim:
 * for `mL in 0..M-1` and `mR in mL..M-1`, computes `scal = 4.0 / numberOfStates` with
 * factors of `0.5` when `mL == 0` and when `mR == 0`, and when `symmetrise == true` sets
 * `tmp = scal * (mu(mL, mR) + conj(mu(mR, mL))) / 2.0`, `mu(mL, mR) = tmp`, and
 * `mu(mR, mL) = conj(tmp)`. Marks `mu.kpmWeighted = true` and `mu.symmetrised = symmetrise`.
 * Complexity is `O(Dims()[0] * Dims()[1])`.
 *
 * @param mu             Rank-2 moment tensor modified in place.
 * @param numberOfStates Number of trace states `N_states > 0` used to accumulate `mu`.
 * @param symmetrise     If `true`, enforces `mu(m, n) = conj(mu(n, m))` on square rank-2 tensors.
 *
 * @warning Throws `std::invalid_argument` if `mu.Rank() != 2`, `numberOfStates <= 0`,
 *          `mu.kpmWeighted` is already `true`, or `symmetrise` is requested on a non-square tensor.
 * @par Pitfall
 *      Must be called at most once after raw rank-2 trace accumulation; rank-1 `dos` plans
 *      fold their weights into `Trace.weights` during accumulation instead.
 */
void ApplyKpmWeights(ChebyshevMoments& mu, int numberOfStates, bool symmetrise);

namespace chebyshev {
using lsquant::ApplyKpmWeights;
using lsquant::ChebyshevMoments;
}  // namespace chebyshev

}  // namespace lsquant

#endif  // LSQUANT_CORE_CHEBYSHEV_CHEBYSHEV_MOMENTS_HPP_
