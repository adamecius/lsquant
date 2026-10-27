/**
 * @file state_generator.hpp
 * @brief Source-state generators for exact and stochastic trace evaluations.
 *
 * Defines the polymorphic `StateGenerator` interface and its concrete implementations:
 * `RandomPhaseStates` (reproducing the legacy Mersenne Twister random-phase convention),
 * `LocalStates` (unit coordinate vectors for local projections or full exact traces),
 * `BlockLocalStates` (per-orbital block-uniform vectors for exact k-space traces in `W` states),
 * `TransformedStates` (applies a unitary `LinearOperator` to an inner state stream),
 * and `UserStates` (caller-provided state vectors).
 */
#ifndef LSQUANT_CORE_ALGEBRA_STATE_GENERATOR_HPP_
#define LSQUANT_CORE_ALGEBRA_STATE_GENERATOR_HPP_

#include <random>
#include <string>
#include <vector>

#include "core/algebra/linear_operator.hpp"
#include "core/algebra/types.hpp"

namespace lsquant {

/**
 * @brief Abstract iterator-style generator yielding quantum state vectors for trace estimation.
 *
 * Solves the problem of streaming source states `|r>` one at a time into the moment engine
 * without storing all `NumberOfStates()` vectors in memory simultaneously. Each call to
 * `Next(out)` populates `out` with a vector of length `Dimension()` and advances the cursor;
 * `Reset()` rewinds the generator to the initial state so multiple plans can be evaluated over
 * the identical state sequence.
 *
 * @warning `Next` mutates the generator's internal state and is not meant to be called
 *          concurrently from multiple threads.
 * @par Pitfall
 *      Always call `Reset()` before reusing a `StateGenerator` across successive `MomentEngine`
 *      runs so that the cursor starts at index 0.
 */
class StateGenerator {
public:
  /**
   * @brief Virtual destructor.
   */
  virtual ~StateGenerator() = default;

  /**
   * @brief Return the Hilbert-space dimension `N` of generated vectors.
   *
   * @return Vector length `N > 0`.
   */
  virtual Index Dimension() const = 0;

  /**
   * @brief Return the total number of state vectors in the sequence.
   *
   * @return Count of states `R >= 0`.
   */
  virtual Index NumberOfStates() const = 0;

  /**
   * @brief Populate `out` with the next state vector in the sequence.
   *
   * @param out Output vector resized/filled with `Dimension()` complex entries.
   * @return    `true` if a state was produced, `false` when the sequence is exhausted.
   */
  virtual bool Next(Vector& out) = 0;

  /**
   * @brief Rewind the generator to its initial state and re-seed any internal PRNG.
   */
  virtual void Reset() = 0;

  /**
   * @brief Return the canonical state-set label used in output filenames and headers.
   *
   * @return One of `"default"`, `"exact"`, `"local"`, or `"user"`.
   */
  virtual std::string Label() const = 0;
};

/**
 * @brief Stochastic-trace generator yielding normalised random-phase vectors.
 *
 * Reproduces the exact legacy random-phase convention (`RANDOM_STATE`):
 * for each component `i in 0..N-1`, draws `u` from `std::uniform_real_distribution<double>(0.0, 1.0)`
 * driven by `std::mt19937(seed)`, sets `phi = 2 * pi * u`, and assigns
 * `x_i = (cos(phi) / sqrt(N), sin(phi) / sqrt(N))`. Both the engine and the distribution
 * object are constructed once per generator instance. Cost per state is `O(N)`.
 *
 * @warning Changing the distribution construction lifetime or PRNG type changes the random
 *          stream and breaks stochastic-trace regression goldens.
 * @par Pitfall
 *      `Reset()` must both re-seed `rng_` with `seed_` and call `uni_.reset()` so that cached
 *      distribution state cannot leak across runs.
 */
class RandomPhaseStates : public StateGenerator {
public:
  /**
   * @brief Construct a random-phase generator for `numberOfStates` vectors of length `dimension`.
   *
   * @param dimension      Hilbert-space dimension `N > 0`.
   * @param numberOfStates Number of random vectors `R >= 0`.
   * @param seed           Deterministic `std::mt19937` seed (default `12345`).
   *
   * @warning Throws `std::invalid_argument` if `dimension <= 0` or `numberOfStates < 0`.
   * @par Pitfall
   *      Pass the same `seed` to reproduce a previous stochastic trace bit for bit.
   */
  RandomPhaseStates(Index dimension, Index numberOfStates, unsigned seed = 12345);

  /**
   * @brief Return the vector dimension `N`.
   *
   * @return Dimension `N > 0`.
   */
  Index Dimension() const override;

  /**
   * @brief Return the number of random vectors `R`.
   *
   * @return State count `R >= 0`.
   */
  Index NumberOfStates() const override;

  /**
   * @brief Generate the next normalised random-phase vector into `out`.
   *
   * @param out Vector overwritten with the next random-phase state.
   * @return    `true` if generated, `false` if all `NumberOfStates()` have been consumed.
   */
  bool Next(Vector& out) override;

  /**
   * @brief Rewind the counter and reset the `std::mt19937` engine to its initial `seed`.
   */
  void Reset() override;

  /**
   * @brief Return `"default"` (matching legacy `gen.StateLabel()`).
   *
   * @return `"default"`.
   */
  std::string Label() const override;

private:
  Index dimension_;
  Index numberOfStates_;
  unsigned seed_;
  Index current_ = 0;
  std::mt19937 rng_;
  std::uniform_real_distribution<double> uni_;
};

/**
 * @brief Deterministic generator yielding unit coordinate vectors `|e_{s}>`.
 *
 * For each index `s` in `sites`, yields a vector of length `dimension` with `1.0` at
 * index `s` and `0.0` elsewhere. Reports label `"exact"` if and only if `sites` is the
 * complete ascending sequence `0, 1, ..., dimension - 1`, and `"local"` otherwise.
 * Cost per state is `O(dimension)`.
 *
 * @warning Throws `std::invalid_argument` if `dimension <= 0` or any site index is outside
 *          `[0, dimension)`.
 * @par Pitfall
 *      Full exact trace over `N` states requires `O(N)` polynomial recursions; use only on
 *      small or medium systems.
 */
class LocalStates : public StateGenerator {
public:
  /**
   * @brief Construct a local-state generator for the specified site indices.
   *
   * @param dimension Hilbert-space dimension `N > 0`.
   * @param sites     Ordered list of site indices in `[0, dimension)`.
   *
   * @warning Throws `std::invalid_argument` if `dimension <= 0` or any site is out of bounds.
   * @par Pitfall
   *      The label is `"exact"` only when `sites` equals `0..dimension-1` in ascending order.
   */
  LocalStates(Index dimension, std::vector<Index> sites);

  /**
   * @brief Return the vector dimension `N`.
   *
   * @return Dimension `N > 0`.
   */
  Index Dimension() const override;

  /**
   * @brief Return `sites.size()`.
   *
   * @return Number of local site states.
   */
  Index NumberOfStates() const override;

  /**
   * @brief Populate `out` with the next unit vector `|e_{sites[current]}>`.
   *
   * @param out Vector overwritten with zeros and `1.0` at `sites[current]`.
   * @return    `true` if a state was written, `false` when `sites` is exhausted.
   */
  bool Next(Vector& out) override;

  /**
   * @brief Rewind the site cursor to `0`.
   */
  void Reset() override;

  /**
   * @brief Return `"exact"` iff `sites == 0..N-1`, else `"local"`.
   *
   * @return `"exact"` or `"local"`.
   */
  std::string Label() const override;

private:
  Index dimension_;
  std::vector<Index> sites_;
  Index current_ = 0;
  std::string label_;
};

/**
 * @brief Deterministic generator yielding `W` block-uniform orbital states for exact block traces.
 *
 * Solves the problem of computing the exact trace `(1 / N_k) sum_k Tr_W[...]` of any
 * block-diagonal operator in `W` state evaluations instead of `N_k * W`. Blocks `ik`
 * are ordered with `m_1` fastest (`ik = m_1 + N_1 * m_2 + N_1 * N_2 * m_3`), matching
 * `wannier2sparse` and `BlochTransform`. For each orbital index `r in 0..W-1`, yields
 * a unit-norm vector of length `N_k * W` with `1 / sqrt(N_k)` at `ik * W + r` for every
 * `ik in 0..N_k-1` and `0` elsewhere. Because `NumberOfStates() == W`, dividing the
 * accumulated trace by `NumberOfStates()` produces `(1 / (N_k * W)) Tr[...]`, identical
 * to a full `N_k * W` coordinate trace.
 *
 * @warning Valid as an exact trace only when all operators in the evaluated word are
 *          block-diagonal with block size `W` (no inter-block mixing).
 * @par Pitfall
 *      Do not use `BlockLocalStates` when real-space disorder couples distinct `k`
 *      blocks; use `RandomPhaseStates` (or `TransformedStates`) in that regime.
 */
class BlockLocalStates : public StateGenerator {
public:
  /**
   * @brief Compute the `n_1`-fastest flat block/cell index `m_1 + N_1 * m_2 + N_1 * N_2 * m_3`.
   */
  static Index BlockIndex(Index m1, Index m2, Index m3, Index n1, Index n2) {
    return m1 + n1 * (m2 + n2 * m3);
  }

  /**
   * @brief Construct a block-uniform state generator for `numBlocks` blocks of size `blockSize`.
   *
   * @param numBlocks Number of diagonal blocks `N_k > 0`.
   * @param blockSize Block dimension (orbital count) `W > 0`.
   *
   * @warning Throws `std::invalid_argument` if `numBlocks <= 0` or `blockSize <= 0`.
   * @par Pitfall
   *      Reports `NumberOfStates() == blockSize` (`W`) and `Dimension() == numBlocks * blockSize`.
   */
  BlockLocalStates(Index numBlocks, Index blockSize);

  /**
   * @brief Construct a block-uniform state generator for an `(n1, n2, n3)` grid with `blockSize` orbitals.
   *
   * @param n1        Grid count `N_1 > 0` (fastest axis).
   * @param n2        Grid count `N_2 > 0` (middle axis).
   * @param n3        Grid count `N_3 > 0` (slowest axis).
   * @param blockSize Block dimension (orbital count) `W > 0`.
   */
  BlockLocalStates(Index n1, Index n2, Index n3, Index blockSize);

  /**
   * @brief Return the total vector dimension `N_k * W`.
   *
   * @return Dimension `N_k * W > 0`.
   */
  Index Dimension() const override;

  /**
   * @brief Return the number of orbital states `W`.
   *
   * @return State count `W > 0`.
   */
  Index NumberOfStates() const override;

  /**
   * @brief Populate `out` with state `r` having `1 / sqrt(N_k)` at `ik * W + r` for all `ik`.
   *
   * @param out Vector of length `N_k * W` overwritten with the `r`-th orbital state.
   * @return    `true` if a state was written, `false` when all `W` states have been yielded.
   */
  bool Next(Vector& out) override;

  /**
   * @brief Rewind the orbital cursor to `0`.
   */
  void Reset() override;

  /**
   * @brief Return `"exact"`.
   *
   * @return `"exact"`.
   */
  std::string Label() const override;

private:
  Index numBlocks_;
  Index blockSize_;
  Index current_ = 0;
};

/**
 * @brief State-generator adapter applying a `LinearOperator` transform to an inner generator.
 *
 * Solves the problem of evaluating stochastic traces in two unitarily equivalent
 * representations (e.g. real space and Bloch space via `BlochTransform`) on the exact
 * same physical state vectors `|r' > = U |r>` so that finite-sample stochastic moments
 * match to machine precision.
 *
 * @warning Both `transform` and `inner` must outlive this `TransformedStates` instance
 *          and have identical `Dimension()`.
 * @par Pitfall
 *      `Reset()` forwards to `inner.Reset()` so repeated runs replay the same transformed
 *      sequence.
 */
class TransformedStates : public StateGenerator {
public:
  /**
   * @brief Construct a generator yielding `transform * |r>` for each state from `inner`.
   *
   * @param transform Unitary transform operator `U` (must match `inner.Dimension()`).
   * @param inner     Wrapped source-state generator.
   */
  TransformedStates(const LinearOperator& transform, StateGenerator& inner);

  /**
   * @brief Return the vector dimension `N`.
   */
  Index Dimension() const override;

  /**
   * @brief Return `inner.NumberOfStates()`.
   */
  Index NumberOfStates() const override;

  /**
   * @brief Draw the next state from `inner` and transform it into `out`.
   */
  bool Next(Vector& out) override;

  /**
   * @brief Reset the inner generator.
   */
  void Reset() override;

  /**
   * @brief Return `inner.Label()`.
   */
  std::string Label() const override;

private:
  const LinearOperator* transform_;
  StateGenerator* inner_;
  Vector scratch_;
};

/**
 * @brief Generator yielding caller-supplied state vectors in order.
 *
 * Stores an explicit list of dense `Vector` states of uniform positive dimension and returns
 * them sequentially on `Next(out)`. Reports label `"user"`. Cost per state is `O(dimension)`.
 *
 * @warning Throws `std::invalid_argument` if the vector list is empty, if vectors have length
 *          `0`, or if vector lengths are inconsistent.
 * @par Pitfall
 *      Caller vectors are not automatically normalised; supply normalised states if unit norm
 *      is expected.
 */
class UserStates : public StateGenerator {
public:
  /**
   * @brief Construct a generator over an explicit list of state vectors.
   *
   * @param states Non-empty list of equal-length non-empty `Vector` instances.
   *
   * @warning Throws `std::invalid_argument` if `states` is empty or dimensions mismatch.
   * @par Pitfall
   *      Moves `states` into internal storage; the caller's vector is consumed when passed as
   *      an rvalue.
   */
  explicit UserStates(std::vector<Vector> states);

  /**
   * @brief Return the common length of the stored user vectors.
   *
   * @return Dimension `N > 0`.
   */
  Index Dimension() const override;

  /**
   * @brief Return the number of stored user vectors.
   *
   * @return State count `R > 0`.
   */
  Index NumberOfStates() const override;

  /**
   * @brief Copy the next user vector into `out`.
   *
   * @param out Output vector receiving the next user state.
   * @return    `true` if a state was copied, `false` when all states have been yielded.
   */
  bool Next(Vector& out) override;

  /**
   * @brief Rewind the cursor to the first user vector.
   */
  void Reset() override;

  /**
   * @brief Return `"user"`.
   *
   * @return `"user"`.
   */
  std::string Label() const override;

private:
  Index dimension_ = 0;
  std::vector<Vector> states_;
  Index current_ = 0;
};

}  // namespace lsquant

#endif  // LSQUANT_CORE_ALGEBRA_STATE_GENERATOR_HPP_
