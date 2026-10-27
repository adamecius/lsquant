/**
 * @file moment_engine.hpp
 * @brief Bytecode execution engine evaluating Chebyshev moment tensors from a Plan.
 *
 * Executes the five-instruction `Plan` state machine over trace vectors yielded by a
 * `StateGenerator`, reproducing the exact recursion, batching, and dot-product accumulation
 * order of legacy `SpectralMoments` and `CorrelationExpansionMoments`.
 */
#ifndef LSQUANT_CORE_CHEBYSHEV_MOMENT_ENGINE_HPP_
#define LSQUANT_CORE_CHEBYSHEV_MOMENT_ENGINE_HPP_

#include <vector>

#include "core/algebra/linear_operator.hpp"
#include "core/algebra/state_generator.hpp"
#include "core/chebyshev/chebyshev_moments.hpp"
#include "core/chebyshev/plan.hpp"

namespace lsquant {

/**
 * @brief Execution engine evaluating `Plan` programs into `ChebyshevMoments` tensors.
 *
 * Solves the problem of computing arbitrary rank-0, rank-1, and rank-2 Chebyshev moment
 * traces with a single representation-independent state machine. For each trace vector `|r>`
 * produced sequentially by `states`, builds the ket and bra sides according to `plan.code`
 * and accumulates `<bra | ket>` inner products in the exact block order of the legacy
 * solvers. Complexity per state is `O(M * cost(H~))` for rank 1 and
 * `O((M / batchSize) * M * cost(H~) + M^2 * N)` for rank 2.
 *
 * @warning All operators and the state generator must share `rescaledH.Dimension()`, and
 *          `batchSize` must be positive.
 * @par Pitfall
 *      `Run` resets `states` at the start of execution so repeated calls on the same engine
 *      traverse the full state sequence deterministically.
 */
class MomentEngine {
public:
  /**
   * @brief Construct a `MomentEngine` bound to `rescaledH`, `operators`, `states`, and `batchSize`.
   *
   * @param rescaledH Rescaled Hamiltonian `H~` with spectrum in `[-1, 1]`.
   * @param operators Positional table of static `LinearOperator` pointers referenced by `ApplyOperator`.
   * @param states    Trace-vector generator yielding source states `|r>`.
   * @param batchSize Number of Chebyshev vectors buffered per batch during `Expand` (`B > 0`).
   *
   * @warning Throws `std::invalid_argument` if `batchSize <= 0`, if any operator pointer is
   *          null, or if any dimension differs from `rescaledH.Dimension()`.
   * @par Pitfall
   *      When `Expand(M)` has `M` not divisible by `batchSize`, `M` is rounded up to the next
   *      multiple of `batchSize` during recursion and truncated on accumulation, matching
   *      legacy `CorrelationExpansionMoments`.
   */
  MomentEngine(const LinearOperator& rescaledH, std::vector<const LinearOperator*> operators,
               StateGenerator& states, int batchSize);

  /**
   * @brief Execute `plan` across all trace states and return the accumulated `ChebyshevMoments`.
   *
   * Validates `plan.code`, rewinds `states` via `Reset()`, and runs the ket/bra state machine
   * sequentially over every trace vector. Populates `systemSize`, `numberOfStates`,
   * `stateLabel`, and `word` on the returned `ChebyshevMoments`; the caller sets `window`.
   *
   * @param plan Instruction sequence ending with `Trace`.
   * @return     Accumulated `ChebyshevMoments` tensor of rank equal to the number of `Expand`
   *             instructions in `plan`.
   *
   * @warning Throws `std::runtime_error` if rank > 2 (a second `Expand` on the same side),
   *          and `std::invalid_argument` if `plan` is malformed.
   * @par Pitfall
   *      Rank-2 plans with `Trace(1)` return raw state sums; call `ApplyKpmWeights` afterwards
   *      when reproducing legacy `CorrelationExpansionMoments`.
   */
  ChebyshevMoments Run(const Plan& plan);

private:
  const LinearOperator* rescaledH_;
  std::vector<const LinearOperator*> operators_;
  StateGenerator* states_;
  int batchSize_;
};

namespace chebyshev {
using lsquant::MomentEngine;
}  // namespace chebyshev

}  // namespace lsquant

#endif  // LSQUANT_CORE_CHEBYSHEV_MOMENT_ENGINE_HPP_
