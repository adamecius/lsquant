/**
 * @file plan.hpp
 * @brief Bytecode representation of a KPM operator word.
 *
 * Defines the physical-name-free instruction set (`ApplyOperator`, `Expand`,
 * `ApplySeries`, `ApplyDerivedSeries`, `MirrorToBra`, `TimeAxis`, `Trace`),
 * register identifiers, register-update step descriptors, and `Plan` container
 * executed by `MomentEngine::Run`.
 */
#ifndef LSQUANT_CORE_CHEBYSHEV_PLAN_HPP_
#define LSQUANT_CORE_CHEBYSHEV_PLAN_HPP_

#include <string>
#include <variant>
#include <vector>

#include "core/algebra/types.hpp"

namespace lsquant {

/**
 * @brief Seed-vector register identifier for `MomentEngine` state machines.
 *
 * `MomentEngine` maintains three vector registers of dimension `N` per trace state:
 * `Ket` (the active ket seed), `Bra` (the active bra seed), and `Aux` (an auxiliary
 * seed register used by coupled register updates between slices of a trailing
 * tensor axis).
 */
enum class Register {
  Ket,
  Bra,
  Aux,
};

/**
 * @brief Instruction applying registered static operator `O_{operatorIndex}` to every live vector.
 *
 * Replaces each live vector `v` on the active side (ket or bra) with `O_k v`.
 *
 * @warning `operatorIndex` must be in `[0, operators.size())` of the executing `MomentEngine`.
 * @par Pitfall
 *      When applied after `Expand`, transforms every vector in the active Chebyshev batch.
 */
struct ApplyOperator {
  int operatorIndex;
};

/**
 * @brief Instruction opening one tensor index of length `numMoments` via Chebyshev recursion.
 *
 * Replaces a live vector `v` on the active side with the sequence `{T_0(H~) v, ..., T_{M-1}(H~) v}`.
 *
 * @warning `numMoments` must be strictly positive; at most one `Expand` per side is supported.
 * @par Pitfall
 *      Each `Expand` increases the rank of the output `ChebyshevMoments` tensor by 1.
 */
struct Expand {
  int numMoments;
};

/**
 * @brief Instruction contracting a Chebyshev series `scale * sum_n coeffs[n] T_n(H~) v` on the fly.
 *
 * Replaces the active seed vector `v` using the forward three-term time-evolution
 * recursion: `cheb0 = v`, `cheb1 = H~ cheb0`, `v = 0`, and for
 * `n = 0 .. coeffs.size() - 1`: `axpy(coeffs[n], cheb0, v)`,
 * `cheb0 <- 2 H~ cheb1 - cheb0`, `swap(cheb0, cheb1)`, followed by a final
 * `scal(scale, v)`. An empty `coeffs` is the zero map (`v <- 0`), which is how a seed
 * register is cleared, and `msd_plan` relies on it.
 *
 * @warning Must appear before `Expand` on the active side.
 * @par Pitfall
 *      Keep the global phase factor in `scale` rather than folding it into `coeffs` when
 *      bit-identical reproduction of the legacy time-evolution recursion is required.
 */
struct ApplySeries {
  std::vector<Scalar> coeffs;
  Scalar scale = Scalar(1.0, 0.0);
};

/**
 * @brief Instruction contracting a derived Chebyshev series `scale * sum_n coeffs[n] D_n v` on the fly.
 *
 * Replaces the active seed vector `v` using the coupled `(T_n, D_n)` commutator-evolution
 * recursion: `D_0 v = 0`, `D_1 v = O_k T_0 v`, and
 * `D_{n+1} v = 2 H~ D_n v - D_{n-1} v + 2 O_k T_n v`, where `k = operatorIndex`,
 * followed by a final `scal(scale, v)`. The `n = 0` term contributes zero regardless of
 * `coeffs[0]`. An empty `coeffs` is the zero map (`v <- 0`), which is how a seed
 * register is cleared, and `msd_plan` relies on it.
 *
 * @warning `operatorIndex` must be in `[0, operators.size())` of the executing `MomentEngine`.
 * @par Pitfall
 *      The recursion order per step (`T` update, `D` update, `+2 O_k T_n` accumulation, then
 *      swaps) must match the legacy commutator-evolution recursion verbatim.
 */
struct ApplyDerivedSeries {
  std::vector<Scalar> coeffs;
  int operatorIndex = 0;
  Scalar scale = Scalar(1.0, 0.0);
};

/**
 * @brief Switch a static word to the bra side on a fresh copy of the same trace state.
 *
 * In a static word, subsequent instructions start from the original `|r>`,
 * independently of ket-side preparation. No new trace state is drawn. Shared
 * transforms must be written explicitly on both sides.
 * The legacy `TimeAxis` register program uses its established copy semantics:
 * seed preparation copies ketSeed; its body copies the current step's ketInit.
 *
 * @warning May appear at most once in seed preparation and at most once inside a `TimeAxis` body.
 * @par Pitfall
 *      If `MirrorToBra` is omitted, the bra seed remains the initial trace state `|r>`.
 */
struct MirrorToBra {};

/**
 * @brief Register-update step between trailing-axis slices applying `ApplySeries` to `target`.
 *
 * Replaces `seed[target]` with `scale * sum_n coeffs[n] T_n(H~) seed[target]`.
 */
struct EvolveRegister {
  Register target = Register::Ket;
  std::vector<Scalar> coeffs;
  Scalar scale = Scalar(1.0, 0.0);
};

/**
 * @brief Register-update step between trailing-axis slices accumulating a derived Chebyshev series from `source` into `target`.
 *
 * Computes `phiT = scale * sum_n coeffs[n] D_n seed[source]` with operator
 * `O_{operatorIndex}` (using the coupled commutator recursion) and accumulates
 * `seed[target] += phiT` via `axpy(1.0, phiT, seed[target])`.
 */
struct EvolveDerived {
  Register target = Register::Ket;
  Register source = Register::Aux;
  std::vector<Scalar> coeffs;
  int operatorIndex = 0;
  Scalar scale = Scalar(1.0, 0.0);
};

/**
 * @brief Register-update step between trailing-axis slices copying one seed register into another (`seed[to] = seed[from]`).
 */
struct CopyRegister {
  Register from = Register::Ket;
  Register to = Register::Bra;
};

/**
 * @brief Variant holding one register-update action executed between slices of a `TimeAxis`.
 */
using AxisStep = std::variant<EvolveRegister, EvolveDerived, CopyRegister>;

/**
 * @brief Instruction opening a trailing tensor axis of length `numSteps` with register updates between slices.
 *
 * May appear at most once in a `Plan`, after all initial seed preparation. For each step
 * `n = 0 .. numSteps - 1`, executes the instructions following `TimeAxis` (the body) on
 * working copies of `seed[Ket]` and `seed[Bra]`, accumulating into slice `mu[..., n]`, and
 * then applies `evolve` in order to update the persistent `seed` registers.
 *
 * @warning `numSteps` must be strictly positive.
 * @par Pitfall
 *      Body instructions act on per-step working copies of `seed[Ket]` and `seed[Bra]` so
 *      that pre-expansion operators inside the body do not corrupt the evolving seeds.
 */
struct TimeAxis {
  int numSteps = 0;
  std::vector<AxisStep> evolve;
};

/**
 * @brief Terminal instruction contracting live bra and ket vectors into the moment tensor.
 *
 * Accumulates `mu[i_bra..., i_ket...] += weights[...] * dot(bra_{i_bra}, ket_{i_ket})` for each
 * trace state. `weights` must have size `1` (uniform scalar weight), match the moment-axis
 * length `M` (for `TimeAxis` plans), or match the total element count of `mu`.
 *
 * @warning Must appear as the final instruction of every `Plan`.
 * @par Pitfall
 *      Real weights are multiplied as `double * Scalar` to preserve IEEE-754 signed zeros in
 *      imaginary components.
 */
struct Trace {
  std::vector<Scalar> weights;
};

/**
 * @brief Variant holding any one of the engine instructions.
 */
using Instruction = std::variant<ApplyOperator, Expand, ApplySeries, ApplyDerivedSeries,
                                 MirrorToBra, TimeAxis, Trace>;

/**
 * @brief Compiled instruction sequence and descriptive word string for `MomentEngine::Run`.
 *
 * Pairs the executable `code` vector with the human-readable `word` label stored in the
 * resulting `ChebyshevMoments` metadata.
 *
 * @warning `code` must end with a single `Trace` instruction.
 * @par Pitfall
 *      Operator indices in `code` refer positionally to the `operators` vector passed to
 *      `MomentEngine`.
 */
struct Plan {
  std::vector<Instruction> code;
  std::string word;
};

}  // namespace lsquant

#endif  // LSQUANT_CORE_CHEBYSHEV_PLAN_HPP_
