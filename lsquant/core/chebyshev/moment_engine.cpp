/**
 * @file moment_engine.cpp
 * @brief Implementation of MomentEngine::Run.
 */
#include "core/chebyshev/moment_engine.hpp"

#include <algorithm>
#include <stdexcept>
#include <variant>
#include <vector>

#include "core/algebra/linalg.hpp"
#include "core/algebra/memory_placement.hpp"
#include "core/chebyshev/chebyshev_recursion.hpp"
#include "core/chebyshev/chebyshev_vectors.hpp"

namespace lsquant {

namespace {

using SeedOp = std::variant<ApplyOperator, ApplySeries, ApplyDerivedSeries>;

struct SideProgram {
  std::vector<SeedOp> preOps;
  bool hasExpand = false;
  int numMoments = 0;
  std::vector<int> postOperators;
};

void ApplySingleVectorOperators(const std::vector<const LinearOperator*>& operators,
                                const std::vector<int>& opIndices, const Vector& init,
                                Vector& out, Vector& scratch) {
  linalg::copy(init, out);
  for (int idx : opIndices) {
    const LinearOperator& op = *operators[static_cast<std::size_t>(idx)];
    if (op.IsIdentity()) {
      continue;
    }
    linalg::copy(out, scratch);
    op.Multiply(Scalar(1.0, 0.0), scratch, Scalar(0.0, 0.0), out);
  }
}

void RunApplySeries(const LinearOperator& rescaledH, const std::vector<Scalar>& coeffs,
                    Scalar scale, Vector& phi, Vector& cheb0, Vector& cheb1) {
  if (coeffs.empty()) {
    linalg::scal(Scalar(0.0, 0.0), phi);
    return;
  }
  linalg::copy(phi, cheb0);
  rescaledH.Multiply(Scalar(1.0, 0.0), cheb0, Scalar(0.0, 0.0), cheb1);
  linalg::scal(Scalar(0.0, 0.0), phi);
  for (std::size_t n = 0; n < coeffs.size(); ++n) {
    linalg::axpy(coeffs[n], cheb0, phi);
    rescaledH.Multiply(Scalar(2.0, 0.0), cheb1, Scalar(-1.0, 0.0), cheb0);
    cheb0.swap(cheb1);
  }
  linalg::scal(scale, phi);
}

void RunApplyDerivedSeries(const LinearOperator& rescaledH, const LinearOperator& commOp,
                           const std::vector<Scalar>& coeffs, Scalar scale, Vector& phi,
                           Vector& cheb0, Vector& cheb1, Vector& chebx0, Vector& chebx1) {
  if (coeffs.empty()) {
    linalg::scal(Scalar(0.0, 0.0), phi);
    return;
  }
  linalg::copy(phi, cheb0);
  linalg::copy(phi, chebx0);
  linalg::scal(Scalar(0.0, 0.0), chebx0);
  rescaledH.Multiply(Scalar(1.0, 0.0), cheb0, Scalar(0.0, 0.0), cheb1);
  commOp.Multiply(Scalar(1.0, 0.0), cheb0, Scalar(0.0, 0.0), chebx1);
  linalg::scal(Scalar(0.0, 0.0), phi);
  for (std::size_t n = 0; n < coeffs.size(); ++n) {
    const Scalar cn = (n == 0) ? Scalar(0.0, 0.0) : coeffs[n];
    linalg::axpy(cn, chebx0, phi);
    rescaledH.Multiply(Scalar(2.0, 0.0), cheb1, Scalar(-1.0, 0.0), cheb0);
    rescaledH.Multiply(Scalar(2.0, 0.0), chebx1, Scalar(-1.0, 0.0), chebx0);
    commOp.Multiply(Scalar(2.0, 0.0), cheb1, Scalar(1.0, 0.0), chebx0);
    cheb0.swap(cheb1);
    chebx0.swap(chebx1);
  }
  linalg::scal(scale, phi);
}

void ApplySeedOps(const LinearOperator& rescaledH,
                  const std::vector<const LinearOperator*>& operators,
                  const std::vector<SeedOp>& ops, const Vector& init, Vector& out,
                  Vector& scratch, Vector& cheb0, Vector& cheb1, Vector& chebx0,
                  Vector& chebx1) {
  linalg::copy(init, out);
  for (const SeedOp& op : ops) {
    if (const auto* app = std::get_if<ApplyOperator>(&op)) {
      const LinearOperator& linOp = *operators[static_cast<std::size_t>(app->operatorIndex)];
      if (linOp.IsIdentity()) {
        continue;
      }
      linalg::copy(out, scratch);
      linOp.Multiply(Scalar(1.0, 0.0), scratch, Scalar(0.0, 0.0), out);
    } else if (const auto* ser = std::get_if<ApplySeries>(&op)) {
      RunApplySeries(rescaledH, ser->coeffs, ser->scale, out, cheb0, cheb1);
    } else if (const auto* dser = std::get_if<ApplyDerivedSeries>(&op)) {
      const LinearOperator& commOp = *operators[static_cast<std::size_t>(dser->operatorIndex)];
      RunApplyDerivedSeries(rescaledH, commOp, dser->coeffs, dser->scale, out, cheb0, cheb1,
                            chebx0, chebx1);
    }
  }
}

void AccumulateWeightedDot(Scalar& target, const Scalar& weight, const Vector& bra,
                           const Vector& ket) {
  const Scalar d = linalg::dot(bra, ket);
  if (weight.imag() == 0.0) {
    target += weight.real() * d;
  } else {
    target += weight * d;
  }
}

int PaddedMoments(int numMoments, int batchSize) {
  return ((numMoments + batchSize - 1) / batchSize) * batchSize;
}

Vector& SelectRegister(Register reg, Vector& ket, Vector& bra, Vector& aux) {
  switch (reg) {
    case Register::Ket:
      return ket;
    case Register::Bra:
      return bra;
    case Register::Aux:
      return aux;
  }
  return ket;
}

}  // namespace

MomentEngine::MomentEngine(const LinearOperator& rescaledH,
                           std::vector<const LinearOperator*> operators,
                           StateGenerator& states, int batchSize)
    : rescaledH_(&rescaledH),
      operators_(std::move(operators)),
      states_(&states),
      batchSize_(batchSize) {
  if (rescaledH.Dimension() <= 0) {
    throw std::invalid_argument("MomentEngine: rescaledH dimension must be positive");
  }
  if (states.Dimension() != rescaledH.Dimension()) {
    throw std::invalid_argument("MomentEngine: states dimension mismatch");
  }
  if (batchSize <= 0) {
    throw std::invalid_argument("MomentEngine: batchSize must be positive");
  }
  for (const LinearOperator* op : operators_) {
    if (op == nullptr) {
      throw std::invalid_argument("MomentEngine: null operator pointer");
    }
    if (op->Dimension() != rescaledH.Dimension()) {
      throw std::invalid_argument("MomentEngine: operator dimension mismatch");
    }
  }
}

ChebyshevMoments MomentEngine::Run(const Plan& plan) {
  if (plan.code.empty()) {
    throw std::invalid_argument("MomentEngine::Run: plan code must not be empty");
  }

  std::size_t timeAxisPos = plan.code.size();
  for (std::size_t pc = 0; pc < plan.code.size(); ++pc) {
    if (std::holds_alternative<TimeAxis>(plan.code[pc])) {
      if (timeAxisPos != plan.code.size()) {
        throw std::invalid_argument("MomentEngine::Run: TimeAxis may appear at most once");
      }
      timeAxisPos = pc;
    }
  }

  const bool hasTimeAxis = (timeAxisPos != plan.code.size());

  std::vector<SeedOp> seedKetOps;
  std::vector<SeedOp> seedBraOps;
  bool seedMirror = false;
  TimeAxis timeAxisSpec;

  if (hasTimeAxis) {
    const auto* ta = std::get_if<TimeAxis>(&plan.code[timeAxisPos]);
    if (ta->numSteps <= 0) {
      throw std::invalid_argument("MomentEngine::Run: TimeAxis numSteps must be positive");
    }
    for (const AxisStep& ev : ta->evolve) {
      if (const auto* ed = std::get_if<EvolveDerived>(&ev)) {
        if (ed->operatorIndex < 0 ||
            ed->operatorIndex >= static_cast<int>(operators_.size())) {
          throw std::invalid_argument(
              "MomentEngine::Run: EvolveDerived operatorIndex out of range");
        }
      }
    }
    timeAxisSpec = *ta;

    std::vector<SeedOp>* activeSeedOps = &seedKetOps;
    for (std::size_t pc = 0; pc < timeAxisPos; ++pc) {
      const Instruction& inst = plan.code[pc];
      if (const auto* app = std::get_if<ApplyOperator>(&inst)) {
        if (app->operatorIndex < 0 ||
            app->operatorIndex >= static_cast<int>(operators_.size())) {
          throw std::invalid_argument("MomentEngine::Run: operatorIndex out of range");
        }
        activeSeedOps->push_back(*app);
      } else if (const auto* ser = std::get_if<ApplySeries>(&inst)) {
        activeSeedOps->push_back(*ser);
      } else if (const auto* dser = std::get_if<ApplyDerivedSeries>(&inst)) {
        if (dser->operatorIndex < 0 ||
            dser->operatorIndex >= static_cast<int>(operators_.size())) {
          throw std::invalid_argument("MomentEngine::Run: operatorIndex out of range");
        }
        activeSeedOps->push_back(*dser);
      } else if (std::holds_alternative<MirrorToBra>(inst)) {
        if (seedMirror) {
          throw std::invalid_argument(
              "MomentEngine::Run: MirrorToBra may appear at most once before TimeAxis");
        }
        seedMirror = true;
        activeSeedOps = &seedBraOps;
      } else {
        throw std::invalid_argument(
            "MomentEngine::Run: Expand and Trace cannot precede TimeAxis");
      }
    }
  }

  SideProgram ketProg;
  SideProgram braProg;
  SideProgram* activeSide = &ketProg;
  bool seenMirror = false;
  bool seenTrace = false;
  std::vector<Scalar> traceWeights;

  const std::size_t bodyStart = hasTimeAxis ? (timeAxisPos + 1) : 0;
  for (std::size_t pc = bodyStart; pc < plan.code.size(); ++pc) {
    if (seenTrace) {
      throw std::invalid_argument("MomentEngine::Run: Trace must be the final instruction");
    }
    const Instruction& inst = plan.code[pc];
    if (const auto* app = std::get_if<ApplyOperator>(&inst)) {
      if (app->operatorIndex < 0 ||
          app->operatorIndex >= static_cast<int>(operators_.size())) {
        throw std::invalid_argument("MomentEngine::Run: operatorIndex out of range");
      }
      if (!activeSide->hasExpand) {
        activeSide->preOps.push_back(*app);
      } else {
        activeSide->postOperators.push_back(app->operatorIndex);
      }
    } else if (const auto* ser = std::get_if<ApplySeries>(&inst)) {
      if (activeSide->hasExpand) {
        throw std::invalid_argument("MomentEngine::Run: ApplySeries must precede Expand");
      }
      activeSide->preOps.push_back(*ser);
    } else if (const auto* dser = std::get_if<ApplyDerivedSeries>(&inst)) {
      if (dser->operatorIndex < 0 ||
          dser->operatorIndex >= static_cast<int>(operators_.size())) {
        throw std::invalid_argument("MomentEngine::Run: operatorIndex out of range");
      }
      if (activeSide->hasExpand) {
        throw std::invalid_argument("MomentEngine::Run: ApplyDerivedSeries must precede Expand");
      }
      activeSide->preOps.push_back(*dser);
    } else if (const auto* exp = std::get_if<Expand>(&inst)) {
      if (exp->numMoments <= 0) {
        throw std::invalid_argument("MomentEngine::Run: Expand numMoments must be positive");
      }
      if (activeSide->hasExpand) {
        throw std::runtime_error("Expand: rank > 2 not implemented");
      }
      activeSide->hasExpand = true;
      activeSide->numMoments = exp->numMoments;
    } else if (std::holds_alternative<MirrorToBra>(inst)) {
      if (seenMirror) {
        throw std::invalid_argument("MomentEngine::Run: MirrorToBra may appear at most once");
      }
      seenMirror = true;
      activeSide = &braProg;
    } else if (const auto* tr = std::get_if<Trace>(&inst)) {
      seenTrace = true;
      traceWeights = tr->weights;
    }
  }

  if (!seenTrace) {
    throw std::invalid_argument("MomentEngine::Run: plan must end with a Trace instruction");
  }

  std::vector<int> dims;
  if (braProg.hasExpand) {
    dims.push_back(braProg.numMoments);
  }
  if (ketProg.hasExpand) {
    dims.push_back(ketProg.numMoments);
  }
  if (hasTimeAxis) {
    if (braProg.hasExpand) {
      throw std::invalid_argument(
          "MomentEngine::Run: TimeAxis body does not support bra Expand");
    }
    dims.push_back(timeAxisSpec.numSteps);
  }

  ChebyshevMoments mu(dims);
  mu.systemSize = rescaledH_->Dimension();
  mu.numberOfStates = static_cast<int>(states_->NumberOfStates());
  mu.stateLabel = states_->Label();
  mu.word = plan.word;
  if (hasTimeAxis) {
    mu.numberOfTimeSteps = timeAxisSpec.numSteps;
  }

  const std::size_t tensorSize = mu.Data().size();
  const bool validMomentWeights =
      hasTimeAxis && ketProg.hasExpand &&
      (traceWeights.size() == static_cast<std::size_t>(ketProg.numMoments));
  if (traceWeights.size() != 1 && traceWeights.size() != tensorSize && !validMomentWeights) {
    throw std::invalid_argument(
        "MomentEngine::Run: Trace weights size must be 1 or match tensor size");
  }
  const bool scalarWeight = (traceWeights.size() == 1);

  const Index dim = rescaledH_->Dimension();
  states_->Reset();

#ifdef LSQUANT_NATIVE_BUILD
  Vector r = MakeZeroVector(dim);
  Vector ketSeed = MakeZeroVector(dim);
  Vector braSeed = MakeZeroVector(dim);
  Vector auxSeed = MakeZeroVector(dim);
  Vector ketInit = MakeZeroVector(dim);
  Vector braInit = MakeZeroVector(dim);
  Vector scratch = MakeZeroVector(dim);
  Vector cheb0 = MakeZeroVector(dim);
  Vector cheb1 = MakeZeroVector(dim);
  Vector chebx0 = MakeZeroVector(dim);
  Vector chebx1 = MakeZeroVector(dim);
  Vector phiT = MakeZeroVector(dim);
  Vector postVec = MakeZeroVector(dim);
#else
  Vector r(static_cast<std::size_t>(dim), Scalar(0.0, 0.0));
  Vector ketSeed(static_cast<std::size_t>(dim), Scalar(0.0, 0.0));
  Vector braSeed(static_cast<std::size_t>(dim), Scalar(0.0, 0.0));
  Vector auxSeed(static_cast<std::size_t>(dim), Scalar(0.0, 0.0));
  Vector ketInit(static_cast<std::size_t>(dim), Scalar(0.0, 0.0));
  Vector braInit(static_cast<std::size_t>(dim), Scalar(0.0, 0.0));
  Vector scratch(static_cast<std::size_t>(dim), Scalar(0.0, 0.0));
  Vector cheb0(static_cast<std::size_t>(dim), Scalar(0.0, 0.0));
  Vector cheb1(static_cast<std::size_t>(dim), Scalar(0.0, 0.0));
  Vector chebx0(static_cast<std::size_t>(dim), Scalar(0.0, 0.0));
  Vector chebx1(static_cast<std::size_t>(dim), Scalar(0.0, 0.0));
  Vector phiT(static_cast<std::size_t>(dim), Scalar(0.0, 0.0));
  Vector postVec(static_cast<std::size_t>(dim), Scalar(0.0, 0.0));
#endif

  if (hasTimeAxis) {
    ChebyshevRecursion rec(*rescaledH_);
    while (states_->Next(r)) {
      linalg::copy(r, auxSeed);
      ApplySeedOps(*rescaledH_, operators_, seedKetOps, r, ketSeed, scratch, cheb0, cheb1,
                   chebx0, chebx1);
      const Vector& braSeedStart = seedMirror ? ketSeed : r;
      ApplySeedOps(*rescaledH_, operators_, seedBraOps, braSeedStart, braSeed, scratch, cheb0,
                   cheb1, chebx0, chebx1);

      for (int stepIdx = 0; stepIdx < timeAxisSpec.numSteps; ++stepIdx) {
        ApplySeedOps(*rescaledH_, operators_, ketProg.preOps, ketSeed, ketInit, scratch, cheb0,
                     cheb1, chebx0, chebx1);
        const Vector& braBodyStart = seenMirror ? ketInit : braSeed;
        ApplySeedOps(*rescaledH_, operators_, braProg.preOps, braBodyStart, braInit, scratch,
                     cheb0, cheb1, chebx0, chebx1);

        if (ketProg.hasExpand) {
          const int numMoms = ketProg.numMoments;
          rec.Start(ketInit);
          for (int m = 0; m < numMoms; ++m) {
            const Vector* ketVec = &rec.Current();
            if (!ketProg.postOperators.empty()) {
              ApplySingleVectorOperators(operators_, ketProg.postOperators, rec.Current(),
                                         postVec, scratch);
              ketVec = &postVec;
            }
            std::size_t wIdx = 0;
            if (validMomentWeights) {
              wIdx = static_cast<std::size_t>(m);
            } else if (!scalarWeight) {
              wIdx = static_cast<std::size_t>(m) *
                         static_cast<std::size_t>(timeAxisSpec.numSteps) +
                     static_cast<std::size_t>(stepIdx);
            }
            AccumulateWeightedDot(mu.At(m, stepIdx), traceWeights[wIdx], braInit, *ketVec);
            rec.Step();
          }
        } else {
          const std::size_t wIdx = scalarWeight ? 0 : static_cast<std::size_t>(stepIdx);
          AccumulateWeightedDot(mu.At(stepIdx), traceWeights[wIdx], braInit, ketInit);
        }

        for (const AxisStep& ev : timeAxisSpec.evolve) {
          if (const auto* er = std::get_if<EvolveRegister>(&ev)) {
            Vector& targetVec = SelectRegister(er->target, ketSeed, braSeed, auxSeed);
            RunApplySeries(*rescaledH_, er->coeffs, er->scale, targetVec, cheb0, cheb1);
          } else if (const auto* ed = std::get_if<EvolveDerived>(&ev)) {
            const Vector& srcVec = SelectRegister(ed->source, ketSeed, braSeed, auxSeed);
            Vector& targetVec = SelectRegister(ed->target, ketSeed, braSeed, auxSeed);
            linalg::copy(srcVec, phiT);
            const LinearOperator& commOp =
                *operators_[static_cast<std::size_t>(ed->operatorIndex)];
            RunApplyDerivedSeries(*rescaledH_, commOp, ed->coeffs, ed->scale, phiT, cheb0,
                                  cheb1, chebx0, chebx1);
            linalg::axpy(Scalar(1.0, 0.0), phiT, targetVec);
          } else if (const auto* cp = std::get_if<CopyRegister>(&ev)) {
            const Vector& fromVec = SelectRegister(cp->from, ketSeed, braSeed, auxSeed);
            Vector& toVec = SelectRegister(cp->to, ketSeed, braSeed, auxSeed);
            linalg::copy(fromVec, toVec);
          }
        }
      }
    }
    return mu;
  }

  if (!braProg.hasExpand && !ketProg.hasExpand) {
    // Rank 0: neither side opens a Chebyshev index.
    while (states_->Next(r)) {
      ApplySeedOps(*rescaledH_, operators_, ketProg.preOps, r, ketInit, scratch, cheb0, cheb1,
                   chebx0, chebx1);
      const Vector& braStart = r;
      ApplySeedOps(*rescaledH_, operators_, braProg.preOps, braStart, braInit, scratch, cheb0,
                   cheb1, chebx0, chebx1);
      AccumulateWeightedDot(mu.Data()[0], traceWeights[0], braInit, ketInit);
    }
    return mu;
  }

  if (!braProg.hasExpand && ketProg.hasExpand) {
    // Rank 1 (ket expanded): matches legacy SpectralMoments when postOperators is empty.
    const int numMoms = ketProg.numMoments;
    if (ketProg.postOperators.empty()) {
      ChebyshevRecursion rec(*rescaledH_);
      while (states_->Next(r)) {
        ApplySeedOps(*rescaledH_, operators_, ketProg.preOps, r, ketInit, scratch, cheb0, cheb1,
                     chebx0, chebx1);
        const Vector& braStart = r;
        ApplySeedOps(*rescaledH_, operators_, braProg.preOps, braStart, braInit, scratch, cheb0,
                     cheb1, chebx0, chebx1);
        rec.Start(ketInit);
        for (int m = 0; m < numMoms; ++m) {
          const Scalar& w = scalarWeight ? traceWeights[0]
                                         : traceWeights[static_cast<std::size_t>(m)];
          AccumulateWeightedDot(mu.At(m), w, braInit, rec.Current());
          rec.Step();
        }
      }
    } else {
      const int effBatch = std::min(batchSize_, numMoms);
      const int paddedM = PaddedMoments(numMoms, effBatch);
      ChebyshevRecursion rec(*rescaledH_);
      ChebyshevVectors vecR(dim, effBatch);
      while (states_->Next(r)) {
        ApplySeedOps(*rescaledH_, operators_, ketProg.preOps, r, ketInit, scratch, cheb0, cheb1,
                     chebx0, chebx1);
        const Vector& braStart = r;
        ApplySeedOps(*rescaledH_, operators_, braProg.preOps, braStart, braInit, scratch, cheb0,
                     cheb1, chebx0, chebx1);
        rec.Start(ketInit);
        for (int mR = 0; mR < paddedM; mR += effBatch) {
          vecR.Fill(rec);
          for (int opIdx : ketProg.postOperators) {
            vecR.ApplyToAll(*operators_[static_cast<std::size_t>(opIdx)]);
          }
          for (int m0 = 0; m0 < effBatch; ++m0) {
            const int m = mR + m0;
            if (m < numMoms) {
              const Scalar& w = scalarWeight ? traceWeights[0]
                                             : traceWeights[static_cast<std::size_t>(m)];
              AccumulateWeightedDot(mu.At(m), w, braInit, vecR.At(m0));
            }
          }
        }
      }
    }
    return mu;
  }

  if (braProg.hasExpand && !ketProg.hasExpand) {
    // Rank 1 (bra expanded).
    const int numMoms = braProg.numMoments;
    const int effBatch = std::min(batchSize_, numMoms);
    const int paddedM = PaddedMoments(numMoms, effBatch);
    ChebyshevRecursion recL(*rescaledH_);
    ChebyshevVectors vecL(dim, effBatch);
    while (states_->Next(r)) {
      ApplySeedOps(*rescaledH_, operators_, ketProg.preOps, r, ketInit, scratch, cheb0, cheb1,
                   chebx0, chebx1);
      const Vector& braStart = r;
      ApplySeedOps(*rescaledH_, operators_, braProg.preOps, braStart, braInit, scratch, cheb0,
                   cheb1, chebx0, chebx1);
      recL.Start(braInit);
      for (int mL = 0; mL < paddedM; mL += effBatch) {
        vecL.Fill(recL);
        for (int opIdx : braProg.postOperators) {
          vecL.ApplyToAll(*operators_[static_cast<std::size_t>(opIdx)]);
        }
        for (int m0 = 0; m0 < effBatch; ++m0) {
          const int m = mL + m0;
          if (m < numMoms) {
            const Scalar& w = scalarWeight ? traceWeights[0]
                                           : traceWeights[static_cast<std::size_t>(m)];
            AccumulateWeightedDot(mu.At(m), w, vecL.At(m0), ketInit);
          }
        }
      }
    }
    return mu;
  }

  // Rank 2 (both bra and ket expanded): matches legacy CorrelationExpansionMoments block order.
  const int numMomsL = braProg.numMoments;
  const int numMomsR = ketProg.numMoments;
  const int effBatchL = std::min(batchSize_, numMomsL);
  const int effBatchR = std::min(batchSize_, numMomsR);
  const int paddedML = PaddedMoments(numMomsL, effBatchL);
  const int paddedMR = PaddedMoments(numMomsR, effBatchR);

  ChebyshevRecursion recL(*rescaledH_);
  ChebyshevRecursion recR(*rescaledH_);
  ChebyshevVectors vecL(dim, effBatchL);
  ChebyshevVectors vecR(dim, effBatchR);

  while (states_->Next(r)) {
    ApplySeedOps(*rescaledH_, operators_, ketProg.preOps, r, ketSeed, scratch, cheb0, cheb1,
                 chebx0, chebx1);
    // Static words build both sides independently from the SAME trace state.
    // Any shared preparation must be explicit on each side of MirrorToBra.
    const Vector& braStart = r;
    ApplySeedOps(*rescaledH_, operators_, braProg.preOps, braStart, braInit, scratch, cheb0,
                 cheb1, chebx0, chebx1);
    recL.Start(braInit);
    for (int mL = 0; mL < paddedML; mL += effBatchL) {
      vecL.Fill(recL);
      for (int opIdx : braProg.postOperators) {
        vecL.ApplyToAll(*operators_[static_cast<std::size_t>(opIdx)]);
      }

      linalg::copy(ketSeed, ketInit);
      recR.Start(ketInit);
      for (int mR = 0; mR < paddedMR; mR += effBatchR) {
        vecR.Fill(recR);
        for (int opIdx : ketProg.postOperators) {
          vecR.ApplyToAll(*operators_[static_cast<std::size_t>(opIdx)]);
        }

#ifdef LSQUANT_NATIVE_BUILD
        std::vector<const Vector*> left, right;
        const int liveL = std::min(effBatchL, numMomsL - mL);
        const int liveR = std::min(effBatchR, numMomsR - mR);
        for (int m = 0; m < liveL; ++m) left.push_back(&vecL.At(m));
        for (int n = 0; n < liveR; ++n) right.push_back(&vecR.At(n));
        const auto table = linalg::MomentTable(left, right);
        for (int m = 0; m < liveL; ++m) {
          for (int n = 0; n < liveR; ++n) {
            const std::size_t index = static_cast<std::size_t>(mL + m) * numMomsR + mR + n;
            const Scalar& weight = scalarWeight ? traceWeights[0] : traceWeights[index];
            const Scalar& value = table[static_cast<std::size_t>(m) * liveR + n];
            if (weight.imag() == 0.0) mu.At(mL + m, mR + n) += weight.real() * value;
            else mu.At(mL + m, mR + n) += weight * value;
          }
        }
#else
        // Each row has one owner; dot products and trace-state additions keep
        // their serial order, so thread count cannot change the moment values.
#pragma omp parallel for schedule(static)
        for (int m0 = 0; m0 < effBatchL; ++m0) {
          const int row = mL + m0;
          if (row >= numMomsL) {
            continue;
          }
          for (int m1 = 0; m1 < effBatchR; ++m1) {
            const int col = mR + m1;
            if (col >= numMomsR) {
              continue;
            }
            const std::size_t wIdx =
                static_cast<std::size_t>(row) * static_cast<std::size_t>(numMomsR) +
                static_cast<std::size_t>(col);
            const Scalar& w = scalarWeight ? traceWeights[0] : traceWeights[wIdx];
            AccumulateWeightedDot(mu.At(row, col), w, vecL.At(m0), vecR.At(m1));
          }
        }
#endif
      }
    }
  }

  return mu;
}

}  // namespace lsquant
