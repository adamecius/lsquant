/**
 * @file stage6_nonorth.cpp
 * @brief Stage-6 gate executable evaluating `dos`, `kubo`, and `kuboxy` moments across
 *        `lowdin`, `generalised`, and `naive` non-orthogonal routes.
 *
 * Parses command-line options specifying the fixture directory, model label, spectral window,
 * non-orthogonal route (`lowdin`, `generalised`, or `naive`), overlap applier (`exact`,
 * `ischo`, or `chebyshev`), Chebyshev applier order, PCG tolerance, trace generator (`exact`
 * or `random`), moment count `M`, observable word (`dos`, `kubo`, or `kuboxy`), and output
 * directory. Note that `--route naive` implements the legacy non-orthogonal estimator
 * (`theory.md` section 4.1): recursion on `S^-1 H~` (`GeneralisedRescaledOperator`) with
 * `S^-1 V'` (`MetricInverseProduct(S^-1, V')`) on both sides and no `S` on the bra; it
 * reproduces the DOS (where no velocity enters) but is wrong for Kubo moments when `S != I`
 * because the symmetric trace estimator is not metric-covariant. Executes the corresponding
 * bytecode `Plan` with `MomentEngine` without modifying the engine or plan instruction set,
 * and writes both the legacy `.chebmom1D`/`.chebmom2D` file and the unified `.chebmom` file.
 * Also supports `--check-appliers` to benchmark the C++ overlap appliers against
 * `ExactInverseOverlap` and `ExactInverseSqrtOverlap` on deterministic random vectors.
 */
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <exception>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include "core/algebra/composite_operators.hpp"
#include "core/algebra/linalg.hpp"
#include "core/algebra/metric_policy.hpp"
#include "core/algebra/overlap_appliers.hpp"
#include "core/algebra/sparse_matrix.hpp"
#include "core/algebra/state_generator.hpp"
#include "core/chebyshev/chebyshev_moments.hpp"
#include "core/chebyshev/moment_engine.hpp"
#include "core/chebyshev/plan.hpp"
#include "core/chebyshev/rescaled_operator.hpp"
#include "core/chebyshev/spectral_window.hpp"
#include "io/chebmom_file.hpp"
#include "io/csr_text.hpp"

namespace {

struct CliOptions {
  std::string fixtureDir;
  std::string label;
  double windowLower = 0.0;
  double windowUpper = 0.0;
  bool hasWindow = false;
  std::string route = "lowdin";
  std::string applier = "exact";
  int order = 60;
  double coeffTol = 1e-12;
  double tol = 1e-13;
  int maxIter = 200;
  lsquant::Index maxDim = 4096;
  std::string traceMode = "exact";
  unsigned seed = 12345;
  bool hasSeed = false;
  lsquant::Index numStates = -1;
  int numMoments = 64;
  std::string word;
  int batchSize = -1;
  std::string outDir = ".";
  bool checkAppliers = false;
  bool showHelp = false;
};

void PrintUsage(std::ostream& os) {
  os << "Usage: stage6_nonorth --fixture <dir> --label <model> --window <lower> <upper>\n"
     << "                      [--route lowdin|generalised|naive]\n"
     << "                      [--applier exact|ischo|chebyshev]\n"
     << "                      [--order <P>] [--coeff-tol <ctol>] [--tol <t>]\n"
     << "                      [--max-iter <n>] [--max-dim <Nmax>]\n"
     << "                      [--trace exact|random] [--seed <s>] [--states <R>]\n"
     << "                      [--moments <M>] [--batch <B>] --word dos|kubo|kuboxy\n"
     << "                      [--out <dir>] [--check-appliers]\n\n"
     << "Routes:\n"
     << "  lowdin       Symmetric Lowdin sandwich H_L = S^-1/2 H' S^-1/2, V_L = S^-1/2 V' S^-1/2.\n"
     << "  generalised  Recursion on S^-1 (a H' + b S') with S^-1 V' on ket and adjoint\n"
     << "               bra S T_m(S^-1 H~) S^-1 V' S^-1 |r> (theory.md Eqs. 2a-2c).\n"
     << "  naive        Legacy estimator (wrong for Kubo when S != I, theory.md section 4.1):\n"
     << "               recursion on S^-1 (a H' + b S') (GeneralisedRescaledOperator) with\n"
     << "               S^-1 V' (MetricInverseProduct(S^-1, V')) on both sides and no S on\n"
     << "               the bra (the symmetric stage-1 plan); passes on dos and fails on\n"
     << "               kubo/kuboxy by O(1) because the trace estimator is not covariant.\n";
}

CliOptions ParseArgs(int argc, char* argv[]) {
  CliOptions opts;

  for (int i = 1; i < argc; ++i) {
    const std::string arg = argv[i];
    if (arg == "--help" || arg == "-h") {
      opts.showHelp = true;
      return opts;
    }
    auto requireNext = [&](const char* flag) -> std::string {
      if (i + 1 >= argc) {
        throw std::invalid_argument(std::string("stage6_nonorth: missing argument for ") + flag);
      }
      return argv[++i];
    };

    if (arg == "--fixture" || arg == "--operators") {
      opts.fixtureDir = requireNext(arg.c_str());
    } else if (arg == "--label") {
      opts.label = requireNext("--label");
    } else if (arg == "--window") {
      opts.windowLower = std::stod(requireNext("--window"));
      opts.windowUpper = std::stod(requireNext("--window"));
      opts.hasWindow = true;
    } else if (arg == "--route") {
      opts.route = requireNext("--route");
    } else if (arg == "--applier") {
      opts.applier = requireNext("--applier");
    } else if (arg == "--order") {
      opts.order = std::stoi(requireNext("--order"));
    } else if (arg == "--coeff-tol") {
      opts.coeffTol = std::stod(requireNext("--coeff-tol"));
    } else if (arg == "--tol") {
      opts.tol = std::stod(requireNext("--tol"));
    } else if (arg == "--max-iter") {
      opts.maxIter = std::stoi(requireNext("--max-iter"));
    } else if (arg == "--max-dim") {
      opts.maxDim = std::stoll(requireNext("--max-dim"));
    } else if (arg == "--trace") {
      opts.traceMode = requireNext("--trace");
    } else if (arg == "--seed") {
      opts.seed = static_cast<unsigned>(std::stoul(requireNext("--seed")));
      opts.hasSeed = true;
    } else if (arg == "--states") {
      opts.numStates = std::stoll(requireNext("--states"));
    } else if (arg == "--moments") {
      opts.numMoments = std::stoi(requireNext("--moments"));
    } else if (arg == "--word") {
      opts.word = requireNext("--word");
    } else if (arg == "--batch") {
      opts.batchSize = std::stoi(requireNext("--batch"));
    } else if (arg == "--out") {
      opts.outDir = requireNext("--out");
    } else if (arg == "--check-appliers") {
      opts.checkAppliers = true;
    } else {
      throw std::invalid_argument("stage6_nonorth: unknown option: " + arg);
    }
  }

  if (opts.fixtureDir.empty()) {
    throw std::invalid_argument("stage6_nonorth: --fixture <dir> is required");
  }
  if (opts.label.empty()) {
    throw std::invalid_argument("stage6_nonorth: --label <label> is required");
  }
  if (opts.coeffTol < 0.0 || !std::isfinite(opts.coeffTol)) {
    throw std::invalid_argument("stage6_nonorth: --coeff-tol must be non-negative and finite");
  }
  if (opts.maxDim <= 0) {
    throw std::invalid_argument("stage6_nonorth: --max-dim must be positive");
  }
  if (opts.checkAppliers) {
    return opts;
  }

  if (!opts.hasWindow || !(opts.windowUpper > opts.windowLower)) {
    throw std::invalid_argument("stage6_nonorth: valid --window <lower> <upper> is required");
  }
  if (opts.route != "lowdin" && opts.route != "generalised" && opts.route != "naive") {
    throw std::invalid_argument(
        "stage6_nonorth: --route must be 'lowdin', 'generalised', or 'naive'");
  }
  if (opts.applier != "exact" && opts.applier != "ischo" && opts.applier != "chebyshev") {
    throw std::invalid_argument(
        "stage6_nonorth: --applier must be 'exact', 'ischo', or 'chebyshev'");
  }
  if (opts.route == "lowdin" && opts.applier == "ischo") {
    throw std::invalid_argument(
        "stage6_nonorth: --applier ischo provides S^-1 only and cannot be used with --route "
        "lowdin (use 'exact' or 'chebyshev')");
  }
  if (opts.traceMode != "exact" && opts.traceMode != "random") {
    throw std::invalid_argument("stage6_nonorth: --trace must be 'exact' or 'random'");
  }
  if (opts.numMoments <= 0) {
    throw std::invalid_argument("stage6_nonorth: --moments must be positive");
  }
  if (opts.word != "dos" && opts.word != "kubo" && opts.word != "kuboxy") {
    throw std::invalid_argument("stage6_nonorth: --word must be 'dos', 'kubo', or 'kuboxy'");
  }
  if (opts.batchSize <= 0) {
    opts.batchSize = opts.numMoments;
  }
  return opts;
}

double RelativeDifference(const lsquant::Vector& a, const lsquant::Vector& b) {
  lsquant::Vector diff = a;
  lsquant::linalg::axpy(lsquant::Scalar(-1.0, 0.0), b, diff);
  const double diffNorm = std::sqrt(lsquant::linalg::norm2(diff));
  const double denom = std::sqrt(lsquant::linalg::norm2(b));
  return denom > 0.0 ? (diffNorm / denom) : diffNorm;
}

lsquant::SpectralBounds ComputePaddedOverlapBounds(const lsquant::LinearOperator& sOp) {
  const lsquant::SpectralBounds ritz = lsquant::LanczosBounds(sOp, 80, 12345u);
  const double span = std::max(1e-6, ritz.upper - ritz.lower);
  const double lower = std::max(0.5 * ritz.lower, ritz.lower - 0.05 * span);
  const double upper = ritz.upper + 0.05 * span;
  return {lower, upper};
}

int RunApplierChecks(const CliOptions& opts) {
  const std::string sPath = opts.fixtureDir + "/" + opts.label + ".S.CSR";
  const lsquant::SparseMatrix sMat = lsquant::ReadCsrText(sPath);
  const lsquant::SparseMatrixOperator sOp(sMat);
  const lsquant::Index dim = sMat.Dimension();
  const lsquant::SpectralBounds boundsOfS = ComputePaddedOverlapBounds(sOp);

  const bool includeExact = (opts.applier == "exact");
  std::unique_ptr<lsquant::ExactInverseOverlap> exactInv;
  std::unique_ptr<lsquant::ExactInverseSqrtOverlap> exactInvSqrt;
  if (includeExact) {
    exactInv = std::make_unique<lsquant::ExactInverseOverlap>(sMat, opts.maxDim);
    exactInvSqrt = std::make_unique<lsquant::ExactInverseSqrtOverlap>(sMat, opts.maxDim);
  }
  const lsquant::IschoInverseOverlap ischoInv(sMat, opts.tol, opts.maxIter);
  const lsquant::ChebyshevInverseOverlap chebInv(sMat, boundsOfS, opts.order, opts.coeffTol);
  const lsquant::ChebyshevInverseSqrtOverlap chebInvSqrt(sMat, boundsOfS, opts.order,
                                                         opts.coeffTol);

  lsquant::RandomPhaseStates rng(dim, 10, opts.seed);
  lsquant::Vector v(static_cast<std::size_t>(dim), lsquant::Scalar(0.0, 0.0));
  lsquant::Vector yExact(static_cast<std::size_t>(dim), lsquant::Scalar(0.0, 0.0));
  lsquant::Vector yIscho(static_cast<std::size_t>(dim), lsquant::Scalar(0.0, 0.0));
  lsquant::Vector yCheb(static_cast<std::size_t>(dim), lsquant::Scalar(0.0, 0.0));
  lsquant::Vector sv(static_cast<std::size_t>(dim), lsquant::Scalar(0.0, 0.0));
  lsquant::Vector tmp(static_cast<std::size_t>(dim), lsquant::Scalar(0.0, 0.0));

  double maxDiffChebVsExact = 0.0;
  double maxDiffIschoVsExact = 0.0;
  double maxDiffSqrtChebVsExact = 0.0;
  double maxResExact = 0.0;
  double maxResIscho = 0.0;
  double maxResCheb = 0.0;
  double maxResSqrtExact = 0.0;
  double maxResSqrtCheb = 0.0;

  while (rng.Next(v)) {
    ischoInv.Multiply(v, yIscho);
    chebInv.Multiply(v, yCheb);

    if (includeExact) {
      exactInv->Multiply(v, yExact);
      maxDiffChebVsExact = std::max(maxDiffChebVsExact, RelativeDifference(yCheb, yExact));
      maxDiffIschoVsExact = std::max(maxDiffIschoVsExact, RelativeDifference(yIscho, yExact));
      sOp.Multiply(yExact, sv);
      maxResExact = std::max(maxResExact, RelativeDifference(sv, v));
    }

    sOp.Multiply(yIscho, sv);
    maxResIscho = std::max(maxResIscho, RelativeDifference(sv, v));

    sOp.Multiply(yCheb, sv);
    maxResCheb = std::max(maxResCheb, RelativeDifference(sv, v));

    chebInvSqrt.Multiply(v, yCheb);
    if (includeExact) {
      exactInvSqrt->Multiply(v, yExact);
      maxDiffSqrtChebVsExact =
          std::max(maxDiffSqrtChebVsExact, RelativeDifference(yCheb, yExact));
    }

    // Check ||(S^-1/2)^2 S v - v|| / ||v||.
    sOp.Multiply(v, sv);
    if (includeExact) {
      exactInvSqrt->Multiply(sv, tmp);
      exactInvSqrt->Multiply(tmp, yExact);
      maxResSqrtExact = std::max(maxResSqrtExact, RelativeDifference(yExact, v));
    }

    chebInvSqrt.Multiply(sv, tmp);
    chebInvSqrt.Multiply(tmp, yCheb);
    maxResSqrtCheb = std::max(maxResSqrtCheb, RelativeDifference(yCheb, v));
  }

  std::cout << std::scientific << std::setprecision(16)
            << "stage6_nonorth: check_appliers label=" << opts.label << " dim=" << dim
            << " order=" << opts.order << " coeff_tol=" << opts.coeffTol << " tol=" << opts.tol
            << " max_iterations=" << opts.maxIter << " max_dim=" << opts.maxDim
            << " max_rel_diff_sinv_cheb_vs_exact=" << maxDiffChebVsExact
            << " max_rel_diff_sinv_ischo_vs_exact=" << maxDiffIschoVsExact
            << " max_rel_diff_sinvsqrt_cheb_vs_exact=" << maxDiffSqrtChebVsExact
            << " max_res_sinv_exact=" << maxResExact
            << " max_res_sinv_ischo=" << maxResIscho
            << " max_res_sinv_cheb=" << maxResCheb
            << " max_res_sinvsqrt_exact=" << maxResSqrtExact
            << " max_res_sinvsqrt_cheb=" << maxResSqrtCheb << "\n";
  return 0;
}

lsquant::Plan BuildDosPlan(int numMoments, lsquant::Index numStates) {
  lsquant::Plan plan;
  plan.word = "1 delta";
  plan.code.push_back(lsquant::Expand{numMoments});
  std::vector<lsquant::Scalar> weights(static_cast<std::size_t>(numMoments));
  for (int m = 0; m < numMoments; ++m) {
    double scal = 2.0 / static_cast<double>(numStates);
    if (m == 0) {
      scal *= 0.5;
    }
    weights[static_cast<std::size_t>(m)] = lsquant::Scalar(scal, 0.0);
  }
  plan.code.push_back(lsquant::Trace{std::move(weights)});
  return plan;
}

lsquant::Plan BuildHermitianKuboPlan(const std::string& wordName, int braOpIdx, int ketOpIdx,
                                     int numMoments) {
  lsquant::Plan plan;
  plan.word = wordName;
  plan.code.push_back(lsquant::Expand{numMoments});
  plan.code.push_back(lsquant::ApplyOperator{ketOpIdx});
  plan.code.push_back(lsquant::MirrorToBra{});
  plan.code.push_back(lsquant::ApplyOperator{braOpIdx});
  plan.code.push_back(lsquant::Expand{numMoments});
  plan.code.push_back(lsquant::Trace{{lsquant::Scalar(1.0, 0.0)}});
  return plan;
}

lsquant::Plan BuildGeneralisedKuboPlan(const std::string& wordName, int sinvBraOpIdx,
                                       int sinvKetOpIdx, int sinvIdx, int sIdx,
                                       int numMoments) {
  lsquant::Plan plan;
  plan.word = wordName;
  plan.code.push_back(lsquant::Expand{numMoments});
  plan.code.push_back(lsquant::ApplyOperator{sinvKetOpIdx});
  plan.code.push_back(lsquant::MirrorToBra{});
  plan.code.push_back(lsquant::ApplyOperator{sinvIdx});
  plan.code.push_back(lsquant::ApplyOperator{sinvBraOpIdx});
  plan.code.push_back(lsquant::Expand{numMoments});
  plan.code.push_back(lsquant::ApplyOperator{sIdx});
  plan.code.push_back(lsquant::Trace{{lsquant::Scalar(1.0, 0.0)}});
  return plan;
}

}  // namespace

int main(int argc, char* argv[]) {
  try {
    const CliOptions opts = ParseArgs(argc, argv);
    if (opts.showHelp) {
      PrintUsage(std::cout);
      return 0;
    }
    if (opts.checkAppliers) {
      return RunApplierChecks(opts);
    }

    const std::string hamPath = opts.fixtureDir + "/" + opts.label + ".HAM.CSR";
    const lsquant::SparseMatrix ham = lsquant::ReadCsrText(hamPath);
    const lsquant::SparseMatrixOperator hamOp(ham);
    const lsquant::SpectralWindow window{opts.windowLower, opts.windowUpper};
    const lsquant::Index dim = ham.Dimension();

    const lsquant::Index numStates =
        (opts.numStates > 0) ? opts.numStates : (opts.traceMode == "exact" ? dim : 1);

    std::unique_ptr<lsquant::StateGenerator> stateGen;
    std::string stateLabel;
    if (opts.traceMode == "exact") {
      if (numStates > dim) {
        throw std::invalid_argument("stage6_nonorth: --states exceeds operator dimension");
      }
      std::vector<lsquant::Index> sites(static_cast<std::size_t>(numStates));
      for (lsquant::Index s = 0; s < numStates; ++s) {
        sites[static_cast<std::size_t>(s)] = s;
      }
      stateGen = std::make_unique<lsquant::LocalStates>(dim, std::move(sites));
      stateLabel = "exact";
    } else {
      stateGen = std::make_unique<lsquant::RandomPhaseStates>(dim, numStates, opts.seed);
      stateLabel = stateGen->Label();
    }

    std::optional<lsquant::SparseMatrix> sMat;
    std::unique_ptr<lsquant::LinearOperator> inverseApplier;
    std::unique_ptr<lsquant::LinearOperator> inverseSqrtApplier;
    std::unique_ptr<lsquant::MetricPolicy> metric;
    const lsquant::IschoInverseOverlap* ischoPtr = nullptr;
    int effectiveOrder = -1;

    const std::string sPath = opts.fixtureDir + "/" + opts.label + ".S.CSR";
    const bool hasOverlapFile = std::ifstream(sPath).good();

    if (opts.route == "naive" && !hasOverlapFile) {
      metric = std::make_unique<lsquant::OrthogonalMetric>(dim);
    } else {
      sMat = lsquant::ReadCsrText(sPath);
      const lsquant::SparseMatrixOperator sOp(*sMat);

      if (opts.route == "lowdin") {
        inverseApplier = std::make_unique<lsquant::IdentityOperator>(dim);
        if (opts.applier == "exact") {
          inverseSqrtApplier = std::make_unique<lsquant::ExactInverseSqrtOverlap>(*sMat, opts.maxDim);
        } else {
          const lsquant::SpectralBounds boundsOfS = ComputePaddedOverlapBounds(sOp);
          auto chebSqrt = std::make_unique<lsquant::ChebyshevInverseSqrtOverlap>(
              *sMat, boundsOfS, opts.order, opts.coeffTol);
          effectiveOrder = chebSqrt->EffectiveOrder();
          inverseSqrtApplier = std::move(chebSqrt);
        }
        metric = std::make_unique<lsquant::NonOrthogonalMetric>(
            *sMat, *inverseApplier, *inverseSqrtApplier);
      } else {
        if (opts.applier == "exact") {
          inverseApplier = std::make_unique<lsquant::ExactInverseOverlap>(*sMat, opts.maxDim);
        } else if (opts.applier == "ischo") {
          auto ischo =
              std::make_unique<lsquant::IschoInverseOverlap>(*sMat, opts.tol, opts.maxIter);
          ischoPtr = ischo.get();
          inverseApplier = std::move(ischo);
        } else {
          const lsquant::SpectralBounds boundsOfS = ComputePaddedOverlapBounds(sOp);
          auto cheb = std::make_unique<lsquant::ChebyshevInverseOverlap>(
              *sMat, boundsOfS, opts.order, opts.coeffTol);
          effectiveOrder = cheb->EffectiveOrder();
          inverseApplier = std::move(cheb);
        }
        inverseSqrtApplier = std::make_unique<lsquant::IdentityOperator>(dim);
        metric = std::make_unique<lsquant::NonOrthogonalMetric>(
            *sMat, *inverseApplier, *inverseSqrtApplier);
      }
    }

    std::unique_ptr<lsquant::LowdinOperator> lowdinH;
    std::unique_ptr<lsquant::LinearOperator> rescaledH;
    if (metric->IsOrthogonal()) {
      rescaledH = std::make_unique<lsquant::RescaledOperator>(ham, window);
    } else if (opts.route == "lowdin") {
      lowdinH = std::make_unique<lsquant::LowdinOperator>(metric->InverseSqrtOverlap(), hamOp);
      rescaledH = std::make_unique<lsquant::RescaledOperator>(*lowdinH, window);
    } else {
      rescaledH = std::make_unique<lsquant::GeneralisedRescaledOperator>(
          ham, *sMat, metric->InverseOverlap(), window);
    }

    lsquant::ChebyshevMoments mu;
    std::string legacyPath;
    std::string newPath;

    if (opts.word == "dos") {
      lsquant::MomentEngine engine(*rescaledH, {}, *stateGen, opts.batchSize);
      const lsquant::Plan plan = BuildDosPlan(opts.numMoments, stateGen->NumberOfStates());
      mu = engine.Run(plan);
      mu.window = window;
      mu.systemSize = dim;
      mu.numberOfStates = static_cast<int>(stateGen->NumberOfStates());
      mu.stateLabel = stateLabel;
      mu.word = plan.word;
      mu.kpmWeighted = true;
      mu.symmetrised = false;

      const std::string stem = opts.outDir + "/SpectralOp1" + opts.label + "KPM_M" +
                               std::to_string(opts.numMoments) + "_state" + stateLabel;
      legacyPath = stem + ".chebmom1D";
      newPath = stem + ".chebmom";
      lsquant::WriteLegacyChebmom1D(legacyPath, mu);
      lsquant::WriteChebmom(newPath, mu);
    } else if (opts.word == "kubo") {
      const std::string vxPath = opts.fixtureDir + "/" + opts.label + ".VX.CSR";
      const lsquant::SparseMatrix vxMat = lsquant::ReadCsrText(vxPath);
      const lsquant::SparseMatrixOperator vxOp(vxMat);

      if (opts.route == "naive") {
        if (metric->IsOrthogonal()) {
          lsquant::MomentEngine engine(*rescaledH, {&vxOp}, *stateGen, opts.batchSize);
          const lsquant::Plan plan =
              BuildHermitianKuboPlan("VX delta VX delta", 0, 0, opts.numMoments);
          mu = engine.Run(plan);
          mu.word = plan.word;
        } else {
          const lsquant::MetricInverseProduct sinvVx(metric->InverseOverlap(), vxOp);
          lsquant::MomentEngine engine(*rescaledH, {&sinvVx}, *stateGen, opts.batchSize);
          const lsquant::Plan plan =
              BuildHermitianKuboPlan("VX delta VX delta", 0, 0, opts.numMoments);
          mu = engine.Run(plan);
          mu.word = plan.word;
        }
      } else if (opts.route == "lowdin") {
        const lsquant::LowdinOperator lowdinVx(metric->InverseSqrtOverlap(), vxOp);
        lsquant::MomentEngine engine(*rescaledH, {&lowdinVx}, *stateGen, opts.batchSize);
        const lsquant::Plan plan =
            BuildHermitianKuboPlan("VX delta VX delta", 0, 0, opts.numMoments);
        mu = engine.Run(plan);
        mu.word = plan.word;
      } else {
        const lsquant::MetricInverseProduct sinvVx(metric->InverseOverlap(), vxOp);
        lsquant::MomentEngine engine(
            *rescaledH, {&sinvVx, &metric->InverseOverlap(), &metric->Overlap()}, *stateGen,
            opts.batchSize);
        const lsquant::Plan plan =
            BuildGeneralisedKuboPlan("VX delta VX delta", 0, 0, 1, 2, opts.numMoments);
        mu = engine.Run(plan);
        mu.word = plan.word;
      }

      mu.window = window;
      mu.systemSize = dim;
      mu.numberOfStates = static_cast<int>(stateGen->NumberOfStates());
      mu.stateLabel = stateLabel;
      lsquant::ApplyKpmWeights(mu, static_cast<int>(stateGen->NumberOfStates()), true);

      const std::string stem = opts.outDir + "/NonEqOpVX-VX" + opts.label + "KPM_M" +
                               std::to_string(opts.numMoments) + "x" +
                               std::to_string(opts.numMoments) + "_state" + stateLabel;
      legacyPath = stem + ".chebmom2D";
      newPath = stem + ".chebmom";
      lsquant::WriteLegacyChebmom2D(legacyPath, mu);
      lsquant::WriteChebmom(newPath, mu);
    } else {
      const std::string vxPath = opts.fixtureDir + "/" + opts.label + ".VX.CSR";
      const std::string vyPath = opts.fixtureDir + "/" + opts.label + ".VY.CSR";
      const lsquant::SparseMatrix vxMat = lsquant::ReadCsrText(vxPath);
      const lsquant::SparseMatrix vyMat = lsquant::ReadCsrText(vyPath);
      const lsquant::SparseMatrixOperator vxOp(vxMat);
      const lsquant::SparseMatrixOperator vyOp(vyMat);

      if (opts.route == "naive") {
        if (metric->IsOrthogonal()) {
          lsquant::MomentEngine engine(*rescaledH, {&vxOp, &vyOp}, *stateGen, opts.batchSize);
          const lsquant::Plan plan =
              BuildHermitianKuboPlan("VX delta VY delta", 0, 1, opts.numMoments);
          mu = engine.Run(plan);
          mu.word = plan.word;
        } else {
          const lsquant::MetricInverseProduct sinvVx(metric->InverseOverlap(), vxOp);
          const lsquant::MetricInverseProduct sinvVy(metric->InverseOverlap(), vyOp);
          lsquant::MomentEngine engine(*rescaledH, {&sinvVx, &sinvVy}, *stateGen, opts.batchSize);
          const lsquant::Plan plan =
              BuildHermitianKuboPlan("VX delta VY delta", 0, 1, opts.numMoments);
          mu = engine.Run(plan);
          mu.word = plan.word;
        }
      } else if (opts.route == "lowdin") {
        const lsquant::LowdinOperator lowdinVx(metric->InverseSqrtOverlap(), vxOp);
        const lsquant::LowdinOperator lowdinVy(metric->InverseSqrtOverlap(), vyOp);
        lsquant::MomentEngine engine(*rescaledH, {&lowdinVx, &lowdinVy}, *stateGen,
                                     opts.batchSize);
        const lsquant::Plan plan =
            BuildHermitianKuboPlan("VX delta VY delta", 0, 1, opts.numMoments);
        mu = engine.Run(plan);
        mu.word = plan.word;
      } else {
        const lsquant::MetricInverseProduct sinvVx(metric->InverseOverlap(), vxOp);
        const lsquant::MetricInverseProduct sinvVy(metric->InverseOverlap(), vyOp);
        lsquant::MomentEngine engine(
            *rescaledH, {&sinvVx, &sinvVy, &metric->InverseOverlap(), &metric->Overlap()},
            *stateGen, opts.batchSize);
        const lsquant::Plan plan =
            BuildGeneralisedKuboPlan("VX delta VY delta", 0, 1, 2, 3, opts.numMoments);
        mu = engine.Run(plan);
        mu.word = plan.word;
      }

      mu.window = window;
      mu.systemSize = dim;
      mu.numberOfStates = static_cast<int>(stateGen->NumberOfStates());
      mu.stateLabel = stateLabel;
      lsquant::ApplyKpmWeights(mu, static_cast<int>(stateGen->NumberOfStates()), true);

      const std::string stem = opts.outDir + "/NonEqOpVX-VY" + opts.label + "KPM_M" +
                               std::to_string(opts.numMoments) + "x" +
                               std::to_string(opts.numMoments) + "_state" + stateLabel;
      legacyPath = stem + ".chebmom2D";
      newPath = stem + ".chebmom";
      lsquant::WriteLegacyChebmom2D(legacyPath, mu);
      lsquant::WriteChebmom(newPath, mu);
    }

    double residual = 0.0;
    if (ischoPtr != nullptr) {
      residual = ischoPtr->LastError();
    } else if (!metric->IsOrthogonal()) {
      lsquant::RandomPhaseStates rng(dim, 1, opts.seed);
      lsquant::Vector v(static_cast<std::size_t>(dim), lsquant::Scalar(0.0, 0.0));
      lsquant::Vector sinvV(static_cast<std::size_t>(dim), lsquant::Scalar(0.0, 0.0));
      lsquant::Vector sSinvV(static_cast<std::size_t>(dim), lsquant::Scalar(0.0, 0.0));
      rng.Next(v);
      if (opts.route == "lowdin") {
        lsquant::Vector tmp(static_cast<std::size_t>(dim), lsquant::Scalar(0.0, 0.0));
        metric->InverseSqrtOverlap().Multiply(v, tmp);
        metric->InverseSqrtOverlap().Multiply(tmp, sinvV);
      } else {
        metric->InverseOverlap().Multiply(v, sinvV);
      }
      metric->Overlap().Multiply(sinvV, sSinvV);
      residual = RelativeDifference(sSinvV, v);
    }

    std::cout << "stage6_nonorth: route=" << opts.route << " applier=" << opts.applier
              << " word=" << opts.word << " dim=" << dim << " states=" << mu.numberOfStates
              << " (" << stateLabel << ") M=" << opts.numMoments << " -> " << legacyPath
              << " order=" << opts.order;
    if (effectiveOrder >= 0) {
      std::cout << " effective_order=" << effectiveOrder;
    }
    std::cout << std::scientific << std::setprecision(6) << " coeff_tol=" << opts.coeffTol
              << " tol=" << opts.tol << std::defaultfloat << " max_iterations=" << opts.maxIter
              << " max_dim=" << opts.maxDim;
    if (ischoPtr != nullptr) {
      std::cout << " iterations=" << ischoPtr->LastIterations();
    }
    std::cout << std::scientific << std::setprecision(6) << " residual=" << residual << "\n";
    return 0;
  } catch (const std::exception& ex) {
    std::cerr << "stage6_nonorth error: " << ex.what() << "\n";
    return 1;
  }
}
