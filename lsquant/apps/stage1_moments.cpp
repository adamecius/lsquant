/**
 * @file stage1_moments.cpp
 * @brief Stage-1 gate executable evaluating `dos` and `kubo` moment plans on CSR operators.
 *
 * Parses command-line arguments specifying the operator directory, model label, spectral
 * window, trace generator (`exact` or `random`), moment count `M`, word (`dos` or `kubo`),
 * batch size, and output directory. Executes the corresponding bytecode `Plan` with
 * `MomentEngine` and writes both the byte-exact legacy `.chebmom1D`/`.chebmom2D` file and
 * the unified `.chebmom` file.
 */
#include <cstdlib>
#include <exception>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

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
  std::string operatorsDir;
  std::string label;
  double windowLower = 0.0;
  double windowUpper = 0.0;
  bool hasWindow = false;
  std::string traceMode = "exact";
  unsigned seed = 12345;
  bool hasSeed = false;
  lsquant::Index numStates = -1;
  int numMoments = 64;
  std::string word;
  int batchSize = -1;
  std::string outDir = ".";
};

CliOptions ParseArgs(int argc, char* argv[]) {
  CliOptions opts;

  for (int i = 1; i < argc; ++i) {
    const std::string arg = argv[i];
    auto requireNext = [&](const char* flag) -> std::string {
      if (i + 1 >= argc) {
        throw std::invalid_argument(std::string("stage1_moments: missing argument for ") + flag);
      }
      return argv[++i];
    };

    if (arg == "--operators") {
      opts.operatorsDir = requireNext("--operators");
    } else if (arg == "--label") {
      opts.label = requireNext("--label");
    } else if (arg == "--window") {
      opts.windowLower = std::stod(requireNext("--window"));
      opts.windowUpper = std::stod(requireNext("--window"));
      opts.hasWindow = true;
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
    } else {
      throw std::invalid_argument("stage1_moments: unknown option: " + arg);
    }
  }

  if (opts.operatorsDir.empty()) {
    throw std::invalid_argument("stage1_moments: --operators <dir> is required");
  }
  if (opts.label.empty()) {
    throw std::invalid_argument("stage1_moments: --label <label> is required");
  }
  if (!opts.hasWindow || !(opts.windowUpper > opts.windowLower)) {
    throw std::invalid_argument("stage1_moments: valid --window <lower> <upper> is required");
  }
  if (opts.traceMode != "exact" && opts.traceMode != "random") {
    throw std::invalid_argument("stage1_moments: --trace must be 'exact' or 'random'");
  }
  if (opts.numMoments <= 0) {
    throw std::invalid_argument("stage1_moments: --moments must be positive");
  }
  if (opts.word != "dos" && opts.word != "kubo") {
    throw std::invalid_argument("stage1_moments: --word must be 'dos' or 'kubo'");
  }
  if (opts.batchSize <= 0) {
    opts.batchSize = opts.numMoments;
  }
  return opts;
}

}  // namespace

int main(int argc, char* argv[]) {
  try {
    const CliOptions opts = ParseArgs(argc, argv);

    const std::string hamPath = opts.operatorsDir + "/" + opts.label + ".HAM.CSR";
    const lsquant::SparseMatrix ham = lsquant::ReadCsrText(hamPath);
    const lsquant::SpectralWindow window{opts.windowLower, opts.windowUpper};
    const lsquant::RescaledOperator rescaledH(ham, window);

    const lsquant::Index dim = ham.Dimension();
    const lsquant::Index numStates =
        (opts.numStates > 0) ? opts.numStates : (opts.traceMode == "exact" ? dim : 1);

    std::unique_ptr<lsquant::StateGenerator> stateGen;
    std::string stateLabel;
    if (opts.traceMode == "exact") {
      if (numStates > dim) {
        throw std::invalid_argument("stage1_moments: --states exceeds operator dimension");
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

    if (opts.word == "dos") {
      lsquant::MomentEngine engine(rescaledH, {}, *stateGen, opts.batchSize);
      lsquant::Plan plan;
      plan.word = "1 delta";
      plan.code.push_back(lsquant::Expand{opts.numMoments});

      std::vector<lsquant::Scalar> weights(static_cast<std::size_t>(opts.numMoments));
      for (int m = 0; m < opts.numMoments; ++m) {
        double scal = 2.0 / static_cast<double>(stateGen->NumberOfStates());
        if (m == 0) {
          scal *= 0.5;
        }
        weights[static_cast<std::size_t>(m)] = lsquant::Scalar(scal, 0.0);
      }
      plan.code.push_back(lsquant::Trace{std::move(weights)});

      lsquant::ChebyshevMoments mu = engine.Run(plan);
      mu.window = window;
      mu.systemSize = dim;
      mu.numberOfStates = static_cast<int>(stateGen->NumberOfStates());
      mu.stateLabel = stateLabel;
      mu.word = plan.word;
      mu.kpmWeighted = true;
      mu.symmetrised = false;

      const std::string stem = opts.outDir + "/SpectralOp1" + opts.label + "KPM_M" +
                               std::to_string(opts.numMoments) + "_state" + stateLabel;
      const std::string legacyPath = stem + ".chebmom1D";
      const std::string newPath = stem + ".chebmom";

      lsquant::WriteLegacyChebmom1D(legacyPath, mu);
      lsquant::WriteChebmom(newPath, mu);
      std::cout << "stage1_moments: word=dos dim=" << dim << " states=" << mu.numberOfStates
                << " (" << stateLabel << ") M=" << opts.numMoments << " -> " << legacyPath
                << "\n";
    } else {
      const std::string vxPath = opts.operatorsDir + "/" + opts.label + ".VX.CSR";
      const lsquant::SparseMatrix vx = lsquant::ReadCsrText(vxPath);
      const lsquant::SparseMatrixOperator vxOp(vx);

      lsquant::MomentEngine engine(rescaledH, {&vxOp}, *stateGen, opts.batchSize);
      lsquant::Plan plan;
      plan.word = "VX delta VX delta";
      plan.code.push_back(lsquant::Expand{opts.numMoments});
      plan.code.push_back(lsquant::ApplyOperator{0});
      plan.code.push_back(lsquant::MirrorToBra{});
      plan.code.push_back(lsquant::ApplyOperator{0});
      plan.code.push_back(lsquant::Expand{opts.numMoments});
      plan.code.push_back(lsquant::Trace{{lsquant::Scalar(1.0, 0.0)}});

      lsquant::ChebyshevMoments mu = engine.Run(plan);
      mu.window = window;
      mu.systemSize = dim;
      mu.numberOfStates = static_cast<int>(stateGen->NumberOfStates());
      mu.stateLabel = stateLabel;
      mu.word = plan.word;
      lsquant::ApplyKpmWeights(mu, static_cast<int>(stateGen->NumberOfStates()), true);

      const std::string stem = opts.outDir + "/NonEqOpVX-VX" + opts.label + "KPM_M" +
                               std::to_string(opts.numMoments) + "x" +
                               std::to_string(opts.numMoments) + "_state" + stateLabel;
      const std::string legacyPath = stem + ".chebmom2D";
      const std::string newPath = stem + ".chebmom";

      lsquant::WriteLegacyChebmom2D(legacyPath, mu);
      lsquant::WriteChebmom(newPath, mu);
      std::cout << "stage1_moments: word=kubo dim=" << dim << " states=" << mu.numberOfStates
                << " (" << stateLabel << ") M=" << opts.numMoments << "x" << opts.numMoments
                << " -> " << legacyPath << "\n";
    }

    return 0;
  } catch (const std::exception& ex) {
    std::cerr << "stage1_moments error: " << ex.what() << "\n";
    return 1;
  }
}
