/**
 * @file stage_1_check.cpp
 * @brief Developer self-check for stage 1: exact-trace `dos` and `kubo` plans vs dense diagonalisation.
 *
 * Constructs a random 64 x 64 complex Hermitian Hamiltonian H and observable A from
 * `std::mt19937(7)` as CSR matrices, evaluates the `dos` and `kubo` plans by exact trace
 * through `MomentEngine`, and compares them against the closed-form eigenbasis sums
 * computed with `Eigen::SelfAdjointEigenSolver`. Also re-executes itself under
 * `OMP_NUM_THREADS=1` and `OMP_NUM_THREADS=8` to verify bit-for-bit thread invariance.
 */
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

#include "core/algebra/sparse_matrix.hpp"
#include "core/algebra/spectral_bounds.hpp"
#include "core/algebra/state_generator.hpp"
#include "core/chebyshev/chebyshev_moments.hpp"
#include "core/chebyshev/moment_engine.hpp"
#include "core/chebyshev/plan.hpp"
#include "core/chebyshev/rescaled_operator.hpp"
#include "core/chebyshev/spectral_window.hpp"

namespace {

struct ModelAndMoments {
  lsquant::ChebyshevMoments dosMoments;
  lsquant::ChebyshevMoments kuboMoments;
  double maxDosError = 0.0;
  double maxKuboError = 0.0;
};

lsquant::SparseMatrix MakeRandomHermitianCsr(lsquant::Index dim, std::mt19937& rng,
                                             Eigen::MatrixXcd& denseOut) {
  std::uniform_real_distribution<double> dist(-1.0, 1.0);
  denseOut = Eigen::MatrixXcd::Zero(dim, dim);
  std::vector<lsquant::Index> rows;
  std::vector<lsquant::Index> cols;
  std::vector<lsquant::Scalar> vals;

  for (lsquant::Index i = 0; i < dim; ++i) {
    const double diagVal = dist(rng);
    denseOut(i, i) = lsquant::Scalar(diagVal, 0.0);
    rows.push_back(i);
    cols.push_back(i);
    vals.push_back(lsquant::Scalar(diagVal, 0.0));
    for (lsquant::Index j = i + 1; j < dim; ++j) {
      if (std::abs(dist(rng)) < 0.35) {
        const lsquant::Scalar v(dist(rng), dist(rng));
        denseOut(i, j) = v;
        denseOut(j, i) = std::conj(v);
        rows.push_back(i);
        cols.push_back(j);
        vals.push_back(v);
        rows.push_back(j);
        cols.push_back(i);
        vals.push_back(std::conj(v));
      }
    }
  }
  return lsquant::SparseMatrix::FromTriplets(dim, rows, cols, vals);
}

ModelAndMoments RunCheck() {
  constexpr lsquant::Index kDim = 64;
  constexpr int kNumMoments = 24;

  std::mt19937 rng(7);
  Eigen::MatrixXcd denseH;
  Eigen::MatrixXcd denseA;
  const lsquant::SparseMatrix csrH = MakeRandomHermitianCsr(kDim, rng, denseH);
  const lsquant::SparseMatrix csrA = MakeRandomHermitianCsr(kDim, rng, denseA);
  const lsquant::SparseMatrixOperator opA(csrA);

  const lsquant::SpectralBounds bounds = lsquant::GershgorinBounds(csrH, 0.10);
  const lsquant::SpectralWindow window{bounds.lower, bounds.upper};
  const lsquant::RescaledOperator rescaledH(csrH, window);

  std::vector<lsquant::Index> allSites(static_cast<std::size_t>(kDim));
  for (lsquant::Index i = 0; i < kDim; ++i) {
    allSites[static_cast<std::size_t>(i)] = i;
  }
  lsquant::LocalStates exactStates(kDim, allSites);

  lsquant::MomentEngine engine(rescaledH, {&opA}, exactStates, kNumMoments);

  // Plan 1: dos
  lsquant::Plan dosPlan;
  dosPlan.word = "1 delta";
  dosPlan.code.push_back(lsquant::Expand{kNumMoments});
  std::vector<lsquant::Scalar> dosWeights(kNumMoments);
  for (int m = 0; m < kNumMoments; ++m) {
    double scal = 2.0 / static_cast<double>(kDim);
    if (m == 0) {
      scal *= 0.5;
    }
    dosWeights[static_cast<std::size_t>(m)] = lsquant::Scalar(scal, 0.0);
  }
  dosPlan.code.push_back(lsquant::Trace{dosWeights});
  lsquant::ChebyshevMoments dosMu = engine.Run(dosPlan);
  dosMu.kpmWeighted = true;

  // Plan 2: kubo
  lsquant::Plan kuboPlan;
  kuboPlan.word = "VX delta VX delta";
  kuboPlan.code.push_back(lsquant::Expand{kNumMoments});
  kuboPlan.code.push_back(lsquant::ApplyOperator{0});
  kuboPlan.code.push_back(lsquant::MirrorToBra{});
  kuboPlan.code.push_back(lsquant::ApplyOperator{0});
  kuboPlan.code.push_back(lsquant::Expand{kNumMoments});
  kuboPlan.code.push_back(lsquant::Trace{{lsquant::Scalar(1.0, 0.0)}});
  lsquant::ChebyshevMoments kuboMu = engine.Run(kuboPlan);
  lsquant::ApplyKpmWeights(kuboMu, static_cast<int>(kDim), true);

  // Dense eigenvalue reference
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> eigSolver(denseH);
  if (eigSolver.info() != Eigen::Success) {
    throw std::runtime_error("stage_1_check: SelfAdjointEigenSolver failed");
  }
  const Eigen::VectorXd evals = eigSolver.eigenvalues();
  const Eigen::MatrixXcd evecs = eigSolver.eigenvectors();
  const Eigen::MatrixXcd aEig = evecs.adjoint() * denseA * evecs;

  // T_m(eps~_k) table of shape (kNumMoments, kDim)
  Eigen::MatrixXd chebTable(kNumMoments, kDim);
  for (lsquant::Index k = 0; k < kDim; ++k) {
    const double x = window.Scale() * evals(k) + window.Shift();
    chebTable(0, k) = 1.0;
    if (kNumMoments > 1) {
      chebTable(1, k) = x;
    }
    for (int m = 2; m < kNumMoments; ++m) {
      chebTable(m, k) = 2.0 * x * chebTable(m - 1, k) - chebTable(m - 2, k);
    }
  }

  double maxDosErr = 0.0;
  for (int m = 0; m < kNumMoments; ++m) {
    double scal = 2.0 / static_cast<double>(kDim);
    if (m == 0) {
      scal *= 0.5;
    }
    double sumT = 0.0;
    for (lsquant::Index k = 0; k < kDim; ++k) {
      sumT += chebTable(m, k);
    }
    const lsquant::Scalar ref(scal * sumT, 0.0);
    maxDosErr = std::max(maxDosErr, std::abs(dosMu.At(m) - ref));
  }

  Eigen::MatrixXcd rawKuboRef = Eigen::MatrixXcd::Zero(kNumMoments, kNumMoments);
  for (int m = 0; m < kNumMoments; ++m) {
    for (int n = 0; n < kNumMoments; ++n) {
      lsquant::Scalar sum(0.0, 0.0);
      for (lsquant::Index k = 0; k < kDim; ++k) {
        for (lsquant::Index l = 0; l < kDim; ++l) {
          sum += aEig(k, l) * chebTable(m, k) * aEig(l, k) * chebTable(n, l);
        }
      }
      rawKuboRef(m, n) = sum;
    }
  }

  double maxKuboErr = 0.0;
  for (int mL = 0; mL < kNumMoments; ++mL) {
    for (int mR = mL; mR < kNumMoments; ++mR) {
      double scal = 4.0 / static_cast<double>(kDim);
      if (mL == 0) {
        scal *= 0.5;
      }
      if (mR == 0) {
        scal *= 0.5;
      }
      const lsquant::Scalar ref =
          scal * (rawKuboRef(mL, mR) + std::conj(rawKuboRef(mR, mL))) / 2.0;
      maxKuboErr = std::max(maxKuboErr, std::abs(kuboMu.At(mL, mR) - ref));
      maxKuboErr = std::max(maxKuboErr, std::abs(kuboMu.At(mR, mL) - std::conj(ref)));
    }
  }

  return ModelAndMoments{std::move(dosMu), std::move(kuboMu), maxDosErr, maxKuboErr};
}

void WriteBinaryDump(const std::string& path, const ModelAndMoments& res) {
  std::ofstream out(path, std::ios::binary);
  if (!out.is_open()) {
    throw std::runtime_error("stage_1_check: cannot open dump file: " + path);
  }
  const auto& d1 = res.dosMoments.Data();
  const auto& d2 = res.kuboMoments.Data();
  out.write(reinterpret_cast<const char*>(d1.data()),
            static_cast<std::streamsize>(d1.size() * sizeof(lsquant::Scalar)));
  out.write(reinterpret_cast<const char*>(d2.data()),
            static_cast<std::streamsize>(d2.size() * sizeof(lsquant::Scalar)));
}

std::vector<char> ReadBinaryDump(const std::string& path) {
  std::ifstream in(path, std::ios::binary | std::ios::ate);
  if (!in.is_open()) {
    throw std::runtime_error("stage_1_check: cannot read dump file: " + path);
  }
  const std::streamsize size = in.tellg();
  in.seekg(0, std::ios::beg);
  std::vector<char> buf(static_cast<std::size_t>(size));
  if (!in.read(buf.data(), size)) {
    throw std::runtime_error("stage_1_check: failed reading dump file: " + path);
  }
  return buf;
}

}  // namespace

int main(int argc, char* argv[]) {
  try {
    const ModelAndMoments res = RunCheck();

    if (argc == 3 && std::string(argv[1]) == "--dump") {
      WriteBinaryDump(argv[2], res);
      return 0;
    }

    const std::string tmp1 = "/tmp/lsquant_stage1_check_t1.bin";
    const std::string tmp8 = "/tmp/lsquant_stage1_check_t8.bin";
    const std::string selfPath = argv[0];

    const std::string cmd1 = "OMP_NUM_THREADS=1 \"" + selfPath + "\" --dump \"" + tmp1 + "\"";
    const std::string cmd8 = "OMP_NUM_THREADS=8 \"" + selfPath + "\" --dump \"" + tmp8 + "\"";

    if (std::system(cmd1.c_str()) != 0 || std::system(cmd8.c_str()) != 0) {
      std::remove(tmp1.c_str());
      std::remove(tmp8.c_str());
      throw std::runtime_error("stage_1_check: sub-invocation for thread comparison failed");
    }

    const std::vector<char> bytes1 = ReadBinaryDump(tmp1);
    const std::vector<char> bytes8 = ReadBinaryDump(tmp8);
    std::remove(tmp1.c_str());
    std::remove(tmp8.c_str());

    const bool bitIdentical =
        (bytes1.size() == bytes8.size()) &&
        (std::memcmp(bytes1.data(), bytes8.data(), bytes1.size()) == 0);

    std::cout << "dos max_abs_error: " << res.maxDosError << "\n";
    std::cout << "kubo max_abs_error: " << res.maxKuboError << "\n";
    std::cout << "OMP_NUM_THREADS=1 vs 8 bit-identical: " << (bitIdentical ? "YES" : "NO")
              << "\n";

    if (res.maxDosError < 1e-11 && res.maxKuboError < 1e-11 && bitIdentical) {
      return 0;
    }
    return 1;
  } catch (const std::exception& ex) {
    std::cerr << "stage_1_check error: " << ex.what() << "\n";
    return 1;
  }
}
