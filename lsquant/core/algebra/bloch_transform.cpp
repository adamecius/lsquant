/**
 * @file bloch_transform.cpp
 * @brief Implementation of BlochTransform using batched 3D FFTW plans.
 */
#include "core/algebra/bloch_transform.hpp"
#include "core/algebra/fftw_planner.hpp"

#include <fftw3.h>

#include <cmath>
#include <complex>
#include <mutex>
#include <stdexcept>
#include <utility>

namespace lsquant {

namespace {

std::mutex& FftwPlannerMutex() {
  return detail::FftwPlannerMutex();
}

}  // namespace

struct BlochTransform::FftwWorkspace {
  fftw_complex* inBuffer = nullptr;
  fftw_complex* outBuffer = nullptr;
  fftw_plan plan = nullptr;

  FftwWorkspace(int n1, int n2, int n3, int howmany, int sign) {
    const std::size_t total = static_cast<std::size_t>(n1) *
                              static_cast<std::size_t>(n2) *
                              static_cast<std::size_t>(n3) *
                              static_cast<std::size_t>(howmany);
    std::lock_guard<std::mutex> lock(FftwPlannerMutex());
    inBuffer = fftw_alloc_complex(total);
    outBuffer = fftw_alloc_complex(total);
    if (inBuffer == nullptr || outBuffer == nullptr) {
      if (inBuffer != nullptr) {
        fftw_free(inBuffer);
      }
      if (outBuffer != nullptr) {
        fftw_free(outBuffer);
      }
      throw std::bad_alloc();
    }
    // wannier2sparse cell order has n1 fastest: iR = n1 + n1_ * n2 + n1_ * n2_ * n3,
    // which corresponds to C row-major order on a 3D grid of shape {n3, n2, n1}.
    int dims[3] = {n3, n2, n1};
    plan = fftw_plan_many_dft(3, dims, howmany, inBuffer, nullptr, howmany, 1,
                              outBuffer, nullptr, howmany, 1, sign, FFTW_ESTIMATE);
    if (plan == nullptr) {
      fftw_free(inBuffer);
      fftw_free(outBuffer);
      throw std::runtime_error("BlochTransform: fftw_plan_many_dft failed");
    }
  }

  ~FftwWorkspace() {
    std::lock_guard<std::mutex> lock(FftwPlannerMutex());
    if (plan != nullptr) {
      fftw_destroy_plan(plan);
    }
    if (inBuffer != nullptr) {
      fftw_free(inBuffer);
    }
    if (outBuffer != nullptr) {
      fftw_free(outBuffer);
    }
  }

  FftwWorkspace(const FftwWorkspace&) = delete;
  FftwWorkspace& operator=(const FftwWorkspace&) = delete;
};

BlochTransform::BlochTransform(Index n1, Index n2, Index n3, Index numOrbitals,
                               std::vector<Scalar> atomPhases, Direction direction)
    : n1_(n1),
      n2_(n2),
      n3_(n3),
      numCells_(n1 * n2 * n3),
      numOrbitals_(numOrbitals),
      dimension_(n1 * n2 * n3 * numOrbitals),
      direction_(direction) {
  if (n1 <= 0 || n2 <= 0 || n3 <= 0 || numOrbitals <= 0) {
    throw std::invalid_argument(
        "BlochTransform: grid dimensions and numOrbitals must be positive");
  }
  const std::size_t dim = static_cast<std::size_t>(dimension_);
  const double invSqrtNk = 1.0 / std::sqrt(static_cast<double>(numCells_));

  std::vector<Scalar> scaled(dim);
  if (atomPhases.empty()) {
    std::fill(scaled.begin(), scaled.end(), Scalar(invSqrtNk, 0.0));
  } else {
    if (atomPhases.size() != dim) {
      throw std::invalid_argument(
          "BlochTransform: atomPhases length must equal n1 * n2 * n3 * numOrbitals");
    }
    for (std::size_t i = 0; i < dim; ++i) {
      scaled[i] = atomPhases[i] * invSqrtNk;
    }
  }
  scaledPhases_ = std::make_shared<const std::vector<Scalar>>(std::move(scaled));
  InitWorkspace();
}

BlochTransform BlochTransform::FromScaledPhases(
    Index n1, Index n2, Index n3, Index numOrbitals,
    std::shared_ptr<const std::vector<Scalar>> scaledPhases,
    Direction direction) {
  if (n1 <= 0 || n2 <= 0 || n3 <= 0 || numOrbitals <= 0) {
    throw std::invalid_argument(
        "BlochTransform::FromScaledPhases: grid dimensions and numOrbitals must be positive");
  }
  const std::size_t dim = static_cast<std::size_t>(n1 * n2 * n3 * numOrbitals);
  if (!scaledPhases || scaledPhases->size() != dim) {
    throw std::invalid_argument(
        "BlochTransform::FromScaledPhases: scaledPhases length must equal n1 * n2 * n3 * numOrbitals");
  }
  return BlochTransform(ScaledPhasesTag{}, n1, n2, n3, numOrbitals, std::move(scaledPhases),
                        direction);
}

BlochTransform::BlochTransform(ScaledPhasesTag, Index n1, Index n2, Index n3, Index numOrbitals,
                               std::shared_ptr<const std::vector<Scalar>> scaledPhases,
                               Direction direction)
    : n1_(n1),
      n2_(n2),
      n3_(n3),
      numCells_(n1 * n2 * n3),
      numOrbitals_(numOrbitals),
      dimension_(n1 * n2 * n3 * numOrbitals),
      direction_(direction),
      scaledPhases_(std::move(scaledPhases)) {
  InitWorkspace();
}

BlochTransform::BlochTransform(const BlochTransform& other)
    : n1_(other.n1_),
      n2_(other.n2_),
      n3_(other.n3_),
      numCells_(other.numCells_),
      numOrbitals_(other.numOrbitals_),
      dimension_(other.dimension_),
      direction_(other.direction_),
      scaledPhases_(other.scaledPhases_) {
  InitWorkspace();
}

BlochTransform& BlochTransform::operator=(const BlochTransform& other) {
  if (this != &other) {
    n1_ = other.n1_;
    n2_ = other.n2_;
    n3_ = other.n3_;
    numCells_ = other.numCells_;
    numOrbitals_ = other.numOrbitals_;
    dimension_ = other.dimension_;
    direction_ = other.direction_;
    scaledPhases_ = other.scaledPhases_;
    InitWorkspace();
  }
  return *this;
}

void BlochTransform::InitWorkspace() {
  const int sign = (direction_ == Direction::RealToBloch) ? FFTW_FORWARD : FFTW_BACKWARD;
  workspace_ = std::make_shared<FftwWorkspace>(
      static_cast<int>(n1_), static_cast<int>(n2_), static_cast<int>(n3_),
      static_cast<int>(numOrbitals_), sign);
}

void BlochTransform::Multiply(Scalar a, const Vector& x, Scalar b, Vector& y) const {
  if (static_cast<Index>(x.size()) != dimension_ ||
      static_cast<Index>(y.size()) != dimension_) {
    throw std::invalid_argument("BlochTransform::Multiply: vector dimension mismatch");
  }

  auto* inPtr = reinterpret_cast<Scalar*>(workspace_->inBuffer);
  const auto* outPtr = reinterpret_cast<const Scalar*>(workspace_->outBuffer);
  const Scalar* phases = scaledPhases_->data();
  const Scalar* xPtr = x.data();
  Scalar* yPtr = y.data();
  const Index dim = dimension_;

  if (direction_ == Direction::RealToBloch) {
    #pragma omp parallel for schedule(static)
    for (Index i = 0; i < dim; ++i) {
      inPtr[i] = xPtr[i];
    }
    fftw_execute(workspace_->plan);
    #pragma omp parallel for schedule(static)
    for (Index i = 0; i < dim; ++i) {
      yPtr[i] = a * (outPtr[i] * phases[i]) + b * yPtr[i];
    }
  } else {
    #pragma omp parallel for schedule(static)
    for (Index i = 0; i < dim; ++i) {
      inPtr[i] = std::conj(phases[i]) * xPtr[i];
    }
    fftw_execute(workspace_->plan);
    #pragma omp parallel for schedule(static)
    for (Index i = 0; i < dim; ++i) {
      yPtr[i] = a * outPtr[i] + b * yPtr[i];
    }
  }
}

BlochTransform BlochTransform::Adjoint() const {
  const Direction adjDir = (direction_ == Direction::RealToBloch)
                               ? Direction::BlochToReal
                               : Direction::RealToBloch;
  return FromScaledPhases(n1_, n2_, n3_, numOrbitals_, scaledPhases_, adjDir);
}

}  // namespace lsquant
