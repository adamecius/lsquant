/**
 * @file chebyshev_moments.cpp
 * @brief Implementation of ChebyshevMoments and ApplyKpmWeights.
 */
#include "core/chebyshev/chebyshev_moments.hpp"

#include <complex>
#include <stdexcept>
#include <utility>

namespace lsquant {

ChebyshevMoments::ChebyshevMoments(std::vector<int> dims) : dims_(std::move(dims)) {
  std::size_t total = 1;
  for (int d : dims_) {
    if (d <= 0) {
      throw std::invalid_argument("ChebyshevMoments: dimensions must be positive");
    }
    total *= static_cast<std::size_t>(d);
  }
  data_.assign(total, Scalar(0.0, 0.0));
}

int ChebyshevMoments::Rank() const { return static_cast<int>(dims_.size()); }

const std::vector<int>& ChebyshevMoments::Dims() const { return dims_; }

std::size_t ChebyshevMoments::FlatIndex(const std::vector<int>& index) const {
  if (index.size() != dims_.size()) {
    throw std::invalid_argument("ChebyshevMoments::At: index rank mismatch");
  }
  std::size_t flat = 0;
  for (std::size_t k = 0; k < dims_.size(); ++k) {
    const int coord = index[k];
    const int bound = dims_[k];
    if (coord < 0 || coord >= bound) {
      throw std::out_of_range("ChebyshevMoments::At: index out of range");
    }
    flat = flat * static_cast<std::size_t>(bound) + static_cast<std::size_t>(coord);
  }
  return flat;
}

Scalar& ChebyshevMoments::At(const std::vector<int>& index) {
  return data_[FlatIndex(index)];
}

const Scalar& ChebyshevMoments::At(const std::vector<int>& index) const {
  return data_[FlatIndex(index)];
}

Scalar& ChebyshevMoments::At(int m) {
  if (dims_.size() != 1) {
    throw std::invalid_argument("ChebyshevMoments::At(m): requires rank 1");
  }
  if (m < 0 || m >= dims_[0]) {
    throw std::out_of_range("ChebyshevMoments::At(m): index out of range");
  }
  return data_[static_cast<std::size_t>(m)];
}

const Scalar& ChebyshevMoments::At(int m) const {
  if (dims_.size() != 1) {
    throw std::invalid_argument("ChebyshevMoments::At(m): requires rank 1");
  }
  if (m < 0 || m >= dims_[0]) {
    throw std::out_of_range("ChebyshevMoments::At(m): index out of range");
  }
  return data_[static_cast<std::size_t>(m)];
}

Scalar& ChebyshevMoments::At(int m, int n) {
  if (dims_.size() != 2) {
    throw std::invalid_argument("ChebyshevMoments::At(m, n): requires rank 2");
  }
  if (m < 0 || m >= dims_[0] || n < 0 || n >= dims_[1]) {
    throw std::out_of_range("ChebyshevMoments::At(m, n): index out of range");
  }
  return data_[static_cast<std::size_t>(m) * static_cast<std::size_t>(dims_[1]) +
               static_cast<std::size_t>(n)];
}

const Scalar& ChebyshevMoments::At(int m, int n) const {
  if (dims_.size() != 2) {
    throw std::invalid_argument("ChebyshevMoments::At(m, n): requires rank 2");
  }
  if (m < 0 || m >= dims_[0] || n < 0 || n >= dims_[1]) {
    throw std::out_of_range("ChebyshevMoments::At(m, n): index out of range");
  }
  return data_[static_cast<std::size_t>(m) * static_cast<std::size_t>(dims_[1]) +
               static_cast<std::size_t>(n)];
}

std::vector<Scalar>& ChebyshevMoments::Data() { return data_; }

const std::vector<Scalar>& ChebyshevMoments::Data() const { return data_; }

void ApplyKpmWeights(ChebyshevMoments& mu, int numberOfStates, bool symmetrise) {
  if (mu.Rank() != 2) {
    throw std::invalid_argument("ApplyKpmWeights: requires a rank-2 ChebyshevMoments tensor");
  }
  if (numberOfStates <= 0) {
    throw std::invalid_argument("ApplyKpmWeights: numberOfStates must be positive");
  }
  if (mu.kpmWeighted) {
    throw std::invalid_argument("ApplyKpmWeights: KPM weights have already been applied");
  }

  const int numMomsL = mu.Dims()[0];
  const int numMomsR = mu.Dims()[1];
  if (symmetrise) {
    if (numMomsL != numMomsR) {
      throw std::invalid_argument("ApplyKpmWeights: symmetrisation requires equal dimensions");
    }
    for (int mL = 0; mL < numMomsL; ++mL) {
      for (int mR = mL; mR < numMomsR; ++mR) {
        double scal = 4.0 / static_cast<double>(numberOfStates);
        if (mL == 0) {
          scal *= 0.5;
        }
        if (mR == 0) {
          scal *= 0.5;
        }
        const Scalar tmp = scal * (mu.At(mL, mR) + std::conj(mu.At(mR, mL))) / 2.0;
        mu.At(mL, mR) = tmp;
        mu.At(mR, mL) = std::conj(tmp);
      }
    }
  } else {
    for (int mL = 0; mL < numMomsL; ++mL) {
      for (int mR = 0; mR < numMomsR; ++mR) {
        double scal = 4.0 / static_cast<double>(numberOfStates);
        if (mL == 0) {
          scal *= 0.5;
        }
        if (mR == 0) {
          scal *= 0.5;
        }
        mu.At(mL, mR) = scal * mu.At(mL, mR);
      }
    }
  }

  mu.kpmWeighted = true;
  mu.symmetrised = symmetrise;
}

}  // namespace lsquant
