/**
 * @file state_generator.cpp
 * @brief Implementation of RandomPhaseStates, LocalStates, BlockLocalStates,
 *        TransformedStates, and UserStates.
 */
#include "core/algebra/state_generator.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <utility>

namespace lsquant {

RandomPhaseStates::RandomPhaseStates(Index dimension, Index numberOfStates, unsigned seed)
    : dimension_(dimension),
      numberOfStates_(numberOfStates),
      seed_(seed),
      current_(0),
      rng_(static_cast<std::mt19937::result_type>(seed)),
      uni_(0.0, 1.0) {
  if (dimension <= 0) {
    throw std::invalid_argument("RandomPhaseStates: dimension must be positive");
  }
  if (numberOfStates < 0) {
    throw std::invalid_argument("RandomPhaseStates: numberOfStates must be non-negative");
  }
}

Index RandomPhaseStates::Dimension() const { return dimension_; }

Index RandomPhaseStates::NumberOfStates() const { return numberOfStates_; }

bool RandomPhaseStates::Next(Vector& out) {
  if (current_ >= numberOfStates_) {
    return false;
  }
  out.resize(static_cast<std::size_t>(dimension_));
  const double norm = std::sqrt(static_cast<double>(out.size()));
  for (auto& x : out) {
    const double phi = 2.0 * M_PI * uni_(rng_);
    x = Scalar(std::cos(phi) / norm, std::sin(phi) / norm);
  }
  ++current_;
  return true;
}

void RandomPhaseStates::Reset() {
  current_ = 0;
  rng_.seed(static_cast<std::mt19937::result_type>(seed_));
  uni_.reset();
}

std::string RandomPhaseStates::Label() const { return "default"; }

LocalStates::LocalStates(Index dimension, std::vector<Index> sites)
    : dimension_(dimension), sites_(std::move(sites)), current_(0) {
  if (dimension <= 0) {
    throw std::invalid_argument("LocalStates: dimension must be positive");
  }
  for (Index s : sites_) {
    if (s < 0 || s >= dimension) {
      throw std::invalid_argument("LocalStates: site index out of range");
    }
  }
  bool isExact = (static_cast<Index>(sites_.size()) == dimension);
  if (isExact) {
    for (Index i = 0; i < dimension; ++i) {
      if (sites_[static_cast<std::size_t>(i)] != i) {
        isExact = false;
        break;
      }
    }
  }
  label_ = isExact ? "exact" : "local";
}

Index LocalStates::Dimension() const { return dimension_; }

Index LocalStates::NumberOfStates() const { return static_cast<Index>(sites_.size()); }

bool LocalStates::Next(Vector& out) {
  if (current_ >= static_cast<Index>(sites_.size())) {
    return false;
  }
  out.assign(static_cast<std::size_t>(dimension_), Scalar(0.0, 0.0));
  const Index site = sites_[static_cast<std::size_t>(current_)];
  out[static_cast<std::size_t>(site)] = Scalar(1.0, 0.0);
  ++current_;
  return true;
}

void LocalStates::Reset() { current_ = 0; }

std::string LocalStates::Label() const { return label_; }

BlockLocalStates::BlockLocalStates(Index numBlocks, Index blockSize)
    : numBlocks_(numBlocks), blockSize_(blockSize), current_(0) {
  if (numBlocks <= 0 || blockSize <= 0) {
    throw std::invalid_argument(
        "BlockLocalStates: numBlocks and blockSize must be positive");
  }
}

BlockLocalStates::BlockLocalStates(Index n1, Index n2, Index n3, Index blockSize)
    : numBlocks_(n1 * n2 * n3), blockSize_(blockSize), current_(0) {
  if (n1 <= 0 || n2 <= 0 || n3 <= 0 || blockSize <= 0) {
    throw std::invalid_argument(
        "BlockLocalStates: grid dimensions and blockSize must be positive");
  }
}

Index BlockLocalStates::Dimension() const { return numBlocks_ * blockSize_; }

Index BlockLocalStates::NumberOfStates() const { return blockSize_; }

bool BlockLocalStates::Next(Vector& out) {
  if (current_ >= blockSize_) {
    return false;
  }
  const Index dim = Dimension();
  out.assign(static_cast<std::size_t>(dim), Scalar(0.0, 0.0));
  const Scalar amp(1.0 / std::sqrt(static_cast<double>(numBlocks_)), 0.0);
  const std::size_t w = static_cast<std::size_t>(blockSize_);
  const std::size_t r = static_cast<std::size_t>(current_);
  for (Index ik = 0; ik < numBlocks_; ++ik) {
    out[static_cast<std::size_t>(ik) * w + r] = amp;
  }
  ++current_;
  return true;
}

void BlockLocalStates::Reset() { current_ = 0; }

std::string BlockLocalStates::Label() const { return "exact"; }

TransformedStates::TransformedStates(const LinearOperator& transform, StateGenerator& inner)
    : transform_(&transform),
      inner_(&inner),
      scratch_(static_cast<std::size_t>(inner.Dimension()), Scalar(0.0, 0.0)) {
  if (transform.Dimension() != inner.Dimension()) {
    throw std::invalid_argument(
        "TransformedStates: transform and inner state dimensions must match");
  }
}

Index TransformedStates::Dimension() const { return transform_->Dimension(); }

Index TransformedStates::NumberOfStates() const { return inner_->NumberOfStates(); }

bool TransformedStates::Next(Vector& out) {
  if (!inner_->Next(scratch_)) {
    return false;
  }
  out.resize(static_cast<std::size_t>(transform_->Dimension()), Scalar(0.0, 0.0));
  transform_->Multiply(Scalar(1.0, 0.0), scratch_, Scalar(0.0, 0.0), out);
  return true;
}

void TransformedStates::Reset() { inner_->Reset(); }

std::string TransformedStates::Label() const { return inner_->Label(); }

UserStates::UserStates(std::vector<Vector> states)
    : states_(std::move(states)), current_(0) {
  if (states_.empty()) {
    throw std::invalid_argument("UserStates: state list must not be empty");
  }
  dimension_ = static_cast<Index>(states_[0].size());
  if (dimension_ <= 0) {
    throw std::invalid_argument("UserStates: state vector dimension must be positive");
  }
  for (const auto& v : states_) {
    if (static_cast<Index>(v.size()) != dimension_) {
      throw std::invalid_argument("UserStates: inconsistent state vector dimensions");
    }
  }
}

Index UserStates::Dimension() const { return dimension_; }

Index UserStates::NumberOfStates() const { return static_cast<Index>(states_.size()); }

bool UserStates::Next(Vector& out) {
  if (current_ >= static_cast<Index>(states_.size())) {
    return false;
  }
  out = states_[static_cast<std::size_t>(current_)];
  ++current_;
  return true;
}

void UserStates::Reset() { current_ = 0; }

std::string UserStates::Label() const { return "user"; }

}  // namespace lsquant
