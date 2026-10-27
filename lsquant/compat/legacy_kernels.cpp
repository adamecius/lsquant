/**
 * @file legacy_kernels.cpp
 * @brief Implementation of gate-only legacy retarded Green's function table evaluators.
 */
#include "compat/legacy_kernels.hpp"

#include <cmath>
#include <gsl/gsl_sf_bessel.h>
#include <stdexcept>

namespace lsquant {
namespace compat {

void ComputeLegacyGreenTables(const double* xGrid, Index numPoints, Index numMoments0,
                              Index numMoments1, Scalar* greenR, Scalar* dgreenR) {
  if (numMoments0 <= 0 || numMoments1 <= 0) {
    throw std::invalid_argument(
        "ComputeLegacyGreenTables: numMoments0 and numMoments1 must be positive");
  }
  if (numPoints < 0) {
    throw std::invalid_argument("ComputeLegacyGreenTables: numPoints must be non-negative");
  }
  if (numPoints > 0 && (xGrid == nullptr || greenR == nullptr || dgreenR == nullptr)) {
    throw std::invalid_argument("ComputeLegacyGreenTables: null buffer pointer");
  }

#pragma omp parallel for schedule(static)
  for (Index i = 0; i < numPoints; ++i) {
    const double x = xGrid[i];
    for (Index m0 = 0; m0 < numMoments0; ++m0) {
      greenR[m0 * numPoints + i] = greenR_chebF(x, static_cast<double>(m0));
    }
    for (Index m1 = 0; m1 < numMoments1; ++m1) {
      dgreenR[m1 * numPoints + i] = DgreenR_chebF(x, static_cast<double>(m1));
    }
  }
}

std::vector<Scalar> ComputeLegacyBesselCoefficients(double x, double tol) {
  if (tol <= 0.0) {
    throw std::invalid_argument("ComputeLegacyBesselCoefficients: tol must be positive");
  }
  const Scalar I(0.0, 1.0);
  double momCutOff = 2.0 * tol;
  int n = 0;
  Scalar nIp(1.0, 0.0);
  double jn = 0.5 * gsl_sf_bessel_Jn(0, x);
  std::vector<Scalar> coeffs;
  while (momCutOff > tol) {
    coeffs.push_back(nIp * Scalar(2.0, 0.0) * jn);
    nIp *= -I;
    ++n;
    jn = gsl_sf_bessel_Jn(n, x);
    momCutOff = std::fabs(jn);
  }
  return coeffs;
}

}  // namespace compat
}  // namespace lsquant
