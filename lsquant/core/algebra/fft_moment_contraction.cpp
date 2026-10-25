#include "core/algebra/fft_moment_contraction.hpp"
#include "core/algebra/fftw_planner.hpp"
#include <fftw3.h>
#include <cmath>
#include <limits>
#include <memory>
#include <stdexcept>

namespace lsquant {
std::vector<Scalar> FftMomentContraction(const Scalar* moments, Index m0, Index m1,
                                        Index points) {
  if (!moments || m0 <= 0 || m1 <= 0 || points <= 0 ||
      points > std::numeric_limits<int>::max()) {
    throw std::invalid_argument("FftMomentContraction: invalid dimensions/buffer");
  }
  const double pi = std::acos(-1.0), offset = pi / (2.0 * points);
  const Index limit = m0 + m1;
  std::vector<Scalar> phases(static_cast<std::size_t>(2 * limit + 1));
  for (Index k = -limit; k <= limit; ++k) phases[k + limit] = std::polar(1.0, k * offset);
  std::vector<Scalar> result(static_cast<std::size_t>(2 * points));
  std::unique_ptr<fftw_complex, decltype(&fftw_free)> in(fftw_alloc_complex(points), &fftw_free);
  std::unique_ptr<fftw_complex, decltype(&fftw_free)> out(fftw_alloc_complex(points), &fftw_free);
  if (!in || !out) throw std::bad_alloc();
  fftw_plan plan;
  {
    std::lock_guard<std::mutex> lock(detail::FftwPlannerMutex());
    plan = fftw_plan_dft_1d(static_cast<int>(points), in.get(), out.get(), FFTW_BACKWARD, FFTW_ESTIMATE);
  }
  if (!plan) throw std::runtime_error("FftMomentContraction: FFTW planning failed");
  // No allocation or exception after planning: destroy before buffer release.
  // cos(m*theta)*exp(i*sign*n*theta) is half the sum of frequencies
  // m+sign*n and -m+sign*n. Fold those frequencies modulo N AFTER applying
  // their offset phases; the positive-sign transform covers folded angles
  // where sin(raw_angle)<0 without conjugating the complex moment tensor.
  for (int weighted = 0; weighted < 2; ++weighted) {
    for (int sign : {-1, 1}) {
      for (Index j = 0; j < points; ++j) { in.get()[j][0] = 0.0; in.get()[j][1] = 0.0; }
      for (Index m = 0; m < m0; ++m) {
        for (Index n = 0; n < m1; ++n) {
          const Scalar value = (weighted ? static_cast<double>(n) : 1.0) * 0.5 * moments[m * m1 + n];
          for (Index k : {m + sign * n, -m + sign * n}) {
            const Index j = (k % points + points) % points;
            const Scalar term = value * phases[k + limit];
            in.get()[j][0] += term.real();
            in.get()[j][1] += term.imag();
          }
        }
      }
      fftw_execute(plan);
      for (Index j = 0; j < points; ++j) {
        const double angle = 2.0 * pi * j / points + offset;
        const int selected = std::sin(angle) >= 0.0 ? -1 : 1;
        if (sign == selected) result[weighted * points + j] = Scalar(out.get()[j][0], out.get()[j][1]);
      }
    }
  }
  {
    std::lock_guard<std::mutex> lock(detail::FftwPlannerMutex());
    fftw_destroy_plan(plan);
  }
  return result;
}
}  // namespace lsquant
