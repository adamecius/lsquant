// Phase R -- synthetic reconstruction oracle (operator-free).
//
// Isolates the KPM *reconstruction* knobs (kernel, moment count M, edge cutoff alpha) from
// all of the operator / Chebyshev-recursion machinery, against an analytically known answer.
//
// A single delta at E0 in (-1,1) has EXACT Chebyshev moments mu_m = T_m(E0) = cos(m*arccos E0).
// We synthesize {mu_m} directly (no Hamiltonian, no recursion), store them the way production
// does -- folding the (2 - delta_{m0}) factor in, exactly like the solver's
// `scal = 2/Nr; if(m==0) scal*=0.5` -- apply the Jackson kernel (same closed form as
// chebyshev::Moments1D::JacksonKernel), and reconstruct with the PRODUCTION kernel
// `delta_chebF` from include/chebyshev_coefficients.hpp:
//
//     rho(x) = sum_m  delta_chebF(x,m) * mu_m ,   delta_chebF(x,m) = T_m(x)/(pi*sqrt(1-x^2)).
//
// Verifications (the "fit it and verify"):
//   * Position : peak at E0 (to grid resolution).
//   * Width    : Jackson FWHM ~ 2.355 * (pi/M) * sqrt(1-E0^2)  (sigma ~ (pi/M)sqrt(1-E0^2));
//                halves when M doubles.
//   * Weight   : integral rho dx = 1 over the interior window [-alpha, alpha].
//   * alpha    : E0 near +-1 -> with alpha=1 the grid touches 1/sqrt(1-x^2) and is NON-finite;
//                with alpha=0.95 it stays finite. This is the Kubo-Bastin edge divergence that
//                motivated alpha, reproduced in a controlled, operator-free setting.
//
// This test GUARDS the production reconstruction kernel (delta_chebF). When Phase 5 extracts a
// unified reconstruction backend, re-point this test at it (it should stay green).
//
// TODO (Phase R extension, per refactor v2 feedback): add a 2-D synthetic density
// (mu_{m,n} = T_m(E0) T_n(E0)) to pin the Kubo-Bastin/Greenwood integrand
// (greenR_chebF / DgreenR_chebF + the cumulative integral), which this 1-D oracle does not cover.
//
// Usage: recon_kernel_oracle_test [TOL]      (TOL default 0.10, relative, for the width law)
#include <iostream>
#include <vector>
#include <cmath>
#include <string>
#include <complex>
#include "observable.hpp"   // lsquant::reconstruct_density_grid (the PRODUCTION 1-D reconstruction)

// Jackson damping factor g_m for an M-moment expansion -- identical closed form to
// chebyshev::Moments1D::JacksonKernel (kept local so the oracle has no object/IO setup).
static double jackson(int m, int M) {
    const double p = M_PI / (M + 1.0);
    return ((M - m + 1) * std::cos(p * m) + std::sin(p * m) / std::tan(p)) * p / M_PI;
}

// rho(x) on a uniform grid over [-alpha, alpha] for a synthetic delta at E0, M moments, Jackson.
// Phase 5: this drives the PRODUCTION primitive lsquant::reconstruct_density_grid (the exact
// routine the DOS / spectral drivers use), so the oracle now guards production -- not a copy.
// We fold the (2 - delta_{m0}) weighting and the Jackson kernel into mu, exactly as the stored,
// kernel-damped moments carry them.
static void reconstruct(double E0, int M, double alpha, int npts,
                        std::vector<double>& x, std::vector<double>& rho) {
    std::vector<double> mu(M);
    for (int m = 0; m < M; ++m)
        mu[m] = (2.0 - (m == 0 ? 1.0 : 0.0)) * std::cos(m * std::acos(E0)) * jackson(m, M);
    const std::vector<std::pair<double,double> > grid =
        lsquant::reconstruct_density_grid(mu, alpha, npts);
    x.assign(npts, 0.0); rho.assign(npts, 0.0);
    for (int i = 0; i < npts; ++i) { x[i] = grid[i].first; rho[i] = grid[i].second; }
}

static double fwhm(const std::vector<double>& x, const std::vector<double>& rho) {
    double mx = 0.0; for (double r : rho) mx = std::max(mx, r);
    const double half = mx / 2.0;
    int lo = -1, hi = -1;
    for (size_t i = 0; i < rho.size(); ++i) if (rho[i] >= half) { if (lo < 0) lo = i; hi = i; }
    return (lo < 0) ? 0.0 : x[hi] - x[lo];
}
static double peak(const std::vector<double>& x, const std::vector<double>& rho) {
    size_t k = 0; for (size_t i = 1; i < rho.size(); ++i) if (rho[i] > rho[k]) k = i;
    return x[k];
}
static double integral(const std::vector<double>& x, const std::vector<double>& rho) {
    double s = 0.0; for (size_t i = 1; i < x.size(); ++i) s += 0.5*(rho[i]+rho[i-1])*(x[i]-x[i-1]);
    return s;
}

int main(int argc, char** argv) {
    const double TOL = (argc > 1) ? std::stod(argv[1]) : 0.10;  // relative tol on the width law
    bool ok = true;

    // --- position: peak at E0 (interior), to grid resolution ---
    std::cout << "[position] peak vs E0 (M=512, alpha=0.999)\n";
    for (double E0 : {-0.5, 0.0, 0.5}) {
        std::vector<double> x, rho; reconstruct(E0, 512, 0.999, 8001, x, rho);
        const double p = peak(x, rho), dx = x[1] - x[0];
        std::cout << "   E0=" << E0 << "  peak=" << p << "  |peak-E0|=" << std::fabs(p-E0)
                  << " (grid dx=" << dx << ")\n";
        if (std::fabs(p - E0) > 5*dx) ok = false;
    }

    // --- width: FWHM ~ 2.355*(pi/M)*sqrt(1-E0^2), and halves when M doubles ---
    std::cout << "[width] Jackson FWHM scaling (E0=0)\n";
    double prev = 0.0;
    for (int M : {128, 256, 512}) {
        std::vector<double> x, rho; reconstruct(0.0, M, 0.999, 16001, x, rho);
        const double w = fwhm(x, rho), pred = 2.355 * M_PI / M;   // sqrt(1-0)=1
        const double rel = std::fabs(w - pred) / pred;
        std::cout << "   M=" << M << "  FWHM=" << w << "  pred~2.355*pi/M=" << pred
                  << "  rel=" << rel*100 << "%";
        if (prev > 0) { double ratio = prev / w; std::cout << "  ratio_vs_prevM=" << ratio;
                        if (std::fabs(ratio - 2.0) > 0.15) ok = false; }
        std::cout << "\n";
        if (rel > TOL) ok = false;
        prev = w;
    }

    // --- weight: integral = 1 over the interior window ---
    std::cout << "[weight] normalization (integral rho dx over [-alpha,alpha])\n";
    for (double E0 : {0.0, 0.5}) {
        std::vector<double> x, rho; reconstruct(E0, 256, 0.99, 16001, x, rho);
        const double I = integral(x, rho);
        std::cout << "   E0=" << E0 << "  integral=" << I << "\n";
        if (std::fabs(I - 1.0) > 0.02) ok = false;
    }

    // --- alpha edge role: E0 near the band edge ---
    std::cout << "[alpha] edge singularity at E0=0.97\n";
    {
        std::vector<double> x1, r1, x2, r2;
        reconstruct(0.97, 256, 1.00, 4001, x1, r1);   // grid reaches x=+-1
        reconstruct(0.97, 256, 0.95, 4001, x2, r2);   // grid stays inside
        bool nan_at_1 = false; for (double r : r1) if (!std::isfinite(r)) nan_at_1 = true;
        bool nan_at_a = false; for (double r : r2) if (!std::isfinite(r)) nan_at_a = true;
        std::cout << "   alpha=1.00 -> non-finite present: " << (nan_at_1 ? "yes" : "no") << "\n";
        std::cout << "   alpha=0.95 -> non-finite present: " << (nan_at_a ? "yes" : "no") << "\n";
        if (!nan_at_1) ok = false;   // alpha=1 MUST hit the 1/sqrt(1-x^2) edge
        if (nan_at_a)  ok = false;   // alpha=0.95 MUST stay finite (the safeguard)
    }

    if (ok) { std::cout << "PASS: synthetic delta oracle -- position, Jackson width law, "
                           "unit weight, and alpha edge behavior all hold\n"; return 0; }
    std::cerr << "FAIL: reconstruction kernel oracle outside tolerance\n";
    return 1;
}
