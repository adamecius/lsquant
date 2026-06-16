// Jackson-kernel application golden-gate (plan Step 0.5c).
//
// Freezes the output of Moments1D::ApplyJacksonKernel and Moments2D::ApplyJacksonKernel on a
// fixed input moment table (all moments = 1, no kernel reduction at this broadening) against
// INDEPENDENT literal golden values (computed once from the damping formula, embedded here so
// the test does not re-derive the in-source formula). This guards B1 (single-source jackson_g)
// and B2 (base retype) against any value drift in the damping coefficient or its application.
//
// 1D:  g_D(m,M) = ((M - m + 1) cos(phi m) + sin(phi m)/tan(phi)) phi/pi,  phi = pi/(M+1)
// 2D:  tensor product g_D(m0,M0) * g_D(m1,M1)
//
// Usage: jackson_kernel_test [TOL]   (default TOL = 1e-12)
#include <iostream>
#include <cmath>
#include <cstdlib>
#include "chebyshev_moments.hpp"

int main(int argc, char *argv[])
{
    const double TOL = (argc > 1) ? std::atof(argv[1]) : 1e-12;
    int failures = 0;

    // --- 1D: M = 8, BandWidth = 1, broadening = 0.5 (eta ~ 1e-3 -> no moment reduction) ----
    // golden: input mu_m = 1 scaled by g_D(m,8)
    const double g1d[8] = {
        1.0,
        0.93969262078590843,
        0.79203950499464715,
        0.59770947128575747,
        0.39710866135994871,
        0.22346048369301835,
        0.097709471285757493,
        0.025995061875669197,
    };
    {
        chebyshev::Moments1D mu(8);
        mu.BandWidth(1.0);
        for (size_t m = 0; m < 8; m++) mu(m) = 1.0;
        mu.ApplyJacksonKernel(0.5);

        if (mu.HighestMomentNumber() != 8) {
            std::cerr << "FAIL: 1D moment count changed (reduction unexpected): "
                      << mu.HighestMomentNumber() << std::endl;
            failures++;
        }
        for (size_t m = 0; m < 8; m++) {
            const double err = std::abs(mu(m).real() - g1d[m]) + std::abs(mu(m).imag());
            if (err > TOL) {
                std::cerr << "FAIL: 1D Jackson m=" << m << " got " << mu(m).real()
                          << " expected " << g1d[m] << " (err " << err << ")" << std::endl;
                failures++;
            }
        }
    }

    // --- 2D: M0 = 4, M1 = 5, BandWidth = 1, broadening = 0.5/0.5 (no reduction) -------------
    const double g0[4] = { 1.0, 0.80901699437494745, 0.44721359549995793, 0.13819660112501056 };
    const double g1[5] = { 1.0, 0.8660254037844386, 0.58333333333333337,
                           0.28867513459481292, 0.08333333333333344 };
    {
        chebyshev::Moments2D mu(4, 5);
        mu.BandWidth(1.0);
        for (int m0 = 0; m0 < 4; m0++)
            for (int m1 = 0; m1 < 5; m1++) mu(m0, m1) = 1.0;
        mu.ApplyJacksonKernel(0.5, 0.5);

        if (mu.HighestMomentNumber(0) != 4 || mu.HighestMomentNumber(1) != 5) {
            std::cerr << "FAIL: 2D moment count changed (reduction unexpected): "
                      << mu.HighestMomentNumber(0) << "x" << mu.HighestMomentNumber(1) << std::endl;
            failures++;
        }
        for (int m0 = 0; m0 < 4; m0++)
            for (int m1 = 0; m1 < 5; m1++) {
                const double expected = g0[m0] * g1[m1];
                const double err = std::abs(mu(m0, m1).real() - expected) + std::abs(mu(m0, m1).imag());
                if (err > TOL) {
                    std::cerr << "FAIL: 2D Jackson (" << m0 << "," << m1 << ") got "
                              << mu(m0, m1).real() << " expected " << expected
                              << " (err " << err << ")" << std::endl;
                    failures++;
                }
            }
    }

    if (failures == 0) {
        std::cout << "PASS: Jackson kernel (1D + 2D) reproduces frozen golden coefficients" << std::endl;
        return 0;
    }
    std::cerr << failures << " Jackson-kernel golden check(s) FAILED" << std::endl;
    return 1;
}
