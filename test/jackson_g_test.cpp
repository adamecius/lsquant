// Isolated scalar unit test for the single-sourced Jackson damping coefficient (plan Step B1).
//
// chebyshev::Moments::JacksonKernel(m, M) is now the ONE place the g_D(m,M) formula lives;
// every ApplyJacksonKernel (1D/2D/TD/Local) and the solver's inline broadening call it. Being
// a pure static formula, it is tested here directly against independent analytic references --
// no Moments object, no operators, no recursion.
//
//   g_D(m, M) = ((M - m + 1) cos(phi m) + sin(phi m)/tan(phi)) phi/pi,   phi = pi/(M+1)
//
// Usage: jackson_g_test [TOL]   (default TOL = 1e-13)
#include <iostream>
#include <cmath>
#include <cstdlib>
#include "chebyshev_moments.hpp"

int main(int argc, char *argv[])
{
    const double TOL = (argc > 1) ? std::atof(argv[1]) : 1e-13;
    int failures = 0;

    // Reference values computed independently from the closed form (not from this build).
    struct Ref { double m, M, g; };
    const Ref refs[] = {
        { 0,   8,   1.0                    },
        { 3,   8,   0.59770947128575747    },
        { 7,   8,   0.025995061875669197   },
        { 0,   64,  1.0                    },
        { 32,  64,  0.33023686812562286    },
        { 63,  64,  7.1821004343993646e-05 },
        { 0,   128, 1.0                    },
        { 1,   128, 0.99970346984513925    },
        { 127, 128, 9.1933702262909864e-06 },
    };
    for (const Ref& r : refs) {
        const double got = chebyshev::Moments::JacksonKernel(r.m, r.M);
        const double err = std::abs(got - r.g);
        if (err > TOL) {
            std::cerr << "FAIL: g_D(" << r.m << "," << r.M << ") got " << got
                      << " expected " << r.g << " (err " << err << ")" << std::endl;
            failures++;
        }
    }

    // Analytic invariants: g_D(0,M) == 1 EXACTLY, and 0 < g_D(m,M) <= 1 (positive, normalised).
    for (double M : {4.0, 16.0, 100.0, 1000.0}) {
        const double g0 = chebyshev::Moments::JacksonKernel(0.0, M);
        if (g0 != 1.0) { std::cerr << "FAIL: g_D(0," << M << ") != 1 exactly: " << g0 << std::endl; failures++; }
        for (double m = 0.0; m < M; m += 1.0) {
            const double g = chebyshev::Moments::JacksonKernel(m, M);
            if (!(g > 0.0) || g > 1.0 + TOL) {
                std::cerr << "FAIL: g_D(" << m << "," << M << ")=" << g << " out of (0,1]" << std::endl;
                failures++;
            }
        }
    }

    if (failures == 0) {
        std::cout << "PASS: Jackson coefficient matches analytic references and invariants" << std::endl;
        return 0;
    }
    std::cerr << failures << " Jackson-coefficient check(s) FAILED" << std::endl;
    return 1;
}
