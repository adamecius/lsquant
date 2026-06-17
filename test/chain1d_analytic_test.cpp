// Analytic golden test (#3): LSQUANT Kubo-Greenwood Chebyshev moment matrix for the
// monoatomic 1D chain vs the closed form (thesis Eq. 4.9, B.15 normalized).
//
// Usage: chain1d_analytic_test <chebmom2D file> <TOL>
//
// The B.15 reference is the closed form proven symbolically in lsquant_derive_core.py
// and exposed by lsquant_chain_reference.kg_moment(); it is embedded here (no Python /
// no data file needed to build or run) -- with gamma=a=hbar=1:
//     M[m,n] = d_mn/(1+d_m0+d_m1) - 1/2 (1-d_mn) d_{|m-n|,2}
// i.e. diagonal {1/2,1/2,1,1,...}, the |m-n|=2 band = -1/2, zero elsewhere.
//
// .chebmom2D layout: line1 "dim bandwidth bandcenter"; line2 "M N"; then M*N
// "(real imag)" pairs, row-major. The stored moments already carry Algorithm-1's
// g_m*g_n folding; the only remaining map to M is one overall scale (the
// t^2 a^2/hbar^2 prefactor) = the m>=2 diagonal plateau, which we divide out.
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <string>

static double Mref(int m, int n) {            // oracle kg_moment, gamma=a=hbar=1
    int dmn = (m == n), dm0 = (m == 0), dm1 = (m == 1), d2 = (std::abs(m - n) == 2);
    return (dmn ? 1.0 / (1 + dm0 + dm1) : 0.0) - 0.5 * (1 - dmn) * d2;
}

int main(int argc, char** argv) {
    if (argc < 3) { std::cerr << "usage: " << argv[0] << " <chebmom2D> <TOL>\n"; return 2; }
    const double TOL = std::stod(argv[2]);
    const int B = 10; // leading block compared

    std::ifstream f(argv[1]);
    if (!f) { std::cerr << "cannot open moments: " << argv[1] << "\n"; return 2; }
    long dim; double bw, bc; int M, N;
    f >> dim >> bw >> bc >> M >> N;
    if (M < B || N < B) { std::cerr << "moment matrix smaller than " << B << "\n"; return 2; }
    std::vector<std::vector<double>> mu(M, std::vector<double>(N, 0.0));
    double maxIm = 0.0;
    for (int m = 0; m < M; ++m)
        for (int n = 0; n < N; ++n) {
            double re, im; f >> re >> im;
            if (!f) { std::cerr << "truncated moment file at (" << m << "," << n << ")\n"; return 2; }
            mu[m][n] = re; maxIm = std::max(maxIm, std::fabs(im));
        }
    std::cout << "moments: dim=" << dim << " bandwidth=" << bw << " bandcenter=" << bc
              << " size=" << M << "x" << N << "  max|Im|=" << maxIm << "\n";

    double plateau = 0.0; int cnt = 0;
    for (int m = 2; m < B; ++m) { plateau += mu[m][m]; ++cnt; }
    plateau /= cnt;
    if (!(plateau > 0.0)) { std::cerr << "non-positive diagonal plateau: " << plateau << "\n"; return 1; }
    std::cout << "diagonal plateau (m=2.." << B-1 << ") = " << plateau << "\n";

    double maxRes = 0.0, maxOff = 0.0; int im_ = 0, jm_ = 0;
    for (int m = 0; m < B; ++m)
        for (int n = 0; n < B; ++n) {
            double comp = mu[m][n] / plateau, ref = Mref(m, n);
            double res = std::fabs(comp - ref);
            if (res > maxRes) { maxRes = res; im_ = m; jm_ = n; }
            if (ref == 0.0) maxOff = std::max(maxOff, std::fabs(comp));
        }
    std::cout << "diag[0,0]=" << mu[0][0]/plateau << " diag[1,1]=" << mu[1][1]/plateau
              << "  band[0,2]=" << mu[0][2]/plateau << " band[2,4]=" << mu[2][4]/plateau << "\n";
    std::cout << "max|computed - M| = " << maxRes << " at (" << im_ << "," << jm_ << ")"
              << "   max|off-structure| = " << maxOff << "   TOL = " << TOL << "\n";
    if (maxRes < TOL) { std::cout << "PASS: chain1d Kubo-Greenwood moments match Eq. (4.9)\n"; return 0; }
    std::cerr << "FAIL: residual " << maxRes << " exceeds TOL " << TOL << "\n";
    return 1;
}
