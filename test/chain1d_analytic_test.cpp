// Analytic golden test: LinQT Kubo-Greenwood Chebyshev moment matrix for the
// monoatomic 1D chain vs the closed form in Angel's thesis Eq. (4.9).
//
// Usage: chain1d_analytic_test <chebmom2D file> <analytic M ref> <TOL>
//
// .chebmom2D layout: line1 "dim bandwidth bandcenter"; line2 "M N";
//                    then M*N "(real imag)" pairs, row-major.
// Normalization ("unit-coefficient", thesis Fig. 4.2): the stored moments already
// carry the Algorithm-1 g_m*g_n folding (CorrelationExpansionMoments applies
// 4/R * (1/2 if m==0) * (1/2 if n==0)); the only remaining map to M[m,n] is one
// overall scale = the t^2 a^2/hbar^2 prefactor, which equals the m>=2 diagonal
// plateau (M[m,m]=1 there). We divide by that plateau, then compare to Eq. (4.9).
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>
#include <string>

int main(int argc, char** argv) {
    if (argc < 4) { std::cerr << "usage: " << argv[0] << " <chebmom2D> <Mref> <TOL>\n"; return 2; }
    const std::string momPath = argv[1], refPath = argv[2];
    const double TOL = std::stod(argv[3]);
    const int B = 10; // leading block compared

    // --- read moments ---
    std::ifstream f(momPath);
    if (!f) { std::cerr << "cannot open moments: " << momPath << "\n"; return 2; }
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

    // --- unit-coefficient normalization: divide by the m>=2 diagonal plateau ---
    double plateau = 0.0; int cnt = 0;
    for (int m = 2; m < B; ++m) { plateau += mu[m][m]; ++cnt; }
    plateau /= cnt;
    if (!(plateau > 0.0)) { std::cerr << "non-positive diagonal plateau: " << plateau << "\n"; return 1; }
    std::cout << "diagonal plateau (m=2.." << B-1 << ") = " << plateau << "\n";

    // --- read analytic reference (skip # comment lines) ---
    std::ifstream r(refPath);
    if (!r) { std::cerr << "cannot open reference: " << refPath << "\n"; return 2; }
    std::vector<std::vector<double>> Mref;
    std::string line;
    while (std::getline(r, line)) {
        size_t p = line.find_first_not_of(" \t");
        if (p == std::string::npos || line[p] == '#') continue;
        std::istringstream ss(line); std::vector<double> row; double v;
        while (ss >> v) row.push_back(v);
        if (!row.empty()) Mref.push_back(row);
    }
    if ((int)Mref.size() < B) { std::cerr << "reference has < " << B << " rows\n"; return 2; }

    // --- compare ---
    double maxRes = 0.0, maxOff = 0.0; int im_=0, jm_=0;
    for (int m = 0; m < B; ++m)
        for (int n = 0; n < B; ++n) {
            double comp = mu[m][n] / plateau;
            double res = std::fabs(comp - Mref[m][n]);
            if (res > maxRes) { maxRes = res; im_ = m; jm_ = n; }
            if (Mref[m][n] == 0.0) maxOff = std::max(maxOff, std::fabs(comp));
        }
    std::cout << "diag[0,0]=" << mu[0][0]/plateau << " diag[1,1]=" << mu[1][1]/plateau
              << "  band[0,2]=" << mu[0][2]/plateau << " band[2,4]=" << mu[2][4]/plateau << "\n";
    std::cout << "max|computed - M| = " << maxRes << " at (" << im_ << "," << jm_ << ")"
              << "   max|off-structure| = " << maxOff << "   TOL = " << TOL << "\n";
    if (maxRes < TOL) { std::cout << "PASS: chain1d Kubo-Greenwood moments match thesis Eq. (4.9)\n"; return 0; }
    std::cerr << "FAIL: residual " << maxRes << " exceeds TOL " << TOL << "\n";
    return 1;
}
