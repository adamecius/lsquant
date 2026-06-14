// DOS-moments golden (#1): 1D-chain DOS Chebyshev moments vs the oracle.
//
// Usage: chain1d_dos_test <chebmom1D file> <DOS reference> <TOL>
//
// The clean 1D chain's rescaled DOS is exactly the Chebyshev weight 1/(pi*sqrt(1-x^2))
// (requires bounds = the true band [-2,2]), so Cbar^DOS_m = (1/N)Tr[T_m(H~)] = delta_{m,0}
// (lsquant_chain_reference.dos_moment, proven in lsquant_derive_core.py Part A). Any weight
// leaking to m>=1 beyond the stochastic band flags a trace/rescale bug.
//
// .chebmom1D layout: line1 "dim bandwidth bandcenter"; line2 "numMoms"; then numMoms
// "(real imag)" pairs. We normalize the computed moments by C_0 and compare to the
// reference (delta_{m,0}); the reference is the oracle's source of truth.
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>
#include <string>

int main(int argc, char** argv) {
    if (argc < 4) { std::cerr << "usage: " << argv[0] << " <chebmom1D> <DOSref> <TOL>\n"; return 2; }
    const double TOL = std::stod(argv[3]);

    std::ifstream f(argv[1]);
    if (!f) { std::cerr << "cannot open moments: " << argv[1] << "\n"; return 2; }
    long dim; double bw, bc; int numMoms;
    f >> dim >> bw >> bc >> numMoms;
    std::vector<double> C; C.reserve(numMoms);
    double maxIm = 0.0;
    for (int m = 0; m < numMoms; ++m) {
        double re, im; f >> re >> im;
        if (!f) { std::cerr << "truncated at moment " << m << "\n"; return 2; }
        C.push_back(re); maxIm = std::max(maxIm, std::fabs(im));
    }
    std::cout << "DOS moments: dim=" << dim << " bandwidth=" << bw << " bandcenter=" << bc
              << " numMoms=" << numMoms << "  max|Im|=" << maxIm << "\n";
    if (!(std::fabs(C[0]) > 0.0)) { std::cerr << "C_0 is zero\n"; return 1; }

    std::ifstream r(argv[2]);
    if (!r) { std::cerr << "cannot open reference: " << argv[2] << "\n"; return 2; }
    std::vector<double> ref; std::string line;
    while (std::getline(r, line)) {
        size_t p = line.find_first_not_of(" \t");
        if (p == std::string::npos || line[p] == '#') continue;
        std::istringstream ss(line); double v; if (ss >> v) ref.push_back(v);
    }
    const int K = std::min((int)ref.size(), numMoms);

    double maxRes = 0.0; int mmax = 0;
    for (int m = 0; m < K; ++m) {
        double cm = C[m] / C[0];                 // normalized: Cbar_0 -> 1
        double res = std::fabs(cm - ref[m]);
        if (res > maxRes) { maxRes = res; mmax = m; }
    }
    std::cout << "max|C_m/C_0 - delta_m0| = " << maxRes << " at m=" << mmax
              << "   (over m=0.." << K-1 << ")   TOL = " << TOL << "\n";
    if (maxRes < TOL) { std::cout << "PASS: 1D-chain DOS moments match delta_{m,0}\n"; return 0; }
    std::cerr << "FAIL: residual " << maxRes << " exceeds TOL " << TOL << "\n";
    return 1;
}
