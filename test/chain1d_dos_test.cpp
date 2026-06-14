// DOS-moments golden (#1): 1D-chain DOS Chebyshev moments vs the closed form.
//
// Usage: chain1d_dos_test <chebmom1D file> <TOL>
//
// The clean 1D chain's rescaled DOS is exactly the Chebyshev weight 1/(pi*sqrt(1-x^2))
// (requires bounds = the true band [-2,2]), so Cbar^DOS_m = (1/N)Tr[T_m(H~)] = delta_{m,0}
// -- proven in lsquant_derive_core.py Part A, exposed by lsquant_chain_reference.dos_moment.
// Reference embedded here (no Python / no data file). Any weight leaking to m>=1 beyond
// the stochastic band flags a trace/rescale bug.
//
// .chebmom1D layout: line1 "dim bandwidth bandcenter"; line2 "numMoms"; then numMoms
// "(real imag)" pairs. We normalize the computed moments by C_0 (=> Cbar_0=1) and
// assert |Cbar_m - delta_{m,0}| < TOL.
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <string>

int main(int argc, char** argv) {
    if (argc < 3) { std::cerr << "usage: " << argv[0] << " <chebmom1D> <TOL>\n"; return 2; }
    const double TOL = std::stod(argv[2]);

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

    double maxRes = 0.0; int mmax = 0;
    for (int m = 0; m < numMoms; ++m) {
        double cm = C[m] / C[0];                 // normalized: Cbar_0 -> 1
        double ref = (m == 0) ? 1.0 : 0.0;       // oracle dos_moment(m) = delta_{m,0}
        double res = std::fabs(cm - ref);
        if (res > maxRes) { maxRes = res; mmax = m; }
    }
    std::cout << "max|Cbar_m - delta_m0| = " << maxRes << " at m=" << mmax
              << "   (over m=0.." << numMoms-1 << ")   TOL = " << TOL << "\n";
    if (maxRes < TOL) { std::cout << "PASS: 1D-chain DOS moments match delta_{m,0}\n"; return 0; }
    std::cerr << "FAIL: residual " << maxRes << " exceeds TOL " << TOL << "\n";
    return 1;
}
