// Generic numerical comparator for two .chebmom moment files (1D or 2D).
//
// Usage: chebmom_compare <fileA> <fileB> <TOL>
//
// Asserts max_k |A_k - B_k| < TOL over the complex moments. Unlike a byte diff this is
// robust to the floating-point reduction order that differs across compiler / Eigen /
// OpenMP versions (~1e-15), while still catching any real physics regression (which a
// src/ refactor bug would produce at O(1e-2) or larger). Use a tolerance well above
// toolchain noise but far below the stochastic-trace error, e.g. 1e-6.
//
// Layout: line1 "dim bandwidth bandcenter"; for 2D line2 "M N" then M*N pairs, for 1D
// line2 "numMoms" then numMoms pairs; each pair is "re im". The header shapes are read
// generically: everything after the first two header lines is parsed as re/im pairs.
#include <iostream>
#include <fstream>
#include <vector>
#include <complex>
#include <cmath>
#include <string>
#include <sstream>

static std::vector<std::complex<double>> load(const std::string& path, long& count) {
    std::ifstream f(path);
    if (!f) { std::cerr << "cannot open " << path << "\n"; std::exit(2); }
    std::string l1, l2; std::getline(f, l1); std::getline(f, l2);   // two header lines
    std::vector<std::complex<double>> v;
    std::string line;
    while (std::getline(f, line)) {
        std::istringstream ss(line);
        double re, im;
        if (ss >> re >> im) v.emplace_back(re, im);
    }
    count = (long)v.size();
    return v;
}

int main(int argc, char** argv) {
    if (argc < 4) { std::cerr << "usage: " << argv[0] << " <A> <B> <TOL>\n"; return 2; }
    const double TOL = std::stod(argv[3]);
    long na, nb;
    auto A = load(argv[1], na);
    auto B = load(argv[2], nb);
    if (na == 0 || na != nb) { std::cerr << "FAIL: moment count mismatch (" << na << " vs " << nb << ")\n"; return 1; }

    double maxd = 0.0; long argmax = 0; bool nan = false;
    for (long k = 0; k < na; ++k) {
        if (!std::isfinite(A[k].real()) || !std::isfinite(A[k].imag())) nan = true;
        double d = std::abs(A[k] - B[k]);
        if (d > maxd) { maxd = d; argmax = k; }
    }
    if (nan) { std::cerr << "FAIL: non-finite moment in " << argv[1] << "\n"; return 1; }
    std::cout << "compared " << na << " moments  max|A-B| = " << maxd
              << " at k=" << argmax << "   TOL = " << TOL << "\n";
    if (maxd < TOL) { std::cout << "PASS: moment files agree within tolerance\n"; return 0; }
    std::cerr << "FAIL: max|A-B| " << maxd << " exceeds TOL " << TOL << "\n";
    return 1;
}
