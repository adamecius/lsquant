// Phase 3 -- manual tight-binding traceability via descriptor basis tags.
//
// The magnetic-chain spin operators SX/SY/SZ are no longer opaque CSR blobs: they are
// DERIVABLE from the basis tags in chain1d_mag.HAM.desc (orbitals_per_cell / orbital_spin /
// num_cells). This test builds each operator from those tags (build_spin_operator) and asserts
// it reproduces the committed CSR EXACTLY -- the conformance gate for the manual-TB path.
// (spin_algebra, which reads the committed CSRs, remains the algebraic check; together they
// show the operators are both correct AND reconstructible from their description.)
//
// Usage: tb_descriptor_test <chain1d_mag.HAM.desc> <SX.CSR> <SY.CSR> <SZ.CSR>
#include <iostream>
#include <fstream>
#include <vector>
#include <complex>
#include <map>
#include <string>
#include "operator_descriptor.hpp"

typedef std::complex<double> cd;
typedef std::map<std::pair<int,int>, cd> SpMap;

// Parse the project CSR format: "dim nnz" / nnz (re im) / nnz columns / (dim+1) rowIndex.
static bool read_csr(const std::string& path, int& dim, SpMap& M) {
    std::ifstream f(path.c_str());
    if (!f) return false;
    int nnz; f >> dim >> nnz;
    std::vector<cd> vals(nnz); std::vector<int> cols(nnz), rowp(dim+1);
    for (int i=0;i<nnz;i++){ double re,im; f>>re>>im; vals[i]=cd(re,im); }
    for (int i=0;i<nnz;i++) f>>cols[i];
    for (int i=0;i<dim+1;i++) f>>rowp[i];
    M.clear();
    for (int r=0;r<dim;r++) for (int k=rowp[r]; k<rowp[r+1]; ++k)
        if (vals[k] != cd(0,0)) M[{r,cols[k]}] = vals[k];
    return true;
}

static int g_fail = 0;
static void compare(const std::string& tag, const SpMap& built, const SpMap& golden) {
    // every nonzero must match (both directions), to machine precision
    double maxdiff = 0.0; int extra = 0, missing = 0;
    for (const auto& kv : built) {
        auto it = golden.find(kv.first);
        if (it == golden.end()) { extra++; continue; }
        maxdiff = std::max(maxdiff, std::abs(kv.second - it->second));
    }
    for (const auto& kv : golden) if (built.find(kv.first) == built.end()) missing++;
    bool ok = (extra==0 && missing==0 && maxdiff < 1e-12 && built.size()==golden.size());
    std::cout << (ok?"  OK   ":"  FAIL ") << tag
              << ": built nnz=" << built.size() << " golden nnz=" << golden.size()
              << " maxdiff=" << maxdiff << " extra=" << extra << " missing=" << missing << "\n";
    if (!ok) g_fail = 1;
}

int main(int argc, char** argv) {
    if (argc < 5) { std::cerr << "usage: " << argv[0] << " <HAM.desc> <SX.CSR> <SY.CSR> <SZ.CSR>\n"; return 2; }
    lsquant::OperatorDescriptor d = lsquant::read_descriptor(argv[1]);
    if (!d.valid || !d.has_basis) { std::cerr << "FAIL: descriptor has no basis block\n"; return 1; }
    std::cout << "basis: orbitals_per_cell=" << d.orbitals_per_cell
              << " num_cells=" << d.num_cells << " spins=";
    for (int s : d.orbital_spin) std::cout << (s>0?'+':'-');
    std::cout << "\n";

    const char comps[3] = {'X','Y','Z'};
    for (int i=0;i<3;i++) {
        std::vector<int> r,c; std::vector<cd> v;
        if (!lsquant::build_spin_operator(d, comps[i], r, c, v)) {
            std::cerr << "FAIL: build_spin_operator(" << comps[i] << ") returned false\n"; return 1;
        }
        SpMap built; for (size_t k=0;k<r.size();++k) if (v[k]!=cd(0,0)) built[{r[k],c[k]}] += v[k];
        int dim; SpMap golden;
        if (!read_csr(argv[2+i], dim, golden)) { std::cerr << "FAIL: cannot read " << argv[2+i] << "\n"; return 1; }
        if (dim != d.orbitals_per_cell*d.num_cells)
            { std::cout << "  FAIL S" << comps[i] << ": dim " << dim << " != basis dim\n"; g_fail=1; }
        compare(std::string("S")+comps[i], built, golden);
    }

    if (g_fail==0) { std::cout << "PASS: SX/SY/SZ reproduced exactly from descriptor basis tags\n"; return 0; }
    std::cerr << "FAIL: built spin operators differ from committed CSRs\n"; return 1;
}
