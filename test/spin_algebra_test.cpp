// Spin-operator algebra golden (#6 prerequisite): the magnetic-chain spin density
// operators S_x, S_y, S_z (built by wannier2sparse from the spin doubling) must form a
// correct spin-1/2 representation, and S_x must NOT commute with H -- the latter is the
// physical source of transverse-spin precession that test #6 will measure. Asserting
// [H,S_x]=0 here would be the physics error (a collinear-but-static reading), so we
// assert it is non-zero.
//
// Usage: spin_algebra_test <HAM.CSR> <SX.CSR> <SY.CSR> <SZ.CSR> <TOL>
//
// CSR text layout (row-major): line1 "dim nnz"; line2 = 2*nnz interleaved "re im"
// values; line3 = nnz column indices; line4 = dim+1 row pointers.
//
// Checks (gamma=hbar=1, S in units of hbar/2 so eigenvalues are +/-1):
//   - each S_a Hermitian and traceless
//   - S_a^2 = I            (involutory: a genuine Pauli operator, eig +/-1)
//   - [S_x,S_y] = 2 i S_z  (and the two cyclic partners)  -> su(2)
//   - S_x != 0  and  [H,S_x] != 0   (precession is reachable; would-be #6 source)
#include <iostream>
#include <fstream>
#include <vector>
#include <complex>
#include <cmath>
#include <string>

using cd = std::complex<double>;
using Mat = std::vector<std::vector<cd>>;

static Mat load_csr(const std::string& path) {
    std::ifstream f(path);
    if (!f) { std::cerr << "cannot open " << path << "\n"; std::exit(2); }
    long dim, nnz; f >> dim >> nnz;
    std::vector<double> re(nnz), im(nnz);
    for (long k = 0; k < nnz; ++k) f >> re[k] >> im[k];
    std::vector<long> col(nnz); for (long k = 0; k < nnz; ++k) f >> col[k];
    std::vector<long> rptr(dim + 1); for (long r = 0; r <= dim; ++r) f >> rptr[r];
    if (!f) { std::cerr << "truncated CSR " << path << "\n"; std::exit(2); }
    Mat M(dim, std::vector<cd>(dim, cd(0, 0)));
    for (long r = 0; r < dim; ++r)
        for (long k = rptr[r]; k < rptr[r + 1]; ++k) M[r][col[k]] = cd(re[k], im[k]);
    return M;
}

static int N(const Mat& A) { return (int)A.size(); }
static Mat mul(const Mat& A, const Mat& B) {
    int n = N(A); Mat C(n, std::vector<cd>(n, cd(0, 0)));
    for (int i = 0; i < n; ++i)
        for (int k = 0; k < n; ++k) { cd a = A[i][k]; if (a == cd(0,0)) continue;
            for (int j = 0; j < n; ++j) C[i][j] += a * B[k][j]; }
    return C;
}
// max_{ij} | A - (s*B + t*C) |   (C optional via t=0)
static double maxdiff(const Mat& A, const Mat& B, cd s = 1, const Mat* C = nullptr, cd t = 0) {
    int n = N(A); double m = 0;
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j) {
            cd v = A[i][j] - s * B[i][j];
            if (C) v -= t * (*C)[i][j];
            m = std::max(m, std::abs(v));
        }
    return m;
}
static double maxabs(const Mat& A) { double m = 0; for (auto& r : A) for (auto& v : r) m = std::max(m, std::abs(v)); return m; }
static Mat adj(const Mat& A) { int n = N(A); Mat T(n, std::vector<cd>(n)); for (int i=0;i<n;++i) for (int j=0;j<n;++j) T[i][j]=std::conj(A[j][i]); return T; }
static cd trace(const Mat& A) { cd t=0; for (int i=0;i<N(A);++i) t+=A[i][i]; return t; }
static Mat identity(int n) { Mat I(n, std::vector<cd>(n, cd(0,0))); for (int i=0;i<n;++i) I[i][i]=cd(1,0); return I; }
static Mat comm(const Mat& A, const Mat& B) { Mat ab=mul(A,B), ba=mul(B,A); int n=N(A); Mat C(n,std::vector<cd>(n)); for(int i=0;i<n;++i) for(int j=0;j<n;++j) C[i][j]=ab[i][j]-ba[i][j]; return C; }

int main(int argc, char** argv) {
    if (argc < 6) { std::cerr << "usage: " << argv[0] << " <HAM> <SX> <SY> <SZ> <TOL>\n"; return 2; }
    const double TOL = std::stod(argv[5]);
    Mat H = load_csr(argv[1]), Sx = load_csr(argv[2]), Sy = load_csr(argv[3]), Sz = load_csr(argv[4]);
    int n = N(Sx);
    if (N(H) != n || N(Sy) != n || N(Sz) != n) { std::cerr << "dimension mismatch\n"; return 2; }
    Mat I = identity(n);
    const cd i_(0, 1);
    int rc = 0;
    auto chk = [&](const std::string& name, double r) {
        std::cout << "  " << name << " = " << r << (r < TOL ? "   ok" : "   FAIL") << "\n";
        if (!(r < TOL)) rc = 1;
    };
    std::cout << "spin algebra (dim " << n << ", TOL " << TOL << "):\n";
    chk("Hermiticity Sx", maxdiff(Sx, adj(Sx)));
    chk("Hermiticity Sy", maxdiff(Sy, adj(Sy)));
    chk("Hermiticity Sz", maxdiff(Sz, adj(Sz)));
    chk("|tr Sx|", std::abs(trace(Sx))); chk("|tr Sy|", std::abs(trace(Sy))); chk("|tr Sz|", std::abs(trace(Sz)));
    chk("|Sx^2 - I|", maxdiff(mul(Sx, Sx), I)); chk("|Sy^2 - I|", maxdiff(mul(Sy, Sy), I)); chk("|Sz^2 - I|", maxdiff(mul(Sz, Sz), I));
    chk("|[Sx,Sy] - 2i Sz|", maxdiff(comm(Sx, Sy), Sz, 2.0 * i_));
    chk("|[Sy,Sz] - 2i Sx|", maxdiff(comm(Sy, Sz), Sx, 2.0 * i_));
    chk("|[Sz,Sx] - 2i Sy|", maxdiff(comm(Sz, Sx), Sy, 2.0 * i_));

    // Physical reachability of precession (#6): S_x nonzero AND [H,S_x] nonzero.
    double sx_norm = maxabs(Sx), hsx = maxabs(comm(H, Sx));
    std::cout << "  |Sx| = " << sx_norm << "  |[H,Sx]| = " << hsx
              << "   (precession source: both must be > " << TOL << ")\n";
    if (!(sx_norm > TOL)) { std::cout << "  FAIL: Sx is (near) zero -- transverse spin unreachable\n"; rc = 1; }
    if (!(hsx     > TOL)) { std::cout << "  FAIL: [H,Sx]=0 -- no precession (collinear-static)\n"; rc = 1; }

    std::cout << (rc ? "FAIL: spin operators are not a correct su(2) representation\n"
                     : "PASS: spin operators form su(2) and [H,Sx]!=0 (precession reachable)\n");
    return rc;
}
