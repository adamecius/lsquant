// spmm_batchmultiply_equiv_test -- gate for the SpMM block kernel (Phase B, W2).
//
// PROBLEM. BatchMultiply (Y = a*M*X + b*Y over R right-hand sides) was a
// std::terminate() stub. Phase B implements it as a row-fused block kernel. It must
// be BIT-IDENTICAL to running R independent SparseMatrixType::Multiply calls -- in
// particular R=1 must equal Multiply exactly -- so that a future multi-vector caller
// cannot perturb the byte-exact moment goldens.
//
// ASSERTS (==, not approx): for R in {1, 8} and several thread counts, column r of
// BatchMultiply(R, a, X, b, Y) equals Multiply(a, X[:,r], b, Y[:,r]) for every r,
// every component, with non-trivial complex a, b and a complex Hermitian operator.
//
// Pitfall: do not relax == to a tolerance; a mismatch means the block kernel diverged
//   from the single-vector kernel (a bug), not a tolerance to absorb.

#include <iostream>
#include <vector>
#include <complex>
#include <eigen3/Eigen/Sparse>
#ifdef _OPENMP
#include <omp.h>
#endif
#include "sparse_matrix.hpp"

using cd = std::complex<double>;

int main() {
    const int N = 777, bw = 4;            // odd N to exercise ragged thread partitions
    std::vector<Eigen::Triplet<cd, indexType>> trip;
    auto re = [](int i,int j){ return std::sin(0.3*i + 0.7*j + 1.0); };
    auto im = [](int i,int j){ return std::cos(0.2*i - 0.5*j + 0.4); };
    for (int i = 0; i < N; ++i) {
        trip.emplace_back(i, i, cd(1.5 + 0.01*i, 0.0));
        for (int d = 1; d <= bw; ++d) { int j=i+d; if(j>=N) continue;
            cd v(re(i,j), im(i,j)); trip.emplace_back(i,j,v); trip.emplace_back(j,i,std::conj(v)); }
    }
    Eigen::SparseMatrix<cd, Eigen::RowMajor, indexType> M(N, N);
    M.setFromTriplets(trip.begin(), trip.end()); M.makeCompressed();
    SparseMatrixType op; op.set_eigen_matrix(M);

    const cd a(0.9, -0.4), b(0.3, 0.7);
    const int Rs[] = {1, 8};
    const int threads[] = {1, 4, 8};
    int failures = 0;

    for (int R : Rs) {
        // column-major blocks X, Y0 (leading dim N): RHS r is the contiguous N-vector at r*N
        std::vector<cd> X(N*R), Y0(N*R);
        for (int r = 0; r < R; ++r)
            for (int i = 0; i < N; ++i) {
                X [r*N+i] = cd(std::sin(0.11*i + 0.2*r), std::cos(0.07*i - 0.3*r));
                Y0[r*N+i] = cd(std::cos(0.05*i + r),     std::sin(0.09*i + 0.1*r));
            }
        // reference: R independent Multiply calls
        std::vector<cd> ref = Y0;
        for (int r = 0; r < R; ++r)
            op.Multiply(a, X.data()+r*N, b, ref.data()+r*N);

        for (int nt : threads) {
#ifdef _OPENMP
            omp_set_num_threads(nt);
#endif
            std::vector<cd> Y = Y0;
            op.BatchMultiply(R, a, X.data(), b, Y.data());
            int mism = 0; double w = 0;
            for (int k = 0; k < N*R; ++k)
                if (Y[k].real()!=ref[k].real() || Y[k].imag()!=ref[k].imag()) {
                    ++mism; double d = std::abs(Y[k]-ref[k]); if (d>w) w=d; }
            std::cout << "R=" << R << " threads=" << nt << " mismatches=" << mism
                      << " max_abs_diff=" << w << std::endl;
            if (mism) ++failures;
        }
    }

    if (failures == 0) {
        std::cout << "PASS: BatchMultiply is bit-identical to R independent Multiply calls "
                     "(R in {1,8})." << std::endl;
        return 0;
    }
    std::cout << "FAIL: BatchMultiply diverged from per-vector Multiply -- bug, do not relax."
              << std::endl;
    return 1;
}
