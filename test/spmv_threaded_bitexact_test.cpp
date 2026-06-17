// spmv_threaded_bitexact_test -- safety net for the threaded SpMV (Phase A, W1).
//
// PROBLEM. The Chebyshev moment recursion's only O(nnz) hot loop is the leaf
// SpMV SparseMatrixType::Multiply (y = a*M*x + b*y). Phase A replaces the stock
// single-threaded Eigen product with a row-parallel OpenMP kernel. The exact-match
// moment goldens (chain1d_analytic, graphene_regression, greenwood_regression,
// hall_haldane, msd_ballistic) are only safe if that replacement does not perturb
// the bytes the recursion produces.
//
// TWO PROPERTIES PIN GOLDEN-SAFETY (this test asserts BOTH, with == not approx):
//
//  (1) THREAD-INVARIANCE for ANY scalars. Each y[i] is an independent reduction
//      over row i's nonzeros in stored order; distinct rows are independent writes.
//      So the result cannot depend on how rows are split across threads. Asserted by
//      comparing 1 vs 2/4/8 threads, bit-for-bit, including non-trivial complex a,b.
//
//  (2) BIT-IDENTITY TO EIGEN on the scalars the recursion actually uses. The kernel
//      computes a*(sum_p v_p x_p) + b*y[i]; Eigen's product folds the scalar
//      differently, so for a GENERAL complex a the two differ by ~1e-15 (IEEE
//      non-associativity -- a*(x+y) != a*x+a*y in floating point). That is HARMLESS
//      here because the KPM recursion only ever calls Multiply with a in {2,1,-1}
//      and b in {0,-1} (Iterate: 2*H~*T1 - 1*T0; SetInitVectors/operator apply:
//      1,0). Multiplying a double by an exact power-of-two / unit / zero is EXACT,
//      so a*(sum) == sum(a*x) bit-for-bit on every recursion path => the goldens
//      cannot move. Asserted by comparing Multiply to the Eigen reference a*M*x+b*y
//      for (a,b) in {(2,-1),(1,0),(1,-1),(-1,2)} at several thread counts.
//
// Pitfall: do NOT "fix" a failure by relaxing == to a tolerance -- that would let the
//   byte-exact goldens drift silently. A property-(1) failure means the kernel
//   introduced thread nondeterminism (a real bug); a property-(2) failure means a
//   recursion-relevant scalar path diverged from Eigen (a real bug). Fix the kernel.
// Note: a general-complex (a,b) is intentionally NOT required to match Eigen here;
//   the recursion never uses one, and demanding it would be demanding a property
//   Eigen's own SIMD product does not even hold under row-blocking.

#include <iostream>
#include <vector>
#include <complex>
#include <eigen3/Eigen/Sparse>
#include <eigen3/Eigen/Core>
#ifdef _OPENMP
#include <omp.h>
#endif
#include "sparse_matrix.hpp"

using cd = std::complex<double>;

static SparseMatrixType make_hermitian(int N, int bw) {
    std::vector<Eigen::Triplet<cd, indexType>> trip;
    auto re = [](int i,int j){ return std::sin(0.3*i + 0.7*j + 1.0); };
    auto im = [](int i,int j){ return std::cos(0.2*i - 0.5*j + 0.4); };
    for (int i = 0; i < N; ++i) {
        trip.emplace_back(i, i, cd(1.5 + 0.01*i, 0.0));
        for (int d = 1; d <= bw; ++d) {
            int j = i + d; if (j >= N) continue;
            cd v(re(i,j), im(i,j));
            trip.emplace_back(i, j, v);
            trip.emplace_back(j, i, std::conj(v));
        }
    }
    Eigen::SparseMatrix<cd, Eigen::RowMajor, indexType> M(N, N);
    M.setFromTriplets(trip.begin(), trip.end());
    M.makeCompressed();
    SparseMatrixType op; op.set_eigen_matrix(M);
    return op;
}

int main() {
    const int N = 1024, bw = 5;
    SparseMatrixType op = make_hermitian(N, bw);
    Eigen::SparseMatrix<cd, Eigen::RowMajor, indexType>& M = op.eigen_matrix();

    std::vector<cd> x(N), y0(N);
    for (int i = 0; i < N; ++i) {
        x[i]  = cd(std::sin(0.11*i + 0.2), std::cos(0.07*i - 0.3));
        y0[i] = cd(std::cos(0.05*i),       std::sin(0.09*i + 0.1));
    }
    Eigen::Map<const Eigen::Vector<cd, -1>> ex(x.data(), N), ey0(y0.data(), N);
    const int threads[] = {1, 2, 4, 8};
    int failures = 0;

    auto run = [&](cd a, cd b, int nt) {
#ifdef _OPENMP
        omp_set_num_threads(nt);
#endif
        std::vector<cd> y = y0;
        op.Multiply(a, x.data(), b, y.data());
        return y;
    };

    // ---- Property (1): thread-invariance for non-trivial complex (a,b) ----
    {
        const cd a(0.9, -0.4), b(0.3, 0.7);
        std::vector<cd> ref = run(a, b, 1);
        for (int nt : threads) {
            std::vector<cd> y = run(a, b, nt);
            int mism = 0;
            for (int i = 0; i < N; ++i)
                if (y[i].real()!=ref[i].real() || y[i].imag()!=ref[i].imag()) ++mism;
            std::cout << "[invariance] complex(a,b) threads=" << nt
                      << " vs 1-thread mismatches=" << mism << std::endl;
            if (mism) ++failures;
        }
    }

    // ---- Property (2): bit-identity to Eigen for recursion scalars ----
    const std::pair<cd,cd> cases[] = {
        {cd(2,0), cd(-1,0)},   // Iterate: 2*H~*T_m - T_{m-1}
        {cd(1,0), cd(0,0)},    // SetInitVectors / operator apply
        {cd(1,0), cd(-1,0)},
        {cd(-1,0), cd(2,0)},
    };
    for (auto& ab : cases) {
        Eigen::Vector<cd, -1> yref = ab.first * M * ex + ab.second * ey0;
        for (int nt : threads) {
            std::vector<cd> y = run(ab.first, ab.second, nt);
            int mism = 0; double w = 0;
            for (int i = 0; i < N; ++i)
                if (y[i].real()!=yref(i).real() || y[i].imag()!=yref(i).imag()) {
                    ++mism; double d = std::abs(y[i]-yref(i)); if (d>w) w=d;
                }
            std::cout << "[vs-Eigen] a=(" << ab.first.real() << ") b=(" << ab.second.real()
                      << ") threads=" << nt << " mismatches=" << mism
                      << " max_abs_diff=" << w << std::endl;
            if (mism) ++failures;
        }
    }

    if (failures == 0) {
        std::cout << "PASS: threaded SpMV is thread-invariant for all scalars and "
                     "bit-identical to Eigen on every recursion scalar (N=" << N << ")." << std::endl;
        return 0;
    }
    std::cout << "FAIL: " << failures << " property violations -- this is a bug, do NOT "
                 "relax to a tolerance; fix the kernel." << std::endl;
    return 1;
}
