#ifndef LSQUANT_CORE_ALGEBRA_FFT_MOMENT_CONTRACTION_HPP_
#define LSQUANT_CORE_ALGEBRA_FFT_MOMENT_CONTRACTION_HPP_
#include "core/algebra/types.hpp"
namespace lsquant {
// Item 5c: evaluate sum_mn mu[m,n] cos(m*theta) exp(-i*n*theta), and
// the same with an extra n, at theta=acos(cos((2*j+0.5)*pi/points)).
// Output is two rows of `points` complex values. Collapse frequencies in fixed
// m,n order, then four single-thread FFTs: O(M0*M1 + N log N) work and
// O(M0+M1+N) workspace (O(N) on the usual N=max(M0,M1) grid).
// Python owns coefficient factors, damping, units and integration. Tester must
// compare against table contraction on the SAME nodes, including rectangular
// tensors and odd/even N; FFT reassociation is not promised byte-equal to scalar.
std::vector<Scalar> FftMomentContraction(const Scalar* moments, Index m0, Index m1,
                                        Index points);
}  // namespace lsquant
#endif
