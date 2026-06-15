// Phase 6A WIRED contract: the recursion's central type is now the operator interface, and the
// scaffolding placeholders stay loud. This replaces phase6_scaffolding (which only checked that
// the headers compiled). It pins three things:
//   (1) COMPILE-TIME: chebyshev::Moments routes through lsquant::HamiltonianOp / MetricPolicy --
//       Moments::op() returns HamiltonianOp&, metric() returns MetricPolicy&, and the live
//       SparseHamiltonianOp implements the interface;
//   (2) OrthogonalMetric is the REAL identity default (no-op S^{-1});
//   (3) the placeholders are LOUD -- MatrixFreeHamiltonianOp / NonOrthogonalMetric THROW.
#include <iostream>
#include <complex>
#include <vector>
#include <stdexcept>
#include <type_traits>
#include "hamiltonian_op.hpp"
#include "metric_policy.hpp"
#include "chebyshev_moments.hpp"   // the migrated core

// (1) The core speaks the interface (a compile-time assertion -- the header is wired in).
static_assert(std::is_same<decltype(std::declval<chebyshev::Moments&>().op()),
                           lsquant::HamiltonianOp&>::value,
              "Moments::op() must return lsquant::HamiltonianOp&");
static_assert(std::is_same<decltype(std::declval<chebyshev::Moments&>().metric()),
                           lsquant::MetricPolicy&>::value,
              "Moments::metric() must return lsquant::MetricPolicy&");
static_assert(std::is_base_of<lsquant::HamiltonianOp, lsquant::SparseHamiltonianOp>::value,
              "SparseHamiltonianOp must implement HamiltonianOp (the live adapter)");

static int g_fail = 0;
static void check(bool ok, const std::string& what) {
    std::cout << (ok ? "  OK   " : "  FAIL ") << what << "\n";
    if (!ok) g_fail = 1;
}
template <class F>
static bool throws(F f) { try { f(); return false; } catch (const std::logic_error&) { return true; } catch (...) { return false; } }

int main() {
    typedef std::complex<double> cd;

    check(true, "Moments::op() returns HamiltonianOp& (compile-time, see static_assert)");

    // (2) OrthogonalMetric: identity S^{-1}, leaves the vector untouched.
    std::vector<cd> v;
    v.push_back(cd(1,2)); v.push_back(cd(-3,4)); v.push_back(cd(0,0));
    const std::vector<cd> before = v;
    lsquant::OrthogonalMetric ortho;
    ortho.apply_inverse_metric(v);
    check(ortho.is_orthogonal(), "OrthogonalMetric.is_orthogonal() == true");
    check(v == before,           "OrthogonalMetric.apply_inverse_metric is the identity");

    // (3a) NonOrthogonalMetric placeholder: reports non-orthogonal AND throws (not a silent no-op).
    lsquant::NonOrthogonalMetric nonortho;
    check(!nonortho.is_orthogonal(), "NonOrthogonalMetric.is_orthogonal() == false");
    std::vector<cd> w(3, cd(1,0));
    check(throws([&]{ nonortho.apply_inverse_metric(w); }),
          "NonOrthogonalMetric.apply_inverse_metric throws (placeholder, not implemented)");

    // (3b) MatrixFreeHamiltonianOp placeholder: throws on use.
    lsquant::MatrixFreeHamiltonianOp mf;
    check(throws([&]{ mf.dim(); }),
          "MatrixFreeHamiltonianOp.dim() throws (placeholder, not implemented)");
    std::vector<cd> x(3, cd(1,0)), y(3, cd(0,0));
    check(throws([&]{ mf.multiply(cd(1,0), x, cd(0,0), y); }),
          "MatrixFreeHamiltonianOp.multiply() throws (placeholder, not implemented)");

    if (g_fail == 0) { std::cout << "PASS: Phase-6A core is wired to HamiltonianOp; placeholders are loud\n"; return 0; }
    std::cerr << "FAIL: Phase-6A wired contract violated\n";
    return 1;
}
