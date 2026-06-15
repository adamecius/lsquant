// Phase 6 scaffolding contract (interfaces only -- the core is not migrated yet).
//
// Pins two things so the placeholders cannot rot or be mistaken for working code:
//   (1) the REAL default behaviours -- OrthogonalMetric is the identity S^{-1};
//   (2) the PLACEHOLDERS are loud -- MatrixFreeHamiltonianOp and NonOrthogonalMetric THROW
//       (a silent no-op would be a correctness trap when the hot loop is migrated in Phase 6).
#include <iostream>
#include <complex>
#include <vector>
#include <stdexcept>
#include "hamiltonian_op.hpp"
#include "metric_policy.hpp"

static int g_fail = 0;
static void check(bool ok, const std::string& what) {
    std::cout << (ok ? "  OK   " : "  FAIL ") << what << "\n";
    if (!ok) g_fail = 1;
}
template <class F>
static bool throws(F f) { try { f(); return false; } catch (const std::logic_error&) { return true; } catch (...) { return false; } }

int main() {
    typedef std::complex<double> cd;

    // (1) OrthogonalMetric: identity S^{-1}, leaves the vector untouched.
    std::vector<cd> v;
    v.push_back(cd(1,2)); v.push_back(cd(-3,4)); v.push_back(cd(0,0));
    const std::vector<cd> before = v;
    lsquant::OrthogonalMetric ortho;
    ortho.apply_inverse_metric(v);
    check(ortho.is_orthogonal(), "OrthogonalMetric.is_orthogonal() == true");
    check(v == before,           "OrthogonalMetric.apply_inverse_metric is the identity");

    // (2a) NonOrthogonalMetric placeholder: reports non-orthogonal AND throws (not a silent no-op).
    lsquant::NonOrthogonalMetric nonortho;
    check(!nonortho.is_orthogonal(), "NonOrthogonalMetric.is_orthogonal() == false");
    std::vector<cd> w(3, cd(1,0));
    check(throws([&]{ nonortho.apply_inverse_metric(w); }),
          "NonOrthogonalMetric.apply_inverse_metric throws (placeholder, not implemented)");

    // (2b) MatrixFreeHamiltonianOp placeholder: throws on use.
    lsquant::MatrixFreeHamiltonianOp mf;
    check(throws([&]{ mf.dim(); }),
          "MatrixFreeHamiltonianOp.dim() throws (placeholder, not implemented)");
    std::vector<cd> x(3, cd(1,0)), y(3, cd(0,0));
    check(throws([&]{ mf.multiply(cd(1,0), x, cd(0,0), y); }),
          "MatrixFreeHamiltonianOp.multiply() throws (placeholder, not implemented)");

    if (g_fail == 0) { std::cout << "PASS: Phase-6 interfaces -- orthogonal default works, placeholders are loud\n"; return 0; }
    std::cerr << "FAIL: Phase-6 scaffolding contract violated\n";
    return 1;
}
