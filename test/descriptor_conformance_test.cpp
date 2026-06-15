// Phase 2 -- descriptor schema conformance.
//
// lsquant's .desc reader (include/operator_descriptor.hpp) is a DUPLICATED-BUT-VERSIONED
// mirror of the wannier2sparse producer schema (no shared submodule -- keeps the build
// decoupled). This test closes the duplication risk: it parses committed, w2s-PRODUCED .desc
// files and asserts every field round-trips, so reader and writer can't drift silently.
//
// Usage: descriptor_conformance_test <HAM.desc> <VX.desc>
#include <iostream>
#include <string>
#include <cmath>
#include "operator_descriptor.hpp"

static int g_fail = 0;
static void check(bool ok, const std::string& what) {
    std::cout << (ok ? "  OK   " : "  FAIL ") << what << "\n";
    if (!ok) g_fail = 1;
}

int main(int argc, char** argv) {
    if (argc < 3) { std::cerr << "usage: " << argv[0] << " <HAM.desc> <VX.desc>\n"; return 2; }

    // Hamiltonian descriptor: identity + Lanczos bounds.
    lsquant::OperatorDescriptor h = lsquant::read_descriptor(argv[1]);
    std::cout << "[HAM] " << argv[1] << "\n";
    check(h.valid,                       "file parsed");
    check(h.observable == "hamiltonian", "observable == hamiltonian (got '" + h.observable + "')");
    check(h.units == "eV",               "units == eV (got '" + h.units + "')");
    check(!h.provenance.empty(),         "provenance non-empty");
    check(h.has_bounds,                  "carries spectral bounds");
    check(std::fabs(h.a + 3.0) < 1e-6,   "spectral_min ~ -3 (graphene bandwidth)");
    check(std::fabs(h.b - 3.0) < 1e-6,   "spectral_max ~ +3 (graphene bandwidth)");

    // Velocity descriptor: identity, NO bounds.
    lsquant::OperatorDescriptor v = lsquant::read_descriptor(argv[2]);
    std::cout << "[VX] " << argv[2] << "\n";
    check(v.valid,                    "file parsed");
    check(v.observable == "velocity", "observable == velocity (got '" + v.observable + "')");
    check(v.component == "X",         "component == X (got '" + v.component + "')");
    check(!v.has_bounds,              "velocity carries NO bounds");

    // Absent file -> valid == false (graceful).
    lsquant::OperatorDescriptor missing = lsquant::read_descriptor("/no/such/file.desc");
    check(!missing.valid, "absent file -> valid == false");

    if (g_fail == 0) { std::cout << "PASS: lsquant reader conforms to the w2s .desc schema (v"
                                 << lsquant::DESCRIPTOR_SCHEMA_VERSION << ")\n"; return 0; }
    std::cerr << "FAIL: descriptor reader/producer schema mismatch\n";
    return 1;
}
