// #2 -- Reconstructed DOS golden: the KPM-reconstructed DOS of the clean 1D chain vs the
// closed form rho(E) = 1/(pi sqrt(4 - E^2)) (oracle lsquant_chain_reference.dos, gamma=1),
// at the verified interior energies, with relative error < 1%.
//
// Usage: dos_reconstructed_test <DOS_reconstructed.dat> <TOL>
//
// The .dat is the per-site DOS averaged over N_r=128 single-vector runs (M=64, exact
// bounds [-2,2], Jackson kernel) -- see dos_reconstructed.sh. Reference comes from the
// oracle closed form, never from another run. Checked away from the band edges (the edges
// diverge and KPM smears them); the verified targets are |E| <= 1.5.
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <string>

int main(int argc, char** argv) {
    if (argc < 3) { std::cerr << "usage: " << argv[0] << " <DOS.dat> <TOL>\n"; return 2; }
    const double TOL = std::stod(argv[2]);

    std::ifstream f(argv[1]);
    if (!f) { std::cerr << "cannot open DOS file: " << argv[1] << "\n"; return 2; }
    std::vector<double> E, rho;
    double e, r;
    while (f >> e >> r) { E.push_back(e); rho.push_back(r); }
    if (E.size() < 2) { std::cerr << "DOS file too short\n"; return 2; }

    auto interp = [&](double x) {           // linear interpolation onto the E grid
        for (size_t i = 1; i < E.size(); ++i)
            if (x <= E[i]) {
                double t = (x - E[i-1]) / (E[i] - E[i-1]);
                return rho[i-1] + t * (rho[i] - rho[i-1]);
            }
        return rho.back();
    };
    auto ref = [](double x) { return 1.0 / (M_PI * std::sqrt(4.0 - x*x)); };

    const double pts[] = {0.0, 0.5, 1.0, 1.5};
    double maxRel = 0.0; double atE = 0.0;
    bool ok = true;
    for (double p : pts) {
        double got = interp(p), want = ref(p), rel = std::fabs(got - want) / want;
        std::cout << "  E=" << p << ": rho=" << got << " ref=" << want
                  << " rel=" << rel*100 << "%\n";
        if (rel > maxRel) { maxRel = rel; atE = p; }
        if (rel >= TOL) ok = false;
    }
    std::cout << "max rel error = " << maxRel*100 << "% at E=" << atE
              << "   TOL = " << TOL*100 << "%\n";
    if (ok) { std::cout << "PASS: reconstructed DOS matches 1/(pi sqrt(4-E^2)) within TOL\n"; return 0; }
    std::cerr << "FAIL: reconstructed DOS exceeds TOL\n";
    return 1;
}
