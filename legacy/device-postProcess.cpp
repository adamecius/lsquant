// [LEGACY -- NOT BUILT] previous-working device-postProcess. Reference only.
// Supported path: lsquant <cmd>  |  port + golden-gate before any reuse.
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <cmath>
#include <iomanip>
#include <map>
#include <stdexcept>


#include <fstream>
#include <sstream>
#include <string>
#include <unordered_map>
#include <stdexcept>
#include <algorithm>
#include <regex>

#include "Graphene_NNN.hpp"

#ifdef _OPENMP
#include <omp.h>
#endif


static constexpr double H_SI    = 6.62607015e-34;   // J·s
static constexpr double HBAR_SI = 1.05457182e-34;   // J·s
static constexpr double Q       = 1.60217663e-19;   // C
static constexpr double PI      = 3.14159265358979;

// hbar in eV·s
static constexpr double HBAR_EV = HBAR_SI / Q;
// Quantum of conductance G0 = 2e²/h
static constexpr double G0      = 2.0 * Q * Q / H_SI;


// ─────────────────────────────────────────────────────────────────────────────
// Low-noise N=11 central-difference differentiator
// (Holoborodko, http://www.holoborodko.com/pavel)
// ─────────────────────────────────────────────────────────────────────────────
static std::vector<double> gradient_lownoise_N11(const std::vector<double>& y, double dx)
{
    int n = static_cast<int>(y.size());
    std::vector<double> d(n, 0.0);
    if (n < 2) return d;

    for (int i = 0; i < n; ++i) {
        if      (i == 0)     d[i] = (y[i+1]-y[i]) / dx;
        else if (i == n-1)   d[i] = (y[i]-y[i-1]) / dx;
        else if (i==1||i==n-2)
            d[i] = (y[i+1]-y[i-1]) / (2*dx);
        else if (i==2||i==n-3)
            d[i] = (2*(y[i+1]-y[i-1]) + (y[i+2]-y[i-2])) / (8*dx);
        else if (i==3||i==n-4)
            d[i] = (5*(y[i+1]-y[i-1]) + 4*(y[i+2]-y[i-2]) + (y[i+3]-y[i-3])) / (32*dx);
        else if (i==4||i==n-5)
            d[i] = (14*(y[i+1]-y[i-1]) + 14*(y[i+2]-y[i-2])
                  + 6*(y[i+3]-y[i-3]) + (y[i+4]-y[i-4])) / (128*dx);
        else
            d[i] = (42*(y[i+1]-y[i-1]) + 48*(y[i+2]-y[i-2])
                  + 27*(y[i+3]-y[i-3]) + 8*(y[i+4]-y[i-4])
                  + (y[i+5]-y[i-5])) / (512*dx);
    }
    return d;
}


// ─────────────────────────────────────────────────────────────────────────────
// Cumulative trapezoidal integration  (matches scipy.integrate.cumtrapz)
// ─────────────────────────────────────────────────────────────────────────────
static std::vector<double> cumtrapz(const std::vector<double>& y,
                                    const std::vector<double>& x)
{
    int n = static_cast<int>(y.size());
    std::vector<double> out(n, 0.0);
    for (int i = 1; i < n; ++i)
        out[i] = out[i-1] + 0.5*(y[i]+y[i-1])*(x[i]-x[i-1]);
    return out;
}


// ─────────────────────────────────────────────────────────────────────────────
// Linear interpolation of src defined on x_src, evaluated at x_dst
// ─────────────────────────────────────────────────────────────────────────────
static std::vector<double> interp(const std::vector<double>& x_dst,
                                  const std::vector<double>& x_src,
                                  const std::vector<double>& y_src)
{
    int n = static_cast<int>(x_dst.size());
    std::vector<double> out(n);
    for (int i = 0; i < n; ++i) {
        double xq = x_dst[i];
        // Find bracket
        auto it = std::lower_bound(x_src.begin(), x_src.end(), xq);
        if (it == x_src.begin()) { out[i] = y_src[0]; continue; }
        if (it == x_src.end())   { out[i] = y_src.back(); continue; }
        int j = static_cast<int>(it - x_src.begin());
        double t = (xq - x_src[j-1]) / (x_src[j] - x_src[j-1]);
        out[i] = y_src[j-1] + t*(y_src[j] - y_src[j-1]);
    }
    return out;
}




// ─────────────────────────────────────────────────────────────────────────────
// Write a 2-column file
// ─────────────────────────────────────────────────────────────────────────────
static void write2col(const std::string& fname,
                      const std::vector<double>& a,
                      const std::vector<double>& b)
{
    std::ofstream f(fname);
    if (!f) throw std::runtime_error("Cannot write: " + fname);
    f << std::scientific << std::setprecision(8);
    for (size_t i = 0; i < a.size(); ++i)
        f << a[i] << "  " << b[i] << "\n";
    std::cout << "  Written " << fname << "\n";
}

// Write a matrix where rows = energies, columns = time steps
// Header line contains the time values
static void writeMatrix(const std::string& fname,
                        const std::vector<double>& E,
                        const std::vector<double>& time,
                        const std::vector<std::vector<double>>& M)
{
    std::ofstream f(fname);
    if (!f) throw std::runtime_error("Cannot write: " + fname);
    f << std::scientific << std::setprecision(8);

    // Header: time values
    f << "#E";
    for (double t : time) f << "  " << t;
    f << "\n";

    // Rows: one per energy
    for (size_t iE = 0; iE < E.size(); ++iE) {
        f << E[iE];
        for (size_t it = 0; it < time.size(); ++it)
            f << "  " << M[iE][it];
        f << "\n";
    }
    std::cout << "  Written " << fname << "\n";
}





// ── Kernels ───────────────────────────────────────────────────────────────────
std::vector<double> JacksonKernel(int N)
{
    std::vector<double> g(N);
    for (int n = 0; n < N; ++n)
        g[n] = ((N - n + 1.0) * std::cos(PI * n / (N + 1.0))
              + std::sin(PI * n / (N + 1.0)) / std::tan(PI / (N + 1.0)))
             / (N + 1.0);
    return g;
}

std::vector<double> LorentzKernel(int N, double lambda)
{
    std::vector<double> g(N);
    for (int n = 0; n < N; ++n)
        g[n] = std::sinh(lambda * (1.0 - double(n) / N)) / std::sinh(lambda);
    return g;
}

// ── Spectral reconstruction ───────────────────────────────────────────────────
// nMom is modified in-place to the truncated value (mirrors Fortran behaviour).
void densite(int& nMom,
             const std::vector<double>& mu,
             double Emin, double dE, int nE,
             double eta, double a, 
             std::vector<double>& dos)
{
    const double lambda  = PI;
    int nMom_eta = static_cast<int>(lambda / (eta / a ));
    if (nMom_eta > nMom) {
        nMom_eta = nMom;
        std::cout << "Warning: not enough moments --> broadening = "
                  << lambda / nMom * a << "\n";
    }
    nMom = nMom_eta;
    std::cout << "Nmom = " << nMom << "\n";

    auto g = JacksonKernel(nMom);

    #pragma omp parallel for schedule(static)
    for (int i = 0; i <= nE; ++i) {
        double E    = (Emin + i * dE) / a ;
        double suma = g[0] * mu[0];
        for (int n = 1; n < nMom; ++n)
            suma += 2.0 * g[n] * mu[n] * std::cos(n * std::acos(E));
        dos[i] = suma / PI / std::sqrt(1.0 - E * E);
    }
}



// ── Main ──────────────────────────────────────────────────────────────────────
int main()
{
    // ── Step 0: parameters ───────────────────────────────────────────────────
  Graphene_NNN graph(0, false, true);

    /*
    auto P = readNamelist("params_cpp.txt");
    
    int    nrecurs = 2215;//get<int>   (P, "nrecurs");
    int    nt      = 10;//get<int>   (P, "nt");
    double tstep   = 100;//get<double>(P, "tstep") / 1000.0;   // fs -> ps
    double Emin    = -1.0;//get<double>(P, "Emin");
    double Emax    = 1.0;//get<double>(P, "Emax");
    double dE      = 0.001;//get<double>(P, "dE");
    
    double eta     = 0.00026;//get<double>(P, "eta");
    double a       = 5;//get<double>(P, "a");
    double teV     = graph.params().teV;//get<double>(P, "teV");
    */

  


  double a;
  double b;
  int    nMom_DOS;
  int    size;
  int    nt;
  double tmax;
  double tstep;
  double Emin;
  double Emax;
  double dE;
  double eta;

  std::ifstream Input;
  Input.open("postProcess.dat");
  std::string dump;

  Input>>dump;
  Input>>a;

  Input>>dump;
  Input>>b;

  Input>>dump;
  Input>>nMom_DOS;

  Input>>dump;
  Input>>size;

  Input>>dump;
  Input>>nt;

  Input>>dump;
  Input>>tmax;

  Input>>dump;
  Input>>tstep;  

  Input>>dump;
  Input>>Emin;

  Input>>dump;
  Input>>Emax;

  Input>>dump;
  Input>>dE;

  Input>>dump;
  Input>>eta;

  
  
    
    std::cout << "tstep (ps) = " << tstep << "\n";

    const int nE = static_cast<int>((Emax - Emin) / dE);

    // ── Step 1: DOS from mu_01.txt ───────────────────────────────────────────
    std::vector<double> mu(nMom_DOS), mux(nMom_DOS), muy(nMom_DOS);

    {
        std::ifstream f("mu_01.txt");
        if (!f) throw std::runtime_error("Cannot open mu_01.txt");
        int idx; double im;
        for (int j = 0; j < nMom_DOS; ++j)
            f >> idx >> mu[j] >> im;
    }

    std::vector<double> dos(nE + 1, 0.0);
    std::vector<double> E(nE + 1, 0.0);
    std::vector<double> time(nt + 1);
    
    for (int i = 0; i <= nt; ++i) time[i] = i * tstep;
    for (int j = 0; j <= nE; ++j) E[j] = (Emin + j * dE);
    
    
    densite(nMom_DOS, mu, Emin, dE, nE, eta, a, dos);

    {
        std::ofstream f("dos.txt");
        for (int j = 0; j <= nE; ++j){
	 
	  f << std::setw(16) << (Emin + j * dE)
              << std::setw(16) << dos[j] / a  << "\n";
	}
    }

    // ── Step 2: MSD for each timestep, stored as columns ────────────────────
    // Layout: all_dx2[j][i] = dx2 at energy index j, timestep i (0-based)
    std::vector<std::vector<double>> all_dx2(nE + 1, std::vector<double>(nt, 0.0));
    std::vector<std::vector<double>> all_dy2(nE + 1, std::vector<double>(nt, 0.0));
    std::vector<std::vector<double>> dL2(nE + 1, std::vector<double>(nt, 0.0));
 
    for (int i = 1; i <= nt; ++i) {
        if (i % 10 == 0)
            std::cout << "time (ps) = " << i * tstep << "\n";

        // Build filename muxy_XXXXX.txt
        std::ostringstream ss;
        ss << "muxy_" << std::setw(5) << std::setfill('0') << i << ".txt";
        std::ifstream f(ss.str());
        if (!f) throw std::runtime_error("Cannot open " + ss.str());

        char   ctemp;
        double normx, normy, r1, r2;
        int    idx;
        f >> ctemp >> normx >> normy;
        for (int j = 0; j < nMom_DOS; ++j)
            f >> idx >> mux[j] >> r1 >> muy[j] >> r2;

        std::vector<double> dx2(nE + 1, 0.0), dy2(nE + 1, 0.0);
        int nMom_x = nMom_DOS, nMom_y = nMom_DOS;
        densite(nMom_x, mux, Emin, dE, nE, eta, a, dx2);
        densite(nMom_y, muy, Emin, dE, nE, eta, a, dy2);

        for (int j = 0; j <= nE; ++j) {
            all_dx2[j][i-1] = dx2[j] * normx * normx / dos[j];
            all_dy2[j][i-1] = dy2[j] * normy * normy / dos[j];
	    dL2[j][i-1] = all_dx2[j][i-1] + all_dy2[j][i-1];
        }
    }

    
    // ── Step 3: write output — one file per component, timesteps as columns ──
    // Header row: E  t=tstep  t=2*tstep  ...
    auto writeOutput = [&](const std::string& fname,
                           const std::vector<std::vector<double>>& data)
    {
        std::ofstream f(fname);
        if (!f) throw std::runtime_error("Cannot open " + fname);
        // Header
        f << std::setw(16) << "# E(eV)";
        for (int i = 1; i <= nt; ++i)
            f << std::setw(16) << (i * tstep);
        f << "\n";
        // Data
        for (int j = 0; j <= nE; ++j) {
            f << std::setw(16) << (Emin + j * dE);
            for (int i = 0; i < nt; ++i)
                f << std::setw(16) << data[j][i];
            f << "\n";
        }
    };

    writeOutput("dX2.txt", all_dx2);
    writeOutput("dY2.txt", all_dy2);

    std::cout << "Done. Output written to dX2.txt and dY2.txt\n";


    // ── COMPUTE PHYSICAL QUANTITIES ─────────────────

    //MISSING: Here we should get the unit cell variables in a generic case somehow    
    double A_m2  = 1.0;
    /*
    double A1x = graph.unit_cell().A1x, A1y = graph.unit_cell().A1y, A2x = graph.unit_cell().A2x, A2y = graph.unit_cell().A2y, acc = graph.unit_cell().acc;
    int nuc = graph.unit_cell().nuc;
    
    double L1   = std::sqrt(A1x*A1x + A1y*A1y);
    double L2   = std::sqrt(A2x*A2x + A2y*A2y);
    // Area per atom [m²]
    double A_m2  = (A1x*A2y - A1y*A2x) * 1e-20 / nuc;
    */
    // Convert DOS from file units to [1/J/m²]: dos = 2*dos_file / A_m2 / Q
    for (size_t i = 0; i < dos.size();     ++i) dos[i]     = 2.0*dos[i]     / A_m2 / Q;



    
    // ── Compute D(t), sigma(t), and semiclassical quantities ─────────────────
    std::vector<std::vector<double>> D    (nE+1, std::vector<double>(nt+1, 0.0));
    std::vector<std::vector<double>> sigma(nE+1, std::vector<double>(nt+1, 0.0));
    std::vector<double> D_sc    (nE+1, 0.0);
    std::vector<double> sigma_sc(nE+1, 0.0);
    std::vector<double> tp      (nE+1, 0.0);

    //MISSING: We need a generic way to compute the fermi velocity. Joaquin's formula??
    // Graphene Fermi velocity  vF = (3/2) * acc[m] * t[eV] / hbar[eV·s]
    double vf = 1;
    //double vf    = 1.5 * (acc * 1e-10) * graph.params().teV / HBAR_EV;


    // Carrier density n2D [1/m²] via cumulative trapezoidal integration
    // Subtract CNP value so n2D = 0 at the Dirac point

    int iCNP = static_cast<int>(
    std::min(dos.begin(), dos.end()) - dos.begin());
    
    auto n2D_raw = cumtrapz(dos, E);   // dosfine in 1/J/m², E in eV → integrate over eV→multiply by Q
    // Actually integrate over eV: n [1/m²] = ∫ dos[1/J/m²] * Q[J/eV] dE[eV]
    std::vector<double> n2D(n2D_raw.size());
    for (size_t i = 0; i < n2D.size(); ++i)
        n2D[i] = (n2D_raw[i] - n2D_raw[iCNP]) * Q;



    
    for (int iE = 0; iE <= nE; ++iE) {
        // Diffusivity D = d(dL2)/dt / 4  [m²/s]
        // Input dL2 in Å², time in fs → convert to m²/s: ×1e-20/1e-15 = ×1e-5
	std::vector<double> dL2_slice(nt+1);
	for (int it = 0; it < nt; ++it) dL2_slice[it] = dL2[iE][it];

	  
        auto grad      = gradient_lownoise_N11(dL2_slice, tstep);
        for (int it = 0; it <= nt; ++it)
            D[iE][it] = grad[it] / 4.0 * 1e-20 / 1e-15;
        D[iE][0] = 0.0;

        // Conductivity σ(t) = e² * DOS * D  [1/Ω]
        for (int it = 0; it <= nt; ++it)
            sigma[iE][it] = Q * Q * dos[iE] * D[iE][it];

        // Semiclassical: max of D and σ over time
        D_sc[iE]     = *std::max(D[iE].begin(),     D[iE].end());
        sigma_sc[iE] = *std::max(sigma[iE].begin(), sigma[iE].end());

        // Momentum relaxation time τ_p = 2*D_sc / vF²  [fs]
        tp[iE] = 2.0 * D_sc[iE] / (vf * vf) * 1e15;
    }

    // ── Mobility [cm²/V·s] from generalised Einstein relation ─────────────────
    std::vector<double> mub(nE+1, 0.0);
    for (int iE = 0; iE <= nE; ++iE) {
        double dndE = dos[iE];   // [1/J/m²]
        double n    = std::abs(n2D[iE]);
        if (n > 0.0)
            mub[iE] = D_sc[iE] * dndE / n * Q * 1e4;  // ×1e4: m²→cm²
    }





    // ── Seebeck coefficient [μV/K] ────────────────────────────────────────────
    // Uses generalised Mahan/BTE transport coefficients:
    //   K_n[E] = (σ_sc ★ E^n·dfdE)[E] · dE       (same-mode convolution)
    //   S[E]   = (1/T) · K1[E]/K0[E] · 1e6
    //
    // dfdE[i] = -df/dE = 1 / (cosh(E[i]/(2kT)))² / (4kT)   (Fermi window)
    // h_Planck in SI = 6.62607015e-34 J·s
    // ─────────────────────────────────────────────────────────────────────────
    const double T_K    = 300.0;                        // temperature [K]
    const double e      = 1.602e-19;                    // electron volt in Joules
    const double kB     = 1.380649e-23 / e;             // Boltzmann [e/K]
    const double kT     = kB * T_K;                     // thermal energy [J]
    const double dE_J   = (E.back() - E.front())  / std::max(nE, 1)       // energy step in eV
     ;        // E array is in eV → ×Q

    // Fermi window: dfdE[i] = 1 / cosh²(E[i]/(2kT)) / (4kT)
    std::vector<double> dfdE(nE+1);
    for (int iE = 0; iE <= nE; ++iE) {
      double E_J  = E[iE];
      double arg  = E_J / (2.0 * kT);
      double ch   = std::cosh(arg);
      dfdE[iE]    = 1.0 / (ch * ch) / (4.0 * kT);
    }


    
    // Same-mode discrete convolution (linear, centred)
    //   K_n[iE] = Σ_j  (E[j]^n · dfdE[j]) · sig_norm[iE - j + nE/2] · dE_J
    // Equivalent to numpy's convolve(...,'same'): output index i gets
    // contributions from kernel index k and signal index i-k+(N-1)/2.
    std::vector<double> K0(nE+1, 0.0), K1(nE+1, 0.0), K2(nE+1, 0.0);
    const int half = nE / 2;
    for (int iE = 0; iE <= nE; ++iE) {
      double k0 = 0.0, k1 = 0.0, k2 = 0.0;
      for (int j = 0; j <= nE; ++j) {
        int si = iE - j + half;                     // signal index
        if (si < 0 || si > nE) continue;
        double E_J  = E[j] ;
        double kern = dfdE[j] * sigma_sc[si] * dE_J;
        k0 += kern;                                 // E^0 · dfdE ★ σ
        k1 += E_J * kern;                           // E^1 · dfdE ★ σ
        k2 += E_J * kern;                           // E^2 · dfdE ★ σ        
      }
      K0[iE] = k0;
      K1[iE] = k1;
      K2[iE] = k2;
    }

    // Seebeck S[iE] = (1/T) · K1/K0 · 1e6  [μV/K]
    std::vector<double> S_seebeck(nE+1, 0.0);
    for (int iE = 0; iE <= nE; ++iE)
      if (std::abs(K0[iE]) > 0.0)
        S_seebeck[iE] = (1.0 / T_K) * (K1[iE] / K0[iE]) * 1e6;



    // Thermal conductivity kappa[iE] = (2/hT) · (K2-K1^2/K0) · 1e6  Units? [μeV/ K]?
    std::vector<double> kappa(nE+1, 0.0);
    for (int iE = 0; iE <= nE; ++iE)
      if (std::abs(K0[iE]) > 0.0)
        kappa[iE] = (1.0 / ( 2 * PI * HBAR_EV * T_K)) * (K2[iE] - K1[iE] * K1[iE] / K0[iE]) * 1e6;










    

    // ── Write output files ────────────────────────────────────────────────────
    std::cout << "\nWriting output files...\n";

    // DOS: E [eV]  dos [1/eV/nm²]
    {
        std::vector<double> dos_plot(nE+1);
        for (int i = 0; i <= nE; ++i)
            dos_plot[i] = dos[i] * Q / 1e18;   // 1/J/m² → 1/eV/nm²
        write2col("dos_plot.txt", E, dos_plot);
    }

    // Semiclassical conductivity [G0]
    {
        std::vector<double> ssc(nE+1);
        for (int i = 0; i <= nE; ++i) ssc[i] = sigma_sc[i] / G0;
        write2col("sigma_sc.txt", E, ssc);
    }

    // Momentum relaxation time [fs]
    write2col("tp.txt", E, tp);

    // Mobility [cm²/V·s]  (also writes n2D in cm⁻²)
    {
        std::ofstream f("mu.txt");
        f << std::scientific << std::setprecision(8);
        for (int i = 0; i <= nE; ++i)
            f << E[i] << "  " << n2D[i]/1e4 << "  " << mu[i] << "\n";
        std::cout << "  Written mu.txt\n";
    }

    
    // Diffusivity vs time [cm²/s]
    {
        std::vector<std::vector<double>> D_cm(nE+1, std::vector<double>(nt+1));
        for (int iE = 0; iE <= nE; ++iE)
            for (int it = 0; it <= nt; ++it)
                D_cm[iE][it] = D[iE][it] * 1e4;   // m²/s → cm²/s
	writeMatrix("D_vs_t.txt", E, time, D_cm);
    }

    // Conductivity vs time [G0]
    {
        std::vector<std::vector<double>> sig_G0(nE+1, std::vector<double>(nt+1));
        for (int iE = 0; iE <= nE; ++iE)
            for (int it = 0; it <= nt; ++it)
                sig_G0[iE][it] = sigma[iE][it] / G0;
        writeMatrix("sigma_vs_t.txt", E, time, sig_G0);
    }



    // MSD vs time: rows=energies, cols=time steps  [µm²]
    {
        std::vector<std::vector<double>> MSD_um(nE, std::vector<double>(nt+1));
        for (int iE = 0; iE < nE; ++iE)
            for (int it = 0; it <= nt; ++it)
                MSD_um[iE][it] = dL2[iE][it] / (1e4*1e4);   // Å² → µm²
        writeMatrix("MSD_vs_t.txt", E, time, MSD_um);
    }



        // Seebeck coefficient [µV/T]  (also writes n2D in cm⁻²)
    {
        std::ofstream f("Seebeck.txt");
        f << std::scientific << std::setprecision(8);
        for (int i = 0; i <= nE; ++i)
            f << E[i] << "  " << S_seebeck[i] << "\n";
        std::cout << "  Written Seebeck.txt\n";
    }



    // Thermal conductivity [Units?]
    {
        std::ofstream f("thermal_conductivity.txt");
        f << std::scientific << std::setprecision(8);
        for (int i = 0; i <= nE; ++i)
            f << E[i] << "  " << kappa[i] << "\n";
        std::cout << "  Written thermal_conductivity.txt\n";
    }








    return 0;
}
