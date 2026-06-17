#!/usr/bin/env python3
"""Validate the random-phase (stochastic-trace) spin relaxation of tutorial 06
against (i) Elliott-Yafet theory and (ii) the exact trace.

Theory. For a weak random spin-flip field eta_i with <eta>=0, <eta^2>=Dsf^2, Fermi's
golden rule gives a spin-flip rate proportional to the perturbation variance:
    1/T1  ~  Dsf^2  (Elliott-Yafet).
So a log-log fit of 1/T1 vs Dsf must have slope ~ 2.

Method. 1/T1 from the exponential decay of S_z(E=2t, t), traced stochastically over
Nr=20 random-phase vectors (rp_spin.rp_spin_vs_time). The stochastic-trace error
scales as ~1/sqrt(Nr*N), so we use a larger lattice than the old exact tutorial.
A small exact-trace point cross-checks convergence.

Run:  LSQUANT_BIN=../../build OMP_NUM_THREADS=64 python3 rp_validate.py
"""
import os, time
import numpy as np
import make_ey, rp_spin as rp

THREADS = int(os.environ.get("OMP_NUM_THREADS", "64"))

def rates_vs_Dsf(L, Ds, Nr=20, M=128, NT=80, TM=260, seed=1):
    N = 2 * L * L
    out = []
    for D in Ds:
        lab = make_ey.build_ey(L, D, seed)
        t0 = time.time()
        t, s = rp.rp_spin_vs_time(lab, N, "SZ", "PZ", M, NT, TM, Nr=Nr, threads=THREADS)
        T1 = rp.fit_rate(t, s)
        out.append((D, T1, 1.0 / T1, time.time() - t0))
        print(f"  Dsf={D:.3f}  T1={T1:7.1f} fs  1/T1={1.0/T1:.4e}  ({time.time()-t0:.1f}s)")
    return out, (t, s)

if __name__ == "__main__":
    print(f"[rp_validate] threads={THREADS}")

    # --- (1) Elliott-Yafet rate law on a LARGER lattice, random-phase Nr=20 -----------
    Lbig, Ds = 50, [0.04, 0.06, 0.08, 0.10, 0.14]      # N = 2*50*50 = 5000
    print(f"\n(1) rate law -- random-phase Nr=20, L={Lbig} (N={2*Lbig*Lbig}):")
    t0 = time.time()
    data, _ = rates_vs_Dsf(Lbig, Ds, Nr=20)
    D = np.array([d for d, *_ in data]); rate = np.array([r for _, _, r, _ in data])
    p, logc = np.polyfit(np.log(D), np.log(rate), 1)
    print(f"  FIT  1/T1 ~ Dsf^p  ->  p = {p:.3f}   (Elliott-Yafet theory: p = 2)")
    print(f"  total wall for the sweep: {time.time()-t0:.1f}s")

    # --- (2) exact-trace cross-check at a small lattice (now cheap; serial guard) ------
    Lsm, Dsm = 8, 0.10                                 # N = 128, exact affordable
    Nsm = 2 * Lsm * Lsm
    lab = make_ey.build_ey(Lsm, Dsm, 1)
    M, NT, TM = 96, 60, 220
    t0 = time.time(); te, se = rp.exact_spin_vs_time(lab, Nsm, "SZ", "PZ", M, NT, TM); tex = time.time()-t0
    t0 = time.time(); tr, sr = rp.rp_spin_vs_time(lab, Nsm, "SZ", "PZ", M, NT, TM, Nr=20, threads=THREADS); trp = time.time()-t0
    T1e, T1r = rp.fit_rate(te, se), rp.fit_rate(tr, sr)
    print(f"\n(2) exact vs random-phase cross-check, L={Lsm} (N={Nsm}), Dsf={Dsm}:")
    print(f"  EXACT  T1={T1e:6.1f} fs  ({Nsm} states, {tex:.1f}s)")
    print(f"  RP(20) T1={T1r:6.1f} fs  (20 vectors, {trp:.1f}s)  ->  rel.diff {abs(T1r-T1e)/T1e*100:.1f}%")
    print(f"  speedup exact/rp = {tex/trp:.1f}x")

    # save the rate-law data for the figure
    np.savetxt("ey_ratelaw.dat", np.column_stack([D, rate]), header="Dsf  1/T1  (p=%.3f)" % p)
    print("\nwrote ey_ratelaw.dat")
