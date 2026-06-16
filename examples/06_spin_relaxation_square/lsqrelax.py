"""
lsqrelax.py -- spin relaxation on the 2D square lattice (illustrative).

Part A (Elliott-Yafet): random onsite spin-flip relaxes S_z. We extract the
energy-resolved S_z(E,t) at E=2t and show (i) it decays exponentially and (ii)
the rate obeys 1/T1 ~ Dsf^2 (the golden rule scaling). Writes fig_ey.png.

Part B (Dyakonov-Perel): uniform Rashba SOC + Anderson disorder. Clean (W=0) gives
precession; disorder drives motional-narrowing relaxation whose lifetime GROWS
with scattering (opposite to EY), and out-of-plane spin relaxes faster than
in-plane. Writes fig_dp.png.

Energy-resolved spin: S_a(E,t) = <phi(t)|S_a delta(E-H)|phi(t)> / <phi|delta(E-H)|phi>.
The denominator is conserved (energy), so we normalise by its time-average. Exact
trace -> deterministic. The binaries are found on PATH or under $LSQUANT_BIN.
Toy-scale and illustrative, not converged production numbers.
"""
import os
import glob
import subprocess
from shutil import which

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import make_ey
import make_rashba

E = 2.0
HBAR = 0.6582119624


def _bin(name):
    b = os.environ.get("LSQUANT_BIN")
    return os.path.join(b, name) if b and os.path.exists(os.path.join(b, name)) else (which(name) or name)


def _run(argv):
    subprocess.run(argv, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)


def spin_vs_time(label, N, op, proj, M, ntime, tmax, eta_mev=140):
    """Return (t, S(E,t)) normalised to S(0)=1, by exact trace at energy E."""
    state = "exact"
    if not os.path.exists(state):
        with open(state, "w") as f:
            f.write("local\n%d\n%s\n" % (N, "\n".join(map(str, range(N)))))
    _run([_bin("inline_compute-kpm-TimeEvProjetedOp"), label, op, proj,
          str(M), str(ntime), str(tmax), state])
    _run([_bin("inline_compute-kpm-TimeEvProjetedOp"), label, "1", proj,
          str(M), str(ntime), str(tmax), state])
    num_mom = sorted(glob.glob("TimeProj%s-%s%s*chebmomTD" % (op, proj, label)))[-1]
    den_mom = sorted(glob.glob("TimeProj1-%s%s*chebmomTD" % (proj, label)))[-1]
    _run([_bin("inline_timeCorrelationsFromChebmom"), num_mom, str(eta_mev), str(E)])
    num_dat = sorted(glob.glob("mean*%s-%s%s*EF%f*dat" % (op, proj, label, E)))[-1]
    _run([_bin("inline_timeCorrelationsFromChebmom"), den_mom, str(eta_mev), str(E)])
    den_dat = sorted(glob.glob("mean*1-%s%s*EF%f*dat" % (proj, label, E)))[-1]
    t, num = np.loadtxt(num_dat, unpack=True)
    _, den = np.loadtxt(den_dat, unpack=True)
    good = t > 0
    s = num[good] / den[good].mean()
    return t[good], s / s[0]


def fit_rate(t, s, lo=0.2, hi=0.85):
    m = (s > lo) & (s < hi)
    if m.sum() < 3:
        return np.nan
    slope = np.polyfit(t[m], np.log(s[m]), 1)[0]
    return -1.0 / slope


def part_A():
    L, seed = 20, 1
    Ds = [0.05, 0.07, 0.10]
    M, NT, TM = 128, 80, 260
    rates, curves = [], {}
    for D in Ds:
        lab = make_ey.build_ey(L, D, seed)
        t, s = spin_vs_time(lab, 2 * L * L, "SZ", "PZ", M, NT, TM)
        curves[D] = (t, s)
        T1 = fit_rate(t, s)
        rates.append(1.0 / T1)
        print("EY  Dsf=%.2f  T1=%.1f fs  1/T1=%.4f" % (D, T1, 1.0 / T1))

    fig, (a, b) = plt.subplots(1, 2, figsize=(9.2, 4.0))
    t, s = curves[0.10]
    T1 = fit_rate(t, s)
    a.plot(t, s, "o", ms=3, color="#1f77b4", label=r"$S_z(E{=}2t,\,t)$")
    a.plot(t, np.exp(-t / T1), "-", color="#1f77b4", alpha=0.6,
           label=r"exponential, $T_1\approx%.0f$ fs" % T1)
    a.axhline(0, color="0.85", lw=0.8)
    a.set_xlabel("time $t$ (fs)"); a.set_ylabel(r"$S_z(E,t)$"); a.legend(fontsize=9)
    a.set_title(r"Elliott-Yafet: $S_z$ decays exponentially ($\Delta_{\rm sf}=0.1$)")
    cols = plt.cm.viridis(np.linspace(0.15, 0.8, len(Ds)))
    for D, c in zip(Ds, cols):
        tt, ss = curves[D]
        b.plot(tt, ss, color=c, lw=1.6, label=r"$\Delta_{\rm sf}=%.2f$" % D)
    b.axhline(0, color="0.85", lw=0.8)
    b.set_xlabel("time $t$ (fs)"); b.set_ylabel(r"$S_z(E,t)$"); b.legend(fontsize=9)
    b.set_title(r"stronger spin-flip $\Rightarrow$ faster relaxation")
    fig.tight_layout(); fig.savefig("fig_ey.png", dpi=150); print("wrote fig_ey.png")


def part_B():
    L, lam, seed = 20, 0.2, 1
    M, NT, TM = 160, 120, 45      # short window: late-time Chebyshev evolution drifts
    curves = {}
    for W in (0.0, 6.0):
        lab = make_rashba.build_rashba(L, lam, W, seed)
        curves[W] = spin_vs_time(lab, 2 * L * L, "SZ", "PZ", M, NT, TM)
    labW = make_rashba.build_rashba(L, lam, 6.0, seed)
    tx, sx = spin_vs_time(labW, 2 * L * L, "SX", "PX", M, NT, TM)
    print("DP  W=0 (precession) and W=6 (relaxation) + in/out-of-plane computed")

    fig, (a, b) = plt.subplots(1, 2, figsize=(9.2, 4.0))
    a.plot(*curves[0.0], color="#9467bd", lw=1.5, label=r"$W=0$ (precession)")
    a.plot(*curves[6.0], color="#2ca02c", lw=1.6, label=r"$W=6$ (relaxes)")
    a.axhline(0, color="0.85", lw=0.8); a.set_ylim(-1.25, 1.25)
    a.set_xlabel("time $t$ (fs)"); a.set_ylabel(r"$S_z(E,t)$"); a.legend(fontsize=9)
    a.set_title("Dyakonov-Perel: clean precesses, disorder relaxes")
    b.plot(*curves[6.0], color="#2ca02c", lw=1.6, label=r"$S_z$ (out of plane)")
    b.plot(tx, sx, color="#d62728", lw=1.6, label=r"$S_x$ (in plane)")
    b.axhline(0, color="0.85", lw=0.8)
    b.set_xlabel("time $t$ (fs)"); b.set_ylabel(r"$S(E,t)$"); b.legend(fontsize=9)
    b.set_title("in-plane spin outlives out-of-plane")
    fig.tight_layout(); fig.savefig("fig_dp.png", dpi=150); print("wrote fig_dp.png")


if __name__ == "__main__":
    part_A()
    part_B()
