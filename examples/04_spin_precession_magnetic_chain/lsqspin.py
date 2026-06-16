"""
lsqspin.py -- Larmor spin precession on the magnetic chain.

For a spin prepared along +x (via the PX = 1/2(I + S_x) projector), this runs the
time-evolved projected-spin driver for S_x, S_y and S_z, reads the m=0 Chebyshev
moment of each (the energy-integrated spin expectation), normalises by <S_x(0)>,
and plots the precession against the analytic Larmor solution. Writes
fig_precession.png. The driver is found on PATH or under $LSQUANT_BIN.
"""
import os
import glob
import subprocess
from shutil import which

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

LABEL = "chain1d_mag"
M, NTIME, TMAX = 32, 96, 95.0     # moments, time steps, max time (fs)
HBAR = 0.6582119624               # eV*fs (LinQT physical units)
J_EX = 0.1                        # eV  (on-site exchange; Zeeman splitting 2*J_ex)


def _bin(name):
    bindir = os.environ.get("LSQUANT_BIN")
    if bindir and os.path.exists(os.path.join(bindir, name)):
        return os.path.join(bindir, name)
    return which(name) or name


def _m0(path):
    # .chebmomTD: line 1 header, line 2 "M Nt", then real/imag pairs, m-major.
    # The first Nt entries are the m=0 row -> the spin time series (real part).
    with open(path) as f:
        f.readline()
        nt = int(f.readline().split()[1])
        return np.array([float(f.readline().split()[0]) for _ in range(nt)]), nt


def main():
    drv = _bin("inline_compute-kpm-TimeEvProjetedOp")
    series = {}
    nt = NTIME
    for op in ("SX", "SY", "SZ"):
        subprocess.run([drv, LABEL, op, "PX", str(M), str(NTIME), str(TMAX), "exact"],
                       check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        mom = sorted(glob.glob("TimeProj%s-PX*stateexact*chebmomTD" % op))[-1]
        series[op], nt = _m0(mom)

    dt = TMAX / (NTIME - 1)
    t = np.arange(nt) * dt
    D = series["SX"][0]                                   # <S_x(0)>
    sx, sy, sz = series["SX"] / D, series["SY"] / D, series["SZ"] / D

    w = 2 * J_EX / HBAR
    tt = np.linspace(0.0, t[-1], 400)
    fig, ax = plt.subplots(figsize=(6.6, 4.0))
    ax.plot(tt, np.cos(w * tt), "-", color="#1f77b4", lw=1.0, alpha=0.5)
    ax.plot(tt, -np.sin(w * tt), "-", color="#d62728", lw=1.0, alpha=0.5)
    ax.plot(t, sx, "o", ms=3.5, color="#1f77b4", label=r"$\langle S_x\rangle$ (KPM)")
    ax.plot(t, sy, "s", ms=3.5, color="#d62728", label=r"$\langle S_y\rangle$ (KPM)")
    ax.plot(t, sz, "^", ms=3.5, color="#2ca02c", label=r"$\langle S_z\rangle$ (KPM)")
    ax.axhline(0.0, color="0.85", lw=0.8)
    ax.set_xlabel("time $t$ (fs)")
    ax.set_ylabel("normalised spin expectation")
    ax.set_title(r"Larmor precession, magnetic chain ($\omega=2J_{\rm ex}/\hbar=%.3f$ rad/fs, period %.1f fs)"
                 % (w, 2 * np.pi / w))
    ax.legend(loc="upper right", ncol=3, fontsize=8.5)
    ax.set_ylim(-1.25, 1.55)
    fig.tight_layout()
    fig.savefig("fig_precession.png", dpi=150)
    print("wrote fig_precession.png")


def figure_zeeman(J=0.1):
    """Opening figure: the exchange-split spin bands."""
    k = np.linspace(-np.pi, np.pi, 500)
    fig, ax = plt.subplots(figsize=(6.4, 3.7))
    ax.plot(k, -2 * np.cos(k) - J, color="#1f77b4", lw=1.8, label=r"spin $\uparrow$")
    ax.plot(k, -2 * np.cos(k) + J, color="#d62728", lw=1.8, label=r"spin $\downarrow$")
    ax.annotate(r"$2J_{\rm ex}$", xy=(0, -2), xytext=(0.8, -1.6),
                arrowprops=dict(arrowstyle="<->"), fontsize=11)
    ax.set_xlabel("$k$"); ax.set_ylabel("energy (eV)"); ax.legend(fontsize=9)
    ax.set_title(r"The exchange field splits $\uparrow$ and $\downarrow$ by $2J_{\rm ex}$")
    fig.tight_layout(); fig.savefig("fig_zeeman.png", dpi=150); print("wrote fig_zeeman.png")


if __name__ == "__main__":
    figure_zeeman()
    main()
