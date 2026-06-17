"""
lsqspin.py -- Larmor spin precession on the magnetic chain.

For a spin prepared along +x (via the PX = 1/2(I + S_x) projector), this runs the
time-evolved projected-spin driver for S_x, S_y and S_z, reads the m=0 Chebyshev
moment of each (the energy-integrated spin expectation), normalises by <S_x(0)>,
and plots the precession against the analytic Larmor solution. Writes
fig_precession.png. The driver is found on PATH or under $LSQUANT_BIN.
"""
import os
import sys
import glob
import subprocess
from shutil import which

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# Shared publication (APS/PRL) plot style: examples/prl_style.py.
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import prl_style as prl

LABEL = "chain1d_mag"
M, NTIME, TMAX = 32, 96, 95.0     # moments, time steps, max time (fs)
HBAR = 0.6582119624               # eV*fs (LSQUANT physical units)
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
    cx, cy, cz = prl.COLORS[1], prl.COLORS[2], prl.COLORS[3]
    prl.use("web", aspect=0.5)
    fig, ax = plt.subplots()
    # analytic Larmor solution (thin reference lines)
    ax.plot(tt, np.cos(w * tt), "-", color=cx, lw=1.0, alpha=0.45)
    ax.plot(tt, -np.sin(w * tt), "--", color=cy, lw=1.0, alpha=0.45)
    # KPM points (redundant: distinct colour AND marker)
    ax.plot(t, sx, "o", ms=4, color=cx, label=r"$\langle S_x\rangle$ (KPM)")
    ax.plot(t, sy, "s", ms=4, color=cy, label=r"$\langle S_y\rangle$ (KPM)")
    ax.plot(t, sz, "^", ms=4, color=cz, label=r"$\langle S_z\rangle$ (KPM)")
    ax.axhline(0.0, color="0.7", lw=0.6, zorder=0)
    ax.set_xlabel(r"time $t$ (fs)")
    ax.set_ylabel(r"normalised spin $\langle S_\alpha(t)\rangle/\langle S_x(0)\rangle$")
    ax.legend(loc="upper right", ncol=3)
    ax.set_ylim(-1.25, 1.6)
    prl.save(fig, "fig_precession")


def figure_zeeman(J=0.1):
    """Opening figure: the exchange-split spin bands."""
    k = np.linspace(-np.pi, np.pi, 500)
    prl.use("web", aspect=0.45)
    fig, ax = plt.subplots()
    ax.plot(k, -2 * np.cos(k) - J, color=prl.COLORS[1], ls="-", lw=1.8, label=r"spin $\uparrow$")
    ax.plot(k, -2 * np.cos(k) + J, color=prl.COLORS[2], ls="--", lw=1.8, label=r"spin $\downarrow$")
    ax.annotate(r"$2J_{\rm ex}$", xy=(0, -2), xytext=(0.8, -1.6),
                arrowprops=dict(arrowstyle="<->"))
    ax.set_xlabel(r"$k$"); ax.set_ylabel(r"energy (eV)"); ax.legend(loc="upper center")
    prl.save(fig, "fig_zeeman")


if __name__ == "__main__":
    figure_zeeman()
    main()
