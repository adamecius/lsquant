#!/usr/bin/env python3
"""Render the four evaluation figures in PRL style (examples/prl_style.py):

    fig_web        a webpage/GitHub figure  (6.8 in, larger type)
    fig_single     a single-column article figure (3.375 in)
    fig_double     a double-column article figure (7.0 in)
    fig_grid       a 2x2 multi-panel, double-column width, shared x, panels (a)-(d)

Data is KPM-flavoured but analytic, so the demo is fast and self-contained --
the point is to judge the STYLE. Outputs go to examples/_style_preview/
(PDF = vector article copy, PNG = on-screen preview).
"""

import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import prl_style as prl

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "_style_preview")
os.makedirs(OUT, exist_ok=True)


# --- analytic KPM-flavoured data ---------------------------------------------
def dos_1d(E, eta):
    """Broadened 1D tight-binding DOS rho(E)=1/(pi sqrt(4-E^2)) (Lorentzian eta)."""
    out = np.zeros_like(E)
    Ep = np.linspace(-1.999, 1.999, 4000)
    rho = 1.0 / (np.pi * np.sqrt(4 - Ep ** 2))
    for i, e in enumerate(E):
        out[i] = np.trapezoid(rho * eta / np.pi / ((e - Ep) ** 2 + eta ** 2), Ep)
    return out


def sigma_xx(E, broad=0.15):
    """A longitudinal-conductivity-like bell, suppressed at band edges."""
    return np.exp(-(E ** 2) / (2 * 0.8 ** 2)) * (1 - 0.3 * np.cos(2 * E)) / (1 + broad)


def sigma_xy(E, gap=0.5, smooth=0.12):
    """A Hall plateau: sigma_xy -> 1 e^2/h for |E| < gap, rounded by `smooth`.

    Difference of two logistic steps -> a plateau of height 1 and width 2*gap."""
    lo = 1.0 / (1.0 + np.exp(-(E + gap) / smooth))
    hi = 1.0 / (1.0 + np.exp(-(E - gap) / smooth))
    return lo - hi


def spin_decay(t, tau, omega=0.0):
    return np.cos(omega * t) * np.exp(-t / tau)


# --- figures ------------------------------------------------------------------
def fig_web():
    prl.use("web", aspect=0.5)
    E = np.linspace(-2.4, 2.4, 600)
    fig, ax = plt.subplots()
    for i, eta in enumerate([0.05, 0.10, 0.20]):
        ax.plot(E, dos_1d(E, eta), label=r"$\eta=%.2f\,t$" % eta, **prl.trace(i))
    ax.set_xlabel(r"$E/t$")
    ax.set_ylabel(r"$\rho(E)\,(1/t)$")
    ax.set_xlim(-2.4, 2.4)
    ax.legend(loc="upper center", ncol=3)
    prl.save(fig, os.path.join(OUT, "fig_web"))
    plt.close(fig)


def fig_single():
    prl.use("single")
    E = np.linspace(-2.4, 2.4, 500)
    fig, ax = plt.subplots()
    for i, eta in enumerate([0.05, 0.15]):       # 2 traces: single column stays legible
        ax.plot(E, dos_1d(E, eta), label=r"$\eta=%.2f\,t$" % eta, **prl.trace(i))
    ax.set_xlabel(r"$E/t$")
    ax.set_ylabel(r"$\rho(E)\,(1/t)$")
    ax.set_xlim(-2.4, 2.4)
    ax.legend(loc="upper right")
    prl.save(fig, os.path.join(OUT, "fig_single"))
    plt.close(fig)


def fig_double():
    prl.use("double", aspect=0.42)
    E = np.linspace(-2.0, 2.0, 600)
    fig, ax = plt.subplots()
    for i, g in enumerate([0.3, 0.5, 0.7, 0.9]):
        ax.plot(E, sigma_xy(E, gap=g), label=r"$\Delta=%.1f\,t$" % g, **prl.trace(i))
    ax.axhline(0.0, color="0.6", lw=0.6, zorder=0)
    ax.set_xlabel(r"$E_F/t$")
    ax.set_ylabel(r"$\sigma_{xy}\,(e^2/h)$")
    ax.set_xlim(-2.0, 2.0)
    ax.legend(loc="upper left", ncol=2)
    prl.save(fig, os.path.join(OUT, "fig_double"))
    plt.close(fig)


def fig_grid():
    # Coherent multi-panel: the LEFT column is energy-resolved and shares its x
    # (the "shared axes where possible" rule); the RIGHT column is a separate
    # observable per panel, so each keeps its own x. Panels tagged (a)-(d).
    prl.use("double", height=4.6)
    E = np.linspace(-2.2, 2.2, 500)
    t = np.linspace(0, 12, 240)
    fig, ax = plt.subplots(2, 2)

    # left column shares the energy axis
    ax[0, 0].sharex(ax[1, 0])
    ax[0, 0].tick_params(labelbottom=False)

    # (a) DOS vs E
    for i, eta in enumerate([0.05, 0.15]):
        ax[0, 0].plot(E, dos_1d(E, eta), label=r"$\eta=%.2f$" % eta, **prl.trace(i))
    ax[0, 0].set_ylabel(r"$\rho(E)\,(1/t)$"); prl.panel(ax[0, 0], "(a)")
    ax[0, 0].legend(loc="upper center", ncol=2)

    # (c) sigma_xx vs E (shares x with a)
    ax[1, 0].plot(E, sigma_xx(E), **prl.trace(2))
    ax[1, 0].set_xlabel(r"$E_F/t$"); ax[1, 0].set_ylabel(r"$\sigma_{xx}$ (arb.)")
    ax[1, 0].set_ylim(0, 0.72); prl.panel(ax[1, 0], "(c)")

    # (b) S(t) spin relaxation, several rates (time axis)
    for i, tau in enumerate([8.0, 3.5, 1.6]):
        ax[0, 1].plot(t, spin_decay(t, tau), label=r"$\tau=%.1f$" % tau, **prl.trace(i))
    ax[0, 1].axhline(0, color="0.6", lw=0.6, zorder=0)
    ax[0, 1].set_xlabel(r"$t\,(\hbar/t)$"); ax[0, 1].set_ylabel(r"$S_z(t)/S_z(0)$")
    ax[0, 1].set_xlim(0, 12); ax[0, 1].set_ylim(-0.03, 1.05)
    prl.panel(ax[0, 1], "(b)", loc="lower left")
    ax[0, 1].legend(loc="upper right")

    # (d) sigma_xy plateau vs E
    for i, g in enumerate([0.4, 0.7]):
        ax[1, 1].plot(E, sigma_xy(E, gap=g), label=r"$\Delta=%.1f$" % g, **prl.trace(i))
    ax[1, 1].axhline(0, color="0.6", lw=0.6, zorder=0)
    ax[1, 1].set_xlabel(r"$E_F/t$"); ax[1, 1].set_ylabel(r"$\sigma_{xy}\,(e^2/h)$")
    ax[1, 1].set_xlim(-2.2, 2.2); ax[1, 1].set_ylim(-0.08, 1.18)
    prl.panel(ax[1, 1], "(d)", loc="upper right")
    ax[1, 1].legend(loc="lower left")

    prl.save(fig, os.path.join(OUT, "fig_grid"))
    plt.close(fig)


if __name__ == "__main__":
    fig_web()
    fig_single()
    fig_double()
    fig_grid()
    print("\npreview figures in:", OUT)
