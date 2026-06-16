#!/usr/bin/env python3
"""lsqloc - localization study on the 1D Anderson chain, over compiled LinQT.

The heavy work stays in C++:

    inline_compute-kpm-MeanSquareDisplacement   the time-dependent Chebyshev
                                                moments of the spreading
    inline_timeCorrelationsFromChebmom          the reconstruction sum that
                                                returns MSD(E, t)

Python builds each disorder realisation, launches the two kernels, averages
MSD(E, t) over realisations, reads the localization length from the late-time
plateau, and plots. Reading and averaging the small reconstructed curves is
light, so that part is in Python; both sums stay compiled.

The clean chain spreads ballistically, MSD ~ t^2. Disorder localizes every
state, so MSD(E, t) rises and then settles on a plateau set by the localization
length xi(E). The relative size sqrt(MSD_plateau) follows the localization
length, so its energy profile and its 1/W^2 scaling can be read directly; the
overall prefactor is fixed once against the clean ballistic reference.

Run from inside the example folder, after LinQT is built and installed
(executables on PATH, or LSQUANT_BIN set to the build directory).
"""

import os
import sys
import glob
import shutil
import tempfile
import subprocess

import numpy as np

import make_disordered_chain as gen

# Shared publication (APS/PRL) plot style: examples/prl_style.py.
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import prl_style as prl


# ----------------------------------------------------------------------
# compiled kernels
# ----------------------------------------------------------------------
def _bin(name):
    root = os.environ.get("LSQUANT_BIN", "")
    if root and os.path.exists(os.path.join(root, name)):
        return os.path.join(root, name)
    found = shutil.which(name)
    if found:
        return found
    raise FileNotFoundError(
        f"could not find '{name}'. Build and install LinQT, or set LSQUANT_BIN.")


def msd_curve(N, W, energies, seed, M=256, num_times=128, tmax=600, broad=10):
    """One disorder realisation. Return {E: (t, MSD)} for each Fermi energy.

    A fresh random vector (KPM_SEED) traces the spreading, so each realisation
    averages over disorder and over the starting vector together.
    """
    out = {}
    with tempfile.TemporaryDirectory() as wd:
        label = gen.build_chain(N, W, seed, outdir=wd)
        env = dict(os.environ, KPM_SEED=str(seed))
        subprocess.run([_bin("inline_compute-kpm-MeanSquareDisplacement"),
                        label, "VX", str(M), str(num_times), str(tmax)],
                       cwd=wd, env=env, check=True,
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        mom = sorted(glob.glob(os.path.join(wd, f"Correlation*{label}*chebmomTD")))
        if not mom:
            raise RuntimeError(f"MSD produced no moments for {label}")
        mom = mom[-1]
        for E in energies:
            subprocess.run([_bin("inline_timeCorrelationsFromChebmom"),
                            os.path.basename(mom), str(broad), str(E)],
                           cwd=wd, check=True,
                           stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            dat = sorted(glob.glob(os.path.join(wd, "mean*EF*.dat")))[-1]
            d = np.loadtxt(dat)
            out[E] = (d[:, 0], d[:, 1])
            os.remove(dat)
    return out


def averaged_msd(N, W, energies, realisations, **kw):
    """Average MSD(E, t) over disorder realisations on a shared time grid."""
    acc, t_axis = {E: None for E in energies}, None
    for r in range(1, realisations + 1):
        curves = msd_curve(N, W, energies, seed=r, **kw)
        for E in energies:
            t, msd = curves[E]
            t_axis = t if t_axis is None else t_axis
            acc[E] = msd.copy() if acc[E] is None else acc[E] + msd
    for E in energies:
        acc[E] /= realisations
    return t_axis, acc


def xi_from_plateau(t, msd, frac=0.3):
    """Relative localization length from the late-time plateau of the MSD.

    xi_rel = sqrt(<MSD>) over the last `frac` of the time window. The constant
    that turns xi_rel into an absolute length is fixed once against the clean
    ballistic reference; energy profiles and W scaling use xi_rel directly.
    """
    cut = t >= (1.0 - frac) * t.max()
    return np.sqrt(np.mean(msd[cut]))


def xi_born(E, W):
    """Weak-disorder (Born) localization length for box disorder, t = 1."""
    return 24.0 * (4.0 - np.asarray(E, float) ** 2) / W ** 2


# ----------------------------------------------------------------------
# figures
# ----------------------------------------------------------------------
def figure_localization(sizes, W, energies, realisations, **kw):
    import matplotlib.pyplot as plt

    # xi(E, N) for every size and energy
    xi = {N: {} for N in sizes}
    for N in sizes:
        t, acc = averaged_msd(N, W, energies, realisations, **kw)
        for E in energies:
            xi[N][E] = xi_from_plateau(t, acc[E])

    # calibrate the single prefactor against Born at the largest size, mid band
    Nref, Eref = max(sizes), 1.0
    cal = xi_born(Eref, W) / xi[Nref][Eref]

    prl.use("web", height=5.8)
    fig, (top, bot) = plt.subplots(2, 1)

    # (a) xi(E) across the band, three sizes, with the Born curve
    egrid = np.linspace(-1.95, 1.95, 200)
    top.plot(egrid, xi_born(egrid, W), color="0.2", ls="--", lw=1.2, label="Born theory")
    for i, N in enumerate(sizes):
        E = sorted(energies)
        top.plot(E, [cal * xi[N][e] for e in E], lw=1.0,
                 label=r"$N=%d$" % N, **prl.trace(i + 1, marker=True))
    top.set_xlabel(r"energy $E$ (eV)")
    top.set_ylabel(r"$\xi(E)$ (sites)")
    prl.panel(top, "(a)")
    top.legend(loc="upper right", title=r"$W=%g$ eV" % W)

    # (b) scaling with size at each marked energy
    for i, E in enumerate(sorted(energies)):
        bot.plot(sizes, [cal * xi[N][E] for N in sizes], lw=1.0,
                 label=r"$E=%+.1f$ eV" % E, **prl.trace(i, marker=True))
        bot.axhline(xi_born(E, W), color="0.6", lw=0.6, ls=":")
    bot.set_xscale("log", base=2)
    bot.set_xlabel(r"system size $N$ (sites)")
    bot.set_ylabel(r"$\xi$ (sites)")
    prl.panel(bot, "(b)", loc="lower left")
    bot.legend(loc="upper center", ncol=2)

    prl.save(fig, "fig_localization")


def demo():
    W = 1.0
    energies = [0.0, 1.0, -1.0, 1.9]      # centre, intermediate +/-, near edge
    sizes = [256, 1024, 4096]
    # At W = 1 the interior localization length is xi ~ 70-100 sites, so the wavepacket needs a
    # long time to spread to xi and saturate: tmax = 600 (with num_times = 128) reaches the
    # plateau for every energy, and N up to 4096 holds the localized state (N >> 2 xi). The
    # original tmax = 40 stopped while the interior MSD was still rising, far below the plateau.
    figure_localization(sizes, W, energies, realisations=16,
                        M=256, num_times=128, tmax=600, broad=10)


def figure_anderson(N=500):
    """Opening figure: a band-centre eigenstate, extended (clean) vs localized (W=4)."""
    import matplotlib.pyplot as plt

    def state(W, seed):
        rng = np.random.default_rng(seed)
        H = np.diag(rng.uniform(-W / 2, W / 2, N)).astype(float)
        idx = np.arange(N); H[idx, (idx + 1) % N] = -1; H[(idx + 1) % N, idx] = -1
        w, v = np.linalg.eigh(H); return np.abs(v[:, np.argmin(np.abs(w))]) ** 2
    prl.use("web", aspect=0.45)
    fig, ax = plt.subplots()
    ax.plot(state(0.0, 3), color=prl.COLORS[1], lw=1.2, ls="-",
            label=r"clean: extended ($|\psi|^2\sim 1/N$)")
    ax.plot(state(4.0, 3), color=prl.COLORS[2], lw=1.2, ls="--", label=r"$W=4$: localized")
    ax.set_xlabel(r"site $i$"); ax.set_ylabel(r"$|\psi_i|^2$ (band-centre state)")
    ax.legend(loc="upper right")
    prl.save(fig, "fig_anderson")


if __name__ == "__main__":
    figure_anderson()
    demo()
