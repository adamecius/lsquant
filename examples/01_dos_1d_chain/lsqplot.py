#!/usr/bin/env python3
"""lsqplot - a thin Python layer over the compiled LinQT kernels.

The heavy work stays in C++. Python orchestrates runs, reads the light moment
files, and plots. It drives the modern unified `lsquant` entry points:

    lsquant compute --config <json>   computes the Chebyshev moments mu_m (mode
                                      spectral; the threaded Hamiltonian recursion)
    lsquant reconstruct <m> dos <eta> performs the reconstruction sum
                                      rho(E) = sum_m g_m mu_m T_m(E)

(Both are byte-identical to the older inline_compute-kpm-spectralOp /
spectralFunctionFromChebmom drivers -- guarded by the config_equivalence test --
but the `lsquant compute|reconstruct` verbs are the supported surface.)

The reconstruction sum runs over a fine energy grid and every moment, so it is
left in C++. Python reads its small text output and draws the figure. Reading
the moments themselves (a few hundred numbers) is light, so that part is done
in Python directly for the moments plot.

The programs are located through the LSQUANT_BIN environment variable, then the
PATH. Run this from inside the example folder, after make_chain.py has built
the operators.
"""

import os
import sys
import glob
import shutil
import tempfile
import subprocess
import json

import numpy as np

# Shared publication (APS/PRL) plot style: examples/prl_style.py.
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import prl_style as prl


# ----------------------------------------------------------------------
# locate and call the compiled kernels
# ----------------------------------------------------------------------
def _bin(name):
    """Resolve a LinQT executable from LSQUANT_BIN or the PATH.

    The build-dir target names carry an 'inline_' prefix for the *FromChebmom drivers
    (e.g. inline_spectralFunctionFromChebmom); an install may expose them without it. Try the
    given name and its 'inline_'-prefixed variant, so the example works against either.
    """
    candidates = [name]
    if not name.startswith("inline_"):
        candidates.append("inline_" + name)
    root = os.environ.get("LSQUANT_BIN", "")
    for nm in candidates:
        if root:
            cand = os.path.join(root, nm)
            if os.path.exists(cand):
                return cand
        found = shutil.which(nm)
        if found:
            return found
    raise FileNotFoundError(
        f"could not find '{name}' (tried {candidates}). Build LinQT and set LSQUANT_BIN "
        f"to the directory that holds the executables."
    )


def _run_in(workdir, argv):
    subprocess.run(argv, cwd=workdir, check=True,
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)


def _stage(workdir, label):
    """Copy the operator and BOUNDS a driver expects into a clean workdir."""
    os.makedirs(os.path.join(workdir, "operators"), exist_ok=True)
    shutil.copy(f"operators/{label}.HAM.CSR",
                os.path.join(workdir, "operators", f"{label}.HAM.CSR"))
    shutil.copy("BOUNDS", os.path.join(workdir, "BOUNDS"))


def compute_moments(label, M, statefile):
    """Run the Chebyshev recursion. Return the path to the .chebmom1D file.

    Operator label 1 selects the identity, so the spectral function is the
    density of states.
    """
    with tempfile.TemporaryDirectory() as wd:
        _stage(wd, label)
        shutil.copy(statefile, os.path.join(wd, os.path.basename(statefile)))
        # Modern entry point: `lsquant compute` (mode spectral, operator 1 = identity).
        # Byte-identical to the legacy `inline_compute-kpm-spectralOp` (guarded by the
        # config_equivalence test); the "state" field is the trace state-vector file.
        with open(os.path.join(wd, "run_dos.json"), "w") as cfg:
            json.dump({"mode": "spectral", "label": label, "operator": "1",
                       "num_moments": M, "state": os.path.basename(statefile)}, cfg)
        _run_in(wd, [_bin("lsquant"), "compute", "--config", "run_dos.json"])
        produced = sorted(glob.glob(os.path.join(wd, "SpectralOp1*.chebmom1D")))
        if not produced:
            raise RuntimeError(f"no moments produced for {label}")
        out = f"moments_{label}_M{M}.chebmom1D"
        shutil.copy(produced[-1], out)
    return out


def read_moments(path):
    """Read a .chebmom1D file. Light: a header line plus M complex moments."""
    with open(path) as f:
        size, band_width, band_center = f.readline().split()[:3]
        M = int(f.readline())
        mu = np.array([complex(float(a), float(b))
                       for a, b in (f.readline().split() for _ in range(M))])
    return {"size": int(size),
            "band_width": float(band_width),
            "band_center": float(band_center),
            "M": M,
            "mu": mu}


def reconstruct_dos(chebmom_path, broadening_meV):
    """Run the reconstruction sum in C++. Return (energy_eV, dos_per_site).

    The kernel truncates to M_eff = ceil(pi * BandWidth / (2 * broadening)), so the broadening
    is set by the number of moments it keeps. spectralFunctionFromChebmom returns the EXTENSIVE
    (total) DOS = N * rho(E); we divide by the system size N (read from the .chebmom1D header) to
    get the intensive per-site DOS that matches rho(E) = 1/(pi sqrt(4 - E^2)) and overlaps across
    sizes.
    """
    cheb_abs = os.path.abspath(chebmom_path)
    size = read_moments(chebmom_path)["size"]
    with tempfile.TemporaryDirectory() as wd:
        shutil.copy(cheb_abs, os.path.join(wd, os.path.basename(cheb_abs)))
        # Modern entry point: `lsquant reconstruct <moments> dos <broadening>`.
        # Byte-identical to the legacy `spectralFunctionFromChebmom` (same kernel).
        _run_in(wd, [_bin("lsquant"), "reconstruct",
                     os.path.basename(cheb_abs), "dos", str(broadening_meV)])
        produced = sorted(glob.glob(os.path.join(wd, "mean*JACKSON.dat")))
        if not produced:
            raise RuntimeError("reconstruction produced no output")
        data = np.loadtxt(produced[-1])
    return data[:, 0], data[:, 1] / size


def m_eff(broadening_meV, band_width_eV):
    """The number of moments the Jackson kernel keeps for a given broadening."""
    return int(np.ceil(np.pi * band_width_eV / (2.0 * broadening_meV / 1000.0)))


# ----------------------------------------------------------------------
# the demo figures
# ----------------------------------------------------------------------
def figure_moments(sizes, M):
    """Overlay the moments of three sizes. They sit on delta_{m,0}."""
    import matplotlib.pyplot as plt
    prl.use("web", aspect=0.42)
    fig, ax = plt.subplots()
    for i, N in enumerate(sizes):
        moments = compute_moments(f"chain1d_N{N}", M, f"exact_chain1d_N{N}")
        d = read_moments(moments)
        ax.plot(np.arange(d["M"]), d["mu"].real, lw=1.0, markersize=3,
                label=r"$N=%d$" % N, **prl.trace(i, marker=True))
    ax.set_xlabel(r"moment index $m$")
    ax.set_ylabel(r"$\mu_m$")
    ax.legend(loc="upper right")
    prl.save(fig, "fig_moments")


def figure_dos(sizes, smallest, M, fixed_broadening_meV, broadenings_meV):
    """Top: DOS vs system size at fixed broadening.
       Bottom: DOS vs broadening at the smallest size."""
    import matplotlib.pyplot as plt
    prl.use("web", height=5.8)
    fig, (top, bot) = plt.subplots(2, 1, sharex=True)

    # (a) vary system size, broadening fixed
    for i, N in enumerate(sizes):
        moments = compute_moments(f"chain1d_N{N}", M, f"exact_chain1d_N{N}")
        E, dos = reconstruct_dos(moments, fixed_broadening_meV)
        top.plot(E, dos, label=r"$N=%d$" % N, **prl.trace(i))
    top.set_ylabel(r"$\rho(E)\;(1/\mathrm{eV})$")
    prl.panel(top, "(a)")
    top.legend(loc="upper center", ncol=3, title=r"$\eta=%d\,$meV fixed" % fixed_broadening_meV)

    # (b) vary broadening, smallest size
    moments = compute_moments(f"chain1d_N{smallest}", M,
                              f"exact_chain1d_N{smallest}")
    d = read_moments(moments)
    for i, eta in enumerate(broadenings_meV):
        E, dos = reconstruct_dos(moments, eta)
        bot.plot(E, dos, label=r"$\eta=%d\,$meV $\to M=%d$" % (eta, m_eff(eta, d['band_width'])),
                 **prl.trace(i))
    bot.set_xlabel(r"energy $E$ (eV)")
    bot.set_ylabel(r"$\rho(E)\;(1/\mathrm{eV})$")
    prl.panel(bot, "(b)")
    bot.legend(loc="upper center", title=r"$N=%d$ fixed" % smallest)

    # dense oscillatory curves -> cap dpi so the committed PNG stays < 200 KB
    prl.save(fig, "fig_dos", dpi=110)


def demo():
    sizes = [256, 1024, 4096]
    M = 512
    figure_moments(sizes, M=512)   # reach m=256 so the smallest chain's revival is in view
    # Exact-trace moments encode a finite chain only through the aliasing revival at m = N,
    # so the discrete spectrum of the small chain becomes visible once the kernel keeps
    # M_eff >= N moments. With N_small = 256 that needs a broadening below ~24 meV
    # (M_eff = ceil(pi*BandWidth/(2*eta))). The top panel fixes such a broadening (20 meV,
    # M_eff ~ 315): the N=256 curve carries discrete structure while the larger chains are
    # smooth bulk. The bottom panel sweeps the broadening across that threshold.
    figure_dos(sizes, smallest=256, M=M,
               fixed_broadening_meV=20,
               broadenings_meV=[80, 24, 10])


def figure_levels(N=64):
    """Opening figure: the discrete chain spectrum under the smooth bulk DOS."""
    import matplotlib.pyplot as plt
    prl.use("web", aspect=0.45)
    Ek = 2 * np.cos(2 * np.pi * np.arange(N) / N)
    E = np.linspace(-1.999, 1.999, 800); rho = 1 / (np.pi * np.sqrt(4 - E ** 2))
    fig, ax = plt.subplots()
    ax.vlines(Ek, 0, 0.6, color="0.55", lw=0.8, label=r"$%d$ chain levels" % N)
    ax.plot(E, rho, color=prl.COLORS[1], lw=1.8,
            label=r"bulk DOS $\;1/\pi\sqrt{4-E^2}$")
    ax.set_xlabel(r"energy $E$ (eV)"); ax.set_ylabel(r"$\rho(E)\;(1/\mathrm{eV})$")
    ax.set_xlim(-2.2, 2.2); ax.set_ylim(0, 0.65); ax.legend(loc="upper center")
    prl.save(fig, "fig_levels")


if __name__ == "__main__":
    figure_levels()
    demo()
