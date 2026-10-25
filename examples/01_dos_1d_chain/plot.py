#!/usr/bin/env python3
"""Generate the Tutorial 01 figures from the output directory of `python -m lsquant run`.

Usage:
    python -m lsquant run examples/01_dos_1d_chain/run_dos.json --out /tmp/tut01
    python examples/01_dos_1d_chain/plot.py /tmp/tut01
"""

from __future__ import annotations

import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
_EX_DIR = os.path.dirname(_HERE)
_REPO_ROOT = os.path.dirname(_EX_DIR)
for _p in (
    _EX_DIR,
    os.path.join(_REPO_ROOT, "python"),
    os.path.join(_REPO_ROOT, "build", "native", "lsquant_core"),
    os.path.join(_REPO_ROOT, "build", "golden", "lsquant_core"),
    os.path.join(_REPO_ROOT, "build", "lsquant_core"),
):
    if os.path.isdir(_p) and _p not in sys.path:
        sys.path.insert(0, _p)

import prl_style as prl
from lsquant.runfile import execute_runfile


def read_moments(path: str) -> dict:
    """Read a .chebmom1D file header and M complex moments."""
    with open(path, encoding="utf-8") as f:
        size, band_width, band_center = f.readline().split()[:3]
        num_moments = int(f.readline())
        mu = np.array(
            [
                complex(float(a), float(b))
                for a, b in (f.readline().split() for _ in range(num_moments))
            ]
        )
    return {
        "size": int(size),
        "band_width": float(band_width),
        "band_center": float(band_center),
        "M": num_moments,
        "mu": mu,
    }


def m_eff(broadening_meV: float, band_width_eV: float) -> int:
    """Number of moments retained by the Jackson kernel for a given broadening."""
    return int(np.ceil(np.pi * band_width_eV / (2.0 * broadening_meV / 1000.0)))


def _ensure_outputs(out_dir: str, sizes: list[int], num_moments: int) -> None:
    expected = os.path.join(
        out_dir, f"moments_chain1d_N{sizes[0]}_M{num_moments}.chebmom1D"
    )
    if not os.path.isfile(expected):
        run_json = os.path.join(_HERE, "run_dos.json")
        execute_runfile(run_json, out_dir_override=out_dir, repo_root=_REPO_ROOT)


def figure_levels(out_dir: str, num_sites: int = 64) -> None:
    """Opening figure: discrete chain levels under the smooth bulk DOS."""
    prl.use("web", aspect=0.45)
    ek = 2.0 * np.cos(2.0 * np.pi * np.arange(num_sites) / num_sites)
    energies = np.linspace(-1.999, 1.999, 800)
    rho = 1.0 / (np.pi * np.sqrt(4.0 - energies**2))
    fig, ax = plt.subplots()
    ax.vlines(ek, 0, 0.6, color="0.55", lw=0.8, label=r"$%d$ chain levels" % num_sites)
    ax.plot(
        energies,
        rho,
        color=prl.COLORS[1],
        lw=1.8,
        label=r"bulk DOS $\;1/\pi\sqrt{4-E^2}$",
    )
    ax.set_xlabel(r"energy $E$ (eV)")
    ax.set_ylabel(r"$\rho(E)\;(1/\mathrm{eV})$")
    ax.set_xlim(-2.2, 2.2)
    ax.set_ylim(0, 0.65)
    ax.legend(loc="upper center")
    prl.save(fig, os.path.join(out_dir, "fig_levels"))


def figure_moments(out_dir: str, sizes: list[int], num_moments: int = 512) -> None:
    """Overlay the moments of three chain sizes showing the revival at m = N."""
    prl.use("web", aspect=0.42)
    fig, ax = plt.subplots()
    for i, n_sites in enumerate(sizes):
        mom_path = os.path.join(
            out_dir, f"moments_chain1d_N{n_sites}_M{num_moments}.chebmom1D"
        )
        d = read_moments(mom_path)
        ax.plot(
            np.arange(d["M"]),
            d["mu"].real,
            lw=1.0,
            markersize=3,
            label=r"$N=%d$" % n_sites,
            **prl.trace(i, marker=True),
        )
    ax.set_xlabel(r"moment index $m$")
    ax.set_ylabel(r"$\mu_m$")
    ax.legend(loc="upper right")
    prl.save(fig, os.path.join(out_dir, "fig_moments"))


def figure_dos(
    out_dir: str,
    sizes: list[int],
    smallest: int = 256,
    num_moments: int = 512,
    fixed_broadening_meV: int = 20,
    broadenings_meV: list[int] = [80, 24, 10],
) -> None:
    """Top: DOS vs system size at fixed broadening. Bottom: DOS vs broadening."""
    prl.use("web", height=5.8)
    fig, (top, bot) = plt.subplots(2, 1, sharex=True)

    for i, n_sites in enumerate(sizes):
        dat_path = os.path.join(
            out_dir,
            f"eta{fixed_broadening_meV}",
            f"meanmoments_chain1d_N{n_sites}_M{num_moments}JACKSON.dat",
        )
        data = np.loadtxt(dat_path)
        top.plot(
            data[:, 0],
            data[:, 1] / float(n_sites),
            label=r"$N=%d$" % n_sites,
            **prl.trace(i),
        )
    top.set_ylabel(r"$\rho(E)\;(1/\mathrm{eV})$")
    prl.panel(top, "(a)")
    top.legend(
        loc="upper center",
        ncol=3,
        title=r"$\eta=%d\,$meV fixed" % fixed_broadening_meV,
    )

    d_small = read_moments(
        os.path.join(out_dir, f"moments_chain1d_N{smallest}_M{num_moments}.chebmom1D")
    )
    for i, eta in enumerate(broadenings_meV):
        dat_path = os.path.join(
            out_dir,
            f"eta{eta}",
            f"meanmoments_chain1d_N{smallest}_M{num_moments}JACKSON.dat",
        )
        data = np.loadtxt(dat_path)
        bot.plot(
            data[:, 0],
            data[:, 1] / float(smallest),
            label=r"$\eta=%d\,$meV $\to M=%d$" % (eta, m_eff(eta, d_small["band_width"])),
            **prl.trace(i),
        )
    bot.set_xlabel(r"energy $E$ (eV)")
    bot.set_ylabel(r"$\rho(E)\;(1/\mathrm{eV})$")
    prl.panel(bot, "(b)")
    bot.legend(loc="upper center", title=r"$N=%d$ fixed" % smallest)
    prl.save(fig, os.path.join(out_dir, "fig_dos"), dpi=110)


def _parse_dirs(argv: list[str] | None, sub_name: str) -> tuple[str, str]:
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("out_dir", nargs="?", default=".")
    parser.add_argument("--run-dir", default=None)
    parser.add_argument("--out", default=None)
    ns = parser.parse_args(argv)
    run_dir = os.path.abspath(ns.run_dir if ns.run_dir is not None else ns.out_dir)
    out_dir = os.path.abspath(ns.out if ns.out is not None else run_dir)
    sub_cand = os.path.join(run_dir, sub_name)
    if os.path.isdir(sub_cand):
        run_dir = sub_cand
    return run_dir, out_dir


def main(argv: list[str] | None = None) -> int:
    import shutil
    run_dir, out_dir = _parse_dirs(argv, "run_dos")
    os.makedirs(run_dir, exist_ok=True)
    os.makedirs(out_dir, exist_ok=True)
    sizes = [256, 1024, 4096]
    num_moments = 512
    _ensure_outputs(run_dir, sizes, num_moments)
    figure_levels(run_dir)
    figure_moments(run_dir, sizes, num_moments=num_moments)
    figure_dos(
        run_dir,
        sizes,
        smallest=256,
        num_moments=num_moments,
        fixed_broadening_meV=20,
        broadenings_meV=[80, 24, 10],
    )
    if os.path.abspath(run_dir) != os.path.abspath(out_dir):
        for fn in ("fig_levels.png", "fig_moments.png", "fig_dos.png"):
            src = os.path.join(run_dir, fn)
            if os.path.isfile(src):
                shutil.copy2(src, os.path.join(out_dir, fn))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
