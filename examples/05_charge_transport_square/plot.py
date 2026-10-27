#!/usr/bin/env python3
"""Generate the Tutorial 05 figures from the output directory of `python -m lsquant run`.

Usage:
    python -m lsquant run examples/05_charge_transport_square/run_transport.json --out /tmp/tut05
    python examples/05_charge_transport_square/plot.py /tmp/tut05
"""

from __future__ import annotations

import glob
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

L_SIDE = 200
E_WORK = 2.0
W_LIST = [0.5, 1.0, 2.0, 4.0]


def fermi_velocity(energy: float = E_WORK) -> float:
    """v_F = <|grad_k E|> over the constant-energy contour of the clean band."""
    k = np.linspace(-np.pi, np.pi, 600, endpoint=False)
    kx, ky = np.meshgrid(k, k)
    ek = -2.0 * (np.cos(kx) + np.cos(ky))
    vmag = 2.0 * np.sqrt(np.sin(kx) ** 2 + np.sin(ky) ** 2)
    shell = np.abs(ek - energy) < 0.05
    return float(vmag[shell].mean())


def _ensure_outputs(out_dir: str) -> None:
    pattern = os.path.join(
        out_dir, f"meanSpectralOp1square_L{L_SIDE}_W{W_LIST[-1]:g}*JACKSON.dat"
    )
    if not glob.glob(pattern):
        run_json = os.path.join(_HERE, "run_transport.json")
        execute_runfile(run_json, out_dir_override=out_dir, repo_root=_REPO_ROOT)


def figure_fermi(out_dir: str, energy: float = E_WORK) -> None:
    """Opening figure: square-lattice band and E = 2t Fermi contour."""
    k = np.linspace(-np.pi, np.pi, 400)
    kx, ky = np.meshgrid(k, k)
    ek = -2.0 * (np.cos(kx) + np.cos(ky))
    prl.use("web", aspect=0.85)
    fig, ax = plt.subplots(constrained_layout=True)
    cf = ax.contourf(kx, ky, ek, levels=20, cmap="coolwarm", alpha=0.9)
    ax.contour(kx, ky, ek, levels=[energy], colors="k", linewidths=2.0)
    ax.set_xlabel(r"$k_x$")
    ax.set_ylabel(r"$k_y$")
    ax.set_aspect("equal")
    fig.colorbar(cf, ax=ax, label=r"$E_{\mathbf{k}}$ (eV)")
    out_path = os.path.join(out_dir, "fig_fermi.png")
    fig.savefig(out_path)
    plt.close(fig)
    print(f"wrote {out_path}")


def figure_diffusion_and_transport(out_dir: str) -> None:
    """Read MSD and DOS curves from out_dir and write fig_diffusion.png and fig_transport.png."""
    vf = fermi_velocity(E_WORK)
    norb = L_SIDE * L_SIDE
    curves = {}
    sigma = []
    ell = []
    for w_val in W_LIST:
        label = f"square_L{L_SIDE}_W{w_val:g}"
        msd_matches = sorted(
            glob.glob(os.path.join(out_dir, f"meanCorrelation*{label}*EF{E_WORK:.6f}*.dat"))
            + glob.glob(os.path.join(out_dir, f"meanMSD*{label}*EF{E_WORK:.6f}*.dat"))
        )
        msd_dat = msd_matches[-1]
        t, dx2 = np.loadtxt(msd_dat, unpack=True)
        good = t > 0
        t_good = t[good]
        d_good = dx2[good] / (2.0 * t_good)
        d_max = float(d_good.max())

        dos_dat = sorted(
            glob.glob(os.path.join(out_dir, f"meanSpectralOp1{label}*JACKSON.dat"))
        )[-1]
        en, dos = np.loadtxt(dos_dat, unpack=True)
        rho = float(np.interp(E_WORK, en, dos)) / float(norb)

        sigma.append(rho * d_max)
        ell.append(2.0 * d_max / vf)
        curves[w_val] = (t_good, d_good)
        print(
            "W=%.1f  D_max=%6.3f  rho(E)=%6.4f  sigma=%6.3f  ell=%6.2f"
            % (w_val, d_max, rho, sigma[-1], ell[-1])
        )

    colors = plt.cm.viridis(np.linspace(0.1, 0.85, len(W_LIST)))
    prl.use("web", height=5.8)
    fig, (ax_a, ax_b) = plt.subplots(2, 1, sharex=True)
    for (w_val, (t_good, d_good)), c in zip(curves.items(), colors):
        ax_a.plot(t_good, d_good, color=c, lw=1.6, label=r"$W=%g$" % w_val)
        ax_b.plot(t_good, d_good / d_good.max(), color=c, lw=1.6)
    ax_a.set_ylabel(r"$D(E,t)=\Delta X^2/2t$")
    prl.panel(ax_a, "(a)")
    ax_a.legend(title=r"disorder $W/t$", ncol=2, loc="upper right")
    ax_b.set_ylabel(r"$D(E,t)/D_{\max}$")
    ax_b.set_xlabel(r"time $t$ (fs)")
    ax_b.axhline(1.0, color="0.7", lw=0.6, zorder=0)
    prl.panel(ax_b, "(b)")
    prl.save(fig, os.path.join(out_dir, "fig_diffusion"))

    c_s, c_l = prl.COLORS[1], prl.COLORS[2]
    prl.use("web", aspect=0.52)
    fig, ax = plt.subplots()
    ax.semilogy(W_LIST, sigma, "o-", color=c_s, label=r"$\sigma=e^2\rho D_{\max}$")
    ax.set_xlabel(r"disorder $W/t$")
    ax.set_ylabel(r"conductivity $\sigma$ (arb.)", color=c_s)
    ax.tick_params(axis="y", labelcolor=c_s)
    ax2 = ax.twinx()
    ax2.semilogy(W_LIST, ell, "s--", color=c_l, label=r"$\ell=2D_{\max}/v_F$")
    ax2.axhline(1.0, color=c_l, lw=0.7, ls=":")
    ax2.set_ylabel(r"mean free path $\ell$ ($a$)", color=c_l)
    ax2.tick_params(axis="y", labelcolor=c_l)
    ax2.spines["top"].set_visible(False)
    prl.save(fig, os.path.join(out_dir, "fig_transport"))


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
    run_dir, out_dir = _parse_dirs(argv, "run_transport")
    os.makedirs(run_dir, exist_ok=True)
    os.makedirs(out_dir, exist_ok=True)
    _ensure_outputs(run_dir)
    figure_fermi(run_dir)
    figure_diffusion_and_transport(run_dir)
    if os.path.abspath(run_dir) != os.path.abspath(out_dir):
        for fn in ("fig_fermi.png", "fig_diffusion.png", "fig_transport.png"):
            src = os.path.join(run_dir, fn)
            if os.path.isfile(src):
                shutil.copy2(src, os.path.join(out_dir, fn))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
