#!/usr/bin/env python3
"""Generate the Tutorial 03 figures from the output directory of `python -m lsquant run`.

Usage:
    python -m lsquant run examples/03_haldane_hall/run_haldane.json --out /tmp/tut03
    python examples/03_haldane_hall/plot.py /tmp/tut03
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

A_CELL = 3.0 * np.sqrt(3.0) / 2.0
E2_OVER_H = 1.0 / (2.0 * np.pi)
DIM_DEFAULT = 512


def bloch(
    k: np.ndarray,
    t1: float = -1.0,
    t2: float = 0.15,
    phi: float = np.pi / 2.0,
    mass: float = 0.0,
) -> np.ndarray:
    """Two-band Bloch Hamiltonian H(k) for the honeycomb Haldane model."""
    a1 = np.array([1.5, np.sqrt(3.0) / 2.0])
    a2 = np.array([1.5, -np.sqrt(3.0) / 2.0])
    nn = [
        np.array([0.5, np.sqrt(3.0) / 2.0]),
        np.array([0.5, -np.sqrt(3.0) / 2.0]),
        np.array([-1.0, 0.0]),
    ]
    nnn = [a1, a2 - a1, -a2]
    h0 = 2.0 * t2 * np.cos(phi) * sum(np.cos(k @ v) for v in nnn)
    hx = t1 * sum(np.cos(k @ v) for v in nn)
    hy = t1 * sum(np.sin(k @ v) for v in nn)
    hz = mass - 2.0 * t2 * np.sin(phi) * sum(np.sin(k @ v) for v in nnn)
    return h0 * np.eye(2) + np.array([[hz, hx - 1j * hy], [hx + 1j * hy, -hz]])


def chern_lower(num_grid: int = 24, **kw) -> float:
    """Chern number of the lower band by the Fukui-Hatsugai lattice method."""
    b1 = 2.0 * np.pi * np.array([1.0 / 3.0, 1.0 / np.sqrt(3.0)])
    b2 = 2.0 * np.pi * np.array([1.0 / 3.0, -1.0 / np.sqrt(3.0)])
    u = np.empty((num_grid, num_grid), object)
    for i in range(num_grid):
        for j in range(num_grid):
            k = (i / num_grid) * b1 + (j / num_grid) * b2
            _, v = np.linalg.eigh(bloch(k, **kw))
            u[i, j] = v[:, 0]
    flux = 0.0
    for i in range(num_grid):
        for j in range(num_grid):
            a = u[i, j]
            b = u[(i + 1) % num_grid, j]
            c = u[(i + 1) % num_grid, (j + 1) % num_grid]
            d = u[i, (j + 1) % num_grid]
            flux += np.angle(
                np.vdot(a, b) * np.vdot(b, c) / np.vdot(d, c) / np.vdot(a, d)
            )
    return float(flux / (2.0 * np.pi))


def chern_from_plateau(
    dim: int, energies: np.ndarray, sigma_xy: np.ndarray, gap_halfwidth: float = 0.5
) -> float:
    """Chern number from the gap plateau: C = (plateau / A) * 2 pi."""
    area = (dim / 2.0) * A_CELL
    plateau = float(np.mean(sigma_xy[np.abs(energies) < gap_halfwidth]))
    return float((plateau / area) * 2.0 * np.pi)


def _ensure_outputs(out_dir: str, label: str = "haldane") -> None:
    pattern = os.path.join(out_dir, f"KuboBastin_NonEqOpVX-VY{label}*.dat")
    if not glob.glob(pattern):
        run_json = os.path.join(_HERE, "run_haldane.json")
        execute_runfile(run_json, out_dir_override=out_dir, repo_root=_REPO_ROOT)


def figure_bands(out_dir: str, num_grid: int = 220) -> None:
    """Opening figure: Haldane density of states showing the bulk gap."""
    b1 = 2.0 * np.pi * np.array([1.0 / 3.0, 1.0 / np.sqrt(3.0)])
    b2 = 2.0 * np.pi * np.array([1.0 / 3.0, -1.0 / np.sqrt(3.0)])
    energies_list: list[float] = []
    for i in range(num_grid):
        for j in range(num_grid):
            k = (i + 0.5) / num_grid * b1 + (j + 0.5) / num_grid * b2
            energies_list.extend(np.linalg.eigvalsh(bloch(k)).tolist())
    es = np.array(energies_list)
    prl.use("web", aspect=0.45)
    fig, ax = plt.subplots()
    ax.hist(es, bins=140, color=prl.COLORS[1], alpha=0.85)
    ax.axvspan(
        float(es[es < 0].max()),
        float(es[es > 0].min()),
        color="0.82",
        label="gap (no states)",
    )
    ax.set_xlabel(r"energy $E$ (eV)")
    ax.set_ylabel(r"density of states (a.u.)")
    ax.legend(loc="upper right")
    prl.save(fig, os.path.join(out_dir, "fig_bands"))


def figure_haldane_hall(
    out_dir: str, label: str = "haldane", dim: int = DIM_DEFAULT
) -> None:
    """Two-panel Hall (Bastin, Sea, Surf) and longitudinal (Greenwood) figure."""
    d_xy = np.loadtxt(
        sorted(glob.glob(os.path.join(out_dir, f"KuboBastin_NonEqOpVX-VY{label}*.dat")))[-1]
    )
    d_sea = np.loadtxt(
        sorted(glob.glob(os.path.join(out_dir, f"KuboBastinSea_NonEqOpVX-VY{label}*.dat")))[-1]
    )
    d_surf = np.loadtxt(
        sorted(glob.glob(os.path.join(out_dir, f"KuboBastinSurf_NonEqOpVX-VY{label}*.dat")))[-1]
    )
    d_xx = np.loadtxt(
        sorted(glob.glob(os.path.join(out_dir, f"KuboGreenWood_NonEqOpVX-VX{label}*.dat")))[-1]
    )

    norm = (dim / 2.0) * A_CELL * E2_OVER_H
    chern = chern_from_plateau(dim, d_xy[:, 0], d_xy[:, 1])

    prl.use("web", height=5.8)
    fig, (top, bot) = plt.subplots(2, 1, sharex=True)

    curves = [
        ("Bastin (total)", d_xy),
        ("Fermi sea", d_sea),
        ("Fermi surface", d_surf),
    ]
    for i, (lab, arr) in enumerate(curves):
        top.plot(arr[:, 0], arr[:, 1] / norm, label=lab, **prl.trace(i))
    top.axhline(1.0, color="0.6", ls=":", lw=0.8)
    top.set_ylabel(r"$\sigma_{xy}$ ($e^2/h$)")
    prl.panel(top, "(a)", loc="upper right")
    top.legend(loc="upper left", title=r"$C=%+.3f$ (gap plateau)" % chern)

    bot.plot(d_xx[:, 0], d_xx[:, 1] / norm, **prl.trace(3))
    bot.set_xlabel(r"energy $E$ (eV)")
    bot.set_ylabel(r"$\sigma_{xx}$ ($e^2/h$)")
    prl.panel(bot, "(b)")

    prl.save(fig, os.path.join(out_dir, "fig_haldane_hall"))
    print(f"Chern from gap plateau: C = {chern:+.3f}")


def figure_phase_diagram(
    out_dir: str, t2: float = 0.15, phi: float = np.pi / 2.0
) -> None:
    """Topological phase transition as a function of sublattice mass M."""
    masses = np.linspace(0.0, 1.6, 17)
    cherns = [round(chern_lower(mass=float(m), t2=t2, phi=phi)) for m in masses]
    mc = 3.0 * np.sqrt(3.0) * t2 * abs(np.sin(phi))
    prl.use("web", aspect=0.42)
    fig, ax = plt.subplots()
    ax.plot(masses, cherns, lw=1.2, **prl.trace(0, marker=True))
    ax.axvline(mc, color="0.5", ls=":", lw=1.0, label=r"gap closes, $M_c=%.3f$" % mc)
    ax.set_xlabel(r"sublattice mass $M$ (eV)")
    ax.set_ylabel(r"Chern number $C$")
    ax.legend(loc="upper right")
    prl.save(fig, os.path.join(out_dir, "fig_phase_diagram"))


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
    run_dir, out_dir = _parse_dirs(argv, "run_haldane")
    os.makedirs(run_dir, exist_ok=True)
    os.makedirs(out_dir, exist_ok=True)
    _ensure_outputs(run_dir)
    figure_bands(run_dir)
    figure_haldane_hall(run_dir)
    figure_phase_diagram(run_dir)
    if os.path.abspath(run_dir) != os.path.abspath(out_dir):
        for fn in ("fig_bands.png", "fig_haldane_hall.png", "fig_phase_diagram.png"):
            src = os.path.join(run_dir, fn)
            if os.path.isfile(src):
                shutil.copy2(src, os.path.join(out_dir, fn))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
