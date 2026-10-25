#!/usr/bin/env python3
"""Generate the Tutorial 04 figures from the output directory of `python -m lsquant run`.

Usage:
    python -m lsquant run examples/04_spin_precession_magnetic_chain/run_spin.json --out /tmp/tut04
    python examples/04_spin_precession_magnetic_chain/plot.py /tmp/tut04
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

NTIME = 96
TMAX = 95.0
HBAR = 0.6582119624
J_EX = 0.1


def _m0(path: str) -> tuple[np.ndarray, int]:
    """Read the m=0 row (energy-integrated spin expectation) from a .chebmomTD file."""
    with open(path, encoding="utf-8") as f:
        f.readline()
        nt = int(f.readline().split()[1])
        return np.array([float(f.readline().split()[0]) for _ in range(nt)]), nt


def _ensure_outputs(out_dir: str) -> None:
    pattern = os.path.join(out_dir, "TimeProjSX-PX*stateexact*chebmomTD")
    if not glob.glob(pattern):
        run_json = os.path.join(_HERE, "run_spin.json")
        execute_runfile(run_json, out_dir_override=out_dir, repo_root=_REPO_ROOT)


def figure_zeeman(out_dir: str, j_ex: float = J_EX) -> None:
    """Opening figure: exchange-split spin bands."""
    k = np.linspace(-np.pi, np.pi, 500)
    prl.use("web", aspect=0.45)
    fig, ax = plt.subplots()
    ax.plot(k, -2.0 * np.cos(k) - j_ex, color=prl.COLORS[1], lw=1.6, label=r"$E_\uparrow(k)$")
    ax.plot(
        k,
        -2.0 * np.cos(k) + j_ex,
        color=prl.COLORS[2],
        lw=1.6,
        ls="--",
        label=r"$E_\downarrow(k)=E_\uparrow+2J_{\rm ex}$",
    )
    ax.set_xlabel(r"wavevector $ka$")
    ax.set_ylabel(r"energy $E$ (eV)")
    ax.legend(loc="upper right")
    prl.save(fig, os.path.join(out_dir, "fig_zeeman"))


def figure_precession(out_dir: str) -> None:
    """Larmor spin precession <S_x>, <S_y>, <S_z> vs analytic cos(wt), -sin(wt), 0."""
    series = {}
    nt = NTIME
    for op in ("SX", "SY", "SZ"):
        mom = sorted(
            glob.glob(os.path.join(out_dir, f"TimeProj{op}-PX*stateexact*chebmomTD"))
        )[-1]
        series[op], nt = _m0(mom)

    dt = TMAX / (NTIME - 1)
    t = np.arange(nt) * dt
    norm = series["SX"][0]
    sx, sy, sz = series["SX"] / norm, series["SY"] / norm, series["SZ"] / norm

    omega = 2.0 * J_EX / HBAR
    tt = np.linspace(0.0, float(t[-1]), 400)
    cx, cy, cz = prl.COLORS[1], prl.COLORS[2], prl.COLORS[3]
    prl.use("web", aspect=0.5)
    fig, ax = plt.subplots()
    ax.plot(tt, np.cos(omega * tt), "-", color=cx, lw=1.0, alpha=0.45)
    ax.plot(tt, -np.sin(omega * tt), "--", color=cy, lw=1.0, alpha=0.45)
    ax.plot(t, sx, "o", ms=4, color=cx, label=r"$\langle S_x\rangle$ (KPM)")
    ax.plot(t, sy, "s", ms=4, color=cy, label=r"$\langle S_y\rangle$ (KPM)")
    ax.plot(t, sz, "^", ms=4, color=cz, label=r"$\langle S_z\rangle$ (KPM)")
    ax.axhline(0.0, color="0.7", lw=0.6, zorder=0)
    ax.set_xlabel(r"time $t$ (fs)")
    ax.set_ylabel(r"normalised spin $\langle S_\alpha(t)\rangle/\langle S_x(0)\rangle$")
    ax.legend(loc="upper right", ncol=3)
    ax.set_ylim(-1.25, 1.6)
    prl.save(fig, os.path.join(out_dir, "fig_precession"))


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
    run_dir, out_dir = _parse_dirs(argv, "run_spin")
    os.makedirs(run_dir, exist_ok=True)
    os.makedirs(out_dir, exist_ok=True)
    _ensure_outputs(run_dir)
    figure_zeeman(run_dir)
    figure_precession(run_dir)
    if os.path.abspath(run_dir) != os.path.abspath(out_dir):
        for fn in ("fig_zeeman.png", "fig_precession.png"):
            src = os.path.join(run_dir, fn)
            if os.path.isfile(src):
                shutil.copy2(src, os.path.join(out_dir, fn))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
