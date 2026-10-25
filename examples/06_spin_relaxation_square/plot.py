#!/usr/bin/env python3
"""Generate the Tutorial 06 figures from the output directory of `python -m lsquant run`.

Usage:
    python -m lsquant run examples/06_spin_relaxation_square/run_relaxation.json --out /tmp/tut06
    python examples/06_spin_relaxation_square/plot.py /tmp/tut06
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

E_WORK = 2.0
EY_L = 20
EY_DS = [0.05, 0.07, 0.10]
DP_L = 20
DP_LAMBDA = 0.2
DP_WS = [0.0, 6.0]


def _load_normalised_spin(
    out_dir: str, label: str, op: str, proj: str, ef: float = E_WORK
) -> tuple[np.ndarray, np.ndarray]:
    num_path = sorted(
        glob.glob(
            os.path.join(
                out_dir, f"meanTimeProj{op}-{proj}{label}*EF{ef:.6f}JACKSON.dat"
            )
        )
    )[-1]
    den_path = sorted(
        glob.glob(
            os.path.join(
                out_dir, f"meanTimeProj1-{proj}{label}*EF{ef:.6f}JACKSON.dat"
            )
        )
    )[-1]
    t, num = np.loadtxt(num_path, unpack=True)
    _, den = np.loadtxt(den_path, unpack=True)
    good = t > 0.0
    s_curve = num[good] / float(np.mean(den[good]))
    return t[good], s_curve / float(s_curve[0])


def fit_rate(
    t: np.ndarray, s: np.ndarray, lo: float = 0.2, hi: float = 0.85
) -> float:
    mask = (s > lo) & (s < hi)
    if int(mask.sum()) < 3:
        return float("nan")
    slope = float(np.polyfit(t[mask], np.log(s[mask]), 1)[0])
    return -1.0 / slope


def _ensure_outputs(out_dir: str) -> None:
    expected = os.path.join(out_dir, "ey_ratelaw.dat")
    if not os.path.isfile(expected):
        run_json = os.path.join(_HERE, "run_relaxation.json")
        execute_runfile(run_json, out_dir_override=out_dir, repo_root=_REPO_ROOT)


def figure_ey(out_dir: str) -> None:
    """Part A: Elliott-Yafet spin relaxation curves and exponential fit."""
    curves = {}
    for d_val in EY_DS:
        label = f"ey_L{EY_L}_D{d_val:g}"
        t, s = _load_normalised_spin(out_dir, label, "SZ", "PZ")
        curves[d_val] = (t, s)
        t1 = fit_rate(t, s)
        print("EY  Dsf=%.2f  T1=%.1f fs  1/T1=%.4f" % (d_val, t1, 1.0 / t1))

    prl.use("web", aspect=0.46)
    fig, (a, b) = plt.subplots(1, 2)
    t, s = curves[0.10]
    t1 = fit_rate(t, s)
    a.plot(t, s, "o", ms=4, color=prl.COLORS[1], label=r"$S_z(E{=}2t,\,t)$")
    a.plot(
        t,
        np.exp(-t / t1),
        "-",
        color=prl.COLORS[1],
        alpha=0.55,
        label=r"$e^{-t/T_1},\ T_1\approx%.0f$ fs" % t1,
    )
    a.axhline(0, color="0.7", lw=0.6, zorder=0)
    a.set_xlabel(r"time $t$ (fs)")
    a.set_ylabel(r"$S_z(E,t)$")
    a.legend(loc="upper right")
    prl.panel(a, "(a)")

    cols = plt.cm.viridis(np.linspace(0.15, 0.8, len(EY_DS)))
    for d_val, c in zip(EY_DS, cols):
        tt, ss = curves[d_val]
        b.plot(tt, ss, color=c, lw=1.6, label=r"$\Delta_{\rm sf}=%.2f$" % d_val)
    b.axhline(0, color="0.7", lw=0.6, zorder=0)
    b.set_xlabel(r"time $t$ (fs)")
    b.set_ylabel(r"$S_z(E,t)$")
    b.legend(loc="upper right")
    prl.panel(b, "(b)")
    prl.save(fig, os.path.join(out_dir, "fig_ey"))


def figure_dp(out_dir: str) -> None:
    """Part B: Dyakonov-Perel precession vs relaxation and spin anisotropy."""
    curves = {}
    for w_val in DP_WS:
        label = f"rashba_L{DP_L}_l{DP_LAMBDA:g}_W{w_val:g}"
        curves[w_val] = _load_normalised_spin(out_dir, label, "SZ", "PZ")
    label_w = f"rashba_L{DP_L}_l{DP_LAMBDA:g}_W{DP_WS[-1]:g}"
    tx, sx = _load_normalised_spin(out_dir, label_w, "SX", "PX")

    prl.use("web", aspect=0.46)
    fig, (a, b) = plt.subplots(1, 2)
    a.plot(*curves[0.0], color=prl.COLORS[4], ls="-", lw=1.6, label=r"$W=0$ (precession)")
    a.plot(*curves[6.0], color=prl.COLORS[3], ls="--", lw=1.6, label=r"$W=6$ (relaxes)")
    a.axhline(0, color="0.7", lw=0.6, zorder=0)
    a.set_ylim(-1.25, 1.25)
    a.set_xlabel(r"time $t$ (fs)")
    a.set_ylabel(r"$S_z(E,t)$")
    a.legend(loc="upper right")
    prl.panel(a, "(a)")

    b.plot(
        *curves[6.0],
        color=prl.COLORS[3],
        ls="-",
        lw=1.6,
        label=r"$S_z$ (out of plane)",
    )
    b.plot(tx, sx, color=prl.COLORS[2], ls="--", lw=1.6, label=r"$S_x$ (in plane)")
    b.axhline(0, color="0.7", lw=0.6, zorder=0)
    b.set_xlabel(r"time $t$ (fs)")
    b.set_ylabel(r"$S(E,t)$")
    b.legend(loc="upper right")
    prl.panel(b, "(b)")
    prl.save(fig, os.path.join(out_dir, "fig_dp"))


def figure_ey_ratelaw(out_dir: str) -> None:
    """Elliott-Yafet rate law 1/T1 vs spin-flip disorder Dsf from random-phase trace."""
    d_arr, rate = np.loadtxt(os.path.join(out_dir, "ey_ratelaw.dat"), unpack=True)
    p, logc = np.polyfit(np.log(d_arr), np.log(rate), 1)

    prl.use("web", aspect=0.62)
    fig, ax = plt.subplots()
    ax.plot(
        d_arr,
        rate,
        label=r"random phase, $N_r{=}20$, $N{=}5000$",
        **prl.trace(0, line=False, marker=True),
    )
    xx = np.linspace(float(d_arr.min()), float(d_arr.max()), 50)
    ax.plot(
        xx,
        np.exp(logc) * xx**p,
        color=prl.COLORS[0],
        ls="-",
        lw=1.2,
        label=r"fit $1/T_1\propto\Delta_{\rm sf}^{%.2f}$" % p,
    )
    c2 = rate[0] / (d_arr[0] ** 2)
    ax.plot(
        xx,
        c2 * xx**2,
        color="0.45",
        ls="--",
        lw=1.1,
        label=r"Elliott-Yafet $\propto\Delta_{\rm sf}^{2}$",
    )
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(r"spin-flip disorder $\Delta_{\rm sf}\,[t]$")
    ax.set_ylabel(r"relaxation rate $1/T_1\,[{\rm fs}^{-1}]$")
    ax.legend(loc="upper left")
    prl.save(fig, os.path.join(out_dir, "fig_ey_ratelaw"))
    print("wrote fig_ey_ratelaw.png (fitted p = %.2f; theory p = 2)" % p)


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
    run_dir, out_dir = _parse_dirs(argv, "run_relaxation")
    os.makedirs(run_dir, exist_ok=True)
    os.makedirs(out_dir, exist_ok=True)
    _ensure_outputs(run_dir)
    figure_ey(run_dir)
    figure_dp(run_dir)
    figure_ey_ratelaw(run_dir)
    if os.path.abspath(run_dir) != os.path.abspath(out_dir):
        for fn in ("fig_ey.png", "fig_dp.png", "fig_ey_ratelaw.png", "ey_ratelaw.dat"):
            src = os.path.join(run_dir, fn)
            if os.path.isfile(src):
                shutil.copy2(src, os.path.join(out_dir, fn))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
