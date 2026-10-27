#!/usr/bin/env python3
"""Generate the Tutorial 02 figures from the output directory of `python -m lsquant run`.

Usage:
    python -m lsquant run examples/02_localization_1d_chain/run_localization.json --out /tmp/tut02
    python examples/02_localization_1d_chain/plot.py /tmp/tut02
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


def xi_from_plateau(t: np.ndarray, msd: np.ndarray, frac: float = 0.3) -> float:
    """Relative localization length from the late-time plateau of the MSD."""
    cut = t >= (1.0 - frac) * float(t.max())
    return float(np.sqrt(np.mean(msd[cut])))


def xi_born(energies: np.ndarray | float, w_disorder: float) -> np.ndarray:
    """Weak-disorder (Born) localization length for box disorder, t = 1."""
    return 24.0 * (4.0 - np.asarray(energies, float) ** 2) / (w_disorder**2)


def _ensure_outputs(out_dir: str, sizes: list[int], realisations: int) -> None:
    r_dir = os.path.join(out_dir, f"r{realisations:02d}")
    pattern = os.path.join(r_dir, f"mean*chain1d_dis_N{sizes[-1]}*.dat")
    if not glob.glob(pattern):
        run_json = os.path.join(_HERE, "run_localization.json")
        execute_runfile(run_json, out_dir_override=out_dir, repo_root=_REPO_ROOT)


def _load_averaged_msd(
    out_dir: str,
    n_sites: int,
    w_disorder: float,
    energies: list[float],
    realisations: int,
) -> tuple[np.ndarray, dict[float, np.ndarray]]:
    label = f"chain1d_dis_N{n_sites}_W{w_disorder:g}"
    acc: dict[float, np.ndarray | None] = {e: None for e in energies}
    t_axis: np.ndarray | None = None
    for r in range(1, realisations + 1):
        r_dir = os.path.join(out_dir, f"r{r:02d}")
        for e in energies:
            matches = sorted(
                glob.glob(os.path.join(r_dir, f"mean*{label}*EF{float(e):.6f}*.dat"))
            )
            if not matches:
                raise FileNotFoundError(
                    f"Missing MSD curve for {label} EF={e:.6f} in {r_dir}"
                )
            d = np.loadtxt(matches[-1])
            t_axis = d[:, 0] if t_axis is None else t_axis
            acc[e] = d[:, 1].copy() if acc[e] is None else acc[e] + d[:, 1]
    assert t_axis is not None
    out: dict[float, np.ndarray] = {}
    for e in energies:
        assert acc[e] is not None
        out[e] = acc[e] / float(realisations)
    return t_axis, out


def figure_anderson(out_dir: str, num_sites: int = 500) -> None:
    """Opening figure: band-centre eigenstate, extended (clean) vs localized (W=4)."""

    def state(w_val: float, seed: int) -> np.ndarray:
        rng = np.random.default_rng(seed)
        ham = np.diag(rng.uniform(-w_val / 2.0, w_val / 2.0, num_sites)).astype(float)
        idx = np.arange(num_sites)
        ham[idx, (idx + 1) % num_sites] = -1.0
        ham[(idx + 1) % num_sites, idx] = -1.0
        w_eig, v_eig = np.linalg.eigh(ham)
        return np.abs(v_eig[:, int(np.argmin(np.abs(w_eig)))]) ** 2

    prl.use("web", aspect=0.45)
    fig, ax = plt.subplots()
    ax.plot(
        state(0.0, 3),
        color=prl.COLORS[1],
        lw=1.2,
        ls="-",
        label=r"clean: extended ($|\psi|^2\sim 1/N$)",
    )
    ax.plot(
        state(4.0, 3),
        color=prl.COLORS[2],
        lw=1.2,
        ls="--",
        label=r"$W=4$: localized",
    )
    ax.set_xlabel(r"site $i$")
    ax.set_ylabel(r"$|\psi_i|^2$ (band-centre state)")
    ax.legend(loc="upper right")
    prl.save(fig, os.path.join(out_dir, "fig_anderson"))


def figure_localization(
    out_dir: str,
    sizes: list[int],
    w_disorder: float,
    energies: list[float],
    realisations: int = 16,
) -> None:
    """Two-panel localization figure from ensemble-averaged MSD plateaus."""
    xi: dict[int, dict[float, float]] = {n_sites: {} for n_sites in sizes}
    for n_sites in sizes:
        t_axis, acc = _load_averaged_msd(
            out_dir, n_sites, w_disorder, energies, realisations
        )
        for e in energies:
            xi[n_sites][e] = xi_from_plateau(t_axis, acc[e])

    n_ref, e_ref = max(sizes), 1.0
    cal = float(xi_born(e_ref, w_disorder)) / xi[n_ref][e_ref]

    prl.use("web", height=5.8)
    fig, (top, bot) = plt.subplots(2, 1)

    egrid = np.linspace(-1.95, 1.95, 200)
    top.plot(
        egrid,
        xi_born(egrid, w_disorder),
        color="0.2",
        ls="--",
        lw=1.2,
        label="Born theory",
    )
    for i, n_sites in enumerate(sizes):
        e_sorted = sorted(energies)
        top.plot(
            e_sorted,
            [cal * xi[n_sites][e] for e in e_sorted],
            lw=1.0,
            label=r"$N=%d$" % n_sites,
            **prl.trace(i + 1, marker=True),
        )
    top.set_xlabel(r"energy $E$ (eV)")
    top.set_ylabel(r"$\xi(E)$ (sites)")
    prl.panel(top, "(a)")
    top.legend(loc="upper right", title=r"$W=%g$ eV" % w_disorder)

    for i, e in enumerate(sorted(energies)):
        bot.plot(
            sizes,
            [cal * xi[n_sites][e] for n_sites in sizes],
            lw=1.0,
            label=r"$E=%+.1f$ eV" % e,
            **prl.trace(i, marker=True),
        )
        bot.axhline(float(xi_born(e, w_disorder)), color="0.6", lw=0.6, ls=":")
    bot.set_xscale("log", base=2)
    bot.set_xlabel(r"system size $N$ (sites)")
    bot.set_ylabel(r"$\xi$ (sites)")
    prl.panel(bot, "(b)", loc="lower left")
    bot.legend(loc="upper center", ncol=2)

    prl.save(fig, os.path.join(out_dir, "fig_localization"))


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
    run_dir, out_dir = _parse_dirs(argv, "run_localization")
    os.makedirs(run_dir, exist_ok=True)
    os.makedirs(out_dir, exist_ok=True)
    sizes = [256, 1024, 4096]
    energies = [0.0, 1.0, -1.0, 1.9]
    w_disorder = 1.0
    realisations = 16
    _ensure_outputs(run_dir, sizes, realisations)
    figure_anderson(run_dir)
    figure_localization(run_dir, sizes, w_disorder, energies, realisations=realisations)
    if os.path.abspath(run_dir) != os.path.abspath(out_dir):
        for fn in ("fig_anderson.png", "fig_localization.png"):
            src = os.path.join(run_dir, fn)
            if os.path.isfile(src):
                shutil.copy2(src, os.path.join(out_dir, fn))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
