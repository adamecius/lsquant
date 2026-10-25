#!/usr/bin/env python3
"""Generate fig_dos_shc.png for Example 08 (pentagonal PdSe2 monolayer spin Hall conductivity).

Usage:
    python -m lsquant run examples/08_spin_hall_pdse2_monolayer/run_dos_quick.json --out /tmp/pdse2
    python -m lsquant run examples/08_spin_hall_pdse2_monolayer/run_shc_quick.json --out /tmp/pdse2
    python examples/08_spin_hall_pdse2_monolayer/plot.py /tmp/pdse2
"""

from __future__ import annotations

import argparse
import json
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


def _resolve_sub_dir(root_dir: str, filename: str, candidates: tuple[str, ...]) -> str:
    if os.path.isfile(os.path.join(root_dir, filename)):
        return root_dir
    for sub in candidates:
        cand = os.path.join(root_dir, sub)
        if os.path.isfile(os.path.join(cand, filename)):
            return cand
    return root_dir


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("out_dir", nargs="?", default=".", help="Directory containing dos.dat and shc.dat")
    parser.add_argument("--run-dir", default=None, help="Run directory containing dos.dat/shc.dat or run_* subdirs")
    parser.add_argument("--dos-dir", default=None, help="Directory containing dos.dat (if separate)")
    parser.add_argument("--shc-dir", default=None, help="Directory containing shc.dat (if separate)")
    parser.add_argument("--out", default=None, help="Output directory for fig_dos_shc.png")
    parser.add_argument("--emin", type=float, default=-1.0, help="Minimum energy (eV) for plot window")
    parser.add_argument("--emax", type=float, default=2.5, help="Maximum energy (eV) for plot window")
    args = parser.parse_args(argv)

    run_root = os.path.abspath(args.run_dir if args.run_dir is not None else args.out_dir)
    base_dir = os.path.abspath(args.out if args.out is not None else run_root)
    dos_dir = (
        os.path.abspath(args.dos_dir)
        if args.dos_dir is not None
        else _resolve_sub_dir(run_root, "dos.dat", ("run_dos_production", "run_dos_quick"))
    )
    shc_dir = (
        os.path.abspath(args.shc_dir)
        if args.shc_dir is not None
        else _resolve_sub_dir(run_root, "shc.dat", ("run_shc_production", "run_shc_quick"))
    )
    os.makedirs(base_dir, exist_ok=True)

    dos_path = os.path.join(dos_dir, "dos.dat")
    shc_path = os.path.join(shc_dir, "shc.dat")
    if not os.path.isfile(dos_path):
        execute_runfile(
            os.path.join(_HERE, "run_dos_quick.json"),
            out_dir_override=dos_dir,
            repo_root=_REPO_ROOT,
        )
    if not os.path.isfile(shc_path):
        execute_runfile(
            os.path.join(_HERE, "run_shc_quick.json"),
            out_dir_override=shc_dir,
            repo_root=_REPO_ROOT,
        )

    e_dos, v_dos = np.loadtxt(dos_path, unpack=True)
    e_shc, v_shc = np.loadtxt(shc_path, unpack=True)

    prl.use("web", height=7.2)
    fig, (top, mid, bot) = plt.subplots(
        3, 1, sharex=True, gridspec_kw={"height_ratios": [1.0, 1.15, 0.65]}
    )

    mask_dos = (e_dos >= args.emin) & (e_dos <= args.emax)
    top.plot(
        e_dos[mask_dos],
        v_dos[mask_dos],
        color=prl.COLORS[1],
        lw=1.5,
        label=r"KPM $\mathrm{Tr}[\delta(E-H)]$",
    )
    top.axvline(0.0, color="0.6", ls=":", lw=0.8)
    top.set_ylabel(r"DOS (states / eV / cell)")
    prl.panel(top, "(a)")
    top.legend(loc="upper right")

    mask_shc = (e_shc >= args.emin) & (e_shc <= args.emax)
    mid.plot(
        e_shc[mask_shc],
        v_shc[mask_shc],
        color=prl.COLORS[0],
        lw=1.6,
        label=r"KPM Kubo-Bastin $\sigma^z_{xy}$",
    )

    oracle_path = os.path.join(shc_dir, "exact_shc.json")
    if os.path.isfile(oracle_path):
        with open(oracle_path, encoding="utf-8") as fh:
            oracle = json.load(fh)
        e_or = np.asarray(oracle["energies"], dtype=float)
        v_or = np.asarray(oracle["values"], dtype=float)
        mask_or = (e_or >= args.emin) & (e_or <= args.emax)
        if np.any(mask_or):
            mid.plot(
                e_or[mask_or],
                v_or[mask_or],
                color=prl.COLORS[2],
                ls="--",
                lw=1.3,
                label="exact eigenbasis oracle",
            )
            v_kpm_on_or = np.interp(e_or[mask_or], e_shc, v_shc)
            bot.plot(
                e_or[mask_or],
                v_kpm_on_or - v_or[mask_or],
                color=prl.COLORS[3],
                lw=1.2,
                label=r"$\sigma^z_{xy,\mathrm{KPM}} - \sigma^z_{xy,\mathrm{oracle}}$",
            )

    mid.axvline(0.0, color="0.6", ls=":", lw=0.8)
    mid.axhline(0.0, color="0.7", lw=0.6, zorder=0)
    mid.set_ylabel(r"$\sigma^z_{xy}(E_F)$ ($e/4\pi$)")
    prl.panel(mid, "(b)")
    mid.legend(loc="upper right")

    bot.axvline(0.0, color="0.6", ls=":", lw=0.8)
    bot.axhline(0.0, color="0.7", lw=0.6, zorder=0)
    bot.set_xlim(args.emin, args.emax)
    bot.set_xlabel(r"Fermi energy $E_F$ (eV)")
    bot.set_ylabel(r"residual ($e/4\pi$)")
    prl.panel(bot, "(c)")
    if os.path.isfile(oracle_path):
        bot.legend(loc="upper right")

    prl.save(fig, os.path.join(base_dir, "fig_dos_shc"))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
