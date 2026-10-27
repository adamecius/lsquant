#!/usr/bin/env python3
"""Stage-3 figure (handoff 0003, codex.md section 3): tutorial 03 reproduced from Python.

Runs `python3 -m lsquant.apps.tutorial03 --out <tmp>` (the Python library over the
pybind11 module lsquant._core, PYTHONPATH=<build>/lsquant_core:python) in a temporary
directory. The app stages the Haldane operators (N = 512, window [-4.5, 4.5] eV), computes
the VX-VY and VX-VX moments by exact trace (M = 128), reconstructs Kubo-Bastin, Fermi sea,
Fermi surface and Kubo-Greenwood at 10 meV and writes the four .dat files with their legacy
names. Each is byte-compared (filecmp, shallow=False: the cmp criterion) with the file of
the same name in goldens/stage0/tutorials/03/, and one figure is drawn in the units of
examples/03_haldane_hall/lsqhall.py (sigma / ((N/2) A_cell e^2/h), A_cell = 3 sqrt(3)/2,
e^2/h = 1/(2 pi)):

  (a) sigma_xy: Bastin total, Fermi sea, Fermi surface; new as lines, golden as markers,
      residual new - golden below; the Chern number C from the gap plateau
      (C = (mean_{|E|<0.5} sigma_xy / A) 2 pi, A = (N/2) A_cell) in the legend.
  (b) sigma_xx (Kubo-Greenwood), same encoding.

Outputs figures/out/stage_3/fig_stage_3.png and fig_stage_3.caption.txt.
Usage (repository root, on TCN):
    MPLBACKEND=Agg python3 figures/stage_3.py [--build build] [--threads 16]
Nothing under goldens/ is written. Seeds and goldens are never changed here.
"""
import argparse
import filecmp
import glob
import math
import os
import re
import shutil
import subprocess
import sys
import tempfile

import numpy as np

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
GOLD = os.path.join(REPO, "goldens", "stage0", "tutorials", "03")
HAM = os.path.join(REPO, "goldens", "stage0", "haldane", "operators", "haldane.HAM.CSR")
M, ETA = 128, 10
XY = "NonEqOpVX-VYhaldaneKPM_M128x128_stateexact"
XX = "NonEqOpVX-VXhaldaneKPM_M128x128_stateexact"
FILES = {"sigma_xy": "KuboBastin_%sJACKSON.dat" % XY,
         "sea": "KuboBastinSea_%sJACKSON.dat" % XY,
         "surf": "KuboBastinSurf_%sJACKSON.dat" % XY,
         "sigma_xx": "KuboGreenWood_%sJACKSON.dat" % XX}
A_CELL = 3.0 * math.sqrt(3.0) / 2.0
E2_OVER_H = 1.0 / (2.0 * math.pi)


def find_pydir(build):
    for sub in ("lsquant_core", "lsquant"):
        d = os.path.join(build, sub)
        if glob.glob(os.path.join(d, "_core*.so")):
            return os.path.abspath(d)
    sys.exit("_core*.so not found under %s/{lsquant_core,lsquant}" % build)


def chern(dim, E, s):               # examples/03_haldane_hall/lsqhall.py chern_from_plateau
    area = (dim / 2) * A_CELL
    return (np.mean(s[np.abs(E) < 0.5]) / area) * 2 * np.pi


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--build", default=os.path.join(REPO, "build"))
    ap.add_argument("--threads", type=int, default=16)
    ap.add_argument("--out", default=os.path.join(REPO, "figures", "out", "stage_3"))
    a = ap.parse_args()
    pydir = find_pydir(a.build)
    os.makedirs(a.out, exist_ok=True)
    with open(HAM) as f:
        dim = int(f.readline().split()[0])

    tmp = tempfile.mkdtemp(prefix="stage3_tut03_")
    env = dict(os.environ, OMP_NUM_THREADS=str(a.threads), KPM_SEED="12345",
               PYTHONPATH=pydir + os.pathsep + os.path.join(REPO, "python"))
    r = subprocess.run([sys.executable, "-m", "lsquant.apps.tutorial03", "--out", tmp],
                       cwd=REPO, env=env, capture_output=True, text=True)
    printed = [l.replace(tmp, "<tmp>") for l in r.stdout.splitlines()]
    if r.returncode != 0:
        sys.exit("tutorial03 exited %d: %s" % (r.returncode, r.stderr.strip()))
    found = [re.search(r"C = \S+ \((\S+)\)", l) for l in printed]
    found = [x for x in found if x]
    c_printed = found[-1].group(1) if found else "not printed"

    # ---- byte comparison, by name
    cmp_lines, ok_all = [], True
    for key in ("sigma_xy", "sea", "surf", "sigma_xx"):
        name = FILES[key]
        new, gold = os.path.join(tmp, name), os.path.join(GOLD, name)
        ok = os.path.isfile(new) and filecmp.cmp(new, gold, shallow=False)
        ok_all = ok_all and ok
        cmp_lines.append("%s: %s" % (name, "IDENTICAL" if ok else "DIFFERS"))
    n_ident = sum(l.endswith("IDENTICAL") for l in cmp_lines)

    norm = (dim / 2) * A_CELL * E2_OVER_H
    data, raw, max_res, max_dE = {}, {}, 0.0, 0.0
    for key, name in FILES.items():
        new = np.loadtxt(os.path.join(tmp, name))
        gold = np.loadtxt(os.path.join(GOLD, name))
        if new.shape != gold.shape:
            sys.exit("line count differs for %s: %d new vs %d golden" % (name, len(new), len(gold)))
        raw[key] = (new, gold)
        data[key] = (new[:, 0], new[:, 1] / norm, gold[:, 0], gold[:, 1] / norm)
        max_res = max(max_res, float(np.max(np.abs(new[:, 1] - gold[:, 1]) / norm)))
        max_dE = max(max_dE, float(np.max(np.abs(new[:, 0] - gold[:, 0]))))
    new_xy, gold_xy = raw["sigma_xy"]
    C_new = chern(dim, new_xy[:, 0], new_xy[:, 1])
    C_gold = chern(dim, gold_xy[:, 0], gold_xy[:, 1])

    # ---- figure
    import matplotlib
    matplotlib.use(os.environ.get("MPLBACKEND", "Agg"))
    import matplotlib.pyplot as plt
    sys.path.insert(0, os.path.join(REPO, "examples"))
    try:
        import prl_style as prl
        prl.use("double", aspect=0.5, usetex=False)
        CO, L, K = prl.COLORS, prl.LINES, prl.MARKERS
        tag = prl.panel
    except ImportError:              # conventions inlined if the style module is absent
        CO, L, K = ["#000000", "#0072B2", "#D55E00", "#009E73"], ["-", "--", "-.", ":"], ["o", "s", "^", "v"]

        def tag(ax, t, loc="upper left", **kw):
            ax.text(0.04, 0.96, t, transform=ax.transAxes, ha="left", va="top", fontweight="bold", **kw)

    fig = plt.figure(figsize=(7.0, 3.6))
    gs = fig.add_gridspec(2, 2, height_ratios=[3.0, 1.0], hspace=0.12, wspace=0.28)
    nmark = 48                        # golden markers: 48 evenly spaced points per curve

    def panel(col, series, ylabel, letter, title, ytop, ncol):
        ax = fig.add_subplot(gs[0, col])
        axr = fig.add_subplot(gs[1, col], sharex=ax)
        zero = True
        for i, (key, lab, ci) in enumerate(series):
            E, y, Eg, yg = data[key]
            c = CO[ci % len(CO)]
            ax.plot(E, y, color=c, lw=0.7, linestyle=L[0], zorder=10 - 2 * i, label=lab + " new")
            k = np.linspace(0, Eg.size - 1, min(nmark, Eg.size)).astype(int)
            ax.plot(Eg[k], yg[k], linestyle="none", marker=K[i % len(K)], ms=2.8, mfc="none",
                    mew=0.6, color=c, zorder=11 - 2 * i, label=lab + " golden")
            axr.plot(E, y - yg, color=c, lw=0.8, linestyle=L[i % len(L)])
            zero = zero and float(np.max(np.abs(y - yg))) == 0.0
        ax.set_xlim(-4.4, 4.4)
        ax.set_ylim(top=ytop)
        ax.set_ylabel(ylabel)
        ax.tick_params(labelbottom=False)
        ax.legend(loc="upper left", ncol=ncol, fontsize=5.5, title=title, title_fontsize=6,
                  handlelength=1.6, columnspacing=0.8, frameon=False)
        tag(ax, letter, loc="upper right")
        axr.axhline(0.0, color="0.7", lw=0.4)
        if zero:
            axr.set_ylim(-1e-15, 1e-15)
            axr.text(0.5, 0.75, "identically 0", transform=axr.transAxes, ha="center", va="center",
                     fontsize=6)
        axr.set_xlabel(r"energy $E$ (eV)")
        axr.set_ylabel(r"new $-$ golden", fontsize=7)
        axr.tick_params(axis="y", labelsize=6)
        axr.yaxis.get_offset_text().set_fontsize(6)
        return ax

    ax = panel(0, [("sigma_xy", "Bastin (total)", 0), ("sea", "Fermi sea", 1), ("surf", "Fermi surface", 2)],
               r"$\sigma_{xy}$ ($e^2/h$)", "(a)", r"$C=%+.3f$ (gap plateau)" % C_new,
               1.62, 2)
    ax.axhline(1.0, color="0.6", ls=":", lw=0.6)
    panel(1, [("sigma_xx", "Kubo-Greenwood", 3)], r"$\sigma_{xx}$ ($e^2/h$)", "(b)",
          r"$\eta=%d$ meV, $M=%d$" % (ETA, M), 17.0, 1)

    png = os.path.join(a.out, "fig_stage_3.png")
    fig.savefig(png, dpi=300, bbox_inches="tight")

    caption = (
        "FIG. 3. Stage 3 (handoff 0003): tutorial 03, the Hall and longitudinal conductivity of the "
        "Haldane model (t1 = -1, t2 = 0.15, phi = pi/2), reproduced from Python by "
        "python3 -m lsquant.apps.tutorial03 (lsquant package over the pybind11 module lsquant._core, "
        "OMP_NUM_THREADS = %d): N = %d orbitals, window [-4.5, 4.5] eV, exact trace over all N site "
        "states, M = %d x %d Chebyshev moments (VX-VY for sigma_xy, VX-VX for sigma_xx), Jackson kernel "
        "at eta = %d meV; conductivities in units of e^2/h as lsqhall.py normalises them (divided by "
        "(N/2) A_cell e^2/h, A_cell = 3 sqrt(3)/2, e = hbar = 1). "
        "(a) sigma_xy: Kubo-Bastin total, Fermi-sea and Fermi-surface parts; new as lines, the stage-0 "
        "golden (goldens/stage0/tutorials/03) as open markers (%d evenly spaced points per curve), "
        "dotted line at 1; residual new - golden below. "
        "(b) sigma_xx (Kubo-Greenwood), same encoding. "
        "Chern number from the gap plateau, C = (mean over |E| < 0.5 eV of sigma_xy / A) 2 pi with "
        "A = (N/2) A_cell: C = %+.12f from the new curve (tutorial03 prints %s), C = %+.12f from the "
        "golden curve. Max |new - golden| over the four curves: %.1e e^2/h, max energy-axis "
        "difference %.1e eV. Byte comparison with goldens/stage0/tutorials/03 by name (cmp criterion), "
        "%d of %d IDENTICAL: %s. %s"
        % (a.threads, dim, M, M, ETA, nmark, C_new, c_printed, C_gold, max_res, max_dE,
           n_ident, len(cmp_lines), "; ".join(cmp_lines),
           "All files IDENTICAL." if ok_all else "AT LEAST ONE FILE DIFFERS."))
    with open(os.path.join(a.out, "fig_stage_3.caption.txt"), "w") as f:
        f.write(caption + "\n")
    print("tutorial03 printed (%d lines):" % len(printed))
    for l in printed:
        print("  " + l)
    print("cmp list:")
    for l in cmp_lines:
        print("  " + l)
    print("C new = %+.12f, C golden = %+.12f" % (C_new, C_gold))
    print("wrote", png)
    print(caption)
    shutil.rmtree(tmp)
    return 0 if ok_all else 1


if __name__ == "__main__":
    sys.exit(main())
