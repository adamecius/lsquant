#!/usr/bin/env python3
"""Stage-4 figure (handoff 0004, codex.md section 3): tutorials 04, 02 and 05 from Python.

Runs, each in its own temporary directory and with its wall time recorded,
`python3 -m lsquant.apps.tutorial04 --out <tmp>`, `... tutorial02 --out <tmp>` and
`... tutorial05 --out <tmp>` (the Python library over the pybind11 module lsquant._core,
PYTHONPATH=<build>/lsquant_core:python, KPM_SEED=12345 in the environment; the apps set
their own trace seeds as the legacy tutorial scripts do). Every file the apps write is
byte-compared (filecmp, shallow=False: the cmp criterion) with goldens/stage0/tutorials/0N:

  * captured/NNNN_<name>: the golden capture with the same NNNN_<name> (same read order
    and same name). Several captures share a name (tutorial 02: 16 disorder realisations
    per size and energy); each is compared with its own counterpart. A DIFFERS capture is
    also checked against every new capture of the same <name> (reported as "matches new
    capture k" if one is identical, which would mean an ordering difference only).
  * the top-level *.dat of the golden tutorial directory: the new file of the same name.

One figure, three columns, curve on top and residual below:
  (a) tutorial 04: <S_x>, <S_y>, <S_z> (m = 0 row of the TimeProj tensors, normalised by
      <S_x(0)> as lsqspin.py does) vs cos(wt), -sin(wt), 0 with w = 2 J_ex / hbar; golden
      as open markers; residual new - analytic (lines) and new - golden (markers).
  (b) tutorial 02: the localisation panel of lsqloc.py, xi(E) = cal sqrt(<MSD plateau>)
      over the last 30 % of the time window, averaged over 16 realisations, for
      N = 256, 1024, 4096, calibrated once against the Born length at N = 4096, E = 1;
      new lines, golden markers, Born curve dashed; residual new - golden.
  (c) tutorial 05: D(E, t) = DeltaX^2 / 2t at E = 2 eV for W = 0.5, 1, 2, 4; new lines,
      golden markers; residual new - golden.

A tutorial app that exits non-zero does not stop the script: its files count as DIFFERS
(not produced) and its panel shows the golden only, with the error in the residual axis.
--reuse NN=DIR takes tutorial NN from an existing output directory instead of running the
app (diagnostics only; written into the caption), --name sets the output basename.

Outputs figures/out/stage_4/fig_stage_4.png, fig_stage_4.caption.txt and fig_stage_4.cmp.txt
(the per-file cmp list).
Usage (repository root, on TCN):
    MPLBACKEND=Agg python3 figures/stage_4.py [--build build] [--threads 16]
Nothing under goldens/ is written. Seeds and goldens are never changed here.
"""
import argparse
import filecmp
import glob
import os
import re
import shutil
import subprocess
import sys
import tempfile
import time

import numpy as np

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
GOLD = os.path.join(REPO, "goldens", "stage0", "tutorials")
HBAR = 0.6582119624          # eV fs (units.hpp; lsqspin.py)
J_EX = 0.1                   # eV (lsqspin.py)
T04 = dict(M=32, NTIME=96, TMAX=95.0)
T02 = dict(W=1.0, SIZES=(256, 1024, 4096), ENERGIES=(0.0, 1.0, -1.0, 1.9), NREAL=16,
           M=256, NTIME=128, TMAX=600.0, ETA=10)
T05 = dict(L=200, W=(0.5, 1.0, 2.0, 4.0), E=2.0, M=256, NTIME=160, TMAX=150.0, ETA=80, MDOS=1024)


def find_pydir(build):
    for sub in ("lsquant_core", "lsquant"):
        d = os.path.join(build, sub)
        if glob.glob(os.path.join(d, "_core*.so")):
            return os.path.abspath(d)
    sys.exit("_core*.so not found under %s/{lsquant_core,lsquant}" % build)


def run_app(app, env, reuse=None):
    """Run one tutorial app into a fresh temporary directory (or take a reused one).

    Returns (dir, wall seconds or None, lines printed, error string or None). A failing app
    does not stop the figure: its files count as DIFFERS (not produced) and its panel shows
    the golden only.
    """
    if reuse:
        return os.path.abspath(reuse), None, 0, None
    tmp = tempfile.mkdtemp(prefix="stage4_%s_" % app)
    t0 = time.time()
    r = subprocess.run([sys.executable, "-m", "lsquant.apps.%s" % app, "--out", tmp],
                       cwd=REPO, env=env, capture_output=True, text=True)
    wall = time.time() - t0
    err = None
    if r.returncode != 0:
        err = "%s exited %d: %s" % (app, r.returncode, r.stderr.strip().splitlines()[-1] if r.stderr.strip() else "")
    return tmp, wall, len(r.stdout.splitlines()), err


def compare_tutorial(num, tmp):
    """cmp every golden capture and top-level .dat of tutorials/<num> with the new files."""
    gdir = os.path.join(GOLD, num)
    lines, n_ok, n_all = [], 0, 0
    newcap = os.path.join(tmp, "captured")
    new_by_name = {}
    for p in sorted(glob.glob(os.path.join(newcap, "[0-9][0-9][0-9][0-9]_*"))):
        new_by_name.setdefault(os.path.basename(p)[5:], []).append(p)
    for g in sorted(glob.glob(os.path.join(gdir, "captured", "[0-9][0-9][0-9][0-9]_*"))):
        base = os.path.basename(g)
        n_all += 1
        new = os.path.join(newcap, base)
        if not os.path.isfile(new):
            lines.append("captured/%s: DIFFERS (not produced)" % base)
            continue
        if filecmp.cmp(new, g, shallow=False):
            n_ok += 1
            lines.append("captured/%s: IDENTICAL" % base)
        else:
            alt = [os.path.basename(p)[:4] for p in new_by_name.get(base[5:], [])
                   if filecmp.cmp(p, g, shallow=False)]
            lines.append("captured/%s: DIFFERS%s" % (base, (" (matches new capture %s)" % ",".join(alt)) if alt else ""))
    for g in sorted(glob.glob(os.path.join(gdir, "*.dat"))):
        base = os.path.basename(g)
        n_all += 1
        new = os.path.join(tmp, base)
        if not os.path.isfile(new):
            lines.append("%s: DIFFERS (not produced)" % base)
        elif filecmp.cmp(new, g, shallow=False):
            n_ok += 1
            lines.append("%s: IDENTICAL" % base)
        else:
            lines.append("%s: DIFFERS" % base)
    n_new = sum(len(v) for v in new_by_name.values())
    return lines, n_ok, n_all, n_new


def m0_row(path):            # lsqspin.py _m0: line 1 header, line 2 "M Nt", m-major re/im pairs
    with open(path) as f:
        f.readline()
        nt = int(f.readline().split()[1])
        return np.array([float(f.readline().split()[0]) for _ in range(nt)])



def xi_born(E, W):
    """Weak-disorder (Born) localisation length for box disorder, t = 1 (lsqloc.py)."""
    return 24.0 * (4.0 - np.asarray(E, float) ** 2) / W ** 2


def xi_table(root):
    """lsqloc.py figure_localization: xi_rel(N, E) from the 16-realisation mean MSD plateau."""
    xi, idx = {}, 0
    for N in T02["SIZES"]:
        acc = {E: None for E in T02["ENERGIES"]}
        tax = None
        for r in range(1, T02["NREAL"] + 1):
            for E in T02["ENERGIES"]:
                name = ("%04d_meanCorrelationVXchain1d_dis_N%d_W1KPM_M256_state600defaultEF%.6fJACKSON.dat"
                        % (idx, N, E))
                d = np.loadtxt(os.path.join(root, "captured", name))
                tax = d[:, 0] if tax is None else tax
                acc[E] = d[:, 1].copy() if acc[E] is None else acc[E] + d[:, 1]
                idx += 1
        for E in T02["ENERGIES"]:
            m = acc[E] / T02["NREAL"]
            cut = tax >= 0.7 * tax.max()
            xi[(N, E)] = float(np.sqrt(np.mean(m[cut])))
    cal = xi_born(1.0, T02["W"]) / xi[(max(T02["SIZES"]), 1.0)]
    return {k: cal * v for k, v in xi.items()}


def d_curves(root):
    """lsqtransport.py diffusion(): D(E, t) = DeltaX^2 / 2t for t > 0, per W."""
    out = {}
    for i, W in enumerate(T05["W"]):
        name = "%04d_meanCorrelationVXsquare_L200_W%gKPM_M256_state150defaultEF2.000000JACKSON.dat" % (2 * i, W)
        d = np.loadtxt(os.path.join(root, "captured", name))
        g = d[:, 0] > 0
        out[W] = (d[g, 0], d[g, 1] / (2.0 * d[g, 0]))
    return out


def spin_series(root):
    """lsqspin.py: m = 0 row of TimeProj{SX,SY,SZ}-PX normalised by <S_x(0)>."""
    ser = {}
    for i, op in enumerate(("SX", "SY", "SZ")):
        name = "%04d_TimeProj%s-PXchain1d_magKPM_M32_stateexact.chebmomTD" % (i, op)
        ser[op] = m0_row(os.path.join(root, "captured", name))
    D = ser["SX"][0]
    return {k: v / D for k, v in ser.items()}


def try_load(fn, root):
    try:
        return fn(root), None
    except Exception as exc:  # a failed app leaves no files: the panel shows the golden only
        return None, "%s: %s" % (type(exc).__name__, exc)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--build", default=os.path.join(REPO, "build"))
    ap.add_argument("--threads", type=int, default=16)
    ap.add_argument("--out", default=os.path.join(REPO, "figures", "out", "stage_4"))
    ap.add_argument("--name", default="fig_stage_4", help="output basename")
    ap.add_argument("--reuse", action="append", default=[], metavar="NN=DIR",
                    help="take tutorial NN's output from DIR instead of running the app "
                         "(diagnostics only; recorded in the caption)")
    a = ap.parse_args()
    pydir = find_pydir(a.build)
    os.makedirs(a.out, exist_ok=True)
    env = dict(os.environ, OMP_NUM_THREADS=str(a.threads), KPM_SEED="12345",
               PYTHONPATH=pydir + os.pathsep + os.path.join(REPO, "python"))
    reuse = dict(r.split("=", 1) for r in a.reuse)

    runs, cmps = {}, {}
    for num, app in (("04", "tutorial04"), ("02", "tutorial02"), ("05", "tutorial05")):
        tmp, wall, nprint, err = run_app(app, env, reuse.get(num))
        runs[num] = (tmp, wall, nprint, err)
        cmps[num] = compare_tutorial(num, tmp)
        print("%s: %s, %d lines printed, %d of %d IDENTICAL (%d new captures)%s"
              % (app, ("reused %s" % tmp) if wall is None else ("wall %.1f s" % wall), nprint,
                 cmps[num][1], cmps[num][2], cmps[num][3], ("; " + err) if err else ""), flush=True)

    gold04 = spin_series(os.path.join(GOLD, "04"))
    gold02 = xi_table(os.path.join(GOLD, "02"))
    gold05 = d_curves(os.path.join(GOLD, "05"))
    new04, lerr04 = try_load(spin_series, runs["04"][0])
    new02, lerr02 = try_load(xi_table, runs["02"][0])
    new05, lerr05 = try_load(d_curves, runs["05"][0])

    nt = gold04["SX"].size
    dt = T04["TMAX"] / (T04["NTIME"] - 1)
    t = np.arange(nt) * dt
    w = 2 * J_EX / HBAR
    ana = {"SX": np.cos(w * t), "SY": -np.sin(w * t), "SZ": np.zeros_like(t)}
    nan = float("nan")
    res04_ana = max(float(np.max(np.abs(gold04[k] - ana[k]))) for k in ana)   # golden vs analytic
    res04_gold = max(float(np.max(np.abs(new04[k] - gold04[k]))) for k in ana) if new04 else nan
    res02 = max(abs(new02[k] - gold02[k]) for k in gold02) if new02 else nan
    res05 = max(float(np.max(np.abs(new05[W][1] - gold05[W][1]))) for W in T05["W"]) if new05 else nan

    # ---- figure
    import matplotlib
    matplotlib.use(os.environ.get("MPLBACKEND", "Agg"))
    import matplotlib.pyplot as plt
    sys.path.insert(0, os.path.join(REPO, "examples"))
    try:
        import prl_style as prl
        prl.use("double", aspect=0.42, usetex=False)
        CO, L, K = prl.COLORS, prl.LINES, prl.MARKERS
    except ImportError:
        CO, L, K = ["#000000", "#0072B2", "#D55E00", "#009E73"], ["-", "--", "-.", ":"], ["o", "s", "^", "v"]

    def tag(ax, s):
        ax.text(0.03, 0.96, s, transform=ax.transAxes, ha="left", va="top", fontweight="bold")

    def residual_axis(axr, res, err):
        axr.axhline(0.0, color="0.7", lw=0.4)
        axr.tick_params(axis="y", labelsize=6)
        axr.yaxis.get_offset_text().set_fontsize(6)
        if err:
            import textwrap
            axr.text(0.5, 0.5, "no new data:\n" + "\n".join(textwrap.wrap(err, 44)[:3]),
                     transform=axr.transAxes, ha="center", va="center", fontsize=4.5, color="#D55E00")
        elif res == 0.0:
            axr.text(0.5, 0.9, r"new $-$ golden identically 0", transform=axr.transAxes,
                     ha="center", va="top", fontsize=5.5,
                     bbox=dict(facecolor="white", edgecolor="none", pad=0.5))

    fig = plt.figure(figsize=(7.0, 3.3))
    gs = fig.add_gridspec(2, 3, height_ratios=[3.0, 1.1], hspace=0.12, wspace=0.36)

    # (a) tutorial 04
    ax = fig.add_subplot(gs[0, 0])
    axr = fig.add_subplot(gs[1, 0], sharex=ax)
    tt = np.linspace(0.0, t[-1], 400)
    labs = {"SX": r"$\langle S_x\rangle$", "SY": r"$\langle S_y\rangle$", "SZ": r"$\langle S_z\rangle$"}
    curves = {"SX": np.cos(w * tt), "SY": -np.sin(w * tt), "SZ": np.zeros_like(tt)}
    for i, op in enumerate(("SX", "SY", "SZ")):
        c = CO[(i + 1) % len(CO)]
        ax.plot(tt, curves[op], color=c, lw=0.6, ls=L[i], alpha=0.6)
        if new04:
            ax.plot(t, new04[op], ls="none", marker=K[i], ms=2.0, color=c, label=labs[op] + " new")
            axr.plot(t, new04[op] - ana[op], color=c, lw=0.7, ls=L[i])
        ax.plot(t[::4], gold04[op][::4], ls="none", marker=K[i], ms=4.2, mfc="none", mew=0.5, color="0.35")
    ax.plot([], [], ls="none", marker="o", mfc="none", color="0.35", label="golden")
    ax.plot([], [], color="0.4", lw=0.6, label=r"$\cos\omega t,\,-\sin\omega t,\,0$")
    ax.set_ylim(-1.25, 2.05)
    ax.set_ylabel(r"$\langle S_\alpha(t)\rangle/\langle S_x(0)\rangle$")
    ax.tick_params(labelbottom=False)
    ax.legend(loc="upper right", ncol=2, fontsize=5, frameon=False, handlelength=1.4, columnspacing=0.6)
    tag(ax, "(a)")
    residual_axis(axr, res04_gold, runs["04"][3] or lerr04)
    axr.set_xlabel(r"time $t$ (fs)")
    axr.set_ylabel(r"new $-$ analytic", fontsize=6.5)

    # (b) tutorial 02
    ax = fig.add_subplot(gs[0, 1])
    axr = fig.add_subplot(gs[1, 1], sharex=ax)
    eg = np.linspace(-1.95, 1.95, 200)
    ax.plot(eg, xi_born(eg, T02["W"]), color="0.2", ls="--", lw=0.8, label="Born")
    Es = sorted(T02["ENERGIES"])
    for i, N in enumerate(T02["SIZES"]):
        c = CO[(i + 1) % len(CO)]
        yg = np.array([gold02[(N, E)] for E in Es])
        ax.plot(Es, yg, ls="none", marker=K[i], ms=3.6, mfc="none", mew=0.6, color=c,
                label=(r"$N=%d$" % N) + ("" if new02 else " golden"))
        if new02:
            yn = np.array([new02[(N, E)] for E in Es])
            ax.plot(Es, yn, color=c, lw=0.8, ls=L[i])
            axr.plot(Es, yn - yg, color=c, lw=0.7, ls=L[i], marker=K[i], ms=2)
    if new02:
        ax.plot([], [], color="0.35", lw=0.8, label="new (lines)")
        ax.plot([], [], ls="none", marker="o", mfc="none", color="0.35", label="golden")
    ax.set_ylabel(r"$\xi(E)$ (sites)")
    ax.tick_params(labelbottom=False)
    ax.set_ylim(0, 112)
    ax.legend(loc="lower center", ncol=3 if new02 else 2, fontsize=5, frameon=False, title=r"$W=1$ eV", title_fontsize=5.5,
              columnspacing=0.8, handlelength=1.6)
    tag(ax, "(b)")
    residual_axis(axr, res02, runs["02"][3] or lerr02)
    axr.set_xlabel(r"energy $E$ (eV)")
    axr.set_ylabel(r"new $-$ golden", fontsize=6.5)

    # (c) tutorial 05
    ax = fig.add_subplot(gs[0, 2])
    axr = fig.add_subplot(gs[1, 2], sharex=ax)
    cols = plt.cm.viridis(np.linspace(0.1, 0.85, len(T05["W"])))
    for i, W in enumerate(T05["W"]):
        tg, Dg5 = gold05[W]
        k = slice(None, None, 8)
        ax.plot(tg[k], Dg5[k], ls="none", marker=K[i], ms=3.2, mfc="none", mew=0.6, color=cols[i],
                label=(r"$W=%g$" % W) + ("" if new05 else " golden"))
        if new05:
            tn, Dn5 = new05[W]
            ax.plot(tn, Dn5, color=cols[i], lw=0.9)
            axr.plot(tn, Dn5 - Dg5, color=cols[i], lw=0.7, ls=L[i])
    if new05:
        ax.plot([], [], color="0.35", lw=0.9, label="new (lines)")
        ax.plot([], [], ls="none", marker="o", mfc="none", color="0.35", label="golden")
    ax.set_ylabel(r"$D(E,t)=\Delta X^2/2t$")
    ax.tick_params(labelbottom=False)
    ax.legend(loc="center right", bbox_to_anchor=(1.0, 0.42), ncol=2, fontsize=5, frameon=False,
              title=r"$E=2$ eV", title_fontsize=5.5, columnspacing=0.8)
    tag(ax, "(c)")
    residual_axis(axr, res05, runs["05"][3] or lerr05)
    axr.set_xlabel(r"time $t$ (fs)")
    axr.set_ylabel(r"new $-$ golden", fontsize=6.5)

    png = os.path.join(a.out, a.name + ".png")
    fig.savefig(png, dpi=300, bbox_inches="tight")

    def cmp_summary(num):
        lines, ok, n, nnew = cmps[num]
        bad = [l for l in lines if "DIFFERS" in l]
        s = "%d of %d IDENTICAL" % (ok, n)
        if runs[num][3]:
            s += " (app failed: %s)" % runs[num][3]
        elif bad:
            s += " (DIFFERS: " + "; ".join(bad[:6]) + ("; ..." if len(bad) > 6 else "") + ")"
        return s

    def wall(num):
        r = runs[num]
        return ("reused output of %s" % r[0]) if r[1] is None else ("%.1f s%s" % (r[1], " (failed)" if r[3] else ""))

    def fmt(x):
        return "n/a (no new data)" if x != x else "%.1e" % x

    all_ok = all(cmps[n][1] == cmps[n][2] for n in cmps)
    reused = ["tutorial %s from %s" % (k, v) for k, v in sorted(reuse.items())]
    caption = (
        "FIG. 4. Stage 4 (handoff 0004): tutorials 04, 02 and 05 reproduced from Python by "
        "python3 -m lsquant.apps.tutorial04|tutorial02|tutorial05 (lsquant package over the pybind11 "
        "module lsquant._core, OMP_NUM_THREADS = %d), each run in a temporary directory%s; time evolution "
        "U(dt) = exp(-i H dt/hbar) by the Bessel-Chebyshev series (GSL J_n via lsquant._core.compat, "
        "|J_n| <= 1e-15 cut, hbar = %.10g eV fs), dt = tmax/(steps - 1). "
        "(a) Tutorial 04, Larmor precession on chain1d_mag (N = 128, t = -1, J_ex = %.1f eV): m = 0 row of "
        "TimeProj{SX,SY,SZ}-PX (projector PX = (1 + S_x)/2, exact trace over the 128 site states, M = %d, "
        "%d time steps, tmax = %g fs, Gershgorin window [-2.31, 2.31] eV), normalised by <S_x(0)> as "
        "lsqspin.py does; new as filled markers, golden as open grey markers (every 4th step), "
        "cos(wt), -sin(wt), 0 with w = 2 J_ex/hbar as thin lines; residual new - analytic below "
        "(golden - analytic max %.2e, round-off level); max |new - golden| = %s. "
        "(b) Tutorial 02, localisation in the 1D Anderson chain (box disorder W = %g eV, disorder seed "
        "and KPM_SEED = r for r = 1..%d, one random-phase vector, window from BOUNDS, M = %d, %d steps, "
        "tmax = %g fs, MSD(E,t) at %d meV): xi(E) = cal sqrt(<MSD>) over the last 30%% of the time window, "
        "averaged over the %d realisations, N = 256, 1024, 4096, E = 0, +-1, 1.9 eV, cal fixed once "
        "against the Born length 24(4 - E^2)/W^2 at N = 4096, E = 1 eV (lsqloc.py); new as lines, golden "
        "as open markers, Born dashed; residual new - golden below (max %s sites). "
        "(c) Tutorial 05, 2D square lattice L = %d (N = 40000), disorder W = 0.5, 1, 2, 4 eV (seed 1, "
        "KPM_SEED = 1, one random-phase vector): D(E,t) = DeltaX^2/(2t) at E = %g eV from the MSD "
        "(M = %d, %d steps, tmax = %g fs, %d meV; the app also writes the DOS, M = %d, %d meV); new as "
        "lines, golden as open markers (every 8th step); residual new - golden below (max %s). "
        "Byte comparison with goldens/stage0/tutorials/0N (captured/NNNN_<name> against the capture "
        "with the same index and name, top-level .dat by name; cmp criterion): tutorial 04 %s; "
        "tutorial 02 %s; tutorial 05 %s. Wall times (TCN, load 60 to 100 on 64 cores): tutorial 04 %s, "
        "tutorial 02 %s, tutorial 05 %s. %s"
        % (a.threads, ("; " + ", ".join(reused) + " (diagnostic)") if reused else "", HBAR, J_EX,
           T04["M"], T04["NTIME"], T04["TMAX"], res04_ana, fmt(res04_gold),
           T02["W"], T02["NREAL"], T02["M"], T02["NTIME"], T02["TMAX"], T02["ETA"], T02["NREAL"], fmt(res02),
           T05["L"], T05["E"], T05["M"], T05["NTIME"], T05["TMAX"], T05["ETA"], T05["MDOS"], T05["ETA"],
           fmt(res05), cmp_summary("04"), cmp_summary("02"), cmp_summary("05"),
           wall("04"), wall("02"), wall("05"),
           "All files IDENTICAL." if all_ok else "AT LEAST ONE FILE DIFFERS OR IS NOT PRODUCED."))
    with open(os.path.join(a.out, a.name + ".caption.txt"), "w") as f:
        f.write(caption + "\n")
    with open(os.path.join(a.out, a.name + ".cmp.txt"), "w") as f:
        for num in ("04", "02", "05"):
            lines, ok, n, nnew = cmps[num]
            f.write("tutorial %s: %d of %d IDENTICAL, %d new captures, %s%s\n"
                    % (num, ok, n, nnew, wall(num), ("; " + runs[num][3]) if runs[num][3] else ""))
            for l in lines:
                f.write("  " + l + "\n")
    for num in ("04", "02", "05"):
        lines, ok, n, nnew = cmps[num]
        bad = [l for l in lines if "DIFFERS" in l]
        print("tutorial %s cmp: %d of %d IDENTICAL%s" % (num, ok, n, (", %d DIFFERS" % len(bad)) if bad else ""))
        for l in bad[:10]:
            print("  " + l)
    print("max |new - golden|: (a) %s  (b) %s  (c) %s; (a) max |golden - analytic| %.3e"
          % (fmt(res04_gold), fmt(res02), fmt(res05), res04_ana))
    print("wrote", png)
    print(caption)
    for num in runs:
        if num not in reuse:
            shutil.rmtree(runs[num][0], ignore_errors=True)
    return 0 if all_ok else 1


if __name__ == "__main__":
    sys.exit(main())
