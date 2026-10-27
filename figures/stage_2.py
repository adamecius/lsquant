#!/usr/bin/env python3
"""Stage-2 figure (handoff 0002, codex.md section 3): tutorial 01 reproduced from Python.

Runs `python3 -m lsquant.apps.tutorial01 --out <tmp>` (the Python library over the
pybind11 module lsquant._core, PYTHONPATH=<build>/lsquant_core:python) in a temporary
directory. The app writes the M = 512 exact-trace moment files with their legacy names at
the top of <tmp> and the reconstructed Jackson curves under <tmp>/eta<meV>/ with their
legacy names. They are matched to goldens/stage0/tutorials/01/captured/ by the tutorial's
order of operations (the captured prefixes follow examples/01_dos_1d_chain/lsqplot.py):

  moment files     by name after the 4-digit prefix (every captured copy of that name)
  eta20/  N = 256, 1024, 4096    ->  0013_, 0015_, 0017_
  eta80/, eta24/, eta10/  N = 256  ->  0019_, 0020_, 0021_

Every pair is byte-compared (filecmp, shallow=False: the cmp criterion) and one figure is
drawn, with the per-site DOS rho(E) = value / N (the normalisation of lsqplot.py):

  (a) rho(E) for N = 256, 1024, 4096 at 20 meV: new as lines, golden as markers, the
      analytic 1/(pi sqrt(4 - E^2)) as a thin reference; residual new - golden below.
  (b) N = 256 at 80, 24, 10 meV, same encoding.

Outputs figures/out/stage_2/fig_stage_2.png and fig_stage_2.caption.txt.
Usage (repository root, on TCN):
    MPLBACKEND=Agg python3 figures/stage_2.py [--build build] [--threads 16]
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

import numpy as np

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
CAPT = os.path.join(REPO, "goldens", "stage0", "tutorials", "01", "captured")
SIZES, M = (256, 1024, 4096), 512
FIXED_ETA, ETAS = 20, (80, 24, 10)
# (eta in meV, N) -> captured prefix, in the order lsqplot.py produced them
CURVES = {(20, 256): "0013", (20, 1024): "0015", (20, 4096): "0017",
          (80, 256): "0019", (24, 256): "0020", (10, 256): "0021"}


def find_pydir(build):
    for sub in ("lsquant_core", "lsquant"):
        d = os.path.join(build, sub)
        if glob.glob(os.path.join(d, "_core*.so")):
            return os.path.abspath(d)
    sys.exit("_core*.so not found under %s/{lsquant_core,lsquant}" % build)


def curve_name(N):
    return "meanmoments_chain1d_N%d_M%dJACKSON.dat" % (N, M)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--build", default=os.path.join(REPO, "build"))
    ap.add_argument("--threads", type=int, default=16)
    ap.add_argument("--out", default=os.path.join(REPO, "figures", "out", "stage_2"))
    a = ap.parse_args()
    pydir = find_pydir(a.build)
    os.makedirs(a.out, exist_ok=True)

    tmp = tempfile.mkdtemp(prefix="stage2_tut01_")
    env = dict(os.environ, OMP_NUM_THREADS=str(a.threads), KPM_SEED="12345",
               PYTHONPATH=pydir + os.pathsep + os.path.join(REPO, "python"))
    r = subprocess.run([sys.executable, "-m", "lsquant.apps.tutorial01", "--out", tmp],
                       cwd=REPO, env=env, capture_output=True, text=True)
    printed = [l.replace(tmp, "<tmp>") for l in r.stdout.splitlines()]
    if r.returncode != 0:
        sys.exit("tutorial01 exited %d: %s" % (r.returncode, r.stderr.strip()))

    # ---- byte comparison
    cmp_lines = []

    def cmp(new, gold):
        ok = os.path.isfile(new) and filecmp.cmp(new, gold, shallow=False)
        cmp_lines.append("%s vs captured/%s: %s" % (os.path.relpath(new, tmp), os.path.basename(gold),
                                                     "IDENTICAL" if ok else "DIFFERS"))
        return ok

    captured = sorted(f for f in os.listdir(CAPT) if re.match(r"^\d{4}_", f))
    for f in captured:
        name = f[5:]
        if name.endswith(".chebmom1D"):
            cmp(os.path.join(tmp, name), os.path.join(CAPT, f))
    for (eta, N), pre in sorted(CURVES.items(), key=lambda kv: kv[1]):
        cmp(os.path.join(tmp, "eta%d" % eta, curve_name(N)), os.path.join(CAPT, pre + "_" + curve_name(N)))
    n_ident = sum(l.endswith("IDENTICAL") for l in cmp_lines)
    all_identical = n_ident == len(cmp_lines)

    def load(eta, N):
        new = np.loadtxt(os.path.join(tmp, "eta%d" % eta, curve_name(N)))
        gold = np.loadtxt(os.path.join(CAPT, CURVES[(eta, N)] + "_" + curve_name(N)))
        return new[:, 0], new[:, 1] / N, gold[:, 0], gold[:, 1] / N

    data_a = {N: load(FIXED_ETA, N) for N in SIZES}
    data_b = {eta: load(eta, SIZES[0]) for eta in ETAS}
    max_res, max_dE, npts = 0.0, 0.0, {}
    for k, (E, y, Eg, yg) in list(data_a.items()) + [(("b", e), v) for e, v in data_b.items()]:
        if E.shape != Eg.shape:
            sys.exit("line count differs for %s: %d new vs %d golden" % (k, E.size, Eg.size))
        max_res = max(max_res, float(np.max(np.abs(y - yg))))
        max_dE = max(max_dE, float(np.max(np.abs(E - Eg))))
        npts[k] = E.size

    # ---- figure
    import matplotlib
    matplotlib.use(os.environ.get("MPLBACKEND", "Agg"))
    import matplotlib.pyplot as plt
    sys.path.insert(0, os.path.join(REPO, "examples"))
    try:
        import prl_style as prl
        prl.use("double", aspect=0.5, usetex=False)
        C, L, K = prl.COLORS, prl.LINES, prl.MARKERS
        tag = prl.panel
    except ImportError:              # conventions inlined if the style module is absent
        C, L, K = ["#000000", "#0072B2", "#D55E00", "#009E73"], ["-", "--", "-.", ":"], ["o", "s", "^", "v"]

        def tag(ax, t, loc="upper left", **kw):
            ax.text(0.04, 0.96, t, transform=ax.transAxes, ha="left", va="top", fontweight="bold", **kw)

    fig = plt.figure(figsize=(7.0, 3.6))
    gs = fig.add_gridspec(2, 2, height_ratios=[3.0, 1.0], hspace=0.12, wspace=0.28)
    Ea = np.linspace(-1.995, 1.995, 2000)
    rho_exact = 1.0 / (np.pi * np.sqrt(4.0 - Ea ** 2))
    nmark = 48                        # golden markers: 48 evenly spaced points per curve (curves have 2370 to 15360)

    def panel(col, series, label_fmt, title, letter):
        ax = fig.add_subplot(gs[0, col])
        axr = fig.add_subplot(gs[1, col], sharex=ax)
        ax.plot(Ea, rho_exact, color="0.55", lw=0.6, zorder=1,
                label=r"$1/\pi\sqrt{4-E^2}$")
        for i, (key, (E, y, Eg, yg)) in enumerate(series):
            c = C[(i + 1) % len(C)]
            ax.plot(E, y, color=c, lw=0.6, linestyle=L[0], zorder=10 - 2 * i, label=label_fmt % key + " new")
            k = np.linspace(0, Eg.size - 1, nmark).astype(int)
            ax.plot(Eg[k], yg[k], linestyle="none", marker=K[i % len(K)], ms=2.8,
                    mfc="none", mew=0.6, color=c, zorder=11 - 2 * i, label=label_fmt % key + " golden")
            axr.plot(E, y - yg, color=c, lw=0.8, linestyle=L[i % len(L)])
        ax.set_ylim(0.0, 0.95)
        ax.set_xlim(-2.05, 2.05)
        ax.set_ylabel(r"$\rho(E)\;(1/\mathrm{eV})$")
        ax.tick_params(labelbottom=False)
        ax.legend(loc="upper center", ncol=2, fontsize=5.5, title=title, title_fontsize=6,
                  handlelength=1.6, columnspacing=0.8, frameon=False)
        tag(ax, letter)
        axr.axhline(0.0, color="0.7", lw=0.4)
        if max(float(np.max(np.abs(y - yg))) for _, (E, y, Eg, yg) in series) == 0.0:
            axr.set_ylim(-1e-15, 1e-15)
            axr.text(0.5, 0.75, "identically 0", transform=axr.transAxes, ha="center", va="center",
                     fontsize=6)
        axr.set_xlabel(r"energy $E$ (eV)")
        axr.set_ylabel(r"new $-$ golden", fontsize=7)
        axr.tick_params(axis="y", labelsize=6)
        axr.yaxis.get_offset_text().set_fontsize(6)

    panel(0, [(N, data_a[N]) for N in SIZES], r"$N=%d$", r"$\eta=%d$ meV, $M=%d$" % (FIXED_ETA, M), "(a)")
    panel(1, [(e, data_b[e]) for e in ETAS], r"$\eta=%d$ meV", r"$N=%d$, $M=%d$" % (SIZES[0], M), "(b)")

    png = os.path.join(a.out, "fig_stage_2.png")
    fig.savefig(png, dpi=300, bbox_inches="tight")

    def m_eff(eta):                   # legacy: M_eff = ceil(pi/eta_r), eta_r = 2 eta / 1000 / BandWidth (4 eV)
        return min(M, int(np.ceil(np.pi / (2.0 * eta / 1000.0 / 4.0))))

    caption = (
        "FIG. 2. Stage 2 (handoff 0002): tutorial 01, the density of states of the periodic "
        "nearest-neighbour chain (t = -1 eV), reproduced from Python by "
        "python3 -m lsquant.apps.tutorial01 (lsquant package over the pybind11 module lsquant._core, "
        "OMP_NUM_THREADS = %d): N = 256, 1024, 4096 sites, window [-2, 2] eV, exact trace over all N "
        "site states, M = %d Chebyshev moments, Jackson kernel, reconstruction on 30*M_eff points over "
        "|x| <= 0.95, per-site DOS rho(E) = value / N. "
        "(a) N = 256, 1024, 4096 at eta = %d meV (M_eff = %d): new as lines, the stage-0 golden "
        "(goldens/stage0/tutorials/01/captured 0013, 0015, 0017) as open markers (%d evenly spaced points per curve), "
        "the analytic bulk DOS 1/(pi sqrt(4 - E^2)) as the thin grey line; residual new - golden below. "
        "(b) N = 256 at eta = 80, 24, 10 meV (M_eff = %d, %d, %d; captured 0019, 0020, 0021), same "
        "encoding. Max |new - golden| over all six curves: %.1e per site, max energy-axis difference "
        "%.1e eV. Byte comparison with the captured goldens (cmp criterion), %d of %d IDENTICAL: %s. %s"
        % (a.threads, M, FIXED_ETA, m_eff(FIXED_ETA), nmark, m_eff(80), m_eff(24), m_eff(10),
           max_res, max_dE, n_ident, len(cmp_lines), "; ".join(cmp_lines),
           "All files IDENTICAL." if all_identical else "AT LEAST ONE FILE DIFFERS."))
    with open(os.path.join(a.out, "fig_stage_2.caption.txt"), "w") as f:
        f.write(caption + "\n")
    print("tutorial01 printed (%d lines):" % len(printed))
    for l in printed:
        print("  " + l)
    print("cmp list:")
    for l in cmp_lines:
        print("  " + l)
    print("wrote", png)
    print(caption)
    shutil.rmtree(tmp)
    return 0 if all_identical else 1


if __name__ == "__main__":
    sys.exit(main())
