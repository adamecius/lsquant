#!/usr/bin/env python3
"""Stage-6 figure (handoff 0006, codex.md section 3): the rotated (non-orthogonal) models
reproduce the orthogonal goldens through both metric routes.

Reads the runs kept by the gate (`S6_KEEP=figures/out/stage_6/data scripts/check_goldens.sh 6
<build>`): figures/out/stage_6/data/{runs/,gate_table.tsv,chern.tsv}. Nothing is recomputed
here except the stage-2 DOS reconstruction of the moment files.

  (a) chain1d DOS (N = 2048, window [-2, 2], exact trace, M = 64), reconstructed with the
      stage-2 lsquant.reconstruct.dos_curve (Jackson, 120 meV as G3): orthogonal golden G1,
      Lowdin and generalised routes with the exact applier on the rotated `complex` fixture;
      residual (route - golden) below for both fixtures (`real`, `complex`).
  (b) Haldane sigma_xy(E) (N = 512, window [-4.5, 4.5], exact trace, M = 128, Bastin 10 meV
      from python3 -m lsquant.apps.stage3_kubo): golden G5 Bastin curve and the two routes
      (exact applier, `complex` fixture), in units of e^2/h per (N/2) A_cell, so the gap
      plateau is the Chern number; residual below.
The caption carries C for every route and applier and the max relative moment difference
per applier (gate_table.tsv).

Outputs figures/out/stage_6/fig_stage_6.png and fig_stage_6.caption.txt.
Usage (repository root, on TCN):
    MPLBACKEND=Agg python3 figures/stage_6.py [--build build] [--data figures/out/stage_6/data]
"""
import argparse
import csv
import glob
import math
import os
import sys

import numpy as np

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
GOLD = os.path.join(REPO, "goldens", "stage0")
G1 = "SpectralOp1chain1dKPM_M64_stateexact.chebmom1D"
XY = "NonEqOpVX-VYhaldaneKPM_M128x128_stateexact"
BASTIN = "KuboBastin_%sJACKSON.dat" % XY
A_CELL = 3.0 * math.sqrt(3.0) / 2.0
ETA_DOS = 120.0


def find_pydir(build):
    for sub in ("lsquant_core", "lsquant"):
        d = os.path.join(build, sub)
        if glob.glob(os.path.join(d, "_core*.so")):
            return os.path.abspath(d)
    sys.exit("_core*.so not found under %s/{lsquant_core,lsquant}" % build)


def chern(dim, E, s):               # test/hall_haldane.sh, tutorial03.chern_from_plateau
    return (np.mean(s[np.abs(E) < 0.5]) / ((dim / 2) * A_CELL)) * 2 * np.pi


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--build", default=os.path.join(REPO, "build"))
    ap.add_argument("--data", default=os.path.join(REPO, "figures", "out", "stage_6", "data"))
    ap.add_argument("--out", default=os.path.join(REPO, "figures", "out", "stage_6"))
    a = ap.parse_args()
    sys.path.insert(0, find_pydir(a.build))
    sys.path.insert(0, os.path.join(REPO, "python"))
    from lsquant import _core
    from lsquant.reconstruct import dos_curve
    runs = os.path.join(a.data, "runs")

    # ---- (a) chain1d DOS
    def dos(path):
        E, y = dos_curve(_core.ReadChebmom(path), ETA_DOS)
        return np.asarray(E), np.asarray(y) / 2048.0      # per site
    Eg, yg = dos(os.path.join(GOLD, "chain1d", G1))
    dosd, dosres = {}, {}
    for v in ("real", "complex"):
        for r in ("lowdin", "generalised"):
            E, y = dos(os.path.join(runs, "chain1d_%s_dos_%s_exact" % (v, r), G1))
            if not np.array_equal(E, Eg):
                sys.exit("energy grid differs for chain1d %s %s" % (v, r))
            dosd[(v, r)] = y
            dosres[(v, r)] = y - yg

    # ---- (b) Haldane sigma_xy
    norm = (512 / 2) * A_CELL / (2 * math.pi)
    gx = np.loadtxt(os.path.join(GOLD, "haldane", BASTIN))
    hal, halres = {}, {}
    for v in ("real", "complex"):
        for r in ("lowdin", "generalised"):
            d = np.loadtxt(os.path.join(runs, "haldane_%s_kuboxy_%s_exact" % (v, r), BASTIN))
            if d.shape != gx.shape or float(np.max(np.abs(d[:, 0] - gx[:, 0]))) > 0:
                sys.exit("energy grid differs for haldane %s %s" % (v, r))
            hal[(v, r)] = d[:, 1] / norm
            halres[(v, r)] = (d[:, 1] - gx[:, 1]) / norm
    C_gold = chern(512, gx[:, 0], gx[:, 1])

    # ---- tables from the gate
    with open(os.path.join(a.data, "gate_table.tsv")) as f:
        rows = list(csv.DictReader(f, delimiter="\t"))
    with open(os.path.join(a.data, "chern.tsv")) as f:
        crow = list(csv.DictReader(f, delimiter="\t"))
    per_app = {}
    for rw in rows:
        if rw["route"] not in ("lowdin", "generalised") or rw["max_rel"] in ("N/A", "nan"):
            continue
        k = (rw["route"], rw["applier"])
        per_app[k] = max(per_app.get(k, 0.0), float(rw["max_rel"]))
    cross = {}
    for rw in rows:
        if rw["route"].startswith("generalised") and "vs" in rw["route"] and rw["max_rel"] != "nan":
            cross[rw["applier"]] = max(cross.get(rw["applier"], 0.0), float(rw["max_rel"]))
    naive = [rw for rw in rows if rw["route"] == "naive"]
    n_pass = sum(rw["verdict"] == "PASS" for rw in rows)
    n_fail = sum(rw["verdict"] == "FAIL" for rw in rows)

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
    except ImportError:
        CO, L, K = ["#000000", "#0072B2", "#D55E00", "#009E73"], ["-", "--", "-.", ":"], ["o", "s", "^", "v"]

        def tag(ax, t, loc="upper left", **kw):
            ax.text(0.04, 0.96, t, transform=ax.transAxes, ha="left", va="top", fontweight="bold", **kw)

    fig = plt.figure(figsize=(3.6, 6.4))
    gs = fig.add_gridspec(5, 1, height_ratios=[3.0, 1.1, 0.55, 3.0, 1.1], hspace=0.12)
    nmark = 40
    lab = {"lowdin": "Löwdin", "generalised": "generalised"}
    sty = {("real", "lowdin"): (1, "-"), ("real", "generalised"): (2, "--"),
           ("complex", "lowdin"): (1, ":"), ("complex", "generalised"): (2, "-.")}

    ax = fig.add_subplot(gs[0])
    axr = fig.add_subplot(gs[1], sharex=ax)
    k = np.linspace(0, Eg.size - 1, nmark).astype(int)
    ax.plot(Eg[k], yg[k], ls="none", marker=K[0], ms=2.8, mfc="none", mew=0.6, color=CO[0],
            label="orthogonal golden G1", zorder=5)
    ax.plot(Eg, dosd[("complex", "lowdin")], color=CO[1], lw=0.9, ls="-", label="Löwdin, exact")
    ax.plot(Eg, dosd[("complex", "generalised")], color=CO[2], lw=0.9, ls="--",
            label="generalised, exact")
    ax.set_ylabel(r"DOS per site (1/eV)")
    ax.set_xlim(-2.0, 2.0)
    ax.set_ylim(0, 1.1 * float(np.max(yg)))
    ax.tick_params(labelbottom=False)
    ax.legend(loc="upper center", fontsize=5.5, frameon=False,
              title="chain1d, rotated (complex T), N = 2048, M = 64", title_fontsize=5.5)
    tag(ax, "(a)", loc="upper left")
    for (v, r), y in dosres.items():
        ci, ls = sty[(v, r)]
        axr.plot(Eg, y, color=CO[ci], lw=0.7, ls=ls, label="%s %s" % (lab[r], v))
    axr.axhline(0, color="0.7", lw=0.4)
    axr.set_ylabel("route $-$ G1", fontsize=6)
    axr.set_xlabel(r"energy $E$ (eV)")
    axr.tick_params(axis="y", labelsize=5.5)
    axr.yaxis.get_offset_text().set_fontsize(5.5)
    axr.legend(fontsize=4.5, ncol=2, frameon=False, loc="lower center")

    bx = fig.add_subplot(gs[3])
    bxr = fig.add_subplot(gs[4], sharex=bx)
    E = gx[:, 0]
    k = np.linspace(0, E.size - 1, nmark).astype(int)
    C = {(v, r): chern(512, E, hal[(v, r)] * norm) for v in ("real", "complex") for r in ("lowdin", "generalised")}
    bx.axvspan(-0.5, 0.5, color="0.92", lw=0)
    bx.axhline(1.0, color="0.6", ls=":", lw=0.6)
    bx.plot(E[k], gx[k, 1] / norm, ls="none", marker=K[0], ms=2.8, mfc="none", mew=0.6, color=CO[0],
            label="golden G5\n$C=%+.4f$" % C_gold, zorder=5)
    bx.plot(E, hal[("complex", "lowdin")], color=CO[1], lw=0.9, ls="-",
            label="Löwdin, exact\n$C=%+.4f$" % C[("complex", "lowdin")])
    bx.plot(E, hal[("complex", "generalised")], color=CO[2], lw=0.9, ls="--",
            label="generalised, exact\n$C=%+.4f$" % C[("complex", "generalised")])
    bx.set_xlim(-4.4, 4.4)
    bx.set_ylabel(r"$\sigma_{xy}$ ($e^2/h$ per $(N/2)A_{\rm cell}$)", fontsize=6.5)
    bx.tick_params(labelbottom=False)
    bx.legend(loc="upper right", bbox_to_anchor=(1.0, 0.8), fontsize=4.6, frameon=False,
              title="Haldane, rotated (complex T)\nN = 512, M = 128, 10 meV", title_fontsize=5.0)
    tag(bx, "(b)", loc="upper left")
    for (v, r), y in halres.items():
        ci, ls = sty[(v, r)]
        bxr.plot(E, y, color=CO[ci], lw=0.7, ls=ls, label="%s %s" % (lab[r], v))
    bxr.axhline(0, color="0.7", lw=0.4)
    bxr.set_ylabel("route $-$ G5", fontsize=6)
    bxr.set_xlabel(r"energy $E$ (eV)")
    bxr.tick_params(axis="y", labelsize=5.5)
    bxr.yaxis.get_offset_text().set_fontsize(5.5)
    bxr.legend(fontsize=4.5, ncol=2, frameon=False, loc="upper center")

    os.makedirs(a.out, exist_ok=True)
    png = os.path.join(a.out, "fig_stage_6.png")
    fig.savefig(png, dpi=300, bbox_inches="tight")

    ctab = "; ".join("%s C = %s (%s)" % (c["source"], c["C"], c["verdict"]) for c in crow)
    atab = "; ".join("%s/%s %.2e" % (r, ap_, per_app[(r, ap_)]) for (r, ap_) in sorted(per_app))
    xtab = "; ".join("%s %.2e" % (ap_, cross[ap_]) for ap_ in sorted(cross))
    ntab = "; ".join("%s/%s %s %s" % (rw["model"], rw["variant"], rw["word"], rw["max_rel"]) for rw in naive)
    mdos = max(float(np.max(np.abs(y))) for y in dosres.values())
    mhal = max(float(np.max(np.abs(y))) for y in halres.values())
    caption = (
        "FIG. 6. Stage 6 (handoff 0006): a non-orthogonal basis through MetricPolicy reproduces the "
        "orthogonal goldens. Fixtures: the committed chain1d and Haldane operators rotated by a local "
        "transform T (H' = T^+ H T, S' = T^+ T, V' = T^+ V T; cond(S') 5.2 to 6.7; variants real and "
        "complex T), goldens/stage0/nonorth. Moments from build/lsquant_core/stage6_nonorth (exact trace, "
        "OMP_NUM_THREADS = 16). (a) chain1d DOS per site, window [-2, 2] eV, M = 64, reconstructed "
        "with the stage-2 lsquant.reconstruct.dos_curve (Jackson, %g meV): orthogonal golden G1 "
        "(open markers, %d points), Lowdin route (S^-1/2 H' S^-1/2) and generalised route "
        "(S^-1 (a H' + b S')) with the exact applier on the complex fixture; residual route - G1 "
        "below for both fixtures, max |residual| = %.1e 1/eV. (b) Haldane sigma_xy(E), "
        "window [-4.5, 4.5] eV, M = 128, Kubo-Bastin at 10 meV by python3 -m lsquant.apps.stage3_kubo, "
        "in units of e^2/h per (N/2) A_cell (A_cell = 3 sqrt(3)/2) so the gap plateau (shaded, "
        "|E| < 0.5 eV) is C; golden G5 as markers; residual below, max |residual| = %.1e. "
        "Chern numbers C = (mean_{|E|<0.5} sigma_xy / A) 2 pi: %s. "
        "Max relative moment difference vs G1, G2, G5 (max_k |mu - mu_gold| / max_k |mu_gold|, over "
        "models, variants and words) per route/applier: %s (tolerances: exact 1e-10, ischo 1e-9, "
        "chebyshev 1e-8; Lowdin/ischo not defined). Generalised vs Lowdin, same applier: %s. "
        "Naive route (legacy estimator, evidence only): %s. Gate rows: %d PASS, %d FAIL."
        % (ETA_DOS, nmark, mdos, mhal, ctab, atab, xtab, ntab, n_pass, n_fail))
    with open(os.path.join(a.out, "fig_stage_6.caption.txt"), "w") as f:
        f.write(caption + "\n")
    print("wrote", png)
    print(caption)
    return 0


if __name__ == "__main__":
    sys.exit(main())
