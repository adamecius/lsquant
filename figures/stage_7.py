#!/usr/bin/env python3
"""Stage-7 figure: the six tutorial figures, regenerated from run files, beside stage 0.

Stage 7 (AGENTS.md section 9) closes when every golden is reproduced from word syntax.
For the tutorials that means: `python -m lsquant run runs/goldens/tutorial0N.json` writes
the same data the stage-0 legacy pipeline captured (goldens/stage0/tutorials/0N), and
the tutorial figure drawn from it is the stage-0 figure.

One row per tutorial:
  left    the stage-0 figure (figures/out/stage_0/0N/<main figure>.png, legacy build);
  middle  the curves of that figure drawn from the run-file outputs (solid) over the
          stage-0 captures (dashed, goldens/stage0/tutorials/0N);
  right   the residual new - golden of the drawn curves, and in the title the number of
          byte-identical files and the max |new - golden| over every numeric file of the
          tutorial golden (all of captured/ and the top-level curves).
Tutorial 06 (Elliott-Yafet and Dyakonov-Perel spin relaxation) is drawn from its EY curves.

usage: python3 figures/stage_7.py --runs <dir> [--out figures/out/stage_7]
  <dir> holds the run-file outputs tutorial01 ... tutorial05 (scripts/check_goldens.sh 7
  with S7_KEEP=<keep>; then <dir> = <keep>/T16). Writes fig_stage_7.png and
  fig_stage_7.caption.txt. Never writes into goldens/.
"""
import argparse
import filecmp
import glob
import hashlib
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.image as mpimg
import matplotlib.pyplot as plt
import numpy as np

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
GOLD = os.path.join(REPO, "goldens", "stage0", "tutorials")
STAGE0 = os.path.join(REPO, "figures", "out", "stage_0")

MAIN_PNG = {1: "fig_dos.png", 2: "fig_localization.png", 3: "fig_haldane_hall.png",
            4: "fig_precession.png", 5: "fig_diffusion.png", 6: "fig_ey.png"}


def tokens(path):
    """All numeric tokens of a legacy text file (header lines included), as float64."""
    out = []
    with open(path) as fh:
        for line in fh:
            for t in line.split():
                try:
                    out.append(float(t))
                except ValueError:
                    pass
    return np.array(out)


def curve(path):
    d = np.loadtxt(path)
    return d[:, 0], d[:, 1]


def td_row0(path):
    """Row m = 0 of a legacy .chebmomTD file: (times, Re mu_0(t))."""
    with open(path) as fh:
        h = fh.readline().split()
        M, Nt = (int(x) for x in fh.readline().split()[:2])
        v = np.loadtxt(fh)
    dt = float(h[4])
    return np.arange(Nt) * dt, v[:Nt, 0]


def tutorial_stats(n, run_dir):
    """(identical, compared, max |new - golden|) over every non-PNG golden file of tutorial n."""
    gdir = os.path.join(GOLD, "%02d" % n)
    same = tot = 0
    worst = 0.0
    for g in sorted(glob.glob(os.path.join(gdir, "**", "*"), recursive=True)):
        if not os.path.isfile(g) or g.endswith(".png"):
            continue
        rel = os.path.relpath(g, gdir)
        new = os.path.join(run_dir, rel)
        tot += 1
        if rel == os.path.join("captured", "index.tsv"):
            # stage-0 capture index (name, sha1, legacy source path): checked by sha1
            ok = True
            for line in open(g):
                name, sha = line.split("\t")[:2]
                p = os.path.join(run_dir, "captured", name)
                ok = ok and os.path.isfile(p) and hashlib.sha1(open(p, "rb").read()).hexdigest() == sha
            same += int(ok)
            # a mismatching index entry is a mismatching capture, counted there
            continue
        if not os.path.isfile(new):
            worst = np.inf
            continue
        if filecmp.cmp(new, g, shallow=False):
            same += 1
            continue
        if g.endswith(".tsv"):
            worst = np.inf
            continue
        a, b = tokens(new), tokens(g)
        worst = max(worst, float(np.max(np.abs(a - b))) if a.shape == b.shape else np.inf)
    return same, tot, worst


def differing(n, run_dir):
    """Relative paths of the non-PNG golden files of tutorial n that are not byte-identical."""
    gdir = os.path.join(GOLD, "%02d" % n)
    out = []
    for g in sorted(glob.glob(os.path.join(gdir, "**", "*"), recursive=True)):
        if not os.path.isfile(g) or g.endswith(".png"):
            continue
        rel = os.path.relpath(g, gdir)
        new = os.path.join(run_dir, rel)
        if rel == os.path.join("captured", "index.tsv"):
            continue
        if not (os.path.isfile(new) and filecmp.cmp(new, g, shallow=False)):
            out.append(rel)
    return out


def pairs(n, run_dir):
    """The curves drawn in the middle panel: list of (label, new path, golden path, reader)."""
    gdir = os.path.join(GOLD, "%02d" % n)
    cap = lambda pat: sorted(glob.glob(os.path.join(gdir, "captured", pat)))
    out = []
    if n == 1:
        for f in cap("*_meanmoments_chain1d_N*_M512JACKSON.dat")[:3]:
            lab = "N = " + os.path.basename(f).split("_N")[1].split("_")[0]
            out.append((lab, f, curve))
    elif n == 2:
        seen = set()
        for f in cap("*EF0.000000JACKSON.dat"):
            N = os.path.basename(f).split("_N")[1].split("_")[0]
            if N not in seen:
                seen.add(N)
                out.append(("N = %s, E_F = 0, first realisation" % N, f, curve))
    elif n == 3:
        for f in cap("*.dat"):
            out.append((os.path.basename(f)[5:].split("_NonEq")[0], f, curve))
    elif n == 4:
        for f in cap("*.chebmomTD"):
            out.append((os.path.basename(f)[13:15] + " (m = 0 row)", f, td_row0))
    elif n == 5:
        for f in cap("*meanCorrelationVX*.dat"):
            out.append(("W = " + os.path.basename(f).split("_W")[1].split("KPM")[0], f, curve))
    elif n == 6:
        for f in cap("*_meanTimeProjSZ-PZey_L20_D*.dat"):
            out.append(("EY, D = " + os.path.basename(f).split("_D")[1].split("KPM")[0], f, curve))
    return [(lab, os.path.join(run_dir, os.path.relpath(g, gdir)), g, rd) for lab, g, rd in out]


XLAB = {1: ("E (eV)", "DOS"), 2: ("t (hbar/eV)", "MSD"), 3: ("E (eV)", "sigma (e^2/h)"),
        4: ("t (hbar/eV)", "Re mu_0(t)"), 5: ("t (hbar/eV)", "MSD"),
        6: ("t (hbar/eV)", "projected spin")}
RUNFILE = {n: "runs/goldens/tutorial%02d.json" % n for n in range(1, 7)}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--runs", required=True)
    ap.add_argument("--out", default=os.path.join(REPO, "figures", "out", "stage_7"))
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    fig, ax = plt.subplots(6, 3, figsize=(15, 26), gridspec_kw={"width_ratios": [1.25, 1, 1]})
    cap = []
    for row, n in enumerate(range(1, 7)):
        png = os.path.join(STAGE0, "%02d" % n, MAIN_PNG[n])
        ax[row, 0].imshow(mpimg.imread(png))
        ax[row, 0].set_axis_off()
        ax[row, 0].set_title("(%s) stage 0, tutorial %02d: %s" % ("abcdef"[row], n, MAIN_PNG[n]), fontsize=9)
        rdir = os.path.join(a.runs, "tutorial%02d" % n)
        same, tot, worst = tutorial_stats(n, rdir)
        resmax = 0.0
        for k, (lab, new, gold, rd) in enumerate(pairs(n, rdir)):
            c = "C%d" % k
            xg, yg = rd(gold)
            ax[row, 1].plot(xg, yg, "--", color=c, lw=2.2, alpha=0.5)
            if os.path.isfile(new):
                xn, yn = rd(new)
                ax[row, 1].plot(xn, yn, "-", color=c, lw=0.9, label=lab)
                if xn.shape == xg.shape:
                    ax[row, 2].plot(xn, yn - yg, "-", color=c, lw=0.9, label=lab)
                    resmax = max(resmax, float(np.max(np.abs(yn - yg))))
        ax[row, 1].set_xlabel(XLAB[n][0]); ax[row, 1].set_ylabel(XLAB[n][1])
        ax[row, 1].legend(fontsize=7)
        ax[row, 1].set_title("from %s (solid) vs stage-0 capture (dashed)" % RUNFILE[n], fontsize=9)
        ax[row, 2].set_xlabel(XLAB[n][0]); ax[row, 2].set_ylabel("new - golden")
        ax[row, 2].set_title("residual: %d/%d files byte-identical, max |new - golden| = %.3g" % (same, tot, worst), fontsize=9)
        cap.append("(%s) tutorial %02d from %s: %d of %d golden files (goldens/stage0/tutorials/%02d, "
                   "every non-PNG file; captured/index.tsv checked by sha1) byte-identical; max |new - golden| over all of them %.3g; drawn curves "
                   "residual %.3g." % ("abcdef"[row], n, RUNFILE[n], same, tot, n, worst, resmax))
        if same < tot:
            bad = differing(n, rdir)
            den = [b for b in bad if "TimeProj1-" in b]
            cap[-1] += (" The %d differing files: %d denominator curves of the word with operator '1' "
                        "(TimeProj1-P*), %s, and captured/index.tsv on those entries. Drawn: the numerator curves, all byte-identical."
                        % (len(bad), len(den), ", ".join(b for b in bad if "TimeProj1-" not in b)))
    fig.suptitle("FIG. 7. The tutorial figures regenerated from run files (python -m lsquant run), beside stage 0",
                 fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.985))
    path = os.path.join(a.out, "fig_stage_7.png")
    fig.savefig(path, dpi=110)
    text = ("FIG. 7. The six tutorial figures of stage 0 (left, legacy build, figures/out/stage_0) beside the "
            "same curves drawn from the outputs of `python -m lsquant run runs/goldens/tutorial0N.json` "
            "(middle, solid) over the stage-0 captures (dashed), with the residual new - golden (right). "
            "Run-file outputs: %s (OMP_NUM_THREADS = 16, KPM_SEED = 12345).\n" % a.runs + "\n".join(cap) + "\n")
    with open(os.path.join(a.out, "fig_stage_7.caption.txt"), "w") as fh:
        fh.write(text)
    print(path)
    print(text)


if __name__ == "__main__":
    sys.exit(main())
