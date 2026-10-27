#!/usr/bin/env python3
"""Stage-1 figure (handoff 0001, codex.md section 4): the new core against G1 and G2.

Runs build/stage1_moments (the gate executable of brief.md section 6) on the committed
golden operators goldens/stage0/chain1d/operators at every requested OpenMP thread count,
into a temporary directory, with the gate arguments

    --label chain1d --window -2 2 --trace exact --states 2048 --moments 64 --word dos|kubo

and KPM_SEED=12345. It byte-compares (filecmp, shallow=False: the same criterion as cmp)
every legacy-layout output with goldens/stage0/chain1d, reads the new tensors from the
rank-generic .chebmom files written alongside, and draws one figure:

  (a) mu_m, m = 0..63 (G1): golden as markers, new as a line, delta_m0 as the analytic
      reference; the residual new - golden on a secondary axis.
  (b) the m, n <= 10 block of mu_mn (G2): |new - golden| as a heat map.
  (c) beside it, the closed-form B.15 block: Mref(m, n) of test/chain1d_analytic_test.cpp
      times the diagonal plateau mean(mu_mm, m = 2..9) of the new tensor.

Outputs figures/out/stage_1/fig_stage_1.png and fig_stage_1.caption.txt.
Usage (repository root, on TCN):
    MPLBACKEND=Agg python3 figures/stage_1.py [--build build] [--threads 1 8 64]
Nothing under goldens/ is written. Seeds and goldens are never changed here.
"""
import argparse
import filecmp
import os
import subprocess
import sys
import tempfile

import numpy as np

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
GOLD = os.path.join(REPO, "goldens", "stage0", "chain1d")
F1 = "SpectralOp1chain1dKPM_M64_stateexact"
F2 = "NonEqOpVX-VXchain1dKPM_M64x64_stateexact"
N_SITES, M, WINDOW, B = 2048, 64, (-2.0, 2.0), 10


def find_exe(build):
    for sub in ("lsquant_core", "lsquant", "."):
        p = os.path.join(build, sub, "stage1_moments")
        if os.access(p, os.X_OK):
            return os.path.abspath(p)
    sys.exit("stage1_moments not found under %s/{lsquant_core,lsquant,.}" % build)


def run(exe, threads, out):
    env = dict(os.environ, OMP_NUM_THREADS=str(threads), KPM_SEED="12345")
    for word in ("dos", "kubo"):
        subprocess.run([exe, "--operators", os.path.join(GOLD, "operators"), "--label", "chain1d",
                        "--window", "-2", "2", "--trace", "exact", "--states", str(N_SITES),
                        "--moments", str(M), "--word", word, "--out", out],
                       env=env, check=True, stdout=subprocess.DEVNULL)


def read_legacy(path):
    """Legacy .chebmom1D (header, M, pairs) or .chebmom2D (header, M N, pairs)."""
    with open(path) as f:
        f.readline()
        dims = [int(x) for x in f.readline().split()]
        v = np.loadtxt(f)
    return (v[:, 0] + 1j * v[:, 1]).reshape(dims)


def read_chebmom(path):
    """New rank-generic .chebmom (brief.md section 7): key lines, 'data', then re im."""
    meta = {}
    with open(path) as f:
        for line in f:
            if line.strip() == "data":
                break
            k, _, rest = line.strip().partition(" ")
            meta[k] = rest
        v = np.loadtxt(f, ndmin=2)
    dims = [int(x) for x in meta["dims"].split()]
    return (v[:, 0] + 1j * v[:, 1]).reshape(dims), meta


def mref(m, n):                     # test/chain1d_analytic_test.cpp, oracle kg_moment
    dmn, dm0, dm1, d2 = m == n, m == 0, m == 1, abs(m - n) == 2
    return (1.0 / (1 + dm0 + dm1) if dmn else 0.0) - 0.5 * (not dmn) * d2


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--build", default=os.path.join(REPO, "build"))
    ap.add_argument("--threads", type=int, nargs="+", default=[1, 8, 64])
    ap.add_argument("--out", default=os.path.join(REPO, "figures", "out", "stage_1"))
    a = ap.parse_args()
    exe = find_exe(a.build)
    os.makedirs(a.out, exist_ok=True)

    tmp = tempfile.mkdtemp(prefix="stage1_fig_")
    cmp_lines, tensors = [], {}
    for T in a.threads:
        d = os.path.join(tmp, "T%d" % T)
        os.makedirs(d)
        run(exe, T, d)
        for stem, ext in ((F1, ".chebmom1D"), (F2, ".chebmom2D")):
            same = filecmp.cmp(os.path.join(d, stem + ext), os.path.join(GOLD, stem + ext), shallow=False)
            cmp_lines.append("T=%d %s%s: %s" % (T, stem, ext, "IDENTICAL" if same else "DIFFERS"))
        tensors[T] = (read_chebmom(os.path.join(d, F1 + ".chebmom"))[0],
                      read_chebmom(os.path.join(d, F2 + ".chebmom"))[0])
    T0 = a.threads[0]
    mu1, mu2 = tensors[T0]
    threads_equal = all(np.array_equal(tensors[T][0], mu1) and np.array_equal(tensors[T][1], mu2)
                        for T in a.threads)
    g1 = read_legacy(os.path.join(GOLD, F1 + ".chebmom1D"))
    g2 = read_legacy(os.path.join(GOLD, F2 + ".chebmom2D"))
    # residuals: the new legacy-layout files (15 digits, same text as the goldens) minus the
    # goldens; the 17-digit .chebmom is checked against them separately (format consistency)
    dT0 = os.path.join(tmp, "T%d" % T0)
    n1 = read_legacy(os.path.join(dT0, F1 + ".chebmom1D"))
    n2 = read_legacy(os.path.join(dT0, F2 + ".chebmom2D"))
    res1 = n1 - g1
    res2 = n2 - g2
    fmt_dev = max(float(np.max(np.abs(mu1 - n1))), float(np.max(np.abs(mu2 - n2))))
    blk = slice(0, B + 1)
    plateau = float(np.mean([mu2[m, m].real for m in range(2, B)]))
    ref = np.array([[mref(m, n) for n in range(B + 1)] for m in range(B + 1)])
    b15_res = float(np.max(np.abs(mu2[:B, :B].real / plateau - ref[:B, :B])))
    max_im = float(np.max(np.abs(mu2[blk, blk].imag)))
    all_identical = all(l.endswith("IDENTICAL") for l in cmp_lines)

    import matplotlib
    matplotlib.use(os.environ.get("MPLBACKEND", "Agg"))
    import matplotlib.pyplot as plt
    sys.path.insert(0, os.path.join(REPO, "examples"))
    try:
        import prl_style as prl
        prl.use("double", aspect=0.36, usetex=False)
        C, L, K = prl.COLORS, prl.LINES, prl.MARKERS
        tag = prl.panel
    except ImportError:              # conventions inlined if the style module is absent
        C, L, K = ["#000000", "#0072B2", "#D55E00"], ["-", "--", "-."], ["o", "s", "^"]

        def tag(ax, t, loc="upper left", **kw):
            up = "upper" in loc
            ax.text(0.04, 0.96 if up else 0.04, t, transform=ax.transAxes, ha="left",
                    va="top" if up else "bottom", fontweight="bold", **kw)

    fig = plt.figure(figsize=(7.0, 2.6))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.3, 1.0, 1.0], wspace=1.0)
    ax = fig.add_subplot(gs[0])
    m = np.arange(M)
    ax.plot(m, g1.real, linestyle="none", marker=K[0], color=C[0], ms=4.5, mfc="none",
            label=r"golden G1")
    ax.plot(m, mu1.real, linestyle=L[1], color=C[1], label=r"new ($T=%s$)" % ",".join(map(str, a.threads)))
    ax.plot(m, (m == 0).astype(float), linestyle="none", marker=K[2], color=C[2], ms=2.5,
            label=r"$\delta_{m0}$")
    ax.set_xlabel(r"$m$")
    ax.set_ylabel(r"$\mu_m$")
    ax.set_ylim(-0.1, 1.25)
    ax.legend(loc="center right", fontsize=7)
    tag(ax, "(a)")
    axr = ax.twinx()
    axr.plot(m, np.abs(res1), linestyle=L[2], color=C[3], lw=0.9)
    axr.set_ylim(-1e-16, 1e-15)
    axr.set_ylabel(r"$|\mu_m^{\rm new}-\mu_m^{\rm golden}|$", color=C[3], fontsize=7)
    axr.tick_params(axis="y", colors=C[3], labelsize=6)
    axr.yaxis.get_offset_text().set_fontsize(6)

    axb = fig.add_subplot(gs[1])
    d = np.abs(res2[blk, blk])
    im = axb.imshow(d, origin="upper", cmap="Greys", vmin=0.0, vmax=max(1e-15, d.max()),
                    extent=(-0.5, B + 0.5, B + 0.5, -0.5))
    cb = fig.colorbar(im, ax=axb, fraction=0.046, pad=0.04)
    cb.set_label(r"$|\mu_{mn}^{\rm new}-\mu_{mn}^{\rm golden}|$", fontsize=7)
    cb.ax.tick_params(labelsize=6)
    cb.ax.yaxis.get_offset_text().set_fontsize(6)
    axb.set_xlabel(r"$n$")
    axb.set_ylabel(r"$m$")
    tag(axb, "(b)", loc="lower left", bbox=dict(fc="white", ec="none", pad=1.0))
    if d.max() == 0.0:
        axb.text(0.5, 0.5, "identically 0", transform=axb.transAxes, ha="center", va="center", fontsize=8)

    axc = fig.add_subplot(gs[2])
    vm = np.max(np.abs(ref * plateau))
    im2 = axc.imshow(ref * plateau, origin="upper", cmap="RdBu_r", vmin=-vm, vmax=vm,
                     extent=(-0.5, B + 0.5, B + 0.5, -0.5))
    cb2 = fig.colorbar(im2, ax=axc, fraction=0.046, pad=0.04)
    cb2.set_label(r"$\bar\mu\,M^{\rm B.15}_{mn}$", fontsize=7)
    cb2.ax.tick_params(labelsize=6)
    axc.set_xlabel(r"$n$")
    axc.set_ylabel(r"$m$")
    tag(axc, "(c)", loc="lower left", bbox=dict(fc="white", ec="none", pad=1.0))
    for x in (axb, axc):
        x.set_xticks(range(0, B + 1, 2))
        x.set_yticks(range(0, B + 1, 2))
        x.minorticks_off()

    png = os.path.join(a.out, "fig_stage_1.png")
    fig.savefig(png, dpi=300, bbox_inches="tight")

    caption = (
        "FIG. 1. Stage 1 of the unified C++ core (handoff 0001): Chebyshev moments of the "
        "periodic tight-binding chain chain1d, N = %d sites, spectral window [-2, 2], exact trace "
        "over all %d site states, M = %d moments, KPM_SEED = 12345, computed by build/stage1_moments "
        "on goldens/stage0/chain1d/operators at OMP_NUM_THREADS = %s. "
        "(a) Rank-1 moments mu_m (G1, word dos): legacy golden (open circles), new core (dashed line) "
        "and the analytic reference delta_m0 (triangles); right axis: |new - golden| (dash-dot), "
        "max %.1e. "
        "(b) |new - golden| over the m, n <= %d block of the rank-2 moments mu_mn (G2, word kubo, "
        "KPM weights and Hermitian symmetrisation applied), max %.1e. "
        "(c) The closed-form block B.15, Mref(m, n) of test/chain1d_analytic_test.cpp times the "
        "diagonal plateau mean(mu_mm, m = 2..%d) = %.6g of the new tensor; max |Re mu_mn / plateau - "
        "Mref| over m, n < %d is %.1e, max |Im mu_mn| over the block is %.1e. "
        "Byte comparison of the legacy-layout files with goldens/stage0/chain1d (cmp criterion): %s. "
        "The new-format .chebmom tensors are %s across the thread counts and agree with the "
        "15-digit legacy-layout files to %.1e. %s"
        % (N_SITES, N_SITES, M, ", ".join(map(str, a.threads)), float(np.max(np.abs(res1))), B,
           float(np.max(d)), B - 1, plateau, B, b15_res, max_im,
           "; ".join(cmp_lines), "bit-identical" if threads_equal else "NOT bit-identical", fmt_dev,
           "All files IDENTICAL." if all_identical else "AT LEAST ONE FILE DIFFERS."))
    with open(os.path.join(a.out, "fig_stage_1.caption.txt"), "w") as f:
        f.write(caption + "\n")
    print("wrote", png)
    print(caption)
    for root, dirs, files in os.walk(tmp, topdown=False):
        for x in files:
            os.remove(os.path.join(root, x))
        for x in dirs:
            os.rmdir(os.path.join(root, x))
    os.rmdir(tmp)
    return 0 if (all_identical and threads_equal) else 1


if __name__ == "__main__":
    sys.exit(main())
