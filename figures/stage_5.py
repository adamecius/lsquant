#!/usr/bin/env python3
"""Stage-5 figure (handoff 0005, codex.md section 3): graphene DOS and sigma_xx in three bases.

Reads the stage5_bases outputs kept by the gate (`S5_KEEP=<dir> scripts/check_goldens.sh 5`,
layout <dir>/T16/graphene/W{0,1}/<basis>/{dos,kuboxx}/, each with its <stem>.basis.json) and
the committed G9 goldens (goldens/stage0/bases/graphene/G9_dos.json, G9_sigma_xx.json:
hr_exactdiag on the Gamma-centred 48 x 48 mesh, Gaussian 40 meV). KPM curves are recomputed
at full precision from the moment files with lsquant.reconstruct at the run's broadening
(from basis.json), as figures/stage_5_table.py does; the eigen exact curve is
EigenBasis.exact_dos / exact_sigma_xx with the G9 Gaussian.

Panels: (a) DOS, W = 0: KPM real (random phase 32), Bloch (random phase 32), eigen (exact
trace) as lines, G9 as markers, EigenBasis exact (Gaussian) as a thin black line; residuals
below: eigen exact - G9 (the gated row; KPM - G9, the kernel shape, is in the caption only);
(b) DOS, W = 1 eV (disorder seed 7): real and Bloch as lines, the clean G9 as markers for
reference, residual real - Bloch below; (c), (d) the same for sigma_xx. Units: DOS per cell;
sigma_xx in e^2/h with the true honeycomb cell A = |a1 x a2| (KPM x 2 pi / (A N_k); G9 x
A_hr / A because hr_exactdiag cmd_sigma divides by A_hr = |a1_x a2_y|, half the cell).

Outputs figures/out/stage_5/fig_stage_5.png and fig_stage_5.caption.txt.
Usage (repository root, on TCN):
    MPLBACKEND=Agg PYTHONPATH=<build>/lsquant_core:python python3 figures/stage_5.py --runs <dir>/T16 \
        [--table <gate table file>]
"""
import argparse
import glob
import json
import math
import os
import re

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from lsquant import _core
from lsquant.bases import make_basis
from lsquant.models.wannier import load_wannier_model
from lsquant.reconstruct import dos_curve, kubo_greenwood

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
GOLD = os.path.join(REPO, "goldens", "stage0", "bases", "graphene")
NK, M = 48 * 48, 256
COL = {"real": "#0072B2", "bloch": "#E69F00", "eigen": "#009E73", "g9": "0.15"}
LBL = {"real": "KPM real space (random 32)", "bloch": "KPM Bloch (random 32)", "eigen": "KPM eigen (exact trace)"}


def curve(d, word):
    meta = json.load(open(glob.glob(os.path.join(d, "*.basis.json"))[0]))
    mu = _core.ReadChebmom(glob.glob(os.path.join(d, "*.chebmom"))[0])
    E, y = (dos_curve if word == "dos" else kubo_greenwood)(mu, float(meta["broadening_meV"]))
    return np.asarray(E), np.asarray(y), meta


def nrel(a, b):
    return float(np.max(np.abs(a - b)) / np.max(np.abs(b)))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--runs", required=True, help="<S5_KEEP>/T16")
    ap.add_argument("--table", default=None, help="gate table output (stage_5_table.py) for the caption")
    ap.add_argument("--out", default=os.path.join(REPO, "figures", "out", "stage_5"))
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    par = {l.split()[0]: l.split()[1] for l in open(os.path.join(GOLD, "G9_params.txt"))}
    eta = float(par["eta_eV"])
    lat = np.array([[float(x) for x in l.split()] for l in open(os.path.join(GOLD, "model", "graphene.uc")) if l.split()][:3])
    A_hr = abs(lat[0, 0] * lat[1, 1])
    A = abs(lat[0, 0] * lat[1, 1] - lat[0, 1] * lat[1, 0])
    unit = {"dos": 1.0 / NK, "kuboxx": 2 * math.pi / (A * NK)}
    j = {w: json.load(open(os.path.join(GOLD, "G9_%s.json" % ("dos" if w == "dos" else "sigma_xx")))) for w in ("dos", "kuboxx")}
    g9 = {"dos": (np.array(j["dos"]["energy_eV"]), np.array(j["dos"]["dos"])),
          "kuboxx": (np.array(j["kuboxx"]["energy_eV"]), np.array(j["kuboxx"]["sigma_xx_e2h"]) * A_hr / A)}
    eb = make_basis("eigen", load_wannier_model(os.path.join(GOLD, "model"), seed="graphene"), (48, 48, 1))
    exact = {"dos": eb.exact_dos(g9["dos"][0], eta), "kuboxx": eb.exact_sigma_xx(g9["kuboxx"][0], eta) * unit["kuboxx"]}

    plt.rcParams.update({"font.size": 8, "axes.linewidth": 0.6, "lines.linewidth": 1.0, "axes.unicode_minus": False})
    fig = plt.figure(figsize=(7.0, 6.6))
    gs = fig.add_gridspec(5, 2, height_ratios=[3, 1.2, 0.45, 3, 1.2], hspace=0.08, wspace=0.30)
    stats, meta = {}, {}
    letters = iter("abcd")
    for r, word in enumerate(("dos", "kuboxx")):
        for c, W in enumerate((0, 1)):
            ax = fig.add_subplot(gs[3 * r, c]); rx = fig.add_subplot(gs[3 * r + 1, c], sharex=ax)
            Eg, yg = g9[word]
            bases = ("real", "bloch", "eigen") if W == 0 else ("real", "bloch")
            C = {}
            for b in bases:
                E, y, meta[(word, W)] = curve(os.path.join(a.runs, "graphene", "W%d" % W, b, word), word)
                C[b] = (E, y)
                ax.plot(E, y * unit[word], color=COL[b], lw=1.6 if b == "real" else 1.0,
                        ls={"real": "-", "bloch": "--", "eigen": ":"}[b], label=LBL[b])
            step = max(1, len(Eg) // 330)
            ax.plot(Eg[::step], yg[::step], "o", ms=1.6, mfc="none", mec=COL["g9"], mew=0.4,
                    label="hr_exactdiag G9, Gaussian" + (" (clean)" if W else ""))
            E = C["real"][0]
            if W == 0:
                ax.plot(Eg, exact[word], color="k", lw=0.4, label="EigenBasis exact, Gaussian")
                rx.plot(Eg, (exact[word] - yg) * 1e14, color="k", lw=0.6, label="exact - G9 (x 1e14)")
                ref = np.interp(E, Eg, yg)
                stats[(word, 0)] = {"exact": nrel(exact[word], yg),
                                    "kpm": {b: nrel(C[b][1] * unit[word], ref) for b in bases}}
                rx.set_ylabel("exact $-$ G9\n($10^{-14}$)", fontsize=6)
            else:
                d = (C["real"][1] - C["bloch"][1]) * unit[word]
                rx.plot(E, d * 1e14, color=COL["real"], lw=0.8)
                stats[(word, 1)] = nrel(C["real"][1], C["bloch"][1])
                rx.set_ylabel("real $-$ Bloch\n($10^{-14}$)", fontsize=6)
            rx.axhline(0, color="0.6", lw=0.5)
            ax.set_xlim(-3.7, 3.7)
            plt.setp(ax.get_xticklabels(), visible=False)
            if r == 1:
                rx.set_xlabel("$E$ (eV)")
            ax.set_ylabel(r"DOS (eV$^{-1}$ cell$^{-1}$)" if word == "dos" else r"$\sigma_{xx}$ ($e^2/h$)")
            ax.set_title("(%s) %s, $W$ = %d eV" % (next(letters), "DOS" if word == "dos" else r"$\sigma_{xx}$", W),
                         loc="left", fontsize=8)
            ax.legend(loc="upper right", fontsize=5.3, frameon=False, handlelength=2.2)
            rx.tick_params(labelsize=6)
    png = os.path.join(a.out, "fig_stage_5.png")
    fig.savefig(png, dpi=300, bbox_inches="tight")

    tab = ""
    if a.table and os.path.isfile(a.table):
        txt = open(a.table).read()

        def num(pat):
            mm = re.search(pat, txt)
            return mm.group(1) if mm else "n/a"
        parts = []
        for w in ("DOS", "sigma_xx"):
            parts.append("%s: real vs Bloch moments %s (W = 0), %s (W = 1), tol 1e-9; Bloch (exact trace) vs eigen curve %s, tol 1e-8; "
                         "EigenBasis exact vs G9 %s, tol 1e-8"
                         % (w, num(r"graphene %s W=0 real vs bloch moments[^:]*: nrel (\S+)" % w),
                            num(r"graphene %s W=1 \(disorder seed 7\) real vs bloch moments[^:]*: nrel (\S+)" % w),
                            num(r"graphene %s W=0 bloch\(exact\) vs eigen curve: nrel (\S+)" % w),
                            num(r"graphene %s W=0 eigen exact [^:]*: nrel (\S+)" % w)))
        tab = " Gate rows (figures/stage_5_table.py, T = 16; identical at T = 1), graphene: " + "; ".join(parts) + "."
    m0, m1 = meta[("dos", 0)], meta[("dos", 1)]
    s0d, s0s = stats[("dos", 0)], stats[("kuboxx", 0)]
    caption = (
        "FIG. 5. Graphene (wannier2sparse examples/gen_models.py, t = -1 eV, nearest-neighbour distance 1 A) on the "
        "48 x 48 x 1 supercell (4608 orbitals; Bloch mesh k = n/48), evaluated by lsquant.apps.stage5_bases in three bases "
        "from one model file, M = 256 Chebyshev moments, Jackson kernel, broadening %s meV. Top: DOS per cell; bottom: "
        "sigma_xx (Kubo-Greenwood) in e^2/h with the true cell area %.4f A^2. Left, W = 0, window %s eV (--window eigen): "
        "real space (the wannier2sparse CSR, random phase, 32 states, seed 12345, solid), Bloch (block-diagonal H(k), the "
        "same 32 vectors mapped by B, dashed), eigen (exact trace in the eigenbasis of each H(k), dotted); markers: G9, "
        "tools/hr_exactdiag.py (wannier2sparse 3f553db) dos and sigma --comp xx on the same Gamma-centred 48 x 48 mesh "
        "(scripts/hr_exactdiag_gamma.py), Gaussian eta = %.0f meV; thin black line: EigenBasis.exact_dos/exact_sigma_xx "
        "with the same Gaussian. Residual below: exact - G9 (units 1e-14; max relative %.1e DOS, %.1e sigma_xx; the gated "
        "row). KPM eigen - G9 (Jackson vs Gaussian kernel shape, report only, not drawn): max relative %.2e DOS, %.2e sigma_xx. "
        "Right, W = 1 eV Anderson disorder (uniform in [-0.5, 0.5] eV, --disorder-seed 7), window %s eV: real space "
        "(CSR + diag V) and Bloch (H(k) + B^+ V B through FFTW); clean G9 markers for reference; residual real - Bloch "
        "(units 1e-14; max relative on the curve %.1e DOS, %.1e sigma_xx); eigen skipped (no per-k eigenbasis with "
        "disorder).%s"
        % (m0["broadening_meV"], A, m0["window"], eta * 1000, s0d["exact"], s0s["exact"], s0d["kpm"]["eigen"], s0s["kpm"]["eigen"],
           m1["window"], stats[("dos", 1)], stats[("kuboxx", 1)], tab))
    with open(os.path.join(a.out, "fig_stage_5.caption.txt"), "w") as f:
        f.write(caption + "\n")
    print(png)
    print(caption)


if __name__ == "__main__":
    main()
