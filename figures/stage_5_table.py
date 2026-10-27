#!/usr/bin/env python3
"""Stage-5 tolerance table (handoff/0005 codex.md section 2, revised after run 1), called by
scripts/check_goldens.sh 5.

usage: figures/stage_5_table.py <run_dir> <goldens/stage0/bases>

<run_dir> holds the stage5_bases outputs of one thread count laid out by the gate as
<model>/W0/{real,bloch,blochx,eigen}/{dos,kuboxx}/ and <model>/W1/{real,bloch}/{dos,kuboxx}/
(real, bloch: random phase, 32 states, seed 12345; blochx, eigen: exact trace; W1: --disorder 1
--disorder-seed 7). The model files are goldens/stage0/bases/<model>/model.

Prints one line per check: "PASS"/"FAIL" for a gated row (counted by the gate), "INFO" for a
report-only row. Definitions:

  nrel(a, b) = max_i |a_i - b_i| / max_i |b_i|      (sup-norm relative; the gated quantity)
  prel(a, b) = max_i |a_i - b_i| / max(|a_i|, |b_i|) (pointwise, as scripts/golden_diff.py)

Gated rows: real vs Bloch moments (W = 0 and W = 1) <= 1e-9; Bloch (exact trace) vs eigen
curve <= 1e-8 (curves recomputed at full precision from the .chebmom files with
lsquant.reconstruct.dos_curve / kubo_greenwood at the run's broadening, read from
<stem>.basis.json, because the written .dat has 6 significant digits); eigen exact
(EigenBasis.exact_dos / exact_sigma_xx, Gaussian eta = G9's) vs G9 on G9's own energy grid
<= 1e-8. Units against G9: exact_dos is per cell like hr_exactdiag; exact_sigma_xx (the
stage-3 kubo_greenwood normalisation, pi Tr[v delta v delta]) is multiplied by
2 pi / (A_hr N_k) with A_hr = |a1_x a2_y|, the area hr_exactdiag cmd_sigma divides by (half
the honeycomb cell |a1 x a2|; so G9's sigma_xx_e2h is twice the e^2/h value).
Report only: KPM (Jackson, M = 256) in every basis vs G9 (linear interpolation onto the KPM
energies; DOS / N_k, sigma x 2 pi / (A_hr N_k)); G9 vs the shifted 240 x 240 mesh (finite size);
the exact Chebyshev moments of the Gamma-mesh system from hr_exactdiag's H(k) vs eigen.
"""
import glob
import json
import math
import os
import sys

import numpy as np

from lsquant import _core
from lsquant.bases import make_basis
from lsquant.models.wannier import load_wannier_model
from lsquant.reconstruct import dos_curve, kubo_greenwood

N = (48, 48, 1)
NK = 48 * 48
M = 256
TOL_MOM, TOL_BE, TOL_G9 = 1e-9, 1e-8, 1e-8
MASK = {"graphene": lambda E: np.abs(E) < 2.5,
        "haldane": lambda E: (np.abs(E) > 0.5) & (np.abs(E) < 2.5)}


def nrel(a, b):
    return float(np.max(np.abs(a - b)) / np.max(np.abs(b)))


def prel(a, b):
    d = np.abs(a - b)
    s = np.maximum(np.abs(a), np.abs(b))
    return float(np.max(np.where(s > 0, d / np.where(s > 0, s, 1), 0.0)))


def row(ok, text):
    print(("PASS  " if ok else "FAIL  ") + text)


def info(text):
    print("INFO  " + text)


def load(d, word):
    """(moments, energies, full-precision curve, written .dat, basis.json) of one run directory."""
    mom = glob.glob(os.path.join(d, "*.chebmom"))
    dat = glob.glob(os.path.join(d, "*.dat"))
    js = glob.glob(os.path.join(d, "*.basis.json"))
    if len(mom) != 1 or len(dat) != 1 or len(js) != 1:
        return None
    meta = json.load(open(js[0]))
    mu = _core.ReadChebmom(mom[0])
    E, y = (dos_curve if word == "dos" else kubo_greenwood)(mu, float(meta["broadening_meV"]))
    return np.asarray(mu.data), np.asarray(E), np.asarray(y), np.loadtxt(dat[0]), meta


def g9(gold, m, word, prefix="G9"):
    f = os.path.join(gold, m, "%s_%s.json" % (prefix, "dos" if word == "dos" else "sigma_xx"))
    d = json.load(open(f))
    return np.array(d["energy_eV"]), np.array(d["dos" if word == "dos" else "sigma_xx_e2h"])


def areas(gold, m):
    lat = np.array([[float(x) for x in l.split()] for l in open(os.path.join(gold, m, "model", m + ".uc")) if l.split()][:3])
    return abs(lat[0, 0] * lat[1, 1]), abs(lat[0, 0] * lat[1, 1] - lat[0, 1] * lat[1, 0])


def g9_eta(gold, m):
    for line in open(os.path.join(gold, m, "G9_params.txt")):
        if line.split()[0] == "eta_eV":
            return float(line.split()[1])
    raise SystemExit("eta_eV missing in G9_params.txt")


def oracle_moments(gold, m, word, lo, hi):
    """Exact Chebyshev moments of the 48 x 48 Gamma-mesh system from hr_exactdiag's read_hr/opk
    H(k) and the wannier2sparse velocity (no lsquant model code), .chebmom layout
    (kpm_weighted 1: (2 - delta_m0) factors, normalised by N_k W)."""
    sys.path.insert(0, os.environ.get("W2S_TOOLS", os.path.expanduser("~/wannier2sparse/tools")))
    from hr_exactdiag import read_hr, opk
    d = os.path.join(gold, m, "model")
    iR, H = read_hr(os.path.join(d, m + "_hr.dat"))
    n = np.arange(48) / 48.0
    K = np.stack(np.meshgrid(n, n, indexing="ij"), -1).reshape(-1, 2)
    kf = np.c_[K, np.zeros(len(K))]
    Hk = opk(iR, H, kf); Hk = 0.5 * (Hk + np.conj(np.transpose(Hk, (0, 2, 1))))
    w, U = np.linalg.eigh(Hk)
    x = (w - 0.5 * (lo + hi)) / (0.5 * (hi - lo))
    T = np.cos(np.arange(M)[None, None, :] * np.arccos(np.clip(x, -1, 1))[:, :, None])
    h = np.where(np.arange(M) == 0, 1.0, 2.0)
    if word == "dos":
        return T.sum((0, 1)) / x.size * h
    iRa, va = read_hr(os.path.join(d, m + "_vx_hr.dat"))
    V = np.einsum("kmi,kmn,knj->kij", np.conj(U), opk(iRa, va, kf), U)
    return np.einsum("kam,kab,kbn->mn", T, np.abs(V) ** 2, T) / x.size * np.outer(h, h)


def main():
    run, gold = sys.argv[1], sys.argv[2]
    for m in ("graphene", "haldane"):
        A_hr, A_true = areas(gold, m)
        eta = g9_eta(gold, m)
        eb = make_basis("eigen", load_wannier_model(os.path.join(gold, m, "model"), seed=m), N)
        for word in ("dos", "kuboxx"):
            tag = "%s %s" % (m, "DOS" if word == "dos" else "sigma_xx")
            R = {b: load(os.path.join(run, m, "W0", b, word), word) for b in ("real", "bloch", "blochx", "eigen")}
            D = {b: load(os.path.join(run, m, "W1", b, word), word) for b in ("real", "bloch")}
            missing = [b for b, v in R.items() if v is None] + ["W1 " + b for b, v in D.items() if v is None]
            if missing:
                row(False, "%s: outputs missing for %s" % (tag, ", ".join(missing)))
                continue
            for b, v in list(R.items()) + [("W1 " + b, v) for b, v in D.items()]:
                mu, E, y, dat, meta = v
                if dat.shape[0] != E.size or np.max(np.abs(dat[:, 0] - E)) > 1e-5 * np.max(np.abs(E)):
                    row(False, "%s %s: written .dat grid differs from the recomputed curve grid" % (tag, b))
                    continue
                info("%s %-8s window %s (%s), broadening %s meV; written .dat vs recomputed curve nrel %.2e (6-digit text)"
                     % (tag, b, meta["window"], meta["window_source"], meta["broadening_meV"], nrel(dat[:, 1], y)))
            # W = 0: real (wannier2sparse CSR) vs Bloch moments, the same random-phase vectors through B
            r = nrel(R["real"][0], R["bloch"][0])
            row(r <= TOL_MOM, "%s W=0 real vs bloch moments (random 32, seed 12345): nrel %.3e (tol %.0e); prel %.3e"
                % (tag, r, TOL_MOM, prel(R["real"][0], R["bloch"][0])))
            info("%s W=0 real vs bloch curve: nrel %.3e" % (tag, nrel(R["real"][2], R["bloch"][2])))
            # Bloch vs eigen, same kernel, exact trace on both sides
            r = nrel(R["blochx"][2], R["eigen"][2])
            row(r <= TOL_BE, "%s W=0 bloch(exact) vs eigen curve: nrel %.3e (tol %.0e); moments nrel %.3e"
                % (tag, r, TOL_BE, nrel(R["blochx"][0], R["eigen"][0])))
            info("%s W=0 bloch(random 32) vs eigen curve: nrel %.3e (stochastic trace error)" % (tag, nrel(R["bloch"][2], R["eigen"][2])))
            # eigen exact (Gaussian) vs G9 on G9's grid
            Eg, yg = g9(gold, m, word)
            if word == "dos":
                ex, fac = eb.exact_dos(Eg, eta), 1.0
                kfac = 1.0 / NK
            else:
                fac = kfac = 2 * math.pi / (A_hr * NK)
                ex = eb.exact_sigma_xx(Eg, eta) * fac
            r = nrel(ex, yg)
            row(r <= TOL_G9, "%s W=0 eigen exact (Gaussian %.0f meV) vs G9 (hr_exactdiag, Gamma mesh 48x48): nrel %.3e (tol %.0e); prel %.3e"
                % (tag, eta * 1000, r, TOL_G9, prel(ex, yg)))
            if word != "dos":
                info("%s units: exact_sigma_xx x 2 pi / (A_hr N_k) = %.6e, A_hr = |a1_x a2_y| = %.6f A^2 (true cell %.6f A^2: G9 is %.1f x the e^2/h value)"
                     % (tag, fac, A_hr, A_true, A_true / A_hr))
            # report only: KPM vs G9 (kernel shape), finite size, independent moments
            E = R["eigen"][1]
            msk = MASK[m](E)
            ref = np.interp(E, Eg, yg)
            for b in ("real", "bloch", "blochx", "eigen"):
                y = R[b][2] * kfac
                info("%s W=0 KPM %-6s vs G9 (report only): nrel %.3e on %s, %.3e on the whole grid"
                     % (tag, b, nrel(y[msk], ref[msk]), "|E|<2.5" if m == "graphene" else "0.5<|E|<2.5", nrel(y, ref)))
            Ef, yf = g9(gold, m, word, "finite_size_shifted240")
            info("%s finite size (report only): G9 (Gamma 48x48) vs hr_exactdiag on its shifted 240x240 mesh: nrel %.3e on the grid, %.3e on the mask"
                 % (tag, nrel(yg, yf), nrel(yg[MASK[m](Eg)], yf[MASK[m](Eg)])))
            lo, hi = R["eigen"][4]["window"]
            mo = oracle_moments(gold, m, word, lo, hi)
            info("%s W=0 independent exact moments (hr_exactdiag H(k) + wannier2sparse VX, Gamma mesh): eigen nrel %.3e, bloch(exact) nrel %.3e"
                 % (tag, nrel(R["eigen"][0].reshape(mo.shape), mo), nrel(R["blochx"][0].reshape(mo.shape), mo)))
            # W = 1: real vs Bloch with the disorder sandwich B^+ V B
            r = nrel(D["real"][0], D["bloch"][0])
            row(r <= TOL_MOM, "%s W=1 (disorder seed 7) real vs bloch moments (random 32, seed 12345): nrel %.3e (tol %.0e); prel %.3e"
                % (tag, r, TOL_MOM, prel(D["real"][0], D["bloch"][0])))
            info("%s W=1 real vs bloch curve: nrel %.3e" % (tag, nrel(D["real"][2], D["bloch"][2])))
            info("%s W=1 eigen: skipped as codex.md 2.4 directs" % tag)


if __name__ == "__main__":
    main()
