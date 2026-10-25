#!/usr/bin/env python
"""Compare Wannier90 (_hr.dat + _wsvec.dat) bands with exact SIESTA bands on Gamma-X-S-Y-Gamma.
usage: compare_w90.py RUNDIR LABEL OUTPNG --model name:seed_prefix:EFshift [--model ...]
       [--spin 0] [--froz LO HI] (relative to E_F)
EFshift: value (eV) subtracted from the model energies (0 if already relative to E_F)."""
import argparse, os, re
import numpy as np, scipy.linalg as sl
import sisl
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt


def read_hr(prefix):
    L = open(prefix + "_hr.dat").read().split("\n")
    nw = int(L[1]); nR = int(L[2])
    deg, i = [], 3
    while len(deg) < nR:
        deg += list(map(int, L[i].split())); i += 1
    data = np.array([l.split() for l in L[i:] if l.strip()], float)
    R = data[:, :3].astype(int); m = data[:, 3].astype(int) - 1; n = data[:, 4].astype(int) - 1
    t = data[:, 5] + 1j * data[:, 6]
    Rlist = [tuple(r) for r in R[::nw * nw]]
    degR = {r: d for r, d in zip(Rlist, deg)}
    ws = {}
    if os.path.exists(prefix + "_wsvec.dat"):
        W = open(prefix + "_wsvec.dat").read().split("\n")[1:]
        j = 0
        while j < len(W) and W[j].strip():
            r1, r2, r3, a, b = map(int, W[j].split()); nd = int(W[j + 1])
            sh = [tuple(map(int, W[j + 2 + q].split())) for q in range(nd)]
            ws[(r1, r2, r3, a - 1, b - 1)] = sh
            j += 2 + nd
    return nw, R, m, n, t, degR, ws


def model_bands(hr, kpts):
    """Vectorised Wannier interpolation with Wigner-Seitz shifts."""
    nw, R, m, n, t, degR, ws = hr
    RR, AA, BB, WW = [], [], [], []
    for r, a, b, val in zip(R, m, n, t):
        rr = tuple(r); w = val / degR[rr]
        shifts = ws.get((*rr, a, b), [(0, 0, 0)])
        for sft in shifts:
            RR.append(np.array(rr) + np.array(sft)); AA.append(a); BB.append(b); WW.append(w / len(shifts))
    RR = np.array(RR, float); AA = np.array(AA); BB = np.array(BB); WW = np.array(WW)
    out = []
    for k in kpts:
        Hk = np.zeros((nw, nw), complex)
        np.add.at(Hk, (AA, BB), WW * np.exp(2j * np.pi * (RR @ np.asarray(k))))
        out.append(np.linalg.eigvalsh(0.5 * (Hk + Hk.conj().T)))
    return np.array(out)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("rundir"); p.add_argument("label"); p.add_argument("png")
    p.add_argument("--model", action="append", required=True)
    p.add_argument("--spin", type=int, default=0)
    p.add_argument("--froz", type=float, nargs=2, default=[-6.4, 1.6])
    a = p.parse_args()
    H = sisl.get_sile(os.path.join(a.rundir, f"{a.label}.HSX")).read_hamiltonian()
    bs = sisl.BandStructure(H, [[0, 0, 0], [.5, 0, 0], [.5, .5, 0], [0, .5, 0], [0, 0, 0]], 160, ["Γ", "X", "S", "Y", "Γ"])
    Eref = []
    for k in bs.k:
        Hk = H.Hk(k, spin=a.spin, format="array") if H.spin.is_polarized else H.Hk(k, format="array")
        Sk = H.Sk(k, format="array")
        if Sk.shape[0] != Hk.shape[0]:
            Sk = np.kron(Sk, np.eye(2))
        Eref.append(sl.eigh(Hk, Sk, eigvals_only=True))
    Eref = np.array(Eref)
    x = bs.lineark(); ticks, labels = bs.lineark(ticks=True)[1:]
    fig, ax = plt.subplots(figsize=(7, 5.5))
    for nb in range(Eref.shape[1]):
        ax.plot(x, Eref[:, nb], c="0.6", lw=0.9, label="SIESTA DZP (exacto)" if nb == 0 else None)
    cols = ["#c0504d", "#1f4e79", "#2e7d32", "#7b3fa0"]
    print(f"{'modelo':28s} {'J':>3s} {'err_max ventana':>16s} {'RMS ventana':>12s} {'err_max (≤F+3)':>15s}")
    for c, spec in zip(cols, a.model):
        name, prefix, shift = spec.split(":"); shift = float(shift)
        hr = read_hr(prefix)
        Em = model_bands(hr, bs.k) - shift
        lo, hi = a.froz
        errs, errs3 = [], []
        for ik in range(len(bs.k)):
            er = Eref[ik]
            sel = er[(er >= lo) & (er <= hi)]
            if len(sel):
                errs.append(np.abs(sel[:, None] - Em[ik][None, :]).min(1))
            sel3 = er[(er >= lo) & (er <= hi + 3)]
            if len(sel3):
                errs3.append(np.abs(sel3[:, None] - Em[ik][None, :]).min(1))
        e = np.concatenate(errs); e3 = np.concatenate(errs3)
        print(f"{name:28s} {hr[0]:3d} {e.max()*1e3:13.1f} meV {np.sqrt((e**2).mean())*1e3:9.1f} meV {e3.max()*1e3:12.1f} meV")
        for nb in range(Em.shape[1]):
            ax.plot(x, Em[:, nb], c=c, lw=1.2, ls="--", label=f"{name} (J={hr[0]})" if nb == 0 else None)
    ax.axhspan(*a.froz, color="#1f4e79", alpha=0.06, lw=0)
    ax.axhline(0, c="0.4", ls=":", lw=0.6)
    ax.set_xticks(ticks); ax.set_xticklabels(labels); [ax.axvline(t, c="0.85", lw=0.5) for t in ticks]
    ax.set_xlim(x[0], x[-1]); ax.set_ylim(-7, 4); ax.set_ylabel("E − E$_F$ (eV)")
    ax.legend(fontsize=7, frameon=False, loc="lower right"); ax.set_title(os.path.basename(a.png).replace(".png", ""))
    fig.tight_layout(); fig.savefig(a.png, dpi=130)


if __name__ == "__main__":
    main()
