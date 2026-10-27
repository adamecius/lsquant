#!/usr/bin/env python
"""
siesta2w90_sisl.py -- write Wannier90 inputs (.win, .eig, .amn, .mmn) from a SIESTA .HSX,
including spin-orbit (spinor) Hamiltonians, which SIESTA's own Wannier90 interface does not
support (siesta2wannier90: "non-collinear spin not implemented yet").

  A_mn(k)   = <psi_mk | g_n>           = (C_k^H S(k))_{m, g_n}     (g = basis orbitals)
  M_mn(k,b) = <u_mk | u_n,k+b>  ~=  C_k^H S(k) diag(exp(-i b.tau_nu)) C_{k'}
      (tight-binding approximation: e^{-i b.r} evaluated at the centre of each orbital;
       k + b = k' + G, cell gauge so C_{k'+G} = C_{k'})

Energies are written relative to E_F (sisl convention), so windows are relative to E_F.
Usage:
  siesta2w90_sisl.py RUNDIR LABEL SEED --mp 8 14 1 --bands 33 70 --proj 10-14,29-33,... \
      --win -6.4 12.5 --froz -6.4 1.6 --num-iter 200 [--spin 0|1]
Then: wannier90.x -pp SEED ; siesta2w90_sisl.py ... --stage mats ; wannier90.x SEED
"""
import argparse, os, subprocess, sys
import numpy as np, scipy.linalg as sl
import sisl

BOHR = 0.529177210903


def parse_ranges(s):
    out = []
    for part in s.split(","):
        a, *b = part.split("-")
        out += list(range(int(a), int(b[0]) + 1)) if b else [int(a)]
    return out


def HSk(H, k, spin):
    Hk = H.Hk(k, spin=spin, format="array") if H.spin.is_polarized else H.Hk(k, format="array")
    Sk = H.Sk(k, format="array")
    if Sk.shape[0] != Hk.shape[0]:
        Sk = np.kron(Sk, np.eye(2))
    return Hk, Sk


def main():
    p = argparse.ArgumentParser()
    p.add_argument("rundir"); p.add_argument("label"); p.add_argument("seed")
    p.add_argument("--mp", type=int, nargs=3, required=True)
    p.add_argument("--bands", type=int, nargs=2, required=True, help="1-based, inclusive (spinor bands if SOC)")
    p.add_argument("--proj", required=True, help="1-based SIESTA orbital indices (spin-less numbering)")
    p.add_argument("--win", type=float, nargs=2, required=True)
    p.add_argument("--froz", type=float, nargs=2, required=True)
    p.add_argument("--num-iter", type=int, default=200)
    p.add_argument("--spin", type=int, default=0)
    p.add_argument("--stage", choices=["win", "mats"], default="win")
    a = p.parse_args()

    H = sisl.get_sile(os.path.join(a.rundir, f"{a.label}.HSX")).read_hamiltonian()
    g = H.geometry
    spinor = H.spin.is_spinorbit or H.spin.is_noncolinear
    ns = 2 if spinor else 1
    proj0 = [i - 1 for i in parse_ranges(a.proj)]
    proj = [ns * i + s for i in proj0 for s in range(ns)]          # sisl interleaves spinor components
    tau = np.repeat(g.xyz[g.o2a(np.arange(H.no))], ns, axis=0)      # orbital centres (Ang)
    b0, b1 = a.bands[0] - 1, a.bands[1]
    nb, nw = b1 - b0, len(proj)
    N = a.mp
    ks = np.array([[i / N[0], j / N[1], l / N[2]] for i in range(N[0]) for j in range(N[1]) for l in range(N[2])])

    if a.stage == "win":
        # projections block: only used by -pp for bookkeeping; the real A comes from our .amn
        species_l = {}
        for i in proj0:
            ia = g.o2a(i); o = g.atoms[ia].orbitals[i - g.a2o(ia)]
            species_l.setdefault((g.atoms[ia].tag, o.l), True)
        lines = []
        lines += [f"num_wann = {nw}", f"num_bands = {nb}", f"num_iter = {a.num_iter}",
                  f"dis_win_min = {a.win[0]}", f"dis_win_max = {a.win[1]}",
                  f"dis_froz_min = {a.froz[0]}", f"dis_froz_max = {a.froz[1]}",
                  "dis_num_iter = 500", "write_hr = true", "use_ws_distance = true",
                  "write_xyz = true", "guiding_centres = true"]
        if spinor:
            lines.append("spinors = true")
        lines += ["begin unit_cell_cart", "ang"] + [" %.8f %.8f %.8f" % tuple(v) for v in g.cell] + ["end unit_cell_cart"]
        lines += ["begin atoms_cart", "ang"] + [f" {at.tag} %.8f %.8f %.8f" % tuple(x) for at, x in zip(g.atoms, g.xyz)] + ["end atoms_cart"]
        lines += ["begin projections"] + [f" {sp}:{'spdf'[l]}" for (sp, l) in species_l] + ["end projections"]
        lines += [f"mp_grid = {N[0]} {N[1]} {N[2]}", "begin kpoints"] + [" %.10f %.10f %.10f" % tuple(k) for k in ks] + ["end kpoints"]
        open(f"{a.seed}.win", "w").write("\n".join(lines) + "\n")
        print(f"{a.seed}.win written: num_wann={nw}, num_bands={nb}, spinor={spinor}")
        return

    # ---- read nnkp
    txt = open(f"{a.seed}.nnkp").read().split("begin nnkpts")[1].split("end nnkpts")[0].split()
    nnb = int(txt[0]); vals = list(map(int, txt[1:]))
    nn = np.array(vals).reshape(-1, 5)                               # k, k', G1, G2, G3 (1-based)
    rc = g.rcell                                                     # 2pi/Ang, rows
    E, C, S = [], [], []
    for k in ks:
        Hk, Sk = HSk(H, k, a.spin)
        e, c = sl.eigh(Hk, Sk)
        E.append(e[b0:b1]); C.append(c[:, b0:b1]); S.append(Sk)
    with open(f"{a.seed}.eig", "w") as fh:
        for ik in range(len(ks)):
            for n in range(nb):
                fh.write(f"{n+1:5d}{ik+1:5d}{E[ik][n]:18.12f}\n")
    with open(f"{a.seed}.amn", "w") as fh:
        fh.write("from siesta2w90_sisl\n")
        fh.write(f"{nb} {len(ks)} {nw}\n")
        for ik in range(len(ks)):
            A = C[ik].conj().T @ S[ik][:, proj]                          # (nb, nw)
            for n in range(nw):
                for m in range(nb):
                    fh.write(f"{m+1:5d}{n+1:5d}{ik+1:5d}{A[m,n].real:18.12f}{A[m,n].imag:18.12f}\n")
    with open(f"{a.seed}.mmn", "w") as fh:
        fh.write("from siesta2w90_sisl (TB approximation)\n")
        fh.write(f"{nb} {len(ks)} {nnb}\n")
        for k1, k2, G1, G2, G3 in nn:
            i1, i2 = k1 - 1, k2 - 1
            b = (ks[i2] + np.array([G1, G2, G3]) - ks[i1]) @ rc            # cartesian 1/Ang
            ph = np.exp(-1j * tau @ b)
            M = C[i1].conj().T @ S[i1] @ (ph[:, None] * C[i2])
            fh.write(f"{k1:5d}{k2:5d}{G1:5d}{G2:5d}{G3:5d}\n")
            for n in range(nb):
                for m in range(nb):
                    fh.write(f"{M[m,n].real:18.12f}{M[m,n].imag:18.12f}\n")
    print(f"{a.seed}: .eig .amn .mmn written ({len(ks)} k, {nnb} neighbours)")


if __name__ == "__main__":
    main()
