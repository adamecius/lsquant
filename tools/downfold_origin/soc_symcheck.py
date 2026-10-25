#!/usr/bin/env python
"""Symmetry checks of spin-orbit Hamiltonians (SIESTA reference and Wannier/TB models).
  * Kramers splitting |E_2n+1 - E_2n| in a window: with P*T symmetry (1T'-WTe2 is centrosymmetric)
    every band must be doubly degenerate at every k; at TRIM points T alone forces it.
  * Time reversal in real space: H(R) = Theta H(R)^* Theta^dagger, Theta = 1 (x) i sigma_y
    (spinor components interleaved: index 2p+s, s=0 up, 1 down).
  * Hermiticity: H(R) = H(-R)^dagger.
usage: soc_symcheck.py RUNDIR LABEL --hr name:prefix:shift ... --dat name:hoppings.dat ...
"""
import argparse, os, sys
import numpy as np, scipy.linalg as sl, sisl
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from compare_w90 import read_hr, model_bands

TRIM = {"Γ": [0, 0, 0], "X": [.5, 0, 0], "S": [.5, .5, 0], "Y": [0, .5, 0]}


def kramers(E, lo, hi):
    """max and rms splitting of consecutive pairs (2n, 2n+1) inside [lo, hi]."""
    sp = []
    for e in E:
        e = np.sort(e)
        for n in range(0, len(e) - 1, 2):
            if lo <= e[n] <= hi:
                sp.append(abs(e[n + 1] - e[n]))
    sp = np.array(sp)
    return sp.max(), np.sqrt((sp ** 2).mean())


def hr_dict(hr):
    nw, R, m, n, t, degR, ws = hr
    D = {}
    for r, a, b, val in zip(R, m, n, t):
        D.setdefault(tuple(r), np.zeros((nw, nw), complex))[a, b] = val
    return D


def dat_dict(fn):
    rows = [l.split() for l in open(fn) if not l.startswith("#") and l.strip()]
    J = max(max(int(r[3]), int(r[4])) for r in rows) + 1
    D = {}
    for r in rows:
        R = tuple(map(int, r[:3])); i, j = int(r[3]), int(r[4])
        D.setdefault(R, np.zeros((J, J), complex))[i, j] += float(r[5]) + 1j * float(r[6])
    return D


def dict_bands(D, kpts):
    out = []
    for k in kpts:
        Hk = sum(M * np.exp(2j * np.pi * np.dot(k, R)) for R, M in D.items())
        out.append(np.linalg.eigvalsh(0.5 * (Hk + Hk.conj().T)))
    return np.array(out)


def realspace_checks(D):
    J = next(iter(D.values())).shape[0]
    Th = np.kron(np.eye(J // 2), np.array([[0, 1], [-1, 0]]))       # i sigma_y
    scale = max(np.abs(M).max() for M in D.values())
    trs = max(np.abs(M - Th @ M.conj() @ Th.T).max() for M in D.values())
    herm = max(np.abs(M - D.get(tuple(-np.array(R)), np.zeros_like(M)).conj().T).max() for R, M in D.items())
    return trs, herm, scale


def main():
    p = argparse.ArgumentParser()
    p.add_argument("rundir"); p.add_argument("label")
    p.add_argument("--hr", action="append", default=[])
    p.add_argument("--dat", action="append", default=[])
    p.add_argument("--win", type=float, nargs=2, default=[-6.4, 1.6])
    a = p.parse_args()
    rng = np.random.default_rng(7)
    kr = rng.random((60, 3)); kr[:, 2] = 0
    H = sisl.get_sile(os.path.join(a.rundir, f"{a.label}.HSX")).read_hamiltonian()

    def ref(kpts):
        out = []
        for k in kpts:
            Hk = H.Hk(k, format="array"); Sk = H.Sk(k, format="array")
            if Sk.shape[0] != Hk.shape[0]:
                Sk = np.kron(Sk, np.eye(2))
            out.append(sl.eigh(Hk, Sk, eigvals_only=True))
        return np.array(out)

    lo, hi = a.win
    print(f"{'modelo':34s} {'Kramers máx (k gen.)':>21s} {'RMS':>9s} {'Kramers máx (TRIM)':>19s} {'TRS H(R)':>10s} {'Hermit.':>9s}")
    E = ref(kr); Et = ref(list(TRIM.values()))
    m1, r1 = kramers(E, lo, hi); m2, _ = kramers(Et, lo, hi)
    print(f"{'SIESTA (referencia)':34s} {m1*1e3:17.3f} meV {r1*1e3:6.3f} meV {m2*1e3:15.4f} meV {'—':>10s} {'—':>9s}")
    for spec in a.hr:
        name, prefix, shift = spec.split(":")
        hr = read_hr(prefix)
        E = model_bands(hr, kr); Et = model_bands(hr, list(TRIM.values()))
        m1, r1 = kramers(E - float(shift), lo, hi); m2, _ = kramers(Et - float(shift), lo, hi)
        trs, herm, sc = realspace_checks(hr_dict(hr))
        print(f"{name:34s} {m1*1e3:17.3f} meV {r1*1e3:6.3f} meV {m2*1e3:15.4f} meV {trs*1e3:7.3f} meV {herm*1e3:6.3f} meV")
    for spec in a.dat:
        name, fn = spec.split(":")
        D = dat_dict(fn)
        E = dict_bands(D, kr); Et = dict_bands(D, list(TRIM.values()))
        m1, r1 = kramers(E, lo, hi); m2, _ = kramers(Et, lo, hi)
        trs, herm, sc = realspace_checks(D)
        print(f"{name:34s} {m1*1e3:17.3f} meV {r1*1e3:6.3f} meV {m2*1e3:15.4f} meV {trs*1e3:7.3f} meV {herm*1e3:6.3f} meV")


if __name__ == "__main__":
    main()
