#!/usr/bin/env python3
"""Rashba (Dyakonov-Perel) model on the 2D square lattice (spinful).

Square lattice with spin-diagonal hopping, a UNIFORM Rashba spin-orbit coupling,
and scalar (Anderson) onsite disorder that sets the momentum-scattering time:

    H = -t sum_<ij> c_i^dag c_j (x) 1_spin
        + lam sum_i [ c_{i+x}^dag (i sigma_y) c_i - c_{i+y}^dag (i sigma_x) c_i ] + h.c.
        + sum_i eps_i (x) 1_spin,     eps_i uniform in [-W/2, W/2].

The Rashba term is an in-plane k-dependent spin field; with disorder it relaxes
the spin by the Dyakonov-Perel mechanism. The spin operators / projectors come
from make_ey (same spin doubling). Operators are written under operators/ at run
time (nothing committed).

    python make_rashba.py <L> <lambda> <W> [seed]
"""

import os
import sys
import random

import make_ey   # reuse _csr and the spin-operator / projector writers


def _write_ham(path, L, lam, W, seed, t=1.0):
    rng = random.Random(seed)
    Ns = L * L
    N = 2 * Ns
    eps = [rng.uniform(-0.5 * W, 0.5 * W) for _ in range(Ns)]
    rows = [dict() for _ in range(N)]

    def add(b, a, M):                       # directed bond a->b (2x2 spin M) + h.c.
        for sp in (0, 1):
            for s in (0, 1):
                v = M[sp][s]
                if v != 0:
                    rb, ca = 2 * b + sp, 2 * a + s
                    rows[rb][ca] = rows[rb].get(ca, 0j) + v
                    rows[ca][rb] = rows[ca].get(rb, 0j) + complex(v).conjugate()

    hop = [[-t + 0j, 0j], [0j, -t + 0j]]                       # -t * I
    rx = [[0j, lam + 0j], [-lam + 0j, 0j]]                     # lam * i sigma_y
    ry = [[0j, complex(0, -lam)], [complex(0, -lam), 0j]]      # -lam * i sigma_x
    for i in range(Ns):
        x, y = i % L, i // L
        jx, jy = y * L + (x + 1) % L, ((y + 1) % L) * L + x
        add(jx, i, hop); add(jx, i, rx)                        # +x bond: hop + Rashba_x
        add(jy, i, hop); add(jy, i, ry)                        # +y bond: hop + Rashba_y
    for i in range(Ns):                                        # onsite disorder
        for s in (0, 1):
            rows[2 * i + s][2 * i + s] = rows[2 * i + s].get(2 * i + s, 0j) + eps[i]

    vals, cols, rowptr = [], [], [0]
    for r in range(N):
        for c in sorted(rows[r]):
            v = rows[r][c]
            cols.append(c)
            vals.append(f"{v.real:.10g} {v.imag:.10g}")
        rowptr.append(len(cols))
    make_ey._csr(path, N, vals, cols, rowptr)


def _write_bounds(path, lam, W):
    edge = 4.0 + 2.0 * abs(lam) + 0.5 * W + 0.2
    with open(path, "w") as f:
        f.write(f"{-edge:g} {edge:g}\n")


def build_rashba(L, lam, W, seed, outdir="."):
    label = f"rashba_L{L}_l{lam:g}_W{W:g}"
    opdir = os.path.join(outdir, "operators")
    os.makedirs(opdir, exist_ok=True)
    _write_ham(os.path.join(opdir, f"{label}.HAM.CSR"), L, lam, W, seed)
    make_ey._write_spin(opdir, label, L)            # SX, SY, SZ, PZ, PX
    _write_bounds(os.path.join(outdir, "BOUNDS"), lam, W)
    with open(os.path.join(opdir, f"{label}.HAM.desc"), "w") as f:
        edge = 4.0 + 2.0 * abs(lam) + 0.5 * W + 0.2
        f.write("observable: hamiltonian\nunits: eV\n")
        f.write(f"system_label: {label}\nsystem_size: {2*L*L}\n")
        f.write("dimensions: 2\nperiodic: true\nspin: spinful\n")
        f.write(f"rashba_lambda_eV: {lam:g}\ndisorder_W_eV: {W:g}\ndisorder_seed: {seed}\n")
        f.write(f"has_bounds: true\nspectral_min: {-edge:g}\nspectral_max: {edge:g}\n")
        f.write("provenance: make_rashba.py: 2D square + uniform Rashba + Anderson disorder\n")
    return label


def main():
    if len(sys.argv) not in (4, 5):
        print("usage: python make_rashba.py <L> <lambda> <W> [seed]")
        sys.exit(1)
    L, lam, W = int(sys.argv[1]), float(sys.argv[2]), float(sys.argv[3])
    seed = int(sys.argv[4]) if len(sys.argv) == 5 else 1
    print("built", build_rashba(L, lam, W, seed), f": {2*L*L} orbitals")


if __name__ == "__main__":
    main()
