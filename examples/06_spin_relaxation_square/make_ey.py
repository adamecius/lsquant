#!/usr/bin/env python3
"""Elliott-Yafet-like spin-flip model on the 2D square lattice (spinful).

Two spin orbitals per site (index 2*site + s, s=0 up, s=1 down) on an L x L
periodic lattice. The clean band acts in each spin channel; a random onsite
spin-flip field relaxes the spin:

    H = H0 (x) 1_spin  +  sum_i ( eta_i c^dag_{i up} c_{i dn} + h.c. ),
    H0 = -t sum_<ij> |i><j|,   eta_i real Gaussian, <eta> = 0, <eta^2> = Dsf^2.

It also writes the spin operators S_x, S_y, S_z and the projector
P_z = 1/2(1 + S_z), which prepares the S_z = +1 initial state.

Operators are written under operators/ at run time (nothing committed). Functions
are importable; the CLI builds one labelled set.

    python make_ey.py <L> <Dsf> [seed]
"""

import os
import sys
import random


def _csr(path, n, vals, cols, rowptr):
    with open(path, "w") as f:
        f.write(f"{n} {len(cols)}\n")
        f.write(" ".join(vals) + " \n")
        f.write(" ".join(map(str, cols)) + " \n")
        f.write(" ".join(map(str, rowptr)) + " \n")


def _write_ham(path, L, Dsf, seed, t=1.0):
    rng = random.Random(seed)
    Ns = L * L
    N = 2 * Ns
    eta = [rng.gauss(0.0, Dsf) for _ in range(Ns)]
    vals, cols, rowptr = [], [], [0]
    for n in range(N):
        site, s = n // 2, n % 2
        x, y = site % L, site // L
        row = {}
        for j in (y * L + (x + 1) % L, y * L + (x - 1) % L,
                  ((y + 1) % L) * L + x, ((y - 1) % L) * L + x):
            row[2 * j + s] = row.get(2 * j + s, 0.0) - t      # spin-diagonal hop
        flip = 2 * site + (1 - s)                              # the opposite spin
        row[flip] = row.get(flip, 0.0) + eta[site]            # eta (real) -> sigma_x field
        for c in sorted(row):
            cols.append(c)
            vals.append(f"{row[c]:.10g} 0")
        rowptr.append(len(cols))
    _csr(path, N, vals, cols, rowptr)


def _write_spin(opdir, label, L):
    Ns = L * L
    N = 2 * Ns
    # S_z = diag(+1 up, -1 dn)
    vals, cols, rowptr = [], [], [0]
    for n in range(N):
        cols.append(n)
        vals.append("1 0" if n % 2 == 0 else "-1 0")
        rowptr.append(len(cols))
    _csr(os.path.join(opdir, f"{label}.SZ.CSR"), N, vals, cols, rowptr)

    # S_x = site sigma_x (up<->dn), S_y = site sigma_y
    for name, up_dn, dn_up in (("SX", "1 0", "1 0"), ("SY", "0 -1", "0 1")):
        vals, cols, rowptr = [], [], [0]
        for n in range(N):
            site, s = n // 2, n % 2
            partner = 2 * site + (1 - s)
            cols.append(partner)
            vals.append(up_dn if s == 0 else dn_up)
            rowptr.append(len(cols))
        _csr(os.path.join(opdir, f"{label}.{name}.CSR"), N, vals, cols, rowptr)

    # P_z = 1/2(1 + S_z) = diag(1 up, 0 dn): keep only the up orbitals
    vals, cols, rowptr = [], [], [0]
    for n in range(N):
        if n % 2 == 0:
            cols.append(n)
            vals.append("1 0")
        rowptr.append(len(cols))
    _csr(os.path.join(opdir, f"{label}.PZ.CSR"), N, vals, cols, rowptr)

    # P_x = 1/2(1 + S_x): site block 1/2[[1,1],[1,1]]
    vals, cols, rowptr = [], [], [0]
    for n in range(N):
        site = n // 2
        for c in (2 * site, 2 * site + 1):
            cols.append(c)
            vals.append("0.5 0")
        rowptr.append(len(cols))
    _csr(os.path.join(opdir, f"{label}.PX.CSR"), N, vals, cols, rowptr)


def _write_bounds(path, Dsf):
    edge = 4.0 + abs(Dsf) + 0.2
    with open(path, "w") as f:
        f.write(f"{-edge:g} {edge:g}\n")


def _write_desc(path, label, L, Dsf, seed):
    edge = 4.0 + abs(Dsf) + 0.2
    with open(path, "w") as f:
        f.write("# LSQUANT operator descriptor (physical sidecar)\n")
        f.write("observable: hamiltonian\nunits: eV\n")
        f.write(f"system_label: {label}\nsystem_size: {2 * L * L}\n")
        f.write("dimensions: 2\nperiodic: true\nspin: spinful\n")
        f.write(f"spin_flip_Dsf_eV: {Dsf:g}\ndisorder_seed: {seed}\n")
        f.write("has_bounds: true\n")
        f.write(f"spectral_min: {-edge:g}\nspectral_max: {edge:g}\n")
        f.write("provenance: make_ey.py: 2D square + random onsite spin-flip\n")


def build_ey(L, Dsf, seed, outdir="."):
    label = f"ey_L{L}_D{Dsf:g}"
    opdir = os.path.join(outdir, "operators")
    os.makedirs(opdir, exist_ok=True)
    _write_ham(os.path.join(opdir, f"{label}.HAM.CSR"), L, Dsf, seed)
    _write_spin(opdir, label, L)
    _write_bounds(os.path.join(outdir, "BOUNDS"), Dsf)
    _write_desc(os.path.join(opdir, f"{label}.HAM.desc"), label, L, Dsf, seed)
    return label


def main():
    if len(sys.argv) not in (3, 4):
        print("usage: python make_ey.py <L> <Dsf> [seed]")
        sys.exit(1)
    L, Dsf = int(sys.argv[1]), float(sys.argv[2])
    seed = int(sys.argv[3]) if len(sys.argv) == 4 else 1
    print("built", build_ey(L, Dsf, seed), f": {2*L*L} orbitals, Dsf={Dsf:g}, seed={seed}")


if __name__ == "__main__":
    main()
