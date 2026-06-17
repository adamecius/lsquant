#!/usr/bin/env python3
"""Build a 2D square-lattice tight-binding model with box Anderson disorder.

One orbital per site on an L x L periodic lattice, hopping amplitude -t (so the
clean band is E_k = -2t(cos kx + cos ky), t = 1, a = 1), with i.i.d. uniform
onsite energies eps_i in [-W/2, +W/2]:

    H = -t sum_<ij> ( |i><j| + h.c. ) + sum_i eps_i |i><i|.

Onsite disorder is diagonal, so the x-velocity operator stays the clean one and
acts only on x-bonds (y-bonds commute with the x-position):

    v_{i,i+x} = +i,   v_{i,i-x} = -i      (units t = a = 1).

The mean-square-displacement driver reads VX, so DeltaX^2(E,t) is the x-projection
of the spread. Gershgorin puts the spectrum inside |E| <= 4t + W/2, so
BOUNDS = ( -(4 + W/2), +(4 + W/2) ); the descriptor records the same.

Operators are written under operators/ at run time (nothing is committed). The
functions are importable: lsqtransport.py calls build_square with a fresh seed per
disorder realisation.

    python make_square.py <L> <W> [seed]
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


def _write_ham(path, L, W, seed, t=1.0):
    rng = random.Random(seed)
    N = L * L
    eps = [rng.uniform(-0.5 * W, 0.5 * W) for _ in range(N)]
    vals, cols, rowptr = [], [], [0]
    for n in range(N):
        x, y = n % L, n // L
        row = {}
        for c in (y * L + (x + 1) % L, y * L + (x - 1) % L,
                  ((y + 1) % L) * L + x, ((y - 1) % L) * L + x):
            row[c] = row.get(c, 0.0) - t            # hopping amplitude -t
        row[n] = row.get(n, 0.0) + eps[n]           # onsite disorder
        for c in sorted(row):
            cols.append(c)
            vals.append(f"{row[c]:.10g} 0")
        rowptr.append(len(cols))
    _csr(path, N, vals, cols, rowptr)


def _write_vx(path, L):
    N = L * L
    vals, cols, rowptr = [], [], [0]
    for n in range(N):
        x, y = n % L, n // L
        row = {y * L + (x + 1) % L: "0 1", y * L + (x - 1) % L: "-0 -1"}
        for c in sorted(row):
            cols.append(c)
            vals.append(row[c])
        rowptr.append(len(cols))
    _csr(path, N, vals, cols, rowptr)


def _write_bounds(path, W):
    edge = 4.0 + 0.5 * W
    with open(path, "w") as f:
        f.write(f"{-edge:g} {edge:g}\n")


def _write_desc(path, label, L, W, seed):
    edge = 4.0 + 0.5 * W
    with open(path, "w") as f:
        f.write("# LSQUANT operator descriptor (physical sidecar)\n")
        f.write("schema: operator_descriptor/v1\n")
        f.write("observable: hamiltonian\n")
        f.write("units: eV\n")
        f.write(f"system_label: {label}\n")
        f.write(f"system_size: {L * L}\n")
        f.write("dimensions: 2\n")
        f.write("periodic: true\n")
        f.write("hopping_t_eV: 1\n")
        f.write("disorder: anderson_box\n")
        f.write(f"disorder_W_eV: {W:g}\n")
        f.write(f"disorder_seed: {seed}\n")
        f.write("has_bounds: true\n")
        f.write(f"spectral_min: {-edge:g}\n")
        f.write(f"spectral_max: {edge:g}\n")
        f.write("basis_tags: site:s\n")
        f.write("spin: spinless\n")
        f.write("provenance: make_square.py: 2D square Anderson lattice, periodic\n")


def build_square(L, W, seed, outdir="."):
    """Write HAM.CSR, VX.CSR, BOUNDS and the descriptor; return the label."""
    label = f"square_L{L}_W{W:g}"
    opdir = os.path.join(outdir, "operators")
    os.makedirs(opdir, exist_ok=True)
    _write_ham(os.path.join(opdir, f"{label}.HAM.CSR"), L, W, seed)
    _write_vx(os.path.join(opdir, f"{label}.VX.CSR"), L)
    _write_bounds(os.path.join(outdir, "BOUNDS"), W)
    _write_desc(os.path.join(opdir, f"{label}.HAM.desc"), label, L, W, seed)
    return label


def main():
    if len(sys.argv) not in (3, 4):
        print("usage: python make_square.py <L> <W> [seed]")
        sys.exit(1)
    L = int(sys.argv[1])
    W = float(sys.argv[2])
    seed = int(sys.argv[3]) if len(sys.argv) == 4 else 1
    label = build_square(L, W, seed)
    edge = 4.0 + 0.5 * W
    print(f"built {label}: {L}x{L}={L*L} sites, t=1, W={W:g}, seed={seed}, "
          f"band [-{edge:g},{edge:g}]")


if __name__ == "__main__":
    main()
