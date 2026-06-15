#!/usr/bin/env python3
"""Build a 1D Anderson chain (box onsite disorder) and its velocity operator.

Hamiltonian, one orbital per site, periodic, hopping t = -1 eV, with i.i.d.
uniform onsite energies:

    H = t * sum_i ( |i><i+1| + |i+1><i| ) + sum_i eps_i |i><i|,
    eps_i uniform in [-W/2, +W/2].

Onsite disorder sits on the diagonal, so it commutes with the position, and the
velocity operator stays the clean one,

    v_{i,i+1} = +i,   v_{i,i-1} = -i      (units t = a = 1).

This VX file is byte-identical to the committed chain1d_ballistic golden, so the
mean-square-displacement driver reads it directly.

Gershgorin puts the spectrum inside |E| <= 2|t| + W/2, so the bounds widen with
disorder: BOUNDS = ( -(2 + W/2), +(2 + W/2) ).

The functions are importable (the lsqloc wrapper calls build_chain with a fresh
seed per disorder realisation). The CLI builds one labelled set in place.

    python make_disordered_chain.py <N> <W> [seed]
"""

import os
import sys
import hashlib


def _write_ham(path, N, W, seed, t=-1.0):
    import random
    rng = random.Random(seed)
    eps = [rng.uniform(-0.5 * W, 0.5 * W) for _ in range(N)]

    vals, cols, rowptr = [], [], [0]
    for i in range(N):
        left, right = (i - 1) % N, (i + 1) % N
        entries = sorted([(left, t), (i, eps[i]), (right, t)])
        for c, v in entries:
            cols.append(c)
            vals.append(f"{v:.10g} 0")
        rowptr.append(len(cols))

    with open(path, "w") as f:
        f.write(f"{N} {3 * N}\n")
        f.write(" ".join(vals) + " \n")
        f.write(" ".join(map(str, cols)) + " \n")
        f.write(" ".join(map(str, rowptr)) + " \n")


def _write_vx(path, N):
    # v_{i,i+1} = +i ("0 1"), v_{i,i-1} = -i ("-0 -1"); matches the golden VX.
    vals, cols, rowptr = [], [], [0]
    for i in range(N):
        left, right = (i - 1) % N, (i + 1) % N
        entries = sorted([(left, "-0 -1"), (right, "0 1")])
        for c, s in entries:
            cols.append(c)
            vals.append(s)
        rowptr.append(len(cols))
    with open(path, "w") as f:
        f.write(f"{N} {2 * N}\n")
        f.write(" ".join(vals) + " \n")
        f.write(" ".join(map(str, cols)) + " \n")
        f.write(" ".join(map(str, rowptr)) + " \n")


def _write_bounds(path, W):
    edge = 2.0 + 0.5 * W
    with open(path, "w") as f:
        f.write(f"{-edge:g} {edge:g}\n")


def _write_desc(path, label, N, W, seed, t=-1.0, csr_path=None):
    edge = 2.0 + 0.5 * W
    csr_hash = ""
    if csr_path and os.path.exists(csr_path):
        with open(csr_path, "rb") as fh:
            csr_hash = hashlib.sha256(fh.read()).hexdigest()[:16]
    with open(path, "w") as f:
        f.write("# LinQT operator descriptor (physical sidecar)\n")
        f.write("schema: operator_descriptor/v1\n")
        f.write("observable: hamiltonian\n")
        f.write("units: eV\n")
        f.write(f"system_label: {label}\n")
        f.write(f"system_size: {N}\n")
        f.write("dimensions: 1\n")
        f.write("periodic: true\n")
        f.write(f"hopping_t_eV: {t:g}\n")
        f.write("disorder: anderson_box\n")
        f.write(f"disorder_W_eV: {W:g}\n")
        f.write(f"disorder_seed: {seed}\n")
        f.write("has_bounds: true\n")
        f.write(f"band_min_eV: {-edge:g}\n")
        f.write(f"band_max_eV: {edge:g}\n")
        f.write("basis_tags: site:s\n")
        f.write("spin: spinless\n")
        f.write("provenance: make_disordered_chain.py: 1D Anderson chain, periodic\n")
        if csr_hash:
            f.write(f"csr_sha256_16: {csr_hash}\n")


def build_chain(N, W, seed, outdir="."):
    """Write HAM.CSR, VX.CSR, BOUNDS, and the descriptor into outdir."""
    label = f"chain1d_dis_N{N}_W{W:g}"
    opdir = os.path.join(outdir, "operators")
    os.makedirs(opdir, exist_ok=True)
    ham = os.path.join(opdir, f"{label}.HAM.CSR")
    vx = os.path.join(opdir, f"{label}.VX.CSR")
    _write_ham(ham, N, W, seed)
    _write_vx(vx, N)
    _write_bounds(os.path.join(outdir, "BOUNDS"), W)
    _write_desc(os.path.join(opdir, f"{label}.HAM.desc"),
                label, N, W, seed, csr_path=ham)
    return label


def main():
    if len(sys.argv) not in (3, 4):
        print("usage: python make_disordered_chain.py <N> <W> [seed]")
        sys.exit(1)
    N = int(sys.argv[1])
    W = float(sys.argv[2])
    seed = int(sys.argv[3]) if len(sys.argv) == 4 else 1
    label = build_chain(N, W, seed)
    print(f"built {label}: {N} sites, t=-1 eV, W={W:g} eV, seed={seed}, "
          f"band [-{2 + 0.5 * W:g},{2 + 0.5 * W:g}] eV")


if __name__ == "__main__":
    main()
