#!/usr/bin/env python3
"""Build a 1D nearest-neighbour tight-binding chain for LinQT.

The chain has N sites, one orbital per site, periodic boundaries, and a single
hopping t. The Hamiltonian is

    H = t * sum_i ( |i><i+1| + |i+1><i| ),     i+1 taken modulo N,

with eigenvalues E_k = 2 t cos(2 pi k / N). The spectrum fills the band
[-2|t|, 2|t|], so the spectral bounds are (a, b) = (-2|t|, 2|t|).

For each size this script writes, inside the current example folder:

    operators/<label>.HAM.CSR   the sparse Hamiltonian (LinQT complex-CSR text)
    operators/<label>.HAM.desc  the physical sidecar (provenance + identity + bounds)
    BOUNDS                      the (a b) line read by the compute driver
    exact_<label>               the exact-trace state set (one local state per basis vector)

The sidecar carries everything that defines the operator as a physical object
(units, observable identity, basis tags, spectral bounds, provenance). The
numerical CSR carries the matrix. They live side by side. The bounds use the
field names the .desc reader parses (`spectral_min`, `spectral_max`), so
`lsquant inspect operators/<label>.HAM.desc` reports them.

The spectral driver this tutorial uses (`inline_compute-kpm-spectralOp`) reads
its bounds from the BOUNDS file, so this script writes BOUNDS as well; the
descriptor travels alongside as the physical record.

Usage:
    python make_chain.py 512
    python make_chain.py 2048
    python make_chain.py 8192
"""

import os
import sys
import hashlib


def write_csr(path, N, t=-1.0):
    """Write the chain Hamiltonian in LinQT complex-CSR text format.

    Format (four lines):
        line 1: DIM NNZ
        line 2: values as 're im' pairs, in row-major CSR order
        line 3: column indices, in the same order (sorted ascending per row)
        line 4: row pointers (N+1 integers)
    """
    nnz = 2 * N
    cols = []
    rowptr = [0]
    for i in range(N):
        left = (i - 1) % N
        right = (i + 1) % N
        cols.extend(sorted((left, right)))
        rowptr.append(len(cols))

    with open(path, "w") as f:
        f.write(f"{N} {nnz}\n")
        f.write(" ".join(f"{t:g} 0" for _ in range(nnz)) + "\n")
        f.write(" ".join(str(c) for c in cols) + "\n")
        f.write(" ".join(str(r) for r in rowptr) + "\n")


def write_desc(path, label, N, t=-1.0, csr_path=None):
    """Write the human-authored physical sidecar for a manually built model.

    The fields match LinQT's OperatorDescriptor: observable identity, component,
    units, spectral bounds (a, b), basis/orbital tags, and provenance. The
    bounds here are exact for the analytic chain; for a generic model they come
    from the producer (wannier2sparse via Lanczos).
    """
    a, b = -2.0 * abs(t), 2.0 * abs(t)
    provenance = "make_chain.py: 1D NN chain, periodic, one orbital/site"
    csr_hash = ""
    if csr_path and os.path.exists(csr_path):
        with open(csr_path, "rb") as fh:
            csr_hash = hashlib.sha256(fh.read()).hexdigest()[:16]

    with open(path, "w") as f:
        f.write("# LinQT operator descriptor (physical sidecar)\n")
        f.write("schema: operator_descriptor/v1\n")
        f.write("observable: hamiltonian\n")
        f.write("component: scalar\n")
        f.write("units: eV\n")
        f.write(f"system_label: {label}\n")
        f.write(f"system_size: {N}\n")
        f.write("dimensions: 1\n")
        f.write("periodic: true\n")
        f.write(f"hopping_t_eV: {t:g}\n")
        f.write("has_bounds: true\n")
        f.write(f"spectral_min: {a:g}\n")          # field names the .desc reader parses
        f.write(f"spectral_max: {b:g}\n")
        f.write(f"band_width_eV: {b - a:g}\n")
        f.write(f"band_center_eV: {0.5 * (a + b):g}\n")
        f.write("basis_tags: site:s\n")          # one s-like orbital per site
        f.write("spin: spinless\n")
        f.write(f"provenance: {provenance}\n")
        if csr_hash:
            f.write(f"csr_sha256_16: {csr_hash}\n")


def write_bounds(path, t=-1.0):
    with open(path, "w") as f:
        f.write(f"{-2.0 * abs(t):g} {2.0 * abs(t):g}\n")


def write_state(path, N):
    """Exact-trace state set: the local states e_0..e_{N-1} over the full basis.

    LinQT evaluates sum_i <e_i|T_m|e_i>/N = Tr[T_m]/N exactly, so the moments are
    deterministic and equal to the analytic trace, with no random vector.
    """
    with open(path, "w") as f:
        f.write("local\n")
        f.write(f"{N}\n")
        f.write("\n".join(str(i) for i in range(N)) + "\n")


def main():
    if len(sys.argv) != 2:
        print("usage: python make_chain.py <N>")
        sys.exit(1)
    N = int(sys.argv[1])
    t = -1.0
    label = f"chain1d_N{N}"

    os.makedirs("operators", exist_ok=True)
    csr = f"operators/{label}.HAM.CSR"
    desc = f"operators/{label}.HAM.desc"
    state = f"exact_{label}"

    write_csr(csr, N, t)
    write_desc(desc, label, N, t, csr_path=csr)
    write_bounds("BOUNDS", t)
    write_state(state, N)

    print(f"built {label}: {N} sites, t={t:g} eV, band [-2,2] eV")
    print(f"  {csr}")
    print(f"  {desc}")
    print(f"  {state}")
    print("  BOUNDS")


if __name__ == "__main__":
    main()
