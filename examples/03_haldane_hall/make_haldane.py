#!/usr/bin/env python3
"""Stage the committed Haldane operators for this tutorial.

The Haldane operators are not generated here: they are the verified Chern-insulator
model committed under test/golden/haldane (t1 = -1, t2 = 0.15, phi = pi/2, C = +1). This
script copies them into the example folder and writes the two sidecar files the run needs:

    operators/haldane.HAM.CSR   the sparse Hamiltonian (copied from the golden)
    operators/haldane.VX.CSR    the x velocity operator (copied)
    operators/haldane.VY.CSR    the y velocity operator (copied)
    operators/haldane.HAM.desc  the physical sidecar (identity + spectral bounds)
    exact                       the exact-trace state set (one local state per basis vector)

The descriptor carries the spectral bounds under the field names the .desc reader parses
(spectral_min, spectral_max), so `lsquant compute` and `lsquant inspect` read them straight
from the sidecar -- this run needs no BOUNDS file. The golden set ships no .desc, so this
script writes it.

    python make_haldane.py
"""

import os
import shutil


GOLDEN = os.path.join("..", "..", "test", "golden", "haldane", "operators")
LABEL = "haldane"
SPECTRAL_MIN, SPECTRAL_MAX = -4.5, 4.5      # enclose the spectrum (|E|max ~ 3)
PHI = 1.5707963267948966                    # pi / 2


def _write_desc(path):
    with open(path, "w") as f:
        f.write("# LSQUANT operator descriptor (physical sidecar)\n")
        f.write("schema: operator_descriptor/v1\n")
        f.write("observable: hamiltonian\n")
        f.write("component: scalar\n")
        f.write("units: eV\n")
        f.write(f"system_label: {LABEL}\n")
        f.write("dimensions: 2\n")
        f.write("hopping_t1_eV: -1\n")
        f.write("hopping_t2_eV: 0.15\n")
        f.write(f"phase_phi: {PHI}\n")
        f.write("has_bounds: true\n")
        f.write(f"spectral_min: {SPECTRAL_MIN:g}\n")    # field names the .desc reader parses
        f.write(f"spectral_max: {SPECTRAL_MAX:g}\n")
        f.write("basis_tags: honeycomb:AB\n")
        f.write("spin: spinless\n")
        f.write("provenance: test/golden/haldane (committed Chern-insulator operators, C=+1); "
                "bounds enclose |E|max~3\n")


def build(outdir="."):
    """Copy the golden operators into outdir/operators and write the sidecars."""
    opdir = os.path.join(outdir, "operators")
    os.makedirs(opdir, exist_ok=True)
    for op in ("HAM", "VX", "VY"):
        shutil.copy(os.path.join(GOLDEN, f"{LABEL}.{op}.CSR"),
                    os.path.join(opdir, f"{LABEL}.{op}.CSR"))
    _write_desc(os.path.join(opdir, f"{LABEL}.HAM.desc"))
    dim = int(open(os.path.join(opdir, f"{LABEL}.HAM.CSR")).readline().split()[0])
    with open(os.path.join(outdir, "exact"), "w") as f:
        f.write("local\n" + f"{dim}\n" + "\n".join(map(str, range(dim))) + "\n")
    return dim


def main():
    dim = build()
    print(f"staged {LABEL}: {dim} sites, t1=-1, t2=0.15, phi=pi/2, "
          f"bounds [{SPECTRAL_MIN:g},{SPECTRAL_MAX:g}] eV")
    print("  operators/haldane.{HAM,VX,VY}.CSR")
    print("  operators/haldane.HAM.desc")
    print("  exact")


if __name__ == "__main__":
    main()
