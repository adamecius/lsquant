#!/usr/bin/env python3
"""Stage the committed Haldane operators for Tutorial 03.

Copies the verified Chern-insulator model (t1 = -1, t2 = 0.15, phi = pi/2, C = +1)
committed under goldens/stage0/haldane/operators into <outdir>/operators and writes
the operator descriptor sidecar and exact-trace state file:

    operators/haldane.HAM.CSR   sparse Hamiltonian
    operators/haldane.VX.CSR    x-velocity operator
    operators/haldane.VY.CSR    y-velocity operator
    operators/haldane.HAM.desc  operator descriptor sidecar (identity + spectral bounds)
    exact                       exact-trace state set

Usage:
    python make_haldane.py
"""

import os
import shutil

_REPO_ROOT = os.path.dirname(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
)
GOLDEN = os.path.join(_REPO_ROOT, "goldens", "stage0", "haldane", "operators")
LABEL = "haldane"
SPECTRAL_MIN, SPECTRAL_MAX = -4.5, 4.5
PHI = 1.5707963267948966


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
        f.write(f"spectral_min: {SPECTRAL_MIN:g}\n")
        f.write(f"spectral_max: {SPECTRAL_MAX:g}\n")
        f.write("basis_tags: honeycomb:AB\n")
        f.write("spin: spinless\n")
        f.write("provenance: goldens/stage0/haldane (committed Chern-insulator operators, C=+1); "
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
