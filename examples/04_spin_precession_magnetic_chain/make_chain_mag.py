#!/usr/bin/env python3
"""Stage or build the 1D magnetic chain operators for Tutorial 04.

Stages the verified magnetic-chain operators (N = 128 sites = 64 spinful unit
cells, hopping gamma = 1 eV, exchange J_ex = 0.1 eV along z) from
goldens/stage0/chain1d_mag/operators into <outdir>/operators:

    operators/chain1d_mag.HAM.CSR   exchange-split 1D chain Hamiltonian
    operators/chain1d_mag.SX.CSR    spin operator S_x
    operators/chain1d_mag.SY.CSR    spin operator S_y
    operators/chain1d_mag.SZ.CSR    spin operator S_z
    operators/chain1d_mag.PX.CSR    initial-state projector P_x = (1 + S_x) / 2

Usage:
    python make_chain_mag.py
"""

from __future__ import annotations

import os
import shutil

_REPO_ROOT = os.path.dirname(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
)
GOLDEN = os.path.join(_REPO_ROOT, "goldens", "stage0", "chain1d_mag", "operators")
LABEL = "chain1d_mag"
OPS = ("HAM", "SX", "SY", "SZ", "PX")


def write_operators(outdir: str = ".") -> str:
    """Copy the committed magnetic-chain operators and sidecars into outdir/operators."""
    opdir = os.path.join(outdir, "operators")
    os.makedirs(opdir, exist_ok=True)
    for op in OPS:
        src_csr = os.path.join(GOLDEN, f"{LABEL}.{op}.CSR")
        dst_csr = os.path.join(opdir, f"{LABEL}.{op}.CSR")
        if os.path.abspath(src_csr) != os.path.abspath(dst_csr):
            shutil.copyfile(src_csr, dst_csr)
        src_desc = os.path.join(GOLDEN, f"{LABEL}.{op}.desc")
        dst_desc = os.path.join(opdir, f"{LABEL}.{op}.desc")
        if os.path.isfile(src_desc) and os.path.abspath(src_desc) != os.path.abspath(dst_desc):
            shutil.copyfile(src_desc, dst_desc)
    return LABEL


def build_chain(outdir: str = ".") -> str:
    """Alias for write_operators(outdir)."""
    return write_operators(outdir)


def build(outdir: str = ".") -> int:
    """Stage operators into outdir/operators and return matrix dimension."""
    write_operators(outdir)
    ham_path = os.path.join(outdir, "operators", f"{LABEL}.HAM.CSR")
    with open(ham_path, encoding="utf-8") as fh:
        return int(fh.readline().split()[0])


def main() -> None:
    dim = build(".")
    print(f"staged {LABEL}: dim={dim} (64 sites x 2 spins), J_ex=0.1 eV")
    print("  operators/chain1d_mag.{HAM,SX,SY,SZ,PX}.CSR")


if __name__ == "__main__":
    main()
