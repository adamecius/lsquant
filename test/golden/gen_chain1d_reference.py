#!/usr/bin/env python3
"""
Generate the committed analytic reference files FROM the verified oracle
(lsquant_chain_reference.py) — the single source of truth. The C++ golden tests
read the generated .txt files, so no Python is needed to RUN ctest; Python is only
needed here, when (re)generating the references after an oracle change.

Run via the project venv:
    .venv/bin/python test/golden/gen_chain1d_reference.py
"""
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))   # repo root (oracle lives there)
sys.path.insert(0, ROOT)
import lsquant_chain_reference as oracle   # noqa: E402

M = 9          # leading (M+1)x(M+1) KG block compared by the test
NDOS = 15      # DOS moments m = 0..NDOS

kg_path  = os.path.join(HERE, "chain1d", "KG_reference.txt")
dos_path = os.path.join(HERE, "chain1d", "DOS_reference.txt")

KG = oracle.kg_moment_matrix(M)            # gamma=a=hbar=1 -> B.15 normalized structural matrix
with open(kg_path, "w") as f:
    f.write("# Kubo-Greenwood Chebyshev moment matrix M[m,n] for the 1D chain.\n")
    f.write("# GENERATED from lsquant_chain_reference.kg_moment_matrix (B.15 normalized).\n")
    f.write("# Do NOT hand-edit; oracle proven in lsquant_derive_core.py. Leading 10x10:\n")
    for row in KG:
        f.write(" ".join(f"{x: .10g}" for x in row) + "\n")

with open(dos_path, "w") as f:
    f.write("# DOS Chebyshev moments Cbar^DOS_m = delta_{m,0} for the 1D chain.\n")
    f.write("# GENERATED from lsquant_chain_reference.dos_moment; do NOT hand-edit.\n")
    for m in range(NDOS + 1):
        f.write(f"{oracle.dos_moment(m): .10g}\n")

print(f"wrote {os.path.relpath(kg_path, ROOT)} and {os.path.relpath(dos_path, ROOT)} from the oracle")
