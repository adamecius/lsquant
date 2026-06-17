#!/usr/bin/env python3
"""Generate a graphene (honeycomb) tight-binding operator at an arbitrary supercell
size, in the CSR text format `lsquant` reads -- for the Phase-F threaded-SpMV scaling
benchmark. Self-contained (no kwant): the official examples/utilities/kwant2wannier
path needs kwant, which is not available offline here.

Model. Honeycomb lattice, 2-atom basis (A,B), L x L cells -> N = 2*L*L sites, periodic.
  - nearest neighbour  : real hopping  t1 = -2.7 eV  (A<->B, 3 bonds/site)
  - next-nearest (Haldane/intrinsic-SOC-like): COMPLEX hopping t2*exp(i*phi)*nu_ij on
    same-sublattice bonds (6 bonds/site), nu = +-1 by chirality -> the operator is a
    genuinely COMPLEX Hermitian matrix (not the real chain), ~9 nonzeros/row.
Hermiticity is exact by construction (each unordered bond is written as (i,j,v) and
(j,i,conj v)); lsquant's ConvertFromCSR fold is then the identity.

Output (under <outdir>/operators and <outdir>/BOUNDS):
  operators/<label>.HAM.CSR   "N nnz" / "re im ..."(2*nnz) / cols(nnz) / rowptr(N+1)
  BOUNDS                       "lo hi"  (Gershgorin enclosure, +2% pad)

Usage:  make_graphene.py <L> [t2] [phi_over_pi] [outdir]
        N = 2*L*L  (L=40 -> 3200 sites; L=400 -> 320k; L=632 -> ~800k)
"""
import sys, os, math, cmath

def build(L, t2=0.10, phi=math.pi/2, t1=-2.7, outdir="."):
    Lx = Ly = int(L)
    N = 2 * Lx * Ly
    def idx(cx, cy, s):
        return 2 * ((cy % Ly) * Lx + (cx % Lx)) + s

    # adjacency: dict per row -> {col: complex value}, summed (handles wrap coincidences)
    rows = [dict() for _ in range(N)]
    def add(i, j, v):
        rows[i][j] = rows[i].get(j, 0j) + v
        rows[j][i] = rows[j].get(i, 0j) + v.conjugate()

    for cx in range(Lx):
        for cy in range(Ly):
            A = idx(cx, cy, 0)
            B = idx(cx, cy, 1)
            # nearest neighbours A(cx,cy) <-> B in cells (0,0),(-1,0),(0,-1)
            add(A, idx(cx,     cy,     1), complex(t1, 0.0))
            add(A, idx(cx - 1, cy,     1), complex(t1, 0.0))
            add(A, idx(cx,     cy - 1, 1), complex(t1, 0.0))
            # next-nearest (same sublattice) Haldane phases; 3 independent +nu vectors
            t2c = t2 * cmath.exp(1j * phi)
            for (dx, dy) in ((1, 0), (-1, 1), (0, -1)):
                add(A, idx(cx + dx, cy + dy, 0),  t2c)            # A: +phase
                add(B, idx(cx + dx, cy + dy, 1),  t2c.conjugate())# B: -phase
    # onsite 0 (clean graphene); ensure a diagonal entry exists for a dense-ish row sum
    nnz = sum(len(r) for r in rows)

    # Gershgorin spectral bound: max over rows of sum |M_ij|
    gmax = 0.0
    for r in rows:
        s = sum(abs(v) for v in r.values())
        if s > gmax:
            gmax = s
    edge = gmax * 1.02   # +2% pad to keep the spectrum strictly inside [-1,1] after rescale

    opdir = os.path.join(outdir, "operators")
    os.makedirs(opdir, exist_ok=True)
    label = f"graphene_L{Lx}"
    vals, cols, rowptr = [], [], [0]
    for i in range(N):
        for j in sorted(rows[i]):
            v = rows[i][j]
            vals.append(f"{v.real:.12g} {v.imag:.12g}")
            cols.append(j)
        rowptr.append(len(cols))
    with open(os.path.join(opdir, f"{label}.HAM.CSR"), "w") as f:
        f.write(f"{N} {len(cols)}\n")
        f.write(" ".join(vals) + " \n")
        f.write(" ".join(map(str, cols)) + " \n")
        f.write(" ".join(map(str, rowptr)) + " \n")
    with open(os.path.join(outdir, "BOUNDS"), "w") as f:
        f.write(f"{-edge:.6f} {edge:.6f}\n")
    return label, N, len(cols), edge

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print(__doc__); sys.exit(1)
    L = int(sys.argv[1])
    t2 = float(sys.argv[2]) if len(sys.argv) > 2 else 0.10
    phiov = float(sys.argv[3]) if len(sys.argv) > 3 else 0.5
    outdir = sys.argv[4] if len(sys.argv) > 4 else "."
    label, N, nnz, edge = build(L, t2=t2, phi=phiov * math.pi, outdir=outdir)
    print(f"built {label}: N={N} sites, nnz={nnz} ({nnz/N:.1f}/row), "
          f"complex Hermitian, bounds=[{-edge:.3f},{edge:.3f}]")
