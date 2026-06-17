#!/usr/bin/env python3
"""
siesta2lsquant_realspace.py

Reads a SIESTA .fdf file, builds a real-space supercell Hamiltonian of
size Nx×Ny unit cells, Löwdin-orthogonalises it, and writes lsquant .CSR files:

    <prefix>.HAM   – Löwdin H_ortho  =  S^{-1/2} H S^{-1/2}
    <prefix>.VX    – velocity Vx = i [H_ortho, X_ortho]
    <prefix>.VY    – velocity Vy = i [H_ortho, Y_ortho]

where X_ortho = S^{-1/2} X S^{-1/2}.

How sisl stores the Hamiltonian in real space
---------------------------------------------
sisl's internal format is a rectangular sparse matrix of shape (no, no*n_s),
where no = orbitals per unit cell and n_s = number of supercell images (R-vectors).
This is NOT directly usable as a square Hamiltonian matrix.

The correct real-space supercell square matrix for an Nx×Ny tiling is obtained by:
    H_super = H.tile(Nx, axis=0).tile(Ny, axis=1)
    H_mat   = H_super.Hk(k=[0,0,0])   # sums all R-images → square (N,N) sparse
    S_mat   = H_super.Sk(k=[0,0,0])

After tiling, all inter-cell hoppings become intra-supercell, so evaluating at
k=0 (which computes sum_R H(R) e^{i·0·R} = sum_R H(R)) gives the full real-space
matrix with no Fourier approximation whatsoever.

Löwdin orthogonalisation
------------------------
    S = U Λ U†   (eigh decomposition)
    S^{-1/2} = U Λ^{-1/2} U†
    H_ortho = S^{-1/2} H S^{-1/2}

S^{-1/2} is dense in general → stored as a dense (N,N) array.
Complexity: O(N³) for the diagonalisation.

Position operator
-----------------
sisl does not expose ⟨μ|r̂|ν⟩ matrix elements directly.  A stub is provided:
  - Default: diagonal approximation  X[μ,μ] = x_μ  (orbital centre coordinate)
  - Replace build_position_operator() with your own off-diagonal implementation.

Usage
-----
    python siesta2lsquant_realspace.py <fdf_file> <Nx> <Ny>
        [--prefix NAME] [--no-vy] [--no-vel] [--eps EPS]

Dependencies: sisl scipy numpy
"""

import argparse
import sys
import time
from pathlib import Path

import numpy as np
import scipy.sparse as sp


# ─────────────────────────────────────────────────────────────────────────────
# lsquant .CSR writer
# ─────────────────────────────────────────────────────────────────────────────

def write_lsquant_csr(mat, path: str):
    """
    Write a sparse or dense matrix as a lsquant .CSR text file.

    Format (4 lines):
        Ndim  nnz
        re0 im0 re1 im1 ...    (2*nnz interleaved floats)
        col0 col1 ...          (nnz column indices)
        ptr0 ptr1 ...          (Ndim+1 row pointers)
    """
    if not sp.issparse(mat):
        mat = sp.csr_matrix(mat)
    mat = mat.tocsr()
    mat.sort_indices()
    mat.eliminate_zeros()

    Ndim = mat.shape[1]
    nnz  = mat.nnz
    data = mat.data.astype(complex)

    interleaved      = np.empty(2 * nnz, dtype=np.float64)
    interleaved[0::2] = data.real
    interleaved[1::2] = data.imag

    with open(path, 'w') as f:
        f.write(f"{Ndim} {nnz}\n")
        f.write(' '.join(f'{v:.22f}' for v in interleaved)); f.write('\n')
        f.write(' '.join(map(str, mat.indices)));            f.write('\n')
        f.write(' '.join(map(str, mat.indptr)));             f.write('\n')

    print(f"  Written: {path}  (Ndim={Ndim}, nnz={nnz})")


# ─────────────────────────────────────────────────────────────────────────────
# S^{-1/2} via eigendecomposition
# ─────────────────────────────────────────────────────────────────────────────

def sqrt_inv_dense(S_dense: np.ndarray, eps: float = 1e-10) -> np.ndarray:
    """
    S^{-1/2} for Hermitian positive-definite S, via eigh.

    Eigenvalues below eps are clamped to avoid divide-by-zero for
    near-linearly-dependent basis sets.
    """
    lam, U = np.linalg.eigh(S_dense)
    lam    = np.maximum(lam, eps)
    return (U * (lam ** -0.5)[np.newaxis, :]) @ U.conj().T


# ─────────────────────────────────────────────────────────────────────────────
# Position operator  (USER STUB — replace off-diagonal part as needed)
# ─────────────────────────────────────────────────────────────────────────────

def build_position_operator(geometry, axis: int) -> sp.csr_matrix:
    """
    Real-space position operator X^axis in the orbital basis, shape (N, N).

    Default: diagonal approximation  X[μ,μ] = r_μ^axis
    where r_μ is the Cartesian coordinate of the atom hosting orbital μ.

    This is exact for atom-centred Wannier functions (Peierls approximation).
    For a full implementation, replace the diagonal assignment below with
    the computed ⟨φ_μ(r)|r^axis|φ_ν(r-R)⟩ matrix elements.

    Parameters
    ----------
    geometry : sisl.Geometry   (supercell)
    axis     : 0=x, 1=y, 2=z
    """
    N        = geometry.no
    atom_idx = np.array([geometry.o2a(io) for io in range(N)])
    coords   = geometry.xyz[atom_idx, axis].astype(complex)   # (N,)

    # ── Diagonal approximation (default) ─────────────────────────────────────
    X = sp.diags(coords, format='csr', dtype=complex)

    # ── Replace above with full off-diagonal implementation if needed ─────────
    # rows, cols, vals = [], [], []
    # for mu in range(N):
    #     for nu, R, x_munu in your_integral_routine(geometry, mu, axis):
    #         rows.append(mu); cols.append(nu); vals.append(x_munu)
    # X = sp.csr_matrix((vals, (rows, cols)), shape=(N, N), dtype=complex)
    # ─────────────────────────────────────────────────────────────────────────

    return X


# ─────────────────────────────────────────────────────────────────────────────
# Main builder
# ─────────────────────────────────────────────────────────────────────────────

def build_realspace_operators(H_source, Nx: int, Ny: int,
                              eps: float    = 1e-10,
                              compute_vel: bool = True):
    """
    Build H_ortho, Vx, Vy for an Nx×Ny supercell.

    Steps
    -----
    1. Tile:   H_super = H.tile(Nx, 0).tile(Ny, 1)
    2. Extract real-space square matrices via Hk/Sk at k=0
       (correct because all hoppings are intra-supercell after tiling)
    3. Löwdin: S^{-1/2} via eigh, H_ortho = S^{-1/2} H S^{-1/2}
    4. Velocity: V^α = i [H_ortho, X_ortho]
                 with X_ortho = S^{-1/2} X S^{-1/2}

    Returns
    -------
    H_ortho : sp.csr_matrix
    Vx      : sp.csr_matrix  (or None)
    Vy      : sp.csr_matrix  (or None)
    """

    # ── Step 1: tile ──────────────────────────────────────────────────────────
    print(f"  Tiling unit cell ({H_source.no} orb) → {Nx}×{Ny} supercell ...",
          end=' ', flush=True)
    t0 = time.time()
    H_super = H_source.tile(Nx, axis=0).tile(Ny, axis=1)
    N       = H_super.no
    print(f"done ({time.time()-t0:.1f}s)  N={N} orbitals")

    # ── Step 2: extract square real-space H and S at k=0 ─────────────────────
    # After tiling, evaluating at k=[0,0,0] gives sum_R H(R) = H_realspace
    # because all inter-cell R vectors are now zero (intra-supercell).
    print("  Extracting H and S at k=0 (real-space square matrices) ...",
          end=' ', flush=True)
    t0 = time.time()
    H_sp = H_super.Hk(k=[0, 0, 0], dtype=np.complex128).tocsr()
    S_sp = H_super.Sk(k=[0, 0, 0], dtype=np.complex128).tocsr()
    print(f"done ({time.time()-t0:.1f}s)"
          f"  H nnz={H_sp.nnz}  S nnz={S_sp.nnz}")

    # ── Step 3: S^{-1/2} and H_ortho ─────────────────────────────────────────
    print(f"  Diagonalising S ({N}×{N}) ...", end=' ', flush=True)
    t0    = time.time()
    Sinvh = sqrt_inv_dense(S_sp.toarray(), eps=eps)     # (N,N) dense
    print(f"done ({time.time()-t0:.1f}s)")

    print("  Computing H_ortho = S^{-1/2} H S^{-1/2} ...", end=' ', flush=True)
    t0 = time.time()
    H_dense  = H_sp.toarray()
    Ho_dense = Sinvh @ H_dense @ Sinvh
    H_ortho  = sp.csr_matrix(Ho_dense)
    H_ortho.eliminate_zeros()
    print(f"done ({time.time()-t0:.1f}s)  nnz={H_ortho.nnz}")

    if not compute_vel:
        return H_ortho, None, None

    # ── Step 4: position operators and velocities ─────────────────────────────
    geom = H_super.geometry
    Vx = Vy = None

    for axis, label in [(0, 'x'), (1, 'y')]:
        print(f"  Building X_{label} ...", end=' ', flush=True)
        t0 = time.time()
        X_sp    = build_position_operator(geom, axis)
        X_dense = X_sp.toarray() if sp.issparse(X_sp) else np.asarray(X_sp)
        print(f"done ({time.time()-t0:.2f}s)")

        print(f"  Computing X_{label}_ortho and V_{label} = i[H, X]_ortho ...",
              end=' ', flush=True)
        t0      = time.time()
        Xo      = Sinvh @ X_dense @ Sinvh          # X_ortho  (N,N) dense
        comm    = Ho_dense @ Xo - Xo @ Ho_dense    # [H_ortho, X_ortho]
        V_sp    = sp.csr_matrix(1j * comm)
        V_sp.eliminate_zeros()
        print(f"done ({time.time()-t0:.1f}s)  nnz={V_sp.nnz}")

        if axis == 0:
            Vx = V_sp
        else:
            Vy = V_sp

    return H_ortho, Vx, Vy


# ─────────────────────────────────────────────────────────────────────────────
# Entry point
# ─────────────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Real-space Löwdin orthogonalisation from SIESTA → lsquant .CSR"
    )
    parser.add_argument("fdf_file",     help="SIESTA .fdf file")
    parser.add_argument("Nx", type=int, help="Supercell tiling along a1")
    parser.add_argument("Ny", type=int, help="Supercell tiling along a2")
    parser.add_argument("--prefix", default=None,
                        help="Output prefix (default: <stem>_NxNy)")
    parser.add_argument("--no-vy",  action="store_true", help="Skip Vy")
    parser.add_argument("--no-vel", action="store_true", help="Skip all velocities")
    parser.add_argument("--eps",    type=float, default=1e-10,
                        help="S^{-1/2} eigenvalue clamp threshold (default: 1e-10)")
    args = parser.parse_args()

    fdf_path = Path(args.fdf_file)
    if not fdf_path.exists():
        sys.exit(f"Error: not found: {fdf_path}")

    prefix = args.prefix or f"{fdf_path.stem}_{args.Nx}x{args.Ny}"

    print(f"Reading {fdf_path} ...")
    try:
        import sisl
        H_src = sisl.get_sile(fdf_path).read_hamiltonian()
    except Exception as e:
        sys.exit(f"sisl read error: {e}")

    print(f"  Unit cell orbitals : {H_src.no}")
    print(f"  Supercell images   : {H_src.n_s}  (R-vectors in _csr)")
    print(f"  Spin               : {H_src.spin}\n")

    H_ortho, Vx, Vy = build_realspace_operators(
        H_src, args.Nx, args.Ny,
        eps=args.eps,
        compute_vel=not args.no_vel,
    )

    print("\nWriting CSR files ...")
    write_lsquant_csr(H_ortho, f"{prefix}.HAM.CSR")
    if Vx is not None:
        write_lsquant_csr(Vx, f"{prefix}.VX.CSR")
    if Vy is not None and not args.no_vy:
        write_lsquant_csr(Vy, f"{prefix}.VY.CSR")

    print("\nAll done.")


if __name__ == "__main__":
    main()
