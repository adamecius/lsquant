#!/usr/bin/env python3
"""
siesta2lsquant.py

Reads a SIESTA .fdf file, builds a k-point grid of size kx × ky × 1,
performs Löwdin orthogonalisation at each k-point, assembles the full
supercell (block-diagonal) operators, and writes them as lsquant .CSR files:

    <prefix>.HAM   – Löwdin-orthogonalised Hamiltonian  H_ortho
    <prefix>.VX    – velocity operator V_x
    <prefix>.VY    – velocity operator V_y

Löwdin velocity formula (from the non-orthogonal Kubo–Bastin literature):

    H_ortho  = S^{-1/2} H S^{-1/2}

    V^α = S^{-1/2} ∂_α H S^{-1/2}
          - ½ H_ortho (∂_α S) S^{-1}
          - ½ (∂_α S) S^{-1} H_ortho

where ∂_α ≡ ∂/∂k_α  (sisl's dHk / dSk, gauge='r').

S^{-1/2} is computed via eigendecomposition of S (Hermitian positive-definite):
    S = U Λ U†  →  S^{-1/2} = U Λ^{-1/2} U†
This is numerically more stable than sqrtm(inv(S)).

Usage
-----
    python siesta2lsquant.py <fdf_file> <kx> <ky> [--prefix NAME]
                           [--kz KZ] [--no-vy] [--gauge GAUGE]

    fdf_file   path to SIESTA RUN.fdf (or equivalent)
    kx, ky     k-point grid dimensions
    --prefix   output file prefix  (default: derived from fdf filename)
    --kz       k-points along z    (default: 1)
    --no-vy    skip V_y output
    --gauge    sisl gauge string   (default: 'r')

Dependencies
------------
    pip install sisl scipy numpy
"""

import argparse
import sys
import time
from pathlib import Path

import numpy as np
import scipy.linalg as sla
from scipy.sparse import csr_matrix, coo_matrix, block_diag


# ─────────────────────────────────────────────────────────────────────────────
# CSR writer  (lsquant format)
# ─────────────────────────────────────────────────────────────────────────────

def write_lsquant_csr(mat, path: str):
    """
    Write a scipy sparse matrix to a lsquant .CSR text file.

    Format (4 lines):
        Ndim  nnz
        re0 im0 re1 im1 ...    (2*nnz floats, interleaved)
        col0 col1 ...          (nnz int column indices)
        ptr0 ptr1 ...          (Ndim+1 int row pointers)
    """
    mat = mat.tocsr()
    mat.sort_indices()

    Ndim = mat.shape[1]
    nnz  = mat.nnz
    data = mat.data  # complex128

    with open(path, 'w') as f:
        f.write(f"{Ndim} {nnz}\n")

        # Interleaved real / imag
        interleaved = np.empty(2 * nnz, dtype=np.float64)
        interleaved[0::2] = data.real
        interleaved[1::2] = data.imag
        f.write(' '.join(f'{v:.22f}' for v in interleaved))
        f.write('\n')

        f.write(' '.join(map(str, mat.indices)))
        f.write('\n')
        f.write(' '.join(map(str, mat.indptr)))
        f.write('\n')

    print(f"  Written: {path}  (Ndim={Ndim}, nnz={nnz})")


# ─────────────────────────────────────────────────────────────────────────────
# Stable S^{-1/2} via eigendecomposition
# ─────────────────────────────────────────────────────────────────────────────

def sqrt_inv(S: np.ndarray, eps: float = 1e-12) -> np.ndarray:
    """
    Compute S^{-1/2} for a Hermitian positive-definite matrix S.

    Uses eigh (stable for Hermitian matrices) rather than sqrtm(inv(S)),
    which can accumulate errors when eigenvalues are small.

    S = U Λ U†  →  S^{-1/2} = U diag(λ^{-1/2}) U†
    """
    lam, U = np.linalg.eigh(S)
    # Guard against near-zero / negative eigenvalues (numerical noise)
    lam = np.maximum(lam, eps)
    return (U * (lam ** -0.5)[np.newaxis, :]) @ U.conj().T


# ─────────────────────────────────────────────────────────────────────────────
# Per-k Löwdin operators
# ─────────────────────────────────────────────────────────────────────────────

def lowdin_operators_k(H_source, k, gauge: str):
    """
    Compute Löwdin-orthogonalised H_ortho, V_x, V_y at a single k-point.

    Returns
    -------
    H_orth : (W,W) complex ndarray
    Vx     : (W,W) complex ndarray
    Vy     : (W,W) complex ndarray
    """
    # Dense matrices at this k
    Hk  = H_source.Hk (k, gauge=gauge).toarray().astype(complex)
    Sk  = H_source.Sk (k, gauge=gauge).toarray().astype(complex)
    dHk = H_source.dHk(k, gauge=gauge)     # list: [x, y, ...]
    dSk = H_source.dSk(k, gauge=gauge)

    dHk_x = dHk[0].toarray().astype(complex)
    dHk_y = dHk[1].toarray().astype(complex)
    dSk_x = dSk[0].toarray().astype(complex)
    dSk_y = dSk[1].toarray().astype(complex)

    # S^{-1/2}  and  S^{-1}
    Sinvh = sqrt_inv(Sk)           # S^{-1/2}
    Sinv  = Sinvh @ Sinvh          # S^{-1}  = (S^{-1/2})²

    # Löwdin Hamiltonian
    H_orth = Sinvh @ Hk @ Sinvh

    # Velocity:
    #   V^α = S^{-1/2} ∂H S^{-1/2}
    #         - ½ H_orth ∂S S^{-1}
    #         - ½ ∂S S^{-1} H_orth
    def velocity(dH, dS):
        return (Sinvh @ dH @ Sinvh
                - 0.5 * H_orth @ dS @ Sinv
                - 0.5 * dS @ Sinv @ H_orth)

    Vx = velocity(dHk_x, dSk_x)
    Vy = velocity(dHk_y, dSk_y)

    return H_orth, Vx, Vy


# ─────────────────────────────────────────────────────────────────────────────
# Build full block-diagonal supercell operators
# ─────────────────────────────────────────────────────────────────────────────

def build_operators(H_source, kx: int, ky: int, kz: int, gauge: str):
    """
    Loop over the kx×ky×kz grid, compute Löwdin operators at each k,
    and return block-diagonal sparse matrices.

    Returns
    -------
    H_blk  : scipy sparse CSR  (Nk*W, Nk*W)
    Vx_blk : scipy sparse CSR
    Vy_blk : scipy sparse CSR  (or None if compute_vy=False)
    """
    K1, K2, K3 = np.mgrid[0:kx, 0:ky, 0:kz]
    Ks = np.column_stack([
        K1.ravel() / kx,
        K2.ravel() / ky,
        K3.ravel() / kz,
    ])
    Nk = len(Ks)

    print(f"  k-grid: {kx}×{ky}×{kz} = {Nk} k-points")

    H_list  = []
    Vx_list = []
    Vy_list = []

    t0 = time.time()
    for ik, k in enumerate(Ks):
        H_orth, Vx, Vy = lowdin_operators_k(H_source, k, gauge)
        H_list .append(csr_matrix(H_orth))
        Vx_list.append(csr_matrix(Vx))
        Vy_list.append(csr_matrix(Vy))

        # Progress on same line
        elapsed = time.time() - t0
        rate    = (ik + 1) / elapsed if elapsed > 0 else 0
        eta     = (Nk - ik - 1) / rate if rate > 0 else 0
        print(f"  k-point {ik+1}/{Nk}  |  {rate:.1f} k/s  |  ETA {eta:.0f}s   ",
              end='\r', flush=True)

    print(f"\n  Done in {time.time()-t0:.1f}s")

    print("  Assembling block-diagonal matrices ...", end=' ', flush=True)
    H_blk  = block_diag(H_list,  format='csr')
    Vx_blk = block_diag(Vx_list, format='csr')
    Vy_blk = block_diag(Vy_list, format='csr') 
    print("done")

    return H_blk, Vx_blk, Vy_blk


# ─────────────────────────────────────────────────────────────────────────────
# Entry point
# ─────────────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Löwdin-orthogonalise a SIESTA Hamiltonian and write "
                    "H, Vx, Vy as lsquant .CSR files."
    )
    parser.add_argument("fdf_file",        help="Path to SIESTA .fdf file")
    parser.add_argument("kx",  type=int,   help="k-points along a1")
    parser.add_argument("ky",  type=int,   help="k-points along a2")
    parser.add_argument("--kz",    type=int, default=1,
                        help="k-points along a3 (default: 1)")
    parser.add_argument("--prefix", default=None,
                        help="Output file prefix (default: fdf stem)")
    parser.add_argument("--gauge", default="r",
                        help="sisl gauge string (default: 'r')")
    args = parser.parse_args()

    fdf_path = Path(args.fdf_file)
    if not fdf_path.exists():
        sys.exit(f"Error: file not found: {fdf_path}")

    prefix = args.prefix if args.prefix else fdf_path.stem

    # ── Read Hamiltonian ──────────────────────────────────────────────────────
    print(f"Reading {fdf_path} ...")
    try:
        import sisl
        fdf   = sisl.get_sile(fdf_path)
        H_src = fdf.read_hamiltonian()
    except Exception as e:
        sys.exit(f"Error reading SIESTA file: {e}")

    print(f"  orbitals per cell (num_wann): {H_src.no}")
    print(f"  spin class: {H_src.spin}")

    # ── Build operators ───────────────────────────────────────────────────────
    print(f"\nBuilding Löwdin operators on {args.kx}×{args.ky}×{args.kz} grid ...")
    H_blk, Vx_blk, Vy_blk = build_operators(
        H_src, args.kx, args.ky, args.kz,
        gauge=args.gauge
    )

    # ── Write outputs ─────────────────────────────────────────────────────────
    print("\nWriting CSR files ...")
    write_lsquant_csr(H_blk,  f"{prefix}.HAM.CSR")
    write_lsquant_csr(Vx_blk, f"{prefix}.VX.CSR")
    write_lsquant_csr(Vy_blk, f"{prefix}.VY.CSR")

    print("\nAll done.")


if __name__ == "__main__":
    main()
