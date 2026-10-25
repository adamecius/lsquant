"""Basis-independent Bloch sums; k in inverse Angstrom, lattice rows in Angstrom."""
from __future__ import annotations

import numpy as np


def reciprocal_lattice(lattice):
    """Return reciprocal rows ``B=2 pi (L^-1)^T``, so ``L B^T=2 pi I``."""
    return 2 * np.pi * np.linalg.inv(lattice).T


def k_mesh(lattice, supercell):
    """Return ``k=(i/N1,j/N2,l/N3) B``, with i varying fastest."""
    if len(supercell) != 3 or any(int(n) != n or n < 1 for n in supercell):
        raise ValueError("supercell must contain three positive integers")
    n1, n2, n3 = map(int, supercell)
    return np.array([[i/n1, j/n2, l/n3] for l in range(n3)
                     for j in range(n2) for i in range(n1)]) @ reciprocal_lattice(lattice)


def bloch_sum(r_cart, centres, matrices, k_points, *, derivative=None, cell=False,
              memory_mb=256.0):
    """Sum ``exp(+ik.(R+tau_j-tau_i)) M_ij(R)`` and optional Cartesian derivative.

    For ``cell=True`` omit tau from the phase and derivative. A derivative
    multiplies by ``+i (R+tau_j-tau_i)_a``. Factor the phase into a lattice
    sum and orbital phases, with O(N_k N_R W^2) work. Chunk k to bound array
    workspace by memory_mb MiB, excluding caller inputs and the final result.
    """
    points = np.asarray(k_points, dtype=float)
    r_cart, centres = np.asarray(r_cart), np.asarray(centres)
    if not np.isfinite(memory_mb) or memory_mb <= 0:
        raise ValueError("Bloch workspace memory_mb must be positive and finite")
    data = np.asarray(matrices, dtype=np.complex128, order="C")
    nr, width, _ = data.shape
    # Reserve block stacks for conversion and derivative weights, plus
    # phase construction, orbital dressing and contraction output per k.
    fixed_bytes = 3*data.nbytes
    bytes_per_k = 32*nr + 48*width + 32*width*width
    chunk = int((memory_mb*1024**2-fixed_bytes)//bytes_per_k)
    if chunk < 1:
        raise ValueError("Bloch workspace limit cannot hold the block stack and one k point")
    if derivative is not None:
        displacement = r_cart[:, derivative, None, None]
        if not cell:
            displacement = displacement + (centres[None, None, :, derivative]
                                           - centres[None, :, None, derivative])
        data = (1j*displacement)*data
        del displacement
    points_2d = np.atleast_2d(points)
    out = np.empty((len(points_2d), width, width), dtype=np.complex128)
    for start in range(0, len(points_2d), chunk):
        stop = min(start+chunk, len(points_2d))
        k = points_2d[start:stop]
        phase = np.asarray(k @ r_cart.T, dtype=np.complex128)
        phase *= 1j
        np.exp(phase, out=phase)
        blocks = out[start:stop]
        np.einsum("kr,rij->kij", phase, data, out=blocks, optimize=True)
        if not cell:
            # The atom-gauge phase is exactly one for i=j. Multiplying two
            # approximate inverse phases would perturb an identity overlap
            # and prevent the exact orthogonal limit of the disorder adapter.
            diagonal = np.arange(width)
            cell_diagonal = blocks[:, diagonal, diagonal].copy()
            orbital = np.asarray(k @ centres.T, dtype=np.complex128)
            orbital *= 1j
            np.exp(orbital, out=orbital)
            blocks *= orbital[:, None, :]
            blocks *= orbital.conj()[:, :, None]
            blocks[:, diagonal, diagonal] = cell_diagonal
            del orbital
        del phase
    return out[0] if points.ndim == 1 else out
