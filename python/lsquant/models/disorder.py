"""Nonorthogonal Anderson potential and direct periodic SIESTA CSR assembly.

No sisl or Wannier dependency. Translations are kept unwrapped until the
velocity derivative has been formed; indices alone are wrapped on the torus.
"""
from __future__ import annotations

import numpy as np
import scipy.sparse as sparse


def checked_csr(matrix, *, name, tolerance=1e-10):
    """Refuse non-Hermitian input before removing accepted round-off."""
    matrix = sparse.csr_matrix(matrix, dtype=np.complex128)
    matrix.sum_duplicates()
    matrix.sort_indices()
    norm = float(sparse.linalg.norm(matrix))
    anti = 0.5*float(sparse.linalg.norm(matrix-matrix.conj().T))
    if not np.isfinite(norm) or not np.isfinite(anti) or anti > tolerance*max(norm, np.finfo(float).tiny):
        raise ValueError(f"{name}: anti-Hermitian norm {anti:g}, operator norm {norm:g}")
    result = (0.5*(matrix+matrix.conj().T)).tocsr()
    result.eliminate_zeros()
    result.sort_indices()
    return result


def anderson_potential(overlap, values):
    """metric_forms.anderson_forms: V_ij=(w_i+w_j)*S_ij/2 as CSR.

    The RNG/layout matches make_anderson_disorder: one draw per real
    coefficient (including spin), n1 fastest, then orbital/spinor index.
    A scalar w gives w*S; using diag(w) alone would give the wrong identity.
    """
    S = checked_csr(overlap, name="Anderson overlap")
    w = np.asarray(values)
    if w.shape != (S.shape[0],) or not np.all(np.isfinite(w)) or np.any(w.imag != 0):
        raise ValueError(
            f"Expected finite real Anderson values of shape {(S.shape[0],)}, got shape {w.shape}"
        )
    diagonal = sparse.diags(w.real, format="csr")
    return checked_csr(0.5*(diagonal @ S+S @ diagonal), name="Anderson potential")


def supercell_csr(model, supercell, matrices, *, derivative=None, cell=False):
    """Assemble M[(c,i),(c+R,j)] += M_ij(R), c+R modulo N.

    ``iR=n1+N1*n2+N1*N2*n3`` matches the BlochTransform cell ordering. Derivatives
    multiply by i*(R.lat+tau_j-tau_i)_a BEFORE wrapping, or i*R_a for
    the bare rung. Together with the midpoint dipole this implements the
    covering-lattice position commutator (theory.md eq. 3), not a commutator
    with discontinuous wrapped coordinates. Duplicate bonds are summed.
    Work/storage are O(N_cells*nnz_R); no Fourier inversion is used.
    """
    if len(supercell) != 3 or any(int(n) != n or n < 1 for n in supercell):
        raise ValueError("supercell must contain three positive integers")
    n1, n2, n3 = map(int, supercell)
    cells = np.array([(i, j, k) for k in range(n3) for j in range(n2) for i in range(n1)])
    W = model.num_wann
    shape = (len(model.r_int), W, W)
    data = np.asarray(matrices, complex)
    if data.shape != shape:
        raise ValueError("real-space blocks must have shape (N_R,W,W)")
    rows, columns, values = [], [], []
    for translation, block in zip(model.r_int, data):
        i, j = np.nonzero(block)
        if i.size == 0:
            continue
        entries = block[i, j]
        if derivative is not None:
            distance = float((translation @ model.uc)[derivative])
            displacement = distance if cell else distance+model.tau[j, derivative]-model.tau[i, derivative]
            entries = 1j*displacement*entries
        target = (cells+translation) % np.array([n1, n2, n3])
        target_flat = target[:, 0]+n1*target[:, 1]+n1*n2*target[:, 2]
        rows.append((np.arange(len(cells))[:, None]*W+i).ravel())
        columns.append((target_flat[:, None]*W+j).ravel())
        values.append(np.broadcast_to(entries, (len(cells), len(i))).ravel())
    dim = len(cells)*W
    if not rows:
        return sparse.csr_matrix((dim, dim), dtype=np.complex128)
    matrix = sparse.coo_matrix((np.concatenate(values), (np.concatenate(rows), np.concatenate(columns))),
                               shape=(dim, dim)).tocsr()
    return checked_csr(matrix, name="SIESTA supercell", tolerance=model.input_tolerance)
