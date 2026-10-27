"""Isolated-window orthogonal downfolder using the origin's polar gauge.

The full origin shell/frozen-window/SCDM searches are also callable through
run_origin_downfolder. This array tool has a narrower, explicit contract:
constant isolated window rank, SCDM trial-column choice, polar gauge, FFT/WS
interpolation, and independent off-grid band validation. No gauge is chosen
from individual eigenvector phases, which are not deterministic at degeneracy.
"""
from __future__ import annotations

from dataclasses import dataclass
import numpy as np
from scipy.linalg import qr

from .origin import origin_module, run_origin_downfolder
from .spectral import spectral_projection, adjoint

__all__ = ["OrthogonalDownfold", "downfold_orthogonal", "validate_bands", "run_origin_downfolder"]


@dataclass
class OrthogonalDownfold:
    """Reduced H, full-to-reduced S-isometry Q and SCDM trial indices."""
    hamiltonian: np.ndarray
    isometry: np.ndarray
    trial_columns: np.ndarray
    minimum_singular_value: float
    original_eigenvalues: np.ndarray
    kept: np.ndarray

    def project_operator(self, operator):
        """O_W=Q^+ O Q, Q^+ S Q=I, for every physical static operator."""
        result = adjoint(self.isometry) @ operator @ self.isometry
        return np.ascontiguousarray((result+adjoint(result))/2)


def downfold_orthogonal(hamiltonian, overlap, window):
    """Reduce a constant isolated manifold with the origin polar helper.

    SCDM-k selects columns by QRCP of mean_k Pt Pt^+, Pt=S^1/2 C_W,
    exactly the isolated-group recipe of downfold_scdm.scdm_columns.
    A=C_W^+ S^1/2[:,cols]; the unchanged origin isolated_U(A) supplies
    U=u v^+; Q=C_W U. Refuse rank loss rather than fill arbitrary vectors.
    """
    projection = spectral_projection(hamiltonian, overlap, window, kind="shift")
    counts = projection.kept.sum(axis=1)
    if np.any(counts != counts[0]) or counts[0] == 0:
        raise ValueError("orthogonal reduction requires constant nonzero window rank")
    rank = int(counts[0])
    C = np.stack([c[:, mask] for c, mask in zip(projection.eigenvectors, projection.kept)])
    E = np.stack([e[mask] for e, mask in zip(projection.original_eigenvalues, projection.kept)])
    s, v = np.linalg.eigh(projection.overlap)
    root = (v*np.sqrt(s)[:, None, :]) @ adjoint(v)
    Pt = root @ C
    rho = np.einsum("kib,kjb->ij", Pt, Pt.conj())/len(Pt)
    _, _, pivots = qr(rho, mode="economic", pivoting=True)
    columns = pivots[:rank]
    A = Pt[:, columns, :].conj().transpose(0, 2, 1)
    U, minimum = origin_module().isolated_U(A)
    smin = float(np.min(minimum))
    if smin < 1e-8:
        raise ValueError(f"orthogonal polar trials lose rank (sigma_min={smin:g}); use the full origin search")
    Q = C @ U
    H = adjoint(U) @ (E[:, :, None]*U)
    return OrthogonalDownfold(np.ascontiguousarray((H+adjoint(H))/2), Q, columns,
                               smin, projection.original_eigenvalues, projection.kept)


def validate_bands(local_hamiltonian, reference, window, lattice, supercell, *, count=40, seed=1):
    """Off-grid band error, independent of training mesh, no fitted alignment.

    reference(cartesian_k) returns generalised eigenvalues. A changed window
    rank is refused: an isolated downfolded manifold must stay isolated away
    from the mesh too. Both eigenvalue lists are sorted with multiplicity.
    """
    rng = np.random.default_rng(seed)
    reduced = rng.random((count, 3))
    reduced[:, np.asarray(supercell) == 1] = 0
    points = reduced @ (2*np.pi*np.linalg.inv(lattice).T)
    maximum = 0.
    for k in points:
        all_e = np.asarray(reference(k))
        selected = np.sort(all_e[(all_e >= window[0]) & (all_e <= window[1])])
        actual = np.linalg.eigvalsh(local_hamiltonian.blocks(k))
        if len(selected) != len(actual):
            raise ValueError("orthogonal window rank changes at an independent off-grid point")
        maximum = max(maximum, float(np.max(np.abs(actual-selected))))
    return dict(max_band_error_eV=maximum, off_grid_points=count, seed=seed)
