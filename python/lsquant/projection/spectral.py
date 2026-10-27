"""S-orthogonal spectral selection without changing coordinates or overlap.

Symbolic source: metric_forms.spectral_projection_forms/proof. All inputs
and outputs called operators are physical (covariant) matrices, not S^-1 O.
The diagonalisation is preprocessing, once per k; the moment engine is unchanged.
"""
from __future__ import annotations

from dataclasses import dataclass
from collections.abc import Mapping
import numpy as np
from scipy.linalg import eigh


def adjoint(matrix):
    """Conjugate transpose on the final two axes."""
    return np.asarray(matrix).swapaxes(-1, -2).conj()


def resolve_projection(config=None):
    """Validate the projection run section before loading any model.

    Windows are absolute stored energies (eV, E_F=0 for SIESTA), inclusive.
    Radius is in Angstrom; None retains every FFT/WS hopping. Orthogonal
    reduction requires constant window rank and an isolated smooth manifold.
    """
    if config is None:
        config = {}
    if not isinstance(config, Mapping):
        raise ValueError("projection must be an object")
    unknown = set(config)-{"kind", "window", "radius"}
    if unknown:
        raise ValueError(f"unknown projection fields: {sorted(unknown)}")
    kind = config.get("kind", "none")
    if kind not in ("none", "hard", "shift", "orthogonal"):
        raise ValueError("projection.kind must be none, hard, shift or orthogonal")
    window = config.get("window")
    if window is not None:
        window = np.asarray(window, float)
        if window.shape != (2,) or not np.all(np.isfinite(window)) or window[0] >= window[1]:
            raise ValueError("projection.window must contain two finite increasing energies")
        window = window.tolist()
    if kind != "none" and window is None:
        raise ValueError("a spectral projection requires projection.window")
    radius = config.get("radius")
    if radius is not None:
        if isinstance(radius, bool) or not np.isfinite(radius) or radius < 0:
            raise ValueError("projection.radius must be finite and nonnegative")
        radius = float(radius)
    if kind == "none" and (window is not None or radius is not None):
        raise ValueError("projection none does not accept window or radius")
    return dict(kind=kind, window=window, radius=radius)


def physical_blocks(matrix, *, name="matrix"):
    """Check finite Hermitian (N_k,W,W) physical blocks, fold round-off only.

    The scale is the largest operator block, so Gamma zeros are harmless.
    Negative/indefinite overlaps are refused by the Hermitian eigensolver.
    """
    data = np.asarray(matrix, dtype=complex)
    if data.ndim != 3 or data.shape[-1] != data.shape[-2] or not data.size:
        raise ValueError(f"{name} must have nonempty shape (N_k,W,W)")
    if not np.all(np.isfinite(data)):
        raise ValueError(f"{name} contains non-finite entries")
    scale = max(float(np.max(np.linalg.norm(data, axis=(-2, -1)))), np.finfo(float).tiny)
    if np.max(np.linalg.norm(data-adjoint(data), axis=(-2, -1))) > 1e-10*scale:
        raise ValueError(f"{name} is not Hermitian")
    return np.ascontiguousarray((data+adjoint(data))/2)


@dataclass
class SpectralProjection:
    """Full-dimensional projected pencil and its original eigenframe.

    eigenvalues remain in ORIGINAL band order; hard and shift introduce
    degeneracies. Do not sort energies without sorting eigenvectors with them.
    ``kept`` is authoritative, including bands exactly on the window edges.
    """
    hamiltonian: np.ndarray
    overlap: np.ndarray
    eigenvalues: np.ndarray
    eigenvectors: np.ndarray
    original_eigenvalues: np.ndarray
    kept: np.ndarray
    projector: np.ndarray
    kind: str
    window: tuple[float, float]

    def project_operator(self, operator):
        """O_P=P^+ O P=S P (S^-1 O) P, for physical velocity/spin/current.

        Project a supplied spin current itself. Rebuilding {V_P,Sigma_P}
        inserts extra projectors and generally defines a different operator.
        """
        data = physical_blocks(operator, name="physical observable")
        if data.shape != self.overlap.shape:
            raise ValueError("operator and projection have different block dimensions")
        projected = adjoint(self.projector) @ data @ self.projector
        return np.ascontiguousarray((projected+adjoint(projected))/2)

    def receipt(self):
        """Finite, JSON-ready diagnostics; bounds exclude engine padding."""
        P, S, C, E = self.projector, self.overlap, self.eigenvectors, self.eigenvalues
        residual = self.hamiltonian @ C-(S @ C)*E[:, None, :]
        return dict(kind=self.kind, window=list(self.window),
                    kept_bands_per_k=self.kept.sum(axis=1).tolist(),
                    original_bounds=[float(self.original_eigenvalues.min()), float(self.original_eigenvalues.max())],
                    spectral_bounds=[float(E.min()), float(E.max())],
                    matrix_dimension=int(S.shape[-1]),
                    complement_energy=(0.0 if self.kind == "hard" else "nearest window edge"),
                    projected_dos="identity DOS retains all complement modes; use a projected identity word for selected DOS",
                    residuals=dict(idempotency=float(np.max(np.abs(P @ P-P))),
                                   s_self_adjoint=float(np.max(np.abs(adjoint(P) @ S-S @ P))),
                                   eigenpairs_eV=float(np.max(np.abs(residual))),
                                   s_normalisation=float(np.max(np.abs(adjoint(C) @ S @ C-np.eye(S.shape[-1]))))))


def spectral_projection(hamiltonian, overlap, window, *, kind="shift"):
    """Build P=C_W C_W^+ S and H_P or H_s, O(N_k W^3) once.

    Hard: H_P=S C diag(kept*E) C^+ S. Complement modes land at zero,
    even when zero lies OUTSIDE the requested window. Shift: clip E to
    [lo,hi], equivalently H+S C_O diag(E_shift-E_O) C_O^+ S.
    Both keep S exactly; all retained eigenpairs are exact. Neither alone
    preserves a word after its static operators have been projected.
    """
    cfg = resolve_projection(dict(kind=kind, window=window))
    if kind not in ("hard", "shift"):
        raise ValueError("spectral_projection requires hard or shift")
    H = physical_blocks(hamiltonian, name="Hamiltonian")
    S = physical_blocks(overlap, name="overlap")
    if H.shape != S.shape:
        raise ValueError("H and S must have identical shapes")
    pairs = [eigh(h, s, check_finite=False) for h, s in zip(H, S)]
    E = np.array([pair[0] for pair in pairs])
    C = np.array([pair[1] for pair in pairs])
    lo, hi = cfg["window"]
    kept = (E >= lo) & (E <= hi)
    if not np.any(kept):
        raise ValueError("projection window contains no bands on this mesh")
    P = (C*kept[:, None, :]) @ adjoint(C) @ S
    mapped = np.where(kept, E, 0.) if kind == "hard" else np.clip(E, lo, hi)
    SC = S @ C
    projected_h = (SC*mapped[:, None, :]) @ adjoint(SC)
    return SpectralProjection(np.ascontiguousarray((projected_h+adjoint(projected_h))/2),
                              S, mapped, C, E, kept, P, kind, (lo, hi))


def word_trace(energies, eigenvectors, operators, functions):
    """Exact per-cell Tr[O0 f0(H) ...], physical operators, array functions.

    Each function is an array of shape (N_k,W) evaluated at the supplied
    energies. This generic diagnostic covers any ordered spectral word;
    the caller supplies the same parameters for the reference and projection.
    """
    E, C = np.asarray(energies), np.asarray(eigenvectors)
    if len(operators) != len(functions):
        raise ValueError("supply one spectral array for each static operator")
    product = np.broadcast_to(np.eye(E.shape[-1], dtype=complex), (*E.shape, E.shape[-1])).copy()
    for operator, function in zip(operators, functions):
        values = np.asarray(function)
        if values.shape != E.shape:
            raise ValueError("spectral values must have shape (N_k,W)")
        product = product @ ((adjoint(C) @ operator @ C)*values[:, None, :])
    return complex(np.trace(product, axis1=-2, axis2=-1).mean())


def word_error(reference, projected):
    """Absolute and relative complex word errors without an acceptance verdict."""
    reference, projected = np.asarray(reference), np.asarray(projected)
    difference = float(np.max(np.abs(projected-reference)))
    scale = float(np.max(np.abs(reference)))
    return dict(max_absolute=difference, max_relative=(difference/scale if scale > 0 else None))
