"""NumPy/SciPy generalised-eigenproblem spin Hall and spin-Chern oracle.

Pauli sigma_z, e=hbar=1 internally; output sigma in e^2/h. The physical
spin-current response need not be quantised when spin is not conserved.
The independent Prodan sector invariant remains integer while both the
energy gap and the projected-spin gap stay open. Never round a conductivity
onto an integer or equate it with the sector invariant by construction.
"""
from __future__ import annotations

import numpy as np
from scipy.special import erfc

UNITS = "e^2/h, sigma_z Pauli"


def adjoint(matrix):
    """Return ``M^+`` on the last two axes."""
    return np.asarray(matrix).swapaxes(-1, -2).conj()


def generalised_eigh(H, S, *, inverse_sqrt=None):
    """Solve ``H C=S C diag(E), C^+ S C=I`` using NumPy Hermitian eigh.

    ``Y=S^-1/2; Y H Y U=U E; C=Y U``. An already cached Y can be supplied
    to avoid a second metric factorisation. Nonpositive overlap is refused.
    """
    if inverse_sqrt is None:
        s, u = np.linalg.eigh(S)
        if np.any(s <= 0):
            raise ValueError("oracle overlap is not positive definite")
        inverse_sqrt = (u * (1/np.sqrt(s))[..., None, :]) @ adjoint(u)
    H_L = inverse_sqrt @ H @ inverse_sqrt
    E, U = np.linalg.eigh(0.5*(H_L+adjoint(H_L)))
    return E, inverse_sqrt @ U


def spin_berry_curvature(energies, eigenvectors, J, V, *, degeneracy_tolerance=1e-6,
                         hamiltonian=None, overlap=None, return_diagnostics=False):
    """Return ``Omega_n=-2 Im sum_(m!=n) J_nm V_mn/(E_n-E_m)^2``.

    Matrix elements are ``C^+ J C, C^+ V C`` with C^+SC=I. Terms within
    degenerate multiplets are excluded: individual band curvatures there
    depend on gauge, but the occupied multiplet sum is invariant if E_F
    lies in a gap. No artificial eta is inserted in interband denominators.
    Adjacent sorted energies separated by at most the tolerance form one
    subspace, including transitive chains. The eV tolerance is at least the
    requested value, ten times the energy residual ||H c-E S c||/||S c||
    (when H and S are supplied), and ten times eps*W*max(bandwidth, |E|).
    Optional diagnostics count unordered intra-subspace pairs, once per k.
    Work O(N_k W^3), scratch O(W^2), result (N_k,W).
    """
    if not np.isfinite(degeneracy_tolerance) or degeneracy_tolerance <= 0:
        raise ValueError("oracle degeneracy tolerance must be positive and finite (eV)")
    if (hamiltonian is None) != (overlap is None):
        raise ValueError("both H and S are required for the eigensolver residual")
    curvature = np.empty_like(energies, dtype=float)
    tolerances, residuals, bandwidths = [], [], []
    merged_pairs = 0
    for k, (E, C) in enumerate(zip(energies, eigenvectors)):
        bandwidth = float(np.ptp(E))
        roundoff = np.finfo(float).eps * len(E) * max(bandwidth, float(np.max(np.abs(E))))
        residual = 0.0
        if hamiltonian is not None:
            sc = overlap[k] @ C
            residual = float(np.max(np.linalg.norm(hamiltonian[k] @ C-sc*E, axis=0) /
                                    np.maximum(np.linalg.norm(sc, axis=0), np.finfo(float).tiny)))
            if not np.isfinite(residual):
                raise ValueError("non-finite oracle eigensolver residual")
        tolerance = max(degeneracy_tolerance, 10*residual, 10*roundoff)
        order = np.argsort(E)
        subspaces = np.empty(len(E), dtype=int)
        subspaces[order] = np.cumsum(np.r_[True, np.diff(E[order]) > tolerance])
        same_subspace = subspaces[:, None] == subspaces[None, :]
        merged_pairs += int(np.count_nonzero(np.triu(same_subspace, k=1)))
        tolerances.append(tolerance)
        residuals.append(residual)
        bandwidths.append(bandwidth)
        j, v = adjoint(C) @ J[k] @ C, adjoint(C) @ V[k] @ C
        gap = E[:, None]-E[None, :]
        inverse_gap = np.zeros_like(gap)
        np.divide(1.0, gap*gap, out=inverse_gap, where=~same_subspace)
        curvature[k] = -2*np.imag(np.sum(j*v.T*inverse_gap, axis=1))
    if return_diagnostics:
        return curvature, dict(requested_tolerance_eV=float(degeneracy_tolerance),
                               tolerance_eV=float(max(tolerances)),
                               tolerance_per_k_eV=[float(t) for t in tolerances],
                               max_eigensolver_residual_eV=max(residuals),
                               max_bandwidth_eV=max(bandwidths), merged_pairs=merged_pairs,
                               pair_count_convention="unordered intra-subspace pairs summed over k")
    return curvature


def eigenphase_trace(blocks, *, states=16, seed=12345):
    """Estimate ``Tr F = mean_r sum_k r_k^+ F(k) r_k`` in the eigenbasis.

    Each r_kn has unit modulus; diagonal F is traced exactly by every sample.
    Only cross terms WITHIN a k block fluctuate. This NumPy helper is also
    useful for oracle diagnostics; KPM uses the engine's seeded generator.
    """
    if states <= 0:
        raise ValueError("states must be positive")
    rng = np.random.default_rng(seed)
    values = []
    for _ in range(states):
        r = np.exp(2j*np.pi*rng.random(np.asarray(blocks).shape[:-1]))
        values.append(np.einsum("ki,kij,kj->", r.conj(), blocks, r, optimize=True))
    return complex(np.mean(values))


def sigma_xy(energies, curvature, fermi_energies, cell_area, *, broadening=0.0):
    """Return ``sigma(E_F)=2 pi/(N_k A_cell) sum_kn f(E_kn) Omega_kn``.

    A Gaussian-broadened step has ``f=erfc((E-E_F)/(sqrt(2)*eta))/2``;
    eta is its standard deviation in eV. eta=0 is the strict occupied sum.
    This is a diagonal eigenbasis trace, therefore ANY eigenphase sample
    gives the same occupied sum. Sigma is in e^2/h with Pauli spin.
    """
    if cell_area <= 0 or broadening < 0:
        raise ValueError("positive area and nonnegative Gaussian width required")
    E, omega = np.asarray(energies), np.asarray(curvature)
    ef = np.atleast_1d(fermi_energies)
    # One array covers all Fermi energies, k points and bands. Reuse it for
    # erfc and the occupied curvature to avoid another array of that size.
    if broadening == 0:
        occupied = (E[None, :, :] < ef[:, None, None]).astype(float)
    else:
        occupied = E[None, :, :]-ef[:, None, None]
        occupied /= np.sqrt(2)*broadening
        erfc(occupied, out=occupied)
        occupied *= 0.5
    occupied *= omega[None, :, :]
    return 2*np.pi*np.sum(occupied, axis=(1, 2))/(len(E)*cell_area)


def _fhs(frames, *, link_tolerance=1e-12):
    """Compute ``C=-sum plaquette arg(Ux Uy_x Ux_y^* Uy^*)/(2 pi)``.

    ``Ua=det(F_k^+ F_(k+a))/abs(det(...))``. The minus sign matches the
    ``-2 Im`` Berry-curvature convention above. Frames live in a common,
    PERIODIC orthonormal representation and have shape (Ny,Nx,W,sector).
    A singular link signals an insufficient mesh or a closed sector gap.
    """
    links = []
    for axis in (1, 0):
        overlap = adjoint(frames) @ np.roll(frames, -1, axis=axis)
        sign, logabs = np.linalg.slogdet(overlap)
        if np.any(np.abs(sign) == 0) or np.any(logabs < np.log(link_tolerance)):
            raise ValueError("singular FHS link: refine the mesh or inspect the spin gap")
        links.append(sign)
    ux, uy = links
    loop = ux * np.roll(uy, -1, axis=1) * np.roll(ux, -1, axis=0).conj() * uy.conj()
    return float(-np.angle(loop).sum()/(2*np.pi))


def spin_chern(energies, eigenvectors, overlap, spin_z, fermi_energy, mesh, *,
               k_points=None, centres=None, spin_gap_tolerance=1e-7, sqrt_overlap=None):
    """Prodan's P sigma_z P sectors, independently integrated by FHS.

    Within occupied C_o, diagonalise ``C_o^+ Sigma_z C_o`` and separate
    positive/negative eigenvalues. Use Lowdin frames ``S^1/2 C_o W_+/-``.
    A cached ``sqrt_overlap=S Y`` avoids another overlap eigensolve.
    Remove the atom-gauge phase with ``exp(+ik.tau)`` BEFORE periodic
    boundary links. Return the difference ``C_spin=C_+-C_-`` (2 for KM),
    the Prodan spin Chern number ``C_s=(C_+-C_-)/2`` (1 for KM), and
    ``Z2 = C_s mod 2`` when C_s is integral. These are the positive/negative
    Prodan projected-spin sectors: the returned spin_chern field is the full
    sector difference (2 for Kane-Mele), and prodan_half_difference is C_s.
    """
    nx, ny, nz = map(int, mesh)
    if nz != 1 or nx < 2 or ny < 2 or nx*ny != len(energies):
        raise ValueError("spin Chern requires a complete n1-fastest 2D mesh")
    occupied = np.sum(energies < fermi_energy, axis=1)
    if np.any(occupied != occupied[0]) or not (0 < occupied[0] < energies.shape[1]):
        raise ValueError("Fermi energy is not in a common occupied-band gap")
    count = int(occupied[0])
    if np.min(energies[:, count])-np.max(energies[:, count-1]) <= 0:
        raise ValueError("no global insulating energy gap")
    if sqrt_overlap is None:
        s, u = np.linalg.eigh(overlap)
        if np.any(s <= 0):
            raise ValueError("spin Chern overlap must be positive definite")
        sqrt_s = (u*np.sqrt(s)[..., None, :]) @ adjoint(u)
    else:
        # The basis already built Y=S^-1/2 once: S^1/2=S Y. Reuse it.
        sqrt_s = np.asarray(sqrt_overlap)
        if sqrt_s.shape != np.asarray(overlap).shape:
            raise ValueError("cached overlap square root has the wrong shape")
    C = eigenvectors[..., :count]
    projected_spin = adjoint(C) @ spin_z @ C
    eigen_spin, w = np.linalg.eigh(0.5*(projected_spin+adjoint(projected_spin)))
    min_gap = float(np.min(np.abs(eigen_spin)))
    if min_gap <= spin_gap_tolerance:
        raise ValueError(f"projected spin gap closed: {min_gap:g}")
    sizes = np.sum(eigen_spin > 0, axis=1)
    if np.any(sizes != sizes[0]) or sizes[0] in (0, count):
        raise ValueError("projected spin sectors do not have fixed positive ranks")
    frames = sqrt_s @ C @ w
    if (k_points is None) != (centres is None):
        raise ValueError("both k_points and centres are needed for atom-gauge sewing")
    if k_points is not None:
        frames = np.exp(1j*(k_points @ centres.T))[..., :, None] * frames
    neg_count = count-int(sizes[0])
    minus = _fhs(frames[..., :neg_count].reshape(ny, nx, frames.shape[-2], neg_count))
    plus = _fhs(frames[..., neg_count:].reshape(ny, nx, frames.shape[-2], count-neg_count))
    difference = plus-minus
    half = 0.5*difference
    return dict(chern_plus=plus, chern_minus=minus, spin_chern=difference,
                prodan_half_difference=half, z2=int(round(half)) % 2 if abs(half-round(half)) < 1e-6 else None,
                projected_spin_gap=min_gap, occupied_bands=count,
                energy_gap=float(np.min(energies[:, count])-np.max(energies[:, count-1])),
                convention="C_plus-C_minus; Z2 uses the half-difference")
