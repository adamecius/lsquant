"""Controlled non-orthogonal basis transformations for tight-binding CSR fixtures.

Provides ``make_local_transform`` to build a sparse, narrow-banded basis change
``R = I + alpha * K`` sharing the off-diagonal sparsity graph of a Hamiltonian ``H``,
and ``rotate_basis`` to map any dictionary of orthogonal CSR operators into the
transformed representation ``S' = R^+ R`` and ``O' = R^+ O R``. Also provides ASCII
``.CSR`` reader and writer routines matching ``lsquant::ReadCsrText`` and
``lsquant::WriteCsrText``.
"""

from __future__ import annotations

from pathlib import Path
import sys
import numpy as np
import scipy.linalg as la
import scipy.sparse as sp


def hermitian_fold(matrix: sp.csr_matrix) -> sp.csr_matrix:
    """Return ``0.5 * (M + M^+)`` in canonical sorted CSR layout.

    Reproduces ``lsquant::SparseMatrix::HermitianFold`` (and the legacy
    CSR reader convention) so that operators read in Python match
    the C++ ``ReadCsrText`` representation down to IEEE-754 rounding.

    Parameters
    ----------
    matrix : scipy.sparse.csr_matrix
        Square complex sparse matrix.

    Returns
    -------
    scipy.sparse.csr_matrix
        Hermitian-folded CSR matrix with sorted column indices.
    """
    csr = matrix.tocsr().astype(np.complex128, copy=False)
    folded = (0.5 * (csr.conj().T + csr)).tocsr()
    folded.sum_duplicates()
    folded.sort_indices()
    return folded


def read_csr_text(path: str | Path, fold_hermitian: bool = True, *, legacy_fold: bool = False) -> sp.csr_matrix:
    """Read an LSQUANT ASCII ``.CSR`` operator file into a ``scipy.sparse.csr_matrix``.

    Parses ``dim nnz``, followed by ``nnz`` pairs of ``(re, im)`` floating-point values,
    ``nnz`` zero-based column indices, and ``dim + 1`` row pointer offsets.

    Parameters
    ----------
    path : str or pathlib.Path
        Path to the ``.CSR`` file.
    fold_hermitian : bool, default True
        When True, applies ``0.5 * (M^+ + M)`` to match ``lsquant::ReadCsrText``.

    legacy_fold : bool, default False
        Permit relative ``||(M-M^+)/2||_F / ||M||_F > 1e-10``. The norm is
        reported even when ``fold_hermitian=False``; that option controls only
        the final fold, not validation.

    Returns
    -------
    scipy.sparse.csr_matrix
        Complex CSR matrix of shape ``(dim, dim)``.
    """
    text = Path(path).read_text(encoding="utf-8").split()
    if len(text) < 2:
        raise ValueError(f"read_csr_text: empty or truncated header in {path}")
    dim = int(text[0])
    nnz = int(text[1])
    if dim < 0 or nnz < 0:
        raise ValueError(f"read_csr_text: negative dimension/count in {path}")
    expected_tokens = 2 + 2 * nnz + nnz + (dim + 1)
    if len(text) < expected_tokens:
        raise ValueError(
            f"read_csr_text: expected {expected_tokens} tokens in {path}, got {len(text)}"
        )

    val_tokens = np.array(text[2 : 2 + 2 * nnz], dtype=np.float64)
    values = val_tokens[0::2] + 1j * val_tokens[1::2]

    col_offset = 2 + 2 * nnz
    columns = np.array(text[col_offset : col_offset + nnz], dtype=np.int64)

    rp_offset = col_offset + nnz
    row_pointers = np.array(text[rp_offset : rp_offset + dim + 1], dtype=np.int64)

    if (row_pointers[0] != 0 or row_pointers[-1] != nnz or
            np.any(np.diff(row_pointers) < 0) or np.any(columns < 0) or np.any(columns >= dim)):
        raise ValueError(f"read_csr_text: invalid CSR indices in {path}")
    mat = sp.csr_matrix((values, columns, row_pointers), shape=(dim, dim), dtype=np.complex128)
    mat.sum_duplicates()
    mat.sort_indices()
    anti_norm = 0.5 * float(sp.linalg.norm(mat - mat.conj().T))
    norm = float(sp.linalg.norm(mat))
    relative = anti_norm / max(norm, np.finfo(float).tiny)
    print(f"read_csr_text: {path} anti_hermitian_norm={anti_norm:.17g} "
          f"relative={relative:.17g} legacy_fold={legacy_fold}", file=sys.stderr)
    if not np.isfinite(norm) or not np.isfinite(relative) or (relative > 1e-10 and not legacy_fold):
        raise ValueError(f"read_csr_text: non-Hermitian operator {path}; "
                         "use legacy_fold only to reproduce legacy inputs")
    if fold_hermitian:
        return hermitian_fold(mat)
    return mat


def write_csr_text(path: str | Path, matrix: sp.csr_matrix) -> None:
    """Write a ``scipy.sparse.csr_matrix`` to an LSQUANT ASCII ``.CSR`` operator file.

    Serialises ``dim nnz``, the ``nnz`` complex pairs ``re im`` at 17-digit precision,
    the ``nnz`` zero-based column indices, and the ``dim + 1`` CSR row pointers.

    Parameters
    ----------
    path : str or pathlib.Path
        Destination ``.CSR`` file path.
    matrix : scipy.sparse.csr_matrix
        Square sparse matrix to serialise.
    """
    csr = matrix.tocsr().astype(np.complex128, copy=False)
    csr.sum_duplicates()
    csr.sort_indices()
    if csr.shape[0] != csr.shape[1]:
        raise ValueError("write_csr_text: matrix must be square")

    dim = int(csr.shape[0])
    nnz = int(csr.nnz)
    out_path = Path(path)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    # Normalise -0.0 to +0.0 for deterministic ASCII output.
    re_vals = np.where(csr.data.real == 0.0, 0.0, csr.data.real)
    im_vals = np.where(csr.data.imag == 0.0, 0.0, csr.data.imag)

    val_strs = [f"{re:.17g} {im:.17g}" for re, im in zip(re_vals, im_vals)]
    col_strs = [str(int(c)) for c in csr.indices]
    rp_strs = [str(int(r)) for r in csr.indptr]

    with out_path.open("w", encoding="utf-8") as fh:
        fh.write(f"{dim} {nnz}\n")
        fh.write(" ".join(val_strs) + (" \n" if val_strs else "\n"))
        fh.write(" ".join(col_strs) + (" \n" if col_strs else "\n"))
        fh.write(" ".join(rp_strs) + " \n")


def make_local_transform(
    H: sp.csr_matrix,
    alpha: float = 0.55,
    complex_entries: bool = True,
    seed: int = 12345,
) -> sp.csr_matrix:
    """Build a local basis transformation ``R = I + alpha * K`` on the sparsity graph of ``H``.

    Off-diagonal entries of ``K`` are drawn deterministically from
    ``numpy.random.default_rng(seed)`` (either real non-symmetric when
    ``complex_entries=False`` or complex non-Hermitian ``K_r + i K_i`` when
    ``complex_entries=True``) and normalised so that ``||K||_2 == 1`` and
    ``||alpha * K||_2 <= 0.55``. For ``N <= 4096``, verifies via dense ``eigvalsh``
    that ``S' = R^+ R`` is positive definite and ``cond(S')`` lies in ``[5.0, 10.0]``,
    raising ``ValueError`` outside that band.

    Parameters
    ----------
    H : scipy.sparse.csr_matrix
        Reference Hamiltonian defining the dimension ``N`` and off-diagonal sparsity graph.
    alpha : float, default 0.55
        Coupling amplitude multiplying ``K`` in ``R = I + alpha * K``.
    complex_entries : bool, default True
        If True, generate complex non-Hermitian off-diagonal entries; if False, generate
        real non-symmetric entries.
    seed : int, default 12345
        Seed for ``numpy.random.default_rng``.

    Returns
    -------
    scipy.sparse.csr_matrix
        Invertible transformation matrix ``R = I + alpha * K`` in CSR format.
    """
    csr = H.tocsr()
    if csr.shape[0] != csr.shape[1]:
        raise ValueError("make_local_transform: H must be square")
    dim = int(csr.shape[0])
    coo = csr.tocoo()
    mask = coo.row != coo.col
    rows = coo.row[mask]
    cols = coo.col[mask]
    nnz_off = int(rows.size)

    rng = np.random.default_rng(seed)
    if complex_entries:
        kr = rng.uniform(-1.0, 1.0, size=nnz_off)
        ki = rng.uniform(-1.0, 1.0, size=nnz_off)
        vals = (kr + 1j * ki).astype(np.complex128)
    else:
        kr = rng.uniform(-1.0, 1.0, size=nnz_off)
        vals = kr.astype(np.complex128)

    K = sp.csr_matrix((vals, (rows, cols)), shape=(dim, dim), dtype=np.complex128)
    K.sum_duplicates()
    K.sort_indices()

    if nnz_off > 0 and alpha != 0.0:
        sigma_max = float(np.linalg.norm(K.toarray(), ord=2))
        if sigma_max > 0.0:
            target_k_norm = min(1.0, 0.55 / abs(float(alpha)))
            K = K * (target_k_norm / sigma_max)

    ident = sp.eye(dim, format="csr", dtype=np.complex128)
    R = (ident + float(alpha) * K).tocsr()
    R.sum_duplicates()
    R.sort_indices()

    if dim <= 4096:
        s_mat = hermitian_fold((R.conj().T @ R).tocsr())
        s_eigs = la.eigvalsh(s_mat.toarray(), driver="evd")
        lmin = float(np.min(s_eigs))
        lmax = float(np.max(s_eigs))
        if lmin <= 0.0:
            raise ValueError(
                f"make_local_transform: S' is not positive definite (lambda_min={lmin})"
            )
        cond_s = lmax / lmin
        if not (5.0 <= cond_s <= 10.0):
            raise ValueError(
                f"make_local_transform: cond(S') = {cond_s:.6f} outside [5.0, 10.0]"
            )

    return R


def rotate_basis(
    ops: dict[str, sp.csr_matrix],
    R: sp.csr_matrix,
) -> dict[str, sp.csr_matrix]:
    """Return ``{'S': R.H @ R, **{k: R.H @ O @ R for k, O in ops.items()}}``.

    Applies the congruent transformation ``O' = R^+ O R`` to every operator in ``ops``
    and forms the overlap metric ``S' = R^+ R``.

    Parameters
    ----------
    ops : dict[str, scipy.sparse.csr_matrix]
        Dictionary mapping operator names (e.g. ``'HAM'``, ``'VX'``, ``'VY'``) to
        orthogonal-basis CSR matrices.
    R : scipy.sparse.csr_matrix
        Basis transformation matrix ``R`` of shape ``(N, N)``.

    Returns
    -------
    dict[str, scipy.sparse.csr_matrix]
        Dictionary containing ``'S'`` and every transformed operator in canonical CSR form.
    """
    r_csr = R.tocsr().astype(np.complex128, copy=False)
    r_dag = r_csr.conj().T.tocsr()

    s_mat = hermitian_fold((r_dag @ r_csr).tocsr())
    result: dict[str, sp.csr_matrix] = {"S": s_mat}
    for key, op in ops.items():
        op_csr = op.tocsr().astype(np.complex128, copy=False)
        rotated = (r_dag @ op_csr @ r_csr).tocsr()
        result[key] = hermitian_fold(rotated)
    return result


def adjoint(matrix: np.ndarray) -> np.ndarray:
    """Return ``M^+`` on the last two axes, preserving any block axes."""
    return np.asarray(matrix).swapaxes(-1, -2).conj()


def checked_hermitian(matrix: np.ndarray, *, name: str, tolerance: float = 1e-10) -> np.ndarray:
    """Check Hermiticity against the largest block norm of this operator.

    Symbolic source: metric_forms.hermiticity. The common operator scale
    accepts absolute round-off in symmetry-zero blocks (for example a Gamma
    velocity), while still refusing anti-Hermitian defects at that scale.
    Folding only removes accepted round-off; non-finite matrices are refused.
    """
    matrix = np.asarray(matrix, dtype=np.complex128)
    if matrix.ndim < 2 or matrix.shape[-2] != matrix.shape[-1]:
        raise ValueError(f"{name}: expected square matrix blocks")
    norm = np.linalg.norm(matrix, axis=(-2, -1))
    error = 0.5 * np.linalg.norm(matrix - adjoint(matrix), axis=(-2, -1))
    relative = error / max(float(np.max(norm)), np.finfo(float).tiny)
    if not np.all(np.isfinite(matrix)) or np.any(relative > tolerance):
        raise ValueError(f"{name}: anti-Hermitian relative norm {np.max(relative):.6g} "
                         f"exceeds {tolerance:g}")
    return np.ascontiguousarray(0.5 * (matrix + adjoint(matrix)))


def exact_block_metric(overlap: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Build ``X=U diag(1/s) U^+`` and ``Y=U diag(1/sqrt(s)) U^+`` once.

    ``S=U diag(s) U^+`` must be positive definite; no eigenvalue clipping or
    iterative solve is allowed. Work is O(N_k W^3), storage O(N_k W^2).
    """
    overlap = checked_hermitian(overlap, name="overlap")
    values, vectors = np.linalg.eigh(overlap)
    if np.any(values <= 0):
        raise ValueError(f"overlap is not positive definite: lambda_min={values.min():.17g}")
    inverse = (vectors * (1.0 / values)[..., None, :]) @ adjoint(vectors)
    inverse_sqrt = (vectors * (1.0 / np.sqrt(values))[..., None, :]) @ adjoint(vectors)
    return np.ascontiguousarray(inverse), np.ascontiguousarray(inverse_sqrt)


def bloch_velocity(H: np.ndarray, S: np.ndarray, dH: np.ndarray, dS: np.ndarray,
                   A: np.ndarray | None = None, mode: str = "berry_connection", *,
                   inverse_overlap: np.ndarray | None = None) -> np.ndarray:
    """Return the Hermitian physical velocity matrix, with hbar=1 (eV Angstrom).

    Symbolic source: metric_forms.velocity_forms. Implements theory.md (4): ``V=dH-(dS X H+H X dS)/2``, ``X=S^-1``.
    For ``covariant`` add ``i(H X A-A^+ X H)``. ``A`` is the Hermitian
    midpoint-centred connection, not the unsymmetrised row-centred dipole:
    ``A_ij(R)=<0i|r-(tau_i+tau_j+R)/2|Rj>``. See siesta.py for the conversion.
    ``bare`` requires cell derivatives expressed in the SAME gauge as H and S;
    the other modes require atom derivatives. The caller supplies the derivative
    pair, so the numerical formula occurs only here. In particular no derivative
    of the Lowdin Hamiltonian is substituted for the velocity.
    """
    if mode not in ("bare", "berry_connection", "covariant"):
        raise ValueError(f"unknown velocity mode {mode!r}")
    if mode == "covariant" and A is None:
        raise ValueError("covariant velocity requires an intra-cell dipole matrix")
    X = exact_block_metric(S)[0] if inverse_overlap is None else inverse_overlap
    velocity = dH - 0.5 * (dS @ X @ H + H @ X @ dS)
    if mode == "covariant":
        A = checked_hermitian(A, name="midpoint connection")
        velocity = velocity + 1j * (H @ X @ A - adjoint(A) @ X @ H)
    return checked_hermitian(velocity, name=f"{mode} velocity")


def applied_spin_operators(size: int, ordering: str = "interleaved") -> dict[str, np.ndarray]:
    """Return Pauli operators ``q_a=I_orb tensor sigma_a`` (eigenvalues +/-1).

    ``interleaved`` orders (orbital up, orbital down); ``blocked`` orders
    (all up, all down). ``size`` includes spin and must be even.
    """
    if size <= 0 or size % 2:
        raise ValueError("spinor dimension must be positive and even")
    pauli = (np.array([[0, 1], [1, 0]]), np.array([[0, -1j], [1j, 0]]),
             np.diag([1, -1]))
    if ordering not in ("interleaved", "blocked"):
        raise ValueError("spin ordering must be 'interleaved' or 'blocked'")
    identity = np.eye(size // 2)
    return {f"S{a}": np.asarray(np.kron(identity, p) if ordering == "interleaved"
                                else np.kron(p, identity), dtype=np.complex128)
            for a, p in zip("XYZ", pauli)}


def spin_operators(S: np.ndarray, ordering: str = "interleaved") -> dict[str, np.ndarray]:
    """Return physical matrices ``Sigma_a=S q_a`` with applied form ``q_a``.

    Symbolic source: metric_forms.spin_forms. S is the full spinor overlap ``S_orb tensor I_2``, not the spatial overlap.
    Refuse overlaps that do not commute with all three Pauli operators: those
    need explicitly transformed spin matrices, not a guessed tensor product.
    """
    applied = applied_spin_operators(S.shape[-1], ordering)
    for name, q in applied.items():
        if not np.allclose(S @ q, q @ S, atol=1e-12, rtol=1e-10):
            raise ValueError(f"overlap is not spin scalar for {name}; provide physical spin matrices")
    return {name: checked_hermitian(S @ q, name=name) for name, q in applied.items()}


def spin_current(V: np.ndarray, sigma: np.ndarray) -> np.ndarray:
    """Return ``J^b_a=(V_a q_b+q_b V_a)/2``, with applied Pauli ``q_b``.

    Symbolic source: metric_forms.spin_current_forms. This is the matrix of the operator anticommutator because the physical
    spin matrix is ``Sigma=S q=q S`` and ``AB`` is ``A S^-1 B``. Pauli units
    are used, without hbar/2. For a general transformed spin use
    ``(V S^-1 Sigma+Sigma S^-1 V)/2`` instead.
    """
    return checked_hermitian(0.5 * (V @ sigma + sigma @ V), name="spin current")
