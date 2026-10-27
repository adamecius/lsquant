"""Wannier90 model reader and k-space / real-space supercell operator builders.

Reads Wannier90 ``<seed>_hr.dat`` (with ``1 / ndegen`` and optional
``<seed>_wsvec.dat`` folded once at load time), ``<seed>.uc`` (the ``3 x 3``
lattice matrix in Angstrom), and ``<seed>.xyz`` (the ``W`` orbital positions in
Angstrom), and constructs:

- ``H(k)`` blocks in the atom gauge:
  ``H(k)_{a,b} = sum_R exp(+i k . (R + tau_b - tau_a)) H_{a,b}(R)``;
- velocity ``v_d(k)`` blocks in the ``berry_connection`` atom gauge matching
  ``wannier2sparse``:
  ``v_d(k)_{a,b} = -i sum_R (R + tau_b - tau_a)_d exp(+i k . (R + tau_b - tau_a)) H_{a,b}(R)``;
- the supercell k-mesh ``k = (m_1 / N_1) b_1 + (m_2 / N_2) b_2 + (m_3 / N_3) b_3``
  and atom-gauge Bloch phase table ``exp(-i k . tau_a)``;
- real-space supercell CSR matrices ``HAM`` and ``VX``/``VY``/``VZ`` in
  ``n_1``-fastest cell order ``iR = n_1 + N_1 * n_2 + N_1 * N_2 * n_3`` matching
  ``wannier2sparse`` and :class:`~lsquant._core.BlochTransform`.
"""

from __future__ import annotations

from dataclasses import dataclass
import math
import os
from typing import Sequence
import numpy as np
import numpy.typing as npt

from .. import _core

_DIRECTION_MAP: dict[str | int, int] = {
    0: 0,
    1: 1,
    2: 2,
    "x": 0,
    "y": 1,
    "z": 2,
    "X": 0,
    "Y": 1,
    "Z": 2,
    "VX": 0,
    "VY": 1,
    "VZ": 2,
    "vx": 0,
    "vy": 1,
    "vz": 2,
}


def _parse_wsvec(path: str) -> dict[tuple[int, int, int, int, int], list[tuple[int, int, int]]]:
    """Parse a Wannier90 ``_wsvec.dat`` file into a lookup table."""
    with open(path, encoding="utf-8") as fh:
        lines = [ln.strip() for ln in fh if ln.strip()]
    idx = 1  # skip comment line
    lut: dict[tuple[int, int, int, int, int], list[tuple[int, int, int]]] = {}
    while idx < len(lines):
        parts = lines[idx].split()
        idx += 1
        if len(parts) < 6:
            break
        r1, r2, r3, ib1, ib2, ndeg = (int(x) for x in parts[:6])
        shifts: list[tuple[int, int, int]] = []
        for _ in range(ndeg):
            t_parts = [int(x) for x in lines[idx].split()[:3]]
            idx += 1
            shifts.append((t_parts[0], t_parts[1], t_parts[2]))
        lut[(r1, r2, r3, ib1 - 1, ib2 - 1)] = shifts
    return lut


def read_hr(
    path: str, wsvec_path: str | None = None
) -> tuple[int, npt.NDArray[np.int64], npt.NDArray[np.complex128]]:
    """Read a Wannier90 ``_hr.dat`` file (sparse or dense) and fold ``1 / ndegen`` into ``H(R)``.

    Parses ``R_1 R_2 R_3 i j re im`` records by content and assigns each unique
    lattice vector ``R = (R_1, R_2, R_3)`` its degeneracy weight ``ndegen[r]``
    by order of first appearance of ``R`` (matching Wannier90 and
    ``hr_exactdiag.py``), supporting both sparse ``_hr.dat`` files (such as
    ``gen_models.py`` that list only non-zero hoppings) and dense ``nrpts * W^2``
    files as a special case.

    Parameters
    ----------
    path : str
        Path to ``<seed>_hr.dat``.
    wsvec_path : str or None, optional
        Optional path to ``<seed>_wsvec.dat``. If ``None`` and a sibling
        ``<seed>_wsvec.dat`` exists next to ``path``, it is loaded
        automatically.

    Returns
    -------
    tuple[int, numpy.ndarray, numpy.ndarray]
        Tuple ``(num_wann, r_int, hoppings)`` where ``num_wann`` is ``W``,
        ``r_int`` has shape ``(n_R, 3)`` with ``int64`` lattice vectors, and
        ``hoppings`` has shape ``(n_R, W, W)`` with ``complex128`` matrices
        after dividing by ``ndegen`` (and splitting across ``wsvec`` shifts when
        present).
    """
    with open(path, encoding="utf-8") as fh:
        raw_lines = fh.read().splitlines()
    if len(raw_lines) < 3:
        raise ValueError(f"read_hr: truncated _hr.dat file '{path}'")

    # Line 0 is the Wannier90 comment/header line; collect remaining non-empty lines.
    lines = [ln.strip() for ln in raw_lines[1:] if ln.strip()]
    if len(lines) < 2:
        raise ValueError(f"read_hr: truncated _hr.dat header in '{path}'")

    num_wann = int(lines[0].split()[0])
    nrpts = int(lines[1].split()[0])
    if num_wann <= 0 or nrpts <= 0:
        raise ValueError(
            f"read_hr: invalid num_wann={num_wann} or nrpts={nrpts} in '{path}'"
        )

    ndegen: list[int] = []
    idx = 2
    while len(ndegen) < nrpts and idx < len(lines):
        ndegen.extend(int(x) for x in lines[idx].split())
        idx += 1
    if len(ndegen) < nrpts:
        raise ValueError(
            f"read_hr: expected {nrpts} ndegen integers in '{path}', found {len(ndegen)}"
        )

    if wsvec_path is None and path.endswith("_hr.dat"):
        candidate = path[: -len("_hr.dat")] + "_wsvec.dat"
        if os.path.isfile(candidate):
            wsvec_path = candidate

    ws_lut = _parse_wsvec(wsvec_path) if (wsvec_path and os.path.isfile(wsvec_path)) else None

    r_order: dict[tuple[int, int, int], int] = {}
    accum: dict[tuple[int, int, int], npt.NDArray[np.complex128]] = {}
    for line in lines[idx:]:
        parts = line.split()
        if len(parts) < 7:
            continue
        r1 = int(parts[0])
        r2 = int(parts[1])
        r3 = int(parts[2])
        a = int(parts[3]) - 1
        b = int(parts[4]) - 1
        if a < 0 or a >= num_wann or b < 0 or b >= num_wann:
            raise ValueError(
                f"read_hr: orbital index ({a + 1}, {b + 1}) out of range 1..{num_wann} in '{path}'"
            )
        raw_r = (r1, r2, r3)
        if raw_r not in r_order:
            r_order[raw_r] = len(r_order)
        r_idx = r_order[raw_r]
        deg = (
            float(ndegen[r_idx])
            if (r_idx < len(ndegen) and ndegen[r_idx] > 0)
            else 1.0
        )
        val = complex(float(parts[5]), float(parts[6])) / deg
        if abs(val) <= np.finfo(np.float64).eps:
            continue
        if ws_lut is not None:
            shifts = ws_lut.get((r1, r2, r3, a, b))
            if shifts:
                sub_val = val / float(len(shifts))
                for t1, t2, t3 in shifts:
                    key = (r1 + t1, r2 + t2, r3 + t3)
                    mat = accum.get(key)
                    if mat is None:
                        mat = np.zeros((num_wann, num_wann), dtype=np.complex128)
                        accum[key] = mat
                    mat[a, b] += sub_val
                continue
        mat = accum.get(raw_r)
        if mat is None:
            mat = np.zeros((num_wann, num_wann), dtype=np.complex128)
            accum[raw_r] = mat
        mat[a, b] += val

    keys = sorted(accum.keys())
    if not keys:
        keys = [(0, 0, 0)]
        accum[(0, 0, 0)] = np.zeros((num_wann, num_wann), dtype=np.complex128)

    r_int = np.array(keys, dtype=np.int64)
    hoppings = np.stack([accum[k] for k in keys], axis=0)
    return num_wann, r_int, hoppings


def read_uc(path: str) -> npt.NDArray[np.float64]:
    """Read the ``3 x 3`` lattice matrix ``A`` (rows = ``a_1, a_2, a_3`` in Angstrom) from ``<seed>.uc``."""
    rows: list[list[float]] = []
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            stripped = line.split("#", 1)[0].strip()
            if not stripped:
                continue
            parts = stripped.split()
            if len(parts) < 3:
                continue
            rows.append([float(parts[0]), float(parts[1]), float(parts[2])])
            if len(rows) == 3:
                break
    if len(rows) != 3:
        raise ValueError(f"read_uc: expected 3 lattice vectors in '{path}'")
    return np.array(rows, dtype=np.float64)


def _try_parse_xyz_coord(line: str) -> list[float] | None:
    """Attempt to parse an orbital line ``<label> <x> <y> <z>``."""
    parts = line.split()
    if len(parts) < 4:
        return None
    try:
        return [float(parts[1]), float(parts[2]), float(parts[3])]
    except ValueError:
        return None


def read_xyz(path: str) -> npt.NDArray[np.float64]:
    """Read the ``W`` orbital positions ``tau_alpha`` in Angstrom from ``<seed>.xyz``.

    Accepts both ``wannier2sparse``'s single-header layout (count line followed
    directly by ``W`` ``<label> <x> <y> <z>`` lines with no comment line) and
    the two-header XYZ layout (count line, comment line, then ``W`` orbital
    lines).
    """
    with open(path, encoding="utf-8") as fh:
        raw_lines = [ln.strip() for ln in fh.read().splitlines()]
    non_empty = [ln for ln in raw_lines if ln]
    if len(non_empty) < 2:
        raise ValueError(f"read_xyz: truncated .xyz file '{path}'")

    num_orbs = int(non_empty[0].split()[0])
    if num_orbs <= 0:
        raise ValueError(f"read_xyz: invalid orbital count {num_orbs} in '{path}'")

    rest = non_empty[1:]
    if len(rest) == num_orbs:
        orb_lines = rest
    elif len(rest) >= num_orbs + 1 and (
        _try_parse_xyz_coord(rest[0]) is None or len(rest) == num_orbs + 1
    ):
        orb_lines = rest[1 : 1 + num_orbs]
    else:
        orb_lines = rest[:num_orbs]

    coords: list[list[float]] = []
    for line in orb_lines:
        parsed = _try_parse_xyz_coord(line)
        if parsed is None:
            raise ValueError(f"read_xyz: malformed orbital line '{line}' in '{path}'")
        coords.append(parsed)
    if len(coords) != num_orbs:
        raise ValueError(
            f"read_xyz: expected {num_orbs} orbital lines in '{path}', found {len(coords)}"
        )
    return np.array(coords, dtype=np.float64)


def reciprocal_lattice(uc: npt.NDArray[np.float64]) -> npt.NDArray[np.float64]:
    """Return the ``3 x 3`` reciprocal lattice matrix ``B = 2 * pi * inv(A)^T`` (rows = ``b_1, b_2, b_3``)."""
    return 2.0 * math.pi * np.linalg.inv(np.asarray(uc, dtype=np.float64)).T


def k_mesh(
    uc: npt.NDArray[np.float64], supercell: Sequence[int]
) -> npt.NDArray[np.float64]:
    """Build the ``N_k = N_1 * N_2 * N_3`` k-points in ``m_1``-fastest order ``(m_1, m_2, m_3)``.

    Each point is ``k = (m_1 / N_1) b_1 + (m_2 / N_2) b_2 + (m_3 / N_3) b_3`` in
    inverse Angstrom, with flat index ``ik = m_1 + N_1 * m_2 + N_1 * N_2 * m_3``
    matching ``wannier2sparse``'s cell numbering.
    """
    if len(supercell) != 3:
        raise ValueError("k_mesh: supercell must have 3 integer dimensions (N1, N2, N3)")
    n1, n2, n3 = (int(supercell[0]), int(supercell[1]), int(supercell[2]))
    if n1 <= 0 or n2 <= 0 or n3 <= 0:
        raise ValueError("k_mesh: supercell dimensions must be positive")

    b_mat = reciprocal_lattice(uc)
    f1 = np.arange(n1, dtype=np.float64) / float(n1)
    f2 = np.arange(n2, dtype=np.float64) / float(n2)
    f3 = np.arange(n3, dtype=np.float64) / float(n3)
    g3, g2, g1 = np.meshgrid(f3, f2, f1, indexing="ij")
    frac = np.stack([g1.ravel(), g2.ravel(), g3.ravel()], axis=-1)
    return np.ascontiguousarray(frac @ b_mat, dtype=np.float64)


def bloch_phases(
    k_points: npt.NDArray[np.float64], tau: npt.NDArray[np.float64]
) -> npt.NDArray[np.complex128]:
    """Compute the atom-gauge phase table ``exp(-i k . tau_alpha)`` of shape ``(N_k, W)``."""
    k_arr = np.atleast_2d(np.asarray(k_points, dtype=np.float64))
    tau_arr = np.asarray(tau, dtype=np.float64)
    phase_arg = k_arr @ tau_arr.T  # shape (N_k, W)
    return np.ascontiguousarray(np.exp(-1j * phase_arg), dtype=np.complex128)


def hk_blocks(
    r_cart: npt.NDArray[np.float64],
    tau: npt.NDArray[np.float64],
    hoppings: npt.NDArray[np.complex128],
    k_points: npt.NDArray[np.float64],
) -> npt.NDArray[np.complex128]:
    """Evaluate ``H(k)_{a,b} = sum_R exp(+i k . (R + tau_b - tau_a)) H_{a,b}(R)`` across ``k_points``.

    Parameters
    ----------
    r_cart : numpy.ndarray
        Cartesian lattice vectors ``R`` of shape ``(n_R, 3)`` in Angstrom.
    tau : numpy.ndarray
        Orbital positions ``tau_alpha`` of shape ``(W, 3)`` in Angstrom.
    hoppings : numpy.ndarray
        Real-space hopping blocks ``H(R)`` of shape ``(n_R, W, W)``.
    k_points : numpy.ndarray
        Wavevectors of shape ``(N_k, 3)`` or ``(3,)`` in inverse Angstrom.

    Returns
    -------
    numpy.ndarray
        Contiguous ``complex128`` array of shape ``(N_k, W, W)`` (or ``(W, W)``
        if a single 1D ``k`` vector was passed).
    """
    k_in = np.asarray(k_points, dtype=np.float64)
    single = k_in.ndim == 1
    k_arr = np.atleast_2d(k_in)

    # Lattice Fourier sum: H_lat(k) = sum_R exp(+i k . R) H(R), shape (N_k, W, W)
    lat_phases = np.exp(1j * (k_arr @ r_cart.T))  # (N_k, n_R)
    h_lat = np.einsum("kr,rab->kab", lat_phases, hoppings, optimize=True)

    # Atom-gauge dressing: exp(+i k . (tau_b - tau_a)) = exp(-i k . tau_a) * exp(+i k . tau_b)
    orb_phase = np.exp(1j * (k_arr @ tau.T))  # (N_k, W)
    blocks = np.conj(orb_phase)[:, :, None] * h_lat * orb_phase[:, None, :]
    blocks = np.ascontiguousarray(blocks, dtype=np.complex128)
    return blocks[0] if single else blocks


def dhk_blocks(
    r_cart: npt.NDArray[np.float64],
    tau: npt.NDArray[np.float64],
    hoppings: npt.NDArray[np.complex128],
    k_points: npt.NDArray[np.float64],
    direction: int | str,
) -> npt.NDArray[np.complex128]:
    """Evaluate the atom-gauge velocity ``v_d(k)_{a,b}`` across ``k_points``.

    Uses the ``wannier2sparse`` ``berry_connection`` sign convention:
    ``v_d(k)_{a,b} = -i sum_R (R + tau_b - tau_a)_d exp(+i k . (R + tau_b - tau_a)) H_{a,b}(R)``
    so that ``B VX_csr B^+ == v_x`` holds to machine precision.
    """
    dir_idx = _DIRECTION_MAP.get(direction, -1)
    if dir_idx not in (0, 1, 2):
        raise ValueError(f"dhk_blocks: unsupported direction '{direction}'")

    k_in = np.asarray(k_points, dtype=np.float64)
    single = k_in.ndim == 1
    k_arr = np.atleast_2d(k_in)

    # Displacement (R + tau_b - tau_a)_d for each (r, a, b), shape (n_R, W, W)
    disp = (
        r_cart[:, dir_idx, None, None]
        + tau[None, None, :, dir_idx]
        - tau[None, :, None, dir_idx]
    )
    weighted_hop = (-1j * disp) * hoppings  # (n_R, W, W)

    lat_phases = np.exp(1j * (k_arr @ r_cart.T))  # (N_k, n_R)
    v_lat = np.einsum("kr,rab->kab", lat_phases, weighted_hop, optimize=True)

    orb_phase = np.exp(1j * (k_arr @ tau.T))  # (N_k, W)
    blocks = np.conj(orb_phase)[:, :, None] * v_lat * orb_phase[:, None, :]
    blocks = np.ascontiguousarray(blocks, dtype=np.complex128)
    return blocks[0] if single else blocks


def _build_supercell_sparse_matrix(
    num_wann: int,
    r_int: npt.NDArray[np.int64],
    hop_blocks: npt.NDArray[np.complex128],
    supercell: Sequence[int],
) -> _core.SparseMatrix:
    """Assemble a periodic supercell :class:`~lsquant._core.SparseMatrix` in ``n_1``-fastest cell order.

    Uses the ``wannier2sparse`` cell index ``iR = n_1 + N_1 * n_2 + N_1 * N_2 * n_3``
    and enforces the minimum-image condition ``N_d >= 2 * max(|R_d|) + 1`` along
    every axis ``d`` that carries non-zero ``R_d``.
    """
    if len(supercell) != 3:
        raise ValueError("supercell must have 3 integer dimensions (N1, N2, N3)")
    n1, n2, n3 = (int(supercell[0]), int(supercell[1]), int(supercell[2]))
    if n1 <= 0 or n2 <= 0 or n3 <= 0:
        raise ValueError("supercell dimensions must be positive")

    dims = (n1, n2, n3)
    for d in range(3):
        max_r = int(np.max(np.abs(r_int[:, d]))) if r_int.shape[0] > 0 else 0
        if max_r > 0 and dims[d] < 2 * max_r + 1:
            raise ValueError(
                f"supercell N={dims[d]} along axis {d} is too small for hopping range "
                f"{max_r} (minimum image requires N >= {2 * max_r + 1})"
            )

    num_cells = n1 * n2 * n3
    total_dim = num_cells * num_wann

    # Extract non-zero primitive hopping entries (dr1, dr2, dr3, a, b, val)
    nz_r, nz_a, nz_b = np.nonzero(np.abs(hop_blocks) > np.finfo(np.float64).eps)
    if nz_r.size == 0:
        rp = np.zeros(total_dim + 1, dtype=np.int64)
        cols = np.zeros(0, dtype=np.int64)
        vals = np.zeros(0, dtype=np.complex128)
        return _core.SparseMatrix.from_csr(total_dim, rp, cols, vals)

    dr1 = r_int[nz_r, 0]
    dr2 = r_int[nz_r, 1]
    dr3 = r_int[nz_r, 2]
    vals_nz = hop_blocks[nz_r, nz_a, nz_b]

    # Grid of source cells (c1, c2, c3) in wannier2sparse order: iR = c1 + n1 * c2 + n1 * n2 * c3
    c3, c2, c1 = np.meshgrid(
        np.arange(n3, dtype=np.int64),
        np.arange(n2, dtype=np.int64),
        np.arange(n1, dtype=np.int64),
        indexing="ij",
    )
    c1_flat = c1.ravel()[:, None]  # (N_k, 1)
    c2_flat = c2.ravel()[:, None]
    c3_flat = c3.ravel()[:, None]
    ir_src = c1_flat + n1 * (c2_flat + n2 * c3_flat)  # (N_k, 1)

    t1 = (c1_flat + dr1[None, :]) % n1
    t2 = (c2_flat + dr2[None, :]) % n2
    t3 = (c3_flat + dr3[None, :]) % n3
    ir_dst = t1 + n1 * (t2 + n2 * t3)  # (N_k, E)

    rows = (ir_src * num_wann + nz_a[None, :]).ravel()
    cols = (ir_dst * num_wann + nz_b[None, :]).ravel()
    data = np.broadcast_to(vals_nz[None, :], ir_dst.shape).ravel()

    # Sort by (row, col) and sum any duplicate entries
    order = np.lexsort((cols, rows))
    rows_s = rows[order]
    cols_s = cols[order]
    data_s = np.ascontiguousarray(data[order], dtype=np.complex128)

    if rows_s.size > 1:
        unique_mask = np.empty(rows_s.size, dtype=bool)
        unique_mask[0] = True
        unique_mask[1:] = (rows_s[1:] != rows_s[:-1]) | (cols_s[1:] != cols_s[:-1])
        if not np.all(unique_mask):
            idx_starts = np.flatnonzero(unique_mask)
            rows_s = rows_s[idx_starts]
            cols_s = cols_s[idx_starts]
            data_s = np.add.reduceat(data_s, idx_starts)

    counts = np.bincount(rows_s, minlength=total_dim)
    rp = np.empty(total_dim + 1, dtype=np.int64)
    rp[0] = 0
    np.cumsum(counts, out=rp[1:])
    cols_out = np.ascontiguousarray(cols_s, dtype=np.int64)
    vals_out = np.ascontiguousarray(data_s, dtype=np.complex128)
    return _core.SparseMatrix.from_csr(total_dim, rp, cols_out, vals_out)


def supercell_csr(
    num_wann: int,
    r_int: npt.NDArray[np.int64],
    hoppings: npt.NDArray[np.complex128],
    supercell: Sequence[int],
) -> _core.SparseMatrix:
    """Assemble the real-space supercell Hamiltonian CSR matrix."""
    return _build_supercell_sparse_matrix(num_wann, r_int, hoppings, supercell)


def supercell_velocity_csr(
    num_wann: int,
    r_int: npt.NDArray[np.int64],
    r_cart: npt.NDArray[np.float64],
    tau: npt.NDArray[np.float64],
    hoppings: npt.NDArray[np.complex128],
    supercell: Sequence[int],
    direction: int | str,
) -> _core.SparseMatrix:
    """Assemble the real-space supercell velocity CSR matrix in the atom gauge."""
    dir_idx = _DIRECTION_MAP.get(direction, -1)
    if dir_idx not in (0, 1, 2):
        raise ValueError(f"supercell_velocity_csr: unsupported direction '{direction}'")

    disp = (
        r_cart[:, dir_idx, None, None]
        + tau[None, None, :, dir_idx]
        - tau[None, :, None, dir_idx]
    )
    vel_blocks = np.ascontiguousarray((-1j * disp) * hoppings, dtype=np.complex128)
    return _build_supercell_sparse_matrix(num_wann, r_int, vel_blocks, supercell)


@dataclass(frozen=True)
class WannierModel:
    """In-memory Wannier90 tight-binding model feeding Real, Bloch, and Eigen bases.

    Attributes
    ----------
    num_wann : int
        Number of Wannier orbitals ``W`` per unit cell.
    uc : numpy.ndarray
        ``3 x 3`` real-space lattice matrix (rows = ``a_1, a_2, a_3`` in Angstrom).
    tau : numpy.ndarray
        ``W x 3`` orbital coordinates in Angstrom.
    r_int : numpy.ndarray
        ``n_R x 3`` integer lattice displacement vectors.
    r_cart : numpy.ndarray
        ``n_R x 3`` Cartesian lattice displacement vectors ``r_int @ uc`` in Angstrom.
    hoppings : numpy.ndarray
        ``n_R x W x W`` complex hopping matrices ``H(R)`` with ``ndegen`` folded.
    seed_name : str
        Model seed prefix (e.g. ``"graphene"`` or ``"haldane"``).
    """

    num_wann: int
    uc: npt.NDArray[np.float64]
    tau: npt.NDArray[np.float64]
    r_int: npt.NDArray[np.int64]
    r_cart: npt.NDArray[np.float64]
    hoppings: npt.NDArray[np.complex128]
    seed_name: str = "model"
    model_dir: str | None = None

    def Hk(self, k_points: npt.NDArray[np.float64]) -> npt.NDArray[np.complex128]:
        """Evaluate ``H(k)`` blocks in the atom gauge for ``k_points``."""
        return hk_blocks(self.r_cart, self.tau, self.hoppings, k_points)

    def dHk(
        self, k_points: npt.NDArray[np.float64], direction: int | str
    ) -> npt.NDArray[np.complex128]:
        """Evaluate atom-gauge velocity ``v_d(k)`` blocks for ``k_points``."""
        return dhk_blocks(self.r_cart, self.tau, self.hoppings, k_points, direction)

    def velocity_blocks(self, k_points, direction, mode="berry_connection"):
        """Build Wannier velocity with the Wannier90 velocity sign convention.

        ``V=-dH`` in the bare/cell or centre/atom convention. Covariant mode
        implements ``V=-dH_cell-i[H,A_full]``, with folded, weighted A loaded
        before the commutator. This preserves the Wannier90 sign convention;
        SIESTA uses the explicit +dH convention. A Hall
        word has two velocities, so their simultaneous reversal cancels.
        """
        from .bloch import bloch_sum
        from .nonorthogonal import bloch_velocity

        axis = _DIRECTION_MAP.get(direction, -1)
        if axis not in (0, 1, 2):
            raise ValueError(f"unknown velocity direction {direction!r}")
        H = self.Hk(k_points)
        S = np.broadcast_to(np.eye(self.num_wann), H.shape)
        # Cell derivatives are dressed into the atom gauge of H.
        dH = bloch_sum(self.r_cart, np.zeros_like(self.tau), self.hoppings,
                      k_points, derivative=axis)
        phases = np.exp(1j * (np.atleast_2d(k_points) @ self.tau.T))
        if np.asarray(k_points).ndim == 1:
            phases = phases[0]
        dH = phases.conj()[..., :, None] * dH * phases[..., None, :]
        if mode == "berry_connection":
            dH = -self.dHk(k_points, direction)
        A = None
        if mode == "covariant":
            if self.model_dir is None:
                raise ValueError("covariant Wannier velocity requires a model directory with _r.dat")
            if not hasattr(self, "_position_cache"):
                seed = os.path.join(self.model_dir, self.seed_name)
                self._position_cache = read_position_matrix(seed + "_r.dat", seed + "_hr.dat")
            r, positions = self._position_cache
            # Here A is the FULL position matrix, so use cell derivatives;
            # atom dressing changes representation but adds no centre term.
            A = bloch_sum(r @ self.uc, self.tau, positions[axis], k_points)
        return -bloch_velocity(H, S, dH, np.zeros_like(dH), A, mode,
                               inverse_overlap=S)

    def k_mesh(self, supercell: Sequence[int]) -> npt.NDArray[np.float64]:
        """Return the ``(N_k, 3)`` k-mesh for ``supercell``."""
        return k_mesh(self.uc, supercell)

    def bloch_phases(
        self, supercell_or_kpoints: Sequence[int] | npt.NDArray[np.float64]
    ) -> npt.NDArray[np.complex128]:
        """Return the ``(N_k, W)`` atom-gauge phase table ``exp(-i k . tau_alpha)``."""
        arr = np.asarray(supercell_or_kpoints)
        if arr.ndim == 1 and arr.size == 3 and np.issubdtype(arr.dtype, np.integer):
            k_pts = self.k_mesh([int(arr[0]), int(arr[1]), int(arr[2])])
        else:
            k_pts = np.asarray(supercell_or_kpoints, dtype=np.float64)
        return bloch_phases(k_pts, self.tau)

    def hamiltonian_csr(self, supercell: Sequence[int]) -> _core.SparseMatrix:
        """Build the real-space supercell Hamiltonian CSR matrix."""
        return supercell_csr(self.num_wann, self.r_int, self.hoppings, supercell)

    def velocity_csr(
        self, supercell: Sequence[int], direction: int | str
    ) -> _core.SparseMatrix:
        """Build the real-space supercell velocity CSR matrix in direction ``direction``."""
        return supercell_velocity_csr(
            self.num_wann,
            self.r_int,
            self.r_cart,
            self.tau,
            self.hoppings,
            supercell,
            direction,
        )


def load_wannier_model(model_dir: str, seed: str | None = None) -> WannierModel:
    """Load ``<seed>_hr.dat``, ``<seed>.uc``, and ``<seed>.xyz`` from ``model_dir``.

    Parameters
    ----------
    model_dir : str
        Directory containing the Wannier90 input triplet (or path to ``<seed>_hr.dat``).
    seed : str or None, optional
        Model prefix (e.g. ``"graphene"``). If ``None``, inferred from the unique
        ``*_hr.dat`` file in ``model_dir``.

    Returns
    -------
    WannierModel
        Loaded tight-binding model with ``ndegen`` folded.
    """
    if os.path.isfile(model_dir) and model_dir.endswith("_hr.dat"):
        hr_path = model_dir
        base_dir = os.path.dirname(model_dir) or "."
        inferred_seed = os.path.basename(model_dir)[: -len("_hr.dat")]
        seed_name = seed if seed is not None else inferred_seed
    else:
        base_dir = model_dir
        if seed is None:
            hr_files = [
                f for f in os.listdir(base_dir) if f.endswith("_hr.dat")
            ]
            if len(hr_files) != 1:
                raise ValueError(
                    f"load_wannier_model: specify seed explicitly in '{base_dir}' "
                    f"(found {hr_files})"
                )
            seed_name = hr_files[0][: -len("_hr.dat")]
        else:
            seed_name = seed
        hr_path = os.path.join(base_dir, f"{seed_name}_hr.dat")

    uc_path = os.path.join(base_dir, f"{seed_name}.uc")
    xyz_path = os.path.join(base_dir, f"{seed_name}.xyz")

    num_wann, r_int, hoppings = read_hr(hr_path)
    uc = read_uc(uc_path)
    tau = read_xyz(xyz_path)
    if tau.shape[0] != num_wann:
        raise ValueError(
            f"load_wannier_model: orbital count mismatch between '{hr_path}' ({num_wann}) "
            f"and '{xyz_path}' ({tau.shape[0]})"
        )
    r_cart = np.ascontiguousarray(r_int.astype(np.float64) @ uc, dtype=np.float64)
    return WannierModel(
        num_wann=num_wann,
        uc=uc,
        tau=tau,
        r_int=r_int,
        r_cart=r_cart,
        hoppings=hoppings,
        seed_name=seed_name,
        model_dir=os.path.abspath(base_dir),
    )


def read_position_matrix(path: str, hr_path: str, wsvec_path: str | None = None
                         ) -> tuple[np.ndarray, np.ndarray]:
    """Read ``_r.dat`` and return Hermitian, ndegen-normalised ``A_a(R)``.

    Implements notes section 11: ``A_a(R) <- [A_a(R)+A_a(-R)^+]/(2 ndegen(R))``
    BEFORE any Fourier sum or commutator. Degeneracies come from the matching
    ``_hr.dat`` (``_r.dat`` has no degeneracy block). Opposite R must have equal
    weights. Optional wsvec image splitting uses the SAME weights as H, after
    this fold; the split result is checked for real-space Hermiticity.
    Returns ``(r_int, A)`` with A shaped ``(3,N_R,W,W)`` in Angstrom.
    """
    from .nonorthogonal import adjoint

    with open(hr_path, encoding="utf-8") as fh:
        lines = [ln.strip() for ln in fh.read().splitlines()[1:] if ln.strip()]
    width, count = int(lines[0]), int(lines[1])
    weights, index = [], 2
    while len(weights) < count:
        weights.extend(map(int, lines[index].split()))
        index += 1
    order = list(dict.fromkeys(tuple(map(int, ln.split()[:3])) for ln in lines[index:]))
    if len(weights) != count or len(order) != count or any(w <= 0 for w in weights):
        raise ValueError("position matrix: invalid matching hr degeneracy table")
    degeneracy = dict(zip(order, weights))
    with open(path, encoding="utf-8") as fh:
        records = [ln.strip() for ln in fh.read().splitlines()[1:] if ln.strip()]
    if int(records[0]) != width or int(records[1]) != count:
        raise ValueError("position matrix dimensions differ from hr")
    raw = {r: np.zeros((3, width, width), complex) for r in order}
    seen = set()
    for line in records[2:]:
        fields = line.split()
        if len(fields) != 11:
            raise ValueError(f"malformed position matrix record: {line}")
        r = tuple(map(int, fields[:3]))
        i, j = int(fields[3])-1, int(fields[4])-1
        if r not in raw or not (0 <= i < width and 0 <= j < width) or (r, i, j) in seen:
            raise ValueError("position matrix has unknown R/index or duplicate record")
        seen.add((r, i, j))
        values = np.array(fields[5:], float).reshape(3, 2)
        raw[r][:, i, j] = values[:, 0] + 1j * values[:, 1]
    if len(seen) != count * width * width:
        raise ValueError("truncated position matrix")
    folded = {}
    for r in order:
        opposite = tuple(-x for x in r)
        if opposite not in raw or degeneracy[opposite] != degeneracy[r]:
            raise ValueError("position matrix requires paired R with equal ndegen")
        folded[r] = 0.5 * (raw[r] + adjoint(raw[opposite])) / degeneracy[r]
    if wsvec_path is None:
        candidate = str(hr_path).removesuffix("_hr.dat") + "_wsvec.dat"
        wsvec_path = candidate if os.path.isfile(candidate) else None
    if wsvec_path:
        lut = _parse_wsvec(wsvec_path)
        split = {}
        for r, matrix in folded.items():
            for i in range(width):
                for j in range(width):
                    shifts = lut.get((*r, i, j), [(0, 0, 0)])
                    for shift in shifts:
                        key = tuple(r[a] + shift[a] for a in range(3))
                        split.setdefault(key, np.zeros_like(matrix))[:, i, j] += matrix[:, i, j] / len(shifts)
        folded = split
    for r, matrix in folded.items():
        opposite = tuple(-x for x in r)
        if opposite not in folded or not np.allclose(matrix, adjoint(folded[opposite]), atol=1e-12, rtol=1e-10):
            raise ValueError("position matrix image splitting violates Hermiticity")
    keys = sorted(folded)
    return np.array(keys, dtype=int), np.stack([folded[r] for r in keys], axis=1)
