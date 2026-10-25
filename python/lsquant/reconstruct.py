"""Rank-1 and rank-2 Chebyshev spectral and conductivity reconstruction routines.

Reconstructs the density of states (:func:`dos_curve`) and the four rank-2
conductivity curves (:func:`kubo_greenwood`, :func:`kubo_bastin`,
:func:`kubo_bastin_sea`, :func:`kubo_bastin_surf`) plus the Chebyshev FFT-node
grid variants (:func:`kubo_fft_grid`) in the exact floating-point operation
order of the legacy C++ reconstruction routines.
"""

from __future__ import annotations

import math
import os
import numpy as np
import numpy.typing as npt

from . import _core
from .kernels import jackson_damping, jackson_damping_2d
from .spectral_functions import legacy_delta_table, legacy_green_tables


def dos_curve(
    moments: _core.ChebyshevMoments,
    broadening_meV: float,
    alpha: float = 0.95,
) -> tuple[npt.NDArray[np.float64], npt.NDArray[np.float64]]:
    """Reconstruct the extensive density of states curve from rank-1 moments.

    Solves the problem of turning a rank-1 Chebyshev moment vector ``mu_m`` and
    a target energy resolution ``broadening_meV`` into the physical energy grid
    ``E`` (in eV) and extensive density of states ``N * rho(E)`` (in 1/eV),
    reproducing the operation order of legacy ``spectralFunctionFromChebmom``
    byte for byte:

    1. ``eta = 2.0 * broadening_meV / 1000 / BandWidth``,
       ``M_eff = min(ceil(pi / eta), M)``.
    2. Truncate moments to ``M_eff``, multiply complex ``mu[m]`` by the Jackson
       kernel ``g(m, M_eff)``, then take the real part ``mu_real[m] = mu[m].real``.
    3. Build ``num_div = 30 * M_eff`` uniform points
       ``x_i = -alpha + i * (2 * alpha) / (num_div - 1)``.
    4. Contract ``rho_i = sum_{m=0}^{M_eff-1} delta_chebF(x_i, m) * mu_real[m]``
       sequentially in ``m`` via :func:`lsquant._core.ContractRank1`.
    5. Map to physical energy ``E_i = (x_i - Shift) / Scale`` and extensive value
       ``value_i = rho_i * (SystemSize * Scale)``.

    Complexity is ``O(M_eff * num_div) = O(30 * M_eff^2)`` time and memory.

    Parameters
    ----------
    moments : lsquant._core.ChebyshevMoments
        Rank-1 Chebyshev moment tensor with populated ``window`` and
        ``system_size``.
    broadening_meV : float
        Target Jackson Gaussian broadening in meV (``broadening_meV > 0``).
    alpha : float, optional
        Dimensionless Chebyshev domain cutoff keeping the grid inside
        ``[-alpha, alpha]`` away from ``x = +-1`` (default ``0.95``).

    Returns
    -------
    tuple[numpy.ndarray, numpy.ndarray]
        Pair ``(energies, values)`` of 1D ``float64`` arrays of length
        ``30 * M_eff``, where ``energies`` is in eV and ``values`` is the
        extensive density of states ``SystemSize * rho(E)`` in 1/eV.

    Warnings
    --------
    Raises ``ValueError`` if ``moments.rank != 1``, ``broadening_meV <= 0``,
    ``!(0 < alpha < 1)``, or the spectral window has non-positive width.

    Notes
    -----
    Pitfall: The returned ``values`` array is the extensive density of states
    ``N * rho(E)`` written by ``spectralFunctionFromChebmom``; divide by
    ``moments.system_size`` when plotting the intensive per-site density of
    states.
    """
    if moments.rank != 1:
        raise ValueError("dos_curve: moments must have rank 1")
    if broadening_meV <= 0.0:
        raise ValueError("dos_curve: broadening_meV must be positive")
    if not (0.0 < alpha < 1.0):
        raise ValueError("dos_curve: alpha must lie strictly between 0 and 1")

    band_width = (moments.window.upper - moments.window.lower) * 1.0
    if band_width <= 0.0:
        raise ValueError("dos_curve: spectral window must have positive width")

    eta = 2.0 * float(broadening_meV) / 1000.0 / band_width
    max_mom = int(math.ceil(math.pi / eta))
    total_mom = int(moments.dims[0])
    if max_mom > total_mom:
        max_mom = total_mom

    g = jackson_damping(max_mom)
    damped = np.empty(max_mom, dtype=np.complex128)
    raw = moments.data
    for m_idx in range(max_mom):
        damped[m_idx] = raw[m_idx] * g[m_idx]

    num_div = 30 * max_mom
    x_grid = [
        -alpha + i * (2.0 * alpha) / (num_div - 1) for i in range(num_div)
    ]

    table = legacy_delta_table(x_grid, max_mom)
    rho = _core.ContractRank1(damped.real, table)

    scale = moments.window.scale
    shift = moments.window.shift
    prefactor = float(moments.system_size) * scale

    energies = np.empty(num_div, dtype=np.float64)
    values = np.empty(num_div, dtype=np.float64)
    for i in range(num_div):
        energies[i] = (x_grid[i] - shift) / scale
        values[i] = float(rho[i]) * prefactor

    return energies, values


def _truncate_and_damp_2d(
    moments: _core.ChebyshevMoments,
    b0_meV: float,
    b1_meV: float,
) -> tuple[npt.NDArray[np.complex128], int, int]:
    """Truncate and apply 2D Jackson damping to ``mu_matrix``."""
    if moments.rank != 2:
        raise ValueError("rank-2 reconstruction requires moments.rank == 2")
    if b0_meV <= 0.0 or b1_meV <= 0.0:
        raise ValueError("broadening_meV must be positive")

    band_width = (moments.window.upper - moments.window.lower) * 1.0
    if band_width <= 0.0:
        raise ValueError("spectral window must have positive width")

    eta0 = 2.0 * float(b0_meV) / 1000.0 / band_width
    eta1 = 2.0 * float(b1_meV) / 1000.0 / band_width

    num_moms0 = int(moments.dims[0])
    num_moms1 = int(moments.dims[1])
    max_mom0 = int(math.ceil(math.pi / eta0))
    max_mom1 = int(math.ceil(math.pi / eta1))
    if max_mom0 > num_moms0 or max_mom0 <= 0:
        max_mom0 = num_moms0
    if max_mom1 > num_moms1 or max_mom1 <= 0:
        max_mom1 = num_moms1

    w2d = jackson_damping_2d(max_mom0, max_mom1)
    raw = moments.data
    damped = np.empty((max_mom0, max_mom1), dtype=np.complex128)
    for m0 in range(max_mom0):
        for m1 in range(max_mom1):
            damped[m0, m1] = raw[m0, m1] * float(w2d[m0, m1])

    return damped, max_mom0, max_mom1


def kubo_greenwood(
    moments: _core.ChebyshevMoments,
    broadening_meV: float,
    alpha: float = 0.95,
) -> tuple[npt.NDArray[np.float64], npt.NDArray[np.float64]]:
    """Reconstruct the Kubo-Greenwood conductivity curve on the uniform grid.

    Solves the problem of evaluating the longitudinal Kubo-Greenwood response
    from a rank-2 moment tensor ``mu_{m0, m1}`` in the exact operation order of
    legacy ``reconstruct_kubo(mu, KUBO_GREENWOOD, GRID_UNIFORM)``:

    1. Truncate and damp ``mu`` with :func:`jackson_damping_2d` to ``(M0, M1)``,
       ``M = max(M0, M1)``.
    2. Build ``num_div = M`` uniform points on ``[-alpha, alpha]``.
    3. Contract ``out_i = sum_{m0, m1} delta(x_i, m0) * imag(greenR(x_i, m1) * mu[m0, m1])``
       via :func:`~lsquant._core.ContractRank2Imag`.
    4. Scale ``kernel_i = -1.0 * out_i * (N * scale * scale)`` and report the
       first ``num_div - 1`` points at ``E_i = (x_i - shift) / scale``.

    Parameters
    ----------
    moments : lsquant._core.ChebyshevMoments
        Rank-2 Chebyshev moment tensor.
    broadening_meV : float
        Jackson Gaussian broadening in meV (``broadening_meV > 0``).
    alpha : float, optional
        Chebyshev domain cutoff ``(0 < alpha < 1)``, default ``0.95``.

    Returns
    -------
    tuple[numpy.ndarray, numpy.ndarray]
        Pair ``(energies, values)`` of 1D ``float64`` arrays of length
        ``max(M0, M1) - 1``.

    Warnings
    --------
    Raises ``ValueError`` if ``moments.rank != 2``, ``broadening_meV <= 0``, or
    ``!(0 < alpha < 1)``.

    Notes
    -----
    Pitfall: Drops the last grid point (``i = num_div - 1``) so the output length
    matches legacy ``accumulate(..., integrate=false)``.
    """
    if not (0.0 < alpha < 1.0):
        raise ValueError("kubo_greenwood: alpha must lie strictly between 0 and 1")

    damped, m0, m1 = _truncate_and_damp_2d(moments, broadening_meV, broadening_meV)
    m_max = m1 if m1 > m0 else m0
    num_div = 1 * m_max
    x_grid = [
        -alpha + i * (2.0 * alpha) / (num_div - 1) for i in range(num_div)
    ]

    delta_tbl = legacy_delta_table(x_grid, m0)
    green_r_tbl, _ = legacy_green_tables(x_grid, m1, m1)
    out = _core.ContractRank2Imag(damped, delta_tbl, green_r_tbl)

    scale = moments.window.scale
    shift = moments.window.shift
    pref_ext = float(moments.system_size) * scale * scale

    n_out = max(0, num_div - 1)
    energies = np.empty(n_out, dtype=np.float64)
    values = np.empty(n_out, dtype=np.float64)
    for i in range(n_out):
        kernel_i = -1.0 * float(out[i]) * pref_ext
        values[i] = kernel_i
        energies[i] = (x_grid[i] - shift) / scale

    return energies, values


def _folded_bastin_contraction(moments, x):
    """Contract on arbitrary real nodes using two real Chebyshev series.

    cos(m theta) exp(-i n theta) = (exp(-i(n-m)theta) +
    exp(-i(n+m)theta))/2. Sum the two frequency diagonals once, also
    weighted by n for the Green-function derivative. Pair +/- frequencies
    into cos(l theta) and sin(l theta)=sin(theta) U_(l-1)(x), then use
    U_j=2(T_j+T_(j-2)+...)-[j even] T_0.

    If F=A+i sin(theta) B and F_n=C+i sin(theta) D, the requested
    imaginary contraction is
    (-x Re(A)+(1-x*x) Re(D)+sin(theta)*(x Im(B)+Im(C))) /
    (pi*(1-x*x)^2). This changes summation order, not nodes or kernels.
    Cost O(M0*M1+(M0+M1)*N), scratch O(M0+M1+N); no M-by-N tables.
    """
    m0, m1 = moments.shape
    offset, highest = m0-1, m0+m1-2
    frequencies = np.zeros((2, offset+highest+1), dtype=np.complex128)
    n = np.arange(m1, dtype=float)
    for m in range(m0):
        row = 0.5*moments[m]
        weighted = n*row
        for start in (offset-m, offset+m):
            frequencies[0, start:start+m1] += row
            frequencies[1, start:start+m1] += weighted
    positive = frequencies[:, offset:].copy()
    negative = np.zeros_like(positive)
    negative[:, 1:offset+1] = frequencies[:, :offset][:, ::-1]
    cosine = positive+negative
    sine = negative-positive

    def u_to_t(coefficients):
        if not len(coefficients):
            return np.zeros(1)
        result = np.empty_like(coefficients)
        for parity in (0, 1):
            result[parity::2] = 2*np.cumsum(coefficients[parity::2][::-1])[::-1]
        if len(result):
            result[0] *= 0.5
        return result

    from numpy.polynomial.chebyshev import chebadd, chebsub, chebmulx, chebval
    b = u_to_t(sine[0, 1:].imag)
    d = u_to_t(sine[1, 1:].real)
    polynomial = chebsub(chebsub(d, chebmulx(chebmulx(d))), chebmulx(cosine[0].real))
    radical = chebadd(chebmulx(b), cosine[1].imag)
    x = np.asarray(x, dtype=float)
    sine_squared = 1-x*x
    return (chebval(x, polynomial)+np.sqrt(sine_squared)*chebval(x, radical)) / sine_squared / sine_squared / math.pi


def kubo_bastin(
    moments: _core.ChebyshevMoments,
    broadening_meV: float,
    alpha: float = 0.95,
    *, backend: str = "tables",
) -> tuple[npt.NDArray[np.float64], npt.NDArray[np.float64]]:
    """Reconstruct the Kubo-Bastin conductivity curve on the uniform grid.

    Solves the problem of evaluating the full Fermi-sea Kubo-Bastin conductivity
    from a rank-2 moment tensor ``mu_{m0, m1}`` in the exact operation order of
    legacy ``reconstruct_kubo(mu, KUBO_BASTIN, GRID_UNIFORM)``:

    1. Truncate and damp ``mu`` with :func:`jackson_damping_2d` to ``(M0, M1)``,
       ``M = max(M0, M1)``.
    2. Build ``num_div = 30 * M`` uniform points on ``[-alpha, alpha]``.
    3. Contract ``out_i = sum_{m0, m1} delta(x_i, m0) * imag(DgreenR(x_i, m1) * mu[m0, m1])``
       via :func:`~lsquant._core.ContractRank2Imag`.
    4. Scale ``kernel_i = -2.0 * out_i * (N * scale * scale)`` and integrate
       cumulatively via the forward rectangle rule
       ``acc += kernel_i * (x_{i+1} - x_i)`` for ``i < num_div - 1``.

    Parameters
    ----------
    moments : lsquant._core.ChebyshevMoments
        Rank-2 Chebyshev moment tensor.
    broadening_meV : float
        Jackson Gaussian broadening in meV (``broadening_meV > 0``).
    alpha : float, optional
        Chebyshev domain cutoff ``(0 < alpha < 1)``, default ``0.95``.
    backend : str, optional
        ``tables`` preserves legacy accumulation; ``folded`` evaluates the
        same finite series and grid with O(M^2+M*N) work and bounded scratch.

    Returns
    -------
    tuple[numpy.ndarray, numpy.ndarray]
        Pair ``(energies, values)`` of 1D ``float64`` arrays of length
        ``30 * max(M0, M1) - 1``.

    Warnings
    --------
    Raises ``ValueError`` if ``moments.rank != 2``, ``broadening_meV <= 0``, or
    ``!(0 < alpha < 1)``.

    Notes
    -----
    Pitfall: Accumulates ``acc += kernel_i * (x_{i+1} - x_i)`` in dimensionless
    ``x`` coordinates while ``pref_ext = N * scale * scale`` carries the second
    factor of ``scale``.
    """
    if not (0.0 < alpha < 1.0):
        raise ValueError("kubo_bastin: alpha must lie strictly between 0 and 1")
    if backend not in ("tables", "folded"):
        raise ValueError("kubo_bastin: backend must be 'tables' or 'folded'")

    damped, m0, m1 = _truncate_and_damp_2d(moments, broadening_meV, broadening_meV)
    m_max = m1 if m1 > m0 else m0
    num_div = 30 * m_max
    x_grid = [
        -alpha + i * (2.0 * alpha) / (num_div - 1) for i in range(num_div)
    ]

    if backend == "folded":
        out = _folded_bastin_contraction(damped, x_grid)
    else:
        delta_tbl = legacy_delta_table(x_grid, m0)
        _, dgreen_r_tbl = legacy_green_tables(x_grid, m0, m1)
        out = _core.ContractRank2Imag(damped, delta_tbl, dgreen_r_tbl)

    scale = moments.window.scale
    shift = moments.window.shift
    pref_ext = float(moments.system_size) * scale * scale
    kernel = [-2.0 * float(out[i]) * pref_ext for i in range(num_div)]

    n_out = max(0, num_div - 1)
    energies = np.empty(n_out, dtype=np.float64)
    values = np.empty(n_out, dtype=np.float64)
    acc = 0.0
    for i in range(n_out):
        acc += kernel[i] * (x_grid[i + 1] - x_grid[i])
        values[i] = acc
        energies[i] = (x_grid[i] - shift) / scale

    return energies, values


def kubo_bastin_sea(
    moments: _core.ChebyshevMoments,
    broadening_meV: float,
    alpha: float = 0.95,
) -> tuple[npt.NDArray[np.float64], npt.NDArray[np.float64]]:
    """Reconstruct the Fermi-sea contribution of Kubo-Bastin conductivity.

    Solves the problem of evaluating the Green-difference Fermi-sea part of
    Kubo-Bastin in the exact operation order of legacy
    ``reconstruct_kubo_bastin_sea``:

    1. Truncate and damp ``mu`` with :func:`jackson_damping_2d` to ``(M0, M1)``,
       ``M = max(M0, M1)``.
    2. Build ``num_div = 30 * M`` uniform points on ``[-alpha, alpha]``.
    3. Build ``A = Gr - Ga`` with ``Ga = conj(Gr)`` and ``B = DGr + DGa`` with
       ``DGr = real(DgreenR)`` and ``DGa = complex(DGr, 0.0)`` (matching C++
       ``std::conj(double)``), then contract
       ``out_i = sum_{m0, m1} real((A[m0, i] * B[m1, i]) * mu[m0, m1])`` via
       :func:`~lsquant._core.ContractRank2Real`.
    4. Set ``kernel_i = -out_i``, multiply by ``N * scale * scale / 2 / pi``,
       and integrate cumulatively with the trapezoid rule
       ``acc += 0.5 * (kernel_i + kernel_{i+1}) * (x_{i+1} - x_i)``.

    Parameters
    ----------
    moments : lsquant._core.ChebyshevMoments
        Rank-2 Chebyshev moment tensor.
    broadening_meV : float
        Jackson Gaussian broadening in meV (``broadening_meV > 0``).
    alpha : float, optional
        Chebyshev domain cutoff ``(0 < alpha < 1)``, default ``0.95``.

    Returns
    -------
    tuple[numpy.ndarray, numpy.ndarray]
        Pair ``(energies, values)`` of 1D ``float64`` arrays of length
        ``30 * max(M0, M1) - 1``.

    Warnings
    --------
    Raises ``ValueError`` if ``moments.rank != 2``, ``broadening_meV <= 0``, or
    ``!(0 < alpha < 1)``.

    Notes
    -----
    Pitfall: Preserve the two sequential divisions ``/ 2.0 / math.pi`` in the
    prefactor to match legacy ``/ 2 / M_PI`` rounding.
    """
    if not (0.0 < alpha < 1.0):
        raise ValueError("kubo_bastin_sea: alpha must lie strictly between 0 and 1")

    damped, m0, m1 = _truncate_and_damp_2d(moments, broadening_meV, broadening_meV)
    m_max = m1 if m1 > m0 else m0
    num_div = 30 * m_max
    x_grid = [
        -alpha + i * (2.0 * alpha) / (num_div - 1) for i in range(num_div)
    ]

    gr_tbl, dgr_tbl = legacy_green_tables(x_grid, m0, m1)
    ga_tbl = np.conj(gr_tbl)
    a_tbl = np.ascontiguousarray(gr_tbl - ga_tbl, dtype=np.complex128)

    dgr_real = dgr_tbl.real
    dga_tbl = np.array(dgr_real, dtype=np.complex128)
    b_tbl = np.ascontiguousarray(dgr_real + dga_tbl, dtype=np.complex128)

    out = _core.ContractRank2Real(damped, a_tbl, b_tbl)

    scale = moments.window.scale
    shift = moments.window.shift
    pref_sea = float(moments.system_size) * scale * scale / 2.0 / math.pi

    kernel = [0.0] * num_div
    for i in range(num_div):
        k_val = -float(out[i])
        k_val *= pref_sea
        kernel[i] = k_val

    n_out = max(0, num_div - 1)
    energies = np.empty(n_out, dtype=np.float64)
    values = np.empty(n_out, dtype=np.float64)
    acc = 0.0
    for i in range(n_out):
        denerg = x_grid[i + 1] - x_grid[i]
        acc += 0.5 * (kernel[i] + kernel[i + 1]) * denerg
        values[i] = acc
        energies[i] = (x_grid[i] - shift) / scale

    return energies, values


def kubo_bastin_surf(
    moments: _core.ChebyshevMoments,
    broadening_meV: float,
    alpha: float = 0.95,
) -> tuple[npt.NDArray[np.float64], npt.NDArray[np.float64]]:
    """Reconstruct the Fermi-surface contribution of Kubo-Bastin conductivity.

    Solves the problem of evaluating the double-delta Fermi-surface part of
    Kubo-Bastin in the exact operation order of legacy
    ``reconstruct_kubo_bastin_surf``:

    1. Truncate and damp ``mu`` with :func:`jackson_damping_2d` to ``(M0, M1)``,
       ``M = max(M0, M1)``.
    2. Build ``num_div = 30 * M`` uniform points on ``[-alpha, alpha]``.
    3. Contract
       ``out_i = sum_{m0, m1} (delta(x_i, m0) * delta(x_i, m1)) * real(mu[m0, m1])``
       via :func:`~lsquant._core.ContractRank2RealReal`.
    4. Scale ``kernel_i = out_i * (pi * N * scale * scale)`` and report the
       midpoint average ``(kernel_{i+1} + kernel_i) / 2.0``.

    Parameters
    ----------
    moments : lsquant._core.ChebyshevMoments
        Rank-2 Chebyshev moment tensor.
    broadening_meV : float
        Jackson Gaussian broadening in meV (``broadening_meV > 0``).
    alpha : float, optional
        Chebyshev domain cutoff ``(0 < alpha < 1)``, default ``0.95``.

    Returns
    -------
    tuple[numpy.ndarray, numpy.ndarray]
        Pair ``(energies, values)`` of 1D ``float64`` arrays of length
        ``30 * max(M0, M1) - 1``.

    Warnings
    --------
    Raises ``ValueError`` if ``moments.rank != 2``, ``broadening_meV <= 0``, or
    ``!(0 < alpha < 1)``.

    Notes
    -----
    Pitfall: Evaluate ``math.pi * float(moments.system_size) * scale * scale``
    left to right before multiplying ``out_i`` to match C++ ``*=`` evaluation.
    """
    if not (0.0 < alpha < 1.0):
        raise ValueError("kubo_bastin_surf: alpha must lie strictly between 0 and 1")

    damped, m0, m1 = _truncate_and_damp_2d(moments, broadening_meV, broadening_meV)
    m_max = m1 if m1 > m0 else m0
    num_div = 30 * m_max
    x_grid = [
        -alpha + i * (2.0 * alpha) / (num_div - 1) for i in range(num_div)
    ]

    delta_0 = legacy_delta_table(x_grid, m0)
    delta_1 = delta_0 if m1 == m0 else legacy_delta_table(x_grid, m1)
    out = _core.ContractRank2RealReal(damped, delta_0, delta_1)

    scale = moments.window.scale
    shift = moments.window.shift
    pref_surf = math.pi * float(moments.system_size) * scale * scale
    kernel = [float(out[i]) * pref_surf for i in range(num_div)]

    n_out = max(0, num_div - 1)
    energies = np.empty(n_out, dtype=np.float64)
    values = np.empty(n_out, dtype=np.float64)
    for i in range(n_out):
        values[i] = (kernel[i + 1] + kernel[i]) / 2.0
        energies[i] = (x_grid[i] - shift) / scale

    return energies, values


def kubo_fft_grid(
    moments: _core.ChebyshevMoments,
    broadening_meV: float,
    formula: str | bool = "bastin",
    *,
    backend: str | None = None,
) -> tuple[npt.NDArray[np.float64], npt.NDArray[np.float64]]:
    """Reconstruct Kubo-Greenwood or Kubo-Bastin on the Chebyshev FFT-node grid.

    Solves the problem of evaluating ``reconstruct_kubo(mu, obs, GRID_FFT_NODES)``
    on the ``M`` Chebyshev nodes ``x_i = cos(pi * (2 * i + 0.5) / M)`` in the
    exact operation order of legacy ``kuboGreenwoodFromChebmom_FFTgrid`` and
    ``kuboBastinFromChebmom_FFTgrid``:

    1. Truncate and damp ``mu`` with :func:`jackson_damping_2d` to ``(M0, M1)``,
       ``M = max(M0, M1)``.
    2. Build ``x_i = cos(pi * (2 * i + 0.5) / M)`` for ``i = 0 .. M - 1``, and
       when integrating (Bastin) rearrange ``x_grid`` to ascending order via
       ``rearrange_crescent_order``.
    3. Contract with :func:`~lsquant._core.ContractRank2Imag` against
       ``greenR`` (Greenwood) or ``DgreenR`` (Bastin).
    4. Scale by ``pref_ext = N / hw / hw`` and ``-1.0`` (Greenwood) or ``-2.0``
       (Bastin), accumulate, and map energies as ``E_i = x_i * hw + bc``.

    Parameters
    ----------
    moments : lsquant._core.ChebyshevMoments
        Rank-2 Chebyshev moment tensor.
    broadening_meV : float
        Jackson Gaussian broadening in meV (``broadening_meV > 0``).
    formula : str or bool, optional
        ``"bastin"`` / ``"bastin_fft"`` / ``True`` for Kubo-Bastin, or
        ``"greenwood"`` / ``"greenwood_fft"`` / ``False`` for Kubo-Greenwood
        (default ``"bastin"``).

    Returns
    -------
    tuple[numpy.ndarray, numpy.ndarray]
        Pair ``(energies, values)`` of 1D ``float64`` arrays of length
        ``max(M0, M1) - 1``.

    Warnings
    --------
    Raises ``ValueError`` if ``moments.rank != 2``, ``broadening_meV <= 0``, or
    ``formula`` is not recognised.

    Notes
    -----
    Pitfall: Use ``N / hw / hw`` (two sequential divisions by ``half_width``)
    and ``E_i = x_i * hw + bc`` on the FFT-node grid rather than ``scale`` and
    ``shift``; although algebraically equal, they differ by 1 ULP for general
    bandwidths.
    """
    # Opt-in FFT polynomial contraction on the same Chebyshev nodes and Python factors.
    backend = backend or os.environ.get("LSQUANT_FFT_BACKEND", "tables")
    if backend not in ("tables", "fft"):
        raise ValueError("kubo_fft_grid: backend must be 'tables' or 'fft'")
    if isinstance(formula, bool):
        integrate = formula
    else:
        key = formula.strip().lower()
        if key in ("bastin", "bastin_fft"):
            integrate = True
        elif key in ("greenwood", "greenwood_fft"):
            integrate = False
        else:
            raise ValueError(
                f"kubo_fft_grid: unknown formula '{formula}' "
                "(expected 'bastin' or 'greenwood')"
            )

    damped, m0, m1 = _truncate_and_damp_2d(moments, broadening_meV, broadening_meV)
    num_div = m1 if m1 > m0 else m0
    disp = 0.5
    x_grid = [
        math.cos(math.pi * (2 * i + disp) / num_div) for i in range(num_div)
    ]
    order = list(range(num_div))

    if integrate:
        orig = list(x_grid)
        orig_order = list(order)
        for k in range(num_div // 2):
            x_grid[2 * k] = orig[k]
            x_grid[2 * k + 1] = orig[num_div - k - 1]
            order[2 * k] = orig_order[k]
            order[2 * k + 1] = orig_order[num_div - k - 1]
        for k in range(num_div // 2):
            t = x_grid[k]
            x_grid[k] = x_grid[num_div - k - 1]
            x_grid[num_div - k - 1] = t
            order[k], order[num_div - k - 1] = order[num_div - k - 1], order[k]

    if backend == "fft":
        sums = _core.FftMomentContraction(np.ascontiguousarray(damped), num_div)
        out = np.empty(num_div, dtype=np.float64)
        for i, j in enumerate(order):
            x = x_grid[i]
            sine = math.sqrt(1.0 - x * x)
            if integrate:
                green = (-1j * x * complex(sums[0, j]) + sine * complex(sums[1, j])) / sine / sine / sine
            else:
                green = -1j * complex(sums[0, j]) / sine
            out[i] = green.imag / sine / math.pi
    else:
        delta_tbl = legacy_delta_table(x_grid, m0)
        green_r_tbl, dgreen_r_tbl = legacy_green_tables(x_grid, m0, m1)
        b_tbl = dgreen_r_tbl if integrate else green_r_tbl
        out = _core.ContractRank2Imag(damped, delta_tbl, b_tbl)

    hw = moments.window.half_width
    bc = moments.window.center
    pref_ext = float(moments.system_size) / hw / hw
    obs_pref = -2.0 if integrate else -1.0
    kernel = [obs_pref * float(out[i]) * pref_ext for i in range(num_div)]

    n_out = max(0, num_div - 1)
    energies = np.empty(n_out, dtype=np.float64)
    values = np.empty(n_out, dtype=np.float64)
    acc = 0.0
    for i in range(n_out):
        if integrate:
            acc += kernel[i] * (x_grid[i + 1] - x_grid[i])
            values[i] = acc
        else:
            values[i] = kernel[i]
        energies[i] = x_grid[i] * hw + bc

    return energies, values


def time_correlation_curve(
    moments: _core.ChebyshevMoments,
    broadening_meV: float,
    fermi_energy: float,
) -> tuple[npt.NDArray[np.float64], npt.NDArray[np.float64]]:
    """Reconstruct a time-dependent correlation curve at fixed ``fermi_energy``.

    Reproduces the operation order of legacy ``timeCorrelationsFromChebmom``
    byte for byte:

    1. ``eta = 2.0 * broadening_meV / 1000 / BandWidth``,
       ``M_eff = min(ceil(pi / eta), M)``.
    2. Multiply ``mu[m, n]`` by the 1D Jackson kernel ``g(m, M_eff)`` along
       ``m = 0 .. M_eff - 1``.
    3. Rescale ``x_F = (fermi_energy - BandCenter) / HalfWidth`` and evaluate
       ``delta_chebF(x_F, m)`` for ``m = 0 .. M_eff - 1``.
    4. For each time step ``n = 0 .. N_t - 1``, accumulate
       ``val_n = (sum_{m=0}^{M_eff-1} delta_chebF(x_F, m) * mu[m, n].real) * ScaleFactor``
       at ``t_n = n * dt``.
    """
    if moments.rank != 2:
        raise ValueError("time_correlation_curve: moments must have rank 2 (M, N_t)")
    if broadening_meV <= 0.0:
        raise ValueError("time_correlation_curve: broadening_meV must be positive")

    band_width = (moments.window.upper - moments.window.lower) * 1.0
    if band_width <= 0.0:
        raise ValueError(
            "time_correlation_curve: spectral window must have positive width"
        )
    half_width = band_width / 2.0
    band_center = (moments.window.upper + moments.window.lower) * 0.5
    scale_factor = 1.0 / half_width

    eta = 2.0 * float(broadening_meV) / 1000.0 / band_width
    max_mom = int(math.ceil(math.pi / eta))
    num_moms = int(moments.dims[0])
    num_times = int(moments.dims[1])
    if max_mom > num_moms:
        max_mom = num_moms

    g = jackson_damping(max_mom)
    raw = moments.data
    damped = np.empty((max_mom, num_times), dtype=np.complex128)
    for m_idx in range(max_mom):
        g_m = float(g[m_idx])
        for n_idx in range(num_times):
            damped[m_idx, n_idx] = raw[m_idx, n_idx] * g_m

    x_fermi = (float(fermi_energy) - band_center) / half_width
    delta_col = legacy_delta_table([x_fermi], max_mom)[:, 0]

    dt = float(moments.time_step)
    times = np.empty(num_times, dtype=np.float64)
    values = np.empty(num_times, dtype=np.float64)
    for n_idx in range(num_times):
        out_val = 0.0
        for m_idx in range(max_mom):
            out_val += float(delta_col[m_idx]) * float(damped[m_idx, n_idx].real)
        out_val *= scale_factor
        times[n_idx] = float(n_idx) * dt
        values[n_idx] = out_val

    return times, values


def shc(moments, broadening_meV, num_cells, cell_area, *, backend="tables"):
    """Return sigma^z_xy in ``e^2/h, sigma_z Pauli`` from a Bastin tensor.

    Symbolic source: metric_forms.spin_hall_normalisation.
    ``sigma=2*pi*KuboBastin(mu)/(N_k*A_cell)``. The Kubo-Bastin reconstruction
    already multiplies the normalised moments by the Hilbert-space size;
    divide by the number of cells and area, not by the number of orbitals.
    The physical word is Tr[JXSZ delta VY dG-/dE], with J in Pauli units.
    """
    if num_cells <= 0 or cell_area <= 0:
        raise ValueError("shc requires positive num_cells and cell_area")
    energies, values = kubo_bastin(moments, broadening_meV, backend=backend)
    return energies, values * (2*np.pi/(num_cells*cell_area))
