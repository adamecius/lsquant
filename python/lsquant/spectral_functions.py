"""SymPy spectral functions and Chebyshev coefficient table evaluators.

Defines the symbolic :class:`SpectralFunction` representation in the rescaled
Chebyshev variable ``x in (-1, 1)`` and moment order ``m``, the Dirac delta
expansion :data:`DeltaFunction`, the retarded Green's function expansion
:data:`GreenRetarded` and its derivative :data:`GreenRetardedDerivative`, and
the bit-exact legacy evaluators :func:`legacy_delta_table` and
:func:`legacy_green_tables`.
"""

from __future__ import annotations

import math
from typing import Sequence
import numpy as np
import numpy.typing as npt
from scipy.special import jv
import sympy as sp

from . import _core


class SpectralFunction:
    """Symbolic Chebyshev expansion coefficient generator for a function of ``H``.

    Solves the problem of turning an analytic spectral kernel ``c_m(x)`` written
    as a SymPy expression in the dimensionless rescaled energy ``x in (-1, 1)``
    and polynomial index ``m >= 0`` into a dense 2D ``(M, N)`` evaluation table
    for moment contraction. Lambdifies the symbolic expression once at
    construction time and broadcasts over ``(m, x)`` in ``O(M * N)`` time and
    memory.

    Parameters
    ----------
    expr : sympy.Expr
        Symbolic expression for the ``m``-th Chebyshev coefficient at rescaled
        coordinate ``x``.
    x : sympy.Symbol, optional
        Symbol representing the rescaled energy variable in ``(-1, 1)``.
    m : sympy.Symbol, optional
        Symbol representing the Chebyshev polynomial order ``0, 1, ..., M - 1``.

    Warnings
    --------
    Evaluating at ``|x| >= 1`` encounters the ``1 / sqrt(1 - x**2)`` Chebyshev
    weight singularity and yields ``inf`` or ``nan``.

    Notes
    -----
    Pitfall: NumPy's vectorised ``cos``, ``arccos``, and ``exp`` inside
    ``lambdify`` agree with ``libm`` and ``libstdc++`` to machine precision
    (``< 1e-12``) but may differ in the last few ULPs; use
    :func:`legacy_delta_table` and :func:`legacy_green_tables` when byte-for-byte
    identity with legacy C++ streams is required.
    """

    def __init__(
        self,
        expr: sp.Expr,
        x: sp.Symbol | None = None,
        m: sp.Symbol | None = None,
    ) -> None:
        self.x: sp.Symbol = x if x is not None else sp.Symbol("x", real=True)
        self.m: sp.Symbol = (
            m if m is not None else sp.Symbol("m", integer=True, nonnegative=True)
        )
        self.expr: sp.Expr = expr
        self._fn = sp.lambdify((self.m, self.x), self.expr, modules="numpy")

    def table(
        self, x_grid: Sequence[float] | npt.NDArray[np.float64], num_moments: int
    ) -> npt.NDArray[np.float64] | npt.NDArray[np.complex128]:
        """Evaluate the ``(num_moments, len(x_grid))`` Chebyshev coefficient table.

        Broadcasts the lambdified SymPy expression over column vector
        ``m = 0 .. num_moments - 1`` and row vector ``x_grid``. Both axes are
        dimensionless. Complexity is ``O(num_moments * len(x_grid))``.

        Parameters
        ----------
        x_grid : Sequence[float] or numpy.ndarray
            Rescaled evaluation points ``x_i in (-1, 1)`` of length ``N``.
        num_moments : int
            Number of Chebyshev moments ``M > 0``.

        Returns
        -------
        numpy.ndarray
            Contiguous 2D ``float64`` (for real kernels) or ``complex128`` (for
            complex kernels) array of shape ``(M, N)``.

        Warnings
        --------
        Raises ``ValueError`` if ``num_moments <= 0``.

        Notes
        -----
        Pitfall: If ``self.expr`` is independent of ``m``, ``lambdify`` returns
        a row array; ``np.broadcast_to`` ensures the returned shape is always
        ``(M, N)``.
        """
        if num_moments <= 0:
            raise ValueError("SpectralFunction.table: num_moments must be positive")
        x_arr = np.asarray(x_grid, dtype=np.float64)[np.newaxis, :]
        m_arr = np.arange(num_moments, dtype=np.float64)[:, np.newaxis]
        evaluated = self._fn(m_arr, x_arr)
        out_dtype = np.complex128 if np.iscomplexobj(evaluated) else np.float64
        broadcasted = np.broadcast_to(evaluated, (num_moments, x_arr.shape[1]))
        return np.ascontiguousarray(broadcasted, dtype=out_dtype)

    def evaluate(
        self,
        x_values: Sequence[float] | npt.NDArray[np.float64] | float,
        param: Sequence[float] | npt.NDArray[np.float64] | float,
    ) -> npt.NDArray[np.float64] | npt.NDArray[np.complex128]:
        """Evaluate the lambdified SymPy expression at ``(param, x_values)``."""
        x_arr = np.asarray(x_values, dtype=np.float64)
        p_arr = np.asarray(param, dtype=np.float64)
        evaluated = self._fn(p_arr, x_arr)
        out_dtype = np.complex128 if np.iscomplexobj(evaluated) else np.float64
        return np.ascontiguousarray(evaluated, dtype=out_dtype)


_x_sym = sp.Symbol("x", real=True)
_m_sym = sp.Symbol("m", integer=True, nonnegative=True)
_t_sym = sp.Symbol("t", real=True)
_eta_sym = sp.Symbol("eta", real=True, positive=True)

GaussianDelta = SpectralFunction(
    sp.exp(-(_x_sym**2) / (2 * _eta_sym**2)) / (sp.sqrt(2 * sp.pi) * _eta_sym),
    x=_x_sym,
    m=_eta_sym,
)

DeltaFunction = SpectralFunction(
    sp.cos(_m_sym * sp.acos(_x_sym)) / (sp.pi * sp.sqrt(1 - _x_sym**2)),
    x=_x_sym,
    m=_m_sym,
)

GreenRetarded = SpectralFunction(
    -sp.I * sp.exp(-sp.I * _m_sym * sp.acos(_x_sym)) / sp.sqrt(1 - _x_sym**2),
    x=_x_sym,
    m=_m_sym,
)

GreenRetardedDerivative = SpectralFunction(
    sp.diff(GreenRetarded.expr, _x_sym),
    x=_x_sym,
    m=_m_sym,
)


class TimeEvolutionFunction(SpectralFunction):
    """Symbolic time-evolution operator ``exp(-I * x * t)`` and Bessel expansion.

    Solves the problem of representing the unitary time-evolution propagator
    ``U(t) = exp(-i H~ t)`` on the rescaled Chebyshev domain ``x in [-1, 1]``
    together with its Jacobi-Anger Bessel series coefficients
    ``c_0 = J_0(t)``, ``c_n = 2 (-i)^n J_n(t)`` evaluated via
    :func:`scipy.special.jv`.
    """

    def __init__(self) -> None:
        self.t: sp.Symbol = _t_sym
        self.kernel_expr: sp.Expr = sp.exp(-sp.I * _x_sym * _t_sym)
        coeff_expr = sp.Piecewise(
            (sp.besselj(_m_sym, _x_sym), sp.Eq(_m_sym, 0)),
            (2 * (-sp.I) ** _m_sym * sp.besselj(_m_sym, _x_sym), True),
        )
        super().__init__(coeff_expr, x=_x_sym, m=_m_sym)
        self.expr = self.kernel_expr

    def table(
        self, x_grid: Sequence[float] | npt.NDArray[np.float64], num_moments: int
    ) -> npt.NDArray[np.complex128]:
        """Evaluate the ``(num_moments, len(x_grid))`` Jacobi-Anger Bessel table."""
        if num_moments <= 0:
            raise ValueError("TimeEvolution.table: num_moments must be positive")
        x_arr = np.asarray(x_grid, dtype=np.float64)[np.newaxis, :]
        m_arr = np.arange(num_moments, dtype=np.float64)[:, np.newaxis]
        j_vals = jv(m_arr, x_arr)
        weights = np.where(m_arr == 0.0, 1.0, 2.0) * ((-1j) ** m_arr)
        return np.ascontiguousarray(weights * j_vals, dtype=np.complex128)

    def series_coefficients(
        self, x: float, tol: float = 1e-15
    ) -> npt.NDArray[np.complex128]:
        """Compute Bessel series coefficients ``c_n(x)`` via :func:`scipy.special.jv`.

        Evaluates ``c_0 = J_0(x)`` and ``c_n = (-i)^n * 2 * J_n(x)`` for ``n >= 1``,
        stopping at the first ``n >= 1`` where ``|J_n(x)| <= tol`` (that term
        excluded).
        """
        if tol <= 0.0:
            raise ValueError("series_coefficients: tol must be positive")
        x_val = float(x)
        tol_val = float(tol)
        mom_cutoff = 2.0 * tol_val
        n = 0
        n_ip = complex(1.0, 0.0)
        minus_i = complex(0.0, -1.0)
        jn = 0.5 * float(jv(0, x_val))
        coeffs: list[complex] = []
        while mom_cutoff > tol_val:
            coeffs.append(n_ip * complex(2.0, 0.0) * jn)
            n_ip *= minus_i
            n += 1
            jn = float(jv(n, x_val))
            mom_cutoff = abs(jn)
        return np.asarray(coeffs, dtype=np.complex128)


TimeEvolution = TimeEvolutionFunction()


def bessel_series_coefficients(
    x: float, tol: float = 1e-15
) -> npt.NDArray[np.complex128]:
    """Compute Jacobi-Anger Bessel coefficients ``c_n(x)`` via SciPy ``jv``."""
    return TimeEvolution.series_coefficients(x, tol=tol)


def legacy_bessel_coefficients(
    x: float, tol: float = 1e-15
) -> npt.NDArray[np.complex128]:
    """Compute Jacobi-Anger Bessel coefficients ``c_n(x)`` via GSL in ``_core.compat``."""
    if tol <= 0.0:
        raise ValueError("legacy_bessel_coefficients: tol must be positive")
    return _core.compat.bessel_series_coefficients(float(x), float(tol))


def legacy_delta_table(
    x: Sequence[float] | npt.NDArray[np.float64], num_moments: int
) -> npt.NDArray[np.float64]:
    """Compute the ``(M_eff, N)`` delta coefficient table in the exact legacy order.

    Solves the problem of reproducing C++ ``delta_chebF(x, m)`` bit for bit
    without SIMD ULP drift by calling Python's ``math.sqrt``, ``math.acos``,
    ``math.cos``, and ``math.pi`` in the exact expression order:

    .. code-block:: python

        f0 = math.sqrt(1.0 - xi * xi)
        fm = math.cos(m * math.acos(xi))
        val = fm / f0 / math.pi

    All inputs and outputs are dimensionless in the rescaled Chebyshev domain.
    Time and space complexity are ``O(num_moments * len(x))``.

    Parameters
    ----------
    x : Sequence[float] or numpy.ndarray
        Rescaled coordinate grid ``x_i in (-1, 1)`` of length ``N``.
    num_moments : int
        Effective number of Chebyshev moments ``M_eff > 0``.

    Returns
    -------
    numpy.ndarray
        2D C-contiguous ``float64`` array of shape ``(num_moments, len(x))``.

    Warnings
    --------
    Raises ``ValueError`` if ``num_moments <= 0`` or if any ``|x_i| >= 1.0``.

    Notes
    -----
    Pitfall: Keep the two sequential divisions ``fm / f0 / math.pi`` rather than
    ``fm / (f0 * math.pi)``; multiplying the denominator first alters rounding
    at the last bit.
    """
    if num_moments <= 0:
        raise ValueError("legacy_delta_table: num_moments must be positive")

    x_list = [float(val) for val in x]
    num_points = len(x_list)
    table = np.empty((num_moments, num_points), dtype=np.float64)
    pi = math.pi
    cos = math.cos
    acos = math.acos
    sqrt = math.sqrt

    for i, xi in enumerate(x_list):
        f0 = sqrt(1.0 - xi * xi)
        acos_x = acos(xi)
        for m_idx in range(num_moments):
            fm = cos(float(m_idx) * acos_x)
            table[m_idx, i] = fm / f0 / pi

    return table


def legacy_green_tables(
    x: Sequence[float] | npt.NDArray[np.float64],
    num_moments_0: int,
    num_moments_1: int | None = None,
) -> tuple[npt.NDArray[np.complex128], npt.NDArray[np.complex128]]:
    """Compute the legacy ``(greenR, DgreenR)`` tables via ``_core.compat``.

    Solves the problem of reproducing C++ ``greenR_chebF(x, m)`` and
    ``DgreenR_chebF(x, m)`` bit for bit by delegating to
    :func:`lsquant._core.compat.legacy_green_tables`, which uses ``libstdc++``
    ``std::pow(std::complex<double>, double)`` and three sequential divisions by
    ``f0 = sqrt(1 - x*x)``. Complexity is
    ``O((num_moments_0 + num_moments_1) * len(x))``.

    Parameters
    ----------
    x : Sequence[float] or numpy.ndarray
        Rescaled coordinate grid ``x_i in (-1, 1)`` of length ``N``.
    num_moments_0 : int
        Number of Chebyshev moments ``M_0 > 0`` for the ``greenR`` table.
    num_moments_1 : int or None, optional
        Number of Chebyshev moments ``M_1 > 0`` for the ``DgreenR`` table
        (defaults to ``num_moments_0``).

    Returns
    -------
    tuple[numpy.ndarray, numpy.ndarray]
        Pair ``(green_r, dgreen_r)`` of 2D C-contiguous ``complex128`` arrays of
        shapes ``(num_moments_0, len(x))`` and ``(num_moments_1, len(x))``.

    Warnings
    --------
    Raises ``ValueError`` if ``num_moments_0 <= 0`` or ``num_moments_1 <= 0``.

    Notes
    -----
    Pitfall: Python's built-in complex ``**`` operator uses repeated
    multiplication for small integer exponents and does not match
    ``std::pow(std::complex<double>, double)`` at the last bit.
    """
    m1 = num_moments_0 if num_moments_1 is None else num_moments_1
    if num_moments_0 <= 0 or m1 <= 0:
        raise ValueError("legacy_green_tables: moment counts must be positive")
    x_arr = np.asarray(x, dtype=np.float64)
    return _core.compat.legacy_green_tables(x_arr, int(num_moments_0), int(m1))
