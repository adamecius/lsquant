"""Chebyshev truncation damping kernels.

Provides 1D and 2D Jackson kernel damping factors computed with standard library
``math`` trigonometric functions in the exact legacy operation order to
guarantee bit-identical moment damping.
"""

from __future__ import annotations

import math
import numpy as np
import numpy.typing as npt


def jackson_damping(num_moments: int) -> npt.NDArray[np.float64]:
    """Compute the Jackson kernel damping weights ``g_m`` for ``m = 0 .. M - 1``.

    Solves the problem of suppressing Gibbs oscillations in a truncated
    Chebyshev series of length ``M = num_moments`` while preserving non-negativity
    of spectral densities. Evaluates the closed-form Jackson weights using
    Python's ``math`` module in the exact operation order of legacy
    ``chebyshev::Moments::JacksonKernel``:

    .. math::

        \\phi_J = \\frac{\\pi}{M + 1}, \\qquad
        g_m = \\left((M - m + 1)\\cos(m\\phi_J)
              + \\frac{\\sin(m\\phi_J)}{\\tan(\\phi_J)}\\right)\\frac{\\phi_J}{\\pi}.

    All quantities are dimensionless and lie in ``(0, 1]`` with ``g_0 = 1``.
    Time and space complexity are ``O(num_moments)``.

    Parameters
    ----------
    num_moments : int
        Effective number of Chebyshev moments ``M > 0``.

    Returns
    -------
    numpy.ndarray
        1D ``float64`` array of length ``num_moments`` containing ``g_m``.

    Warnings
    --------
    Raises ``ValueError`` if ``num_moments <= 0``.

    Notes
    -----
    Pitfall: Do not replace ``math.cos``, ``math.sin``, and ``math.tan`` with
    vectorised NumPy ufuncs; SIMD trigonometric implementations can differ by
    1 ULP and break byte-for-byte golden curve reproduction.
    """
    if num_moments <= 0:
        raise ValueError("jackson_damping: num_moments must be positive")

    mom = float(num_moments)
    phi = math.pi / (mom + 1.0)
    weights = np.empty(num_moments, dtype=np.float64)
    for m_idx in range(num_moments):
        m = float(m_idx)
        weights[m_idx] = (
            ((mom - m + 1.0) * math.cos(phi * m) + math.sin(phi * m) / math.tan(phi))
            * phi
            / math.pi
        )
    return weights


def jackson_damping_2d(
    num_moments_0: int, num_moments_1: int
) -> npt.NDArray[np.float64]:
    """Compute the 2D Jackson damping matrix ``g(m0, M0) * g(m1, M1)``.

    Solves the problem of damping a truncated ``M0 x M1`` rank-2 Chebyshev
    moment tensor in the exact floating-point operation order of the legacy 2D
    Jackson damping loop, where the two ``double``
    damping weights ``g(m0, M0)`` and ``g(m1, M1)`` are multiplied together
    before scaling the complex moment ``mu[m0, m1]``. Complexity is
    ``O(num_moments_0 * num_moments_1)`` time and space.

    Parameters
    ----------
    num_moments_0 : int
        Effective number of Chebyshev moments ``M0 > 0`` along axis 0.
    num_moments_1 : int
        Effective number of Chebyshev moments ``M1 > 0`` along axis 1.

    Returns
    -------
    numpy.ndarray
        2D ``float64`` array of shape ``(num_moments_0, num_moments_1)`` with
        entries ``g(m0, M0) * g(m1, M1)``.

    Warnings
    --------
    Raises ``ValueError`` if ``num_moments_0 <= 0`` or ``num_moments_1 <= 0``.

    Notes
    -----
    Pitfall: Multiply ``g0[m0] * g1[m1]`` in ``float64`` first and then multiply
    ``mu[m0, m1]`` by that product; scaling ``(mu[m0, m1] * g0[m0]) * g1[m1]``
    sequentially alters rounding in the last bit.
    """
    if num_moments_0 <= 0 or num_moments_1 <= 0:
        raise ValueError("jackson_damping_2d: moment counts must be positive")

    g0 = jackson_damping(num_moments_0)
    g1 = jackson_damping(num_moments_1) if num_moments_1 != num_moments_0 else g0
    weights = np.empty((num_moments_0, num_moments_1), dtype=np.float64)
    for m0 in range(num_moments_0):
        g_m0 = float(g0[m0])
        for m1 in range(num_moments_1):
            weights[m0, m1] = g_m0 * float(g1[m1])
    return weights
