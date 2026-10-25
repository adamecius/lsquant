"""Bytecode plan constructors, spectral window helpers, and execution routines.

Defines the rank-1 density-of-states plan (:func:`dos_plan`) and generalised
rank-2 Kubo correlation plan (:func:`kubo_plan`) as data over the
five-instruction C++ engine grammar, together with
:func:`window_from_gershgorin` and :func:`run_plan`.
"""

from __future__ import annotations

import cmath
from typing import Sequence
import numpy as np
import numpy.typing as npt

from . import _core
from .spectral_functions import (
    bessel_series_coefficients,
    legacy_bessel_coefficients,
)

HBAR: float = 0.6582119624


def _kpm_weights_1d(num_moments: int, num_states: int) -> list[complex]:
    """Build the 1D KPM ``(2 - delta_{m,0}) / num_states`` weight vector."""
    weights: list[complex] = []
    base = 2.0 / float(num_states)
    for m in range(num_moments):
        scal = base * 0.5 if m == 0 else base
        weights.append(complex(scal, 0.0))
    return weights


def time_step_from_tmax(num_steps: int, tmax: float) -> float:
    """Compute the legacy uniform time step ``dt = tmax / (num_steps - 1)``."""
    if num_steps <= 1:
        raise ValueError("time_step_from_tmax: num_steps must be > 1")
    return float(tmax) / float(num_steps - 1)


def time_evolution_series(
    window: _core.SpectralWindow,
    dt: float,
    tol: float = 1e-15,
    use_legacy_bessel: bool = False,
) -> tuple[npt.NDArray[np.complex128], complex]:
    """Compute Bessel series coefficients and band-center phase for ``U(dt)``.

    Reproduces the ``MomentsTD::Evolve`` parameters verbatim:
    ``x = dt * (BandWidth / 2.0 / 1.0 / HBAR)``,
    ``exp0 = exp(-i * (BandCenter / HBAR) * dt)``.
    """
    band_width = (window.upper - window.lower) * 1.0
    if band_width <= 0.0:
        raise ValueError("time_evolution_series: window must have positive width")
    half_width = band_width / 2.0
    band_center = (window.upper + window.lower) * 0.5
    cheb_freq = half_width / 1.0 / HBAR
    cheb_freq_0 = band_center / HBAR
    x = float(dt) * cheb_freq
    coeffs = (
        legacy_bessel_coefficients(x, tol=tol)
        if use_legacy_bessel
        else bessel_series_coefficients(x, tol=tol)
    )
    phase = cmath.exp(complex(0.0, -1.0) * cheb_freq_0 * float(dt))
    return coeffs, phase


def dos_plan(
    num_moments: int,
    num_states: int,
    word: str,
    *,
    operator_index: int | None = None,
) -> _core.Plan:
    """Build the rank-1 density-of-states bytecode :class:`~lsquant._core.Plan`.

    Solves the problem of expressing ``Tr[A * delta(E - H)]`` as a ket-only
    Chebyshev expansion ``Expand(M)`` (optionally followed by
    ``ApplyOperator(operator_index)`` when ``A != 1``) and ``Trace(w)`` with the
    legacy per-state KPM weights ``w[0] = 1 / num_states`` and
    ``w[m] = 2 / num_states`` for ``m > 0``. Complexity is ``O(num_moments)``.

    Parameters
    ----------
    num_moments : int
        Number of Chebyshev moments ``M > 0``.
    num_states : int
        Number of trace vectors ``N_states > 0`` used for weight normalisation.
    word : str
        Symbolic correlator word recorded in the plan and moment tensor metadata
        (e.g. ``"1 delta"``).
    operator_index : int or None, optional
        Optional post-expansion operator index for ``Tr[A delta(E-H)]``.

    Returns
    -------
    lsquant._core.Plan
        Compiled bytecode plan with ``plan.word == word``.

    Warnings
    --------
    Raises ``ValueError`` if ``num_moments <= 0``, ``num_states <= 0``,
    ``operator_index < 0``, or ``word`` is empty.

    Notes
    -----
    Pitfall: Because the ``(2 - delta_{m,0}) / num_states`` weights are folded
    into :class:`~lsquant._core.Trace` per state, do not apply KPM weights a
    second time after running this plan.
    """
    if num_moments <= 0:
        raise ValueError("dos_plan: num_moments must be positive")
    if num_states <= 0:
        raise ValueError("dos_plan: num_states must be positive")
    if not word:
        raise ValueError("dos_plan: word must be non-empty")
    if operator_index is not None and operator_index < 0:
        raise ValueError("dos_plan: operator_index must be non-negative")

    weights: list[complex] = []
    base = 2.0 / float(num_states)
    for m in range(num_moments):
        scal = base * 0.5 if m == 0 else base
        weights.append(complex(scal, 0.0))

    code: list[object] = [_core.Expand(num_moments)]
    if operator_index is not None:
        code.append(_core.ApplyOperator(operator_index))
    code.append(_core.Trace(weights))
    return _core.Plan(code, word)


def kubo_plan(
    num_moments: int,
    *,
    left_index: int,
    right_index: int,
    word: str,
) -> _core.Plan:
    """Build the rank-2 Kubo correlation bytecode :class:`~lsquant._core.Plan`.

    Solves the problem of expressing the two-operator correlation trace
    ``mu_{m,n} = Tr[O_left T_m(H~) O_right T_n(H~)]`` in the exact ket/bra
    evaluation order of legacy ``CorrelationExpansionMoments(OPL, OPR)``:

    - ket: ``Expand(M)``, ``ApplyOperator(right_index)``
      (computes ``O_right T_n |r>``);
    - ``MirrorToBra()``;
    - bra: ``ApplyOperator(left_index)``, ``Expand(M)``
      (computes ``T_m O_left |r>``);
    - ``Trace([1])`` (accumulates raw inner products across states).

    Complexity is ``O(1)``.

    Parameters
    ----------
    num_moments : int
        Number of Chebyshev moments ``M > 0`` along each of the two axes.
    left_index : int
        Keyword-only index of the bra-side operator ``O_left`` in the engine's
        operator list.
    right_index : int
        Keyword-only index of the ket-side operator ``O_right`` in the engine's
        operator list.
    word : str
        Keyword-only correlator word built by the caller from its operator names
        (e.g. ``f"{left} delta {right} delta"``).

    Returns
    -------
    lsquant._core.Plan
        Compiled bytecode plan with ``plan.word == word``.

    Warnings
    --------
    Raises ``ValueError`` if ``num_moments <= 0``, operator indices are
    negative, or ``word`` is empty.

    Notes
    -----
    In ``Tr[OP_L delta(E-H) OP_R delta(E-H)]``, ``left_index`` seeds the bra
    recursion and ``right_index`` is applied to the ket batch.
    """
    if num_moments <= 0:
        raise ValueError("kubo_plan: num_moments must be positive")
    if left_index < 0 or right_index < 0:
        raise ValueError("kubo_plan: operator indices must be non-negative")
    if not word:
        raise ValueError("kubo_plan: word must be non-empty")

    code = [
        _core.Expand(num_moments),
        _core.ApplyOperator(right_index),
        _core.MirrorToBra(),
        _core.ApplyOperator(left_index),
        _core.Expand(num_moments),
        _core.Trace([complex(1.0, 0.0)]),
    ]
    return _core.Plan(code, word)


def msd_plan(
    num_moments: int,
    num_steps: int,
    num_states: int,
    coeffs: Sequence[complex] | npt.NDArray[np.complex128],
    *,
    velocity_index: int = 0,
    phase: complex = 1.0 + 0.0j,
    word: str = "[X, U(t)]^+ delta [X, U(t)]",
) -> _core.Plan:
    """Build the mean-square displacement bytecode :class:`~lsquant._core.Plan`.

    Solves the problem of evaluating ``Tr[[X, U(t)]^+ delta(E - H) [X, U(t)]]``
    via the coupled ``(Ket = MSDWF, Aux = PhiR)`` time-evolution recursion of
    legacy ``MeanSquareDisplacement`` and ``MomentsTD::MSD_Evolve``:

    - Seed preparation: ``ApplySeries([])`` zeros ``seed[Ket]`` (``MSDWF``) while
      ``seed[Aux]`` starts as ``|r>``, followed by ``MirrorToBra()``. An empty
      ``coeffs`` is the zero map (``v <- 0``), which is how a seed register is
      cleared, and ``msd_plan`` relies on it.
    - ``TimeAxis(num_steps, evolve=[EvolveRegister(Ket, coeffs, phase),
      EvolveDerived(Ket, Aux, coeffs, velocity_index, phase),
      EvolveRegister(Aux, coeffs, phase)])``.
    - Body: ``Expand(num_moments)``, ``MirrorToBra()`` (sets ``workBra = workKet``
      for that step), ``Trace(w)``.
    """
    if num_moments <= 0:
        raise ValueError("msd_plan: num_moments must be positive")
    if num_steps <= 0:
        raise ValueError("msd_plan: num_steps must be positive")
    if num_states <= 0:
        raise ValueError("msd_plan: num_states must be positive")
    if velocity_index < 0:
        raise ValueError("msd_plan: velocity_index must be non-negative")
    if not word:
        raise ValueError("msd_plan: word must be non-empty")

    weights = _kpm_weights_1d(num_moments, num_states)
    evolve = [
        _core.EvolveRegister(_core.Register.Ket, coeffs, phase),
        _core.EvolveDerived(
            _core.Register.Ket,
            _core.Register.Aux,
            coeffs,
            velocity_index,
            phase,
        ),
        _core.EvolveRegister(_core.Register.Aux, coeffs, phase),
    ]
    code = [
        _core.ApplySeries([], complex(1.0, 0.0)),
        _core.MirrorToBra(),
        _core.TimeAxis(num_steps, evolve),
        _core.Expand(num_moments),
        _core.MirrorToBra(),
        _core.Trace(weights),
    ]
    return _core.Plan(code, word)


def projected_evolution_plan(
    num_moments: int,
    num_steps: int,
    num_states: int,
    coeffs: Sequence[complex] | npt.NDArray[np.complex128],
    *,
    operator_index: int,
    projector_index: int,
    phase: complex = 1.0 + 0.0j,
    word: str = "P U(t)^+ O U(t) delta P",
) -> _core.Plan:
    """Build the time-evolved projected operator bytecode :class:`~lsquant._core.Plan`.

    Solves the problem of evaluating ``Tr[P U(t)^+ O U(t) delta(E - H) P]`` in
    the exact evaluation order of legacy ``TimeEvolvedProjectedOperator(OP, P)``:

    - Seed preparation: ``ApplyOperator(projector_index)`` sets
      ``seed[Ket] = P |r>``, followed by ``MirrorToBra()`` so
      ``seed[Bra] = P |r>``.
    - ``TimeAxis(num_steps, evolve=[EvolveRegister(Bra, coeffs, phase),
      EvolveRegister(Ket, coeffs, phase)])``.
    - Body: ``Expand(num_moments)``, ``ApplyOperator(operator_index)`` (post),
      ``Trace(w)``.
    """
    if num_moments <= 0:
        raise ValueError("projected_evolution_plan: num_moments must be positive")
    if num_steps <= 0:
        raise ValueError("projected_evolution_plan: num_steps must be positive")
    if num_states <= 0:
        raise ValueError("projected_evolution_plan: num_states must be positive")
    if operator_index < 0 or projector_index < 0:
        raise ValueError(
            "projected_evolution_plan: operator indices must be non-negative"
        )
    if not word:
        raise ValueError("projected_evolution_plan: word must be non-empty")

    weights = _kpm_weights_1d(num_moments, num_states)
    evolve = [
        _core.EvolveRegister(_core.Register.Bra, coeffs, phase),
        _core.EvolveRegister(_core.Register.Ket, coeffs, phase),
    ]
    code = [
        _core.ApplyOperator(projector_index),
        _core.MirrorToBra(),
        _core.TimeAxis(num_steps, evolve),
        _core.Expand(num_moments),
        _core.ApplyOperator(operator_index),
        _core.Trace(weights),
    ]
    return _core.Plan(code, word)


def time_correlation_plan(
    num_moments: int,
    num_steps: int,
    num_states: int,
    coeffs: Sequence[complex] | npt.NDArray[np.complex128],
    *,
    left_index: int,
    right_index: int,
    phase: complex = 1.0 + 0.0j,
    word: str = "A delta U(t)^+ B U(t)",
) -> _core.Plan:
    """Build the time-dependent two-operator correlation :class:`~lsquant._core.Plan`.

    Solves the problem of evaluating ``Tr[O_left U(t) delta(E - H) O_right U(t)^+]``
    in the exact evaluation order of legacy ``TimeDependentCorrelations(OPL, OPR)``:

    - Seed preparation: ``MirrorToBra()``, ``ApplyOperator(left_index)`` so
      ``seed[Bra] = O_left |r>`` while ``seed[Ket] = |r>``.
    - ``TimeAxis(num_steps, evolve=[EvolveRegister(Bra, coeffs, phase),
      EvolveRegister(Ket, coeffs, phase)])``.
    - Body: ``ApplyOperator(right_index)`` (pre, on a copy of ``seed[Ket]``),
      ``Expand(num_moments)``, ``Trace(w)``.
    """
    if num_moments <= 0:
        raise ValueError("time_correlation_plan: num_moments must be positive")
    if num_steps <= 0:
        raise ValueError("time_correlation_plan: num_steps must be positive")
    if num_states <= 0:
        raise ValueError("time_correlation_plan: num_states must be positive")
    if left_index < 0 or right_index < 0:
        raise ValueError("time_correlation_plan: operator indices must be non-negative")
    if not word:
        raise ValueError("time_correlation_plan: word must be non-empty")

    weights = _kpm_weights_1d(num_moments, num_states)
    evolve = [
        _core.EvolveRegister(_core.Register.Bra, coeffs, phase),
        _core.EvolveRegister(_core.Register.Ket, coeffs, phase),
    ]
    code = [
        _core.MirrorToBra(),
        _core.ApplyOperator(left_index),
        _core.TimeAxis(num_steps, evolve),
        _core.ApplyOperator(right_index),
        _core.Expand(num_moments),
        _core.Trace(weights),
    ]
    return _core.Plan(code, word)


def window_from_gershgorin(
    hamiltonian: _core.SparseMatrix, pad: float = 0.10
) -> _core.SpectralWindow:
    """Compute the symmetric Gershgorin :class:`~lsquant._core.SpectralWindow` for ``hamiltonian``.

    Solves the problem of enclosing the spectrum of ``hamiltonian`` when no
    explicit ``BOUNDS`` or ``.desc`` sidecar is supplied, reproducing the legacy
    ``utility::SpectralBounds`` fallback ``[-rho * (1 + pad), +rho * (1 + pad)]``
    where ``rho = max_i sum_j |H_{ij}|`` is the infinity norm of ``hamiltonian``.
    Complexity is ``O(dim + nnz)``.

    Parameters
    ----------
    hamiltonian : lsquant._core.SparseMatrix
        Compressed sparse Hamiltonian matrix.
    pad : float, optional
        Relative safety padding fraction (default ``0.10``).

    Returns
    -------
    lsquant._core.SpectralWindow
        Symmetric spectral window ``[-rho * (1 + pad), +rho * (1 + pad)]``.

    Warnings
    --------
    Raises ``ValueError`` if ``pad < 0``.

    Notes
    -----
    Pitfall: Ensure ``hamiltonian`` has undergone :meth:`hermitian_fold` (as
    :func:`~lsquant._core.ReadCsrText` does) before computing Gershgorin row
    sums so that the resulting ``BandWidth`` matches legacy files to the last
    bit.
    """
    bounds = _core.GershgorinBounds(hamiltonian, float(pad))
    return _core.SpectralWindow(bounds.lower, bounds.upper)


def run_plan(
    plan: _core.Plan,
    hamiltonian: _core.SparseMatrix | _core.LinearOperator,
    window: _core.SpectralWindow,
    states: _core.StateGenerator,
    operators: Sequence[_core.LinearOperator] = (),
    batch_size: int | None = None,
    symmetrise: bool = True,
    dt: float | None = None,
) -> _core.ChebyshevMoments:
    """Execute a bytecode :class:`~lsquant._core.Plan` and populate tensor metadata.

    Rescales ``hamiltonian`` to ``window``, runs ``plan`` through
    :class:`~lsquant._core.MomentEngine`, and records ``window``,
    ``kpm_weighted``, and ``symmetrised`` on the returned
    :class:`~lsquant._core.ChebyshevMoments`. For rank-2 plans with unit trace
    weight, applies :func:`lsquant._core.ApplyKpmWeights`.

    Parameters
    ----------
    plan : lsquant._core.Plan
        Compiled bytecode instruction sequence.
    hamiltonian : lsquant._core.SparseMatrix
        Physical tight-binding Hamiltonian in CSR format.
    window : lsquant._core.SpectralWindow
        Spectral window ``[lower, upper]`` enclosing the spectrum of ``hamiltonian``.
    states : lsquant._core.StateGenerator
        Source-state generator for trace evaluation.
    operators : Sequence[lsquant._core.LinearOperator], optional
        Positional list of static operators referenced by ``ApplyOperator``.
    batch_size : int or None, optional
        Chebyshev batch size (defaults to the first ``Expand`` moment count).
    symmetrise : bool, optional
        Whether to Hermitian-symmetrise rank-2 moment tensors (default ``True``).

    Returns
    -------
    lsquant._core.ChebyshevMoments
        Accumulated and weighted Chebyshev moment tensor.

    Warnings
    --------
    Raises ``ValueError`` or ``RuntimeError`` if dimensions mismatch or ``plan``
    is invalid.

    Notes
    -----
    Pitfall: ``states`` is automatically reset by ``MomentEngine.run`` before
    traversing the state sequence.
    """
    rescaled_h = _core.RescaledOperator(hamiltonian, window)
    if batch_size is None:
        batch_size = 1
        for instr in plan.code:
            if isinstance(instr, _core.Expand):
                batch_size = instr.num_moments
                break

    engine = _core.MomentEngine(rescaled_h, list(operators), states, batch_size)
    mu = engine.run(plan)
    mu.window = window
    mu.system_size = hamiltonian.dimension
    mu.number_of_states = int(states.number_of_states)
    mu.state_label = states.label
    mu.word = plan.word
    if dt is not None:
        mu.time_step = float(dt)

    if mu.rank == 1 or mu.number_of_time_steps > 0:
        mu.kpm_weighted = True
        mu.symmetrised = False
    elif mu.rank == 2 and not mu.kpm_weighted:
        _core.ApplyKpmWeights(mu, int(states.number_of_states), symmetrise)

    return mu
