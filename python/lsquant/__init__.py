"""LSQUANT Python library for Kernel Polynomial Method calculations.

Exposes the C++ execution engine via :mod:`lsquant._core` alongside SymPy
spectral functions, damping kernels, moment plan builders, and rank-1/rank-2
curve reconstructions.
"""

from __future__ import annotations

import sys as _sys

from .runtime import configure_threads as _configure_threads

_configure_threads()

try:
    from . import _core
except ImportError:
    import _core as _core  # type: ignore[no-redef]

    _sys.modules[f"{__name__}._core"] = _core

from .runtime import configure_kernels as _configure_kernels

_configure_kernels(_core)

from . import bases, cli, kernels, models, moments, reconstruct, runfile, spectral_functions, words
from .bases import BlochMetricBasis, EigenMetricBasis, BlochBasis, EigenBasis, RealBasis, make_basis
from .runfile import RunFile, execute_runfile, inspect_target, reconstruct_target
from .kernels import jackson_damping, jackson_damping_2d
from .moments import (
    HBAR,
    dos_plan,
    kubo_plan,
    msd_plan,
    projected_evolution_plan,
    run_plan,
    time_correlation_plan,
    time_evolution_series,
    time_step_from_tmax,
    window_from_gershgorin,
)
from .reconstruct import (
    dos_curve,
    shc,
    kubo_bastin,
    kubo_bastin_sea,
    kubo_bastin_surf,
    kubo_fft_grid,
    kubo_greenwood,
    time_correlation_curve,
)
from .spectral_functions import (
    DeltaFunction,
    GaussianDelta,
    GreenRetarded,
    GreenRetardedDerivative,
    SpectralFunction,
    TimeEvolution,
    TimeEvolutionFunction,
    bessel_series_coefficients,
    legacy_bessel_coefficients,
    legacy_delta_table,
    legacy_green_tables,
)

__all__ = [
    "_core",
    "bases",
    "kernels",
    "models",
    "moments",
    "reconstruct",
    "runfile",
    "spectral_functions",
    "words",
    "cli",
    "RunFile",
    "execute_runfile",
    "reconstruct_target",
    "inspect_target",
    "RealBasis",
    "BlochBasis",
    "BlochMetricBasis",
    "EigenMetricBasis",
    "shc",
    "EigenBasis",
    "make_basis",
    "jackson_damping",
    "jackson_damping_2d",
    "HBAR",
    "dos_plan",
    "kubo_plan",
    "msd_plan",
    "projected_evolution_plan",
    "time_correlation_plan",
    "time_evolution_series",
    "time_step_from_tmax",
    "window_from_gershgorin",
    "run_plan",
    "dos_curve",
    "kubo_greenwood",
    "kubo_bastin",
    "kubo_bastin_sea",
    "kubo_bastin_surf",
    "kubo_fft_grid",
    "time_correlation_curve",
    "SpectralFunction",
    "DeltaFunction",
    "GreenRetarded",
    "GreenRetardedDerivative",
    "TimeEvolutionFunction",
    "TimeEvolution",
    "bessel_series_coefficients",
    "legacy_bessel_coefficients",
    "legacy_delta_table",
    "legacy_green_tables",
]
