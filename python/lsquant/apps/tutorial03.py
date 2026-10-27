"""Reproduce Tutorial 03 (Haldane anomalous Hall conductivity) through the Python API.

Stages the committed Haldane operators by importing
``examples/03_haldane_hall/make_haldane.py``, computes the ``VX-VY`` and
``VX-VX`` rank-2 Chebyshev moments (``M = 128``, exact trace) through
:class:`~lsquant._core.MomentEngine`, reconstructs the Kubo-Bastin, Sea, Surf,
and Kubo-Greenwood curves at ``10`` meV, writes the four ``.dat`` files with
their legacy names into ``<out>``, and prints the Chern number extracted from
the gap plateau (``|E| < 0.5`` eV).
"""

from __future__ import annotations

import argparse
import importlib.util
import math
import os
import shutil
import sys
import tempfile
from types import ModuleType
import numpy as np
import numpy.typing as npt

from .. import _core
from ..moments import kubo_plan, run_plan
from ..reconstruct import (
    kubo_bastin,
    kubo_bastin_sea,
    kubo_bastin_surf,
    kubo_greenwood,
)

A_CELL: float = 3.0 * math.sqrt(3.0) / 2.0


def _load_make_haldane() -> tuple[ModuleType, str]:
    """Import ``examples/03_haldane_hall/make_haldane.py`` from the repository tree.

    Returns
    -------
    tuple[types.ModuleType, str]
        Pair ``(module, repo_root)``.

    Warnings
    --------
    Raises ``FileNotFoundError`` if ``examples/03_haldane_hall/make_haldane.py``
    cannot be located.

    Notes
    -----
    Pitfall: ``make_haldane`` is imported only for its descriptor writer
    (``_write_desc``) and exact-state writer (``build``); its ``GOLDEN`` path is
    overridden by the caller with ``goldens/stage0/haldane/operators`` (or
    ``--operators``).
    """
    repo_root = os.path.abspath(
        os.path.join(os.path.dirname(__file__), "..", "..", "..")
    )
    script_path = os.path.join(repo_root, "examples", "03_haldane_hall", "make_haldane.py")
    if not os.path.isfile(script_path):
        raise FileNotFoundError(f"tutorial03: make_haldane.py not found at {script_path}")

    spec = importlib.util.spec_from_file_location("lsquant_tutorial03_make_haldane", script_path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"tutorial03: could not load spec for {script_path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module, repo_root


def _read_desc_window(desc_path: str) -> _core.SpectralWindow:
    """Parse ``spectral_min`` and ``spectral_max`` from an operator ``.desc`` sidecar."""
    s_min: float | None = None
    s_max: float | None = None
    with open(desc_path, encoding="utf-8") as fh:
        for raw_line in fh:
            line = raw_line.strip()
            if not line or line.startswith("#") or ":" not in line:
                continue
            key, val = line.split(":", 1)
            key = key.strip()
            val = val.strip()
            if key == "spectral_min":
                s_min = float(val)
            elif key == "spectral_max":
                s_max = float(val)
    if s_min is None or s_max is None or not (s_max > s_min):
        raise ValueError(f"tutorial03: missing valid spectral bounds in {desc_path}")
    return _core.SpectralWindow(s_min, s_max)


def chern_from_plateau(
    dim: int,
    energies: npt.NDArray[np.float64],
    sigma_xy: npt.NDArray[np.float64],
    gap_halfwidth: float = 0.5,
) -> float:
    """Compute the Chern number from the gap plateau ``|E| < gap_halfwidth``.

    Solves the problem of reading the integer topological invariant off the
    extensive Kubo-Bastin Hall conductivity ``sigma_xy(E)`` in units ``e = hbar = 1``
    using the exact formula of ``examples/03_haldane_hall/lsqhall.py``:

    .. math::

        C = \\frac{\\langle \\sigma_{xy} \\rangle_{|E| < 0.5}}{A} \\, 2\\pi, \\qquad
        A = \\frac{\\text{dim}}{2} \\frac{3\\sqrt{3}}{2}.

    Parameters
    ----------
    dim : int
        Hilbert-space dimension ``N`` (two orbitals per honeycomb unit cell).
    energies : numpy.ndarray
        1D energy grid in eV.
    sigma_xy : numpy.ndarray
        1D extensive Hall conductivity curve matching ``energies``.
    gap_halfwidth : float, optional
        Half-width of the plateau averaging window around ``E = 0`` in eV
        (default ``0.5``).

    Returns
    -------
    float
        Estimated Chern number ``C``.

    Warnings
    --------
    Raises ``ValueError`` if no energy points fall inside ``|E| < gap_halfwidth``.

    Notes
    -----
    Pitfall: ``sigma_xy`` returned by :func:`~lsquant.reconstruct.kubo_bastin`
    is extensive (proportional to total system area ``A``).
    """
    mask = np.abs(energies) < gap_halfwidth
    if not np.any(mask):
        raise ValueError("chern_from_plateau: no energy points inside gap window")
    area = (float(dim) / 2.0) * A_CELL
    plateau = float(np.mean(sigma_xy[mask]))
    return (plateau / area) * 2.0 * np.pi


def run_tutorial03(
    out_dir: str,
    operators_dir: str | None = None,
    num_moments: int = 128,
    broadening_meV: float = 10.0,
) -> tuple[list[str], float]:
    """Stage Haldane operators, compute moments, write the four curves, and return ``C``.

    Parameters
    ----------
    out_dir : str
        Destination directory for the four reconstructed ``.dat`` files.
    operators_dir : str or None, optional
        Directory containing ``haldane.{HAM,VX,VY}.CSR`` (defaults to
        ``goldens/stage0/haldane/operators`` relative to the repository root).
    num_moments : int, optional
        Chebyshev expansion order along both axes (default ``128``).
    broadening_meV : float, optional
        Jackson broadening in meV (default ``10.0``).

    Returns
    -------
    tuple[list[str], float]
        Pair ``(written_paths, chern_number)``.

    Warnings
    --------
    Overwrites any existing output ``.dat`` files in ``out_dir``.

    Notes
    -----
    Pitfall: Evaluates ``chern_from_plateau`` on the formatted ``KuboBastin``
    ``.dat`` file read back via ``numpy.loadtxt`` so the reported ``C`` matches
    ``lsqhall.py`` to the last digit.
    """
    make_haldane, repo_root = _load_make_haldane()
    resolved_op_dir = (
        os.path.abspath(operators_dir)
        if operators_dir is not None
        else os.path.join(repo_root, "goldens", "stage0", "haldane", "operators")
    )
    make_haldane.GOLDEN = resolved_op_dir
    os.makedirs(out_dir, exist_ok=True)
    label = "haldane"

    with tempfile.TemporaryDirectory() as stage_dir:
        dim = int(make_haldane.build(stage_dir))
        staged_op_dir = os.path.join(stage_dir, "operators")
        window = _read_desc_window(os.path.join(staged_op_dir, f"{label}.HAM.desc"))
        with open(os.path.join(stage_dir, "exact"), encoding="utf-8") as fh:
            state_lines = [ln.strip() for ln in fh if ln.strip()]
        num_local = int(state_lines[1])
        local_indices = [int(x) for x in state_lines[2 : 2 + num_local]]

    ham = _core.ReadCsrText(os.path.join(resolved_op_dir, f"{label}.HAM.CSR"))
    vx = _core.ReadCsrText(os.path.join(resolved_op_dir, f"{label}.VX.CSR"))
    vy = _core.ReadCsrText(os.path.join(resolved_op_dir, f"{label}.VY.CSR"))

    vx_op = _core.SparseMatrixOperator(vx)
    vy_op = _core.SparseMatrixOperator(vy)
    states = _core.LocalStates(dim, local_indices)

    # Hall response: VX (bra) - VY (ket)
    plan_xy = kubo_plan(
        num_moments, left_index=0, right_index=1, word="VX delta VY delta"
    )
    mu_xy = run_plan(
        plan_xy,
        ham,
        window,
        states,
        operators=[vx_op, vy_op],
        batch_size=num_moments,
        symmetrise=True,
    )
    mu_xy.state_label = "exact"
    stem_xy = f"NonEqOpVX-VY{label}KPM_M{num_moments}x{num_moments}_stateexact"

    # Longitudinal response: VX (bra) - VX (ket)
    plan_xx = kubo_plan(
        num_moments, left_index=0, right_index=0, word="VX delta VX delta"
    )
    mu_xx = run_plan(
        plan_xx,
        ham,
        window,
        states,
        operators=[vx_op],
        batch_size=num_moments,
        symmetrise=True,
    )
    mu_xx.state_label = "exact"
    stem_xx = f"NonEqOpVX-VX{label}KPM_M{num_moments}x{num_moments}_stateexact"

    # Legacy pipeline artefact: legacy lsqhall.py wrote the .chebmom2D file at
    # 15-digit ASCII precision (digits10) and reconstructed from the re-read
    # file; the text round-trip drops the last two digits of the in-memory
    # double-precision moments. The library itself keeps full precision; the
    # tutorial app round-trips through WriteLegacyChebmom2D / ReadChebmom so
    # its four output curves are byte-identical to the stage-0 capture.
    xy_mom_path = os.path.join(out_dir, f"{stem_xy}.chebmom2D")
    xx_mom_path = os.path.join(out_dir, f"{stem_xx}.chebmom2D")
    _core.WriteLegacyChebmom2D(xy_mom_path, mu_xy)
    _core.WriteLegacyChebmom2D(xx_mom_path, mu_xx)
    mu_xy_disk = _core.ReadChebmom(xy_mom_path)
    mu_xx_disk = _core.ReadChebmom(xx_mom_path)

    written: list[str] = []
    captured_dir = os.path.join(out_dir, "captured")
    os.makedirs(captured_dir, exist_ok=True)

    e_bastin, s_bastin = kubo_bastin(mu_xy_disk, broadening_meV)
    bastin_path = os.path.join(out_dir, f"KuboBastin_{stem_xy}JACKSON.dat")
    _core.WriteCurveText(bastin_path, e_bastin, s_bastin)
    written.append(bastin_path)

    e_sea, s_sea = kubo_bastin_sea(mu_xy_disk, broadening_meV)
    sea_path = os.path.join(out_dir, f"KuboBastinSea_{stem_xy}JACKSON.dat")
    _core.WriteCurveText(sea_path, e_sea, s_sea)
    written.append(sea_path)

    e_surf, s_surf = kubo_bastin_surf(mu_xy_disk, broadening_meV)
    surf_path = os.path.join(out_dir, f"KuboBastinSurf_{stem_xy}JACKSON.dat")
    _core.WriteCurveText(surf_path, e_surf, s_surf)
    written.append(surf_path)

    e_gw, s_gw = kubo_greenwood(mu_xx_disk, broadening_meV)
    gw_path = os.path.join(out_dir, f"KuboGreenWood_{stem_xx}JACKSON.dat")
    _core.WriteCurveText(gw_path, e_gw, s_gw)
    written.append(gw_path)

    for cap_idx, src_path in enumerate((bastin_path, sea_path, surf_path, gw_path)):
        shutil.copyfile(
            src_path,
            os.path.join(captured_dir, f"{cap_idx:04d}_{os.path.basename(src_path)}"),
        )

    loaded_bastin = np.loadtxt(bastin_path)
    c_val = chern_from_plateau(dim, loaded_bastin[:, 0], loaded_bastin[:, 1])
    return written, c_val


def main(argv: list[str] | None = None) -> int:
    """Command-line entry point for ``lsquant.apps.tutorial03``.

    Parameters
    ----------
    argv : list[str] or None, optional
        Command-line arguments.

    Returns
    -------
    int
        Exit status (``0`` on success, ``1`` on error).

    Warnings
    --------
    Prints any uncaught exception to ``stderr`` and returns ``1``.

    Notes
    -----
    Pitfall: Prints each written ``.dat`` path followed by the plateau Chern
    number line.
    """
    try:
        parser = argparse.ArgumentParser(prog="lsquant.apps.tutorial03")
        parser.add_argument("pos_out", nargs="?", default=None, help="Output directory")
        parser.add_argument("--out", "-o", dest="out", default=None, help="Output directory")
        parser.add_argument(
            "--operators",
            default=None,
            help="Directory containing haldane.{HAM,VX,VY}.CSR (default: goldens/stage0/haldane/operators)",
        )
        args = parser.parse_args(argv)
        out_dir = (
            args.out
            if args.out is not None
            else (args.pos_out if args.pos_out is not None else ".")
        )
        paths, c_val = run_tutorial03(out_dir, operators_dir=args.operators)
        for p in paths:
            print(p)
        print(f"Chern from gap plateau: C = {c_val:+.3f} ({c_val:+.12f})")
        return 0
    except Exception as exc:  # noqa: BLE001
        print(f"tutorial03 error: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
