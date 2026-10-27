"""Reproduce Tutorial 04 (Larmor spin precession in a magnetic 1D chain) through the Python API.

Evaluates the three projected time-evolution tensors ``TimeProj{SX,SY,SZ}-PX``
on ``goldens/stage0/chain1d_mag/operators`` (``N = 128``, ``M = 32``, ``96`` time
steps, ``tmax = 95`` fs, exact trace, Gershgorin window ``[-2.31, 2.31]``),
writes the legacy ``.chebmomTD`` files into ``<out>`` (and ``<out>/captured/``),
and extracts the ``m = 0`` row ``S_a(t) = Re mu_0(t)``.
"""

from __future__ import annotations

import argparse
import os
import shutil
import sys
import numpy as np
import numpy.typing as npt

from .. import _core
from ..moments import (
    projected_evolution_plan,
    run_plan,
    time_evolution_series,
    time_step_from_tmax,
    window_from_gershgorin,
)

LABEL: str = "chain1d_mag"
SPIN_OPS: tuple[str, ...] = ("SX", "SY", "SZ")
PROJECTOR: str = "PX"
M_MOMENTS: int = 32
NTIME: int = 96
TMAX: float = 95.0


def run_tutorial04(
    out_dir: str,
    operators_dir: str | None = None,
) -> tuple[list[str], dict[str, npt.NDArray[np.float64]]]:
    """Run Tutorial 04 and write the three ``TimeProj*.chebmomTD`` files into ``out_dir``."""
    repo_root = os.path.abspath(
        os.path.join(os.path.dirname(__file__), "..", "..", "..")
    )
    resolved_op_dir = (
        os.path.abspath(operators_dir)
        if operators_dir is not None
        else os.path.join(repo_root, "goldens", "stage0", "chain1d_mag", "operators")
    )
    os.makedirs(out_dir, exist_ok=True)
    captured_dir = os.path.join(out_dir, "captured")
    os.makedirs(captured_dir, exist_ok=True)

    ham = _core.ReadCsrText(os.path.join(resolved_op_dir, f"{LABEL}.HAM.CSR"))
    dim = ham.dimension
    window = window_from_gershgorin(ham)
    dt = time_step_from_tmax(NTIME, TMAX)
    # Byte identity with the legacy captures requires the GSL Bessel values.
    coeffs, phase = time_evolution_series(window, dt, use_legacy_bessel=True)

    proj_mat = _core.ReadCsrText(
        os.path.join(resolved_op_dir, f"{LABEL}.{PROJECTOR}.CSR")
    )
    proj_op = _core.SparseMatrixOperator(proj_mat)
    states = _core.LocalStates(dim, list(range(dim)))

    written: list[str] = []
    m0_rows: dict[str, npt.NDArray[np.float64]] = {}

    for idx, op_name in enumerate(SPIN_OPS):
        op_mat = _core.ReadCsrText(
            os.path.join(resolved_op_dir, f"{LABEL}.{op_name}.CSR")
        )
        op_lin = _core.SparseMatrixOperator(op_mat)
        word = f"{PROJECTOR} U(t)^+ {op_name} U(t) delta {PROJECTOR}"
        plan = projected_evolution_plan(
            M_MOMENTS,
            NTIME,
            dim,
            coeffs,
            operator_index=0,
            projector_index=1,
            phase=phase,
            word=word,
        )
        mu = run_plan(
            plan,
            ham,
            window,
            states,
            operators=[op_lin, proj_op],
            batch_size=M_MOMENTS,
            symmetrise=False,
            dt=dt,
        )
        mu.state_label = "exact"

        fname = f"TimeProj{op_name}-{PROJECTOR}{LABEL}KPM_M{M_MOMENTS}_stateexact.chebmomTD"
        out_path = os.path.join(out_dir, fname)
        _core.WriteLegacyChebmomTD(out_path, mu)
        written.append(out_path)

        cap_path = os.path.join(captured_dir, f"{idx:04d}_{fname}")
        shutil.copyfile(out_path, cap_path)

        mu_disk = _core.ReadChebmom(out_path)
        m0_rows[op_name] = np.array(mu_disk.data[0, :].real, dtype=np.float64)

    return written, m0_rows


def main(argv: list[str] | None = None) -> int:
    """Command-line entry point for ``lsquant.apps.tutorial04``."""
    try:
        parser = argparse.ArgumentParser(prog="lsquant.apps.tutorial04")
        parser.add_argument("pos_out", nargs="?", default=None, help="Output directory")
        parser.add_argument(
            "--out", "-o", dest="out", default=None, help="Output directory"
        )
        parser.add_argument(
            "--operators",
            default=None,
            help="Directory containing chain1d_mag.{HAM,SX,SY,SZ,PX}.CSR (default: goldens/stage0/chain1d_mag/operators)",
        )
        args = parser.parse_args(argv)
        out_dir = (
            args.out
            if args.out is not None
            else (args.pos_out if args.pos_out is not None else ".")
        )
        paths, _ = run_tutorial04(out_dir, operators_dir=args.operators)
        for p in paths:
            print(p)
        return 0
    except Exception as exc:  # noqa: BLE001
        print(f"tutorial04 error: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
