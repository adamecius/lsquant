"""Stage-4 gate program for time-evolution moment computation and curve reconstruction.

Supports two invocation forms:

1. Reconstruct from an existing time-dependent moment file (``.chebmomTD`` or ``.chebmom``):
   ``python -m lsquant.apps.stage4_time --moments <file.chebmomTD> --broadening-meV 10 --fermi-energy <E> --out <dir>``
2. Compute time-dependent moments from CSR operators with :class:`~lsquant._core.MomentEngine`:
   - MSD (G6):
     ``python -m lsquant.apps.stage4_time --operators <dir> --label chain1d --velocity VX --window -2 2 --trace exact --states 256 --moments 256 --steps 64 --tmax 20 --word msd --out <dir>``
   - Projected operator time evolution (G7):
     ``python -m lsquant.apps.stage4_time --operators <dir> --label chain1d_mag --operator SX --projector PX --window <lo> <hi> --trace exact --states 128 --moments 32 --steps 96 --tmax 95 --word projected --out <dir>``
   - Two-operator time correlation (VAC):
     ``python -m lsquant.apps.stage4_time --operators <dir> --label chain1d --left VX --right VX --window -2 2 --trace exact --states 256 --moments 64 --steps 16 --tmax 5 --word correlation --out <dir>``
"""

from __future__ import annotations

import argparse
import os
import sys

from .. import _core
from ..moments import (
    msd_plan,
    projected_evolution_plan,
    run_plan,
    time_correlation_plan,
    time_evolution_series,
    time_step_from_tmax,
    window_from_gershgorin,
)
from ..reconstruct import time_correlation_curve

_WORDS = ("msd", "projected", "correlation", "corr")


def _parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    """Parse command-line arguments for ``stage4_time``."""
    parser = argparse.ArgumentParser(prog="lsquant.apps.stage4_time")
    parser.add_argument(
        "--moments",
        required=True,
        help="Path to a .chebmomTD/.chebmom file, or moment count M when --operators is set",
    )
    parser.add_argument(
        "--broadening-meV",
        type=float,
        default=None,
        dest="broadening_meV",
        help="Jackson Gaussian broadening in meV",
    )
    parser.add_argument(
        "--fermi-energy",
        type=float,
        nargs="+",
        default=None,
        dest="fermi_energy",
        help="Fermi energy (or energies) in eV for time-correlation reconstruction",
    )
    parser.add_argument(
        "--out",
        default=".",
        help="Output directory",
    )
    parser.add_argument(
        "--operators",
        default=None,
        help="Directory containing <label>.HAM.CSR and operator CSR files",
    )
    parser.add_argument(
        "--label",
        default=None,
        help="System label when computing from --operators",
    )
    parser.add_argument(
        "--word",
        choices=_WORDS,
        default=None,
        help="Time-evolution word preset ('msd', 'projected', or 'correlation')",
    )
    parser.add_argument(
        "--velocity",
        default=None,
        help="Velocity / commutator operator name for --word msd (e.g. VX)",
    )
    parser.add_argument(
        "--operator",
        default=None,
        help="Observable operator name for --word projected (e.g. SX)",
    )
    parser.add_argument(
        "--projector",
        default=None,
        help="Projector operator name for --word projected (e.g. PX)",
    )
    parser.add_argument(
        "--left",
        default=None,
        help="Bra-side operator name for --word correlation (e.g. VX)",
    )
    parser.add_argument(
        "--right",
        default=None,
        help="Ket-side operator name for --word correlation (e.g. VX)",
    )
    parser.add_argument(
        "--window",
        nargs="+",
        default=None,
        help="Spectral window: '<lower> <upper>' in eV or 'gershgorin'",
    )
    parser.add_argument(
        "--trace",
        choices=("exact", "random"),
        default="exact",
        help="Trace evaluation mode ('exact' or 'random')",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=12345,
        help="Random-phase PRNG seed (default 12345)",
    )
    parser.add_argument(
        "--states",
        type=int,
        default=None,
        help="Number of trace states",
    )
    parser.add_argument(
        "--steps",
        type=int,
        default=None,
        help="Number of time steps",
    )
    parser.add_argument(
        "--tmax",
        type=float,
        default=None,
        help="Maximum evolution time in fs",
    )
    parser.add_argument(
        "--batch",
        type=int,
        default=None,
        help="Chebyshev batch size (default M; does not apply to the time-axis path)",
    )

    args = parser.parse_args(argv)

    if args.operators is None:
        if args.broadening_meV is None or args.broadening_meV <= 0.0:
            raise ValueError(
                "stage4_time: positive --broadening-meV is required when reconstructing"
            )
        if not args.fermi_energy:
            raise ValueError(
                "stage4_time: --fermi-energy is required when reconstructing"
            )
    else:
        if not args.label:
            raise ValueError("stage4_time: --label is required with --operators")
        if not args.word:
            raise ValueError("stage4_time: --word is required with --operators")
        if args.steps is None or args.steps <= 1:
            raise ValueError("stage4_time: --steps (> 1) is required with --operators")
        if args.tmax is None or args.tmax < 0.0:
            raise ValueError("stage4_time: --tmax is required with --operators")
        try:
            args.num_moments = int(args.moments)
        except ValueError as exc:
            raise ValueError(
                "stage4_time: --moments must be an integer when --operators is provided"
            ) from exc
        if args.num_moments <= 0:
            raise ValueError("stage4_time: --moments must be positive")

        if args.word == "msd":
            if not args.velocity and not args.operator:
                raise ValueError(
                    "stage4_time: --velocity is required with --word msd"
                )
        elif args.word == "projected":
            if not args.operator or not args.projector:
                raise ValueError(
                    "stage4_time: --operator and --projector are required with --word projected"
                )
        elif args.word in ("correlation", "corr"):
            if not args.left or not args.right:
                raise ValueError(
                    "stage4_time: --left and --right are required with --word correlation"
                )

    return args


def _resolve_window(
    window_tokens: list[str] | None, ham: _core.SparseMatrix
) -> _core.SpectralWindow:
    """Resolve ``--window <lo> <hi>`` or ``--window gershgorin``."""
    if window_tokens is None or (
        len(window_tokens) == 1 and window_tokens[0].strip().lower() == "gershgorin"
    ):
        return window_from_gershgorin(ham)
    if len(window_tokens) == 2:
        lo = float(window_tokens[0])
        hi = float(window_tokens[1])
        if not (hi > lo):
            raise ValueError(
                "stage4_time: --window upper bound must exceed lower bound"
            )
        return _core.SpectralWindow(lo, hi)
    raise ValueError(
        "stage4_time: --window must be either '<lo> <hi>' or 'gershgorin'"
    )


def _load_static_op(
    operators_dir: str, label: str, name: str, dim: int
) -> tuple[_core.SparseMatrix | None, _core.LinearOperator]:
    """Load a static CSR operator (or identity when ``name == '1'``)."""
    if name == "1":
        return None, _core.IdentityOperator(dim)
    mat = _core.ReadCsrText(os.path.join(operators_dir, f"{label}.{name}.CSR"))
    return mat, _core.SparseMatrixOperator(mat)


def _reconstruct_one(
    mu: _core.ChebyshevMoments,
    stem: str,
    broadening_meV: float,
    fermi_energy: float,
    out_dir: str,
) -> str:
    """Reconstruct one time-correlation curve at ``fermi_energy`` and write ``.dat``."""
    times, values = time_correlation_curve(mu, broadening_meV, fermi_energy)
    ef_token = f"{float(fermi_energy):.6f}"
    dat_name = f"mean{stem}EF{ef_token}JACKSON.dat"
    dat_path = os.path.join(out_dir, dat_name)
    _core.WriteCurveText(dat_path, times, values)
    return dat_path


def main(argv: list[str] | None = None) -> int:
    """Run ``stage4_time`` from command-line arguments."""
    try:
        args = _parse_args(argv)
        os.makedirs(args.out, exist_ok=True)

        if args.operators is None:
            moments_path = args.moments
            mu = _core.ReadChebmom(moments_path)
            stem = os.path.splitext(os.path.basename(moments_path))[0]
            for ef in args.fermi_energy:
                dat_path = _reconstruct_one(
                    mu, stem, float(args.broadening_meV), float(ef), args.out
                )
                print(dat_path)
            return 0

        ham_path = os.path.join(args.operators, f"{args.label}.HAM.CSR")
        ham = _core.ReadCsrText(ham_path)
        dim = ham.dimension
        window = _resolve_window(args.window, ham)

        if args.trace == "exact":
            num_states = args.states if args.states is not None else dim
            if num_states <= 0 or num_states > dim:
                raise ValueError(
                    "stage4_time: --states must be in 1..dimension for exact trace"
                )
            state_gen: _core.StateGenerator = _core.LocalStates(
                dim, list(range(num_states))
            )
            state_label = "exact"
        else:
            num_states = args.states if args.states is not None else 1
            if num_states <= 0:
                raise ValueError("stage4_time: --states must be positive")
            state_gen = _core.RandomPhaseStates(dim, num_states, int(args.seed))
            state_label = state_gen.label

        num_moments = int(args.num_moments)
        num_steps = int(args.steps)
        tmax = float(args.tmax)
        dt = time_step_from_tmax(num_steps, tmax)
        # Byte identity with the legacy captures requires the GSL Bessel values.
        coeffs, phase = time_evolution_series(window, dt, use_legacy_bessel=True)
        batch_size = int(args.batch) if args.batch is not None else num_moments

        kept_matrices: list[_core.SparseMatrix] = []
        operators: list[_core.LinearOperator] = []

        if args.word == "msd":
            vel_name = args.velocity if args.velocity is not None else args.operator
            if vel_name == "1":
                vel_op: _core.LinearOperator = _core.IdentityOperator(dim)
                prefix = "Correlations1"
            else:
                vel_mat = _core.ReadCsrText(
                    os.path.join(args.operators, f"{args.label}.{vel_name}.CSR")
                )
                kept_matrices.append(vel_mat)
                # Legacy RescaleCon2ChebyshevDomain rescales the commutator with
                # both ScaleFactor and ShiftFactor.
                vel_op = _core.RescaledOperator(vel_mat, window)
                prefix = f"Correlation{vel_name}"
            operators.append(vel_op)
            word_str = f"[{vel_name}, U(t)]^+ delta [{vel_name}, U(t)]"
            plan = msd_plan(
                num_moments,
                num_steps,
                num_states,
                coeffs,
                velocity_index=0,
                phase=phase,
                word=word_str,
            )
            tmax_str = f"{tmax:g}"
            stem = (
                f"{prefix}{args.label}KPM_M{num_moments}_state{tmax_str}{state_label}"
            )
        elif args.word == "projected":
            op_name = args.operator
            proj_name = args.projector
            op_mat, op_lin = _load_static_op(
                args.operators, args.label, op_name, dim
            )
            if op_mat is not None:
                kept_matrices.append(op_mat)
            operators.append(op_lin)

            if proj_name == op_name:
                proj_idx = 0
            else:
                proj_mat, proj_lin = _load_static_op(
                    args.operators, args.label, proj_name, dim
                )
                if proj_mat is not None:
                    kept_matrices.append(proj_mat)
                operators.append(proj_lin)
                proj_idx = 1

            word_str = f"{proj_name} U(t)^+ {op_name} U(t) delta {proj_name}"
            plan = projected_evolution_plan(
                num_moments,
                num_steps,
                num_states,
                coeffs,
                operator_index=0,
                projector_index=proj_idx,
                phase=phase,
                word=word_str,
            )
            stem = (
                f"TimeProj{op_name}-{proj_name}{args.label}KPM_M"
                f"{num_moments}_state{state_label}"
            )
        else:
            left_name = args.left
            right_name = args.right
            left_mat, left_lin = _load_static_op(
                args.operators, args.label, left_name, dim
            )
            if left_mat is not None:
                kept_matrices.append(left_mat)
            operators.append(left_lin)

            if right_name == left_name:
                right_idx = 0
            else:
                right_mat, right_lin = _load_static_op(
                    args.operators, args.label, right_name, dim
                )
                if right_mat is not None:
                    kept_matrices.append(right_mat)
                operators.append(right_lin)
                right_idx = 1

            word_str = f"{left_name} delta U(t)^+ {right_name} U(t)"
            plan = time_correlation_plan(
                num_moments,
                num_steps,
                num_states,
                coeffs,
                left_index=0,
                right_index=right_idx,
                phase=phase,
                word=word_str,
            )
            stem = (
                f"TimeCorr{left_name}-{right_name}{args.label}KPM_M"
                f"{num_moments}_state{state_label}"
            )

        mu = run_plan(
            plan,
            ham,
            window,
            state_gen,
            operators=operators,
            batch_size=batch_size,
            symmetrise=False,
            dt=dt,
        )
        mu.state_label = state_label

        leg_path = os.path.join(args.out, f"{stem}.chebmomTD")
        new_path = os.path.join(args.out, f"{stem}.chebmom")
        _core.WriteLegacyChebmomTD(leg_path, mu)
        _core.WriteChebmom(new_path, mu)
        print(leg_path)
        print(new_path)

        if args.broadening_meV is not None and args.fermi_energy:
            mu_disk = _core.ReadChebmom(leg_path)
            for ef in args.fermi_energy:
                dat_path = _reconstruct_one(
                    mu_disk, stem, float(args.broadening_meV), float(ef), args.out
                )
                print(dat_path)

        return 0
    except Exception as exc:  # noqa: BLE001
        print(f"stage4_time error: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
