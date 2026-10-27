"""Stage-2 gate program for rank-1 density-of-states computation and reconstruction.

Supports two invocation forms:

1. Reconstruct from an existing moment file (``.chebmom1D`` or ``.chebmom``):
   ``python -m lsquant.apps.stage2_dos --moments <file> --broadening-meV 120 --out <dir>``
2. Compute moments from CSR operators with :class:`~lsquant._core.MomentEngine`
   and reconstruct:
   ``python -m lsquant.apps.stage2_dos --operators <dir> --label chain1d --window -2 2 --trace exact --states 2048 --moments 64 --broadening-meV 120 --out <dir>``
"""

from __future__ import annotations

import argparse
import os
import sys

from .. import _core
from ..moments import dos_plan, run_plan
from ..reconstruct import dos_curve


def _parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    """Parse command-line options for the two ``stage2_dos`` forms.

    Parameters
    ----------
    argv : list[str] or None, optional
        Argument list (defaults to ``sys.argv[1:]``).

    Returns
    -------
    argparse.Namespace
        Parsed options namespace.

    Warnings
    --------
    Raises ``ValueError`` if required arguments for the selected form are
    missing or invalid.

    Notes
    -----
    Pitfall: ``--moments`` denotes a file path in the reconstruction-only form
    and an integer moment count ``M`` when ``--operators`` is given.
    """
    parser = argparse.ArgumentParser(prog="lsquant.apps.stage2_dos")
    parser.add_argument(
        "--moments",
        required=True,
        help="Path to a .chebmom1D/.chebmom file, or moment count M when --operators is set",
    )
    parser.add_argument(
        "--broadening-meV",
        type=float,
        required=True,
        dest="broadening_meV",
        help="Jackson Gaussian broadening in meV",
    )
    parser.add_argument(
        "--out",
        default=".",
        help="Output directory",
    )
    parser.add_argument(
        "--operators",
        default=None,
        help="Directory containing <label>.HAM.CSR",
    )
    parser.add_argument(
        "--label",
        default=None,
        help="System label when computing from --operators",
    )
    parser.add_argument(
        "--window",
        type=float,
        nargs=2,
        default=None,
        metavar=("LOWER", "UPPER"),
        help="Spectral window lower and upper bounds in eV",
    )
    parser.add_argument(
        "--trace",
        choices=("exact", "random"),
        default="exact",
        help="Trace evaluation mode ('exact' or 'random')",
    )
    parser.add_argument(
        "--states",
        type=int,
        default=None,
        help="Number of trace states",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=12345,
        help="Random-phase PRNG seed (default 12345)",
    )
    parser.add_argument(
        "--batch",
        type=int,
        default=None,
        help="Chebyshev batch size (default M)",
    )

    args = parser.parse_args(argv)
    if args.broadening_meV <= 0.0:
        raise ValueError("stage2_dos: --broadening-meV must be positive")

    if args.operators is not None:
        if not args.label:
            raise ValueError("stage2_dos: --label is required with --operators")
        if args.window is None or not (args.window[1] > args.window[0]):
            raise ValueError(
                "stage2_dos: valid --window <lower> <upper> is required with --operators"
            )
        try:
            args.num_moments = int(args.moments)
        except ValueError as exc:
            raise ValueError(
                "stage2_dos: --moments must be an integer when --operators is provided"
            ) from exc
        if args.num_moments <= 0:
            raise ValueError("stage2_dos: --moments must be positive")

    return args


def main(argv: list[str] | None = None) -> int:
    """Run ``stage2_dos`` from command-line arguments.

    Parameters
    ----------
    argv : list[str] or None, optional
        Command-line arguments.

    Returns
    -------
    int
        Process exit status (``0`` on success, ``1`` on error).

    Warnings
    --------
    Prints any uncaught exception message to ``sys.stderr`` and returns ``1``.

    Notes
    -----
    Pitfall: Prints exactly one line to ``stdout`` per written output file.
    """
    try:
        args = _parse_args(argv)
        os.makedirs(args.out, exist_ok=True)

        if args.operators is None:
            moments_path = args.moments
            mu = _core.ReadChebmom(moments_path)
            stem = os.path.splitext(os.path.basename(moments_path))[0]
            energies, values = dos_curve(mu, args.broadening_meV)
            dat_path = os.path.join(args.out, f"mean{stem}JACKSON.dat")
            _core.WriteCurveText(dat_path, energies, values)
            print(dat_path)
        else:
            ham_path = os.path.join(args.operators, f"{args.label}.HAM.CSR")
            ham = _core.ReadCsrText(ham_path)
            dim = ham.dimension
            window = _core.SpectralWindow(float(args.window[0]), float(args.window[1]))

            if args.trace == "exact":
                num_states = args.states if args.states is not None else dim
                if num_states <= 0 or num_states > dim:
                    raise ValueError(
                        "stage2_dos: --states must be in 1..dimension for exact trace"
                    )
                state_gen: _core.StateGenerator = _core.LocalStates(
                    dim, list(range(num_states))
                )
                state_label = "exact"
            else:
                num_states = args.states if args.states is not None else 1
                if num_states <= 0:
                    raise ValueError("stage2_dos: --states must be positive")
                state_gen = _core.RandomPhaseStates(dim, num_states, int(args.seed))
                state_label = state_gen.label

            num_moments = int(args.num_moments)
            batch_size = int(args.batch) if args.batch is not None else num_moments
            plan = dos_plan(num_moments, num_states, word="1 delta")
            mu = run_plan(plan, ham, window, state_gen, batch_size=batch_size)
            mu.state_label = state_label

            stem = f"SpectralOp1{args.label}KPM_M{num_moments}_state{state_label}"
            mom_path = os.path.join(args.out, f"{stem}.chebmom1D")
            _core.WriteLegacyChebmom1D(mom_path, mu)
            print(mom_path)

            energies, values = dos_curve(mu, args.broadening_meV)
            dat_path = os.path.join(args.out, f"mean{stem}JACKSON.dat")
            _core.WriteCurveText(dat_path, energies, values)
            print(dat_path)

        return 0
    except Exception as exc:  # noqa: BLE001
        print(f"stage2_dos error: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
