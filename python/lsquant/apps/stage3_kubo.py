"""Stage-3 gate program for rank-2 Kubo moment computation and curve reconstruction.

Supports two invocation forms:

1. Reconstruct from an existing rank-2 moment file (``.chebmom2D`` or ``.chebmom``):
   ``python -m lsquant.apps.stage3_kubo --moments <file.chebmom2D> --formula greenwood|bastin|sea|surf|greenwood_fft|bastin_fft --broadening-meV <b> --out <dir>``
2. Compute rank-2 moments from CSR operators with :class:`~lsquant._core.MomentEngine`:
   ``python -m lsquant.apps.stage3_kubo --operators <dir> --label <label> --left VX --right VY --window <lo> <hi> | --window gershgorin --trace exact|random --seed 12345 --states <N> --moments <M> --out <dir>``
"""

from __future__ import annotations

import argparse
import os
import sys
import numpy as np
import numpy.typing as npt

from .. import _core
from ..moments import kubo_plan, run_plan, window_from_gershgorin
from ..reconstruct import (
    kubo_bastin,
    kubo_bastin_sea,
    kubo_bastin_surf,
    kubo_fft_grid,
    kubo_greenwood,
)

_FORMULAS = (
    "greenwood",
    "bastin",
    "sea",
    "surf",
    "greenwood_fft",
    "bastin_fft",
)


def _parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    """Parse command-line options for the two ``stage3_kubo`` forms.

    Parameters
    ----------
    argv : list[str] or None, optional
        Argument list (defaults to ``sys.argv[1:]``).

    Returns
    -------
    argparse.Namespace
        Parsed and validated options namespace.

    Warnings
    --------
    Raises ``ValueError`` if required arguments for the selected form are
    missing or invalid.

    Notes
    -----
    Pitfall: ``--moments`` denotes a file path in the reconstruction form and
    an integer moment count ``M`` when ``--operators`` is given.
    """
    parser = argparse.ArgumentParser(prog="lsquant.apps.stage3_kubo")
    parser.add_argument(
        "--moments",
        required=True,
        help="Path to a .chebmom2D/.chebmom file, or moment count M when --operators is set",
    )
    parser.add_argument(
        "--formula",
        choices=_FORMULAS,
        default=None,
        help="Conductivity reconstruction formula",
    )
    parser.add_argument(
        "--broadening-meV",
        type=float,
        default=None,
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
        help="Directory containing <label>.HAM.CSR and operator CSR files",
    )
    parser.add_argument(
        "--label",
        default=None,
        help="System label when computing from --operators",
    )
    parser.add_argument(
        "--left",
        default=None,
        help="Bra-side operator name (first legacy CLI operator, e.g. VX)",
    )
    parser.add_argument(
        "--right",
        default=None,
        help="Ket-side operator name (second legacy CLI operator, e.g. VY)",
    )
    parser.add_argument(
        "--window",
        nargs="+",
        default=None,
        help="Spectral window: either '<lower> <upper>' in eV or 'gershgorin'",
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
        "--batch",
        type=int,
        default=None,
        help="Chebyshev batch size (default M)",
    )

    args = parser.parse_args(argv)

    if args.operators is None:
        if args.formula is None:
            raise ValueError("stage3_kubo: --formula is required when reconstructing from --moments")
        if args.broadening_meV is None or args.broadening_meV <= 0.0:
            raise ValueError(
                "stage3_kubo: positive --broadening-meV is required when reconstructing"
            )
    else:
        if not args.label:
            raise ValueError("stage3_kubo: --label is required with --operators")
        if not args.left or not args.right:
            raise ValueError("stage3_kubo: --left and --right are required with --operators")
        if not args.window:
            raise ValueError("stage3_kubo: --window is required with --operators")
        try:
            args.num_moments = int(args.moments)
        except ValueError as exc:
            raise ValueError(
                "stage3_kubo: --moments must be an integer when --operators is provided"
            ) from exc
        if args.num_moments <= 0:
            raise ValueError("stage3_kubo: --moments must be positive")
        if args.formula is not None and (
            args.broadening_meV is None or args.broadening_meV <= 0.0
        ):
            raise ValueError(
                "stage3_kubo: positive --broadening-meV is required when --formula is set"
            )

    return args


def _reconstruct_one(
    mu: _core.ChebyshevMoments,
    stem: str,
    formula: str,
    broadening_meV: float,
    out_dir: str,
) -> str:
    """Run one rank-2 reconstruction formula and write the legacy ``.dat`` file."""
    energies: npt.NDArray[np.float64]
    values: npt.NDArray[np.float64]
    if formula == "greenwood":
        energies, values = kubo_greenwood(mu, broadening_meV)
        dat_name = f"KuboGreenWood_{stem}JACKSON.dat"
    elif formula == "bastin":
        energies, values = kubo_bastin(mu, broadening_meV)
        dat_name = f"KuboBastin_{stem}JACKSON.dat"
    elif formula == "sea":
        energies, values = kubo_bastin_sea(mu, broadening_meV)
        dat_name = f"KuboBastinSea_{stem}JACKSON.dat"
    elif formula == "surf":
        energies, values = kubo_bastin_surf(mu, broadening_meV)
        dat_name = f"KuboBastinSurf_{stem}JACKSON.dat"
    elif formula == "greenwood_fft":
        energies, values = kubo_fft_grid(mu, broadening_meV, "greenwood")
        dat_name = f"KuboGreenWood_{stem}JACKSON_FFTgrid.dat"
    elif formula == "bastin_fft":
        energies, values = kubo_fft_grid(mu, broadening_meV, "bastin")
        dat_name = f"KuboBastin_{stem}JACKSON_FFTgrid.dat"
    else:
        raise ValueError(f"stage3_kubo: unsupported formula '{formula}'")

    dat_path = os.path.join(out_dir, dat_name)
    _core.WriteCurveText(dat_path, energies, values)
    return dat_path


def _resolve_window(
    window_tokens: list[str], ham: _core.SparseMatrix
) -> _core.SpectralWindow:
    """Resolve ``--window <lo> <hi>`` or ``--window gershgorin``."""
    if len(window_tokens) == 1 and window_tokens[0].strip().lower() == "gershgorin":
        return window_from_gershgorin(ham)
    if len(window_tokens) == 2:
        lo = float(window_tokens[0])
        hi = float(window_tokens[1])
        if not (hi > lo):
            raise ValueError("stage3_kubo: --window upper bound must exceed lower bound")
        return _core.SpectralWindow(lo, hi)
    raise ValueError(
        "stage3_kubo: --window must be either '<lo> <hi>' or 'gershgorin'"
    )


def main(argv: list[str] | None = None) -> int:
    """Run ``stage3_kubo`` from command-line arguments.

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
    Pitfall: Prints one line to ``stdout`` per written output file.
    """
    try:
        args = _parse_args(argv)
        os.makedirs(args.out, exist_ok=True)

        if args.operators is None:
            moments_path = args.moments
            mu = _core.ReadChebmom(moments_path)
            stem = os.path.splitext(os.path.basename(moments_path))[0]
            dat_path = _reconstruct_one(
                mu, stem, args.formula, float(args.broadening_meV), args.out
            )
            print(dat_path)
        else:
            ham_path = os.path.join(args.operators, f"{args.label}.HAM.CSR")
            ham = _core.ReadCsrText(ham_path)
            dim = ham.dimension
            window = _resolve_window(args.window, ham)

            left_mat: _core.SparseMatrix | None = None
            right_mat: _core.SparseMatrix | None = None
            operators: list[_core.LinearOperator] = []

            if args.left == "1":
                op_left: _core.LinearOperator = _core.IdentityOperator(dim)
            else:
                left_mat = _core.ReadCsrText(
                    os.path.join(args.operators, f"{args.label}.{args.left}.CSR")
                )
                op_left = _core.SparseMatrixOperator(left_mat)
            operators.append(op_left)

            if args.right == args.left:
                left_idx = 0
                right_idx = 0
            else:
                if args.right == "1":
                    op_right: _core.LinearOperator = _core.IdentityOperator(dim)
                else:
                    right_mat = _core.ReadCsrText(
                        os.path.join(args.operators, f"{args.label}.{args.right}.CSR")
                    )
                    op_right = _core.SparseMatrixOperator(right_mat)
                operators.append(op_right)
                left_idx = 0
                right_idx = 1

            if args.trace == "exact":
                num_states = args.states if args.states is not None else dim
                if num_states <= 0 or num_states > dim:
                    raise ValueError(
                        "stage3_kubo: --states must be in 1..dimension for exact trace"
                    )
                state_gen: _core.StateGenerator = _core.LocalStates(
                    dim, list(range(num_states))
                )
                state_label = "exact"
            else:
                num_states = args.states if args.states is not None else 1
                if num_states <= 0:
                    raise ValueError("stage3_kubo: --states must be positive")
                state_gen = _core.RandomPhaseStates(dim, num_states, int(args.seed))
                state_label = state_gen.label

            num_moments = int(args.num_moments)
            batch_size = int(args.batch) if args.batch is not None else num_moments
            word = f"{args.left} delta {args.right} delta"
            plan = kubo_plan(
                num_moments,
                left_index=left_idx,
                right_index=right_idx,
                word=word,
            )
            mu = run_plan(
                plan,
                ham,
                window,
                state_gen,
                operators=operators,
                batch_size=batch_size,
                symmetrise=True,
            )
            mu.state_label = state_label

            stem = (
                f"NonEqOp{args.left}-{args.right}{args.label}KPM_M"
                f"{num_moments}x{num_moments}_state{state_label}"
            )
            leg_path = os.path.join(args.out, f"{stem}.chebmom2D")
            new_path = os.path.join(args.out, f"{stem}.chebmom")
            _core.WriteLegacyChebmom2D(leg_path, mu)
            _core.WriteChebmom(new_path, mu)
            print(leg_path)
            print(new_path)

            if args.formula is not None:
                # Re-read the written .chebmom2D (15-digit legacy ASCII precision)
                # so form 2 with --formula agrees byte-for-byte with form 1.
                mu_disk = _core.ReadChebmom(leg_path)
                dat_path = _reconstruct_one(
                    mu_disk, stem, args.formula, float(args.broadening_meV), args.out
                )
                print(dat_path)

        return 0
    except Exception as exc:  # noqa: BLE001
        print(f"stage3_kubo error: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
