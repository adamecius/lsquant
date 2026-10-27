"""Command-line interface for LSQUANT (``run``, ``reconstruct``, ``inspect``, ``words``)."""

from __future__ import annotations

import argparse
import json
import sys
from typing import Sequence

from .runfile import execute_runfile, inspect_target, reconstruct_target
from .words import describe, list_catalogue


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="lsquant",
        description="LSQUANT unified Kernel Polynomial Method CLI.",
    )
    sub = parser.add_subparsers(dest="command", required=True)

    # 1. run
    p_run = sub.add_parser(
        "run", help="Validate and execute a JSON run file."
    )
    p_run.add_argument("runfile", help="Path to the JSON run file.")
    p_run.add_argument(
        "--out",
        "-o",
        default=None,
        help="Override the output directory specified in the run file.",
    )
    p_run.add_argument(
        "--json",
        action="store_true",
        help="Print execution summary as JSON.",
    )

    # 2. reconstruct
    p_rec = sub.add_parser(
        "reconstruct",
        help="Reconstruct curves from a moment file or run output directory.",
    )
    p_rec.add_argument(
        "target", help="Path to a moment file or run output directory."
    )
    p_rec.add_argument(
        "pos_formula",
        nargs="?",
        default=None,
        help="Optional positional reconstruction formula.",
    )
    p_rec.add_argument(
        "pos_broadening",
        nargs="?",
        type=float,
        default=None,
        help="Optional positional broadening in meV.",
    )
    p_rec.add_argument(
        "--formula",
        "-f",
        default=None,
        help="Reconstruction formula ('shc', 'dos', 'greenwood', 'bastin', 'sea', 'surf', 'greenwood_fft', 'bastin_fft', 'time_correlation').",
    )
    p_rec.add_argument(
        "--broadening-meV",
        "-b",
        type=float,
        default=None,
        dest="broadening_meV",
        help="Gaussian broadening in meV.",
    )
    p_rec.add_argument(
        "--fermi-energies",
        "--fermi-energy",
        "-e",
        nargs="+",
        type=float,
        default=None,
        dest="fermi_energies",
        help="Fermi energy values in eV (for time correlations).",
    )
    p_rec.add_argument(
        "--out",
        "-o",
        default=None,
        help="Destination directory for reconstructed curve files.",
    )
    p_rec.add_argument(
        "--json",
        action="store_true",
        help="Print reconstruction summary as JSON.",
    )

    # 3. inspect
    p_ins = sub.add_parser(
        "inspect",
        help="Inspect a moment file, run file, or output directory and print JSON metadata.",
    )
    p_ins.add_argument(
        "path", help="Path to a moment file, run file, or output directory."
    )

    # 4. words
    p_words = sub.add_parser(
        "words",
        help="List recognised observable word shapes or describe a specific word.",
    )
    p_words.add_argument(
        "action",
        nargs="?",
        default="list",
        help="'list', 'describe', or an observable word string.",
    )
    p_words.add_argument(
        "word",
        nargs="?",
        default=None,
        help="Observable word string when action is 'describe'.",
    )

    return parser


def main(argv: Sequence[str] | None = None) -> int:
    """Entry point for ``python -m lsquant`` and the ``lsquant`` console script."""
    parser = _build_parser()
    try:
        args = parser.parse_args(list(argv) if argv is not None else None)
        cmd = args.command

        if cmd == "run":
            res = execute_runfile(args.runfile, out_dir_override=args.out)
            if args.json:
                print(json.dumps(res, indent=2))
            else:
                for path in res["written_files"]:
                    print(path)
            return 0

        if cmd == "reconstruct":
            formula = (
                args.formula
                if args.formula is not None
                else args.pos_formula
            )
            broadening = (
                args.broadening_meV
                if args.broadening_meV is not None
                else args.pos_broadening
            )
            res = reconstruct_target(
                args.target,
                formula=formula,
                broadening_meV=broadening,
                fermi_energies=args.fermi_energies,
                out_dir=args.out,
            )
            if args.json:
                print(json.dumps(res, indent=2))
            else:
                for path in res["curve_files"]:
                    print(path)
            return 0

        if cmd == "inspect":
            info = inspect_target(args.path)
            print(json.dumps(info, indent=2))
            return 0

        if cmd == "words":
            action = args.action
            if action == "list" and args.word is None:
                print(json.dumps(list_catalogue(), indent=2))
                return 0
            if action == "describe":
                if not args.word:
                    raise ValueError(
                        "lsquant words describe requires a '<word>' argument"
                    )
                print(json.dumps(describe(args.word), indent=2))
                return 0
            # Treat action itself as the word string if not 'list'/'describe'
            print(json.dumps(describe(action), indent=2))
            return 0

        raise ValueError(f"Unknown command '{cmd}'")
    except Exception as exc:  # noqa: BLE001
        print(f"lsquant error: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
