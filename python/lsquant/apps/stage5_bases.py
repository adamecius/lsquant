"""Stage-5 gate program evaluating DOS and Kubo words across Real, Bloch, and Eigen bases.

Invocation contract:

.. code-block:: bash

   python -m lsquant.apps.stage5_bases \\
     --model <dir with <seed>_hr.dat, <seed>.uc, <seed>.xyz> \\
     --seed graphene \\
     --supercell 48 48 1 \\
     --basis real|bloch|eigen \\
     --word dos|kuboxx \\
     --moments 256 \\
     --trace exact|random --states 32 --seed 12345 \\
     [--disorder W --disorder-seed 7] \\
     [--window eigen|gershgorin|<lo> <hi>] \\
     --broadening-meV 20 \\
     --out <dir>

The default spectral window is ``--window eigen`` for every basis: the app
computes the eigenvalue bounds once from the ``H(k)`` blocks (widened by ``W/2``
when ``--disorder W > 0`` and padded by ``0.10`` of the spectral width) so that
``real``, ``bloch``, and ``eigen`` runs of the same comparison share the
identical window. Gershgorin is not used by default on Wannier models because
spread-out hoppings make the Gershgorin enclosure several times wider than the
true spectrum and waste Chebyshev moments.

In all three bases, ``stage5_bases`` runs the identical :class:`~lsquant._core.Plan`
(:func:`~lsquant.moments.dos_plan` or :func:`~lsquant.moments.kubo_plan`) through
:class:`~lsquant._core.MomentEngine` and reconstructs with
:func:`~lsquant.reconstruct.dos_curve` or
:func:`~lsquant.reconstruct.kubo_greenwood`, writing the legacy-layout moment
file (``.chebmom1D`` or ``.chebmom2D``), the new ``.chebmom`` file, the curve
with the legacy filename (``mean<stem>JACKSON.dat`` or
``KuboGreenWood_<stem>JACKSON.dat``), and ``<stem>.basis.json`` to ``--out``.
Output filenames use ``<seed>_<basis>`` as the label part of the legacy stem
(e.g. ``SpectralOp1graphene_bloch_KPM_M256_stateexact.chebmom1D``) so that the
three bases written to one directory do not overwrite each other; the legacy
stem is kept otherwise.
"""

from __future__ import annotations

import argparse
import json
import os
import sys
from typing import Sequence

from .. import _core
from ..bases import BlochBasis, EigenBasis, RealBasis, make_basis
from ..models.wannier import WannierModel, load_wannier_model
from ..moments import dos_plan, kubo_plan, run_plan
from ..reconstruct import dos_curve, kubo_greenwood


def _resolve_seeds(
    model_dir: str,
    seed_tokens: list[str] | None,
    explicit_label: str | None,
    explicit_rng_seed: int | None,
) -> tuple[str | None, int]:
    """Disambiguate ``--seed <model_name>`` and ``--seed <int_rng_seed>`` when both are passed."""
    model_seed: str | None = explicit_label
    rng_seed: int | None = explicit_rng_seed

    for tok in seed_tokens or []:
        try:
            val = int(tok)
            if not os.path.isfile(os.path.join(model_dir, f"{tok}_hr.dat")):
                rng_seed = val
                continue
        except ValueError:
            pass
        model_seed = tok

    if rng_seed is None:
        rng_seed = 12345
    return model_seed, int(rng_seed)


def _parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    """Parse and validate command-line options for ``stage5_bases``."""
    parser = argparse.ArgumentParser(prog="lsquant.apps.stage5_bases")
    parser.add_argument(
        "--model",
        "--model-dir",
        required=True,
        dest="model_dir",
        help="Directory containing <seed>_hr.dat, <seed>.uc, <seed>.xyz",
    )
    parser.add_argument(
        "--seed",
        action="append",
        default=[],
        help="Model seed prefix (e.g. 'graphene') and/or PRNG integer seed (e.g. 12345)",
    )
    parser.add_argument(
        "--label",
        default=None,
        help="Explicit model seed prefix if --seed is used only for the PRNG seed",
    )
    parser.add_argument(
        "--rng-seed",
        type=int,
        default=None,
        dest="rng_seed",
        help="Explicit integer state PRNG seed",
    )
    parser.add_argument(
        "--basis",
        required=True,
        choices=("real", "bloch", "eigen"),
        help="Basis adapter ('real', 'bloch', or 'eigen')",
    )
    parser.add_argument(
        "--supercell",
        nargs=3,
        type=int,
        required=True,
        metavar=("N1", "N2", "N3"),
        help="Supercell dimensions N1 N2 N3",
    )
    parser.add_argument(
        "--disorder",
        "--disorder-w",
        type=float,
        default=0.0,
        dest="disorder_w",
        help="Real-space Anderson box disorder strength W in eV (default 0.0)",
    )
    parser.add_argument(
        "--disorder-seed",
        type=int,
        default=7,
        dest="disorder_seed",
        help="PRNG seed for real-space Anderson disorder (default 7, independent of state seed)",
    )
    parser.add_argument(
        "--csr-dir",
        default=None,
        dest="csr_dir",
        help="Optional directory containing wannier2sparse <seed>.HAM.CSR and <seed>.VX.CSR",
    )
    parser.add_argument(
        "--word",
        required=True,
        choices=("dos", "kuboxx", "kuboxy"),
        help="Observable word ('dos', 'kuboxx', or 'kuboxy')",
    )
    parser.add_argument(
        "--moments",
        type=int,
        required=True,
        help="Number of Chebyshev moments M",
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
        help="Number of random-phase trace vectors (when --trace random)",
    )
    parser.add_argument(
        "--window",
        nargs="+",
        default=["eigen"],
        help="Spectral window: 'eigen' (default), 'gershgorin', or '<lower> <upper>' in eV",
    )
    parser.add_argument(
        "--broadening-meV",
        type=float,
        required=True,
        dest="broadening_meV",
        help="Jackson Gaussian broadening in meV",
    )
    parser.add_argument(
        "--batch",
        type=int,
        default=None,
        help="Chebyshev batch size (default M)",
    )
    parser.add_argument(
        "--out",
        default=".",
        help="Output directory",
    )

    args = parser.parse_args(argv)
    if args.moments <= 0:
        raise ValueError("stage5_bases: --moments must be positive")
    if args.broadening_meV <= 0.0:
        raise ValueError("stage5_bases: --broadening-meV must be positive")
    if args.disorder_w < 0.0:
        raise ValueError("stage5_bases: --disorder must be non-negative")
    return args


def _resolve_window(
    window_tokens: list[str],
    model: WannierModel,
    supercell: Sequence[int],
    basis_obj: RealBasis | BlochBasis | EigenBasis,
) -> tuple[_core.SpectralWindow, str]:
    """Resolve ``--window eigen``, ``--window gershgorin``, or ``--window <lo> <hi>``.

    Defaults to ``eigen`` across all three bases so that ``real``, ``bloch``,
    and ``eigen`` runs of one comparison share the exact same window computed
    from the ``H(k)`` block eigenvalues (widened by ``W/2`` before the 10% pad
    when ``disorder_w > 0``). Gershgorin is not used by default on Wannier
    models because spread-out hoppings make the Gershgorin enclosure several
    times wider than the spectrum.
    """
    if len(window_tokens) == 1:
        key = window_tokens[0].strip().lower()
        if key == "eigen":
            if isinstance(basis_obj, (BlochBasis, EigenBasis)):
                return basis_obj.window(), "eigen"
            clean_win = EigenBasis(model, supercell).window(pad=0.0)
            emin = float(clean_win.lower) - 0.5 * float(basis_obj.disorder_w)
            emax = float(clean_win.upper) + 0.5 * float(basis_obj.disorder_w)
            width = emax - emin
            if width <= 0.0:
                emin -= 1.0
                emax += 1.0
                width = 2.0
            margin = 0.10 * width
            return _core.SpectralWindow(emin - margin, emax + margin), "eigen"
        if key == "gershgorin":
            if isinstance(basis_obj, RealBasis):
                return basis_obj.window(), "gershgorin"
            return (
                RealBasis(
                    model,
                    supercell,
                    disorder_w=basis_obj.disorder_w,
                    seed=basis_obj.seed,
                    disorder_seed=getattr(basis_obj, "disorder_seed", basis_obj.seed),
                ).window(),
                "gershgorin",
            )
    if len(window_tokens) == 2:
        lo = float(window_tokens[0])
        hi = float(window_tokens[1])
        if not (hi > lo):
            raise ValueError("stage5_bases: --window upper bound must exceed lower bound")
        return _core.SpectralWindow(lo, hi), "explicit"
    raise ValueError(
        "stage5_bases: --window must be 'eigen', 'gershgorin', or '<lo> <hi>'"
    )


def main(argv: list[str] | None = None) -> int:
    """Run ``stage5_bases`` from command-line arguments."""
    try:
        args = _parse_args(argv)
        os.makedirs(args.out, exist_ok=True)

        model_seed, rng_seed = _resolve_seeds(
            args.model_dir, args.seed, args.label, args.rng_seed
        )
        model = load_wannier_model(args.model_dir, seed=model_seed)
        basis = make_basis(
            args.basis,
            model,
            args.supercell,
            disorder_w=float(args.disorder_w),
            seed=rng_seed,
            disorder_seed=int(args.disorder_seed),
            csr_dir=args.csr_dir,
        )

        window, window_source = _resolve_window(
            args.window, model, args.supercell, basis
        )
        state_gen = basis.states(
            trace_mode=args.trace,
            num_states=args.states,
            seed=rng_seed,
        )

        num_moments = int(args.moments)
        batch_size = int(args.batch) if args.batch is not None else num_moments
        word_key = args.word.strip().lower()

        if word_key == "dos":
            used_ops = ["HAM"]
            plan = dos_plan(
                num_moments,
                int(state_gen.number_of_states),
                word="1 delta",
            )
            mu = run_plan(
                plan,
                basis.hamiltonian(),
                window,
                state_gen,
                operators=[],
                batch_size=batch_size,
                symmetrise=False,
            )
            if args.trace == "exact":
                mu.number_of_states = basis.dimension
            energies, values = dos_curve(mu, float(args.broadening_meV))
            stem = (
                f"SpectralOp1{model.seed_name}_{args.basis}_KPM_M{num_moments}"
                f"_state{mu.state_label}"
            )
            legacy_cheb_path = os.path.join(args.out, f"{stem}.chebmom1D")
            cheb_path = os.path.join(args.out, f"{stem}.chebmom")
            dat_path = os.path.join(args.out, f"mean{stem}JACKSON.dat")
            _core.WriteLegacyChebmom1D(legacy_cheb_path, mu)
        elif word_key in ("kuboxx", "kuboxy"):
            op_vx = basis.operator("VX")
            if word_key == "kuboxx":
                used_ops = ["HAM", "VX"]
                operators = [op_vx]
                left_idx, right_idx = 0, 0
                word_str = "VX delta VX delta"
                op_tag = "VX-VX"
            else:
                used_ops = ["HAM", "VX", "VY"]
                op_vy = basis.operator("VY")
                operators = [op_vx, op_vy]
                left_idx, right_idx = 0, 1
                word_str = "VX delta VY delta"
                op_tag = "VX-VY"
            plan = kubo_plan(
                num_moments,
                left_index=left_idx,
                right_index=right_idx,
                word=word_str,
            )
            mu = run_plan(
                plan,
                basis.hamiltonian(),
                window,
                state_gen,
                operators=operators,
                batch_size=batch_size,
                symmetrise=(word_key == "kuboxx"),
            )
            if args.trace == "exact":
                mu.number_of_states = basis.dimension
            energies, values = kubo_greenwood(mu, float(args.broadening_meV))
            stem = (
                f"NonEqOp{op_tag}{model.seed_name}_{args.basis}_KPM_M"
                f"{num_moments}x{num_moments}_state{mu.state_label}"
            )
            legacy_cheb_path = os.path.join(args.out, f"{stem}.chebmom2D")
            cheb_path = os.path.join(args.out, f"{stem}.chebmom")
            dat_path = os.path.join(args.out, f"KuboGreenWood_{stem}JACKSON.dat")
            _core.WriteLegacyChebmom2D(legacy_cheb_path, mu)
        else:
            raise ValueError(f"stage5_bases: unsupported word '{args.word}'")

        json_path = os.path.join(args.out, f"{stem}.basis.json")
        basis_meta = {
            "basis": args.basis,
            "seed": model.seed_name,
            "supercell": list(basis.supercell),
            "mesh": list(basis.supercell),
            "dimension": int(basis.dimension),
            "trace": args.trace,
            "states": int(mu.number_of_states),
            "state_seed": int(rng_seed),
            "window_source": window_source,
            "window": [float(window.lower), float(window.upper)],
            "disorder_w": float(args.disorder_w),
            "disorder_seed": (
                int(args.disorder_seed) if float(args.disorder_w) > 0.0 else None
            ),
            "moments": int(num_moments),
            "broadening_meV": float(args.broadening_meV),
            "word": word_key,
            "operators": {
                op_name: {"source": basis.operator_source(op_name)}
                for op_name in used_ops
            },
        }
        _core.WriteChebmom(cheb_path, mu)
        _core.WriteCurveText(dat_path, energies, values)
        with open(json_path, "w", encoding="utf-8") as fh:
            json.dump(basis_meta, fh, indent=2)
            fh.write("\n")
        print(legacy_cheb_path)
        print(cheb_path)
        print(dat_path)
        print(json_path)
        return 0
    except Exception as exc:  # noqa: BLE001
        print(f"stage5_bases error: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
