"""Reproduce Tutorial 05 (diffusive charge transport on a 2D square lattice) through the Python API.

Imports ``examples/05_charge_transport_square/make_square.py``, generates the
``L = 200`` (``N = 40000``) disordered square-lattice operators for
``W in (0.5, 1.0, 2.0, 4.0)`` eV with disorder seed ``1`` and ``KPM_SEED = 1``,
computes the MSD moments (``M = 256``, ``160`` steps, ``tmax = 150`` fs) and
reconstructs at ``E_F = 2.0`` eV (``80`` meV broadening), and computes the DOS
moments (``M = 1024``) and reconstructs at ``80`` meV broadening.
"""

from __future__ import annotations

import argparse
import importlib.util
import os
import shutil
import sys
import tempfile
from collections.abc import Sequence
from types import ModuleType

from .. import _core
from ..moments import run_plan
from ..reconstruct import dos_curve, time_correlation_curve
from ..words import compile as compile_word

L_SIDE: int = 200
W_LIST: tuple[float, ...] = (0.5, 1.0, 2.0, 4.0)
DISORDER_SEED: int = 1
KPM_SEED: int = 1
M_MSD: int = 256
NTIME: int = 160
TMAX: float = 150.0
EF_TARGET: float = 2.0
M_DOS: int = 1024
BROAD_MEV: float = 80.0


def _load_make_square() -> ModuleType:
    """Import ``examples/05_charge_transport_square/make_square.py``."""
    repo_root = os.path.abspath(
        os.path.join(os.path.dirname(__file__), "..", "..", "..")
    )
    script_path = os.path.join(
        repo_root, "examples", "05_charge_transport_square", "make_square.py"
    )
    if not os.path.isfile(script_path):
        raise FileNotFoundError(
            f"tutorial05: make_square.py not found at {script_path}"
        )
    spec = importlib.util.spec_from_file_location(
        "lsquant_tutorial05_make_square", script_path
    )
    if spec is None or spec.loader is None:
        raise RuntimeError(f"tutorial05: could not load spec for {script_path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _read_bounds(bounds_path: str) -> _core.SpectralWindow:
    """Read ``lo hi`` from a legacy ``BOUNDS`` file."""
    with open(bounds_path, encoding="utf-8") as fh:
        tokens = fh.read().split()
    return _core.SpectralWindow(float(tokens[0]), float(tokens[1]))


def run_tutorial05(
    out_dir: str,
    *,
    l_side: int = L_SIDE,
    w_list: Sequence[float] = W_LIST,
    disorder_seed: int = DISORDER_SEED,
    dos_moments: int = M_DOS,
    kpm_seed: int = KPM_SEED,
    num_moments: int = M_MSD,
    num_states: int = 1,
    window: tuple[float, float] | list[float] | _core.SpectralWindow | str | None = None,
    ntime: int = NTIME,
    tmax: float = TMAX,
    bessel: str = "legacy",
    words: Sequence[str] = (
        "Tr[[X,U(t)]^+ delta(E-H) [X,U(t)]]",
        "Tr[delta(E-H)]",
    ),
    broadening_meV: float = BROAD_MEV,
    fermi_energies: Sequence[float] = (EF_TARGET,),
    legacy_text_round_trip: bool = True,
) -> list[str]:
    """Run Tutorial 05 and write the eight MSD and DOS ``.dat`` curves into ``out_dir``."""
    make_square = _load_make_square()
    os.makedirs(out_dir, exist_ok=True)
    captured_dir = os.path.join(out_dir, "captured")
    os.makedirs(captured_dir, exist_ok=True)

    win_override: _core.SpectralWindow | None = None
    if isinstance(window, _core.SpectralWindow):
        win_override = window
    elif isinstance(window, (list, tuple)) and len(window) == 2:
        win_override = _core.SpectralWindow(float(window[0]), float(window[1]))

    written: list[str] = []
    cap_idx = 0
    tmax_str = f"{float(tmax):g}"
    ef_list = tuple(float(e) for e in fermi_energies) or (EF_TARGET,)

    for w_val in w_list:
        with tempfile.TemporaryDirectory() as stage_dir:
            label = str(
                make_square.build_square(
                    int(l_side),
                    float(w_val),
                    int(disorder_seed),
                    outdir=stage_dir,
                )
            )
            op_dir = os.path.join(stage_dir, "operators")
            win = (
                win_override
                if win_override is not None
                else _read_bounds(os.path.join(stage_dir, "BOUNDS"))
            )

            ham = _core.ReadCsrText(os.path.join(op_dir, f"{label}.HAM.CSR"))
            vx_mat = _core.ReadCsrText(os.path.join(op_dir, f"{label}.VX.CSR"))
            vx_op = _core.RescaledOperator(vx_mat, win)
            dim = ham.dimension

            for w_str in words:
                compiled_probe = compile_word(
                    w_str,
                    basis={"name": "real", "operators": {"VX": vx_op}},
                    params={
                        "expansion": {
                            "moments": int(num_moments),
                            "states": int(num_states),
                            "window": win,
                        },
                        "time": {
                            "steps": int(ntime),
                            "tmax": float(tmax),
                            "bessel": str(bessel),
                        },
                    },
                )
                if compiled_probe.word.is_time_dependent:
                    m_cur = int(num_moments)
                    states_msd = _core.RandomPhaseStates(
                        dim, int(num_states), int(kpm_seed)
                    )
                    mu_msd = run_plan(
                        compiled_probe.plan,
                        ham,
                        win,
                        states_msd,
                        operators=compiled_probe.operators,
                        batch_size=m_cur,
                        symmetrise=False,
                        dt=compiled_probe.dt,
                    )
                    prefix = compiled_probe.reconstruction.legacy_moment_prefix
                    msd_stem = (
                        f"{prefix}{label}KPM_M{m_cur}_state{tmax_str}{states_msd.label}"
                    )
                    msd_mom_path = os.path.join(stage_dir, f"{msd_stem}.chebmomTD")
                    _core.WriteLegacyChebmomTD(msd_mom_path, mu_msd)
                    mu_msd_disk = (
                        _core.ReadChebmom(msd_mom_path)
                        if legacy_text_round_trip
                        else mu_msd
                    )

                    for ef in ef_list:
                        times, msd_vals = time_correlation_curve(
                            mu_msd_disk, float(broadening_meV), ef
                        )
                        msd_dat_name = f"mean{msd_stem}EF{ef:.6f}JACKSON.dat"
                        msd_dat_path = os.path.join(out_dir, msd_dat_name)
                        _core.WriteCurveText(msd_dat_path, times, msd_vals)
                        written.append(msd_dat_path)

                        msd_cap_path = os.path.join(
                            captured_dir, f"{cap_idx:04d}_{msd_dat_name}"
                        )
                        shutil.copyfile(msd_dat_path, msd_cap_path)
                        cap_idx += 1
                else:
                    m_dos = int(dos_moments)
                    compiled_dos = compile_word(
                        w_str,
                        basis={"name": "real", "operators": {"VX": vx_op}},
                        params={
                            "expansion": {
                                "moments": m_dos,
                                "states": int(num_states),
                                "window": win,
                            },
                        },
                    )
                    states_dos = _core.RandomPhaseStates(
                        dim, int(num_states), int(kpm_seed)
                    )
                    mu_dos = run_plan(
                        compiled_dos.plan,
                        ham,
                        win,
                        states_dos,
                        operators=compiled_dos.operators,
                        batch_size=m_dos,
                    )
                    prefix = compiled_dos.reconstruction.legacy_moment_prefix
                    dos_stem = f"{prefix}{label}KPM_M{m_dos}_state{states_dos.label}"
                    dos_mom_path = os.path.join(stage_dir, f"{dos_stem}.chebmom1D")
                    _core.WriteLegacyChebmom1D(dos_mom_path, mu_dos)
                    mu_dos_disk = (
                        _core.ReadChebmom(dos_mom_path)
                        if legacy_text_round_trip
                        else mu_dos
                    )

                    energies, dos_vals = dos_curve(
                        mu_dos_disk, float(broadening_meV)
                    )
                    dos_dat_name = f"mean{dos_stem}JACKSON.dat"
                    dos_dat_path = os.path.join(out_dir, dos_dat_name)
                    _core.WriteCurveText(dos_dat_path, energies, dos_vals)
                    written.append(dos_dat_path)

                    dos_cap_path = os.path.join(
                        captured_dir, f"{cap_idx:04d}_{dos_dat_name}"
                    )
                    shutil.copyfile(dos_dat_path, dos_cap_path)
                    cap_idx += 1

    return written


def main(argv: list[str] | None = None) -> int:
    """Command-line entry point for ``lsquant.apps.tutorial05``."""
    try:
        parser = argparse.ArgumentParser(prog="lsquant.apps.tutorial05")
        parser.add_argument("pos_out", nargs="?", default=None, help="Output directory")
        parser.add_argument(
            "--out", "-o", dest="out", default=None, help="Output directory"
        )
        args = parser.parse_args(argv)
        out_dir = (
            args.out
            if args.out is not None
            else (args.pos_out if args.pos_out is not None else ".")
        )
        paths = run_tutorial05(out_dir)
        for p in paths:
            print(p)
        return 0
    except Exception as exc:  # noqa: BLE001
        print(f"tutorial05 error: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
