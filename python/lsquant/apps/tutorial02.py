"""Reproduce Tutorial 02 (Anderson localisation in a 1D disordered chain) through the Python API.

Imports ``examples/02_localization_1d_chain/make_disordered_chain.py``, runs the
README ``N = 512`` single-chain MSD calculation (``W = 1.0``, disorder seed ``1``,
``KPM_SEED = 12345``, ``M = 256``, ``128`` steps, ``tmax = 600`` fs, reconstructed
at ``10`` meV and ``E = 0``) and the ``lsqloc.py`` ensemble over
``N in (256, 1024, 4096)``, realisations ``r = 1 .. 16`` (with ``KPM_SEED = r``),
and energies ``E in (0.0, 1.0, -1.0, 1.9)`` eV.
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
from ..reconstruct import time_correlation_curve
from ..words import compile as compile_word

W_DISORDER: float = 1.0
SIZES: tuple[int, ...] = (256, 1024, 4096)
N_REAL: int = 16
ENERGIES: tuple[float, ...] = (0.0, 1.0, -1.0, 1.9)
M_MOMENTS: int = 256
NTIME: int = 128
TMAX: float = 600.0
BROAD_MEV: float = 10.0


def _load_make_disordered_chain() -> ModuleType:
    """Import ``examples/02_localization_1d_chain/make_disordered_chain.py``."""
    repo_root = os.path.abspath(
        os.path.join(os.path.dirname(__file__), "..", "..", "..")
    )
    script_path = os.path.join(
        repo_root, "examples", "02_localization_1d_chain", "make_disordered_chain.py"
    )
    if not os.path.isfile(script_path):
        raise FileNotFoundError(
            f"tutorial02: make_disordered_chain.py not found at {script_path}"
        )
    spec = importlib.util.spec_from_file_location(
        "lsquant_tutorial02_make_disordered_chain", script_path
    )
    if spec is None or spec.loader is None:
        raise RuntimeError(f"tutorial02: could not load spec for {script_path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _read_bounds(bounds_path: str) -> _core.SpectralWindow:
    """Read ``lo hi`` from a legacy ``BOUNDS`` file."""
    with open(bounds_path, encoding="utf-8") as fh:
        tokens = fh.read().split()
    return _core.SpectralWindow(float(tokens[0]), float(tokens[1]))


def _run_one_chain(
    make_mod: ModuleType,
    n_sites: int,
    w_val: float,
    disorder_seed: int,
    kpm_seed: int,
    energies: tuple[float, ...],
    stage_dir: str,
    dest_dir: str,
    *,
    num_moments: int = M_MOMENTS,
    num_states: int = 1,
    window_override: _core.SpectralWindow | None = None,
    ntime: int = NTIME,
    tmax: float = TMAX,
    bessel: str = "legacy",
    word: str = "Tr[[X,U(t)]^+ delta(E-H) [X,U(t)]]",
    broadening_meV: float = BROAD_MEV,
    legacy_text_round_trip: bool = True,
) -> list[tuple[str, str]]:
    """Build one disordered chain, compute MSD moments, and reconstruct ``energies``."""
    label = str(make_mod.build_chain(n_sites, w_val, disorder_seed, outdir=stage_dir))
    op_dir = os.path.join(stage_dir, "operators")
    window = (
        window_override
        if window_override is not None
        else _read_bounds(os.path.join(stage_dir, "BOUNDS"))
    )

    ham = _core.ReadCsrText(os.path.join(op_dir, f"{label}.HAM.CSR"))
    vx_mat = _core.ReadCsrText(os.path.join(op_dir, f"{label}.VX.CSR"))
    vx_op = _core.RescaledOperator(vx_mat, window)
    states = _core.RandomPhaseStates(n_sites, int(num_states), int(kpm_seed))

    compiled = compile_word(
        word,
        basis={"name": "real", "operators": {"VX": vx_op}},
        params={
            "expansion": {
                "moments": int(num_moments),
                "states": int(num_states),
                "window": window,
            },
            "time": {
                "steps": int(ntime),
                "tmax": float(tmax),
                "bessel": str(bessel),
            },
            "reconstruction": {
                "broadening_meV": float(broadening_meV),
                "fermi_energies": list(energies),
                "legacy_text_round_trip": bool(legacy_text_round_trip),
            },
        },
    )
    mu = run_plan(
        compiled.plan,
        ham,
        window,
        states,
        operators=compiled.operators,
        batch_size=int(num_moments),
        symmetrise=False,
        dt=compiled.dt,
    )
    tmax_str = f"{float(tmax):g}"
    prefix = compiled.reconstruction.legacy_moment_prefix
    stem = f"{prefix}{label}KPM_M{int(num_moments)}_state{tmax_str}{states.label}"
    mom_path = os.path.join(stage_dir, f"{stem}.chebmomTD")
    _core.WriteLegacyChebmomTD(mom_path, mu)
    mu_disk = _core.ReadChebmom(mom_path) if legacy_text_round_trip else mu

    os.makedirs(dest_dir, exist_ok=True)
    results: list[tuple[str, str]] = []
    for ef in energies:
        times, values = time_correlation_curve(mu_disk, float(broadening_meV), ef)
        dat_name = f"mean{stem}EF{ef:.6f}JACKSON.dat"
        dat_path = os.path.join(dest_dir, dat_name)
        _core.WriteCurveText(dat_path, times, values)
        results.append((dat_name, dat_path))
    return results


def run_tutorial02(
    out_dir: str,
    *,
    sizes: Sequence[int] = SIZES,
    w_disorder: float = W_DISORDER,
    realisations: int = N_REAL,
    readme_size: int = 512,
    readme_disorder_seed: int = 1,
    kpm_seed: int = 12345,
    num_moments: int = M_MOMENTS,
    num_states: int = 1,
    window: tuple[float, float] | list[float] | _core.SpectralWindow | str | None = None,
    ntime: int = NTIME,
    tmax: float = TMAX,
    bessel: str = "legacy",
    word: str = "Tr[[X,U(t)]^+ delta(E-H) [X,U(t)]]",
    broadening_meV: float = BROAD_MEV,
    energies: Sequence[float] = ENERGIES,
    legacy_text_round_trip: bool = True,
) -> list[str]:
    """Run Tutorial 02 and write the top-level and captured ``.dat`` files into ``out_dir``."""
    make_mod = _load_make_disordered_chain()
    os.makedirs(out_dir, exist_ok=True)
    captured_dir = os.path.join(out_dir, "captured")
    os.makedirs(captured_dir, exist_ok=True)

    win_override: _core.SpectralWindow | None = None
    if isinstance(window, _core.SpectralWindow):
        win_override = window
    elif isinstance(window, (list, tuple)) and len(window) == 2:
        win_override = _core.SpectralWindow(float(window[0]), float(window[1]))

    sizes_tup = tuple(int(s) for s in sizes)
    energies_tup = tuple(float(e) for e in energies)
    written: list[str] = []

    # 1. README single-chain run: N = 512, W = 1.0, disorder seed = 1, KPM_SEED = 12345, E = 0.0
    readme_ef = (energies_tup[0],) if energies_tup else (0.0,)
    with tempfile.TemporaryDirectory() as stage_dir:
        res_512 = _run_one_chain(
            make_mod,
            int(readme_size),
            float(w_disorder),
            int(readme_disorder_seed),
            int(kpm_seed),
            readme_ef,
            stage_dir,
            out_dir,
            num_moments=int(num_moments),
            num_states=int(num_states),
            window_override=win_override,
            ntime=int(ntime),
            tmax=float(tmax),
            bessel=str(bessel),
            word=str(word),
            broadening_meV=float(broadening_meV),
            legacy_text_round_trip=bool(legacy_text_round_trip),
        )
        for _, path in res_512:
            written.append(path)

    # 2. lsqloc.py ensemble: N in (256, 1024, 4096), r in 1..16 (KPM_SEED = r), E in (0, 1, -1, 1.9)
    cap_idx = 0
    for n_sites in sizes_tup:
        for r in range(1, int(realisations) + 1):
            r_dir = os.path.join(out_dir, f"r{r:02d}")
            with tempfile.TemporaryDirectory() as stage_dir:
                res_r = _run_one_chain(
                    make_mod,
                    n_sites,
                    float(w_disorder),
                    r,
                    r,
                    energies_tup,
                    stage_dir,
                    r_dir,
                    num_moments=int(num_moments),
                    num_states=int(num_states),
                    window_override=win_override,
                    ntime=int(ntime),
                    tmax=float(tmax),
                    bessel=str(bessel),
                    word=str(word),
                    broadening_meV=float(broadening_meV),
                    legacy_text_round_trip=bool(legacy_text_round_trip),
                )
                for dat_name, dat_path in res_r:
                    cap_path = os.path.join(
                        captured_dir, f"{cap_idx:04d}_{dat_name}"
                    )
                    shutil.copyfile(dat_path, cap_path)
                    written.append(cap_path)
                    cap_idx += 1

    return written


def main(argv: list[str] | None = None) -> int:
    """Command-line entry point for ``lsquant.apps.tutorial02``."""
    try:
        parser = argparse.ArgumentParser(prog="lsquant.apps.tutorial02")
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
        paths = run_tutorial02(out_dir)
        for p in paths:
            print(p)
        return 0
    except Exception as exc:  # noqa: BLE001
        print(f"tutorial02 error: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
