"""Reproduce Tutorial 01 (1D chain density of states) through the Python API.

Builds the periodic 1D tight-binding chains for ``N = 256, 1024, 4096`` by
importing ``examples/01_dos_1d_chain/make_chain.py``, evaluates ``M = 512``
exact-trace Chebyshev moments through :class:`~lsquant._core.MomentEngine`,
writes the legacy ``.chebmom1D`` files at the top of the output directory, and
reconstructs the density of states at ``20, 80, 24, 10`` meV under
``<out>/eta<meV>/`` with :func:`~lsquant.reconstruct.dos_curve`.
"""

from __future__ import annotations

import argparse
import importlib.util
import os
import shutil
import sys
from collections.abc import Sequence
from types import ModuleType

from .. import _core
from ..moments import run_plan, window_from_gershgorin
from ..reconstruct import dos_curve
from ..words import compile as compile_word


def _load_make_chain() -> ModuleType:
    """Import ``examples/01_dos_1d_chain/make_chain.py`` from the repository tree.

    Returns
    -------
    types.ModuleType
        Loaded ``make_chain`` module.

    Warnings
    --------
    Raises ``FileNotFoundError`` if ``examples/01_dos_1d_chain/make_chain.py``
    cannot be located relative to the package path.

    Notes
    -----
    Pitfall: Loads ``make_chain.py`` via :mod:`importlib.util` so that no copy
    of the model generator exists in the library.
    """
    repo_root = os.path.abspath(
        os.path.join(os.path.dirname(__file__), "..", "..", "..")
    )
    script_path = os.path.join(repo_root, "examples", "01_dos_1d_chain", "make_chain.py")
    if not os.path.isfile(script_path):
        raise FileNotFoundError(f"tutorial01: make_chain.py not found at {script_path}")

    spec = importlib.util.spec_from_file_location("lsquant_tutorial01_make_chain", script_path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"tutorial01: could not load spec for {script_path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def run_tutorial01(
    out_dir: str,
    *,
    sizes: Sequence[int] = (256, 1024, 4096),
    hopping_t: float = -1.0,
    num_moments: int = 512,
    window: tuple[float, float] | list[float] | _core.SpectralWindow | str | None = None,
    word: str = "Tr[delta(E-H)]",
    broadening_meV: float = 20.0,
) -> list[str]:
    """Generate the Tutorial 01 operators, exact moments, and DOS curves in ``out_dir``.

    Solves the problem of reproducing the moment files and reconstructed density
    of states curves of ``examples/01_dos_1d_chain`` using the new Python API.
    Writes the ``.chebmom1D`` files with their legacy names at the top of
    ``out_dir`` and the reconstructed curves under ``<out_dir>/eta<meV>/`` with
    their legacy names (``meanmoments_chain1d_N<N>_M512JACKSON.dat``).

    Parameters
    ----------
    out_dir : str
        Destination directory for generated operators, moment files, and curves.

    Returns
    -------
    list[str]
        Ordered list of output file paths written.

    Warnings
    --------
    Overwrites existing tutorial output files in ``out_dir``.

    Notes
    -----
    Pitfall: Because ``N = 256`` is reconstructed at four different broadenings
    (``20, 80, 24, 10`` meV) with the same legacy filename
    ``meanmoments_chain1d_N256_M512JACKSON.dat``, each broadening is placed in
    its own ``eta<meV>/`` subdirectory.
    """
    make_chain = _load_make_chain()
    os.makedirs(out_dir, exist_ok=True)
    op_dir = os.path.join(out_dir, "operators")
    os.makedirs(op_dir, exist_ok=True)
    captured_dir = os.path.join(out_dir, "captured")
    os.makedirs(captured_dir, exist_ok=True)

    sizes_tup = tuple(int(s) for s in sizes)
    num_moments = int(num_moments)
    hopping_t = float(hopping_t)
    broadening_meV = float(broadening_meV)
    if isinstance(window, _core.SpectralWindow):
        default_win = window
    elif isinstance(window, (list, tuple)) and len(window) == 2:
        default_win = _core.SpectralWindow(float(window[0]), float(window[1]))
    else:
        default_win = _core.SpectralWindow(-2.0 * abs(hopping_t), 2.0 * abs(hopping_t))

    moments_by_size: dict[int, _core.ChebyshevMoments] = {}
    spec_paths: dict[int, str] = {}
    short_paths: dict[int, str] = {}
    written: list[str] = []

    for n_sites in sizes_tup:
        label = f"chain1d_N{n_sites}"
        csr_path = os.path.join(op_dir, f"{label}.HAM.CSR")
        desc_path = os.path.join(op_dir, f"{label}.HAM.desc")
        bounds_path = os.path.join(out_dir, "BOUNDS")
        state_path = os.path.join(out_dir, f"exact_{label}")

        make_chain.write_csr(csr_path, n_sites, hopping_t)
        make_chain.write_desc(desc_path, label, n_sites, hopping_t, csr_path=csr_path)
        make_chain.write_bounds(bounds_path, hopping_t)
        make_chain.write_state(state_path, n_sites)

        ham = _core.ReadCsrText(csr_path)
        win_n = (
            window_from_gershgorin(ham)
            if window == "gershgorin"
            else default_win
        )
        states = _core.LocalStates(n_sites, list(range(n_sites)))
        compiled = compile_word(
            word,
            basis={"name": "real"},
            params={"moments": num_moments, "states": n_sites, "window": win_n},
        )
        mu = run_plan(
            compiled.plan,
            ham,
            win_n,
            states,
            operators=compiled.operators,
            batch_size=num_moments,
        )
        mu.state_label = f"exact_{label}"
        moments_by_size[n_sites] = mu

        spec_name = f"SpectralOp1{label}KPM_M{num_moments}_stateexact_{label}.chebmom1D"
        short_name = f"moments_{label}_M{num_moments}.chebmom1D"
        spec_path = os.path.join(out_dir, spec_name)
        short_path = os.path.join(out_dir, short_name)

        _core.WriteLegacyChebmom1D(spec_path, mu)
        shutil.copyfile(spec_path, short_path)

        spec_paths[n_sites] = spec_path
        short_paths[n_sites] = short_path
        written.append(spec_path)
        written.append(short_path)

    primary_size = 256 if 256 in moments_by_size else sizes_tup[0]
    # README step 3: reconstruct primary_size (256) at 120 meV at the top level of out_dir
    e_120, v_120 = dos_curve(moments_by_size[primary_size], 120.0)
    top_120_path = os.path.join(
        out_dir,
        f"meanSpectralOp1chain1d_N{primary_size}KPM_M{num_moments}_stateexact_chain1d_N{primary_size}JACKSON.dat",
    )
    _core.WriteCurveText(top_120_path, e_120, v_120)
    written.append(top_120_path)

    # Captured 0000..0011: figure_moments reads (twice: standalone step + lsqplot.py)
    cap_idx = 0
    for _ in range(2):
        for n_sites in sizes_tup:
            s_path = spec_paths[n_sites]
            sh_path = short_paths[n_sites]
            shutil.copyfile(
                s_path,
                os.path.join(captured_dir, f"{cap_idx:04d}_{os.path.basename(s_path)}"),
            )
            cap_idx += 1
            shutil.copyfile(
                sh_path,
                os.path.join(captured_dir, f"{cap_idx:04d}_{os.path.basename(sh_path)}"),
            )
            cap_idx += 1

    # Panel (a): fixed broadening (default 20 meV) across sizes (captures 0012..0017)
    for n_sites in sizes_tup:
        eta = broadening_meV
        energies, values = dos_curve(moments_by_size[n_sites], eta)
        eta_dir = os.path.join(out_dir, f"eta{int(eta)}")
        os.makedirs(eta_dir, exist_ok=True)
        dat_name = f"meanmoments_chain1d_N{n_sites}_M{num_moments}JACKSON.dat"
        dat_path = os.path.join(eta_dir, dat_name)
        _core.WriteCurveText(dat_path, energies, values)
        written.append(dat_path)

        s_path = spec_paths[n_sites]
        shutil.copyfile(
            s_path,
            os.path.join(captured_dir, f"{cap_idx:04d}_{os.path.basename(s_path)}"),
        )
        cap_idx += 1
        shutil.copyfile(
            dat_path,
            os.path.join(captured_dir, f"{cap_idx:04d}_{dat_name}"),
        )
        cap_idx += 1

    # Panel (b): primary_size (256) at broadenings 80, 24, 10 meV (captures 0018..0021)
    s_path_256 = spec_paths[primary_size]
    shutil.copyfile(
        s_path_256,
        os.path.join(captured_dir, f"{cap_idx:04d}_{os.path.basename(s_path_256)}"),
    )
    cap_idx += 1
    for eta in (80.0, 24.0, 10.0):
        energies, values = dos_curve(moments_by_size[primary_size], eta)
        eta_dir = os.path.join(out_dir, f"eta{int(eta)}")
        os.makedirs(eta_dir, exist_ok=True)
        dat_name = f"meanmoments_chain1d_N{primary_size}_M{num_moments}JACKSON.dat"
        dat_path = os.path.join(eta_dir, dat_name)
        _core.WriteCurveText(dat_path, energies, values)
        written.append(dat_path)
        shutil.copyfile(
            dat_path,
            os.path.join(captured_dir, f"{cap_idx:04d}_{dat_name}"),
        )
        cap_idx += 1

    return written


def main(argv: list[str] | None = None) -> int:
    """Command-line entry point for ``lsquant.apps.tutorial01``.

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
    Pitfall: Accepts both ``--out <dir>`` and a positional ``<dir>`` argument,
    defaulting to the current working directory ``.``.
    """
    try:
        parser = argparse.ArgumentParser(prog="lsquant.apps.tutorial01")
        parser.add_argument("pos_out", nargs="?", default=None, help="Output directory")
        parser.add_argument("--out", "-o", dest="out", default=None, help="Output directory")
        args = parser.parse_args(argv)
        out_dir = args.out if args.out is not None else (args.pos_out if args.pos_out is not None else ".")
        paths = run_tutorial01(out_dir)
        for p in paths:
            print(p)
        return 0
    except Exception as exc:  # noqa: BLE001
        print(f"tutorial01 error: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
