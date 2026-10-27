"""Narrow, lazy access to the unchanged origin downfolder.

Only the origin's pure array helpers and FFT/WS Model are reused in-process.
The full automatic shell and SCDM searches remain callable as subprocesses;
no source rewriting, import-path mutation, or global warning suppression leaks.
"""
from __future__ import annotations

from functools import lru_cache
import importlib.util
from pathlib import Path
import subprocess
import sys
import warnings


def origin_directory():
    """Locate the vendored tools in a source checkout; never fetch a substitute."""
    directory = Path(__file__).resolve().parents[3]/"tools"/"downfold_origin"
    if not (directory/"downfold.py").is_file():
        raise FileNotFoundError(
            f"Expected downfold.py under {str(directory)!r}; origin helpers require tools/downfold_origin"
        )
    return directory


@lru_cache(maxsize=1)
def origin_module():
    """Import the immutable source once (requires the pinned sisl environment)."""
    path = origin_directory()/"downfold.py"
    spec = importlib.util.spec_from_file_location("_lsquant_downfold_origin", path)
    module = importlib.util.module_from_spec(spec)
    with warnings.catch_warnings():
        spec.loader.exec_module(module)
    return module


def run_origin_downfolder(rundir, label, *, output, method="polar", arguments=()):
    """Run the complete original automatic reducer, including off-grid validation.

    method='polar' includes isolated/frozen-window searches; 'scdm' selects
    the original SCDM candidate scan. Output files and exit failures belong
    to the origin program, whose main guard is invoked unchanged.
    """
    if method not in ("polar", "scdm"):
        raise ValueError("orthogonal origin method must be polar or scdm")
    script = "downfold.py" if method == "polar" else "downfold_scdm.py"
    return subprocess.run([sys.executable, str(origin_directory()/script),
                           str(rundir), str(label), "--out", str(output),
                           *map(str, arguments)], check=True)
