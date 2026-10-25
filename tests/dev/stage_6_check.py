"""Developer self-check for Stage 6 (non-orthogonal metric and overlap appliers).

Executes the four checks specified in Section 7 of ``handoff/0006-nonorthogonal-metric/brief.md``
on the rotated ``chain1d`` (``N = 2048``) ``real`` and ``complex`` fixtures:
1. Generalised spectrum invariance ``max |eigvalsh(H', S') - eigvalsh(H)| < 1e-12``.
2. ``stage6_nonorth`` exact-trace DOS moments (``M = 64``, ``window = [-2, 2]``) via
   ``lowdin/exact`` and ``generalised/exact`` compared against G1 (``mu_m = delta_m0``)
   to ``< 1e-12``.
3. Approximate appliers ``ischo`` and ``chebyshev`` (``M_S = 60``, plus error vs ``M_S``
   and ``tol``) reporting applier residuals ``||S (S^-1 v) - v|| / ||v||`` and
   ``||(S^-1/2)^2 S v - v|| / ||v||`` alongside moment errors.
4. Direct C++ comparison of ``ExactInverseOverlap`` vs ``ChebyshevInverseOverlap``
   (``M_S = 60``) on ten random vectors:
   ``max ||S^-1_cheb v - S^-1_exact v|| / ||S^-1_exact v|| < 1e-10``.
"""

from __future__ import annotations

import argparse
import os
from pathlib import Path
import subprocess
import sys
import tempfile

os.environ.setdefault("OPENBLAS_NUM_THREADS", "4")
os.environ.setdefault("MKL_NUM_THREADS", "4")
os.environ.setdefault("OMP_NUM_THREADS", "8")

import numpy as np
import scipy.linalg as la

_REPO_ROOT = Path(__file__).resolve().parents[2]
_RAW_BUILD_DIR = os.environ.get("LSQUANT_BUILD_DIR", "build")
_BUILD_DIR = (
    Path(_RAW_BUILD_DIR).resolve()
    if os.path.isabs(_RAW_BUILD_DIR)
    else (_REPO_ROOT / _RAW_BUILD_DIR).resolve()
)
if not (_BUILD_DIR / "lsquant_core").is_dir() and not (_BUILD_DIR / "stage6_nonorth").is_file():
    for _cand in (_REPO_ROOT / "build" / "golden", _REPO_ROOT / "build" / "native"):
        if (_cand / "lsquant_core").is_dir():
            _BUILD_DIR = _cand.resolve()
            break
for _p in (_BUILD_DIR / "lsquant_core", _BUILD_DIR, _REPO_ROOT / "python"):
    if _p.is_dir() and str(_p) not in sys.path:
        sys.path.insert(0, str(_p))

from lsquant.models.nonorthogonal import read_csr_text  # noqa: E402


def read_legacy_chebmom1d(path: Path) -> np.ndarray:
    """Parse a legacy ``.chebmom1D`` file and return complex moments ``mu[0..M-1]``."""
    lines = path.read_text(encoding="utf-8").strip().splitlines()
    if len(lines) < 2:
        raise ValueError(f"read_legacy_chebmom1d: truncated file {path}")
    m_count = int(lines[1].split()[0])
    moments = np.zeros(m_count, dtype=np.complex128)
    for idx in range(m_count):
        parts = lines[2 + idx].split()
        moments[idx] = float(parts[0]) + 1j * float(parts[1])
    return moments


def parse_kv_line(line: str) -> dict[str, str]:
    """Extract ``key=value`` tokens from a single summary line."""
    out: dict[str, str] = {}
    for token in line.strip().split():
        if "=" in token:
            k, v = token.split("=", 1)
            out[k] = v
    return out


def run_dos(
    exe: Path,
    fixture_dir: Path,
    route: str,
    applier: str,
    out_dir: Path,
    order: int = 60,
    tol: float = 1e-13,
) -> tuple[np.ndarray, str]:
    """Run ``stage6_nonorth`` for ``word=dos`` on ``chain1d`` and return moments + stdout."""
    cmd = [
        str(exe),
        "--fixture",
        str(fixture_dir),
        "--label",
        "chain1d",
        "--window",
        "-2",
        "2",
        "--route",
        route,
        "--applier",
        applier,
        "--order",
        str(order),
        "--tol",
        f"{tol:.1e}",
        "--trace",
        "exact",
        "--moments",
        "64",
        "--word",
        "dos",
        "--out",
        str(out_dir),
    ]
    res = subprocess.run(cmd, check=True, capture_output=True, text=True)
    chebmom_path = out_dir / "SpectralOp1chain1dKPM_M64_stateexact.chebmom1D"
    return read_legacy_chebmom1d(chebmom_path), res.stdout.strip()


def run_check_appliers(
    exe: Path,
    fixture_dir: Path,
    order: int = 60,
    tol: float = 1e-13,
    applier: str = "exact",
) -> dict[str, float]:
    """Run ``stage6_nonorth --check-appliers`` and return parsed floating-point metrics."""
    cmd = [
        str(exe),
        "--fixture",
        str(fixture_dir),
        "--label",
        "chain1d",
        "--applier",
        applier,
        "--order",
        str(order),
        "--tol",
        f"{tol:.1e}",
        "--check-appliers",
    ]
    res = subprocess.run(cmd, check=True, capture_output=True, text=True)
    kv = parse_kv_line(res.stdout)
    metrics: dict[str, float] = {}
    for key, val in kv.items():
        if key.startswith("max_"):
            metrics[key] = float(val)
    return metrics


def main() -> int:
    parser = argparse.ArgumentParser(description="Run Stage-6 developer self-check.")
    parser.add_argument(
        "--root",
        type=Path,
        default=_REPO_ROOT,
        help="Repository root directory.",
    )
    default_exe = (
        _BUILD_DIR / "stage6_nonorth"
        if (_BUILD_DIR / "stage6_nonorth").is_file()
        else _BUILD_DIR / "lsquant_core" / "stage6_nonorth"
    )
    parser.add_argument(
        "--exe",
        type=Path,
        default=default_exe,
        help="Path to the stage6_nonorth executable.",
    )
    args = parser.parse_args()

    root = args.root.resolve()
    exe = args.exe.resolve()
    if not exe.is_file():
        raise FileNotFoundError(f"stage6_nonorth executable not found at {exe}")

    orth_dir = root / "goldens" / "stage0" / "chain1d" / "operators"
    g1_path = root / "goldens" / "stage0" / "chain1d" / "SpectralOp1chain1dKPM_M64_stateexact.chebmom1D"
    mu_g1 = read_legacy_chebmom1d(g1_path)
    delta_m0 = np.zeros_like(mu_g1)
    delta_m0[0] = 1.0 + 0.0j

    h_orth = read_csr_text(orth_dir / "chain1d.HAM.CSR").toarray()
    eigs_orth = la.eigvalsh(h_orth, driver="evd")

    print("=== Stage 6 Self-Check: Rotated chain1d (N = 2048) ===", flush=True)

    with tempfile.TemporaryDirectory(prefix="stage6_check_") as tmp:
        tmp_dir = Path(tmp)

        for variant in ("real", "complex"):
            fix_dir = root / "goldens" / "stage0" / "nonorth" / "chain1d" / variant
            h_prime = read_csr_text(fix_dir / "chain1d.HAM.CSR").toarray()
            s_prime = read_csr_text(fix_dir / "chain1d.S.CSR").toarray()

            # 1. Generalised eigenvalue check.
            eigs_gen = la.eigvalsh(h_prime, s_prime, driver="gvd")
            max_eig_diff = float(np.max(np.abs(eigs_gen - eigs_orth)))
            print(
                f"[1] {variant}: max |eigvalsh(H', S') - eigvalsh(H)| = {max_eig_diff:.6e}",
                flush=True,
            )
            if not (max_eig_diff < 1e-12):
                raise AssertionError(
                    f"Generalised eigenvalue check failed for {variant}: {max_eig_diff}"
                )

            # 4. Direct applier check on 10 random vectors (M_S = 60, tol = 1e-13).
            metrics_60 = run_check_appliers(
                exe, fix_dir, order=60, tol=1e-13, applier="exact"
            )
            rel_cheb_vs_exact = metrics_60["max_rel_diff_sinv_cheb_vs_exact"]
            print(
                f"[4] {variant}: max ||S^-1_cheb v - S^-1_exact v|| / ||S^-1_exact v|| (M_S=60) = "
                f"{rel_cheb_vs_exact:.6e}",
                flush=True,
            )
            if not (rel_cheb_vs_exact < 1e-10):
                raise AssertionError(
                    f"ChebyshevInverseOverlap vs ExactInverseOverlap failed for {variant}: "
                    f"{rel_cheb_vs_exact}"
                )

            # 2. Exact appliers through lowdin and generalised routes vs G1 (delta_m0).
            for route in ("lowdin", "generalised"):
                run_out = tmp_dir / f"{variant}_{route}_exact"
                run_out.mkdir(parents=True, exist_ok=True)
                mu, _ = run_dos(exe, fix_dir, route, "exact", run_out)
                err_g1 = float(np.max(np.abs(mu - mu_g1)))
                err_delta = float(np.max(np.abs(mu - delta_m0)))
                res_key = (
                    "max_res_sinvsqrt_exact"
                    if route == "lowdin"
                    else "max_res_sinv_exact"
                )
                print(
                    f"[2] {variant} route={route:<11} applier=exact     "
                    f"applier_res={metrics_60[res_key]:.6e} "
                    f"max|mu - G1|={err_g1:.6e} max|mu - delta_m0|={err_delta:.6e}",
                    flush=True,
                )
                if not (err_g1 < 1e-12 and err_delta < 1e-12):
                    raise AssertionError(
                        f"Exact DOS moment check failed for {variant}/{route}: "
                        f"err_g1={err_g1}, err_delta={err_delta}"
                    )

            # 3. Approximate appliers (ischo and chebyshev, M_S = 60, plus M_S=12 / tol=1e-6 comparison).
            orders_to_run = (12, 60) if variant == "real" else (60,)
            for order in orders_to_run:
                m_ord = (
                    metrics_60
                    if order == 60
                    else run_check_appliers(
                        exe, fix_dir, order=order, tol=1e-13, applier="chebyshev"
                    )
                )
                for route, res_key in (
                    ("lowdin", "max_res_sinvsqrt_cheb"),
                    ("generalised", "max_res_sinv_cheb"),
                ):
                    run_out = tmp_dir / f"{variant}_{route}_cheb_{order}"
                    run_out.mkdir(parents=True, exist_ok=True)
                    mu, out_line = run_dos(exe, fix_dir, route, "chebyshev", run_out, order=order)
                    err_delta = float(np.max(np.abs(mu - delta_m0)))
                    kv = parse_kv_line(out_line)
                    eff_str = (
                        f", eff_order={kv['effective_order']}"
                        if "effective_order" in kv
                        else ""
                    )
                    print(
                        f"[3] {variant} route={route:<11} applier=chebyshev (M_S={order}{eff_str}) "
                        f"applier_res={m_ord[res_key]:.6e} max|mu - delta_m0|={err_delta:.6e}",
                        flush=True,
                    )

            tols_to_run = (1e-6, 1e-13) if variant == "real" else (1e-13,)
            for tol in tols_to_run:
                m_tol = (
                    metrics_60
                    if tol == 1e-13
                    else run_check_appliers(
                        exe, fix_dir, order=60, tol=tol, applier="ischo"
                    )
                )
                run_out = tmp_dir / f"{variant}_generalised_ischo_{tol:.0e}"
                run_out.mkdir(parents=True, exist_ok=True)
                mu, _ = run_dos(exe, fix_dir, "generalised", "ischo", run_out, tol=tol)
                err_delta = float(np.max(np.abs(mu - delta_m0)))
                print(
                    f"[3] {variant} route=generalised applier=ischo     (tol={tol:.0e}) "
                    f"applier_res={m_tol['max_res_sinv_ischo']:.6e} "
                    f"max|mu - delta_m0|={err_delta:.6e}",
                    flush=True,
                )

    print("=== Stage 6 Self-Check PASSED ===", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
