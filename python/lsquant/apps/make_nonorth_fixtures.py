"""CLI entry point generating the Stage-6 non-orthogonal CSR fixtures and manifests.

Generates ``real`` and ``complex`` rotated basis fixtures for ``chain1d`` (``N = 2048``)
and ``haldane`` (``N = 512``) under ``goldens/stage0/nonorth/<model>/<variant>/``,
verifying generalised spectrum invariance ``max |eigvalsh(H', S') - eigvalsh(H)| < 1e-12``
and ``cond(S') in [5.0, 10.0]``.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import sys

os.environ.setdefault("OPENBLAS_NUM_THREADS", "2")
os.environ.setdefault("MKL_NUM_THREADS", "2")
os.environ.setdefault("OMP_NUM_THREADS", "2")

import numpy as np
import scipy.linalg as la

# Ensure build/lsquant_core and python/ are importable when invoked as a script.
_REPO_ROOT = Path(__file__).resolve().parents[3]
for _p in (_REPO_ROOT / "build" / "lsquant_core", _REPO_ROOT / "python"):
    if _p.is_dir() and str(_p) not in sys.path:
        sys.path.insert(0, str(_p))

from lsquant.models.nonorthogonal import (  # noqa: E402
    make_local_transform,
    read_csr_text,
    rotate_basis,
    write_csr_text,
)


def _sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def build_fixture_variant(
    model: str,
    src_dir: Path,
    out_dir: Path,
    op_names: list[str],
    window: list[float],
    alpha: float,
    complex_entries: bool,
    seed: int,
) -> dict:
    """Build and write one non-orthogonal fixture variant and its ``manifest.json``."""
    ops = {name: read_csr_text(src_dir / f"{model}.{name}.CSR") for name in op_names}
    h_orth = ops["HAM"]
    dim = int(h_orth.shape[0])

    r_mat = make_local_transform(
        h_orth, alpha=alpha, complex_entries=complex_entries, seed=seed
    )
    rotated = rotate_basis(ops, r_mat)

    out_dir.mkdir(parents=True, exist_ok=True)
    written_files = [f"{model}.S.CSR"] + [f"{model}.{name}.CSR" for name in op_names]
    write_csr_text(out_dir / f"{model}.S.CSR", rotated["S"])
    for name in op_names:
        write_csr_text(out_dir / f"{model}.{name}.CSR", rotated[name])
    sha256_map = {fname: _sha256_file(out_dir / fname) for fname in written_files}

    s_dense = rotated["S"].toarray()
    h_prime_dense = rotated["HAM"].toarray()
    h_orth_dense = h_orth.toarray()

    s_eigs = la.eigvalsh(s_dense, driver="evd")
    lambda_min_s = float(np.min(s_eigs))
    lambda_max_s = float(np.max(s_eigs))
    if lambda_min_s <= 0.0:
        raise RuntimeError(
            f"{model}/{out_dir.name}: S' is not positive definite (lambda_min={lambda_min_s})"
        )
    cond_s = float(lambda_max_s / lambda_min_s)

    orth_eigs = la.eigvalsh(h_orth_dense, driver="evd")
    gen_eigs = la.eigvalsh(h_prime_dense, s_dense, driver="gvd")
    max_eig_diff = float(np.max(np.abs(gen_eigs - orth_eigs)))

    if not (5.0 <= cond_s <= 10.0):
        raise RuntimeError(
            f"{model}/{out_dir.name}: cond(S') = {cond_s:.6f} outside [5.0, 10.0]"
        )
    if max_eig_diff >= 1e-12:
        raise RuntimeError(
            f"{model}/{out_dir.name}: max |eigvalsh(H', S') - eigvalsh(H)| = "
            f"{max_eig_diff:.3e} >= 1e-12"
        )

    manifest = {
        "model": model,
        "variant": "complex" if complex_entries else "real",
        "dimension": dim,
        "seed": int(seed),
        "alpha": float(alpha),
        "window": [float(window[0]), float(window[1])],
        "lambda_min(S')": lambda_min_s,
        "lambda_max(S')": lambda_max_s,
        "cond(S')": cond_s,
        "nnz(S')": int(rotated["S"].nnz),
        "nnz(H')": int(rotated["HAM"].nnz),
        "max |eigvalsh(H', S') - eigvalsh(H)|": max_eig_diff,
        "sha256": sha256_map,
    }
    manifest_path = out_dir / "manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    return manifest


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Generate Stage-6 non-orthogonal CSR fixtures and manifests."
    )
    parser.add_argument(
        "--root",
        type=Path,
        default=_REPO_ROOT,
        help="Repository root directory (default: auto-detected).",
    )
    parser.add_argument("--alpha", type=float, default=0.55, help="Coupling amplitude alpha.")
    parser.add_argument("--seed", type=int, default=12345, help="Deterministic RNG seed.")
    args = parser.parse_args()

    root = args.root.resolve()
    specs = [
        (
            "chain1d",
            root / "goldens" / "stage0" / "chain1d" / "operators",
            ["HAM", "VX"],
            [-2.0, 2.0],
        ),
        (
            "haldane",
            root / "goldens" / "stage0" / "haldane" / "operators",
            ["HAM", "VX", "VY"],
            [-4.5, 4.5],
        ),
    ]

    for model, src_dir, op_names, window in specs:
        for variant, is_complex in (("real", False), ("complex", True)):
            out_dir = root / "goldens" / "stage0" / "nonorth" / model / variant
            manifest = build_fixture_variant(
                model=model,
                src_dir=src_dir,
                out_dir=out_dir,
                op_names=op_names,
                window=window,
                alpha=args.alpha,
                complex_entries=is_complex,
                seed=args.seed,
            )
            lmin = manifest["lambda_min(S')"]
            lmax = manifest["lambda_max(S')"]
            cond = manifest["cond(S')"]
            nnz_s = manifest["nnz(S')"]
            nnz_h = manifest["nnz(H')"]
            ediff = manifest["max |eigvalsh(H', S') - eigvalsh(H)|"]
            print(
                f"make_nonorth_fixtures: {model}/{variant} N={manifest['dimension']} "
                f"spec(S')=[{lmin:.6f}, {lmax:.6f}] "
                f"cond(S')={cond:.6f} "
                f"nnz(S')={nnz_s} nnz(H')={nnz_h} "
                f"max_eig_diff={ediff:.3e}",
                flush=True,
            )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
