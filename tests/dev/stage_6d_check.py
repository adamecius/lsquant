#!/usr/bin/env python3
"""The one stage-6d developer self-check, intended for the Aron tester.

PYTHONPATH=python:build:build/lsquant_core python tests/dev/stage_6d_check.py
Defaults deliberately retain the strict requested 1e-3 plateau comparison.
No developer result is claimed until the tester executes this file.
"""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import sys
import numpy as np


def augmented_kane_mele():
    """KM plus TWO spinful high-energy orbitals per site, then local rotation.

    The added orbitals are flat, decoupled and at +30/+60 eV. A fixed
    invertible spin-scalar intra-site T mixes all three zetas, making the
    physical S nontrivial while preserving the exact low-energy response.
    Each site has [low up/down, extra1 up/down, extra2 up/down]. This toy
    certifies a decoupled complement; it cannot establish Hall invariance
    after removing virtual transitions of a real material.
    """
    from stage_6b_check import kane_mele_arrays
    from lsquant.models.siesta import SiestaModel
    lattice, centres, R, H4, S4 = kane_mele_arrays()
    H = np.zeros((len(R), 12, 12), complex)
    S = np.zeros_like(H)
    low = np.array([0, 1, 6, 7])
    H[:, low[:, None], low[None, :]] = H4
    zero = np.flatnonzero(np.all(R == 0, axis=1))[0]
    H[zero, [2, 3, 4, 5, 8, 9, 10, 11], [2, 3, 4, 5, 8, 9, 10, 11]] = [30, 30, 60, 60, 30, 30, 60, 60]
    S[zero] = np.eye(12)
    local = np.kron(np.array([[1.3, .2, .1], [.1, .8, .15], [0., .1, 1.1]]), np.eye(2))
    T = np.kron(np.eye(2), local)
    H, S = T.conj().T @ H @ T, T.conj().T @ S @ T
    tau = np.repeat(centres[[0, 2]], 6, axis=0)
    A = np.zeros((3, *H.shape), complex)
    return SiestaModel(lattice, tau, R, H, S, A, "KM_extra_rotated")


def check_symbolic():
    """(a) Arbitrary-dimension NC SymPy identities and Fourier gauge."""
    from lsquant.metric_forms import (spectral_projection_forms, spectral_projection_proof,
                                      orthogonal_projection_forms)
    proofs = spectral_projection_proof()
    print("(a) symbolic residuals:", {key: str(value) for key, value in proofs.items()}, flush=True)
    assert all(value == 0 for value in proofs.values()), proofs
    print("(a) symbolic forms:", spectral_projection_forms(), flush=True)
    assert orthogonal_projection_forms()["cell_phase_residual"] == 0


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mesh", type=int, default=12)
    parser.add_argument("--moments", type=int, default=1024)
    parser.add_argument("--broadening-meV", type=float, default=10.)
    parser.add_argument("--backend", choices=("tables", "folded"), default="tables",
                        help="uniform-grid Bastin contraction; the accuracy assertions are unchanged")
    args = parser.parse_args()
    if min(args.mesh, args.moments) < 2 or args.broadening_meV <= 0:
        parser.error("mesh/moments >=2 and positive broadening required")
    root = Path(__file__).resolve().parents[2]
    build = Path(os.environ.get("LSQUANT_BUILD_DIR", root/"build"))
    for path in (root/"python", build, build/"lsquant_core"):
        sys.path.insert(0, str(path))
    from scipy.linalg import eigvalsh
    from lsquant.bases import BlochMetricBasis
    from lsquant.projection import spectral_projection, downfold_orthogonal, LocalMatrix, validate_bands, word_trace
    from lsquant.projection.bases import ProjectedBlochBasis
    from lsquant.projection.runfile import exact_word_curve
    from lsquant.projection.words import compute_projection_moments
    from lsquant.projection.spectral import adjoint
    from lsquant.apps.stage6b_siesta import compute_moments, SHC_WORD
    from lsquant.reconstruct import shc
    from lsquant.runfile import RunFile
    check_symbolic()
    model = augmented_kane_mele()
    source = BlochMetricBasis(model, [args.mesh, args.mesh, 1], velocity_mode="covariant")
    window = [-3.1, 3.1]
    area = float(np.linalg.norm(np.cross(model.uc[0], model.uc[1])))
    fermi_grid = np.array([-.1, 0., .1])
    original_oracle, _ = exact_word_curve(source, SHC_WORD, fermi_grid,
                                          broadening=args.broadening_meV/1000, cell_area=area)
    # Compare both exact metric routes against the same finite-mesh oracle.
    results = []
    for route in ("lowdin", "generalised"):
        source = BlochMetricBasis(model, [args.mesh, args.mesh, 1], velocity_mode="covariant", route=route)
        for kind in ("hard", "shift"):
            projection = spectral_projection(source.hk_array, source.sk_array, window, kind=kind)
            receipt = projection.receipt()
            assert max(receipt["residuals"].values()) < 1e-10
            assert np.array_equal(projection.overlap, source.sk_array)
            E = projection.eigenvalues
            assert E.min() >= window[0]-1e-10 and E.max() <= window[1]+1e-10
            outside = ~projection.kept
            if kind == "hard":
                assert np.all(E[outside] == 0)
            else:
                assert np.max(E) <= window[1]
            # Every physical operator, including current, is projected itself.
            for name in ("VX", "VY", "SZ", "JXSZ"):
                projected = projection.project_operator(source.operator_blocks(name))
                C = projection.eigenvectors
                original = adjoint(C) @ source.operator_blocks(name) @ C
                selected = projection.kept[:, :, None]*original*projection.kept[:, None, :]
                assert np.max(np.abs(adjoint(C) @ projected @ C-selected)) < 1e-10
            # Generic word diagnostic on a smooth spectral word, independent
            # of SHC-specific reconstruction or oracle normalisation.
            op = source.operator_blocks("SZ")
            first = np.exp(-source.eigenvalues**2)
            before = word_trace(source.eigenvalues, source.eigenvectors, [op], [first])
            after = word_trace(E, projection.eigenvectors, [projection.project_operator(op)], [np.exp(-E**2)])
            assert abs(before-after) < 1e-10
            # Non-centred hard projection must include complement zero modes.
            offset = spectral_projection(source.hk_array, source.sk_array, [20., 40.], kind="hard")
            assert offset.eigenvalues.min() == 0.
            selected = ProjectedBlochBasis(source, dict(kind=kind, window=window))
            curve, _ = exact_word_curve(selected, SHC_WORD, fermi_grid,
                                        broadening=args.broadening_meV/1000, cell_area=area)
            oracle_error = float(np.max(np.abs(curve-original_oracle)))
            assert oracle_error < 1e-3
            results.append(dict(route=route, kind=kind, oracle=curve.tolist(), oracle_error=oracle_error))
        shifted = ProjectedBlochBasis(source, dict(kind="shift", window=window))
        # Exercise the actual nonorthogonal compiler plans as well as the
        # diagonal-H eigen calculation used for the longer plateau check.
        for word in (SHC_WORD, "Tr[1 delta(E-H)]"):
            left, _ = compute_projection_moments(shifted, word, moments=8, trace="exact")
            right, _ = compute_projection_moments(shifted.eigenbasis(), word, moments=8, trace="exact")
            assert np.max(np.abs(left.data-right.data)) < 1e-10
            if word != SHC_WORD:
                # Four selected states per cell; literal 1 must survive the
                # catalogue's usual simplification of identity operators.
                assert abs(float(np.asarray(left.data).ravel()[0].real)-4/source.num_wann) < 1e-9
        # Equal M for the wide and shifted spectra. Exact block trace is W
        # vectors, not one vector per (k,band); print oracle/KPM separately.
        kpm_values = []
        for tag, basis in (("unprojected", source.eigenbasis()), ("shift", shifted.eigenbasis())):
            print(f"(b) {route} {tag} exact-trace KPM M={args.moments} starting", flush=True)
            mu, _ = compute_moments(basis, SHC_WORD, moments=args.moments, trace="exact")
            energy, values = shc(mu, args.broadening_meV, source.num_cells, area, backend=args.backend)
            plateau = np.interp(fermi_grid, energy, values)
            kpm_values.append(plateau)
            print(f"(b) {route} {tag} KPM plateau={plateau.tolist()}, oracle={original_oracle.tolist()}", flush=True)
            assert np.max(np.abs(plateau-original_oracle)) < 1e-3
        assert np.max(np.abs(kpm_values[0]-kpm_values[1])) < 1e-3
        wide_range, shifted_range = np.ptp(source.eigenvalues), np.ptp(shifted.eigenvalues)
        print(f"(b) {route} spectral ranges {wide_range:.12g} -> {shifted_range:.12g} eV; "
              f"nominal fixed-M resolution gain={wide_range/shifted_range:.12g}; "
              f"pi*range/(2M)={np.pi*wide_range/(2*args.moments):.12g} -> "
              f"{np.pi*shifted_range/(2*args.moments):.12g} eV", flush=True)
    print("(b) oracle comparisons:", json.dumps(results, indent=2), flush=True)
    reduced = downfold_orthogonal(source.hk_array, source.sk_array, window)
    expected = np.stack([e[mask] for e, mask in zip(reduced.original_eigenvalues, reduced.kept)])
    assert np.max(np.abs(np.linalg.eigvalsh(reduced.hamiltonian)-expected)) < 1e-10
    assert np.max(np.abs(adjoint(reduced.isometry) @ source.sk_array @ reduced.isometry-np.eye(4))) < 1e-10
    local = LocalMatrix(reduced.hamiltonian, model.uc, model.tau[reduced.trial_columns],
                         source.supercell, source.k_points)
    band_check = validate_bands(local, lambda k: eigvalsh(model.Hk(k), model.Sk(k)),
                                window, model.uc, source.supercell)
    print("(c) orthogonal off-grid band validation:", band_check, flush=True)
    assert band_check["max_band_error_eV"] < 1e-6
    # Local per-k projector/physical operator export and rectangular mesh
    # exercise origin ordering, WS ties, gauge conversion and radius errors.
    small = BlochMetricBasis(model, [4, 3, 1], velocity_mode="covariant")
    projected = ProjectedBlochBasis(small, dict(kind="shift", window=window))
    local_h = projected.local_matrix("H")
    local_p = projected.local_matrix("1")  # physical S P; coefficient projector is S^-1 (S P)
    assert local_h.receipt()["interpolation_max_frobenius"] < 1e-9
    assert local_p.receipt()["interpolation_max_frobenius"] < 1e-9
    csr = projected.supercell_operators(["1", "SZ", "JXSZ"])
    for name, matrix in csr.items():
        assert np.max(np.abs((matrix-matrix.getH()).toarray())) < 1e-9, name
    print("local untruncated/truncated:", local_h.receipt(),
          LocalMatrix(projected.hk_array, model.uc, model.tau, small.supercell,
                       small.k_points, radius=0.).receipt(), flush=True)
    # The strict run-file round trip must retain the independent projection
    # window, rather than confuse it with the padded expansion window.
    raw = dict(model=dict(kind="siesta", dir="/tmp/siesta", label="toy"),
               basis=dict(name="eigen", supercell=[12, 12, 1]),
               metric=dict(policy="nonorthogonal", route="generalised", applier="exact"),
               projection=dict(kind="shift", window=window, radius=12.),
               word=SHC_WORD, expansion=dict(moments=256, window="eigen"),
               reconstruction=dict(broadening_meV=10.))
    rf = RunFile.from_dict(raw)
    assert RunFile.from_dict(rf.to_resolved_dict()).projection == rf.projection
    print("stage6d developer self-check reached all assertions; this is not a golden gate", flush=True)


if __name__ == "__main__":
    main()
