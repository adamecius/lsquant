#!/usr/bin/env python3
"""The single stage-6c developer self-check, exactly brief section 1 (a)-(c).

On Aron, after building the extension:
  PYTHONPATH=python:build:build/lsquant_core python tests/dev/stage_6c_check.py

No SIESTA/material data needed. This does not run the Bi gate or certify a
stage. Mesh=12x12, W=1 eV and disorder seed=7 are fixed, not CLI knobs.
"""
from __future__ import annotations

import os
from pathlib import Path
import sys
from types import SimpleNamespace
import numpy as np


def compare(label, actual, expected, tolerance, *, absolute=False):
    """Relative max error uses a tensor norm, so symmetry zeros stay meaningful."""
    scale = max(float(np.max(np.abs(expected))), np.finfo(float).tiny)
    error = float(np.max(np.abs(actual-expected)))
    print(f"{label}: max_abs={error:.17g} relative_max={error/scale:.17g} tolerance={tolerance:g}", flush=True)
    assert np.all(np.isfinite(actual)) and (error if absolute else error/scale) <= tolerance, label


def fixtures():
    """Use precisely 6b's rotated Kane-Mele fixture, including its full dipole."""
    from stage_6b_check import kane_mele_arrays
    from lsquant.models.siesta import SiestaModel, midpoint_connection
    from lsquant.models.nonorthogonal import adjoint
    lattice, centres, r_int, H, S = kane_mele_arrays()
    zero = int(np.flatnonzero(np.all(r_int == 0, axis=1))[0])
    T = np.kron(np.diag([1.5, .6]) @
                np.array([[np.cos(.37), np.sin(.37)], [-np.sin(.37), np.cos(.37)]]), np.eye(2))
    Hp, Sp = adjoint(T) @ H @ T, adjoint(T) @ S @ T
    Q = np.zeros((3, *H.shape), complex)
    for axis in range(3):
        Q[axis, zero] = adjoint(T) @ np.diag(centres[:, axis]) @ T
    A = midpoint_connection(Q, Sp, centres, r_int @ lattice)
    original = SiestaModel(lattice, centres, r_int, H, S, np.zeros_like(A), "KM")
    rotated = SiestaModel(lattice, centres, r_int, Hp, Sp, A, "KM_rotated")
    return original, rotated


def main():
    root = Path(__file__).resolve().parents[2]
    build = Path(os.environ.get("LSQUANT_BUILD_DIR", root/"build"))
    for path in (root/"python", build, build/"lsquant_core"):
        sys.path.insert(0, str(path))
    from lsquant import _core
    from lsquant.metric_forms import (anderson_forms, disorder_velocity_forms,
                                     periodic_position_velocity_proof)
    from lsquant.bases import BlochBasis, BlochMetricBasis, make_anderson_disorder
    from lsquant.disorder_bases import DisorderedBlochMetricBasis, RealSpaceMetricBasis
    from lsquant.apps.stage6b_siesta import compute_moments, DOS_WORD, SHC_WORD
    from lsquant.models.disorder import anderson_potential
    from lsquant.moments import dos_plan, run_plan

    symbolic = anderson_forms()
    print("Anderson sympy derivation:", symbolic, flush=True)
    assert symbolic["constant_residual"] == 0
    assert symbolic["orthogonal_residual"] == 0
    assert symbolic["hermiticity"]["physical"]["hermitian"]
    position_proof = periodic_position_velocity_proof()
    print("periodic position proof:", position_proof, flush=True)
    assert position_proof["residual"] == 0
    for mode in ("bare", "berry_connection", "covariant"):
        derivation = disorder_velocity_forms(mode)
        print("disorder velocity", mode, derivation, flush=True)
        assert derivation["substitution_residual"] == 0 and derivation["hermitian"]

    original, model = fixtures()
    mesh, M, R, phase_seed = [12, 12, 1], 64, 2, 12345
    rng = np.random.default_rng(31)
    for route in ("lowdin", "generalised"):
        for width in (0., 1.):
            tag = "(a) W=0" if width == 0 else "(b) W=1 eV seed=7 12x12"
            bloch = DisorderedBlochMetricBasis(model, mesh, disorder_w=width, disorder_seed=7,
                                               velocity_mode="covariant", route=route)
            real = RealSpaceMetricBasis(model, mesh, disorder=bloch.disorder_values,
                                         disorder_w=width, disorder_seed=7,
                                         velocity_mode="covariant", route=route, applier="exact")
            window = bloch.window()
            rb, rr = bloch.states("random", R, phase_seed), real.states("random", R, phase_seed)
            # These are the streams used by compute_moments, not independent
            # NumPy phases or separately sampled bra/ket vectors.
            rb.reset()
            rr.reset()
            for _ in range(R):
                b, r = rb.next(), rr.next()
                compare(tag+" mapped phases "+route, bloch.transform.multiply(r), b, 1e-12)
            clean = (BlochMetricBasis(model, mesh, velocity_mode="covariant", route=route)
                     if width == 0 else None)
            for word in (DOS_WORD, SHC_WORD):
                options = dict(moments=M, states=R, seed=phase_seed, window=window, batch=17)
                mu_b, _ = compute_moments(bloch, word, **options)
                mu_r, _ = compute_moments(real, word, **options)
                compare(tag+" real/Bloch "+route+" "+word, mu_r.data, mu_b.data,
                        1e-10 if width == 0 else 1e-9, absolute=width == 0)
                if clean is not None:
                    mu_6b, _ = compute_moments(clean, word, **options)
                    compare("(a) clean Bloch/6b "+route+" "+word, mu_b.data, mu_6b.data, 1e-12)
            # (c) Hermiticity and the explicit Fourier image of the same V.
            V = real.potential_csr
            anti = V-V.conj().T
            anti_norm = float(np.max(np.abs(anti.data))) if anti.nnz else 0.
            print(f"(c) V Hermitian {route} W={width}: max_abs={anti_norm:.17g}", flush=True)
            assert anti_norm < 1e-13
            for _ in range(2):
                x = rng.normal(size=bloch.dimension)+1j*rng.normal(size=bloch.dimension)
                y = rng.normal(size=bloch.dimension)+1j*rng.normal(size=bloch.dimension)
                bx, by = bloch.potential().multiply(x), bloch.potential().multiply(y)
                hermitian_residual = abs(np.vdot(x, by)-np.vdot(bx, y))
                assert hermitian_residual <= 1e-12*max(np.linalg.norm(x)*np.linalg.norm(by), 1.)
                reference = bloch.transform.multiply(real.potential().multiply(bloch.transform.adjoint().multiply(x)))
                if width > 0:
                    compare("(c) FFT/CSR potential "+route, bx, reference, 1e-12)
                else:
                    assert np.array_equal(bx, np.zeros_like(bx))

    # (c) Orthogonal limit: identical RNG and identical stage-5 D, including
    # bytes. Give both engines the SAME orbital-Bloch phase stream and the
    # SAME window: stage 5's default B*r_real is a different sample than the
    # 6b orbital-Bloch sampler required here, although both are unbiased.
    orth = DisorderedBlochMetricBasis(original, mesh, disorder_w=1., disorder_seed=7,
                                      velocity_mode="covariant", route="lowdin")
    from scipy.sparse import eye
    w = make_anderson_disorder(orth.dimension, 1., 7)
    assert orth.disorder_values.tobytes() == w.tobytes()
    csr = anderson_potential(eye(orth.dimension, format="csr"), w)
    assert csr.diagonal().tobytes() == w.tobytes()
    assert csr.nnz == orth.dimension
    stage5_model = SimpleNamespace(num_wann=original.num_wann, uc=original.uc, tau=original.tau,
                                   Hk=original.Hk, k_mesh=original.k_mesh,
                                   bloch_phases=lambda k: np.exp(-1j*(k @ original.tau.T)))
    old = BlochBasis(stage5_model, mesh, disorder_w=1., disorder_seed=7)
    window = orth.window()
    # Compare actions as well as moments: an average of two half-potentials
    # can otherwise change bytes when it is accumulated onto the clean H.
    x = rng.normal(size=orth.dimension)+1j*rng.normal(size=orth.dimension)
    actual, expected = orth.potential().multiply(x), old._sandwich_v.multiply(x)
    assert actual.tobytes() == expected.tobytes()
    phase_generator = _core.RandomPhaseStates(orth.dimension, R, phase_seed)
    reference = run_plan(dos_plan(M, R, DOS_WORD), old.hamiltonian(), window,
                          phase_generator, batch_size=17)
    actual, _ = compute_moments(orth, DOS_WORD, moments=M, states=R, seed=phase_seed,
                                window=window, batch=17)
    print(f"(c) orthogonal stage-5 DOS bytes equal: {actual.data.tobytes() == reference.data.tobytes()}", flush=True)
    assert actual.data.tobytes() == reference.data.tobytes()
    print("stage_6c_check: (a)-(c) passed; Bi gate and thread-count checks remain with the tester", flush=True)


if __name__ == "__main__":
    main()
