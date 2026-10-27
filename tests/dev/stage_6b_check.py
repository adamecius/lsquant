#!/usr/bin/env python3
"""The single stage-6b developer self-check (brief section 8), for Aron.

Run after building the extension:
  PYTHONPATH=python:build/lsquant_core python tests/dev/stage_6b_check.py
No SIESTA, Wannier, or material files are needed. This is not a golden gate.
"""
from __future__ import annotations

import argparse
import os
from pathlib import Path
import sys
import numpy as np


def kane_mele_arrays():
    """Return honeycomb KM hoppings, t1=-1, lambda_SO=0.15, no Rashba.

    ``H_AB=-sum_(d=0,-a1,-a2) exp(ik.d) I_spin``;
    ``H_AA=2 lambda [sin(q1)-sin(q2)+sin(q2-q1)] sigma_z``;
    ``H_BB=-H_AA``. Cell rows (sqrt(3)/2,3/2),(-sqrt(3)/2,3/2),
    bond length 1 Angstrom. Spin is interleaved, Pauli units. The orientation
    is pinned by the independent -2 Im curvature and FHS computations, with
    the positive spin sector carrying C=+1 (not by taking abs of a result).
    """
    lattice = np.array([[np.sqrt(3)/2, 1.5, 0], [-np.sqrt(3)/2, 1.5, 0], [0, 0, 10.]])
    centres = np.repeat(np.array([[0., 0., 0.], [0., 1., 0.]]), 2, axis=0)
    hoppings = {}
    def add(r, block):
        r = tuple(r)
        hoppings[r] = hoppings.get(r, np.zeros((4, 4), complex)) + block
    for r in [(0, 0, 0), (-1, 0, 0), (0, -1, 0)]:
        block = np.zeros((4, 4), complex)
        block[:2, 2:] = -np.eye(2)
        add(r, block)
        add(tuple(-x for x in r), block.conj().T)
    for r, orientation in [((1, 0, 0), 1), ((0, 1, 0), -1), ((-1, 1, 0), 1)]:
        # i*h*exp(iq)+h.c.=-2*h*sin(q).
        block = -1j*0.15*orientation*np.diag([1, -1, -1, 1])
        add(r, block)
        add(tuple(-x for x in r), block.conj().T)
    vectors = sorted(hoppings)
    H = np.array([hoppings[r] for r in vectors])
    S = np.zeros_like(H)
    S[vectors.index((0, 0, 0))] = np.eye(4)
    return lattice, centres, np.array(vectors), H, S


def check_one_vector(original, rotated):
    """(g) Compare the actual engine stream with one NumPy quadratic form.

    Inspect raw moments BEFORE any trace-cyclic Hermitian projection, then
    check the unsymmetrised KPM weights separately. No independent random
    generator is used for the NumPy sample: engine.run resets this stream.
    """
    from lsquant import _core
    from lsquant.apps.stage6b_siesta import prepare_moments, compute_moments, SHC_WORD
    from lsquant.bases import BlochMetricBasis
    from lsquant.models.nonorthogonal import adjoint

    M = 8
    for model in (original, rotated):
        for route in ("generalised", "lowdin"):
            bloch = BlochMetricBasis(model, [4, 4, 1], velocity_mode="covariant", route=route)
            Y, C = bloch.inverse_sqrt_array, bloch.eigenvectors
            U = bloch.sk_array @ Y @ C
            A = adjoint(C) @ bloch.operator_blocks("JXSZ") @ C
            B = adjoint(C) @ bloch.operator_blocks("VY") @ C
            for basis in (bloch, bloch.eigenbasis()):
                prepared = prepare_moments(basis, SHC_WORD, moments=M, states=1, seed=12345,
                                           batch=3)
                engine, plan, _, generator, window, _, _ = prepared
                generator.reset()
                r = generator.next()  # Same vector that run() receives after its Reset().
                assert r is not None and generator.next() is None
                r = r.reshape(bloch.num_cells, bloch.num_wann)
                z = ((adjoint(U) @ r[..., None])[..., 0] if basis.name == "bloch" else r)
                x = window.scale*bloch.eigenvalues + window.shift
                polynomials = np.empty((M, *x.shape))
                polynomials[0], polynomials[1] = 1., x
                for m in range(2, M):
                    polynomials[m] = 2*x*polynomials[m-1]-polynomials[m-2]
                direct = np.empty((M, M), complex)
                for m in range(M):
                    for n in range(M):
                        right = (B @ (polynomials[n]*z)[..., None])[..., 0]
                        word_r = (A @ (polynomials[m]*right)[..., None])[..., 0]
                        direct[m, n] = np.vdot(z, word_r)
                mu = engine.run(plan)
                residual = float(np.max(np.abs(mu.data-direct)))
                print(f"(g) {model.seed_name} {route} {basis.name} one-vector quadratic-form residual={residual:.17g}", flush=True)
                assert residual < 1e-10
                weights = np.full(M, 2.)
                weights[0] = 1.
                _core.ApplyKpmWeights(mu, 1, False)
                assert np.max(np.abs(mu.data-direct*weights[:, None]*weights[None, :])) < 1e-10
                projected = 0.5*(direct+direct.conj().T)*weights[:, None]*weights[None, :]
                production, _ = compute_moments(basis, SHC_WORD, _prepared=prepared)
                assert production.symmetrised
                assert np.max(np.abs(production.data-projected)) < 1e-10


def main():
    """Run (a)-(g), exact rank-2 comparisons and both phase noise diagnostics."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--noise-samples", type=int, default=20)
    parser.add_argument("--noise-states", type=int, default=16)
    parser.add_argument("--noise-mesh", type=int, default=60)
    parser.add_argument("--noise-moments", type=int, default=512)
    args = parser.parse_args()
    if args.noise_samples < 2 or min(args.noise_states, args.noise_mesh, args.noise_moments) < 1:
        parser.error("noise requires at least two samples and positive states/mesh/moments")
    root = Path(__file__).resolve().parents[2]
    build = Path(os.environ.get("LSQUANT_BUILD_DIR", root/"build"))
    for path in (root/"python", build, build/"lsquant_core"):
        sys.path.insert(0, str(path))
    import json
    import scipy.sparse as sparse
    from lsquant.metric_forms import describe, trace_equivalence
    print(json.dumps(describe(), indent=2))
    proof = trace_equivalence()
    print("symbolic trace proofs:", proof)
    assert proof["hamiltonian_similarity"]
    for name, item in proof.items():
        if name not in ("hamiltonian_similarity", "anticommutator"):
            assert item["equal"], (name, item)
    assert proof["anticommutator"]["single_operator_equal"]
    assert not proof["anticommutator"]["generic_products_equal"]
    from lsquant.models.bloch import bloch_sum
    from lsquant.models.siesta import SiestaModel, midpoint_connection
    from lsquant.models.nonorthogonal import adjoint, rotate_basis, bloch_velocity, spin_current
    from lsquant.bases import BlochMetricBasis
    from lsquant.oracle.spin_hall import spin_berry_curvature, sigma_xy, spin_chern, eigenphase_trace
    from lsquant.apps.stage6b_siesta import compute_moments, SHC_WORD, DOS_WORD
    from lsquant.reconstruct import shc

    lattice, centres, r_int, H, S = kane_mele_arrays()
    zero = int(np.flatnonzero(np.all(r_int == 0, axis=1))[0])
    r_cart = r_int @ lattice
    # A fixed local, spin-independent, nonsymmetric T. Its singular values
    # are 1.5 and 0.6, hence cond(S')=6.25, in the stage-6 fixture range.
    T = np.kron(np.array([[1.5, 0.], [0., 0.6]]) @
                np.array([[np.cos(.37), np.sin(.37)], [-np.sin(.37), np.cos(.37)]]), np.eye(2))
    rotated_H, rotated_S = adjoint(T) @ H @ T, adjoint(T) @ S @ T
    Q = np.zeros((3, *H.shape), complex)
    for axis in range(3):
        Q[axis, zero] = adjoint(T) @ np.diag(centres[:, axis]) @ T
    A = midpoint_connection(Q, rotated_S, centres, r_cart)
    original = SiestaModel(lattice, centres, r_int, H, S, np.zeros_like(A), "KM")
    rotated = SiestaModel(lattice, centres, r_int, rotated_H, rotated_S, A, "KM_rotated")
    check_one_vector(original, rotated)
    rng = np.random.default_rng(12345)
    points = rng.normal(size=(5, 3)) * np.array([1, 1, 0])
    maximum_covariance, maximum_hf, maximum_hermiticity = 0., 0., 0.
    for k in points:
        phase = np.diag(np.exp(1j*(centres @ k)))
        Tk = adjoint(phase) @ T @ phase
        Hp, Sp = rotated.Hk(k), rotated.Sk(k)
        assert np.max(np.abs(Hp-adjoint(Tk) @ original.Hk(k) @ Tk)) < 1e-12
        assert np.max(np.abs(Sp-adjoint(Tk) @ original.Sk(k) @ Tk)) < 1e-12
        # Use the stage-6 rotate_basis implementation on Hermitian k blocks;
        # applying its fold separately to non-Hermitian H(R) would be wrong.
        old = rotate_basis({"H": sparse.csr_matrix(bloch_sum(r_cart, centres, H, k, cell=True))}, sparse.csr_matrix(T))
        assert np.max(np.abs(old["H"].toarray()-bloch_sum(r_cart, centres, rotated_H, k, cell=True))) < 1e-12
        assert np.max(np.abs(old["S"].toarray()-bloch_sum(r_cart, centres, rotated_S, k, cell=True))) < 1e-12
        from lsquant.oracle.spin_hall import generalised_eigh
        E, C = generalised_eigh(Hp, Sp)
        X = np.linalg.inv(Sp)
        for axis in (0, 1):
            dH, dS = rotated.derivatives(k, axis)
            berry = bloch_velocity(Hp, Sp, dH, dS, mode="berry_connection", inverse_overlap=X)
            connection = rotated.Ak(k, axis)
            covariant = bloch_velocity(Hp, Sp, dH, dS, connection, "covariant", inverse_overlap=X)
            dipole_term = 1j*(Hp @ X @ connection-adjoint(connection) @ X @ Hp)
            reference = adjoint(Tk) @ original.derivatives(k, axis)[0] @ Tk
            error = float(np.max(np.abs(covariant-reference)))
            maximum_covariance = max(maximum_covariance, error)
            assert error < 1e-10
            assert np.max(np.abs(covariant-berry-dipole_term)) < 1e-12
            print(f"dipole axis={axis} norm={np.linalg.norm(dipole_term):.17g}")
            for mode in ("bare", "berry_connection", "covariant"):
                dh, ds = rotated.derivatives(k, axis, cell=mode == "bare")
                velocity = bloch_velocity(Hp, Sp, dh, ds, connection if mode == "covariant" else None,
                                          mode, inverse_overlap=X)
                J = spin_current(velocity, np.diag([1, -1, 1, -1]))
                herm = max(np.max(np.abs(velocity-adjoint(velocity))), np.max(np.abs(J-adjoint(J))))
                maximum_hermiticity = max(maximum_hermiticity, float(herm))
                assert herm < 1e-13
            # Sorted bands remain spin-degenerate; both members have the same
            # intraband derivative. Never compare a gauge-dependent interband
            # element to a finite difference of eigenvalues.
            step = np.eye(3)[axis]*1e-5
            plus = generalised_eigh(rotated.Hk(k+step), rotated.Sk(k+step))[0]
            minus = generalised_eigh(rotated.Hk(k-step), rotated.Sk(k-step))[0]
            diagonal = np.diag(adjoint(C) @ covariant @ C).real
            hf = float(np.max(np.abs(diagonal-(plus-minus)/(2e-5))))
            maximum_hf = max(maximum_hf, hf)
            assert hf < 1e-6
    try:
        bloch_velocity(Hp, Sp, dH, dS, mode="covariant")
    except ValueError:
        pass
    else:
        raise AssertionError("covariant velocity accepted missing dipoles")
    print(f"(a)-(d) covariance={maximum_covariance:.17g} Hermiticity={maximum_hermiticity:.17g} HF={maximum_hf:.17g}")

    area = float(np.linalg.norm(np.cross(lattice[0], lattice[1])))
    # Dense mesh for the continuum integer check; the deliberately short
    # 12x12 KPM below is compared to its OWN finite-mesh oracle.
    dense = BlochMetricBasis(rotated, [96, 96, 1], velocity_mode="covariant")
    curvature = spin_berry_curvature(dense.eigenvalues, dense.eigenvectors,
                                     dense.operator_blocks("JXSZ"), dense.operator_blocks("VY"))
    plateau = float(sigma_xy(dense.eigenvalues, curvature, [0.], area)[0])
    topology = spin_chern(dense.eigenvalues, dense.eigenvectors, dense.sk_array,
                          dense.operator_blocks("SZ"), 0., dense.supercell,
                          k_points=dense.k_points, centres=centres,
                          sqrt_overlap=dense.sk_array @ dense.inverse_sqrt_array)
    print(f"(e) sigma_xy={plateau:.12f} spin_chern={topology['spin_chern']:.12f} Z2={topology['z2']}")
    Y = dense.inverse_sqrt_array
    U = Y @ dense.sk_array @ dense.physical_eigenvectors
    lowdin_curvature = spin_berry_curvature(
        dense.physical_eigenvalues, U, Y @ dense.operator_blocks("JXSZ") @ Y,
        Y @ dense.operator_blocks("VY") @ Y)
    lowdin_plateau = float(sigma_xy(dense.physical_eigenvalues, lowdin_curvature, [0.], area)[0])
    print(f"(e) explicit Lowdin oracle plateau={lowdin_plateau:.12f}")
    assert abs(lowdin_plateau-2.) < 1e-6
    assert abs(lowdin_plateau-plateau) < 1e-10
    assert abs(plateau-2.) < 1e-6
    assert abs(topology["spin_chern"]-2.) < 1e-6
    diagonal = np.zeros_like(dense.hk_array)
    idx = np.arange(4)
    diagonal[:, idx, idx] = curvature*(dense.eigenvalues < 0)
    random_trace = eigenphase_trace(diagonal, states=16, seed=12345).real
    exact_trace = float(np.sum(curvature*(dense.eigenvalues < 0)))
    print(f"oracle eigenphase-minus-exact={random_trace-exact_trace:.17g}")
    assert abs(random_trace-exact_trace) < 1e-10

    small = BlochMetricBasis(rotated, [12, 12, 1], velocity_mode="covariant", route="generalised")
    small_curvature = spin_berry_curvature(small.eigenvalues, small.eigenvectors,
                                           small.operator_blocks("JXSZ"), small.operator_blocks("VY"))
    small_oracle = float(sigma_xy(small.eigenvalues, small_curvature, [0.], area)[0])
    exact, _ = compute_moments(small.eigenbasis(), SHC_WORD, moments=256, trace="exact")
    bloch_exact, _ = compute_moments(small, SHC_WORD, moments=256, trace="exact")
    exact_difference = float(np.max(np.abs(bloch_exact.data-exact.data)))
    print(f"SHC exact block Bloch-minus-eigen moments={exact_difference:.17g}")
    assert exact_difference < 1e-10
    energies, curve = shc(exact, 10., small.num_cells, area)
    value = float(np.interp(0., energies, curve))
    print(f"(f) generalised exact KPM M=256 Nk=12x12 plateau={value:.12f} oracle={small_oracle:.12f}")
    assert abs(value-small_oracle) < 1e-2
    random, _ = compute_moments(small, SHC_WORD, moments=256, states=16, seed=12345)
    eigen_random, _ = compute_moments(small.eigenbasis(), SHC_WORD, moments=256, states=16, seed=12345)
    residual = float(np.max(np.abs(random.data-eigen_random.data)))
    print(f"SHC orbital-phase same-sample Bloch-minus-eigen moments={residual:.17g}; random-minus-exact={np.max(np.abs(random.data-exact.data)):.17g} (diagnostic only)")
    lowdin = BlochMetricBasis(rotated, [12, 12, 1], velocity_mode="covariant", route="lowdin")
    lowdin_exact, _ = compute_moments(lowdin, SHC_WORD, moments=256, trace="exact")
    assert np.max(np.abs(lowdin_exact.data-exact.data)) < 1e-10
    el, vl = shc(lowdin_exact, 10., small.num_cells, area)
    print(f"(f) Lowdin exact KPM M=256 Nk=12x12 plateau={np.interp(0.,el,vl):.12f}")
    assert abs(float(np.interp(0., el, vl))-small_oracle) < 1e-2
    approximate = BlochMetricBasis(rotated, [12, 12, 1], velocity_mode="covariant", route="anticommutator")
    approx_mu, _ = compute_moments(approximate.eigenbasis(), SHC_WORD, moments=256, trace="exact")
    ea, va = shc(approx_mu, 10., approximate.num_cells, area)
    print(f"anticommutator exact plateau={np.interp(0.,ea,va):.12f}; equality is NOT required by the derivation")

    # A direct eigenphase still traces a diagonal DOS exactly. This diagnostic
    # deliberately bypasses the production orbital-phase state adapter.
    from lsquant import _core
    from lsquant.moments import dos_plan, run_plan
    eigen_basis = small.eigenbasis()
    dos_eigenphase = run_plan(dos_plan(64, 1, DOS_WORD), eigen_basis.hamiltonian(),
                              small.window(), _core.RandomPhaseStates(small.dimension, 1, 12345))
    dos_exact, _ = compute_moments(eigen_basis, DOS_WORD, moments=64, trace="exact")
    dos_difference = float(np.max(np.abs(dos_eigenphase.data-dos_exact.data)))
    print(f"DOS single eigenphase-minus-exact={dos_difference:.17g}")
    assert dos_difference < 1e-12

    # Review 2: same ORTHOGONAL KM model/settings as the tester's stage-3 run.
    # Both estimates include the same Hermitian projection and reconstruction.
    # The reference spread below is quoted evidence, not a new local measurement.
    gate_mesh = [args.noise_mesh, args.noise_mesh, 1]
    noise_basis = BlochMetricBasis(original, gate_mesh, velocity_mode="covariant", route="lowdin")
    eigen_noise_basis = noise_basis.eigenbasis()
    reference, _ = compute_moments(eigen_noise_basis, SHC_WORD,
                                    moments=args.noise_moments, trace="exact")
    en, vn = shc(reference, 30., noise_basis.num_cells, area)
    exact_plateau = float(np.interp(0., en, vn))
    from lsquant.apps.stage6b_siesta import prepare_moments
    stage3_std = 0.0041
    print("stage-3 real-space baseline (tester review-2 addendum, 10 seeds): "
          "orthogonal KM mesh=60x60 M=512 broadening_meV=30 states=16 "
          f"mean=2.0027 std={stage3_std:.4f}; seed=12345 plateau=1.9957", flush=True)
    for phase_mode in ("eigenbasis", "orbital-Bloch"):
        samples = []
        for sample in range(args.noise_samples):
            # Direct phases in eigen coordinates versus phases in orthonormal
            # orbital coordinates. Never confuse U^+r with direct eigen phases.
            basis = eigen_noise_basis if phase_mode == "eigenbasis" else noise_basis
            prepared = prepare_moments(basis, SHC_WORD, moments=args.noise_moments,
                                       trace="random", states=args.noise_states, seed=12345+sample)
            if phase_mode == "eigenbasis":
                prepared = list(prepared)
                prepared[3] = _core.RandomPhaseStates(basis.dimension, args.noise_states, 12345+sample)
                prepared[0] = _core.MomentEngine(prepared[6], prepared[5], prepared[3], args.noise_moments)
            mu, _ = compute_moments(basis, SHC_WORD, _prepared=prepared)
            en, vn = shc(mu, 30., noise_basis.num_cells, area)
            samples.append(float(np.interp(0., en, vn)))
            print(f"SHC {phase_mode} noise sample={sample+1}/{args.noise_samples} seed={12345+sample} plateau={samples[-1]:.12f}", flush=True)
        deviation = float(np.std(samples, ddof=1))
        print(f"SHC {phase_mode} noise orthogonal_KM mesh={args.noise_mesh}x{args.noise_mesh} M={args.noise_moments} broadening_meV=30 states={args.noise_states} seeds={args.noise_samples} mean={np.mean(samples):.12f} std={deviation:.12f} standard_error={deviation/np.sqrt(args.noise_samples):.12f} exact_KPM={exact_plateau:.12f} stage3_real_space_std={stage3_std:.4f} (quoted 10-seed baseline) [e^2/h, sigma_z Pauli]", flush=True)
    print("stage_6b_check: developer self-check passed (not a golden gate)")


if __name__ == "__main__":
    main()
