#!/usr/bin/env python3
"""Developer self-check for stage 4 (time evolution, MSD, projected evolution, G6 & G7).

Performs the stage-4 checks without a test framework:
  (a) G7: runs ``projected_evolution_plan`` for ``SX``, ``SY``, ``SZ`` with
      projector ``PX`` on ``chain1d_mag`` and asserts ``numpy.array_equal``
      against the committed G7 ``.chebmomTD`` goldens.
  (b) G6: runs ``msd_plan`` on ``chain1d_ballistic`` ``VX`` and asserts
      ``numpy.array_equal`` against the committed G6 ``.chebmomTD`` golden, and
      reconstructs the four ``EF`` curves (``0.0, 0.5, 1.0, 1.5`` eV at ``10``
      meV) asserting 0 differing lines against the committed golden ``.dat``
      files.
  (c) Compares ``legacy_bessel_coefficients`` (GSL ``gsl_sf_bessel_Jn``) against
      ``bessel_series_coefficients`` (SciPy ``scipy.special.jv``) and asserts
      max abs diff < 1e-13.
  (d) Runs ``time_correlation_plan`` on ``chain1d`` ``VX, VX`` (``M = 64``,
      ``16`` steps, ``tmax = 5``, exact trace) and asserts ``numpy.array_equal``
      against the reference SHA-256 of the ``Tr[VX delta(E-H) U(t)^+ VX U(t)]``
      capture.
"""

from __future__ import annotations

import hashlib
import os
import sys
import tempfile
import numpy as np

from lsquant import (
    _core,
    bessel_series_coefficients,
    legacy_bessel_coefficients,
    msd_plan,
    projected_evolution_plan,
    run_plan,
    time_correlation_curve,
    time_correlation_plan,
    time_evolution_series,
    time_step_from_tmax,
    window_from_gershgorin,
)


def _roundtrip_legacy_td(mu: _core.ChebyshevMoments) -> _core.ChebyshevMoments:
    """Write ``mu`` with ``WriteLegacyChebmomTD`` and read it back via ``ReadChebmom``."""
    with tempfile.NamedTemporaryFile(suffix=".chebmomTD", delete=False) as tmp:
        tmp_path = tmp.name
    try:
        _core.WriteLegacyChebmomTD(tmp_path, mu)
        return _core.ReadChebmom(tmp_path)
    finally:
        if os.path.exists(tmp_path):
            os.remove(tmp_path)


def _diff_lines(
    times: np.ndarray, values: np.ndarray, gold_path: str
) -> tuple[int, int]:
    """Write ``(times, values)`` with ``WriteCurveText`` and count line diffs vs ``gold_path``."""
    with tempfile.NamedTemporaryFile(mode="r", suffix=".dat", delete=False) as tmp:
        tmp_path = tmp.name
    try:
        _core.WriteCurveText(tmp_path, times, values)
        with open(tmp_path, encoding="utf-8") as fh_new:
            new_lines = fh_new.readlines()
        with open(gold_path, encoding="utf-8") as fh_gold:
            gold_lines = fh_gold.readlines()
    finally:
        if os.path.exists(tmp_path):
            os.remove(tmp_path)

    diffs = abs(len(new_lines) - len(gold_lines)) + sum(
        1 for a, b in zip(new_lines, gold_lines) if a != b
    )
    return diffs, len(gold_lines)


# SHA-256 of TimeCorrVX-VXchain1dKPM_M64_stateexact.chebmomTD for
# Tr[VX delta(E-H) U(t)^+ VX U(t)] on chain1d (M=64, 16 steps, tmax=5, exact).
_LEGACY_CORR_SHA256 = (
    "84047151f40de0933e725fbfbf56627bc1e5b33f1bba7ec9671da1908e661445"
)


def main() -> int:
    repo_root = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))

    # (a) G7: TimeProj{SX,SY,SZ}-PX on chain1d_mag vs goldens/stage0/chain1d_mag
    g7_dir = os.path.join(repo_root, "goldens", "stage0", "chain1d_mag")
    g7_op_dir = os.path.join(g7_dir, "operators")
    mag_ham = _core.ReadCsrText(os.path.join(g7_op_dir, "chain1d_mag.HAM.CSR"))
    mag_px = _core.ReadCsrText(os.path.join(g7_op_dir, "chain1d_mag.PX.CSR"))
    mag_px_op = _core.SparseMatrixOperator(mag_px)
    mag_dim = mag_ham.dimension
    mag_win = window_from_gershgorin(mag_ham)
    mag_bw = (mag_win.upper - mag_win.lower) * 1.0
    mag_dt = time_step_from_tmax(96, 95.0)
    # Byte identity with the legacy captures requires the GSL Bessel values.
    mag_coeffs, mag_phase = time_evolution_series(
        mag_win, mag_dt, use_legacy_bessel=True
    )
    mag_states = _core.LocalStates(mag_dim, list(range(mag_dim)))

    print(
        f"(a) G7 chain1d_mag Gershgorin window = [{mag_win.lower:.2f}, {mag_win.upper:.2f}], "
        f"BandWidth = {mag_bw:.2f}"
    )
    assert f"{mag_bw:.2f}" == "4.62", f"Expected G7 BandWidth 4.62, got {mag_bw}"

    for op_name in ("SX", "SY", "SZ"):
        op_mat = _core.ReadCsrText(
            os.path.join(g7_op_dir, f"chain1d_mag.{op_name}.CSR")
        )
        op_lin = _core.SparseMatrixOperator(op_mat)
        plan_g7 = projected_evolution_plan(
            32,
            96,
            mag_dim,
            mag_coeffs,
            operator_index=0,
            projector_index=1,
            phase=mag_phase,
            word=f"PX U(t)^+ {op_name} U(t) delta PX",
        )
        mu_g7 = run_plan(
            plan_g7,
            mag_ham,
            mag_win,
            mag_states,
            operators=[op_lin, mag_px_op],
            batch_size=32,
            symmetrise=False,
            dt=mag_dt,
        )
        ref_g7 = _core.ReadChebmom(
            os.path.join(
                g7_dir,
                f"TimeProj{op_name}-PXchain1d_magKPM_M32_stateexact.chebmomTD",
            )
        )
        g7_equal = bool(
            np.array_equal(_roundtrip_legacy_td(mu_g7).data, ref_g7.data)
        )
        print(f"(a) G7 TimeProj{op_name}-PX array_equal vs golden: {g7_equal}")
        assert g7_equal, f"G7 TimeProj{op_name}-PX moments do not equal golden"

    # (b) G6: CorrelationVX on chain1d_ballistic vs goldens/stage0/chain1d_ballistic + 4 EF curves
    g6_dir = os.path.join(repo_root, "goldens", "stage0", "chain1d_ballistic")
    g6_op_dir = os.path.join(g6_dir, "operators")
    b_ham = _core.ReadCsrText(os.path.join(g6_op_dir, "chain1d.HAM.CSR"))
    b_vx = _core.ReadCsrText(os.path.join(g6_op_dir, "chain1d.VX.CSR"))
    b_dim = b_ham.dimension
    b_win = _core.SpectralWindow(-2.0, 2.0)
    b_vx_op = _core.RescaledOperator(b_vx, b_win)
    b_dt = time_step_from_tmax(64, 20.0)
    # Byte identity with the legacy captures requires the GSL Bessel values.
    b_coeffs, b_phase = time_evolution_series(
        b_win, b_dt, use_legacy_bessel=True
    )
    b_states = _core.LocalStates(b_dim, list(range(b_dim)))

    plan_g6 = msd_plan(
        256,
        64,
        b_dim,
        b_coeffs,
        velocity_index=0,
        phase=b_phase,
        word="[VX, U(t)]^+ delta [VX, U(t)]",
    )
    mu_g6 = run_plan(
        plan_g6,
        b_ham,
        b_win,
        b_states,
        operators=[b_vx_op],
        batch_size=256,
        symmetrise=False,
        dt=b_dt,
    )
    ref_g6 = _core.ReadChebmom(
        os.path.join(
            g6_dir, "CorrelationVXchain1dKPM_M256_state20exact.chebmomTD"
        )
    )
    mu_g6_disk = _roundtrip_legacy_td(mu_g6)
    g6_equal = bool(np.array_equal(mu_g6_disk.data, ref_g6.data))
    print(f"(b) G6 CorrelationVXchain1d array_equal vs golden: {g6_equal}")
    assert g6_equal, "G6 CorrelationVX moments do not equal golden"

    for ef in (0.0, 0.5, 1.0, 1.5):
        t_arr, v_arr = time_correlation_curve(mu_g6_disk, 10.0, ef)
        gold_dat = os.path.join(
            g6_dir,
            f"meanCorrelationVXchain1dKPM_M256_state20exactEF{ef:.6f}JACKSON.dat",
        )
        diffs, total = _diff_lines(t_arr, v_arr, gold_dat)
        print(
            f"(b) G6 EF={ef:+.1f} curve differing lines: {diffs} (total lines: {total})"
        )
        assert diffs == 0, f"G6 EF={ef} curve has {diffs} differing lines"

    # (c) compat GSL Bessel coefficients vs SciPy Bessel coefficients
    max_bessel_diff = 0.0
    for x_test in (0.5, 1.0, 3.5, 12.0, 25.0):
        c_gsl = legacy_bessel_coefficients(x_test)
        c_scipy = bessel_series_coefficients(x_test)
        assert c_gsl.shape == c_scipy.shape, (
            f"Bessel series length mismatch at x={x_test}: "
            f"{c_gsl.shape} vs {c_scipy.shape}"
        )
        diff_x = float(np.max(np.abs(c_gsl - c_scipy)))
        if diff_x > max_bessel_diff:
            max_bessel_diff = diff_x
    print(f"(c) compat GSL Bessel vs SciPy jv max_abs_diff = {max_bessel_diff:.6e}")
    assert max_bessel_diff < 1e-13, (
        f"Bessel coefficient diff {max_bessel_diff} >= 1e-13"
    )

    # (d) time_correlation_plan on chain1d VX, VX (M=64, 16 steps, tmax=5, exact)
    #     vs exact SHA-256 reference of Tr[VX delta(E-H) U(t)^+ VX U(t)]
    corr_dt = time_step_from_tmax(16, 5.0)
    # Byte identity with the legacy captures requires the GSL Bessel values.
    corr_coeffs, corr_phase = time_evolution_series(
        b_win, corr_dt, use_legacy_bessel=True
    )
    b_vx_static = _core.SparseMatrixOperator(b_vx)
    plan_corr = time_correlation_plan(
        64,
        16,
        b_dim,
        corr_coeffs,
        left_index=0,
        right_index=0,
        phase=corr_phase,
        word="VX delta U(t)^+ VX U(t)",
    )
    mu_corr = run_plan(
        plan_corr,
        b_ham,
        b_win,
        b_states,
        operators=[b_vx_static],
        batch_size=64,
        symmetrise=False,
        dt=corr_dt,
    )
    with tempfile.NamedTemporaryFile(suffix=".chebmomTD", delete=False) as tmp:
        tmp_corr_path = tmp.name
    try:
        _core.WriteLegacyChebmomTD(tmp_corr_path, mu_corr)
        with open(tmp_corr_path, "rb") as fh_corr:
            corr_sha256 = hashlib.sha256(fh_corr.read()).hexdigest()
    finally:
        if os.path.exists(tmp_corr_path):
            os.remove(tmp_corr_path)

    corr_equal = corr_sha256 == _LEGACY_CORR_SHA256
    print(
        f"(d) time_correlation_plan chain1d VX-VX (M=64, steps=16, tmax=5, exact) "
        f"array_equal vs reference SHA-256: {corr_equal}"
    )
    assert corr_equal, "time_correlation_plan moments do not match reference SHA-256"

    return 0


if __name__ == "__main__":
    sys.exit(main())
