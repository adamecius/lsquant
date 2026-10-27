#!/usr/bin/env python3
"""Developer self-check for stage 3 (from_csr fix, rank-2 Kubo, G4 and G5).

Performs the stage-3 checks without a test framework:
  (0) Loads ``chain1d.HAM.CSR`` with ``scipy.sparse``, passes its ``indptr``,
      ``indices``, ``data`` (``int64`` / ``complex128``) to
      ``SparseMatrix.from_csr``, folds with ``hermitian_fold()``, and asserts
      the ``dos`` moments equal those from ``ReadCsrText`` (plus checks that
      non-``int64`` or non-contiguous inputs raise ``TypeError`` naming the
      offending argument).
  (a) Runs ``kubo_plan`` on Haldane ``VX, VY`` and ``VX, VX`` by exact trace
      and asserts ``numpy.array_equal`` against the committed G5 ``.chebmom2D``
      files.
  (b) Runs graphene random-phase trace with seed 12345 and
      ``window_from_gershgorin`` and asserts ``numpy.array_equal`` against G4.
  (c) Compares ``GreenRetarded.table`` and ``GreenRetardedDerivative.table``
      against ``legacy_green_tables`` and asserts max abs diff < 1e-12.
  (d) Runs the reconstructions of G5 and G4 and asserts 0 differing lines
      against all committed golden ``.dat`` files.
"""

from __future__ import annotations

import os
import sys
import tempfile
import numpy as np
import scipy.sparse as sp_sparse

from lsquant import (
    GreenRetarded,
    GreenRetardedDerivative,
    _core,
    dos_plan,
    kubo_bastin,
    kubo_bastin_sea,
    kubo_bastin_surf,
    kubo_fft_grid,
    kubo_greenwood,
    kubo_plan,
    legacy_green_tables,
    run_plan,
    window_from_gershgorin,
)


def _load_csr_as_scipy(path: str) -> sp_sparse.csr_matrix:
    """Parse an LSQUANT ASCII ``.CSR`` file into a ``scipy.sparse.csr_matrix``."""
    with open(path, encoding="utf-8") as fh:
        tokens = fh.read().split()
    dim = int(tokens[0])
    nnz = int(tokens[1])
    val_tokens = tokens[2 : 2 + 2 * nnz]
    col_tokens = tokens[2 + 2 * nnz : 2 + 3 * nnz]
    rp_tokens = tokens[2 + 3 * nnz : 2 + 3 * nnz + dim + 1]

    vals = np.empty(nnz, dtype=np.complex128)
    for i in range(nnz):
        vals[i] = complex(float(val_tokens[2 * i]), float(val_tokens[2 * i + 1]))
    cols = np.array([int(c) for c in col_tokens], dtype=np.int64)
    indptr = np.array([int(r) for r in rp_tokens], dtype=np.int64)
    return sp_sparse.csr_matrix((vals, cols, indptr), shape=(dim, dim))


def _roundtrip_legacy_2d(mu: _core.ChebyshevMoments) -> _core.ChebyshevMoments:
    """Write ``mu`` with ``WriteLegacyChebmom2D`` and read it back via ``ReadChebmom``."""
    with tempfile.NamedTemporaryFile(suffix=".chebmom2D", delete=False) as tmp:
        tmp_path = tmp.name
    try:
        _core.WriteLegacyChebmom2D(tmp_path, mu)
        return _core.ReadChebmom(tmp_path)
    finally:
        if os.path.exists(tmp_path):
            os.remove(tmp_path)


def _diff_lines(
    energies: np.ndarray, values: np.ndarray, gold_path: str
) -> tuple[int, int]:
    """Write ``(energies, values)`` with ``WriteCurveText`` and count line diffs vs ``gold_path``."""
    with tempfile.NamedTemporaryFile(mode="r", suffix=".dat", delete=False) as tmp:
        tmp_path = tmp.name
    try:
        _core.WriteCurveText(tmp_path, energies, values)
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


def main() -> int:
    repo_root = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))

    # (0) Item 0: SparseMatrix.from_csr from scipy CSR arrays without intermediate copies
    chain_csr_path = os.path.join(
        repo_root, "goldens", "stage0", "chain1d", "operators", "chain1d.HAM.CSR"
    )
    sp_ham = _load_csr_as_scipy(chain_csr_path)
    indptr_64 = np.ascontiguousarray(sp_ham.indptr, dtype=np.int64)
    indices_64 = np.ascontiguousarray(sp_ham.indices, dtype=np.int64)
    data_128 = np.ascontiguousarray(sp_ham.data, dtype=np.complex128)

    ham_from_scipy = _core.SparseMatrix.from_csr(
        sp_ham.shape[0], indptr_64, indices_64, data_128
    ).hermitian_fold()
    ham_from_text = _core.ReadCsrText(chain_csr_path)
    dim_chain = ham_from_text.dimension
    win_chain = _core.SpectralWindow(-2.0, 2.0)
    states_chain = _core.LocalStates(dim_chain, list(range(dim_chain)))
    mu_scipy = run_plan(
        dos_plan(64, dim_chain, word="1 delta"),
        ham_from_scipy,
        win_chain,
        states_chain,
        batch_size=64,
    )
    mu_text = run_plan(
        dos_plan(64, dim_chain, word="1 delta"),
        ham_from_text,
        win_chain,
        states_chain,
        batch_size=64,
    )
    csr_dos_equal = bool(np.array_equal(mu_scipy.data, mu_text.data))
    print(f"(0) SparseMatrix.from_csr (scipy) dos array_equal vs ReadCsrText: {csr_dos_equal}")
    assert csr_dos_equal, "SparseMatrix.from_csr dos moments differ from ReadCsrText"

    # Verify forcecast is off and TypeError names the offending argument
    for bad_arg, kwargs in (
        ("row_pointers", {"row_pointers": indptr_64.astype(np.int32)}),
        ("columns", {"columns": indices_64[::2]}),
        ("values", {"values": data_128.real}),
    ):
        call_args = {
            "dim": sp_ham.shape[0],
            "row_pointers": indptr_64,
            "columns": indices_64,
            "values": data_128,
        }
        call_args.update(kwargs)
        try:
            _core.SparseMatrix.from_csr(**call_args)
            raise AssertionError(f"Expected TypeError for bad {bad_arg}")
        except TypeError as exc:
            assert bad_arg in str(exc), f"TypeError did not name {bad_arg}: {exc}"

    # (a) kubo_plan on Haldane VX-VY and VX-VX by exact trace vs G5 .chebmom2D
    haldane_dir = os.path.join(repo_root, "goldens", "stage0", "haldane")
    haldane_op_dir = os.path.join(haldane_dir, "operators")
    h_ham = _core.ReadCsrText(os.path.join(haldane_op_dir, "haldane.HAM.CSR"))
    h_vx = _core.ReadCsrText(os.path.join(haldane_op_dir, "haldane.VX.CSR"))
    h_vy = _core.ReadCsrText(os.path.join(haldane_op_dir, "haldane.VY.CSR"))
    h_vx_op = _core.SparseMatrixOperator(h_vx)
    h_vy_op = _core.SparseMatrixOperator(h_vy)
    h_dim = h_ham.dimension
    h_win = _core.SpectralWindow(-4.5, 4.5)
    h_states = _core.LocalStates(h_dim, list(range(h_dim)))

    mu_h_xy = run_plan(
        kubo_plan(128, left_index=0, right_index=1, word="VX delta VY delta"),
        h_ham,
        h_win,
        h_states,
        operators=[h_vx_op, h_vy_op],
        batch_size=128,
        symmetrise=True,
    )
    mu_h_xx = run_plan(
        kubo_plan(128, left_index=0, right_index=0, word="VX delta VX delta"),
        h_ham,
        h_win,
        h_states,
        operators=[h_vx_op],
        batch_size=128,
        symmetrise=True,
    )

    ref_h_xy = _core.ReadChebmom(
        os.path.join(haldane_dir, "NonEqOpVX-VYhaldaneKPM_M128x128_stateexact.chebmom2D")
    )
    ref_h_xx = _core.ReadChebmom(
        os.path.join(haldane_dir, "NonEqOpVX-VXhaldaneKPM_M128x128_stateexact.chebmom2D")
    )
    g5_xy_equal = bool(np.array_equal(_roundtrip_legacy_2d(mu_h_xy).data, ref_h_xy.data))
    g5_xx_equal = bool(np.array_equal(_roundtrip_legacy_2d(mu_h_xx).data, ref_h_xx.data))
    print(f"(a) Haldane VX-VY array_equal vs G5: {g5_xy_equal}")
    print(f"(a) Haldane VX-VX array_equal vs G5: {g5_xx_equal}")
    assert g5_xy_equal, "Haldane VX-VY moments do not equal G5"
    assert g5_xx_equal, "Haldane VX-VX moments do not equal G5"

    # (b) Graphene random phase, seed 12345, Gershgorin window vs G4 .chebmom2D
    graphene_dir = os.path.join(repo_root, "goldens", "stage0", "graphene")
    graphene_op_dir = os.path.join(graphene_dir, "operators")
    g_ham = _core.ReadCsrText(os.path.join(graphene_op_dir, "graphene_golden.HAM.CSR"))
    g_vx = _core.ReadCsrText(os.path.join(graphene_op_dir, "graphene_golden.VX.CSR"))
    g_vx_op = _core.SparseMatrixOperator(g_vx)
    g_dim = g_ham.dimension
    g_win = window_from_gershgorin(g_ham)
    g_bw = (g_win.upper - g_win.lower) * 1.0
    g_states = _core.RandomPhaseStates(g_dim, 1, 12345)

    mu_g4 = run_plan(
        kubo_plan(64, left_index=0, right_index=0, word="VX delta VX delta"),
        g_ham,
        g_win,
        g_states,
        operators=[g_vx_op],
        batch_size=64,
        symmetrise=True,
    )
    ref_g4 = _core.ReadChebmom(
        os.path.join(
            graphene_dir, "NonEqOpVX-VXgraphene_goldenKPM_M64x64_statedefault.chebmom2D"
        )
    )
    g4_equal = bool(np.array_equal(_roundtrip_legacy_2d(mu_g4).data, ref_g4.data))
    print(f"(b) Graphene Gershgorin BandWidth = {g_bw:.13f}, array_equal vs G4: {g4_equal}")
    assert f"{g_bw:.13f}" == "15.0366847652085", f"Unexpected BandWidth {g_bw:.13f}"
    assert g4_equal, "Graphene VX-VX moments do not equal G4"

    # (c) SymPy GreenRetarded and GreenRetardedDerivative vs compat legacy_green_tables
    alpha_c = 0.80
    m_test = 16
    n_test = 30 * m_test
    x_grid = [-alpha_c + i * (2.0 * alpha_c) / (n_test - 1) for i in range(n_test)]
    sym_gr = GreenRetarded.table(x_grid, m_test)
    sym_dgr = GreenRetardedDerivative.table(x_grid, m_test)
    compat_gr, compat_dgr = legacy_green_tables(x_grid, m_test)
    diff_gr = float(np.max(np.abs(sym_gr - compat_gr)))
    diff_dgr = float(np.max(np.abs(sym_dgr - compat_dgr)))
    rel_dgr = float(np.max(np.abs(sym_dgr - compat_dgr) / np.abs(compat_dgr)))
    print(f"(c) GreenRetarded.table vs compat max_abs_diff           = {diff_gr:.6e}")
    print(
        f"(c) GreenRetardedDerivative.table vs compat max_abs_diff = {diff_dgr:.6e} "
        f"(max_rel_diff = {rel_dgr:.6e})"
    )
    assert diff_gr < 1e-12, f"GreenRetarded table diff {diff_gr} >= 1e-12"
    assert diff_dgr < 1e-12, f"GreenRetardedDerivative table diff {diff_dgr} >= 1e-12"

    # (d) Reconstructions of G5 and G4 vs golden .dat files read as text
    checks = [
        (
            "G5 Bastin (10 meV)",
            *kubo_bastin(ref_h_xy, 10.0),
            os.path.join(
                haldane_dir,
                "KuboBastin_NonEqOpVX-VYhaldaneKPM_M128x128_stateexactJACKSON.dat",
            ),
        ),
        (
            "G5 Sea (10 meV)",
            *kubo_bastin_sea(ref_h_xy, 10.0),
            os.path.join(
                haldane_dir,
                "KuboBastinSea_NonEqOpVX-VYhaldaneKPM_M128x128_stateexactJACKSON.dat",
            ),
        ),
        (
            "G5 Surf (10 meV)",
            *kubo_bastin_surf(ref_h_xy, 10.0),
            os.path.join(
                haldane_dir,
                "KuboBastinSurf_NonEqOpVX-VYhaldaneKPM_M128x128_stateexactJACKSON.dat",
            ),
        ),
        (
            "G5 Greenwood (10 meV)",
            *kubo_greenwood(ref_h_xx, 10.0),
            os.path.join(
                haldane_dir,
                "KuboGreenWood_NonEqOpVX-VXhaldaneKPM_M128x128_stateexactJACKSON.dat",
            ),
        ),
        (
            "G4 Bastin (64 meV)",
            *kubo_bastin(ref_g4, 64.0),
            os.path.join(
                graphene_dir,
                "KuboBastin_NonEqOpVX-VXgraphene_goldenKPM_M64x64_statedefaultJACKSON.dat",
            ),
        ),
        (
            "G4 Sea (10 meV)",
            *kubo_bastin_sea(ref_g4, 10.0),
            os.path.join(
                graphene_dir,
                "sea_surf",
                "KuboBastinSea_NonEqOpVX-VXgraphene_goldenKPM_M64x64_statedefaultJACKSON.dat",
            ),
        ),
        (
            "G4 Surf (10 meV)",
            *kubo_bastin_surf(ref_g4, 10.0),
            os.path.join(
                graphene_dir,
                "sea_surf",
                "KuboBastinSurf_NonEqOpVX-VXgraphene_goldenKPM_M64x64_statedefaultJACKSON.dat",
            ),
        ),
        (
            "G4 Bastin FFTgrid (10 meV)",
            *kubo_fft_grid(ref_g4, 10.0, "bastin"),
            os.path.join(
                graphene_dir,
                "fft_grid",
                "KuboBastin_NonEqOpVX-VXgraphene_goldenKPM_M64x64_statedefaultJACKSON_FFTgrid.dat",
            ),
        ),
        (
            "G4 Greenwood FFTgrid (10 meV)",
            *kubo_fft_grid(ref_g4, 10.0, "greenwood"),
            os.path.join(
                graphene_dir,
                "fft_grid",
                "KuboGreenWood_NonEqOpVX-VXgraphene_goldenKPM_M64x64_statedefaultJACKSON_FFTgrid.dat",
            ),
        ),
    ]

    for label, e_arr, v_arr, gold_path in checks:
        diffs, total = _diff_lines(e_arr, v_arr, gold_path)
        print(f"(d) {label:28s} differing lines: {diffs} (total lines: {total})")
        assert diffs == 0, f"{label}: expected 0 differing lines, got {diffs}"

    return 0


if __name__ == "__main__":
    sys.exit(main())
