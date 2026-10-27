#!/usr/bin/env python3
"""Developer self-check for stage 2 (pybind11 bindings, SymPy table, rank-1 DOS).

Performs three checks without a test framework:
  (a) Builds the stage-1 ``dos`` and ``kubo`` plans from ``lsquant.moments``,
      runs them on ``goldens/stage0/chain1d/operators`` through ``lsquant._core``,
      and asserts bit-exact array equality (``numpy.array_equal``) with
      ``ReadChebmom`` of the committed legacy goldens.
  (b) Compares ``DeltaFunction.table`` against ``legacy_delta_table`` on the
      uniform grid and prints the maximum absolute difference (expected < 1e-14).
  (c) Runs ``dos_curve`` on G1 at 120 meV, writes the formatted text curve, and
      prints the number of lines differing from the golden G3 file (expected 0).
"""

from __future__ import annotations

import os
import sys
import tempfile
import numpy as np

from lsquant import (
    DeltaFunction,
    _core,
    dos_curve,
    dos_plan,
    kubo_plan,
    legacy_delta_table,
    run_plan,
)


def main() -> int:
    repo_root = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
    op_dir = os.path.join(repo_root, "goldens", "stage0", "chain1d", "operators")
    g1_path = os.path.join(
        repo_root,
        "goldens",
        "stage0",
        "chain1d",
        "SpectralOp1chain1dKPM_M64_stateexact.chebmom1D",
    )
    g2_path = os.path.join(
        repo_root,
        "goldens",
        "stage0",
        "chain1d",
        "NonEqOpVX-VXchain1dKPM_M64x64_stateexact.chebmom2D",
    )
    g3_path = os.path.join(
        repo_root,
        "goldens",
        "stage0",
        "chain1d",
        "meanSpectralOp1chain1dKPM_M64_stateexactJACKSON.dat",
    )

    # (a) Stage-1 dos and kubo plans through _core vs ReadChebmom of goldens
    ham = _core.ReadCsrText(os.path.join(op_dir, "chain1d.HAM.CSR"))
    vx = _core.ReadCsrText(os.path.join(op_dir, "chain1d.VX.CSR"))
    vx_op = _core.SparseMatrixOperator(vx)
    dim = ham.dimension
    window = _core.SpectralWindow(-2.0, 2.0)
    states = _core.LocalStates(dim, list(range(dim)))
    num_moments = 64

    mu_dos = run_plan(
        dos_plan(num_moments, dim, word="1 delta"),
        ham,
        window,
        states,
        batch_size=num_moments,
    )
    mu_kubo = run_plan(
        kubo_plan(num_moments, left_index=0, right_index=0, word="VX delta VX delta"),
        ham,
        window,
        states,
        operators=[vx_op],
        batch_size=num_moments,
        symmetrise=True,
    )

    ref_dos = _core.ReadChebmom(g1_path)
    ref_kubo = _core.ReadChebmom(g2_path)

    dos_equal = bool(np.array_equal(mu_dos.data, ref_dos.data))
    kubo_equal = bool(np.array_equal(mu_kubo.data, ref_kubo.data))
    print(f"(a) dos_plan array_equal vs G1:  {dos_equal}")
    print(f"(a) kubo_plan array_equal vs G2: {kubo_equal}")
    assert dos_equal, "dos_plan tensor does not equal G1"
    assert kubo_equal, "kubo_plan tensor does not equal G2"

    # (b) DeltaFunction.table vs legacy_delta_table max abs diff
    alpha = 0.95
    m_eff = 53
    num_div = 30 * m_eff
    x_grid = [-alpha + i * (2.0 * alpha) / (num_div - 1) for i in range(num_div)]
    sym_table = DeltaFunction.table(x_grid, m_eff)
    leg_table = legacy_delta_table(x_grid, m_eff)
    max_table_diff = float(np.max(np.abs(sym_table - leg_table)))
    print(f"(b) DeltaFunction.table vs legacy_delta_table max_abs_diff = {max_table_diff:.6e}")
    assert max_table_diff < 1e-14, f"table difference {max_table_diff} >= 1e-14"

    # (c) dos_curve on G1 at 120 meV vs G3 file read as text
    energies, values = dos_curve(ref_dos, 120.0)
    with tempfile.NamedTemporaryFile(mode="r", suffix=".dat", delete=False) as tmp:
        tmp_path = tmp.name
    try:
        _core.WriteCurveText(tmp_path, energies, values)
        with open(tmp_path, encoding="utf-8") as fh_new:
            new_lines = fh_new.readlines()
        with open(g3_path, encoding="utf-8") as fh_gold:
            gold_lines = fh_gold.readlines()
    finally:
        if os.path.exists(tmp_path):
            os.remove(tmp_path)

    diff_lines = abs(len(new_lines) - len(gold_lines)) + sum(
        1 for a, b in zip(new_lines, gold_lines) if a != b
    )
    print(f"(c) dos_curve vs G3 differing lines: {diff_lines} (total lines: {len(gold_lines)})")
    assert diff_lines == 0, f"expected 0 differing lines against G3, got {diff_lines}"

    return 0


if __name__ == "__main__":
    sys.exit(main())
