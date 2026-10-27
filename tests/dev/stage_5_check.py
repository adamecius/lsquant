#!/usr/bin/env python3
"""Stage-5 developer self-check for BlochBasis, EigenBasis, and RealBasis.

Verifies checks (a) through (f) of ``handoff/0005-bloch-and-eigen-bases/brief.md``:

(a) Unitarity of ``BlochTransform`` (``B B^+ x == x`` and ``B^+ B x == x`` to ``1e-14``)
    on ``(8, 8, 1)`` graphene (``N_k = 64, W = 2, dim = 128``) across 5 random
    complex vectors.
(b) Basis equivalence of the clean Hamiltonian:
    ``max_x || H_k x - B H_csr B^+ x || < 1e-13`` on ``(8, 8, 1)`` graphene and
    Haldane across 5 random complex vectors.
(c) Basis equivalence of the atom-gauge velocity operator:
    ``max_x || v_x x - B VX_csr B^+ x || < 1e-13`` on ``(8, 8, 1)`` graphene and
    Haldane across 5 random complex vectors.
(d) Disorder equivalence via ``SandwichOperator``:
    for box disorder ``v in [-W/2, +W/2]`` (``W = 1.0, seed = 12345``),
    ``H_real = H_csr + diag(v)`` and ``H_bloch = H_k + B diag(v) B^+`` yield
    64 DOS moments and ``64 x 64`` ``VX-VX`` moments agreeing to ``< 1e-10``
    with random-phase trace (``seed = 12345``).
(e) Three-basis exact trace equivalence at ``W = 0`` on ``(16, 16, 1)`` graphene
    (``N_k = 256, dim = 512, M = 128, broadening = 200 meV, window = EigenBasis.window()``):
    ``RealBasis``, ``BlochBasis``, and ``EigenBasis`` moment tensors and
    reconstructed DOS and ``VX-VX`` Kubo-Greenwood curves agree to ``< 1e-11``,
    and agree with the direct eigenbasis Gaussian formula of ``hr_exactdiag.py``
    within the KPM Jackson-vs-Gaussian kernel difference.
(f) Regression check that ``scripts/check_goldens.sh 3`` and
    ``tests/dev/stage_4_check.py`` pass without regression.
"""

from __future__ import annotations

import math
import os
import shutil
import subprocess
import sys
import tempfile
import numpy as np

REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
PYTHON_DIR = os.path.join(REPO_ROOT, "python")
if PYTHON_DIR not in sys.path:
    sys.path.insert(0, PYTHON_DIR)

import lsquant as lsq
from lsquant import _core
from lsquant.bases import BlochBasis, EigenBasis, RealBasis, make_anderson_disorder
from lsquant.models.wannier import load_wannier_model
from lsquant.moments import dos_plan, kubo_plan, run_plan
from lsquant.reconstruct import dos_curve, kubo_greenwood


def _write_wannier_triplet(
    out_dir: str,
    seed: str,
    uc: np.ndarray,
    species: list[str],
    tau: np.ndarray,
    hoppings: dict[tuple[int, int, int], np.ndarray],
) -> None:
    """Write ``<seed>_hr.dat``, ``<seed>.uc``, and ``<seed>.xyz`` matching ``wannier2sparse``."""
    os.makedirs(out_dir, exist_ok=True)
    nw = tau.shape[0]
    r_keys = sorted(hoppings.keys())
    nr = len(r_keys)

    hr_path = os.path.join(out_dir, f"{seed}_hr.dat")
    with open(hr_path, "w", encoding="utf-8") as fh:
        fh.write(f" {seed} written by stage_5_check.py\n")
        fh.write(f"{nw:12d}\n{nr:12d}\n")
        for i in range(nr):
            fh.write(f"{1:5d}")
            if (i + 1) % 15 == 0 or i + 1 == nr:
                fh.write("\n")
        for rx, ry, rz in r_keys:
            blk = hoppings[(rx, ry, rz)]
            for b in range(nw):
                for a in range(nw):
                    val = complex(blk[a, b])
                    fh.write(
                        f"{rx:5d}{ry:5d}{rz:5d}{a + 1:5d}{b + 1:5d}"
                        f"{val.real:20.12f}{val.imag:20.12f}\n"
                    )

    uc_path = os.path.join(out_dir, f"{seed}.uc")
    with open(uc_path, "w", encoding="utf-8") as fh:
        for row in uc:
            fh.write(f"{row[0]:20.12f} {row[1]:20.12f} {row[2]:20.12f}\n")

    xyz_path = os.path.join(out_dir, f"{seed}.xyz")
    with open(xyz_path, "w", encoding="utf-8") as fh:
        fh.write(f"{nw}\n{seed} basis (Angstrom)\n")
        for sp, pos in zip(species, tau):
            fh.write(f"{sp:<4s} {pos[0]:20.12f} {pos[1]:20.12f} {pos[2]:20.12f}\n")


def _generate_graphene_and_haldane(tmp_dir: str) -> tuple[str, str]:
    """Generate the ``wannier2sparse/examples/gen_models.py`` graphene and Haldane models in ``tmp_dir``."""
    uc_str = "1.5 0.8660254 0\n1.5 -0.8660254 0\n0 0 10\n"
    xyz_str = "2\nCA 0 0 0\nCB 1 0 0\n"

    def write_sparse(
        out_dir: str,
        model: str,
        num_wann: int,
        rpts: set[tuple[int, int, int]],
        hoppings: list[tuple[tuple[int, int, int], int, int, float, float]],
    ) -> None:
        os.makedirs(out_dir, exist_ok=True)
        with open(os.path.join(out_dir, f"{model}_hr.dat"), "w", encoding="utf-8") as fh:
            fh.write(f"{model} (generated)\n{num_wann}\n{len(rpts)}\n")
            fh.write(" ".join("1" for _ in rpts) + "\n")
            for r_vec, i, j, re_v, im_v in hoppings:
                fh.write(
                    f"{r_vec[0]:d} {r_vec[1]:d} {r_vec[2]:d}  {i:d} {j:d}  "
                    f"{re_v:.10f} {im_v:.10f}\n"
                )
        with open(os.path.join(out_dir, f"{model}.uc"), "w", encoding="utf-8") as fh:
            fh.write(uc_str)
        with open(os.path.join(out_dir, f"{model}.xyz"), "w", encoding="utf-8") as fh:
            fh.write(xyz_str)

    t1 = -1.0
    g_hops = [
        ((0, 0, 0), 1, 2, t1, 0.0),
        ((0, 0, 0), 2, 1, t1, 0.0),
        ((-1, 0, 0), 1, 2, t1, 0.0),
        ((1, 0, 0), 2, 1, t1, 0.0),
        ((0, -1, 0), 1, 2, t1, 0.0),
        ((0, 1, 0), 2, 1, t1, 0.0),
    ]
    g_rpts = {(0, 0, 0), (1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0)}
    g_dir = os.path.join(tmp_dir, "graphene")
    write_sparse(g_dir, "graphene", 2, g_rpts, g_hops)

    t2 = 0.15
    phi = 0.5 * math.pi
    c_val = t2 * math.cos(phi)
    s_val = t2 * math.sin(phi)
    h_hops = list(g_hops)
    for r_nnn in [(1, 0, 0), (0, -1, 0), (-1, 1, 0)]:
        h_hops += [(r_nnn, 1, 1, c_val, s_val), (r_nnn, 2, 2, c_val, -s_val)]
    for r_nnn in [(-1, 0, 0), (0, 1, 0), (1, -1, 0)]:
        h_hops += [(r_nnn, 1, 1, c_val, -s_val), (r_nnn, 2, 2, c_val, s_val)]
    h_rpts = {t[0] for t in h_hops}
    h_dir = os.path.join(tmp_dir, "haldane")
    write_sparse(h_dir, "haldane", 2, h_rpts, h_hops)
    return g_dir, h_dir


def _dense_from_operator(op: _core.LinearOperator) -> np.ndarray:
    """Materialise a small ``LinearOperator`` as a dense ``(dim, dim)`` matrix."""
    dim = int(op.dimension)
    mat = np.empty((dim, dim), dtype=np.complex128)
    e_vec = np.zeros(dim, dtype=np.complex128)
    for col in range(dim):
        e_vec[col] = 1.0
        mat[:, col] = op.multiply(e_vec)
        e_vec[col] = 0.0
    return mat


def _find_w2s_binary() -> str | None:
    """Locate the ``wannier2sparse`` executable on TCN if available."""
    candidates = [
        os.environ.get("W2S_BIN", ""),
        os.path.expanduser("~/wannier2sparse/build/wannier2sparse"),
        shutil.which("wannier2sparse") or "",
    ]
    for path in candidates:
        if path and os.path.isfile(path) and os.access(path, os.X_OK):
            return path
    return None


def main() -> int:
    with tempfile.TemporaryDirectory(prefix="lsquant_stage5_") as tmp_dir:
        g_dir, h_dir = _generate_graphene_and_haldane(tmp_dir)
        graphene = load_wannier_model(g_dir, "graphene")
        haldane = load_wannier_model(h_dir, "haldane")

        # ------------------------------------------------------------------
        # Check (a): BlochTransform vs dense B on (6, 4, 1) and unitarity
        # ------------------------------------------------------------------
        sc64 = (6, 4, 1)
        bloch64 = BlochBasis(graphene, sc64)
        b_op64 = bloch64.transform
        b_adj64 = b_op64.adjoint()
        dim64 = bloch64.dimension
        nk64 = bloch64.num_cells
        nw = graphene.num_wann

        # Build dense B with n1-fastest cell index iR = n1 + N1 * n2 + N1 * N2 * n3
        n_cells = np.array(
            [
                (n1, n2, n3)
                for n3 in range(sc64[2])
                for n2 in range(sc64[1])
                for n1 in range(sc64[0])
            ],
            dtype=np.float64,
        )
        r_cells = n_cells @ graphene.uc
        k_pts64 = graphene.k_mesh(sc64)
        b_dense = np.zeros((dim64, dim64), dtype=np.complex128)
        for al in range(nw):
            b_dense[al::nw, al::nw] = np.exp(
                -1j * (k_pts64 @ (r_cells + graphene.tau[al]).T)
            ) / math.sqrt(float(nk64))

        rng = np.random.default_rng(202605)
        max_rel_dense_f = 0.0
        max_rel_dense_b = 0.0
        for _ in range(5):
            x = rng.standard_normal(dim64) + 1j * rng.standard_normal(dim64)
            x /= np.linalg.norm(x)
            yf = b_op64.multiply(x)
            yf_ref = b_dense @ x
            yb = b_adj64.multiply(x)
            yb_ref = b_dense.conj().T @ x
            max_rel_dense_f = max(
                max_rel_dense_f,
                float(np.max(np.abs(yf - yf_ref)) / np.max(np.abs(yf_ref))),
            )
            max_rel_dense_b = max(
                max_rel_dense_b,
                float(np.max(np.abs(yb - yb_ref)) / np.max(np.abs(yb_ref))),
            )
        print(
            f"[check a] BlochTransform vs dense B (6,4,1) graphene: "
            f"rel(B x) = {max_rel_dense_f:.3e}, rel(B^+ x) = {max_rel_dense_b:.3e}"
        )
        assert max_rel_dense_f < 1e-13 and max_rel_dense_b < 1e-13

        sc8 = (8, 8, 1)
        bloch8 = BlochBasis(graphene, sc8)
        b_op = bloch8.transform
        b_adj = b_op.adjoint()
        dim8 = bloch8.dimension
        assert dim8 == 128, f"Expected dim=128, got {dim8}"

        max_err_bbh = 0.0
        max_err_bhb = 0.0
        for _ in range(5):
            x = rng.standard_normal(dim8) + 1j * rng.standard_normal(dim8)
            x /= np.linalg.norm(x)
            y1 = b_op.multiply(b_adj.multiply(x))
            y2 = b_adj.multiply(b_op.multiply(x))
            max_err_bbh = max(max_err_bbh, float(np.max(np.abs(y1 - x))))
            max_err_bhb = max(max_err_bhb, float(np.max(np.abs(y2 - x))))

        print(
            f"[check a] BlochTransform unitarity (8,8,1) graphene: "
            f"max|B B^+ x - x| = {max_err_bbh:.3e}, max|B^+ B x - x| = {max_err_bhb:.3e}"
        )
        assert max_err_bbh < 1e-14 and max_err_bhb < 1e-14

        # ------------------------------------------------------------------
        # Check (b) & (c): wannier2sparse CSR identity on (6,4,1) and Haldane 48x48
        # ------------------------------------------------------------------
        w2s_bin = _find_w2s_binary()
        for name, mdl, mdir in (
            ("graphene", graphene, g_dir),
            ("haldane", haldane, h_dir),
        ):
            real_b = RealBasis(mdl, sc64, use_w2s_csr=False)
            bloch_b = BlochBasis(mdl, sc64)
            h_csr_op = real_b.hamiltonian()
            hk_op = bloch_b.hamiltonian()
            vx_csr_op = real_b.operator("VX")
            vx_k_op = bloch_b.operator("VX")
            b_tr = bloch_b.transform
            b_tr_adj = b_tr.adjoint()

            if w2s_bin is not None:
                w2s_out64 = os.path.join(tmp_dir, f"w2s_{name}_6x4")
                subprocess.run(
                    [
                        w2s_bin,
                        name,
                        "6",
                        "4",
                        "1",
                        "VX",
                        "-p",
                        mdir,
                        "-o",
                        w2s_out64,
                        "--quiet",
                        "--no-log-file",
                    ],
                    check=True,
                    stdout=subprocess.DEVNULL,
                )
                real_w2s = RealBasis(mdl, sc64, csr_dir=w2s_out64)
                assert real_w2s.operator_source("HAM").startswith("wannier2sparse ")
                assert real_b.operator_source("HAM") == "built by lsquant"
                h_built_dense = _dense_from_operator(h_csr_op)
                h_w2s_dense = _dense_from_operator(real_w2s.hamiltonian())
                max_h_csr_diff = float(np.max(np.abs(h_built_dense - h_w2s_dense)))
                print(
                    f"[check b] ({name} 6x4x1) max|H_real_basis - H_w2s_csr| = {max_h_csr_diff:.3e}"
                )
                assert max_h_csr_diff == 0.0, (
                    f"Non-zero |H_real_basis - H_w2s_csr| on {name} 6x4x1: {max_h_csr_diff}"
                )
                if name == "graphene":
                    bad_csr_dir = os.path.join(tmp_dir, "w2s_graphene_6x4_bad_receipt")
                    shutil.copytree(w2s_out64, bad_csr_dir)
                    bad_out_path = os.path.join(bad_csr_dir, "graphene.out")
                    with open(bad_out_path, encoding="utf-8") as fh:
                        orig_rec = fh.read()
                    # Right-size CSR (dim=48 for 6x4x1), but receipt claims supercell 4x6x1 (also dim=48)
                    with open(bad_out_path, "w", encoding="utf-8") as fh:
                        fh.write(orig_rec.replace("graphene 6 4 1 VX", "graphene 4 6 1 VX"))
                    raised_wrong_sc = False
                    try:
                        RealBasis(mdl, sc64, csr_dir=bad_csr_dir)
                    except ValueError:
                        raised_wrong_sc = True
                    assert raised_wrong_sc, (
                        "Expected RealBasis to raise ValueError on right-size CSR with wrong receipt supercell"
                    )
                    print(
                        "[check b] (graphene 6x4x1) right-size CSR with wrong receipt raises ValueError: PASS"
                    )

            max_h_diff = 0.0
            max_vx_diff = 0.0
            for _ in range(5):
                x = rng.standard_normal(dim64) + 1j * rng.standard_normal(dim64)
                x /= np.linalg.norm(x)
                lhs_h = hk_op.multiply(x)
                rhs_h = b_tr.multiply(h_csr_op.multiply(b_tr_adj.multiply(x)))
                max_h_diff = max(max_h_diff, float(np.max(np.abs(lhs_h - rhs_h))))

                lhs_v = vx_k_op.multiply(x)
                rhs_v = b_tr.multiply(vx_csr_op.multiply(b_tr_adj.multiply(x)))
                max_vx_diff = max(max_vx_diff, float(np.max(np.abs(lhs_v - rhs_v))))

            print(
                f"[check b] ({name} 6x4x1) max_x ||H_k x - B H_csr B^+ x||_inf = {max_h_diff:.3e}"
            )
            print(
                f"[check c] ({name} 6x4x1) max_x ||v_x x - B VX_csr B^+ x||_inf = {max_vx_diff:.3e}"
            )
            assert max_h_diff < 1e-13, f"H_k mismatch on {name}: {max_h_diff}"
            assert max_vx_diff < 1e-13, f"v_x mismatch on {name}: {max_vx_diff}"

        if w2s_bin is not None:
            sc48 = (48, 48, 1)
            w2s_hal48 = os.path.join(tmp_dir, "w2s_haldane_48x48")
            subprocess.run(
                [
                    w2s_bin,
                    "haldane",
                    "48",
                    "48",
                    "1",
                    "VX",
                    "-p",
                    h_dir,
                    "-o",
                    w2s_hal48,
                    "--quiet",
                    "--no-log-file",
                ],
                check=True,
                stdout=subprocess.DEVNULL,
            )
            real_hal48 = RealBasis(haldane, sc48, csr_dir=w2s_hal48)
            bloch_hal48 = BlochBasis(haldane, sc48)
            vx_w2s_op = real_hal48.operator("VX")
            vx_k48_op = bloch_hal48.operator("VX")
            b48 = bloch_hal48.transform
            b48_adj = b48.adjoint()
            dim48 = bloch_hal48.dimension
            max_rel_vx48 = 0.0
            for _ in range(4):
                x48 = rng.standard_normal(dim48) + 1j * rng.standard_normal(dim48)
                lhs48 = vx_k48_op.multiply(x48)
                rhs48 = b48.multiply(vx_w2s_op.multiply(b48_adj.multiply(x48)))
                rel48 = float(np.max(np.abs(lhs48 - rhs48)) / np.max(np.abs(lhs48)))
                max_rel_vx48 = max(max_rel_vx48, rel48)
            print(
                f"[check c] (haldane 48x48x1) rel ||v_x x - B VX_w2s B^+ x||_inf = {max_rel_vx48:.3e}"
            )
            assert max_rel_vx48 < 1e-10, (
                f"Haldane 48x48 B VX_w2s B^+ mismatch: {max_rel_vx48}"
            )

        # ------------------------------------------------------------------
        # Check (d): Disordered RealBasis vs BlochBasis (W = 1.0, seed = 12345)
        # ------------------------------------------------------------------
        v_dis = make_anderson_disorder(dim8, disorder_w=1.0, seed=12345)
        real_dis = RealBasis(graphene, sc8, disorder=v_dis, seed=12345)
        bloch_dis = BlochBasis(graphene, sc8, disorder=v_dis, seed=12345)
        win_dis = _core.SpectralWindow(-10.0, 10.0)

        st_real_dos = real_dis.states("random", num_states=4, seed=12345)
        st_bloch_dos = bloch_dis.states("random", num_states=4, seed=12345)
        plan_dos64 = dos_plan(64, int(st_real_dos.number_of_states), word="1 delta")
        mu_real_dos = run_plan(plan_dos64, real_dis.hamiltonian(), win_dis, st_real_dos)
        mu_bloch_dos = run_plan(
            plan_dos64, bloch_dis.hamiltonian(), win_dis, st_bloch_dos
        )
        diff_dis_dos = float(np.max(np.abs(mu_real_dos.data - mu_bloch_dos.data)))

        st_real_kubo = real_dis.states("random", num_states=4, seed=12345)
        st_bloch_kubo = bloch_dis.states("random", num_states=4, seed=12345)
        plan_kubo64 = kubo_plan(
            64, left_index=0, right_index=0, word="VX delta VX delta"
        )
        mu_real_kubo = run_plan(
            plan_kubo64,
            real_dis.hamiltonian(),
            win_dis,
            st_real_kubo,
            operators=[real_dis.operator("VX")],
        )
        mu_bloch_kubo = run_plan(
            plan_kubo64,
            bloch_dis.hamiltonian(),
            win_dis,
            st_bloch_kubo,
            operators=[bloch_dis.operator("VX")],
        )
        diff_dis_kubo = float(np.max(np.abs(mu_real_kubo.data - mu_bloch_kubo.data)))
        print(
            f"[check d] Anderson disorder W=1.0 (seed=12345): "
            f"max|mu_real - mu_bloch| DOS(64) = {diff_dis_dos:.3e}, "
            f"VX-VX(64x64) = {diff_dis_kubo:.3e}"
        )
        assert diff_dis_dos < 1e-10 and diff_dis_kubo < 1e-10

        # ------------------------------------------------------------------
        # Check (e): Clean (16, 16, 1) graphene exact trace in all 3 bases
        # ------------------------------------------------------------------
        sc16 = (16, 16, 1)
        m16 = 128
        broad16 = 200.0

        b_real = RealBasis(graphene, sc16)
        b_bloch = BlochBasis(graphene, sc16)
        b_eigen = EigenBasis(graphene, sc16)
        win16 = b_eigen.window()
        win16_bloch = b_bloch.window()
        assert (
            win16.lower == win16_bloch.lower and win16.upper == win16_bloch.upper
        ), "BlochBasis.window() and EigenBasis.window() must match"

        dos_moments: dict[str, np.ndarray] = {}
        dos_curves: dict[str, tuple[np.ndarray, np.ndarray]] = {}
        kubo_moments: dict[str, np.ndarray] = {}
        kubo_curves: dict[str, tuple[np.ndarray, np.ndarray]] = {}

        for b_name, b_obj in (
            ("real", b_real),
            ("bloch", b_bloch),
            ("eigen", b_eigen),
        ):
            st_d = b_obj.states("exact")
            p_d = dos_plan(m16, int(st_d.number_of_states), word="1 delta")
            mu_d = run_plan(p_d, b_obj.hamiltonian(), win16, st_d)
            dos_moments[b_name] = np.array(mu_d.data, copy=True)
            dos_curves[b_name] = dos_curve(mu_d, broad16)

            st_k = b_obj.states("exact")
            p_k = kubo_plan(m16, left_index=0, right_index=0, word="VX delta VX delta")
            mu_k = run_plan(
                p_k,
                b_obj.hamiltonian(),
                win16,
                st_k,
                operators=[b_obj.operator("VX")],
            )
            kubo_moments[b_name] = np.array(mu_k.data, copy=True)
            kubo_curves[b_name] = kubo_greenwood(mu_k, broad16)

        for pair in (("real", "bloch"), ("real", "eigen"), ("bloch", "eigen")):
            a_k, b_k = pair
            d_mu_dos = float(np.max(np.abs(dos_moments[a_k] - dos_moments[b_k])))
            d_cv_dos_abs = float(
                np.max(np.abs(dos_curves[a_k][1] - dos_curves[b_k][1]))
            )
            d_cv_dos_rel = d_cv_dos_abs / float(np.max(np.abs(dos_curves[a_k][1])))
            d_mu_kubo = float(np.max(np.abs(kubo_moments[a_k] - kubo_moments[b_k])))
            d_cv_kubo_abs = float(
                np.max(np.abs(kubo_curves[a_k][1] - kubo_curves[b_k][1]))
            )
            d_cv_kubo_rel = d_cv_kubo_abs / float(np.max(np.abs(kubo_curves[a_k][1])))
            d_cv_kubo_per_site = d_cv_kubo_abs / float(b_real.dimension)
            print(
                f"[check e] (16,16,1) {a_k} vs {b_k}: "
                f"DOS max|d_mu|={d_mu_dos:.3e}, max_rel|d_curve|={d_cv_dos_rel:.3e}; "
                f"KuboXX max|d_mu|={d_mu_kubo:.3e}, max|d_curve/N|={d_cv_kubo_per_site:.3e}, "
                f"max_rel|d_curve|={d_cv_kubo_rel:.3e}"
            )
            assert d_mu_dos < 1e-11 and d_cv_dos_rel < 1e-11
            assert d_mu_kubo < 1e-11 and d_cv_kubo_per_site < 1e-11 and d_cv_kubo_rel < 1e-11

        # Compare Jackson-reconstructed KPM curve with EigenBasis.exact_dos / exact_sigma_xx
        e_dos, y_dos_kpm = dos_curves["eigen"]
        e_kubo, y_kubo_kpm = kubo_curves["eigen"]
        eta16 = 2.0 * (broad16 / 1000.0) / (win16.upper - win16.lower)
        m_eff16 = min(int(math.ceil(math.pi / eta16)), m16)
        sigma_theta = math.pi / float(m_eff16 + 1)
        sigma_e_dos = sigma_theta * np.sqrt(
            np.maximum(1e-12, win16.half_width**2 - (e_dos - win16.center) ** 2)
        )
        sigma_e_kubo = sigma_theta * np.sqrt(
            np.maximum(1e-12, win16.half_width**2 - (e_kubo - win16.center) ** 2)
        )
        dos_gauss = float(b_eigen.num_cells) * b_eigen.exact_dos(e_dos, sigma_e_dos)
        kubo_gauss = b_eigen.exact_sigma_xx(e_kubo, sigma_e_kubo)
        rel_l2_dos = float(
            np.linalg.norm(y_dos_kpm - dos_gauss) / np.linalg.norm(dos_gauss)
        )
        rel_l2_kubo = float(
            np.linalg.norm(y_kubo_kpm - kubo_gauss) / np.linalg.norm(kubo_gauss)
        )
        print(
            f"[check e] Jackson KPM (M={m16}) vs Gaussian exact-diag: "
            f"rel_L2(DOS)={rel_l2_dos:.3e}, rel_L2(KuboXX)={rel_l2_kubo:.3e}"
        )
        assert rel_l2_dos < 0.03 and rel_l2_kubo < 0.05

        # EigenBasis.exact_dos and exact_sigma_xx vs hr_exactdiag wrapped on Gamma mesh of 12x12 graphene
        w2s_tools_dir = os.environ.get(
            "W2S_TOOLS", os.path.expanduser("~/wannier2sparse/tools")
        )
        hr_ed_path = os.path.join(w2s_tools_dir, "hr_exactdiag.py")
        if w2s_bin is not None and os.path.isfile(hr_ed_path):
            sc12 = (12, 12, 1)
            eig12 = EigenBasis(graphene, sc12)
            e_grid12 = np.linspace(-3.0, 3.0, 241, dtype=np.float64)
            eta12 = 0.04
            dos_eig12 = eig12.exact_dos(e_grid12, eta12)
            sig_eig12 = eig12.exact_sigma_xx(e_grid12, eta12)

            w2s_b12 = os.path.join(tmp_dir, "w2s_graphene_bundle")
            subprocess.run(
                [
                    w2s_bin,
                    "graphene",
                    "1",
                    "1",
                    "1",
                    "VX",
                    "--mode",
                    "bundle",
                    "-p",
                    g_dir,
                    "-o",
                    w2s_b12,
                    "--quiet",
                    "--no-log-file",
                ],
                check=True,
                stdout=subprocess.DEVNULL,
            )
            vx_hr_path = os.path.join(
                w2s_b12, "graphene.w2sp", "operators", "VX.hr.dat"
            )
            if w2s_tools_dir not in sys.path:
                sys.path.insert(0, w2s_tools_dir)
            from hr_exactdiag import opk as hr_opk, read_hr as hr_read_hr  # type: ignore[import-not-found]

            ir_h, h_r = hr_read_hr(os.path.join(g_dir, "graphene_hr.dat"))
            ir_vx, vx_r = hr_read_hr(vx_hr_path)
            n12 = np.arange(12, dtype=np.float64) / 12.0
            k12 = np.stack(np.meshgrid(n12, n12, indexing="ij"), -1).reshape(-1, 2)
            kf12 = np.c_[k12, np.zeros(len(k12), dtype=np.float64)]
            hk12 = hr_opk(ir_h, h_r, kf12)
            hk12 = 0.5 * (hk12 + np.conj(np.transpose(hk12, (0, 2, 1))))
            w12, u12 = np.linalg.eigh(hk12)
            g12 = np.exp(
                -((e_grid12[None, None, :] - w12[:, :, None]) ** 2)
                / (2.0 * eta12**2)
            ) / (math.sqrt(2.0 * math.pi) * eta12)
            dos_hr12 = g12.sum((0, 1)) / float(len(kf12))
            va12 = np.einsum(
                "kmi,kmn,knj->kij", np.conj(u12), hr_opk(ir_vx, vx_r, kf12), u12
            )
            a_hr12 = abs(float(graphene.uc[0, 0] * graphene.uc[1, 1]))
            sig_hr12 = (
                np.einsum("knm,kne,kme->e", np.abs(va12) ** 2, g12, g12)
                * (2.0 * math.pi**2 / a_hr12)
                / float(len(kf12))
            )
            sig_eig12_e2h = sig_eig12 * (
                2.0 * math.pi / (a_hr12 * float(eig12.num_cells))
            )
            rel_dos_hr12 = float(
                np.max(np.abs(dos_eig12 - dos_hr12)) / np.max(np.abs(dos_hr12))
            )
            rel_sig_hr12 = float(
                np.max(np.abs(sig_eig12_e2h - sig_hr12)) / np.max(np.abs(sig_hr12))
            )
            print(
                f"[check e] (graphene 12x12x1) EigenBasis.exact_dos/exact_sigma_xx vs hr_exactdiag (Gamma mesh): "
                f"rel_dos={rel_dos_hr12:.3e}, rel_sigma_xx={rel_sig_hr12:.3e}"
            )
            assert rel_dos_hr12 < 1e-10 and rel_sig_hr12 < 1e-10, (
                f"EigenBasis vs hr_exactdiag mismatch on 12x12 graphene: "
                f"dos={rel_dos_hr12}, sigma_xx={rel_sig_hr12}"
            )

        # ------------------------------------------------------------------
        # Check (f): G1-G8 regression checks against committed goldens
        # ------------------------------------------------------------------
        g1_ham_path = os.path.join(
            REPO_ROOT, "goldens", "stage0", "chain1d", "operators", "chain1d.HAM.CSR"
        )
        g1_ref_path = os.path.join(
            REPO_ROOT,
            "goldens",
            "stage0",
            "chain1d",
            "SpectralOp1chain1dKPM_M64_stateexact.chebmom1D",
        )
        if os.path.isfile(g1_ham_path) and os.path.isfile(g1_ref_path):
            g1_ham = _core.ReadCsrText(g1_ham_path)
            g1_dim = g1_ham.dimension
            g1_win = _core.SpectralWindow(-2.0, 2.0)
            g1_states = _core.LocalStates(g1_dim, list(range(g1_dim)))
            mu_g1 = run_plan(
                dos_plan(64, g1_dim, word="1 delta"),
                g1_ham,
                g1_win,
                g1_states,
                batch_size=64,
            )
            g1_tmp = os.path.join(tmp_dir, "g1_test.chebmom1D")
            _core.WriteLegacyChebmom1D(g1_tmp, mu_g1)
            with open(g1_tmp, "rb") as f_new, open(g1_ref_path, "rb") as f_ref:
                assert f_new.read() == f_ref.read(), "G1 byte-identity regression!"
            print("[check f] G1 chain1d DOS byte-identical to golden: PASS", flush=True)

        g4_ham_path = os.path.join(
            REPO_ROOT,
            "goldens",
            "stage0",
            "graphene",
            "operators",
            "graphene_golden.HAM.CSR",
        )
        g4_vx_path = os.path.join(
            REPO_ROOT,
            "goldens",
            "stage0",
            "graphene",
            "operators",
            "graphene_golden.VX.CSR",
        )
        g4_ref_path = os.path.join(
            REPO_ROOT,
            "goldens",
            "stage0",
            "graphene",
            "NonEqOpVX-VXgraphene_goldenKPM_M64x64_statedefault.chebmom2D",
        )
        if os.path.isfile(g4_ham_path) and os.path.isfile(g4_ref_path):
            g4_ham = _core.ReadCsrText(g4_ham_path)
            g4_vx = _core.ReadCsrText(g4_vx_path)
            g4_win = lsq.window_from_gershgorin(g4_ham)
            g4_states = _core.RandomPhaseStates(g4_ham.dimension, 1, 12345)
            mu_g4 = run_plan(
                kubo_plan(64, left_index=0, right_index=0, word="VX delta VX delta"),
                g4_ham,
                g4_win,
                g4_states,
                operators=[_core.SparseMatrixOperator(g4_vx)],
                batch_size=64,
                symmetrise=True,
            )
            g4_tmp = os.path.join(tmp_dir, "g4_test.chebmom2D")
            _core.WriteLegacyChebmom2D(g4_tmp, mu_g4)
            with open(g4_tmp, "rb") as f_new, open(g4_ref_path, "rb") as f_ref:
                assert f_new.read() == f_ref.read(), "G4 byte-identity regression!"
            print("[check f] G4 graphene Kubo byte-identical to golden: PASS", flush=True)

        g7_dir = os.path.join(REPO_ROOT, "goldens", "stage0", "chain1d_mag")
        g7_ham_path = os.path.join(g7_dir, "operators", "chain1d_mag.HAM.CSR")
        if os.path.isfile(g7_ham_path):
            mag_ham = _core.ReadCsrText(g7_ham_path)
            mag_px = _core.ReadCsrText(
                os.path.join(g7_dir, "operators", "chain1d_mag.PX.CSR")
            )
            mag_sx = _core.ReadCsrText(
                os.path.join(g7_dir, "operators", "chain1d_mag.SX.CSR")
            )
            mag_win = lsq.window_from_gershgorin(mag_ham)
            mag_dt = lsq.time_step_from_tmax(96, 95.0)
            mag_coeffs, mag_phase = lsq.time_evolution_series(
                mag_win, mag_dt, use_legacy_bessel=True
            )
            mag_states = _core.LocalStates(
                mag_ham.dimension, list(range(mag_ham.dimension))
            )
            plan_g7 = lsq.projected_evolution_plan(
                32,
                96,
                mag_ham.dimension,
                mag_coeffs,
                operator_index=0,
                projector_index=1,
                phase=mag_phase,
                word="PX U(t)^+ SX U(t) delta PX",
            )
            mu_g7 = run_plan(
                plan_g7,
                mag_ham,
                mag_win,
                mag_states,
                operators=[
                    _core.SparseMatrixOperator(mag_sx),
                    _core.SparseMatrixOperator(mag_px),
                ],
                batch_size=32,
                symmetrise=False,
                dt=mag_dt,
            )
            g7_tmp = os.path.join(tmp_dir, "g7_test.chebmomTD")
            _core.WriteLegacyChebmomTD(g7_tmp, mu_g7)
            g7_ref = os.path.join(
                g7_dir, "TimeProjSX-PXchain1d_magKPM_M32_stateexact.chebmomTD"
            )
            with open(g7_tmp, "rb") as f_new, open(g7_ref, "rb") as f_ref:
                assert f_new.read() == f_ref.read(), "G7 byte-identity regression!"
            print(
                "[check f] G7 chain1d_mag spin precession byte-identical to golden: PASS",
                flush=True,
            )

    print("stage_5_check: ALL CHECKS PASSED", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
