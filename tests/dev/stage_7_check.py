#!/usr/bin/env python3
"""Stage-7 developer self-check (handoff/0007-words-cli-mcp/brief.md section 7).

Checks:
(a) Every word of the catalogue parses and prints back to the same canonical text;
(b) Five malformed or unsupported inputs are refused with the expected messages
    (unknown operator, rank 3, mixed basis without transform, U(t) without time,
    unknown key in the run file);
(c) runs/goldens/g1.json and runs/goldens/g5.json through the CLI give files
    byte-identical to G1 and G5 (cmp);
(d) The MCP server started in a subprocess answers words_describe and inspect
    with the same JSON as the CLI.
"""

from __future__ import annotations

import asyncio
import filecmp
import glob
import json
import os
import subprocess
import sys
import tempfile

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
_RAW_BUILD_DIR = os.environ.get("LSQUANT_BUILD_DIR", "build")
BUILD_DIR = (
    os.path.abspath(_RAW_BUILD_DIR)
    if os.path.isabs(_RAW_BUILD_DIR)
    else os.path.abspath(os.path.join(REPO, _RAW_BUILD_DIR))
)
if not os.path.isdir(os.path.join(BUILD_DIR, "lsquant_core")):
    for _cand in (os.path.join(REPO, "build", "golden"), os.path.join(REPO, "build", "native")):
        if os.path.isdir(os.path.join(_cand, "lsquant_core")):
            BUILD_DIR = _cand
            break
_PY_PATHS = [
    os.path.join(REPO, "python"),
    os.path.join(BUILD_DIR, "lsquant_core"),
    BUILD_DIR,
]
for _p in _PY_PATHS:
    if os.path.isdir(_p) and _p not in sys.path:
        sys.path.insert(0, _p)
_existing_pypath = os.environ.get("PYTHONPATH", "")
os.environ["PYTHONPATH"] = os.pathsep.join(
    [p for p in _PY_PATHS if os.path.isdir(p)]
    + ([_existing_pypath] if _existing_pypath else [])
)
os.environ.setdefault("OMP_NUM_THREADS", "1")

from lsquant.runfile import RunFile  # noqa: E402
from lsquant.words import CATALOGUE, compile as compile_word, describe, parse  # noqa: E402


def check_a_catalogue_roundtrip() -> None:
    print("== (a) Catalogue canonical parse/print round-trip ==")
    words_to_test = [entry.canonical_example for entry in CATALOGUE] + [
        "Tr[sea[VX delta(E-H) VY dG-(E)/dE]]",
        "Tr[P:PX U(t)^+ SX U(t) delta(E-H) P:PX]",
        "Tr[B^+ VX@bloch B delta(E-H)]",
    ]
    for w_text in words_to_test:
        parsed = parse(w_text)
        canonical = str(parsed)
        if canonical != w_text:
            raise AssertionError(
                f"Round-trip mismatch: input '{w_text}' -> '{canonical}'"
            )
        info = describe(w_text)
        assert info["word"]["canonical"] == w_text
        print(
            f"  PASS  {w_text:<42} -> shape={info['catalogue']['shape_id']:<20} "
            f"plan_builder={info['catalogue']['plan_builder']}"
        )


def check_b_refusals() -> None:
    print("\n== (b) Refusal of 5 malformed / unsupported inputs ==")

    cases = [
        (
            "unknown_operator",
            lambda: compile_word(
                "Tr[FOO delta(E-H)]",
                basis={"name": "real", "operators": ["VX", "VY"]},
                params={"moments": 64, "states": 1},
            ),
            "Unknown operator 'FOO'",
        ),
        (
            "rank_3",
            lambda: compile_word(
                "Tr[VX delta(E-H) VY delta(E-H) VZ delta(E-H)]",
                basis="real",
                params={"moments": 64, "states": 1},
            ),
            "rank 3",
        ),
        (
            "mixed_basis_without_transform",
            lambda: compile_word(
                "Tr[VX@bloch delta(E-H) VX@real delta(E-H)]",
                basis="real",
                params={"moments": 64, "states": 1},
            ),
            "Mixed basis without explicit transform",
        ),
        (
            "ut_without_time",
            lambda: compile_word(
                "Tr[A delta(E-H) U(t)^+ B U(t)]",
                basis="real",
                params={"moments": 64, "states": 1},
            ),
            "without 'time' parameters",
        ),
        (
            "unknown_key_in_runfile",
            lambda: RunFile.from_dict(
                {
                    "model": {
                        "kind": "csr",
                        "dir": "goldens/stage0/chain1d/operators",
                        "label": "chain1d",
                    },
                    "word": "Tr[delta(E-H)]",
                    "expansion": {"moments": 64},
                    "unknown_extra_key": 123,
                },
                repo_root=REPO,
            ),
            "Unknown key 'unknown_extra_key' in run file",
        ),
        (
            "unsupported_green_word",
            lambda: compile_word(
                "Tr[VX delta(E-H) VY G+(E)]",
                basis="real",
                params={"moments": 64, "states": 1},
            ),
            "kubo_bastin",
        ),
        (
            "eigen_window_on_csr",
            lambda: RunFile.from_dict(
                {
                    "model": {
                        "kind": "csr",
                        "dir": "goldens/stage0/chain1d/operators",
                        "label": "chain1d",
                    },
                    "word": "Tr[delta(E-H)]",
                    "expansion": {"moments": 64, "window": "eigen"},
                },
                repo_root=REPO,
            ),
            "expansion.window='eigen' is not supported for model.kind='csr'",
        ),
        (
            "unmatched_formula_for_word",
            lambda: RunFile.from_dict(
                {
                    "model": {
                        "kind": "csr",
                        "dir": "goldens/stage0/haldane/operators",
                        "label": "haldane",
                        "operators": ["VX", "VY"],
                    },
                    "word": "Tr[VX delta(E-H) VY dG-(E)/dE]",
                    "expansion": {"moments": 128, "window": [-4.5, 4.5]},
                    "reconstruction": {
                        "broadening_meV": 10.0,
                        "formulas": ["bastin", "greenwood"],
                    },
                },
                repo_root=REPO,
            ),
            "does not match any word",
        ),
        (
            "unknown_generator_param",
            lambda: RunFile.from_dict(
                {
                    "model": {
                        "kind": "generator",
                        "name": "make_chain",
                        "params": {"bogus_param": 1},
                    },
                    "word": "Tr[delta(E-H)]",
                    "expansion": {"moments": 64},
                },
                repo_root=REPO,
            ),
            "Unknown key 'bogus_param'",
        ),
    ]

    for label, fn, expected_substr in cases:
        try:
            fn()
        except ValueError as exc:
            msg = str(exc)
            if expected_substr not in msg:
                raise AssertionError(
                    f"Case '{label}' raised ValueError('{msg}'), "
                    f"missing expected substring '{expected_substr}'"
                ) from exc
            print(f"  PASS  [{label}] refused with: {msg}")
        else:
            raise AssertionError(
                f"Case '{label}' did NOT raise ValueError (expected '{expected_substr}')"
            )

    # Verify unknown generator name is refused
    try:
        RunFile.from_dict(
            {
                "model": {"kind": "generator", "name": "evil_script"},
                "word": "Tr[delta(E-H)]",
                "expansion": {"moments": 64},
            },
            repo_root=REPO,
        )
    except ValueError as exc:
        assert "Unknown generator name 'evil_script'" in str(exc)
        print(f"  PASS  [unknown_generator] refused with: {exc}")
    else:
        raise AssertionError("Unknown generator 'evil_script' was NOT refused")

    # Verify basis-dependent default window ("gershgorin" on real, "eigen" on bloch/eigen)
    for b_name, exp_win in (("real", "gershgorin"), ("bloch", "eigen"), ("eigen", "eigen")):
        rf_win = RunFile.from_dict(
            {
                "model": {
                    "kind": "wannier",
                    "dir": "goldens/stage0/bases/graphene/model",
                    "seed": "graphene",
                },
                "basis": {"name": b_name, "supercell": [48, 48, 1]},
                "word": "Tr[delta(E-H)]",
                "expansion": {"moments": 64},
            },
            repo_root=REPO,
        )
        assert rf_win.to_resolved_dict()["expansion"]["window"] == exp_win
    print("  PASS  [window_default_by_basis] real->gershgorin, bloch->eigen, eigen->eigen")

    # Also verify all committed run files under runs/goldens/*.json validate cleanly
    rf_paths = sorted(glob.glob(os.path.join(REPO, "runs", "goldens", "*.json")))
    assert rf_paths, "No run files found in runs/goldens/"
    for p in rf_paths:
        RunFile.from_file(p, repo_root=REPO)
    print(f"  PASS  [runfile_schema] all {len(rf_paths)} files in runs/goldens/*.json validated")


def _run_cli(args: list[str]) -> subprocess.CompletedProcess[str]:
    env = dict(os.environ, KPM_SEED="12345")
    return subprocess.run(
        [sys.executable, "-m", "lsquant", *args],
        cwd=REPO,
        env=env,
        capture_output=True,
        text=True,
        check=False,
    )


def _cmp_file(new_path: str, gold_path: str, tag: str) -> None:
    if not os.path.isfile(new_path):
        raise AssertionError(f"DIFFERS {tag}: missing produced file {new_path}")
    if not os.path.isfile(gold_path):
        raise AssertionError(f"DIFFERS {tag}: missing golden file {gold_path}")
    r = subprocess.run(["cmp", "-s", new_path, gold_path], check=False)
    ok = (r.returncode == 0) and filecmp.cmp(new_path, gold_path, shallow=False)
    if not ok:
        raise AssertionError(f"DIFFERS {tag}: {new_path} != {gold_path}")
    print(f"  PASS  cmp IDENTICAL: {tag}")


def check_c_cli_goldens(tmp_root: str) -> str:
    print(
        "\n== (c) CLI reproduction of G1, G5, G9 (g9_graphene_eigen.json), and G10 (g10_chain1d_real_generalised.json) =="
    )
    out_g1 = os.path.join(tmp_root, "g1")
    r1 = _run_cli(["run", "runs/goldens/g1.json", "--out", out_g1])
    if r1.returncode != 0:
        raise RuntimeError(f"CLI g1.json failed ({r1.returncode}): {r1.stderr}")

    g1_name = "SpectralOp1chain1dKPM_M64_stateexact.chebmom1D"
    g1_new = os.path.join(out_g1, g1_name)
    g1_gold = os.path.join(REPO, "goldens", "stage0", "chain1d", g1_name)
    _cmp_file(g1_new, g1_gold, f"G1 ({g1_name})")

    out_g5 = os.path.join(tmp_root, "g5")
    r5 = _run_cli(["run", "runs/goldens/g5.json", "--out", out_g5])
    if r5.returncode != 0:
        raise RuntimeError(f"CLI g5.json failed ({r5.returncode}): {r5.stderr}")

    g5_files = [
        "NonEqOpVX-VYhaldaneKPM_M128x128_stateexact.chebmom2D",
        "NonEqOpVX-VXhaldaneKPM_M128x128_stateexact.chebmom2D",
        "KuboBastin_NonEqOpVX-VYhaldaneKPM_M128x128_stateexactJACKSON.dat",
        "KuboBastinSea_NonEqOpVX-VYhaldaneKPM_M128x128_stateexactJACKSON.dat",
        "KuboBastinSurf_NonEqOpVX-VYhaldaneKPM_M128x128_stateexactJACKSON.dat",
        "KuboGreenWood_NonEqOpVX-VXhaldaneKPM_M128x128_stateexactJACKSON.dat",
    ]
    for fname in g5_files:
        _cmp_file(
            os.path.join(out_g5, fname),
            os.path.join(REPO, "goldens", "stage0", "haldane", fname),
            f"G5 ({fname})",
        )

    # Exact-trace G9 run file: runs/goldens/g9_graphene_eigen.json vs stage5_bases
    out_g9 = os.path.join(tmp_root, "g9_graphene_eigen")
    ref_g9 = os.path.join(tmp_root, "ref_g9_graphene_eigen")
    os.makedirs(ref_g9, exist_ok=True)
    r9 = _run_cli(
        ["run", "runs/goldens/g9_graphene_eigen.json", "--out", out_g9]
    )
    if r9.returncode != 0:
        raise RuntimeError(
            f"CLI g9_graphene_eigen.json failed ({r9.returncode}): {r9.stderr}"
        )
    for w_tag in ("dos", "kuboxx"):
        subprocess.run(
            [
                sys.executable,
                "-m",
                "lsquant.apps.stage5_bases",
                "--model",
                os.path.join(REPO, "goldens", "stage0", "bases", "graphene", "model"),
                "--seed",
                "graphene",
                "--basis",
                "eigen",
                "--supercell",
                "48",
                "48",
                "1",
                "--disorder",
                "0.0",
                "--disorder-seed",
                "7",
                "--moments",
                "256",
                "--trace",
                "exact",
                "--window",
                "eigen",
                "--broadening-meV",
                "40.0",
                "--word",
                w_tag,
                "--out",
                ref_g9,
            ],
            cwd=REPO,
            check=True,
            stdout=subprocess.DEVNULL,
        )
    for ref_f in sorted(os.listdir(ref_g9)):
        _cmp_file(
            os.path.join(out_g9, ref_f),
            os.path.join(ref_g9, ref_f),
            f"G9 g9_graphene_eigen ({ref_f})",
        )

    # G10 run file: runs/goldens/g10_chain1d_real_generalised.json vs stage6_nonorth
    out_g10 = os.path.join(tmp_root, "g10_chain1d_real_generalised")
    ref_g10 = os.path.join(tmp_root, "ref_g10_chain1d_real_generalised")
    os.makedirs(ref_g10, exist_ok=True)
    r10 = _run_cli(
        [
            "run",
            "runs/goldens/g10_chain1d_real_generalised.json",
            "--out",
            out_g10,
        ]
    )
    if r10.returncode != 0:
        raise RuntimeError(
            f"CLI g10_chain1d_real_generalised.json failed ({r10.returncode}): {r10.stderr}"
        )
    stage6_bin = os.path.join(BUILD_DIR, "stage6_nonorth")
    if not os.path.isfile(stage6_bin):
        stage6_bin = os.path.join(BUILD_DIR, "lsquant_core", "stage6_nonorth")
    for w_tag in ("dos", "kubo"):
        subprocess.run(
            [
                stage6_bin,
                "--fixture",
                os.path.join(REPO, "goldens", "stage0", "nonorth", "chain1d", "real"),
                "--label",
                "chain1d",
                "--route",
                "generalised",
                "--applier",
                "exact",
                "--order",
                "60",
                "--tol",
                "1e-13",
                "--window",
                "-2.0",
                "2.0",
                "--trace",
                "exact",
                "--states",
                "2048",
                "--moments",
                "64",
                "--word",
                w_tag,
                "--out",
                ref_g10,
            ],
            cwd=REPO,
            check=True,
            stdout=subprocess.DEVNULL,
        )
    subprocess.run(
        [
            sys.executable,
            "-m",
            "lsquant.apps.stage2_dos",
            "--moments",
            os.path.join(ref_g10, "SpectralOp1chain1dKPM_M64_stateexact.chebmom1D"),
            "--broadening-meV",
            "120.0",
            "--out",
            ref_g10,
        ],
        cwd=REPO,
        check=True,
        stdout=subprocess.DEVNULL,
    )
    subprocess.run(
        [
            sys.executable,
            "-m",
            "lsquant.apps.stage3_kubo",
            "--moments",
            os.path.join(
                ref_g10, "NonEqOpVX-VXchain1dKPM_M64x64_stateexact.chebmom2D"
            ),
            "--formula",
            "greenwood",
            "--broadening-meV",
            "120.0",
            "--out",
            ref_g10,
        ],
        cwd=REPO,
        check=True,
        stdout=subprocess.DEVNULL,
    )
    for ref_f in sorted(os.listdir(ref_g10)):
        _cmp_file(
            os.path.join(out_g10, ref_f),
            os.path.join(ref_g10, ref_f),
            f"G10 g10_chain1d_real_generalised ({ref_f})",
        )

    return g1_new


async def _run_mcp_check_async(g1_gold_path: str) -> None:
    from mcp import ClientSession, StdioServerParameters
    from mcp.client.stdio import stdio_client

    word_query = "Tr[VX delta(E-H) VY dG-(E)/dE]"

    cli_desc_proc = _run_cli(["words", "describe", word_query])
    if cli_desc_proc.returncode != 0:
        raise RuntimeError(f"CLI words describe failed: {cli_desc_proc.stderr}")
    cli_desc_json = json.loads(cli_desc_proc.stdout)

    cli_insp_proc = _run_cli(["inspect", g1_gold_path])
    if cli_insp_proc.returncode != 0:
        raise RuntimeError(f"CLI inspect failed: {cli_insp_proc.stderr}")
    cli_insp_json = json.loads(cli_insp_proc.stdout)

    server_params = StdioServerParameters(
        command=sys.executable,
        args=["-m", "lsquant.mcp_server", "--root", REPO],
        env=dict(os.environ, LSQUANT_MCP_ROOT=REPO),
    )

    async with stdio_client(server_params) as (read_stream, write_stream):
        async with ClientSession(read_stream, write_stream) as session:
            await session.initialize()
            tools_resp = await session.list_tools()
            tool_names = [t.name for t in tools_resp.tools]
            print(f"  MCP tools registered: {tool_names}")

            # 1. Call words_describe
            req_1 = {"tool": "words_describe", "arguments": {"word": word_query}}
            print(f"  MCP -> call_tool: {json.dumps(req_1)}")
            res_1 = await session.call_tool("words_describe", {"word": word_query})
            mcp_desc_text = res_1.content[0].text
            mcp_desc_json = json.loads(mcp_desc_text)
            print("  MCP <- words_describe response:")
            for line in mcp_desc_text.splitlines():
                print(f"    {line}")
            if mcp_desc_json != cli_desc_json:
                raise AssertionError(
                    "MCP words_describe JSON differs from CLI words describe JSON"
                )
            print("  PASS  MCP words_describe JSON == CLI words describe JSON")

            # 2. Call inspect
            rel_g1 = os.path.relpath(g1_gold_path, REPO)
            req_2 = {"tool": "inspect", "arguments": {"path": rel_g1}}
            print(f"  MCP -> call_tool: {json.dumps(req_2)}")
            res_2 = await session.call_tool("inspect", {"path": rel_g1})
            mcp_insp_text = res_2.content[0].text
            mcp_insp_json = json.loads(mcp_insp_text)
            print("  MCP <- inspect response:")
            for line in mcp_insp_text.splitlines():
                print(f"    {line}")
            if mcp_insp_json != cli_insp_json:
                raise AssertionError(
                    "MCP inspect JSON differs from CLI inspect JSON"
                )
            print("  PASS  MCP inspect JSON == CLI inspect JSON")

            # 3. Verify MCP error responses carry their reason in content[0].text
            res_err_word = await session.call_tool(
                "words_describe", {"word": "Tr[VX delta(E-H) VY G+(E)]"}
            )
            is_err_word = bool(
                getattr(res_err_word, "isError", getattr(res_err_word, "is_error", False))
            )
            assert is_err_word, "Expected isError=True for unsupported G+(E) word"
            assert "kubo_bastin" in res_err_word.content[0].text, (
                f"Expected error reason in MCP response, got: {res_err_word.content[0].text}"
            )
            res_err_root = await session.call_tool(
                "inspect", {"path": "/etc/passwd"}
            )
            is_err_root = bool(
                getattr(res_err_root, "isError", getattr(res_err_root, "is_error", False))
            )
            assert is_err_root, "Expected isError=True for path outside root"
            assert "outside configured root" in res_err_root.content[0].text, (
                f"Expected confinement error reason in MCP response, got: {res_err_root.content[0].text}"
            )
            print(
                f"  PASS  MCP error responses preserve reason: '{res_err_word.content[0].text}'"
            )


def check_d_mcp_subprocess() -> None:
    print("\n== (d) MCP server subprocess vs CLI JSON equivalence ==")
    g1_gold = os.path.join(
        REPO,
        "goldens",
        "stage0",
        "chain1d",
        "SpectralOp1chain1dKPM_M64_stateexact.chebmom1D",
    )
    asyncio.run(_run_mcp_check_async(g1_gold))


def main() -> int:
    with tempfile.TemporaryDirectory(prefix="stage7_check_") as tmp_root:
        check_a_catalogue_roundtrip()
        check_b_refusals()
        check_c_cli_goldens(tmp_root)
        check_d_mcp_subprocess()
    print("\nSTAGE 7 SELF-CHECK PASS")
    return 0


if __name__ == "__main__":
    sys.exit(main())
