#!/usr/bin/env python3
"""Item 6: prepare fixed performance run files without generating operators.

The tester supplies graphene CSR from pinned wannier2sparse and the local 6b Bi
run. Measure the resulting files unchanged at every revision/profile/thread count.
This stdlib-only driver also works before stage 8 and never modifies its inputs.
"""
from __future__ import annotations

import argparse
import copy
import json
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", type=Path, default=Path.cwd())
    parser.add_argument("--graphene-csr", required=True, type=Path)
    parser.add_argument("--bi-run", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    args = parser.parse_args()
    root, out = args.repo.resolve(), args.out.resolve()
    out.mkdir(parents=True, exist_ok=True)
    csr = args.graphene_csr.resolve()
    with (csr / "graphene.HAM.CSR").open() as handle:
        dimension, nnz = map(int, handle.readline().split())
    if dimension != 2 * 2000 * 2000:
        parser.error("graphene CSR must contain 2000 x 2000 cells with two orbitals")
    if not (csr / "graphene.out").is_file():
        parser.error("graphene CSR needs the wannier2sparse receipt graphene.out")
    graphene = {
        "model": {"kind": "csr", "dir": str(csr), "label": "graphene", "operators": []},
        "basis": {"name": "real"}, "metric": {"policy": "orthogonal"},
        "word": "Tr[delta(E-H)]",
        "expansion": {"moments": 2048, "window": [-3.2, 3.2], "trace": "random", "states": 1, "seed": 12345},
        "reconstruction": {"formulas": [], "broadening_meV": None},
        "output": {"dir": str(out / "graphene-output"), "legacy_layout": False, "write_chebmom": True}}
    # CSR run avoids constructing a clean-system eigen reference for eight million
    # sites. The input Hamiltonian still comes exclusively from wannier2sparse.
    haldane = json.loads((root / "runs/goldens/g9_haldane_bloch_exact.json").read_text())
    haldane["model"]["dir"] = str(root / "goldens/stage0/bases/haldane/model")
    haldane["word"] = "Tr[VX delta(E-H) VX delta(E-H)]"
    haldane.pop("words", None)
    haldane["expansion"].update(trace="random", states=16, seed=12345, window=[-4.5, 4.5], moments=256)
    haldane["reconstruction"] = {"kernel": "jackson", "broadening_meV": 10.0,
                                "formulas": ["greenwood", "greenwood_fft"]}
    haldane["output"] = {"dir": str(out / "haldane-output"), "legacy_layout": False, "write_chebmom": True}
    bi_source = json.loads(args.bi_run.read_text())
    bi = copy.deepcopy(bi_source)
    if bi.get("basis", {}).get("name") != "bloch":
        parser.error("Bi benchmark requires the 6b production Bloch/block spin Hall run")
    if bi.get("basis", {}).get("supercell") != [48, 48, 1]:
        parser.error("Bi benchmark requires 48 x 48 x 1 (2304 blocks)")
    bi["expansion"].update(moments=1024, states=16, trace="random", seed=12345)
    # Resolve model inputs before moving the JSON away from its source directory.
    # Unknown 6b fields are preserved; the owning branch validates their schema.
    for key in ("dir", "csr_dir"):
        value = bi.get("model", {}).get(key)
        if value and not Path(value).is_absolute():
            repo_path = root / value
            source_path = args.bi_run.resolve().parent / value
            bi["model"][key] = str((repo_path if repo_path.exists() else source_path).resolve())
    bi.pop("build_profile", None)
    bi.setdefault("output", {})["dir"] = str(out / "bi-output")
    for name, run in (("graphene", graphene), ("haldane", haldane), ("bi", bi)):
        (out / f"{name}.json").write_text(json.dumps(run, indent=2) + "\n")
    (out / "inputs.json").write_text(json.dumps({"graphene_dimension": dimension,
        "graphene_input_nnz": nnz, "graphene_producer_receipt": str(csr / "graphene.out"),
        "bi_source": str(args.bi_run.resolve()), "bi_source_config": bi_source,
        "note": "Bi reference/oracle settings preserved; use the same settings before/after."}, indent=2) + "\n")


if __name__ == "__main__":
    main()
