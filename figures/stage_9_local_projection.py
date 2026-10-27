"""Wrap-up evidence for the local projection adapter and run-file replay.

Use the existing stage-6d toy and stage-6c W=0 comparison limit. This is
stage evidence, not a new unit-test suite or a material accuracy claim.
"""
import argparse
import json
from pathlib import Path
import sys
import numpy as np

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--out", type=Path, required=True)
parser.add_argument("--repo", type=Path, default=Path(__file__).resolve().parents[1])
args = parser.parse_args()
args.out.mkdir(parents=True, exist_ok=True)
sys.path.insert(0, str(args.repo/"tests/dev"))
from stage_6d_check import augmented_kane_mele
from lsquant import _core
from lsquant.bases import BlochMetricBasis
from lsquant.apps.stage6b_siesta import SHC_WORD
from lsquant.projection.bases import ProjectedBlochBasis, ProjectedRealSpaceBasis
from lsquant.projection.words import compute_projection_moments
from lsquant.runfile import execute_runfile, RunFile
import lsquant.projection.runfile as projection_runfile

model = augmented_kane_mele()
records = []
curves = []
for route in ("lowdin", "generalised"):
    source = BlochMetricBasis(model, [4, 3, 1], velocity_mode="covariant", route=route)
    for kind in ("hard", "shift"):
        for radius in (None, 1.2):
            selected = ProjectedBlochBasis(source, dict(kind=kind, window=[-3.1, 3.1], radius=radius))
            real = ProjectedRealSpaceBasis(selected)
            for word in ("Tr[1 delta(E-H)]", SHC_WORD):
                options = dict(moments=32, states=2, seed=12345, batch=11, window=selected.window())
                b, _ = compute_projection_moments(selected, word, **options)
                r, _ = compute_projection_moments(real, word, **options)
                error = float(np.max(np.abs(b.data-r.data)))
                row = dict(route=route, kind=kind, radius=radius, word=word,
                           max_absolute_moment_difference=error, limit=1e-10,
                           passed=bool(np.isfinite(error) and error <= 1e-10))
                records.append(row)
                print(row, flush=True)
                (args.out/"local-projection.json").write_text(json.dumps(records, indent=2)+"\n")
                assert row["passed"], row
            # Exercise the actual dispatcher and its resolved-file replay.
            # The input loader alone is replaced by the existing in-memory
            # toy; parsing, adapter, compiler, engine and writers are real.
            projection_runfile.load_siesta_model = lambda *a, **k: model
            config = dict(
                model=dict(kind="siesta", dir=str(args.out), label="toy",
                           params=dict(velocity_mode="covariant")),
                basis=dict(name="real", supercell=[4, 3, 1]),
                metric=dict(policy="nonorthogonal", route=route, applier="exact"),
                projection=dict(kind=kind, window=[-3.1, 3.1], radius=radius),
                word="Tr[1 delta(E-H)]",
                words=["Tr[1 delta(E-H)]", SHC_WORD],
                expansion=dict(moments=32, states=2, seed=12345, batch=11,
                               trace="random", window="eigen"),
                reconstruction=dict(broadening_meV=10., energies=[-.1, 0., .1],
                                    bastin_backend="folded"),
                output=dict(legacy_layout=False))
            label=f"{route}-{kind}-r{radius}"
            result=execute_runfile(config, out_dir_override=str(args.out/label))
            saved=_core.ReadChebmom(result["moment_files"][0])
            expected, _=compute_projection_moments(real, config["word"], **options)
            assert np.array_equal(saved.data, expected.data), label
            replay=execute_runfile(RunFile.from_file(result["resolved_runfile"]),
                                   out_dir_override=str(args.out/(label+"-replay")))
            assert Path(result["moment_files"][0]).read_bytes() == Path(replay["moment_files"][0]).read_bytes()
            assert Path(result["moment_files"][1]).read_bytes() == Path(replay["moment_files"][1]).read_bytes()
            for filename in ("word_0.dat", "word_1.dat"):
                assert (args.out/label/filename).read_bytes() == (args.out/(label+"-replay")/filename).read_bytes()
            from lsquant.reconstruct import shc
            mu = _core.ReadChebmom(result["moment_files"][1])
            area = float(np.linalg.norm(np.cross(model.uc[0], model.uc[1])))
            energy, table = shc(mu, 10., real.num_cells, area)
            folded = np.loadtxt(args.out/label/"word_1.dat")
            assert np.array_equal(energy, folded[:, 0])
            difference = float(np.max(np.abs(table-folded[:, 1])))
            assert difference <= 1e-10, (label, difference)
            curves.append((label, np.loadtxt(args.out/label/"word_0.dat")))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
fig, axes=plt.subplots(2, 1, figsize=(8, 7))
for label, curve in curves:
    if label.startswith("lowdin-shift"):
        axes[0].plot(curve[:, 0], curve[:, 1], label=label)
axes[0].set(xlabel="Energy (eV)", ylabel="Selected DOS (states/eV/cell)",
            title="Local projected CSR: radius changes the model")
axes[0].legend()
axes[1].semilogy(range(len(records)), [max(r["max_absolute_moment_difference"], 1e-18) for r in records], "o")
axes[1].axhline(1e-10, color="black", linestyle="--", label="Existing W=0 comparison limit")
axes[1].set(xlabel="Layout/route/kind/radius/word comparison", ylabel="Maximum moment residual")
axes[1].legend()
fig.tight_layout()
fig.savefig(args.out/"local-projection.png", dpi=180)
(args.out/"complete.json").write_text(json.dumps(dict(
    comparisons=len(records), runfile_replays=len(curves), passed=True,
    replay_words=["selected DOS", "SHC"], bastin_backend="folded",
    scope="clean toy, hard/shift, both metric routes; radius 1.2 Angstrom is a truncation diagnostic, not an accepted model"), indent=2))
