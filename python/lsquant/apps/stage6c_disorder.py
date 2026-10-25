"""Stage-6c run-file executor: nonorthogonal Anderson disorder in two layouts.

The physical DOS/SHC words and all reconstruction stay in the existing
catalogue. Receipts include phase wall times and the disorder hash. Material
measurements and plateau comparisons belong to the tester on Aron.
"""
from __future__ import annotations

import argparse
from contextlib import contextmanager
import hashlib
import json
from pathlib import Path
from time import perf_counter
import numpy as np

from .. import _core
from ..disorder_bases import DisorderedBlochMetricBasis, RealSpaceMetricBasis
from ..models.siesta import load_siesta_model
from ..reconstruct import dos_curve, shc
from ..words import parse as parse_word
from .stage6b_siesta import compute_moments, prepare_moments, DOS_WORD, SHC_WORD, _write_json


class PhaseReceipt:
    """Persist timings at phase boundaries, including an interrupted/failed phase."""

    def __init__(self, path, resolved):
        self.path, self.resolved = path, resolved
        resolved["receipt"] = dict(executor="stage6c_disorder", phases=[],
                                   reference="separate real-space run with the same disorder and B-mapped phases",
                                   per_k_oracle="not applicable to disorder")

    @contextmanager
    def phase(self, name, **details):
        entry = dict(phase=name, status="running", **details)
        self.resolved["receipt"]["phases"].append(entry)
        _write_json(self.path, self.resolved)
        print(f"stage6c: {name} started {details}", flush=True)
        started = perf_counter()
        try:
            yield
        except BaseException:
            entry["status"] = "failed"
            raise
        else:
            entry["status"] = "completed"
        finally:
            entry["wall_seconds"] = perf_counter()-started
            _write_json(self.path, self.resolved)
            print(f"stage6c: {name} {entry['status']} wall_seconds={entry['wall_seconds']:.6f} {details}", flush=True)


def execute_disorder_runfile(rf, out_dir):
    """Run both physical words through the same compiler and stage-6b sampler.

    ``model.params.nonorthogonal_disorder=true`` opts W=0 Bloch files into
    this executor too. Real-space SIESTA or W>0 also selects it. The basis
    choice and metric route are independent; the per-k metric stays exact
    for Bloch, while the CSR reference uses stage-6 exact/Chebyshev appliers.
    """
    if rf.basis["name"] not in ("real", "bloch"):
        raise ValueError("stage6c supports real/bloch; disorder has no per-k eigenbasis")
    if rf.metric["policy"] != "nonorthogonal" or rf.metric["route"] not in ("lowdin", "generalised"):
        raise ValueError("stage6c needs nonorthogonal lowdin or generalised metric")
    if rf.basis["name"] == "bloch" and rf.metric["applier"] != "exact":
        raise ValueError("stage6c Bloch metric is exact per block")
    if rf.basis["name"] == "real" and rf.metric["applier"] not in ("exact", "chebyshev"):
        raise ValueError("stage6c real metric must use exact or chebyshev appliers")
    params = rf.model["params"]
    if params.get("oracle", False):
        raise ValueError("a per-k oracle is invalid for disorder; run the real-space reference")
    supported = {parse_word(DOS_WORD).canonical(), parse_word(SHC_WORD).canonical()}
    if any(parse_word(word).canonical() not in supported for word in rf.words):
        raise ValueError("stage6c supports the identity DOS and JXSZ/VY spin Hall words")
    if not set(rf.reconstruction["formulas"] or []) <= {"dos", "dos_curve", "shc"}:
        raise ValueError("stage6c reconstruction supports dos and shc")
    if rf.reconstruction["targets"]:
        raise ValueError("stage6c does not accept legacy tutorial targets")
    broadening = rf.reconstruction["broadening_meV"]
    if broadening is None or not np.isfinite(broadening) or broadening <= 0:
        raise ValueError("stage6c needs positive reconstruction.broadening_meV")
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    resolved = rf.to_resolved_dict()
    resolved["output"]["dir"] = str(out)
    # Mark resolved W=0 files so a replay cannot dispatch them to stage 6b.
    resolved["model"]["params"]["nonorthogonal_disorder"] = True
    timing = PhaseReceipt(out/"run.resolved.json", resolved)
    started = perf_counter()
    mode = params.get("velocity_mode", "covariant")
    with timing.phase("model load"):
        model = load_siesta_model(rf.model["dir"], rf.model["label"], dipoles=mode == "covariant",
                                  quadrature_order=int(params.get("quadrature_order", 24)),
                                  quadrature_max_order=int(params.get("quadrature_max_order", 384)),
                                  quadrature_tolerance=float(params.get("quadrature_tolerance", 1e-5)),
                                  phase=timing.phase)
    options = dict(disorder_w=rf.basis["disorder"]["W"], disorder_seed=rf.basis["disorder"]["seed"],
                   velocity_mode=mode, route=rf.metric["route"], seed=rf.expansion["seed"])
    with timing.phase("metric and Hamiltonian build", layout=rf.basis["name"]):
        if rf.basis["name"] == "bloch":
            basis = DisorderedBlochMetricBasis(model, rf.basis["supercell"], **options)
        else:
            basis = RealSpaceMetricBasis(model, rf.basis["supercell"], **options,
                                          applier=rf.metric["applier"], order=rf.metric["order"],
                                          coeff_tol=rf.metric["coeff_tol"], max_dim=rf.metric["max_dim"],
                                          tol=rf.metric["tol"])
    spec = rf.expansion["window"]
    if spec == "gershgorin":
        raise ValueError("stage6c uses the clean spectrum plus disorder norm enclosure, not Gershgorin")
    window = basis.window() if spec == "eigen" else _core.SpectralWindow(*spec)
    enclosure = basis.window(pad=0.)
    if window.lower > enclosure.lower or window.upper < enclosure.upper:
        raise ValueError("explicit window must enclose the clean spectrum plus sqrt(cond(S))*max|w|")
    area = float(np.linalg.norm(np.cross(model.uc[0], model.uc[1])))
    provenance = dict(model.metadata, layout=basis.name, route=basis.route,
                      metric_applier=rf.metric["applier"], metric_residuals=basis.metric_residuals,
                      disorder=dict(W_eV=basis.disorder_w, seed=basis.disorder_seed,
                                    sha256=hashlib.sha256(basis.disorder_values.astype("<c16").tobytes()).hexdigest(),
                                    distribution="numpy.default_rng uniform(-W/2,W/2), one draw per spinor coefficient",
                                    potential="(D S+S D)/2", overlap_unchanged=True,
                                    layout="n1 fastest, cell major, orbital/spinor inside cell"),
                      window=[float(window.lower), float(window.upper)],
                      window_source="clean generalised eigenvalues plus sqrt(cond(S))*max|w| enclosure",
                      overlap_bounds=list(basis.overlap_bounds), num_cells=basis.num_cells, cell_area=area,
                      trace=rf.expansion["trace"],
                      trace_basis="orthonormal orbital Bloch phases; real vectors mapped by B adjoint",
                      units="e^2/h, sigma_z Pauli")
    resolved["expansion"]["window"] = provenance["window"]
    resolved["reconstruction"].update(num_cells=basis.num_cells, cell_area=area)
    resolved["receipt"]["basis"] = provenance
    _write_json(out/"basis.json", provenance)
    np.save(out/"disorder.npy", basis.disorder_values, allow_pickle=False)
    # Compile/validate every operator before starting the first word's moments.
    with timing.phase("operators and plans"):
        prepared = [(word, prepare_moments(
            basis, word, moments=rf.expansion["moments"], trace=rf.expansion["trace"],
            states=rf.expansion["states"] or 16, seed=rf.expansion["seed"],
            window=window, batch=rf.expansion["batch"])) for word in rf.words]
    moment_files, curve_files, plateaus = [], [], []
    for word, plan in prepared:
        with timing.phase("moments", word=word):
            mu, compiled = compute_moments(basis, word, _prepared=plan)
        shape = str(rf.expansion["moments"])
        if mu.rank == 2:
            shape += "x"+shape
        stem = f"{compiled.reconstruction.legacy_moment_prefix}{model.seed_name}KPM_M{shape}_state{mu.state_label}"
        path = out/(stem+".chebmom")
        _core.WriteChebmom(str(path), mu)
        moment_files.append(str(path))
        if rf.output["legacy_layout"] or rf.reconstruction["legacy_text_round_trip"]:
            path = out/(stem+compiled.reconstruction.legacy_moment_ext)
            writer = _core.WriteLegacyChebmom2D if mu.rank == 2 else _core.WriteLegacyChebmom1D
            writer(str(path), mu)
            moment_files.append(str(path))
        reconstruction_mu = _core.ReadChebmom(str(path)) if rf.reconstruction["legacy_text_round_trip"] else mu
        with timing.phase("reconstruction", word=word):
            if mu.rank == 2:
                energies, values = shc(reconstruction_mu, broadening, basis.num_cells, area,
                                       backend=rf.reconstruction.get("bastin_backend", "tables"))
                label, header = "shc", "energy_eV sigma_z_xy [e^2/h, sigma_z Pauli]"
            else:
                energies, values = dos_curve(reconstruction_mu, broadening)
                values = values/basis.num_cells
                label, header = "dos", "energy_eV DOS [states/eV/cell]"
            curve = out/(label+".dat")
            np.savetxt(curve, np.column_stack((energies, values)), fmt="%.17g", header=header)
            curve_files.append(str(curve))
            for ef in rf.reconstruction["fermi_energies"] or [0.]:
                if not energies[0] <= ef <= energies[-1]:
                    raise ValueError("requested Fermi energy lies outside the reconstructed curve")
                plateaus.append(dict(word=word, fermi_energy_eV=ef,
                                     value=float(np.interp(ef, energies, values)),
                                     units=header.split("[", 1)[1].rstrip("]")))
    resolved["receipt"].update(total_wall_seconds=perf_counter()-started, fermi_values=plateaus)
    _write_json(out/"run.resolved.json", resolved)
    _write_json(out/"fermi_values.json", plateaus)
    return dict(output_dir=str(out), resolved_runfile=str(out/"run.resolved.json"),
                moment_files=moment_files, curve_files=curve_files,
                written_files=[str(path) for path in sorted(out.iterdir()) if path.is_file()])


def main(argv=None):
    """The CLI constructs a RunFile; it has no separate numerical path."""
    from ..runfile import execute_runfile
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--siesta", required=True)
    parser.add_argument("--label", required=True)
    parser.add_argument("--supercell", nargs=3, type=int, default=[24, 24, 1])
    parser.add_argument("--basis", choices=["bloch", "real"], default="bloch")
    parser.add_argument("--route", choices=["lowdin", "generalised"], default="lowdin")
    parser.add_argument("--applier", choices=["exact", "chebyshev"], default="exact")
    parser.add_argument("--velocity-mode", choices=["bare", "berry_connection", "covariant"], default="covariant")
    parser.add_argument("--disorder", type=float, default=0.)
    parser.add_argument("--disorder-seed", type=int, default=7)
    parser.add_argument("--word", choices=["dos", "shc", "both"], default="both")
    parser.add_argument("--moments", type=int, default=1024)
    parser.add_argument("--states", type=int, default=16)
    parser.add_argument("--seed", type=int, default=12345)
    parser.add_argument("--metric-order", type=int, default=1024)
    parser.add_argument("--broadening-meV", type=float, default=10.)
    parser.add_argument("--out", required=True)
    args = parser.parse_args(argv)
    words = [DOS_WORD, SHC_WORD] if args.word == "both" else [DOS_WORD if args.word == "dos" else SHC_WORD]
    raw = dict(model=dict(kind="siesta", dir=str(Path(args.siesta).expanduser().resolve()), label=args.label,
                          params=dict(nonorthogonal_disorder=True, velocity_mode=args.velocity_mode)),
               basis=dict(name=args.basis, supercell=args.supercell,
                          disorder=dict(W=args.disorder, seed=args.disorder_seed)),
               metric=dict(policy="nonorthogonal", route=args.route, applier=args.applier,
                           order=args.metric_order, coeff_tol=1e-14, tol=1e-10),
               word=words[0], words=words,
               expansion=dict(moments=args.moments, trace="random", states=args.states,
                              seed=args.seed, window="eigen"),
               reconstruction=dict(formulas=["dos", "shc"] if args.word == "both" else [args.word],
                                   broadening_meV=args.broadening_meV, fermi_energies=[0.]),
               output=dict(dir=args.out, write_chebmom=True, legacy_layout=True))
    print(json.dumps(execute_runfile(raw), indent=2))


if __name__ == "__main__":
    main()
