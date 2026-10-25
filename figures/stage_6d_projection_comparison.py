#!/usr/bin/env python3
"""Tester study: saved stage-6b references, then independent projection errors.

Run on Aron with the isolated stage-6d build in PYTHONPATH. Reference logs
are copied and hashed, never regenerated. All energies use the model's
existing absolute energy convention, matching runs/projection/*.json.
"""
import argparse
import gc
import hashlib
import json
from pathlib import Path
import re
import shutil
import time

import numpy as np

FROZEN_MODELS = None


def freeze_models(out, materials):
    """Run with the saved stage-6b checkout/build in PYTHONPATH.

    The old loader validates its cache key against the current SIESTA input
    files, radial inputs and quadrature settings, then validates the stored
    dipole checksum and certificate. No cache is converted or re-blessed.
    """
    from lsquant.models.siesta import load_siesta_model
    for material in materials:
        reference = Path.home()/f"t6b2/pd/deg_{material}_covariant.log"
        expected = re.search(r"dipole cache hit: ([^;]+);",reference.read_text()).group(1)
        model = load_siesta_model(str(Path.home()/f"from_fatnode/shc_{material}/run"),material,
                                 dipoles=True,quadrature_tolerance=5e-5)
        cache = model.metadata["dipole_cache"]
        if not cache["hit"] or cache["path"] != expected:
            raise ValueError("Reference loader did not reuse the exact cache named in the golden log")
        target = out/f"{material}.npz"
        if target.exists():
            raise FileExistsError(f"Refusing to overwrite frozen input {target}")
        np.savez_compressed(target,uc=model.uc,tau=model.tau,r_int=model.r_int,
                            hoppings=model.hoppings,overlaps=model.overlaps,dipoles=model.dipoles,
                            metadata=json.dumps(model.metadata),input_tolerance=model.input_tolerance)
        print("FROZEN",material,hashlib.sha256(target.read_bytes()).hexdigest(),cache,flush=True)


def load_model(material, mode):
    """Consume matrices exported by the saved reference's validated loader.

    This is intentionally separate from current dipole recertification: the
    experiment changes projection, not the physical input matrices.
    """
    from lsquant.models.siesta import SiestaModel
    path = FROZEN_MODELS / f"{material}.npz"
    with np.load(path, allow_pickle=False) as arrays:
        metadata = json.loads(str(arrays["metadata"]))
        metadata["frozen_model_sha256"] = hashlib.sha256(path.read_bytes()).hexdigest()
        return SiestaModel(arrays["uc"], arrays["tau"], arrays["r_int"],
                           arrays["hoppings"], arrays["overlaps"],
                           arrays["dipoles"] if mode == "covariant" else None,
                           seed_name=material, metadata=metadata,
                           input_tolerance=float(arrays["input_tolerance"]))


def save(path, value):
    path.write_text(json.dumps(value, indent=2, default=lambda x: x.tolist()) + "\n")


def saved_reference(path, mesh):
    section = path.read_text().split(f"{mesh}x{mesh}:", 1)[1]
    line = next(line for line in section.splitlines() if "degeneracy_tolerance=1e-06:" in line)
    pairs = re.findall(r"E=([+\-\d.]+):([+\-\d.]+)", line)
    return np.array([float(x) for x, _ in pairs]), np.array([float(y) for _, y in pairs])


def k_space(out, meshes, materials, modes):
    from lsquant.bases import BlochMetricBasis
    from lsquant.apps.stage6b_siesta import SHC_WORD
    from lsquant.projection.bases import ProjectedBlochBasis
    from lsquant.projection.runfile import exact_word_curve
    refs = out / "references"
    refs.mkdir(exist_ok=True)
    provenance = {}
    for material in ("bi", "pdse2"):
        for mode in ("berry_connection", "covariant"):
            src = Path.home() / f"t6b2/pd/deg_{material}_{mode}.log"
            dest = refs / src.name
            if not dest.exists():
                shutil.copy2(src, dest)
            provenance[src.name] = dict(source=str(src), sha256=hashlib.sha256(dest.read_bytes()).hexdigest())
    save(out / "reference_provenance.json", provenance)
    for material in materials:
        for mode in modes:
            print(f"LOAD {material} {mode}", flush=True)
            model = load_model(material, mode)
            area = float(np.linalg.norm(np.cross(model.uc[0], model.uc[1])))
            for mesh in meshes:
                prefix = f"{material}_{mode}_{mesh}"
                energies, reference = saved_reference(refs/f"deg_{material}_{mode}.log", mesh)
                start = time.monotonic()
                source = BlochMetricBasis(model, [mesh, mesh, 1], velocity_mode=mode, route="lowdin", seed=12345)
                baseline, diagnostics = exact_word_curve(source, SHC_WORD, energies, broadening=0.,
                                                          cell_area=area, degeneracy_tolerance=1e-6)
                record = dict(material=material, velocity_mode=mode, mesh=[mesh,mesh,1], energies=energies,
                              model_sha256=model.metadata["frozen_model_sha256"],
                              saved_reference=reference, unprojected=baseline,
                              baseline_minus_saved=baseline-reference, diagnostics=diagnostics,
                              seconds=time.monotonic()-start, units="e^2/h, sigma_z Pauli",
                              original_bounds=[source.eigenvalues.min(), source.eigenvalues.max()],
                              reference_precision="five printed decimal places", projections=[])
                print(f"BASELINE {prefix}: saved={reference} current={baseline}", flush=True)
                save(out/f"{prefix}.json", record)
                for width in (3., 5.):
                    for kind in ("shift", "hard"):
                        start = time.monotonic()
                        projected = ProjectedBlochBasis(source, dict(kind=kind, window=[-width, width]))
                        values, diag = exact_word_curve(projected, SHC_WORD, energies, broadening=0.,
                                                       cell_area=area, degeneracy_tolerance=1e-6)
                        row = dict(kind=kind, window=[-width,width], projected=values,
                                   projected_minus_unprojected=values-baseline,
                                   projected_minus_saved=values-reference,
                                   relative_to_unprojected=(values-baseline)/np.maximum(np.abs(baseline),1e-15),
                                   receipt=projected.projection_receipt(), diagnostics=diag,
                                   seconds=time.monotonic()-start)
                        record["projections"].append(row)
                        save(out/f"{prefix}.json", record)
                        print(f"PROJECTED {prefix} {kind} +/-{width}: {values}; error={values-baseline}", flush=True)
                        del projected
                        gc.collect()
                del source
                gc.collect()
            del model
            gc.collect()


def real_space(out, materials):
    """Validate actual CSR exports on a tractable identical finite mesh.

    The 4x4 mesh tests representation equivalence, not bulk convergence.
    S is exported unchanged. Its inverse square root is applied by the exact
    Bloch FFT composite, while H and currents use the exported CSR matrices.
    Identical transformed random vectors isolate representation error.
    """
    from lsquant import _core
    from lsquant.bases import BlochMetricBasis
    from lsquant.apps.stage6b_siesta import SHC_WORD
    from lsquant.projection.bases import ProjectedBlochBasis
    from lsquant.projection.runfile import exact_word_curve
    from lsquant.projection.real_space import cell_gauge
    from lsquant.words import compile as compile_word

    def sparse_operator(matrix):
        return _core.SparseMatrixOperator(_core.SparseMatrix.from_csr(
            matrix.shape[0], matrix.indptr.astype(np.int64), matrix.indices.astype(np.int64),
            matrix.data.astype(np.complex128)))

    def run_moments(name, h, ops, vectors, window):
        states = _core.UserStates(vectors)
        metric = dict(policy="orthogonal", route="lowdin")
        compiled = compile_word(SHC_WORD, basis=dict(name=name, operators=ops, metric=metric),
                                params=dict(expansion=dict(moments=64, states=len(vectors), window=window), metric=metric))
        scaled = _core.RescaledOperator(h, window)
        engine = _core.MomentEngine(scaled, compiled.operators, states, 64)
        mu = engine.run(compiled.plan)
        return np.asarray(mu.data).copy()

    for material in materials:
        prior = [out/f"{material}_{kind}_{width}.json" for kind in ("shift","hard") for width in (3,5)]
        if all(path.exists() and json.loads(path.read_text()).get("passed") for path in prior):
            print(f"REAL_ALREADY_COMPLETE {material}", flush=True)
            continue
        print(f"REAL_LOAD {material}", flush=True)
        model = load_model(material, "covariant")
        source = BlochMetricBasis(model, [4,4,1], velocity_mode="covariant", route="lowdin", seed=12345)
        area = float(np.linalg.norm(np.cross(model.uc[0],model.uc[1])))
        energies = [0.] if material == "bi" else [1.052]
        baseline, _ = exact_word_curve(source, SHC_WORD, energies, broadening=0., cell_area=area)
        phases = np.ascontiguousarray(np.exp(-1j*(source.k_points @ model.tau.T))).ravel()
        B = _core.BlochTransform(4,4,1,source.num_wann,phases,
                                 _core.BlochDirection.RealToBloch)
        Bd = B.adjoint()
        Yk = _core.BlockDiagonalOperator(source.inverse_sqrt_array)
        Yr = _core.SandwichOperator(Bd,Yk)
        rng = np.random.default_rng(12345)
        vectors = [np.ascontiguousarray(np.exp(2j*np.pi*rng.random(source.dimension))) for _ in range(2)]
        real_vectors = [Bd.multiply(v) for v in vectors]
        for width in (3.,5.):
            for kind in ("shift","hard"):
                start = time.monotonic()
                proj = ProjectedBlochBasis(source, dict(kind=kind,window=[-width,width]))
                oracle, _ = exact_word_curve(proj, SHC_WORD, energies, broadening=0., cell_area=area)
                print(f"REAL_EXPORT {material} {kind} +/-{width}", flush=True)
                matrices = proj.supercell_operators(["JXSZ","VY"])
                operators = {key:sparse_operator(value) for key,value in matrices.items()}
                errors = {}
                for key, mat in matrices.items():
                    blocks = proj.hk_array if key == "H" else proj.operator_blocks(key)
                    ok = _core.BlockDiagonalOperator(blocks)
                    residual = max(np.linalg.norm(B.multiply(mat @ r)-ok.multiply(k))/
                                   max(np.linalg.norm(ok.multiply(k)),1e-15)
                                   for r,k in zip(real_vectors,vectors))
                    hermiticity = np.linalg.norm((mat-mat.getH()).data)/max(np.linalg.norm(mat.data),1e-15)
                    # Independent full-matrix check using a direct discrete Fourier
                    # sum, without the origin WS exporter or the C++ FFT operator.
                    width_orbitals = source.num_wann
                    cell_blocks = cell_gauge(blocks, source.k_points, model.tau)
                    hoppings = np.fft.fftn(cell_blocks.reshape(1,4,4,width_orbitals,width_orbitals),
                                          axes=(0,1,2)).reshape(16,width_orbitals,width_orbitals)/16
                    coordinates = np.array(np.unravel_index(np.arange(16),(4,4,1),order="F")).T
                    squared_error, squared_scale, max_element = 0., 0., 0.
                    for cell,coordinate in enumerate(coordinates):
                        delta_cells = (coordinates-coordinate) % np.array([4,4,1])
                        indices = np.ravel_multi_index(delta_cells.T,(4,4,1),order="F")
                        expected = np.concatenate(hoppings[indices],axis=1)
                        actual = mat[cell*width_orbitals:(cell+1)*width_orbitals].toarray()
                        difference = actual-expected
                        squared_error += float(np.linalg.norm(difference)**2)
                        squared_scale += float(np.linalg.norm(expected)**2)
                        max_element = max(max_element,float(np.max(np.abs(difference))))
                    errors[key] = dict(action_relative=float(residual), hermiticity_relative=float(hermiticity),
                                       full_matrix_relative_frobenius=float(np.sqrt(squared_error/squared_scale)),
                                       full_matrix_max_element=max_element,
                                       nnz=int(mat.nnz))
                window = proj.window()
                hk = _core.LowdinOperator(Yk, _core.BlockDiagonalOperator(proj.hk_array))
                hr = _core.LowdinOperator(Yr, operators["H"])
                opsk = {key:_core.LowdinOperator(Yk,proj.operator(key)) for key in ("JXSZ","VY")}
                opsr = {key:_core.LowdinOperator(Yr,operators[key]) for key in ("JXSZ","VY")}
                muk = run_moments("bloch",hk,opsk,vectors,window)
                mur = run_moments("real",hr,opsr,real_vectors,window)
                delta = np.abs(mur-muk)
                record = dict(material=material, mesh=[4,4,1], dimension=source.dimension,
                              model_sha256=model.metadata["frozen_model_sha256"],
                              kind=kind, window=[-width,width], energies=energies,
                              finite_mesh_unprojected=baseline, finite_mesh_projected=oracle,
                              matrices=errors, moments=64, trace_vectors=2, seed=12345,
                              raw_moments_max_absolute=float(delta.max()),
                              raw_moments_relative_frobenius=float(np.linalg.norm(mur-muk)/np.linalg.norm(muk)),
                              seconds=time.monotonic()-start, radius=None,
                              scope="4x4 finite-system CSR/FFT equivalence; not a converged material golden",
                              metric="original exported S; exact FFT application of S^-1/2",
                              passed=bool(max(row["action_relative"] for row in errors.values()) < 1e-10 and
                                          max(row["hermiticity_relative"] for row in errors.values()) < 1e-10 and
                                          max(row["full_matrix_relative_frobenius"] for row in errors.values()) < 1e-10 and
                                          np.linalg.norm(mur-muk)/np.linalg.norm(muk) < 1e-10))
                stem=f"{material}_{kind}_{width:g}"
                save(out/f"{stem}.json", record)
                np.savez(out/f"{stem}_moments.npz", bloch=muk,real=mur)
                print(f"REAL_RESULT {stem} pass={record['passed']} moments_relative={record['raw_moments_relative_frobenius']:.3g}", flush=True)
                del proj, matrices, operators, hk,hr,opsk,opsr
                gc.collect()
    failures = [str(out/f"{material}_{kind}_{width}.json")
                for material in materials for kind in ("shift","hard") for width in (3,5)
                if not json.loads((out/f"{material}_{kind}_{width}.json").read_text())["passed"]]
    if failures:
        raise RuntimeError(f"Real-space equivalence failed: {failures}")


def plot_results(out):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    records = []
    for path in sorted((out/"kspace").glob("*.json")):
        if path.name == "reference_provenance.json":
            continue
        record = json.loads(path.read_text())
        if len(record.get("projections",[])) != 4:
            continue
        # Use the implementation's published convention; no numerical rescaling.
        record["units"] = "e^2/h, sigma_z Pauli"
        save(path, record)
        records.append(record)
    fig, axes = plt.subplots(2, 4, figsize=(16,7), sharex=True,
                             gridspec_kw=dict(height_ratios=[2,1]))
    for column,(material,mode) in enumerate((m,v) for m in ("bi","pdse2")
                                          for v in ("covariant","berry_connection")):
        ax, residual = axes[:,column]
        for mesh,color in ((24,"#176b9a"),(48,"#d96c26")):
            matching = [r for r in records if r["material"] == material and
                        r["velocity_mode"] == mode and r["mesh"][0] == mesh]
            if not matching:
                continue
            record = matching[0]
            index = int(np.argmin(np.abs(np.asarray(record["energies"])-(0 if material == "bi" else 1.052))))
            reference = record["saved_reference"][index]
            ax.axhline(reference,color=color,linestyle=":",label=f"{mesh}x{mesh} saved reference")
            for kind,marker,line in (("shift","o","-"),("hard","x","--")):
                rows = [r for r in record["projections"] if r["kind"] == kind]
                x = [r["window"][1] for r in rows]
                y = [r["projected"][index] for r in rows]
                ax.plot(x,y,marker=marker,linestyle=line,color=color,label=f"{mesh}x{mesh} {kind}")
                residual.plot(x,100*(np.array(y)-reference)/abs(reference),marker=marker,linestyle=line,color=color)
        ax.set_title(f"{'Bi' if material == 'bi' else 'PdSe2'} / {mode.replace('_',' ')}")
        ax.set_ylabel("Exact SHC [e^2/h, Pauli spin]")
        ax.grid(alpha=.2)
        ax.legend(fontsize=7)
        residual.axhline(0,color="black",linewidth=.7)
        residual.set_ylabel("Error vs saved [%]")
        residual.set_xlabel("Projection window half-width [eV]")
        residual.set_xticks([3,5])
        residual.grid(alpha=.2)
    pending = len(records) < 8
    fig.suptitle("Projection versus frozen stage-6b SIESTA references - exact k-space calculation" +
                 ("\nPdSe2 48x48 still running; completed meshes shown" if pending else ""))
    fig.tight_layout()
    fig.savefig(out/"projection-vs-golden.png",dpi=170)
    plt.close(fig)
    (out/"projection-vs-golden-caption.txt").write_text(
        "Exact clean-system spin Hall response at E=0 eV (Bi) and E=1.052 eV (PdSe2), "
        "in the implementation convention e^2/h with Pauli sigma_z. Dotted lines are "
        "frozen stage-6b SIESTA references on each matching mesh; markers are the "
        "delivered hard/shift projection. Residuals use the saved five-decimal reference. "
        "Every calculation uses the exact cached dipoles named by the reference logs. "
        "Windows are [-3,3] and [-5,5] eV relative to the existing model energy zero. "
        "The separate untruncated 4x4 real-space check validates full exported CSR matrices "
        "and identical-vector M=64 moments, not bulk mesh convergence. " +
        ("PdSe2 48x48 remains pending.\n" if pending else "All requested exact k-space meshes are complete.\n"))


def kpm_comparison(out, materials):
    """Separate finite-M error from the already measured projection error."""
    from lsquant import _core
    from lsquant.bases import BlochMetricBasis
    from lsquant.apps.stage6b_siesta import SHC_WORD
    from lsquant.projection.bases import ProjectedBlochBasis
    from lsquant.projection.runfile import exact_word_curve
    from lsquant.projection.words import compute_projection_moments
    from lsquant.reconstruct import shc
    for material in materials:
        model = load_model(material, "covariant")
        source = BlochMetricBasis(model,[24,24,1],velocity_mode="covariant",route="lowdin",seed=12345)
        area = float(np.linalg.norm(np.cross(model.uc[0],model.uc[1])))
        energies = np.array([0.]) if material == "bi" else np.array([1.052])
        for width in (3.,5.):
            proj = ProjectedBlochBasis(source,dict(kind="shift",window=[-width,width]))
            exact, _ = exact_word_curve(proj,SHC_WORD,energies,broadening=.01,cell_area=area)
            for moments in (256,1024):
                print(f"KPM_START {material} shift +/-{width} M={moments} R=16",flush=True)
                start = time.monotonic()
                mu, _ = compute_projection_moments(proj.eigenbasis(),SHC_WORD,moments=moments,
                                                   trace="random",states=16,seed=12345)
                moment_seconds = time.monotonic()-start
                stem=f"{material}_shift_{width:g}_M{moments}"
                _core.WriteChebmom(str(out/f"{stem}.chebmom"),mu)
                start = time.monotonic()
                energy,curve = shc(mu,10.,proj.num_cells,area)
                order = np.argsort(energy)
                values = np.interp(energies,energy[order],curve[order])
                np.savetxt(out/f"{stem}.dat",np.column_stack((energy,curve)))
                record = dict(material=material,mesh=[24,24,1],kind="shift",window=[-width,width],
                              model_sha256=model.metadata["frozen_model_sha256"],
                              moments=moments,states=16,seed=12345,broadening_meV=10.,energies=energies,
                              projected_exact=exact,kpm=values,kpm_minus_projected_exact=values-exact,
                              moment_seconds=moment_seconds,reconstruction_seconds=time.monotonic()-start,
                              note="KPM uses delivered reconstruction kernel; oracle uses Gaussian occupancy; probes are inside gaps")
                save(out/f"{stem}.json",record)
                print(f"KPM_RESULT {stem}: kpm={values} projected_exact={exact}; moments_seconds={moment_seconds:.3f}",flush=True)


def main():
    global FROZEN_MODELS
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--meshes", type=int, nargs="+", default=[24,48])
    parser.add_argument("--phase", choices=["freeze", "kspace", "realspace", "kpm", "plot"], default="kspace")
    parser.add_argument("--materials", choices=["bi", "pdse2"], nargs="+", default=["bi","pdse2"])
    parser.add_argument("--modes", choices=["covariant", "berry_connection"], nargs="+", default=["covariant","berry_connection"])
    parser.add_argument("--frozen-models", type=Path)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    FROZEN_MODELS = args.frozen_models
    if args.phase not in ("freeze", "plot") and FROZEN_MODELS is None:
        parser.error("--frozen-models is required to preserve reference inputs")
    if args.phase == "freeze":
        freeze_models(args.out,args.materials)
    elif args.phase == "kspace":
        k_space(args.out, args.meshes, args.materials, args.modes)
    elif args.phase == "realspace":
        real_space(args.out, args.materials)
    elif args.phase == "kpm":
        kpm_comparison(args.out, args.materials)
    else:
        plot_results(args.out)
    print(f"{args.phase.upper()}_COMPLETE", flush=True)


if __name__ == "__main__":
    main()
