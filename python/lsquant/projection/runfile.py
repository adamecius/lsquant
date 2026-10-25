"""Stage-6d run executor, using the 6b word compiler/engine and physical oracle.

The none path remains the original executor. Every projected run compares
its exact word oracle with the UNPROJECTED matrices on the same mesh and
energy grid, even if an external material reference is not available.
"""
from __future__ import annotations

from pathlib import Path
import numpy as np

from .. import _core
from ..apps.stage6b_siesta import (PhaseReceipt, _write_json, oracle_energies, SHC_WORD)
from ..bases import BlochMetricBasis
from ..models.siesta import load_siesta_model
from ..oracle.spin_hall import spin_berry_curvature, sigma_xy, UNITS
from ..reconstruct import dos_curve, shc
from ..words import parse as parse_word
from ..words.grammar import OperatorSlot
from .bases import ProjectedBlochBasis, ProjectedRealSpaceBasis
from .spectral import adjoint, word_error
from .words import prepare_projection_moments, compute_projection_moments


def exact_word_curve(basis, word, energies, *, broadening, cell_area,
                     degeneracy_tolerance=1e-6):
    """Exact SHC or one-spectral-slot Gaussian curve on physical matrices.

    The oracle compares the FULL word after consistent static projection,
    including lost virtual transitions; no retained-band-only reference is
    substituted. The Gaussian width is its standard deviation in eV.
    """
    parsed = parse_word(word)
    E, C = basis.physical_eigenvalues, basis.physical_eigenvectors
    if parsed.canonical() == parse_word(SHC_WORD).canonical():
        omega, diagnostic = spin_berry_curvature(
            E, C, basis.operator_blocks("JXSZ"), basis.operator_blocks("VY"),
            hamiltonian=basis.hk_array, overlap=basis.sk_array,
            degeneracy_tolerance=degeneracy_tolerance, return_diagnostics=True)
        return sigma_xy(E, omega, energies, cell_area, broadening=broadening), diagnostic
    operators = [slot.operator_key for slot in parsed.slots if isinstance(slot, OperatorSlot)]
    if len(operators) > 1:
        raise ValueError("stage6d oracle supports SHC and single-slot spectral words")
    operator = None
    if operators:
        key = operators[0].upper()
        operator = (basis.sk_array if key in ("1", "WINDOW") and not isinstance(basis, ProjectedBlochBasis)
                    else basis.operator_blocks(key))
    weights = (np.ones_like(E) if operator is None else
               np.diagonal(adjoint(C) @ operator @ C, axis1=-2, axis2=-1).real)
    values = []
    # Bound scratch by 16*N_k*W real numbers, including on large PdSe2 meshes.
    for start in range(0, len(energies), 16):
        diff = (np.asarray(energies[start:start+16])[:, None, None]-E)/broadening
        gaussian = np.exp(-0.5*diff**2)/(np.sqrt(2*np.pi)*broadening)
        values.extend(np.sum(gaussian*weights, axis=(1, 2))/basis.num_cells)
    return np.asarray(values), dict(kernel="Gaussian", standard_deviation_eV=broadening)


def off_grid_errors(projected, operator_names, *, count=16, seed=1):
    """Compare local H/O(R) against fresh per-k projections off the FFT mesh.

    The physical source velocity ladder is evaluated on a tiny independent
    mesh; no reference is interpolated from the projected training values.
    For the optional reduced gauge, validate bands instead (basis receipt).
    """
    if projected.config["kind"] == "orthogonal":
        return {"status": "band comparison in band_validation; reduced operator gauge differs off-grid"}
    if projected.config["radius"] is None:
        return {"status": "per-k projection; local interpolation not requested (set projection.radius)"}
    from .spectral import spectral_projection
    rng = np.random.default_rng(seed)
    reduced = rng.random((count, 3))
    reduced[:, np.asarray(projected.supercell) == 1] = 0
    points = reduced @ (2*np.pi*np.linalg.inv(projected.model.uc).T)
    # Use the existing physical ladder without constructing a second engine.
    # A lightweight model proxy supplies the independent points to the 6b adapter.
    class PointModel:
        def __getattr__(self, name):
            return getattr(projected.source.model, name)
        def k_mesh(self, _):
            return points
    source = BlochMetricBasis(PointModel(), [count, 1, 1],
                               velocity_mode=projected.velocity_mode, route=projected.route,
                               seed=seed)
    exact = spectral_projection(source.hk_array, source.sk_array,
                                projected.config["window"], kind=projected.config["kind"])
    result = dict(points=count, seed=seed, matrices={})
    targets = {"H": exact.hamiltonian}
    for name in operator_names:
        original = source.sk_array if name in ("1", "I", "IDENTITY", "WINDOW") else source.operator_blocks(name)
        targets[name] = exact.project_operator(original)
    for name, target in targets.items():
        matrix = projected.local_matrices[name]
        full = matrix.blocks(points, radius=np.inf)
        cut = matrix.blocks(points)
        result["matrices"][name] = dict(
            interpolation_max_frobenius=float(np.linalg.norm(full-target, axis=(-2, -1)).max()),
            total_with_truncation_max_frobenius=float(np.linalg.norm(cut-target, axis=(-2, -1)).max()))
    return result


def execute_projection_runfile(rf, out_dir):
    """Execute clean projected SIESTA words and persist reproducible receipts.

    Hard/shift support Lowdin and generalised; orthogonal is an explicit
    dimension-reducing tool. Projection in the presence of disorder needs a
    separate 6c composition and error measurement; it is refused here.
    """
    if rf.basis["name"] not in ("real", "bloch", "eigen"):
        raise ValueError("stage6d run files require real, bloch or eigen")
    appliers = ("exact", "chebyshev") if rf.basis["name"] == "real" else ("exact",)
    if (rf.metric["policy"] != "nonorthogonal" or rf.metric["applier"] not in appliers
            or rf.metric["route"] not in ("lowdin", "generalised")):
        raise ValueError("stage6d needs Lowdin/generalised with exact metric (or real-space Chebyshev)")
    if rf.basis["name"] == "real" and rf.projection["kind"] == "orthogonal":
        raise ValueError("stage6d local nonorthogonal adapter requires hard or shift projection")
    if rf.basis["disorder"]["W"] != 0:
        raise ValueError("clean spectral projection with disorder must be composed and measured by stage6c")
    broadening_meV = rf.reconstruction["broadening_meV"]
    if broadening_meV is None or not np.isfinite(broadening_meV) or broadening_meV <= 0:
        raise ValueError("stage6d requires positive reconstruction.broadening_meV")
    if rf.reconstruction["targets"] or rf.time is not None:
        raise ValueError("stage6d run executor supports spectral and SHC words only")
    formulas = rf.reconstruction["formulas"]
    if formulas and not set(formulas) <= {"shc", "dos", "dos_curve"}:
        raise ValueError("stage6d supports dos and shc reconstruction")
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    resolved = rf.to_resolved_dict()
    resolved["output"]["dir"] = str(out)
    timing = PhaseReceipt(out/"run.resolved.json", resolved)
    params = rf.model["params"]
    mode = params.get("velocity_mode", "berry_connection")
    tolerance = float(params.get("oracle_degeneracy_tolerance", 1e-6))
    if not np.isfinite(tolerance) or tolerance <= 0:
        raise ValueError("oracle_degeneracy_tolerance must be positive and finite")
    model = load_siesta_model(rf.model["dir"], rf.model["label"], dipoles=mode == "covariant",
                              quadrature_order=int(params.get("quadrature_order", 24)),
                              quadrature_tolerance=float(params.get("quadrature_tolerance", 1e-5)),
                              quadrature_max_order=int(params.get("quadrature_max_order", 384)),
                              phase=timing.phase)
    with timing.phase("metric build"):
        original = BlochMetricBasis(model, rf.basis["supercell"], velocity_mode=mode,
                                    route=rf.metric["route"], seed=rf.expansion["seed"])
    with timing.phase("spectral projection"):
        projected = ProjectedBlochBasis(original, rf.projection)
    if rf.basis["name"] == "real":
        with timing.phase("local CSR and metric build"):
            basis = ProjectedRealSpaceBasis(
                projected, applier=rf.metric["applier"], order=rf.metric["order"],
                coeff_tol=rf.metric["coeff_tol"], max_dim=rf.metric["max_dim"], tol=rf.metric["tol"])
    else:
        basis = projected if rf.basis["name"] == "bloch" else projected.eigenbasis()
    requested = rf.expansion["window"]
    if requested == "gershgorin":
        raise ValueError("projected spectral windows must be eigen or an enclosing explicit interval")
    window = projected.window() if requested == "eigen" else _core.SpectralWindow(*requested)
    if window.lower > projected.eigenvalues.min() or window.upper < projected.eigenvalues.max():
        raise ValueError("expansion.window must enclose the actual projected spectrum, including zero modes")
    area = float(np.linalg.norm(np.cross(model.uc[0], model.uc[1])))
    resolved["expansion"]["window"] = [float(window.lower), float(window.upper)]
    resolved["reconstruction"].update(num_cells=projected.num_cells, cell_area=area)
    resolved["receipt"]["trace_basis"] = (
        "full trace with B-adjoint mapped local states" if rf.expansion["trace"] == "exact" and basis.name == "real" else
        "orbital Bloch phases mapped by B adjoint" if basis.name == "real" else
        "exact block trace" if rf.expansion["trace"] == "exact" else
                                             "independent eigenbasis phases" if rf.basis["name"] == "eigen" else
                                             "orthonormal orbital Bloch phases")
    resolved["receipt"]["metric_residuals"] = basis.metric_residuals if basis.name == "real" else projected.metric_residuals
    options = dict(moments=rf.expansion["moments"], trace=rf.expansion["trace"],
                   states=rf.expansion["states"] or 16, seed=rf.expansion["seed"],
                   batch=rf.expansion["batch"], window=window)
    names = sorted({slot.operator_key.upper() for word in rf.words for slot in parse_word(word).slots
                    if isinstance(slot, OperatorSlot)})
    with timing.phase("operators"):
        prepared = []
        for word in rf.words:
            production = prepare_projection_moments(basis, word, **options)
            if (production[2].match.entry.shape_id == "kubo_bastin" and
                    parse_word(word).canonical() != parse_word(SHC_WORD).canonical()):
                raise ValueError("stage6d rank-two oracle requires the physical JXSZ/VY SHC word")
            # Validate source operators before any oracle or moments run.
            for slot in parse_word(word).slots:
                if isinstance(slot, OperatorSlot):
                    if slot.operator_key.upper() != "WINDOW":
                        original.operator(slot.operator_key)
            prepared.append((word, production))
    with timing.phase("projection diagnostics"):
        for name in names:
            projected.operator_blocks(name)
        receipt = projected.projection_receipt()
        receipt["off_grid"] = off_grid_errors(projected, names)
        receipt["words"] = []
        resolved["receipt"]["projection"] = receipt
        _write_json(out/"projection.json", receipt)
    broadening = broadening_meV/1000
    energy_grid = oracle_energies(rf.reconstruction, original.eigenvalues, broadening_meV)
    moments, curves = [], []
    for index, (word, plan) in enumerate(prepared):
        with timing.phase("oracle", word=word):
            before, before_diag = exact_word_curve(original, word, energy_grid,
                                                   broadening=broadening, cell_area=area,
                                                   degeneracy_tolerance=tolerance)
            after, after_diag = exact_word_curve(projected, word, energy_grid,
                                                 broadening=broadening, cell_area=area,
                                                 degeneracy_tolerance=tolerance)
            error = dict(word=word, window=rf.projection["window"], oracle_error=word_error(before, after),
                         unprojected_diagnostics=before_diag, projected_diagnostics=after_diag)
            receipt["words"].append(error)
            path = out/f"word_{index}_oracle.json"
            _write_json(path, dict(energies=energy_grid.tolist(), unprojected=before.tolist(),
                                   projected=after.tolist(), **error))
            curves.append(str(path))
        with timing.phase("moments", word=word):
            mu, compiled = compute_projection_moments(basis, word, _prepared=plan)
            path = out/f"word_{index}.chebmom"
            _core.WriteChebmom(str(path), mu)
            moments.append(str(path))
            if rf.output["legacy_layout"]:
                path = out/(f"word_{index}"+compiled.reconstruction.legacy_moment_ext)
                writer = _core.WriteLegacyChebmom2D if mu.rank == 2 else _core.WriteLegacyChebmom1D
                writer(str(path), mu)
                moments.append(str(path))
            reconstruction_mu = _core.ReadChebmom(str(path)) if rf.reconstruction["legacy_text_round_trip"] else mu
        with timing.phase("reconstruction", word=word):
            if mu.rank == 2:
                E, values = shc(reconstruction_mu, broadening_meV, projected.num_cells, area,
                                backend=rf.reconstruction.get("bastin_backend", "tables"))
                units = UNITS
            else:
                E, values = dos_curve(reconstruction_mu, broadening_meV)
                values /= projected.num_cells
                units = "states/eV/cell"
            path = out/f"word_{index}.dat"
            np.savetxt(path, np.column_stack((E, values)), fmt="%.17g", header=f"energy_eV value [{units}]")
            curves.append(str(path))
            inside = (energy_grid >= E.min()) & (energy_grid <= E.max())
            if np.any(inside):
                error["kpm_vs_projected_oracle"] = word_error(after[inside], np.interp(energy_grid[inside], E, values))
                error["kpm_comparison"] = "KPM/Jackson versus Gaussian oracle; kernels differ"
            else:
                error["kpm_comparison"] = "oracle grid outside reconstruction range"
        _write_json(out/"projection.json", receipt)
        _write_json(out/"run.resolved.json", resolved)
    return dict(output_dir=str(out), resolved_runfile=str(out/"run.resolved.json"),
                moment_files=moments, curve_files=curves,
                written_files=[str(p) for p in sorted(out.iterdir()) if p.is_file()])
