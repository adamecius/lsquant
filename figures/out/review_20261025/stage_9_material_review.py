"""Material review panels from frozen SIESTA models; no golden is modified.

Top: full/selected bands, Gaussian DOS, broadened Kubo-Greenwood charge
conductivity, exact spin Hall response and the saved projected KPM curve.
Bottom: upper-valence spin sectors. Degenerate eigenvectors are resolved
by diagonalising physical sigma_z inside the valence multiplet, not by
plotting the arbitrary vectors returned by an eigensolver.
"""
import argparse
import hashlib
import importlib.util
import json
from pathlib import Path
import sys
import time

import numpy as np
from scipy.linalg import eigvalsh
from scipy.spatial import Voronoi
from scipy.special import erfc

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--root", type=Path, required=True)
parser.add_argument("--repo", type=Path, required=True)
parser.add_argument("--out", type=Path, required=True)
parser.add_argument("--material", choices=("bi", "pdse2"), required=True)
parser.add_argument("--plot-only", action="store_true")
args = parser.parse_args()
args.out.mkdir(parents=True, exist_ok=True)
material = args.material
reference_mu = 0.0 if material == "bi" else 1.052
label = r"$\mathrm{Bi}$" if material == "bi" else r"$\mathrm{PdSe}_2$"
data_path = args.out/(material+"-review-data.npz")
receipt_path = args.out/(material+"-review.json")
started = time.monotonic()


def phase(message):
    print(json.dumps(dict(material=material, phase=message,
                         wall_seconds=time.monotonic()-started)), flush=True)


def k_path(lattice):
    vertices = (np.array([[0, 0, 0], [.5, 0, 0], [2/3, 1/3, 0], [0, 0, 0]])
                if material == "bi" else
                np.array([[0, 0, 0], [.5, 0, 0], [.5, .5, 0], [0, .5, 0], [0, 0, 0]]))
    labels = ([r"$\Gamma$", r"$M$", r"$K$", r"$\Gamma$"] if material == "bi"
              else [r"$\Gamma$", r"$X$", r"$S$", r"$Y$", r"$\Gamma$"])
    reciprocal = 2*np.pi*np.linalg.inv(lattice).T
    points = np.concatenate([np.linspace(a, b, 41)[:-1] for a, b in zip(vertices[:-1], vertices[1:])]
                            + [vertices[-1:]]) @ reciprocal
    distance = np.r_[0, np.cumsum(np.linalg.norm(np.diff(points, axis=0), axis=1))]
    return points, distance, distance[np.arange(len(vertices))*40], labels


def gaussian_dos(energies, mesh_eigenvalues, kept=None, eta=.05):
    result = []
    for e in energies:
        weight = np.exp(-.5*((e-mesh_eigenvalues)/eta)**2)/(np.sqrt(2*np.pi)*eta)
        if kept is not None:
            weight *= kept
        result.append(weight.sum()/len(mesh_eigenvalues))
    return np.array(result)


if not args.plot_only:
    from lsquant.bases import BlochMetricBasis
    from lsquant.projection.bases import ProjectedBlochBasis
    from lsquant.projection.spectral import spectral_projection
    from lsquant.projection.runfile import exact_word_curve
    from lsquant.apps.stage6b_siesta import SHC_WORD
    from lsquant.models.nonorthogonal import applied_spin_operators
    spec = importlib.util.spec_from_file_location("frozen_study", args.root/"projection-comparison.py")
    harness = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(harness)
    harness.FROZEN_MODELS = args.root/"projection-reference-study/frozen-models"
    phase("load frozen SIESTA reference")
    model = harness.load_model(material, "covariant")
    source = BlochMetricBasis(model, [24, 24, 1], velocity_mode="covariant", route="lowdin")
    phase("clean mesh prepared")
    projected = ProjectedBlochBasis(source, dict(kind="shift", window=[-5., 5.]))
    grid = np.linspace(reference_mu-3, reference_mu+3, 601)
    area = float(np.linalg.norm(np.cross(model.uc[0], model.uc[1])))
    full_dos = gaussian_dos(grid, source.eigenvalues)
    selected_dos = gaussian_dos(grid, source.eigenvalues,
                                (source.eigenvalues >= -5)&(source.eigenvalues <= 5))
    phase("DOS complete; spin Hall oracle starting")
    full_shc, diag_full = exact_word_curve(source, SHC_WORD, grid, broadening=.01, cell_area=area)
    selected_shc, diag_selected = exact_word_curve(projected, SHC_WORD, grid,
                                                   broadening=.01, cell_area=area)
    phase("spin Hall oracle complete; charge conductivity starting")
    vx = source.operator_blocks("VX")
    charge = np.zeros((2, len(grid)))
    # Kubo-Greenwood in e^2/h: 2*pi^2/(Nk*A) Tr[Vx delta_eta Vx delta_eta].
    # Vx is the physical covariant derivative (eV Angstrom, hbar=1).
    # Gaussian eta makes a finite clean-crystal diagnostic, not a lifetime
    # or a converged disordered DC conductivity. It includes intraband terms.
    for k, (e, c) in enumerate(zip(source.eigenvalues, source.eigenvectors)):
        v = c.conj().T @ vx[k] @ c
        pairs = np.abs(v)**2
        delta = np.exp(-.5*((grid[:, None]-e)/.05)**2)/(np.sqrt(2*np.pi)*.05)
        charge[0] += np.einsum("en,nm,em->e", delta, pairs, delta, optimize=True)
        selected = delta*((e >= -5)&(e <= 5))
        charge[1] += np.einsum("en,nm,em->e", selected, pairs, selected, optimize=True)
    charge *= 2*np.pi**2/(source.num_cells*area)
    if np.min(charge) < -1e-12:
        raise RuntimeError("negative longitudinal Gaussian trace")
    phase("charge conductivity complete; valence spin sectors starting")
    q = applied_spin_operators(model.num_wann)
    spin = np.full((2, source.num_cells, 3), np.nan)
    spin_trace = np.zeros((source.num_cells, 3))
    multiplicity, valence_energy, sector_gap = [], [], []
    spin_projection_error = 0.
    occupied = (source.eigenvalues < reference_mu).sum(axis=1)
    if not np.all(occupied == occupied[0]):
        raise RuntimeError("review chemical potential does not select a fixed valence manifold")
    band = int(occupied[0])-1
    for k, (e, c, s) in enumerate(zip(source.eigenvalues, source.eigenvectors, source.sk_array)):
        indices = np.flatnonzero((np.abs(e-e[band]) <= 1e-6)&(e < reference_mu))
        if np.any((e[indices] < -5)|(e[indices] > 5)):
            raise RuntimeError("valence multiplet outside retained window")
        frame = c[:, indices]
        small = [frame.conj().T @ s @ q["S"+axis] @ frame for axis in "XYZ"]
        small = [.5*(a+a.conj().T) for a in small]
        z, u = np.linalg.eigh(small[2])
        spin_trace[k] = [np.trace(a).real/len(indices) for a in small]
        multiplicity.append(len(indices))
        valence_energy.append(e[band])
        sector_gap.append(float(np.min(np.abs(z))))
        # A nondegenerate band needs no spin-sector convention. If it is
        # degenerate, select the unique maximal-sigma_z eigenstate; never
        # resolve an additional degeneracy by an arbitrary array column.
        resolved = len(z) == 1 or np.sum(np.abs(z-z[-1]) < 1e-8) == 1
        if resolved:
            vector = u[:, -1]
            spin[0, k] = [float(np.vdot(vector, a @ vector).real) for a in small]
            selected_state = projected.projection.projector[k] @ (frame @ vector)
            normalisation = float(np.vdot(selected_state, s @ selected_state).real)
            spin[1, k] = [float(np.vdot(selected_state, s @ q["S"+a] @ selected_state).real
                                / normalisation) for a in "XYZ"]
        pframe = projected.projection.projector[k] @ frame
        for a in "XYZ":
            original = frame.conj().T @ s @ q["S"+a] @ frame
            selected = pframe.conj().T @ s @ q["S"+a] @ pframe
            spin_projection_error = max(spin_projection_error, float(np.max(np.abs(original-selected))))
    if np.nanmax(np.linalg.norm(spin, axis=-1)) > 1+1e-8:
        raise RuntimeError("spin expectation exceeds Pauli bound")
    phase("valence spin sectors complete; band path starting")
    points, distance, ticks, ticklabels = k_path(model.uc)
    path_h, path_s = model.Hk(points), model.Sk(points)
    path_projection = spectral_projection(path_h, path_s, [-5., 5.], kind="shift")
    original_bands = path_projection.original_eigenvalues
    projected_bands = np.array([np.diag(c.conj().T @ h @ c).real for c, h in
                                zip(path_projection.eigenvectors, path_projection.hamiltonian)])
    kept = path_projection.kept
    deviation = (projected_bands-original_bands)[kept]
    sorted_projected = np.array([eigvalsh(h, s) for h, s in
                                 zip(path_projection.hamiltonian, path_s)])
    spectral_error = float(np.max(np.abs(sorted_projected-np.sort(path_projection.eigenvalues, axis=1))))
    selected_bands = np.where(kept, projected_bands, np.nan)
    old = args.root/"reduced-validation/materials"/(material+"-covariant-k24-shift5")
    sampled = np.loadtxt(old/"M1024.dat")
    gap = [float(source.eigenvalues[:, band].max()), float(source.eigenvalues[:, band+1].min())]
    receipt = dict(material=material, source_tree="bf70c9b8c8b8ec66df329265d6c170a1cb315e0e",
                   model_sha256=model.metadata["frozen_model_sha256"],
                   material_model="SIESTA; covariant velocity; exact Lowdin metric; W=0",
                   mesh=[24, 24, 1], projection=dict(kind="shift", window_eV=[-5, 5], radius=None),
                   reference_mu_eV=reference_mu, stored_fermi_level_eV=0.,
                   energy_reference="Reference chemical potential relative to stored SIESTA Fermi zero",
                   gaussian_DOS_KG_eta_eV=.05, gaussian_SHC_occupation_eta_eV=.01,
                   charge_units="e^2/h", spin_units="Pauli", shc_units="e^2/h, Pauli spin current",
                   longitudinal_formula="2*pi^2/(Nk*A) sum_knm |C_n^+ VX C_m|^2 delta_eta(E-e_n) delta_eta(E-e_m)",
                   longitudinal_scope="Finite Gaussian width and k mesh; clean broadened diagnostic, not converged DC transport",
                   projected_DOS="Selected physical identity only; shifted complement spikes excluded",
                   valence_band_index_zero_based=band, gap_eV=gap, valence_multiplicity=sorted(set(multiplicity)),
                   spin_texture="Upper valence band; if energy-degenerate, unique maximal sigma_z eigenstate within that multiplet. Original versus explicitly projected state.",
                   min_abs_projected_spin_eigenvalue=float(min(sector_gap)),
                   max_multiplet_average_spin=float(np.max(np.linalg.norm(spin_trace, axis=1))),
                   spin_projection_max_error=spin_projection_error,
                   spin_texture_max_difference=float(np.nanmax(np.abs(spin[1]-spin[0]))),
                   unresolved_spin_points=int(np.sum(~np.all(np.isfinite(spin[0]), axis=1))),
                   band_path_labels=ticklabels, path_MAE_meV=float(np.mean(np.abs(deviation))*1000),
                   path_max_band_error_meV=float(np.max(np.abs(deviation))*1000),
                   projected_path_eigenvalue_max_error_eV=spectral_error,
                   shc_at_reference=dict(full=float(np.interp(reference_mu, grid, full_shc)),
                                         projected=float(np.interp(reference_mu, grid, selected_shc)),
                                         saved_KPM_M1024_R16=float(np.interp(reference_mu, sampled[:, 0], sampled[:, 1]))),
                   max_SHC_multiplet_tolerance_eV=max(diag_full["tolerance_eV"], diag_selected["tolerance_eV"]),
                   saved_KPM_file=str(old/"M1024.dat"),
                   saved_KPM_sha256=hashlib.sha256((old/"M1024.dat").read_bytes()).hexdigest(),
                   note="Review figure; does not close any material gate.")
    if receipt["path_MAE_meV"] >= 10:
        raise RuntimeError("retained band MAE fails user criterion")
    np.savez_compressed(data_path, lattice=model.uc, kpoints=source.k_points, energy=grid,
                        dos=np.array([full_dos, selected_dos]), charge=charge,
                        shc=np.array([full_shc, selected_shc]), sampled=sampled,
                        path_distance=distance, path_ticks=ticks,
                        bands=original_bands, projected_bands=selected_bands,
                        spin=spin, valence_energy=valence_energy)
    receipt_path.write_text(json.dumps(receipt, indent=2)+"\n")
    phase("numerical data saved")

receipt = json.loads(receipt_path.read_text())
with np.load(data_path) as saved:
    data = {key: saved[key] for key in saved.files}
receipt["spin_texture_max_difference"] = float(np.nanmax(np.abs(data["spin"][1]-data["spin"][0])))
receipt["unresolved_spin_points"] = int(np.sum(~np.all(np.isfinite(data["spin"][0]), axis=1)))
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.tri import Triangulation
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.path import Path as PlotPath
from matplotlib.lines import Line2D
from matplotlib.ticker import MaxNLocator
sys.path.insert(0, str(args.repo/"examples"))
import prl_style

prl_style.use("double", height=5.8)
blue, gray, orange = "#0072B2", "#999999", "#D55E00"
fig = plt.figure()
gs = fig.add_gridspec(2, 4, width_ratios=[2.5, 1., 1.15, 1.25],
                     height_ratios=[1.45, 1.], left=.09, right=.91,
                     bottom=.095, top=.93, hspace=.32, wspace=.17)
axes = [fig.add_subplot(gs[0, 0])]
axes += [fig.add_subplot(gs[0, i], sharey=axes[0]) for i in range(1, 4)]
energy = data["energy"]-reference_mu
axes[0].plot(data["path_distance"], data["bands"]-reference_mu, color=gray, lw=.85)
axes[0].plot(data["path_distance"], data["projected_bands"]-reference_mu,
             color=blue, ls="--", lw=.65)
axes[0].set_xticks(data["path_ticks"], receipt["band_path_labels"])
axes[0].set_xlim(data["path_distance"][[0, -1]])
axes[0].set_ylabel(r"$(E-\mu_{\mathrm{ref}})\ [\mathrm{eV}]$")
axes[0].set_xlabel(r"$\mathrm{Wave\ vector}$")
for tick in data["path_ticks"][1:-1]:
    axes[0].axvline(tick, color=".8", lw=.5, zorder=0)
for index, values in enumerate((data["dos"], data["charge"], data["shc"]), 1):
    for j, (color, line) in enumerate(((gray, "-"), (blue, "--"))):
        axes[index].plot(values[j], energy, color=color, ls=line, lw=1.15)
    axes[index].tick_params(labelleft=False)
    axes[index].xaxis.set_major_locator(MaxNLocator(2))
axes[1].set_xlabel(r"$\rho\ [\mathrm{eV}^{-1}\mathrm{cell}^{-1}]$")
axes[2].set_xlabel(r"$\sigma_{xx}\ [e^2/h]$")
axes[3].set_xlabel(r"$\sigma^{z}_{xy}\ [e^2/h]$")
sampled = data["sampled"]
for i, ax in enumerate(axes):
    ax.set_ylim(-3, 3)
    ax.axhline(0, color="black", lw=.5, ls=":")
    prl_style.panel(ax, "$("+chr(97+i)+")$")
axes[0].text(.97, .96, label, transform=axes[0].transAxes, ha="right", va="top")
fig.legend([Line2D([], [], color=gray, lw=1.4),
            Line2D([], [], color=blue, ls="--", lw=1.4),
            Line2D([], [], color=orange, ls=":", lw=1.4)],
           [r"$\mathrm{Full}$", r"$\mathrm{Projected}$", r"$\mathrm{KPM}\ M=1024\ (d\ \mathrm{inset})$"],
           loc="upper center", ncol=3, bbox_to_anchor=(.53, 1.005))

# Full original energy range remains visible despite the gap-centred main view.
inset = axes[0].inset_axes([.07, .045, .28, .31])
inset.plot(data["path_distance"], data["bands"]-reference_mu, color=gray, lw=.5)
inset.plot(data["path_distance"], data["projected_bands"]-reference_mu, color=blue, ls="--", lw=.5)
inset.set_xticks([])
inset.yaxis.set_major_locator(MaxNLocator(2))
inset.tick_params(labelsize=8)
inset.set_facecolor("white")

# The overview covers the full energy-dependent response. A separate gap
# inset exposes the plateau discrepancy without letting noisy KPM tails
# dictate the scale of the physical projection comparison.
detail = axes[3].inset_axes([.16, .39, .78, .30])
for j, (color, line) in enumerate(((gray, "-"), (blue, "--"))):
    detail.plot(data["shc"][j], energy, color=color, ls=line, lw=1.)
detail.plot(sampled[:, 1], sampled[:, 0]-reference_mu, color=orange, ls=":", lw=1.)
gap_half = min(.15, .8*min(reference_mu-receipt["gap_eV"][0],
                         receipt["gap_eV"][1]-reference_mu))
gap_half = max(gap_half, .01)
values = np.array(list(receipt["shc_at_reference"].values()))
margin = max(.08, .35*float(np.ptp(values)))
detail.set_xlim(values.min()-margin, values.max()+margin)
detail.set_ylim(-gap_half, gap_half)
detail.xaxis.set_major_locator(MaxNLocator(2))
detail.yaxis.set_major_locator(MaxNLocator(2))
detail.tick_params(labelsize=8)
detail.set_facecolor("white")

# Tile the reciprocal primitive mesh and retain points in its Wigner-Seitz cell.
reciprocal = 2*np.pi*np.linalg.inv(data["lattice"]).T
sites = np.array([i*reciprocal[0, :2]+j*reciprocal[1, :2]
                  for i in range(-2, 3) for j in range(-2, 3)])
vor = Voronoi(sites)
centre = np.argmin(np.linalg.norm(sites, axis=1))
polygon = vor.vertices[vor.regions[vor.point_region[centre]]]
order = np.argsort(np.arctan2(polygon[:, 1], polygon[:, 0]))
polygon = polygon[order]
boundary = np.vstack([polygon, polygon[:1]])
path = PlotPath(boundary)
points, indices = [], []
for i in (-1, 0, 1):
    for j in (-1, 0, 1):
        p = data["kpoints"][:, :2]+i*reciprocal[0, :2]+j*reciprocal[1, :2]
        choose = path.contains_points(p, radius=1e-8)
        points.extend(p[choose]); indices.extend(np.flatnonzero(choose))
points, indices = np.array(points), np.array(indices)
_, unique = np.unique(np.round(points, 10), axis=0, return_index=True)
points, indices = points[unique], indices[unique]
spin_axes = [fig.add_subplot(gs[1, :2]), fig.add_subplot(gs[1, 2:])]
cmap = LinearSegmentedColormap.from_list("spin_blue_white_orange", [blue, "#ffffff", orange])
for branch, ax in enumerate(spin_axes):
    values = data["spin"][branch, indices]
    finite = np.all(np.isfinite(values), axis=1)
    if np.count_nonzero(finite) >= 3:
        triangles = Triangulation(points[:, 0], points[:, 1])
        triangles.set_mask(np.any(~finite[triangles.triangles], axis=1))
        field = ax.tripcolor(triangles, np.nan_to_num(values[:, 2]),
                             shading="gouraud", cmap=cmap, vmin=-1, vmax=1, rasterized=True)
        choose = finite & (indices % 2 == 0) & ((indices//24) % 2 == 0)
        ax.quiver(points[choose, 0], points[choose, 1], values[choose, 0], values[choose, 1],
                  color="black", angles="xy", scale_units="xy", scale=8, width=.005)
    ax.plot(boundary[:, 0], boundary[:, 1], color="black", lw=.7)
    ax.set_aspect("equal")
    ax.set_xlabel(r"$k_x\ [\mathrm{\AA}^{-1}]$")
    ax.set_ylabel(r"$k_y\ [\mathrm{\AA}^{-1}]$")
    ax.xaxis.set_major_locator(MaxNLocator(3))
    ax.yaxis.set_major_locator(MaxNLocator(3))
    text = r"$(e)\quad\mathrm{Full\ valence}$" if branch == 0 else r"$(f)\quad\mathrm{Projected\ valence}$"
    prl_style.panel(ax, text, bbox=dict(facecolor="white", edgecolor="none", alpha=.85, pad=1))
cax = fig.add_axes([.93, .13, .014, .24])
bar = fig.colorbar(field, cax=cax, ticks=[-1, 0, 1])
bar.set_label(r"$\langle\sigma_z\rangle$", labelpad=1)
stem = args.out/("fig_"+material+"_material_review_double")
for ext in ("png", "svg"):
    fig.savefig(str(stem)+"."+ext, dpi=600, bbox_inches="tight", pad_inches=0,
                bbox_extra_artists=[fig.patch])
from PIL import Image
with Image.open(str(stem)+".png") as raster:
    export_pixels = list(raster.size)
if abs(export_pixels[0]-4200) > 1:
    raise RuntimeError(f"figure export width must be seven inches at 600 dpi: {export_pixels}")
caption = (
    f"FIG. {1 if material == 'bi' else 2}. {material.upper()} clean SIESTA material review. "
    "Solid gray denotes the full Hamiltonian and operators; dashed blue denotes the retained "
    "spectral subspace of the implemented shift projection, with absolute window [-5,5] eV. "
    "(a) Bands on the indicated reciprocal-space path, main view within 3 eV of the reference; "
    "the inset contains the full original energy range. Retained bands overlap the gray curves. "
    "(b) Exact Gaussian DOS per cell with projected identity, excluding shifted complement states. "
    "(c) Exact finite-width Kubo-Greenwood longitudinal conductivity, including intraband terms. "
    "(d) Exact spin Hall response with Pauli spin, in e^2/h; the inset zooms into the gap and "
    "includes the dotted-orange saved projected M=1024, R=16, seed=12345, Jackson 10 meV KPM "
    "result, which uses a different kernel. "
    f"Band energies and response chemical potentials are measured from mu_ref={reference_mu:g} eV "
    "relative to stored SIESTA Fermi zero; for PdSe2 this selects the previously studied upper gap. "
    "(e,f) Spin texture of the highest occupied valence band below mu_ref, original and explicitly "
    "projected, respectively. Where energy-degenerate, the unique maximal eigenstate of physical "
    "sigma_z within the multiplet is selected; an unresolved spin eigenspace is left blank. "
    "Blue-to-orange fill encodes sigma_z from -1 to +1; black arrows encode "
    "(sigma_x,sigma_y), with one Pauli unit equal to 0.125 inverse Angstrom on the wave-vector "
    "axes. The black outline is the first Brillouin zone. Degenerate spin-resolved states "
    "depend on the stated spin quantisation convention, not an arbitrary energy eigenvector. "
    "The retained projected spin agrees with the original within the error in the JSON receipt. "
    "Parameters: 24x24 Gamma-centred mesh, W=0, exact Lowdin metric, covariant physical velocity, "
    "Gaussian standard deviations 50 meV for DOS/longitudinal conductivity and 10 meV for the "
    "spin Hall occupation. The charge curve depends on artificial width and mesh and is not a "
    "converged DC transport result. This is a review figure, not a stage-closure certificate."
)
Path(str(stem)+".caption.txt").write_text(caption+"\n")
receipt["figure_style"] = dict(width_inches=float(fig.get_figwidth()), height_inches=float(fig.get_figheight()),
                               export_pixels=export_pixels,
                               dpi=600, renderer="latex" if matplotlib.rcParams["text.usetex"] else "mathtext-cm",
                               tick_font_points=matplotlib.rcParams["xtick.labelsize"],
                               minimum_line_points=.5, redundant_curves="gray solid / blue dashed / orange dotted",
                               caption=True, no_titles=all(not a.get_title() for a in fig.axes),
                               scientific_scope="diagnostic; no transport convergence claim")
receipt["figure_script_sha256"] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
receipt_path.write_text(json.dumps(receipt, indent=2)+"\n")
phase("figure exported")
