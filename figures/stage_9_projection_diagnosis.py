"""Diagnose the existing stage-6d failure using saved and dense exact moments.

This is a tester figure/diagnostic, not a replacement gate or a new self-check.
It keeps the production source, existing assertion and golden data unchanged.
"""
import argparse
import json
from pathlib import Path
import sys
import time
import numpy as np
from lsquant import _core
from lsquant.bases import BlochMetricBasis
from lsquant.apps.stage6b_siesta import SHC_WORD
from lsquant.projection.bases import ProjectedBlochBasis
from lsquant.projection.runfile import exact_word_curve
from lsquant.kernels import jackson_damping
from lsquant.spectral_functions import legacy_delta_table, legacy_green_tables

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--repo", type=Path, required=True)
parser.add_argument("--capture", type=Path, required=True)
parser.add_argument("--out", type=Path, required=True)
args = parser.parse_args()
sys.path.insert(0, str(args.repo/"tests"/"dev"))
from stage_6d_check import augmented_kane_mele
args.out.mkdir(parents=True, exist_ok=True)
records = []
def persist():
    (args.out/"diagnosis.json").write_text(json.dumps(records, indent=2)+"\n")
def record(**values):
    records.append(values)
    persist()
    print(json.dumps(values), flush=True)

def direct_moments(basis, count):
    """Evaluate Tr[J T_m(H) V T_n(H)] directly in the eigenbasis."""
    eigenvalues, vectors = basis.eigenvalues, basis.eigenvectors
    adjoint = np.swapaxes(vectors.conj(), -1, -2)
    current = adjoint @ basis.operator_blocks("JXSZ") @ vectors
    velocity = adjoint @ basis.operator_blocks("VY") @ vectors
    pair = current * np.swapaxes(velocity, -1, -2)
    window = basis.window()
    x = eigenvalues*window.scale + window.shift
    cheb = np.empty((count, *x.shape))
    cheb[0] = 1.
    cheb[1] = x
    for m in range(2, count):
        cheb[m] = 2*x*cheb[m-1]-cheb[m-2]
    right = np.einsum("kab,mka->mkb", pair, cheb, optimize=True)
    raw = cheb.reshape(count,-1) @ right.reshape(count,-1).T
    weights = np.full(count, 2.)
    weights[0] = 1.
    weighted = raw/basis.dimension * weights[:,None]*weights[None,:]
    return weighted, window

def plateau(mu, window, dimension, cells, area, count, grid_factor=30):
    """Use delivered coefficient tables/contraction on the needed grid prefix.

    The rectangle rule and energy labels are those of kubo_bastin. Points
    above the largest requested energy cannot affect this cumulative integral.
    """
    probes=np.array([-.1,0.,.1])
    full_grid=np.array([-.95+i*1.9/(grid_factor*count-1)
                        for i in range(grid_factor*count)])
    energy=(full_grid-window.shift)/window.scale
    end=min(len(full_grid),int(np.searchsorted(energy,probes.max()))+3)
    x=full_grid[:end]
    kernel=jackson_damping(count)
    damped=np.ascontiguousarray(mu[:count,:count]*kernel[:,None]*kernel[None,:])
    delta=legacy_delta_table(x,count)
    _, derivative=legacy_green_tables(x,count,count)
    contracted=_core.ContractRank2Imag(damped,delta,derivative)
    integrand=-2*np.asarray(contracted)*dimension*window.scale**2
    curve=np.cumsum(integrand[:-1]*np.diff(x))*(2*np.pi/(cells*area))
    return np.interp(probes,energy[:end-1],curve).tolist()

source=BlochMetricBasis(augmented_kane_mele(),[12,12,1],
                        velocity_mode="covariant",route="lowdin")
area=float(np.linalg.norm(np.cross(source.model.uc[0],source.model.uc[1])))
probes=np.array([-.1,0.,.1])
oracle,_=exact_word_curve(source,SHC_WORD,probes,broadening=.01,cell_area=area)
saved=_core.ReadChebmom(str(args.capture/"moments-05.chebmom"))
direct64,window=direct_moments(source,64)
record(phase="engine_vs_direct",moments=64,
       max_abs=float(np.max(np.abs(direct64-saved.data[:64,:64]))),
       max_abs_transposed=float(np.max(np.abs(direct64.T-saved.data[:64,:64]))),
       window=[window.lower,window.upper],oracle=oracle.tolist())
for count in (256,512,1024):
    begin=time.monotonic()
    values=plateau(saved.data,saved.window,saved.system_size,source.num_cells,area,count)
    record(phase="saved_wide_convergence",moments=count,plateau=values,
           max_error=float(np.max(np.abs(np.array(values)-oracle))),
           nominal_resolution_meV=float(np.pi*(window.upper-window.lower)/(2*count)*1000),
           seconds=time.monotonic()-begin)
for tag,basis in (("wide",source),
                  ("shift",ProjectedBlochBasis(source,dict(kind="shift",window=[-3.1,3.1])))):
    begin=time.monotonic()
    raw,window=direct_moments(basis,2048 if tag=="wide" else 1024)
    record(phase="direct_moments",model=tag,moments=len(raw),
           seconds=time.monotonic()-begin,
           window=[window.lower,window.upper])
    for count,factor in ((1024,30),(1024,60),(2048,30)) if tag=="wide" else ((256,30),(512,30),(1024,30)):
        begin=time.monotonic()
        values=plateau(raw,window,basis.dimension,basis.num_cells,area,count,factor)
        record(phase="direct_convergence",model=tag,moments=count,grid_factor=factor,
               plateau=values,max_error=float(np.max(np.abs(np.array(values)-oracle))),
               seconds=time.monotonic()-begin)

# Independently execute the existing orthogonal-tool checks after the
# baseline failure. Their observations do not make the original self-check pass.
from scipy.linalg import eigvalsh
from lsquant.projection import downfold_orthogonal,LocalMatrix,validate_bands
reduced=downfold_orthogonal(source.hk_array,source.sk_array,[-3.1,3.1])
expected=np.stack([e[mask] for e,mask in zip(reduced.original_eigenvalues,reduced.kept)])
record(phase="orthogonal_training_bands",
       max_error=float(np.max(np.abs(np.linalg.eigvalsh(reduced.hamiltonian)-expected))))
local=LocalMatrix(reduced.hamiltonian,source.model.uc,source.model.tau[reduced.trial_columns],
                  source.supercell,source.k_points)
check=validate_bands(local,lambda k:eigvalsh(source.model.Hk(k),source.model.Sk(k)),
                     [-3.1,3.1],source.model.uc,source.supercell)
record(phase="orthogonal_offgrid",result=check)

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
fig,ax=plt.subplots(figsize=(7,4.5))
for tag,label in (("wide","Original spectrum"),("shift","Shifted spectrum")):
    selected=[r for r in records if r["phase"]=="direct_convergence" and r["model"]==tag and r["grid_factor"]==30]
    if tag=="wide":
        selected=[r for r in records if r["phase"]=="saved_wide_convergence" and r["moments"]<1024]+selected
    ax.loglog([r["moments"] for r in selected],[r["max_error"] for r in selected],"o-",label=label)
ax.axhline(1e-3,color="black",linestyle="--",label="Existing self-check limit")
ax.set(xlabel="Chebyshev moments",ylabel="Maximum SHC error at the three gap probes",
       title="Stage 6d: exact-trace convergence diagnosis")
ax.grid(alpha=.2)
ax.legend()
fig.tight_layout()
fig.savefig(args.out/"selfcheck-convergence.png",dpi=180)
record(phase="complete",status="diagnosis only; original gate unchanged")
