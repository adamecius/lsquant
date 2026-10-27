"""Separate finite-order and stochastic material errors from saved runs.

This tester diagnostic sums the exact projected trace in retained eigenstates.
Zero operator rows in the discarded complement contribute nothing; normalization
still uses the original Hilbert dimension. It does not replace the production
model, generate a golden, or change a tolerance.
"""
import argparse
import importlib.util
import json
from pathlib import Path
import time
import numpy as np
from scipy.special import erfc
from lsquant import _core
from lsquant.bases import BlochMetricBasis
from lsquant.kernels import jackson_damping
from lsquant.spectral_functions import legacy_delta_table, legacy_green_tables

parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument("--root",type=Path,required=True)
parser.add_argument("--material",choices=("bi","pdse2"),required=True)
parser.add_argument("--out",type=Path,required=True)
parser.add_argument("--max-moments",type=int,default=1024)
parser.add_argument("--folded-convergence",action="store_true")
args=parser.parse_args()
if args.max_moments < 1024:
    parser.error("--max-moments must be at least 1024")
args.root=args.root.resolve()
args.out.mkdir(parents=True,exist_ok=True)
spec=importlib.util.spec_from_file_location("projection_study",args.root/"projection-comparison.py")
harness=importlib.util.module_from_spec(spec)
spec.loader.exec_module(harness)
harness.FROZEN_MODELS=args.root/"projection-reference-study"/"frozen-models"
old_dir=args.root/"reduced-validation"/"materials"/(args.material+"-covariant-k24-shift5")
old=json.loads((old_dir/"run.json").read_text())
records=[]
def record(**values):
    records.append(values)
    (args.out/"error-budget.json").write_text(json.dumps(records,indent=2)+"\n")
    print(json.dumps(values),flush=True)
model=harness.load_model(args.material,"covariant")
if model.metadata["frozen_model_sha256"]!=old["model_sha256"]:
    raise RuntimeError("Frozen input checksum mismatch")
record(phase="source_started",material=args.material,model_sha256=old["model_sha256"])
source=BlochMetricBasis(model,[24,24,1],velocity_mode="covariant",route="lowdin",seed=12345)
area=float(np.linalg.norm(np.cross(model.uc[0],model.uc[1])))
probes=np.array(old["energies"])
current=source.operator_blocks("JXSZ")
velocity=source.operator_blocks("VY")
blocks=[]
curvatures=[]
oracle=np.zeros(len(probes))
for k,(energies,vectors) in enumerate(zip(source.eigenvalues,source.eigenvectors)):
    selected=(energies>=-5)&(energies<=5)
    e=energies[selected]
    c=vectors[:,selected]
    j=c.conj().T@current[k]@c
    v=c.conj().T@velocity[k]@c
    j=.5*(j+j.conj().T)
    v=.5*(v+v.conj().T)
    pair=j*v.T
    groups=np.cumsum(np.r_[True,np.diff(e)>1e-6])
    gaps=e[:,None]-e[None,:]
    inverse=np.zeros_like(gaps)
    np.divide(1.,gaps*gaps,out=inverse,where=groups[:,None]!=groups[None,:])
    omega=-2*np.imag(np.sum(pair*inverse,axis=1))
    curvatures.append((e,omega))
    oracle+=np.sum(.5*erfc((e[None,:]-probes[:,None])/(np.sqrt(2)*.01))*omega,axis=1)
    blocks.append((e,pair))
oracle*=2*np.pi/(source.num_cells*area)
record(phase="oracle_crosscheck",exact_selected_trace=oracle.tolist(),
       previous_projected_oracle=old["projected_oracle_10meV"],
       max_difference=float(np.max(np.abs(oracle-np.array(old["projected_oracle_10meV"])))),
       min_kept=min(len(e) for e,p in blocks),max_kept=max(len(e) for e,p in blocks))
saved=_core.ReadChebmom(str(old_dir/"M1024.chebmom"))
window=saved.window
begin=time.monotonic()
left=[]
right=[]
for e,pair in blocks:
    x=e*window.scale+window.shift
    cheb=np.empty((args.max_moments,len(e)))
    cheb[0]=1.
    cheb[1]=x
    for m in range(2,args.max_moments):
        cheb[m]=2*x*cheb[m-1]-cheb[m-2]
    left.append(cheb)
    right.append(cheb@pair)
left=np.concatenate(left,axis=1)
right=np.concatenate(right,axis=1)
mu=left@right.T/source.dimension
weights=np.full(args.max_moments,2.)
weights[0]=1.
mu*=weights[:,None]*weights[None,:]
record(phase="exact_moments",seconds=time.monotonic()-begin,
       weighted_hermiticity_max=float(np.max(np.abs(mu-mu.conj().T))),
       original_dimension=source.dimension)
np.savez_compressed(args.out/"exact-moments.npz",moments=mu,
                     window=[window.lower,window.upper],dimension=source.dimension)
if args.folded_convergence:
    from types import SimpleNamespace
    from lsquant.reconstruct import shc
    for count in (256,1024,2048,4096,8192):
        if count > args.max_moments:
            continue
        for broadening in (10.,1.):
            begin=time.monotonic()
            target=np.zeros(len(probes))
            for e,omega in curvatures:
                target+=np.sum(.5*erfc((e[None,:]-probes[:,None])/
                                      (np.sqrt(2)*broadening/1000))*omega,axis=1)
            target*=2*np.pi/(source.num_cells*area)
            tensor=SimpleNamespace(rank=2,dims=[count,count],data=mu[:count,:count],
                                   window=window,system_size=source.dimension)
            energy,curve=shc(tensor,broadening,source.num_cells,area,backend="folded")
            values=np.interp(probes,energy,curve)
            effective=min(count,int(np.ceil(np.pi*(window.upper-window.lower)/(2*broadening/1000))))
            record(phase="folded_convergence",moments=count,effective_moments=effective,
                   broadening_meV=broadening,energies=probes.tolist(),
                   exact_trace_kpm=values.tolist(),gaussian_oracle=target.tolist(),
                   difference=(values-target).tolist(),seconds=time.monotonic()-begin)
    record(phase="complete",status="convergence diagnostic; kernels differ; no gate changed")
    raise SystemExit(0)
for count in (256,1024):
    begin=time.monotonic()
    xfull=np.array([-.95+i*1.9/(30*count-1) for i in range(30*count)])
    efull=(xfull-window.shift)/window.scale
    end=min(len(xfull),int(np.searchsorted(efull,probes.max()))+3)
    x=xfull[:end]
    g=jackson_damping(count)
    damped=np.ascontiguousarray(mu[:count,:count]*g[:,None]*g[None,:])
    delta=legacy_delta_table(x,count)
    _,derivative=legacy_green_tables(x,count,count)
    out=_core.ContractRank2Imag(damped,delta,derivative)
    kernel=-2*np.asarray(out)*source.dimension*window.scale**2
    curve=np.cumsum(kernel[:-1]*np.diff(x))*(2*np.pi/(source.num_cells*area))
    exact=np.interp(probes,efull[:end-1],curve)
    sampled=np.array(next(r["kpm"] for r in old["results"] if r["moments"]==count))
    record(phase="error_budget",moments=count,energies=probes.tolist(),
           exact_trace_kpm=exact.tolist(),sampled_R16_kpm=sampled.tolist(),
           projected_oracle=oracle.tolist(),
           finite_order_error=(exact-oracle).tolist(),
           sample_error_at_same_order=(sampled-exact).tolist(),
           total_error=(sampled-oracle).tolist(),seconds=time.monotonic()-begin)
record(phase="complete",status="diagnosis; no gate or tolerance changed")
