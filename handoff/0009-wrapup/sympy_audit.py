"""User-requested symbolic audit; diagnostic only, not a new golden gate.

Run on TCN with candidate-v5/python and its compiled extension on PYTHONPATH.
No production sources, reference data, seeds or tolerances are modified.
"""
import argparse
import hashlib
import inspect
import json
from pathlib import Path
import sympy as sp
import numpy as np
from lsquant import metric_forms as mf
from lsquant import reconstruct

parser = argparse.ArgumentParser()
parser.add_argument("--out", type=Path, required=True)
args = parser.parse_args()
args.out.mkdir(parents=True, exist_ok=True)
results = {}

def record(name, residual, scope):
    value = sp.simplify(residual)
    assert value == 0, (name, value)
    results[name] = dict(residual=str(value), scope=scope)

a, b = sp.symbols("a b", real=True)
record("chebyshev_product_all_angles",
       2*sp.cos(a)*sp.cos(b)-sp.cos(a+b)-sp.cos(a-b),
       "For all real angles; substitute a=m theta,b=n theta and x=cos(theta).")
x = sp.symbols("x", real=True)
for m in range(1, 13):
    record(f"dos_even_{m}", sp.chebyshevt(2*m,x)-2*sp.chebyshevt(m,x)**2+1,
           "Polynomial corroboration; same-seed Hermitian identity-DOS shortcut.")
    record(f"dos_odd_{m}", sp.chebyshevt(2*m+1,x)-2*sp.chebyshevt(m,x)*sp.chebyshevt(m+1,x)+x,
           "Polynomial corroboration; raw unweighted moments only.")

# Generic formal trace identities: no commutation of static/spectral slots.
A, B = mf.A, mf.B
X, Y = sp.symbols("F_m F_n", commutative=False)
record("trace_adjoint_exchange", mf.normal_form(Y*B*X*A-A*Y*B*X, trace=True),
       "A,B,F_m,F_n Hermitian; conjugation reverses factors, trace is cyclic.")
record("equal_static_trace_real", mf.normal_form(Y*A*X*A-A*X*A*Y, trace=True),
       "Equal Hermitian statics: exact-trace moments are real; not generic SHC.")

sx=sp.Matrix([[0,1],[1,0]])
sy=sp.Matrix([[0,-sp.I],[sp.I,0]])
sz=sp.diag(1,-1)
I=sp.eye(2)
H=sz/2
AA=I+sx
BB=I+sy
r=sp.Matrix([1,sp.I])/sp.sqrt(2)
def mean(matrix):
    return sp.simplify((r.conjugate().T*matrix*r)[0])
delta=sp.simplify(mean(AA*BB*H)-sp.conjugate(mean(AA*H*BB)))
assert delta != 0
results["stochastic_triangle_counterexample"] = dict(
    A=str(AA), B=str(BB), H=str(H), r=str(r),
    mu01_minus_conjugate_mu10=str(delta),
    implication="Exact trace Hermiticity is not a per-vector identity. Mirroring one sampled triangle changes the estimator.")
bad=sp.trace(sx*H*sx*H)-sp.trace(sx*sx*H*H)
assert bad != 0
results["illegal_rank2_collapse"] = dict(residual=str(bad),
    implication="Do not merge Chebyshev factors across a static operator that does not commute with H.")

# General retarded derivative and folded expression, interior of Chebyshev domain.
theta=sp.symbols("theta", positive=True)
n=sp.symbols("n", integer=True, nonnegative=True)
s=sp.sin(theta)
g=-sp.I*sp.exp(-sp.I*n*theta)/s
dg=-sp.I*sp.exp(-sp.I*n*theta)*(sp.cos(theta)+sp.I*n*s)/s**3
record("retarded_derivative", -sp.diff(g,theta)/s-dg,
       "0<theta<pi; d/dx=-(1/sin(theta))*d/dtheta.")
record("folded_frequency_identity",
       sp.expand_complex(sp.cos(a)*sp.exp(-sp.I*b)-
       (sp.exp(-sp.I*(b-a))+sp.exp(-sp.I*(b+a)))/2),
       "All real a,b: combines frequencies n-m and n+m.")
ar,ai,br,bi,cr,ci,dr,di,q = sp.symbols("ar ai br bi cr ci dr di q", real=True)
F=ar+sp.I*ai+sp.I*q*(br+sp.I*bi)
Fn=cr+sp.I*ci+sp.I*q*(dr+sp.I*di)
record("folded_imaginary_numerator",
       sp.im(-sp.I*(x*F+sp.I*q*Fn))-(-x*ar+q*q*dr+q*(x*bi+ci)),
       "q=sqrt(1-x^2); divide by pi*q^4. Arbitrary complex moments.")
for j in range(13):
    series=2*sum(sp.chebyshevt(k,x) for k in range(j,-1,-2))-(1 if j%2==0 else 0)
    record(f"U_to_T_{j}", sp.chebyshevu(j,x)-series,
           "Finite-degree corroboration of parity-tail U-to-T conversion.")
record("U_to_T_inductive_step",
       sp.trigsimp((sp.sin((n+1)*theta)-sp.sin((n-1)*theta))/sp.sin(theta)-2*sp.cos(n*theta)),
       "General induction step U_n=2*T_n+U_(n-2); bases U_0=T_0,U_1=2*T_1.")

# Record the audited production implementation.
source=inspect.getsource(reconstruct._folded_bastin_contraction)
results["folded_source_sha256"]=hashlib.sha256(source.encode()).hexdigest()
# Independent numerical corroboration of the actual implementation on complex
# rectangular tensors; symbolic identities above establish the algebra.
from lsquant import _core
from lsquant.spectral_functions import legacy_delta_table, legacy_green_tables
rng=np.random.default_rng(20261025)
errors=[]
for shape in [(1,1),(3,5),(5,3),(8,8)]:
    mu=rng.normal(size=shape)+1j*rng.normal(size=shape)
    nodes=np.array([-.93,-.5,0.,.21,.91])
    delta_table=legacy_delta_table(nodes,shape[0])
    _,derivative=legacy_green_tables(nodes,shape[0],shape[1])
    reference=_core.ContractRank2Imag(mu,delta_table,derivative)
    folded=reconstruct._folded_bastin_contraction(mu,nodes)
    errors.append(dict(shape=shape,max_abs=float(np.max(np.abs(reference-folded))),
                       relative_to_peak=float(np.max(np.abs(reference-folded))/max(1.,np.max(np.abs(reference))))))
results["actual_folded_numeric_corroboration"]=errors

# Exact symbolic reduction of a spectral selected block; all matrix entries of
# the two retained operators are independent, not assumed commuting.
e0,e1,ec=sp.symbols("e0 e1 ec", real=True)
f0,f1,fc,g0,g1,gc=sp.symbols("f0 f1 fc g0 g1 gc")
P=sp.diag(1,1,0)
Q=sp.Matrix([[1,0],[0,1],[0,0]])
O=sp.Matrix(3,3,sp.symbols("o0:9"))
V=sp.Matrix(3,3,sp.symbols("v0:9"))
Fdiag=sp.diag(f0,f1,fc)
Gdiag=sp.diag(g0,g1,gc)
record("projected_trace_compression",
       sp.trace(P*O*P*Fdiag*P*V*P*Gdiag)-
       sp.trace((Q.T*O*Q)*sp.diag(f0,f1)*(Q.T*V*Q)*sp.diag(g0,g1)),
       "Generic retained 2x2 block embedded in 3x3; proof generalizes by block multiplication. Independent complement functions cancel.")
record("projected_identity_compression",
       sp.trace(P*Fdiag)-(f0+f1),
       "Explicit projected identity required; unprojected DOS includes fc.")
rv=sp.Matrix(sp.symbols("r0:3", complex=True))
record("mapped_sample_compression",
       (rv.conjugate().T*P*O*P*Fdiag*P*V*P*Gdiag*rv)[0]-
       ((Q.T*rv).conjugate().T*(Q.T*O*Q)*sp.diag(f0,f1)*(Q.T*V*Q)*sp.diag(g0,g1)*(Q.T*rv))[0],
       "Same mapped vector, no renormalization or fresh random stream; original trace prefactor unchanged.")
spin_loss=P*(O*V+V*O)*P/2-((P*O*P)*(P*V*P)+(P*V*P)*(P*O*P))/2
expected=(P*O*(sp.eye(3)-P)*V*P+P*V*(sp.eye(3)-P)*O*P)/2
assert sp.simplify(spin_loss-expected)==sp.zeros(3)
results["projected_spin_current_missing_term"]=str(expected)
D=sp.Matrix([[0,0,1],[0,0,0],[1,0,0]])
assert P*D-D*P != sp.zeros(3)
results["disorder_commutator_counterexample"]=str(P*D-D*P)

# Existing arbitrary-word metric proofs, interpreted with their documented
# assumptions (positive overlap; spin-scalar overlap where spin is commuted).
metric=mf.trace_equivalence()
assert metric["hamiltonian_similarity"]
for name,item in metric.items():
    if isinstance(item,dict) and "equal" in item:
        assert item["equal"], name
results["existing_metric_trace_proof"]={k:str(v) for k,v in metric.items()}
projection=mf.spectral_projection_proof()
assert all(sp.simplify(value)==0 for value in projection.values())
results["existing_spectral_projection_proof"]={k:str(v) for k,v in projection.items()}

# Exact eigenbasis factorization of rank-2 moments.
t0,t1,u0,u1=sp.symbols("t0 t1 u0 u1",real=True)
AE=sp.Matrix(2,2,sp.symbols("a0:4"))
BE=sp.Matrix(2,2,sp.symbols("b0:4"))
direct=sp.trace(AE*sp.diag(t0,t1)*BE*sp.diag(u0,u1))
factor=sum(AE[i,j]*BE[j,i]*[t0,t1][j]*[u0,u1][i] for i in range(2) for j in range(2))
record("exact_eigen_pair_factorization",direct-factor,
       "mu_mn=sum_ij A_ij B_ji T_m(e_j) T_n(e_i); no stochastic estimator required for this exact block route.")
results["sympy_version"]=sp.__version__
results["audit_sha256"]=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
results["scope"]="Symbolic audit plus small deterministic corroborations; not a golden gate or performance measurement."
(args.out/"sympy-results.json").write_text(json.dumps(results,indent=2)+"\n")
print(json.dumps(dict(symbolic_zero_checks=sum(isinstance(v,dict) and v.get("residual")== "0" for v in results.values()),
    folded_numeric=errors,stochastic_counterexample=str(delta),output=str(args.out/"sympy-results.json")),indent=2))



