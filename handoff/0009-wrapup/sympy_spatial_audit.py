"""SymPy spatial-symmetry audit requested by the user; no material symmetry assumed."""
import argparse
import hashlib
import json
from pathlib import Path
import sympy as sp

parser=argparse.ArgumentParser()
parser.add_argument("--out",type=Path,required=True)
args=parser.parse_args()
args.out.mkdir(parents=True,exist_ok=True)
out={"sympy_version":sp.__version__}
R=sp.Matrix([[0,-1,0],[1,0,0],[0,0,1]])
mirror=sp.diag(-1,1,1)
inversion=-sp.eye(3)
M=sp.Matrix(3,3,sp.symbols("m0:9"))
rotations=[R**k for k in range(4)]
avg=sp.simplify(sum((r*M*r.T for r in rotations),sp.zeros(3))/4)
assert sp.simplify(R*avg*R.T-avg)==sp.zeros(3)
assert avg[0,0]==avg[1,1] and sp.simplify(avg[0,1]+avg[1,0])==0
out["C4_charge_group_projection"]=str(avg)
out["C4_projection_idempotent"]=sp.simplify(sum((r*avg*r.T for r in rotations),sp.zeros(3))/4-avg)==sp.zeros(3)
mirrored=sp.simplify((M+mirror*M*mirror.T)/2)
assert mirrored[0,1]==0 and mirrored[1,0]==0
out["mirror_x_charge_projection"]=str(mirrored)

def component_signs(r):
    polar=list(r.diagonal())
    axial=[r.det()*v for v in polar]
    return {"Vx":str(polar[0]),"Vy":str(polar[1]),"Sz":str(axial[2]),
            "Jxz":str(polar[0]*axial[2]),
            "charge_xy":str(polar[0]*polar[1]),
            "spin_hall_xzy":str(polar[0]*axial[2]*polar[1])}
out["mirror_x_parities"]=component_signs(mirror)
out["inversion_parities"]=component_signs(inversion)
assert out["mirror_x_parities"]["charge_xy"]=="-1"
assert out["mirror_x_parities"]["spin_hall_xzy"]=="1"
assert out["inversion_parities"]["spin_hall_xzy"]=="1"
assert out["inversion_parities"]["charge_xy"]=="1"

# Antiunitary selection expressed in real and imaginary components.
real,imag=sp.symbols("real imag",real=True)
z=real+sp.I*imag
out["TR_charge_constraint"]=str(sp.simplify(z-sp.conjugate(z)))
out["TR_spin_hall_constraint"]=str(sp.simplify(z+sp.conjugate(z)))
assert sp.simplify(z-sp.conjugate(z))==2*sp.I*imag
assert sp.simplify(z+sp.conjugate(z))==2*real

# Counterexample: exact mirror-odd trace vanishes; a finite phase vector does not.
U=sp.diag(1,-1)
H=sp.diag(sp.Rational(1,3),sp.Rational(2,3))
A=sp.Matrix([[0,1],[1,0]])
r=sp.Matrix([1,1])/sp.sqrt(2)
assert U*H*U.H==H and U*A*U.H==-A
raw=(r.H*A*r)[0]
assert sp.trace(A)==0 and raw==1
orbit=(r*r.H+U*r*r.H*U.H)/2
assert sp.trace(orbit*A)==0
out["finite_random_state_counterexample"]={"H":str(H),"U":str(U),"A":str(A),
 "r":str(r),"exact_trace":str(sp.trace(A)),"single_sample":str(raw),
 "orbit_averaged_sample":str(sp.trace(orbit*A))}

# Spatial covariance preserves polynomial index. Chiral antisymmetry is distinct.
x=sp.symbols("x",real=True)
for n in range(13):
    assert sp.expand(sp.chebyshevt(n,-x)-(-1)**n*sp.chebyshevt(n,x))==0
out["chiral_parity_polynomial_checks"]="T_n(-x)=(-1)^n T_n(x), n=0..12; general parity follows from Chebyshev recurrence."
center,E0,width,h=sp.symbols("center E0 width h",real=True,nonzero=True)
out["chiral_offcenter_residual"]=str(sp.simplify((2*E0-h-center)/width+(h-center)/width))
assert sp.simplify((2*E0-h-center)/width+(h-center)/width-2*(E0-center)/width)==0

# Rotational tensor averaging includes components, not just scalar orbit weights.
xy_orbit=[sp.expand((r*M*r.T)[0,1]) for r in rotations]
out["C4_xy_orbit"]=list(map(str,xy_orbit))
out["C4_xy_sum"]=str(sp.expand(sum(xy_orbit)))
assert sp.expand(sum(xy_orbit))==2*(M[0,1]-M[1,0])
out["scope"]="Group/tensor identities and counterexamples only; no symmetry certification of the frozen Bi/PdSe2 models."
out["audit_sha256"]=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
(args.out/"sympy-spatial-results.json").write_text(json.dumps(out,indent=2)+"\n")
print(json.dumps(out,indent=2))

