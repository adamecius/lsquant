import sympy as sp

th, a, g, hbar, t, Jex, EF = sp.symbols('theta a gamma hbar t J_ex E_F', real=True)

print("="*70)
print("C') KG matrix: confirm thesis B.15 INCLUDING (1+d_m0)(1+d_n0) norm")
print("="*70)
V2 = (2*g*a/hbar)**2 * sp.sin(th)**2
def Rraw(mm,nn):
    integ = V2*sp.chebyshevt(mm,-sp.cos(th))*sp.chebyshevt(nn,-sp.cos(th))
    return sp.simplify(sp.integrate(integ,(th,-sp.pi,sp.pi))/(2*sp.pi))
pref = g**2*a**2/hbar**2
def closed(mm,nn):
    d=lambda i,j:1 if i==j else 0
    return sp.Rational(d(mm,nn),1+d(mm,0)+d(mm,1)) - sp.Rational(1,2)*(1-d(mm,nn))*d(abs(mm-nn),2)
ok=True
N=6
for mm in range(N+1):
    for nn in range(N+1):
        norm=(1+(1 if mm==0 else 0))*(1+(1 if nn==0 else 0))
        Cnorm = sp.simplify(Rraw(mm,nn)/pref/norm)
        if sp.simplify(Cnorm-closed(mm,nn))!=0:
            ok=False; print("  mismatch",mm,nn,Cnorm,closed(mm,nn))
print("  thesis B.15 (normalized) matches symbolic integral:", ok)

print()
print("="*70)
print("E) Spin precession, exchange chain  H_k = eps(k) I - J_ex sigma_z")
print("   Heisenberg s(t)=U^dag s U ; orbital phase cancels (k-independent)")
print("="*70)
I2=sp.eye(2)
sx=sp.Matrix([[0,1],[1,0]]); sy=sp.Matrix([[0,-sp.I],[sp.I,0]]); sz=sp.Matrix([[1,0],[0,-1]])
# U_ex(t) = exp(-i H_ex t/hbar), H_ex=-J_ex sigma_z  -> exp(+i J_ex t sigma_z/hbar)
phi = Jex*t/hbar
Uex = sp.cos(phi)*I2 + sp.I*sp.sin(phi)*sz       # = exp(i phi sigma_z)
sx_t = sp.simplify(Uex.H*sx*Uex)
sy_t = sp.simplify(Uex.H*sy*Uex)
print("  sigma_x(t) =", sp.simplify(sx_t))
print("  sigma_y(t) =", sp.simplify(sy_t))
# initial +x state
xp = sp.Matrix([1,1])/sp.sqrt(2)
Sx = sp.simplify((xp.H*sx_t*xp)[0]);  Sy = sp.simplify((xp.H*sy_t*xp)[0])
print("  <sigma_x>(t) =", Sx, "   expected cos(2 J_ex t/hbar)")
print("  <sigma_y>(t) =", Sy, "   expected -sin(2 J_ex t/hbar)")

print()
print("="*70)
print("E2) Equivalence of the two flavors for the exchange chain")
print("    (k-independent spin block => both reduce to same single-mode result)")
print("="*70)
# Density-projection: rho ~ P_x delta(H-EF) P_x ; trace over the 2-dim spin block
# at the selected k. P_x = |+x><+x|.  Spin block at any selected k is identical.
Px = xp*xp.H
# numerator Tr[U^dag sx U Px]/Tr[Px] within block, orbital phase cancels:
num_dens = sp.simplify(sp.trace(Uex.H*sx*Uex*Px))
den_dens = sp.simplify(sp.trace(Px))
Sx_dens = sp.simplify(num_dens/den_dens)
print("  density-projection <sigma_x>(t) =", Sx_dens)
# State-projection: |phi(0)> = (1/2)(I + x.s)|phi_r>; for the spin block the spin
# texture picks +x; numerator = Re<phi(t)| sx delta |phi(t)> / <phi(t)|delta|phi(t)>
# in the block this is Re<+x| U^dag sx U |+x> = <sigma_x>(t)  (already real)
Sx_state = sp.re(Sx)
print("  state-projection   <sigma_x>(t) =", sp.simplify(Sx_state))
print("  identical:", sp.simplify(Sx_dens - Sx_state)==0)

print()
print("="*70)
print("F) TOY k-dependent field: H_k = eps(k) I + B(k) sigma_z, B(k)=J_ex+delta*cos(ka)")
print("   precession freq omega(k)=2 B(k)/hbar ; isolates flavor difference")
print("="*70)
delta, EFs, w = sp.symbols('delta E_F omega', real=True)
ka = sp.symbols('ka', real=True)
Bk = Jex + delta*sp.cos(ka)
# DENSITY (sharp E_F): shell {+k0,-k0}, eps even in k => B(+k0)=B(-k0) -> single freq.
# cos(k0 a) = -E_F/(2 gamma)  => B = J_ex - delta*E_F/(2 gamma)
B_at_EF = Jex + delta*(-EFs/(2*g))
omega_dens = sp.simplify(2*B_at_EF/hbar)
print("  DENSITY: coherent precession, omega(E_F) =", omega_dens)
print("           <sigma_x>(t) = cos(omega(E_F) t)  (NO dephasing in 1D)")
# STATE (energy window): average cos(2 B(k) t/hbar) over the band with weight rho(E)dE.
# Map: E=-2 gamma cos(ka). Average of cos(2B t/hbar) over uniform-in-k (schematic) shows decay.
# Demonstrate decay numerically below.
print("  STATE: averages cos(2 B(k) t/hbar) over an energy window -> dephasing")

print()
print("="*70)
print("F2) Numerical demo of flavor difference (delta != 0)")
print("="*70)
import numpy as np
gN, JexN, dN, hbarN = 1.0, 0.1, 0.05, 1.0    # eV, eV, eV, (hbar=1 units, t in 1/eV)
Nk = 200000
kk = np.linspace(-np.pi, np.pi, Nk, endpoint=False)
epsk = -2*gN*np.cos(kk)
Bk_n = JexN + dN*np.cos(kk)
tt = np.linspace(0, 200, 5)   # few sample times
# density at sharp E_F=0: shell cos(ka)=0 -> B=J_ex, single freq
omega0 = 2*JexN/hbarN
dens = np.cos(omega0*tt)
# state: energy window around E_F=0 of half-width eta=0.3 eV, weight ~ within window
eta = 0.3
mask = np.abs(epsk-0.0) < eta
state = np.array([np.mean(np.cos(2*Bk_n[mask]*tm/hbarN)) for tm in tt])
print("  times:", tt)
print("  DENSITY <sx>(t):", np.round(dens,4))
print("  STATE   <sx>(t):", np.round(state,4), " (decays/dephases)")
