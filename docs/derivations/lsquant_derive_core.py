import sympy as sp

th, m, n, k, a, hbar, t, Jex, EF, E = sp.symbols(
    'theta m n k a hbar t J_ex E_F E', real=True)
g = sp.symbols('gamma', positive=True)  # hopping magnitude > 0 (lets sympy prove |g|/g = 1)

print("="*70)
print("A) DOS Chebyshev moments  Cbar^DOS_m = (1/2pi) int T_m(-cos th) dth")
print("="*70)
# tilde eps = -cos(theta), theta = k a in [-pi,pi]; Delta E = 2 gamma, Ebar=0
for mm in range(0, 9):
    Tm = sp.chebyshevt(mm, -sp.cos(th))
    val = sp.integrate(Tm, (th, -sp.pi, sp.pi)) / (2*sp.pi)
    print(f"  m={mm}:  Cbar^DOS_m = {sp.simplify(val)}")

print()
print("="*70)
print("B) DOS reconstruction check (only m=0 survives) -> closed-form DOS")
print("="*70)
Et = sp.symbols('Etil', real=True)  # tilde E = E/DeltaE
# rho(Et) = (1/(DeltaE*Omega)) sum_m Cdos_m g_m delta_m ;  only m=0, g0=1, Cdos0=1
# delta_0(Et) = (2-1) T_0 / (pi sqrt(1-Et^2)) = 1/(pi sqrt(1-Et^2))
DeltaE = 2*g
rho_tilde = sp.Rational(1,1) * (1/sp.pi/sp.sqrt(1-Et**2))
rho_of_E = (1/DeltaE) * rho_tilde.subs(Et, E/DeltaE)
print("  rho(E) =", sp.simplify(rho_of_E), "  (per site)")
print("  expected 1/(pi sqrt(4 gamma^2 - E^2)) ->",
      sp.simplify(rho_of_E - 1/(sp.pi*sp.sqrt(4*g**2 - E**2))) == 0)

print()
print("="*70)
print("C) Velocity matrix element & KG moment matrix C^KG_mn")
print("="*70)
# v(theta) = (2 gamma a/hbar) sin(theta);  V^2 = (2 g a/hbar)^2 sin^2
# Cbar^KG_mn = (1/2pi) int V^2 T_m(-cos) T_n(-cos) dth
V2 = (2*g*a/hbar)**2 * sp.sin(th)**2
def CKG(mm, nn):
    integrand = V2 * sp.chebyshevt(mm, -sp.cos(th)) * sp.chebyshevt(nn, -sp.cos(th))
    return sp.simplify(sp.integrate(integrand, (th, -sp.pi, sp.pi)) / (2*sp.pi))

N = 6
print("  C^KG_mn / (gamma^2 a^2 / hbar^2):")
pref = g**2 * a**2 / hbar**2
for mm in range(N+1):
    row = []
    for nn in range(N+1):
        row.append(sp.nsimplify(sp.simplify(CKG(mm, nn)/pref)))
    print(f"   m={mm}: ", row)

print()
print("  Compare to thesis closed form B.15:")
print("   C^KG_mn = (g^2 a^2/hbar^2)[ delta_mn/(1+delta_m0+delta_m1)"
      " - 1/2 (1-delta_mn) delta_{|m-n|,2} ]")
def closed(mm, nn):
    d_mn = 1 if mm==nn else 0
    d_m0 = 1 if mm==0 else 0
    d_m1 = 1 if mm==1 else 0
    d_abs2 = 1 if abs(mm-nn)==2 else 0
    return sp.Rational(d_mn, 1+d_m0+d_m1) - sp.Rational(1,2)*(1-d_mn)*d_abs2
# The raw integral Tr[T_m V T_n V]/N gives the UN-normalized moments (C00=2, C02=-1);
# the B.15 closed form is the NORMALIZED convention (C00=1/2, C02=-1/2). They reconcile
# via the per-index Chebyshev factor (1+delta_m0)(1+delta_n0): raw / factor = normalized
# (4 at (0,0), 2 at (0,2), 1 elsewhere). Divide the raw by it before comparing.
def norm_factor(mm, nn):
    return (1 + (1 if mm==0 else 0)) * (1 + (1 if nn==0 else 0))
ok = True
for mm in range(N+1):
    for nn in range(N+1):
        normalized = sp.nsimplify(CKG(mm,nn)/pref / norm_factor(mm,nn))
        diff = sp.simplify(normalized - closed(mm,nn))
        if diff != 0:
            ok = False
            print(f"   MISMATCH m={mm} n={nn}: raw/factor={normalized} closed={closed(mm,nn)}")
print("  raw/((1+d_m0)(1+d_n0)) matches B.15 normalized closed form:", ok)
