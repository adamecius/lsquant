#!/usr/bin/env python3
"""
Symbolic derivation of the clean 1D-chain transport closed forms used by test #4
(VAC <-> MSD equivalence). Documentation/proof only -- nothing in the build/test
pipeline runs this; the C++ comparators embed the resulting closed form. Mirrors the
derive_core.py / derive_spin.py pattern.

Result (gamma=a=hbar=1, chain E(k)=2 gamma cos(ka)):
    v^2(E)      = (a/hbar)^2 (4 gamma^2 - E^2)          group velocity squared
    VAC(E,t)    = v^2(E)                                 constant (ballistic, no scattering)
    MSD(E,t)    = v^2(E) t^2                             = (4 gamma^2 - E^2) t^2
    sigma(E,t)  = e^2 rho(E) v^2(E) t                    -> slope (e^2 a^2/pi hbar^2) sqrt(4 gamma^2 - E^2)
VAC and MSD give the SAME sigma (MSD = double time-integral of VAC), so #4 has an exact
analytic target -- no fit. This matches lsquant_chain_reference.{vac,msd,sigma_t}.
"""
import sympy as sp

k, a, g, hbar, t, E, e = sp.symbols('k a gamma hbar t E e', positive=True)

Ek  = 2*g*sp.cos(k*a)                         # dispersion
v   = sp.diff(Ek, k)/hbar                      # group velocity v = (1/hbar) dE/dk
v2  = sp.simplify(v**2)                         # (2 g a sin(ka)/hbar)^2
v2E = sp.simplify((2*g*a/hbar)**2 * (1-(E/(2*g))**2))   # in terms of E via sin^2 = 1-(E/2g)^2

rho   = 1/(sp.pi*sp.sqrt(4*g**2 - E**2))       # chain DOS
MSD   = v2E*t**2                               # ballistic
slope = sp.simplify(e**2 * rho * v2E)          # d sigma / dt

ok_v   = sp.simplify(v2E - a**2*(4*g**2 - E**2)/hbar**2) == 0
ok_sig = sp.simplify(slope.subs({a:1, hbar:1}) - e**2*sp.sqrt(4*g**2 - E**2)/sp.pi) == 0

print("v^2(E)        =", v2E)
print("MSD(E,t)      =", sp.simplify(MSD))
print("sigma slope   =", sp.simplify(slope.subs({a:1, hbar:1})), "(a=hbar=1)")
print("v^2 closed form matches (a/hbar)^2(4g^2-E^2):", ok_v)
print("sigma slope matches (e^2 a^2/pi hbar^2) sqrt(4g^2-E^2):", ok_sig)
for Ev in (0, 1):
    print(f"  slope(E={Ev}) =", float(slope.subs({a:1, hbar:1, g:1, e:1, E:Ev})))
