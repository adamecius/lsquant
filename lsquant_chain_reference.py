"""
chain_reference.py
Closed-form analytics of the clean 1D monatomic tight-binding chain.
Verified symbolically with sympy in derive_core.py / derive_spin.py.

Conventions
-----------
hopping        gamma > 0           (the thesis' 't')
lattice        a
dispersion     eps(k) = -2 gamma cos(k a),  k in (-pi/a, pi/a]
half-bandwidth DeltaE = 2 gamma,  band centre Ebar = 0
rescaled       eps~(k) = eps/DeltaE = -cos(k a)
velocity       v(k)   = (1/hbar) d eps/dk = (2 gamma a / hbar) sin(k a)

All quantities are PER SITE.  Hopping and time are 'gamma' and 't' (no name clash).
"""
import numpy as np

# ----------------------------------------------------------------------
# Spectral / band quantities  (charge sector, spinless)
# ----------------------------------------------------------------------
def dos(E, gamma=1.0):
    """DOS per site: 1/(pi sqrt(4 gamma^2 - E^2)), |E|<2 gamma, else 0."""
    E = np.asarray(E, float)
    out = np.zeros_like(E)
    inside = np.abs(E) < 2*gamma
    out[inside] = 1.0/(np.pi*np.sqrt(4*gamma**2 - E[inside]**2))
    return out

def v2(E, gamma=1.0, a=1.0, hbar=1.0):
    """Group velocity squared on the energy shell: (a^2/hbar^2)(4 gamma^2 - E^2)."""
    E = np.asarray(E, float)
    return (a**2/hbar**2)*np.clip(4*gamma**2 - E**2, 0.0, None)

# ----------------------------------------------------------------------
# Chebyshev moments  (exact)
# ----------------------------------------------------------------------
def dos_moment(m):
    """Cbar^DOS_m = (1/N) Tr[T_m(H~)] = delta_{m,0}."""
    return 1.0 if m == 0 else 0.0

def kg_moment(m, n, gamma=1.0, a=1.0, hbar=1.0):
    """
    Kubo-Greenwood Chebyshev moment, thesis B.15 normalization
        C^KG_mn = Tr[T_m V T_n V] / ( N (1+d_m0)(1+d_n0) ).
    Closed form:
        (gamma^2 a^2/hbar^2)[ d_mn/(1+d_m0+d_m1) - 1/2 (1-d_mn) d_{|m-n|,2} ].
    """
    d = lambda i, j: 1.0 if i == j else 0.0
    pref = gamma**2 * a**2 / hbar**2
    val = d(m, n)/(1 + d(m, 0) + d(m, 1)) - 0.5*(1 - d(m, n))*d(abs(m-n), 2)
    return pref*val

def kg_moment_matrix(M, **kw):
    return np.array([[kg_moment(m, n, **kw) for n in range(M+1)] for m in range(M+1)])

# ----------------------------------------------------------------------
# Ballistic time-dependent transport  (clean chain: [H,V]=0)
# ----------------------------------------------------------------------
def vac(E, t, gamma=1.0, a=1.0, hbar=1.0):
    """Velocity autocorrelation C_vv(E,t) = v^2(E), constant in time (ballistic)."""
    return np.broadcast_to(v2(E, gamma, a, hbar), np.shape(np.broadcast_arrays(E, t)[0])).copy()

def msd(E, t, gamma=1.0, a=1.0, hbar=1.0):
    """Mean-square displacement DX^2(E,t) = v^2(E) t^2 (ballistic)."""
    return v2(E, gamma, a, hbar)*np.asarray(t, float)**2

def sigma_t(E, t, gamma=1.0, a=1.0, hbar=1.0, e=1.0):
    """
    Time-dependent conductivity (ballistic), equal from VAC and MSD:
        sigma(E,t) = e^2 rho(E) v^2(E) t = (e^2 a^2/pi hbar^2) sqrt(4 gamma^2-E^2) * t.
    Diverges linearly in t: no DC limit for the clean chain.
    """
    return (e**2 * a**2/(np.pi*hbar**2))*np.sqrt(np.clip(4*gamma**2 - E**2, 0, None))*np.asarray(t, float)

# ----------------------------------------------------------------------
# Spin dynamics  (magnetic chain  H = H0 - J_ex sigma_z, initial spin +x)
# Both flavors (state-projection / density-projection) give this, E_F-independent.
# ----------------------------------------------------------------------
def spin_precession(t, Jex=0.1, hbar=1.0):
    """Returns (<Sx>, <Sy>, <Sz>)/(hbar/2) = (cos w t, -sin w t, 0), w=2 J_ex/hbar."""
    t = np.asarray(t, float)
    w = 2*Jex/hbar
    return np.cos(w*t), -np.sin(w*t), np.zeros_like(t)

def precession_frequency(Jex=0.1, hbar=1.0):
    """|f| = J_ex/(pi hbar)  (thesis Fig 4.3)."""
    return Jex/(np.pi*hbar)

# ----------------------------------------------------------------------
# Self-test
# ----------------------------------------------------------------------
if __name__ == "__main__":
    print("DOS moments m=0..4:", [dos_moment(m) for m in range(5)])
    print("KG matrix 5x5 (units gamma^2 a^2/hbar^2):")
    print(kg_moment_matrix(4))
    print("rho(0)=", dos(0.0), " expected 1/(2 pi)=", 1/(2*np.pi))
    print("v2(0)=", v2(0.0), " expected 4 (gamma=a=hbar=1)")
    print("MSD/t^2 at E=0:", msd(0.0, 3.0)/9, " expected v2(0)=4")
    Sx, Sy, Sz = spin_precession(np.array([0.0, np.pi/(2*0.1*2)]), Jex=0.1)
    print("precession Sx at t=0 and quarter period:", Sx)
    print("|f| =", precession_frequency(0.1), " expected", 0.1/np.pi)
