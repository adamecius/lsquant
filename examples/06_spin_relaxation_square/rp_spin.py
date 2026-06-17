"""Random-phase (stochastic-trace) spin-relaxation evaluation for tutorial 06,
and its validation against the exact trace and against Elliott-Yafet theory.

Method. S_a(E,t) = Tr[S_a delta(E-H) U(t) P U(t)^dag] / Tr[delta(E-H) U(t)P U(t)^dag].
The exact trace loops over all N basis states (O(N) recursions -> slow). The
stochastic trace replaces Tr[A] by (N/Nr) sum_r <r|A|r> over Nr random-phase
vectors |r>; the error scales as ~1/sqrt(Nr*N), so a LARGER system needs FEWER
vectors. Here each random vector is one driver run with no state file
(KPM_SEED=r selects the vector); we accumulate the numerator and denominator
moments across r and take the ratio (the S(0)=1 normalisation cancels the trace
prefactor). No lsquant source change -- pure orchestration.
"""
import os, glob, subprocess
import numpy as np
from shutil import which

E = 2.0
DEVNULL = subprocess.DEVNULL

def _bin(name):
    b = os.environ.get("LSQUANT_BIN")
    return os.path.join(b, name) if b and os.path.exists(os.path.join(b, name)) else (which(name) or name)

def _te(env, label, op, proj, M, ntime, tmax):
    for f in glob.glob("TimeProj%s-%s%s*chebmomTD" % (op, proj, label)): os.remove(f)
    subprocess.run([_bin("inline_compute-kpm-TimeEvProjetedOp"), label, op, proj,
                    str(M), str(ntime), str(tmax)], env=env, check=True, stdout=DEVNULL, stderr=DEVNULL)
    return glob.glob("TimeProj%s-%s%s*chebmomTD" % (op, proj, label))[0]

def _recon(mom, op, proj, label, eta_mev):
    for f in glob.glob("mean*%s-%s%s*EF*dat" % (op, proj, label)): os.remove(f)
    subprocess.run([_bin("inline_timeCorrelationsFromChebmom"), mom, str(eta_mev), str(E)],
                   check=True, stdout=DEVNULL, stderr=DEVNULL)
    dat = glob.glob("mean*%s-%s%s*EF%f*dat" % (op, proj, label, E))[0]
    t, v = np.loadtxt(dat, unpack=True)
    return t, v

def exact_spin_vs_time(label, N, op, proj, M, ntime, tmax, eta_mev=140):
    """Reference: exact trace over all N basis states (one driver call, 'local' file)."""
    state = "exact_%d" % N
    with open(state, "w") as f:
        f.write("local\n%d\n%s\n" % (N, "\n".join(map(str, range(N)))))
    env = dict(os.environ)
    def one(o, p):
        for f in glob.glob("TimeProj%s-%s%s*chebmomTD" % (o, p, label)): os.remove(f)
        subprocess.run([_bin("inline_compute-kpm-TimeEvProjetedOp"), label, o, p,
                        str(M), str(ntime), str(tmax), state], env=env, check=True, stdout=DEVNULL, stderr=DEVNULL)
        mom = glob.glob("TimeProj%s-%s%s*chebmomTD" % (o, p, label))[0]
        return _recon(mom, o, p, label, eta_mev)
    t, num = one(op, proj)
    _, den = one("1", proj)
    good = t > 0
    s = num[good] / den[good].mean()
    return t[good], s / s[0]

def rp_spin_vs_time(label, N, op, proj, M, ntime, tmax, Nr=20, eta_mev=140, threads=16):
    """Stochastic trace over Nr random-phase vectors (KPM_SEED = 1..Nr)."""
    NUM = DEN = t = None
    for r in range(1, Nr + 1):
        env = dict(os.environ, KPM_SEED=str(r), OMP_NUM_THREADS=str(threads),
                   OMP_PROC_BIND="close", OMP_PLACES="cores")
        mom = _te(env, label, op, proj, M, ntime, tmax)
        tt, num = _recon(mom, op, proj, label, eta_mev); os.remove(mom)
        mom = _te(env, label, "1", proj, M, ntime, tmax)
        _, den = _recon(mom, "1", proj, label, eta_mev); os.remove(mom)
        NUM = num if NUM is None else NUM + num
        DEN = den if DEN is None else DEN + den
        t = tt
    good = t > 0
    s = NUM[good] / DEN[good].mean()
    return t[good], s / s[0]

def fit_rate(t, s, lo=0.2, hi=0.85):
    m = (s > lo) & (s < hi)
    if m.sum() < 3:
        return np.nan
    slope = np.polyfit(t[m], np.log(s[m]), 1)[0]
    return -1.0 / slope
