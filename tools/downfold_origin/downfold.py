#!/usr/bin/env python
"""
downfold.py -- deterministic, automatic reduction of a SIESTA Hamiltonian (non-orthogonal
multi-zeta LCAO basis, with its spurious high-energy states) to a small orthogonal
tight-binding model whose cost for the kernel polynomial method (KPM) is minimal.

Procedure (closed form, no minimisation):
  1. Read H(R), S(R) from <label>.HSX and the basis (n, l, zeta) from *.ion.xml; valence,
     semicore and polarisation shells from the PSML valence configuration.
  2. 2D systems with a horizontal mirror: split into mirror-even / mirror-odd blocks.
  3. Solve H(k)c = E S(k)c on a Gamma-centred grid; find groups of bands separated by gaps.
  4. Isolated candidates: unions of consecutive groups per block, at most one per block,
     whose union contains every state of the Fermi window [E_VBM - w, E_CBM + w].  For each,
     every combination of atomic shells with the right count is tried; the one with the
     largest min_k singular value of A(k) = <psi|S|g> is kept, then the zetas of each shell
     are contracted to maximise the projection.  Gauge: U = u w^H from the SVD of A
     (Loewdin / polar decomposition, as in Souza et al. 2001 Eq. 23, FPLO, Ozaki 2024).
  5. Frozen-window candidates for every valence(+polarisation) shell set: window states are
     kept exactly and the remaining J - n_F states are those of Souza et al. 2001, Eq. 27
     (outer states spanning the trial directions not covered by the window), here with an
     extra energy weight exp(-(E - W_top)/tau) and used as the final model (no iteration).
  6. H_W(R) by FFT with Wigner-Seitz minimal images; validation on random off-grid k-points
     plus the high-symmetry path.
  7. Choice among all candidates: fewest non-zero elements per cell (J x hoppings) needed to
     reach --tol, i.e. the cost of a KPM matrix-vector product.
  8. Write hoppings, summary JSON and a band comparison plot.
"""
import argparse, itertools, json, os, sys, warnings
import numpy as np, scipy.linalg as sl
from scipy.special import erfc
import sisl
warnings.filterwarnings("ignore")

LSYM = "spdf"


# ----------------------------------------------------------------------------- input
def load(rundir, label):
    sile = sisl.get_sile(os.path.join(rundir, f"{label}.HSX"))
    H = sile.read_hamiltonian()
    geom = H.geometry
    spinor = H.spin.is_noncolinear or H.spin.is_spinorbit
    orbinfo = []          # per basis function (spinor components interleaved, as in sisl)
    for ia in range(geom.na):
        atom = geom.atoms[ia]
        basis = sisl.get_sile(os.path.join(rundir, f"{atom.tag}.ion.xml")).read_basis()
        val = psml_valence(os.path.join(rundir, f"{atom.tag}.psml"))
        if basis.no != atom.no:
            raise RuntimeError(f"basis mismatch for {atom.tag}: {basis.no} vs {atom.no}")
        nmax = max(n for n, l in val) if val else None
        for j, (o, ob) in enumerate(zip(atom.orbitals, basis.orbitals)):
            if (o.l, o.m, o.zeta) != (ob.l, ob.m, ob.zeta):
                raise RuntimeError(f"orbital order mismatch {atom.tag} #{j}")
            if val is None:
                pol, semi = bool(ob.P), False
            else:
                pol = (ob.n, ob.l) not in val
                # semicore: a filled valence-configuration shell with another shell of larger n
                semi = (not pol) and val[(ob.n, ob.l)] >= 2 * (2 * ob.l + 1) - 1e-6 and ob.n < nmax
            for sp in ((0, 1) if spinor else (None,)):
                orbinfo.append(dict(atom=ia, species=atom.tag, n=ob.n, l=ob.l, m=ob.m,
                                    zeta=ob.zeta, pol=pol, semi=semi, spin=sp))
    nspin = 2 if H.spin.is_polarized else 1
    return H, orbinfo, nspin


def psml_valence(fn):
    """{(n, l): occupation} of the valence configuration of a PSML pseudopotential."""
    import re
    if not os.path.exists(fn):
        return None
    txt = open(fn).read()
    m = re.search(r"<valence-configuration.*?</valence-configuration>", txt, re.S)
    sh = re.findall(r'<shell\s+n="(\d)"\s+l="([spdf])"\s+occupation="([\d.]+)"', m.group(0) if m else txt)
    return {(int(n), "spdf".index(l)): float(o) for n, l, o in sh}


def is_2d(geom):
    return geom.nsc[2] == 1


# ----------------------------------------------------------------------------- symmetry
def mirror_blocks(geom, orbinfo, tol=1e-3):
    """Return list of (name, B) with B (norb x nb) orthonormal basis of each mirror block,
    or [("full", I)] if there is no horizontal mirror."""
    no = len(orbinfo)
    full = [("full", np.eye(no))]
    if not is_2d(geom) or any(o["spin"] is not None for o in orbinfo):
        return full
    z = geom.xyz[:, 2]
    z0 = 0.5 * (z.max() + z.min())
    amap = {}
    fxy = geom.fxyz[:, :2] % 1.0
    for ia in range(geom.na):
        found = None
        for ja in range(geom.na):
            if geom.atoms[ja].tag != geom.atoms[ia].tag:
                continue
            dxy = (fxy[ja] - fxy[ia] + 0.5) % 1.0 - 0.5
            if np.abs(dxy).max() < tol and abs(z[ja] - (2 * z0 - z[ia])) < tol:
                found = ja
        if found is None:
            return full
        amap[ia] = found
    first = np.cumsum([0] + [geom.atoms[i].no for i in range(geom.na)])
    M = np.zeros((no, no))
    for ia, ja in amap.items():
        for j in range(geom.atoms[ia].no):
            o = orbinfo[first[ia] + j]
            M[first[ja] + j, first[ia] + j] = (-1) ** (o["l"] + o["m"])
    w, v = np.linalg.eigh(M)
    return [("mirror-even", v[:, w > 0]), ("mirror-odd", v[:, w < 0])]


# ----------------------------------------------------------------------------- k grids
def kgrid(geom, n2d, n3d):
    N = (n2d, n2d, 1) if is_2d(geom) else (n3d, n3d, n3d)
    ks = np.array([[i / N[0], j / N[1], l / N[2]] for i in range(N[0])
                   for j in range(N[1]) for l in range(N[2])])
    return N, ks


def validation_kpoints(geom, N, nval):
    """Random off-grid points plus the high-symmetry path (it crosses Dirac points / Fermi surfaces)."""
    rng = np.random.default_rng(1)
    kr = rng.random((nval, 3)); kr[:, [d for d in range(3) if N[d] == 1]] = 0
    kp = kpath(geom, 80).k
    return np.concatenate([kr, kp])


def kpath(geom, npts=120):
    """High-symmetry path in reduced coordinates of geom.rcell (simple lattice detection)."""
    A = geom.cell
    def frac(kcart_2pi_over_a, a):          # cartesian (2pi/a units) -> reduced
        return np.array(kcart_2pi_over_a) @ A.T / a
    L = np.linalg.norm(A, axis=1)
    cos = lambda i, j: A[i] @ A[j] / (L[i] * L[j])
    if is_2d(geom) and abs(L[0] - L[1]) < 1e-3 and abs(abs(cos(0, 1)) - 0.5) < 1e-3:
        b = geom.rcell
        K = None
        for cand in ([1/3, 1/3, 0], [2/3, 1/3, 0], [1/3, 2/3, 0]):   # pick the true K
            kc = np.array(cand) @ b
            if abs(np.linalg.norm(kc) - 4 * np.pi / (3 * L[0])) < 1e-6:
                K = cand; break
        pts, lab = [[0, 0, 0], [0.5, 0, 0], K, [0, 0, 0]], ["Γ", "M", "K", "Γ"]
    elif not is_2d(geom) and np.allclose(L, L[0], atol=1e-3):
        c = cos(0, 1)
        if abs(c - 0.5) < 1e-3:                       # fcc primitive
            a = L[0] * np.sqrt(2)
            P = dict(Γ=[0, 0, 0], X=[0, 1, 0], W=[.5, 1, 0], L=[.5, .5, .5], K=[.75, .75, 0])
            lab = ["Γ", "X", "W", "L", "Γ", "K"]
        elif abs(abs(c) - 1/3) < 1e-3:                # bcc primitive
            a = L[0] * 2 / np.sqrt(3)
            P = dict(Γ=[0, 0, 0], H=[0, 0, 1], N=[.5, .5, 0], P=[.5, .5, .5])
            lab = ["Γ", "H", "N", "Γ", "P", "H"]
        else:
            P = None
        if P is not None:
            pts = [frac(P[l], a) for l in lab]
    else:
        P = None
    if "pts" not in locals():
        pts = [[0, 0, 0], [.5, 0, 0], [.5, .5, 0], [0, 0, 0]]
        lab = ["Γ", "X1", "X12", "Γ"]
    bs = sisl.BandStructure(geom, pts, npts, lab)
    return bs


# ----------------------------------------------------------------------------- spectra
class Spectrum:
    """Eigen-decomposition of every block on a set of k-points for one spin."""
    def __init__(self, H, spin, blocks, ks):
        self.ks = ks
        self.blocks = blocks
        self.E = {b: [] for b, _ in blocks}        # (nk, nb)
        self.V = {b: [] for b, _ in blocks}        # (nk, norb, nb) full-space S-normalised
        self.S = []
        for k in ks:
            Hk, Sk = HSk(H, spin, k)
            self.S.append(Sk)
            for b, B in blocks:
                e, c = sl.eigh(B.T @ Hk @ B, B.T @ Sk @ B)
                self.E[b].append(e)
                self.V[b].append(B @ c)
        for b, _ in blocks:
            self.E[b] = np.array(self.E[b])
            self.V[b] = np.array(self.V[b])
        self.S = np.array(self.S)

    def all_E(self):
        return np.sort(np.concatenate([self.E[b] for b, _ in self.blocks], axis=1), axis=1)


def band_groups(E, gap_tol):
    """Groups of consecutive bands separated by a global gap > gap_tol."""
    groups, start = [], 0
    nb = E.shape[1]
    for n in range(nb - 1):
        if E[:, n + 1].min() - E[:, n].max() > gap_tol:
            groups.append((start, n)); start = n + 1
    groups.append((start, nb - 1))
    return [dict(lo=a, hi=b, J=b - a + 1, Emin=E[:, a].min(), Emax=E[:, b].max()) for a, b in groups]


def fermi_window(Eall, w):
    occ = Eall[Eall <= 0]; emp = Eall[Eall > 0]
    vbm, cbm = occ.max(), emp.min()
    metal = cbm - vbm < 0.05 or np.any((Eall.min(0) < 0) & (Eall.max(0) > 0))
    return vbm, cbm, metal, (vbm - w, cbm + w)


# ----------------------------------------------------------------------------- shells / trials
def shells(orbinfo):
    """Atomic shells: (species, n, l, pol) -> list of zeta-1 orbital indices (all atoms)."""
    sh = {}
    for i, o in enumerate(orbinfo):
        if o["zeta"] != 1:
            continue
        key = (o["species"], o["n"], o["l"], o["pol"])
        sh.setdefault(key, []).append(i)
    return sh



def zeta_rows(orbinfo):
    """(species,n,l,pol) -> {(atom,m): [orbital indices ordered by zeta]}"""
    out = {}
    for i, o in enumerate(orbinfo):
        key = (o["species"], o["n"], o["l"], o["pol"])
        out.setdefault(key, {}).setdefault((o["atom"], o["m"], o["spin"]), []).append((o["zeta"], i))
    return {k: {am: [i for _, i in sorted(v)] for am, v in d.items()} for k, d in out.items()}


def contracted_trials(combo, ZR, V, S, S0, f=None):
    """Trial orbitals: for every shell the combination of its zetas that maximises the
    projection onto the target subspace sum_k S V f V^H S (same coefficients for all m and
    all atoms of the species, so the symmetry is kept). Returns norb x J (S0-normalised)."""
    norb = S.shape[1]
    W = np.einsum("kop,kpb->kob", S, V)                      # S(k) V(k)
    cols = []
    for key in combo:
        groups = ZR[key]
        nz = len(next(iter(groups.values())))
        if nz == 1:
            c = np.array([1.0])
        else:
            MP = np.zeros((nz, nz), complex); MS = np.zeros((nz, nz))
            for rows in groups.values():
                Wr = W[:, rows, :]
                Wf = Wr * (f[:, None, :] if f is not None else 1.0)
                MP += np.einsum("kib,kjb->ij", Wf, Wr.conj())
                MS += S0[np.ix_(rows, rows)] * len(W)
            w, v = sl.eigh(MP.real, MS)
            c = v[:, -1]; c = c * np.sign(c[0])
        for rows in groups.values():
            g = np.zeros(norb); g[rows] = c
            g /= np.sqrt(g[rows] @ S0[np.ix_(rows, rows)] @ g[rows])
            cols.append(g)
    return np.array(cols).T


def shell_name(key):
    sp, n, l, pol = key
    return f"{sp}-{n}{LSYM[l]}" + ("(pol)" if pol else "")


def trials_in_block(idx, B, norb):
    """Combinations of the selected atomic orbitals (indices, or an norb x n matrix of
    orbital vectors) that live in block B."""
    if isinstance(idx, np.ndarray) and idx.ndim == 2:
        G = idx
    else:
        G = np.zeros((norb, len(idx))); G[idx, range(len(idx))] = 1.0
    if B.shape[1] == norb and np.allclose(B, np.eye(norb)):
        return G
    C = B.T @ G                                   # block coordinates
    u, s, _ = np.linalg.svd(C, full_matrices=False)
    r = int((s > 1e-6).sum())
    return B @ u[:, :r]                           # norb x r, each a localised combination


def projection_matrix(V, S, G):
    """A(k) = V^H S G for every k: (nk, nb, nt)."""
    return np.einsum("kob,kop,pt->kbt", V.conj(), S, G)


# ----------------------------------------------------------------------------- models
def isolated_U(A):
    u, s, wh = np.linalg.svd(A)
    return u @ wh, s.min(axis=1)


class Model:
    def __init__(self, geom, centers, N, HWk, ks):
        """HWk: (nk, J, J) on the regular Gamma-centred grid ks of size N."""
        self.geom = geom; self.centers = centers
        J = HWk.shape[1]; self.J = J
        N = np.array(N)
        HR = np.fft.fftn(HWk.reshape(*N, J, J), axes=(0, 1, 2)) / N.prod()   # 1/Nk sum_k e^{-2pi i k.R} H(k)
        A = geom.cell
        m = np.array(list(itertools.product(*[range(n) for n in N])))          # (nm,3)
        m = np.where(m >= N / 2, m - N, m)
        Ls = np.array(list(itertools.product(*[(-1, 0, 1) if n > 1 else (0,) for n in N])))
        Rimg = m[:, None, :] + Ls[None, :, :] * N                               # (nm,nL,3)
        Rcart = Rimg @ A                                                         # (nm,nL,3)
        dvec = centers[None, None, None, :, :] + Rcart[:, :, None, None, :] - centers[None, None, :, None, :]
        d = np.linalg.norm(dvec, axis=-1)                                        # (nm,nL,J,J)
        dmin = d.min(axis=1)
        isbest = d < dmin[:, None] + 1e-4
        deg = isbest.sum(axis=1)                                                 # (nm,J,J)
        t = HR[m[:, 0] % N[0], m[:, 1] % N[1], m[:, 2] % N[2]]                   # (nm,J,J)
        im, il, i, j = np.nonzero(isbest)
        self.R = Rimg[im, il]
        self.i, self.j = i, j
        self.t = t[im, i, j] / deg[im, i, j]
        self.d = dmin[im, i, j]

    def Hk(self, k, rcut=np.inf, tmin=0.0):
        sel = (self.d <= rcut + 1e-6) & (np.abs(self.t) >= tmin)
        ph = np.exp(2j * np.pi * (self.R[sel] @ np.asarray(k)))
        Hk = np.zeros((self.J, self.J), complex)
        np.add.at(Hk, (self.i[sel], self.j[sel]), self.t[sel] * ph)
        return 0.5 * (Hk + Hk.conj().T)

    @property
    def hops(self):
        return [(tuple(R), i, j, t, d) for R, i, j, t, d in zip(self.R, self.i, self.j, self.t, self.d)]

    def shells(self):
        return sorted(set(np.round(self.d, 3)))


# ----------------------------------------------------------------------------- driver
def analyse_spin(H, orbinfo, spin, args, out, tag):
    geom = H.geometry
    norb = len(orbinfo)
    blocks = mirror_blocks(geom, orbinfo)
    N, ks = kgrid(geom, args.n2d, args.n3d)
    kval = validation_kpoints(geom, N, args.nval)
    log = []
    say = lambda s: (print(s), log.append(s))
    say(f"\n=== {tag}: {norb} orbitales, bloques: {[b for b, _ in blocks]}, malla {N}")

    spec = Spectrum(H, spin, blocks, ks)
    full_spec = spec if blocks[0][0] == "full" else Spectrum(H, spin, [("full", np.eye(norb))], ks)
    Eall = full_spec.E["full"]
    vbm, cbm, metal, Wf = fermi_window(Eall, args.w)
    say(f"  E_F=0  VBM={vbm:+.3f}  CBM={cbm:+.3f}  {'metal' if metal else 'gap=%.3f eV' % (cbm - vbm)}"
        f"  ventana Fermi=[{Wf[0]:+.2f},{Wf[1]:+.2f}]  E_max espectro={Eall.max():.1f} eV")
    n_in_W = ((Eall >= Wf[0]) & (Eall <= Wf[1])).sum(1)            # per k
    SH = shells(orbinfo)
    ZR = zeta_rows(orbinfo)
    S0 = H.tocsr(H.S_idx)[:, :H.no].toarray()                    # S(R=0)
    if S0.shape[0] != norb:
        S0 = np.kron(S0, np.eye(2))
    # semicore: a shell with the same (species, l) as another non-polarisation shell of larger n
    semicore = {(o["species"], o["n"], o["l"], o["pol"]) for o in orbinfo if o["semi"]}
    if semicore:
        say("  semicore: " + ", ".join(shell_name(k) for k in sorted(semicore)))

    # reference eigenvalues on the validation k-points (cached)
    Eval_full = [ref_block(H, spin, np.eye(norb), k) for k in kval]
    Eval_blk = {b: [ref_block(H, spin, B, k) for k in kval] for b, B in blocks}

    # ---------------- isolated candidates, per block
    cands = {b: [] for b, _ in blocks}
    for b, B in blocks:
        E = spec.E[b]
        grps = [g for g in band_groups(E, args.gap) if g["Emin"] > -args.emax_below and g["Emax"] < args.emax_above]
        say(f"  [{b}] grupos aislados: " + ", ".join(f"{g['J']}b[{g['Emin']:+.1f},{g['Emax']:+.1f}]" for g in grps))
        for a in range(len(grps)):
            for c in range(a, len(grps)):
                if any(grps[x + 1]["lo"] != grps[x]["hi"] + 1 for x in range(a, c)):
                    continue
                lo, hi = grps[a]["lo"], grps[c]["hi"]; J = hi - lo + 1
                if J > args.jmax:
                    continue
                Ec = E[:, lo:hi + 1]
                V = spec.V[b][:, :, lo:hi + 1]
                best = None
                keys = list(SH)
                for r in range(1, len(keys) + 1):
                    for combo in itertools.combinations(keys, r):
                        idx = sum((SH[k] for k in combo), [])
                        G = trials_in_block(idx, B, norb)
                        if G.shape[1] != J:
                            continue
                        A = projection_matrix(V, spec.S, G)
                        smin = np.linalg.svd(A, compute_uv=False).min()
                        score = (smin, -sum(k[3] for k in combo))
                        if best is None or score > best[0]:
                            best = (score, combo, G, A)
                if best is None:
                    continue
                # refine: optimal zeta contraction of the chosen shells
                Gc = trials_in_block(contracted_trials(best[1], ZR, V, spec.S, S0), B, norb)
                if Gc.shape[1] == J:
                    Ac = projection_matrix(V, spec.S, Gc)
                    sc = np.linalg.svd(Ac, compute_uv=False).min()
                    if sc > best[0][0]:
                        best = ((sc, best[0][1]), best[1], Gc, Ac)
                cands[b].append(dict(block=b, lo=lo, hi=hi, J=J, smin=best[0][0], npol=-best[0][1],
                                     combo=best[1], G=best[2], A=best[3], E=Ec,
                                     inW=((Ec >= Wf[0]) & (Ec <= Wf[1])).sum(1)))
                c_ = cands[b][-1]
                say(f"    cand. {b:11s} J={J:2d} [{Ec.min():+6.1f},{Ec.max():+6.1f}] eV  σ_min={c_['smin']:.3f}  "
                    f"capas={'+'.join(shell_name(k) for k in c_['combo'])}")
    # combine at most one candidate per block; the union must contain every state in the Fermi window
    combos = []
    for pick in itertools.product(*[[None] + cands[b] for b, _ in blocks]):
        pick = [c for c in pick if c is not None]
        if not pick:
            continue
        if not np.all(sum(c["inW"] for c in pick) == n_in_W):
            continue
        smin = min(c["smin"] for c in pick)
        if smin < args.smin:
            continue
        combos.append((sum(c["J"] for c in pick), -smin, sum(c["npol"] for c in pick), pick))
    options = []          # (choice, HWk, G, errfn)
    for J, msmin, _, pick in sorted(combos, key=lambda x: x[:3])[:3]:
        Js = [c["J"] for c in pick]; off = np.cumsum([0] + Js)
        HWk = np.zeros((len(ks), J, J), complex)
        for c, o in zip(pick, off):
            U, _ = isolated_U(c["A"])
            HWk[:, o:o + c["J"], o:o + c["J"]] = np.einsum("kbi,kb,kbj->kij", U.conj(), c["E"], U)
        G = np.concatenate([c["G"] for c in pick], axis=1)
        ch = dict(kind="aislado", block="+".join(f"{c['block']}({c['J']})" for c in pick), J=J,
                  combo=sorted(set(sum((list(c["combo"]) for c in pick), [])), key=str), smin=-msmin)
        Eref_val = [np.sort(np.concatenate([Eval_blk[c["block"]][ik][c["lo"]:c["hi"] + 1] for c in pick]))
                    for ik in range(len(kval))]
        errfn = (lambda Ev: (lambda m, rc: max(np.abs(np.linalg.eigvalsh(m.Hk(k, rc)) - Er).max()
                                               for k, Er in zip(kval, Ev))))(Eref_val)
        options.append((ch, HWk, G, errfn))

    # ---------------- frozen-window projection (Souza, Marzari, Vanderbilt, PRB 65, 035109, Eq. 27; concept: Thygesen et al. 2005)
    valence = [k for k in SH if not k[3] and k not in semicore]
    pols = [k for k in SH if k[3]]
    sets = [valence] + [valence + list(p) for r in range(1, len(pols) + 1) for p in itertools.combinations(pols, r)]
    V = full_spec.V["full"]; Sf = full_spec.S
    werr = lambda m, rc: window_error(m, kval, Eval_full, Wf, rc)
    p_low = None
    for combo in sets:
        G0 = trials_in_block(sum((SH[k] for k in combo), []), np.eye(norb), norb)
        J = G0.shape[1]
        if J < n_in_W.max():
            continue
        if p_low is None:                                   # semicore cut: groups dominated by semicore shells
            p_low = -np.inf
            if semicore:
                Gs = trials_in_block(sum((SH[k] for k in semicore), []), np.eye(norb), norb)
                As = projection_matrix(V, Sf, Gs); p = np.empty(Eall.shape)
                for ik in range(len(ks)):
                    X = sl.fractional_matrix_power(Gs.T @ Sf[ik] @ Gs, -0.5)
                    p[ik] = np.sum(np.abs(As[ik] @ X) ** 2, axis=1).real
                for g in band_groups(Eall, args.gap):
                    if g["Emax"] < Wf[0] and p[:, g["lo"]:g["hi"] + 1].mean() > 0.5:
                        p_low = max(p_low, g["Emax"])
                say(f"  estados de semicore excluidos: E <= {p_low:+.1f} eV")
        frozen = (Eall >= Wf[0]) & (Eall <= Wf[1])
        outer = (~frozen) & (Eall > p_low) & (Eall < args.outer_max)
        for tau in args.tau_grid:
            # weight of outer states: 1 below the window, exp(-(E-Wtop)/tau) above it
            wout = np.where(Eall > Wf[1], np.exp(-np.clip(Eall - Wf[1], 0, None) / tau), 1.0) * outer
            Gc = contracted_trials(combo, ZR, V, Sf, S0, np.where(frozen, 1.0, wout))
            if np.any(frozen.sum(1) > J) or np.any(frozen.sum(1) + outer.sum(1) < J):
                say(f"    (descartado {'+'.join(shell_name(k) for k in combo)}: J={J} no cabe en los estados razonables)")
                break
            HWk = np.zeros((len(ks), J, J), complex); smin = np.inf
            for ik in range(len(ks)):
                F = np.where(frozen[ik])[0]; O = np.where(outer[ik])[0]
                AF = V[ik][:, F].conj().T @ Sf[ik] @ Gc
                AO = V[ik][:, O].conj().T @ Sf[ik] @ Gc
                nx = J - len(F)
                # trial directions already described by the frozen states are projected out, so the
                # extra states cover what the window misses (Souza et al. 2001, Eq. 27; tau-weighted here)
                if len(F):
                    _, sF, whF = np.linalg.svd(AF, full_matrices=False)
                    Wr = whF[sF > 1e-6 * max(sF.max(), 1)]
                    Pp = np.eye(J) - Wr.conj().T @ Wr
                else:
                    Pp = np.eye(J)
                u, sv, _ = np.linalg.svd(wout[ik, O][:, None] * (AO @ Pp), full_matrices=True)
                C = u[:, :nx]
                Hsub = sl.block_diag(np.diag(Eall[ik, F]), C.conj().T @ np.diag(Eall[ik, O]) @ C)
                Asub = np.vstack([AF, C.conj().T @ AO])
                uu, ss, wh = np.linalg.svd(Asub)
                smin = min(smin, ss.min())
                Uk = uu @ wh
                HWk[ik] = Uk.conj().T @ Hsub @ Uk
            m = Model(geom, centers_of(Gc, geom, orbinfo), N, HWk, ks)
            ch = dict(kind="ventana congelada", block="full", J=J, combo=list(combo), smin=smin, tau=tau)
            options.append((ch, HWk, Gc, werr, m))
    choice, HWk, G, errfn = select_option(options, geom, orbinfo, N, ks, args, say)

    window = Wf if choice["kind"] != "aislado" else (-np.inf, np.inf)
    say(f"  >>> ELEGIDO: {choice['kind']} {choice['block']} J={choice['J']} capas="
        f"{'+'.join(shell_name(k) for k in choice['combo'])}  σ_min={choice['smin']:.3f}")
    model = Model(geom, centers_of(G, geom, orbinfo), N, HWk, ks)

    # ---------------- truncation radius (bisection over neighbour shells) and error/range trade-off
    nh_of = lambda r: ((model.d > 1e-3) & (model.d <= r + 1e-6)).sum() / model.J
    table = []
    for mark in sorted({0.2, 0.1, 0.05, 0.02, args.tol}, reverse=True):
        r, e = rcut_for(model, errfn, mark)
        if r is not None:
            table.append((mark, r, nh_of(r), e))
    say("  compromiso error / alcance:  " + "   ".join(f"<{m*1e3:.0f} meV: r={d:.2f} Å, {nh:.0f} hop/orb"
                                                     for m, d, nh, e in table))
    rcut, err_final = rcut_for(model, errfn, args.tol)
    if rcut is None:
        rcut = model.shells()[-1]
    nhop = nh_of(rcut)
    say(f"  >>> r_cut={rcut:.2f} Å  ({nhop:.0f} hoppings/orbital)  error_max={err_final*1e3:.1f} meV")

    # ---------------- outputs
    bs = kpath(geom)
    Em = np.array([np.linalg.eigvalsh(model.Hk(k, rcut)) for k in bs.k])
    Er = np.array([ref_block(H, spin, np.eye(norb), k) for k in bs.k])
    width = Em.max() - Em.min()
    say(f"  ancho del espectro: SIESTA {Er.max()-Er.min():.1f} eV  ->  modelo {width:.2f} eV")
    write_outputs(out, tag, model, rcut, choice, bs, Er, Em, window, log, geom, orbinfo)
    return dict(tag=tag, kind=choice["kind"], block=choice["block"], J=choice["J"],
                shells=[shell_name(k) for k in choice["combo"]], smin=float(choice["smin"]),
                rcut=float(rcut), err_meV=float(err_final * 1e3), width_model=float(width),
                width_siesta=float(Er.max() - Er.min()), metal=bool(metal), gap=float(max(cbm - vbm, 0)))


def HSk(H, spin, k):
    Hk = H.Hk(k, spin=spin, format="array") if H.spin.is_polarized else H.Hk(k, format="array")
    Sk = H.Sk(k, format="array")
    if Sk.shape[0] != Hk.shape[0]:                      # spinor H with spin-less S
        Sk = np.kron(Sk, np.eye(2))
    return Hk, Sk


def rcut_for(model, errfn, tol, rmax=None):
    """Smallest neighbour-shell radius with error <= tol (bisection; error ~ monotonic)."""
    sh = [d for d in model.shells() if d > 1e-3 and (rmax is None or d <= rmax)]
    if errfn(model, sh[-1]) > tol:
        return None, errfn(model, sh[-1])
    lo, hi = 0, len(sh) - 1
    while lo < hi:
        mid = (lo + hi) // 2
        if errfn(model, sh[mid]) <= tol:
            hi = mid
        else:
            lo = mid + 1
    return sh[lo], errfn(model, sh[lo])


def select_option(options, geom, orbinfo, N, ks, args, say):
    """Pick the option that reaches --tol with the fewest non-zero elements per cell
    (J x hoppings per orbital): that is the cost of a KPM matrix-vector product.
    Options that cannot reach --tol are ranked after, by their untruncated error."""
    scored = []
    for opt in options:
        ch, HWk, G, ef = opt[:4]
        m = opt[4] if len(opt) > 4 else Model(geom, centers_of(G, geom, orbinfo), N, HWk, ks)
        e_inf = ef(m, np.inf)
        rc, e = rcut_for(m, ef, args.tol) if e_inf <= args.tol else (None, e_inf)
        nnz = ((m.d <= rc + 1e-6).sum()) if rc is not None else np.inf
        scored.append((nnz, e_inf, ch["J"], ch, HWk, G, ef))
        say(f"    opción {ch['kind']:17s} {ch['block']:28s} J={ch['J']:2d} "
            f"capas={'+'.join(shell_name(k) for k in ch['combo']):30s}"
            + (f" τ={ch['tau']:4.1f}" if "tau" in ch else "       ")
            + f"  σ_min={ch['smin']:.3f}  err(sin truncar)={e_inf*1e3:6.1f} meV  "
            + (f"r({args.tol*1e3:.0f} meV)={rc:5.2f} Å  nnz/celda={nnz:6d}" if rc is not None else "no llega a la tolerancia"))
    nnz, e_inf, J, choice, HWk, G, errfn = sorted(scored, key=lambda x: (x[0], x[1], x[2]))[0]
    return choice, HWk, G, errfn


def ref_block(H, spin, B, k):
    Hk, Sk = HSk(H, spin, k)
    return sl.eigh(B.T @ Hk @ B, B.T @ Sk @ B, eigvals_only=True)


def window_error(model, kval, Eref, W, rcut):
    e = 0.0
    for k, Er_all in zip(kval, Eref):
        Er = Er_all[(Er_all >= W[0]) & (Er_all <= W[1])]
        Em = np.linalg.eigvalsh(model.Hk(k, rcut))
        if len(Er):
            e = max(e, np.abs(Er[:, None] - Em[None, :]).min(1).max())
        Ein = Em[(Em >= W[0] + 0.05) & (Em <= W[1] - 0.05)]        # spurious model states in the window
        if len(Ein):
            e = max(e, np.abs(Ein[:, None] - Er_all[None, :]).min(1).max())
    return e


def centers_of(G, geom, orbinfo):
    at = np.array([o["atom"] for o in orbinfo])
    C = []
    for t in G.T:
        w = np.zeros(geom.na)
        np.add.at(w, at, np.abs(t) ** 2)
        C.append((w[:, None] * geom.xyz).sum(0) / w.sum())
    return np.array(C)


def write_outputs(out, tag, model, rcut, choice, bs, Er, Em, window, log, geom, orbinfo):
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    os.makedirs(out, exist_ok=True)
    with open(os.path.join(out, f"{tag}_hoppings.dat"), "w") as fh:
        fh.write(f"# {tag}: orthogonal TB model, J={model.J}, r_cut={rcut:.3f} Ang, energies in eV (E_F = 0)\n")
        fh.write("# cell (Ang):\n")
        for v in geom.cell:
            fh.write("#   %12.6f %12.6f %12.6f\n" % tuple(v))
        fh.write("# orbital centres (Ang):\n")
        for i, c in enumerate(model.centers):
            fh.write("#   %3d %12.6f %12.6f %12.6f\n" % (i, *c))
        fh.write("# R1 R2 R3  i  j   Re(t)   Im(t)   dist\n")
        for R, i, j, t, d in model.hops:
            if d <= rcut + 1e-6 and abs(t) > 1e-6:
                fh.write(f"{R[0]:3d} {R[1]:3d} {R[2]:3d} {i:3d} {j:3d} {t.real:+.8f} {t.imag:+.8f} {d:.4f}\n")
    with open(os.path.join(out, f"{tag}_log.txt"), "w") as fh:
        fh.write("\n".join(log) + "\n")
    x = bs.lineark(); ticks, labels = bs.lineark(ticks=True)[1:]
    fig, ax = plt.subplots(1, 2, figsize=(11, 5), gridspec_kw=dict(width_ratios=[1, 1.4]))
    for a, (lo, hi) in zip(ax, [(Er.min() - 2, min(Er.max(), 300)), (min(Em.min(), -12) - 1, max(Em.max(), 8) + 1)]):
        for n in range(Er.shape[1]):
            a.plot(x, Er[:, n], c="0.65", lw=0.8, label="SIESTA DZP" if n == 0 else None)
        for n in range(Em.shape[1]):
            a.plot(x, Em[:, n], c="#c0504d", lw=1.3, ls="--", label=f"modelo J={model.J}" if n == 0 else None)
        if np.isfinite(window[0]):
            a.axhspan(*window, color="#1f4e79", alpha=0.08, lw=0, label="ventana de Fermi")
        a.set_xticks(ticks); a.set_xticklabels(labels); a.set_xlim(x[0], x[-1]); a.set_ylim(lo, hi)
        [a.axvline(t, c="0.85", lw=0.5) for t in ticks]; a.axhline(0, c="0.4", ls=":", lw=0.6)
        a.set_ylabel("E − E$_F$ (eV)")
    ax[0].set_title(f"{tag}: espectro completo"); ax[1].set_title("zoom: modelo reducido")
    ax[1].legend(fontsize=8, frameon=False, loc="lower right")
    fig.suptitle(f"{tag} — {choice['kind']} {choice['block']}: " + "+".join(shell_name(k) for k in choice.get("combo", []))
                 + (" columnas: " + ", ".join(choice["centres"]) if "centres" in choice else ""), fontsize=10)
    fig.tight_layout(); fig.savefig(os.path.join(out, f"{tag}_bands.png"), dpi=120); plt.close(fig)


def build_parser(doc=__doc__):
    p = argparse.ArgumentParser(description=doc, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("rundir"); p.add_argument("label")
    p.add_argument("--out", default="downfold_out")
    p.add_argument("--n2d", type=int, default=24, help="k-grid per in-plane direction (2D)")
    p.add_argument("--n3d", type=int, default=10, help="k-grid per direction (3D)")
    p.add_argument("--w", type=float, default=1.0, help="Fermi window margin below VBM / above CBM (eV)")
    p.add_argument("--gap", type=float, default=0.05, help="minimum gap between groups (eV)")
    p.add_argument("--smin", type=float, default=0.15, help="minimum singular value of the projection")
    p.add_argument("--tol", type=float, default=0.010, help="target max band error after truncation (eV)")
    p.add_argument("--tol-ent", type=float, default=0.050, help="acceptable SCDM error in the Fermi window (eV)")
    p.add_argument("--jmax", type=int, default=40)
    p.add_argument("--emax-above", type=float, default=25.0, help="ignore groups reaching above E_F + this")
    p.add_argument("--emax-below", type=float, default=30.0, help="ignore groups reaching below E_F - this")
    p.add_argument("--outer-max", type=float, default=60.0,
                   help="highest energy (eV) of the extra states used to complete an entangled model")
    p.add_argument("--tau-grid", type=float, nargs="+", default=[1e9, 8.0, 4.0, 1.5],
                   help="decay (eV) of the weight of outer states above the Fermi window")
    p.add_argument("--mu-grid", type=float, nargs="+", default=[0.5, 1, 2, 3, 5, 8],
                   help="SCDM: erfc centre above the top of the Fermi window (eV)")
    p.add_argument("--sigma-grid", type=float, nargs="+", default=[0.25, 0.5, 1, 2], help="SCDM: erfc widths (eV)")
    p.add_argument("--nval", type=int, default=40, help="random validation k-points")
    return p


def main():
    p = build_parser()
    a = p.parse_args()
    H, orbinfo, nspin = load(a.rundir, a.label)
    res = [analyse_spin(H, orbinfo, s, a, a.out, a.label + (f"_spin{'UD'[s]}" if nspin == 2 else ""))
           for s in range(nspin)]
    with open(os.path.join(a.out, f"{a.label}_summary.json"), "w") as fh:
        json.dump(res, fh, indent=1, ensure_ascii=False)


if __name__ == "__main__":
    main()
