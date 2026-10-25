#!/usr/bin/env python
"""
downfold_scdm.py -- same goal as downfold.py, but the localised functions are chosen by
SCDM-k (Damle, Lin & Ying, J. Chem. Theory Comput. 11, 1463 (2015); Multiscale Model. Simul.
16, 1392 (2018)) instead of projections on atomic shells.  No trial orbitals are given.

In the LCAO basis the "real-space grid" of SCDM is the set of Loewdin-orthonormalised
orbitals: psi~(k) = S(k)^{1/2} c(k).  QR with column pivoting of the real-space density matrix
rho(R=0) = <psi~ f^2 psi~^H>_k selects the J orbitals that best span the target subspace; their Bloch sums are the trial functions
at every k:  A_mn(k) = f(E_mk) conj(psi~_mk[col_n]),  U(k) = A (A^H A)^{-1/2},
H_W(k) = U^H diag(E) U.

  * isolated manifolds (per mirror block, combined as in downfold.py): f = 1 on the group.
  * entangled: f = erfc((E-mu)/sigma)/2, with E_semicore < E < --outer-max; mu, sigma and J are
    scanned on a grid (deterministic) and the result with the lowest KPM cost is kept.
Everything else (symmetry blocks, Fourier/Wigner-Seitz, truncation, validation, output, the
selection by non-zeros per cell) is shared with downfold.py.
"""
import itertools, json, os, sys
import numpy as np, scipy.linalg as sl
from scipy.special import erfc
import downfold as D


def sqrtm_h(S):
    w, v = np.linalg.eigh(S)
    return (v * np.sqrt(np.clip(w, 0, None))) @ v.conj().T


def lowdin_coeffs(spec, B, Vfull):
    """psi~ in the block's Loewdin basis for every k: (nk, nb_block, nstates)."""
    out = []
    for ik in range(len(spec.S)):
        Sb = B.T @ spec.S[ik] @ B
        c = B.T @ Vfull[ik]
        out.append(sqrtm_h(Sb) @ c)
    return np.array(out)


def scdm_columns(Pt, f, J):
    """QRCP of the real-space density matrix rho(R=0) = <P~(k)>_k in the Loewdin basis
    (k-averaged SCDM; more robust than Gamma alone for small cells / metals)."""
    rho0 = np.einsum("kob,kb,kpb->op", Pt, f ** 2, Pt.conj()) / len(Pt)
    _, _, piv = sl.qr(rho0, mode="economic", pivoting=True)
    return piv[:J]


def scdm_HW(Pt, E, f, cols):
    A = f[:, :, None] * Pt[:, cols, :].conj().transpose(0, 2, 1)       # (nk, nstates, J)
    u, s, wh = np.linalg.svd(A, full_matrices=False)
    U = u @ wh
    return np.einsum("kbi,kb,kbj->kij", U.conj(), E, U), s.min()


def analyse_spin(H, orbinfo, spin, args, out, tag):
    geom = H.geometry
    norb = len(orbinfo)
    blocks = D.mirror_blocks(geom, orbinfo)
    N, ks = D.kgrid(geom, args.n2d, args.n3d)
    kval = D.validation_kpoints(geom, N, args.nval)
    log = []
    say = lambda s: (print(s), log.append(s))
    say(f"\n=== {tag} [SCDM]: {norb} funciones de base, bloques: {[b for b, _ in blocks]}, malla {N}")
    spec = D.Spectrum(H, spin, blocks, ks)
    full = spec if blocks[0][0] == "full" else D.Spectrum(H, spin, [("full", np.eye(norb))], ks)
    Eall = full.E["full"]
    vbm, cbm, metal, Wf = D.fermi_window(Eall, args.w)
    n_in_W = ((Eall >= Wf[0]) & (Eall <= Wf[1])).sum(1)
    say(f"  VBM={vbm:+.3f} CBM={cbm:+.3f} {'metal' if metal else 'gap=%.3f eV' % (cbm - vbm)}  "
        f"ventana=[{Wf[0]:+.2f},{Wf[1]:+.2f}]  E_max espectro={Eall.max():.1f} eV")
    Eval_full = [D.ref_block(H, spin, np.eye(norb), k) for k in kval]
    Eval_blk = {b: [D.ref_block(H, spin, B, k) for k in kval] for b, B in blocks}
    options = []

    # ---------------- isolated groups: f = indicator of the group
    cands = {b: [] for b, _ in blocks}
    for b, B in blocks:
        E = spec.E[b]
        Pt = lowdin_coeffs(spec, B, spec.V[b])
        grps = [g for g in D.band_groups(E, args.gap) if g["Emin"] > -args.emax_below and g["Emax"] < args.emax_above]
        for a in range(len(grps)):
            for c in range(a, len(grps)):
                if any(grps[x + 1]["lo"] != grps[x]["hi"] + 1 for x in range(a, c)):
                    continue
                lo, hi = grps[a]["lo"], grps[c]["hi"]; J = hi - lo + 1
                if J > args.jmax:
                    continue
                Ec = E[:, lo:hi + 1]
                f = np.ones_like(Ec)
                cols = scdm_columns(Pt[:, :, lo:hi + 1], f, J)
                HWk, smin = scdm_HW(Pt[:, :, lo:hi + 1], Ec, f, cols)
                cands[b].append(dict(block=b, lo=lo, hi=hi, J=J, smin=smin, HWk=HWk, G=B[:, cols],
                                     inW=((Ec >= Wf[0]) & (Ec <= Wf[1])).sum(1)))
    for pick in itertools.product(*[[None] + cands[b] for b, _ in blocks]):
        pick = [c for c in pick if c is not None]
        if not pick or not np.all(sum(c["inW"] for c in pick) == n_in_W):
            continue
        smin = min(c["smin"] for c in pick)
        if smin < 1e-3:
            continue
        J = sum(c["J"] for c in pick); off = np.cumsum([0] + [c["J"] for c in pick])
        HWk = np.zeros((len(ks), J, J), complex)
        for c, o in zip(pick, off):
            HWk[:, o:o + c["J"], o:o + c["J"]] = c["HWk"]
        G = np.concatenate([c["G"] for c in pick], axis=1)
        Ev = [np.sort(np.concatenate([Eval_blk[c["block"]][ik][c["lo"]:c["hi"] + 1] for c in pick]))
              for ik in range(len(kval))]
        ef = (lambda Ev: (lambda m, rc: max(np.abs(np.linalg.eigvalsh(m.Hk(k, rc)) - Er).max()
                                            for k, Er in zip(kval, Ev))))(Ev)
        ch = dict(kind="SCDM aislado", block="+".join(f"{c['block']}({c['J']})" for c in pick), J=J,
                  combo=[], smin=smin)
        options.append((ch, HWk, G, ef))

    # ---------------- entangled: erfc window, scan mu, sigma, J
    SH = D.shells(orbinfo)
    semicore = {(o["species"], o["n"], o["l"], o["pol"]) for o in orbinfo if o["semi"]}
    Elow = -np.inf
    V = full.V["full"]
    if semicore:
        Gs = D.trials_in_block(sum((SH[k] for k in semicore), []), np.eye(norb), norb)
        As = D.projection_matrix(V, full.S, Gs); p = np.empty(Eall.shape)
        for ik in range(len(ks)):
            X = sl.fractional_matrix_power(Gs.T @ full.S[ik] @ Gs, -0.5)
            p[ik] = np.sum(np.abs(As[ik] @ X) ** 2, axis=1).real
        for g in D.band_groups(Eall, args.gap):
            if g["Emax"] < Wf[0] and p[:, g["lo"]:g["hi"] + 1].mean() > 0.5:
                Elow = max(Elow, g["Emax"])
        say(f"  semicore excluido: E <= {Elow:+.1f} eV")
    Pt = lowdin_coeffs(full, np.eye(norb), V)
    werr = lambda m, rc: D.window_error(m, kval, Eval_full, Wf, rc)
    step = 2 if any(o["spin"] is not None for o in orbinfo) else 1
    Jmin = int(n_in_W.max()); Jmin += (-Jmin) % step
    Js = list(range(max(Jmin, step), Jmin + args.jextra + 1, step))
    ent = []
    for dmu in args.mu_grid:
        for sig in args.sigma_grid:
            mu = Wf[1] + dmu
            f = 0.5 * erfc((Eall - mu) / sig) * (Eall > Elow) * (Eall < args.outer_max)
            piv = scdm_columns(Pt, f, max(Js))
            for J in Js:
                cols = piv[:J]
                HWk, smin = scdm_HW(Pt, Eall, f, cols)
                if smin < 1e-3:
                    continue
                m = D.Model(geom, D.centers_of(np.eye(norb)[:, cols], geom, orbinfo), N, HWk, ks)
                e = werr(m, np.inf)
                ent.append((e, J, mu, sig, smin, HWk, cols, m))
    ent.sort(key=lambda x: (x[0], x[1]))
    say(f"  SCDM con ventana: {len(ent)} combinaciones (μ, σ, J) evaluadas; las mejores:")
    for e, J, mu, sig, smin, HWk, cols, m in ent[:args.keep]:
        lab = sorted(set(f"{orbinfo[c]['species']}-{orbinfo[c]['n']}{D.LSYM[orbinfo[c]['l']]}" for c in cols))
        ch = dict(kind="SCDM ventana", block=f"μ=EF{mu:+.1f},σ={sig:g}", J=J, combo=[], smin=smin,
                  cols=[int(c) for c in cols], centres=lab)
        say(f"    J={J:2d} μ=E_F{mu:+.2f} σ={sig:<4g} σ_min={smin:.3f} err(sin truncar)={e*1e3:7.1f} meV  "
            f"columnas en: {', '.join(lab)}")
        options.append((ch, HWk, np.eye(norb)[:, cols], werr, m))

    if not options:
        say("  >>> SCDM no encuentra ningún modelo válido para este sistema")
        return dict(tag=tag, kind="ninguno", block="", J=0, smin=0.0, rcut=0.0, err_meV=float("nan"),
                    width_model=float("nan"), width_siesta=float(Eall.max() - Eall.min()))
    choice, HWk, G, errfn = D.select_option(options, geom, orbinfo, N, ks, args, say)
    window = Wf if choice["kind"] == "SCDM ventana" else (-np.inf, np.inf)
    say(f"  >>> ELEGIDO: {choice['kind']} {choice['block']} J={choice['J']}  σ_min={choice['smin']:.3f}")
    model = D.Model(geom, D.centers_of(G, geom, orbinfo), N, HWk, ks)
    rc, e = D.rcut_for(model, errfn, args.tol)
    if rc is None:
        rc = model.shells()[-1]; e = errfn(model, rc)
    nhop = ((model.d > 1e-3) & (model.d <= rc + 1e-6)).sum() / model.J
    say(f"  >>> r_cut={rc:.2f} Å  ({nhop:.0f} hoppings/orbital)  error_max={e*1e3:.1f} meV")
    bs = D.kpath(geom)
    Em = np.array([np.linalg.eigvalsh(model.Hk(k, rc)) for k in bs.k])
    Er = np.array([D.ref_block(H, spin, np.eye(norb), k) for k in bs.k])
    say(f"  ancho del espectro: SIESTA {Er.max()-Er.min():.1f} eV  ->  modelo {Em.max()-Em.min():.2f} eV")
    choice.setdefault("combo", [])
    D.write_outputs(out, tag, model, rc, choice, bs, Er, Em, window, log, geom, orbinfo)
    return dict(tag=tag, kind=choice["kind"], block=choice["block"], J=int(choice["J"]), smin=float(choice["smin"]),
                rcut=float(rc), err_meV=float(e * 1e3), width_model=float(Em.max() - Em.min()),
                width_siesta=float(Er.max() - Er.min()))


def main():
    p = D.build_parser(__doc__)
    p.add_argument("--jextra", type=int, default=12, help="SCDM: scan J from the window count up to this many more")
    p.add_argument("--keep", type=int, default=4, help="SCDM: windowed candidates passed to the final selection")
    a = p.parse_args()
    if a.out == "downfold_out":
        a.out = "downfold_scdm_out"
    H, orbinfo, nspin = D.load(a.rundir, a.label)
    res = [analyse_spin(H, orbinfo, s, a, a.out, a.label + "_scdm" + (f"_spin{'UD'[s]}" if nspin == 2 else ""))
           for s in range(nspin)]
    with open(os.path.join(a.out, f"{a.label}_scdm_summary.json"), "w") as fh:
        json.dump(res, fh, indent=1, ensure_ascii=False)


if __name__ == "__main__":
    main()
