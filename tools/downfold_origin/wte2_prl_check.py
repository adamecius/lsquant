#!/usr/bin/env python
"""Contrast the claims of Garcia et al. PRL 125, 256603 (2020) with SIESTA PBE+SOC (DZP).
Axes: paper x (zigzag, a_x) = our y (b = 3.496 A); paper y (a_y) = our x (a = 6.282 A); z = z.
  1. Z2 (Fu-Kane parities) and band inversion at Gamma, on the P2_1/m (symmetrised) geometry.
  2. CBM position Q.
  3. Spin texture of the lowest conduction band around Q:
     - 1Td-like geometry (Pm, inversion broken by ~1 mA): direct <sigma> of each band.
     - 1T' (P2_1/m): Kramers pairs are degenerate -> axis n that diagonalises the projected
       spin P sigma P inside the pair (the conserved spin axis if it exists).
"""
import numpy as np, scipy.linalg as sl, sisl, spglib, warnings
warnings.filterwarnings("ignore")
NOCC = 92                         # 2 W x 14 + 4 Te x 16 = 92 electrons -> 92 occupied spinor bands
SIG = [np.array([[0, 1], [1, 0]]), np.array([[0, -1j], [1j, 0]]), np.array([[1, 0], [0, -1]])]


def load(d):
    return sisl.get_sile(f"{d}/WTe2.HSX").read_hamiltonian()


def eig(H, k):
    Hk = H.Hk(k, format="array"); Sk = H.Sk(k, format="array")
    if Sk.shape[0] != Hk.shape[0]:
        Sk = np.kron(Sk, np.eye(2))
    e, c = sl.eigh(Hk, Sk)
    return e, c, Sk


def spin_mats(H, k):
    S = H.Sk(k, format="array")
    S = S[0::2, 0::2] if S.shape[0] == 2 * H.no else S          # sisl returns S (x) 1, spins interleaved
    return [np.kron(S, s) for s in SIG]


def orb_weights(H, c, Sk, sel):
    """Mulliken weight of states c on the orbital subset sel (spinor indices)."""
    w = np.real(np.conj(c) * (Sk @ c))
    return w[sel].sum(0)


def parity_matrix(H, k, center):
    g = H.geometry; no = H.no
    P = np.zeros((2 * no, 2 * no), complex)
    fx = g.fxyz; cf = center
    for ia in range(g.na):
        target = (2 * cf - fx[ia])
        for ja in range(g.na):
            if g.atoms[ja].Z != g.atoms[ia].Z: continue
            d = target - fx[ja]
            if np.allclose(d - np.round(d), 0, atol=2e-3):
                R = np.round(d)                         # lattice vector (fractional)
                break
        else:
            raise RuntimeError("no inversion image")
        ph = np.exp(2j * np.pi * np.dot(k, R))
        for io, o in enumerate(g.atoms[ia].orbitals):
            i = g.a2o(ia) + io; j = g.a2o(ja) + io
            for s in range(2):
                P[2 * j + s, 2 * i + s] = (-1) ** o.l * ph
    return P


def main():
    # ---------------- 1. Z2 and inversion at Gamma on the symmetric geometry
    H = load("D_soc_sym"); g = H.geometry
    ds = spglib.get_symmetry_dataset((g.cell, g.fxyz, [a.Z for a in g.atoms]), symprec=1e-4)
    inv = [t for r, t in zip(ds.rotations, ds.translations) if np.allclose(r, -np.eye(3))]
    center = inv[0] / 2
    print(f"[1] grupo {ds.international}; centro de inversión (frac.) = {np.round(center, 4)}")
    TRIM = {"Γ": [0, 0, 0], "X(nuestro)": [.5, 0, 0], "Y(nuestro)": [0, .5, 0], "S": [.5, .5, 0]}
    Wd = [2 * i + s for i in range(H.no) for s in range(2)
          if g.atoms[g.o2a(i)].Z == 74 and g.atoms[g.o2a(i)].orbitals[i - g.a2o(g.o2a(i))].l == 2]
    Tep = [2 * i + s for i in range(H.no) for s in range(2)
           if g.atoms[g.o2a(i)].Z == 52 and g.atoms[g.o2a(i)].orbitals[i - g.a2o(g.o2a(i))].l == 1]
    prod_all = 1
    for name, k in TRIM.items():
        e, c, Sk = eig(H, k)
        P = parity_matrix(H, k, center)
        # parity eigenvalues inside each degenerate block (non-symmorphic sticking at zone boundary)
        xi = np.zeros(len(e)); b = 0
        while b < len(e):
            nb = 1
            while b + nb < len(e) and abs(e[b + nb] - e[b]) < 1e-4:
                nb += 1
            Cb = c[:, b:b + nb]
            Pb = Cb.conj().T @ Sk @ P @ Cb
            xi[b:b + nb] = np.sort(np.real(np.linalg.eigvals(Pb)))
            if b < NOCC < b + nb:
                print(f"    aviso: bloque degenerado cortado por NOCC en {name}")
            b += nb
        prod = np.prod(np.sign(xi[:NOCC]))                   # each Kramers pair counted twice -> use sqrt
        npos = (xi[:NOCC] > 0).sum(); nneg = (xi[:NOCC] < 0).sum()
        prod = (-1) ** (nneg // 2); prod_all *= prod
        print(f"    {name:11s}: ξ>0: {npos}, ξ<0: {nneg} (estados, contando ambos de cada par)  -> δ = {int(prod):+d}  (|ξ| mín = {np.abs(xi[:NOCC]).min():.3f})")
        if name == "Γ":
            print("    estados en Γ cerca de E_F (E-EF, paridad, peso W-d, peso Te-p):")
            wd = orb_weights(H, c, Sk, Wd); wp = orb_weights(H, c, Sk, Tep)
            for b in range(NOCC - 6, NOCC + 4, 2):
                tag = "ocup." if b < NOCC else "vacío"
                print(f"      banda {b+1:3d}-{b+2:3d} {tag}: {e[b]:+.3f} eV  ξ={xi[b]:+.2f}  W-d={wd[b]:.2f}  Te-p={wp[b]:.2f}")
    print(f"    ν(Z2) = {0 if prod_all > 0 else 1}  (producto total {int(prod_all):+d})")

    # ---------------- 2. CBM position
    for d in ["D_soc_sym", "C_soc"]:
        Hd = load(d)
        n1, n2 = 24, 120
        best = (9, None)
        for i in range(-n1 // 2, n1 // 2 + 1):
            for j in range(-n2 // 2, n2 // 2 + 1):
                k = [i / n1, j / n2, 0]
                e = eig(Hd, k)[0][NOCC]
                if e < best[0]: best = (e, k)
        # refine along both axes
        from scipy.optimize import minimize
        f = lambda kk: eig(Hd, [kk[0], kk[1], 0])[0][NOCC]
        r = minimize(f, best[1][:2], method="Nelder-Mead", options=dict(xatol=1e-5, fatol=1e-6))
        kq = r.x
        print(f"[2] {d}: CBM en k = ({kq[0]:+.4f}, {kq[1]:+.4f}) (frac. nuestros); E_CBM-E_F = {r.fun:+.3f} eV")
        print(f"    en unidades del PRL (π/a_x a lo largo del zigzag = nuestro b): q_x = {2*abs(kq[1]):.3f},"
              f" q_y = {2*abs(kq[0]):.3f}   [PRL: Q = (0.266, 0.0)]")
        if d == "C_soc":
            kQ_C = kq
        else:
            kQ_D = kq

    # ---------------- 3. spin texture around Q
    def texture(d, kq, nondeg):
        Hd = load(d); out = []
        for dx in np.linspace(-0.02, 0.02, 5):
            for dy in np.linspace(-0.03, 0.03, 7):
                k = [kq[0] + dx, kq[1] + dy, 0]
                e, c, Sk = eig(Hd, k); Sm = spin_mats(Hd, k)
                if nondeg:
                    v = np.array([np.real(c[:, NOCC].conj() @ M @ c[:, NOCC]) for M in Sm])
                else:
                    C2 = c[:, NOCC:NOCC + 2]
                    S2 = [C2.conj().T @ M @ C2 for M in Sm]          # projected spin in the Kramers pair
                    # axis n maximising the spin polarisation: largest eigvec of sum_a S2_a S2_a^dagger
                    Tm = np.array([[np.real(np.trace(S2[a] @ S2[b])) / 2 for b in range(3)] for a in range(3)])
                    w, V = np.linalg.eigh(Tm); v = V[:, -1] * np.sqrt(max(w[-1], 0)); v *= np.sign(v[2] or 1)
                out.append((e[NOCC], v))
        return out

    for d, kq, nondeg in [("C_soc", kQ_C, True), ("D_soc_sym", kQ_D, False)]:
        T = texture(d, kq, nondeg)
        V = np.array([v for _, v in T])
        # paper axes: x_p = our y, y_p = our x, z_p = z
        Vp = np.stack([V[:, 1], V[:, 0], V[:, 2]], axis=1)
        n = Vp / np.linalg.norm(Vp, axis=1)[:, None]
        n *= np.sign(n[:, 2])[:, None]
        theta = np.degrees(np.arctan(n[:, 2] / n[:, 1]))          # axis angle from y_p in the y_p-z_p plane, in (-90, 90]
        kind = "<σ> banda de conducción más baja" if nondeg else "eje de espín conservado en el par de Kramers"
        print(f"[3] {d}: {kind}, {len(T)} puntos alrededor de Q")
        print(f"    componentes medias (ejes PRL) |n_x|={np.abs(n[:,0]).mean():.3f}  n_y={n[:,1].mean():+.3f}  n_z={n[:,2].mean():+.3f};"
              f"  |<σ>| medio = {np.linalg.norm(V,axis=1).mean():.3f}")
        print(f"    ángulo respecto a y (plano yz): {theta.mean():+.1f}° ± {theta.std():.1f}°  "
              f"(máx. dispersión {np.ptp(theta):.1f}°)   [PRL: −56°; experimentos ~50–59°]")


if __name__ == "__main__":
    main()
