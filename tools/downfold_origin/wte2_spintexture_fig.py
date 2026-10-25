#!/usr/bin/env python
"""Spin texture of the lowest conduction band around Q (cf. Fig. 1(b) inset of PRL 125, 256603).
Paper axes: x = zigzag (our y, b), y = our x (a), z = out of plane.  k in units of pi/a_x, centred at Q.
(a) 1T' (P2_1/m): Kramers pairs are degenerate -> inside each pair take the state with spin +n
    along the locally conserved axis (largest eigenvector of the projected spin).
(b) 1Td-like (Pm): spin expectation of the lowest conduction band itself.
"""
import numpy as np, scipy.linalg as sl, sisl, warnings
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
warnings.filterwarnings("ignore")
NOCC = 92
SIG = [np.array([[0, 1], [1, 0]]), np.array([[0, -1j], [1j, 0]]), np.array([[1, 0], [0, -1]])]


def solve(H, k):
    Hk = H.Hk(k, format="array"); Sx = H.Sk(k, format="array")
    S = Sx[0::2, 0::2]
    e, c = sl.eigh(Hk, Sx)
    return e, c, [np.kron(S, s) for s in SIG]


def texture(run, kQ, nx, ny, half_x, half_y, mode):
    H = sisl.get_sile(f"{run}/WTe2.HSX").read_hamiltonian()
    a, b = H.geometry.cell[0, 0], H.geometry.cell[1, 1]
    # paper units pi/a_x with a_x = b:  q_x (zigzag) <-> our fractional ky = q_x/2 ; q_y <-> our kx = q_y*a/(2b)
    qx = np.linspace(-half_x, half_x, nx); qy = np.linspace(-half_y, half_y, ny)
    E = np.zeros((ny, nx)); Sv = np.zeros((ny, nx, 3))
    for iy, y in enumerate(qy):
        for ix, x in enumerate(qx):
            k = [kQ[0] + y * a / (2 * b), kQ[1] + x / 2, 0]
            e, c, Sm = solve(H, k)
            if mode == "pair":
                C2 = c[:, NOCC:NOCC + 2]
                S2 = [C2.conj().T @ M @ C2 for M in Sm]
                T = np.array([[np.real(np.trace(S2[i] @ S2[j])) / 2 for j in range(3)] for i in range(3)])
                n = np.linalg.eigh(T)[1][:, -1]
                if n[2] < 0: n = -n                                  # fix the global sign (+ state)
                w, v = np.linalg.eigh(sum(n[i] * S2[i] for i in range(3)))
                psi = C2 @ v[:, -1]
                E[iy, ix] = e[NOCC:NOCC + 2].mean()
            else:
                psi = c[:, NOCC]; E[iy, ix] = e[NOCC]
            s = np.array([np.real(psi.conj() @ M @ psi) for M in Sm])  # our axes
            Sv[iy, ix] = [s[1], s[0], s[2]]                            # -> paper axes (x_p=y, y_p=x, z)
    return qx, qy, E, Sv


def main():
    kQ = {"D_soc_sym": (0.0, 0.1633), "C_soc": (0.0, 0.1632)}
    panels = [("D_soc_sym", "pair", "(a) 1T′ (P2$_1$/m): estado +n de cada par de Kramers"),
              ("C_soc", "band", "(b) 1T$_d$-like (Pm): banda de conducción más baja")]
    fig, axs = plt.subplots(1, 2, figsize=(11.5, 4.6), constrained_layout=True)
    for ax, (run, mode, title) in zip(axs, panels):
        qx, qy, E, Sv = texture(run, kQ[run], 61, 33, 0.16, 0.08, mode)
        Emin = E.min()
        im = ax.pcolormesh(qx, qy, (E - Emin) * 1e3, cmap="viridis", vmin=0, vmax=80, shading="auto")
        ax.contour(qx, qy, (E - Emin) * 1e3, levels=[10, 20, 40, 60], colors="w", linewidths=0.6, alpha=0.7)
        st = (slice(1, None, 4), slice(2, None, 5))
        X, Y = np.meshgrid(qx, qy)
        ax.quiver(X[st], Y[st], Sv[..., 1][st], Sv[..., 2][st], color="#ff8c1a", scale=18, width=0.004)
        ax.plot(0, 0, "o", color="w", ms=5)
        ax.set_aspect("equal"); ax.set_xlabel("q$_x$ (π/a$_x$)  [zigzag]"); ax.set_ylabel("q$_y$ (π/a$_x$)")
        ax.set_title(title, fontsize=10)
        # statistics inside the pocket (E - Emin < 40 meV)
        m = (E - Emin) < 0.040
        ang = np.degrees(np.arctan2(Sv[..., 2][m], Sv[..., 1][m]))
        pol = np.linalg.norm(Sv[m], axis=1)
        print(f"{run} [{mode}]: ángulo en el plano yz (desde +y hacia +z) = {np.mean(ang):+.1f}° ± {np.std(ang):.1f}°, "
              f"|S_x| medio = {np.abs(Sv[...,0][m]).mean():.3f}, |S| medio = {pol.mean():.3f}, puntos = {m.sum()}")
        cb = fig.colorbar(im, ax=ax, shrink=0.85); cb.set_label("E − E$_{CBM}$ (meV)")
    fig.suptitle("WTe$_2$ monocapa, SIESTA PBE+SOC (DZP): textura de espín S$_y$, S$_z$ alrededor de Q", fontsize=11)
    fig.savefig("wte2_spin_texture_Q.png", dpi=150)


if __name__ == "__main__":
    main()
