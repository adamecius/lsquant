#!/usr/bin/env python
"""Spin texture around Q: SIESTA (exact, P2_1/m) vs the 44-spinor Wannier model built from it.
Same Kramers-pair projection in both (state +n of each degenerate pair).
Wannier model: spin taken as 1 (x) sigma in the Wannier basis (projection-only WFs close to the
atomic spinor trials; the exact operator would need the .spn matrix elements)."""
import os, sys, numpy as np, scipy.linalg as sl, sisl, warnings
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
sys.path.insert(0, os.path.expanduser("~/wte2"))
from compare_w90 import read_hr
warnings.filterwarnings("ignore")
NOCC_SIESTA, NOCC_MODEL = 92, 92 - 64
SIG = [np.array([[0, 1], [1, 0]]), np.array([[0, -1j], [1j, 0]]), np.array([[1, 0], [0, -1]])]
kQ = (0.0, 0.1632)

RUN = sys.argv[1] if len(sys.argv) > 1 else "C_soc"
DAT = sys.argv[2] if len(sys.argv) > 2 else "df_C_soc/WTe2_hoppings.dat"
TAG = sys.argv[3] if len(sys.argv) > 3 else "Pm"
H = sisl.get_sile(os.path.expanduser(f"~/wte2/{RUN}/WTe2.HSX")).read_hamiltonian()
a, b = H.geometry.cell[0, 0], H.geometry.cell[1, 1]
rows = np.array([l.split() for l in open(os.path.expanduser(f"~/wte2/{DAT}")) if not l.startswith("#") and l.strip()], float)
RR = rows[:, :3]; AA = rows[:, 3].astype(int); BB = rows[:, 4].astype(int); WW = rows[:, 5] + 1j * rows[:, 6]
nw = int(max(AA.max(), BB.max())) + 1
Sw = [np.kron(np.eye(nw // 2), s) for s in SIG]


def pair_state(e, v, n0, Sm):
    C2 = v[:, n0:n0 + 2]
    S2 = [C2.conj().T @ M @ C2 for M in Sm]
    T = np.array([[np.real(np.trace(S2[i] @ S2[j])) / 2 for j in range(3)] for i in range(3)])
    nn = np.linalg.eigh(T)[1][:, -1]
    if nn[2] < 0: nn = -nn
    w, u = np.linalg.eigh(sum(nn[i] * S2[i] for i in range(3)))
    psi = C2 @ u[:, -1]
    s = np.array([np.real(psi.conj() @ M @ psi) for M in Sm])
    return e[n0:n0 + 2].mean(), np.array([s[1], s[0], s[2]])            # paper axes (x_p = our y)


def siesta(k):
    Hk = H.Hk(k, format="array"); Sx = H.Sk(k, format="array"); S = Sx[0::2, 0::2]
    e, v = sl.eigh(Hk, Sx)
    return pair_state(e, v, NOCC_SIESTA, [np.kron(S, s) for s in SIG])


def model_eig(k):
    Hk = np.zeros((nw, nw), complex)
    np.add.at(Hk, (AA, BB), WW * np.exp(2j * np.pi * (RR @ np.asarray(k))))
    return np.linalg.eigh(0.5 * (Hk + Hk.conj().T))


def model(k):
    e, v = model_eig(k)
    return pair_state(e, v, NOCC_MODEL, Sw)


kq0 = [kQ[0], kQ[1], 0]
eS = sl.eigh(H.Hk(kq0, format="array"), H.Sk(kq0, format="array"), eigvals_only=True)[NOCC_SIESTA]
eM = model_eig(kq0)[0]; NOCC_MODEL = int(np.argmin(np.abs(eM - eS)))
NOCC_MODEL -= NOCC_MODEL % 2                                               # start of the Kramers pair
print(f"modelo downfold: J = {nw}, par de conducción = bandas {NOCC_MODEL}-{NOCC_MODEL+1}; E(Q): SIESTA {eS*1e3:+.1f}, modelo {eM[NOCC_MODEL]*1e3:+.1f} meV")
qx = np.linspace(-0.16, 0.16, 61); qy = np.linspace(-0.08, 0.08, 33)
E = {"s": np.zeros((33, 61)), "m": np.zeros((33, 61))}; S = {"s": np.zeros((33, 61, 3)), "m": np.zeros((33, 61, 3))}
for iy, y in enumerate(qy):
    for ix, x in enumerate(qx):
        k = [kQ[0] + y * a / (2 * b), kQ[1] + x / 2, 0]
        E["s"][iy, ix], S["s"][iy, ix] = siesta(k)
        E["m"][iy, ix], S["m"][iy, ix] = model(k)

ang = {k: np.degrees(np.arctan2(S[k][..., 2], S[k][..., 1])) for k in S}
dth = (ang["m"] - ang["s"] + 90) % 180 - 90                       # axis difference in (-90, 90]
Emin = E["s"].min(); pocket = (E["s"] - Emin) < 0.040
dE = (E["m"] - E["s"]) * 1e3
print(f"bolsillo (E-E_CBM < 40 meV): {pocket.sum()} puntos")
for k, name in [("s", "SIESTA"), ("m", "downfold.py")]:
    print(f"  {name:11s}: ángulo = {ang[k][pocket].mean():+.1f}° ± {ang[k][pocket].std():.1f}°, "
          f"|S_x| = {np.abs(S[k][...,0][pocket]).mean():.3f}, |S| = {np.linalg.norm(S[k][pocket],axis=1).mean():.3f}, "
          f"E_CBM-E_F = {E[k].min()*1e3:+.1f} meV")
print(f"  Δángulo (downfold−SIESTA): medio {np.abs(dth[pocket]).mean():.2f}°, máx {np.abs(dth[pocket]).max():.2f}°")
print(f"  ΔE en el bolsillo: RMS {np.sqrt((dE[pocket]**2).mean()):.2f} meV, máx {np.abs(dE[pocket]).max():.2f} meV")

fig, axs = plt.subplots(1, 3, figsize=(16, 4.6), constrained_layout=True)
X, Y = np.meshgrid(qx, qy); st = (slice(1, None, 4), slice(2, None, 5))
for ax, k, title in [(axs[0], "s", "(a) SIESTA (PBE+SOC, DZP)"), (axs[1], "m", f"(b) Modelo downfold.py, J = {nw} espinores")]:
    im = ax.pcolormesh(qx, qy, (E[k] - Emin) * 1e3, cmap="viridis", vmin=0, vmax=80, shading="auto")
    ax.contour(qx, qy, (E[k] - Emin) * 1e3, levels=[10, 20, 40, 60], colors="w", linewidths=0.6, alpha=0.7)
    ax.quiver(X[st], Y[st], S[k][..., 1][st], S[k][..., 2][st], color="#ff8c1a", scale=18, width=0.004)
    ax.plot(0, 0, "o", color="w", ms=5); ax.set_aspect("equal"); ax.set_title(title, fontsize=10)
    ax.set_xlabel("q$_x$ (π/a$_x$)  [zigzag]"); ax.set_ylabel("q$_y$ (π/a$_x$)")
    fig.colorbar(im, ax=ax, shrink=0.8).set_label("E − E$_{CBM}$ (meV)")
im = axs[2].pcolormesh(qx, qy, np.where(pocket, dth, np.nan), cmap="RdBu_r", vmin=-2, vmax=2, shading="auto")
axs[2].contour(qx, qy, (E["s"] - Emin) * 1e3, levels=[10, 20, 40], colors="k", linewidths=0.5, alpha=0.5)
axs[2].set_aspect("equal"); axs[2].set_title("(c) Δθ = θ$_{downfold}$ − θ$_{SIESTA}$ en el bolsillo", fontsize=10)
axs[2].set_xlabel("q$_x$ (π/a$_x$)  [zigzag]"); axs[2].set_ylabel("q$_y$ (π/a$_x$)")
fig.colorbar(im, ax=axs[2], shrink=0.8).set_label("Δθ (°)")
fig.suptitle(f"WTe$_2$ monocapa ({TAG}): textura de espín S$_y$, S$_z$ alrededor de Q — SIESTA (referencia) frente a downfold.py", fontsize=11)
fig.savefig(os.path.expanduser(f"~/wte2/wte2_spintex_siesta_vs_downfold_{TAG}.png"), dpi=150)

np.savez(os.path.expanduser(f"~/wte2/spintex_grid_{TAG}.npz"), qx=qx, qy=qy, Es=E["s"], Em=E["m"], Ss=S["s"], Sm=S["m"])
