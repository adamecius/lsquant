#!/usr/bin/env python3
"""Bulk Haldane model in k-space: bands and Chern number.

This companion explains the quantized plateau that the real-space transport
calculation produces. It builds the two-band Bloch Hamiltonian of the Haldane
model and reads the Chern number of the lower band off a k-grid (Fukui-Hatsugai
gauge-invariant sum). It then sweeps the sublattice mass to show the topological
transition: the Chern number jumps from a quantized value to zero when the gap
closes at M = 3*sqrt(3)*t2*sin(phi).

The real-space LinQT calculation (see lsqhall.py / README) reproduces the same
quantized number as a Hall conductivity plateau in the gap.

    python haldane_kspace.py
"""

import numpy as np


def bloch(k, t1=-1.0, t2=0.15, phi=np.pi / 2, mass=0.0):
    """Two-band Bloch Hamiltonian H(k) for the honeycomb Haldane model."""
    a1 = np.array([1.5, np.sqrt(3) / 2])
    a2 = np.array([1.5, -np.sqrt(3) / 2])
    nn = [np.array([0.5, np.sqrt(3) / 2]),
          np.array([0.5, -np.sqrt(3) / 2]),
          np.array([-1.0, 0.0])]                      # A -> B
    nnn = [a1, a2 - a1, -a2]                           # same-sublattice, one chirality
    h0 = 2 * t2 * np.cos(phi) * sum(np.cos(k @ v) for v in nnn)
    hx = t1 * sum(np.cos(k @ v) for v in nn)
    hy = t1 * sum(np.sin(k @ v) for v in nn)
    hz = mass - 2 * t2 * np.sin(phi) * sum(np.sin(k @ v) for v in nnn)
    return h0 * np.eye(2) + np.array([[hz, hx - 1j * hy],
                                      [hx + 1j * hy, -hz]])


def chern_lower(Ng=24, **kw):
    """Chern number of the lower band by the Fukui-Hatsugai lattice method."""
    b1 = 2 * np.pi * np.array([1 / 3, 1 / np.sqrt(3)])
    b2 = 2 * np.pi * np.array([1 / 3, -1 / np.sqrt(3)])
    u = np.empty((Ng, Ng), object)
    for i in range(Ng):
        for j in range(Ng):
            k = (i / Ng) * b1 + (j / Ng) * b2
            _, v = np.linalg.eigh(bloch(k, **kw))
            u[i, j] = v[:, 0]
    ov = lambda a, b: np.vdot(a, b)
    F = 0.0
    for i in range(Ng):
        for j in range(Ng):
            a, b = u[i, j], u[(i + 1) % Ng, j]
            c, d = u[(i + 1) % Ng, (j + 1) % Ng], u[i, (j + 1) % Ng]
            F += np.angle(ov(a, b) * ov(b, c) / ov(d, c) / ov(a, d))
    return F / (2 * np.pi)


def gap(Ng=48, **kw):
    """Minimum direct gap over the Brillouin zone."""
    b1 = 2 * np.pi * np.array([1 / 3, 1 / np.sqrt(3)])
    b2 = 2 * np.pi * np.array([1 / 3, -1 / np.sqrt(3)])
    g = np.inf
    for i in range(Ng):
        for j in range(Ng):
            k = (i / Ng) * b1 + (j / Ng) * b2
            w = np.linalg.eigvalsh(bloch(k, **kw))
            g = min(g, w[1] - w[0])
    return g


def figure_phase_diagram(t2=0.15, phi=np.pi / 2):
    import matplotlib.pyplot as plt
    masses = np.linspace(0.0, 1.6, 17)
    Cs = [round(chern_lower(mass=m, t2=t2, phi=phi)) for m in masses]
    Mc = 3 * np.sqrt(3) * t2 * abs(np.sin(phi))      # gap-closing mass
    fig, ax = plt.subplots(figsize=(7, 3.2))
    ax.plot(masses, Cs, "o-", lw=1.0)
    ax.axvline(Mc, color="grey", ls=":", lw=1.0,
               label=f"gap closes at M = {Mc:.3f}")
    ax.set_xlabel("sublattice mass  M  (eV)")
    ax.set_ylabel("Chern number  C")
    ax.set_title(f"Haldane phase diagram  (t2 = {t2:g}, phi = pi/2)")
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig("fig_phase_diagram.png", dpi=150)
    print("wrote fig_phase_diagram.png")


if __name__ == "__main__":
    print(f"Chern (t2=0.15, phi=pi/2, M=0):    {chern_lower():+.3f}  "
          f"(topological, |C| = 1)")
    print(f"Chern (phi=0, trivial):            {chern_lower(phi=0.0):+.3f}")
    print(f"gap at golden parameters:          {gap():.3f} eV")
    figure_phase_diagram()
