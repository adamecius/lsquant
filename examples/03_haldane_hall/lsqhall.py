#!/usr/bin/env python3
"""lsqhall - Hall and longitudinal conductivity of the Haldane model, over LinQT.

The heavy work stays in C++, driven through the single `lsquant` binary:

    lsquant compute --config run.json the two-operator Chebyshev moments
                                      (VX-VY for the Hall response, VX-VX for the
                                      longitudinal response); bounds come from
                                      the operators/haldane.HAM.desc sidecar
    lsquant reconstruct ... bastin    the Kubo-Bastin sum    -> sigma_xy(E)
    lsquant reconstruct ... greenwood the Kubo-Greenwood sum -> sigma_xx(E)
    inline_kuboBastinSea/SurfFromChebmom   the Fermi-sea and Fermi-surface parts

Python stages the operators, launches the kernels, reads the small reconstructed
curves, reads the Chern number off the gap plateau, and plots.

The Haldane model (t1 = -1, t2 = 0.15, phi = pi/2) is a Chern insulator with
C = +1. The committed operators in test/golden/haldane are the verified model;
copy them into ./operators (see README). Units follow the golden: e = hbar = 1,
so e^2/h = 1/(2 pi).

Run from inside the example folder, after LinQT is built and installed
(executables on PATH, or LSQUANT_BIN set to the build directory).
"""

import os
import json
import glob
import shutil
import tempfile
import subprocess

import numpy as np

A_CELL = 3 * np.sqrt(3) / 2          # honeycomb unit-cell area
E2_OVER_H = 1.0 / (2 * np.pi)        # conductance quantum in e = hbar = 1 units


def _bin(name):
    root = os.environ.get("LSQUANT_BIN", "")
    if root and os.path.exists(os.path.join(root, name)):
        return os.path.join(root, name)
    found = shutil.which(name)
    if found:
        return found
    raise FileNotFoundError(
        f"could not find '{name}'. Build and install LinQT, or set LSQUANT_BIN.")


def _stage(wd, opdir, label, ops):
    """Stage the operators and the exact-trace state into a run directory.

    The Hamiltonian descriptor (operators/<label>.HAM.desc) travels with the operators, so
    `lsquant compute` reads the spectral bounds straight from the sidecar -- no BOUNDS file.
    """
    os.makedirs(os.path.join(wd, "operators"), exist_ok=True)
    for op in ops:
        shutil.copy(os.path.join(opdir, f"{label}.{op}.CSR"),
                    os.path.join(wd, "operators", f"{label}.{op}.CSR"))
    shutil.copy(os.path.join(opdir, f"{label}.HAM.desc"),
                os.path.join(wd, "operators", f"{label}.HAM.desc"))
    dim = int(open(os.path.join(opdir, f"{label}.HAM.CSR")).readline().split()[0])
    with open(os.path.join(wd, "exact"), "w") as f:        # exact trace
        f.write("local\n" + f"{dim}\n" + "\n".join(map(str, range(dim))) + "\n")
    return dim


def _moments(wd, label, opr, opl, M):
    """Compute the two-operator moments through the single binary: lsquant compute --config.

    The run.json carries the run; compute reads the spectral bounds from the staged .desc. This
    is byte-identical to the legacy `inline_compute-kpm-nonEqOp <label> <opr> <opl> <M> exact`.
    """
    cfg = os.path.join(wd, f"run_{opr}_{opl}.json")
    with open(cfg, "w") as f:
        json.dump({"label": label, "operator_right": opr, "operator_left": opl,
                   "num_moments": M, "kernel": "jackson", "state": "exact"}, f, indent=2)
    subprocess.run([_bin("lsquant"), "compute", "--config", cfg],
                   cwd=wd, check=True,
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    return os.path.basename(sorted(
        glob.glob(os.path.join(wd, f"NonEqOp{opr}-{opl}{label}*chebmom2D")))[-1])


def _reconstruct(wd, mom, argv, out_glob, kernel=None):
    """Run a reconstruction kernel. `kernel` is set for the single-binary path
    (lsquant reconstruct <mom> <kernel> <broad>); the inline drivers take
    <mom> <broad> only."""
    full = list(argv) + ([mom, kernel, "10"] if kernel else [mom, "10"])
    subprocess.run(full, cwd=wd, check=True,
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    dat = sorted(glob.glob(os.path.join(wd, out_glob)))[-1]
    d = np.loadtxt(dat)
    os.remove(dat)
    return d[:, 0], d[:, 1]


def haldane_conductivities(opdir, label="haldane", M=128):
    """Return dim and the conductivity curves as {name: (E, sigma)}.

    Each curve keeps its OWN energy grid: the Kubo-Bastin route (sigma_xy, sea, surf) uses a
    fine 30*M grid, while Kubo-Greenwood (sigma_xx) uses an M-point grid, so they are not the
    same length.
    """
    lsq = [_bin("lsquant"), "reconstruct"]
    with tempfile.TemporaryDirectory() as wd:
        dim = _stage(wd, opdir, label, ["HAM", "VX", "VY"])
        xy = _moments(wd, label, "VX", "VY", M)                       # Hall moments
        Exy, sxy = _reconstruct(wd, xy, lsq, "KuboBastin*haldane*.dat", "bastin")
        Esea, sea = _reconstruct(wd, xy, [_bin("inline_kuboBastinSeaFromChebmom")],
                                 "KuboBastin*haldane*.dat")
        Esurf, surf = _reconstruct(wd, xy, [_bin("inline_kuboBastinSurfFromChebmom")],
                                   "KuboBastin*haldane*.dat")

        _stage(wd, opdir, label, ["HAM", "VX"])                       # longitudinal
        xx = _moments(wd, label, "VX", "VX", M)
        Exx, sxx = _reconstruct(wd, xx, lsq, "KuboGreenWood*haldane*.dat", "greenwood")
    return dim, {"sigma_xy": (Exy, sxy), "sea": (Esea, sea),
                 "surf": (Esurf, surf), "sigma_xx": (Exx, sxx)}


def chern_from_plateau(dim, E, sigma_xy, gap_halfwidth=0.5):
    """Chern number from the gap plateau: C = (plateau / A) * 2 pi.

    Matches hall_haldane.sh: A = (dim/2) * A_cell, e = hbar = 1. The Kubo-Bastin sigma_xy is
    area-extensive, so the plateau divided by A and multiplied by 2 pi reads off the integer.
    """
    area = (dim / 2) * A_CELL
    plateau = np.mean(sigma_xy[np.abs(E) < gap_halfwidth])
    return (plateau / area) * 2 * np.pi


def figure(opdir, label="haldane", M=128):
    import matplotlib.pyplot as plt
    dim, c = haldane_conductivities(opdir, label, M)
    norm = (dim / 2) * A_CELL * E2_OVER_H          # divide to read off e^2/h units
    Exy, sxy = c["sigma_xy"]
    C = chern_from_plateau(dim, Exy, sxy)

    fig, (top, bot) = plt.subplots(2, 1, figsize=(7, 6.6), sharex=True)

    for key, lab in [("sigma_xy", "Bastin (total)"), ("sea", "Fermi sea"),
                     ("surf", "Fermi surface")]:
        E, s = c[key]
        top.plot(E, s / norm, lw=1.0 if key != "sigma_xy" else 1.2, label=lab)
    top.axhline(1.0, color="grey", ls=":", lw=0.8)
    top.set_ylabel(r"$\sigma_{xy}$  ($e^2/h$)")
    top.set_title(f"Hall conductivity: the gap plateau gives C = {C:+.3f}")
    top.legend(frameon=False)

    Exx, sxx = c["sigma_xx"]
    bot.plot(Exx, sxx / norm, lw=1.2, color="C3")
    bot.set_xlabel("energy  (eV)")
    bot.set_ylabel(r"$\sigma_{xx}$  ($e^2/h$)")
    bot.set_title("longitudinal conductivity: zero in the gap, finite in the bands")

    fig.tight_layout()
    fig.savefig("fig_haldane_hall.png", dpi=150)
    print(f"Chern from gap plateau: C = {C:+.3f}")
    print("wrote fig_haldane_hall.png")


def demo():
    opdir = os.environ.get("HALDANE_OPDIR", "operators")
    if not os.path.exists(os.path.join(opdir, "haldane.HAM.desc")):
        import make_haldane                                 # stage the golden operators + .desc
        make_haldane.build()
    figure(opdir, label="haldane", M=128)


if __name__ == "__main__":
    demo()
