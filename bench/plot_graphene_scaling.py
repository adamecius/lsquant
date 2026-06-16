#!/usr/bin/env python3
"""Phase-F figures: threaded-SpMV scaling of lsquant on graphene. PRL/prl_style.
Reads graphene_threads.csv + graphene_size.csv (from run_graphene_scaling.sh),
writes fig_scaling_threads, fig_scaling_size, fig_roofline (png+pdf) + captions.
"""
import csv, os, sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "examples"))
import prl_style as prl
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(__file__)
ROOF_SOCKET, ROOF_DUAL = 348.2, 700.8

def load(path):
    rows = []
    with open(os.path.join(HERE, path)) as f:
        for r in csv.DictReader(f):
            rows.append(r)
    return rows

def fnum(x):
    try: return float(x)
    except: return None

T = load("graphene_threads.csv")
S = load("graphene_size.csv")
prl.use("web")
USETEX = matplotlib.rcParams["text.usetex"]
N_fixed = T[0]["N"] if T else "?"

# ---- FIG 1: thread scaling (SpMV kernel vs full DOS recursion), log-log + ideal ----
thr = [int(r["threads"]) for r in T]
spmv = [fnum(r["spmv_ms"]) for r in T]
recur = [fnum(r["recursion_ms"]) for r in T]
fig, ax = plt.subplots()
ax.plot(thr, spmv,  label=r"SpMV kernel (threaded)", **prl.trace(0, line=True, marker=True))
ax.plot(thr, recur, label=r"full DOS recursion",    **prl.trace(1, line=True, marker=True))
ax.plot(thr, [spmv[0]/t for t in thr], color="0.6", lw=0.8, ls=(0,(1,1)))
ax.text(0.04, 0.08, r"dotted grey: ideal $\propto 1/p$", transform=ax.transAxes, fontsize=8, color="0.4")
ax.set_xscale("log", base=2); ax.set_yscale("log")
ax.set_xlabel(r"threads $p$")
ax.set_ylabel(r"time for $M{=}512$ moments [ms]")
ax.legend(loc="upper right")
prl.save(fig, os.path.join(HERE, "fig_scaling_threads"), formats=("png","pdf")); plt.close(fig)

# ---- FIG 2: size scaling, spmv total time vs N (log-log) + linear-in-N ref ----
Ns = [int(r["N"]) for r in S]
sp_ms = [fnum(r["spmv_ms"]) for r in S]
fig, ax = plt.subplots()
ax.plot(Ns, sp_ms, label=r"SpMV total ($M{=}256$, 128 threads)", **prl.trace(0, line=True, marker=True))
# linear-in-N reference anchored at the largest point
k = sp_ms[-1]/Ns[-1]
ax.plot(Ns, [k*n for n in Ns], color="0.6", lw=0.8, ls=(0,(1,1)))
ax.text(0.04, 0.86, r"dotted grey: $\propto N$ (linear KPM)", transform=ax.transAxes, fontsize=8, color="0.4")
ax.set_xscale("log"); ax.set_yscale("log")
ax.set_xlabel(r"system size $N$ (sites)")
ax.set_ylabel(r"SpMV time, $M{=}256$ [ms]")
ax.legend(loc="upper left")
prl.save(fig, os.path.join(HERE, "fig_scaling_size"), formats=("png","pdf")); plt.close(fig)

# ---- FIG 3: roofline -- SpMV GB/s vs threads + measured STREAM read ceilings ----
g = [fnum(r["gbs"]) for r in T]
fig, ax = plt.subplots()
ax.plot(thr, g, label=r"graphene SpMV (matrix-stream)", **prl.trace(0, line=True, marker=True))
ax.axhline(ROOF_SOCKET, color="#000000", lw=1.0, ls="--")
ax.text(1.1, ROOF_SOCKET*1.02, "1-socket read roofline 348 GB/s", fontsize=7.5)
ax.axhline(ROOF_DUAL, color="0.45", lw=1.0, ls=":")
ax.text(1.1, ROOF_DUAL*1.02, "2-socket read roofline 700 GB/s", fontsize=7.5, color="0.3")
ax.set_xscale("log", base=2)
ax.set_xlabel(r"threads $p$")
ax.set_ylabel(r"SpMV matrix-stream [GB/s]")
ax.legend(loc="upper left")
prl.save(fig, os.path.join(HERE, "fig_roofline"), formats=("png","pdf")); plt.close(fig)

# ---- captions ----
peak_g = max([x for x in g if x is not None], default=0)
peak_pct = 100*peak_g/ROOF_DUAL
with open(os.path.join(HERE, "CAPTIONS.md"), "w") as f:
    f.write(f"""# Phase-F figure captions (PRL/prl_style; usetex={USETEX})

FIG. 1. Strong-thread scaling of the threaded lsquant Chebyshev DOS recursion on a
graphene (complex Hermitian) Hamiltonian of N={N_fixed} sites (honeycomb, 2-atom basis,
nearest-neighbour t1=-2.7 eV plus complex Haldane next-nearest hopping t2=0.10 eV,
phi=pi/2; ~9 nonzeros/row), M=512 Chebyshev moments, vs thread count p (log-log).
Solid blue circles: the SpMV kernel region ("spmv" timer); dashed orange squares: the
full DOS recursion ("compute_spectral/SpectralMoments"); dotted grey: ideal 1/p. The
SpMV scales toward the bandwidth knee; the full recursion saturates earlier because the
moment dot-products / AXPY are deliberately left serial (golden-safe). Pinned
OMP_PROC_BIND=close OMP_PLACES=cores, numactl per socket. dual-EPYC-9754.

FIG. 2. Size scaling: SpMV total time for M=256 moments at 128 threads vs system size N
(log-log); dotted grey is a proportional-to-N reference. The KPM cost is linear in N as
expected. Same graphene model and pinning as FIG. 1.

FIG. 3. SpMV matrix-stream bandwidth (GB/s) vs thread count for the FIG.-1 graphene
system; dashed black and dotted grey are the measured STREAM read-only roofline for one
socket (348 GB/s) and two sockets (700 GB/s). Peak achieved {peak_g:.0f} GB/s
(~{peak_pct:.0f}% of the 2-socket roofline). Bytes counted: values(16 B)+index(4 B) per
nonzero + outer(8 B) per row, times the number of SpMV calls.
""")
print("wrote figs + CAPTIONS.md; peak GB/s =", round(peak_g,1), "(", round(peak_pct,1), "% of dual roofline)")
