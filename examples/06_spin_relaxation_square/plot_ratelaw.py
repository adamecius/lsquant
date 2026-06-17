#!/usr/bin/env python3
"""Elliott-Yafet rate law from the random-phase (stochastic-trace) spin relaxation:
1/T1 vs spin-flip disorder Dsf, with the golden-rule theory slope. Reads
ey_ratelaw.dat (written by rp_validate.py). PRL/prl_style.
"""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ))
import prl_style as prl
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

D, rate = np.loadtxt("ey_ratelaw.dat", unpack=True)
p, logc = np.polyfit(np.log(D), np.log(rate), 1)

prl.use("web", aspect=0.62)
fig, ax = plt.subplots()
ax.plot(D, rate, label=r"random phase, $N_r{=}20$, $N{=}5000$", **prl.trace(0, line=False, marker=True))
xx = np.linspace(D.min(), D.max(), 50)
ax.plot(xx, np.exp(logc) * xx**p, color=prl.COLORS[0], ls="-", lw=1.2,
        label=r"fit $1/T_1\propto\Delta_{\rm sf}^{%.2f}$" % p)
# Elliott-Yafet golden-rule reference: slope 2, anchored at the smallest-Dsf point
c2 = rate[0] / D[0]**2
ax.plot(xx, c2 * xx**2, color="0.45", ls="--", lw=1.1,
        label=r"Elliott-Yafet $\propto\Delta_{\rm sf}^{2}$")
ax.set_xscale("log"); ax.set_yscale("log")
ax.set_xlabel(r"spin-flip disorder $\Delta_{\rm sf}\,[t]$")
ax.set_ylabel(r"relaxation rate $1/T_1\,[{\rm fs}^{-1}]$")
ax.legend(loc="upper left")
prl.save(fig, "fig_ey_ratelaw", formats=("png", "pdf"))
print("wrote fig_ey_ratelaw.png/pdf  (fitted p = %.2f; theory p = 2)" % p)
