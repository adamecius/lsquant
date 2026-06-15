# Tutorial 2 - Localization on a 1D Anderson chain

This tutorial adds onsite disorder to the chain of Tutorial 1 and measures the
localization length from how a wavepacket spreads. It shows two things: the
localization length varies across the band the way weak-disorder theory
predicts, and a finite sample reaches that theory once its length passes a few
localization lengths.

This tutorial assumes you are inside the repository, in this example folder,
with LinQT built and installed (executables on `PATH`, or `LSQUANT_BIN` set to
the build directory). It reuses Tutorial 1; read that one first. The Python
helpers need `numpy` and `matplotlib` (`pip install -r requirements.txt`).

## The physics

The Anderson chain keeps the nearest-neighbour hopping and adds an independent
random onsite energy on each site:

```
H = t * sum_i ( |i><i+1| + |i+1><i| ) + sum_i eps_i |i><i|,
eps_i uniform in [-W/2, +W/2],   t = -1 eV.
```

For any nonzero `W` every eigenstate is exponentially localized, so 1D stays
insulating at all energies; the one-parameter scaling flow runs to the localized
fixed point for all conductances [4]. The localization length is the inverse
Lyapunov exponent, `xi(E) = 1 / gamma(E)`. In the weak-disorder (Born) regime,
away from the band centre and edges, box disorder gives [3]

```
xi(E) ~ 24 (4 - E^2) / W^2.
```

Two parts of the band carry their own behaviour. At the centre `E = 0` the
Kappus-Wegner anomaly raises the value to `xi(0) ~ 105 / W^2`, above the Born
estimate `96 / W^2` [3]. Near the edges `E -> ±2` the Lyapunov exponent follows
a separate `gamma ~ W^(2/3)` law, so the Born curve gives way there [3].

LinQT measures this through the mean-square displacement (MSD) of a wavepacket
[1, 2]. The clean chain spreads ballistically, `MSD ~ t^2`. Disorder localizes
the motion, so `MSD(E, t)` rises and then settles on a plateau set by `xi(E)`:
the relative size `sqrt(MSD_plateau)` tracks the localization length. The
moments of the spreading come from the same Chebyshev recursion as the DOS [1].

## Step 1 - build a disordered chain

The generator writes the disordered Hamiltonian, the velocity operator, the
widened bounds, and the sidecar. Onsite disorder sits on the diagonal, so the
velocity operator stays the clean one and the bounds widen to
`±(2 + W/2)` (Gershgorin).

```bash
python make_disordered_chain.py 512 1.0 1
```

This writes, for label `chain1d_dis_N512_W1`:

```
operators/chain1d_dis_N512_W1.HAM.CSR    Hamiltonian with box disorder (seed 1)
operators/chain1d_dis_N512_W1.VX.CSR     velocity operator v = i[H,X]
operators/chain1d_dis_N512_W1.HAM.desc   sidecar (records W, seed, bounds)
BOUNDS                                    (a b) = (-2.5 2.5) for W = 1
```

The sidecar records the disorder strength and the seed, so a result stays
traceable to the exact sample that produced it.

## Step 2 - the spreading moments and MSD

Compute the time-dependent moments with the MSD driver, then reconstruct
`MSD(E, t)` at a Fermi energy. Both steps are the heavy ones, so they run in
C++.

```bash
inline_compute-kpm-MeanSquareDisplacement chain1d_dis_N512_W1 VX 256 128 600
inline_timeCorrelationsFromChebmom Correlation*chain1d_dis_N512_W1*chebmomTD 10 0.0
```

The driver reads `operators/chain1d_dis_N512_W1.{HAM,VX}.CSR` and `BOUNDS`,
runs the recursion with 256 moments over 128 time steps up to `t = 600`, and
writes a time-dependent moment file. The reconstruction returns
`mean...EF0.000000...JACKSON.dat` with two columns, time and `MSD(E, t)`. The
window reaches `t = 600` because at `W = 1` the interior localization length is
`xi ~ 70-100 sites`: the wavepacket spreads ballistically up to `xi` and then
settles, and that settling takes a long time, so a short window (e.g. `t = 40`)
stops while the MSD is still rising, far below its plateau. A random starting
vector seeds the trace, so setting `KPM_SEED` makes each run reproducible; one
realisation per disorder sample keeps the average honest.

## Step 3 - localization length across the band

A single realisation is noisy, so the study averages `MSD(E, t)` over disorder
samples, reads `sqrt(MSD_plateau)` at each energy, and fixes one overall
constant against the clean ballistic reference. The wrapper does this and draws
the figure:

```bash
python lsqloc.py
```

This writes `fig_localization.png` with two panels. It runs three sizes at four
energies, `E = 0` (centre), `E = ±1` (intermediate), and `E = 1.9` (near the
edge), averaged over disorder realisations at `W = 1 eV`.

**Top panel - `xi(E)` across the band, three sizes, with the Born curve.** The
largest size follows `xi(E) ~ 24 (4 - E^2) / W^2` through the band interior: the
`E = ±1` points land on the curve (and agree with each other, confirming the
symmetry of the clean dispersion), and the centre `E = 0` approaches the Born
value as the size and the averaging grow. The near-edge point `E = 1.9` departs
from the Born curve, where the `W^(2/3)` edge law takes over [3]. The centre
carries the longest `xi`, so a small sample falls below the curve there first,
and the Born interior form is recovered with the largest size.

The plateau gives a localization-length proxy `sqrt(MSD_plateau)`, fixed to an
absolute scale by one global constant (set at `E = 1`). That proxy reproduces
the Born `(4 - E^2)` interior trend and the near-edge departure; the small
Kappus-Wegner enhancement at the very centre (`xi(0) ~ 105/W^2`, a few percent
above the Born `96/W^2` [3]) sits within the disorder noise of this study, so
the centre point lands on the Born curve rather than visibly above it.

**Bottom panel - scaling with size, at each marked energy.** Each curve rises
with `N` while the sample is shorter than `xi`, then settles on the
disorder-limited value once `N` passes a few `xi`. The near-edge energy has the
shortest `xi`, so it plateaus at the smallest `N`. The centre has the longest
`xi`, so it needs the largest `N`. Reaching the theory asks for two choices
together: `N` beyond the plateau onset (the sample holds the localized state),
and `W` strong enough to keep `xi` inside reachable `N` while staying in the
weak-disorder regime where the Born formula applies. At `W = 1 eV` the
centre value `xi(0) ~ 100 sites` sets the largest `N` the study needs.

## What to take away

- Disorder localizes every state in 1D; the chain stays insulating across the
  whole band [4].
- The localization length follows `xi(E) ~ 24 (4 - E^2) / W^2` in the band
  interior and changes law near the edges [3]; the centre carries a small
  Kappus-Wegner enhancement that needs heavy averaging to resolve.
- The MSD plateau gives a localization-length proxy, and its `1/W^2` and
  `(4 - E^2)` trends read straight from the spreading (halving `W` roughly
  quadruples `xi`).
- A finite sample reproduces the theory once its length passes a few `xi`, which
  fixes how large `N` and how strong `W` the study needs.

## References

[1] A. Weiße, G. Wellein, A. Alvermann, H. Fehske,
*The kernel polynomial method*, Rev. Mod. Phys. **78**, 275 (2006).
https://journals.aps.org/rmp/abstract/10.1103/RevModPhys.78.275

[2] J. H. García et al., *LinQT / linear-scaling quantum transport*,
arXiv:1811.07387. https://arxiv.org/abs/1811.07387

[3] Weak-disorder localization length and band anomalies (1D Anderson model),
arXiv:cond-mat/9801226. https://arxiv.org/abs/cond-mat/9801226

[4] E. Abrahams, P. W. Anderson, D. C. Licciardello, T. V. Ramakrishnan,
*Scaling theory of localization*, Phys. Rev. Lett. **42**, 673 (1979).
