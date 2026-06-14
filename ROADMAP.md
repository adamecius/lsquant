# lsquant — consolidation roadmap

Trunk `main` tracks the most advanced line
(SantiagoGimenezDC/lsquant-official @ main). Every branch that ever existed is
preserved as a git tag `archive/<name>` and can be restored with
`git checkout archive/<name>` — nothing was lost in consolidation.

## Status of the feature evaluation

Each unique commit on every kept branch was classified against `main`:

- **A — Superseded:** `main` already implements equivalent or better code.
- **B — Out of scope:** optimisation/infrastructure (memory, GPU, parallelism,
  build/toolchain) — the owner will handle these later.
- **C — Functionality `main` lacks:** evaluated for merge.

**Outcome:** `main` already contains essentially all of the KPM physics carried
by the old branches — optical/AC conductivity, local DOS / local densities,
mean-square displacement, torque, disorder, current operators, spin-projection
and projected time-evolution operators, and the non-orthogonal solver (further
developed than any branch). No feature met the bar to merge, so nothing was
merged into `main`; the genuinely-open items below are kept for review.

> **Verification note.** This evaluation could not build/run the library in its
> environment (no Intel MKL and no GSL available; could not be installed). Per
> the operating rule *"never merge code you cannot verify; doubt → roadmap"*,
> any feature that would require a build-and-smoke-test to confirm is listed
> here rather than merged. Re-run the merge evaluation in an environment with
> the GCC+Eigen+FFTW+GSL toolchain (`toolchain-gcc.cmake`) to close them.

## Dropped during this pass (recoverable via `archive/<name>`)

- `kpm_cuda` (GPU/CUDA) and `TB-QuantumTransp` (domain decomposition) — branches
  **deleted**; out of scope (optimisation). Tags `archive/kpm_cuda`,
  `archive/TB-QuantumTransp` retain the code if GPU/parallel work resumes.
- `master` — branch **deleted**; its only source work (2× "Refined nonOrth",
  2024-04) is **superseded** — `main` has the non-orthogonal solver in many more
  commits (`Non-orthogonal solver tested`, `ICHOL preconditioning`,
  `noOrth precision test`, `Flexible bandwidth shift`, …). Other commits were
  `.gitignore` only. Tag `archive/master`.
- `stable` — branch **deleted**; 2019 legacy snapshot, fully superseded. Tag
  `archive/stable`.
- `new_currents` — branch **deleted**; a strict subset of `alpha_release`
  (0 unique commits). Tag `archive/new_currents`.

## Branches kept — open items to review against `main` and reference

### `New_functions` — `7fa8ba5` (2025-02-20)

The 2025 observables work, marked "to be checked" by its author.

| Feature | Class | Evidence |
|---|---|---|
| Optical / AC conductivity | A | `main`: `opticalConductivityFromChebmom.cpp`, `chebyshev_momentsAC.cpp`, `optical_conductivity_tools.hpp` |
| Local DOS / local densities | A | `main`: `inline_compute-kpm-LocalDOS.cpp`, `localDensitiesFromChebmom.cpp`, `chebyshev_momentsLocal.cpp` |
| Mean-square displacement | A | `main`: `MSD_Evolve`, `MeanSquareDisplacement_Graphene_NNN`, `inline_compute-kpm-MeanSquareDisplacement_DeviceSim.cpp` |
| Time-evolution inline variants (`TimeEvOp`, `TimeEvProjetedOp`, `GetDeltaPhi`, `spacialresolution`) | A? | `main` has projected time-evolution (`inline_compute-kpm-TimeEvProjetedOp-Projecting_WF.cpp`); confirm these specific entry points are covered |
| **Quantum metric** | **C** | absent from `main` — **but no identifiable `metric`/geometric-tensor source exists on this branch either** (named only in the commit message) |

- [ ] **Quantum metric:** confirm whether a real implementation exists (here, or
      in the author's reference). If it exists, port to `main`'s interfaces,
      verify (positive semi-definite output on a known model), then merge. If it
      does not, it is a *new feature to implement*, not a merge.
- [ ] Confirm the time-evolution inline entry points are functionally covered by
      `main`; if any is genuinely new and verifiable, merge it; else drop.

### `python_supported` — `d49bc58` (2022-03-27)

| Feature | Class | Evidence |
|---|---|---|
| Spin-projection operator | A | `main`: `spin_project(...)` in `Device/Graphene_NNN.cpp`, `Device/Graphene_KaneMele.cpp` |
| Projected time-evolution operator | A | `main`: `inline_compute-kpm-TimeEvProjetedOp-Projecting_WF.cpp`, `time_evolution` namespace |
| Python interface / `wannier2sparse` bindings / `band_processing.py` / plotting | C (interface) | 2020–22 layer; will not build against current `main` as-is |

- [ ] **Restore Python interface:** decide whether to revive the pybind11
      bindings and tooling against current `main`. This is an *interface* task,
      not a numerics merge — port and confirm it builds before adopting; do not
      merge broken bindings.

### `alpha_release` — `2bc3df6` (2020-02-12)  ·  `beta_release` — `b2eeb8e` (2020-03-10)

2020 feature branches. Headline physics is already in `main`:

| Feature | Class | Evidence |
|---|---|---|
| Torque support (`beta_release`) | A | `main`: `createTorqueMatrix`, `createHoppingTorqueDensity_list` in `utilities/wannier2sparse/src/tbmodel.cpp` |
| Disorder (`beta_release`) | A | `main`: Anderson disorder + electron-hole puddles in `Graphene_NNN` |
| New current operators (`alpha_release`) | A | `main` has the full Kubo–Bastin/Greenwood current machinery (`currentTm`, `resetCurrentTm`, FFT solvers) |
| Electric-field input option (`alpha_release`) | C? | not clearly present in `main` as a discrete input flag — **only open item** |
| Fat-node / large-memory / high-memory mode | B | optimisation — out of scope |
| Parallel Kubo-Bastin, cmake/README fixes | B | infra — out of scope |

- [ ] Confirm the **electric-field input option** is covered by `main`'s
      conductivity machinery. If genuinely missing and worth having, port +
      verify + merge; otherwise these two branches are fully superseded and may
      be deleted (`archive/alpha_release`, `archive/beta_release` retain them).
