# lsquant — consolidation roadmap

Trunk `main` now tracks the most advanced line
(SantiagoGimenezDC/lsquant-official @ main).
Branches deleted during consolidation are preserved as git tags
`archive/<name>` and can be restored with `git checkout archive/<name>`.

## Highlighted directions
- **GPU / CUDA KPM backend** — branch `kpm_cuda` (PhD-era GPU implementation).
- **Parallelism / domain decomposition** — branch `TB-QuantumTransp`.
- **Python interface & tooling** — branch `python_supported`.
- **New observables (to be verified)** — branch `New_functions`:
  optical conductivity, local densities, mean-square displacement, quantum metric.

## Branches kept for review against `main` and reference

### `master` — e8ff89e (2024-05-01), 4 commit(s) not in `main`

- [ ] Review against `main`
- [ ] Validate against reference implementation/results

<details><summary>unique commits</summary>

- \`928d17d\` 2024-04-30 — Refined nonOrth
- \`0bbf919\` 2024-04-30 — Refined nonOrth
- \`2e3e478\` 2024-05-01 — Ignoring build files
- \`e8ff89e\` 2024-05-01 — Ignoring build files

</details>

### `New_functions` — 7fa8ba5 (2025-02-20), 6 commit(s) not in `main`

- [ ] Review against `main`
- [ ] Validate against reference implementation/results

<details><summary>unique commits</summary>

- \`928d17d\` 2024-04-30 — Refined nonOrth
- \`0bbf919\` 2024-04-30 — Refined nonOrth
- \`2e3e478\` 2024-05-01 — Ignoring build files
- \`e8ff89e\` 2024-05-01 — Ignoring build files
- \`1b1f9f3\` 2025-02-20 — New functions implemented (to be checked) in optical conductivity, local densities, mean square displacement, quantum metric...
- \`7fa8ba5\` 2025-02-20 — The version with new functions clean

</details>

### `python_supported` — d49bc58 (2022-03-27), 11 commit(s) not in `main`

- [ ] Review against `main`
- [ ] Validate against reference implementation/results

<details><summary>unique commits</summary>

- \`14f461a\` 2020-07-27 — In this branch we included some python supported modules
- \`2b6e855\` 2020-07-27 — fix the python2.X problem and add example of wannier2sparse module
- \`e033a0b\` 2020-08-26 — Added onsite_spinor function in wannier2sparse
- \`465883f\` 2020-09-17 — Added the spin projection operator to the operators to be obtained from the wannier system
- \`09da47e\` 2020-09-18 — Merge pull request #3 from AngelDPS/python_supported
- \`f9d571a\` 2020-09-26 — Implemented a program for time evolution projected operator exp val
- \`e18b184\` 2020-10-05 — Rearanged the order of the operators for the expression <Psi|PRJ U OP delta U PRJ|Psi>
- \`182db75\` 2020-10-28 — Merge pull request #4 from AngelDPS/python_supported
- \`bfad86d\` 2020-11-24 — improve the plotting scripts
- \`eb96997\` 2021-01-30 — added band_processing.py
- \`d49bc58\` 2022-03-27 — Hotfix for the spin operators, and update on changed files

</details>

### `kpm_cuda` — d590ed7 (2019-10-23), 1 commit(s) not in `main`

- [ ] Review against `main`
- [ ] Validate against reference implementation/results

<details><summary>unique commits</summary>

- \`d590ed7\` 2019-10-23 — This is the sofware I used for my PhD

</details>

### `TB-QuantumTransp` — 6a6c976 (2019-10-23), 1 commit(s) not in `main`

- [ ] Review against `main`
- [ ] Validate against reference implementation/results

<details><summary>unique commits</summary>

- \`6a6c976\` 2019-10-23 — This is a old code where I try to implement domain decomposition

</details>

### `alpha_release` — 2bc3df6 (2020-02-12), 6 commit(s) not in `main`

- [ ] Review against `main`
- [ ] Validate against reference implementation/results

<details><summary>unique commits</summary>

- \`1abfeed\` 2020-01-18 — include new current operators, include a testing vector os kubobastinfromchebmom parallel
- \`d5f214e\` 2020-01-18 — add option for changn the electric field
- \`5678cb8\` 2020-01-20 — add support for high-memory calculations
- \`42baee3\` 2020-01-29 — Fix cmake file
- \`5f48a88\` 2020-02-12 — Create README.md
- \`2bc3df6\` 2020-02-12 — Update README.md

</details>

### `beta_release` — b2eeb8e (2020-03-10), 5 commit(s) not in `main`

- [ ] Review against `main`
- [ ] Validate against reference implementation/results

<details><summary>unique commits</summary>

- \`fe5b00f\` 2020-02-20 — include support of torque
- \`8e455e5\` 2020-02-21 — Fix large memory problems
- \`eac8ef6\` 2020-02-27 — added full fatnode support
- \`708f42a\` 2020-02-28 — added disorder support
- \`b2eeb8e\` 2020-03-10 — Fix compilation issues

</details>

### `new_currents` — 5678cb8 (2020-01-20), 3 commit(s) not in `main`

- [ ] Review against `main`
- [ ] Validate against reference implementation/results

<details><summary>unique commits</summary>

- \`1abfeed\` 2020-01-18 — include new current operators, include a testing vector os kubobastinfromchebmom parallel
- \`d5f214e\` 2020-01-18 — add option for changn the electric field
- \`5678cb8\` 2020-01-20 — add support for high-memory calculations

</details>

### `stable` — 0210d1a (2019-10-23), 1 commit(s) not in `main`

- [ ] Review against `main`
- [ ] Validate against reference implementation/results

<details><summary>unique commits</summary>

- \`0210d1a\` 2019-10-23 — The stable version of my code 0.0.0

</details>

