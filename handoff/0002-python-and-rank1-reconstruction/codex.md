# Handoff 0002 (tester) — stage 2 gate: G3 byte for byte, tutorial 01 from Python

Read `AGENTS.md` sections 8, 9 and 12 and `brief.md` of this handoff. Everything on TCN,
in `~/lsquant-unified`, branch `handoff/0002-codex` cut from `cpp-unified` after merging
`origin/handoff/0002-agi`. You never edit `lsquant/`, `python/lsquant/` (except nothing),
the legacy tree, a golden, a seed or a tolerance.

## 1. Build

`pip install --user pybind11 sympy` if the developer's report says so (record versions);
fresh build with `-DLSQUANT_BUILD_NEW_CORE=ON -DLSQUANT_BUILD_PYTHON=ON` and the
`-Dpybind11_DIR=` line from the report; `PYTHONPATH=build/lsquant_core:python` for every
Python run. Run `python tests/dev/stage_2_check.py` and record its output.

## 2. The gate (`scripts/check_goldens.sh 2`)

Extend the harness with stage 2:

1. `python -m lsquant.apps.stage2_dos --moments goldens/stage0/chain1d/SpectralOp1chain1dKPM_M64_stateexact.chebmom1D --broadening-meV 120 --out $TMP/a`
   then `cmp $TMP/a/meanSpectralOp1chain1dKPM_M64_stateexactJACKSON.dat goldens/stage0/chain1d/meanSpectralOp1chain1dKPM_M64_stateexactJACKSON.dat` (G3).
2. `python -m lsquant.apps.stage2_dos --operators goldens/stage0/chain1d/operators --label chain1d --window -2 2 --trace exact --states 2048 --moments 64 --broadening-meV 120 --out $TMP/b`
   then `cmp` of both the `.chebmom1D` (against G1) and the `JACKSON.dat` (against G3).
3. Both at `OMP_NUM_THREADS=1` and `16`.
4. Stage 1 must still pass: run `check_goldens.sh 1` as part of stage 2.
5. `grep -rn "delta\|DOS\|Kubo\|basis\|std::cout\|getenv" lsquant/core lsquant/io` on
   **code lines only** (exclude lines whose first non-blank characters are `//`, `*`
   or `/*`) must be empty; the mathematical notation `delta_{m0}` in stage-1 comments
   is not a finding.

A `DIFFERS` on the `.dat` is localised with `scripts/golden_diff.py`: first differing
line, both values, and whether the difference is in the 6th significant digit only (a
last-bit difference in the table or the sum order) or larger (grid, kernel cutoff,
prefactor, axis). Classify by symptom: all lines differ by a scale factor (prefactor
`N*scale`), a shift (axis `x/scale + shift`), a different line count (`num_div`,
`M_eff`), differences growing toward the band edge (`alpha`), isolated last-digit
flips (summation order or libm). Do not fix.

## 3. The figure (`figures/stage_2.py`)

Tutorial 01 reproduced from Python: run `python -m lsquant.apps.tutorial01 --out <tmp>`;
it writes the moment files with their legacy names at the top of `<tmp>` and the curves
under `<tmp>/eta<meV>/`. Match them to `goldens/stage0/tutorials/01/captured/` by the
tutorial's order of operations (the captured prefixes follow `lsqplot.py`): the moment
files by name after the prefix; the curves at 20 meV for `N = 256, 1024, 4096` are
`0013_`, `0015_`, `0017_`; the `N = 256` curves at 80, 24, 10 meV are `0019_`, `0020_`,
`0021_`. `cmp` each pair and draw: (a) `rho(E)` for `N = 256, 1024, 4096` at 20 meV,
new as lines, golden as markers, the analytic `1/(pi sqrt(4 - E^2))` as a thin
reference, residual new minus golden below; (b) `N = 256` at 80, 24, 10 meV, same
encoding. Caption `FIG. 2.` with sizes, `M = 512`, broadenings, the `cmp` results.
Save `figures/out/stage_2/fig_stage_2.png` and `.caption.txt`.

## 4. Report

`handoff/0002-python-and-rank1-reconstruction/tests.md`: versions, build, self-check
output, gate output verbatim (stages 1 and 2), the tutorial-01 `cmp` list, the figure
path and caption, smells (numpy used where `math` was required, a physics name in C++,
a table computed in C++, a copy in the bindings, a `print` in the package). Commit on
`handoff/0002-codex`, push, report to Fable.
