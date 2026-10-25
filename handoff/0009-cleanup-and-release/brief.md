# Handoff 0009 (Agi) — clean-up, build, README, examples, and the release pull request

Roles: Agi implements (maximum reasoning effort), the tester (`codex.md`) runs the gates and
regenerates every example figure on TCN, Fable reviews and opens the pull request. If Agi runs out
of credits or quota, Codex (reasoning effort `medium`) takes over this brief as developer from Agi's
last commit, with the same rules; the tester role then goes to an Opus subagent driven by Fable.

Read `AGENTS.md` (sections 3, 8, 10, 11 and 12 are binding) and this brief completely. Nothing is
optional; deviations go in `report.md`.

Starting point: branch `handoff/0009-agi` from `origin/handoff/0009-agi`. It is
`handoff/0006b-dev` (`2c790f1`, stages 1 to 7 plus the SIESTA spin Hall work) plus this brief and
the SIESTA inputs of the two new examples. `origin/handoff/0008-dev` (`10d02e4`, the stage-8
performance work) is merged in step 1.

## 0. Scope

The code works and its goldens pass. This handoff makes it a release:

1. one branch that holds stages 1 to 8 (step 1);
2. no file, function, option or comment left over from the legacy code or from the stage
   scaffolding that a user does not need (steps 2 and 3);
3. a CMake build a user can read and drive from options or presets (step 4);
4. a README that says what the code is now, written like the old one (step 5);
5. the six tutorials driven by the new code, plus two new SIESTA examples: the Bi(111) bilayer
   and the PdSe2 monolayer, both with spin-orbit coupling (step 6).

Out of scope: new physics, new spectral functions, any change of a numerical result, the
performance items still open from 0008 (`MomentTable` by rows, rank 2 with stored right vectors,
k/-k pairing), and the numerical difference measured on Aron (Bi 1.135 against 1.155, PdSe2 2.52
against 2.34). Fable diagnoses that difference separately; do not touch the SIESTA, oracle or
velocity numerics while it is open. If Fable's diagnosis needs a code change, it arrives as an
amendment to this brief.

## 1. Restrictions

- Every golden stays byte-identical: `scripts/check_goldens.sh 7` must give the same numeric
  result it gives today (70 numeric files identical) and pass the text checks.
- Do not edit `goldens/`, `scripts/check_goldens.sh`, `figures/` or `runs/goldens/`: they belong
  to the tester. If a cleanup would require a change there (a path moves, a file is renamed),
  do not make it; list it in `report.md` under "needs harness change" and Fable decides.
- `lsquant/compat/` stays as it is (AGENTS.md section 3).
- The C++ and Python entry points the gate calls keep their names in this handoff:
  `stage1_moments`, `stage6_nonorth`, `stage_1_check`, `python -m lsquant.apps.stage*`,
  `python -m lsquant.apps.tutorial0*`. The README calls them what they are (golden reproduction
  drivers). Renaming them is a follow-up with a harness change.
- No new tests. No new dependencies beyond what `pyproject.toml` and the CMake already require,
  except `sisl` as an optional extra for the SIESTA examples (it is already imported by
  `python/lsquant/models/siesta.py`).
- Everything in English. Commits small and by topic (merge, deletions, code fixes, CMake, README,
  each example). Never push.

## 2. Step 1 — merge stage 8

`git merge origin/handoff/0008-dev`. Three files conflict (checked by Fable on a trial merge), one
hunk each:

- `lsquant/core/chebyshev/moment_engine.cpp`: 0006b parallelised the rank-2 rows with a static
  schedule (`d1b4202`); 0008 added a `LSQUANT_NATIVE_BUILD` path through `linalg::MomentTable`.
  Keep both: the golden profile runs 0006b's row-owned loop unchanged, the native profile runs the
  `MomentTable` block. Both branches must accumulate into `mu` in the same order as before.
- `lsquant/core/algebra/block_diagonal_operator.cpp` and `python/lsquant/runfile.py`: union of
  both sides; neither side removes what the other added.

Check: the tree builds with `LSQUANT_BUILD_PROFILE=golden` and `=native`, zero warnings at
`-Wall -Wextra`; `tests/dev/stage_1_check` passes. Record both build logs in `report.md`.

## 3. Step 2 — remove dead files

For every file below: run `git grep -n <basename>` outside `handoff/` and `goldens/`, paste the
result in `report.md`, then delete, move or keep as the table says. A file that the gate or a
figure script references (`scripts/check_goldens.sh`, `figures/*.py`) is never deleted; it goes to
"needs harness change".

| Path | Action | Reason |
|---|---|---|
| `TODO.md`, `infrastructure_plan.md`, `reporting_system.md`, `unified_io.md` | delete | plans for the legacy reporting system; superseded by AGENTS.md and the run file |
| `GOLDEN_PROVENANCE.txt` | delete after copying anything not already in `goldens/stage0/PROVENANCE.md` into `report.md` (Fable moves it) | duplicate provenance |
| `install.py` | delete; fix the reference in `docs/main.md` | legacy installer |
| `lsquant_derive_core.py`, `lsquant_derive_spin.py`, `lsquant_derive_transport.py` | move to `docs/derivations/` with a one-line README each | sympy derivations kept as documentation, unreferenced from code |
| `lsquant_chain_reference.py` | move to `python/lsquant/oracle/chain.py`, keep a CLI `python -m lsquant.oracle.chain` | G2 oracle (AGENTS.md section 8); update the AGENTS.md mention in `report.md`, Fable edits AGENTS.md |
| `utilities/python/` (all: `bin/*chebmom2D.py`, `lib/sciplots.py`, `lsquant_report.py`, `test/`, `.ipynb_checkpoints/`, every `__pycache__/*.pyc`) | delete | legacy `.chebmom2D` tools, Python 3.6 bytecode, notebook checkpoints |
| `utilities/postprocessing/average_moments1D.cpp`, `utilities/CMakeLists.txt` | delete | legacy format, never built by the top-level CMake |
| `utilities/siesta2sparse/` | delete if `python/lsquant/models/siesta.py` covers both scripts (state which function replaces each); otherwise keep and say why | superseded by the SIESTA model loader |
| `examples/*/lsq*.py`, `examples/06_spin_relaxation_square/rp_*.py`, `plot_ratelaw.py` | replace in step 6 | they drive the legacy `lsquant compute` CLI that no longer exists |
| `examples/style_demo.py` | delete | demo of the plotting style |
| `examples/prl_style.py` | keep | `figures/stage_*.py` import it |
| `examples/utilities/kwant2wannier/` | keep only if `bench/make_graphene.py` still needs it; delete the notebook | |
| `bench/run_graphene_scaling.sh`, `bench/make_graphene.py`, `bench/plot_graphene_scaling.py`, `bench/CAPTIONS.md` | delete if they drive the legacy binaries; the 0008 `bench/*.py` replace them | |
| `tests/block_domains/` (untracked, Fable's local-domain benchmark) | move to `bench/block_domains/` without `BRIEFING_LSQUANT_SESSION.md`; Fable hands you the files | evidence for the README performance section |
| `docs/Doxyfile.in`, `docs/main.cpp`, `docs/main.md`, `docs/CMakeLists.txt`, `docs/requirements.txt` | update to the new tree (inputs `lsquant/` and `python/lsquant/`), or delete if Read the Docs no longer builds from them; say which | |
| `docs/*.md` | remove every mention of `.chebmom1D/2D/TD` as a user format and of `inline_*`, `lsquant compute`, `lsquant reconstruct <file> dos` | |
| `.gitignore`, `examples/.gitignore` | drop the `test/golden`, `test/models` and legacy `bench/gr` lines | the paths no longer exist |

Then search the whole tree for what the table does not cover: `git ls-files` for `*.pyc`,
`.ipynb_checkpoints`, `*.orig`, editor files, files larger than 1 MB outside `goldens/` and
`examples/*/siesta/`; any Python module under `python/lsquant/` that nothing imports (vulture or
`python -m pyflakes` plus a grep); any C++ function under `lsquant/` with no caller. List every
hit in `report.md` with the action taken.

The `inline_*` text check of gate 7 fails today on six files: `README.md`,
`agent_template/tutorial_style.md`, `agent_template/tutorial_template.md`,
`python/lsquant/moments.py`, `tests/dev/stage_4_check.py` and `AGENTS.md`. Remove the mention from
the first five (the legacy name may be replaced by the word it computes). `AGENTS.md` is the
project record; the tester adds it to the record exclusion of the harness.

## 4. Step 3 — code review pass

Over `lsquant/` and `python/lsquant/`:

1. **Names.** Every identifier follows AGENTS.md section 11. List every name that is an
   abbreviation not in the originals, a single letter outside a loop index or a formula symbol, a
   leftover of the legacy name (`chebmom1D` in a variable that now holds a rank-generic tensor),
   or a name that says the stage rather than the job (`stage6b_*` helpers inside library modules).
   Rename inside the library; keep the gate entry points (section 1).
2. **Comments.** Remove comments that refer to handoffs, reviews, testers, items or stages
   (`// Item 3:`, `# Review 2`, `Tester: compare ...`): the reason for the code stays, its history
   goes to git. Keep the comments that state the problem, the convention, the units, the
   assumption or the legacy defect reproduced on purpose (AGENTS.md section 10).
3. **Dead code.** Functions, branches, parameters and options with no caller; diagnostic code
   added for one review (`d3680c9`, `3f3bfa5` printed spreads and NumPy comparisons) unless a run
   file or the gate uses it.
4. **Warnings.** Zero warnings at `-Wall -Wextra` in both profiles; `python -m pyflakes
   python/lsquant` clean.
5. **Errors.** Every `throw`/`raise` message says what was expected and what was received.
   Hard-coded numerical constants that are not conventions (tolerances, block sizes, thread
   counts) become named constants or run-file fields; list each with its value.

Every change in this step is behaviour-preserving. Paste the before/after list in `report.md`.

## 5. Step 4 — CMake

One top-level `CMakeLists.txt` that a user reads in a minute:

- options: `LSQUANT_BUILD_PROFILE` (`golden` | `native`, default `golden`),
  `LSQUANT_BUILD_PYTHON` (ON), `LSQUANT_BUILD_BENCHMARKS` (OFF), `LSQUANT_BUILD_DEV_CHECKS` (ON,
  builds `stage_1_check`), `LSQUANT_ENABLE_FFT_CONTRACTION` if 0008 made it optional. Remove
  `LSQUANT_BUILD_NEW_CORE` (there is no old core) and `LSQUANT_BUILD_INPUT_TOOLS` unless
  `cmake/FetchWannier2Sparse.cmake` exists and works; `scripts/fetch_wannier2sparse.sh` is the
  supported route.
- dependencies found once, at the top, each with the target it serves: Eigen3, OpenMP and FFTW3
  for `lsquant_core`; GSL and pybind11 only when `LSQUANT_BUILD_PYTHON` (GSL is used only by
  `compat/`).
- `install()` rules: `lsquant_core` and its headers, the two apps, and `_core` into the Python
  package directory, so that `cmake --install` followed by `pip install ./python` gives a working
  `lsquant` command.
- `CMakePresets.json` with `golden` and `native` configure presets and matching build presets.
- `compile_commands.json` exported.
- No per-target repetition of `-Wall -Wextra`: one interface target or one function.

Check: fresh configure and build from each preset; `cmake --install` into a temporary prefix,
then `python -c "import lsquant, lsquant._core"` from outside the source tree. Paste the commands
and the output.

## 6. Step 5 — README

Rewrite `README.md` with the structure of the current one (what it computes, references,
performance, installation, interface, tutorials, how it works, inputs, verification,
documentation, repository layout, licence and citation) so a reader of the old README finds
everything where it was. The content describes the code as it is now:

- every observable is a word `Tr[A_0 f_1(H) A_1 ... f_p(H) A_p]`; the word language, the
  catalogue (`python -m lsquant words`), a short table of the words that replace the old drivers;
- the run file (one JSON: model, basis, metric, word, expansion, reconstruction, output), the
  CLI (`run`, `reconstruct`, `inspect`, `words`) and the MCP server, with one complete example;
- the C++ engine and its five instructions; reconstruction in Python with sympy coefficient
  tables;
- the three bases (real space, Bloch with real-space disorder through FFT, eigenbasis as the exact
  oracle) and the non-orthogonal metric (Lowdin, generalised, anticommutator); SIESTA models
  through `sisl`; Wannier models through wannier2sparse;
- what is new with respect to the previous release, as a list: unified engine, words, bases,
  non-orthogonal metric, SIESTA spin Hall, build profiles, FFT contraction, the corrected legacy
  defects that matter to users (AGENTS.md section 10: the batched rank-2 shift, the off-centre
  energy axis, the identity operator of tutorial 06);
- performance: only numbers that a committed file supports, each with its source
  (`bench/block_domains/RESULTS.md`, the 0008 report); no number from memory;
- verification: the golden masters, `scripts/check_goldens.sh`, the stage figures;
- installation from source with the presets; Python install; optional `sisl`.

`python/lsquant/README.md` becomes a short API page that points to the main README. No marketing
language; the register of the current README.

## 7. Step 6 — examples

Each example directory holds: `README.md` (physics, model, the word, the run file, how to run,
the figure and what it shows), the model generator (`make_*.py`, unchanged where it already
exists), one or more run files (`run_*.json`), and one `plot.py` that reads the run output and
draws the figure with `examples/prl_style.py`. The figure files (`fig_*.png`) are regenerated by
the tester on TCN and committed by the tester; Agi does not commit figures.

1. **Tutorials 01 to 06.** Same physics, same figures as before. The run files reproduce the
   tutorial captures already committed as goldens (`runs/goldens/tutorial0*`): reuse those run
   files' parameters; the example run files differ only in paths. The legacy drivers (`lsq*.py`,
   `rp_*.py`, `plot_ratelaw.py`) are deleted. Tutorial 06 uses the identity operator explicitly
   (AGENTS.md section 10).
2. **`examples/07_spin_hall_bi_bilayer/`.** Bi(111) bilayer with SOC. Inputs in `siesta/`
   (`bi.fdf`, `bi.TSHS`, `Bi.ion.xml`, `Bi.psml`; committed by Fable, `SHA256SUMS` alongside).
   Run files: DOS and spin Hall conductivity `Tr[JXSZ delta(E-H) VY dG-(E)/dE]` with the
   covariant velocity and the Lowdin metric, as in `runs/goldens/g11_bi_covariant_lowdin_shc.json`
   but with `model.dir` relative to the example. Two sizes: `quick` (24 x 24 k-points, `M = 256`,
   minutes on a laptop) and `production` (48 x 48, `M = 1024`, the G11 setting). The figure:
   `sigma^z_xy(E_F)` from KPM against the eigenbasis oracle, residual panel. The README states
   the production value at `E_F = 0` only after Fable's diagnosis of the Aron difference;
   until then it says "value pending" and the tester's number goes in `tests.md`.
3. **`examples/08_spin_hall_pdse2_monolayer/`.** Same structure for PdSe2 (`pdse2.fdf`,
   `pdse2.TSHS`, `Pd.ion.xml`, `Se.ion.xml`, `Pd.psml`, `Se.psml`), the parameters of
   `runs/goldens/g11_pdse2_covariant_lowdin_shc.json`, `quick` at 12 x 12.
4. The dipole cache that `models/siesta.py` writes to disk goes to the run output directory or a
   cache directory, never into `siesta/`; add it to `examples/.gitignore`.

## 8. Report

`handoff/0009-cleanup-and-release/report.md`: the merge resolution; per step the files changed,
the grep evidence and the check output verbatim; the "needs harness change" list; the list of
renamed identifiers and named constants; deviations; open questions. Agi does not declare the
item done. Fable reviews, the tester runs `codex.md`, and after a clean review Fable merges into
`cpp-unified` and opens the pull request `cpp-unified -> main`.
