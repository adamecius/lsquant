# Handoff 0007 (Agi) — stage 7: the word language, the run file, the CLI, the MCP tool, and the end of the legacy tree

Roles: Agi implements, the tester (`codex.md`) reproduces every golden from words on TCN,
Fable reviews. Read `AGENTS.md` (sections 2, 3, 8, 9, 11 are binding) and this brief
completely. Nothing is optional; deviations go in `report.md`. Starting point: branch
`handoff/0007-agi` from `origin/cpp-unified` **after the stage-6 merge** (Fable tells you the
commit); it contains stages 1 to 6.

## 0. Scope

Everything the code can compute is a word, `Tr[A_0 f_1(H) A_1 ... f_p(H) A_p]` (AGENTS.md
section 2). Until now every word was a Python function (`dos_plan`, `kubo_plan`, `msd_plan`,
`projected_evolution_plan`, `time_correlation_plan`) called from a stage app. Stage 7 makes
the word a text and the calculation a file:

1. a **word language** and its compiler to `Plan`s (plus the reconstruction that goes with
   each spectral function);
2. a **run file** (JSON) that names the model, the basis, the word and the parameters;
3. a **CLI** `lsquant` (`python -m lsquant`) with `run`, `reconstruct`, `inspect`, `words`;
4. an **MCP tool** that exposes the same four operations to an agent;
5. the **deletion of the legacy tree** (`src/`, `include/`, `test/`, the legacy CMake
   targets, the legacy shell scripts under `examples/`), once every golden is reproduced
   from words. `lsquant/compat/` **stays** (it holds the two legacy-exact kernels that byte
   identity with the stage-0 goldens needs; AGENTS.md section 3 is amended at closure).

Out of scope: new physics, new spectral functions, rank 3, optical conductivity, the
metric study (0006c), performance (stage 8).

## 1. Where things go

```
python/lsquant/words/
  grammar.py        tokenizer + parser: text -> Word (AST: operators, spectral slots, modes)
  compiler.py       Word + basis + parameters -> (Plan, operator list, reconstruction spec); refuses what it cannot compile
  catalogue.py      the table of recognised word shapes (section 3) with their legacy names and output names
  __init__.py       parse(text), compile(word, basis, params), describe(word)
python/lsquant/runfile.py     RunFile dataclass, JSON schema (section 4), validation, resolution of paths and defaults
python/lsquant/cli.py         argparse: run | reconstruct | inspect | words ; `python -m lsquant`; console script `lsquant` (pyproject)
python/lsquant/mcp_server.py  the MCP server (section 6)
python/lsquant/__main__.py
runs/goldens/*.json           one run file per stage-0 golden and per tutorial (section 5), written by you, used by the tester
tests/dev/stage_7_check.py    the single self-check (section 7)
```
No C++ changes except what the compiler needs to express plans it cannot build today
(there should be none; if there is, stop and report).

## 2. Restrictions

All of stages 1 to 6. In addition: the parser is a hand-written recursive descent over a
tokenizer (no parser generator, no `eval`); a word that the catalogue does not recognise
is refused with a message naming the nearest recognised shape, never silently approximated;
mixed bases are refused unless the transform operators are explicit (AGENTS.md section 6);
the run file is validated before anything is loaded; the MCP server only calls the same
Python functions as the CLI and never blocks for more than a dispatch (long runs return a
run id and a log path); nothing under `lsquant/core` learns about words.

## 3. The word language

Tokens: `Tr[` `]`; operator names (`1`, `VX`, `VY`, `VZ`, `SX`, `SY`, `SZ`, `VXSZ`, ...,
any name the basis registers, and `P:<name>` for a projector given as an operator file);
spectral functions `delta(E-H)`, `G+(E)`, `G-(E)`, `dG-(E)/dE`, `U(t)`, `U(t)^+`,
`[X,U(t)]`, `[X,U(t)]^+`; the Fermi-sea marker `sea[...]` around a product; `*` or
whitespace as the product. Case-sensitive. Examples and what they compile to:

| word | plan | reconstruction | legacy name |
|---|---|---|---|
| `Tr[delta(E-H)]` | `dos_plan` | `dos_curve` | `SpectralOp1` |
| `Tr[A delta(E-H)]` | `dos_plan` with `A` as the post operator | `dos_curve` | `SpectralOp<A>` |
| `Tr[VX delta(E-H) VX delta(E-H)]` | `kubo_plan` | `kubo_greenwood` | `NonEqOpVX-VX`, `KuboGreenWood` |
| `Tr[VX delta(E-H) VY dG-(E)/dE]` + h.c. (Bastin) | `kubo_plan` | `kubo_bastin`, `sea`, `surf` | `NonEqOpVX-VY`, `KuboBastin*` |
| `Tr[[X,U(t)]^+ delta(E-H) [X,U(t)]]` | `msd_plan` | `time_correlation_curve` | `CorrelationVX`, `mean...EF` |
| `Tr[P U(t)^+ SX U(t) delta(E-H) P]` | `projected_evolution_plan` | row `m = 0` | `TimeProjSX-PX` |
| `Tr[A delta(E-H) U(t)^+ B U(t)]` | `time_correlation_plan` | `time_correlation_curve` | `TimeCorrA-B` |

The compiler decides which slots are expanded and which are applied from the catalogue
(section 2 of AGENTS.md: `delta` and `G` expanded, `U(t)` applied); it emits the `Plan`,
the list of operators to register (in the engine's index order), and a reconstruction spec
(formula name, kernel, broadening, Fermi energies) for `reconstruct`. On a non-orthogonal
basis (stage 6) the compiler inserts the metric itself from the basis policy: the user
writes the physical word; Löwdin or generalised is a run-file choice (`"route"`), never a
word change. Operator Hermiticity is checked where the plan assumes it.

`words describe "<word>"` prints the parse, the plan as instructions, the operators and
the reconstruction spec; this is what the MCP `describe_word` returns.

## 4. The run file

```json
{
  "model":   {"kind": "wannier", "dir": "goldens/stage0/haldane/operators", "seed": "haldane"}
             or {"kind": "csr", "dir": "...", "label": "haldane", "operators": ["VX", "VY"]},
  "basis":   {"name": "real" | "bloch" | "eigen", "supercell": [48, 48, 1], "disorder": {"W": 1.0, "seed": 7}},
  "metric":  {"policy": "orthogonal" | "nonorthogonal", "route": "lowdin" | "generalised", "applier": "exact" | "chebyshev", "order": 60, "tol": 1e-13},
  "word":    "Tr[VX delta(E-H) VY dG-(E)/dE]",
  "expansion": {"moments": 128, "window": "eigen" | "gershgorin" | [-4.5, 4.5], "trace": "exact" | "random", "states": 512, "seed": 12345, "batch": 128},
  "time":    {"steps": 64, "tmax": 20.0},
  "reconstruction": {"kernel": "jackson", "broadening_meV": 10, "fermi_energies": [0.0, 0.5], "formulas": ["bastin", "sea", "surf"]},
  "output":  {"dir": "out/haldane_hall", "legacy_names": true, "legacy_layout": true}
}
```
Unknown keys are errors; defaults are listed in `runfile.py` and printed by `inspect`.
`legacy_names`/`legacy_layout` reproduce the stage-0 file names and text layouts (this is
how the goldens are matched); with `false` the new `.chebmom` and a single JSON curve file
are written. `inspect <file>` prints the header and dimensions of any moment file and
the resolved run file of a directory (`run.resolved.json` is written by every run).

## 5. Reproducing the goldens from words

Write one run file per stage-0 golden and per tutorial under `runs/goldens/` (G1, G2, G3,
G4, G4c, G5, G6, G7, G9 both models and bases, G10 the four fixtures on both routes;
tutorials 01 to 05 with their exact parameters), such that `python -m lsquant run
runs/goldens/<name>.json` writes files byte-identical to the golden (the tester's gate).
The text round trip of tutorial 03 and the legacy Bessel path of stage 4 are run-file
options (`"legacy_text_round_trip": true`, `"bessel": "legacy"`), documented as such.

## 6. The MCP tool

`python -m lsquant.mcp_server` serves four tools over stdio with the `mcp` package
(FastMCP; `pip install mcp` on TCN): `words_describe(word)`, `run(runfile_json) -> {run_id,
log, output_dir}` (runs in a subprocess, returns at once), `status(run_id)`,
`reconstruct(output_dir, formula, broadening_meV, fermi_energies)`, `inspect(path)`. Every
tool validates its input with the same code as the CLI and returns JSON; no tool reads
files outside the configured root. A `README` under `python/lsquant/` shows the Claude
Code `mcp.json` entry to use it.

## 7. The one self-check (`tests/dev/stage_7_check.py`)

(a) every word of the catalogue parses and prints back to the same canonical text;
(b) five malformed or unsupported words are refused with the expected messages
(unknown operator, rank 3, mixed basis without transform, `U(t)` without `time`, unknown
key in the run file); (c) `runs/goldens/g1.json` and `g5.json` through the CLI give
files byte-identical to G1 and G5 (two `cmp`); (d) the MCP server started in a
subprocess answers `words_describe` and `inspect` with the same JSON as the CLI.

## 8. Deleting the legacy tree

Only after (c) passes: tag the last commit that contains the legacy tree `legacy-final`
(Fable does this at closure), then delete `src/`, `include/`, `test/`, the legacy
`CMakeLists.txt` targets, `scripts/regenerate_goldens.sh` (it becomes a note in
`goldens/stage0/PROVENANCE.md` saying that goldens are regenerated from the `legacy-final`
worktree), the legacy shell scripts under `examples/*/` (keep the Python model generators
the tutorials use). The repository must build and the gate must pass without them.

## 9. Report

`handoff/0007-words-cli-mcp/report.md`: files; the grammar (as a short EBNF); the
catalogue; the run-file schema; the list of run files and which golden each reproduces;
the self-check output; the MCP transcript of (d); deviations; open questions. No "done",
no gate comparison beyond the self-check, no push, no deletion before (c).
