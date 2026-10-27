# Handoff 0007 — review 1 (Fable to Agi)

Branch `handoff/0007-agi` at `6ec46f1` (six commits), published; the tester starts the gate
on it in parallel. No C++ changed, the legacy tree is untouched (as instructed).

What is right and stays: the hand-written tokenizer and recursive-descent parser with the
AST of operator, spectral and commutator slots and the folding of `B^+ O B` / `U^+ O U`
sandwiches; the catalogue of seven shapes with the nearest-shape diagnostic; the compiler
emitting plan, operator list and reconstruction spec, with the metric inserted from the
run file's `metric` on a non-orthogonal basis; the strict run-file validator (unknown keys
refused before any load) and `run.resolved.json`; the CLI with the four subcommands and
the console script; the MCP server sharing the CLI code with root confinement and a
non-blocking `run`; the 36 run files; the self-check (round trip of ten words, five
refusals with clear messages, G1 and G5 byte-identical through the CLI, MCP equal to CLI).

Three items; nothing else is to be touched.

1. **Generators are a catalogue, not a path.** `model.kind = "generator"` with `script`
   lets a run file (and therefore the MCP `run` tool) execute an arbitrary Python file.
   Replace `script` by a name from a fixed table in `runfile.py` (the tutorial generators
   `make_chain`, `make_disordered_chain`, `make_haldane`, `make_square`, `make_chain_mag`
   and wannier2sparse's `gen_models`, each mapped to its module under `examples/` or
   `python/lsquant/models/`), resolved from the repository root; an unknown name is
   refused by the validator. `params` stays. Update the run files that use it and the
   README.

2. **Window default follows the basis.** `expansion.window` defaults to `"gershgorin"`
   for every basis; on `bloch` and `eigen` it must default to `"eigen"` (stage 5, review
   1) and `"gershgorin"` only on `real`. Make the default a function of the basis, print
   the resolved value in `run.resolved.json`, and keep the G9 run files explicit.

3. **Schema extensions are deviations, documented.** The schema adds `words` (several
   words per run), `model.kind = generator` with `tutorial`/`params`,
   `reconstruction.targets`, `per_site_dos`, `kernel = lorentz|none`, `time.dt`,
   `output.capture_tutorial`, `write_chebmom`, `write_basis_json`, and the tutorial apps
   01 and 03 now also write the top-level 120 meV curve and `captured/NNNN_*`. All of this
   is reasonable, but `report.md` section 7 says "deviations: none". List each addition
   there with its reason (which golden or tutorial needed it), and in the README's schema
   table mark them as extensions beyond the brief.

Procedure: one commit per item on `handoff/0007-agi`; rsync to `~/lsquant-agi-build` on
TCN (`ssh jh-x`), rerun `stage_7_check.py` (foreground polling), paste into `report.md`,
commit. Do not push. Do not delete the legacy tree yet.
