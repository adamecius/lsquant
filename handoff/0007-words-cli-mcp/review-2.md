# Handoff 0007 — review 2 (Fable to Agi): the tester's run on `9f8ccc8`

Tester (`handoff/0007-codex` at `65f8510`): every golden a run file produced is
byte-identical through `python -m lsquant run` at 16 and 1 threads: G1, G2, G3, G4, G4c, G5,
G6, G7, tutorials 01 to 05 (194 and 17 files for 02 and 05), the eight random-phase G9
files; MCP `run`/`status` of g1 and g5 identical, out-of-root refused; six CLI refusals
with exit 1 before any load. Two run-time defects stop the rest, and eight smells.

Items, one commit each; nothing else is to be touched. Items 1 and 2 are defects.

1. **D1: the metric accessors are properties, the Python calls them as methods.**
   `bindings.cpp` exposes `overlap`, `inverse_overlap`, `inverse_sqrt_overlap` as
   read-only properties; `runfile.py` and `words/compiler.py` call them (`...()`), and
   `_core.pyi` declares methods. Every G10 run file dies with `object is not callable`.
   Decide one convention (properties, as the binding has) and fix the eleven call sites
   and the stub; the tester's scratch fix gives byte identity with `stage6_nonorth` on
   all four fixtures.

2. **D2: `expansion.states` absent means `None`, not the default.** `compiler.py` reads
   `exp_cfg.get("states", default)` on a dict whose key is present with `None`, so every
   exact-trace G9 run file dies at `int(None)`. Treat `None` as absent (and audit the
   other `.get(key, default)` reads of the resolved run file for the same pattern).

3. **A word is refused, not approximated** (brief section 2). `Tr[VX delta(E-H) VY G+(E)]`,
   `... G-(E)` and `... dG+(E)/dE` currently compile to Kubo-Bastin. Only the catalogue's
   `dG-(E)/dE` form compiles; the others are refused with the nearest-shape message
   (`G+`, `G-` and `dG+/dE` are tokens of the grammar for future shapes, not aliases).
   Likewise `window: "eigen"` on a `csr` model is refused with a message ("no k blocks;
   use gershgorin or an explicit window"), not silently replaced.

4. **Generator run files honour their parameters, or refuse them.** For tutorials 01,
   02, 05 the run file's word, expansion and params are validated and then ignored by
   `run_tutorial0N(out_dir)`. Either the generator path passes them through (preferred
   when the tutorial app already takes them), or the run file for a tutorial declares
   only what the app uses and the validator refuses the rest. No field may be accepted
   and ignored.

5. **No implicit words.** When a VX-VY run asks for `greenwood`, `execute_runfile`
   appends the Greenwood word silently. Require it in `words` (the run files that need
   it list both words) and refuse a formula whose word is not present.

6. **G9 real-basis run files use the wannier2sparse CSR** (`csr_dir`) like the stage-5
   gate, so the configuration they reproduce is the gated one.

7. **MCP errors carry their reason** in the tool result (`isError` with the message),
   not only on the server's stderr.

8. **Coverage.** Add `runs/goldens/tutorial06.json` if `goldens/stage0/tutorials/06/`
   exists (G8 lists tutorial 06; if the stage-0 capture does not exist, say so in the
   report and skip); add G10 run files for the `chebyshev` and `ischo` appliers on one
   fixture each; the self-check (c) runs, in addition to g1 and g5, one exact-trace G9
   file and one G10 file (the two classes that reached the tester).

Procedure: after the legacy-deletion commit (same branch), one commit per item; rsync to
`~/lsquant-agi-build` on TCN, rerun `stage_7_check.py` and the two newly covered run files
with `cmp` against their goldens, paste into `report.md`, commit. Do not push.
