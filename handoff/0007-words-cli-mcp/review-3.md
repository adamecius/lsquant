# Handoff 0007 — review 3 (Fable to Agi): the tester's closure run on `5e8186d`

Tester (`handoff/0007-codex` at `22f92a8`, fat node, 16 threads, all 40 run files): every
golden reproduces byte for byte from its run file: G1, G2, G3, G4, G4c, G5, G6, G7,
tutorials 01 to 05, G9 (96 files, and the CSRs the run files generate equal a fresh
wannier2sparse run), G10 (44 files, 11 run files including the new chebyshev and ischo
ones; worst moment difference 9.5e-12 with the stage-6 measure `max|Δ|/max|gold|`), MCP
(7/7: run/status of g1 and g5, out-of-root refused with the reason in the result). The 11
CLI refusals exit 1 with a message and no output directory. The deleted-tree checks (18)
pass: no `src/ include/ test/`, no legacy CMake targets, `legacy-final` builds the six
legacy drivers and reproduces G1, G2, G5 from that worktree. The 1-thread repetition runs
on Aron (machine-local comparison) and is reported separately.

One golden does not reproduce, and it is the golden's fault, not the code's.

## Decision: tutorial 06 is re-blessed from the new code (legacy defect D5)

142 of 281 tutorial-06 files differ: the 139 curves of the identity-operator word
(`TimeProj1-P*`), `ey_ratelaw.dat` (1.8e-16 relative) and `captured/index.tsv`. Cause,
localised by the tester: the stage-0 pipeline ran the legacy driver with operator `1`; the
driver loads `operators/<label>.1.CSR`, which the generators never write; the legacy loader
(`BuildOPFromCSRFile`, `legacy-final:include/sparse_matrix.hpp:132`) returns silently on a
missing file, so the golden "identity" denominators came from an operator that was never
built (values ~1e-6, noisy). The new code applies the identity and gives the constant the
word has at `t = 0`. The tutorial normalises by this curve, so the stage-0 figure was right
by cancellation. As with G4c: the legacy curves move to
`goldens/stage0/legacy_defects/tutorial06_identity_operator/`, the new-code files become
the golden, PROVENANCE says why, AGENTS.md section 10 gets entry D5 ("silent loader on a
missing operator file"). The tester does the mechanics on `handoff/0007-codex` (it holds the
new-code outputs); nothing else under `goldens/` changes.

## Items for Agi, one commit each, on `handoff/0007-agi` after `5e8186d`

1. **Tutorial 06 honours the run file or refuses it** (same rule as review-2 item 4):
   `apps/tutorial06.py:125-127` builds its own word strings and only validates the run
   file's `words`. Pass them through (the app takes the words, the expansion, the time and
   the reconstruction from the run file) or declare in `tutorial06.json` only what the app
   uses and refuse the rest. No field accepted and ignored.
2. **No implicit generation: remove `ensure_w2s_csr_dir`** (review-2 item 6 went too far).
   It runs an external binary silently, writes into the repository tree (`build/w2s_csr`),
   depends on a path that exists on one machine, skips silently when the binary is absent
   and then fails with a message that does not name the cause, and never regenerates a
   stale directory. Replace it with a refusal: when `csr_dir` lacks `<seed>.HAM.CSR` (or its
   receipt disagrees with the supercell), the run stops before any load with a message that
   names the directory, the mismatch, and the exact `wannier2sparse` command that produces
   it. The gate (`scripts/check_goldens.sh`) generates the CSRs explicitly, as the stage-5
   gate did. The G9 real run files keep `csr_dir`.
3. **`tests/dev/stage_7_check.py` takes the build directory from the environment**
   (`LSQUANT_BUILD_DIR`, default `build/`), like the other checks; no hard-coded
   `REPO/build` (lines 370-372).
4. **The rest of the legacy material goes** (brief section 8, "the end of the legacy
   tree"): `legacy/` (19 tracked driver sources: the `legacy-final` tag holds them),
   every mention of `inline_*` drivers in `docs/` and `examples/`, and
   `examples/06_spin_relaxation_square/lsqrelax.py` and `rp_spin.py`, which call drivers
   that no longer exist: rewrite them as run files (`runs/goldens/tutorial06.json` is the
   model) or delete them with a line in PROVENANCE. `lsquant/compat/` stays.
5. **Report.** Append "Review 3 applied" to `report.md`: one line per item; the
   `stage_7_check.py` output on the fat node; the tutorial-06 run file through
   `python -m lsquant run` with `cmp` against the re-blessed golden (fetch
   `origin/handoff/0007-codex` for it after the tester's commit). Do not push.

Closure after this: the tester reruns the 16-thread gate on the fat node (the byte verdict,
with the re-blessed tutorial 06, expected 281/281) and the 1-thread repetition on Aron;
then stage 7 merges into `cpp-unified`, AGENTS.md section 3 is amended (compat stays, legacy
tree gone, `legacy-final` is the provenance), and the final validation run follows.

Noted, no action in this stage: the dense exact Löwdin chain1d run files take 17 to 42
minutes at 16 threads and tutorial 06 about two hours on the fat node (stage 8).
