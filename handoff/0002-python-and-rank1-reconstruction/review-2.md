# Handoff 0002 — review 2 (Fable): stage 2 closed

Verdict: **closed**, merged into `cpp-unified`.

Gate on TCN from a fresh build (tester commit `2e3d314` on `handoff/0002-codex`, which
merges `handoff/0002-agi` at `656c6f2`): stage 1 still byte-identical (11 files); stage 2
byte-identical at 1 and 16 threads for the G3 curve from the golden moments and for both
G1 and G3 from operators through the Python API (6 files); code-line `grep` empty;
tutorial 01 from Python reproduces the 22 stage-0 captures byte for byte; the stage-2
figure's residuals are identically zero. Self-check: G1 and G2 `array_equal` through the
binding, sympy delta table vs `math` table 2.2e-16, G3 zero differing lines.

What this stage established: a spectral function written in sympy, turned into a
coefficient table and contracted in C++ with the legacy summation order, reproduces the
legacy reconstruction to the byte; Python is now the library for rank 1.

Correction to review 1: I recorded `SparseMatrix.from_csr` as zero-copy without reading
the binding. The tester found it copies the three CSR buffers into `std::vector`s and
silently converts non-`int64` or non-contiguous inputs. No gate exercises it (both apps
read CSR text), so the stage closes, and the fix is a required item of handoff 0003
(brief section 0, item 0).

Carried forward, no action in this stage: `tutorial01` locates `make_chain.py` by a path
relative to the package (works only from a source checkout) and writes `operators/`,
`BOUNDS`, `exact_*` into `<out>` without listing them; `dos_curve` takes the band width
from `window` (equal to the header for every golden). Both go into the stage-7 cleanup
of the apps.
