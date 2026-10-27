# Handoff 0002 — review 1 (Fable to Agi)

Result of the implementation (commits `c01b0f4` to `6397f79`): correct and complete.
The self-check reproduces G1 and G2 to the bit through the binding, the sympy delta
table agrees with the `math` table to 2.2e-16, and `dos_curve` on G1 gives the G3 curve
with zero differing lines; Fable's own `cmp` of the form-A output on TCN against G3 is
byte-identical. The binding is zero-copy where the brief asked (`from_csr` buffers,
`data` view with the parent kept alive), the Python package has no `print`, the C++
additions name no physics, and both declared deviations are accepted (the two
`ContractRank1` overloads; the `const` accessors from stage 1).

Two changes before the handoff closes; nothing else is to be touched.

1. **`python/lsquant/apps/tutorial01.py` knows the golden capture layout.** It writes
   files named `0013_meanmoments_...`, `0019_...`, `0020_...`, `0021_...` and a
   `captured/` sequence that exist only because the stage-0 tester numbered the files
   it captured. A library application must not encode the layout of a test fixture.
   Fix: keep exactly two kinds of output, the moment files with their legacy names at
   the top of `<out>` and the reconstructed curves under `<out>/eta<meV>/` with their
   legacy names (`meanmoments_chain1d_N<N>_M512JACKSON.dat`); delete the prefixed
   copies and the `captured/` sequence. The tester's figure script does the matching
   (its instructions are being updated accordingly).

2. **`SpectralFunction.__call__` returning `self`** serves nothing and invites misuse
   (`DeltaFunction()` looks like a constructor). Remove it.

Accepted as is, for the record: the extended raw-pointer `ContractRank1` signature with
strides (it is what lets `.real` of a complex array be read in place); the Python-side
loops over the grid in `dos_curve` (they reproduce the legacy operation order, speed
is irrelevant); `stage_1` comments containing the notation `delta_{m0}` (mathematical
notation in a comment is not a physics name in code; the tester's `grep` is being
narrowed to code lines).

Procedure: two commits on `handoff/0002-agi`; rsync to `~/lsquant-agi-build`, rebuild
(no C++ change expected), rerun `python -m lsquant.apps.tutorial01 --out /tmp/stage2_tut01`
and paste the new file list in `report.md` (section 5); commit the report. Do not push.
