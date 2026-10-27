# Handoff 0003 — review 1 (Fable to Agi)

Result (commits `64beba1` to `c01b026`): the rank-2 reconstruction is correct and
complete. Self-check: G5 (`VX-VY`, `VX-VX`) and G4 (random phase, Gershgorin window,
`BandWidth = 15.0366847652085`) `array_equal` through the Python API; sympy Green tables
vs compat tables 1e-14 / 7e-15 relative; all nine golden curves (G5 Bastin, Sea, Surf,
Greenwood; G4 Bastin 64 meV, Sea, Surf, Bastin FFT-grid, Greenwood FFT-grid) with zero
differing lines; tutorial 03 gives `C = +1.017`, the legacy value. Item 0 (`from_csr`
without staging copies, strict dtypes) is done. Names, `compat/` isolation and the
absence of physics in `core/` are as required. The tester is running the gate on this
branch in parallel.

One change before the handoff closes; nothing else is to be touched.

1. **`moments.kubo_plan` guesses its argument order and hard-codes the word.**
   `if num_moments <= 1 and right_index >= 2: swap(...)` is a heuristic that makes a
   wrong call silently "work"; and `word = "VX delta VX delta" if left == right else
   "VX delta VY delta"` writes physics names the function cannot know (the indices
   could be `SZ` and `VY`) into the tensor's metadata. Fix: signature
   `kubo_plan(num_moments: int, *, left_index: int, right_index: int, word: str)`,
   all keyword-only after `num_moments`, no normalisation, no default indices; the
   callers (`stage3_kubo.py`, `tutorial03.py`, `stage_3_check.py`, and the stage-2
   `dos`/`kubo` users if any) pass the word built from the operator names they know
   (`f"{left} delta {right} delta"`). `dos_plan` likewise takes `word: str` from the
   caller (`"1 delta"` today). Outputs do not change; the self-check must still show
   `array_equal` and the gate program the same files.

Procedure: one commit on `handoff/0003-agi`; rsync, rebuild, rerun
`tests/dev/stage_3_check.py` and one `stage3_kubo` form-2 run (Haldane `VX VY`);
update `report.md` sections 4 and 8; commit. Do not push.
