# Stage 6c tester checkpoint — 2026-10-23

Tested delivery: `b2ab659af4bba4e64fecded4c3d579c28023b8c3`, recovered from
`/root/lsquant-codex-6c` on Ionos. No source changes, commits, pushes, golden
changes or tolerance changes were made by this verification session.

Machine: `tcn-aronnode`; system GCC 9.4.0, C++17 Release, Eigen from
`~/lsq-env`, Python 3.10.22, NumPy 2.2.6, SciPy 1.15.2, sisl 0.16.4.
Fresh configure and build both exited 0. The isolated checkout and build are
`~/lsquant-verify-20261023-1100/stage6c` and `build6c-golden`.

## Observed result

The existing `tests/dev/stage_6c_check.py` passed without argument overrides,
with `KPM_SEED=12345`, `OMP_NUM_THREADS=48`, static core binding and one BLAS
thread. Wall time: 5.077 s; maximum RSS: 130784 KiB; exit status: 0.

- Symbolic potential, constant-shift, orthogonal-limit and periodic velocity
  residuals are zero; all three velocity substitutions pass.
- Clean Bloch moments equal the stage-6b moments exactly, for DOS and SHC,
  under both Lowdin and generalised metrics.
- W=0 real/Bloch moment maximum absolute difference is at most 7.665e-15.
- W=1, seed 7, 12x12: DOS/SHC relative maximum difference is at most
  1.553e-15, below the prescribed 1e-9.
- The orthogonal stage-5 DOS moment byte comparison passes.

Full log: `figures/out/verification_20261023/stage6c-selfcheck.log` locally;
`~/lsquant-verify-20261023-1100/aron-results/stage6c-selfcheck.log` on Aron.

## Pending

The twelve committed Bi disorder run files are queued after the existing
self-checks and stage-8 standalone measurements. The queue stops that material
sweep on its first nonzero exit. Material results, thread-count comparisons,
the stage figure and reviewer acceptance remain outstanding. This self-check
does not close the stage or certify a golden.

## Low CPU investigation (user-requested stop)

The Aron controller and its current material child were suspended by their
recorded PIDs (642161, 648134, 648135). Jose's independent gate was left running.
The stop receipt is `aron-results/user-stop.json` on Aron.

The earlier reported 1.3 cores came from `ps %CPU`, a lifetime average, not an
interval measurement of the current engine phase. The first Bi run spent
533.356 s certifying dipoles (largest quadrature order 192, margin 10.0042).
Its certified cache was written successfully. A separate diagnostic reused
that same entry in 0.010875 s; no cache bypass or tolerance change was made.
All 48 worker threads have distinct CPU bindings (cores 0 through 47).

A diagnostic with the same Bi model and metric, M=256 and one trace vector,
completed its SHC moments in 0.730358 s at 48 threads. Full-order interval CPU
measurements and a 1-versus-48-thread output comparison are in progress under
`aron-results/diagnostic/`. Production remains suspended until verification.

### Verification and resumption

The subsequent M=1024, two-state comparison completed its moment phases:
SHC 406.085441 s at one thread and 19.456420 s at 48 threads (20.87x).
Interval CPU samples during the 48-thread SHC phase average 46.766 cores,
maximum 47.567. The single-thread diagnostic used core 0; the 48-thread
diagnostic used cores 48-95, sequentially for the measured moment phases.
All four produced DOS/SHC moment files (native and legacy layouts) are
byte-identical. Both diagnostics were stopped during subsequent SHC
reconstruction; no completed curve comparison is claimed for them.

Production was resumed only after those comparisons. The low lifetime CPU
average was dominated by the first serial dipole certification, whose result
is now cached and verified reusable. No numerical source or tolerance change
was required. This does not optimise cold certification for a different input,
nor guarantee full-core use during reconstruction. Resume receipt:
`aron-results/resumed-after-diagnosis.json`. Use interval CPU measurements,
not process-lifetime `ps %CPU`, when assessing an active engine phase.

## Preservation checkpoint, 2026-10-23

This branch preserves the existing tester report and available artifacts at the
maintainer's request. Earlier statements about running jobs describe the time of
the original observation; this snapshot does not certify their present status.
The development tip is recorded in `checkpoint-20261023.json`. No integration
merge, implementation change, reference change or stage closure is made here.

The complete original uncommitted files, including raw numerical receipts and
logs, are retained outside Git in
`C:/Users/adame/LSQUANT/checkpoints/20261023/local-uncommitted.tar.gz`, with
`local-uncommitted-manifest.json` recording their SHA-256s. Selected figures and
small summaries are preserved in this branch. Material arrays, generated moment
files, build directories and independent nested repositories are not staged.

Before this checkpoint commit, a fresh Release build and the existing
`KPM_SEED=12345 bash scripts/check_goldens.sh 3 <build>` ran on TCN
(`fatnode-tacn`). This checks historical stages 1-3, including G1-G5 and G4c,
without changing seeds or tolerances. Result: `GATE PASS (stage 3)`;
`== 47 files compared; failures: 0`. Exit status: 0.
The transcript is `checkpoint-gate-stage3.log`. The three independent gates
ran concurrently with passive OpenMP waiting, spread/core binding and one
BLAS thread. Their wall times are not performance measurements.
This limited regression does not close this stage or supersede its specific
pending or failed validation above.
