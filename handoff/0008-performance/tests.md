# Stage 8 tester launch checkpoint — 2026-10-23

Recovered delivery: `10d02e43309d40292afc656834ca09975ef4b012` from
`/root/lsquant-codex-8` on Ionos. No implementation edits or integration merge
were made by the tester.

Fresh GCC 9.4, C++17 builds of both golden and native profiles passed on Aron,
including the optional `stage8_kernels` target. Logs and actual compile commands
are under `~/lsquant-verify-20261023-1100/` in `aron-results/` and
`build8-{golden,native}/`. The fat-node GCC 11.4 golden build and historical
stage-7 gate are queued after the stage-6b gate.

Aron's persistent verification controller (PID 642161) queues G5 profiling for
both profiles, then the existing SpMV, SpMM, block, dot and contraction
benchmarks at 1/24/48/96 threads, three repetitions each. Measurements run
sequentially after builds and the new self-checks, before material sweeps.
Each invocation records command, profile, thread count, PID, elapsed time and
exit status. No performance number or speedup is claimed yet.

The delivery does not contain the stage-6b SIESTA implementation. Its requested
Bi performance and per-milestone before/after study therefore require the
integration identified in the developer report. They are not included in the
standalone benchmark queue. No unreviewed integration was invented. The large
graphene workload, NUMA/socket and crossover sweeps, thread-count output hash
comparisons, FFT residual analysis and stage figure remain pending.

The historical gate, benchmarks and stage figure have not passed merely because
the build succeeded. Current receipts are in `aron-results/checks-status.jsonl`
and `fat-results/status.tsv`; this checkpoint is not a closure verdict.

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
