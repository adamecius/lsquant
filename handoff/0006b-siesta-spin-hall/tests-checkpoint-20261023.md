# Stage 6b tester launch checkpoint — 2026-10-23

Recovered the preceding developer sessions from Ionos's Codex JSONL history.
The review-3 and review-4 prompts instructed the developer to edit and commit
without running builds or calculations; numerical verification was deferred.
The current delivery is `2c790f19320ad41abfee3043b00ee379f2fe35ac` in
`/root/lsquant-codex`, branch `handoff/0006b-dev`.

## Running work

Isolated root on each node: `~/lsquant-verify-20261023-1100`.
Each delivery was transferred as a Git bundle and checked out detached.
No existing source checkout, running user job or golden was changed.

| Node | Controller PID | Work |
|---|---:|---|
| Jose, `fatnode-tacn` | 689669 | Fresh GCC 11.4 build, historical stage-7 golden gate; stage-8 golden build/gate queued afterward |
| Aron | 638510 | Build controller; all five requested builds completed successfully |
| Aron | 642161 | Existing 6c, 6d, 6b self-checks; standalone stage-8 profiling/kernels; passing branches' material sweeps |

Jose's build passed. Its historical gate uses `S7_SET=list` with every committed
G1-G10/tutorial run, at 16 and 1 threads. New G11 SIESTA runs are excluded from
this historical comparison and queued on Aron. `S7_LEGACY_WT` and all retained
outputs are private to this session. The `legacy-final` tag was fetched locally
for the gate's provenance check.

An initial gate launch used the default all-files selection. Its own launcher
tree was stopped before restarting with the historical selection; the recorded
PIDs are in `fat-results/stopped-pids.txt`. The incomplete log is retained as
`fat-results/stage6b.gate7-initial-incomplete.log` and is not a gate result.

Aron's build passed under GCC 9.4. The default 6b self-check is running,
including the prescribed 20 seeds, 16 states, 60x60 mesh and M=512 noise study.
No completed self-check or golden verdict is claimed at this checkpoint.

## Receipts and remaining scope

Jose: `fat-results/controller.log`, `status.tsv`, `stage6b.gate7.log`,
`controller.pid` under the isolated root. Aron: `aron-results/` contains
`build-status.tsv`, `checks-status.jsonl`, `checks-controller.log`, individual
logs and per-job PID files. The background controllers survive disconnection.

The material queues invoke committed run files, with only output directories
redirected. They stop a branch's sweep on its first nonzero exit; 6d is withheld
after its baseline self-check failure. See the sibling 6c and 6d tester reports.
The 6b material queue covers Bi/PdSe2, the three velocity modes and Lowdin and
generalised routes. Timings, reference-policy comparisons, degeneracy/chunking
diagnostics, complete byte verdicts and figures remain to be reviewed as results
arrive; launching the queue is not stage closure.

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
`KPM_SEED=12345 bash scripts/check_goldens.sh 3 <build>` passed on TCN
(`fatnode-tacn`). This runs the historical stages 1-3 regression checks,
including G1-G5 and G4c, without changing seeds or tolerances. Result:
`GATE PASS (stage 3)`. The complete transcript is
`checkpoint-gate-stage3.log`. OpenMP uses passive waiting and spread/core
binding; no performance conclusion is drawn. This limited regression does
not close stage 6b or supersede the pending/failed material checks above.
