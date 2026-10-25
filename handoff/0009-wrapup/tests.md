# Stage 9 tester execution report

Checkpoint: 2026-10-24. Status: in progress, not a closure verdict.

## Revisions and integration preflight

The four delivery heads still match plan.md. No integrated candidate exists
on origin/cpp-unified. Read-only merge-tree checks identified conflicts in
runfile.py for each combination, plus block_diagonal_operator.cpp and
moment_engine.cpp for 6b + 8. The working tree and delivery source files
were not modified. See [integration-preflight.md](integration-preflight.md)
for exact trees and required semantic preservation.

The final golden/material rerun cannot certify the separate branch heads
as one implementation.

## W02: fresh reproduction and localisation

On Aron, stage-6d revision
22965bfe40f75a6d479f25b906c5ccc069ea8827, build6d-golden, the existing
tests/dev/stage_6d_check.py was rerun with its defaults unchanged. A
tester wrapper captured returned moment tensors without altering them,
then allowed the original assertion to fail normally.

Result: FAIL at line 147, unprojected Lowdin exact trace, M=1024.
Original oracle: 1.9999237004848083 at the three gap energies.
The source, assertion and 0.001 limit are unchanged.

The tester diagnostic figures/stage_9_projection_diagnosis.py independently
formed the exact moments from the eigenvalues and physical current/velocity
matrix elements, using polynomial recursion and dense NumPy contraction.
It did not use MomentEngine for that independent sum.

- Leading 64x64 engine-versus-dense moment discrepancy: 6.73397e-16 maximum.
- Independent dense M=1024 reconstruction reproduces the failed plateau
  within about 1e-14.
- Doubling the reconstruction grid from 30*M to 60*M changes the maximum
  error from 0.0170613418543 to 0.0170612241207.
- The source window is [-9.3, 66.3] eV, width 75.6 eV. Nominal resolution
  pi*width/(2*M) is 115.97 meV at M=1024. Requesting 10 meV broadening does
  not create missing moments.

| Model | M | Maximum SHC error at gap probes |
|---|---:|---:|
| Original | 256 | 0.29863251 |
| Original | 512 | 0.06824509 |
| Original | 1024 | 0.01706134 |
| Original | 2048 | 0.00429533 |
| Shift +/-3.1 eV | 256 | 0.00577751 |
| Shift +/-3.1 eV | 512 | 0.00145836 |
| Shift +/-3.1 eV | 1024 | 0.00036661 |

Errors are in e^2/h with sigma_z Pauli. Exact trace removes stochastic
error in this diagnosis. The observed failure is a finite-order
resolution requirement in the wide baseline, not a discrepancy between
the engine and the independent moment evaluation. This finding does not
waive the delivered self-check or establish every subsequent assertion.

The developer/reviewer must provide a convergence-valid baseline comparison
while preserving the accuracy requirement. Proposed choices are recorded
in integration-preflight.md; none has been applied by the tester.

Figure: figures/out/verification_20261024_wrapup/diagnosis/selfcheck-convergence.png.
Caption: Exact-trace SHC error on the existing augmented Kane-Mele fixture,
original and shifted spectra, versus moment count. The line is the existing
0.001 limit. The shifted model reaches it at M=1024 while the original
wide spectrum remains under-resolved. This is a diagnostic, not closure.

## W07: orthogonal toy checks previously hidden by early failure

Executed the existing toy's orthogonal-tool calculations separately:

- Training-band maximum error: 3.9968e-15 eV.
- Off-grid maximum band error: 2.62013e-14 eV, 40 points, seed 1.
- These observations satisfy the toy's stated 1e-6 band limit.
- Material off-grid/radius validity and full self-check completion remain open.

## W08: FFT execution defect

The existing bench/compare_contraction.py was run on committed G4/G5
moments, Greenwood and Bastin, golden/native builds, 1/48 threads.
All 16 invocations failed in the FFT backend at:

    python/lsquant/reconstruct.py:578
    out[i] = green.imag() / sine / math.pi
    TypeError: 'float' object is not callable

The table backend executes first. The proposed correction is the property
access green.imag. It is recorded in integration-preflight.md but has not
been applied. FFT arithmetic, residuals and speedup remain unverified.

These are accuracy diagnostics, run alongside other diagnostic work;
their timing must not be reported as isolated performance measurements.

## W03: material error decomposition

Completed for the covariant 24x24 cases on the same frozen SIESTA inputs
as the accepted projection study. The diagnostic
figures/stage_9_material_error_budget.py forms the
exact projected trace in retained eigenstates, preserving the original
Hilbert-space normalization and the production spectral window. The
discarded complement has zero projected current/velocity matrix elements.

The comparisons separate exact-trace KPM/oracle discrepancy from the
saved R=16 sampled-minus-exact-trace difference at the same M. This is not
a dimensional-reduction change to production or a replacement material gate.

| Material/probe | M | Exact oracle | Exact-trace KPM | Saved R=16 KPM | Exact KPM minus oracle | R=16 minus exact KPM |
|---|---:|---:|---:|---:|---:|---:|
| Bi, E=0 eV | 256 | 1.07975855 | 1.16184298 | 1.21716593 | +0.08208443 | +0.05532295 |
| Bi, E=0 eV | 1024 | 1.07975855 | 1.08559928 | 1.13534094 | +0.00584072 | +0.04974167 |
| PdSe2, E=1.052 eV | 256 | 1.98650292 | 1.98025082 | 2.27634536 | -0.00625211 | +0.29609455 |
| PdSe2, E=1.052 eV | 1024 | 1.98650292 | 2.12607643 | 2.52038425 | +0.13957350 | +0.39430783 |

All values are SHC in e^2/h with sigma_z Pauli, shift +/-5 eV. The
independent selected-state oracle agrees with the prior full projected
oracle within 9.15e-7 for Bi and 4.04e-9 for PdSe2. That residual is
recorded, not silently removed. Original input hashes were checked.

Bi's exact-trace discrepancy decreases with M, while the 16-vector sample
remains offset. PdSe2's exact-trace discrepancy also worsens at M=1024,
so increasing random-vector count alone cannot solve its result. Further
moment/kernel and integration convergence is required before attributing
the remaining exact-trace discrepancy to a specific formula defect.
The oracle uses Gaussian occupations; KPM uses its delivered Jackson
reconstruction, so these are not identical finite-width kernels.

This decomposition covers two configurations, not the whole 16-run sweep,
and is not a multi-seed uncertainty estimate. The small-basis toy
independently verifies the dense moment formula against the engine; no
full exact-trace material MomentEngine pass was substituted or claimed.
Exact projected traces here were evaluated by dense retained-state sums.

Local compact receipt:
figures/out/verification_20261024_wrapup/material-error-budget.json.

## Receipts and process ownership

Remote root: ~/lsquant-verify-20261023-1100/wrapup-20261024/.
The self-check, diagnostic and FFT launchers each recorded their PIDs.
Only threads of the recorded material diagnostic child were adjusted for
affinity; the change is recorded in diagnostic-affinity.json. It is not
a performance benchmark.

Local receipts: figures/out/verification_20261024_wrapup/.
Captured tensors and material exact tensors remain on Aron.
Existing full-model queues remain paused. No multi-day validation,
commit, push, golden replacement or tolerance change was performed.

## Authorized integration and candidate v1

The user explicitly overrode the tester-only restriction on 2026-10-24
("Overwrittrn"). Codex is implementing and validating the integration in
the isolated `codex/lsquant-wrapup` worktree. No commit or push was made.

Candidate v1 has staged source tree
`20f8f795552e60d196408257ab04c40dfc09b2d3`, based on 6b
`2c790f19320ad41abfee3043b00ee379f2fe35ac`, with the 6c, 6d and 8
delivery diffs applied. The source archive was transferred to fresh
`~/lsquant-wrapup-20261024/candidate-v1` checkouts on Jose and Aron.
Reports added afterward are not part of this source-tree identifier.

Conflict decisions preserve both run-file parameter sets, reference and
projection fields, the 6b Hermiticity checks, native first-touch allocation,
and the 6b fixed-order parallel golden accumulation. Projection dispatch
precedes disorder dispatch so a requested projection cannot be ignored.
The FFT Python complex-property fix (`green.imag`) is applied.

Jose's golden build passes. Aron's first build selected system Eigen 3.3
and failed because the existing source uses the Eigen 3.4 Vector alias.
Reconfiguration with the previously validated environment's Eigen headers
fixes this dependency mismatch; golden and native builds pass. The original
failed logs are retained. No source workaround for the old headers was made.

The first Jose Haldane comparisons now complete on both FFT and tables.
Their grids are identical. Maximum absolute differences are 1.66e-9 for
Greenwood and 5.90e-11 for Bastin; relative to the curve maximum these are
1.30e-12 and 5.47e-13. These are numerical comparisons, not byte passes.
Pointwise relative differences near zero are larger and remain reportable.

Aron is running the full 16-case FFT comparison and the existing 6b/6c
self-checks on the combined candidate. The existing 6d self-check is queued
with its supported `--moments 8192` option, retaining equal M for both spectra
and every existing assertion and tolerance. This tests convergence without
weakening the contract or replacing its engine calculation with a dense oracle.
Results are pending; M=1024 remains a documented failing invocation.

## Integrated validation checkpoint: candidate v5

This section supersedes the pending statements above. Source tree
bf70c9b8c8b8ec66df329265d6c170a1cb315e0e is reproduced on Aron and Jose
under ~/lsquant-wrapup-20261024/candidate-v5. It contains the integrated
deliveries, FFT crash correction, exact-overlap correction, opt-in folded
Bastin contraction, clean projected CSR adapter, and SIESTA run-file
backend selection. The implementation report records the differences.
No golden, seed or numerical gate tolerance was changed.

### Existing checks

| Check | Result and scope |
|---|---|
| 6b, v1 | Existing developer check passes, including sampler study |
| 6c, v1 | Fails the orthogonal-limit byte assertion; the phase-factored identity overlap was rounded |
| 6c, v2 golden and native | All original assertions pass after restoring exact diagonal Bloch entries; stage-5 orthogonal DOS byte assertion is true |
| 6d, v3 native, M=8192 folded | All original assertions pass, Lowdin and generalised, original and projected exact-trace KPM at the same M; symbolic and off-grid orthogonal checks reached and passed |
| 6d, default M=1024 | Remains underresolved and fails the original plateau limit; not reclassified as a pass |
| Clean local projection, v5 golden | 16 Bloch/CSR comparisons pass; eight run-file/resolved replay pairs pass for selected DOS and SHC; moments and saved curves replay byte-identically |
| 6b/6c, v5 native | Both complete with exit 0; all existing assertions pass, including the 20-seed study and orthogonal-limit bytes |
| 6d, v5 native | Fresh M=8192 folded rerun completes with exit 0; every assertion passes, both metric routes, original and shifted, symbolic/local/orthogonal/run-file checks |
| Historical G1-G10/tutorial inventory, v5 Jose | Fresh golden build passes; full 64-run, 16/1-thread sweep is running; no final verdict yet |

The 6d unprojected maximum plateau error is about 0.000270887;
shifted maximum error about 0.000290853, below the unchanged 0.001
limit. The off-grid orthogonal toy band maximum is 2.31e-14 eV.
Its previous M=8192 table attempt was paused during the enormous serial
coefficient-table fill, then explicitly cancelled after the folded
replacement passed. Its pause and cancellation receipts are retained.

Clean local comparison: hard and shift, both metric routes, 4x3 mesh,
M=32, R=2, seed 12345, batch 11, full radius and 1.2 Angstrom radius.
Maximum moment residual is 2.382e-15 against the existing 1e-10 W=0
limit. This is a basis/layout check, not a transport-convergence pass.
The input loader alone supplies the existing in-memory toy in the run-file
harness; dispatch, parser, compiler, adapters, engine and writers are real.
The final harness also compares folded and table SHC on identical grids.

Radius 0/0.8 Angstrom hard-projection attempts were refused by relative
Hermiticity validation: norm 5.47687e-14, anti-Hermitian norm 2.55228e-17.
Those cutoffs remove every physical hopping of this toy. The limitation
is retained and logged, rather than weakening validation. Radius 1.2
retains nearest-neighbour bonds; its changed DOS is visible in the figure
and is not accepted as an accurate material cutoff.

Figure: figures/out/verification_20261024_wrapup/integrated/checks-v5/local/local-projection.png.
Caption: selected toy DOS with all hoppings versus a 1.2 Angstrom cutoff;
the lower panel compares Bloch and CSR moments for the same model in
each of 16 configurations. Layout error is near round-off, while radius
truncation visibly changes the model. Figure visually inspected.

### Reconstruction comparison

All 16 original FFT invocations complete (G4/G5, Greenwood/Bastin,
golden/native, 1/48 threads). This removes the crash; it is not a byte
acceptance. Jose's G5 FFT/table peak-relative discrepancies are
1.30e-12 Greenwood and 5.47e-13 Bastin. Pointwise relative errors near
zeros are larger. The stage-8 FFT golden verdict remains open.

| Uniform-grid Bastin comparison | Maximum absolute difference | Relative to peak | Tables median | Folded median |
|---|---:|---:|---:|---:|
| G4 graphene | 6.7502e-14 | 3.8909e-15 | 0.04779 s | 0.006286 s |
| G5 Haldane | 4.8479e-13 | 4.5015e-15 | 0.16777 s | 0.018717 s |

These are diagnostic measurements, about 7.6x and 9.0x respectively,
not the isolated full-workload performance campaign. Grid arrays are
byte-identical; curve arrays are not. Peak-relative error must not be
substituted for the project's pointwise or byte criteria.

Saved material M=1024 tensors give folded/table maximum differences
1.6862e-12 for Bi and 1.2799e-12 for PdSe2 (peak-relative 3.89e-14 and
3.50e-14); the energy grids are identical. Folded reconstruction took
about 0.92 s in each diagnostic. Defaults retain the table path.

### Material convergence diagnosis

Frozen input SHA values are unchanged:
Bi e17355e699678ea4f1f6228a1ce436c5e83d6cd90b4d7f61942945f2f170f6bc;
PdSe2 588c6364cb34282d9bb7e1343895c02743bb88651dd7533336a2eb8f71d82d27.
These diagnostics evaluate exact retained-state sums, without stochastic
trace error; they do not replace the material MomentEngine gate.
Configurations: covariant, 24x24, shift +/-5 eV.

| Material / energy | M | Broadening | Exact projected KPM | Gaussian projected oracle | Difference |
|---|---:|---:|---:|---:|---:|
| Bi / 0 eV | 4096 | 1 meV | 1.08013473411 | 1.07975855226 | +0.000376182 |
| PdSe2 / 1.052 eV | 1024 | 1 meV | 2.12607642788 | 1.98650998840 | +0.139566439 |
| PdSe2 / 1.052 eV | 2048 | 1 meV | 2.01973608429 | 1.98650998840 | +0.033226096 |
| PdSe2 / 1.052 eV | 4096 | 1 meV | 1.99490085407 | 1.98650998840 | +0.008390866 |

The Bi discrepancy is about 0.0348%; PdSe2 decreases from 7.03% through
1.67% to 0.422%. At 10 meV, the existing broadening rule caps PdSe2 at
1885 effective moments: M=2048 and M=4096 therefore both give
2.02597118614. This explains why increasing stored M alone cannot resolve
that residual. KPM/Jackson and Gaussian oracle kernels differ. No new
SHC acceptance limit is inferred from the user's 10 meV band criterion.
The accidental closeness at M=256 is not evidence of convergence.

Receipts: integrated/checks-v4/{bi,pdse2}-convergence/error-budget.json
and their logs under the figure output root above. Exact M=4096 tensors
remain on Aron. Retained-band acceptance, projection transport bias,
finite-M/kernel error, and stochastic error remain separate.

### Historical regression and remaining work

The v1 historical sweep was intentionally cancelled on its recorded
process tree after the final source changed. Its completed generalised
case and interrupted Lowdin outputs were preserved. A brief v5 launch
with the unchanged v1 C++ binary was also replaced: the gate requires
CMAKE_HOME_DIRECTORY to identify the tested checkout, so v5 gets a
fresh build. These attempts have no final gate verdict.

Current Jose controller: historical-v5/controller.pid (2069007 at launch).
The fresh build log is logs/build-v5-golden.log. Gate log:
historical-v5/gate7.log. A recorded follow-on controller (2070033)
runs the separate 6only gate only if stage 7 passes. The G11 material
files are explicitly excluded from this historical inventory and remain
assigned to the separate material gate, not counted as passes.

Completed Aron controller: checks-v5/selfchecks-controller.pid (697946).
checks-v5/selfchecks.exit is 0. Local copies of the final logs and
receipts are saved under the integrated/checks-v5 output directory.
Only owned, recorded processes have been stopped. Original full-material
and disorder queues remain paused.

Outstanding closure items are W03 acceptance, W04-W09 remaining material,
disorder, projection and performance requirements, W10 completed historical
and full-model long validation, and W11 written review. Clean local
run-file support is implemented, but W06 projected disorder is still
unsupported. The project remains open.

### Live checkpoint, 2026-10-24 14:26 CEST

Jose's v5 controller has run for 4 h 37 min. The 16-thread pass has
completed its production/comparison phase: 703 files reported IDENTICAL
and 33 numerical/reference checks reported PASS. No DIFFERS, FAIL,
ERROR, or nonzero run-exit line is present so far. Informational G10/G5
surface-curve discrepancies remain in the log and are not byte verdicts.
The one-thread pass is active in g10_chain1d_real_lowdin_chebyshev:
about 16 minutes at 99.8% CPU (one requested core), with four preceding
one-thread runs completed. No final gate exit file exists yet.

The harness defaults omit the two dense exact Lowdin chain cases in
its second thread pass. To avoid carrying that historical exception
into final coverage, a recorded follow-on controller (2553875) queues
both at one thread with S7_T1_SKIP empty after the main gate passes.
Its script is run-v5-one-thread-completion.sh and its outputs will be
under historical-v5/one-thread-completion. The separate 6only queue is
unchanged. Aron self-checks and local replay checks are all complete;
the original long material queue remains paused.

### Subsequent checkpoint: historical stage-7 gate finished

The v5 stage-7 gate has finished with exit 1: 1542 comparisons/checks
in its final tally, 1413 IDENTICAL lines, 82 PASS lines, and one FAIL.
No numerical DIFFERS or failed run-exit line was reported. The sole
failure is the legacy-name structural check: old executable names remain
in AGENTS.md, README.md, two tutorial templates, moments.py, and the
stage-4 developer check. All six files are byte-identical to the 6b
base commit, so this is an inherited unresolved structural failure,
not a new numerical integration regression. It is not waived.

The guarded follow-on queues correctly did not run after that failure.
After localising it, independent launchers explicitly resumed 6only
(controller 3111292) and the two missing one-thread exact Lowdin cases
(controller 3111293). Their outputs use the original gate scripts and
settings; the failed stage-7 result is preserved. Reason and file
comparison evidence are in historical-v5/independent-continuation.json.
The earlier .blocked receipts describe the initial guarded queues;
subsequent .exit files from these independent runs determine completion.

Aron still has completed passing 6b, 6c, 6d and local-replay checks;
the original full-material process remains paused. Project closure
requires resolving the structural failure and completing the remaining
coverage already listed above.

### Checkpoint, 2026-10-25: independent numerical coverage completed

Both independent controllers have finished. Stage 6only returns exit 0:
98 comparisons/checks, zero failures (24 IDENTICAL and 74 PASS lines).
The two missing exact Lowdin one-thread cases both completed, taking
9117.69 s (complex) and 8641.15 s (real), about 4 h 56 min combined.
Their stage-7 wrapper reports 70 comparisons/checks and one failure:
the same inherited legacy-executable-name structural check in six files.
It reports 27 IDENTICAL and 28 PASS lines, with no numerical difference
or failed production run. The wrapper's exit 1 is preserved.

The v5 historical numerical coverage and stage-6-specific validation
are now finished; the historical gate still lacks a passing structural
verdict. The recorded Jose controllers are no longer running. Aron
6b/6c/6d and local checks remain complete, and the original material
process remains stopped (T state). No new production material campaign
or long full-model run has been started by this checkpoint.

### Requested material review figures, 2026-10-25

The user requested two figures, confirmed as Bi and PdSe2: bands,
DOS, longitudinal conductivity, spin Hall response, and valence spin
texture. The reproducible figure harness is
figures/stage_9_material_review.py. Final PNG/SVG, plotted data,
captions and numerical provenance are in figures/out/review_20261025/.
REVIEW.md records the numerical values and eight plotting-style checks.
Both final PNGs were visually inspected.

The figures compare the full frozen SIESTA models with the implemented
shift +/-5 eV projection at 24x24, covariant velocity, exact Lowdin,
W=0. They show original/projected textures rather than inventing two
spin partners for Bi's nondegenerate upper valence band. PdSe2's
degenerate points use the stated projected-sigma_z convention.
Retained-band path MAE is about 2.5e-12 meV and original/projected
spin-texture differences are below 1e-14 in Pauli units.

The SHC inset preserves the outstanding KPM discrepancy: PdSe2 gives
2.52038425 for the saved projected M=1024, R=16 run versus
1.98650293 for the projected exact Gaussian oracle. The new longitudinal
curve is an exact finite-width Gaussian Kubo-Greenwood diagnostic,
not a converged DC material-gate result. No acceptance threshold,
golden or production source was changed for these figures.

## 2026-10-25 profiling and explicit closure hold

The user requested profiling to identify improvements after the material review
figures, and explicitly reserved project closure for a later authorization.
Project closure therefore remains on hold regardless of subsequent gate results.

Initial profiling is running sequentially on Aron against candidate-v5, using
the existing frozen Bi and PdSe2 reduced-material workloads: 24x24 mesh,
covariant velocity, Lowdin, shift window [-5,5] eV, M=256/1024, R=16,
seed=12345, original table reconstruction. Eight threads are pinned to CPUs
0-7 with BLAS limited to one thread. The wrapper records cProfile cumulative
and self times, phase timings, CPU time, peak RSS, affinity and source paths.
It does not alter production sources, goldens, numerical parameters or limits.

Remote artifacts: ~/lsquant-wrapup-20261024/profiling-20261025/.
Controller PID: 734194. Per-material child PIDs and exit files are recorded.
An initial missing optional threadpoolctl import was fixed in the diagnostic
wrapper; the original import failure is preserved. No production run occurred
in that failed attempt. Hardware counters are unavailable (perf is absent;
perf_event_paranoid=4).

Still required after these profiles: analyze the dominant functions; measure
unprofiled scaling on identical inputs and compare numerical outputs; prioritize
improvements by measured end-to-end impact. Profiled elapsed times must not be
reported as unbiased benchmarks or as an improvement. This initial pair does
not complete W08/W09 or authorize final project closure.

### Profiling results retrieved: Bi and PdSe2

Both sequential Aron profiles completed with exit 0. These are profiled totals
for M256 plus M1024 at R16, not unbiased benchmark or speedup measurements.

| Material | Total s | Projection moments s | SHC reconstruction s | Other s | Peak RSS GiB | CPU seconds / wall seconds |
|---|---:|---:|---:|---:|---:|---:|
| Bi | 304.753 | 207.748 | 85.793 | 11.212 | 2.791 | 7.258 |
| PdSe2 | 1070.146 | 859.114 | 84.509 | 126.523 | 17.004 | 6.817 |

Moments account for 68.2% and 80.3% respectively, making them the primary
optimization target. Table reconstruction accounts for 28.2% and 7.9%.
Within reconstruction, ContractRank2Imag takes 66.844/65.343 seconds and
legacy_delta_table takes 15.206/15.399 seconds (Bi/PdSe2). These nested
costs must not be added again to reconstruction totals. Preparation and
reference calculations are secondary on this workload. cProfile does not
resolve internal native engine kernels, so it does not yet establish whether
recursion, operator application or moment-table contraction dominates there.

No single-core utilization bottleneck is supported by these whole-run CPU
ratios. The receipt affinity [0] is the calling thread after OpenMP binding;
it is not evidence that all workers ran on CPU 0. Worker placement still
needs per-thread verification during unprofiled scaling.

Numerical comparison remains necessary: candidate-v5 M1024 reference-point
values are 1.1550715643600988 (Bi) and 2.3441496803075585 (PdSe2 upper gap),
whereas earlier saved reduced runs gave 1.1353409434012425 and
2.520384253866164. This is not established as floating-point reordering or
an optimization improvement; investigate source differences and trace-vector
identity before making before/after claims. No golden or threshold is changed.

Next priorities: native engine attribution, unprofiled matched-input scaling,
then a same-moment comparison of existing folded versus table reconstruction.
The project remains open pending the user's explicit closure authorization.

### Targeted integrated kernel scaling confirmation

Recompiled the existing bench/kernel_bench.cpp against the integrated native
library, using candidate-v5 headers. A recursive comparison confirmed that
candidate-v1 and candidate-v5 lsquant/core sources are identical; the existing
build-v1-native library is therefore applicable to these unchanged kernels.
Executable/library SHA256 and source comparison are saved remotely under
~/lsquant-wrapup-20261024/profiling-20261025/scaling/.

Sequential unprofiled measurements on Aron: 1/8/24/48/96 threads, taskset to
0 through T-1, OMP_PLACES=cores, close binding, dynamic teams disabled,
BLAS one thread. Each point is the median of three fresh processes, each
with one excluded warm-up and three timed operations. Setup is excluded.
The two workloads are the existing stage-8 sizes, not material-run timings.

| Native function | 1 core ms | 8 cores ms | 24 cores ms | 48 cores ms | 96 cores ms |
|---|---:|---:|---:|---:|---:|
| linalg::MomentTable, dimension 165888, batch 128 | 1835.867 | 290.141 | 147.750 | 133.418 | 174.241 |
| ContractRank2Imag, M256, 7680 points | 181.227 | 23.438 | 9.050 | 6.475 | 4.606 |

MomentTable accelerates 12.43x/13.76x/10.54x at 24/48/96 cores. It is
30.6% slower at 96 than 48 cores. ContractRank2Imag accelerates
20.03x/27.99x/39.35x and continues improving across the second socket.
All output hashes agree across all thread counts and repetitions for each
kernel (MomentTable 8221316592892497927, contraction 11479475180656984162).
This is a checksum comparison, not a new material or stage golden verdict.

Code inspection provides one concrete limit to investigate: MomentTable's
fixed 16x16 result tiles generate only 64 multiplication tasks at batch 128,
so that phase cannot occupy all 96 workers. Packing/copying, allocations,
NUMA and per-tile GEMM costs may also contribute; their separate times have
not been measured. Do not generalize this 64-task limit to M1024/batch1024.

Priorities: (1) native MomentTable scaling and packing/tile investigation;
(2) engine-internal attribution for the material workload; (3) legacy_delta_table
Python generation/caching and existing folded reconstruction comparison.
The previous material cProfile used the golden build, whose rank-2 engine
uses a different accumulation path. Its 68/80% engine shares must not be
attributed directly to native MomentTable. No source optimization was made
and no end-to-end speedup is claimed. Closure remains explicitly withheld.

### Consolidated profiling report

See [profiling-results.md](profiling-results.md) for measured costs, native
scaling, the batch-specific 64-task limit, conventional optimization options,
provenance and unresolved numerical comparisons. The measured plateau is not
an established hardware limit. No optimization or closure is claimed.

### SymPy mathematical and spatial-symmetry audit

User requested a mathematical audit and explicit spatial-symmetry relations
for moments. See [mathematical-symmetry-audit.md](mathematical-symmetry-audit.md).
Aron execution with SymPy 1.14.0 completed: 48 zero-residual records in
sympy-results.json, the existing metric trace identities and 12 exact spectral
projection identities, plus the spatial tensor/projection counterexamples in
sympy-spatial-results.json. Finite-degree corroborations are distinguished
from general symbolic trace, trigonometric and block identities.

Findings: folded reconstruction is already mathematically justified but was
not used in the profile; shift keeps 72/220 orbitals although 12/36-40 are
selected in Bi/PdSe2; clean projected traces admit selected-block compression
with consistent metric, mapped samples and normalization. Ordinary spatial
symmetry constrains operator components and does not by itself imply odd-order
zeros. Axial spin, antiunitary time reversal and chiral spectral reflection
have separate rules. Finite random samples need not obey exact-trace symmetry.
No actual space group was certified for either frozen material. No production
source, golden, seed, numerical tolerance or project closure status changed.
