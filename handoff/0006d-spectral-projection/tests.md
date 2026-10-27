# Stage 6d tester checkpoint — 2026-10-23

Tested delivery: `22965bfe40f75a6d479f25b906c5ccc069ea8827`, recovered from
`/root/lsquant-codex-6d` on Ionos. Source and tolerances were left unchanged.

Machine: `tcn-aronnode`; GCC 9.4.0, C++17 Release, Python 3.10.22,
NumPy 2.2.6, SciPy 1.15.2, sisl 0.16.4. Fresh configure and build passed.
The isolated checkout/build are under `~/lsquant-verify-20261023-1100/`.

## Failure

The existing `tests/dev/stage_6d_check.py`, with its default 12x12 mesh,
M=1024 and 10 meV broadening, exited 1 after 30.910 s. Environment:
`KPM_SEED=12345`, 48 OpenMP threads, core binding, one BLAS thread.
Maximum RSS: 1387672 KiB.

Every printed symbolic residual is zero. The first numerical failure is at
`tests/dev/stage_6d_check.py:147`: the unprojected Lowdin exact-trace KPM
baseline is outside the required 1e-3 plateau tolerance:

| Point | KPM | Unprojected oracle | Difference |
|---|---:|---:|---:|
| 1 | 2.016985042339149 | 1.9999237004848083 | 0.017061341854341 |
| 2 | 2.016965157027824 | 1.9999237004848083 | 0.017041456543016 |
| 3 | 2.0169409507273683 | 1.9999237004848083 | 0.017017250242560 |

The failure precedes the hard/shift comparison, the generalised route and
the orthogonal/off-grid checks. Those assertions were not reached. This is
evidence against convergence of the prescribed unprojected baseline at these
settings; it does not establish a defect in the projection operators. The
developer report explicitly identified this baseline-convergence risk.

Full log: `figures/out/verification_20261023/stage6d-selfcheck.log` locally;
`~/lsquant-verify-20261023-1100/aron-results/stage6d-selfcheck.log` on Aron.

## Disposition

The original controller withheld the 6d material sweep because its required
self-check failed. The user subsequently requested an independent projection
comparison, first in k space and then real space, and explicitly selected the
saved stage-6b exact SIESTA references. That study is authorised independently
of the failed self-check; it does not close the stage or relax its gate.

## User-requested projection study: completed 24x24 and real-space checkpoint

The previous Aron and Jose production queues were suspended by their recorded
PIDs. Receipts are `aron-results/projection-pivot-stop.json` and
`fat-results/projection-pivot-stop.json` under the isolated verification root.
They remain paused. The reference-matched study runs separately on Aron in
`~/lsquant-verify-20261023-1100/projection-reference-study/`.

Harness: `figures/stage_6d_projection_comparison.py`. References are the four
`~/t6b2/pd/deg_{bi,pdse2}_{berry_connection,covariant}.log` files, copied with
SHA-256 receipts. These contain values printed to five decimal places, not
full-precision byte-comparable curve goldens. The requested degeneracy tolerance
is 1e-6 eV, broadening is zero, and each comparison matches velocity mode and
k mesh. Baseline drift from the saved reference is reported separately from
the exact projected-minus-unprojected difference. The numerical convention is
the oracle's `e^2/h, sigma_z Pauli`; no prefactor was adjusted.

The current loader initially recertified dipoles under a newer policy. Those
preliminary results are preserved in `projection-study/`, and only that
study's recorded PIDs were terminated when this provenance difference was
identified. They are superseded by the comparison below. The old stage-6b
loader at commit `99cbda7b9fdb2cb777e87ae3d11da5cfd55755ec` was used to validate
and reuse the exact dipole caches named in the saved logs (quadrature tolerance
5e-5, maximum order 384). Its input key, checksum and certificate checks ran
before the physical H, S, A and geometry arrays were exported. The new
projection code consumes those frozen arrays; no cache or certificate was
rewritten. `study_provenance.json` records full-file and per-array SHA-256s.
The harness's `--phase freeze`, run with that reference checkout/build on
PYTHONPATH, reproduces the export and refuses to overwrite an existing one.

Windows are [-3,3] and [-5,5] eV in the existing model energy convention,
matching the delivered projection run files. They are not silently recentered
on the PdSe2 gap probe at 1.052 eV. Bi is probed at 0, -0.05 and 0.1 eV;
PdSe2 at -0.3844 and 1.052 eV, exactly as in the saved logs.

Completed exact shift results, at E=0 for Bi and E=1.052 for PdSe2:

| Material / velocity | Mesh | Saved | Unprojected | Shift +/-3 | Error | Shift +/-5 | Error |
|---|---|---:|---:|---:|---:|---:|---:|
| Bi / covariant | 24x24 | 1.14415 | 1.14415315 | 1.28938803 | +12.69% | 1.07975947 | -5.63% |
| Bi / berry_connection | 24x24 | 1.51609 | 1.51608670 | 1.39798944 | -7.79% | 1.28960785 | -14.94% |
| PdSe2 / covariant | 24x24 | 1.97911 | 1.97910576 | 1.96773956 | -0.574% | 1.98650999 | +0.374% |
| PdSe2 / berry_connection | 24x24 | 2.08124 | 2.08123511 | 1.95809886 | -5.92% | 1.97684882 | -5.02% |
| Bi / covariant | 48x48 | 1.17473 | 1.17472608 | 1.31794294 | +12.19% | 1.11033507 | -5.48% |
| Bi / berry_connection | 48x48 | 1.55269 | 1.55268433 | 1.43496931 | -7.58% | 1.32639323 | -14.57% |

Errors in the table are against the fresh unprojected oracle on identical
inputs. All 24x24 baselines round to the saved values. Bi berry_connection
at 48x48 differs by -5.67e-6 from the rounded reference, one final printed
digit; byte identity is not claimed. The saved oracle excluded pairs by an
energy-difference threshold, while the current oracle groups transitive
degenerate subspaces. That baseline difference is retained in the receipts
and is much smaller than the measured projection changes.

Hard projection is nearly identical to shift at these probes. PdSe2's
second saved probe, E=-0.3844 eV, is also important: the covariant reference
is +0.01051, but the +/-3 and +/-5 projections give -0.00968587 and
-0.00113794, respectively. Matching the gap plateau alone is insufficient.

These are exact-oracle errors before KPM truncation.
The current implementation projects both currents as well as H, removing
virtual transitions involving discarded states. Preserving the retained
eigenpairs alone does not certify preservation of the Hall response.

Real-space checks ran after both materials' 24x24 comparisons. They export actual
CSR H, S, JXSZ and VY on a common 4x4 periodic mesh without radius truncation,
compare every matrix entry against an independent NumPy discrete Fourier
construction and their action under the atom-gauge C++ FFT, and compare raw M=64 moment
tensors using two identical transformed random vectors. The original S is
retained; its inverse square root is applied by an exact FFT composite.
This checks finite-system representation equivalence, not convergence of a
4x4 material result to the 24x24/48x48 references. All eight material/kind/window
cases pass: the largest full-matrix relative Frobenius difference is
1.51e-16, and the largest raw-moment relative Frobenius difference is
5.69e-14. This supports the Fourier/CSR implementation; it does not restore
the Hall contributions removed by projection.

Local evidence: `figures/out/verification_20261023/projection-reference-study/`.
`projection-vs-golden.png` overlays the completed exact projected values and
saved reference plateaus, with percentage residuals. Its adjacent caption
states the meshes and limits. The JSON receipts and real-space moment NPZs
are included; frozen material arrays remain on Aron. The figure was rendered
on Aron and visually inspected locally.

Still running: both PdSe2 48x48 velocity cases, with PIDs recorded in the
study directory. `finish-study.py` waits for them, validates that all eight
material/mode/mesh records contain four projections, checks the completed
real-space receipts, and refreshes the remote figure. The local snapshot
explicitly marks these larger-mesh results pending.

The full M=256/1024, 16-vector material gate, finite-radius convergence and
off-grid validation remain unclosed. The preliminary KPM exploration was
stopped during the provenance restart and is not used as validation. No
source, acceptance tolerance, reference or golden was changed. No stage
closure is claimed.

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

The earlier reference to included JSON/NPZ receipts refers to the original
working snapshot. The checkpoint Git tree contains the figure, caption and
provenance record; raw local JSON receipts are in the archive above. Frozen
arrays and remote-only NPZ receipts remain at their recorded remote locations.

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
