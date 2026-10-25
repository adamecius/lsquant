# Stage 9: wrap-up and final validation

Prepared by the tester on 2026-10-24 at the user's request.
Status: execution started; see tests.md for diagnosis and integration findings.

This stage collects the remaining work from 6b, 6c, 6d and 8, then validates
one integrated revision. The original wrap-up adds no new physics scope;
the user-requested symmetry-input extension is tracked separately as W12 below.
Stages 0-7 have a
historically closed baseline; their current-delivery regression is unfinished.
A completed process is not a numerical pass, and a passed self-check is not
a stage closure.

The user authorized Codex to implement source corrections as well as test
them on 2026-10-24, explicitly overriding the tester-only role restriction.
The existing briefs, accuracy criteria and golden provenance remain the
contracts. Implementation and validation evidence are recorded separately;
no new unit-test suite is required.

## 1. Starting evidence

| Area | Evidence available | What it does not establish |
|---|---|---|
| 6b SIESTA | Build and existing developer self-check passed, including the 20-seed sampler study | Completed material gate, all review items, or historical byte regression |
| 6c disorder | Existing self-check passed; finite-system real/Bloch agreement and orthogonal limit checked | Complete Bi W sweep, all metric routes, or integrated projection with disorder |
| 6d projection | All 16 retained-band checks passed the user's MAE <10 meV criterion; 32 exact projection cases completed; eight untruncated 4x4 real-space comparisons passed | Transport preservation, finite-radius/off-grid validity, orthogonal downfolding, or production real-space closure |
| Reduced KPM | 16 calculations finished in about 1 h 57 min on Aron, shift +/-5 eV, W=0 | Numerical acceptance: at M=1024, projected-oracle plateau errors are 3.7-7.9% for Bi and 10.3-27.2% for PdSe2 |
| 8 performance | Both builds passed; 120 kernel measurements completed, hashes stable across thread counts within each profile | Full-workload scaling, cross-profile accuracy, FFT acceptance, or completed golden gate |
| Historical regression | Partial 16-thread sweep completed on Jose | Final byte verdict: paused during tutorial06; complete 16/1-thread gate unfinished |

The independent 6d self-check fails before its projected KPM assertions:
the unprojected Lowdin exact-trace M=1024 baseline differs by about 0.017
from its oracle against the existing 0.001 limit. It must be investigated
separately from stochastic errors in the material runs.

The retained-band criterion accepts the reduced model for test execution.
It does not re-bless spin Hall goldens. Projection also changes the exact
transport observable by removing transitions through discarded states;
this model error must be reported separately from KPM error.

The original full-material/disorder queues and Jose regression were paused
during the projection pivot. Their partial receipts remain historical
evidence. The final run below uses fresh outputs on the integrated revision,
rather than treating those partial runs as a final certificate.

Evidence sources:
- [6b tester report](../0006b-siesta-spin-hall/tests.md)
- [6c tester report](../0006c-nonorth-disorder/tests.md)
- [6d tester report](../0006d-spectral-projection/tests.md)
- [8 tester report and scaling measurements](../0008-performance/tests.md)
- Reduced numerical receipt:
  `figures/out/verification_20261023/reduced-validation/completion-snapshot/summary.json`
- Exact projection comparison and figure:
  `figures/out/verification_20261023/projection-reference-study/`

## 2. Pending work, with closure evidence

| ID | Stage | Required action | Owner | Completion evidence |
|---|---|---|---|---|
| W01 | Integration | Combine the reviewed 6b improvements, 6c, 6d and 8 into one candidate; resolve overlapping changes in metric, sampler, operators and reconstruction | Developer and Fable | Candidate commit, component ancestry, conflict decisions, complete developer reports |
| W02 | 6d | Localise and fix the failed unprojected exact-trace baseline before interpreting projected KPM results | Tester diagnosis; developer fix | Existing self-check passes unchanged, including assertions previously not reached |
| W03 | 6b/6d | Explain reduced-material SHC residuals, including PdSe2 getting worse from M=256 to 1024 | Tester diagnosis; developer fix if needed | Comparison separating projection bias, finite-M/kernel error, stochastic error and implementation error; applicable existing criteria met |
| W04 | 6b | Verify review 1-4 changes on the integrated candidate and finish the SIESTA material gate | Tester | Velocity/metric comparisons, sampler receipts, dipole certificates/cache provenance, degeneracy and chunking diagnostics, oracle and KPM outputs |
| W05 | 6c | Complete Bi disorder validation in Bloch and real space, both metric routes | Tester | DOS/SHC at W=0, 0.1, 0.3 eV, 24x24; same mapped random phases and disorder; residuals and W-sweep figure |
| W06 | 6d/6c | Complete the local projected model and its integration with disorder and run files | Developer; tester validates | Consistent projected operators, unchanged overlap, CSR/run-file equivalence, radius convergence and off-grid errors; projected disorder supported and compared between spaces |
| W07 | 6d | Finish original projection scope beyond the chosen shift +/-5 eV subset | Tester | Full hard/shift, +/-3 and +/-5 eV, 24/48 mesh, Bi/PdSe2 comparison; orthogonal tool provenance and validation |
| W08 | 8 | Validate native/golden and FFT/table numerical behaviour on the integrated code | Tester | Thread comparisons, same-grid FFT residuals and golden verdicts, profile provenance, no hidden tolerance change |
| W09 | 8 | Finish performance measurements and resolve or document regressions | Tester; developer if changes needed | Full-workload timings, 24/48/96-core and per-socket measurements, crossover/batch study, before/after evidence, stage figure |
| W10 | 0-8 | Run the complete final regression and long full-model validation on one frozen revision | Tester | Per-case verdicts, logs, figures and complete inventory; no unexplained missing or skipped case |
| W11 | Closure | Reconcile reports and review all residual findings | Fable with developer/tester evidence | Written stage verdicts and final validation record tied to the tested revision |
| W12 | User-requested input extension | Add global spatial-symmetry declarations and an independent time-reversal section to the run input; see section 8 | Developer; tester validates | Documented schema and transformations, model/sampler validation, unchanged behavior when omitted, comparison against the unreduced calculation |

W06 is needed to use the reduced model for the intended disorder/real-space
workflow. The present executor refuses nonzero disorder and requires the
CSR/API path for real-space work. Passing a small API comparison does not
close this interface or the radius-truncated model.

Do not rebuild velocities by differentiating the projected Hamiltonian as
a substitute for consistently projecting the physical operators. That would
change the observable under test.

For W07, the original KPM matrix contains 64 configurations:
2 materials x 2 velocity modes x 2 windows x 2 projection kinds x 2 meshes
x 2 moment counts. The completed 16 configurations cover only shift +/-5.
The existing 32 exact cases cover the full corresponding oracle matrix.
Reuse their provenance where unchanged; repeat affected comparisons on the
final revision. The user's choice of a fast default does not silently remove
the other brief requirements.

The orthogonal tool must retain origin provenance, be callable, pass the
existing toy check at 1e-6, and have the specified off-grid band comparison.
Its status is unverified in this campaign, rather than assumed failed.

## 3. Acceptance rules

1. Preserve committed goldens, seeds, normalisations and prescribed tolerances.
   A numerical reorder is reported with the first difference and maximum
   absolute/relative error. Any acceptance under the project's round-off
   exception needs Fable's written verdict.
2. Keep the user's retained-band MAE threshold strictly below 10 meV,
   evaluated separately by material/mesh/window/kind. Discarded states are
   excluded explicitly. Record maximum error too.
3. Report both SHC differences: KPM minus the projected exact oracle, and
   projected exact minus the saved unprojected SIESTA reference. Preserve
   velocity mode, mesh, energy, broadening, units and degeneracy policy.
4. Use the saved stage-6b SIESTA references selected by the user. Their
   printed precision is not a byte-level oracle. Preserve input checksums,
   frozen matrices, quadrature certificates and cache keys; make any new
   certification policy a separately reported comparison.
5. Use existing quantitative gates where specified, including the 6d
   self-check's 1e-3 limit. Where a brief requests an error report without a
   numerical acceptance limit, record that explicitly and obtain a written
   review verdict; do not invent a pass threshold or infer one from meV.
6. Keep spin Chern/Z2 checks distinct from SHC. Do not force real-material
   SHC to an integer or substitute one observable for another. Apply the
   later 6b clarifications: the SIESTA route stands on its own and Wannier
   comparisons are optional where valid inputs exist.
7. Byte verdicts against the historical goldens belong on Jose/fatnode-tacn.
   Aron comparisons are machine-local and profile-local. Matching kernel
   hashes do not prove cross-profile accuracy or complete-run equality.

## 4. Execution order for the final rerun

### A. Resolve blockers and freeze the candidate

Reproduce W02 using its existing self-check. Use the existing material
harness and saved tensors for W03 where possible before spending more
large-model time. Record the first failing quantity and separate exact
trace, stochastic sampling and reconstruction effects.

Integrate W01, implement the missing requirements, and rerun the affected
existing checks after each fix. Once those pass, freeze one candidate SHA.
Record source/diff state, toolchains, dependencies, run-file and input hashes,
build flags and reference provenance. Use isolated builds and output roots.
Changes after this freeze require affected checks to be repeated and a final
full regression on the resulting candidate.

The separate delivered commits tested so far are starting evidence only:
6b `2c790f19320ad41abfee3043b00ee379f2fe35ac`;
6c `b2ab659af4bba4e64fecded4c3d579c28023b8c3`;
6d `22965bfe40f75a6d479f25b906c5ccc069ea8827`;
8 `10d02e43309d40292afc656834ca09975ef4b012`.
They do not collectively certify an integrated build.

### B. Run the historical regression first

On Jose, fresh golden build, `KPM_SEED=12345`:

- Run the existing `scripts/check_goldens.sh 7 <build>` for the complete
  historical G1-G10/tutorial inventory at 16 and 1 threads, including
  tutorial06, CLI/MCP/refusal checks and structural checks.
- Use an explicit historical manifest if newer G11 runs make the default
  all-files selection include material workloads. Record every exclusion
  and its dedicated material-gate destination; no reduced quick subset
  counts as the final historical gate.
- Run `scripts/check_goldens.sh 6only <build>` for the stage-6-specific
  assertions not replaced by running its JSON files.
- Run the existing prescribed developer self-checks, including 6b/6c/6d.
  They supplement the golden gate, not replace it.
- Reproduce the historical figures and residuals from the same candidate.
  Investigate any historical expected failure against the current code;
  old exceptions are not carried forward automatically.

Record per-file byte verdicts and localisation for each mismatch. Stop
dependent expensive work on a real failure; independent diagnostics may
continue.

### C. Validate reduced materials, k-space then real-space

On Aron, run W04-W07 with explicit resolved configurations:

- Clean k-space first: Bi and PdSe2, matching saved reference modes/meshes,
  full projection matrix, M=256/1024, R=16, seed 12345 and 10 meV broadening
  where prescribed. Report band and transport errors separately.
- Compare Lowdin and generalised routes and the production Bloch path with
  the eigenbasis reference using the prescribed common sampling.
- Treat anticommutator cases according to the existing symbolic validity
  conditions: diagnostics or expected refusals where equivalence is not
  established, never silently promote them to exact routes.
- Real space second: reproduce the exact finite-mesh FFT/CSR comparison,
  then test finite radius and off-grid behaviour. A 4x4 no-cutoff pass
  does not substitute for those checks or the 24/48 production meshes.
- Complete the 6c Bi 24x24 W sweep, both spaces and metric routes, with
  matched disorder realisations and mapped phases. Validate the integrated
  reduced-disorder path before using it for accelerated production.

The fast reduced campaign must meet numerical acceptance, not just return
exit code zero, before the deferred multi-day full-model work starts.

### D. Complete performance and reproducibility

On the same candidate, run the stage-8 prescribed workloads: graphene
2000x2000 DOS M=2048, Haldane 48x48 rank-2 M=256, and Bi block SHC with
the specified material settings. Use golden and native profiles,
24/48/96 cores, three repetitions, and one-process-per-socket comparisons.
Keep before/after input and reference settings identical.

Record total time, preparation/certification, moments, reconstruction and
reference costs separately. Include median/range, CPU occupancy, memory,
affinity/NUMA placement and output comparisons. Report modeled throughput
as modeled, not measured DRAM traffic. Archive unavailable counters as
unavailable. Audit:

- moment-table scaling flattening at 24-48 cores and slowing at 96;
- native single-core sparse-batch regression against golden;
- the sharp sparse-kernel change between 48 and 96 cores;
- measured sparse parallel crossover and batch-size settings;
- FFT versus table/GEMM contraction on identical nodes for G4/G5 and
  prescribed workloads, with residuals and timing.

The saved 120 kernel runs remain a baseline. Repeat affected measurements
after integration/changes, rather than rerunning unchanged diagnostics
without a reason. The optional oneAPI study is not a closure blocker.

### E. Run the deferred full-model validation

After B-C pass, run the original unprojected SIESTA material matrix with
the validated implementation: Bi/PdSe2, bare/berry_connection/covariant,
Lowdin/generalised, original prescribed settings and explicit reference
policy. Include any remaining original full-model 6c requirements.

Keep the user-requested long validation separate from the projected tests.
Record time and numerical comparisons against both original references
and the reduced results. Do not substitute narrower spectra, fewer states,
or relaxed convergence settings to make this final run finish faster.

The long-run duration is unknown until the short integrated profiles are
available. Do not promise a completion date from kernel speedups alone.

### F. Review and close

Publish the final per-case inventory, all unresolved differences, stage
figures/captions, provenance and exact candidate SHA. Reconcile developer
reports with tester results. Fable issues the closure verdict only when
all required items have evidence and any deviations are explicitly resolved.

No source commit, push, release or project closure is implied by preparing
this plan. Follow the user's commit/push instructions when that step is
requested.

## 5. Final deliverable and completion checklist

The final tester report belongs at `handoff/0009-wrapup/tests.md`. It must
include a row for every required case: completed/pass/fail/blocked,
configuration, source/build SHA, reference, result path, comparison metric,
limit or review verdict, first difference, timing and figure.

One final overview figure should show material oracle/KPM overlays and
residuals, with reduced and full models clearly labelled, plus the reviewed
performance comparison. It supplements the required per-stage figures and
uses saved results; it is generated after the numerical rerun.

- [ ] Integrated candidate and developer reports complete (W01).
- [x] Existing 6d self-check failure diagnosed and passes at equal M=8192 with folded reconstruction (W02); M=1024 remains underresolved.
- [ ] Reduced transport residuals explained and accepted (W03).
- [ ] 6b material/review requirements complete (W04).
- [ ] 6c disorder matrix complete (W05).
- [ ] Local projection, disorder and run-file integration complete (W06).
- [ ] Full 6d projection and orthogonal-tool requirements complete (W07).
- [ ] Profile/FFT numerical validation complete (W08).
- [ ] Stage-8 performance evidence and figure complete (W09).
- [ ] Historical regression and long full-model validation pass (W10).
- [ ] Final figures, provenance, per-case report and written review complete (W11).
- [ ] Global spatial-symmetry input and separate time-reversal input implemented and validated (W12; TODO only at present).

## 6. Scope sources

- [AGENTS.md](../../AGENTS.md), especially verification and roles.
- [6b brief](../0006b-siesta-spin-hall/brief.md) and reviews 1-4.
- [6c brief](../0006c-nonorth-disorder/brief.md).
- [6d brief](../0006d-spectral-projection/brief.md).
- [8 brief](../0008-performance/brief.md).
- [Historical final-validation procedure](../0007-words-cli-mcp/final-validation.md),
  used for coverage, not copied as an acceptance of its old exceptions.
- User instructions in this chat: saved exact SIESTA references, retained
  band MAE <10 meV, reduced validation first, multi-day full-model test
  afterward, and this wrap-up stage before final closure.

## 7. User-directed profiling and closure hold (2026-10-25)

After the material review figures, profile the integrated reduced Bi/PdSe2
workloads, identify preparation/projection, moments and reconstruction costs,
and validate candidate improvements using unprofiled timings and unchanged
numerical comparisons. Preserve the W08/W09 scaling requirements above.

The user explicitly requires a further authorization before project closure.
Passing tests or completing profiling does not itself authorize closure.

## 8. TODO: global symmetry input (user-requested W12)

Status: recorded for future implementation. This change updates the TODO only;
the following draft fields are not yet accepted by the run-file parser.
Mathematical basis: [mathematical-symmetry-audit.md](mathematical-symmetry-audit.md).

Add one global symmetry configuration for the calculation, with two distinct
parts. Spatial symmetries must not contain time reversal as a spatial entry.

### Spatial symmetries

- Inversion about a specified origin.
- Rotation by 90 degrees about a specified axis and origin; expose the common
  z-axis case explicitly, without assuming every material has that symmetry.
- Mirror X, defined as x -> -x relative to the specified origin (normal along x).
- Mirror Y, defined as y -> -y relative to the specified origin (normal along y).

Use explicit Cartesian axes and origin units in the documented input contract.
Resolve orbital/atom mappings and spin transformations from the model and the
declared operation; a name alone is not sufficient to certify a symmetry.
Treat supplied operations as generators and account for their group closure
when building equivalent moment components or k-point orbits.

### Time reversal

Provide an independent `time_reversal` section. Disabled/omitted means that
no time-reversal reduction is requested; it does not assert that the model
breaks time reversal. The enabled case must use the model's antiunitary
representation, including spin where applicable. It must not be treated as
another real-space reflection or as ordinary complex conjugation in an
arbitrary orbital/spin basis.

### Draft input shape

This illustrates the requested interface, not a certified Bi/PdSe2 input:

```json
{
  "symmetries": {
    "spatial": {
      "frame": "cartesian",
      "origin": [0.0, 0.0, 0.0],
      "origin_units": "angstrom",
      "operations": [
        {"type": "inversion"},
        {"type": "rotation", "axis": [0.0, 0.0, 1.0], "angle_degrees": 90},
        {"type": "mirror_x"},
        {"type": "mirror_y"}
      ]
    },
    "time_reversal": {"enabled": false}
  }
}
```

Defaults: no declared spatial operations and no time-reversal reduction.
Existing run files must retain their behavior when `symmetries` is absent.
Finalize field names during implementation; document the coordinate and
spin conventions together with a minimal example for each operation.

### Implementation and acceptance TODO

- [ ] Parse, validate and preserve both sections in resolved run metadata.
- [ ] Build metric-aware model representations and verify covariance of H, S,
  each physical observable, projection and the actual mesh/boundaries/cutoff.
- [ ] Reject invalid or unsupported declarations clearly before applying
  reductions; never modify the supplied Hamiltonian to force a symmetry.
- [ ] Generate moment identities with the correct polar velocity / axial spin
  transformations and separate antiunitary conjugation rules.
- [ ] Apply only validated reductions: forbidden tensor components, equivalent
  components and, where supported, irreducible k-point orbits with tensor
  transformations and proper weights.
- [ ] Define mapped or symmetry-adapted stochastic states explicitly; distinguish
  exact-trace identities from finite-sample estimates. Preserve the reference
  path for comparison and record estimator changes.
- [ ] Require symmetry of the actual disorder realization for per-run reduction;
  distribution/ensemble symmetry alone must not be used to reduce that run.
- [ ] Validate with the existing appropriate harness/oracles and record residuals
  and timings against unreduced inputs; preserve seeds, units and golden rules.

No production implementation, new simulation or closure is authorized merely
by recording this TODO. Final closure continues to require the user's explicit
authorization.
