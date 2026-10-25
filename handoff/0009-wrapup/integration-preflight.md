# Integration preflight and proposed corrections

Initial tester checkpoint, 2026-10-24. The proposals below are historical:
the user subsequently authorized source implementation, and candidate v1
applies the conflict resolutions and FFT correction. See tests.md for the
source-tree identifier, builds and ongoing validation.

## Available deliveries

Remote heads still match the four revisions recorded in plan.md. There is
no combined candidate on origin/cpp-unified. Local working-tree files were
not used as an integration target.

Read-only git merge-tree comparisons (the tree objects are diagnostic,
not commits or buildable candidates):

| Pair | Diagnostic tree | Conflicted paths |
|---|---|---|
| 6b + 6c | 8dca642a81bc075e0b007045e8cdfec74554a626 | python/lsquant/runfile.py |
| 6b + 6d | b6504b4387d4b70a02aab5f9bddbd8a2715c0882 | python/lsquant/runfile.py |
| 6b + 8 | e3a7e6c4d455675f0efbb63f0ba0a221e472f288 | runfile.py, block_diagonal_operator.cpp, moment_engine.cpp |

These pairwise checks do not prove a combined four-branch merge is clean.

## Required conflict resolutions

1. Preserve both bloch_memory_mb and nonorthogonal_disorder parameters and
   their validation. Keep the reference option and projection configuration
   together in validation, RunFile fields, parsing and serialization.
   Audit dispatch so projected disorder cannot silently take a clean path.
2. Keep the 6b relative Hermiticity/finite-value checks in
   BlockDiagonalOperator while retaining native first-touch placement.
   Choosing the older stage-8 boolean/absolute test would revert a 6b fix.
3. Preserve the native MomentTable path. In the golden fallback retain
   the 6b static outer-loop parallelization and fixed per-element sum order.
   Do not apply both accumulation paths to the same tensor.
4. Audit auto-merged binding declarations, reconstruction, metric and
   sampling paths as well as explicit conflicts. Require the existing
   self-checks and same-machine thread comparison after the merge.

## FFT defect with a concrete correction

Observed in stage-8 commit 10d02e43309d40292afc656834ca09975ef4b012,
python/lsquant/reconstruct.py:578, inside kubo_fft_grid, backend="fft":

    out[i] = green.imag() / sine / math.pi

green is a Python complex number. Its imag attribute is a float, so every
invocation raises TypeError: 'float' object is not callable. Proposed fix:

    out[i] = green.imag / sine / math.pi

Reproduced through the existing bench/compare_contraction.py on G4 and G5,
Greenwood and Bastin, golden and native profiles, 1 and 48 threads:
16/16 processes fail at this same line. The table backend runs before the
failure. No FFT numerical residual or speedup can be claimed from these
failed attempts.

After correction, rerun those exact comparisons and inspect identical
energy grids, absolute/relative residuals and thread-count hashes.
Do not assume the one-character correction proves the FFT arithmetic.

## Existing 6d assertion is under-resolved

The default self-check was rerun unchanged and failed at line 147.
Independent dense eigenbasis moments agree with the engine's leading
64x64 block to 6.74e-16 maximum absolute error. Independently generated
M=1024 moments reproduce the failed plateau to about 1e-14.
Doubling integration-grid density changes the worst error by about 1.2e-7.

Wide-spectrum exact-trace errors at M=512,1024,2048 are
0.0682451, 0.0170613, 0.00429533. The shifted model's error at M=1024 is
0.000366605, below the existing 0.001 limit. The wide window spans
75.6 eV, so M=1024 has nominal resolution about 116 meV despite the
10 meV reconstruction request.

Required developer/reviewer decision: make the existing self-check
distinguish finite-M resolution gain from converged exact agreement.
Keep the 0.001 criterion; provide sufficient order for a converged wide
baseline, or use an independently verified finite-M reference for the
engine comparison while retaining a separate converged-oracle assertion.
The tester has not changed moments, assertions or tolerances in the
delivered self-check.

The orthogonal toy's existing training/off-grid calculations were also
executed independently of the failing assertion. Maximum errors:
4.0e-15 eV on the training mesh and 2.62e-14 eV on 40 off-grid points.
This does not establish material or finite-radius validity.
