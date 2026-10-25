# Authorized integration implementation

2026-10-24. The user explicitly authorized Codex to integrate the four
deliveries, fix source, and validate, overriding the tester-only restriction.
No commit or push has been made. This is an implementation report, not a
stage closure.

## Candidate and provenance

Worktree: codex/lsquant-wrapup, based on 6b commit
2c790f19320ad41abfee3043b00ee379f2fe35ac. Applied delivery diffs:
6c b2ab659af4bba4e64fecded4c3d579c28023b8c3,
6d 22965bfe40f75a6d479f25b906c5ccc069ea8827,
8 10d02e43309d40292afc656834ca09975ef4b012.

Candidate v5 source tree: bf70c9b8c8b8ec66df329265d6c170a1cb315e0e.
This is a Git tree, not a commit. Reports written afterward are outside
that identifier. The LF archive reproduces this tree on Aron.
Remote source: ~/lsquant-wrapup-20261024/candidate-v5.

Integration retains both 6b/6c run-file parameter sets, reference and
projection fields, relative Hermiticity checks, native first-touch
allocation, and fixed-order golden accumulation. Projection dispatch
precedes disorder dispatch to prevent silently ignoring a projection.
The implementation still refuses projected disorder until its composition
has been implemented and validated.

## Corrections

1. FFT Bastin contraction accessed Python complex.imag as a function.
   The property access fixes the crash. FFT/table results are numerically
   close, but this does not satisfy the required byte gate.
2. Factoring orbital Bloch phases multiplied diagonal entries by a phase
   and its conjugate. Their rounded product perturbed exact identity
   overlap, disabling the stage-5 orthogonal-limit path. Preserve diagonal
   entries before phase dressing and restore them afterward; compute
   tau_j-tau_i before adding R. The existing 6c byte assertion now passes.
3. Large uniform-grid Bastin reconstruction allocated and filled O(M*N)
   coefficient tables and contracted them at O(M*M*N) cost. The opt-in
   folded backend uses
   cos(m theta) exp(-i n theta) =
   (exp(-i(n-m)theta)+exp(-i(n+m)theta))/2.
   Fold moment entries onto frequency diagonals, pair positive/negative
   frequencies, and evaluate the resulting real Chebyshev series on the
   original uniform nodes. Work is O(M*M+M*N), with O(M+N) scratch beyond
   the moments and existing damped moment array. Broadening truncation,
   Jackson coefficients, integration, prefactors, and energy grid stay
   the same. Floating-point operation ordering changes.
4. The local projected matrices lacked a production real-space adapter.
   ProjectedRealSpaceBasis applies the existing LocalMatrix CSR exports
   and unchanged original overlap through exact or Chebyshev metric
   appliers. It reuses the 6c mapped random stream and metric trace.
   Velocities and currents are projected physical operators, never
   differentiated from the projected Hamiltonian. Existing checked_csr
   validates relative Hermiticity before folding accepted round-off.

The folded backend is available through kubo_bastin/shc backend="folded"
and SIESTA run-file reconstruction.bastin_backend="folded", including
projection and disorder executors. It is saved in resolved run files.
The default remains tables; unsupported backend names and a folded
FFT-grid request are refused. Non-SIESTA run files cannot opt in yet.

## Build and execution

Jose golden and Aron golden/native builds succeed. On Aron, select Eigen
from ~/lsq-env explicitly with CMAKE_DISABLE_FIND_PACKAGE_Eigen3=ON and
EIGEN3_INCLUDE_DIR=$HOME/lsq-env/include/eigen3; the system Eigen 3.3
cannot compile an existing Eigen::Vector alias. Build logs, including the
initial dependency failure, remain under the remote logs directory.
No C++ code changed after these successful builds.

The first source archives inherited Windows CRLF from Git attributes.
Those historical-gate attempts failed during shell setup, before a gate
verdict. Archives now use git -c core.autocrlf=false archive; the unpacked
index tree is checked against the local tree. Both failed packaging logs
remain available.

Run the existing developer check with:

    python tests/dev/stage_6d_check.py --moments 8192 --backend folded

It passes every original assertion on v3. Both original and projected
calculations use M=8192 and exact trace. The M=1024 invocation remains
underresolved for the original wide-spectrum model; its failure was not
hidden or relaxed. Candidate v5 repeats 6b, 6c and this 6d invocation;
all three complete with exit 0 and every original assertion passes.

## Validation and limits

Detailed receipts and numerical comparisons are in tests.md. Candidate
v5's local figure harness passes 16 API comparisons and eight resolved
run-file replays (selected DOS and SHC), with folded/table comparison on
the same grid. Source and numerical inputs are preserved remotely.

Radius 0 and 0.8 Angstrom on the toy remove every physical hopping from
the hard-projected Hamiltonian, leaving numerical cancellation residue
near 5.5e-14. Relative Hermiticity validation refuses this residue. That
limitation is retained; no tolerance was changed. The finite-cutoff
comparison uses 1.2 Angstrom, retaining the 1 Angstrom nearest-neighbour
bonds and omitting longer bonds. Its altered DOS is an explicitly
reported truncation error, not an accepted physical approximation.

Open work: projected disorder composition; production material and
radius/off-grid validation; remaining projection matrix and orthogonal
tool scope; complete FFT/golden verdict; full-workload performance;
historical regression on the final candidate; deferred full-model long
validation and written review. The project is not closed.
