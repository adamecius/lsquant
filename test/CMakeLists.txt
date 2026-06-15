# Golden-master harness for LinQT.
#
# Goldens (committed under test/golden/) are produced by the reference pipeline
# external wannier2sparse -> inline_compute-kpm-nonEqOp -> inline_kuboBastinFromChebmom,
# pinned by KPM_SEED (see scripts/regenerate_goldens.sh). The human refactors src/
# and diffs the regenerated outputs against the committed goldens.

set(GOLDEN_DIR    ${CMAKE_CURRENT_SOURCE_DIR}/golden)
set(REGEN_SCRIPT  ${CMAKE_SOURCE_DIR}/scripts/regenerate_goldens.sh)
set(GOLDEN_STAMP  ${GOLDEN_DIR}/operators/graphene_golden.HAM.CSR)

# On-demand: `cmake --build <dir> --target regenerate-goldens`
add_custom_target(regenerate-goldens
    COMMAND bash ${REGEN_SCRIPT} ${CMAKE_BINARY_DIR}
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    COMMENT "Regenerating golden masters (KPM_SEED-pinned)")
add_dependencies(regenerate-goldens
    inline_compute-kpm-nonEqOp inline_kuboBastinFromChebmom)

# First-build bootstrap: generate goldens as part of ALL when they are absent
# (or forced via -DLINQT_REGEN_GOLDENS=ON) AND the external input tool is enabled.
# The EXISTS check is evaluated at configure time, so committed goldens make this
# a no-op on subsequent configures (idempotent).
if(LINQT_BUILD_INPUT_TOOLS AND (LINQT_REGEN_GOLDENS OR NOT EXISTS ${GOLDEN_STAMP}))
    add_custom_target(bootstrap-goldens ALL
        COMMAND bash ${REGEN_SCRIPT} ${CMAKE_BINARY_DIR}
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
        COMMENT "Bootstrapping golden masters on first build")
    add_dependencies(bootstrap-goldens
        inline_compute-kpm-nonEqOp inline_kuboBastinFromChebmom)
endif()

# --- chain1d ANALYTIC goldens vs the closed form (thesis Eq. 4.9 / B.15) ---------
# Native CTest (already include()'d at the root) -- not a second framework. The C++
# comparators embed the closed-form reference (the proven oracle formula in
# lsquant_chain_reference.py), so there is NO Python / data-file dependency to build
# or run ctest. The oracle .py remain committed only as the symbolic derivation.
set(CHAIN1D_KG_MOM  ${GOLDEN_DIR}/chain1d/NonEqOpVX-VXchain1dKPM_M64x64_stateexact.chebmom2D)
set(CHAIN1D_DOS_MOM ${GOLDEN_DIR}/chain1d/SpectralOp1chain1dKPM_M64_stateexact.chebmom1D)
# TOL justified by the single random-vector (R=1, N=2048) stochastic-trace error
# (~1/sqrt(N)): observed max residual ~0.05 (KG) / ~0.066 (DOS); 0.10 leaves margin
# while the structure is resolved to ~0.05. Deterministic under KPM_SEED=12345.
set(CHAIN1D_TOL 1e-9)   # exact trace -> machine precision (was 0.10 stochastic)

# #3 Kubo-Greenwood moment matrix vs B.15 normalized form.
add_executable(chain1d_analytic_test chain1d_analytic_test.cpp)
add_test(NAME chain1d_analytic
         COMMAND chain1d_analytic_test ${CHAIN1D_KG_MOM} ${CHAIN1D_TOL})

# #1 DOS Chebyshev moments vs delta_{m,0} (cheapest, most sensitive trace/rescale check).
add_executable(chain1d_dos_test chain1d_dos_test.cpp)
add_test(NAME chain1d_dos
         COMMAND chain1d_dos_test ${CHAIN1D_DOS_MOM} ${CHAIN1D_TOL})

# #2 Reconstructed DOS vs the closed form 1/(pi sqrt(4-E^2)) (oracle), rel err < 1% at the
# verified interior energies |E|<=1.5. The committed DOS_reconstructed.dat is the per-site
# DOS averaged over N_r=128 single-vector runs (M=64, exact bounds); regenerate with
# test/dos_reconstructed.sh. Reference is the oracle closed form (not another run).
set(CHAIN1D_DOS_RECON ${GOLDEN_DIR}/chain1d/DOS_reconstructed.dat)
set(DOS_RECON_TOL 0.005)   # exact trace -> kernel-limited; 0.5% (was 1% stochastic)
add_executable(dos_reconstructed_test dos_reconstructed_test.cpp)
add_test(NAME dos_reconstructed
         COMMAND dos_reconstructed_test ${CHAIN1D_DOS_RECON} ${DOS_RECON_TOL})

# Reproducibility + fresh-run physics. Same KPM_SEED => byte-identical moments across
# two runs (strict determinism gate), AND the fresh moments satisfy the analytic
# reference within CHAIN1D_TOL (checked by the C++ comparators, no Python). Does NOT
# byte-compare the fresh run to the committed golden (toolchain-dependent to ~1e-15);
# the committed golden's physics is guarded with tolerance by chain1d_analytic/_dos.
add_test(NAME chain1d_reproducible
         COMMAND bash ${CMAKE_CURRENT_SOURCE_DIR}/chain1d_reproduce.sh ${CMAKE_BINARY_DIR} ${CHAIN1D_TOL})

# --- spin-operator algebra (G2 #6 prerequisite) ---------------------------------
# The magnetic-chain spin operators (built by wannier2sparse from the spin doubling)
# must form su(2) and have [H,Sx]!=0 (precession reachable). Reads the committed CSRs;
# integer-valued Pauli operators, so a tight tolerance is appropriate.
set(SPINOPS_DIR ${CMAKE_CURRENT_SOURCE_DIR}/models/chain1d_mag/operators)
set(SPIN_TOL 1e-9)
add_executable(spin_algebra_test spin_algebra_test.cpp)
add_test(NAME spin_algebra
         COMMAND spin_algebra_test
                 ${SPINOPS_DIR}/chain1d_mag.HAM.CSR
                 ${SPINOPS_DIR}/chain1d_mag.SX.CSR
                 ${SPINOPS_DIR}/chain1d_mag.SY.CSR
                 ${SPINOPS_DIR}/chain1d_mag.SZ.CSR
                 ${SPIN_TOL})

# --- graphene regression + spectral-bounds safeguard ----------------------------
# Regenerates the graphene VX.VX moments from the committed operators (no BOUNDS file,
# so the Gershgorin+10% safeguard is exercised) and guards: determinism, regression vs
# the committed golden (numerical TOL, toolchain-robust), and that the bounds fallback
# is the Gershgorin estimate rather than the old [-100,100].
add_executable(chebmom_compare chebmom_compare.cpp)
set(GRAPHENE_REGR_TOL 1e-6)
add_test(NAME graphene_regression
         COMMAND bash ${CMAKE_CURRENT_SOURCE_DIR}/graphene_reproduce.sh ${CMAKE_BINARY_DIR} ${GRAPHENE_REGR_TOL})

# --- reconstruction NaN/Inf guard (alpha safeguard, Phase 1) ---------------------
# The Kubo-Bastin kernels carry 1/sqrt(1-x^2), NaN at the band edges x=+-1; the
# cumulative integral then poisons the whole curve (the historical all-NaN .dat).
# KPM_ALPHA=0.95 keeps the reconstruction grid off the singularity. This asserts the
# Kubo-Bastin reconstruction of the committed graphene moments is finite everywhere.
add_test(NAME reconstruction_nan_guard
         COMMAND bash ${CMAKE_CURRENT_SOURCE_DIR}/reconstruction_nan_guard.sh ${CMAKE_BINARY_DIR})

# --- ballistic MSD transport on the clean chain (Phase 4 #4) ---------------------
# Exact-trace mean-square displacement: assert MSD(E,t) ~ t^2 (ballistic) and that its
# energy-resolved coefficient tracks the conductivity slope sqrt(4-E^2)/2 (oracle
# dsigma/dt = (1/pi)sqrt(4-E^2); 0.6366 at E=0, 0.5513 at E=1). Deterministic.
set(MSD_M 256)
add_test(NAME msd_ballistic
         COMMAND bash ${CMAKE_CURRENT_SOURCE_DIR}/msd_ballistic.sh ${CMAKE_BINARY_DIR} ${MSD_M})

# --- Hall conductivity / Chern number on Haldane (Phase 4 #5) --------------------
# The Haldane Chern insulator (t1=-1,t2=0.15,phi=pi/2; 16x16) has a quantized Hall
# plateau in the gap: sigma_xy = C*e^2/h with C=+1. Regenerates VX-VY moments by exact
# trace from the committed operators, reconstructs Kubo-Bastin (NaN-free via alpha), and
# asserts (plateau/A)*2pi = +1 within finite-size tolerance. This is LinQT's quantized
# "conductance quantum" test. Reference = the topological invariant (oracle), not a run.
set(HALL_M 128)         # 8s exact-trace; C=+1.017 at this M (1.7% finite-M error)
set(HALL_TOL 0.05)
add_test(NAME hall_haldane
         COMMAND bash ${CMAKE_CURRENT_SOURCE_DIR}/hall_haldane.sh ${CMAKE_BINARY_DIR} ${HALL_M} ${HALL_TOL})

# --- spin precession on the magnetic chain (Phase 4 #6) --------------------------
# A +x-polarised spin in the Zeeman field H_hop - J_ex*sigma_z precesses:
# <S(t)> = cos(wt)x - sin(wt)y, S_z=0, w = 2*J_ex/hbar (physical hbar=0.6582 eV*fs;
# period 20.68 fs), E_F-independent. Density projection path, exact trace -> machine
# precision. Validates LinQT's time-evolution + spin-projection machinery end-to-end.
set(PREC_M 32)
add_test(NAME spin_precession
         COMMAND bash ${CMAKE_CURRENT_SOURCE_DIR}/spin_precession.sh ${CMAKE_BINARY_DIR} ${PREC_M})

# --- synthetic reconstruction oracle (Phase R, operator-free) --------------------
# Feeds analytically-known delta moments (mu_m = T_m(E0)) through the PRODUCTION kernel
# delta_chebF and checks peak position, Jackson width law (~pi/M), unit weight, and the
# alpha edge singularity -- isolating the reconstruction knobs from operators/recursion.
add_executable(recon_kernel_oracle_test recon_kernel_oracle_test.cpp)
target_include_directories(recon_kernel_oracle_test PRIVATE ${CMAKE_SOURCE_DIR}/include)
target_link_libraries(recon_kernel_oracle_test PRIVATE kpm_lib)
add_test(NAME recon_kernel_oracle COMMAND recon_kernel_oracle_test 0.10)

# --- Kubo-Greenwood reconstruction regression (Bug 1, Phase 0) -------------------
# Pins PR #19: kuboGreenwoodFromChebmom multiplied by ShiftFactor() (==0 for centred bands)
# so sigma_KG was identically 0. chain1d_analytic/graphene_regression guard MOMENTS, not the
# reconstruction prefactor -- a regression to sigma_KG==0 would pass CI silently. This asserts
# the reconstructed sigma_KG(E) from the committed graphene moments is finite and non-zero.
add_test(NAME greenwood_regression
         COMMAND bash ${CMAKE_CURRENT_SOURCE_DIR}/greenwood_regression.sh ${CMAKE_BINARY_DIR})

# --- Phase 1 safety-factor invariants (alpha decoupled, HBAR single, alpha runtime) ---
# Guards: (1) ScaleFactor()/ShiftFactor() use CUTOFF only -- alpha must NEVER enter the moment
# rescale (it would break the analytic moment identities); (2) the HBAR literal lives only in
# units.hpp; (3) LSQUANT_RECON_ALPHA changes the reconstruction grid at runtime.
add_test(NAME safety_factors_guard
         COMMAND bash ${CMAKE_CURRENT_SOURCE_DIR}/safety_factors_guard.sh ${CMAKE_BINARY_DIR})

# --- descriptor schema conformance (Phase 2) ------------------------------------
# lsquant's .desc reader mirrors the w2s producer schema (duplicated-but-versioned, no
# submodule). This parses committed w2s-produced .desc files and asserts every field
# round-trips, so reader and writer cannot drift silently.
add_executable(descriptor_conformance_test descriptor_conformance_test.cpp
               ${CMAKE_SOURCE_DIR}/src/operator_descriptor.cpp)
target_include_directories(descriptor_conformance_test PRIVATE ${CMAKE_SOURCE_DIR}/include)
add_test(NAME descriptor_conformance
         COMMAND descriptor_conformance_test
                 ${GOLDEN_DIR}/operators/graphene_golden.HAM.desc
                 ${GOLDEN_DIR}/operators/graphene_golden.VX.desc)

# --- spectral bounds ingested from .desc (Phase 2) ------------------------------
# Reconstructs the clean-chain DOS with NO BOUNDS file: the window (a,b) must come from the
# committed chain1d.HAM.desc (Lanczos [-2,2]) and the DOS must match 1/(pi sqrt(4-E^2)).
add_test(NAME descriptor_bounds
         COMMAND bash ${CMAKE_CURRENT_SOURCE_DIR}/descriptor_bounds.sh ${CMAKE_BINARY_DIR})

# --- manual tight-binding traceability via descriptor basis tags (Phase 3) ------
# Builds SX/SY/SZ from the basis tags in chain1d_mag.HAM.desc and asserts they reproduce the
# committed CSRs exactly -- the conformance gate for the manual-TB descriptor path.
add_executable(tb_descriptor_test tb_descriptor_test.cpp
               ${CMAKE_SOURCE_DIR}/src/operator_descriptor.cpp)
target_include_directories(tb_descriptor_test PRIVATE ${CMAKE_SOURCE_DIR}/include)
add_test(NAME tb_descriptor
         COMMAND tb_descriptor_test
                 ${SPINOPS_DIR}/chain1d_mag.HAM.desc
                 ${SPINOPS_DIR}/chain1d_mag.SX.CSR
                 ${SPINOPS_DIR}/chain1d_mag.SY.CSR
                 ${SPINOPS_DIR}/chain1d_mag.SZ.CSR)

# --- JSON run-config equivalence + inspect (Phase 4) ----------------------------
# A run.json reproduces the argv path byte-for-byte (same computation, new description), and
# `lsquant inspect` reports the run + resolved numerical provenance (alpha, bounds source, hbar).
add_test(NAME config_equivalence
         COMMAND bash ${CMAKE_CURRENT_SOURCE_DIR}/config_equivalence.sh ${CMAKE_BINARY_DIR})

# --- CorrOp state-file regression (Phase 4 bug fix) -----------------------------
# Pins the fix to the atoi(argv[7]) bug: the VAC driver now loads its 7th arg as a state file
# (was misparsed as a vector count, so "exact" -> 0 -> default). Asserts the exact-trace spec
# is honored (output tagged stateexact) and deterministic.
add_test(NAME corrop_state_honored
         COMMAND bash ${CMAKE_CURRENT_SOURCE_DIR}/corrop_state_honored.sh ${CMAKE_BINARY_DIR})

# --- FFT-node reconstruction backend regression (Phase 5 slice 3) ---------------
# Pins the FFT Chebyshev-node Kubo reconstruction (previously untested) by byte-comparing the
# FFT-grid Bastin/Greenwood of the committed graphene moments to committed goldens. Safety net
# for the GRID_FFT_NODES backend in reconstruct_kubo.
add_test(NAME fft_grid_regression
         COMMAND bash ${CMAKE_CURRENT_SOURCE_DIR}/fft_grid_regression.sh ${CMAKE_BINARY_DIR})
