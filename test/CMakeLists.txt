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
