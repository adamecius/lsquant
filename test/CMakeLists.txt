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

# --- chain1d ANALYTIC golden (Kubo-Greenwood moments vs thesis Eq. 4.9) -------
# Native CTest (already include()'d at the root) -- not a second framework.
# Dependency-free C++ comparator (parses .chebmom2D + the checked-in reference);
# no Python/numpy needed at test time.
set(CHAIN1D_MOM ${GOLDEN_DIR}/chain1d/NonEqOpVX-VXchain1dKPM_M64x64_statedefault.chebmom2D)
set(CHAIN1D_REF ${GOLDEN_DIR}/chain1d_analytic_M.txt)
# TOL justified by the single random-vector (R=1, N=2048) stochastic-trace error:
# observed max off-structure residual ~0.049 (~2 sigma, sigma~1/sqrt(N)); 0.10 is
# ~2x margin while the diagonal (0.5/1.0) and |m-n|=2 band (-0.5) are resolved to
# ~0.05. Deterministic under KPM_SEED=12345, so the residual is fixed run-to-run.
set(CHAIN1D_TOL 0.10)

add_executable(chain1d_analytic_test chain1d_analytic_test.cpp)
add_test(NAME chain1d_analytic
         COMMAND chain1d_analytic_test ${CHAIN1D_MOM} ${CHAIN1D_REF} ${CHAIN1D_TOL})

# Reproducibility: same KPM_SEED => byte-identical moments across two runs and
# identical to the committed golden (mandatory determinism gate). Recomputes from
# the committed CSRs, so it needs no external tool at test time.
add_test(NAME chain1d_reproducible
         COMMAND bash ${CMAKE_CURRENT_SOURCE_DIR}/chain1d_reproduce.sh ${CMAKE_BINARY_DIR})
