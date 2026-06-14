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
set(CHAIN1D_KG_MOM  ${GOLDEN_DIR}/chain1d/NonEqOpVX-VXchain1dKPM_M64x64_statedefault.chebmom2D)
set(CHAIN1D_DOS_MOM ${GOLDEN_DIR}/chain1d/SpectralOp1chain1dKPM_M64_statedefault.chebmom1D)
# TOL justified by the single random-vector (R=1, N=2048) stochastic-trace error
# (~1/sqrt(N)): observed max residual ~0.05 (KG) / ~0.066 (DOS); 0.10 leaves margin
# while the structure is resolved to ~0.05. Deterministic under KPM_SEED=12345.
set(CHAIN1D_TOL 0.10)

# #3 Kubo-Greenwood moment matrix vs B.15 normalized form.
add_executable(chain1d_analytic_test chain1d_analytic_test.cpp)
add_test(NAME chain1d_analytic
         COMMAND chain1d_analytic_test ${CHAIN1D_KG_MOM} ${CHAIN1D_TOL})

# #1 DOS Chebyshev moments vs delta_{m,0} (cheapest, most sensitive trace/rescale check).
add_executable(chain1d_dos_test chain1d_dos_test.cpp)
add_test(NAME chain1d_dos
         COMMAND chain1d_dos_test ${CHAIN1D_DOS_MOM} ${CHAIN1D_TOL})

# Reproducibility: same KPM_SEED => byte-identical moments across two runs and
# identical to the committed goldens (mandatory determinism gate). Recomputes from
# the committed CSRs, so it needs no external tool at test time.
add_test(NAME chain1d_reproducible
         COMMAND bash ${CMAKE_CURRENT_SOURCE_DIR}/chain1d_reproduce.sh ${CMAKE_BINARY_DIR})
