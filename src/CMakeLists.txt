# @author  Jose H. Garcia
#
# KPM core library + command-line drivers.
# Single Eigen numeric backend; the dead INTEL_MKL branch (which compiled the
# Eigen sources anyway) and the mkl_*.cpp / *.bak / *.save files have been removed.

set(KPM_LIB_SOURCES
    chebyshev_solver.cpp
    quantum_states.cpp
    chebyshev_vectors.cpp
    chebyshev_moments.cpp
    chebyshev_moments2D.cpp
    chebyshev_momentsAC.cpp
    chebyshev_moments1D.cpp
    chebyshev_momentsTD.cpp
    chebyshev_momentsLocal.cpp
    sparse_matrix.cpp
    operator_descriptor.cpp
    run_config.cpp
    inspect.cpp
    reconstruction.cpp
    special_functions.cpp
    Kubo_solver_FFT.cpp
    Kubo_solver_FFT_postProcess.cpp
    Greenwood_FFTs.cpp
    Bastin_FFTs.cpp
    Device/Device.cpp
    Device/Graphene_KaneMele.cpp
    chebyshev_solver_Device.cpp
    chebyshev_solver_Graphene_NNN.cpp
    chebyshev_vectors_Device.cpp
    chebyshev_moments_Graphene_NNN.cpp
    chebyshev_moments_Device.cpp
    chebyshev_moments1D_Graphene_NNN.cpp
    chebyshev_vectors_Graphene_NNN.cpp
    chebyshev_momentsTD_Device.cpp
    chebyshev_momentsTD_Graphene_NNN.cpp
    eigen_sparse_matrix.cpp
    eigen_linear_algebra.cpp
)

add_library(kpm_lib STATIC ${KPM_LIB_SOURCES})

target_include_directories(kpm_lib PUBLIC ${PROJECT_SOURCE_DIR}/include)

target_link_libraries(kpm_lib PUBLIC
    Eigen3::Eigen
    GSL::gsl
    PkgConfig::FFTW3
    ${FFTW3_OMP_LIBRARY}
    OpenMP::OpenMP_CXX
    Threads::Threads
    ${CMAKE_DL_LIBS}
)

# --- Command-line drivers ----------------------------------------------------
# "<target-name>|<source-file>" — names preserved verbatim from the legacy file
# (several intentionally differ from the source stem). The replicated
# *FromChebmom drivers are still built individually here; collapsing them into a
# single parametrised driver is a separate source-level refactor (Plan 2), not a
# build-file change.
set(LINQT_APPS
    "inline_compute-kpm-nonEqOp|inline_compute-kpm-nonEqOp.cpp"
    "inline_compute-kpm-CorrOp|inline_compute-kpm-CorrOp.cpp"
    "inline_compute-kpm-spectralOp|inline_compute-kpm-spectralOp.cpp"
    "inline_kpm_S-inversion|inline_kpm-S-inversion.cpp"
    "inline_compute-kpm-GetDeltaPhi|inline_compute-kpm-GetDeltaPhi.cpp"
    "inline_compute-kpm-LocalDOS|inline_compute-kpm-LocalDOS.cpp"
    "inline_compute-kpm-spacialresolution|inline_compute-kpm-spacialresolution.cpp"
    "inline_spectralFunctionFromChebmom|spectralFunctionFromChebmom.cpp"
    "inline_localDensitiesFromChebmom|localDensitiesFromChebmom.cpp"
    "inline_kuboGreenwoodFromChebmom|kuboGreenwoodFromChebmom.cpp"
    "inline_kuboBastinIFromChebmom|kuboBastinIFromChebmom.cpp"
    "inline_kuboBastinIIFromChebmom|kuboBastinIIFromChebmom.cpp"
    "inline_kuboBastinFromChebmom|kuboBastinFromChebmom.cpp"
    "inline_kuboBastinSeaFromChebmom|kuboBastinSeaFromChebmom.cpp"
    "inline_kuboBastinSurfFromChebmom|kuboBastinSurfFromChebmom.cpp"
    "inline_timeCorrelationsFromChebmom|timeCorrelationsFromChebmom.cpp"
    "inline_compute-Kubo-kpm-FFT|inline_compute-Kubo-kpm-FFT.cpp"
    "inline_opticalConductivityFromChebmom|opticalConductivityFromChebmom.cpp"
    "inline_spectralFunctionFromChebmom_FFTgrid|spectralFunctionFromChebmom_FFTgrid.cpp"
    "inline_compute-kpm-TimeEvProjetedOpWF|inline_compute-kpm-TimeEvProjetedOp-Projecting_WF.cpp"
    "inline_spectralFunctionFromChebmom_Lorentz|spectralFunctionFromChebmom_Lorentz.cpp"
    "inline_kuboGreenwoodFromChebmom_FFTgrid|kuboGreenwoodFromChebmom_FFTgrid.cpp"
    "inline_compute-kpm-MeanSquareDisplacement|inline_compute-kpm-MeanSquareDisplacement.cpp"
    "inline_compute-kpm-MeanSquareDisplacement_DeviceSim|inline_compute-kpm-MeanSquareDisplacement_DeviceSim.cpp"
    "inline_compute-kpm-MeanSquareDisplacement_Graphene_NNN_2|inline_compute-kpm-MeanSquareDisplacement_Graphene_NNN_2.cpp"
    "inline_kuboBastinFromChebmom_FFTgrid|kuboBastinFromChebmom_FFTgrid.cpp"
    "inline_compute-kpm-TimeEvProjetedOp|inline_compute-kpm-TimeEvProjetedOp.cpp"
    "inline_compute-kpm-spectralOp-nonOrth|inline_compute-kpm-spectralOp-nonOrth.cpp"
    "inline_compute-kpm-TimeEvOp|inline_compute-kpm-TimeEvOp.cpp"
    "inline_compute-kpm-TimeEvOp_Device|inline_compute-kpm-TimeEvOp_Device.cpp"
    "kpm-nonOrth-precision-test|kpm-nonOrth-precision-test.cpp"
    "inline_kpm-DOS-standalone|inline_kpm-DOS-standalone.cpp"
    "inline_kpm-DOS-standalone_Graphene_NNN|inline_kpm-DOS-standalone_Graphene_NNN.cpp"
    "device-postProcess|device-postProcess.cpp"
    "device-postProcess_Graphene_NNN|device-postProcess_Graphene_NNN.cpp"
    "inline_compute-Kubo-kpm-FFT-nonOrth|inline_compute-Kubo-kpm-FFT-nonOrth.cpp"
    "inline_compute-kpm-TimeEvOpwithWF|inline_compute-kpm-TimeEvOpwithWF.cpp"
    "lsquant_inspect|lsquant_inspect.cpp"
    "lsquant|lsquant.cpp"
)

foreach(entry IN LISTS LINQT_APPS)
    string(REPLACE "|" ";" parts "${entry}")
    list(GET parts 0 app_name)
    list(GET parts 1 app_src)
    add_executable(${app_name} ${app_src})
    target_link_libraries(${app_name} PRIVATE kpm_lib)
endforeach()
