set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

option(LSQUANT_BUILD_PYTHON "Build the pybind11 module lsquant._core" ON)
option(LSQUANT_BUILD_BENCHMARKS "Build stage-8 measurement drivers (not gates)" OFF)

add_library(lsquant_core STATIC
    core/algebra/sparse_matrix.cpp
    core/algebra/block_diagonal_operator.cpp
    core/algebra/bloch_transform.cpp
    core/algebra/composite_operators.cpp
    core/algebra/linalg.cpp
    core/algebra/state_generator.cpp
    core/algebra/spectral_bounds.cpp
    core/algebra/metric_policy.cpp
    core/algebra/overlap_appliers.cpp
    core/chebyshev/rescaled_operator.cpp
    core/chebyshev/chebyshev_recursion.cpp
    core/chebyshev/chebyshev_vectors.cpp
    core/chebyshev/chebyshev_moments.cpp
    core/chebyshev/moment_engine.cpp
    core/chebyshev/contraction.cpp
    core/algebra/fft_moment_contraction.cpp
    io/csr_text.cpp
    io/chebmom_file.cpp
    io/curve_text.cpp
)

set_target_properties(lsquant_core PROPERTIES POSITION_INDEPENDENT_CODE ON)
target_include_directories(lsquant_core PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")
target_compile_features(lsquant_core PUBLIC cxx_std_17)
target_compile_options(lsquant_core PUBLIC -Wall -Wextra)
target_link_libraries(lsquant_core PUBLIC Eigen3::Eigen OpenMP::OpenMP_CXX PkgConfig::FFTW3)
if(LSQUANT_BUILD_PROFILE STREQUAL "native")
    target_compile_options(lsquant_core PUBLIC -O3 -march=native -mtune=native -fno-math-errno)
    target_compile_definitions(lsquant_core PUBLIC LSQUANT_NATIVE_BUILD=1)
endif()

add_executable(stage1_moments apps/stage1_moments.cpp)
set_target_properties(stage1_moments PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}")
target_compile_options(stage1_moments PRIVATE -Wall -Wextra)
target_link_libraries(stage1_moments PRIVATE lsquant_core)

add_executable(stage6_nonorth apps/stage6_nonorth.cpp)
set_target_properties(stage6_nonorth PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}")
target_compile_options(stage6_nonorth PRIVATE -Wall -Wextra)
target_link_libraries(stage6_nonorth PRIVATE lsquant_core)

add_executable(stage_1_check "${PROJECT_SOURCE_DIR}/tests/dev/stage_1_check.cpp")
set_target_properties(stage_1_check PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}")
target_compile_options(stage_1_check PRIVATE -Wall -Wextra)
target_link_libraries(stage_1_check PRIVATE lsquant_core)

if(LSQUANT_BUILD_PYTHON)
    find_package(pybind11 CONFIG REQUIRED)
    pybind11_add_module(_core MODULE
        api/bindings.cpp
        compat/legacy_kernels.cpp
    )
    set_target_properties(_core PROPERTIES
        LIBRARY_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
        RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
    )
    target_compile_options(_core PRIVATE -Wall -Wextra)
    target_link_libraries(_core PRIVATE lsquant_core GSL::gsl)
endif()

if(LSQUANT_BUILD_BENCHMARKS)
    add_executable(stage8_kernels "${PROJECT_SOURCE_DIR}/bench/kernel_bench.cpp")
    set_target_properties(stage8_kernels PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}")
    target_link_libraries(stage8_kernels PRIVATE lsquant_core)
endif()
