# LSQUANT — Linear Quantum Transport (KPM core)
#
# Unified C++17 core + pybind11 bindings, with Eigen3, FFTW3, GSL and OpenMP.
cmake_minimum_required(VERSION 3.16)

project(LSQUANT VERSION 2.0.0 LANGUAGES CXX)

include(GNUInstallDirs)

# --- Language / build type ---------------------------------------------------
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type (Debug, Release, RelWithDebInfo)" FORCE)
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS Debug Release RelWithDebInfo)
endif()

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}")
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")

# --- Options -----------------------------------------------------------------
set(LSQUANT_BUILD_PROFILE "golden" CACHE STRING
    "Numerical build profile: 'golden' (portable IEEE, no -march=native or -ffast-math) or 'native' (-O3 -march=native -mtune=native -fno-math-errno, defines LSQUANT_NATIVE_BUILD=1)")
set_property(CACHE LSQUANT_BUILD_PROFILE PROPERTY STRINGS golden native)
if(NOT LSQUANT_BUILD_PROFILE MATCHES "^(golden|native)$")
    message(FATAL_ERROR "LSQUANT_BUILD_PROFILE must be 'golden' or 'native' (got '${LSQUANT_BUILD_PROFILE}')")
endif()
message(STATUS "LSQUANT build profile: ${LSQUANT_BUILD_PROFILE}")

option(LSQUANT_BUILD_PYTHON     "Build the pybind11 extension module lsquant._core" ON)
option(LSQUANT_BUILD_BENCHMARKS "Build the C++ kernel benchmark executable stage8_kernels" OFF)
option(LSQUANT_BUILD_DEV_CHECKS "Build the C++ developer self-check executable stage_1_check" ON)

# --- Shared warning flags ----------------------------------------------------
add_library(lsquant_warnings INTERFACE)
target_compile_options(lsquant_warnings INTERFACE -Wall -Wextra -Wpedantic)

# --- Dependencies ------------------------------------------------------------
# Eigen3: prefer a CMake config package (Eigen3::Eigen); fall back to a header dir.
find_package(Eigen3 3.3 QUIET NO_MODULE)
if(NOT TARGET Eigen3::Eigen)
    find_path(EIGEN3_INCLUDE_DIR
        NAMES signature_of_eigen3_matrix_library
        HINTS $ENV{CONDA_PREFIX}/include/eigen3 /usr/include/eigen3 /usr/local/include/eigen3
        DOC "Directory containing the Eigen/ header folder")
    if(NOT EIGEN3_INCLUDE_DIR)
        message(FATAL_ERROR
            "Eigen3 (>= 3.3) not found. Install libeigen3-dev or pass -DEIGEN3_INCLUDE_DIR=/path/to/eigen3.")
    endif()
    add_library(Eigen3::Eigen INTERFACE IMPORTED)
    set_target_properties(Eigen3::Eigen PROPERTIES
        INTERFACE_INCLUDE_DIRECTORIES "${EIGEN3_INCLUDE_DIR}")
    message(STATUS "Eigen3 headers: ${EIGEN3_INCLUDE_DIR}")
endif()

# OpenMP
find_package(OpenMP QUIET)
if(NOT OpenMP_CXX_FOUND)
    message(FATAL_ERROR
        "OpenMP for C++ not found. Install libomp-dev / gcc OpenMP support or pass -DOpenMP_CXX_FLAGS=...")
endif()

# FFTW3 (double precision).
find_package(PkgConfig QUIET)
if(NOT PkgConfig_FOUND)
    message(FATAL_ERROR "pkg-config not found; required to locate FFTW3.")
endif()
pkg_check_modules(FFTW3 QUIET IMPORTED_TARGET fftw3)
if(NOT FFTW3_FOUND)
    message(FATAL_ERROR
        "FFTW3 (double precision, 'fftw3') not found via pkg-config. Install libfftw3-dev.")
endif()

# GSL and pybind11 (required when LSQUANT_BUILD_PYTHON=ON; GSL is used only by compat/).
if(LSQUANT_BUILD_PYTHON)
    find_package(GSL QUIET)
    if(NOT GSL_FOUND)
        message(FATAL_ERROR
            "GSL (GNU Scientific Library) not found. Install libgsl-dev or pass -DGSL_ROOT_DIR=/path/to/gsl.")
    endif()
    if(NOT pybind11_DIR)
        find_program(LSQUANT_PYTHON_EXECUTABLE NAMES python3 python)
        if(LSQUANT_PYTHON_EXECUTABLE)
            execute_process(
                COMMAND "${LSQUANT_PYTHON_EXECUTABLE}" -m pybind11 --cmakedir
                OUTPUT_VARIABLE _pybind11_cmakedir
                OUTPUT_STRIP_TRAILING_WHITESPACE
                ERROR_QUIET
                RESULT_VARIABLE _pybind11_rc
            )
            if(_pybind11_rc EQUAL 0 AND EXISTS "${_pybind11_cmakedir}")
                list(APPEND CMAKE_PREFIX_PATH "${_pybind11_cmakedir}")
            endif()
        endif()
    endif()
    find_package(pybind11 CONFIG QUIET)
    if(NOT pybind11_FOUND)
        message(FATAL_ERROR
            "pybind11 not found (required when LSQUANT_BUILD_PYTHON=ON). "
            "Install pybind11 ('pip install pybind11' or 'apt install pybind11-dev') "
            "or pass -Dpybind11_DIR=$(python3 -m pybind11 --cmakedir).")
    endif()
endif()

# --- Targets -----------------------------------------------------------------
add_subdirectory(lsquant lsquant_core)
