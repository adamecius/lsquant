# LSQUANT — Linear Quantum Transport (KPM core)
#
# Reflects the build as it actually is: a single Eigen numeric backend plus
# FFTW3, GSL and OpenMP. The former INTEL_MKL / EIGEN / BLAS_LIB switches and the
# (nonexistent) icc+mkl toolchain are gone; dependencies are resolved with
# find_package and consumed as imported targets, so no hand-passed toolchain file
# is required:  cmake -B build && cmake --build build -j
cmake_minimum_required(VERSION 3.16)

project(LSQUANT VERSION 2.0.0 LANGUAGES CXX)

# --- Language / build type ---------------------------------------------------
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS Debug Release RelWithDebInfo)
endif()

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}")
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")

# --- Options -----------------------------------------------------------------
option(LSQUANT_BUILD_UTILITIES   "Build the Python helper utilities" OFF)
option(LSQUANT_BUILD_INPUT_TOOLS "Fetch + build the external wannier2sparse for tests" ON)

# --- Dependencies ------------------------------------------------------------
# Eigen3: prefer a real config-package (Eigen3::Eigen); fall back to a header dir.
find_package(Eigen3 3.3 QUIET NO_MODULE)
if(NOT TARGET Eigen3::Eigen)
    find_path(EIGEN3_INCLUDE_DIR
        NAMES signature_of_eigen3_matrix_library
        HINTS $ENV{CONDA_PREFIX}/include/eigen3 /usr/include/eigen3 /usr/local/include/eigen3
        DOC "Directory containing the Eigen/ folder")
    if(NOT EIGEN3_INCLUDE_DIR)
        message(FATAL_ERROR
            "Eigen3 not found. Install libeigen3-dev or pass -DEIGEN3_INCLUDE_DIR=/path/to/eigen3")
    endif()
    add_library(Eigen3::Eigen INTERFACE IMPORTED)
    set_target_properties(Eigen3::Eigen PROPERTIES
        INTERFACE_INCLUDE_DIRECTORIES "${EIGEN3_INCLUDE_DIR}")
    message(STATUS "Eigen3 headers: ${EIGEN3_INCLUDE_DIR}")
endif()

find_package(OpenMP REQUIRED)
find_package(GSL    REQUIRED)   # -> GSL::gsl (carries gslcblas)
find_package(Threads REQUIRED)

# FFTW3 (single precision not used). The OpenMP-threaded variant lives in a
# separate library (libfftw3_omp); link it when present.
find_package(PkgConfig REQUIRED)
pkg_check_modules(FFTW3 REQUIRED IMPORTED_TARGET fftw3)
find_library(FFTW3_OMP_LIBRARY fftw3_omp)
if(NOT FFTW3_OMP_LIBRARY)
    message(STATUS "libfftw3_omp not found; building without threaded FFTW")
    set(FFTW3_OMP_LIBRARY "")
endif()

include_directories(${PROJECT_SOURCE_DIR}/include)

# --- Targets -----------------------------------------------------------------
add_subdirectory(src)

if(LSQUANT_BUILD_UTILITIES)
    add_subdirectory(utilities)
endif()

# Golden-input generation depends on the external wannier2sparse, which is NOT
# vendored or distributed. It is fetched + built only on demand for tests.
# (cmake/FetchWannier2Sparse.cmake wraps scripts/fetch_wannier2sparse.sh, pinned
#  to a commit of github.com/adamecius/wannier2sparse.)
if(LSQUANT_BUILD_INPUT_TOOLS)
    include(cmake/FetchWannier2Sparse.cmake OPTIONAL)
endif()

include(CTest)
if(BUILD_TESTING AND EXISTS "${PROJECT_SOURCE_DIR}/test/CMakeLists.txt")
    add_subdirectory(test)
endif()
