# docs/CMakeLists.txt — Doxygen documentation target (opt-in).
#
# Included from the top-level CMakeLists only when BUILD_DOCS=ON, so a normal
# build and the test suite never require Doxygen or Graphviz and are completely
# unaffected. The target configures docs/Doxyfile.in with the CMake-derived
# project name/version and runs Doxygen against it, writing HTML + XML under
# ${CMAKE_BINARY_DIR}/docs (outside the source tree).
#
# Pitfall: dot graphs need Graphviz at build time. find_package locates it; if
# absent we keep building (HAVE_DOT=NO) rather than failing, and say so.
#
# Usage:
#   cmake -B build -DBUILD_DOCS=ON
#   cmake --build build --target docs
#   # -> build/docs/html/index.html , build/docs/xml/ , build/docs/doxygen-warnings.log

find_package(Doxygen REQUIRED)

# --- Values fed into Doxyfile.in (@VAR@ slots) -------------------------------
set(DOXYGEN_PROJECT_NAME    "${PROJECT_NAME}")     # = LSQUANT (project(), not a literal)
set(DOXYGEN_PROJECT_NUMBER  "${PROJECT_VERSION}")  # = 2.0.0   (single source of truth)
set(DOXYGEN_OUTPUT_DIR      "${CMAKE_BINARY_DIR}/docs")
set(DOXYGEN_WARN_LOGFILE    "${DOXYGEN_OUTPUT_DIR}/doxygen-warnings.log")
set(DOXYGEN_STRIP_FROM_PATH "${CMAKE_SOURCE_DIR}")

# Narrative pages first (order fixes the page sequence), then the API headers and
# sources. Absolute paths: the target runs Doxygen from the build tree.
set(DOXYGEN_INPUT
    "\"${CMAKE_SOURCE_DIR}/docs/main.md\" \
\"${CMAKE_SOURCE_DIR}/docs/spectral_scales.md\" \
\"${CMAKE_SOURCE_DIR}/docs/reconstruction.md\" \
\"${CMAKE_SOURCE_DIR}/docs/time_evolution.md\" \
\"${CMAKE_SOURCE_DIR}/docs/spinrelaxation.md\" \
\"${CMAKE_SOURCE_DIR}/include\" \
\"${CMAKE_SOURCE_DIR}/src\"")
set(DOXYGEN_MAINPAGE     "${CMAKE_SOURCE_DIR}/docs/main.md")
set(DOXYGEN_EXAMPLE_PATH "${CMAKE_SOURCE_DIR}/examples")
set(DOXYGEN_IMAGE_PATH   "${CMAKE_SOURCE_DIR}/docs")

# Graphviz: degrade gracefully rather than fail the configure step.
if(DOXYGEN_DOT_FOUND)
    set(DOXYGEN_HAVE_DOT "YES")
    get_filename_component(DOXYGEN_DOT_PATH "${DOXYGEN_DOT_EXECUTABLE}" DIRECTORY)
    message(STATUS "LSQUANT docs: Graphviz dot found (${DOXYGEN_DOT_EXECUTABLE}); graphs enabled")
else()
    set(DOXYGEN_HAVE_DOT "NO")
    set(DOXYGEN_DOT_PATH "")
    message(WARNING "LSQUANT docs: Graphviz not found; building without dot graphs")
endif()

set(DOXYFILE_IN  "${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in")
set(DOXYFILE_OUT "${DOXYGEN_OUTPUT_DIR}/Doxyfile")

file(MAKE_DIRECTORY "${DOXYGEN_OUTPUT_DIR}")
configure_file("${DOXYFILE_IN}" "${DOXYFILE_OUT}" @ONLY)

add_custom_target(docs
    COMMAND ${DOXYGEN_EXECUTABLE} "${DOXYFILE_OUT}"
    WORKING_DIRECTORY "${CMAKE_BINARY_DIR}"
    COMMENT "Generating LSQUANT API + narrative docs with Doxygen (HTML + XML)"
    VERBATIM)
