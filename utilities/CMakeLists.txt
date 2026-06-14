# Helper utilities. The vendored wannier2sparse and wannier2sparse2 generators
# have been removed: input generation now uses the standalone, externally-built
# github.com/adamecius/wannier2sparse (see LINQT_BUILD_INPUT_TOOLS). siesta2sparse
# is pure Python and self-contained.
add_subdirectory(python)
