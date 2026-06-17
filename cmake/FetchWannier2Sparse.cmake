execute_process(
  COMMAND bash ${CMAKE_SOURCE_DIR}/scripts/fetch_wannier2sparse.sh
  OUTPUT_VARIABLE W2S_BIN OUTPUT_STRIP_TRAILING_WHITESPACE
  RESULT_VARIABLE W2S_RC)
if(W2S_RC EQUAL 0 AND EXISTS "${W2S_BIN}")
  message(STATUS "wannier2sparse: ${W2S_BIN}")
  set(LSQUANT_W2S_BIN "${W2S_BIN}" CACHE INTERNAL "external wannier2sparse binary")
else()
  message(WARNING "wannier2sparse unavailable (offline?). Goldens will not be (re)generated; "
                  "committed goldens remain usable.")
endif()
