include(test_macros)
execute_process(COMMAND ${ELMERGRID_BIN} 1 2 radiation2)
RUN_ELMER_TEST()
