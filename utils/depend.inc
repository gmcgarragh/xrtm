callxrtm.o: callxrtm.c input_util.h
input_util.o: input_util.c input_local.h input_util.h input_yylex.h \
 parse_data.h
input_yylex.o: input_yylex.c input_local.h input_util.h
test_core.o: test_core.c test.h input_util.h test_result.h \
 prototypes/test_result_p.h test_xrtm.h prototypes/test_xrtm_p.h \
 ../utils2/test2.h test_core.h prototypes/test_core_p.h test_loop.h \
 prototypes/test_loop_p.h test_macros.h
test_execute.o: test_execute.c input_util.h test.h test_result.h \
 prototypes/test_result_p.h test_xrtm.h prototypes/test_xrtm_p.h \
 ../utils2/test2.h test_util.h prototypes/test_util_p.h test_write.h \
 prototypes/test_write_p.h
test_loop.o: test_loop.c test.h input_util.h test_result.h \
 prototypes/test_result_p.h test_xrtm.h prototypes/test_xrtm_p.h \
 ../utils2/test2.h test_execute.h prototypes/test_execute_p.h test_loop.h \
 prototypes/test_loop_p.h test_macros.h test_util.h \
 prototypes/test_util_p.h
test_result.o: test_result.c test_result.h prototypes/test_result_p.h
test_util.o: test_util.c test.h input_util.h test_result.h \
 prototypes/test_result_p.h test_xrtm.h prototypes/test_xrtm_p.h \
 ../utils2/test2.h test_util.h prototypes/test_util_p.h
test_write.o: test_write.c test.h input_util.h test_result.h \
 prototypes/test_result_p.h test_xrtm.h prototypes/test_xrtm_p.h \
 ../utils2/test2.h
test_xrtm.o: test_xrtm.c test.h input_util.h test_result.h \
 prototypes/test_result_p.h test_xrtm.h prototypes/test_xrtm_p.h \
 ../utils2/test2.h
testxrtm.o: testxrtm.c input_util.h test.h test_result.h \
 prototypes/test_result_p.h test_xrtm.h prototypes/test_xrtm_p.h \
 ../utils2/test2.h test_core.h prototypes/test_core_p.h test_macros.h \
 test_util.h prototypes/test_util_p.h test_write.h \
 prototypes/test_write_p.h
