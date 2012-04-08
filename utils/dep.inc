callxrtm.o: callxrtm.c input_util.h
input_util.o: input_util.c input_local.h input_util.h input_yylex.h \
 parse_data.h
input_yylex.o: input_yylex.c input_local.h input_util.h
test_core.o: test_core.c test.h input_util.h ../config.h test_macros.h \
 test_util.h
test_errors.o: test_errors.c test.h input_util.h ../config.h \
 test_macros.h
test_execute.o: test_execute.c input_util.h test.h ../config.h
test_util.o: test_util.c test.h input_util.h ../config.h test_util.h \
 test_macros.h
testxrtm.o: testxrtm.c input_util.h test.h ../config.h test_macros.h
