include ../make.inc
include ../common.inc

SUBDIRS =

OBJECTS_XRTM_INPUT = input_util.o \
                     input_yylex.o

OBJECTS_XRTM_TEST  = test_core.o \
                     test_execute.o \
                     test_loop.o \
                     test_result.o \
                     test_util.o \
                     test_write.o \
                     test_xrtm.o

OBJECTS  = $(OBJECTS_XRTM_INPUT) $(OBJECTS_XRTM_TEST)

BINARIES = callxrtm testxrtm

PRODUCTS =

all: libxrtm_input.a $(BINARIES)

libxrtm_input.a: $(OBJECTS_XRTM_INPUT)
	ar -rs libxrtm_input.a $(OBJECTS_XRTM_INPUT)

callxrtm: callxrtm.o libxrtm_input.a ../src/libxrtm.a ../misc/libxrtm_misc.a
	$(CFMIX_LD) $(CFMIX_LDFLAGS) -o callxrtm callxrtm.o libxrtm_input.a ../src/libxrtm.a ../misc/libxrtm_misc.a \
        $(LIBDIRS) $(LINKS)

testxrtm: testxrtm.o ${OBJECTS_XRTM_TEST} libxrtm_input.a ../src/libxrtm.a ../misc/libxrtm_misc.a
	$(CFMIX_LD) $(CFMIX_LDFLAGS) -o testxrtm testxrtm.o ${OBJECTS_XRTM_TEST} libxrtm_input.a ../src/libxrtm.a ../misc/libxrtm_misc.a \
        $(LIBDIRS) $(LINKS)

input_yylex.o: input_yylex.c
	$(CC) $(CCFLAGS) -Wno-unused-function -c $(INCDIRS) -o input_yylex.o input_yylex.c

input_yylex.c input_yylex.h: input_yylex.l
	flex --header-file=input_yylex.h -oinput_yylex.c input_yylex.l


include depend.inc
