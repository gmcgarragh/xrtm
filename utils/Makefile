include ../make.inc

BINARIES = callxrtm testxrtm

OBJECTS = input_util.o \
          input_yylex.o

all: libxrtm_input.a $(BINARIES)

libxrtm_input.a: $(OBJECTS)
	ar -rs libxrtm_input.a $(OBJECTS)

callxrtm: callxrtm.o libxrtm_input.a ../src/libxrtm.a ../misc/libxrtm_misc.a
	$(CFMIX_LD) $(CFMIX_LDFLAGS) -o callxrtm callxrtm.o libxrtm_input.a ../src/libxrtm.a ../misc/libxrtm_misc.a \
        $(LIBDIRS) $(LINKS)

testxrtm: testxrtm.o test_core.o test_errors.o test_execute.o test_util.o libxrtm_input.a ../src/libxrtm.a ../misc/libxrtm_misc.a
	$(CFMIX_LD) $(CFMIX_LDFLAGS) -o testxrtm testxrtm.o test_core.o test_errors.o test_execute.o test_util.o libxrtm_input.a ../src/libxrtm.a ../misc/libxrtm_misc.a \
        $(LIBDIRS) $(LINKS)

input_yylex.o: input_yylex.c
	$(CC) $(CCFLAGS) -Wno-unused-function -c $(INCDIRS) -o input_yylex.o input_yylex.c

input_yylex.c input_yylex.h: input_yylex.l
	flex --header-file=input_yylex.h -oinput_yylex.c input_yylex.l


include dep.inc
