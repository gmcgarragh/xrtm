#*******************************************************************************
#
#*******************************************************************************
CCDEFINES += -DGBYTE_ORDER=BYTE_ORDER_LE -DNATIVE_SIZE=NATIVE_SIZE_64_BIT \
             -DPLATFORM=LINUX_GCC -DUNIX


MAKE_SUBDIRS  = + \
                @target=`echo $@ | sed s/-subdirs//`; \
                list='$(SUBDIRS)'; \
                for subdir in $$list; do \
                     echo "Making $$target in $$subdir"; \
                     (cd $$subdir && $(MAKE) $$target) || exit 1; \
                done;


all: all-subdirs

all-subdirs:
	$(MAKE_SUBDIRS)


clean: clean-subdirs
	rm -f *.a *.mod *.o *.so lib*.so.* $(OBJECTS) $(BINARIES) $(PRODUCTS) $(EXTRA_CLEANS)

clean-subdirs:
	$(MAKE_SUBDIRS)


README: readme_source.txt
	fold --spaces --width=80 readme_source.txt > README
	sed -i 's/[ \t]*$$//' README


# Requires gcc
depend: depend-subdirs
	@files=`find . -maxdepth 1 -name "*.c" | sort`; \
        if test $${#files} != 0; then \
             if (eval gcc -v 1> /dev/null 2>&1); then \
                  echo gcc -MM -w ${CCDEFINES} $$files "> depend.inc"; \
                  gcc -MM -w ${CCDEFINES} $$files > depend.inc; \
             else \
                  echo makedepend -f- -Y -- ${CCDEFINES} -- $$files "> depend.inc"; \
                  makedepend -f- -Y -- ${CCDEFINES} -- $$files > depend.inc; \
             fi \
        else \
             echo -n '' > depend.inc; \
        fi
	@files=`find . -maxdepth 1 -name "*.cpp" | sort`; \
        if test $${#files} != 0; then \
             if (eval gcc -v 1> /dev/null 2>&1); then \
                  echo gcc -MM -std=c++0x -w ${CCDEFINES} $$files ">> depend.inc"; \
                  gcc -MM -std=c++0x -w ${CCDEFINES} $$files >> depend.inc; \
             else \
                  echo makedepend -f- -Y -- ${CCDEFINES} -- $$files ">> depend.inc"; \
                  makedepend -f- -Y -- ${CCDEFINES} -- $$files >> depend.inc; \
             fi \
        else \
             echo -n '' >> depend.inc; \
        fi


# Requires cproto
depend-subdirs:
	${MAKE_SUBDIRS}

proto_types: proto_types-subdirs
	@if [ -e prototypes ]; then \
             files='$(OBJECTS)'; \
	     for object in $$files; do \
                  source=`echo $$object | sed 's|\.o|\.c|'`; \
                  if [ -e $$source ]; then \
                       header=`echo $$object | sed 's|\.o|_p\.h|'`; \
                       echo "cproto ${CCDEFINES} ${INCDIRS} $$source > prototypes/$$header"; \
                       cproto ${CCDEFINES} -DCPROTO ${INCDIRS} $$source > prototypes/$$header; \
                  fi \
	     done \
        fi

proto_types-subdirs:
	${MAKE_SUBDIRS}


.SUFFIXES: .c .cpp .f .f90 .o

.c.o:
	$(CC) $(CCFLAGS) -c $(CCDEFINES) $(INCDIRS) -o $*.o $<
.cpp.o:
	$(CXX) $(CXXFLAGS) -c $(CCDEFINES) $(INCDIRS) -o $*.o $<
.f.o:
	$(F77) $(F77FLAGS) -c $(INCDIRS) -o $*.o $<
.f90.o:
	$(F90) $(F90FLAGS) -c $(F90DEFINES) $(INCDIRS) -o $*.o $<
%.o : %.mod


clean_doc: clean_doc-subdirs

clean_doc-subdirs:
	$(MAKE_SUBDIRS)
