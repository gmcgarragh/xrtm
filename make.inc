#*******************************************************************************
# You should only need to modify the first two sections on this file
#*******************************************************************************


#*******************************************************************************
# Compiler and linker settings
#*******************************************************************************
CC            = gcc
CXX           = g++
F77           = gfortran
F90           = gfortran

CCFLAGS       = -O2 -Wall -Werror -Wno-uninitialized
CXXFLAGS      = $(CCFLAGS)
F77FLAGS      = -O2 -Wall -Werror -Wno-uninitialized -fno-second-underscore \
                -ffixed-line-length-none
F90FLAGS      = -O2 -Wall -Werror -Wno-uninitialized -fno-second-underscore \
                -ffree-line-length-0

# Link program that should be used for a mixed C/Fortran 
CFMIX_LD      = $(F90)
CFMIX_LDFLAGS = 


#*******************************************************************************
# BLAS and LAPACK settings.  This includes whatever needs to be on the command
# line to link to your BLAS and LAPACK.  If these are installed in standard
# places (/usr/include, /usr/lib) then these variable definitions may be left
# blank.
#*******************************************************************************
LIB_BLAS      = $(HOME)/opt/intel/mkl/10.2.2.025/lib/em64t/libmkl_intel_lp64.a \
                $(HOME)/opt/intel/mkl/10.2.2.025/lib/em64t/libmkl_intel_thread.a \
                $(HOME)/opt/intel/mkl/10.2.2.025/lib/em64t/libmkl_core.a \
                $(HOME)/opt/intel/mkl/10.2.2.025/lib/em64t/libguide.a \
                -lpthread
LIB_LAPACK    = $(HOME)/opt/lapack/liblapack.a


#*******************************************************************************
# 
#*******************************************************************************
VERSION       = 0.91


#*******************************************************************************
# 
#*******************************************************************************
CCDEFINES     = -D'VERSION="$(VERSION)"'

INCDIRS       = -I../src -I../interfaces -I../misc
LIBDIRS       = 

LINKS         = $(LIB_LAPACK) $(LIB_BLAS)


#*******************************************************************************
# 
#*******************************************************************************
MAKE_SUBDIRS  = + \
                @target=`echo $@ | sed s/-subdirs//`; \
                list='$(SUBDIRS)'; \
                for subdir in $$list; do \
                     echo "Making $$target in $$subdir"; \
                     (cd $$subdir && $(MAKE) $$target) || exit 1; \
                done;


#*******************************************************************************
# 
#*******************************************************************************
all: all-subdirs

all-subdirs:
	$(MAKE_SUBDIRS)


clean: clean-subdirs
	rm -f *.a *.o *.mod $(OBJECTS) $(BINARIES) $(PRODUCTS)

clean-subdirs:
	$(MAKE_SUBDIRS)


#*******************************************************************************
# 
#*******************************************************************************
.SUFFIXES: .c .f .f90 .o

.c.o:
	$(CC) $(CCFLAGS) $(CCDEFINES) -c $(INCDIRS) -o $*.o $<

.cpp.o:
	$(CXX) $(CXXFLAGS) $(CCDEFINES) -c $(INCDIRS) -o $*.o $<

.f.o:
	$(F77) $(F77FLAGS) -c $(INCDIRS) -o $*.o $<

.f90.o:
	$(F90) $(F90FLAGS) -c $(F90DEFINES) $(INCDIRS) -o $*.o $<
