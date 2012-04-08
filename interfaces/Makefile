include ../make.inc

OBJECTS = xrtm_int_cpp.o      \
          xrtm_int_f77.o      \
          xrtm_int_f77_util.o \
          xrtm_int_f90.o      \
          xrtm_int_f90_util.o

all: libxrtm_interfaces.a

libxrtm_interfaces.a: $(OBJECTS)
	ar -rs libxrtm_interfaces.a $(OBJECTS)

include dep.inc
