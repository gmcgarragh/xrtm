include ../make.inc

F90FLAGS += -fno-underscoring

BINARIES = example_c example_cpp example_f77 example_f90

all: $(BINARIES)

example_c: example_c.o ../src/libxrtm.a ../misc/libxrtm_misc.a
	$(CFMIX_LD) $(CFMIX_LDFLAGS) -o example_c example_c.o ../src/libxrtm.a ../misc/libxrtm_misc.a \
        $(LIBDIRS) $(INCDIRS) $(LINKS) $(CFMIX_LDLINKS)

example_cpp: example_cpp.o ../interfaces/libxrtm_interfaces.a ../src/libxrtm.a ../misc/libxrtm_misc.a
	$(CFMIX_LD) $(CFMIX_LDFLAGS) -o example_cpp example_cpp.o ../interfaces/libxrtm_interfaces.a ../src/libxrtm.a ../misc/libxrtm_misc.a \
        $(LIBDIRS) $(INCDIRS) $(LINKS) $(CFMIX_LDLINKS) -lstdc++

example_f77: example_f77.o ../interfaces/libxrtm_interfaces.a ../src/libxrtm.a ../misc/libxrtm_misc.a
	$(CFMIX_LD) $(CFMIX_LDFLAGS) -o example_f77 example_f77.o ../interfaces/libxrtm_interfaces.a ../src/libxrtm.a ../misc/libxrtm_misc.a \
        $(LIBDIRS) $(INCDIRS) $(LINKS) $(CFMIX_LDLINKS)

example_f90: example_f90.o ../interfaces/libxrtm_interfaces.a ../src/libxrtm.a ../misc/libxrtm_misc.a
	$(CFMIX_LD) $(CFMIX_LDFLAGS) -o example_f90 example_f90.o ../interfaces/libxrtm_interfaces.a ../src/libxrtm.a ../misc/libxrtm_misc.a \
        $(LIBDIRS) $(INCDIRS) $(LINKS) $(CFMIX_LDLINKS)


include dep.inc
