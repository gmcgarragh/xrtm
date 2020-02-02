include ../make.inc
include ../common.inc

SUBDIRS =

OBJECTS = dasymtx.o           \
          pentdag.o           \
          xrtm_adding.o       \
          xrtm_adding_a.o     \
          xrtm_brdf.o         \
          xrtm_brdf_a.o       \
          xrtm_derivs.o       \
          xrtm_doub_rts.o     \
          xrtm_doubling.o     \
          xrtm_eig_bvp.o      \
          xrtm_eig_bvp_a.o    \
          xrtm_eig_rts.o      \
          xrtm_eig_rts_a.o    \
          xrtm_eig_util.o     \
          xrtm_eig_util_a.o   \
          xrtm_fd_interface.o \
          xrtm_four_stream.o  \
          xrtm_interface.o    \
          xrtm_matrix.o       \
          xrtm_mem_bvp.o      \
          xrtm_model.o        \
          xrtm_model_a.o      \
          xrtm_pade_rts.o     \
          xrtm_radiance.o     \
          xrtm_radiance_a.o   \
          xrtm_save_tree.o    \
          xrtm_scatter.o      \
          xrtm_scatter_a.o    \
          xrtm_sfi.o          \
          xrtm_single.o       \
          xrtm_single_a.o     \
          xrtm_six_stream.o   \
          xrtm_sos.o          \
          xrtm_source.o       \
          xrtm_source_a.o     \
          xrtm_stacks.o       \
          xrtm_support.o      \
          xrtm_two_os.o       \
          xrtm_two_stream.o   \
          xrtm_utility.o      \
          xrtm_utility_a.o    \
          xrtm_work.o         \
          xrtm_external.o     \
          rdi1mach.o          \
          version.o

ifdef make_dev_build
OBJECTS += ../src2/dsyev2.o      \
           ../src2/rs2.o         \
           ../src2/ext_2stream.o \
           ../src2/ext_lidort.o  \
           ../src2/ext_lrad.o    \
           ../src2/ext_polrad.o  \
           ../src2/ext_radiant.o \
           ../src2/ext_soi.o     \
           ../src2/ext_vlidort.o
endif

BINARIES =

PRODUCTS = version.c

all: libxrtm.a libxrtm.so

libxrtm.a: $(OBJECTS)
	ar -rs libxrtm.a $(OBJECTS)

libxrtm.so: $(OBJECTS)
	gcc -shared -Wl,-soname,libxrtm.so.1 -o libxrtm.so.1.0 $(OBJECTS)

version.o: version.c version.h

# Requires UNIX/LINUX date command and git
version.c:
	echo "#include \"version.h\"" > version.c; \
        echo "" >> version.c; \
        echo "const char *build_date  = \"`date +"%F %X %z"`\";" >> version.c; \
        echo "const char *build_sha_1 = \"`git describe --abbrev=7 --dirty --always --tags`\";" >> version.c

include depend.inc
