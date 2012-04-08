include ../make.inc

OBJECTS = dasymtx.o         \
          xrtm_adding.o     \
          xrtm_adding_a.o   \
          xrtm_brdf.o       \
          xrtm_brdf_a.o     \
          xrtm_doub_rts.o   \
          xrtm_doubling.o   \
          xrtm_eig_bvp.o    \
          xrtm_eig_bvp_a.o  \
          xrtm_eig_rts.o    \
          xrtm_eig_rts_a.o  \
          xrtm_eig_util.o   \
          xrtm_eig_util_a.o \
          xrtm_fd.o         \
          xrtm_interface.o  \
          xrtm_matrix.o     \
          xrtm_mem_bvp.o    \
          xrtm_model.o      \
          xrtm_model_a.o    \
          xrtm_pade_rts.o   \
          xrtm_radiance.o   \
          xrtm_radiance_a.o \
          xrtm_save_tree.o  \
          xrtm_scatter.o    \
          xrtm_scatter_a.o  \
          xrtm_single.o     \
          xrtm_single_a.o   \
          xrtm_sos.o        \
          xrtm_source.o     \
          xrtm_source_a.o   \
          xrtm_stacks.o     \
          xrtm_support.o    \
          xrtm_two_os.o     \
          xrtm_utility.o    \
          xrtm_utility_a.o  \
          xrtm_work.o       \
          xrtm_external.o   \
          rdi1mach.o

all: libxrtm.a

libxrtm.a: $(OBJECTS)
	ar -rs libxrtm.a $(OBJECTS)

include dep.inc
