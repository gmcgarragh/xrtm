include ../make.inc
include ../common.inc

SUBDIRS =

OBJECTS = gindex_name_value.o \
          glist.o             \
          gmath_mat2d3d.o     \
          gmath_matrix.o      \
          gmath_vector.o      \
          gstruct.o           \
          gutil_alloc.o       \
          gutil_cmp.o         \
          gutil_error.o       \
          gutil_parse.o       \
          gutil_sort.o        \
          gutil_string.o      \
          rtutil_math.o       \
          rtutil_scat_fft.o   \
          rtutil_scat_io.o    \
          rtutil_support.o    \
          zeroin.o

BINARIES =

PRODUCTS =

all: libxrtm_misc.a

libxrtm_misc.a: $(OBJECTS)
	ar -rs libxrtm_misc.a $(OBJECTS)

include depend.inc
