include ../make.inc

OBJECTS = gmath_matrix.o   \
          gmath_vector.o   \
          gstruct.o        \
          gutil_alloc.o    \
          gutil_error.o    \
          gutil_parse.o    \
          rtutil_math.o    \
          rtutil_scat_io.o \
          rtutil_support.o

all: libxrtm_misc.a

libxrtm_misc.a: $(OBJECTS)
	ar -rs libxrtm_misc.a $(OBJECTS)

include dep.inc
