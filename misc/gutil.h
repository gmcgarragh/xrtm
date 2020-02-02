/*******************************************************************************
**
**    Copyright (C) 1998-2020 Greg McGarragh <greg.mcgarragh@colostate.edu>
**
**    This source code is licensed under the GNU General Public License (GPL),
**    Version 3.  See the file COPYING for more details.
**
*******************************************************************************/

#ifndef GUTIL_H
#define GUTIL_H

#include <assert.h>
#ifdef __cplusplus
#include <complex>
#else
#include <complex.h>
#ifdef I
#undef I
#endif
#endif
#include <ctype.h>
#include <errno.h>
#include <float.h>
#include <limits.h>
#include <locale.h>
#include <math.h>
#include <setjmp.h>
#include <signal.h>
#include <stdarg.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#ifdef UNIX
#include <sys/stat.h>
#include <sys/types.h>
#endif

#ifdef __cplusplus
extern "C" {
#endif


/* Byte Order */
#define BYTE_ORDER_LE	0
#define BYTE_ORDER_BE	1

#define BYTE_ORDER_END	2


/* Native Size */
#define NATIVE_SIZE_16_BIT	1
#define NATIVE_SIZE_32_BIT	2
#define NATIVE_SIZE_64_BIT	3


/* Platforms */
#define AIX_CC		1
#define AIX_GCC		2
#define CYGWIN_GCC      3
#define DEC_UNIX_CC	4
#define DEC_UNIX_GCC	5
#define HP_CC		6
#define HP_GCC		7
#define IRIX_CC		8
#define IRIX_GCC	9
#define LINUX_GCC	10
#define LINUX_ICC	11
#define LINUX_PGCC	12
#define MACOSX_GCC	13
#define MACOSX_XLC	14
#define SOLARIS_CC	15
#define SOLARIS_GCC	16
#define WIN32_MSVC	17


/* Standard string sizes  */
#define NM	128
#define FN	256
#define LN	1024


/* Data structures stuff */
#define DS_STATIC	0
#define DS_DYNAMIC	1

#define MIN_QUEUE_SIZE	8
#define MIN_STACK_SIZE	8


/* Quicksort cutoff */
#define CUTOFF	3


/* Floating point coordinate precision */
#define COORD_PRECISION 0.0001


#if PLATFORM == WIN32_MSVC
#define finite		_finite
#define isinf		_isinf
#define isnan		_isnan
#define snprintf	_snprintf
#define strdup		_strdup

#define strtok_r(a, b, c) (strtok(a, b))
#endif


/* Data types */

/*
Type		ILP64	LP64	LLP64
char		8	8	8
short		16	16	16
int		64	32	32
long		64	64	32
long long	64	64	64
pointer		64	64	64
*/

typedef unsigned char		uchar;
typedef unsigned short		ushort;
typedef unsigned int		uint;
typedef unsigned long		ulong;
typedef unsigned long long	ulonglong;


#if PLATFORM == WIN32_MSVC
#define LONG_LONG signed __int64
#define ULONG_LONG unsigned __int64
#else
#define LONG_LONG signed long long
#define ULONG_LONG unsigned long long
#endif


#ifndef __cplusplus
typedef float _Complex 		fcomplex;
typedef double _Complex 	dcomplex;
#else
typedef std::complex<float> 	fcomplex;
typedef std::complex<double> 	dcomplex;

#define _Complex_I		std::complex<double>(0., 1.)

#define cabs			abs
#define cexp			exp
#define cpow			pow
#define creal			real
#define cimag			imag
#define csqrt			sqrt
#endif


typedef signed char		int8_t;
typedef signed short		int16_t;
typedef unsigned char		uint8_t;
typedef unsigned short 		uint16_t;
#if   (NATIVE_SIZE == NATIVE_SIZE_16_BIT)
typedef signed long		int32_t;
typedef unsigned long		uint32_t;
#elif (NATIVE_SIZE == NATIVE_SIZE_32_BIT)
typedef signed int		int32_t;
typedef unsigned int		uint32_t;

#if   (PLATFORM    == WIN32_MSVC)
typedef signed __int64		int64_t;
typedef unsigned __int64	uint64_t;
#else
typedef signed long long	int64_t;
typedef unsigned long long	uint64_t;
#endif

#else

typedef signed int		int32_t;
typedef unsigned int		uint32_t;
#if   (PLATFORM    == WIN32_MSVC)
typedef signed __int64		int64_t;
typedef unsigned __int64	uint64_t;
#else
typedef signed long		int64_t;
typedef unsigned long		uint64_t;
#endif

#endif

typedef float			float32_t;
typedef double			float64_t;


/* Data types ranges */

#ifndef INT8_MIN
#define INT8_MIN		CHAR_MIN
#define INT8_MAX		CHAR_MAX
#define UINT8_MAX		UCHAR_MAX
#endif

#ifndef INT16_MIN
#define INT16_MIN		SHRT_MIN
#define INT16_MAX		SHRT_MAX
#define UINT16_MAX		USHRT_MAX
#endif

#if   (NATIVE_SIZE == NATIVE_SIZE_16_BIT)
#ifndef INT32_MIN
#define INT32_MIN		LONG_MIN
#define INT32_MAX		LONG_MAX
#define UINT32_MAX		ULONG_MAX
#endif

#elif (NATIVE_SIZE == NATIVE_SIZE_32_BIT)
#ifndef INT32_MIN
#define INT32_MIN		INT_MIN
#define INT32_MAX		INT_MAX
#define UINT32_MAX		UINT_MAX
#endif

#if   (PLATFORM    == WIN32_MSVC)
#ifndef INT64_MIN
#define INT64_MAX		9223372036854775807i64
#define INT64_MIN		(-INT64_MAX - 1i64)
#define UINT64_MAX		18446744073709551615ui64
#endif
#else
#ifndef INT64_MIN
#define INT64_MAX		9223372036854775807LL
#define INT64_MIN		(-INT64_MAX - 1LL)
#define UINT64_MAX		18446744073709551615ULL
#endif
#endif

#else

#ifndef INT32_MIN
#define INT32_MIN		INT_MIN
#define INT32_MAX		INT_MAX
#define UINT32_MAX		UINT_MAX
#endif
#if   (PLATFORM    == WIN32_MSVC)
#ifndef INT64_MIN
#define INT64_MIN		_I64_MIN
#define INT64_MAX		_I64_MAX
#define UINT64_MAX		_UI64_MAX
#endif
#else
#ifndef INT64_MIN
#define INT64_MIN		LONG_MIN
#define INT64_MAX		LONG_MAX
#define UINT64_MAX		ULONG_MAX
#endif
#endif

#endif

#define FLOAT32_MIN		-FLT_MAX
#define FLOAT32_MAX		FLT_MAX
#define FLOAT64_MIN		-DBL_MAX
#define FLOAT64_MAX		DBL_MAX


/* Concatenation macros */
#ifdef CAT
#undef CAT
#endif
#define CAT(a, b) a##b

#ifdef XCAT
#undef XCAT
#endif
#define XCAT(x, y) CAT(x, y)


#ifdef ABS
#undef ABS
#endif
#define ABS(a) ((a) < (0) ? (-(a)) : (a))


/* Min/Max macros */
#ifdef MIN
#undef MIN
#endif
#define MIN(a,b) ((a) < (b) ? (a) : (b))

#ifdef MAX
#undef MAX
#endif
#define MAX(a,b) ((a) > (b) ? (a) : (b))


/* Array length macro */
#define LENGTH(NAME) (sizeof(NAME) / sizeof(NAME[0]))


/* Debug Macros */
#define debug_info(TEXT) do {					\
    fprintf(stderr, "ERROR: file = %s, line = %d: %s\n",	\
            __FILE__, __LINE__, TEXT);				\
} while (0)


#define debug_info_exit(TEXT) do {				\
    fprintf(stderr, "ERROR: file = %s, line = %d: %s\n",	\
            __FILE__, __LINE__, TEXT);				\
    exit(1);							\
} while (0)


/* Powers of two */
#define N0_16 0
#define N1_16 1

#if (NATIVE_SIZE == NATIVE_SIZE_64_BIT)
#define N0_32 0
#define N1_32 1
#else
#define N0_32 0UL
#define N1_32 1UL
#endif

#if (NATIVE_SIZE == NATIVE_SIZE_64_BIT)
#if (PLATFORM == WIN32_MSVC)
#define N0_64 0ui64
#define N1_64 1ui64
#else
#define N0_64 0UL
#define N1_64 1UL
#endif
#else
#if (PLATFORM == WIN32_MSVC)
#define N0_64 0ui64
#define N1_64 1ui64
#else
#define N0_64 0ULL
#define N1_64 1ULL
#endif

#endif

#define TWO00	(N1_16<<0)
#define TWO01	(N1_16<<1)
#define TWO02	(N1_16<<2)
#define TWO03	(N1_16<<3)
#define TWO04	(N1_16<<4)
#define TWO05	(N1_16<<5)
#define TWO06	(N1_16<<6)
#define TWO07	(N1_16<<7)
#define TWO08	(N1_16<<8)
#define TWO09	(N1_16<<9)
#define TWO10	(N1_16<<10)
#define TWO11	(N1_16<<11)
#define TWO12	(N1_16<<12)
#define TWO13	(N1_16<<13)
#define TWO14	(N1_16<<14)
#define TWO15	(N1_16<<15)

#define TWO16m1	(~N0_16)

#define TWO16	(N1_32<<16)
#define TWO17	(N1_32<<17)
#define TWO18	(N1_32<<18)
#define TWO19	(N1_32<<19)
#define TWO20	(N1_32<<20)
#define TWO21	(N1_32<<21)
#define TWO22	(N1_32<<22)
#define TWO23	(N1_32<<23)
#define TWO24	(N1_32<<24)
#define TWO25	(N1_32<<25)
#define TWO26	(N1_32<<26)
#define TWO27	(N1_32<<27)
#define TWO28	(N1_32<<28)
#define TWO29	(N1_32<<29)
#define TWO30	(N1_32<<30)
#define TWO31	(N1_32<<31)

#define TWO32m1	(~N0_32)

#define TWO32   (N1_64<<32)
#define TWO33   (N1_64<<33)
#define TWO34   (N1_64<<34)
#define TWO35   (N1_64<<35)
#define TWO36   (N1_64<<36)
#define TWO37   (N1_64<<37)
#define TWO38   (N1_64<<38)
#define TWO39   (N1_64<<39)
#define TWO40   (N1_64<<40)
#define TWO41   (N1_64<<41)
#define TWO42   (N1_64<<42)
#define TWO43   (N1_64<<43)
#define TWO44   (N1_64<<44)
#define TWO45   (N1_64<<45)
#define TWO46   (N1_64<<46)
#define TWO47   (N1_64<<47)
#define TWO48   (N1_64<<48)
#define TWO49   (N1_64<<49)
#define TWO50   (N1_64<<50)
#define TWO51   (N1_64<<51)
#define TWO52   (N1_64<<52)
#define TWO53   (N1_64<<53)
#define TWO54   (N1_64<<54)
#define TWO55   (N1_64<<55)
#define TWO56   (N1_64<<56)
#define TWO57   (N1_64<<57)
#define TWO58   (N1_64<<58)
#define TWO59   (N1_64<<59)
#define TWO60   (N1_64<<60)
#define TWO61   (N1_64<<61)
#define TWO62   (N1_64<<62)
#define TWO63   (N1_64<<63)

#define TWO64m1	(~N0_64)


/* Math constants */
#ifdef M_PI
#define PI	M_PI
#else
#define PI	3.14159265358979323846264338327
#endif


/* Physical constants */

/* http://physics.nist.gov/cuu/Constants/ */

#define SPEED_OF_LIGHT			2.99792458e+08
#define GRAV_ACCEL_EARTH		9.80665			/* m / s^2	*/

#define AVOGADRO_CONSTANT		6.02214179e23		/* 1 / mol	*/
#define LOSCHMIDT_CONSTANT		2.6867774e25		/* 1 / m3	*/
#define UNIVERAL_GAS_CONSTANT		8.314472		/* J / K mol	*/

#define BOLTZMAN_CONSTANT		1.3806503e-23		/* J / K	*/
#define PLANCK_CONSTANT			6.62606876e-34		/* J s		*/
#define STEFAN_BOLTZMAN_CONSTANT	5.670400E-8		/* W m-2 K-4	*/


/* Derived constants */
#define MEAN_MOLAR_MASS_DRY_AIR		28.9644			/* g mol-1 */

/* HITRAN 2004 */
#define MEAN_MOLAR_MASS_H2O		1.801526e+01		/* g mol-1 */

#define SPEC_GAS_CONST_DRY_AIR		(UNIVERAL_GAS_CONSTANT / (MEAN_MOLAR_MASS_DRY_AIR * 1.e-3))	/* J / kg K	*/
#define SPEC_GAS_CONST_H2O		(UNIVERAL_GAS_CONSTANT / (MEAN_MOLAR_MASS_H2O     * 1.e-3))	/* J / kg K	*/


/* Unitary conversion macros */
#define D2R	(PI/180.0)
#define R2D	(180.0/PI)

#define FT2M	(12.0*2.54/100.0)
#define M2FT	(100/(12.0*2.54))


/* For read_line(), gutil_fio.c */
#define NOBLANKS	(1<<0)
#define COMMENTS	(1<<1)



/* Structure types */

typedef struct {
     const char **comments;
     uint options;
     int count;
     int size;
} read_line_data;

/*
typedef struct {
     const char **comments;
     char **mul_str;
     char **mul_end;
     uint options;
     int mul_count;
     int count;
     int size;
} read_line_data;
*/
typedef struct {
     char comment;
     int count;
     int newline_found;
} gparse_data;

typedef struct {
     int year;
     int month;
     double day;
     double hour;
     double min;
     double sec;
} gtime_data;


typedef struct {
     time_t time1;
     time_t time2;
} timing_data;


typedef struct {
     int offset;
     int size;
     int n;
     int size_lst;
     int perc_lst;
} progress_data;


typedef struct {
     int offset;
     int size;
     int n;
     int size_lst;
     int perc_lst;
     int current;
} progress_ll_data;



/* Function prototypes */

/* **** gutil_alloc.c **** */

size_t array_mem_length(int rank, size_t *dimen, size_t size);
void *alloc_array(int, size_t *, int);
void *realloc_array(void *, int, size_t *, int);
void *array_from_mem(void *, int, size_t *, int, int);
void free_array(void *, int);
void copy_array(void *, void *, int, size_t *, int);
void *dup_array(void *, int, size_t *, int);

size_t array_mem_length1(size_t m, size_t size);
void *alloc_array1(size_t, int);
void *realloc_array1(void *, size_t, int);
void *array_from_mem1(void *, size_t, int, int);
void free_array1(void *);
void copy_array1(void *, void *, size_t, int);
void *dup_array1(void *, size_t, int);

size_t array_mem_length2(size_t m, size_t n, size_t size);
void **alloc_array2(size_t, size_t, int);
void **realloc_array2(void **, size_t, size_t, int);
void **array_from_mem2(void *, size_t, size_t, int, int);
void free_array2(void **);
void copy_array2(void **, void **, size_t, size_t, int);
void **dup_array2(void **, size_t, size_t, int);

size_t array_mem_length3(size_t m, size_t n, size_t o, size_t size);
void ***alloc_array3(size_t, size_t, size_t, int);
void ***realloc_array3(void ***, size_t, size_t, size_t, int);
void ***array_from_mem3(void *, size_t, size_t, size_t, int, int);
void free_array3(void ***);
void copy_array3(void ***, void ***, size_t, size_t, size_t, int);
void ***dup_array3(void ***, size_t, size_t, size_t, int);

size_t array_mem_length4(size_t m, size_t n, size_t o, size_t p, size_t size);
void ****alloc_array4(size_t, size_t, size_t, size_t, int);
void ****realloc_array4(void ****, size_t, size_t, size_t, size_t, int);
void ****array_from_mem4(void *, size_t, size_t, size_t, size_t, int, int);
void free_array4(void ****);
void copy_array4(void ****, void ****, size_t, size_t, size_t, size_t, int);
void ****dup_array4(void ****, size_t, size_t, size_t, size_t, int);

size_t array_mem_length5(size_t m, size_t n, size_t o, size_t p, size_t q, size_t size);
void *****alloc_array5(size_t, size_t, size_t, size_t, size_t, int);
void *****realloc_array5(void *****, size_t, size_t, size_t, size_t, size_t, int);
void *****array_from_mem5(void *, size_t, size_t, size_t, size_t, size_t, int, int);
void free_array5(void *****);
void copy_array5(void *****, void *****, size_t, size_t, size_t, size_t, size_t, int);
void *****dup_array5(void *****, size_t, size_t, size_t, size_t, size_t, int);

size_t array_mem_length6(size_t m, size_t n, size_t o, size_t p, size_t q, size_t r, size_t size);
void ******alloc_array6(size_t, size_t, size_t, size_t, size_t, size_t, int);
void ******realloc_array6(void ******, size_t, size_t, size_t, size_t, size_t, size_t, int);
void ******array_from_mem6(void *, size_t, size_t, size_t, size_t, size_t, size_t, int, int);
void free_array6(void ******);
void copy_array6(void ******, void ******, size_t, size_t, size_t, size_t, size_t, size_t, int);
void ******dup_array6(void ******, size_t, size_t, size_t, size_t, size_t, size_t, int);

size_t array_mem_length7(size_t m, size_t n, size_t o, size_t p, size_t q, size_t r, size_t s, size_t size);
void *******alloc_array7(size_t m, size_t n, size_t o, size_t p, size_t q, size_t r, size_t s, int size);
void *******realloc_array7(void *******a, size_t m, size_t n, size_t o, size_t p, size_t q, size_t r, size_t s, int size);
void *******array_from_mem7(void *a, size_t m, size_t n, size_t o, size_t p, size_t q, size_t r, size_t s, int size, int flag);
void free_array7(void *******a);
void copy_array7(void *******, void ******, size_t, size_t, size_t, size_t, size_t, size_t, size_t, int);
void ******dup_array7(void *******, size_t, size_t, size_t, size_t, size_t, size_t, size_t, int);

#define all_alloc_proto(postfix_, type_)							\
type_ *XCAT(alloc_array1_, postfix_)(size_t);							\
type_ *XCAT(realloc_array1_, postfix_)(type_ *, size_t);						\
type_ *XCAT(array_from_mem1_, postfix_)(type_ *, size_t);						\
void XCAT(free_array1_, postfix_)(type_ *);							\
void XCAT(init_array1_, postfix_)(type_ *, size_t, type_);					\
void XCAT(copy_array1_, postfix_)(type_ *, const type_ *, size_t);				\
type_ *XCAT(dup_array1_, postfix_)(type_ *, size_t);						\
\
type_ **XCAT(alloc_array2_, postfix_)(size_t, size_t);						\
type_ **XCAT(realloc_array2_, postfix_)(type_ **, size_t, size_t);					\
type_ **XCAT(array_from_mem2_, postfix_)(type_ *, size_t, size_t);					\
void XCAT(free_array2_, postfix_)(type_ **);							\
void XCAT(init_array2_, postfix_)(type_ **, size_t, size_t, type_);					\
void XCAT(copy_array2_, postfix_)(type_ **, type_ **, size_t, size_t);				\
type_ **XCAT(dup_array2_, postfix_)(type_ **, size_t, size_t);					\
\
type_ ***XCAT(alloc_array3_, postfix_)(size_t, size_t, size_t);					\
type_ ***XCAT(realloc_array3_, postfix_)(type_ ***, size_t, size_t, size_t);				\
type_ ***XCAT(array_from_mem3_, postfix_)(type_ *, size_t, size_t, size_t);				\
void XCAT(free_array3_, postfix_)(type_ ***);							\
void XCAT(init_array3_, postfix_)(type_ ***, size_t, size_t, size_t, type_);				\
void XCAT(copy_array3_, postfix_)(type_ ***, type_ ***, size_t, size_t, size_t);			\
type_ ***XCAT(dup_array3_, postfix_)(type_ ***, size_t, size_t, size_t);				\
\
type_ ****XCAT(alloc_array4_, postfix_)(size_t, size_t, size_t, size_t);				\
type_ ****XCAT(realloc_array4_, postfix_)(type_ ****, size_t, size_t, size_t, size_t);			\
type_ ****XCAT(array_from_mem4_, postfix_)(type_ *, size_t, size_t, size_t, size_t);			\
void XCAT(free_array4_, postfix_)(type_ ****);							\
void XCAT(init_array4_, postfix_)(type_ ****, size_t, size_t, size_t, size_t, type_);			\
void XCAT(copy_array4_, postfix_)(type_ ****, type_ ****, size_t, size_t, size_t, size_t);		\
type_ ****XCAT(dup_array4_, postfix_)(type_ ****, size_t, size_t, size_t, size_t);			\
\
type_ *****XCAT(alloc_array5_, postfix_)(size_t, size_t, size_t, size_t, size_t);				\
type_ *****XCAT(realloc_array5_, postfix_)(type_ *****, size_t, size_t, size_t, size_t, size_t);		\
type_ *****XCAT(array_from_mem5_, postfix_)(type_ *, size_t, size_t, size_t, size_t, size_t);		\
void XCAT(free_array5_, postfix_)(type_ *****);							\
void XCAT(init_array5_, postfix_)(type_ *****, size_t, size_t, size_t, size_t, size_t, type_);		\
void XCAT(copy_array5_, postfix_)(type_ *****, type_ *****, size_t, size_t, size_t, size_t, size_t);	\
type_ *****XCAT(dup_array5_, postfix_)(type_ *****, size_t, size_t, size_t, size_t, size_t);		\
\
type_ ******XCAT(alloc_array6_, postfix_)(size_t, size_t, size_t, size_t, size_t, size_t);			\
type_ ******XCAT(realloc_array6_, postfix_)(type_ ******, size_t, size_t, size_t, size_t, size_t, size_t);	\
type_ ******XCAT(array_from_mem6_, postfix_)(type_ *, size_t, size_t, size_t, size_t, size_t, size_t);	\
void XCAT(free_array6_, postfix_)(type_ ******);						\
void XCAT(init_array6_, postfix_)(type_ ******, size_t, size_t, size_t, size_t, size_t, size_t, type_);	\
void XCAT(copy_array6_, postfix_)(type_ ******, type_ ******, size_t, size_t, size_t, size_t, size_t, size_t);		\
type_ ******XCAT(dup_array6_, postfix_)(type_ ******, size_t, size_t, size_t, size_t, size_t, size_t);			\
\
type_ *******XCAT(alloc_array7_, postfix_)(size_t, size_t, size_t, size_t, size_t, size_t, size_t);				\
type_ *******XCAT(realloc_array7_, postfix_)(type_ *******, size_t, size_t, size_t, size_t, size_t, size_t, size_t);		\
type_ *******XCAT(array_from_mem7_, postfix_)(type_ *, size_t, size_t, size_t, size_t, size_t, size_t, size_t);		\
void XCAT(free_array7_, postfix_)(type_ *******);								\
void XCAT(init_array7_, postfix_)(type_ *******, size_t, size_t, size_t, size_t, size_t, size_t, size_t, type_);		\
void XCAT(copy_array7_, postfix_)(type_ *******, type_ *******, size_t, size_t, size_t, size_t, size_t, size_t, size_t);	\
type_ *******XCAT(dup_array7_, postfix_)(type_ *******, size_t, size_t, size_t, size_t, size_t, size_t, size_t);		\

all_alloc_proto(uc,  uchar)
all_alloc_proto(c,   char)
all_alloc_proto(us,  ushort)
all_alloc_proto(s,   short)
all_alloc_proto(ui,  uint)
all_alloc_proto(i,   int)
all_alloc_proto(ul,  ulong)
all_alloc_proto(l,   long)
all_alloc_proto(ull, uint64_t)
all_alloc_proto(ll,  int64_t)
all_alloc_proto(f,   float)
all_alloc_proto(fc,  fcomplex)
all_alloc_proto(d,   double)
all_alloc_proto(dc,  dcomplex)

void nr_to_c_darray2(double ***, int);
void c_to_nr_darray2(double ***, int);


/* **** gutil_bio.c **** */
/*
void read_line_data_init(read_line_data *);
*/
int bgets(char *, int, const char *);
char *bread_line(const char *, char *, char **, read_line_data *);
char **bread_lines(char *, char **, char **, int, int, read_line_data *);
char **bread_lines_all(char *, char **, char **, int, read_line_data *);
char *bskip_lines(char *, char *, char **, int, read_line_data *);
int bcount_lines(char *, char *, char **, read_line_data *);

size_t bread_order(void *, size_t, size_t, uchar *, int);
size_t bwrite_order(void *, size_t, size_t, uchar *, int);


/* **** gutil_cmp.c **** */

#define all_cmp_proto(postfix_, type_)			\
int XCAT(cmp_, postfix_)(type_, type_);			\
int XCAT(cmp_p_, postfix_)(type_ *, type_ *);		\
int XCAT(cmp_n_, postfix_)(type_ *, type_ *, long);

all_cmp_proto(uc,  uchar)
all_cmp_proto(c,   char)
all_cmp_proto(us,  ushort)
all_cmp_proto(s,   short)
all_cmp_proto(ui,  uint)
all_cmp_proto(i,   int)
all_cmp_proto(ul,  ulong)
all_cmp_proto(l,   long)
all_cmp_proto(ull, uint64_t)
all_cmp_proto(ll,  int64_t)
all_cmp_proto(f,   float)
all_cmp_proto(d,   double)


#define all_cmp_dec_proto(postfix_, type_)			\
int XCAT(cmp_dec_, postfix_)(type_, type_);			\
int XCAT(cmp_dec_p_, postfix_)(type_ *, type_ *);		\
int XCAT(cmp_dec_n_, postfix_)(type_ *, type_ *, long);

all_cmp_dec_proto(uc,  uchar)
all_cmp_dec_proto(c,   char)
all_cmp_dec_proto(us,  ushort)
all_cmp_dec_proto(s,   short)
all_cmp_dec_proto(ui,  uint)
all_cmp_dec_proto(i,   int)
all_cmp_dec_proto(ul,  ulong)
all_cmp_dec_proto(l,   long)
all_cmp_dec_proto(ull, uint64_t)
all_cmp_dec_proto(ll,  int64_t)
all_cmp_dec_proto(f,   float)
all_cmp_dec_proto(d,   double)


/* **** gutil_error.c **** */

void check_arg_count(int i, int argc, int n, const char *s);

int strtoi_errmsg(const char *string, const char *name, int *value);
int strtoi_errmsg_exit(const char *string, const char *name);
long strtol_errmsg(const char *string, const char *name, long *value);
long strtol_errmsg_exit(const char *string, const char *name);
double strtod_errmsg(const char *string, const char *name, double *value);
double strtod_errmsg_exit(const char *string, const char *name);

int system_errmsg(const char *, int);
void system_errmsg_exit(const char *, int);


/* **** gutil_file.c **** */

int fend(FILE *);
int fexists(const char *);
FILE *fopen_errmsg_exit(const char *, const char *);
FILE *freopen_errmsg_exit(const char *, const char *, FILE *);
int fflush_errmsg_exit(FILE *);
int fclose_errmsg_exit(FILE *);
int fseek_errmsg_exit(FILE *, int, int);
int fremove_errmsg_exit(const char *);
int frename_errmsg_exit(const char *, const char *);
int file_copy(FILE *old_fp, FILE *new_fp);
int copy_file_errmsg_exit(const char *, const char *);
int fread_errmsg_exit(void *, size_t, size_t, FILE *);
int fwrite_errmsg_exit(void *, size_t, size_t, FILE *);
#ifdef UNIX
int chdir_errmsg_exit(const char *);
int mkdir_errmsg_exit(const char *, mode_t);
#endif
int system_utils_match_file(const char *dir_name, const char *file_regex, char *file_name, int n);


/* **** gutil_fio.c **** */

void read_line_data_init(read_line_data *);

char *fread_line(FILE *, char *, read_line_data *);
char **fread_lines(FILE *, char **, int, int, read_line_data *);
char **fread_lines_all(FILE *, char **, int, read_line_data *);
char *fskip_lines(FILE *, char *, int, read_line_data *);
int fcount_lines(FILE *, char *, read_line_data *);

size_t fread_order(void *, size_t, size_t, FILE *, int);
size_t fwrite_order(void *, size_t, size_t, FILE *, int);


/* **** gutil_mem.c **** */

void swap_bits(void *, int, int);
void swap_bits1(void *, int);
void swap_bits2(void *, int);
void swap_bits4(void *, int);
void swap_bits8(void *, int);
void swap_bits16(void *, int);
void swap_bits32(void *, int);
void swap_bitsNN(void *, int, int);
void swap_bytes(void *, int, int);
void swap_bytes2(void *, int);
void swap_bytes4(void *, int);
void swap_bytes8(void *, int);
void swap_bytes16(void *, int);
void swap_bytes32(void *, int);
void swap_bytesNN(void *, long, long);
void swap_words(void *, long, long);
void copy_swap_bytes(void *, void *, int, int);
void copy_swap_bytes2(void *, void *, int, int);
void copy_swap_bytes3(void *, void *, int, int);
void copy_swap_bytes4(void *, void *, int, int);
void copy_swap_bytes8(void *, void *, int, int);
void copy_swap_bytes16(void *, void *, int, int);
void copy_swap_bytes32(void *, void *, int, int);
void copy_swap_bytesNN(void *, void *, int, int);
void copy_swap_words(void *, void *, int, int);


/* **** gutil_parse.c **** */

int gpass_newline(FILE *fp);
int gcount_lines(FILE *fp);
int gparse_comment(FILE *fp, char comment);
void gparse_comments(FILE *fp, char comment);
void trimm_string(char *s);

void gparse_init(gparse_data *d, char comment);
char gparse_string(gparse_data *d, FILE *fp, int n, char *s, const char *delims);
int gparse_newline(gparse_data *d, FILE *fp);
char gparse_char(gparse_data *d, FILE *fp, const char *delims);
char *gparse_quoted(gparse_data *d, FILE *fp, int n, char *s, const char *delims);
int gparse_int(gparse_data *d, FILE *fp, const char *delims);
double gparse_double(gparse_data *d, FILE *fp, const char *delims);


/* **** gutil_search.c **** */

long bin_search_inc_uc(const uchar *a, long n, long accel, uchar x);
long bin_search_inc_c(const char *a, long n, long accel, char x);
long bin_search_inc_us(const ushort *a, long n, long accel, ushort x);
long bin_search_inc_s(const short *a, long n, long accel, short x);
long bin_search_inc_ui(const uint *a, long n, long accel, uint x);
long bin_search_inc_i(const int *a, long n, long accel, int x);
long bin_search_inc_ul(const ulong *a, long n, long accel, ulong x);
long bin_search_inc_l(const long *a, long n, long accel, long x);
long bin_search_inc_ull(const uint64_t *a, long n, long accel, uint64_t x);
long bin_search_inc_ll(const int64_t *a, long n, long accel, int64_t x);
long bin_search_inc_f(const float *a, long n, long accel, float x);
long bin_search_inc_d(const double *a, long n, long accel, double x);
long bin_search_dec_uc(const uchar *a, long n, long accel, uchar x);
long bin_search_dec_c(const char *a, long n, long accel, char x);
long bin_search_dec_us(const ushort *a, long n, long accel, ushort x);
long bin_search_dec_s(const short *a, long n, long accel, short x);
long bin_search_dec_ui(const uint *a, long n, long accel, uint x);
long bin_search_dec_i(const int *a, long n, long accel, int x);
long bin_search_dec_ul(const ulong *a, long n, long accel, ulong x);
long bin_search_dec_l(const long *a, long n, long accel, long x);
long bin_search_dec_ull(const uint64_t *a, long n, long accel, uint64_t x);
long bin_search_dec_ll(const int64_t *a, long n, long accel, int64_t x);
long bin_search_dec_f(const float *a, long n, long accel, float x);
long bin_search_dec_d(const double *a, long n, long accel, double x);

long bin_search_data(const void *, long, const void *, long, int (*) (const void *, const void *));
long bin_search_index(const void *, long *, long, const void *, long, int (*) (const void *, const void *));
long bin_search_ptr(const void **, long, const void *, long, int (*) (const void *, const void *));

long bin_search_int_inc_uc(const uchar *a, long n, long accel, uchar x);
long bin_search_int_inc_c(const char *a, long n, long accel, char x);
long bin_search_int_inc_us(const ushort *a, long n, long accel, ushort x);
long bin_search_int_inc_s(const short *a, long n, long accel, short x);
long bin_search_int_inc_ui(const uint *a, long n, long accel, uint x);
long bin_search_int_inc_i(const int *a, long n, long accel, int x);
long bin_search_int_inc_ul(const ulong *a, long n, long accel, ulong x);
long bin_search_int_inc_l(const long *a, long n, long accel, long x);
long bin_search_int_inc_ull(const uint64_t *a, long n, long accel, uint64_t x);
long bin_search_int_inc_ll(const int64_t *a, long n, long accel, int64_t x);
long bin_search_int_inc_f(const float *a, long n, long accel, float x);
long bin_search_int_inc_d(const double *a, long n, long accel, double x);
long bin_search_int_dec_uc(const uchar *a, long n, long accel, uchar x);
long bin_search_int_dec_c(const char *a, long n, long accel, char x);
long bin_search_int_dec_us(const ushort *a, long n, long accel, ushort x);
long bin_search_int_dec_s(const short *a, long n, long accel, short x);
long bin_search_int_dec_ui(const uint *a, long n, long accel, uint x);
long bin_search_int_dec_i(const int *a, long n, long accel, int x);
long bin_search_int_dec_ul(const ulong *a, long n, long accel, ulong x);
long bin_search_int_dec_l(const long *a, long n, long accel, long x);
long bin_search_int_dec_ull(const uint64_t *a, long n, long accel, uint64_t x);
long bin_search_int_dec_ll(const int64_t *a, long n, long accel, int64_t x);
long bin_search_int_dec_f(const float *a, long n, long accel, float x);
long bin_search_int_dec_d(const double *a, long n, long accel, double x);

long bin_search_int_data(const void *a, long n, const void *x, long size, int (*comp) (const void *, const void *));
long bin_search_int_index(const void *a, long *index, long n, const void *x, long size, int (*comp) (const void *, const void *));
long bin_search_int_ptr(const void **a, long n, const void *x, long size, int (*comp) (const void *, const void *));


/* **** gutil_sort.c **** */

void swap_data(void *, void *, void *, size_t size);
void swap_index(int *, int *);
void swap_pointer(uchar **, uchar **);
int insertion_sort_data(void *, int, size_t, int (*) (const void *, const void *));
void insertion_sort_index(void *, int *, int, size_t, int (*) (const void *, const void *));
void insertion_sort_pointer(void **, int, size_t, int (*) (const void *, const void *));
void *median_3_data(void *, void *, int, int, size_t, int (*) (const void *, const void *));
void q_sort_data(void *, void *, int, int, size_t, int (*) (const void *, const void *));
int quick_sort_data(void *, int, size_t, int (*) (const void *, const void *));
void *median_3_index(void *, int *, int, int, size_t, int (*) (const void *, const void *));
void q_sort_index(void *, int *, int, int, size_t, int (*) (const void *, const void *));
void quick_sort_index(void *, int *, int, size_t, int (*) (const void *, const void *));
void **median_3_pointer(void **, int, int, size_t, int (*) (const void *, const void *));
void q_sort_pointer(void **, int, int, size_t, int (*) (const void *, const void *));
void quick_sort_pointer(void **, int, size_t, int (*) (const void *, const void *));


/* **** gutil_string.c **** */

void str_replace_char(char *s, char c1, char c2);
int gstrtob(const char *, int *);
const char *gbtostr(int);
int gstrtous(const char *, ushort *);
int gstrtos(const char *, short *);
int gstrtoui(const char *, uint *);
int gstrtoi(const char *, int *);
int gstrtoul(const char *, ulong *);
int gstrtol(const char *, long *);
int gstrtoull(const char *, uint64_t *);
int gstrtoll(const char *, int64_t *);
int gstrtof(const char *, float *);
int gstrtod(const char *, double *);
char *gtrtok_r(char *, const char *, char **);
char *get_token(int, char *, char **);
char *get_token_init(char *, int, char *, char **);
int count_tokens(char *, char *);
int tokenize(char *, char **, char *);
char *remove_pad(char *, int);
char *int_to_string(int, char *);
char *uint_to_string(uint, char *);
char *long_to_string(long, char *);
char *ulong_to_string(ulong, char *);
char *float_to_string(float, char *);
char *double_to_string(double, char *);
int strcmp2(const void *, const void *);
int strcmp3(const void *, const void *);
int mincmp(const void *, const void *);
int minncmp(const void *, const void *, int);
int ncscmp(const void *, const void *);
int ncsncmp(const void *, const void *, int);
void strtolower(const char *, char *);
void strtoupper(const char *, char *);
void extract_array_1d(uchar *, void *, int, int, int);
void extract_array_2d(uchar *, void **, int, int, int);
void extract_c_array_1d(uchar *, char *, int, int, int);
void extract_c_array_2d(uchar *, char **, int, int, int);
void extract_uc_array_1d(uchar *, uchar *, int, int, int);
void extract_uc_array_2d(uchar *, uchar **, int, int, int);
void extract_int(uchar *, int *, int, int, int);
void extract_long(uchar *, long *, int, int, int);
void extract_float(uchar *, float *, int, int, int);
void extract_double(uchar *, double *, int, int, int);
void insert_array_1d(uchar *, void *, int, int, int);
void insert_array_2d(uchar *, void **, int, int, int);
void insert_c_array_1d(uchar *s, char *, int j, int, int);
void insert_c_array_2d(uchar *s, char **, int, int, int);
void insert_uc_array_1d(uchar *s, uchar *, int, int, int);
void insert_uc_array_2d(uchar *s, uchar **, int, int, int);
void insert_int(uchar *, int *, int, int, int);
void insert_long(uchar *, long *, int, int, int);
void insert_float(uchar *, float *, int, int, int, int);
void insert_double(uchar *, double *, int, int, int, int);
void insert_scientific_e(uchar *, double *, int, int, int, int);
void insert_scientific_E(uchar *, double *, int, int, int, int);


/* **** gutil_struct.c **** */

void stdfree(void *);


/* **** gutil_time.c **** */

int get_abrv_month_num(const char *);
char *get_abrv_month_name(char *, int);
int get_full_month_num(const char *);
char *get_full_month_name(char *, int);
void jul_to_cal_date(long, int *, int *, int *);
long cal_to_jul_day(int, int, int);
void gtime_init(gtime_data *);
double gtime_to_jtime(gtime_data *);
void jtime_to_gtime(double jtime, gtime_data *);
void convert_gtime(gtime_data *, gtime_data *);
double calender_time_to_jtime(int year, int month, double day, double hour, double minute, double second);
void jtime_to_calender_time(double jtime, int *year, int *month, double *day, double *hour, double *minute, double *second);
void convert_calender_time(int year1, int month1, double day1, double hour1, double minute1, double second1, int *year2, int *month2, double *day2, double *hour2, double *minute2, double *second2);
int is_leap_year(int year);
int get_days_in_year(int);
int get_days_in_month(int, int);
int month_day_to_doy(int, int, int);
void doy_to_month_day(int, int, int *, int *);
long year_doy_to_jul(int, int);
void jul_to_year_doy(long, int *, int *);
long year_month_day_to_jul(int, int, int);
void jul_to_year_month_day(long, int *, int *, int *);


/* **** gutil_timing.c **** */

void start_timing(timing_data *);
void stop_timing(timing_data *);
void diff_timing(timing_data *);


/* **** gutil_verb.c **** */

void errmsg(const char *msg, int flag);
void errmsg_v(const char *fmt, ...);
void errmsg_exit(const char *msg);
void errmsg_exit_v(const char *fmt, ...);

void init_progress(progress_data *, int, int, int);
void show_progress(progress_data *, int);
void init_progress_ll(progress_ll_data *, int, int, int);
void show_progress_ll(progress_ll_data *, int, int *);


#ifdef __cplusplus
}
#endif

#endif /* GUTIL_H */
