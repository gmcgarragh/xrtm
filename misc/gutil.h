/******************************************************************************%
**
**    Copyright (C) 1998-2012 Greg McGarragh <gregm@atmos.colostate.edu>
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

#ifdef __cplusplus
extern "C" {
#endif


#define WIN32_MSVC	16


#if PLATFORM == WIN32_MSVC
#define finite		_finite
#define isinf		_isinf
#define isnan		_isnan
#define snprintf	_snprintf
#define strdup		_strdup

#define strtok_r(a, b, c) (strtok(a, b))
#endif


typedef unsigned char	uchar;
typedef unsigned short	ushort;
typedef unsigned int	uint;
typedef unsigned long	ulong;


#ifndef __cplusplus
typedef float _Complex 		fcomplex;
typedef double _Complex 	dcomplex;
#else
#define _Complex_I		std::complex<double>(0., 1.)
#define cabs			abs
#define cexp			exp
#define cpow			pow
#define creal			real
#define cimag			imag
#define csqrt			sqrt
typedef std::complex<float> 	fcomplex;
typedef std::complex<double> 	dcomplex;
#endif


#ifdef CAT
#undef CAT
#endif
#define CAT(a, b) a##b

#ifdef XCAT
#undef XCAT
#endif
#define XCAT(x, y) CAT(x, y)


#ifdef MIN
#undef MIN
#endif
#define MIN(a,b) ((a) < (b) ? (a) : (b)) 

#ifdef MAX
#undef MAX
#endif
#define MAX(a,b) ((a) > (b) ? (a) : (b)) 


#define LENGTH(NAME) (sizeof(NAME) / sizeof(NAME[0]))


#ifdef M_PI
#define PI	M_PI
#else
#define PI	3.14159265358979323846264338327
#endif


#define D2R	(PI/180.0)
#define R2D	(180.0/PI)


/* Structure types */

typedef struct {
     int offset;
     int size;
     int n;
     int size_lst;
     int perc_lst;
} progress;


typedef struct {
     int offset;
     int size;
     int n;
     int size_lst;
     int perc_lst;
     int current;
} progress_ll;



/* Function prototypes */

/* **** gutil_alloc.c **** */

size_t array_mem_length(int rank, long *dimen, size_t size);
void *alloc_array(int, long *, int);
void *realloc_array(void *, int, long *, int);
void *array_from_mem(void *, int, long *, int, int);
void free_array(void *, int);

size_t array_mem_length1(long m, size_t size);
void *alloc_array1(long, int);
void *realloc_array1(void *, long, int);
void *array_from_mem1(void *, long, int, int);
void free_array1(void *);
size_t array_mem_length2(long m, long n, size_t size);
void **alloc_array2(long, long, int);
void **realloc_array2(void **, long, long, int);
void **array_from_mem2(void *, long, long, int, int);
void free_array2(void **);
size_t array_mem_length3(long m, long n, long o, size_t size);
void ***alloc_array3(long, long, long, int);
void ***realloc_array3(void ***, long, long, long, int);
void ***array_from_mem3(void *, long, long, long, int, int);
void free_array3(void ***);
size_t array_mem_length4(long m, long n, long o, long p, size_t size);
void ****alloc_array4(long, long, long, long, int);
void ****realloc_array4(void ****, long, long, long, long, int);
void ****array_from_mem4(void *, long, long, long, long, int, int);
void free_array4(void ****);
size_t array_mem_length5(long m, long n, long o, long p, long q, size_t size);
void *****alloc_array5(long, long, long, long, long, int);
void *****realloc_array5(void *****, long, long, long, long, long, int);
void *****array_from_mem5(void *, long, long, long, long, long, int, int);
void free_array5(void *****);
size_t array_mem_length6(long m, long n, long o, long p, long q, long r, size_t size);
void ******alloc_array6(long, long, long, long, long, long, int);
void ******realloc_array6(void ******, long, long, long, long, long, long, int);
void ******array_from_mem6(void *, long, long, long, long, long, long, int, int);
void free_array6(void ******);

#define all_alloc_proto(postfix_, type_)							\
type_ *XCAT(alloc_array1_, postfix_)(long);							\
type_ *XCAT(realloc_array1_, postfix_)(type_ *, long);						\
type_ *XCAT(array_from_mem1_, postfix_)(type_ *, long);						\
void XCAT(free_array1_, postfix_)(type_ *);							\
void XCAT(init_array1_, postfix_)(type_ *, long, type_);					\
void XCAT(copy_array1_, postfix_)(type_ *, type_ *, long);					\
type_ **XCAT(alloc_array2_, postfix_)(long, long);						\
type_ **XCAT(realloc_array2_, postfix_)(type_ **, long, long);					\
type_ **XCAT(array_from_mem2_, postfix_)(type_ *, long, long);					\
void XCAT(free_array2_, postfix_)(type_ **);							\
void XCAT(init_array2_, postfix_)(type_ **, long, long, type_);					\
void XCAT(copy_array2_, postfix_)(type_ **, type_ **, long, long);				\
type_ ***XCAT(alloc_array3_, postfix_)(long, long, long);					\
type_ ***XCAT(realloc_array3_, postfix_)(type_ ***, long, long, long);				\
type_ ***XCAT(array_from_mem3_, postfix_)(type_ *, long, long, long);				\
void XCAT(free_array3_, postfix_)(type_ ***);							\
void XCAT(init_array3_, postfix_)(type_ ***, long, long, long, type_);				\
void XCAT(copy_array3_, postfix_)(type_ ***, type_ ***, long, long, long);			\
type_ ****XCAT(alloc_array4_, postfix_)(long, long, long, long);				\
type_ ****XCAT(realloc_array4_, postfix_)(type_ ****, long, long, long, long);			\
type_ ****XCAT(array_from_mem4_, postfix_)(type_ *, long, long, long, long);			\
void XCAT(free_array4_, postfix_)(type_ ****);							\
void XCAT(init_array4_, postfix_)(type_ ****, long, long, long, long, type_);			\
void XCAT(copy_array4_, postfix_)(type_ ****, type_ ****, long, long, long, long);		\
type_ *****XCAT(alloc_array5_, postfix_)(long, long, long, long, long);				\
type_ *****XCAT(realloc_array5_, postfix_)(type_ *****, long, long, long, long, long);		\
type_ *****XCAT(array_from_mem5_, postfix_)(type_ *, long, long, long, long, long);		\
void XCAT(free_array5_, postfix_)(type_ *****);							\
void XCAT(init_array5_, postfix_)(type_ *****, long, long, long, long, long, type_);		\
void XCAT(copy_array5_, postfix_)(type_ *****, type_ *****, long, long, long, long, long);	\
type_ ******XCAT(alloc_array6_, postfix_)(long, long, long, long, long, long);			\
type_ ******XCAT(realloc_array6_, postfix_)(type_ ******, long, long, long, long, long, long);	\
type_ ******XCAT(array_from_mem6_, postfix_)(type_ *, long, long, long, long, long, long);	\
void XCAT(free_array6_, postfix_)(type_ ******);						\
void XCAT(init_array6_, postfix_)(type_ ******, long, long, long, long, long, long, type_);	\
void XCAT(copy_array6_, postfix_)(type_ ******, type_ ******, long, long, long, long, long, long);

all_alloc_proto(uc,  uchar)
all_alloc_proto(c,   char)
all_alloc_proto(us,  ushort)
all_alloc_proto(s,   short)
all_alloc_proto(ui,  uint)
all_alloc_proto(i,   int)
all_alloc_proto(ul,  ulong)
all_alloc_proto(l,   long)
/*
all_alloc_proto(ull, uint64)
all_alloc_proto(ll,  int64)
*/
all_alloc_proto(f,   float)
all_alloc_proto(fc,  fcomplex)
all_alloc_proto(d,   double)
all_alloc_proto(dc,  dcomplex)

void nr_to_c_dmatrix(double ***, int);
void c_to_nr_dmatrix(double ***, int);


/* **** gutil_error.c **** */

void eprintf(const char *fmt, ...);

void check_arg_count(int i, int argc, int n, const char *s);

int strtoi_errmsg(const char *string, const char *name, int *value);
int strtoi_errmsg_exit(const char *string, const char *name);
long strtol_errmsg(const char *string, const char *name, long *value);
long strtol_errmsg_exit(const char *string, const char *name);
double strtod_errmsg(const char *string, const char *name, double *value);
double strtod_errmsg_exit(const char *string, const char *name);

int gsystem_errmsg(const char *, int);
int gsystem_errmsg_exit(const char *, int);



/* **** gutil_parse.c **** */

void read_newline(FILE *fp);
int read_comment(FILE *fp, char comment);
void read_comments(FILE *fp, char comment);


#ifdef __cplusplus
}
#endif

#endif /* GUTIL_H */
