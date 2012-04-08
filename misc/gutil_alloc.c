/******************************************************************************%
**
**    Copyright (C) 1998-2012 Greg McGarragh <gregm@atmos.colostate.edu>
**
**    This source code is licensed under the GNU General Public License (GPL),
**    Version 3.  See the file COPYING for more details.
**
*******************************************************************************/

#include "gutil.h"


#define MAX_DIMENS 6


/*******************************************************************************
 *
 ******************************************************************************/
size_t array_mem_length(int rank, long *dimen, size_t size) {

     int i;

     size_t length;
     size_t stride;

     length = 0;
     stride = 1;
     for (i = 0; i < rank - 1; ++i) {
          stride *= dimen[i];
          length += stride;
     }

     return length * sizeof(void *) + stride * dimen[i] * size;
}



/*******************************************************************************
 *
 ******************************************************************************/
static void **ptr_array_init(void **a, int rank, long *dimen, int length) {

     if (rank > 0)
          a[0] = ptr_array_init(a + length, rank - 1, dimen + 1, length * dimen[0]);

     return a;
}



/*******************************************************************************
 *
 ******************************************************************************/
static void *ptr_array_alloc(int rank, long *dimen, int size, int depth) {

     long length;

     void *a;

     length = array_mem_length(rank, dimen, size);

     a = (void *) malloc(length);
     if (a == NULL) {
          fprintf(stderr, "ERROR: memory allocation failed: size = %ld bytes\n", length);
          return NULL;
     }

     return ptr_array_init((void **) a, rank - 1, dimen + 1, dimen[0]);
}



/*******************************************************************************
 *
 ******************************************************************************/
static void *ptr_array_realloc(void **a1, int rank, long *dimen, int size, int depth) {

     long length;

     void *a2;

     length = array_mem_length(rank, dimen, size);

     a2 = (void *) realloc(a1, length);
     if (a2 == NULL) {
          fprintf(stderr, "ERROR: memory allocation failed: size = %ld bytes\n", length);
          return NULL;
     }

     return ptr_array_init((void **) a2, rank - 1, dimen + 1, dimen[0]);
}



/*******************************************************************************
 *
 ******************************************************************************/
static void *from_mem_ptr_array(void **a, void *a2, int rank, long *dimen, int length) {

     if (rank > 1)
          a[0] = from_mem_ptr_array(a + length, a2, rank - 1, dimen + 1, length * dimen[0]);
     else
     if (rank > 0)
          a[0] = a2;

     return (void *) a;
}



static void *ptr_array_from_mem(void *a2, int rank, long *dimen, int depth, int flag) {

     long length;

     void *a;

     length = array_mem_length(rank, dimen, 0);

     if (! flag) {
          a = a2;

          a2 = (uchar *) a2 + length;
     }
     else {
          a = (void *) malloc(length);
          if (a == NULL) {
               fprintf(stderr, "ERROR: memory allocation failed: size = %ld bytes\n", length);
               return NULL;
          }
     }

     return from_mem_ptr_array((void **) a, a2, rank - 1, dimen + 1, dimen[0]);
}



/*******************************************************************************
 *
 ******************************************************************************/
static void ptr_array_assign(void **a0, void **a, int rank,
                             long *dimen, int size, long *index, int depth) {

     int i;
     int j;
     int n;

     int size2 = sizeof(void *);

     int rank_m_2;
     int depth_p_1;

     long offset;
     long offset2;
     long offset3;

     rank_m_2  = rank  - 2;
     depth_p_1 = depth + 1;

     if (depth == rank_m_2)
          size2 = size;

     offset  = 0;

     offset2 = dimen[depth_p_1];

     offset3 = offset2 * dimen[depth];
     for (j = depth - 1; j >= 0; --j) {
          offset  += index[j] * offset3;
          offset3 *= dimen[j];
     }

     offset  *= size2;

     offset2 *= size2;

     n = dimen[depth];

     if (depth >= rank_m_2) {
          for (i = 0; i < n; ++i) {
               a[i] = (uchar *) a0[0] + offset;

               offset += offset2;
          }
     }
     else {
          for (i = 0; i < n; ++i) {
               a[i] = (uchar *) a0[0] + offset;

               index[depth] = i;
               ptr_array_assign((void **) a0[0], (void **) a[i], rank, dimen, size, index, depth_p_1);

               offset += offset2;
          }
     }
}



/*******************************************************************************
 *
 ******************************************************************************/
static void ptr_array_free(void **a, int rank, int depth) {

     free(a);
}



/*******************************************************************************
 *
 ******************************************************************************/
void *alloc_array(int rank, long *dimen, int size) {

     long index[MAX_DIMENS];

     void *a;

     a = ptr_array_alloc(rank, dimen, size, 0);

     if (a && rank > 1)
          ptr_array_assign((void **) a, (void **) a, rank, dimen, size, index, 0);

     return a;
}



void *realloc_array(void *a1, int rank, long *dimen, int size) {

     long index[MAX_DIMENS];

     void *a2;

     a2 = ptr_array_realloc((void **) a1, rank, dimen, size, 0);

     if (a2 && rank > 1)
          ptr_array_assign((void **) a2, (void **) a2, rank, dimen, size, index, 0);

     return a2;
}



void *array_from_mem(void *a2, int rank, long *dimen, int size, int flag) {

     long index[MAX_DIMENS];

     void *a;

     a = ptr_array_from_mem((void **) a2, rank, dimen, 0, flag);

     if (a && rank > 1)
          ptr_array_assign((void **) a, (void **) a, rank, dimen, size, index, 0);

     return a;
}



void free_array(void *a, int rank) {

     if (a)
          ptr_array_free((void **) a, rank, 0);
}



/*******************************************************************************
 *
 ******************************************************************************/
size_t array_mem_length1(long m, size_t size) {

     long dimen[1];

     dimen[0] = m;

     return array_mem_length(1, dimen, size);
}



void *alloc_array1(long m, int size) {

     long dimen[1];

     dimen[0] = m;

     return alloc_array(1, dimen, size);
}



void *realloc_array1(void *a, long m, int size) {

     long dimen[1];

     dimen[0] = m;

     return realloc_array(a, 1, dimen, size);
}



void *array_from_mem1(void *a, long m, int size, int flag) {

     long dimen[1];

     dimen[0] = m;

     return array_from_mem(a, 1, dimen, size, flag);
}



void free_array1(void *a) {

     free_array(a, 1);
}



/*******************************************************************************
 *
 ******************************************************************************/
size_t array_mem_length2(long m, long n, size_t size) {

     long dimen[2];

     dimen[0] = m;
     dimen[1] = n;

     return array_mem_length(2, dimen, size);
}



void **alloc_array2(long m, long n, int size) {

     long dimen[2];

     dimen[0] = m;
     dimen[1] = n;

     return (void **) alloc_array(2, dimen, size);
}



void **realloc_array2(void **a, long m, long n, int size) {

     long dimen[2];

     dimen[0] = m;
     dimen[1] = n;

     return (void **) realloc_array(a, 2, dimen, size);
}



void **array_from_mem2(void *a, long m, long n, int size, int flag) {

     long dimen[2];

     dimen[0] = m;
     dimen[1] = n;

     return (void **) array_from_mem(a, 2, dimen, size, flag);
}



void free_array2(void **a) {

     free_array(a, 2);
}



/*******************************************************************************
 *
 ******************************************************************************/
size_t array_mem_length3(long m, long n, long o, size_t size) {

     long dimen[3];

     dimen[0] = m;
     dimen[1] = n;
     dimen[2] = o;

     return array_mem_length(3, dimen, size);
}



void ***alloc_array3(long m, long n, long o, int size) {

     long dimen[3];

     dimen[0] = m;
     dimen[1] = n;
     dimen[2] = o;

     return (void ***) alloc_array(3, dimen, size);
}



void ***realloc_array3(void ***a, long m, long n, long o, int size) {

     long dimen[3];

     dimen[0] = m;
     dimen[1] = n;
     dimen[2] = o;

     return (void ***) realloc_array(a, 3, dimen, size);
}



void ***array_from_mem3(void *a, long m, long n, long o, int size, int flag) {

     long dimen[3];

     dimen[0] = m;
     dimen[1] = n;
     dimen[2] = o;

     return (void ***) array_from_mem(a, 3, dimen, size, flag);
}



void free_array3(void ***a) {

     free_array(a, 3);
}



/*******************************************************************************
 *
 ******************************************************************************/
size_t array_mem_length4(long m, long n, long o, long p, size_t size) {

     long dimen[4];

     dimen[0] = m;
     dimen[1] = n;
     dimen[2] = o;
     dimen[3] = p;

     return array_mem_length(4, dimen, size);
}



void ****alloc_array4(long m, long n, long o, long p, int size) {

     long dimen[4];

     dimen[0] = m;
     dimen[1] = n;
     dimen[2] = o;
     dimen[3] = p;

     return (void ****) alloc_array(4, dimen, size);
}



void ****realloc_array4(void ****a, long m, long n, long o, long p, int size) {

     long dimen[4];

     dimen[0] = m;
     dimen[1] = n;
     dimen[2] = o;
     dimen[3] = p;

     return (void ****) realloc_array(a, 4, dimen, size);
}



void ****array_from_mem4(void *a, long m, long n, long o, long p, int size, int flag) {

     long dimen[4];

     dimen[0] = m;
     dimen[1] = n;
     dimen[2] = o;
     dimen[3] = p;

     return (void ****) array_from_mem(a, 4, dimen, size, flag);
}



void free_array4(void ****a) {

     free_array(a, 4);
}



/*******************************************************************************
 *
 ******************************************************************************/
size_t array_mem_length5(long m, long n, long o, long p, long q, size_t size) {

     long dimen[5];

     dimen[0] = m;
     dimen[1] = n;
     dimen[2] = o;
     dimen[3] = p;
     dimen[4] = q;

     return array_mem_length(5, dimen, size);
}



void *****alloc_array5(long m, long n, long o, long p, long q, int size) {

     long dimen[5];

     dimen[0] = m;
     dimen[1] = n;
     dimen[2] = o;
     dimen[3] = p;
     dimen[4] = q;

     return (void *****) alloc_array(5, dimen, size);
}



void *****realloc_array5(void *****a, long m, long n, long o, long p, long q, int size) {

     long dimen[5];

     dimen[0] = m;
     dimen[1] = n;
     dimen[2] = o;
     dimen[3] = p;
     dimen[3] = q;

     return (void *****) realloc_array(a, 5, dimen, size);
}



void *****array_from_mem5(void *a, long m, long n, long o, long p, long q, int size, int flag) {

     long dimen[5];

     dimen[0] = m;
     dimen[1] = n;
     dimen[2] = o;
     dimen[3] = p;
     dimen[4] = q;

     return (void *****) array_from_mem(a, 5, dimen, size, flag);
}



void free_array5(void *****a) {

     free_array(a, 5);
}



/*******************************************************************************
 *
 ******************************************************************************/
size_t array_mem_length6(long m, long n, long o, long p, long q, long r, size_t size) {

     long dimen[6];

     dimen[0] = m;
     dimen[1] = n;
     dimen[2] = o;
     dimen[3] = p;
     dimen[4] = q;
     dimen[5] = r;

     return array_mem_length(6, dimen, size);
}



void ******alloc_array6(long m, long n, long o, long p, long q, long r, int size) {

     long dimen[6];

     dimen[0] = m;
     dimen[1] = n;
     dimen[2] = o;
     dimen[3] = p;
     dimen[4] = q;
     dimen[5] = r;

     return (void ******) alloc_array(6, dimen, size);
}



void ******realloc_array6(void ******a, long m, long n, long o, long p, long q, long r, int size) {

     long dimen[6];

     dimen[0] = m;
     dimen[1] = n;
     dimen[2] = o;
     dimen[3] = p;
     dimen[3] = q;
     dimen[5] = r;

     return (void ******) realloc_array(a, 6, dimen, size);
}



void ******array_from_mem6(void *a, long m, long n, long o, long p, long q, long r, int size, int flag) {

     long dimen[6];

     dimen[0] = m;
     dimen[1] = n;
     dimen[2] = o;
     dimen[3] = p;
     dimen[4] = q;
     dimen[5] = r;

     return (void ******) array_from_mem(a, 6, dimen, size, flag);
}



void free_array6(void ******a) {

     free_array(a, 6);
}



/*******************************************************************************
 *
 ******************************************************************************/
#define postfix_ uc
#define type_ uchar
#include "galloc.h"
#undef postfix_ 
#undef type_


#define postfix_ c
#define type_ char
#include "galloc.h"
#undef postfix_ 
#undef type_


#define postfix_ us
#define type_ ushort
#include "galloc.h"
#undef postfix_ 
#undef type_


#define postfix_ s
#define type_ short
#include "galloc.h"
#undef postfix_ 
#undef type_


#define postfix_ ui
#define type_ uint
#include "galloc.h"
#undef postfix_ 
#undef type_


#define postfix_ i
#define type_ int
#include "galloc.h"
#undef postfix_ 
#undef type_


#define postfix_ ul
#define type_ ulong
#include "galloc.h"
#undef postfix_ 
#undef type_


#define postfix_ l
#define type_ long
#include "galloc.h"
#undef postfix_ 
#undef type_

/*
#define postfix_ ull
#define type_ uint64
#include "galloc.h"
#undef postfix_ 
#undef type_


#define postfix_ ll
#define type_ int64
#include "galloc.h"
#undef postfix_ 
#undef type_
*/

#define postfix_ f
#define type_ float
#include "galloc.h"
#undef postfix_ 
#undef type_


#define postfix_ fc
#define type_ fcomplex
#include "galloc.h"
#undef postfix_ 
#undef type_


#define postfix_ d
#define type_ double
#include "galloc.h"
#undef postfix_ 
#undef type_


#define postfix_ dc
#define type_ dcomplex
#include "galloc.h"
#undef postfix_ 
#undef type_



/*******************************************************************************
 *
 ******************************************************************************/
void nr_to_c_darray2(double ***a, int n) {

     int i;

     ++(*a);

     for (i = 0; i < n; ++i)
          ++((*a)[i]);
}



void c_to_nr_darray2(double ***a, int n) {

     int i;

     for (i = 0; i < n; ++i)
          --((*a)[i]);

     --(*a);
}
