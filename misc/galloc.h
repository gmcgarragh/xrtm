/******************************************************************************%
**
**    Copyright (C) 1998-2012 Greg McGarragh <gregm@atmos.colostate.edu>
**
**    This source code is licensed under the GNU General Public License (GPL),
**    Version 3.  See the file COPYING for more details.
**
*******************************************************************************/

#include "gutil.h"



/*******************************************************************************
 *
 ******************************************************************************/
type_ *XCAT(alloc_array1_, postfix_)(long m) {

     long dimen[1];

     dimen[0] = m;

     return (type_ *) alloc_array(1, dimen, sizeof(type_));
}



type_ *XCAT(realloc_array1_, postfix_)(type_ *a, long m) {

     long dimen[1];

     dimen[0] = m;

     return (type_ *) realloc_array(a, 1, dimen, sizeof(type_));
}



type_ *XCAT(array_from_mem1_, postfix_)(type_ *a, long m) {

     long dimen[1];

     dimen[0] = m;

     return (type_ *) array_from_mem(a, 1, dimen, sizeof(type_), 1);
}



void XCAT(free_array1_, postfix_)(type_ *a) {

     free_array(a, 1);
}



void XCAT(init_array1_, postfix_)(type_ *a, long m, type_ x) {

     long i;

     for (i = 0; i < m; ++i)
          a[i] = x;
}



void XCAT(copy_array1_, postfix_)(type_ *a2, type_ *a1, long m) {

     long i;

     for (i = 0; i < m; ++i)
          a2[i] = a1[i];
}



/*******************************************************************************
 *
 ******************************************************************************/
type_ **XCAT(alloc_array2_, postfix_)(long m, long n) {

     long dimen[2];

     dimen[0] = m;
     dimen[1] = n;

     return (type_ **) alloc_array(2, dimen, sizeof(type_));
}



type_ **XCAT(realloc_array2_, postfix_)(type_ **a, long m, long n) {

     long dimen[2];

     dimen[0] = m;
     dimen[1] = n;

     return (type_ **) realloc_array(a, 2, dimen, sizeof(type_));
}



type_ **XCAT(array_from_mem2_, postfix_)(type_ *a, long m, long n) {

     long dimen[2];

     dimen[0] = m;
     dimen[1] = n;

     return (type_ **) array_from_mem(a, 2, dimen, sizeof(type_), 1);
}



void XCAT(free_array2_, postfix_)(type_ **a) {

     free_array(a, 2);
}



void XCAT(init_array2_, postfix_)(type_ **a, long m, long n, type_ x) {

     long i;

     long size;

     type_ *a2 = *a;

     size = m * n;

     for (i = 0; i < size; ++i)
          a2[i] = x;
}



void XCAT(copy_array2_, postfix_)(type_ **a2, type_ **a1, long m, long n) {

     long i;

     long size;

     type_ *a12 = *a1;
     type_ *a22 = *a2;

     size = m * n;

     for (i = 0; i < size; ++i)
          a22[i] = a12[i];
}



/*******************************************************************************
 *
 ******************************************************************************/
type_ ***XCAT(alloc_array3_, postfix_)(long m, long n, long o) {

     long dimen[3];

     dimen[0] = m;
     dimen[1] = n;
     dimen[2] = o;

     return (type_ ***) alloc_array(3, dimen, sizeof(type_));
}



type_ ***XCAT(realloc_array3_, postfix_)(type_ ***a, long m, long n, long o) {

     long dimen[3];

     dimen[0] = m;
     dimen[1] = n;
     dimen[2] = o;

     return (type_ ***) realloc_array(a, 3, dimen, sizeof(type_));
}



type_ ***XCAT(array_from_mem3_, postfix_)(type_ *a, long m, long n, long o) {

     long dimen[3];

     dimen[0] = m;
     dimen[1] = n;
     dimen[2] = o;

     return (type_ ***) array_from_mem(a, 3, dimen, sizeof(type_), 1);
}



void XCAT(free_array3_, postfix_)(type_ ***a) {

     free_array(a, 3);
}



void XCAT(init_array3_, postfix_)(type_ ***a, long m, long n, long o, type_ x) {

     long i;

     long size;

     type_ *a2 = **a;

     size = m * n * o;

     for (i = 0; i < size; ++i)
          a2[i] = x;
}



void XCAT(copy_array3_, postfix_)(type_ ***a2, type_ ***a1, long m, long n, long o) {

     long i;

     long size;

     type_ *a12 = **a1;
     type_ *a22 = **a2;

     size = m * n * o;

     for (i = 0; i < size; ++i)
          a22[i] = a12[i];
}



/*******************************************************************************
 *
 ******************************************************************************/
type_ ****XCAT(alloc_array4_, postfix_)(long m, long n, long o, long p) {

     long dimen[4];

     dimen[0] = m;
     dimen[1] = n;
     dimen[2] = o;
     dimen[3] = p;

     return (type_ ****) alloc_array(4, dimen, sizeof(type_));
}



type_ ****XCAT(realloc_array4_, postfix_)(type_ ****a, long m, long n, long o, long p) {

     long dimen[4];

     dimen[0] = m;
     dimen[1] = n;
     dimen[2] = o;
     dimen[3] = p;

     return (type_ ****) realloc_array(a, 4, dimen, sizeof(type_));
}



type_ ****XCAT(array_from_mem4_, postfix_)(type_ *a, long m, long n, long o, long p) {

     long dimen[4];

     dimen[0] = m;
     dimen[1] = n;
     dimen[2] = o;
     dimen[3] = p;

     return (type_ ****) array_from_mem(a, 4, dimen, sizeof(type_), 1);
}



void XCAT(free_array4_, postfix_)(type_ ****a) {

     free_array(a, 4);
}



void XCAT(init_array4_, postfix_)(type_ ****a, long m, long n, long o, long p, type_ x) {

     long i;

     long size;

     type_ *a2 = ***a;

     size = m * n * o * p;

     for (i = 0; i < size; ++i)
          a2[i] = x;
}



void XCAT(copy_array4_, postfix_)(type_ ****a2, type_ ****a1, long m, long n, long o, long p) {

     long i;

     long size;

     type_ *a12 = ***a1;
     type_ *a22 = ***a2;

     size = m * n * o * p;

     for (i = 0; i < size; ++i)
          a22[i] = a12[i];
}



/*******************************************************************************
 *
 ******************************************************************************/
type_ *****XCAT(alloc_array5_, postfix_)(long m, long n, long o, long p, long q) {

     long dimen[5];

     dimen[0] = m;
     dimen[1] = n;
     dimen[2] = o;
     dimen[3] = p;
     dimen[4] = q;

     return (type_ *****) alloc_array(5, dimen, sizeof(type_));
}



type_ *****XCAT(realloc_array5_, postfix_)(type_ *****a, long m, long n, long o, long p, long q) {

     long dimen[5];

     dimen[0] = m;
     dimen[1] = n;
     dimen[2] = o;
     dimen[3] = p;
     dimen[4] = q;

     return (type_ *****) realloc_array(a, 5, dimen, sizeof(type_));
}



type_ *****XCAT(array_from_mem5_, postfix_)(type_ *a, long m, long n, long o, long p, long q) {

     long dimen[5];

     dimen[0] = m;
     dimen[1] = n;
     dimen[2] = o;
     dimen[3] = p;
     dimen[4] = q;

     return (type_ *****) array_from_mem(a, 5, dimen, sizeof(type_), 1);
}



void XCAT(free_array5_, postfix_)(type_ *****a) {

     free_array(a, 5);
}



void XCAT(init_array5_, postfix_)(type_ *****a, long m, long n, long o, long p, long q, type_ x) {

     long i;

     long size;

     type_ *a2 = ****a;

     size = m * n * o * p * q;

     for (i = 0; i < size; ++i)
          a2[i] = x;
}



void XCAT(copy_array5_, postfix_)(type_ *****a2, type_ *****a1, long m, long n, long o, long p, long q) {

     long i;

     long size;

     type_ *a12 = ****a1;
     type_ *a22 = ****a2;

     size = m * n * o * p * q;

     for (i = 0; i < size; ++i)
          a22[i] = a12[i];
}



/*******************************************************************************
 *
 ******************************************************************************/
type_ ******XCAT(alloc_array6_, postfix_)(long m, long n, long o, long p, long q, long r) {

     long dimen[6];

     dimen[0] = m;
     dimen[1] = n;
     dimen[2] = o;
     dimen[3] = p;
     dimen[4] = q;
     dimen[5] = r;

     return (type_ ******) alloc_array(6, dimen, sizeof(type_));
}



type_ ******XCAT(realloc_array6_, postfix_)(type_ ******a, long m, long n, long o, long p, long q, long r) {

     long dimen[6];

     dimen[0] = m;
     dimen[1] = n;
     dimen[2] = o;
     dimen[3] = p;
     dimen[4] = q;
     dimen[5] = r;

     return (type_ ******) realloc_array(a, 6, dimen, sizeof(type_));
}



type_ ******XCAT(array_from_mem6_, postfix_)(type_ *a, long m, long n, long o, long p, long q, long r) {

     long dimen[6];

     dimen[0] = m;
     dimen[1] = n;
     dimen[2] = o;
     dimen[3] = p;
     dimen[4] = q;
     dimen[5] = r;

     return (type_ ******) array_from_mem(a, 6, dimen, sizeof(type_), 1);
}



void XCAT(free_array6_, postfix_)(type_ ******a) {

     free_array(a, 6);
}



void XCAT(init_array6_, postfix_)(type_ ******a, long m, long n, long o, long p, long q, long r, type_ x) {

     long i;

     long size;

     type_ *a2 = *****a;

     size = m * n * o * p * q * r;

     for (i = 0; i < size; ++i)
          a2[i] = x;
}



void XCAT(copy_array6_, postfix_)(type_ ******a2, type_ ******a1, long m, long n, long o, long p, long q, long r) {

     long i;

     long size;

     type_ *a12 = *****a1;
     type_ *a22 = *****a2;

     size = m * n * o * p * q * r;

     for (i = 0; i < size; ++i)
          a22[i] = a12[i];
}
