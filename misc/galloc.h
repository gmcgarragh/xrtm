/*******************************************************************************
**
**    Copyright (C) 1998-2020 Greg McGarragh <greg.mcgarragh@colostate.edu>
**
**    This source code is licensed under the GNU General Public License (GPL),
**    Version 3.  See the file COPYING for more details.
**
*******************************************************************************/

#include "gutil.h"



/*******************************************************************************
 *
 ******************************************************************************/
type_ *XCAT(alloc_array1_, postfix_)(size_t m) {

     size_t dimen[1];

     dimen[0] = m;

     return (type_ *) alloc_array(1, dimen, sizeof(type_));
}



type_ *XCAT(realloc_array1_, postfix_)(type_ *a, size_t m) {

     size_t dimen[1];

     dimen[0] = m;

     return (type_ *) realloc_array(a, 1, dimen, sizeof(type_));
}



type_ *XCAT(array_from_mem1_, postfix_)(type_ *a, size_t m) {

     size_t dimen[1];

     dimen[0] = m;

     return (type_ *) array_from_mem(a, 1, dimen, sizeof(type_), 1);
}



void XCAT(free_array1_, postfix_)(type_ *a) {

     free_array(a, 1);
}



void XCAT(init_array1_, postfix_)(type_ *a, size_t m, type_ x) {

     size_t i;

     for (i = 0; i < m; ++i)
          a[i] = x;
}



void XCAT(copy_array1_, postfix_)(type_ *a2, const type_ *a1, size_t m) {

     size_t i;

     for (i = 0; i < m; ++i)
          a2[i] = a1[i];
}



type_ *XCAT(dup_array1_, postfix_)(type_ *a1, size_t m) {

     type_ *a2;

     a2 = XCAT(alloc_array1_, postfix_)(m);
     XCAT(copy_array1_, postfix_)(a2, a1, m);

     return a2;
}



/*******************************************************************************
 *
 ******************************************************************************/
type_ **XCAT(alloc_array2_, postfix_)(size_t m, size_t n) {

     size_t dimen[2];

     dimen[0] = m;
     dimen[1] = n;

     return (type_ **) alloc_array(2, dimen, sizeof(type_));
}



type_ **XCAT(realloc_array2_, postfix_)(type_ **a, size_t m, size_t n) {

     size_t dimen[2];

     dimen[0] = m;
     dimen[1] = n;

     return (type_ **) realloc_array(a, 2, dimen, sizeof(type_));
}



type_ **XCAT(array_from_mem2_, postfix_)(type_ *a, size_t m, size_t n) {

     size_t dimen[2];

     dimen[0] = m;
     dimen[1] = n;

     return (type_ **) array_from_mem(a, 2, dimen, sizeof(type_), 1);
}



void XCAT(free_array2_, postfix_)(type_ **a) {

     free_array(a, 2);
}



void XCAT(init_array2_, postfix_)(type_ **a, size_t m, size_t n, type_ x) {
/*
     size_t i;

     size_t size;

     type_ *a2 = *a;

     size = m * n;

     for (i = 0; i < size; ++i)
          a2[i] = x;
*/
     size_t i;

     for (i = 0; i < m; ++i)
          XCAT(init_array1_, postfix_)(a[i], n, x);
}



void XCAT(copy_array2_, postfix_)(type_ **a2, type_ **a1, size_t m, size_t n) {
/*
     size_t i;

     size_t size;

     type_ *a12 = *a1;
     type_ *a22 = *a2;

     size = m * n;

     for (i = 0; i < size; ++i)
          a22[i] = a12[i];
*/
     size_t i;

     for (i = 0; i < m; ++i)
          XCAT(copy_array1_, postfix_)(a2[i], a1[i], n);
}



type_ **XCAT(dup_array2_, postfix_)(type_ **a1, size_t m, size_t n) {

     type_ **a2;

     a2 = XCAT(alloc_array2_, postfix_)(m, n);
     XCAT(copy_array2_, postfix_)(a2, a1, m, n);

     return a2;
}



/*******************************************************************************
 *
 ******************************************************************************/
type_ ***XCAT(alloc_array3_, postfix_)(size_t m, size_t n, size_t o) {

     size_t dimen[3];

     dimen[0] = m;
     dimen[1] = n;
     dimen[2] = o;

     return (type_ ***) alloc_array(3, dimen, sizeof(type_));
}



type_ ***XCAT(realloc_array3_, postfix_)(type_ ***a, size_t m, size_t n, size_t o) {

     size_t dimen[3];

     dimen[0] = m;
     dimen[1] = n;
     dimen[2] = o;

     return (type_ ***) realloc_array(a, 3, dimen, sizeof(type_));
}



type_ ***XCAT(array_from_mem3_, postfix_)(type_ *a, size_t m, size_t n, size_t o) {

     size_t dimen[3];

     dimen[0] = m;
     dimen[1] = n;
     dimen[2] = o;

     return (type_ ***) array_from_mem(a, 3, dimen, sizeof(type_), 1);
}



void XCAT(free_array3_, postfix_)(type_ ***a) {

     free_array(a, 3);
}



void XCAT(init_array3_, postfix_)(type_ ***a, size_t m, size_t n, size_t o, type_ x) {
/*
     size_t i;

     size_t size;

     type_ *a2 = **a;

     size = m * n * o;

     for (i = 0; i < size; ++i)
          a2[i] = x;
*/
     size_t i;

     for (i = 0; i < m; ++i)
          XCAT(init_array2_, postfix_)(a[i], n, o, x);
}



void XCAT(copy_array3_, postfix_)(type_ ***a2, type_ ***a1, size_t m, size_t n, size_t o) {
/*
     size_t i;

     size_t size;

     type_ *a12 = **a1;
     type_ *a22 = **a2;

     size = m * n * o;

     for (i = 0; i < size; ++i)
          a22[i] = a12[i];
*/
     size_t i;

     for (i = 0; i < m; ++i)
          XCAT(copy_array2_, postfix_)(a2[i], a1[i], n, o);
}



type_ ***XCAT(dup_array3_, postfix_)(type_ ***a1, size_t m, size_t n, size_t o) {

     type_ ***a2;

     a2 = XCAT(alloc_array3_, postfix_)(m, n, o);
     XCAT(copy_array3_, postfix_)(a2, a1, m, n, o);

     return a2;
}



/*******************************************************************************
 *
 ******************************************************************************/
type_ ****XCAT(alloc_array4_, postfix_)(size_t m, size_t n, size_t o, size_t p) {

     size_t dimen[4];

     dimen[0] = m;
     dimen[1] = n;
     dimen[2] = o;
     dimen[3] = p;

     return (type_ ****) alloc_array(4, dimen, sizeof(type_));
}



type_ ****XCAT(realloc_array4_, postfix_)(type_ ****a, size_t m, size_t n, size_t o, size_t p) {

     size_t dimen[4];

     dimen[0] = m;
     dimen[1] = n;
     dimen[2] = o;
     dimen[3] = p;

     return (type_ ****) realloc_array(a, 4, dimen, sizeof(type_));
}



type_ ****XCAT(array_from_mem4_, postfix_)(type_ *a, size_t m, size_t n, size_t o, size_t p) {

     size_t dimen[4];

     dimen[0] = m;
     dimen[1] = n;
     dimen[2] = o;
     dimen[3] = p;

     return (type_ ****) array_from_mem(a, 4, dimen, sizeof(type_), 1);
}



void XCAT(free_array4_, postfix_)(type_ ****a) {

     free_array(a, 4);
}



void XCAT(init_array4_, postfix_)(type_ ****a, size_t m, size_t n, size_t o, size_t p, type_ x) {
/*
     size_t i;

     size_t size;

     type_ *a2 = ***a;

     size = m * n * o * p;

     for (i = 0; i < size; ++i)
          a2[i] = x;
*/
     size_t i;

     for (i = 0; i < m; ++i)
          XCAT(init_array3_, postfix_)(a[i], n, o, p, x);
}



void XCAT(copy_array4_, postfix_)(type_ ****a2, type_ ****a1, size_t m, size_t n, size_t o, size_t p) {
/*
     size_t i;

     size_t size;

     type_ *a12 = ***a1;
     type_ *a22 = ***a2;

     size = m * n * o * p;

     for (i = 0; i < size; ++i)
          a22[i] = a12[i];
*/
     size_t i;

     for (i = 0; i < m; ++i)
          XCAT(copy_array3_, postfix_)(a2[i], a1[i], n, o, p);
}



/*******************************************************************************
 *
 ******************************************************************************/
type_ *****XCAT(alloc_array5_, postfix_)(size_t m, size_t n, size_t o, size_t p, size_t q) {

     size_t dimen[5];

     dimen[0] = m;
     dimen[1] = n;
     dimen[2] = o;
     dimen[3] = p;
     dimen[4] = q;

     return (type_ *****) alloc_array(5, dimen, sizeof(type_));
}



type_ *****XCAT(realloc_array5_, postfix_)(type_ *****a, size_t m, size_t n, size_t o, size_t p, size_t q) {

     size_t dimen[5];

     dimen[0] = m;
     dimen[1] = n;
     dimen[2] = o;
     dimen[3] = p;
     dimen[4] = q;

     return (type_ *****) realloc_array(a, 5, dimen, sizeof(type_));
}



type_ *****XCAT(array_from_mem5_, postfix_)(type_ *a, size_t m, size_t n, size_t o, size_t p, size_t q) {

     size_t dimen[5];

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



void XCAT(init_array5_, postfix_)(type_ *****a, size_t m, size_t n, size_t o, size_t p, size_t q, type_ x) {
/*
     size_t i;

     size_t size;

     type_ *a2 = ****a;

     size = m * n * o * p * q;

     for (i = 0; i < size; ++i)
          a2[i] = x;
*/
     size_t i;

     for (i = 0; i < m; ++i)
          XCAT(init_array4_, postfix_)(a[i], n, o, p, q, x);
}



void XCAT(copy_array5_, postfix_)(type_ *****a2, type_ *****a1, size_t m, size_t n, size_t o, size_t p, size_t q) {
/*
     size_t i;

     size_t size;

     type_ *a12 = ****a1;
     type_ *a22 = ****a2;

     size = m * n * o * p * q;

     for (i = 0; i < size; ++i)
          a22[i] = a12[i];
*/
     size_t i;

     for (i = 0; i < m; ++i)
          XCAT(copy_array4_, postfix_)(a2[i], a1[i], n, o, p, q);
}



/*******************************************************************************
 *
 ******************************************************************************/
type_ ******XCAT(alloc_array6_, postfix_)(size_t m, size_t n, size_t o, size_t p, size_t q, size_t r) {

     size_t dimen[6];

     dimen[0] = m;
     dimen[1] = n;
     dimen[2] = o;
     dimen[3] = p;
     dimen[4] = q;
     dimen[5] = r;

     return (type_ ******) alloc_array(6, dimen, sizeof(type_));
}



type_ ******XCAT(realloc_array6_, postfix_)(type_ ******a, size_t m, size_t n, size_t o, size_t p, size_t q, size_t r) {

     size_t dimen[6];

     dimen[0] = m;
     dimen[1] = n;
     dimen[2] = o;
     dimen[3] = p;
     dimen[4] = q;
     dimen[5] = r;

     return (type_ ******) realloc_array(a, 6, dimen, sizeof(type_));
}



type_ ******XCAT(array_from_mem6_, postfix_)(type_ *a, size_t m, size_t n, size_t o, size_t p, size_t q, size_t r) {

     size_t dimen[6];

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



void XCAT(init_array6_, postfix_)(type_ ******a, size_t m, size_t n, size_t o, size_t p, size_t q, size_t r, type_ x) {
/*
     size_t i;

     size_t size;

     type_ *a2 = *****a;

     size = m * n * o * p * q * r;

     for (i = 0; i < size; ++i)
          a2[i] = x;
*/
     size_t i;

     for (i = 0; i < m; ++i)
          XCAT(init_array5_, postfix_)(a[i], n, o, p, q, r, x);
}



void XCAT(copy_array6_, postfix_)(type_ ******a2, type_ ******a1, size_t m, size_t n, size_t o, size_t p, size_t q, size_t r) {
/*
     size_t i;

     size_t size;

     type_ *a12 = *****a1;
     type_ *a22 = *****a2;

     size = m * n * o * p * q * r;

     for (i = 0; i < size; ++i)
          a22[i] = a12[i];
*/
     size_t i;

     for (i = 0; i < m; ++i)
          XCAT(copy_array5_, postfix_)(a2[i], a1[i], n, o, p, q, r);
}



/*******************************************************************************
 *
 ******************************************************************************/
type_ *******XCAT(alloc_array7_, postfix_)(size_t m, size_t n, size_t o, size_t p, size_t q, size_t r, size_t s) {

     size_t dimen[7];

     dimen[0] = m;
     dimen[1] = n;
     dimen[2] = o;
     dimen[3] = p;
     dimen[4] = q;
     dimen[5] = r;
     dimen[6] = s;

     return (type_ *******) alloc_array(7, dimen, sizeof(type_));
}



type_ *******XCAT(realloc_array7_, postfix_)(type_ *******a, size_t m, size_t n, size_t o, size_t p, size_t q, size_t r, size_t s) {

     size_t dimen[7];

     dimen[0] = m;
     dimen[1] = n;
     dimen[2] = o;
     dimen[3] = p;
     dimen[4] = q;
     dimen[5] = r;
     dimen[6] = s;

     return (type_ *******) realloc_array(a, 7, dimen, sizeof(type_));
}



type_ *******XCAT(array_from_mem7_, postfix_)(type_ *a, size_t m, size_t n, size_t o, size_t p, size_t q, size_t r, size_t s) {

     size_t dimen[7];

     dimen[0] = m;
     dimen[1] = n;
     dimen[2] = o;
     dimen[3] = p;
     dimen[4] = q;
     dimen[5] = r;
     dimen[6] = s;

     return (type_ *******) array_from_mem(a, 7, dimen, sizeof(type_), 1);
}



void XCAT(free_array7_, postfix_)(type_ *******a) {

     free_array(a, 7);
}



void XCAT(init_array7_, postfix_)(type_ *******a, size_t m, size_t n, size_t o, size_t p, size_t q, size_t r, size_t s, type_ x) {
/*
     size_t i;

     size_t size;

     type_ *a2 = ******a;

     size = m * n * o * p * q * r * s;

     for (i = 0; i < size; ++i)
          a2[i] = x;
*/
     size_t i;

     for (i = 0; i < m; ++i)
          XCAT(init_array6_, postfix_)(a[i], n, o, p, q, r, s, x);
}



void XCAT(copy_array7_, postfix_)(type_ *******a2, type_ *******a1, size_t m, size_t n, size_t o, size_t p, size_t q, size_t r, size_t s) {
/*
     size_t i;

     size_t size;

     type_ *a12 = ******a1;
     type_ *a22 = ******a2;

     size = m * n * o * p * q * r * s;

     for (i = 0; i < size; ++i)
          a22[i] = a12[i];
*/
     size_t i;

     for (i = 0; i < m; ++i)
          XCAT(copy_array6_, postfix_)(a2[i], a1[i], n, o, p, q, r, s);
}
