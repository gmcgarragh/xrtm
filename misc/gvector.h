/*******************************************************************************
**
**    Copyright (C) 1998-2020 Greg McGarragh <greg.mcgarragh@colostate.edu>
**
**    This source code is licensed under the GNU General Public License (GPL),
**    Version 3.  See the file COPYING for more details.
**
*******************************************************************************/

#include "gutil.h"


void XCAT(prefix_, vec_print)(type_ *a, long n) {

     long i;

     for (i = 0; i < n; ++i)
#if is_complex_ == 0
          printf("%e\n", a[i]);
#else
          printf("(%e, %e)\n", creal(a[i]), cimag(a[i]));
#endif
}



void XCAT(prefix_, vec_zero)(type_ *a, long n) {

     long i;

     for (i = 0; i < n; ++i)
          a[i] = 0.0;
}



void XCAT(prefix_, vec_init)(type_ *a, type_ alpha, long n) {

     long i;

     for (i = 0; i < n; ++i)
          a[i] = alpha;
}



void XCAT(prefix_, vec_copy)(type_ *a2, const type_ *a1, long n) {

     long i;

     for (i = 0; i < n; ++i)
          a2[i] = a1[i];
}



void XCAT(prefix_, vec_scale)(type_ s, const type_ *a, type_ *b, long n) {

     long i;

     if (a == b) {
          if (s == -1.) {
               for (i = 0; i < n; ++i)
                    b[i] = -a[i];
          }
          else {
               for (i = 0; i < n; ++i)
                    b[i] *= s;
          }
     }
     else {
          if (s ==  1.) {
               for (i = 0; i < n; ++i)
                    b[i] =  a[i];
          }
          else
          if (s == -1.) {
               for (i = 0; i < n; ++i)
                    b[i] = -a[i];
          }
          else {
               for (i = 0; i < n; ++i)
                    b[i] = a[i] * s;
          }
     }
}



void XCAT(prefix_, vec_add)(const type_ *a, const type_ *b, type_ *c, long n) {

     long i;

     for (i = 0; i < n; ++i)
          c[i] = a[i] + b[i];
}



void XCAT(prefix_, vec_sub)(const type_ *a, const type_ *b, type_ *c, long n) {

     long i;

     for (i = 0; i < n; ++i)
          c[i] = a[i] - b[i];
}



void XCAT(prefix_, vec_lin_cmb)(const type_ *a1, type_ r1,
                                const type_ *a2, type_ r2, type_ *b, long n) {

     long i;

     for (i = 0; i < n; ++i)
          b[i] = r1 * a1[i] + r2 * a2[i];
}



type_ XCAT(prefix_, vec_dot)(const type_ *a, const type_ *b, long n) {
     long i;

     type_ x = 0.;

     for (i = 0; i < n; ++i)
          x += a[i]*b[i];

     return x;
}



type_ XCAT(prefix_, vec_mag)(const type_ *a, long n) {

     return xsqrt_(XCAT(prefix_, vec_dot)(a, a, n));
}



void XCAT(prefix_, vec_unit)(const type_ *a, type_ *u, long n) {

     long i;

     type_ mag;

     mag = XCAT(prefix_, vec_mag)(a, n);

     for (i = 0; i < n; ++i)
          u[i] = a[i] / mag;
}



void XCAT(prefix_, vec_inv_elem)(const type_ *a, type_ *b, long n) {

     long i;

     for (i = 0; i < n; ++i)
          b[i] = 1. / a[i];
}



void XCAT(prefix_, vec_mul_elem)(const type_ *a, const type_ *b, type_ *c, long n) {

     long i;

     for (i = 0; i < n; ++i)
          c[i] = a[i] * b[i];
}



void XCAT(prefix_, vec_div_elem)(const type_ *a, const type_ *b, type_ *c, long n) {

     long i;

     for (i = 0; i < n; ++i)
          c[i] = a[i] / b[i];
}
