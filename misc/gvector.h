/******************************************************************************%
**
**    Copyright (C) 1998-2012 Greg McGarragh <gregm@atmos.colostate.edu>
**
**    This source code is licensed under the GNU General Public License (GPL),
**    Version 3.  See the file COPYING for more details.
**
*******************************************************************************/

#include "gutil.h"


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



void XCAT(prefix_, vec_copy)(type_ *a2, type_ *a1, long n) {

     long i;

     for (i = 0; i < n; ++i)
          a2[i] = a1[i];
}



void XCAT(prefix_, vec_scale)(type_ s, type_ *a, type_ *b, long n) {

     long i;

     if (a == b) {
          if (s == -1.) {
               for (i = 0; i < n; ++i)
                    b[i] = -a[i];
          }
          else {
               for (i = 0; i < n; ++i)
                    a[i] *= s;
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



void XCAT(prefix_, vec_add)(type_ *a, type_ *b, type_ *c, long n) {

     long i;

     for (i = 0; i < n; ++i)
          c[i] = a[i] + b[i];
}



void XCAT(prefix_, vec_sub)(type_ *a, type_ *b, type_ *c, long n) {

     long i;

     for (i = 0; i < n; ++i)
          c[i] = a[i] - b[i];
}



void XCAT(prefix_, vec_comb)(type_ *a1, type_ r1,
                             type_ *a2, type_ r2, type_ *b, long n) {

     long i;

     for (i = 0; i < n; ++i)
          b[i] = r1 * a1[i] + r2 * a2[i];
}



type_ XCAT(prefix_, vec_dot)(type_ *a, type_ *b, long n) {
     long i;

     type_ x = 0.;

     for (i = 0; i < n; ++i)
          x += a[i]*b[i];

     return x;
}



type_ XCAT(prefix_, vec_mag)(type_ *a, long n) {

     return sqrt(XCAT(prefix_, vec_dot)(a, a, n));
}



void XCAT(prefix_, vec_unit)(type_ *a, type_ *u, long n) {

     long i;

     type_ mag;

     mag = XCAT(prefix_, vec_mag)(a, n);

     for (i = 0; i < n; ++i)
          u[i] = a[i] / mag;
}



void XCAT(prefix_, vec_inv)(type_ *a, type_ *b, long n) {

     long i;

     for (i = 0; i < n; ++i)
          b[i] = 1. / a[i];
}
