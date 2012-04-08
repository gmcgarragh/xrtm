/******************************************************************************%
**
**    Copyright (C) 2007-2012 Greg McGarragh <gregm@atmos.colostate.edu>
**
**    This source code is licensed under the GNU General Public License (GPL),
**    Version 3.  See the file COPYING for more details.
**
*******************************************************************************/

#include <gutil.h>
#include <gmath_matrix.h>

#include "xrtm.h"
#include "xrtm_matrix.h"


#ifdef __cplusplus
extern "C" {
#endif

void dgemm_(const char *, const char *, long *, long *, long *, double *, double *, long *, double *, long *, double *, double *, long *);

void dsyev_(const char *, const char *, int *, double *, int *, double *, double *, int *, int *);

#ifdef __cplusplus
}
#endif


/*******************************************************************************
 *
 ******************************************************************************/
void dvec_scale_add(double *a, double *b, double *c, double alpha, long n) {

     long i;

     for (i = 0; i < n; ++i)
          c[i] = a[i] + b[i] * alpha;
}



/*******************************************************************************
 *
 ******************************************************************************/
void dmat_scale_trans(double **a, double **b,
                      double alpha, long m, long n) {

     long i;
     long j;

     if (alpha ==  1.) {
          for (i = 0; i < m; ++i) {
               for (j = 0; j < n; ++j) {
                    b[j][i] =  a[i][j];
               }
          }
     }
     else
     if (alpha == -1.) {
          for (i = 0; i < m; ++i) {
               for (j = 0; j < n; ++j) {
                    b[j][i] = -a[i][j];
               }
          }
     }
     else {
          for (i = 0; i < m; ++i) {
               for (j = 0; j < n; ++j) {
                    b[j][i] = a[i][j] * alpha;
               }
          }
     }
}



/*******************************************************************************
 *
 ******************************************************************************/
void dmat_scale_add(double **a, double **b,
                    double **c, double alpha, long m, long n) {

     long i;
     long j;

     if (alpha ==  1.) {
          for (i = 0; i < m; ++i) {
               for (j = 0; j < n; ++j) {
                    c[i][j] = a[i][j] + b[i][j];
               }
          }
     }
     else
     if (alpha == -1.) {
          for (i = 0; i < m; ++i) {
               for (j = 0; j < n; ++j) {
                    c[i][j] = a[i][j] - b[i][j];
               }
          }
     }
     else {
          for (i = 0; i < m; ++i) {
               for (j = 0; j < n; ++j) {
                    c[i][j] = a[i][j] + b[i][j] * alpha;
               }
          }
     }
}



/*******************************************************************************
 *
 ******************************************************************************/
double dmat_p_one_norm_A(double **r, double **t, long n) {

     long i;
     long j;

     double sum;
     double max = -DBL_MAX;

     for (j = 0; j < n; ++j) {
          sum = 0.;
          for (i = 0; i < n; ++i) {
               sum += fabs(t[i][j]) + fabs(r[i][j]);
          }
          if (sum > max)
               max = sum;
     }

     return max;
}



/*******************************************************************************
 *
 ******************************************************************************/
double dmat_p_inf_norm_A(double **r, double **t, long n) {

     long i;
     long j;

     double sum;
     double max = -DBL_MAX;

     for (i = 0; i < n; ++i) {
          sum = 0.;
          for (j = 0; j < n; ++j) {
               sum += fabs(t[i][j]) + fabs(r[i][j]);
          }
          if (sum > max)
               max = sum;
     }

     return max;
}



/*******************************************************************************
 *
 ******************************************************************************/
double dmat_frob_norm_A(double **r, double **t, long n) {

     long i;
     long j;

     double sum;

     sum = 0.;
     for (i = 0; i < n; ++i) {
          for (j = 0; j < n; ++j) {
               sum += r[i][j] * r[i][j] + t[i][j] * t[i][j];
          }
     }

     return sqrt(2. * sum);
}



/*******************************************************************************
 *
 ******************************************************************************/
int positive_definite(double **a, int n) {

     int i;

     int lwork = 4 * n;
     int info;

     double **b;

     double *evals;

     double *dwork;

     b = alloc_array2_d(n, n);

     evals = alloc_array1_d(n);

     dwork = alloc_array1_d(4 * n);

     dmat_copy(b, a, n, n);

     dsyev_("N", "L", &n, *b, &n, evals, dwork, &lwork, &info);
     if (info) {
          eprintf("ERROR: dsyev() info = %d\n", info);
          exit(1);
     }

     for (i = 0; i < n; ++i) {
          if (evals[i] < 0.) {
               return 0;
          }
     }

     free_array2_d(b);

     free_array1_d(evals);

     free_array1_d(dwork);

     return 1;
}



/*******************************************************************************
 *
 ******************************************************************************/
void dmat_sym(double **a, long m, double alpha) {

     long i;
     long j;
/*
     if (alpha ==  1.)
          for (i = 0; i < m; ++i) {
               for (j = i + 1; j < m; ++j) {
                    a[j][i] =  a[i][j];
               }
          }
     }
     else
     if (alpha == -1.) {
          for (i = 0; i < m; ++i) {
               for (j = i + 1; j < m; ++j) {
                    a[j][i] = -a[i][j];
               }
          }
     }
     else {
          for (i = 0; i < m; ++i) {
               for (j = i + 1; j < m; ++j) {
                    a[j][i] =  a[i][j] * alpha;
               }
          }
     }
*/
     if (alpha ==  1.) {
          for (i = 0; i < m; ++i) {
               for (j = 0; j < i + 1; ++j) {
                    a[i][j] =  a[j][i];
               }
          }
     }
     else
     if (alpha == -1.) {
          for (i = 0; i < m; ++i) {
               for (j = 0; j < i + 1; ++j) {
                    a[i][j] = -a[j][i];
               }
          }
     }
     else {
          for (i = 0; i < m; ++i) {
               for (j = 0; j < i + 1; ++j) {
                    a[i][j] = a[j][i] * alpha;
               }
          }
     }
}



/*******************************************************************************
 *
 ******************************************************************************/
void dsym_add(double **a, double **b, double **c, long m, int flag) {

     long i;
     long j;

     if (! flag) {
          for (i = 0; i < m; ++i) {
               c[i][i] = a[i][i] + b[i][i];
               for (j = i + 1; j < m; ++j) {
                    c[i][j] = a[i][j] + b[i][j];
               }
          }
     }
     else {
          for (i = 0; i < m; ++i) {
               c[i][i] = a[i][i] + b[i][i];
               for (j = i + 1; j < m; ++j) {
                    c[j][i] = c[i][j] = a[i][j] + b[i][j];
               }
          }
     }
}



/*******************************************************************************
 *
 ******************************************************************************/
void dsym_sub(double **a, double **b, double **c, long m, int flag) {

     long i;
     long j;

     if (! flag) {
          for (i = 0; i < m; ++i) {
               c[i][i] = a[i][i] - b[i][i];
               for (j = i + 1; j < m; ++j) {
                    c[i][j] = a[i][j] - b[i][j];
               }
          }
     }
     else {
          for (i = 0; i < m; ++i) {
               c[i][i] = a[i][i] - b[i][i];
               for (j = i + 1; j < m; ++j) {
                    c[j][i] = c[i][j] = a[i][j] - b[i][j];
               }
          }
     }
}



/*******************************************************************************
 *
 ******************************************************************************/
#ifdef CRAP
static void dsym_mul_kernel0(double **a, double **b, double alpha, long m, long o,
                             double **c, double beta, int flag) {

     long i;
     long j;
     long k;

     double x;

     if (alpha == 1.) {
          if (beta == 0.) {
               for (i = 0; i < m; ++i) {
                    for (j = i; j < m; ++j) {

                         x = 0.;
                         for (k = 0; k < o; ++k) {
                              x += a[i][k]*b[k][j];
                         }
                         c[i][j] = x;
                    }
               }
          }
          else
          if (beta == 1.) {
               for (i = 0; i < m; ++i) {
                    for (j = i; j < m; ++j) {
                         x = 0.;
                         for (k = 0; k < o; ++k) {
                              x += a[i][k]*b[k][j];
                         }
                         c[i][j] += x;
                    }
               }
          }
          else
          if (beta == -1.) {
               for (i = 0; i < m; ++i) {
                    for (j = i; j < m; ++j) {
                         x = 0.;
                         for (k = 0; k < o; ++k) {
                              x += a[i][k]*b[k][j];
                         }
                         c[i][j] = -c[i][j] + x;
                    }
               }
          }
          else {
               for (i = 0; i < m; ++i) {
                    for (j = i; j < m; ++j) {
                         x = 0.;
                         for (k = 0; k < o; ++k) {
                              x += a[i][k]*b[k][j];
                         }
                         c[i][j] = c[i][j] * beta + x;
                    }
               }
          }
     }
     else
     if (alpha == -1.) {
          if (beta == 0.) {
               for (i = 0; i < m; ++i) {
                    for (j = i; j < m; ++j) {
                         x = 0.;
                         for (k = 0; k < o; ++k) {
                              x += a[i][k]*b[k][j];
                         }
                         c[i][j] = -x;
                    }
               }
          }
          else
          if (beta == 1.) {
               for (i = 0; i < m; ++i) {
                    for (j = i; j < m; ++j) {
                         x = 0.;
                         for (k = 0; k < o; ++k) {
                              x += a[i][k]*b[k][j];
                         }
                         c[i][j] -= x;
                    }
               }
          }
          else
          if (beta == -1.) {
               for (i = 0; i < m; ++i) {
                    for (j = i; j < m; ++j) {
                         x = 0.;
                         for (k = 0; k < o; ++k) {
                              x += a[i][k]*b[k][j];
                         }
                         c[i][j] = -c[i][j] - x;
                    }
               }
          }
          else {
               for (i = 0; i < m; ++i) {
                    for (j = i; j < m; ++j) {
                         x = 0.;
                         for (k = 0; k < o; ++k) {
                              x += a[i][k]*b[k][j];
                         }
                         c[i][j] = c[i][j] * beta - x;
                    }
               }
          }
     }
     else {
          if (beta == 0.) {
               for (i = 0; i < m; ++i) {
                    for (j = i; j < m; ++j) {
                         x = 0.;
                         for (k = 0; k < o; ++k) {
                              x += a[i][k]*b[k][j];
                         }
                         c[i][j] = x * alpha;
                    }
               }
          }
          else
          if (beta == 1.) {
               for (i = 0; i < m; ++i) {
                    for (j = i; j < m; ++j) {
                         x = 0.;
                         for (k = 0; k < o; ++k) {
                              x += a[i][k]*b[k][j];
                         }
                         c[i][j] += x * alpha;
                    }
               }
          }
          else
          if (beta == -1.) {
               for (i = 0; i < m; ++i) {
                    for (j = i; j < m; ++j) {
                         x = 0.;
                         for (k = 0; k < o; ++k) {
                              x += a[i][k]*b[k][j];
                         }
                         c[i][j] = -c[i][j] + x * alpha;
                    }
               }
          }
          else {
               for (i = 0; i < m; ++i) {
                    for (j = i; j < m; ++j) {
                         x = 0.;
                         for (k = 0; k < o; ++k) {
                              x += a[i][k]*b[k][j];
                         }
                         c[i][j] = c[i][j] * beta + x * alpha;
                    }
               }
          }
     }

     if (flag)
          dmat_sym(c, m, 1.);
}



static void dsym_mul_kernel1(double **a, double **b, double alpha, long m, long o,
                             double **c, double beta, int flag, int thresh) {

     long io;
     long jo;
     long mo;
     long me;

     mo = (m + 1) / 2;

     if (mo <= thresh) {
          dsym_mul_kernel0(a, b, alpha, m, o, c, beta, flag);
     }
     else {
          me = m - mo;

          io = mo * o;
          jo = io + mo;

          dgemm_("n", "n", &mo, &mo, &o,
                 &alpha, *b+me, &m, *a, &o, &beta, *c+me, &m);

          dgemm_("n", "n", &me, &me, &o,
                 &alpha, *b,    &m, *a,    &o, &beta, *c,    &m);
          dgemm_("n", "n", &me, &me, &o,
                 &alpha, *b+mo, &m, *a+io, &o, &beta, *c+jo, &m);
/*
          dsym_mul_kernel1(a,    b,    alpha, me, o, c,    beta, 0);
          dsym_mul_kernel1(a+io, b+mo, alpha, me, o, c+jo, beta, 0);
*/
     }

     if (flag)
          dmat_sym(c, m, 1.);
}
#endif


void dsym_mul_kernel(double **a, double **b, double alpha, long m, long o,
                     double **c, double beta, int flag, int thresh) {

     dgemm_("n", "n", &m, &o, &o, &alpha, *b, &m, *a, &o, &beta, *c, &m);
/*
     dsym_mul_kernel0(a, b, alpha, m, o, c, beta, 0);

     dsym_mul_kernel1(a, b, alpha, m, o, c, beta, 0, thresh);

     if (flag)
          dmat_sym(c, m, 1.);
*/
}



/*******************************************************************************
 *
 ******************************************************************************/
void dsym_mul(double **a, double **b,
              long m, long o, double **c, int flag) {

     dsym_mul_kernel(a, b, 1., m, o, c, 0., flag, 0);
}



/*******************************************************************************
 *
 ******************************************************************************/
void dsym_mma(double **a1, double **b1, double **a2, double **b2,
              long m, long o, double **c, int flag) {

     dsym_mul_kernel(a1, b1, 1., m, o, c, 0., 0,    0);
     dsym_mul_kernel(a2, b2, 1., m, o, c, 1., flag, 0);
}



/*******************************************************************************
 *
 ******************************************************************************/
void dsym_mms(double **a1, double **b1, double **a2, double **b2,
              long m, long o, double **c, int flag) {

     dsym_mul_kernel(a1, b1,  1., m, o, c, 0., 0,    0);
     dsym_mul_kernel(a2, b2, -1., m, o, c, 1., flag, 0);
}



/*******************************************************************************
 *
 ******************************************************************************/
void dsym_m_a(double **a1, double **b1, double **a2,
              long m, long o, double **c, int flag) {

     dsym_mul_kernel(a1, b1, 1., m, o, c, 0., 0, 0);

     dsym_add(c, a2, c, m, flag);
}



/*******************************************************************************
 *
 ******************************************************************************/
void dsym_m_s(double **a1, double **b1, double **a2,
              long m, long o, double **c, int flag) {

     dsym_mul_kernel(a1, b1, 1., m, o, c, 0., 0, 0);

     dsym_sub(c, a2, c, m, flag);
}



/*******************************************************************************
 *
 ******************************************************************************/
void dsym_mta(double **a, double **b, long m, long o, double **c) {

     dmat_mul(a, b, m, m, o, c);
     dmat_add_trans(c, c, m, m);
}



/*******************************************************************************
 *
 ******************************************************************************/
void dsym_mts(double **a, double **b, long m, long o, double **c) {

     dmat_mul(a, b, m, m, o, c);
     dmat_sub_trans(c, c, m, m);
}



/*******************************************************************************
 *
 ******************************************************************************/
void zmat_real(dcomplex **z, double **d, long m, long n) {

     long i;
     long j;

     for (i = 0; i < m; ++i) {
          for (j = 0; j < n; ++j) {
               d[i][j] = creal(z[i][j]);
          }
     }
}



void zmat_imag(dcomplex **z, double **d, long m, long n) {

     long i;
     long j;

     for (i = 0; i < m; ++i) {
          for (j = 0; j < n; ++j) {
               d[i][j] = cimag(z[i][j]);
          }
     }
}



void dmat_complex(double **dr, double **di, dcomplex **z, long m, long n) {

     long i;
     long j;

     for (i = 0; i < m; ++i) {
          for (j = 0; j < n; ++j) {
               z[i][j] = dr[i][j] + _Complex_I * di[i][j];
          }
     }
}



void dzmat_diag_mul(double *a, dcomplex **b, dcomplex **c, long m, long n) {

     long i;
     long j;

     for (i = 0; i < m; ++i) {
          for (j = 0; j < n; ++j) {
               c[i][j] = a[i] * b[i][j];
          }
     }
}



void dzmat_mul(double **a, dcomplex **b, long m, long n, long o, dcomplex **c) {

     long i;
     long j;
     long k;

     dcomplex x;

     for (i = 0; i < m; ++i) {
          for (j = 0; j < n; ++j) {
               x = 0.;
               for (k = 0; k < o; ++k) {
                    x += a[i][k]*b[k][j];
               }
               c[i][j] = x;
          }
     }
}



void dzmat_mul2(double **a, dcomplex **b, long m, long n, long o, dcomplex **c,
                double **don, double **dmn1, double **dmn2) {

     zmat_real(b, don, o, n);
     dmat_mul(a, don, m, n, o, dmn1);

     zmat_imag(b, don, o, n);
     dmat_mul(a, don, m, n, o, dmn2);

     dmat_complex(dmn1, dmn2, c, m, n);
}



/*******************************************************************************
 *
 ******************************************************************************/
void copy_to_band_storage_d(double **ab, double **a, double alpha,
                            int m, int n, int kl, int ku, int i1, int j1) {

     int ii;
     int i;
     int j;
     int jj;
     int k;

     k = kl + ku + 1;

     ii = i1 + k - 1;

     if (alpha ==  1.) {
          for (i = 0; i < m; ++i) {
               jj = j1;
               for (j = 0; j < n; ++j) {
                    ab[jj][ii - jj] =  a[i][j];
                    jj++;
               }
               ii++;
          }
     }
     else
     if (alpha == -1.) {
          for (i = 0; i < m; ++i) {
               jj = j1;
               for (j = 0; j < n; ++j) {
                    ab[jj][ii - jj] = -a[i][j];
                    jj++;
               }
               ii++;
          }
     }
     else {
          for (i = 0; i < m; ++i) {
               jj = j1;
               for (j = 0; j < n; ++j) {
                    ab[jj][ii - jj] =  a[i][j] * alpha;
                    jj++;
               }
               ii++;
          }
     }
}



void copy_to_band_storage2_d(double **ab, double **a, double alpha, double *beta,
                             int m, int n, int kl, int ku, int i1, int j1) {

     int ii;
     int i;
     int j;
     int jj;
     int k;

     k = kl + ku + 1;

     ii = i1 + k - 1;

     if (alpha ==  1.) {
          for (i = 0; i < m; ++i) {
               jj = j1;
               for (j = 0; j < n; ++j) {
                    ab[jj][ii - jj] =  a[i][j] * beta[j];
                    jj++;
               }
               ii++;
          }
     }
     else
     if (alpha == -1.) {
          for (i = 0; i < m; ++i) {
               jj = j1;
               for (j = 0; j < n; ++j) {
                    ab[jj][ii - jj] = -a[i][j] * beta[j];
                    jj++;
               }
               ii++;
          }
     }
     else {
          for (i = 0; i < m; ++i) {
               jj = j1;
               for (j = 0; j < n; ++j) {
                    ab[jj][ii - jj] =  a[i][j] * alpha * beta[j];
                    jj++;
               }
               ii++;
          }
     }
}



/*******************************************************************************
 *
 ******************************************************************************/
void copy_to_band_storage_dc(dcomplex **ab, dcomplex **a, dcomplex alpha,
                             int m, int n, int kl, int ku, int i1, int j1) {

     int ii;
     int i;
     int j;
     int jj;
     int k;

     k = kl + ku + 1;

     ii = i1 + k - 1;

     if (alpha ==  1.) {
          for (i = 0; i < m; ++i) {
               jj = j1;
               for (j = 0; j < n; ++j) {
                    ab[jj][ii - jj] =  a[i][j];
                    jj++;
               }
               ii++;
          }
     }
     else
     if (alpha == -1.) {
          for (i = 0; i < m; ++i) {
               jj = j1;
               for (j = 0; j < n; ++j) {
                    ab[jj][ii - jj] = -a[i][j];
                    jj++;
               }
               ii++;
          }
     }
     else {
          for (i = 0; i < m; ++i) {
               jj = j1;
               for (j = 0; j < n; ++j) {
                    ab[jj][ii - jj] =  a[i][j] * alpha;
                    jj++;
               }
               ii++;
          }
     }
}



void copy_to_band_storage2_dc(dcomplex **ab, dcomplex **a,
                              dcomplex alpha, dcomplex *beta,
                              int m, int n, int kl, int ku, int i1, int j1) {

     int ii;
     int i;
     int j;
     int jj;
     int k;

     k = kl + ku + 1;

     ii = i1 + k - 1;

     if (alpha ==  1.) {
          for (i = 0; i < m; ++i) {
               jj = j1;
               for (j = 0; j < n; ++j) {
                    ab[jj][ii - jj] =  a[i][j] * beta[j];
                    jj++;
               }
               ii++;
          }
     }
     else
     if (alpha == -1.) {
          for (i = 0; i < m; ++i) {
               jj = j1;
               for (j = 0; j < n; ++j) {
                    ab[jj][ii - jj] = -a[i][j] * beta[j];
                    jj++;
               }
               ii++;
          }
     }
     else {
          for (i = 0; i < m; ++i) {
               jj = j1;
               for (j = 0; j < n; ++j) {
                    ab[jj][ii - jj] =  a[i][j] * alpha * beta[j];
                    jj++;
               }
               ii++;
          }
     }
}



/*******************************************************************************
 *
 ******************************************************************************/
#ifdef INCLUDE_DEV_SOURCE
#include "../src2/xrtm_matrix2.c"
#endif
