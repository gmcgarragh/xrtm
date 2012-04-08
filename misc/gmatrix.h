/******************************************************************************%
**
**    Copyright (C) 1998-2012 Greg McGarragh <gregm@atmos.colostate.edu>
**
**    This source code is licensed under the GNU General Public License (GPL),
**    Version 3.  See the file COPYING for more details.
**
*******************************************************************************/

#include "gutil.h"


#ifdef USE_BLAS

#ifdef __cplusplus
extern "C" {
#endif

void XCAT(prefix_, gemv_)(const char *, long *, long *, type_ *, type_ *,
                          long *, type_ *, long *, type_ *, type_ *, long *);
void XCAT(prefix_, gemm_)(const char *, const char *, long *, long *, long *,
                          type_ *, type_ *, long *, type_ *, long *, type_ *,
                          type_ *, long *);

#ifdef __cplusplus
}
#endif

#endif

void XCAT(prefix_, mat_zero)(type_ **a, long m, long n) {
/*
     long i;
     long j;

     for (i = 0; i < m; ++i) {
          for (j = 0; j < n; ++j) {
               a[i][j] = 0.0;
          }
     }
*/
     long i;

     size_t size = n * sizeof(type_);

     for (i = 0; i < m; ++i)
          memset(*a+i*n, 0, size);
}



void XCAT(prefix_, mat_init)(type_ **a, double alpha, long m, long n) {

     long o;

     long i;
     long j;

     o = MIN(m, n);

     for (i = 0; i < o; ++i) {
          for (j = 0; j < n; ++j) {
               a[i][j] = 0.0;
          }
          a[i][i] = alpha;
     }

     for (     ; i < m; ++i) {
          for (j = 0; j < n; ++j) {
               a[i][j] = 0.0;
          }
     }
}



void XCAT(prefix_, mat_copy)(type_ **a2, type_ **a1, long m, long n) {
/*
     long i;
     long j;

     for (i = 0; i < m; ++i) {
          for (j = 0; j < n; ++j) {
               a2[i][j] = a1[i][j];
          }
     }
*/
     long i;

     size_t size = n * sizeof(type_);

     for (i = 0; i < m; ++i)
          memcpy(*a2+i*n, *a1+i*n, size);
}



void XCAT(prefix_, mat_trans)(type_ **a, type_ **b, long m, long n) {

     long i;
     long j;

     for (i = 0; i < m; ++i) {
          for (j = 0; j < n; ++j) {
               b[j][i] = a[i][j];
          }
     }
}



double XCAT(prefix_, mat_p_one_norm)(type_ **a, long m, long n) {

     long i;
     long j;

     double sum;
     double max = -DBL_MAX;

     for (j = 0; j < n; ++j) {
          sum = 0.;
          for (i = 0; i < m; ++i) {
               sum += XCAT(cprefix_, abs)(a[i][j]);
          }
          if (sum > max)
               max = sum;
     }

     return max;
}



double XCAT(prefix_, mat_p_inf_norm)(type_ **a, long m, long n) {

     long i;
     long j;

     double sum;
     double max = -DBL_MAX;

     for (i = 0; i < n; ++i) {
          sum = 0.;
          for (j = 0; j < m; ++j) {
               sum += XCAT(cprefix_, abs)(a[i][j]);
          }
          if (sum > max)
               max = sum;
     }

     return max;
}



double XCAT(prefix_, mat_frob_norm)(type_ **a, long m, long n) {

     long i;
     long j;

     double sum;

     sum = 0.;
     for (i = 0; i < m; ++i) {
          for (j = 0; j < n; ++j) {
               sum += (double) xreal(a[i][j]) * (double) xreal(a[i][j]);
          }
     }

     return sqrt(sum);
}



void XCAT(prefix_, mat_add)(type_ **a, type_ **b, type_ **c, long m, long n) {

     long i;
     long j;

     for (i = 0; i < m; ++i) {
          for (j = 0; j < n; ++j) {
               c[i][j] = a[i][j] + b[i][j];
          }
     }
}



void XCAT(prefix_, mat_add_trans)(type_ **a, type_ **c, long m, long n) {

     long i;
     long j;

     for (i = 0; i < m; ++i) {
          c[i][i] += c[i][i];

          for (j = i + 1; j < n; ++j) {
               c[i][j] = c[j][i] = a[i][j] + a[j][i];
          }
     }
}



void XCAT(prefix_, mat_add_diag)(type_ **a, type_ *b, type_ **c, long m) {

     long i;
     long j;

     for (i = 0; i < m; ++i) {
          for (j = 0; j < m; ++j) {
               c[i][j] = a[i][j];
          }

          c[i][i] += b[i];
     }
}



void XCAT(prefix_, mat_sub)(type_ **a, type_ **b, type_ **c, long m, long n) {

     long i;
     long j;

     for (i = 0; i < m; ++i) {
          for (j = 0; j < n; ++j) {
               c[i][j] = a[i][j] - b[i][j];
          }
     }
}



void XCAT(prefix_, mat_sub_trans)(type_ **a, type_ **c, long m, long n) {

     long i;
     long j;

     for (i = 0; i < m; ++i) {
          c[i][i] -= c[i][i];

          for (j = i + 1; j < n; ++j) {
               c[i][j] = c[j][i] = a[i][j] - a[j][i];
          }
     }
}



void XCAT(prefix_, mat_sub_diag)(type_ **a, type_ *b, type_ **c, long m) {

     long i;
     long j;

     for (i = 0; i < m; ++i) {
          for (j = 0; j < m; ++j) {
               c[i][j] = a[i][j];
          }

          c[i][i] -= b[i];
     }
}



void XCAT(prefix_, mat_diag_sub)(type_ *a, type_ **b, type_ **c, long m) {

     long i;
     long j;

     for (i = 0; i < m; ++i) {
          for (j = 0; j < m; ++j) {
               c[i][j] = -b[i][j];
          }

          c[i][i] += a[i];
    }
}



void XCAT(prefix_, mat_i_sub)(type_ **b, type_ **c, long m) {

     long i;
     long j;

     for (i = 0; i < m; ++i) {
          for (j = 0; j < m; ++j) {
               c[i][j] = -b[i][j];
          }

          c[i][i] += 1.;
     }
}



void XCAT(prefix_, mat_scale)(type_ alpha, type_ **b, type_ **c, long m, long n) {

     long i;
     long j;

     for (i = 0; i < m; ++i) {
          for (j = 0; j < n; ++j) {
               c[i][j] = b[i][j] * alpha;
          }
     }
}



void XCAT(prefix_, m_v_mul)(type_ **a, type_ *b, long m, long o, type_ *c) {
#ifdef USE_BLAS
     long incb = 1;
     long incy = 1;

     type_ alpha = 1.;
     type_ beta  = 0.;

     XCAT(prefix_, gemv_)("t", &o, &m, &alpha, *a, &o, b, &incb, &beta, c, &incy);
#else
     long i;
     long k;

     type_ x;

     for (i = 0; i < m; ++i) {
               x = 0.;
               for (k = 0; k < o; ++k) {
                    x += a[i][k]*b[k];
               }
               c[i] = x;
     }
#endif
}



void XCAT(prefix_, mat_gxvxmx)(int trans_a, type_ **a, type_ *b, type_ alpha, type_ *c, type_ beta, long m, long o) {
#ifdef USE_BLAS
     char transa = trans_a ? 'n' : 't';

     long incb = 1;
     long incy = 1;

     XCAT(prefix_, gemv_)(&transa, &o, &m, &alpha, *a, &o, b, &incb, &beta, c, &incy);
#else

#endif
}



void XCAT(prefix_, v_m_mul)(type_ *a, type_ **b, long m, type_ **c) {

     XCAT(prefix_, mat_mul)(&a, b, m, m, 1, c);
}



void XCAT(prefix_, m_v_diag_mul)(type_ *a, type_ *b, type_ *c, long m) {

     long i;

     for (i = 0; i < m; ++i) {
          c[i] = a[i] * b[i];
     }
}



void XCAT(prefix_, m_v_dinv_mul)(type_ *a, type_ *b, type_ *c, long m) {

     long i;

     for (i = 0; i < m; ++i) {
          c[i] = b[i] / a[i];
     }
}



void XCAT(prefix_, mat_mul)(type_ **a, type_ **b,
                            long m, long n, long o, type_ **c) {

#ifdef USE_BLAS
     type_ alpha = 1.;
     type_ beta  = 0.;

     XCAT(prefix_, gemm_)("n", "n", &n, &m, &o, &alpha,
                          (type_ *) *b, &n, (type_ *) *a, &o, &beta, (type_ *) *c, &n);
#else
     long i;
     long j;
     long k;

     type_ x;

     type_ (*a2)[lda] = (type_ (*)[lda]) a;
     type_ (*b2)[ldb] = (type_ (*)[ldb]) b;
     type_ (*c2)[ldc] = (type_ (*)[ldc]) c;

     for (i = 0; i < m; ++i) {
          for (j = 0; j < n; ++j) {
               x = 0.;
               for (k = 0; k < o; ++k) {
                    x += a2[i][k]*b2[k][j];
               }
               c2[i][j] = x;
          }
     }
#endif
}



void XCAT(prefix_, mat_gxgxmx)(int trans_a, type_ **a, int trans_b, type_ **b,
                               type_ alpha, type_ **c, type_ beta, long m, long n, long o) {

#ifdef USE_BLAS
     char transa = trans_a ? 't' : 'n';
     char transb = trans_b ? 't' : 'n';

     XCAT(prefix_, gemm_)(&transb, &transa, &n, &m, &o, &alpha,
                          (type_ *) *b, &n, (type_ *) *a, &o, &beta, (type_ *) *c, &n);
#else
     long i;
     long j;
     long k;

     type_ x;

     type_ (*a2)[lda] = (type_ (*)[lda]) a;
     type_ (*b2)[ldb] = (type_ (*)[ldb]) b;
     type_ (*c2)[ldc] = (type_ (*)[ldc]) c;
/*
     for (i = 0; i < m; ++i) {
          for (j = 0; j < n; ++j) {
               x = 0.;
               for (k = 0; k < o; ++k) {
                    x += a2[i][k]*b2[k][j];
               }
               c2[i][j] = x;
          }
     }
*/
     if (! trans_a) {
          if (! trans_b) {
               for (i = 0; i < m; ++i) {
                    for (j = 0; j < n; ++j) {
                         x = 0.;
                         for (k = 0; k < o; ++k) {
                              x += a2[i][k]*b2[k][j];
                         }
                         c2[i][j] = beta * c2[i][j] + alpha * x;
                    }
               }
          }
          else {
               for (i = 0; i < m; ++i) {
                    for (j = 0; j < n; ++j) {
                         x = 0.;
                         for (k = 0; k < o; ++k) {
                              x += a2[i][k]*b2[j][k];
                         }
                         c2[i][j] = beta * c2[i][j] + alpha * x;
                    }
               }
          }
     }
     else {
          if (! trans_b) {
               for (i = 0; i < m; ++i) {
                    for (j = 0; j < n; ++j) {
                         x = 0.;
                         for (k = 0; k < o; ++k) {
                              x += a2[k][i]*b2[k][j];
                         }
                         c2[i][j] = beta * c2[i][j] + alpha * x;
                    }
               }
          }
          else {
               for (i = 0; i < m; ++i) {
                    for (j = 0; j < n; ++j) {
                         x = 0.;
                         for (k = 0; k < o; ++k) {
                              x += a2[k][i]*b2[j][k];
                         }
                         c2[i][j] = beta * c2[i][j] + alpha * x;
                    }
               }
          }
     }
#endif
}



void XCAT(prefix_, mat_mul_diag)(type_ **a, type_ *b, type_ **c,
                                 long m, long n) {

     long i;
     long j;

     for (i = 0; i < m; ++i) {
          for (j = 0; j < n; ++j) {
               c[i][j] = a[i][j] * b[j];
          }
     }
}



void XCAT(prefix_, mat_gxdxmx)(int trans_a, type_ **a, type_ *b, type_ alpha, type_ **c, type_ beta, long m, long n) {

     long i;
     long j;

     if (! trans_a) {
/*
          if (alpha ==  1.) {
               if (beta ==  0.) {
               }
               else
               if (beta ==  1.) {
               }
               else
               if (beta == -1.) {
               }
               else {

               }
          }
          else
          if (alpha == -1.) {
               if (beta ==  0.) {
               }
               else
               if (beta ==  1.) {
               }
               else
               if (beta == -1.) {
               }
               else {

               }
          }
          else {
               if (beta ==  0.) {
               }
               else
               if (beta ==  1.) {
               }
               else
               if (beta == -1.) {
               }
               else {
*/
                    for (i = 0; i < m; ++i) {
                         for (j = 0; j < n; ++j) {
                              c[i][j] = beta * c[i][j] + alpha * a[i][j] * b[j];
                         }
                    }
/*
               }
          }
*/
     }
     else {
/*
          if (alpha ==  1.) {
               if (beta ==  0.) {
               }
               else
               if (beta ==  1.) {
               }
               else
               if (beta == -1.) {
               }
               else {

               }
          }
          else
          if (alpha == -1.) {
               if (beta ==  0.) {
               }
               else
               if (beta ==  1.) {
               }
               else
               if (beta == -1.) {
               }
               else {

               }
          }
          else {
               if (beta ==  0.) {
               }
               else
               if (beta ==  1.) {
               }
               else
               if (beta == -1.) {
               }
               else {
*/
                    for (i = 0; i < m; ++i) {
                         for (j = 0; j < n; ++j) {
                              c[i][j] = beta * c[i][j] + alpha * a[j][i] * b[j];
                         }
                    }
/*
               }
          }
*/
     }
}



void XCAT(prefix_, mat_mul_dinv)(type_ **a, type_ *b, type_ **c,
                                 long m, long n) {

     long i;
     long j;

     for (i = 0; i < m; ++i) {
          for (j = 0; j < n; ++j) {
               c[i][j] = a[i][j] / b[j];
          }
     }
}



void XCAT(prefix_, mat_diag_mul)(type_ *a, type_ **b, type_ **c,
                                 long m, long n) {

     long i;
     long j;

     for (i = 0; i < m; ++i) {
          for (j = 0; j < n; ++j) {
               c[i][j] = a[i] * b[i][j];
          }
     }
}



void XCAT(prefix_, mat_dinv_mul)(type_ *a, type_ **b, type_ **c,
                                 long m, long n) {

     long i;
     long j;

     for (i = 0; i < m; ++i) {
          for (j = 0; j < n; ++j) {
               c[i][j] = b[i][j] / a[i];
          }
     }
}



void XCAT(prefix_, mat_vxvtmx)(type_ *a, type_ *b, double alpha,
                               type_ **c, double beta, long m, long n) {

     int i;
     int j;

     type_ aa;

     if (alpha ==  1.) {
          if (beta ==  0.) {
               for (i = 0; i < m; ++i) {
                    for (j = 0; j < n; ++j) {
                         c[i][j] = a[i] * b[j];
                    }
               }
          }
          else
          if (beta ==  1.) {
               for (i = 0; i < m; ++i) {
                    for (j = 0; j < n; ++j) {
                         c[i][j] += a[i] * b[j];
                    }
               }
          }
          else
          if (beta == -1.) {
               for (i = 0; i < m; ++i) {
                    for (j = 0; j < n; ++j) {
                         c[i][j] = -c[i][j] + a[i] * b[j];
                    }
               }
          }
          else {
               for (i = 0; i < m; ++i) {
                    for (j = 0; j < n; ++j) {
                         c[i][j] = beta * c[i][j] + a[i] * b[j];
                    }
               }
          }
     }
     else
     if (alpha == -1.) {
          if (beta ==  0.) {
               for (i = 0; i < m; ++i) {
                    for (j = 0; j < n; ++j) {
                         c[i][j] = -a[i] * b[j];
                    }
               }
          }
          else
          if (beta ==  1.) {
               for (i = 0; i < m; ++i) {
                    for (j = 0; j < n; ++j) {
                         c[i][j] -= a[i] * b[j];
                    }
               }
          }
          else
          if (beta == -1.) {
               for (i = 0; i < m; ++i) {
                    for (j = 0; j < n; ++j) {
                         c[i][j] = -c[i][j] - a[i] * b[j];
                    }
               }
          }
          else {
               for (i = 0; i < m; ++i) {
                    for (j = 0; j < n; ++j) {
                         c[i][j] = beta * c[i][j] - a[i] * b[j];
                    }
               }
          }
     }
     else {
          if (beta ==  0.) {
               for (i = 0; i < m; ++i) {
                    aa = alpha * a[i];
                    for (j = 0; j < n; ++j) {
                         c[i][j] = aa * b[j];
                    }
               }
          }
          else
          if (beta ==  1.) {
               for (i = 0; i < m; ++i) {
                    aa = alpha * a[i];
                    for (j = 0; j < n; ++j) {
                         c[i][j] = c[i][j] + aa * b[j];
                    }
               }
          }
          else
          if (beta == -1.) {
               for (i = 0; i < m; ++i) {
                    aa = alpha * a[i];
                    for (j = 0; j < n; ++j) {
                         c[i][j] = -c[i][j] + aa * b[j];
                    }
               }
          }
          else {
               for (i = 0; i < m; ++i) {
                    aa = alpha * a[i];
                    for (j = 0; j < n; ++j) {
                         c[i][j] = beta * c[i][j] + aa * b[j];
                    }
               }
          }
     }
}



int XCAT(prefix_, mat_pow_count)(int s) {

     int p;
     int t;

     int count;

     t = max_pow_of_two[s];

     s -= powers_of_two[t];
     count = t;

     while (s > 0) {
          p = max_pow_of_two[s];
          s -= powers_of_two[p];
          count++;
     }

     return count;
}



void XCAT(prefix_, mat_pow)(type_ **a, int s, type_ **c,
                            long m, long n, type_ **w1, type_ **w2) {

     int i;
     int j;

     int p;
     int t;

     char b[127];

     type_ **w[] = {w1, w2};

     if (s == 0) {
          XCAT(prefix_, mat_init)(c, 1., m, n);
          return;
     }

     t = max_pow_of_two[s];

     memset(b, 0, t + 1);

     b[t] = 1;
     s -= powers_of_two[t];

     while (s > 0) {
          p = max_pow_of_two[s];
          b[p] = 1;
          s -= powers_of_two[p];
     }

     XCAT(prefix_, mat_copy)(w[0], a, m, n);

     j = 0;
     for (i  = 0; i <= t; ++i) {
          if (b[i] != 0)
               break;

          XCAT(prefix_, mat_mul)(w[j], w[j], m, n, n,  w[1 - j]);

          j = 1 - j;
     }

     XCAT(prefix_, mat_copy)(c, w[j], m, n);

     for (i += 1; i <= t; ++i) {
          XCAT(prefix_, mat_mul)(w[j], w[j], m, n, n,  w[1 - j]);

          if (b[i]  != 0) {
               XCAT(prefix_, mat_mul)(c, w[1 - j], m, n, n,  w[j]);
               XCAT(prefix_, mat_copy)(c, w[j], m, n);
          }

          j = 1 - j;
     }
}



#ifdef USE_LAPACK

#ifdef __cplusplus
extern "C" {
#endif

void XCAT(prefix_, pocon_)(const char *uplo, int *n, type_ *a, int *lda, double *anorm, double *rcond, double *work, int *iwork, int *info);
void XCAT(prefix_, gecon_)(char *norm, int *n, type_ *a, int *lda, double *anorm, double *rcond, double *work, int *iwork, int *info);
void XCAT(prefix_, potrf_)(const char *, int *, type_ *, int *, int *);
void XCAT(prefix_, potri_)(const char *, int *, type_ *, int *, int *);
void XCAT(prefix_, potrs_)(const char *, int *, int *, type_ *, int *, type_ *, int *, int *);
void XCAT(prefix_, getrf_)(int *, int *, type_ *, int *, int *, int *);
void XCAT(prefix_, getri_)(int *, type_ *, int *, int *, type_ *, int *, int *);
void XCAT(prefix_, getrs_)(const char *, int *, int *, type_ *, int *, int *, type_ *, int *, int *);
void XCAT(prefix_, gttrf_)(int *, type_ *, type_ *, type_ *, type_ *, int *, int *);
void XCAT(prefix_, gttrs_)(const char *, int *, int *, type_ *, type_ *, type_ *, type_ *, int *, type_ *, int *, int *);

#ifdef __cplusplus
}
#endif

double XCAT(prefix_, mat_pocon)(type_ **a, int n, double anorm, int *iwork, double *dwork) {

     int info;

     double rcond;

     XCAT(prefix_, pocon_)("l", &n, *a, &n, &anorm, &rcond, dwork, iwork, &info);
     if (info) {
          printf("ERROR: dpocon() info = %d\n", info);
          exit(1);
     }

     return rcond;
}



double XCAT(prefix_, mat_gecon)(char norm, type_ **a, int n, double anorm, int *iwork, double *dwork) {

     int info;

     double rcond;

     XCAT(prefix_, gecon_)(&norm, &n, *a, &n, &anorm, &rcond, dwork, iwork, &info);
     if (info) {
          printf("ERROR: dgecon() info = %d\n", info);
          exit(1);
     }

     return rcond;
}



void XCAT(prefix_, mat_potrf)(type_ **a, int n) {

     int info;

     XCAT(prefix_, potrf_)("l", &n, *a, &n, &info);
     if (info) {
          printf("ERROR: dpotrf() info = %d\n", info);
          exit(1);
     }
}



void XCAT(prefix_, mat_potri)(type_ **a, int n) {

     int info;

     XCAT(prefix_, potri_)("l", &n, *a, &n, &info);
     if (info) {
          printf("ERROR: dpotri() info = %d\n", info);
          exit(1);
     }
}



void XCAT(prefix_, mat_potrs)(type_ **a, type_ **b, int n, int nrhs) {

     int info;

     XCAT(prefix_, potrs_)("l", &n, &nrhs, *a, &n, *b, &nrhs, &info);

     if (info) {
          printf("ERROR: dpotrs() info = %d\n", info);
          exit(1);
     }
}



void XCAT(prefix_, mat_getrf)(type_ **a, int m, int n, int *ipiv) {

     int info;

     XCAT(prefix_, getrf_)(&m, &n, *a, &n, ipiv, &info);
     if (info) {
          printf("ERROR: dgetrf() info = %d\n", info);
          exit(1);
     }
}



void XCAT(prefix_, mat_getri)(type_ **a, int n, int *ipiv) {

     int lwork;

     int info;

     double size[2];

     type_ *twork;

     lwork = -1;
     XCAT(prefix_, getri_)(&n, *a, &n, ipiv, (type_ *) size, &lwork, &info);

     lwork = (int) size[0];
     twork = (type_ *) malloc((size_t) size[0] * sizeof(type_));
     XCAT(prefix_, getri_)(&n, *a, &n, ipiv, twork, &lwork, &info);
     if (info) {
          printf("ERROR: dgetri() info = %d\n", info);
          exit(1);
     }

     free(twork);
}



void XCAT(prefix_, mat_getrs)(type_ **a, type_ **b, int n, int nrhs, int *ipiv) {

     XCAT(prefix_, mat_getrs2)('n', a, b, n, nrhs, ipiv);
}



void XCAT(prefix_, mat_getrs2)(char trans, type_ **a, type_ **b, int n, int nrhs, int *ipiv) {

     int info;

     trans = trans == 'n' ? 't' : 'n';

     XCAT(prefix_, getrs_)(&trans, &n, &nrhs, *a, &n, ipiv, *b, &n, &info);

     if (info) {
          printf("ERROR: dgetrs() info = %d\n", info);
          exit(1);
     }
}



void XCAT(prefix_, mat_gttrf)(type_ *dl, type_ *d, type_ *du, type_ *du2, int n, int *ipiv) {

     int info;

     XCAT(prefix_, gttrf_)(&n, dl, d, du, du2, ipiv, &info);
     if (info) {
          printf("ERROR: dgttrf() info = %d\n", info);
          exit(1);
     }
}



void XCAT(prefix_, mat_gttrs)(type_ *dl, type_ *d, type_ *du, type_ *du2, type_ **b, int n, int nrhs, int *ipiv) {

     int info;

     XCAT(prefix_, gttrs_)("n", &n, &nrhs, dl, d, du, du2, ipiv, *b, &n, &info);

     if (info) {
          printf("ERROR: dgttrs() info = %d\n", info);
          exit(1);
     }
}
#endif
