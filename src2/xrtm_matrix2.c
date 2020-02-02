/*******************************************************************************
**
**    Copyright (C) 2007-2020 Greg McGarragh <greg.mcgarragh@colostate.edu>
**
**    This source code is licensed under the GNU General Public License (GPL),
**    Version 3.  See the file COPYING for more details.
**
*******************************************************************************/


#ifdef __cplusplus
extern "C" {
#endif

void dscal_(long *, double *, double *, long *);
void daxpy_(long *, double *, double *, long *, double *, long *);

void dtrmm_(const char *, const char *, const char *, const char *, long *, long *, double *, double *, long *, double *, long *);

void dgemv_(const char *, long *, long *, double *, double *, long *, double *, long *, double *, double *, long *);

#ifdef __cplusplus
}
#endif


/*******************************************************************************
 *
 ******************************************************************************/
void dmat_dtutumx(double **a, double **b, double **c, long n) {

}



/*******************************************************************************
 *
 ******************************************************************************/
void dmat_dtutlmx(double **a, double **b, double **c, long n) {

     int i;
     int j;
     int k;

     double x;

     for (i = 0; i < n; ++i) {
          for (j = 0; j < i; ++j) {
               x = 0.;
               for (k = i; k < n; ++k) {
                    x += a[i][k]*b[k][j];
               }
               c[i][j] = x;
          }
          for (j = i; j < n; ++j) {
               x = 0.;
               for (k = j; k < n; ++k) {
                    x += a[i][k]*b[k][j];
               }
               c[i][j] = x;
          }
     }
}



void dmat_dtltumx(double **a, double **b, double **c, long n) {

     int i;
     int j;
     int k;

     double x;

     for (i = 0; i < n; ++i) {
          for (j = 0; j < i; ++j) {
               x = 0.;
               for (k = 0; k <= j; ++k) {
                    x += a[i][k]*b[k][j];
               }
               c[i][j] = x;
          }
          for (j = i; j < n; ++j) {
               x = 0.;
               for (k = 0; k <= i; ++k) {
                    x += a[i][k]*b[k][j];
               }
               c[i][j] = x;
          }
     }
}



/*******************************************************************************
 *
 ******************************************************************************/
void dmat_dgxtumx(double **a, double **b, double **c, long n) {
#ifdef USE_BLAS
     double alpha = 1.;

     dmat_copy(c, a, n, n);

     dtrmm_("l", "l", "n", "n", &n, &n, &alpha, *b, &n, *c, &n);
#else
     int i;
     int j;
     int k;

     double x;

     for (i = 0; i < n; ++i) {
          for (j = 0; j < n; ++j) {
               x = 0.;
               for (k = 0; k <= j; ++k) {
                    x += a[i][k]*b[k][j];
               }
               c[i][j] = x;
          }
     }
#endif
}



void dmat_dgxtlmx(double **a, double **b, double **c, long n) {
#ifdef USE_BLAS
     double alpha = 1.;

     dmat_copy(c, a, n, n);

     dtrmm_("l", "u", "n", "n", &n, &n, &alpha, *b, &n, *c, &n);
#else
     int i;
     int j;
     int k;

     double x;

     for (i = 0; i < n; ++i) {
          for (j = 0; j < n; ++j) {
               x = 0.;
               for (k = j; k < n; ++k) {
                    x += a[i][k]*b[k][j];
               }

               c[i][j] = x;
          }
     }
#endif
}



void dmat_dtugxmx(double **a, double **b, double **c, long n) {
#ifdef USE_BLAS
     double alpha = 1.;

     dmat_copy(c, b, n, n);

     dtrmm_("r", "l", "n", "n", &n, &n, &alpha, *a, &n, *c, &n);
#else
     int i;
     int j;
     int k;

     double x;

     for (i = 0; i < n; ++i) {
          for (j = 0; j < n; ++j) {
               x = 0.;
               for (k = i; k < n; ++k) {
                    x += a[i][k]*b[k][j];
               }

               c[i][j] = x;
          }
     }
#endif
}



void dmat_dtlgxmx(double **a, double **b, double **c, long n) {
#ifdef USE_BLAS
     double alpha = 1.;

     dmat_copy(c, b, n, n);

     dtrmm_("r", "u", "n", "n", &n, &n, &alpha, *a, &n, *c, &n);
#else
     int i;
     int j;
     int k;

     double x;

     for (i = 0; i < n; ++i) {
          for (j = 0; j < n; ++j) {
               x = 0.;
               for (k = 0; k <= i; ++k) {
                    x += a[i][k]*b[k][j];
               }
               c[i][j] = x;
          }
     }
#endif
}



/*******************************************************************************
 *
 ******************************************************************************/
void dmat_dtutlms(double **a, double **b, double **c, int n, int flag) {

     int i;
     int j;
     int k;

     double x;

     for (i = 0; i < n; ++i) {
          for (j = i; j < n; ++j) {
               x = 0.;
               for (k = j; k < n; ++k) {
                    x += a[i][k]*b[k][j];
               }
               c[i][j] = x;
               if (flag)
                    c[j][i] = x;
          }
     }
}



void dmat_dtltums(double **a, double **b, double **c, int n, int flag) {

     int i;
     int j;
     int k;

     double x;

     for (i = 0; i < n; ++i) {
          for (j = i; j < n; ++j) {
               x = 0.;
               for (k = 0; k <= i; ++k) {
                    x += a[i][k]*b[k][j];
               }
               c[i][j] = x;
               if (flag)
                    c[j][i] = x;
          }
     }
}



/*******************************************************************************
 *
 ******************************************************************************/
void dmat_dgxtums(double **a, double **b, double **c, int n, int flag) {

     int i;
     int j;
     int k;

     double x;

     for (i = 0; i < n; ++i) {
          for (j = i; j < n; ++j) {
               x = 0.;
               for (k = 0; k <= j; ++k) {
                    x += a[i][k]*b[k][j];
               }
               c[i][j] = x;
               if (flag)
                    c[j][i] = x;
          }
     }
}



void dmat_dgxtlms(double **a, double **b, double **c, int n, int flag) {

     int i;
     int j;
     int k;

     double x;

     for (i = 0; i < n; ++i) {
          for (j = i; j < n; ++j) {
               x = 0.;
               for (k = j; k < n; ++k) {
                    x += a[i][k]*b[k][j];
               }

               c[i][j] = x;
               if (flag)
                    c[j][i] = x;
          }
     }
}



void dmat_dtugxms(double **a, double **b, double **c, int n, int flag) {

     int i;
     int j;
     int k;

     double x;

     for (i = 0; i < n; ++i) {
          for (j = i; j < n; ++j) {
               x = 0.;
               for (k = i; k < n; ++k) {
                    x += a[i][k]*b[k][j];
               }

               c[i][j] = x;
               if (flag)
                    c[j][i] = x;
          }
     }
}



void dmat_dtlgxms(double **a, double **b, double **c, int n, int flag) {

     int i;
     int j;
     int k;

     double x;

     for (i = 0; i < n; ++i) {
          for (j = i; j < n; ++j) {
               x = 0.;
               for (k = 0; k <= i; ++k) {
                    x += a[i][k]*b[k][j];
               }
               c[i][j] = x;
               if (flag)
                    c[j][i] = x;
          }
     }
}



/*******************************************************************************
 *
 ******************************************************************************/
void dmat_dtugxtlms(double **a, double **b,
                    double **c, double **d, double **w, int n) {

     int i;
     int j;
     int k;

     double x;

     for (i = 0; i < n; ++i) {
          for (j = 0; j <= i; ++j) {
               x = 0.;
               for (k = j; k < n; ++k) {
                    x += b[i][k]*c[k][j];
               }
               w[i][j] = x;
          }
     }

     for (i = 0; i < n; ++i) {
          for (j = 0; j <= i; ++j) {
               x = 0.;
               for (k = i; k < n; ++k) {
                    x += a[i][k]*w[k][j];
               }
               d[j][i] = x;
          }
     }
}



void dmat_dtugxtlms_l(double **a, double **b, double **c,
                      double **d, double **e, double **f,
                      double **g, double **h, double **w1, int n) {

     int i;
     int j;
     int k;

     double x;
     double y;
/*
     for (i = 0; i < n; ++i) {
          for (j = 0; j <= i; ++j) {
               x = 0.;
               for (k = j; k < n; ++k) {
                    x += f[i][k]*c[k][j];
               }
               w1[i][j]  = x;
          }
     }

     for (i = 0; i < n; ++i) {
          for (j = 0; j <= i; ++j) {
               x = 0.;
               for (k = j; k < n; ++k) {
                    x += b[i][k]*g[k][j];
               }
               w1[i][j] += x;
          }
     }
*/
     for (i = 0; i < n; ++i) {
          for (j = 0; j <= i; ++j) {
               x = 0.;
               for (k = j; k < n; ++k) {
                    x += f[i][k]*c[k][j];
                    x += b[i][k]*g[k][j];
               }
               w1[i][j]  = x;
          }
     }

     for (i = 0; i < n; ++i) {
          for (j = 0; j <= i; ++j) {
               x = 0.;
               for (k = i; k < n; ++k) {
                    x += a[i][k]*w1[k][j];
               }
               h[j][i] = x;
          }
     }

     for (i = 0; i < n; ++i) {
          for (j = 0; j <= i; ++j) {
               x = 0.;
               for (k = j; k < n; ++k) {
                    x += b[i][k]*c[k][j];
               }
               w1[i][j] = x;
          }
     }
/*
     for (i = 0; i < n; ++i) {
          for (j = 0; j <= i; ++j) {
               x = 0.;
               for (k = i; k < n; ++k) {
                    x += e[i][k]*w1[k][j];
               }
               h[j][i] += x;
          }
     }

     for (i = 0; i < n; ++i) {
          for (j = 0; j <= i; ++j) {
               x = 0.;
               for (k = i; k < n; ++k) {
                    x += a[i][k]*w1[k][j];
               }
               d[j][i]  = x;
          }
     }
*/
     for (i = 0; i < n; ++i) {
          for (j = 0; j <= i; ++j) {
               x = 0.;
               y = 0.;
               for (k = i; k < n; ++k) {
                    x += e[i][k]*w1[k][j];
                    y += a[i][k]*w1[k][j];
               }
               h[j][i] += x;
               d[j][i]  = y;
          }
     }
}



/*******************************************************************************
 *
 ******************************************************************************/
void dmat_chol(double **a, long n) {

     long i;

     long mm;
     long nn;

     long incy = 1;

     double alpha = -1.;
     double beta  =  1.;

     double da;

     for (i = 0; i < n; ++i) {
          mm = n - i;

          if (i > 0) {
               nn = i;

               dgemv_("t", &nn, &mm, &alpha, &a[i][0], &n,
                      &a[i][0], &incy, &beta, &a[i][i], &n);
          }

          da = 1. / sqrt(a[i][i]);

          dscal_(&mm, &da, &a[i][i], &n);
     }
}



void dmat_chol_l(double **a, double **b, long n) {

     long i;

     long mm;
     long nn;

     long incy = 1;

     double alpha = -1.;
     double beta  =  1.;

     double da;
     double db;

     for (i = 0; i < n; ++i) {
          mm = n - i;

          if (i > 0) {
               nn = i;

               dgemv_("t", &nn, &mm, &alpha, &b[i][0], &n,
                      &a[i][0], &incy, &beta, &b[i][i], &n);

               dgemv_("t", &nn, &mm, &alpha, &a[i][0], &n,
                      &b[i][0], &incy, &beta, &b[i][i], &n);


               dgemv_("t", &nn, &mm, &alpha, &a[i][0], &n,
                      &a[i][0], &incy, &beta, &a[i][i], &n);
          }

          da = 1. / sqrt(a[i][i]);
          db = da * -b[i][i] / (2. * a[i][i]);


          dscal_(&mm, &da, &b[i][i], &n);

          daxpy_(&mm, &db, &a[i][i], &n, &b[i][i], &n);

          dscal_(&mm, &da, &a[i][i], &n);
     }
}



/*******************************************************************************
 *
 ******************************************************************************/
void dmat_mul_brow(double **ax1, double **ax2, double **b11, double **b12,
                   double **b21, double **b22, long n,
                   double **cx1, double **cx2) {

     double alpha = 1.;
     double beta  = 0.;

     beta  = 0.;
     dgemm_("n", "n", &n, &n, &n, &alpha, *b11, &n, *ax1, &n, &beta, *cx1, &n);
     beta  = 1.;
     dgemm_("n", "n", &n, &n, &n, &alpha, *b21, &n, *ax2, &n, &beta, *cx1, &n);
     beta  = 0.;
     dgemm_("n", "n", &n, &n, &n, &alpha, *b12, &n, *ax1, &n, &beta, *cx2, &n);
     beta  = 1.;
     dgemm_("n", "n", &n, &n, &n, &alpha, *b22, &n, *ax2, &n, &beta, *cx2, &n);
}



/*******************************************************************************
 *
 ******************************************************************************/
void dmat_mul_bcol(double **a11, double **a12, double **a21, double **a22,
                   double **b1x, double **b2x, long n,
                   double **c1x, double **c2x) {

     double alpha = 1.;
     double beta  = 0.;

     beta  = 0.;
     dgemm_("n", "n", &n, &n, &n, &alpha, *b1x, &n, *a11, &n, &beta, *c1x, &n);
     beta  = 1.;
     dgemm_("n", "n", &n, &n, &n, &alpha, *b2x, &n, *a12, &n, &beta, *c1x, &n);
     beta  = 0.;
     dgemm_("n", "n", &n, &n, &n, &alpha, *b1x, &n, *a21, &n, &beta, *c2x, &n);
     beta  = 1.;
     dgemm_("n", "n", &n, &n, &n, &alpha, *b2x, &n, *a22, &n, &beta, *c2x, &n);
}



/*******************************************************************************
 *
 ******************************************************************************/
void dmat_copy_1x2(double **b11, double **b12,
                  double **a11, double **a12, int n) {

     int i;
     int j;

     for (i = 0; i < n; ++i) {
          for (j = 0; j < n; ++j) {
               b11[i][j] = a11[i][j];
               b12[i][j] = a12[i][j];
          }
     }
}



/*******************************************************************************
 *
 ******************************************************************************/
void dmat_copy_2x2(double **b11, double **b12, double **b21, double **b22,
                  double **a11, double **a12, double **a21, double **a22,
                  int n) {

     int i;
     int j;

     for (i = 0; i < n; ++i) {
          for (j = 0; j < n; ++j) {
               b11[i][j] = a11[i][j];
               b12[i][j] = a12[i][j];
               b21[i][j] = a21[i][j];
               b22[i][j] = a22[i][j];
          }
     }
}



/*******************************************************************************
 *
 ******************************************************************************/
void dmat_init_1x2(double **a11, double **a12, double alpha, int n) {

     int i;
     int j;

     for (i = 0; i < n; ++i) {
          for (j = 0; j < n; ++j) {
               a11[i][j] = 0.;
               a12[i][j] = 0.;
          }

          a11[i][i] = alpha;
     }
}



/*******************************************************************************
 *
 ******************************************************************************/
void dmat_init_2x2(double **a11, double **a12, double **a21, double **a22,
                  double alpha, int n) {

     int i;
     int j;

     for (i = 0; i < n; ++i) {
          for (j = 0; j < n; ++j) {
               a11[i][j] = 0.;
               a12[i][j] = 0.;
               a21[i][j] = 0.;
               a22[i][j] = 0.;
          }

          a11[i][i] = alpha;
          a22[i][i] = alpha;
     }
}



/*******************************************************************************
 *
 ******************************************************************************/
void dmat_scale_add_1x2(double **a11, double **a12, double **b11, double **b12,
                       double **P11, double **P12, double alpha, int n) {

     int i;
     int j;

     if (alpha ==  1.) {
          for (i = 0; i < n; ++i) {
               for (j = 0; j < n; ++j) {
                    P11[i][j] = a11[i][j] + b11[i][j];
                    P12[i][j] = a12[i][j] + b12[i][j];
               }
          }
     }
     else
     if (alpha == -1.) {
          for (i = 0; i < n; ++i) {
               for (j = 0; j < n; ++j) {
                    P11[i][j] = a11[i][j] - b11[i][j];
                    P12[i][j] = a12[i][j] - b12[i][j];
               }
          }
     }
     else {
          for (i = 0; i < n; ++i) {
               for (j = 0; j < n; ++j) {
                    P11[i][j] = a11[i][j] + b11[i][j] * alpha;
                    P12[i][j] = a12[i][j] + b12[i][j] * alpha;
               }
          }
     }
}



/*******************************************************************************
 *
 ******************************************************************************/
void dmat_scale_add_2x2(double **a11, double **a12, double **a21, double **a22,
                       double **b11, double **b12, double **b21, double **b22,
                       double **P11, double **P12, double **P21, double **P22,
                       double alpha, int n) {

     int i;
     int j;

     if (alpha ==  1.) {
          for (i = 0; i < n; ++i) {
               for (j = 0; j < n; ++j) {
                    P11[i][j] = a11[i][j] + b11[i][j];
                    P12[i][j] = a12[i][j] + b12[i][j];
                    P21[i][j] = a21[i][j] + b21[i][j];
                    P22[i][j] = a22[i][j] + b22[i][j];
               }
          }
     }
     else
     if (alpha == -1.) {
          for (i = 0; i < n; ++i) {
               for (j = 0; j < n; ++j) {
                    P11[i][j] = a11[i][j] - b11[i][j];
                    P12[i][j] = a12[i][j] - b12[i][j];
                    P21[i][j] = a21[i][j] - b21[i][j];
                    P22[i][j] = a22[i][j] - b22[i][j];
               }
          }
     }
     else {
          for (i = 0; i < n; ++i) {
               for (j = 0; j < n; ++j) {
                    P11[i][j] = a11[i][j] + b11[i][j] * alpha;
                    P12[i][j] = a12[i][j] + b12[i][j] * alpha;
                    P21[i][j] = a21[i][j] + b21[i][j] * alpha;
                    P22[i][j] = a22[i][j] + b22[i][j] * alpha;
               }
          }
     }
}



/*******************************************************************************
 *
 ******************************************************************************/
void dmat_mul_2x2(double **a11, double **a12, double **a21, double **a22,
                 double **b11, double **b12, double **b21, double **b22,
                 long n, double **P11, double **P12, double **P21, double **P22) {

     dmat_mul_brow(a11, a12, b11, b12, b21, b22, n, P11, P12);
     dmat_mul_brow(a21, a22, b11, b12, b21, b22, n, P21, P22);
}
