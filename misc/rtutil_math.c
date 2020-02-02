/*******************************************************************************
**
**    Copyright (C) 2007-2020 Greg McGarragh <greg.mcgarragh@colostate.edu>
**
**    This source code is licensed under the GNU General Public License (GPL),
**    Version 3.  See the file COPYING for more details.
**
*******************************************************************************/

#include <gutil.h>

#include "rtutil_math.h"


/*******************************************************************************
 *
 ******************************************************************************/
double factorial(int x) {

     return factorial_evaluate(1., 2, x);
}



double factorial_evaluate(double y, int x1, int x2) {

     int i;

     for (i = x1; i <= x2; ++i)
          y *= i;

     return y;
}



/*******************************************************************************
 *
 ******************************************************************************/
void extended_trapezoidal_quad(int n, double a, double b, double *x, double *w) {

     /* Mishchenko, 2006, page 413 */

     int i;

     double dx;

     if (n == 1) {
          x[0] = (a + b) / 2.;
          w[0] = 1.;
     }
     else {
          dx = (b - a) / (n - 1);

          x[0] = a;
          w[0] = dx / 2.;

          for (i = 1; i < n-1; ++i) {
               x[i] = a + (i - 1) * dx;
               w[i] = dx;
          }

          x[i] = b;
          w[i] = dx / 2.;
     }
}



/*******************************************************************************
 *
 ******************************************************************************/
static void gauss_leg_point(int n, int i, int ia, int ib,
                            double da, double *y1_, double *p2_) {

     /* Davis 1975. page 95 */

     int j;

     double db;

     double y2;
     double y1;

     double p3;
     double p2;
     double p1;

     double p1p;

     double p1pp;

/*
     y1 = cos(PI * (i + .75) / (n + .5));
*/
     db = (4. * i + 3.) / ia * PI;
     y1 = cos(db + da / tan(db));

     do {
          p2 = 0.;
          p1 = 1.;
          for (j = 0; j < n; ++j) {
               p3 = p2;
               p2 = p1;
/*
               p1 = ((2 * j + 1) * y1 * p2 - j * p3) / (j + 1);
*/
               db = y1 * p2;
               p1 = j / (j + 1.) * (db - p3) + db;
          }
/*
          p1p = n * (p2 - y1 * p1) / (1.0 - y1*y1);
          y2  = y1;
          y1  = y2 - p1 / p1p;
*/
          db = 1. - y1*y1;
          p1p = n * (p2 - y1 * p1) / db;

          p1pp = (2. * y1 * p1p - ib * p1) / db;

          y2 = y1;
          db = p1 / p1p;
          y1 = y2 - db * (1. + db * p1pp / (2. * p1p));

     } while (fabs(y1 - y2) > DBL_EPSILON);

     *y1_ = y1;
     *p2_ = p2;
}



void gauss_leg_quad(int n, double *x, double *w) {

     /* Davis 1975. page 95 */

     int i;
     int ii;

     int nn;

     int ia;
     int ib;

     double da;

     double dd;

     double y1;

     double p2;

     nn = 2 * n;

     ia = 4 * nn + 2;
     ib = nn * (nn + 1);

     da = (nn - 1) / (8. * pow((double) nn, 3));

     ii = 0;
     for (i = n; i < nn; ++i) {
          gauss_leg_point(nn, i, ia, ib, da, &y1, &p2);

          x[ii] = -y1;

          dd = nn * p2;
          w[ii] = 2. *  (1. - y1*y1) / (dd * dd);

          ii++;
     }
}



void gauss_leg_quad2(int n, double *x, double *w) {

     /* Davis 1975. page 95 */

     int i;
     int ii;

     int nn;

     int ia;
     int ib;

     double da;
     double db;
     double dc;

     double y1;

     double p2;

     nn = (n + 1) / 2;

     ia = 4 * n + 2;
     ib = n * (n + 1);

     da = (n - 1) / (8. * pow((double) n, 3));

     ii = n - 1;
     for (i = 0; i < nn; ++i) {
          gauss_leg_point(n, i, ia, ib, da, &y1, &p2);

          db = .5 * y1;
          x[ i] = .5 - db;
          x[ii] = .5 + db;
/*
          db = n * p2;
          w[ i] = (1. - y1*y1) / (db * db);
*/
          dc = n * p2;
          w[ i] = (1. - 2. * db*y1) / (dc * dc);

          w[ii] = w[ i];

          ii--;
     }
}



void gauss_leg_quadx(int n, double x1, double x2, double *x, double *w) {

     gauss_leg_quadx_l(n, 0, x1, x2, NULL, NULL, x, w, NULL, NULL);
}



void gauss_leg_quadx_l(int n, int n_derivs,
                       double x1, double x2, double *x1_l, double *x2_l,
                       double *x, double *w, double **x_l, double **w_l) {

     /* Davis 1975. page 95 */

     int i;
     int ii;
     int j;

     int nn;

     int ia;
     int ib;

     double da;
     double db;
     double *db_l;
     double dc;
     double *dc_l;
     double dd;
     double dd_l;
     double de;

     double y1;

     double p2;

     if (n_derivs > 0) {
          db_l = alloc_array1_d(n_derivs);
          dc_l = alloc_array1_d(n_derivs);
     }

     nn = (n + 1) / 2;

     ia = 4 * n + 2;
     ib = n * (n + 1);

     da = (n - 1) / (8. * pow((double) n, 3));

     db = (x1 + x2) / 2.;
     dc = (x2 - x1) / 2.;

     for (i = 0; i < n_derivs; ++i) {
          db_l[i] = (x1_l[i] + x2_l[i]) / 2.;
          dc_l[i] = (x2_l[i] - x1_l[i]) / 2.;
     }

     ii = n - 1;
     for (i = 0; i < nn; ++i) {
          gauss_leg_point(n, i, ia, ib, da, &y1, &p2);

          dd = dc * y1;
          x[ i] = db - dd;
          x[ii] = db + dd;
/*
          dd = n * p2;
          dd = 2. / (dd * dd);
          w[ i] = dc * (1. - y1*y1) * dd;
*/
          de = n * p2;
          de = 2. / (de * de);
          w[ i] = (dc - dd*y1) * de;

          w[ii] = w[ i];

          for (j = 0; j < n_derivs; ++j) {
               dd_l = dc_l[j] * y1;
               x_l[ i][j] = db_l[j] - dd_l;
               x_l[ii][j] = db_l[j] + dd_l;
/*
               w_l[ i][j] = dc_l[j] * (1. - y1*y1) * dd;
*/
               w_l[ i][j] = (dc_l[j] - dd_l*y1) * de;

               w_l[ii][j] = w_l[i][j];
          }

          ii--;
     }

     if (n_derivs > 0) {
          free_array1_d(db_l);
          free_array1_d(dc_l);
     }
}



void gauss_leg_quad_fit(int n, double x1, double x2,
                        double *xa, double *wa, double *xb, double *wb) {

     gauss_leg_quad_fit_l(n, 0, x1, x2, NULL, NULL,
                          xa, wa, xb, wb, NULL, NULL, NULL, NULL);
}



void gauss_leg_quad_fit_l(int n, int n_derivs,
                          double x1, double x2, double *x1_l, double *x2_l,
                          double *xa, double *wa, double *xb, double *wb,
                          double **xa_l, double **wa_l, double **xb_l, double **wb_l) {

     int i;
     int ii;
     int j;

     int nn;

     double db;
     double *db_l;
     double dc;
     double *dc_l;
     double dd;
     double dd_l;

     if (n_derivs > 0) {
          db_l = alloc_array1_d(n_derivs);
          dc_l = alloc_array1_d(n_derivs);
     }

     nn = (n + 1) / 2;

     db = (x1 + x2) / 2.;
     dc = (x2 - x1) / 2.;

     for (i = 0; i < n_derivs; ++i) {
          db_l[i] = (x1_l[i] + x2_l[i]) / 2.;
          dc_l[i] = (x2_l[i] - x1_l[i]) / 2.;
     }

     ii = n - 1;
     for (i = 0; i < nn; ++i) {
          dd = dc * xa[ii];
          xb[ i] = db - dd;
          xb[ii] = db + dd;

          wb[ i] = wa[ i] * dc;
          wb[ii] = wb[ i];

          for (j = 0; j < n_derivs; ++j) {
               dd_l = dc_l[j] * xa[ii] + dc * xa_l[ii][j];
               xb_l[ i][j] = db_l[j] - dd_l;
               xb_l[ii][j] = db_l[j] + dd_l;

               wb_l[ i][j] = wa_l[i][j] * dc + wa[i] * dc_l[j];
               wb_l[ii][j] = wb_l[i][j];
          }

          ii--;
     }

     if (n_derivs > 0) {
          free_array1_d(db_l);
          free_array1_d(dc_l);
     }
}



/*******************************************************************************
 *
 ******************************************************************************/
int simpsons_rule_quad(int n, double a, double b, double *x, double *w) {

     int i;

     double f;

     if (n % 2 == 0) {
          fprintf(stderr, "ERROR: number of points for Simpson's rule must be odd\n");
          return -1;
     }

     f = (b - a) / (n - 1);

     for (i = 0; i < n; ++i)
          x[i] = a + i * f;

     w[0    ] = f / 3.;
     w[n - 1] = f / 3.;

     for (i = 1; i < n - 1; i += 2)
          w[i] = 4. * f / 3.;

     for (i = 2; i < n - 2; i += 2)
          w[i] = 2. * f / 3.;

     return 0;
}



/*******************************************************************************
 *
 ******************************************************************************/
void trapezoidal_quad(int n, double a, double b, double *x, double *w) {

     int i;

     double f;

     f = (b - a) / (n - 1);

     for (i = 0; i < n; ++i)
          x[i] = a + i * f;

     w[0    ] = .5 * f;
     w[n - 1] = .5 * f;

     for (i = 1; i < n - 1; ++i)
          w[i] = f;
}



/*******************************************************************************
 *
 ******************************************************************************/
void leg_poly(int n_l, double mu, double *p) {

     leg_poly2(1, n_l, &mu, &p);
/*
     double cn;
     double dn;

     leg_poly_assoc(0, n_l, mu, &cn, &dn, p, 1);
*/
/*
     int i;

     dcomplex *gsf;

     gsf = alloc_array1_dc(n_l);

     gen_spher_funcs(0, 0, n_l, mu, gsf);

     for (i = 0; i < n_l; ++i)
          p[i] = creal(gsf[i]);

     free_array1_dc(gsf);
*/
}



void leg_poly2(int n_mu, int n_l, double *mu, double **p) {

     /* Davis 1975. page 95 */

     int i;
     int j;
/*
     int jp1;
     int twojp1;
*/
     double a;
     double b;

     for (i = 0; i < n_mu; ++i) {
          p[i][0] = 1.;

          if (n_l > 1)
               p[i][1] = mu[i];
     }
/*
     for (j = 1; j < n_l - 1; ++j) {
          jp1    = j + 1;
          twojp1 = 2 * j + 1;

          for (i = 0; i < n_mu; ++i) {
               p[i][j+1] = (twojp1 * mu[i] * p[i][j] - j * p[i][j-1]) / jp1;
          }
     }
*/
     for (j = 1; j < n_l - 1; ++j) {
          a = j / (j + 1.);

          for (i = 0; i < n_mu; ++i) {
               b = mu[i] * p[i][j];
               p[i][j+1] = a * (b - p[i][j-1]) + b;
          }
     }
/*
     double cn;
     double dn;

     leg_poly_assoc2(0, n_mu, n_l, mu, &cn, &dn, p, 1);
*/
/*
     int i;
     int j;

     dcomplex **gsf;

     gsf = alloc_array2_dc(n_mu, n_l);

     gen_spher_funcs2(0, 0, n_mu, n_l, mu, gsf);

     for (i = 0; i < n_mu; ++i) {
          for (j = 0; j < n_l; ++j) {
               p[i][j] = creal(gsf[i][j]);
          }
     }

     free_array2_dc(gsf);
*/
}



/*******************************************************************************
 *
 ******************************************************************************/
void leg_poly_assoc(int m, int n_l, double mu,
                    double *cn, double *dn, double *y, int flag) {

     leg_poly_assoc2(m, 1, n_l, &mu, cn, dn, &y, flag);
}



void leg_poly_assoc2(int m, int n_mu, int n_l, double *mu,
                     double *cn, double *dn, double **y, int flag) {

     /* Dave 1970 */

     int i;
     int j;

     int twojp1;

     double a;
     double b;

     if (flag) {
          if (m == 0) {
               *cn = 1.;
               *dn = 1.;
          }
          else {
               a = 2. * m;
               *cn *= -sqrt((a - 1.) / a);
               *dn *= -sqrt((a + 1.) / a);
          }
     }

     for (i = 0; i < n_mu; ++i) {
          if (m == 0) {
               y[i][0] = 1.;

               if (n_l > 1) {
                    y[i][1] = mu[i];
               }
          }
          else {
               for (j = 0; j < m; ++j) {
                    y[i][j] = 0.;
               }

               a = pow((1. - mu[i]*mu[i]), m / 2.);

               y[i][m+0] = *cn * a;

               if (n_l > m + 1) {
                    y[i][m+1] = *dn * mu[i] * a;
               }
          }
     }

     for (j = m + 1; j < n_l - 1; ++j) {
          twojp1 = 2 * j + 1;
          b = sqrt((double) ((j - m + 1) * (j + m + 1)));
          a = sqrt((double) ((j + m + 0) * (j - m + 0)));
          for (i = 0; i < n_mu; ++i) {
               y[i][j+1] = (twojp1 * mu[i] * y[i][j] - a * y[i][j-1]) / b;
          }
     }
}



/*******************************************************************************
 *
 ******************************************************************************/
void gen_spher_funcs(int p, int q, int n_l, double mu, dcomplex *gsf) {

     gen_spher_funcs2(p, q, 1, n_l, &mu, &gsf);
}



void gen_spher_funcs2(int p, int q, int n_mu, int n_l, double *mu, dcomplex **gsf) {

     /* Mischenko 1991 */

     int i;
     int j;

     int r;

     int p_2;
     int q_2;
     int pq;

     int apmq;
     int appq;
/*
     int j_2;
     int jp1;
     int jjp1;
     int jp1_2;
     int twojp1;
*/
     double j_2;
     double jp1;
     double jjp1;
     double jp1_2;
     double twojp1;

     double da;
     double db;

     double apmqd2;
     double appqd2;

     r = MAX(abs(p), abs(q));

     p_2 = p * p;
     q_2 = q * q;
     pq  = p * q;

     apmq = abs(p - q);
     appq = abs(p + q);

     apmqd2 = apmq / 2.;
     appqd2 = appq / 2.;

     da = creal(cpow(-_Complex_I, apmq) / pow(2., r) *
                sqrt(factorial(2 * r) / (factorial(apmq) * factorial(appq))));

     for (i = 0; i < n_mu; ++i) {
          if (r == 0) {
               gsf[i][0] = 1.;

               if (n_l > 1) {
                    gsf[i][1] = mu[i];
               }

               j = 1;
          }
          else {
               for (j = 0; j < r; ++j) {
                    gsf[i][j] = 0.;
               }

               gsf[i][r] = da * pow(1. - mu[i], apmqd2) *
                                pow(1. + mu[i], appqd2);

               j = r;
          }
     }

     for ( ; j < n_l - 1; ++j) {
          j_2    = j * j;
          jp1    = j + 1;
          jjp1   = j * jp1;
          jp1_2  = jp1 * jp1;
          twojp1 = 2 * j + 1;

          da = jp1 * sqrt(j_2 - p_2) * sqrt(j_2 - q_2);
          db = j * sqrt(jp1_2 - p_2) * sqrt(jp1_2 - q_2);

          for (i = 0; i < n_mu; ++i) {
               gsf[i][j+1] = (twojp1 * (jjp1 * mu[i] - pq) * gsf[i][j] -
                              da * gsf[i][j-1]) / db;
          }
     }
}



void gen_spher_funcs_0q(int q, int n_l, double mu, dcomplex *gsf) {

     gen_spher_funcs_0q_2(q, 1, n_l, &mu, &gsf);
}



void gen_spher_funcs_0q_2(int q, int n_mu, int n_l, double *mu, dcomplex **gsf) {

     int i;
     int j;

     int q_2;

     int j_2;
     int jp1;
     int jp1_2;
     int twojp1;

     double da;
     double db;

     double qd2;

     double fq;

     qd2 = q / 2.;

     q_2 = q * q;

     fq  = factorial(q);

     da = creal(cpow(-_Complex_I, q) / pow(2., q) * sqrt(factorial(2 * q) / (fq * fq)));

     for (i = 0; i < n_mu; ++i) {
          if (q == 0) {
               gsf[i][0] = 1.;

               if (n_l > 1) {
                    gsf[i][1] = mu[i];
               }

               j = 1;
          }
          else {
               for (j = 0; j < q; ++j) {
                    gsf[i][j] = 0.;
               }

               gsf[i][q] = da * pow(1. - mu[i], qd2) *
                                pow(1. + mu[i], qd2);

               j = q;
          }
     }

     for ( ; j < n_l - 1; ++j) {
          j_2    = j * j;
          jp1    = j + 1;
          jp1_2  = jp1 * jp1;
          twojp1 = 2 * j + 1;

          da = sqrt((double) (j_2 - q_2));
          db = sqrt((double) (jp1_2 - q_2));

          for (i = 0; i < n_mu; ++i) {
               gsf[i][j+1] = (twojp1 * mu[i] * gsf[i][j] - da * gsf[i][j-1]) / db;
          }
     }
}



void gen_spher_funcs_p0(int p, int n_l, double mu, dcomplex *gsf) {

     gen_spher_funcs_p0_2(p, 1, n_l, &mu, &gsf);
}



void gen_spher_funcs_p0_2(int p, int n_mu, int n_l, double *mu, dcomplex **gsf) {

     int i;
     int j;

     int p_2;

     int j_2;
     int jp1;
     int jp1_2;
     int twojp1;

     double da;
     double db;

     double pd2;

     double fp;

     pd2 = p / 2.;

     p_2 = p * p;

     fp  = factorial(p);

     da = creal(cpow(-_Complex_I, p) / pow(2., p) * sqrt(factorial(2 * p) / (fp * fp)));

     for (i = 0; i < n_mu; ++i) {
          if (p == 0) {
               gsf[i][0] = 1.;

               if (n_l > 1) {
                    gsf[i][1] = mu[i];
               }

               j = 1;
          }
          else {
               for (j = 0; j < p; ++j) {
                    gsf[i][j] = 0.;
               }

               gsf[i][p] = da * pow(1. - mu[i], pd2) *
                                pow(1. + mu[i], pd2);

               j = p;
          }
     }

     for ( ; j < n_l - 1; ++j) {
          j_2    = j * j;
          jp1    = j + 1;
          jp1_2  = jp1 * jp1;
          twojp1 = 2 * j + 1;

          da = sqrt((double) (j_2 - p_2));
          db = sqrt((double) (jp1_2 - p_2));

          for (i = 0; i < n_mu; ++i) {
               gsf[i][j+1] = (twojp1 * mu[i] * gsf[i][j] - da * gsf[i][j-1]) / db;
          }
     }
}



void gen_spher_funcs_prt(int n_l, double mu,
                         double *p00, double *p0p2, double *p2p2, double *p2m2) {

     gen_spher_funcs_prt2(1, n_l, &mu, &p00, &p0p2, &p2p2, &p2m2);
}



void gen_spher_funcs_prt2(int n_mu, int n_l, double *mu,
                          double **p00, double **p0p2, double **p2p2, double **p2m2) {
#ifdef CRAP
     int i;
     int j;
/*
     double cn;
     double dn;
*/
     dcomplex **gsf;

     gsf = alloc_array2_dc(n_mu, n_l);
/*
     leg_poly2(n_mu, n_l, mu, p00);
*/
/*
     leg_poly_assoc2(0, n_mu, n_l, mu, &cn, &dn, p00, 1);
*/
/*
     gen_spher_funcs_0q_2(0, n_mu, n_l, mu, gsf);
     for (i = 0; i < n_mu; ++i) {
          for (j = 0; j < n_l; ++j) {
               p00[i][j] = creal(gsf[i][j]);
          }
     }
*/
     gen_spher_funcs2(0,  0, n_mu, n_l, mu, gsf);
     for (i = 0; i < n_mu; ++i) {
          for (j = 0; j < n_l; ++j) {
               p00[i][j] = creal(gsf[i][j]);
          }
     }
/*
     gen_spher_funcs_0q_2(2, n_mu, n_l, mu, gsf);
     for (i = 0; i < n_mu; ++i) {
          for (j = 0; j < n_l; ++j) {
               p0p2[i][j] = creal(gsf[i][j]);
          }
     }
*/
     gen_spher_funcs2(0,  2, n_mu, n_l, mu, gsf);
     for (i = 0; i < n_mu; ++i) {
          for (j = 0; j < n_l; ++j) {
               p0p2[i][j] = creal(gsf[i][j]);
          }
     }

     gen_spher_funcs2(2,  2, n_mu, n_l, mu, gsf);
     for (i = 0; i < n_mu; ++i) {
          for (j = 0; j < n_l; ++j) {
               p2p2[i][j] = creal(gsf[i][j]);
          }
     }

     gen_spher_funcs2(2, -2, n_mu, n_l, mu, gsf);
     for (i = 0; i < n_mu; ++i) {
          for (j = 0; j < n_l; ++j) {
               p2m2[i][j] = creal(gsf[i][j]);
          }
     }

     free_array2_dc(gsf);
#endif
     int i;
     int j;
/*
     int j_2;
     int jp1;
     int jjp1;
     int jp1_2;
     int twojp1;
*/
     double j_2;
     double jp1;
     double jjp1;
     double jp1_2;
     double twojp1;

     double da;

     double mu2;

     double da_0p2;
     double db_0p2;

     double da_2x2;
     double db_2x2;
     double dc_2x2;

     da_0p2 = .25 * sqrt(6.);

     for (i = 0; i < n_mu; ++i) {
          mu2 = mu[i]*mu[i];

          if (n_l > 0) {
               p00[i][0] = 1.;
               p0p2[i][0] = 0.;
               p2p2[i][0] = 0.;
               p2m2[i][0] = 0.;
          }

          if (n_l > 1) {
               p00[i][1] = mu[i];
               p0p2[i][1] = 0.;
               p2p2[i][1] = 0.;
               p2m2[i][1] = 0.;
          }

          if (n_l > 2) {
               p00[i][2] = (3. * mu2 - 1.) / 2.;
               p0p2[i][2] = da_0p2 * (mu2 - 1.);

               da = (1. + mu[i]);
               p2p2[i][2] = 0.25 * da*da;

               da = (1. - mu[i]);
               p2m2[i][2] = 0.25 * da*da;
          }
     }

     for (j = 2; j < n_l - 1; ++j) {
          j_2    = j * j;
          jp1    = j + 1;
          jjp1   = j * jp1;
          jp1_2  = jp1 * jp1;
          twojp1 = 2 * j + 1;

          da_0p2 = sqrt(j_2 - 4);
          db_0p2 = sqrt(jp1_2 - 4);

          da_2x2 = twojp1 * 4;
          db_2x2 = jp1 * (j_2 - 4);
          dc_2x2 = j * (jp1_2 - 4);

          for (i = 0; i < n_mu; ++i) {
               da = twojp1 * mu[i];
               p00 [i][j+1] = (da * p00[i][j] -
                               j * p00[i][j-1]) / jp1;

               p0p2[i][j+1] = (da * p0p2[i][j] -
                               da_0p2 * p0p2[i][j-1]) / db_0p2;

               da = twojp1 * jjp1 * mu[i];
               p2p2[i][j+1] = ((da - da_2x2) * p2p2[i][j] -
                               db_2x2 * p2p2[i][j-1]) / dc_2x2;

               p2m2[i][j+1] = ((da + da_2x2) * p2m2[i][j] -
                               db_2x2 * p2m2[i][j-1]) / dc_2x2;
          }
     }
}
