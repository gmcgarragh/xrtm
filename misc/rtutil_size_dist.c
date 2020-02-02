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
#include "rtutil_size_dist.h"
#include "zeroin.h"


/*******************************************************************************
 *
 ******************************************************************************/
int get_particle_dist_n(enum size_dist_type dist_type, int n_int1, int n_int2, int n_quad) {

     if (dist_type == SIZE_DIST_MONO)
          return 1;
     else
     if (dist_type == SIZE_DIST_MODIFIED_GAMMA)
          return n_int2 * n_quad;
     else
     if (dist_type == SIZE_DIST_LOG_NORMAL)
          return n_int2 * n_quad;
     else
     if (dist_type == SIZE_DIST_POWER_LAW)
          return n_int2 * n_quad;
     else
     if (dist_type == SIZE_DIST_GAMMA)
          return n_int2 * n_quad;
     else
     if (dist_type == SIZE_DIST_MODIFIED_POWER_LAW)
          return n_int1 * n_quad + n_int2 * n_quad;
     else
     if (dist_type == SIZE_DIST_MODIFIED_BIMODAL_LOG_NORMAL)
          return n_int2 * n_quad;
     else
     if (dist_type == SIZE_DIST_EXPONENTIAL)
          return n_int2 * n_quad;
     else {
          fprintf(stderr, "ERROR: Invalid distribution type: %d\n", dist_type);
          return -1;
     }
}



/*******************************************************************************
 *
 ******************************************************************************/
/*
static double zeroin_func(double *r1, double *data) {

     double r2;

     r2 = (1. + data[1]) * 2. * data[0] - *r1;

     return (r2 - *r1) / log(r2 / *r1) - data[0];
}



#ifdef __cplusplus
extern "C" {
#endif
     double zeroin_(double *ax, double *bx, double (*f)(double *, double *), double *data, double *tol);
#ifdef __cplusplus
}
#endif

void get_power_law_range(double a1, double b1, double *r1, double *r2) {

     double a;
     double b;

     double data[2];

     data[0] = a1;
     data[1] = b1;

     a = DBL_EPSILON;
     b = 0.;
     *r1 = zeroin_(&a, &a1, zeroin_func, data, &b);

     *r2 = (1. + b1) * 2. * a1 - *r1;
}
*/


static double zeroin_func(double r1, double *data) {

     double r2;

     r2 = (1. + data[1]) * 2. * data[0] - r1;

     return (r2 - r1) / log(r2 / r1) - data[0];
}



void get_power_law_range(double a1, double a2, double *r1, double *r2) {

     double a;
     double b;

     double data[2];

     data[0] = a1;
     data[1] = a2;

     a = DBL_EPSILON;
     b = 0.;
     *r1 = zeroin(a, a1, zeroin_func, data, b);

     *r2 = (1. + a2) * 2. * a1 - *r1;
}



/*******************************************************************************
 *
 ******************************************************************************/
static void get_particle_dist_quad(int n_int, int n_quad, int n_derivs,
                         double r1, double r2, double *r1_l, double *r2_l,
                         double *qx, double *qw, double **qx_l, double **qw_l) {

     int i;
     int ii;
     int j;

     double a;
     double b;

     double *a_l;
     double *b_l;
     double *c_l;

     double *tx;
     double *tw;
     double **tx_l;
     double **tw_l;

     tx = alloc_array1_d(n_quad);
     tw = alloc_array1_d(n_quad);

     if (n_derivs > 0) {
          a_l = alloc_array1_d(n_derivs);
          b_l = alloc_array1_d(n_derivs);
          c_l = alloc_array1_d(n_derivs);

          tx_l = alloc_array2_d(n_quad, n_derivs);
          tw_l = alloc_array2_d(n_quad, n_derivs);
     }

     for (i = 0; i < n_derivs; ++i) {
          a_l[i] = 0.;
          b_l[i] = 0.;
     }

     gauss_leg_quadx_l(n_quad, n_derivs, -1., 1., a_l, b_l, tx, tw, tx_l, tw_l);

     a = (r2 - r1) / n_int;
     for (i = 0; i < n_derivs; ++i)
          a_l[i] = (r2_l[i] - (r1_l ? r1_l[i] : 0.)) / n_int;

     for (i = 0; i < n_int; ++i) {
          ii = i * n_quad;
          b = r1 + i * a;
          for (j = 0; j < n_derivs; ++j) {
               b_l[j] = (r1_l ? r1_l[j] : 0.) + i * a_l[j];
               c_l[j] = b_l[j] + a_l[j];
          }
          gauss_leg_quad_fit_l(n_quad, n_derivs, b, b + a, b_l, c_l, tx, tw,
                               qx+ii, qw+ii, tx_l, tw_l, qx_l+ii, qw_l+ii);
     }

     free_array1_d(tx);
     free_array1_d(tw);

     if (n_derivs > 0) {
          free_array1_d(a_l);
          free_array1_d(b_l);
          free_array1_d(c_l);

          free_array2_d(tx_l);
          free_array2_d(tw_l);
     }
}



/*******************************************************************************
 *
 ******************************************************************************/
int get_particle_dist(enum size_dist_type dist_type,
                      int n_int1, int n_int2, int n_quad, int n_derivs,
                      double a1, double a2, double a3, double a4,
                      double a5, double r1, double r2,
                      double *a1_l, double *a2_l, double *a3_l, double *a4_l,
                      double *a5_l, double *r1_l, double *r2_l,
                      double *qx, double *qw, double *nr,
                      double *r1_, double *r2_,
                      double **qx_l, double **qw_l,double **nr_l,
                      double *r1_l_, double *r2_l_,
                      double *norm, double *norm_l) {

     int i;
     int j;

     int n_qsiz;
     int n_qsiz1;
     int n_qsiz2;

     double a;
     double b;
     double c;
     double d;
     double e;
     double f;
     double f2;
     double g;
     double h;
     double o;
     double p;

     double b_1;
     double c_1;
     double d_1;
     double e_1;
     double h_1;

     double b_2;
     double c_2;
     double d_2;
     double e_2;
     double h_2;
     double o_1;
     double o_2;

     double *a_l;
     double *b_l;
     double *c_l;

     double *b_1_l;
     double *c_1_l;

     double *b_2_l;
     double *c_2_l;

     if (n_derivs > 0) {
          a_l   = alloc_array1_d(n_derivs);
          b_l   = alloc_array1_d(n_derivs);
          c_l   = alloc_array1_d(n_derivs);

          b_1_l = alloc_array1_d(n_derivs);
          c_1_l = alloc_array1_d(n_derivs);

          b_2_l = alloc_array1_d(n_derivs);
          c_2_l = alloc_array1_d(n_derivs);
     }

     if (dist_type == SIZE_DIST_MONO) {
          n_qsiz = 1;

          *r1_ = a1;
          *r2_ = a1;
          for (i = 0; i < n_derivs; ++i) {
               r1_l_[i] = 0.;
               r2_l_[i] = 0.;
          }

          qx[0] = a1;
          qw[0] = 1.;

          for (i = 0; i < n_derivs; ++i) {
               qx_l[0][i] = a1_l[i];
               qw_l[0][i] = 0.;
          }

          a = 1.;

          for (i = 0; i < n_derivs; ++i)
               a_l[i] = 0.;

          nr[0] = 1.;

          for (i = 0; i < n_derivs; ++i)
               nr_l[0][i] = 0.;
     }
     else
     if (dist_type == SIZE_DIST_MODIFIED_GAMMA) {
          n_qsiz = n_int2 * n_quad;

          *r1_ = r1;
          *r2_ = r2;
          for (i = 0; i < n_derivs; ++i) {
               r1_l_[i] = r1_l[i];
               r2_l_[i] = r2_l[i];
          }

          get_particle_dist_quad(n_int2, n_quad, n_derivs, *r1_, *r2_,
                                 r1_l_, r2_l_, qx, qw, qx_l, qw_l);

          a = 0.;

          b = pow(a2, a3);
          c = a3 * b;
          d = b / c;
          e = a3 / a2;
          f = log(a2);
          g = a1 / c;

          for (i = 0; i < n_derivs; ++i) {
               a_l[i] = 0.;

               c_l[i] = a3_l[i] * d + (a2_l[i] * e + a3_l[i] * f);
          }

          for (i = 0; i < n_qsiz; ++i) {
               h = pow(qx[i], a3);

               nr[i] = pow(qx[i], a1) * exp(-g * h);

               a += nr[i] * qw[i];

               if (n_derivs > 0)
                    o = log(qx[i]);

               for (j = 0; j < n_derivs; ++j) {
                    p = qx_l[i][j] / qx[i];

                    nr_l[i][j] = nr[i] * ((p * a1 + a1_l[j] * o) + h * ((-a1_l[j] +
                                 a1 * (-(p * a3 + a3_l[j] * o) + c_l[j])) / c));
                    a_l[j] += nr_l[i][j] * qw[i] + nr[i] * qw_l[i][j];
               }
          }
     }
     else
     if (dist_type == SIZE_DIST_LOG_NORMAL) {
          n_qsiz = n_int2 * n_quad;

          *r1_ = r1;
          *r2_ = r2;
          for (i = 0; i < n_derivs; ++i) {
               r1_l_[i] = r1_l[i];
               r2_l_[i] = r2_l[i];
          }

          get_particle_dist_quad(n_int2, n_quad, n_derivs, *r1_, *r2_,
                                 r1_l_, r2_l_, qx, qw, qx_l, qw_l);

          a = 0.;

          b = log(a1);
          c = 2. * a2;

          for (i = 0; i < n_derivs; ++i) {
               a_l[i] = 0.;

               b_l[i] = a1_l[i] / a1;
               c_l[i] = 2. * a2_l[i];
          }

          for (i = 0; i < n_qsiz; ++i) {
               d = log(qx[i]) - b;
               e = d * d / c;

               nr[i] = exp(-e) / qx[i];
               a += nr[i] * qw[i];

               for (j = 0; j < n_derivs; ++j) {
                    f = qx_l[i][j] / qx[i];
                    nr_l[i][j] = nr[i] * (-f - ((2. * d *
                                 (f - b_l[j]) - e * c_l[j]) / c));
                    a_l[j] += nr_l[i][j] * qw[i] + nr[i] * qw_l[i][j];
               }
          }
     }
     else
     if (dist_type == SIZE_DIST_POWER_LAW) {
          n_qsiz = n_int2 * n_quad;

          get_power_law_range(a1, a2, r1_, r2_);

          for (i = 0; i < n_derivs; ++i) {
               /* from the system of the 2 equations in the function f */
               a = log(*r2_) - log(*r1_);
               b = -a + (*r2_ - *r1_) / *r1_;
               c =  a - (*r2_ - *r1_) / *r2_;
               a = a1_l[i] * a * a;

               d = a2_l[i] * 2. * a1 + (1. + a2) * 2. * a1_l[i];

               r1_l_[i] = (c * d - a) / (c - b);
               r2_l_[i] = (b * d - a) / (b - c);
          }

          get_particle_dist_quad(n_int2, n_quad, n_derivs, *r1_, *r2_,
                                 r1_l_, r2_l_, qx, qw, qx_l, qw_l);

          a = 0.;

          for (i = 0; i < n_derivs; ++i)
               a_l[i] = 0.;

          for (i = 0; i < n_qsiz; ++i) {
               nr[i] = 1. / pow(qx[i], 3);
               a += nr[i] * qw[i];

               for (j = 0; j < n_derivs; ++j) {
                    nr_l[i][j] = -nr[i] * 3. * qx_l[i][j] / qx[i];
                    a_l[j] += nr_l[i][j] * qw[i] + nr[i] * qw_l[i][j];
               }
          }
     }
     else
     if (dist_type == SIZE_DIST_GAMMA) {
          n_qsiz = n_int2 * n_quad;

          *r1_ = r1;
          *r2_ = r2;
          for (i = 0; i < n_derivs; ++i) {
               r1_l_[i] = r1_l[i];
               r2_l_[i] = r2_l[i];
          }

          get_particle_dist_quad(n_int2, n_quad, n_derivs, *r1_, *r2_,
                                 r1_l_, r2_l_, qx, qw, qx_l, qw_l);

          a = 0.;

          b = a1 * a2;
          c = (1. - 3. * a2) / a2;

          for (i = 0; i < n_derivs; ++i) {
               a_l[i] = 0.;

               b_l[i] =  a1_l[i] * a2 + a1 * a2_l[i];
               c_l[i] = (-3. * a2_l[i] - c * a2_l[i]) / a2;
          }

          for (i = 0; i < n_qsiz; ++i) {
               d = pow(qx[i],  c);
               e = qx[i] / b;
               f = exp(-e);

               nr[i] = d * f;
               a += nr[i] * qw[i];

               if (n_derivs > 0) {
                    g = c / qx[i];
                    h = log(qx[i]);
               }

               for (j = 0; j < n_derivs; ++j) {
                    nr_l[i][j] = d * ((qx_l[i][j] * g + c_l[j] * h) +
                                 (-qx_l[i][j] + e * b_l[j]) / b) * f;
                    a_l[j] += nr_l[i][j] * qw[i] + nr[i] * qw_l[i][j];
               }
          }
     }
     else
     if (dist_type == SIZE_DIST_MODIFIED_POWER_LAW) {
          n_qsiz1 = n_int1 * n_quad;
          n_qsiz2 = n_int2 * n_quad;

          n_qsiz  = n_qsiz1 + n_qsiz2;

          *r1_ = r1;
          *r2_ = r2;
          for (i = 0; i < n_derivs; ++i) {
               r1_l_[i] = r1_l[i];
               r2_l_[i] = r2_l[i];
          }

          get_particle_dist_quad(n_int1, n_quad, n_derivs, 0, *r1_,
                                 NULL,  r1_l_, qx, qw, qx_l, qw_l);
          get_particle_dist_quad(n_int2, n_quad, n_derivs, *r1_, *r2_,
                                 r1_l_, r2_l_, qx   + n_qsiz1, qw   + n_qsiz1,
                                               qx_l + n_qsiz1, qw_l + n_qsiz1);

          a = 0.;

          for (i = 0; i < n_derivs; ++i)
               a_l[i] = 0.;

          for (i = 0; i < n_qsiz1; ++i) {
               nr[i] = 1.;
               a += nr[i] * qw[i];

               for (j = 0; j < n_derivs; ++j) {
                    nr_l[i][j] = 0.;
                    a_l[j] += nr_l[i][j] * qw[i] + nr[i] * qw_l[i][j];
               }
          }

          for (     ; i < n_qsiz ; ++i) {
               b = qx[i] / *r1_;

               nr[i] = pow(b, a1);
               a += nr[i] * qw[i];

               for (j = 0; j < n_derivs; ++j) {
                    nr_l[i][j] = nr[i] * ((qx_l[i][j] - b * r1_l_[j]) / *r1_ *
                                 a1 / b + a1_l[j] * log(b));
                    a_l[j] += nr_l[i][j] * qw[i] + nr[i] * qw_l[i][j];
               }
          }
     }
     else
     if (dist_type == SIZE_DIST_MODIFIED_BIMODAL_LOG_NORMAL) {
          n_qsiz = n_int2 * n_quad;

          *r1_ = r1;
          *r2_ = r2;
          for (i = 0; i < n_derivs; ++i) {
               r1_l_[i] = r1_l[i];
               r2_l_[i] = r2_l[i];
          }

          get_particle_dist_quad(n_int2, n_quad, n_derivs, *r1_, *r2_,
                                 r1_l_, r2_l_, qx, qw, qx_l, qw_l);

          a = 0.;

          b_1 = log(a1);
          c_1 = 2. * a2;

          b_2 = log(a3);
          c_2 = 2. * a4;

          for (i = 0; i < n_derivs; ++i) {
               a_l[i] = 0.;

               b_1_l[i] = a1_l[i] / a1;
               c_1_l[i] = 2. * a2_l[i];

               b_2_l[i] = a3_l[i] / a3;
               c_2_l[i] = 2. * a4_l[i];
          }

          for (i = 0; i < n_qsiz; ++i) {
               d_1 = log(qx[i]) - b_1;
               e_1 = d_1 * d_1 / c_1;

               d_2 = log(qx[i]) - b_2;
               e_2 = d_2 * d_2 / c_2;

               o_1 = exp(-e_1);
               o_2 = exp(-e_2);

               f  = qx[i] * qx[i];
               f2 = f * f;

               nr[i] = (o_1 + a5 * o_2) / f2;
               a += nr[i] * qw[i];

               d_1 *= 2.;
               d_2 *= 2.;

               o_1 /= c_1;
               o_2 /= c_2;

               p = nr[i] * 4. * qx[i] * f;

               for (j = 0; j < n_derivs; ++j) {
                    g = qx_l[i][j] / qx[i];

                    h_1 = ((-g + b_1_l[j]) * d_1 + e_1 * c_1_l[j]) * o_1;
                    h_2 = ((-g + b_2_l[j]) * d_2 + e_2 * c_2_l[j]) * o_2;

                    nr_l[i][j] = (h_1 + a5_l[j] * exp(-e_2) + a5 * h_2 -
                                  p * qx_l[i][j]) / f2;

                    a_l[j] += nr_l[i][j] * qw[i] + nr[i] * qw_l[i][j];
               }
          }
     }
     else
     if (dist_type == SIZE_DIST_EXPONENTIAL) {
          n_qsiz = n_int2 * n_quad;

          *r1_ = r1;
          *r2_ = r2;
          for (i = 0; i < n_derivs; ++i) {
               r1_l_[i] = r1_l[i];
               r2_l_[i] = r2_l[i];
          }

          a = 0.;

          get_particle_dist_quad(n_int2, n_quad, n_derivs, *r1_, *r2_,
                                 r1_l_, r2_l_, qx, qw, qx_l, qw_l);

          for (i = 0; i < n_qsiz; ++i) {
               /* n(D) = exp(-lambda * D); */

               nr[i] = 2. * exp(-a1 * 2. * qx[i]);

               a += nr[i] * qw[i];

               for (j = 0; j < n_derivs; ++j) {
/*
                    nr_l[i][j] =

                    a_l[j] += nr_l[i][j] * qw[i] + nr[i] * qw_l[i][j];
*/
               }
          }
     }
     else {
          fprintf(stderr, "ERROR: Invalid distribution type: %d\n", dist_type);
          return -1;
     }

     *norm = a;
     for (i = 0; i < n_derivs; ++i)
          norm_l[i] = a_l[i];

     for (i = 0; i < n_qsiz; ++i) {
          nr[i] /= a;

          for (j = 0; j < n_derivs; ++j) {
               nr_l[i][j] = (nr_l[i][j] - nr[i] * a_l[j]) / a;
          }
     }

     if (n_derivs > 0) {
          free_array1_d(a_l);
          free_array1_d(b_l);
          free_array1_d(c_l);

          free_array1_d(b_1_l);
          free_array1_d(c_1_l);

          free_array1_d(b_2_l);
          free_array1_d(c_2_l);
     }

     return 0;
}



/*******************************************************************************
 *
 ******************************************************************************/
int get_dist_parameters(int n_quad, int n_derivs,
                        double *qx,     double *qw,     double *nr,
                        double **qx_l,  double **qw_l,  double **nr_l,
                        double *reff,   double *veff,   double *gavg,
                        double *vavg,   double *ravg,   double *rvw,
                        double *reff_l, double *veff_l, double *gavg_l,
                        double *vavg_l, double *ravg_l, double *rvw_l) {

     int i;
     int j;

     double a;
     double b;
     double c;
     double d;
     double e;
     double f;
     double g;
     double h;
     double p;
     double q;

     double r[4];

     double *a_l;
     double *b_l;
     double *c_l;
     double *d_l;

     if (n_derivs > 0) {
          a_l = alloc_array1_d(n_derivs);
          b_l = alloc_array1_d(n_derivs);
          c_l = alloc_array1_d(n_derivs);
          d_l = alloc_array1_d(n_derivs);
     }

     a = 0.;
     b = 0.;
     c = 0.;
     d = 0.;

     for (i = 0; i < n_derivs; ++i) {
          a_l[i] = 0.;
          b_l[i] = 0.;
          c_l[i] = 0.;
          d_l[i] = 0.;
     }

     for (i = 0; i < n_quad; ++i) {
          e = nr[i] * qx[i] * qw[i];
          a += e;
          e *= qx[i];
          b += e;
          e *= qx[i];
          c += e;
          e *= qx[i];
          d += e;

          r[0] = qx[i];
          r[1] = r[0] * qx[i];
          r[2] = r[1] * qx[i];

          e = qx[i] * qw[i];
          f = nr[i] * qw[i];
          g = nr[i] * qx[i];

          for (j = 0; j < n_derivs; ++j) {
               h       = e * nr_l[i][j];
               p       = f * qx_l[i][j];
               q       = g * qw_l[i][j];
               a_l[j] +=  h +      p + q;
               b_l[j] += (h + 2. * p + q) * r[0];
               c_l[j] += (h + 3. * p + q) * r[1];
               d_l[j] += (h + 4. * p + q) * r[2];
          }
     }

     *gavg = b * PI;
     *vavg = 4. / 3. * PI * c;
     *ravg = a;
     *rvw  = 4. / 3. * PI * d / *vavg;

     *reff = c / b;

     e = *vavg * *vavg;
     f = b * b;

     for (i = 0; i < n_derivs; ++i) {
          gavg_l[i] = b_l[i] * PI;
          vavg_l[i] = 4. / 3. * PI * c_l[i];
          ravg_l[i] = a_l[i];
          rvw_l [i] = 4. / 3. * PI * (d_l[i] / *vavg - d * vavg_l[i] / e);

          reff_l[i] = c_l[i] / b - c * b_l[i] / f;
     }

     a = 0.;

     for (i = 0; i < n_derivs; ++i)
          a_l[i] = 0.;

     for (i = 0; i < n_quad; ++i) {
          b = (qx[i] - *reff);
          c = b * b;
          d = qx[i] * qx[i];
          e = d * qw[i];
          f = c * e;

          a += nr[i] * f;

          g = 2. * b * e;
          h = 2. * qx[i] * qw[i];

          for (j = 0; j < n_derivs; ++j) {
               p = (qx_l[i][j] - reff_l[j]);

               a_l[j] += nr_l[i][j] * f + nr[i] *
                         (g * p + c * (h * qx_l[i][j] + d * qw_l[i][j]));
          }
     }

     b = *gavg * *reff * *reff;

     *veff = a * PI / b;

     c = PI / b;
     d = *veff * *reff / b;
     e = d * *reff;
     f = d * *gavg * 2.;

     for (i = 0; i < n_derivs; ++i)
          veff_l[i] = a_l[i] * c - (gavg_l[i] * e + f * reff_l[i]);

     if (n_derivs > 0) {
          free_array1_d(a_l);
          free_array1_d(b_l);
          free_array1_d(c_l);
          free_array1_d(d_l);
     }

     return 0;
}
