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

#include <rtutil_math.h>
#include <rtutil_support.h>

#include "xrtm.h"
#include "xrtm_adding.h"
#include "xrtm_utility.h"
#include "xrtm_utility_a.h"


/*******************************************************************************
 *
 ******************************************************************************/
void delta_m_coef_a(int n_coef, double f, double *f_a,
                    double **coef, double **coef_a, double **coef_prime_a,
                    int allsix, int coef_type) {

     int j;

     double a;
     double b;
     double c;

     double t;

     a = 1. - f;

     for (j = 0; j < n_coef; ++j) {
          b = 2 * j + 1;
          c = b * f;

          t = coef_prime_a[0][j] / a;
          coef_a[0][j] += t;
          *f_a += t * (-b + (coef[0][j] - c) / a);

          if (allsix) {
               t = coef_prime_a[1][j] / a;
               coef_a[1][j] += t;
               *f_a += t * (-b + (coef[1][j] - c) / a);

               t = coef_prime_a[2][j] / a;
               coef_a[2][j] += t;
               *f_a += t * (-b + (coef[2][j] - c) / a);

               t = coef_prime_a[3][j] / a;
               coef_a[3][j] += t;
               *f_a += t * (-b + (coef[3][j] - c) / a);

               if (coef_type == SCAT_COEF_GC) {
                    t = coef_prime_a[4][j] / a;
                    coef_a[4][j] += t;
                    *f_a += t * coef[4][j] / a;

                    t = coef_prime_a[5][j] / a;
                    coef_a[5][j] += t;
                    *f_a += t * coef[5][j] / a;
               }
               else {
                    coef_a[4][j] += coef_prime_a[4][j];
                    coef_a[5][j] += coef_prime_a[5][j];
               }
          }
     }

     for (j = 0; j < n_coef; ++j) {
          coef_prime_a[0][j] = 0.;

          if (allsix) {
               coef_prime_a[1][j] = 0.;
               coef_prime_a[2][j] = 0.;
               coef_prime_a[3][j] = 0.;
               coef_prime_a[4][j] = 0.;
               coef_prime_a[5][j] = 0.;
          }
     }
}



void delta_m_omega_a(double f, double *f_a,
                     double omega, double *omega_a, double *omega_prime_a) {

     double a;
     double b;

     double t;

     a = 1. - f;
     b = 1. - omega * f;

     t = *omega_prime_a / b;
     *omega_a += t * a;
     *f_a     -= t * omega;

     t *= omega * a / b;
     *omega_a += t * f;
     *f_a     += t * omega;

     *omega_prime_a = 0.;
}



void delta_m_ltau_a(double f, double *f_a,
                    double omega, double *omega_a,
                    double ltau, double *ltau_a, double *ltau_prime_a) {

     double a;

     double t;

     a = 1. - omega * f;

     t = - *ltau_prime_a * ltau;

     *omega_a += t * f;
     *f_a     += t * omega;
     *ltau_a  += *ltau_prime_a * a;

     *ltau_prime_a = 0.;
}



void delta_m_a(int n_coef, double f, double *f_a,
               double **coef, double omega, double ltau,
               double **coef_a, double **coef_prime_a,
               double *omega_a, double *omega_prime_a,
               double *ltau_a, double *ltau_prime_a,
               int allsix, int coef_type) {

     delta_m_ltau_a(f, f_a, omega, omega_a,
                    ltau, ltau_a, ltau_prime_a);

     delta_m_omega_a(f, f_a,
                     omega, omega_a, omega_prime_a);

     delta_m_coef_a(n_coef, f, f_a,
                    coef, coef_a, coef_prime_a, allsix, coef_type);
}



/*******************************************************************************
 *
 ******************************************************************************/
void n_t_tms_scaling_a(double f, double *f_a, double omega0, double omega_tms, double *omega0_a, double omega_tms_a) {

     double t;

     t = omega_tms_a / (1. - f * omega0);

     *omega0_a += t;

     *f_a += t * omega_tms * omega0;

     *omega0_a += t * omega_tms * f;
}



/*******************************************************************************
 *
 ******************************************************************************/
static void build_local_r_or_t_mat_a(int i_four, int n_quad,
                                     double *v1, double *qw_v, double a, 
                                     double omega, double *omega_a,
                                     double **P_m, double **P_m_a,
                                     double **r, double **r_a,
                                     work_data work) {

     int i;
     int j;

     double b;

     double **w1;

     w1 = get_work1(&work, WORK_DXX);

     b = a * (1. + (i_four == 0 ? 1. : 0.)) / 4.;

     dmat_diag_mul(v1, r_a, w1, n_quad, n_quad);
     dmat_mul_diag(w1, qw_v, w1, n_quad, n_quad);
     dmat_scale(b, w1, w1, n_quad, n_quad);

     for (i = 0; i < n_quad; ++i) {
          for (j = 0; j < n_quad; ++j) {
               *omega_a += w1[i][j] * P_m[i][j];
          }
     }

     dmat_scale(omega, w1, w1, n_quad, n_quad);
     dmat_add(P_m_a, w1, P_m_a, n_quad, n_quad);

     dmat_zero(r_a, n_quad, n_quad);
}



void build_local_r_and_t_a(int i_four, int n_quad,
                           double *qx_v, double *qw_v, double omega, double *omega_a,
                           double **P_p, double **P_m, double **r, double **t,
                           double **P_p_a, double **P_m_a, double **r_a, double **t_a,
                           work_data work) {

     double *v1;

     v1 = get_work1(&work, WORK_DX);

     dvec_inv(qx_v, v1, n_quad);

     build_local_r_or_t_mat_a(i_four, n_quad, v1, qw_v,  1.,
                              omega, omega_a, P_p, P_p_a, t, t_a, work);

     build_local_r_or_t_mat_a(i_four, n_quad, v1, qw_v, -1.,
                              omega, omega_a, P_m, P_m_a, r, r_a, work);
}



/*******************************************************************************
 *
 ******************************************************************************/
void build_txr_a(int n_quad, int n_stokes,
                 double **r_p, double **t_p, double **tpr, double **tmr,
                 double **r_p_a, double **t_p_a, double **tpr_a, double **tmr_a,
                 work_data work) {

     int n_quad_v;

     double **w1;

     n_quad_v = n_quad * n_stokes;

     w1 = get_work1(&work, WORK_DXX);

     dmat_add(tpr_a, tmr_a, w1, n_quad_v, n_quad_v);
     dmat_add(t_p_a, w1, t_p_a, n_quad_v, n_quad_v);

     dmat_sub(tpr_a, tmr_a, w1, n_quad_v, n_quad_v);
     dmat_add(r_p_a, w1, r_p_a, n_quad_v, n_quad_v);

     if (n_stokes > 2)
          dmat_mul_D_A(n_quad, n_stokes, r_p_a, r_p_a);

     dmat_zero(tpr_a, n_quad_v, n_quad_v);
     dmat_zero(tmr_a, n_quad_v, n_quad_v);
}



/*******************************************************************************
 *
 ******************************************************************************/
void build_gamma_a(int n_quad,
                   double **tpr, double **tmr, double **gamma,
                   double **tpr_a, double **tmr_a, double **gamma_a,
                   work_data work) {
/*
     double **w1;
     double **w2;

     w1 = get_work1(&work, WORK_DXX);
     w2 = get_work1(&work, WORK_DXX);

     dmat_trans(tpr, w1, n_quad, n_quad);
     dmat_mul(gamma_a, w1, n_quad, n_quad, n_quad, w2);
     dmat_add(tmr_a, w2, tmr_a, n_quad, n_quad);
*/
     dmat_gxgxmx(0, gamma_a, 1, tpr, 1., tmr_a, 1., n_quad, n_quad, n_quad);
/*
     dmat_trans(tmr, w1, n_quad, n_quad);
     dmat_mul(w1, gamma_a, n_quad, n_quad, n_quad, w2);
     dmat_add(tpr_a, w2, tpr_a, n_quad, n_quad);
*/
     dmat_gxgxmx(1, tmr, 0, gamma_a, 1., tpr_a, 1., n_quad, n_quad, n_quad);

     dmat_zero(gamma_a, n_quad, n_quad);
}



/*******************************************************************************
 *
 ******************************************************************************/
#ifdef INCLUDE_DEV_SOURCE
#include "../src2/xrtm_utility_a2.c"
#endif
