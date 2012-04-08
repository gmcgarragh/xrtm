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
void calc_quad_and_umus(int n_quad, int n_umus, int n_stokes, int *n_quad_x_, int *n_quad_d_, int *n_quad_v, int *n_quad_v_x, int *n_quad_v_d, int *n_umus_v, int flag) {

     int n_quad_x;
     int n_quad_d;

     if (! flag)
          n_quad_x = n_quad + n_umus;
     else
          n_quad_x = n_quad;

     n_quad_d = n_quad + n_umus;

     if (n_quad_x_)  *n_quad_x_  = n_quad_x;
     if (n_quad_d_)  *n_quad_d_  = n_quad_d;

     if (n_quad_v)   *n_quad_v   = n_quad   * n_stokes;
     if (n_quad_v_x) *n_quad_v_x = n_quad_x * n_stokes;
     if (n_quad_v_d) *n_quad_v_d = n_quad_d * n_stokes;

     if (n_umus_v)   *n_umus_v   = n_umus   * n_stokes;
}



/*******************************************************************************
 *
 ******************************************************************************/
typedef struct {
     uchar or_;
     uchar and_;
     int count;
} flags_meta_data;


uchar *flags_alloc(int n) {

     uchar *a;

     a = (uchar *) malloc(sizeof(flags_meta_data) +
         n * sizeof(uchar)) + sizeof(flags_meta_data);
     if (a == NULL) {
          eprintf("ERROR: memory allocation failed\n");
          return NULL;
     }

     return a;
}



void flags_free(uchar *a) {

     free((uchar *) a - sizeof(flags_meta_data));
}



uchar **flags_alloc2(int m, int n) {

     int i;

     uchar **a;

     a = (uchar **)
         ((uchar *) malloc(sizeof(flags_meta_data) +
         m * sizeof(uchar *)) + sizeof(flags_meta_data));
     if (a == NULL) {
          eprintf("ERROR: memory allocation failed\n");
          return NULL;
     }

     a[0] = (uchar *) malloc(m * (sizeof(flags_meta_data) +
            n * sizeof(uchar))) + sizeof(flags_meta_data);
     if (a[0] == NULL) {
          eprintf("ERROR: memory allocation failed\n");
          return NULL;
     }

     for (i = 1; i < m; ++i)
          a[i] = a[0] + i * (sizeof(flags_meta_data) + n);

     return a;
}



void flags_free2(uchar **a) {

     free((uchar *) a[0] - sizeof(flags_meta_data));
     free((uchar *) a    - sizeof(flags_meta_data));
}



int flags_fill_meta(uchar *flags, int n) {

     int i;

     flags_meta_data *meta;

     meta = (flags_meta_data *) ((uchar *) flags - sizeof(flags_meta_data));

     meta->or_   = 0;
     meta->and_  = 0;
     meta->count = 0;
     for (i = 0; i < n; ++i) {
          if (flags[i]) {
               meta->or_  = 1;
               meta->count++;
          }
          else
               meta->and_ = 0;
     }

     return meta->count;
}



int flags_fill_meta2(uchar **flags, int m, int n) {

     int i;

     flags_meta_data *meta;

     for (i = 0; i < m; ++i)
          flags_fill_meta(flags[i], n);

     meta = (flags_meta_data *) ((uchar *) flags - sizeof(flags_meta_data));

     meta->or_   = 0;
     meta->and_  = 0;
     meta->count = 0;
     for (i = 0; i < m; ++i) {
          if (flags_or(flags[i], n)) {
               meta->or_  = 1;
               meta->count += flags_count(flags[i], n);
          }
          else
               meta->and_ = 0;
     }

     return meta->count;
}


int flags_count(void *flags, int n) {

     if (n == 0)
          return 0;

     return ((flags_meta_data *) ((uchar *) flags - sizeof(flags_meta_data)))->count;
}



int flags_or(void *flags, int n) {

     if (n == 0)
          return 0;

     return ((flags_meta_data *) ((uchar *) flags - sizeof(flags_meta_data)))->or_;
}



int flags_and(void *flags, int n) {

     if (n == 0)
          return 0;

     return ((flags_meta_data *) ((uchar *) flags - sizeof(flags_meta_data)))->and_;
}



int flags_count2(uchar **flags, int m, int n) {

     if (m == 0 || n == 0)
          return 0;

     return ((flags_meta_data *) ((uchar *) flags - sizeof(flags_meta_data)))->count;
}



int flags_or2(uchar **flags, int m, int n) {

     if (m == 0 || n == 0)
          return 0;

     return ((flags_meta_data *) ((uchar *) flags - sizeof(flags_meta_data)))->or_;
}



int flags_and2(uchar **flags, int m, int n) {

     if (m == 0 || n == 0)
          return 0;

     return ((flags_meta_data *) ((uchar *) flags - sizeof(flags_meta_data)))->and_;
}



/*******************************************************************************
 *
 ******************************************************************************/
void derivs_merge_s(int n_derivs, uchar *derivs_s1, uchar *derivs_sm, uchar *derivs_s) {

     int i;

     for (i = 0; i < n_derivs; ++i) {
          if (! derivs_s1[i] && ! derivs_sm[i])
               derivs_s[i] = ADDING_U_U;
          else
          if (! derivs_s1[i] &&   derivs_sm[i])
               derivs_s[i] = ADDING_U_P;
          else
          if (  derivs_s1[i] &&   derivs_sm[i])
               derivs_s[i] = ADDING_P_P;
#ifdef DEBUG
          else {
               eprintf("ERROR: invalid layer adding combination: %d %d\n",
                       derivs_s1[i], derivs_sm[i]);
               exit(1);
          }
#endif
     }

     flags_fill_meta(derivs_s, n_derivs);
}



uchar derivs_union_d(int n_derivs, uchar *derivs_d1) {

     int i;

     uchar derivs_d2;

     derivs_d2 = ADDING_U_U;

     for (i = 0; i < n_derivs; ++i)
          derivs_d2 |= derivs_d1[i];

     return derivs_d2;
}



void derivs_union_h2(int n_layers, int n_derivs, uchar **derivs_h1, uchar *derivs_h2) {

     int i;

     for (i = 0; i < n_layers; ++i)
          derivs_h2[i] = flags_or(derivs_h1[i], n_derivs);
}



void derivs_union_d2(int n_layers, int n_derivs, uchar **derivs_d1, uchar *derivs_d2) {

     int i;

     for (i = 0; i < n_layers; ++i)
          derivs_d2[i] = derivs_union_d(n_derivs, derivs_d1[i]);
}



/*******************************************************************************
 *
 ******************************************************************************/
double singularity_adjust_up(double base, double value, double epsilon) {

     if (base - value < epsilon)
          return base + epsilon;
     else
          return value;
}



double singularity_adjust_down(double base, double value, double epsilon) {

     if (base - value < epsilon)
          return base - epsilon;
     else
          return value;
}



double singularity_adjust_auto(double base, double value, double epsilon) {

     double a;

     a = base - value;

     if (a <  0. && a > -epsilon)
          return base + epsilon;
     else
     if (a >= 0. && a <  epsilon)
          return base - epsilon;
     else
          return value;
}



/*******************************************************************************
 *
 ******************************************************************************/
int check_phase_vecs_norm(int n_quad, int n_stokes,
                          double *qx, double *qw, double *P_p, double *P_m) {

     int i;
     int ii;

     double a;

     a = 0.;

     ii = 0;
     for (i = 0; i < n_quad; ++i) {
          a += (P_p[ii] + P_m[ii]) * qw[i];
          ii += n_stokes;
     }

     if (fabs(2. - a) > 10.e-8) {
          eprintf("ERROR: phase vectors not normalized: %.16f != 2.0\n", a);
          return 1;
     }

     return 0;
}



int check_phase_mats_norm(int n_quad1, int n_quad2, int n_stokes,
                          double *qx, double *qw, double **P_p, double **P_m) {

     int i;
     int ii;
     int j;
     int jj;

     double a;

     ii = 0;
     for (i = 0; i < n_quad1; ++i) {
          a = 0.;

          jj = 0;
          for (j = 0; j < n_quad2; ++j) {
               a += (P_p[ii][jj] + P_m[ii][jj]) * qw[j];
               jj += n_stokes;
          }

          if (fabs(2. - a) > 10.e-8) {
               eprintf("ERROR: phase matrices not normalized: %.16f != 2.0\n", a);
               return 1;
          }

          ii += n_stokes;
     }

     return 0;
}



int check_phase_vecs_norm_l(int n_quad, int n_stokes,
                            double *qx, double *qw, double *P_p_l, double *P_m_l) {

     int i;
     int ii;

     double a;

     a = 0.;

     ii = 0;
     for (i = 0; i < n_quad; ++i) {
          a += (P_p_l[ii] + P_m_l[ii]) * qw[i];
          ii += n_stokes;
     }

     if (fabs(0. - a) > 10.e-8) {
          eprintf("ERROR: integration of the columns of the derivatives of "
                 "the phase matrices not equal to zero: %.16f != 0.0\n", a);
          return 1;
     }

     return 0;
}



int check_phase_mats_norm_l(int n_quad1, int n_quad2, int n_stokes,
                            double *qx, double *qw, double **P_p_l, double **P_m_l) {

     int i;
     int ii;
     int j;
     int jj;

     double a;

     ii = 0;
     for (i = 0; i < n_quad1; ++i) {
          a = 0.;

          jj = 0;
          for (j = 0; j < n_quad2; ++j) {
               a += (P_p_l[ii][jj] + P_m_l[ii][jj]) * qw[j];
               jj += n_stokes;
          }

          if (fabs(0. - a) > 10.e-8) {
               eprintf("ERROR: integration of the columns of the derivatives of "
                      "the phase matrices not equal to zero: %.16f != 0.0\n", a);
               return 1;
          }

          ii += n_stokes;
     }

     return 0;
}



int check_R_and_T_norm(int n_quad, int n_stokes, double **R, double **T) {

     int i;
     int ii;
     int j;
     int jj;

     double a;

     ii = 0;
     for (i = 0; i < n_quad; ++i) {
          a = 0.;

          jj = 0;
          for (j = 0; j < n_quad; ++j) {
               a += R[ii][jj] + T[ii][jj];
               jj += n_stokes;
          }

          if (fabs(1. - a) > 10.e-8) {
               eprintf("ERROR: R and T not normalized: %.16f != 1.0\n", a);
               return 1;
          }

          ii += n_stokes;
     }

     return 0;
}



int check_R_and_T_norm_l(int n_quad, int n_stokes, double **R_l, double **T_l) {

     int i;
     int ii;
     int j;
     int jj;

     double a;

     ii = 0;
     for (i = 0; i < n_quad; ++i) {
          a = 0.;

          jj = 0;
          for (j = 0; j < n_quad; ++j) {
               a += R_l[ii][jj] + T_l[ii][jj];
               jj += n_stokes;
          }

          if (fabs(0. - a) > 10.e-8) {
               eprintf("ERROR: summation of the columns of the derivatives of "
                      "the R and T matrics not equal to zero: %.16f != 0.0\n", a);
               return 1;
          }

          ii += n_stokes;
     }

     return 0;
}



int check_conserve_energy(int n_quad, int n_stokes, double *qx,
                          double *qw, double F_0, double mu_0,
                          double btran, double *I_p, double *I_m) {

     int i;
     int ii;

     double a;

     a = 0.;

     ii = 0;
     for (i = 0; i < n_quad; ++i) {
          a += (I_p[ii] + I_m[ii]) * qx[i] * qw[i];
          ii += n_stokes;
     }

     a = (2. * PI * a + F_0 * mu_0 * btran) / mu_0;

     if (fabs(1. - a) > 10.e-8) {
          eprintf("ERROR: model does not conserve energy: %.16f != 1.0\n", a);
          return 1;
     }

     return 0;
}



int check_conserve_energy_l(int n_quad, int n_stokes, double *qx,
                            double *qw, double F_0, double mu_0,
                            double btran_l, double *I_p_l, double *I_m_l) {

     int i;
     int ii;

     double a;

     a = 0.;

     ii = 0;
     for (i = 0; i < n_quad; ++i) {
          a += (I_p_l[ii] + I_m_l[ii]) * qx[i] * qw[i];
          ii += n_stokes;
     }

     a = (2. * PI * a + F_0 * mu_0 * btran_l) / mu_0;

     if (fabs(0. - a) > 10.e-8) {
          eprintf("ERROR: model does not conserve energy: %.16f != 1.0\n", a);
          return 1;
     }

     return 0;
}



/*******************************************************************************
 *
 ******************************************************************************/
void delta_m_g(int n_derivs, double f, double *f_l,
               double g, double *g_prime,
               double *g_l, double *g_prime_l, int derivs) {

     double a;

     a = 1. - f;

     *g_prime = (g - f) / a;

     if (derivs)
          delta_m_g_l(n_derivs, f, f_l, g, g_l, g_prime_l);
}



void delta_m_g_l(int n_derivs, double f, double *f_l,
                 double g, double *g_l, double *g_prime_l) {

     int i;

     double a;
     double b;

     a = 1. - f;
     b = g  - f;

     for (i = 0; i < n_derivs; ++i)
          g_prime_l[i] = (g_l[i] - f_l[i] + b * f_l[i] / a) / a;
}



void delta_m_coef(int n_coef, int n_derivs, double f, double *f_l,
                  double **coef, double **coef_prime,
                  double ***coef_l, double ***coef_prime_l,
                  int allsix, int coef_type, int derivs) {

     int i;

     double a;
     double b;

     a = 1. - f;

     for (i = 0; i < n_coef; ++i) {
          b = (2 * i + 1) * f;
          coef_prime[0][i] = (coef[0][i] - b) / a;
          if (allsix) {
               coef_prime[1][i] = (coef[1][i] - b) / a;
               coef_prime[2][i] = (coef[2][i] - b) / a;
               coef_prime[3][i] = (coef[3][i] - b) / a;

               if (coef_type == SCAT_COEF_GC) {
                    coef_prime[4][i] = coef[4][i] / a;
                    coef_prime[5][i] = coef[5][i] / a;
               }
               else {
                    coef_prime[4][i] = coef[4][i];
                    coef_prime[5][i] = coef[5][i];
               }
          }
     }

     if (derivs)
          delta_m_coef_l(n_coef, n_derivs, f, f_l,
                         coef, coef_l, coef_prime_l, allsix, coef_type);
}



void delta_m_coef_l(int n_coef, int n_derivs, double f, double *f_l,
                    double **coef, double ***coef_l, double ***coef_prime_l,
                    int allsix, int coef_type) {

     int i;
     int j;

     double a;
     double b;
     double c;
     double d;
     double e;

     a = 1. - f;

     for (i = 0; i < n_derivs; ++i) {
          for (j = 0; j < n_coef; ++j) {
               b = 2 * j + 1;
               c = b * f;
               d = b * f_l[i];
               e = f_l[i] / a;
               coef_prime_l[i][0][j] =
                    ((coef_l[i][0][j] - d) + (coef[0][j] - c) * e) / a;
               if (allsix) {
                    coef_prime_l[i][1][j] =
                         ((coef_l[i][1][j] - d) + (coef[1][j] - c) * e) / a;
                    coef_prime_l[i][2][j] =
                         ((coef_l[i][2][j] - d) + (coef[2][j] - c) * e) / a;
                    coef_prime_l[i][3][j] =
                         ((coef_l[i][3][j] - d) + (coef[3][j] - c) * e) / a;
                    if (coef_type == SCAT_COEF_GC) {
                         coef_prime_l[i][4][j] =
                              (coef_l[i][4][j] + coef[4][j] * e) / a;
                         coef_prime_l[i][5][j] =
                              (coef_l[i][5][j] + coef[5][j] * e) / a;
                    }
                    else {
                         coef_prime_l[i][4][j] = coef_l[i][4][j];
                         coef_prime_l[i][5][j] = coef_l[i][5][j];
                    }
               }
          }
     }
}



void delta_m_omega(int n_derivs, double f, double *f_l,
                   double omega, double *omega_prime,
                   double *omega_l, double *omega_prime_l, int derivs) {

     double a;
     double b;

     a = 1. - f;
     b = 1. - omega * f;

     *omega_prime = a / b * omega;

     if (derivs)
          delta_m_omega_l(n_derivs, f, f_l, omega, omega_l, omega_prime_l);
}



void delta_m_omega_l(int n_derivs, double f, double *f_l,
                     double omega, double *omega_l, double *omega_prime_l) {

     int i;

     double a;
     double b;
     double c;

     a = 1. - f;
     b = 1. - omega * f;

     for (i = 0; i < n_derivs; ++i) {
          c = omega * f_l[i];

          omega_prime_l[i] = (omega_l[i] * a - c) / b +
                              omega * a * (omega_l[i] * f + c) / (b * b);
     }
}



void delta_m_ltau(int n_derivs, double f, double *f_l,
                  double omega, double *omega_l,
                  double ltau, double *ltau_prime,
                  double *ltau_l, double *ltau_prime_l, int derivs) {

     double a;

     a = 1. - omega * f;

     *ltau_prime = a * ltau;

     if (derivs)
          delta_m_ltau_l(n_derivs, f, f_l,
                         omega, omega_l, ltau, ltau_l, ltau_prime_l);
}



void delta_m_ltau_l(int n_derivs, double f, double *f_l,
                    double omega, double *omega_l,
                    double ltau, double *ltau_l, double *ltau_prime_l) {

     int i;

     double a;

     a = 1. - omega * f;

     for (i = 0; i < n_derivs; ++i)
          ltau_prime_l[i] = -(omega_l[i] * f + omega * f_l[i]) * ltau +
                             a * ltau_l[i];
}



void delta_m(int n_coef, int n_derivs, double f, double *f_l,
             double **coef, double **coef_prime,
             double omega, double *omega_prime,          
             double ltau, double *ltau_prime,
             double ***coef_l, double ***coef_prime_l,
             double *omega_l, double *omega_prime_l,
             double *ltau_l, double *ltau_prime_l,
             int allsix, int coef_type, int derivs) {

     delta_m_coef(n_coef, n_derivs, f, f_l,
                  coef, coef_prime,
                  coef_l, coef_prime_l, allsix, coef_type, derivs);

     delta_m_omega(n_derivs, f, f_l,
                   omega, omega_prime,
                   omega_l, omega_prime_l, derivs);

     delta_m_ltau(n_derivs, f, f_l, omega, omega_l,
                  ltau, ltau_prime, ltau_l, ltau_prime_l, derivs);
}



void delta_m_l(int n_coef, int n_derivs, double f, double *f_l,
               double **coef, double omega, double ltau,
               double ***coef_l, double ***coef_prime_l,
               double *omega_l, double *omega_prime_l,
               double *ltau_l, double *ltau_prime_l,
               int allsix, int coef_type) {

     delta_m_coef_l(n_coef, n_derivs, f, f_l,
                   coef, coef_l, coef_prime_l, allsix, coef_type);

     delta_m_omega_l(n_derivs, f, f_l,
                     omega, omega_l, omega_prime_l);

     delta_m_ltau_l(n_derivs, f, f_l, omega, omega_l,
                    ltau, ltau_l, ltau_prime_l);
}



/*******************************************************************************
 *
 ******************************************************************************/
void n_t_tms_scaling(int n_derivs, double f, double *f_l, double omega0, double *omega_tms, double *omega0_l, double *omega_tms_l) {

     int i;

     double a;

     a = 1. - f * omega0;

     *omega_tms = omega0 / a;

     for (i = 0; i < n_derivs; ++i)
          omega_tms_l[i] = (omega0_l[i] - *omega_tms * (-f_l[i] * omega0 - f * omega0_l[i])) / a;
}



/*******************************************************************************
 *
 ******************************************************************************/
void chapman_functions(int n_layers, double mu_0, double z0, double *z, double **chap) {

     int i;
     int j;

     double a;
     double b;

     double r1;
     double r2;

     a = sin(acos(mu_0));

     a *= a;

     for (i = 0; i < n_layers; ++i) {
          b = z0 + z[i+1];
          b = b * b * a;
          for (j = 0; j <= i; ++j) {
               r1 = z0 + z[j];
               r2 = z0 + z[j+1];
               chap[i][j] = (sqrt(r1*r1 - b) - sqrt(r2*r2 - b)) / (r1 - r2);
          }
     }
}



/*******************************************************************************
 *
 ******************************************************************************/
static void build_local_r_or_t_mat(int i_four, int n_quad, int n_derivs,
                                   double *v1, double *qw_v, double a, 
                                   double omega, double *omega_l,
                                   double **P_m, double ***P_m_l,
                                   double **r, double ***r_l, uchar *derivs,
                                   work_data work) {
if (0) {
     int i;

     double b;
     double c;
     double d;

     double **w1;
     double **w2;

     w1 = get_work1(&work, WORK_DXX);

     b = a * (1. + (i_four == 0 ? 1. : 0.)) / 4.;

     c = b * omega;

     dmat_diag_mul(v1, P_m, w1, n_quad, n_quad);
     dmat_mul_diag(w1, qw_v, w1, n_quad, n_quad);
     dmat_scale(c, w1, r, n_quad, n_quad);

     if (flags_or(derivs, n_derivs)) {
          w2 = get_work1(&work, WORK_DXX);

          for (i = 0; i < n_derivs; ++i) {
               if (! derivs[i])
                    continue;

               d = b * omega_l[i];
               dmat_scale(d, w1, r_l[i], n_quad, n_quad);

               dmat_diag_mul(v1, P_m_l[i], w2, n_quad, n_quad);
               dmat_mul_diag(w2, qw_v, w2, n_quad, n_quad);
               dmat_scale(c, w2, w2, n_quad, n_quad);

               dmat_add(r_l[i], w2, r_l[i], n_quad, n_quad);
          }
     }
}
else {
     int i;

     double b;
     double c;
     double d;

     double **w1;
     double **w2;

     w1 = get_work1(&work, WORK_DXX);

     b = a * (1. + (i_four == 0 ? 1. : 0.)) / 4.;

     c = b * omega;

     dmat_diag_mul(v1, P_m, w1, n_quad, n_quad);
     dmat_mul_diag(w1, qw_v, r, n_quad, n_quad);
     dmat_scale(c, r, r, n_quad, n_quad);

     if (flags_or(derivs, n_derivs)) {
          w2 = get_work1(&work, WORK_DXX);

          for (i = 0; i < n_derivs; ++i) {
               if (! derivs[i])
                    continue;

               d = b * omega_l[i];
               dmat_scale(d, P_m, w1, n_quad, n_quad);
               dmat_scale(c, P_m_l[i], w2, n_quad, n_quad);
               dmat_add(w1, w2, r_l[i], n_quad, n_quad);

               dmat_diag_mul(v1, r_l[i], r_l[i], n_quad, n_quad);
               dmat_mul_diag(r_l[i], qw_v, r_l[i], n_quad, n_quad);
          }
     }
}
}



void build_local_r_and_t(int i_four, int n_quad, int n_derivs,
                         double *qx_v, double *qw_v, double omega, double *omega_l,
                         double **P_p, double **P_m, double **r, double **t,
                         double ***P_p_l, double ***P_m_l, double ***r_l, double ***t_l,
                         uchar *derivs, work_data work) {

     double *v1;

     v1 = get_work1(&work, WORK_DX);

     dvec_inv(qx_v, v1, n_quad);

     build_local_r_or_t_mat(i_four, n_quad, n_derivs, v1, qw_v, -1.,
                            omega, omega_l, P_m, P_m_l, r, r_l, derivs, work);

     build_local_r_or_t_mat(i_four, n_quad, n_derivs, v1, qw_v,  1.,
                            omega, omega_l, P_p, P_p_l, t, t_l, derivs,  work);
     dmat_sub_diag(t, v1, t, n_quad);
#ifdef USE_AD_FOR_TL_CALC_BUILD_LOCAL_R_AND_T
     build_local_r_and_t_tl_with_ad(i_four, n_quad, n_derivs, qx_v, qw_v, omega, omega_l, P_p, P_m, r, t, P_p_l, P_m_l, r_l, t_l, derivs, work);
#endif
}



/*******************************************************************************
 *
 ******************************************************************************/
void build_txr(int n_quad, int n_stokes, int n_derivs,
               double **r_p, double **t_p, double **tpr, double **tmr,
               double ***r_p_l, double ***t_p_l, double ***tpr_l, double ***tmr_l,
               uchar *derivs, work_data work) {

     int i;

     int n_quad_v;

     double **w1;

     n_quad_v = n_quad * n_stokes;

     if (n_stokes <= 2) {
          dmat_add(t_p, r_p, tpr, n_quad_v, n_quad_v);
          dmat_sub(t_p, r_p, tmr, n_quad_v, n_quad_v);

          for (i = 0; i < n_derivs; ++i) {
               if (! derivs[i])
                    continue;

               dmat_add(t_p_l[i], r_p_l[i], tpr_l[i], n_quad_v, n_quad_v);
               dmat_sub(t_p_l[i], r_p_l[i], tmr_l[i], n_quad_v, n_quad_v);
          }
     }
     else {
          w1 = get_work1(&work, WORK_DXX);

          dmat_mul_D_A(n_quad, n_stokes, r_p, w1);

          dmat_add(t_p, w1, tpr, n_quad_v, n_quad_v);
          dmat_sub(t_p, w1, tmr, n_quad_v, n_quad_v);

          for (i = 0; i < n_derivs; ++i) {
               if (! derivs[i])
                    continue;

               dmat_mul_D_A(n_quad, n_stokes, r_p_l[i], w1);

               dmat_add(t_p_l[i], w1, tpr_l[i], n_quad_v, n_quad_v);
               dmat_sub(t_p_l[i], w1, tmr_l[i], n_quad_v, n_quad_v);
          }
     }
#ifdef USE_AD_FOR_TL_CALC_BUILD_TXR
     build_txr_tl_with_ad(n_quad, n_stokes, n_derivs, r_p, t_p, tpr, tmr, r_p_l, t_p_l, tpr_l, tmr_l, derivs, work);
#endif
}



/*******************************************************************************
 *
 ******************************************************************************/
void build_gamma(int n_quad, int n_derivs,
                 double **tpr, double **tmr, double **gamma,
                 double ***tpr_l, double ***tmr_l, double ***gamma_l,
                 uchar *derivs, work_data work) {

     dmat_mul(tmr, tpr, n_quad, n_quad, n_quad, gamma);

     if (flags_or(derivs, n_derivs)) {
          int i;

          double **w1;
          double **w2;

          w1 = get_work1(&work, WORK_DXX);
          w2 = get_work1(&work, WORK_DXX);

          for (i = 0; i < n_derivs; ++i) {
               if (! derivs[i])
                    continue;

               dmat_mul(tmr, tpr_l[i], n_quad, n_quad, n_quad, w1);
               dmat_mul(tmr_l[i], tpr, n_quad, n_quad, n_quad, w2);
               dmat_add(w1, w2, gamma_l[i], n_quad, n_quad);
          }
     }
#ifdef USE_AD_FOR_TL_CALC_BUILD_GAMMA
     build_gamma_tl_with_ad(n_quad, n_derivs, tpr, tmr, gamma, tpr_l, tmr_l, gamma_l, derivs, work);
#endif
}



/*******************************************************************************
 *
 ******************************************************************************/
void build_sim_vectors(int n_quad, int n_umus,
                       double *qx_v, double *qw_v, double *v1, double *v2) {

     int i;

     double a;

     for (i = 0; i < n_quad; ++i) {
          a = sqrt(qx_v[i] * qw_v[i]);
          v1[i] = a;
          a = 1. / a;
          v2[i] = a;
     }

     for ( ; i < n_quad + n_umus; ++i) {
          a = sqrt(qx_v[i]);
          v1[i] = a;
          a = 1. / a;
          v2[i] = a;
     }
}



void vec_sim_trans(int m, double *a, double *s) {

     dm_v_diag_mul(s, a, a, m);
}


void mat_sim_trans(int m, int n, double **a, double *s1, double *s2) {

     dmat_diag_mul(s1, a, a, m, n);
     dmat_mul_diag(a, s2, a, m, n);
}



void vec_sim_trans2(int m, double *a1, double *a2, double *s1) {

     dm_v_diag_mul(s1, a1, a1, m);

     dm_v_diag_mul(s1, a2, a2, m);
}



void mat_sim_trans2(int m, int n, double **a1, double **a2, double *s1, double *s2) {

     dmat_diag_mul(s1, a1, a1, m, n);
     dmat_mul_diag(a1, s2, a1, m, n);

     dmat_diag_mul(s1, a2, a2, m, n);
     dmat_mul_diag(a2, s2, a2, m, n);
}



/*******************************************************************************
 *
 ******************************************************************************/
void dm_v_mul_D_A(int n_quad, int n_stokes, double *a, double *b) {

     int i;
     int ii;
     int j;
     int jj;

     int n;

     if (b == a) {
          ii = 0;
          for (i = 0; i < n_quad; ++i) {
               jj = ii + 2;

               for (j = 2; j < n_stokes; ++j) {
                    a[jj] = -a[jj];
                    jj++;
               }
               ii += n_stokes;
          }
     }
     else {
          n = MIN(2, n_stokes);

          ii = 0;
          for (i = 0; i < n_quad; ++i) {
               jj = ii;
               for (j = 0; j < n; ++j) {
                    b[jj] =  a[jj];
                    jj++;
               }

               for (j = 2; j < n_stokes; ++j) {
                    b[jj] = -a[jj];
                    jj++;
               }
               ii += n_stokes;
          }
     }
}



void dmat_mul_D_A(int n_quad, int n_stokes, double **a, double **b) {

     dmat_mul_D_A2(n_quad, n_stokes, n_quad, n_stokes, a, b);
}



void dmat_mul_D_A2(int n_quad1, int n_stokes1, int n_quad2, int n_stokes2, double **a, double **b) {

     int i;
     int ii;
     int j;
     int jj;
     int k;

     int n;

     int n_quad_v2;

     n_quad_v2 = n_quad2 * n_stokes2;

     if (b == a) {
          ii = 0;
          for (i = 0; i < n_quad1; ++i) {
               jj = ii + 2;

               for (j = 2; j < n_stokes1; ++j) {
                    for (k = 0; k < n_quad_v2; ++k) {
                         a[jj][k] = -a[jj][k];
                    }
                    jj++;
               }
               ii += n_stokes1;
          }
     }
     else {
          n = MIN(2, n_stokes1);
          ii = 0;
          for (i = 0; i < n_quad1; ++i) {
               jj = ii;
               for (j = 0; j < n; ++j) {
                    for (k = 0; k < n_quad_v2; ++k) {
                         b[jj][k] =  a[jj][k];
                    }
                    jj++;
               }

               for (j = 2; j < n_stokes1; ++j) {
                    for (k = 0; k < n_quad_v2; ++k) {
                         b[jj][k] = -a[jj][k];
                    }
                    jj++;
               }
               ii += n_stokes1;
          }
     }
}



/*******************************************************************************
 *
 ******************************************************************************/
/*
void dm_v_mul_D_A_2(int n_quad, int n_stokes, double *a, double *b) {

     int i;
     int ii;
     int j;
     int jj;

     int n;

     if (b == a) {
          ii = 0;
          for (i = 0; i < n_quad; ++i) {
               jj = ii + 2;

               for (j = 2; j < n_stokes; ++j) {
                    a[jj] = -a[jj];
                    jj++;
               }
               ii += n_stokes;
          }
     }
     else {
          n = MIN(2, n_stokes);

          ii = 0;
          for (i = 0; i < n_quad; ++i) {
               jj = ii;
               for (j = 0; j < n; ++j) {
                    b[jj] =  a[jj];
                    jj++;
               }

               for (j = 2; j < n_stokes; ++j) {
                    b[jj] = -a[jj];
                    jj++;
               }
               ii += n_stokes;
          }
     }
}



void dmat_mul_D_A_2(int n_quad, int n_stokes, double **a, double **b) {

     dmat_mul_D_A2_2(n_quad, n_stokes, n_quad, n_stokes, a, b);
}



void dmat_mul_D_A2_2(int n_quad1, int n_stokes1, int n_quad2, int n_stokes2, double **a, double **b) {

     int i;
     int ii;
     int j;
     int jj;
     int k;

     int n;

     int n_quad_v2;

     n_quad_v2 = n_quad2 * n_stokes2;

     if (b == a) {
          ii = 0;
          for (i = 0; i < n_quad1; ++i) {
               jj = ii + 2;

               for (j = 2; j < n_stokes1; ++j) {
                    for (k = 0; k < n_quad_v2; ++k) {
                         a[jj][k] = -a[jj][k];
                    }
                    jj++;
               }
               ii += n_stokes1;
          }
     }
     else {
          n = MIN(2, n_stokes1);
          ii = 0;
          for (i = 0; i < n_quad1; ++i) {
               jj = ii;
               for (j = 0; j < n; ++j) {
                    for (k = 0; k < n_quad_v2; ++k) {
                         b[jj][k] =  a[jj][k];
                    }
                    jj++;
               }

               for (j = 2; j < n_stokes1; ++j) {
                    for (k = 0; k < n_quad_v2; ++k) {
                         b[jj][k] = -a[jj][k];
                    }
                    jj++;
               }
               ii += n_stokes1;
          }
     }
}
*/


/*******************************************************************************
 *
 ******************************************************************************/
void no_scatter_r_t_s(int n_quad, int n_derivs,
                      double *qx_v, double *qw_v, double ltau, double *ltau_l,
                      double **R_p, double **T_p, double **R_m, double **T_m,
                      double *S_p, double *S_m,
                      double ***R_p_l, double ***T_p_l, double ***R_m_l, double ***T_m_l,
                      double **S_p_l, double **S_m_l, uchar *derivs_h, uchar *derivs_p) {
     int i;
     int j;

     dmat_zero(R_p, n_quad, n_quad);
     dmat_zero(T_p, n_quad, n_quad);
     dmat_zero(R_m, n_quad, n_quad);
     dmat_zero(T_m, n_quad, n_quad);
     dvec_zero(S_p, n_quad);
     dvec_zero(S_m, n_quad);

     for (i = 0; i < n_quad; ++i)
          T_p[i][i] = T_m[i][i] = exp(-ltau / qx_v[i]);

     if (flags_or(derivs_h, n_derivs)) {
          for (i = 0; i < n_derivs; ++i) {
               if (! derivs_h[i])
                    continue;

               dmat_zero(R_p_l[i], n_quad, n_quad);
               dmat_zero(T_p_l[i], n_quad, n_quad);
               dmat_zero(R_m_l[i], n_quad, n_quad);
               dmat_zero(T_m_l[i], n_quad, n_quad);

               for (j = 0; j < n_quad; ++j)
                    T_p_l[i][j][j] = T_m_l[i][j][j] = -ltau_l[i] / qx_v[j] * T_p[j][j];
          }
     }

     if (flags_or(derivs_p, n_derivs)) {
          for (i = 0; i < n_derivs; ++i) {
               if (! derivs_p[i])
                    continue;

               dvec_zero(S_p_l[i], n_quad);
               dvec_zero(S_m_l[i], n_quad);
          }
     }
}
