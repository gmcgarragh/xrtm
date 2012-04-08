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

#include "xrtm.h"
#include "xrtm_scatter.h"
#include "xrtm_scatter_a.h"
#include "xrtm_utility.h"


/*******************************************************************************
 *
 ******************************************************************************/
void phase_func_a(int n_coef, double *p, double *chi_a, double P_a) {

     int i;

     for (i = 0; i < n_coef; ++i)
          chi_a[i] += p[i] * P_a;
}



/*******************************************************************************
 *
 ******************************************************************************/
void build_phase_vecs_scalar_a(int i_four, int n_coef, int n_mus, double **Y1, double *Y2, int lda, double *chi_a, double *P_pp_a, double *P_mp_a) {

     int i;
     int k;

     double a;
     double b;
     double c;

     a = 2. - (i_four == 0 ? 1. : 0.);
     for (i = 0; i < n_mus; ++i) {
          b = P_pp_a[i] * a;
          c = P_mp_a[i] * a;
          P_pp_a[i] = b + c;
          P_mp_a[i] = b - c;
     }

     for (k = i_four    ; k < n_coef; k += 2) {
          a = 0.;
          for (i = 0; i < n_mus; ++i) {
               a += P_pp_a[i] * Y1[i][k];
          }
          chi_a[k] += a * Y2[k];
     }

     for (k = i_four + 1; k < n_coef; k += 2) {
          a = 0.;
          for (i = 0; i < n_mus; ++i) {
               a += P_mp_a[i] * Y1[i][k];
          }
          chi_a[k] += a * Y2[k];
     }

     for (i = 0; i < n_mus; ++i) {
          P_pp_a[i] = 0.;
          P_mp_a[i] = 0.;
     }
}



/*******************************************************************************
 *
 ******************************************************************************/
void build_phase_mats_scalar_a(int i_four, int n_coef, int n_mus1, int n_mus2, double **Y1, double **Y2, int lda, double *chi_a, double **P_pp_a, double **P_mp_a) {

     int i;
     int j;
     int k;

     double a;
     double b;
     double c;

     a = 2. - (i_four == 0 ? 1. : 0.);
     for (i = 0; i < n_mus1; ++i) {
          for (j = 0; j < n_mus2; ++j) {
               b = P_pp_a[i][j] * a;
               c = P_mp_a[i][j] * a;
               P_pp_a[i][j] = b + c;
               P_mp_a[i][j] = b - c;
          }
     }

     for (k = i_four    ; k < n_coef; k += 2) {
          for (i = 0; i < n_mus1; ++i) {
               a = 0.;
               for (j = 0; j < n_mus2; ++j) {
                    a += P_pp_a[i][j] * Y2[j][k];
               }
               chi_a[k] += a * Y1[i][k];
          }
     }

     for (k = i_four + 1; k < n_coef; k += 2) {
          for (i = 0; i < n_mus1; ++i) {
               a = 0.;
               for (j = 0; j < n_mus2; ++j) {
                    a += P_mp_a[i][j] * Y2[j][k];
               }
               chi_a[k] += a * Y1[i][k];
          }
     }

     for (i = 0; i < n_mus1; ++i) {
          for (j = 0; j < n_mus2; ++j) {
               P_pp_a[i][j] = 0.;
               P_mp_a[i][j] = 0.;
          }
     }
}



/*******************************************************************************
 *
 ******************************************************************************/
void build_scat_vector_gc_a(int n_coef, int n_stokes, double **gsf, double **chi_a, double *F_a) {

     int i;
     int ii;

     ii = 0;

     for (i = 0; i < n_coef; ++i)
          chi_a[0][i] += F_a[ii] * gsf[0][i];
     ii++;

     if (n_stokes > 1) {
          for (i = 2; i < n_coef; ++i)
               chi_a[4][i] += F_a[ii] * gsf[1][i];
          ii++;
     }
}



/*******************************************************************************
 *
 ******************************************************************************/
static void gc_to_B_a(int i, int n_stokes, double **gc_a, double **B_a, double c) {

     gc_a[0][i] += c * B_a[0][0];
     if (n_stokes >= 3) {
          gc_a[1][i] += c *   B_a[1][1];
          gc_a[4][i] += c * (-B_a[0][1] - B_a[1][0]);
          gc_a[2][i] += c *   B_a[2][2];
     }
     if (n_stokes >= 4) {
          gc_a[3][i] += c *   B_a[3][3];
          gc_a[5][i] += c * ( B_a[2][3] - B_a[3][2]);
     }
}



/*******************************************************************************
 *
 ******************************************************************************/
void build_phase_vecs_vector_gc_a(int i_four, int n_coef, int n_mus, int n_stokes, double *mu, double mu_0, double **gc_a, double *P_pp_a, double *P_mp_a, work_data work) {

     int i;
     int ii;
     int j;
     int k;

     int n_mus_v;

     int n_stokes1;
     int n_stokes2;

     double a;
     double b;
     double c;

     double D[] = {1., 1., -1., -1.};

     double *v1;

     double **w1;
     double **w2;

     double **l0;

     double **B;

     double ***A;
     double ***ll;

     if (n_stokes == 2) {
          n_stokes1 = 2;
          n_stokes2 = 3;
     }
     else {
          n_stokes1 = n_stokes;
          n_stokes2 = n_stokes;
     }

     n_mus_v = n_mus * n_stokes1;

     if (i_four >= n_coef) {
          init_array2_d(gc_a, 6, n_coef, 0.);

          return;
     }

     v1 = get_work_d1(&work, n_stokes2);

     w1 = get_work_d2(&work, n_stokes2, n_mus_v);
     w2 = get_work_d2(&work, n_stokes2, n_stokes2);

     l0 = get_work_d2(&work, n_coef, n_stokes2);

     B  = get_work_d2(&work, n_stokes2, n_stokes2);

     A  = get_work_d3(&work, n_coef, n_stokes2, n_stokes2);

     ll = get_work_d3(&work, n_coef, n_mus_v,    n_stokes2);

     ii = 0;
     for (i = 0; i < n_mus; ++i) {
          basic_matrix(i_four, n_coef, n_stokes2, mu[i], A);

          for (j = i_four; j < n_coef; ++j) {
               dmat_copy(&ll[j][ii], A[j], n_stokes1, n_stokes2);
          }

          ii += n_stokes1;
     }

     basic_matrix(i_four, n_coef, 1, mu_0, A);
     for (i = i_four; i < n_coef; ++i) {
          l0[i][0] = A[i][0][0];
          for (j = 1; j < n_stokes2; ++j) {
               l0[i][j] = 0.;
          }
     }

     dm_v_mul_D_A(n_mus, n_stokes1, P_mp_a, P_mp_a);

     a = 1.;
     b = 1. / factorial(2 * i_four);
     c = 2. - (i_four == 0 ? 1. : 0.);
     for (i = i_four; i < n_coef; ++i) {
          dmat_trans(ll[i], w1, n_mus_v, n_stokes2);

          dm_v_mul(w1, P_pp_a, n_stokes2, n_mus_v, v1);
          for (j = 0; j < n_stokes2; ++j) {
               for (k = 0; k < n_stokes2; ++k) {
                    B[j][k] = v1[j] * l0[i][k];
               }
          }

          dm_v_mul(w1, P_mp_a, n_stokes2, n_mus_v, v1);
          for (j = 0; j < n_stokes2; ++j) {
               for (k = 0; k < n_stokes2; ++k) {
                    w2[j][k]  = v1[j] * l0[i][k];
               }
          }
          dmat_mul_diag(w2, D, w2, n_stokes2, n_stokes2);
          if (a > 0.)
               dmat_add(B, w2, B, n_stokes2, n_stokes2);
          else
               dmat_sub(B, w2, B, n_stokes2, n_stokes2);

          gc_to_B_a(i, n_stokes2, gc_a, B, b * c);

          a *= -1.;

          b *= (i + 1 - i_four) / (i + 1. + i_four);
     }

     dvec_zero(P_pp_a, n_mus_v);
     dvec_zero(P_mp_a, n_mus_v);
}



/*******************************************************************************
 *
 ******************************************************************************/
void build_phase_mats_vector_gc_a(int i_four, int n_coef, int n_mus1, int n_mus2, int n_stokes, double *mu1, double *mu2, double **gc_a, double **P_pp_a, double **P_mp_a, double **P_mm_a, double **P_pm_a, int vector, work_data work) {

     int i;
     int ii;
     int j;

     int n_mus_v1;
     int n_mus_v2;

     int n_stokes1;
     int n_stokes2;

     double a;
     double b;
     double c;

     double D[] = {1., 1., -1., -1.};

     double **w1;
     double **w2;
     double **w3;

     double **B;

     double ***A;
     double ***ll;
     double ***ll2;

     if (n_stokes == 2) {
          n_stokes1 = 2;
          n_stokes2 = 3;
     }
     else {
          n_stokes1 = n_stokes;
          n_stokes2 = n_stokes;
     }

     n_mus_v1 = n_mus1 * n_stokes1;
     n_mus_v2 = n_mus2 * n_stokes1;

     w1  = get_work_d2(&work, n_stokes2, n_mus_v1);
     w2  = get_work_d2(&work, n_stokes2, n_mus_v2);
     w3  = get_work_d2(&work, n_stokes2, n_stokes2);

     B   = get_work_d2(&work, n_stokes2, n_stokes2);

     A   = get_work_d3(&work, n_coef, n_stokes2, n_stokes2);

     ll  = get_work_d3(&work, n_coef, n_mus_v1, n_stokes2);
     ll2 = get_work_d3(&work, n_coef, n_mus_v2, n_stokes2);

     ii = 0;
     for (i = 0; i < n_mus1; ++i) {
          basic_matrix(i_four, n_coef, n_stokes2, mu1[i], A);

          for (j = i_four; j < n_coef; ++j) {
               dmat_copy(&ll[j][ii], A[j], n_stokes1, n_stokes2);
          }

          ii += n_stokes1;
     }

     ii = 0;
     for (i = 0; i < n_mus2; ++i) {
          basic_matrix(i_four, n_coef, n_stokes2, mu2[i], A);

          for (j = i_four; j < n_coef; ++j) {
               dmat_copy(&ll2[j][ii], A[j], n_stokes1, n_stokes2);
          }

          ii += n_stokes1;
     }

     if (vector)
          phase_matrix_symmetry_a3(n_mus1, n_stokes, n_mus2, n_stokes, P_pp_a, P_mp_a, P_mm_a, P_pm_a, 1.);

     dmat_mul_D_A2(n_mus1, n_stokes, n_mus2, n_stokes, P_mp_a, P_mp_a);

     a =  1.;
     b =  1. / factorial(2 * i_four);
     c =  2. - (i_four == 0 ? 1. : 0.);
     for (i = i_four; i < n_coef; ++i) {
          dmat_zero(B, n_stokes2, n_stokes2);

          dmat_trans(ll[i], w1, n_mus_v1, n_stokes2);

          dmat_mul(w1, P_pp_a, n_stokes2, n_mus_v2, n_mus_v1, w2);
          dmat_mul(w2, ll2[i], n_stokes2, n_stokes2, n_mus_v2, w3);
          dmat_add(B, w3, B, n_stokes2, n_stokes2);

          dmat_mul(w1, P_mp_a, n_stokes2, n_mus_v2, n_mus_v1, w2);
          dmat_mul(w2, ll2[i], n_stokes2, n_stokes2, n_mus_v2, w3);
          dmat_mul_diag(w3, D, w3, n_stokes2, n_stokes2);

          if (a > 0.)
               dmat_add(B, w3, B, n_stokes2, n_stokes2);
          else
               dmat_sub(B, w3, B, n_stokes2, n_stokes2);

          gc_to_B_a(i, n_stokes2, gc_a, B, b * c);

          a *= -1.;

          b *= (i + 1 - i_four) / (i + 1. + i_four);
     }

     dmat_zero(P_pp_a, n_mus_v1, n_mus_v2);
     dmat_zero(P_mp_a, n_mus_v1, n_mus_v2);
     dmat_zero(P_mm_a, n_mus_v1, n_mus_v2);
     dmat_zero(P_pm_a, n_mus_v1, n_mus_v2);
}



/*******************************************************************************
 *
 ******************************************************************************/
void phase_matrix_symmetry_a3(int n_quad1, int n_stokes1, int n_quad2, int n_stokes2,
                              double **P_pp_a, double **P_mp_a,
                              double **P_mm_a, double **P_pm_a, double f) {

     phase_matrix_symmetry_a_ldx3(n_quad1, n_stokes1, n_quad2, n_stokes2,
                                  *P_pp_a, *P_mp_a, *P_mm_a, *P_pm_a, n_quad2 * n_stokes2, f);
}



void phase_matrix_symmetry_a_ldx3(int n_quad1, int n_stokes1, int n_quad2, int n_stokes2,
                                  double *P_pp_a, double *P_mp_a,
                                  double *P_mm_a, double *P_pm_a, int ldx, double f) {

     int i;
     int ii;
     int j;
     int jj;

     int n_quad_v1;
     int n_quad_v2;

     double a;
     double b;

     n_quad_v1 = n_quad1 * n_stokes1;
     n_quad_v2 = n_quad2 * n_stokes2;
/*
     if (n_stokes1 != 4 || n_stokes2 != 4) {
*/
          ii = 0;
          for (i = 0; i < n_quad_v1; ++i) {
               a = f * (i % n_stokes1 >= 2 ? -1. : 1.);
               jj = ii;
               for (j = 0; j < n_quad_v2; ++j) {
                    b = a * (j % n_stokes2 >= 2 ? -1. : 1.);
                    P_pp_a[jj] += b * P_mm_a[jj];
                    P_mp_a[jj] += b * P_pm_a[jj];
                    jj++;
               }
               ii += ldx;
          }
/*
     }
     else {

     }
*/
}



/*******************************************************************************
 *
 ******************************************************************************/
#ifdef INCLUDE_DEV_SOURCE
#include "../src2/xrtm_scatter_a2.c"
#endif
