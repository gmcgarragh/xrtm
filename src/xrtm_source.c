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
#include "xrtm_save_tree.h"
#include "xrtm_source.h"
#include "xrtm_source_a.h"
#include "xrtm_utility.h"


/*******************************************************************************
 *
 ******************************************************************************/
void build_source_vectors_1n(int n_quad, int n_stokes, int n_derivs, double *qx_v,
                             double F_0, double omega, double *omega_l,
                             double as_0, double *as_0_l,
                             double  *P_q0_mm, double  *P_q0_pm,
                             double  **tpr, double  **tmr, double  **gamma,
                             double  *F_p, double  *F_m,
                             double **P_q0_mm_l, double **P_q0_pm_l,
                             double ***tpr_l, double ***tmr_l, double ***gamma_l,
                             double **F_p_l, double **F_m_l,
                             uchar *derivs_h, uchar *derivs_p,
                             save_tree_data save_tree, work_data work) {

     int i;
     int j;

     int n_quad_v;

     int *i1;

     double aa;
     double bb;
     double cc;

     double *v1;
     double *v2;

     double *d;
     double *e;
     double *h;
     double *p;

     double **f;

     double *d_l;
     double *e_l;
     double *h_l;
     double *p_l;

     forward_save_build_source_vectors_1n_data *save;


     n_quad_v = n_quad * n_stokes;


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (save_tree.t) {
          save_tree_encode_s(&save_tree, "build_source_vectors_1n");

          if (save_tree_retrieve_data(&save_tree, forward_save_build_source_vectors_1n_data, &save))
               forward_save_build_source_vectors_1n_alloc(save, n_quad_v);
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     i1  = (int *) get_work1(&work, WORK_IX);

     v1  = get_work1(&work, WORK_DX);
     v2  = get_work1(&work, WORK_DX);

     d   = get_work1(&work, WORK_DX);
     e   = get_work1(&work, WORK_DX);
     p   = get_work1(&work, WORK_DX);
     h   = get_work1(&work, WORK_DX);

     f   = get_work1(&work, WORK_DXX);


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     dmat_copy(f, gamma, n_quad_v, n_quad_v);
     aa = as_0*as_0;
     for (i = 0; i < n_quad_v; ++i)
          f[i][i] -= aa;

     dmat_getrf(f, n_quad_v, n_quad_v, i1);

     if (save_tree.t) {
          copy_array1_i(save->ip, i1, n_quad_v);
          dmat_copy(save->f, f, n_quad_v, n_quad_v);
     }

     aa = F_0 / (4. * PI);

     bb = aa * omega;

     for (i = 0; i < n_quad_v; ++i) {
          cc = bb / qx_v[i];
          v1[i] = cc * P_q0_pm[i];
          v2[i] = cc * P_q0_mm[i];
     }

     dm_v_mul_D_A(n_quad, n_stokes, v2, v2);

     dvec_add(v1, v2, d, n_quad_v);
     dvec_sub(v1, v2, e, n_quad_v);

     dm_v_mul(tmr, e, n_quad_v, n_quad_v, v1);
     dvec_scale(-as_0, d, v2, n_quad_v);
     dvec_sub(v2, v1, p, n_quad_v);

     dmat_getrs(f, &p, n_quad_v, 1, i1);

     dm_v_mul(tpr, p, n_quad_v, n_quad_v, v1);
     dvec_add(v1, e, h, n_quad_v);
     dvec_scale(1. / as_0, h, v1, n_quad_v);

     dvec_add(v1, p, v1, n_quad_v);
     dvec_scale(.5, v1, F_p, n_quad_v);

     dvec_sub(F_p, p, F_m, n_quad_v);

     dm_v_mul_D_A(n_quad, n_stokes, F_m, F_m);

     if (save_tree.t) {
          dvec_copy(save->d, d, n_quad_v);
          dvec_copy(save->e, e, n_quad_v);
          dvec_copy(save->p, p, n_quad_v);
          dvec_copy(save->h, h, n_quad_v);
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (flags_or(derivs_h, n_derivs)) {
          d_l = get_work1(&work, WORK_DX);
          e_l = get_work1(&work, WORK_DX);
     }

     if (flags_or(derivs_p, n_derivs)) {
          h_l = get_work1(&work, WORK_DX);
          p_l = get_work1(&work, WORK_DX);
     }

     for (i = 0; i < n_derivs; ++i) {
          if (derivs_h[i]) {
               cc = aa * omega_l[i];

               for (j = 0; j < n_quad_v; ++j) {
                    v1[j] = (cc * P_q0_pm[j] + bb * P_q0_pm_l[i][j]) / qx_v[j];
                    v2[j] = (cc * P_q0_mm[j] + bb * P_q0_mm_l[i][j]) / qx_v[j];
               }

               dm_v_mul_D_A(n_quad, n_stokes, v2, v2);

               dvec_add(v1, v2, d_l, n_quad_v);
               dvec_sub(v1, v2, e_l, n_quad_v);

               dm_v_mul(gamma_l[i], p, n_quad_v, n_quad_v, v1);
               dvec_scale(2. * as_0 * as_0_l[i], p, v2, n_quad_v);
               dvec_sub(v2, v1, v1, n_quad_v);

               dm_v_mul(tmr_l[i], e, n_quad_v, n_quad_v, v2);
               dvec_sub(v1, v2, v1, n_quad_v);

               dm_v_mul(tmr, e_l, n_quad_v, n_quad_v, v2);
               dvec_sub(v1, v2, v1, n_quad_v);

               dvec_scale(as_0_l[i], d, v2, n_quad_v);
               dvec_sub(v1, v2, v1, n_quad_v);

               dvec_scale(as_0, d_l, v2, n_quad_v);
               dvec_sub(v1, v2, p_l, n_quad_v);

               dmat_getrs(f, &p_l, n_quad_v, 1, i1);

               dm_v_mul(tpr_l[i], p, n_quad_v, n_quad_v, v1);
               dm_v_mul(tpr, p_l, n_quad_v, n_quad_v, v2);
               dvec_add(v1, v2, v1, n_quad_v);
               dvec_add(v1, e_l, h_l, n_quad_v);

               dvec_scale(-as_0_l[i] / (as_0 * as_0), h, v1, n_quad_v);
               dvec_scale(1. / as_0, h_l, v2, n_quad_v);
               dvec_add(v1, v2, v1, n_quad_v);
               dvec_add(v1, p_l, v1, n_quad_v);
               dvec_scale(.5, v1, F_p_l[i], n_quad_v);

               dvec_sub(F_p_l[i], p_l, F_m_l[i], n_quad_v);

               dm_v_mul_D_A(n_quad, n_stokes, F_m_l[i], F_m_l[i]);
          }
          else
          if (derivs_p[i] && as_0_l[i] != 0.) {
               dvec_scale(2. * as_0 * as_0_l[i], p, v1, n_quad_v);
               dvec_scale(as_0_l[i], d, v2, n_quad_v);
               dvec_sub(v1, v2, p_l, n_quad_v);

               dmat_getrs(f, &p_l, n_quad_v, 1, i1);

               dm_v_mul(tpr, p_l, n_quad_v, n_quad_v, h_l);

               dvec_scale(-as_0_l[i] / (as_0 * as_0), h, v1, n_quad_v);
               dvec_scale(1. / as_0, h_l, v2, n_quad_v);
               dvec_add(v1, v2, v1, n_quad_v);
               dvec_add(v1, p_l, v1, n_quad_v);
               dvec_scale(.5, v1, F_p_l[i], n_quad_v);

               dvec_sub(F_p_l[i], p_l, F_m_l[i], n_quad_v);

               dm_v_mul_D_A(n_quad, n_stokes, F_m_l[i], F_m_l[i]);
          }
          else
          if (derivs_p[i] && as_0_l[i] == 0.) {
               dvec_zero(F_p_l[i], n_quad_v);
               dvec_zero(F_m_l[i], n_quad_v);
          }
     }
#ifdef USE_AD_FOR_TL_CALC_BUILD_SOURCE_VECTORS_1N
     build_source_vectors_1n_tl_with_ad(n_quad, n_stokes, n_derivs, qx_v,
                                     F_0, omega, omega_l,
                                     as_0, as_0_l,
                                     P_q0_mm, P_q0_pm,
                                     tpr, tmr, gamma,
                                     F_p, F_m,
                                     P_q0_mm_l, P_q0_pm_l,
                                     tpr_l, tmr_l, gamma_l,
                                     F_p_l, F_m_l, derivs_h, derivs_p, save_tree, work);
#endif
}



/*******************************************************************************
 *
 ******************************************************************************/
void build_source_vectors_2n(int n_quad_v, int n_derivs, double *qx_v, double F_0,
                        double omega, double *omega_l, double as_0, double *as_0_l,
                        double  *P_q0_mm, double  *P_q0_pm,
                        double  **r, double  **t, double  *F_p, double  *F_m,
                        double **P_q0_mm_l, double **P_q0_pm_l,
                        double ***r_l, double ***t_l, double **F_p_l, double **F_m_l,
                        uchar *derivs_h, uchar *derivs_p, work_data work) {

     build_source_vectors_2n2(n_quad_v, n_derivs, qx_v, F_0, omega, omega_l, as_0, as_0_l, P_q0_mm, P_q0_pm, r, t, r, t, F_p, F_m, P_q0_mm_l, P_q0_pm_l, r_l, t_l, r_l, t_l, F_p_l, F_m_l, derivs_h, derivs_p, work);
}



void build_source_vectors_2n2(int n_quad_v, int n_derivs, double *qx_v, double F_0,
                         double omega, double *omega_l, double as_0, double *as_0_l,
                         double  *P_q0_mm, double  *P_q0_pm,
                         double  **r_p, double  **t_p, double  **r_m, double  **t_m,
                         double  *F_p, double  *F_m,
                         double **P_q0_mm_l, double **P_q0_pm_l,
                         double ***r_p_l, double ***t_p_l, double ***r_m_l, double ***t_m_l,
                         double **F_p_l, double **F_m_l,
                         uchar *derivs_h, uchar *derivs_p, work_data work) {

     int i;
     int ii;
     int j;
     int jj;
     int k;
     int kk;

     int *i1;

     int n_quad_v2;

     double a;
     double b;
     double c;

     double *v1;
     double *v2;

     double **w1;
     double **w2;

     n_quad_v2 = n_quad_v * 2;

     i1 = (int *) get_work1(&work, WORK_IX2);

     v1 = get_work1(&work, WORK_DX2);
     if (flags_or(derivs_h, n_derivs))
          v2 = get_work1(&work, WORK_DX2);

     w1 = get_work1(&work, WORK_DXX2);

     ii = n_quad_v;
     for (i = 0; i < n_quad_v; ++i) {
          jj = n_quad_v;
          for (j = 0; j < n_quad_v; ++j) {
               w1[i ][j ] = -t_p[i][j];
               w1[i ][jj] =  r_m[i][j];
               w1[ii][j ] = -r_p[i][j];
               w1[ii][jj] =  t_m[i][j];

               if (i == j) {
                    w1[i ][j ] += as_0;
                    w1[ii][jj] += as_0;
               }

               jj++;
          }
          ii++;
     }

     dmat_getrf(w1, n_quad_v2, n_quad_v2, i1);

     a = F_0 / (4. * PI);

     b = a * omega;

     ii = n_quad_v;
     for (i = 0; i < n_quad_v; ++i) {
          c = b / qx_v[i];
          v1[i ] =  c * P_q0_pm[i];
          v1[ii] = -c * P_q0_mm[i];
          ii++;
     }

     dmat_getrs(w1, &v1, n_quad_v2, 1, i1);

     dvec_copy(F_p, v1,        n_quad_v);
     dvec_copy(F_m, v1+n_quad_v, n_quad_v);

     if (flags_or(derivs_h, n_derivs))
          w2 = get_work1(&work, WORK_DXX2);

     for (i = 0; i < n_derivs; ++i) {
          if (derivs_h[i]) {
               jj = n_quad_v;
               for (j = 0; j < n_quad_v; ++j) {
                    kk = n_quad_v;
                    for (k = 0; k < n_quad_v; ++k) {
                         w2[j ][k ] =  t_p_l[i][j][k];
                         w2[j ][kk] = -r_m_l[i][j][k];
                         w2[jj][k ] =  r_p_l[i][j][k];
                         w2[jj][kk] = -t_m_l[i][j][k];

                         if (j == k) {
                              w2[j ][k ] -= as_0_l[i];
                              w2[jj][kk] -= as_0_l[i];
                         }

                        kk++;
                    }
                    jj++;
               }

               dm_v_mul(w2, v1, n_quad_v2, n_quad_v2, v2);

               c = omega_l[i] * a;

               jj = n_quad_v;
               for (j = 0; j < n_quad_v; ++j) {
                    v2[j ] += (c * P_q0_pm[j] + b * P_q0_pm_l[i][j]) / qx_v[j];
                    v2[jj] -= (c * P_q0_mm[j] + b * P_q0_mm_l[i][j]) / qx_v[j];
                    jj++;
               }

               dmat_getrs(w1, &v2, n_quad_v2, 1, i1);

               dvec_copy(F_p_l[i], v2,        n_quad_v);
               dvec_copy(F_m_l[i], v2+n_quad_v, n_quad_v);
          }
          else
          if (derivs_p[i]) {

          }
     }
}



/*******************************************************************************
 *
 ******************************************************************************/
void build_source_vectors_thermal(int n_quad, int n_stokes, int n_derivs, double *qx_v,
                             double F_0, double planck0, double planck1,
                             double omega, double *omega_l, double ltau, double *ltau_l,
                             double as_0, double *as_0_l,
                             double  *P_q0_mm, double  *P_q0_pm,
                             double  **r_p, double  **t_p, double  **r_m, double  **t_m,
                             double  *F0_p, double  *F0_m, double  *F1_p, double  *F1_m,
                             double **P_q0_mm_l, double **P_q0_pm_l,
                             double ***r_p_l, double ***t_p_l, double ***r_m_l, double ***t_m_l,
                             double  *F0_p_l, double  *F0_m_l, double  *F1_p_l, double  *F1_m_l,
                             int thermal, uchar *derivs_h, uchar *derivs_p, work_data work) {

     int i;
     int ii;
     int j;
     int jj;

     int *i1;

     int n_quad_v;
     int n_quad_v2;

     double b;
     double c;
     double d;

     double *v1;
     double *v2;

     double **w1;
/*
     double **w2;
*/
     n_quad_v  = n_quad * n_stokes;
     n_quad_v2 = n_quad_v * 2;

     i1 = (int *) get_work1(&work, WORK_IX2);

     v1 = get_work1(&work, WORK_DX2);
     if (thermal || flags_or(derivs_h, n_derivs))
          v2 = get_work1(&work, WORK_DX2);

     w1 = get_work1(&work, WORK_DXX2);

     ii = n_quad_v;
     for (i = 0; i < n_quad_v; ++i) {
          jj = n_quad_v;
          for (j = 0; j < n_quad_v; ++j) {
               w1[i ][j ] = -t_p[i][j];
               w1[i ][jj] =  r_m[i][j];
               w1[ii][j ] = -r_p[i][j];
               w1[ii][jj] =  t_m[i][j];
               jj++;
          }
          ii++;
     }

     dmat_getrf(w1, n_quad_v2, n_quad_v2, i1);

     b = 1. - omega;

     ii = n_quad_v;
     for (i = 0; i < n_quad_v; ++i) {
          v2[i ] = v2[ii] = 0.;
          if (i % n_stokes == 0) {
               c = b / qx_v[i];
               v2[i ] =  c;
               v2[ii] = -c;
          }
          ii++;
     }

     dmat_getrs(w1, &v2, n_quad_v2, 1, i1);
     dvec_copy(v1, v2, n_quad_v2);
     dmat_getrs(w1, &v2, n_quad_v2, 1, i1);

     b = (planck0 - planck1) / ltau;

     ii = n_quad_v;
     for (i = 0; i < n_quad_v; ++i) {
          c = b * v2[i ];
          d = b * v2[ii];
          F0_p[i] = planck0 * v1[i ] - c;
          F0_m[i] = planck0 * v1[ii] - d;
          F1_p[i] = planck1 * v1[i ] - c;
          F1_m[i] = planck1 * v1[ii] - d;
          ii++;
     }
/*
     if (flags_or(derivs_h, n_derivs))
          w2 = get_work1(&work, WORK_DXX2);

     for (i = 0; i < n_derivs; ++i) {
          if (derivs_h[i]) {

          }
          else
          if (derivs_p[i]) {

          }
     }
*/
}



/*******************************************************************************
 *
 ******************************************************************************/
void build_global_source(int n_quad_v, int n_derivs, double atran, double *atran_l,
                    double  **R, double  **T,
                    double  *F_p, double  *F_m, double  *S_p, double  *S_m,
                    double ***U, double ***W,
                    double **F_p_l, double **F_m_l, double **S_p_l, double **S_m_l,
                    uchar *derivs_h, uchar *derivs_p, work_data work) {

     build_global_source2(n_quad_v, n_derivs, atran, atran_l,
                          R, T, R, T, F_p, F_m, S_p, S_m,
                          U, W, U, W, F_p_l, F_m_l, S_p_l, S_m_l,
                          derivs_h, derivs_p, work);
}



void build_global_source2(int n_quad_v, int n_derivs, double atran, double *atran_l,
                     double  **R_p, double  **T_p, double  **R_m, double  **T_m,
                     double  *F_p, double  *F_m, double  *S_p, double  *S_m,
                     double ***R_p_l, double ***T_p_l, double ***R_m_l, double ***T_m_l,
                     double **F_p_l, double **F_m_l, double **S_p_l, double **S_m_l,
                     uchar *derivs_h, uchar *derivs_p, work_data work) {

     int i;

     double *A;

     double *v1;
     double *v2;

     v1 = get_work1(&work, WORK_DX);
     v2 = get_work1(&work, WORK_DX);

     dm_v_mul(T_p, F_p, n_quad_v, n_quad_v, v1);
     dvec_scale(atran, v1, v1, n_quad_v);
     dvec_sub(F_p, v1, v1, n_quad_v);
     dm_v_mul(R_m, F_m, n_quad_v, n_quad_v, v2);
     dvec_sub(v1, v2, S_p, n_quad_v);

     dm_v_mul(R_p, F_p, n_quad_v, n_quad_v, v1);
     dvec_sub(F_m, v1, v1, n_quad_v);
     dvec_scale(atran, v1, v1, n_quad_v);
     dm_v_mul(T_m, F_m, n_quad_v, n_quad_v, v2);
     dvec_sub(v1, v2, S_m, n_quad_v);

     if (flags_or(derivs_p, n_derivs)) {
          A = get_work1(&work, WORK_DX);

          for (i = 0; i < n_derivs; ++i) {
               if (derivs_h[i]) {
                    dvec_scale(atran, F_p_l[i], v1, n_quad_v);
                    dvec_scale(atran_l[i], F_p, v2, n_quad_v);
                    dvec_add(v1, v2, A, n_quad_v);

                    dm_v_mul(T_p_l[i], F_p, n_quad_v, n_quad_v, v1);
                    dvec_scale(atran, v1, v1, n_quad_v);
                    dvec_sub(F_p_l[i], v1, v1, n_quad_v);

                    dm_v_mul(T_p, A, n_quad_v, n_quad_v, v2);
                    dvec_sub(v1, v2, v1, n_quad_v);

                    dm_v_mul(R_m_l[i], F_m, n_quad_v, n_quad_v, v2);
                    dvec_sub(v1, v2, v1, n_quad_v);

                    dm_v_mul(R_m, F_m_l[i], n_quad_v, n_quad_v, v2);
                    dvec_sub(v1, v2, S_p_l[i], n_quad_v);


                    dvec_scale(atran, F_m_l[i], v1, n_quad_v);

                    dvec_scale(atran_l[i], F_m, v2, n_quad_v);
                    dvec_add(v1, v2, v1, n_quad_v);

                    dm_v_mul(T_m_l[i], F_m, n_quad_v, n_quad_v, v2);
                    dvec_sub(v1, v2, v1, n_quad_v);

                    dm_v_mul(T_m, F_m_l[i], n_quad_v, n_quad_v, v2);
                    dvec_sub(v1, v2, v1, n_quad_v);

                    dm_v_mul(R_p_l[i], F_p, n_quad_v, n_quad_v, v2);
                    dvec_scale(atran, v2, v2, n_quad_v);
                    dvec_sub(v1, v2, v1, n_quad_v);

                    dm_v_mul(R_p, A, n_quad_v, n_quad_v, v2);
                    dvec_sub(v1, v2, S_m_l[i], n_quad_v);
               }
               else
               if (derivs_p[i] && atran_l[i] != 0.) {
                    dvec_scale(atran, F_p_l[i], v1, n_quad_v);
                    dvec_scale(atran_l[i], F_p, v2, n_quad_v);
                    dvec_add(v1, v2, A, n_quad_v);

                    dm_v_mul(T_p, A, n_quad_v, n_quad_v, v2);
                    dvec_sub(F_p_l[i], v2, v1, n_quad_v);

                    dm_v_mul(R_m, F_m_l[i], n_quad_v, n_quad_v, v2);
                    dvec_sub(v1, v2, S_p_l[i], n_quad_v);


                    dvec_scale(atran, F_m_l[i], v1, n_quad_v);

                    dvec_scale(atran_l[i], F_m, v2, n_quad_v);
                    dvec_add(v1, v2, v1, n_quad_v);

                    dm_v_mul(T_m, F_m_l[i], n_quad_v, n_quad_v, v2);
                    dvec_sub(v1, v2, v1, n_quad_v);

                    dm_v_mul(R_p, A, n_quad_v, n_quad_v, v2);
                    dvec_sub(v1, v2, S_m_l[i], n_quad_v);
               }
               else
               if (derivs_p[i] && atran_l[i] == 0.) {
                    dvec_zero(S_p_l[i], n_quad_v);
                    dvec_zero(S_m_l[i], n_quad_v);
               }
          }
     }
#ifdef USE_AD_FOR_TL_CALC_BUILD_GLOBAL_SOURCE
     build_global_source_tl_with_ad2(n_quad_v, n_derivs, atran, atran_l,
                                     R_p, T_p, R_m, T_m,
                                     F_p, F_m, S_p, S_m,
                                     R_p_l, T_p_l, R_m_l, T_m_l,
                                     F_p_l, F_m_l, S_p_l, S_m_l,
                                     derivs_h, derivs_p, work);
#endif
}



/*******************************************************************************
 *
 ******************************************************************************/
void build_global_source_thermal(int n_quad_v, int n_derivs,
                       double atran, double *atran_l,
                       double **R_p, double **T_p, double **R_m, double **T_m,
                       double *F0_p, double *F0_m, double *F1_p, double *F1_m,
                       double *S_p, double *S_m,
                       double ***R_p_l, double ***T_p_l, double ***R_m_l, double ***T_m_l,
                       double **F0_p_l, double **F0_m_l, double **F1_p_l, double **F1_m_l,
                       double **S_p_l, double **S_m_l,
                       uchar *derivs_h, uchar *derivs_p, work_data work) {

     double *v1;
     double *v2;

     v1 = get_work1(&work, WORK_DX);
     v2 = get_work1(&work, WORK_DX);

     dm_v_mul(T_p, F1_p, n_quad_v, n_quad_v, v1);
     dvec_sub(F0_p, v1, v1, n_quad_v);
     dm_v_mul(R_m, F0_m, n_quad_v, n_quad_v, v2);
     dvec_sub(v1, v2, S_p, n_quad_v);

     dm_v_mul(R_p, F1_p, n_quad_v, n_quad_v, v1);
     dvec_sub(F1_m, v1, v1, n_quad_v);
     dm_v_mul(T_m, F0_m, n_quad_v, n_quad_v, v2);
     dvec_sub(v1, v2, S_m, n_quad_v);

     if (flags_or(derivs_p, n_derivs)) {

     }
}



/*******************************************************************************
 *
 ******************************************************************************/
void scale_source_vectors(int n_quad_v, int n_derivs,
                          double btran, double *btran_l,
                          double  *F_p1, double  *F_m1,
                          double  *F_p2, double  *F_m2,
                          double **F_p_l1, double **F_m_l1,
                          double **F_p_l2, double **F_m_l2,
                          uchar *derivs_h, uchar *derivs_p, work_data work) {

     int i;

     double *v1;
     double *v2;

     if (n_derivs > 0) {
          v1 = get_work1(&work, WORK_DX);
          v2 = get_work1(&work, WORK_DX);
     }

     for (i = 0; i < n_derivs; ++i) {
          if (derivs_p[i]) {
/*
          if (derivs_h[i]) {
*/
               dvec_scale(btran_l[i], F_p1, v1, n_quad_v);
               dvec_scale(btran, F_p_l1[i], v2, n_quad_v);
               dvec_add(v1, v2, F_p_l2[i], n_quad_v);

               dvec_scale(btran_l[i], F_m1, v1, n_quad_v);
               dvec_scale(btran, F_m_l1[i], v2, n_quad_v);
               dvec_add(v1, v2, F_m_l2[i], n_quad_v);
          }
/*
          else
          if (derivs_p[i]) {
               dvec_scale(btran_l[i], F_p1, F_p_l2[i], n_quad_v);
               dvec_scale(btran_l[i], F_m1, F_m_l2[i], n_quad_v);
          }
*/
     }

     dvec_scale(btran, F_p1, F_p2, n_quad_v);
     dvec_scale(btran, F_m1, F_m2, n_quad_v);
}
