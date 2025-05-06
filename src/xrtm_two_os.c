/*******************************************************************************
**
**    Copyright (C) 2007-2020 Greg McGarragh <greg.mcgarragh@colostate.edu>
**
**    This source code is licensed under the GNU General Public License (GPL),
**    Version 3.  See the file COPYING for more details.
**
*******************************************************************************/

#include <gutil.h>
#include <gmath_matrix.h>

#define NO_WORK_CAST

#include "xrtm.h"
#include "xrtm_derivs.h"
#include "xrtm_sos.h"
#include "xrtm_two_os.h"
#include "xrtm_utility.h"


/*******************************************************************************
 *
 ******************************************************************************/
int rtm_two_os(int i_four,
               int n_quad, int n_stokes, int n_derivs, int n_layers,
               double qf, double *qx_v, double *qw_v, double F_0, double mu_0,
               int n_ulevels, int *ulevels, double *utaus,
               int n_umus, double *umus_v,
               double *omega, double **omega_l, double *ltau, double **ltau_l,
               double *Rs_q0, double **Rs_q0_l, double **Rs_qq, double ***Rs_qq_l,
               double *btran, double **btran_l,
               double *as_0, double **as_0_l, double *atran, double **atran_l,
               double **P_q0_mm, double **P_q0_pm,
               double ***P_qq_pp, double ***P_qq_mp, double ***P_qq_mm, double ***P_qq_pm,
               double ***P_q0_mm_l, double ***P_q0_pm_l,
               double ****P_qq_pp_l, double ****P_qq_mp_l, double ****P_qq_mm_l, double ****P_qq_pm_l,
               double **I_p, double **I_m, double ***I_p_l, double ***I_m_l,
               int sfi, int surface, int utau_output, int vector, int flag,
               derivs_data *derivs, work_data work) {

     int i;
     int j;
     int k;
     int l;
     int m;

     int n_quad_v_x;
     int n_quad_v_d;

     int n_umus_v;

     int n_xmus_v;

     enum work_type work_dx;
     enum work_type work_dx2;
     enum work_type work_dxx;

     int offset;

     double a;
     double b;
     double c;
     double d;
     double e;

     double *xmus_v;

     double *v1;
     double *v2;
     double *v3;
     double *v4;
     double *v5;
     double *v6;
     double *v7;

     double **w1;
     double **w2;
     double **w3;
     double **w4;

     double *x;
     double *y;

     double **alpha;
     double **beta;

     double **r_0;
     double **t_0;
     double ***r_0_l;
     double ***t_0_l;
/*
     double ***r_p;
*/
     double ***t_p;
     double ***r_m;
/*
     double ***t_m;

     double ****r_p_l;
*/
     double ****t_p_l;
     double ****r_m_l;
/*
     double ****t_m_l;
*/
     double **Psi1;
     double **Psi2;
     double ***Psi3;

     double ***Psi1_l;
     double ***Psi2_l;
     double ****Psi3_l;

     double **Phi1;
     double ***Phi2;
     double ***Phi3;

     double ***Phi1_l;
     double ****Phi2_l;
     double ****Phi3_l;

     double ***R_p;
/*
     double ***R_m;
*/
     double ****R_p_l;
/*
     double ****R_m_l;
*/
     double **i_p[2];
     double **i_m[2];

     double ***i_p_l[2];
     double ***i_m_l[2];


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     calc_quad_and_umus(n_quad, n_umus, n_stokes, NULL, NULL, NULL, &n_quad_v_x, &n_quad_v_d, &n_umus_v, flag);

     if (! flag) {
          n_xmus_v   = n_quad_v_x; xmus_v = qx_v;

          offset   = 0;

          work_dx  = WORK_DX;
          work_dx2 = WORK_DX;
          work_dxx = WORK_DXX;
     }
     else {
          n_xmus_v   = n_umus_v; xmus_v = umus_v;

          offset   = n_quad_v_x;

          work_dx  = WORK_DU;
          work_dx2 = WORK_DD;
          work_dxx = WORK_DUX;
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     v1     = (double *)   get_work1(&work, WORK_DD);
     v2     = (double *)   get_work1(&work, WORK_DD);
     v3     = (double *)   get_work1(&work, WORK_DD);
     v4     = (double *)   get_work1(&work, WORK_DD);
     v5     = (double *)   get_work1(&work, WORK_DD);
     v6     = (double *)   get_work1(&work, WORK_DD);
     v7     = (double *)   get_work1(&work, WORK_DD);

     w1     = (double **)  get_work1(&work, work_dxx);
     w2     = (double **)  get_work1(&work, work_dxx);
     w3     = (double **)  get_work1(&work, work_dxx);
     w4     = (double **)  get_work1(&work, work_dxx);

     x      = (double *)   get_work1(&work, WORK_DX);
     y      = (double *)   get_work1(&work, work_dx2);

     alpha  = (double **)  get_work2(&work, WORK_DX,  WORK_LAYERS_V, NULL);
     beta   = (double **)  get_work2(&work, work_dx,  WORK_LAYERS_V, NULL);

     r_0    = (double **)  get_work2(&work, WORK_DX,  WORK_LAYERS_V, NULL);
     t_0    = (double **)  get_work2(&work, WORK_DX,  WORK_LAYERS_V, NULL);

     r_m    = (double ***) get_work2(&work, work_dxx, WORK_LAYERS_V, NULL);
/*
     r_p    = (double ***) get_work2(&work, work_dxx, WORK_LAYERS_V, NULL);
*/
     t_p    = (double ***) get_work2(&work, work_dxx, WORK_LAYERS_V, NULL);

     Psi1   = (double **)  get_work2(&work, work_dx2, WORK_LAYERS_V, NULL);
     Psi2   = (double **)  get_work2(&work, work_dx2, WORK_LAYERS_V, NULL);
     Psi3   = (double ***) get_work2(&work, work_dxx, WORK_LAYERS_V, NULL);
if (0)
     Phi1   = (double **)  get_work2(&work, WORK_DX,  WORK_LAYERS_V, NULL);
     Phi2   = (double ***) get_work2(&work, work_dxx, WORK_LAYERS_V, NULL);
     Phi3   = (double ***) get_work2(&work, work_dxx, WORK_LAYERS_V, NULL);
/*
     R_p    = (double ***) get_work2(&work, work_dxx, WORK_LAYERS_V, NULL);
     R_m    = (double ***) get_work2(&work, work_dxx, WORK_LAYERS_V, NULL);

     i_p[0] = (double **)  get_work2(&work, WORK_DX,  WORK_LAYERS_V, NULL);
     i_p[1] = (double **)  get_work2(&work, work_dx,  WORK_LAYERS_V, NULL);
     i_m[0] = (double **)  get_work2(&work, WORK_DX,  WORK_LAYERS_V, NULL);
     i_m[1] = (double **)  get_work2(&work, work_dx,  WORK_LAYERS_V, NULL);
*/
     R_p    = (double ***) get_work_d3(&work, n_layers + 1, n_xmus_v, n_quad_v_x);
/*
     R_m    = (double ***) get_work_d3(&work, n_layers + 1, n_xmus_v, n_quad_v_x);
*/
     i_p[0] = (double **)  get_work_d2(&work, n_layers + 1, n_quad_v_x);
     i_p[1] = (double **)  get_work_d2(&work, n_layers + 1, n_xmus_v);
     i_m[0] = (double **)  get_work_d2(&work, n_layers + 1, n_quad_v_x);
     i_m[1] = (double **)  get_work_d2(&work, n_layers + 1, n_xmus_v);
/*
     if (! vector) {
          t_m = t_p;
          r_m = r_p;
     }
     else {
          r_m = (double **)  get_work2(&work, WORK_DXX, WORK_LAYERS_V, NULL);
          t_m = (double **)  get_work2(&work, WORK_DXX, WORK_LAYERS_V, NULL);
     }
*/
     if (n_derivs > 0) {
          r_0_l    = (double ***)  get_work3(&work, WORK_DX,  WORK_BOTH_V, NULL);
          t_0_l    = (double ***)  get_work3(&work, WORK_DX,  WORK_BOTH_V, NULL);

          r_m_l    = (double ****) get_work3(&work, work_dxx, WORK_BOTH_V, NULL);
/*
          r_p_l    = (double ****) get_work3(&work, work_dxx, WORK_BOTH_V, NULL);
*/
          t_p_l    = (double ****) get_work3(&work, work_dxx, WORK_BOTH_V, NULL);

          Psi1_l = (double ***)  get_work3(&work, work_dx2, WORK_BOTH_V, NULL);
          Psi2_l = (double ***)  get_work3(&work, work_dx2, WORK_BOTH_V, NULL);
          Psi3_l = (double ****) get_work3(&work, work_dxx, WORK_BOTH_V, NULL);
if (0)
          Phi1_l = (double ***)  get_work3(&work, WORK_DX,  WORK_BOTH_V, NULL);
          Phi2_l = (double ****) get_work3(&work, work_dxx, WORK_BOTH_V, NULL);
          Phi3_l = (double ****) get_work3(&work, work_dxx, WORK_BOTH_V, NULL);
/*
          R_p_l    = (double ***)  get_work3(&work, work_dxx, WORK_BOTH_V, NULL);
          R_m_l    = (double ***)  get_work3(&work, work_dxx, WORK_BOTH_V, NULL);

          i_p_l[0] = (double ***)  get_work3(&work, WORK_DX,  WORK_BOTH_V, NULL);
          i_p_l[1] = (double ***)  get_work3(&work, work_dx,  WORK_BOTH_V, NULL);
          i_m_l[0] = (double ***)  get_work3(&work, WORK_DX,  WORK_BOTH_V, NULL);
          i_m_l[1] = (double ***)  get_work3(&work, work_dx,  WORK_BOTH_V, NULL);
*/
          R_p_l    = (double ****) get_work_d4(&work, n_layers + 1, n_derivs, n_xmus_v, n_quad_v_x);
/*
          R_m_l    = (double ****) get_work_d4(&work, n_layers + 1, n_derivs, n_xmus_v, n_quad_v_x);
*/
          i_p_l[0] = (double ***)  get_work_d3(&work, n_layers + 1, n_derivs, n_quad_v_x);
          i_p_l[1] = (double ***)  get_work_d3(&work, n_layers + 1, n_derivs, n_xmus_v);
          i_m_l[0] = (double ***)  get_work_d3(&work, n_layers + 1, n_derivs, n_quad_v_x);
          i_m_l[1] = (double ***)  get_work_d3(&work, n_layers + 1, n_derivs, n_xmus_v);
/*
          if (! vector) {
               t_m_l = t_p_l;
               r_m_l = r_p_l;
          }
          else {
               r_m_l = get_work3(&work, WORK_DXX, WORK_BOTH_V, NULL);
               t_m_l = get_work3(&work, WORK_DXX, WORK_BOTH_V, NULL);
          }
*/
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     init_array1_d(I_p[0], n_quad_v_d, 0.);
     init_array1_d(I_m[0], n_quad_v_d, 0.);

     if (n_ulevels == 2) {
          init_array1_d(I_p[1], n_quad_v_d, 0.);
          init_array1_d(I_m[1], n_quad_v_d, 0.);
     }
     if (n_derivs > 0) {
          init_array2_d(I_p_l[0], n_derivs, n_quad_v_d, 0.);
          init_array2_d(I_m_l[0], n_derivs, n_quad_v_d, 0.);

          if (n_ulevels == 2) {
               init_array2_d(I_p_l[1], n_derivs, n_quad_v_d, 0.);
               init_array2_d(I_m_l[1], n_derivs, n_quad_v_d, 0.);
          }
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     for (i = 0; i < n_layers; ++i) {
          for (k = 0; k < n_quad_v_x; ++k) {
               r_0[i][k] = omega[i] * P_q0_pm[i][k];
               t_0[i][k] = omega[i] * P_q0_mm[i][k];
          }

          for (j = 0; j < n_derivs; ++j) {
               for (k = 0; k < n_quad_v_x; ++k) {
                    r_0_l[i][j][k] = omega_l[i][j] * P_q0_pm[i][k];
                    t_0_l[i][j][k] = omega_l[i][j] * P_q0_mm[i][k];

                    if (derivs->layers[i][j]) {
                         r_0_l[i][j][k] += omega[i] * P_q0_pm_l[i][j][k];
                         t_0_l[i][j][k] += omega[i] * P_q0_mm_l[i][j][k];
                    }
               }
          }
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     for (i = 0; i < n_quad_v_x; ++i)
          x[i] = 1. / qx_v[i];

     for (i = 0; i < n_xmus_v; ++i)
          y[i] = 1. / xmus_v[i];

     for (i = 0; i < n_layers; ++i) {
          for (j = 0; j < n_quad_v_x; ++j)
               alpha[i][j] = exp(-ltau[i] * x[j]);

          for (j = 0; j < n_xmus_v; ++j)
               beta [i][j] = exp(-ltau[i] * y[j]);
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     for (i = 0; i < n_layers; ++i) {
          for (j = 0; j < n_quad_v_x; ++j) {
               Psi1[i][j] = alpha[i][j] * atran[i];
               Psi2[i][j] = alpha[i][j];
if (0)
               Phi1[i][j] = alpha[i][j] - atran[i];

               a = (x[j] + as_0[i]);
               b = x[j] * Psi2[i][j];

               for (k = 0; k < n_derivs; ++k) {
                    Psi1_l[i][k][j] = (-ltau_l[i][k] * a - ltau[i] * as_0_l[i][k]) * Psi1[i][j];
                    Psi2_l[i][k][j] =  -ltau_l[i][k] * b;
if (0)
                    Phi1_l[i][k][j] = Psi2_l[i][k][j] - (-ltau_l[i][k] * as_0[i] - ltau[i] * as_0_l[i][k]) * atran[i];
               }
          }
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     a = 1. / 4.;

     if (! surface) {
          init_array1_d(i_p[0][n_layers], n_quad_v_x, 0.);
          for (j = 0; j < n_derivs; ++j)
               init_array1_d(i_p_l[0][n_layers][j], n_quad_v_x, 0.);
     }
     else {
          dvec_copy(i_p[0][n_layers], Rs_q0, n_quad_v_x);
          for (j = 0; j < n_derivs; ++j) {
               if (! derivs->layers[n_layers][j])
                    init_array1_d(i_p_l[0][n_layers][j], n_quad_v_x, 0.);
               else
                    dvec_copy(i_p_l[0][n_layers][j], Rs_q0_l[j], n_quad_v_x);
          }
     }

     for (j = n_layers - 1; j >= 0; --j) {
          for (k = 0; k < n_quad_v_x; ++k) {
/*
               v1[k] = x[k] * as_0[j] / (x[k] + as_0[j]);
*/
               v1[k] = x[k] * 1. / mu_0 / (x[k] + as_0[j]);

               v2[k] = v1[k] * (1. - Psi1[j][k]);
               v3[k] = (1. - Psi1[j][k]) * r_0[j][k];
               v4[k] = v1[k] * r_0[j][k];

               i_p[0][j][k] = Psi1[j][k] * i_p[0][j+1][k] + a * v2[k] * r_0[j][k];
          }

          for (l = 0; l < n_derivs; ++l) {
               for (k = 0; k < n_quad_v_x; ++k) {
/*
                    b = (((x[k] - v1[k]) * as_0_l[j][l]) / (x[k] + as_0[j])) * v3[k] + v4[k] * -Psi1_l[j][l][k];
*/
                    b = (-v1[k] * as_0_l[j][l] / (x[k] + as_0[j])) * v3[k] + v4[k] * -Psi1_l[j][l][k];

                    if (derivs->layers[j][l])
                         b += v2[k] * r_0_l[j][l][k];

                    i_p_l[0][j][l][k] = Psi1_l[j][l][k] * i_p[0][j+1][k] + Psi1[j][k] * i_p_l[0][j+1][l][k] + a * b;
               }
          }
     }
/*
     dvec_copy(I_p[0], i_p[0][0], n_quad_v_x);
     for (j = 0; j < n_derivs; ++j)
          dvec_copy(I_p_l[0][j], i_p_l[0][0][j], n_quad_v_x);
*/

     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
if (0) {
     init_array1_d(i_m[0][0], n_quad_v_x, 0.);
     for (j = 0; j < n_derivs; ++j)
          init_array1_d(i_m_l[0][0][j], n_quad_v_x, 0.);

     for (j = 0; j < n_layers; ++j) {

          b = a * btran[j];

          for (k = 0; k < n_quad_v_x; ++k) {
/*
               v1[k] = x[k] * as_0[j] / (as_0[j] - x[k]);
*/
               v1[k] = x[k] * 1. / mu_0 / (as_0[j] - x[k]);

               v2[k] = b * v1[k];
               v3[k] = v2[k] * Phi1[j][k];
               v4[k] = Phi1[j][k] * t_0[j][k];
               v5[k] = v1[k] * v4[k];
               v6[k] = b * v4[k];
               v7[k] = v2[k] * t_0[j][k];

               i_m[0][j+1][k] = Psi2[j][k] * i_m[0][j][k] + v3[k] * t_0[j][k];
          }

          for (l = 0; l < n_derivs; ++l) {
               c = a * btran_l[j][l];

               for (k = 0; k < n_quad_v_x; ++k) {
/*
                    i_m_l[0][j+1][l][k] = Psi2_l[j][l][k] * i_m[0][j][k] + Psi2[j][k] * i_m_l[0][j][l][k] +
                                        c * v5[k] + v6[k] * ((x[k] * as_0_l[j][l] - v1[k] * as_0_l[j][l]) / (as_0[j] - x[k])) + v7[k] * Phi1_l[j][l][k];
*/
                    i_m_l[0][j+1][l][k] = Psi2_l[j][l][k] * i_m[0][j][k] + Psi2[j][k] * i_m_l[0][j][l][k] +
                                        c * v5[k] + v6[k] * (-v1[k] * as_0_l[j][l] / (as_0[j] - x[k])) + v7[k] * Phi1_l[j][l][k];

                    if (derivs->layers[j][l])
                         i_m_l[0][j+1][l][k] += v3[k] * t_0_l[j][l][k];
               }
          }
     }
/*
     dvec_copy(I_m[1], i_m[0][n_layers], n_quad_v_x);
     for (j = 0; j < n_derivs; ++j)
          dvec_copy(I_m_l[1][j], i_m_l[0][n_layers][j], n_quad_v_x);
*/
}

     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (! vector) {
          P_qq_mm = P_qq_pp;
          P_qq_pm = P_qq_mp;

          if (n_derivs > 0) {
               P_qq_mm_l = P_qq_pp_l;
               P_qq_pm_l = P_qq_mp_l;
          }
     }

     a = (1. + (i_four == 0 ? 1. : 0.)) / 4.;

     for (i = 0; i < n_layers; ++i) {

          dvec_scale(a * omega[i], qw_v, v1, n_quad_v_x);

          dmat_mul_diag(P_qq_pm[i], v1, r_m[i], n_xmus_v, n_quad_v_x);
/*
          dmat_mul_diag(P_qq_mp[i], v1, r_p[i], n_xmus_v, n_quad_v_x);
*/
          dmat_mul_diag(P_qq_pp[i], v1, t_p[i], n_xmus_v, n_quad_v_x);

          for (j = 0; j < n_derivs; ++j) {
               dvec_scale(a * omega_l[i][j], qw_v, v2, n_quad_v_x);

               dmat_mul_diag(P_qq_pm[i], v2, r_m_l[i][j], n_xmus_v, n_quad_v_x);
               if (derivs->layers[i][j]) {
                    dmat_mul_diag(P_qq_pm_l[i][j], v1, w1, n_xmus_v, n_quad_v_x);
                    dmat_add(r_m_l[i][j], w1, r_m_l[i][j], n_xmus_v, n_quad_v_x);
               }
/*
               dmat_mul_diag(P_qq_mp[i], v2, r_p_l[i][j], n_xmus_v, n_quad_v_x);
               if (derivs->layers[i][j]) {
                    dmat_mul_diag(P_qq_mp_l[i][j], v1, w1, n_xmus_v, n_quad_v_x);
                    dmat_add(r_p_l[i][j], w1, r_p_l[i][j], n_xmus_v, n_quad_v_x);
               }
*/
               dmat_mul_diag(P_qq_pp[i], v2, t_p_l[i][j], n_xmus_v, n_quad_v_x);
               if (derivs->layers[i][j]) {
                    dmat_mul_diag(P_qq_pp_l[i][j], v1, w1, n_xmus_v, n_quad_v_x);
                    dmat_add(t_p_l[i][j], w1, t_p_l[i][j], n_xmus_v, n_quad_v_x);
               }
          }
/*
          if (vector) {
               dmat_mul_diag(P_qq_pm[i], v1, r_m[i], n_xmus_v, n_quad_v_x);

               dmat_mul_diag(P_qq_mm[i], v1, t_m[i], n_xmus_v, n_quad_v_x);

               for (j = 0; j < n_derivs; ++j) {
                    dvec_scale(a * omega_l[i][j], qw_v, v2, n_quad_v_x);

                    dmat_mul_diag(P_qq_pm[i], v2, r_m_l[i][j], n_xmus_v, n_quad_v_x);
                    if (derivs->layers[i][j]) {
                         dmat_mul_diag(P_qq_pm_l[i][j], v1, w1, n_xmus_v, n_quad_v_x);
                         dmat_add(r_m_l[i][j], w1, r_m_l[i][j], n_xmus_v, n_quad_v_x);
                    }

                    dmat_mul_diag(P_qq_mm[i], v2, t_m_l[i][j], n_xmus_v, n_quad_v_x);
                    if (derivs->layers[i][j]) {
                         dmat_mul_diag(P_qq_mm_l[i][j], v1, w1, n_xmus_v, n_quad_v_x);
                         dmat_add(t_m_l[i][j], w1, t_m_l[i][j], n_xmus_v, n_quad_v_x);
                    }
               }
          }
*/
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     for (i = 0; i < n_layers; ++i) {
          for (j = 0; j < n_xmus_v; ++j) {
               Psi1[i][j] = beta[i][j] * atran[i];
               Psi2[i][j] = beta[i][j];

               for (k = 0; k < n_quad_v_x; ++k)
                    Psi3[i][j][k] = beta[i][j] * alpha[i][k];

               a = (y[j] + as_0[i]);
               b = y[j] * Psi2[i][j];

               for (k = 0; k < n_derivs; ++k) {
                    Psi1_l[i][k][j] = (-ltau_l[i][k] * a - ltau[i] * as_0_l[i][k]) * Psi1[i][j];
                    Psi2_l[i][k][j] =  -ltau_l[i][k] * b;

                    for (l = 0; l < n_quad_v_x; ++l)
                         Psi3_l[i][k][j][l] = -ltau_l[i][k] * (y[j] + x[l]) * Psi3[i][j][l];
               }
          }
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (! surface) {
          dmat_zero(R_p[n_layers], n_xmus_v, n_quad_v_x);
          for (j = 0; j < n_derivs; ++j) {
               dmat_zero(R_p_l[n_layers][j], n_xmus_v, n_quad_v_x);
          }
     }
     else {
          dmat_copy(R_p[n_layers], Rs_qq, n_xmus_v, n_quad_v_x);
          dmat_mul_diag(R_p[n_layers], x, R_p[n_layers], n_xmus_v, n_quad_v_x);

          for (j = 0; j < n_derivs; ++j) {
               if (! derivs->layers[n_layers][j]) {
                    dmat_zero(R_p_l[n_layers][j], n_xmus_v, n_quad_v_x);
               }
               else {
                    dmat_copy(R_p_l[n_layers][j], Rs_qq_l[j], n_xmus_v, n_quad_v_x);
                    dmat_mul_diag(R_p_l[n_layers][j], x, R_p_l[n_layers][j], n_xmus_v, n_quad_v_x);
               }
          }
     }

     for (j = n_layers - 1; j >= 0; --j) {
          for (k = 0; k < n_xmus_v; ++k) {
               for (l = 0; l < n_quad_v_x; ++l) {
                    w1[k][l] = y[k] * x[l] / (y[k] + x[l]);

                    R_p[j][k][l] = R_p[j+1][k][l] * Psi3[j][k][l] + w1[k][l] * (1. - Psi3[j][k][l]) * r_m[j][k][l];
               }
          }

          for (m = 0; m < n_derivs; ++m) {
               for (k = 0; k < n_xmus_v; ++k) {
                    for (l = 0; l < n_quad_v_x; ++l) {
                         R_p_l[j][m][k][l] = R_p_l[j+1][m][k][l] * Psi3[j][k][l] + R_p[j+1][k][l] * Psi3_l[j][m][k][l] +
                                           w1[k][l] * (-Psi3_l[j][m][k][l] * r_m[j][k][l] + (1. - Psi3[j][k][l]) * r_m_l[j][m][k][l]);
                    }
               }
          }
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     i = 1;
/*
     if (! surface) {
*/
          init_array1_d(i_p[i][n_layers], n_xmus_v, 0.);
          for (j = 0; j < n_derivs; ++j)
               init_array1_d(i_p_l[i][n_layers][j], n_xmus_v, 0.);
/*
     }
     else {
          dm_v_mul(Rs_qq, i_m[1-i][n_layers], n_xmus_v, n_quad_v_x, i_p[i][n_layers]);

          for (j = 0; j < n_derivs; ++j) {
               dm_v_mul(Rs_qq, i_m_l[1-i][n_layers][j], n_xmus_v, n_quad_v_x, i_p_l[i][n_layers][j]);

               if (derivs->layers[n_layers][j]) {
                    dm_v_mul(Rs_qq_l[j], i_m[1-i][n_layers], n_xmus_v, n_quad_v_x, v1);
                    dvec_add(i_p_l[i][n_layers][j], v1, i_p_l[i][n_layers][j], n_quad_v_x);
               }
          }
     }
*/


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     for (j = 0; j < n_layers; ++j) {
          a = atran[j];

          for (k = 0; k < n_xmus_v; ++k) {
               b = beta[j][k] * a;
               e = a * beta[j][k];

               for (l = 0; l < n_quad_v_x; ++l) {
                    c = alpha[j][l] * a;
                    d = beta[j][k] * alpha[j][l];

                    if (y[k] == x[l])
                         Phi2[j][k][l] = ltau[j] * b;
                    else
                         Phi2[j][k][l] = (b - c) / (x[l] - y[k]);
                    if (as_0[j] == x[l])
                         Phi3[j][k][l] = ltau[j] * d;
                    else
                         Phi3[j][k][l] = (d - e) / (as_0[j] - x[l]);

                    for (m = 0; m < n_derivs; ++m) {
                         if (y[k] == x[l])
                              Phi2_l[j][m][k][l] = (ltau_l[j][m] + ltau[j] * ((-ltau_l[j][m] * (y[k] + as_0[j]) - as_0_l[j][m]))) * b;
                         else
                              Phi2_l[j][m][k][l] = ((-ltau_l[j][m] * (y[k] + as_0[j]) - ltau[j] * as_0_l[j][m]) * b - (-ltau_l[j][m] * (x[l] + as_0[j]) - ltau[j] * as_0_l[j][m]) * c) / (x[l] - y[k]);
                         if (as_0[j] == x[l])
                              Phi3_l[j][m][k][l] = (ltau_l[j][m] + ltau[j] * -ltau[j] * (x[l] + y[k])) * d;
                         else
                              Phi3_l[j][m][k][l] = ((-ltau_l[j][m] * (x[l] + y[k]) * d - (-ltau_l[j][m] * (as_0[j] + y[k]) - ltau[j] * as_0_l[j][m]) * e) - Phi3[j][k][l] * as_0_l[j][m]) / (as_0[j] - x[l]);
                    }
               }
          }
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     i = 1;

     for (j = n_layers - 1; j >= 0; --j) {

          dm_v_diag_mul(Psi1[j], i_p[i][j+1], i_p[i][j], n_xmus_v);

          for (k = 0; k < n_xmus_v; ++k) {
               for (l = 0; l < n_quad_v_x; ++l) {
/*
                    w1[k][l] = x[l] * as_0[j] / (4. * (x[l] + as_0[j]) * (y[k] + as_0[j]));
*/
                    w1[k][l] = x[l] * 1. / mu_0 / (4. * (x[l] + as_0[j]) * (y[k] + as_0[j]));
                    w2[k][l] = 1. - Psi1[j][k] - (y[k] + as_0[j]) * Phi2[j][k][l];
                    w3[k][l] = r_0[j][l] * w1[k][l];
                    w4[k][l] = w1 [k][l] * w2[k][l];
               }

               a = 0.;
               for (l = 0; l < n_quad_v_x; ++l)
                    a += t_p[j][k][l] * (Phi2[j][k][l] * i_p[1-i][j+1][l] + r_0[j][l] * w4[k][l]);

               i_p[i][j][k] += y[k] * a;
          }

          for (m = 0; m < n_derivs; ++m) {
               dm_v_diag_mul(Psi1_l[j][m], i_p[i][j+1], v1, n_xmus_v);
               dm_v_diag_mul(Psi1[j], i_p_l  [i][j+1][m], v2, n_xmus_v);
               dvec_add(v1, v2, i_p_l[i][j][m], n_xmus_v);

               for (k = 0; k < n_xmus_v; ++k) {
                    a = 0.;

                    for (l = 0; l < n_quad_v_x; ++l) {
/*
                         b = Phi2_l[j][m][k][l] * i_p[1-i][j+1][l] + Phi2[j][k][l] * i_p_l[1-i][j+1][m][l] +
                             r_0[j][l] * ((x[l] * as_0_l[j][m] - w1[k][l] * 4. * (as_0_l[j][m] * (y[k] + as_0[j]) + (x[l] + as_0[j]) * as_0_l[j][m])) / (4. * (x[l] + as_0[j]) * (y[k] + as_0[j]))) * w2[k][l] + w3[k][l] * (- Psi1_l[j][m][k] - as_0_l[j][m] * Phi2[j][k][l]  - (y[k] + as_0[j]) * Phi2_l[j][m][k][l]);
*/
                         b = Phi2_l[j][m][k][l] * i_p[1-i][j+1][l] + Phi2[j][k][l] * i_p_l[1-i][j+1][m][l] +
                             r_0[j][l] * (-w1[k][l] * (4. * (as_0_l[j][m] * (y[k] + as_0[j]) + (x[l] + as_0[j]) * as_0_l[j][m])) / (4. * (x[l] + as_0[j]) * (y[k] + as_0[j]))) * w2[k][l] + w3[k][l] * (- Psi1_l[j][m][k] - as_0_l[j][m] * Phi2[j][k][l]  - (y[k] + as_0[j]) * Phi2_l[j][m][k][l]);

                         if (derivs->layers[j][m])
                              b += r_0_l[j][m][l] * w4[k][l];

                         a += t_p_l[j][m][k][l] * (Phi2[j][k][l] * i_p[1-i][j+1][l] + r_0[j][l] * w4[k][l]) + t_p[j][k][l] * b;

                    }

                    i_p_l[i][j][m][k] += y[k] * a;
               }
          }

          for (k = 0; k < n_xmus_v; ++k) {
               for (l = 0; l < n_quad_v_x; ++l) {
                    w1[k][l] = x[l] * y[k] / ((y[k] + as_0[j]) * (x[l] + y[k]));
                    w2[k][l] = 1. - Psi1[j][k] - (y[k] + as_0[j]) * Phi3[j][k][l];
                    w3[k][l] = r_m[j][k][l] * w1[k][l];
                    w4[k][l] = w1[k][l] * w2[k][l];
               }

               v1[k] = 0.;
               for (l = 0; l < n_quad_v_x; ++l)
                    v1[k] += (Phi3[j][k][l] * R_p[j+1][k][l] + r_m[j][k][l] * w4[k][l]) * t_0[j][l];

               i_p[i][j][k] += as_0[j] / 4. * v1[k];
          }

          for (m = 0; m < n_derivs; ++m) {
               for (k = 0; k < n_xmus_v; ++k) {
                    a = 0.;
                    for (l = 0; l < n_quad_v_x; ++l) {
                         a += (Phi3_l[j][m][k][l] * R_p[j+1][k][l] + Phi3[j][k][l] * R_p_l[j+1][m][k][l] +
                               r_m_l[j][m][k][l] * w4[k][l] + r_m[j][k][l] * w1[k][l] * as_0_l[j][m] * (x[l] + y[k]) / ((y[k] + as_0[j]) * (x[l] + y[k])) * w2[k][l] + w3[k][l] * (- Psi1_l[j][m][k] - as_0_l[j][m] * Phi3[j][k][l] - (y[k] + as_0[j]) * Phi3_l[j][m][k][l])) * t_0[j][l];

                         if (derivs->layers[j][m])
                              a += (Phi3[j][k][l] * R_p[j+1][k][l] + r_m[j][k][l] * w4[k][l]) * t_0_l[j][m][l];
                    }

                    i_p_l[i][j][m][k] += (as_0_l[j][m] * v1[k] + as_0[j] * a) / 4.;
               }
          }
     }

     dvec_copy(I_p[0]+offset, i_p[i][0], n_xmus_v);
     for (j = 0; j < n_derivs; ++j)
          dvec_copy(I_p_l[0][j]+offset, i_p_l[i][0][j], n_xmus_v);
/*
     dvec_add(I_p[0]+offset, i_p[i][0], I_p[0]+n_quad_v_x, n_xmus_v);
     for (j = 0; j < n_derivs; ++j)
          dvec_add(I_p_l[0][j]+offset, i_p_l[i][0][j], I_p_l[0][j]+n_quad_v_x, n_xmus_v);
*/

     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     a = F_0 * mu_0 / PI;

     for (j = offset; j < offset + n_xmus_v; ++j) {
          I_p[0][j] *= a;
/*
          I_m[1][j] *= a;
*/
     }

     for (k = 0; k < n_derivs; ++k) {
          for (j = offset; j < offset + n_xmus_v; ++j) {
               I_p_l[0][k][j] *= a;
/*
               I_m_l[1][k][j] *= a;
*/
          }
     }


     return 0;
}
