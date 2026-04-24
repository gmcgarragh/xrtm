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

#include "xrtm.h"
#include "xrtm_derivs.h"
#include "xrtm_sos.h"


/*******************************************************************************
 *
 ******************************************************************************/
int rtm_sos(int i_four,
            int n_quad, int n_stokes, int n_derivs, int n_layers,
            double qf, double *qx_v, double *qw_v, double F_0, double mu_0,
            int n_ulevels, int *ulevels, double *utaus,
            double *omega, double **omega_l, double *ltau, double **ltau_l,
            double *btran, double **btran_l,
            double *as_0, double **as_0_l, double *atran, double **atran_l,
            double **P_q0_mm, double **P_q0_pm,
            double ***P_qq_pp, double ***P_qq_mp, double ***P_qq_mm, double ***P_qq_pm,
            double ***P_q0_mm_l, double ***P_q0_pm_l,
            double ****P_qq_pp_l, double ****P_qq_mp_l, double ****P_qq_mm_l, double ****P_qq_pm_l,
            double **Rs_qq, double ***Rs_qq_l,
            double *I1_m, double **I1_m_l, double *In_p, double **In_p_l,
            double **I_p, double **I_m, double ***I_p_l, double ***I_m_l,
            int max_os, double max_tau, double sos_tol,
            int sfi, int surface, int utau_output, int vector,
            derivs_data *derivs, work_data work) {

     int i;
     int ii;
     int j;
     int jj;
     int k;
     int l;
     int n;

     int n_quad_v;

     int n_layers2;

     int i_last_ulevel;

     int i_os;

     int *index;

     double a;
     double b;
     double b_l;
     double c;

     double *aa_l;

     double *v1;
     double *v2;
     double *v3;
     double *v4;
     double *v5;
     double *v6;

     double **w1;

     double *ltau2;
     double *btran2;

     double **ltau_l2;
     double **btran_l2;

     double **T;
     double **E;

     double ***W;
     double ***F;

     double ***r_p;
     double ***t_p;
     double ***r_m;
     double ***t_m;

     double ****t_p_l;
     double ****r_p_l;
     double ****t_m_l;
     double ****r_m_l;

     double ***i_p;
     double ***i_m;

     double ****i_p_l;
     double ****i_m_l;


     n_quad_v = n_quad * n_stokes;


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (n_derivs > 0)
          aa_l = get_work1(&work, WORK_DDERIVS);


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     n = 0;
     for (i = 0; i < n_layers; ++i) {
         a = ltau[i] / max_tau;

         n++;
         if (a > 1.)
              n += (int) a;
     }

     n_layers2 = n;


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     index  = get_work_i1(&work, n_layers2);
     ltau2  = get_work_d1(&work, n_layers2);
     btran2 = get_work_d1(&work, n_layers2 + 1);

     if (n_derivs > 0) {
          ltau_l2  = get_work_d2(&work, n_layers2,     n_derivs);
          btran_l2 = get_work_d2(&work, n_layers2 + 1, n_derivs);
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     ii = 0;
     for (i = 0; i < n_layers; ++i) {
         a = ltau[i] / max_tau;

         n = 1;
         if (a > 1.)
              n += (int) a;

         a = ltau[i] / n;
         for (j = 0; j < n_derivs; ++j)
              aa_l[j] = ltau_l[i][j] / n;

         for (j = 0; j < n; ++j) {
              index [ii] = i;
              ltau2 [ii] = a;
              b = exp(j * -a * as_0[i]);
              btran2[ii] = btran[i] * b;
              for (k = 0; k < n_derivs; ++k) {
                   ltau_l2 [ii][k] = aa_l[k];
                   btran_l2[ii][k] = (btran_l[i][k] + btran[i] * j * (-aa_l[k] * as_0[i] - a * as_0_l[i][k])) * b;
              }
              ii++;
         }
     }

     btran2[ii] = btran[n_layers];
     for (i = 0; i < n_derivs; ++i)
          btran_l2[ii][i] = btran_l[n_layers][i];


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     v1  = get_work1(&work, WORK_DX);
     v2  = get_work1(&work, WORK_DX);
     v3  = get_work1(&work, WORK_DX);
     v4  = get_work1(&work, WORK_DX);

     w1  = get_work1(&work, WORK_DXX);

     T   = get_work_d2(&work, n_layers2, n_quad_v);
     E   = get_work_d2(&work, n_layers2, n_quad_v);

     i_p = get_work_d3(&work, 2, n_layers2 + 1, n_quad_v);
     i_m = get_work_d3(&work, 2, n_layers2 + 1, n_quad_v);

     if (n_derivs > 0) {
          v5  = get_work1(&work, WORK_DX);
          v6  = get_work1(&work, WORK_DX);

          W   = get_work_d3(&work, n_layers2, n_derivs, n_quad_v);
          F   = get_work_d3(&work, n_layers2, n_derivs, n_quad_v);

          i_p_l = get_work_d4(&work, 2, n_layers2 + 1, n_derivs, n_quad_v);
          i_m_l = get_work_d4(&work, 2, n_layers2 + 1, n_derivs, n_quad_v);
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     for (i = 0; i < n_layers2; ++i) {
          for (j = 0; j < n_quad_v; ++j) {
               T[i][j] = exp(-ltau2[i] / qx_v[j]);
               E[i][j] = 1. - T[i][j];
          }

          for (j = 0; j < n_derivs; ++j) {
               for (k = 0; k < n_quad_v; ++k) {
                    W[i][j][k] = -ltau_l2[i][j] / qx_v[k] * T[i][k];
                    F[i][j][k] = - W[i][j][k];
               }
          }
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
i_last_ulevel = n_ulevels - 1;

if (ulevels[0] == 0)
     init_array1_d(I_m[0], n_quad_v, 0.);
if (ulevels[i_last_ulevel] == n_layers)
     init_array1_d(I_p[i_last_ulevel], n_quad_v, 0.);

     if (n_derivs > 0) {
if (ulevels[0            ] == 0)
          init_array2_d(I_m_l[0], n_derivs, n_quad_v, 0.);
if (ulevels[i_last_ulevel] == n_layers)
          init_array2_d(I_p_l[i_last_ulevel], n_derivs, n_quad_v, 0.);
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     a = F_0 / (4. * PI) / 2.;

     dvec_copy(i_p[0][n_layers2], In_p, n_quad_v);
     if (n_derivs > 0)
          dmat_copy(i_p_l[0][n_layers2], In_p_l, n_derivs, n_quad_v);

     for (j = n_layers2 - 1; j >= 0; --j) {
          jj = index[j];

          b = a * omega[jj] * (btran2[j] + btran2[j+1]);

          for (k = 0; k < n_quad_v; ++k)
               i_p[0][j][k] = T[j][k] * i_p[0][j+1][k] + b * E[j][k] * P_q0_pm[jj][k];

          for (k = 0; k < n_derivs; ++k) {
               b_l = a * (omega_l[jj][k] * (btran2[j] + btran2[j+1]) + omega[jj] * (btran_l2[j][k] + btran_l2[j+1][k]));

               for (l = 0; l < n_quad_v; ++l) {
                    i_p_l[0][j][k][l] = W[j][k][l] * i_p[0][j+1][l] + T[j][l] * i_p_l[0][j+1][k][l] + (b_l * E[j][l] + b * F[j][k][l]) * P_q0_pm[jj][l];

                    if (derivs->layers[jj][k])
                         i_p_l[0][j][k][l] += b * E[j][l] * P_q0_pm_l[jj][k][l];
               }
          }
     }
if (ulevels[0] == 0) {
     dvec_copy(I_p[0], i_p[0][0], n_quad_v);
     if (n_derivs > 0)
          dmat_copy(I_p_l[0], i_p_l[0][0], n_derivs, n_quad_v);
}

     dvec_copy(i_m[0][0], I1_m, n_quad_v);
     if (n_derivs > 0)
          dmat_copy(i_m_l[0][0], I1_m_l, n_derivs, n_quad_v);

     for (j = 0; j < n_layers2; ++j) {
          jj = index[j];

          b = a * omega[jj] * (btran2[j] + btran2[j+1]);

          for (k = 0; k < n_quad_v; ++k)
               i_m[0][j+1][k] = T[j][k] * i_m[0][j][k] + b * E[j][k] * P_q0_mm[jj][k];

          for (k = 0; k < n_derivs; ++k) {
               b_l = a * (omega_l[jj][k] * (btran2[j] + btran2[j+1]) + omega[jj] * (btran_l2[j][k] + btran_l2[j+1][k]));

               for (l = 0; l < n_quad_v; ++l) {
                    i_m_l[0][j+1][k][l] = W[j][k][l] * i_m[0][j][l] + T[j][l] * i_m_l[0][j][k][l] + (b_l * E[j][l] + b * F[j][k][l]) * P_q0_mm[jj][l];

                    if (derivs->layers[jj][k])
                         i_m_l[0][j+1][k][l] += b * E[j][l] * P_q0_mm_l[jj][k][l];
               }
          }
     }
if (ulevels[i_last_ulevel] == n_layers) {
     dvec_copy(I_m[i_last_ulevel], i_m[0][n_layers2], n_quad_v);
     if (n_derivs > 0)
          dmat_copy(I_m_l[i_last_ulevel], i_m_l[0][n_layers2], n_derivs, n_quad_v);
}

     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     r_p = get_work_d3(&work, n_layers2, n_quad_v, n_quad_v);
     t_p = get_work_d3(&work, n_layers2, n_quad_v, n_quad_v);

     if (! vector) {
          r_m = r_p;
          t_m = t_p;
     }
     else {
          r_m = get_work_d3(&work, n_layers2, n_quad_v, n_quad_v);
          t_m = get_work_d3(&work, n_layers2, n_quad_v, n_quad_v);
     }

     if (n_derivs > 0) {
          r_p_l = get_work_d4(&work, n_layers2, n_derivs, n_quad_v, n_quad_v);
          t_p_l = get_work_d4(&work, n_layers2, n_derivs, n_quad_v, n_quad_v);

          if (! vector) {
               r_m_l = r_p_l;
               t_m_l = t_p_l;
          }
          else {
               r_m_l = get_work_d4(&work, n_layers2, n_derivs, n_quad_v, n_quad_v);
               t_m_l = get_work_d4(&work, n_layers2, n_derivs, n_quad_v, n_quad_v);
          }
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     /* Includes division by 2. from os > 1 loops. */
     a = (1. + (i_four == 0 ? 1. : 0.)) / 4. / 2.;

     for (i = 0; i < n_layers2; ++i) {
          ii = index[i];

          b = a * omega[ii];

          dvec_scale(b, qw_v, v1, n_quad_v);

          dmat_mul_diag(P_qq_mp[ii], v1, r_p[i], n_quad_v, n_quad_v);

          dmat_mul_diag(P_qq_pp[ii], v1, t_p[i], n_quad_v, n_quad_v);

          for (j = 0; j < n_derivs; ++j) {
               b_l = a * omega_l[ii][j];

               dvec_scale(b_l, qw_v, v2, n_quad_v);

               dmat_mul_diag(P_qq_mp[ii], v2, r_p_l[i][j], n_quad_v, n_quad_v);
               if (derivs->layers[ii][j]) {
                    dmat_mul_diag(P_qq_mp_l[ii][j], v1, w1, n_quad_v, n_quad_v);
                    dmat_add(r_p_l[i][j], w1, r_p_l[i][j], n_quad_v, n_quad_v);
               }

               dmat_mul_diag(P_qq_pp[ii], v2, t_p_l[i][j], n_quad_v, n_quad_v);
               if (derivs->layers[ii][j]) {
                    dmat_mul_diag(P_qq_pp_l[ii][j], v1, w1, n_quad_v, n_quad_v);
                    dmat_add(t_p_l[i][j], w1, t_p_l[i][j], n_quad_v, n_quad_v);
               }
          }

          if (vector) {
               dmat_mul_diag(P_qq_pm[ii], v1, r_m[i], n_quad_v, n_quad_v);

               dmat_mul_diag(P_qq_mm[ii], v1, t_m[i], n_quad_v, n_quad_v);

               for (j = 0; j < n_derivs; ++j) {
                    b_l = a * omega_l[ii][j];

                    dvec_scale(b_l, qw_v, v2, n_quad_v);

                    dmat_mul_diag(P_qq_pm[ii], v2, r_m_l[i][j], n_quad_v, n_quad_v);
                    if (derivs->layers[ii][j]) {
                         dmat_mul_diag(P_qq_pm_l[ii][j], v1, w1, n_quad_v, n_quad_v);
                         dmat_add(r_m_l[i][j], w1, r_m_l[i][j], n_quad_v, n_quad_v);
                    }

                    dmat_mul_diag(P_qq_mm[ii], v2, t_m_l[i][j], n_quad_v, n_quad_v);
                    if (derivs->layers[ii][j]) {
                         dmat_mul_diag(P_qq_mm_l[ii][j], v1, w1, n_quad_v, n_quad_v);
                         dmat_add(t_m_l[i][j], w1, t_m_l[i][j], n_quad_v, n_quad_v);
                    }
               }
          }
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     a = DBL_MAX;

     i = 1;
     for (i_os = 1; i_os < max_os; ++i_os) {
          if (! surface) {
               init_array1_d(i_p[i][n_layers2], n_quad_v, 0.);
               if (n_derivs > 0)
                    init_array2_d(i_p_l[i][n_layers2], n_derivs, n_quad_v, 0.);
          }
          else {
               dm_v_mul(Rs_qq, i_m[1-i][n_layers2], n_quad_v, n_quad_v, i_p[i][n_layers2]);
               for (j = 0; j < n_derivs; ++j) {
                    dm_v_mul(Rs_qq, i_m_l[1-i][n_layers2][j], n_quad_v, n_quad_v, i_p_l[i][n_layers2][j]);

                    if (derivs->layers[n_layers][j]) {
                         dm_v_mul(Rs_qq_l[j], i_m[1-i][n_layers2], n_quad_v, n_quad_v, v1);
                         dvec_add(i_p_l[i][n_layers2][j], v1, i_p_l[i][n_layers2][j], n_quad_v);
                    }
               }
          }

          for (j = n_layers2 - 1; j >= 0; --j) {
               dm_v_diag_mul(T[j], i_p[i][j+1], i_p[i][j], n_quad_v);

               dvec_add(i_p[1-i][j], i_p[1-i][j+1], v1, n_quad_v);
               dm_v_mul(t_p[j], v1, n_quad_v, n_quad_v, v3);

               dvec_add(i_m[1-i][j], i_m[1-i][j+1], v2, n_quad_v);
               dm_v_mul(r_m[j], v2, n_quad_v, n_quad_v, v4);

               dvec_add(v3, v4, v3, n_quad_v);
               dm_v_diag_mul(E[j], v3, v4, n_quad_v);
               dvec_add(i_p[i][j], v4, i_p[i][j], n_quad_v);

               for (k = 0; k < n_derivs; ++k) {
                    dm_v_diag_mul(W[j][k], i_p[i][j+1], v4, n_quad_v);
                    dm_v_diag_mul(T[j], i_p_l[i][j+1][k], v5, n_quad_v);
                    dvec_add(v4, v5, i_p_l[i][j][k], n_quad_v);

                    dm_v_mul(t_p_l[j][k], v1, n_quad_v, n_quad_v, v4);
                    dvec_add(i_p_l[1-i][j][k], i_p_l[1-i][j+1][k], v5, n_quad_v);
                    dm_v_mul(t_p[j], v5, n_quad_v, n_quad_v, v6);
                    dvec_add(v4, v6, v5, n_quad_v);

                    dm_v_mul(r_m_l[j][k], v2, n_quad_v, n_quad_v, v4);
                    dvec_add(v5, v4, v5, n_quad_v);

                    dvec_add(i_m_l[1-i][j][k], i_m_l[1-i][j+1][k], v4, n_quad_v);
                    dm_v_mul(r_m[j], v4, n_quad_v, n_quad_v, v6);
                    dvec_add(v5, v6, v5, n_quad_v);

                    dm_v_diag_mul(F[j][k], v3, v4, n_quad_v);
                    dvec_add(i_p_l[i][j][k], v4, i_p_l[i][j][k], n_quad_v);
                    dm_v_diag_mul(E[j], v5, v4, n_quad_v);
                    dvec_add(i_p_l[i][j][k], v4, i_p_l[i][j][k], n_quad_v);
               }
          }

          b = i_p[i][0][0];

          jj = n_stokes;
          for (j = 1; j < n_quad; ++j) {
               if (i_p[i][0][jj] > b)
                    b = i_p[i][0][jj];
               jj += n_stokes;
          }
if (ulevels[0] == 0) {
/*
          dvec_copy(I_p[0], i_p[i][0], n_quad_v);
          if (n_derivs > 0)
               dmat_copy(I_p_l[0], i_p_l[i][0], n_derivs, n_quad_v);
*/
          dvec_add(I_p[0], i_p[i][0], I_p[0], n_quad_v);
          if (n_derivs > 0)
               dmat_add(I_p_l[0], i_p_l[i][0], I_p_l[0], n_derivs, n_quad_v);
}
          dvec_copy(i_m[i][0], I1_m, n_quad_v);
          if (n_derivs > 0)
               dmat_copy(i_m_l[i][0], I1_m_l, n_derivs, n_quad_v);

          for (j = 0; j < n_layers2; ++j) {
               dm_v_diag_mul(T[j], i_m[i][j], i_m[i][j+1], n_quad_v);

               dvec_add(i_p[1-i][j], i_p[1-i][j+1], v1, n_quad_v);
               dm_v_mul(r_p[j], v1, n_quad_v, n_quad_v, v3);

               dvec_add(i_m[1-i][j], i_m[1-i][j+1], v2, n_quad_v);
               dm_v_mul(t_m[j], v2, n_quad_v, n_quad_v, v4);

               dvec_add(v3, v4, v3, n_quad_v);
               dm_v_diag_mul(E[j], v3, v4, n_quad_v);
               dvec_add(i_m[i][j+1], v4, i_m[i][j+1], n_quad_v);

               for (k = 0; k < n_derivs; ++k) {
                    dm_v_diag_mul(W[j][k], i_m[i][j], v4, n_quad_v);
                    dm_v_diag_mul(T[j], i_m_l[i][j][k], v5, n_quad_v);
                    dvec_add(v4, v5, i_m_l[i][j+1][k], n_quad_v);

                    dm_v_mul(r_p_l[j][k], v1, n_quad_v, n_quad_v, v4);
                    dvec_add(i_p_l[1-i][j][k], i_p_l[1-i][j+1][k], v5, n_quad_v);
                    dm_v_mul(r_p[j], v5, n_quad_v, n_quad_v, v6);
                    dvec_add(v4, v6, v5, n_quad_v);

                    dm_v_mul(t_m_l[j][k], v2, n_quad_v, n_quad_v, v4);
                    dvec_add(v5, v4, v5, n_quad_v);

                    dvec_add(i_m_l[1-i][j][k], i_m_l[1-i][j+1][k], v4, n_quad_v);
                    dm_v_mul(t_m[j], v4, n_quad_v, n_quad_v, v6);
                    dvec_add(v5, v6, v5, n_quad_v);

                    dm_v_diag_mul(F[j][k], v3, v4, n_quad_v);
                    dvec_add(i_m_l[i][j+1][k], v4, i_m_l[i][j+1][k], n_quad_v);
                    dm_v_diag_mul(E[j], v5, v4, n_quad_v);
                    dvec_add(i_m_l[i][j+1][k], v4, i_m_l[i][j+1][k], n_quad_v);
               }
          }

          c = i_m[i][n_layers2][0];

          jj = n_stokes;
          for (j = 1; j < n_quad; ++j) {
               if (i_m[i][n_layers2][jj] > c)
                    c = i_m[i][n_layers2][jj];
               jj += n_stokes;
          }

          a = MIN(a, MAX(b, c));
if (ulevels[i_last_ulevel] == n_layers) {
/*
          dvec_copy(I_m[i_last_ulevel], i_m[i][n_layers2], n_quad_v);
          if (n_derivs > 0)
               dmat_copy(I_m_l[i_last_ulevel], i_m_l[i][n_layers2], n_derivs, n_quad_v);
*/
          dvec_add(I_m[i_last_ulevel], i_m[i][n_layers2], I_m[i_last_ulevel], n_quad_v);
          if (n_derivs > 0)
               dmat_add(I_m_l[i_last_ulevel], i_m_l[i][n_layers2], I_m_l[i_last_ulevel], n_derivs, n_quad_v);
}
          if (a < sos_tol)
               break;

          i = 1 - i;
     }
/*
     if (i == max_os) {
          fprintf(stderr, "ERROR: order of scattering not converging on %e\n", 1.e-4);
          return -1;
     }
*/

     return 0;
}
