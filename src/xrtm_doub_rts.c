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
#include "xrtm_doub_rts.h"
#include "xrtm_doubling.h"
#include "xrtm_save_tree.h"
#include "xrtm_scatter.h"
#include "xrtm_source.h"
#include "xrtm_utility.h"


/*******************************************************************************
 *
 ******************************************************************************/
static void calc_doub_params(int n_derivs, double ltau, double *ltau_l, int *n_doub,
                             double *d_tau, double *d_tau_l, uchar *derivs) {

     int i;

     *d_tau = MIN(ltau, *d_tau);

     *n_doub = (int) (log(ltau / *d_tau) / log(2.) + 1.);

     *d_tau  = ltau / pow(2., *n_doub);

     for (i = 0; i < n_derivs; ++i) {
          if (! derivs[i])
               continue;

          d_tau_l[i] = ltau_l[i] / pow(2., *n_doub);
     }
}



/*******************************************************************************
 *
 ******************************************************************************/
void rtm_doub_rts(int n_quad, int n_stokes, int n_derivs, double F_0,
                  double *qx_v, double *qw_v,
                  double planck0, double planck1, double *planck0_l, double *planck1_l,
                  double omega, double *omega_l, double ltau, double *ltau_l,
                  double as_0, double *as_0_l,
                  double *P_q0_mm, double *P_q0_pm,
                  double **r_p, double **t_p, double **r_m, double **t_m,
                  double **R_p, double **T_p, double **R_m, double **T_m,
                  double *S_p, double *S_m, double *Sl_p, double *Sl_m,
                  double **P_q0_mm_l, double **P_q0_pm_l,
                  double ***r_p_l, double ***t_p_l, double ***r_m_l, double ***t_m_l,
                  double ***R_p_l, double ***T_p_l, double ***R_m_l, double ***T_m_l,
                  double **S_p_l, double **S_m_l, double **Sl_p_l, double **Sl_m_l,
                  double d_tau, int symmetric, int solar, int thermal, int vector,
                  uchar *derivs_layers, uchar *derivs_beam, uchar *derivs_thermal,
                  save_tree_data save_tree, work_data work) {

     int i;
     int j;

     int n_quad_v;

     int n_doub;

     double a;
     double b;
     double c;
     double d;

     double d_tran;

     double lin_fac;

     double *v1;

     double *d_tau_l;

     double *d_path_l;
     double *d_tran_l;

     double *lin_fac_l;

     double **w1;

     double *L_p;
     double *L_m;

     double **L_p_l;
     double **L_m_l;


     n_quad_v = n_quad * n_stokes;


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     v1        = get_work1(&work, WORK_DX);

     d_tau_l   = get_work1(&work, WORK_DDERIVS);

     d_path_l  = get_work1(&work, WORK_DDERIVS);
     d_tran_l  = get_work1(&work, WORK_DDERIVS);

     lin_fac_l = get_work1(&work, WORK_DDERIVS);

     if (thermal) {
          L_p = get_work1(&work, WORK_DX);
          L_m = get_work1(&work, WORK_DX);
     }

     if (flags_or(derivs_layers, n_derivs))
          w1 = get_work1(&work, WORK_DXX);

     if (thermal && flags_or(derivs_thermal, n_derivs)) {
          L_p_l = get_work2(&work, WORK_DX, WORK_DERIVS_V, NULL);
          L_m_l = get_work2(&work, WORK_DX, WORK_DERIVS_V, NULL);
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     calc_doub_params(n_derivs, ltau, ltau_l, &n_doub, &d_tau, d_tau_l, derivs_layers);


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     dmat_scale(-d_tau, r_p, R_p, n_quad_v, n_quad_v);

     dmat_scale(-d_tau, t_p, T_p, n_quad_v, n_quad_v);
     dmat_i_sub(T_p, T_p, n_quad_v);

     if (vector)
          dmat_mul_D_A(n_quad, n_stokes, R_p, R_p);

     if (flags_or(derivs_layers, n_derivs)) {
          for (i = 0; i < n_derivs; ++i) {
               if (! derivs_layers[i])
                    continue;

               dmat_scale(-d_tau, r_p_l[i], R_p_l[i] , n_quad_v, n_quad_v);
               dmat_scale(-d_tau_l[i], r_p, w1, n_quad_v, n_quad_v);
               dmat_add(R_p_l[i], w1, R_p_l[i], n_quad_v, n_quad_v);

               dmat_scale( d_tau, t_p_l[i], T_p_l[i] , n_quad_v, n_quad_v);
               dmat_scale( d_tau_l[i], t_p, w1, n_quad_v, n_quad_v);
               dmat_add(T_p_l[i], w1, T_p_l[i], n_quad_v, n_quad_v);

               if (vector)
                    dmat_mul_D_A(n_quad, n_stokes, R_p_l[i], R_p_l[i]);
          }
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (solar) {
          d_tran = exp(-d_tau * as_0);

          if (flags_or(derivs_beam, n_derivs)) {
               for (i = 0; i < n_derivs; ++i) {
                    if (derivs_layers[i]) {
                         d_path_l[i] = -d_tau_l[i] * as_0 - d_tau * as_0_l[i];
                         d_tran_l[i] = d_path_l[i] * d_tran;
                    }
                    else
                    if (derivs_beam[i]) {
                         d_path_l[i] =                    - d_tau * as_0_l[i];
                         d_tran_l[i] = d_path_l[i] * d_tran;
                    }
               }
          }

          a = F_0 / (4. * PI);

          b = a * omega * d_tran * d_tau;

          dm_v_dinv_mul(qx_v, P_q0_pm, S_p, n_quad_v);
          dvec_scale(b, S_p, S_p, n_quad_v);

          dm_v_dinv_mul(qx_v, P_q0_mm, S_m, n_quad_v);
          dvec_scale(b, S_m, S_m, n_quad_v);

          if (vector)
               dm_v_mul_D_A(n_quad, n_stokes, S_m, S_m);
     }

     if (thermal) {
          if (planck0 == 0.)
               lin_fac = 0.;
          else
               lin_fac = (planck1 / planck0 - 1.) / pow(2., n_doub);

          if (flags_or(derivs_thermal, n_derivs)) {
               for (i = 0; i < n_derivs; ++i) {
                    if (derivs_thermal[i]) {
                         if (planck0 == 0.)
                              lin_fac_l[i] = 0.;
                         else
                              lin_fac_l[i] = (planck1_l[i] - planck1 * planck0_l[i] / planck0) / planck0 / pow(2., n_doub);
                    }
               }
          }

          c = (1. - omega) * planck0 * d_tau;
          for (i = 0; i < n_quad_v; ++i)
               v1[i] = c;
          dm_v_dinv_mul(qx_v, v1, v1, n_quad_v);
          dvec_copy(Sl_p, v1, n_quad_v);
          dvec_copy(Sl_m, v1, n_quad_v);
     }
if (solar) {
     if (flags_or(derivs_beam, n_derivs)) {
          for (i = 0; i < n_derivs; ++i) {
               if (derivs_layers[i]) {
                    c = a * (omega_l[i] * d_tran      * d_tau +
                             omega      * d_tran_l[i] * d_tau +
                             omega      * d_tran      * d_tau_l[i]);

                    dvec_scale(c, P_q0_pm, S_p_l[i], n_quad_v);
                    dvec_scale(b, P_q0_pm_l[i], v1, n_quad_v);
                    dvec_add(S_p_l[i], v1, S_p_l[i], n_quad_v);
                    dm_v_dinv_mul(qx_v, S_p_l[i], S_p_l[i], n_quad_v);

                    dvec_scale(c, P_q0_mm, S_m_l[i], n_quad_v);
                    dvec_scale(b, P_q0_mm_l[i], v1, n_quad_v);
                    dvec_add(S_m_l[i], v1, S_m_l[i], n_quad_v);
                    dm_v_dinv_mul(qx_v, S_m_l[i], S_m_l[i], n_quad_v);

                    if (vector)
                         dm_v_mul_D_A(n_quad, n_stokes, S_m_l[i], S_m_l[i]);
               }
               else
               if (derivs_beam[i]) {
                    dvec_zero(S_p_l[i], n_quad_v);
                    dvec_zero(S_m_l[i], n_quad_v);
               }
          }
     }
}
if (thermal) {
     if (flags_or(derivs_thermal, n_derivs)) {
          for (i = 0; i < n_derivs; ++i) {
               if (derivs_layers[i]) {
                    d =     - omega_l[i] * planck0      * d_tau +
                        (1. - omega)     * planck0_l[i] * d_tau +
                        (1. - omega)     * planck0      * d_tau_l[i];
                    for (j = 0; j < n_quad_v; ++j)
                         v1[j] = d;
                    dm_v_dinv_mul(qx_v, v1, v1, n_quad_v);
                    dvec_copy(Sl_p_l[i], v1, n_quad_v);
                    dvec_copy(Sl_m_l[i], v1, n_quad_v);
               }
               else
               if (derivs_thermal[i]) {
                    d = (1. - omega)     * planck0_l[i] * d_tau;
                    for (j = 0; j < n_quad_v; ++j)
                         v1[j] = d;
                    dm_v_dinv_mul(qx_v, v1, v1, n_quad_v);
                    dvec_copy(Sl_p_l[i], v1, n_quad_v);
                    dvec_copy(Sl_m_l[i], v1, n_quad_v);
               }
          }
     }
}

     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     c = 1.;
     for (i = 0; i < n_doub; ++i) {
          if (! symmetric) {
               if (! flags_or(derivs_layers, n_derivs) && ! (solar && flags_or(derivs_beam, n_derivs)) && ! (thermal && flags_or(derivs_thermal, n_derivs)))
                    layer_double  (R_p, T_p, S_m, S_p, Sl_m, Sl_p, L_m, L_p, n_quad_v, d_tran, lin_fac, ! DOUB_RTS_USE_PARTICULAR_SOLUTION && solar, ! DOUB_RTS_USE_PARTICULAR_SOLUTION && thermal, i == 0, work);
               else
                    layer_double_l(R_p, T_p, S_m, S_p, Sl_m, Sl_p, L_m, L_p, R_p_l, T_p_l, S_m_l, S_p_l, Sl_m_l, Sl_p_l, L_m_l, L_p_l, n_quad_v, n_derivs, d_tran, d_tran_l, lin_fac, lin_fac_l, ! DOUB_RTS_USE_PARTICULAR_SOLUTION && solar, ! DOUB_RTS_USE_PARTICULAR_SOLUTION && thermal, i == 0, derivs_layers, derivs_beam, derivs_thermal, work);
          }
          else {
               if (! flags_or(derivs_layers, n_derivs) && ! (solar && flags_or(derivs_beam, n_derivs)) && ! (thermal && flags_or(derivs_thermal, n_derivs)))
                    layer_double_s  (R_p, T_p, S_m, S_p, Sl_m, Sl_p, L_m, L_p, n_quad_v, d_tran, lin_fac, ! DOUB_RTS_USE_PARTICULAR_SOLUTION && solar, ! DOUB_RTS_USE_PARTICULAR_SOLUTION && thermal, i == 0, work);
               else
                    layer_double_s_l(R_p, T_p, S_m, S_p, Sl_m, Sl_p, L_m, L_p, R_p_l, T_p_l, S_m_l, S_p_l, Sl_m_l, Sl_p_l, L_m_l, L_p_l, n_quad_v, n_derivs, d_tran, d_tran_l, lin_fac, lin_fac_l, ! DOUB_RTS_USE_PARTICULAR_SOLUTION && solar, ! DOUB_RTS_USE_PARTICULAR_SOLUTION && thermal, i == 0, derivs_layers, derivs_beam, derivs_thermal, work);
          }

          if (solar) {
               d_tran *= d_tran;

               if (flags_or(derivs_beam, n_derivs)) {
                    c *= 2.;

                    for (j = 0; j < n_derivs; ++j) {
                         if (! derivs_beam[j])
                              continue;

                         d_tran_l[j] = d_path_l[j] * d_tran * c;
                    }
               }
          }

          if (thermal) {
               lin_fac *= 2.;

               if (flags_or(derivs_thermal, n_derivs)) {
                    for (j = 0; j < n_derivs; ++j) {
                         if (! derivs_thermal[j])
                              continue;

                         lin_fac_l[j] *= 2.;
                    }
               }
          }
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     dmat_mul_D_A(n_quad, n_stokes, R_p, R_p);
     phase_matrix_symmetry2(n_quad, n_stokes, T_p, R_p, T_m, R_m, 1.);

     dm_v_mul_D_A(n_quad, n_stokes, S_m, S_m);

     if (flags_or(derivs_layers, n_derivs)) {
          for (i = 0; i < n_derivs; ++i) {
               if (derivs_layers[i]) {
                    dmat_mul_D_A(n_quad, n_stokes, R_p_l[i], R_p_l[i]);
                    phase_matrix_symmetry2(n_quad, n_stokes, T_p_l[i], R_p_l[i], T_m_l[i], R_m_l[i], 1.);
               }
          }
     }

     if (solar && flags_or(derivs_beam, n_derivs)) {
          for (i = 0; i < n_derivs; ++i) {
               if (derivs_beam[i])
                    dm_v_mul_D_A(n_quad, n_stokes, S_m_l[i], S_m_l[i]);
          }
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (DOUB_RTS_USE_PARTICULAR_SOLUTION) {
          double atran;
          double *atran_l;

          double *Fs_p;
          double *Fs_m;

          double *Ft0_p;
          double *Ft0_m;
          double *Ft1_p;
          double *Ft1_m;

          double **Fs_p_l;
          double **Fs_m_l;

          double **Ft0_p_l;
          double **Ft0_m_l;
          double **Ft1_p_l;
          double **Ft1_m_l;

          if (solar) {
               Fs_p = get_work1(&work, WORK_DX);
               Fs_m = get_work1(&work, WORK_DX);
          }

          if (thermal) {
               Ft0_p = get_work1(&work, WORK_DX);
               Ft0_m = get_work1(&work, WORK_DX);
               Ft1_p = get_work1(&work, WORK_DX);
               Ft1_m = get_work1(&work, WORK_DX);
          }

          if (solar) {
               if (flags_or(derivs_beam, n_derivs)) {
                    atran_l = get_work1(&work, WORK_DDERIVS);

                    Fs_p_l   = get_work2(&work, WORK_DX, WORK_DERIVS_V, derivs_beam);
                    Fs_m_l   = get_work2(&work, WORK_DX, WORK_DERIVS_V, derivs_beam);
               }
          }

          if (thermal) {
               if (flags_or(derivs_thermal, n_derivs)) {
                    Ft0_p_l = get_work2(&work, WORK_DX, WORK_DERIVS_V, derivs_thermal);
                    Ft0_m_l = get_work2(&work, WORK_DX, WORK_DERIVS_V, derivs_thermal);
                    Ft1_p_l = get_work2(&work, WORK_DX, WORK_DERIVS_V, derivs_thermal);
                    Ft1_m_l = get_work2(&work, WORK_DX, WORK_DERIVS_V, derivs_thermal);
               }
          }

          if (! solar) {
/*
               dvec_zero(S_p, n_quad_v);
               dvec_zero(S_m, n_quad_v);
*/
          }
          else {
               if (! CLASSICAL_PARTICULAR_SOLUTION_USE_2D) {
                    double **tpr;
                    double **tmr;
                    double **gamma;
                    double ***tpr_l;
                    double ***tmr_l;
                    double ***gamma_l;

                    tpr   = get_work1(&work, WORK_DXX);
                    tmr   = get_work1(&work, WORK_DXX);
                    gamma = get_work1(&work, WORK_DXX);

                    if (flags_or(derivs_layers, n_derivs)) {
                         tpr_l   = get_work2(&work, WORK_DXX, WORK_DERIVS_V, derivs_layers);
                         tmr_l   = get_work2(&work, WORK_DXX, WORK_DERIVS_V, derivs_layers);
                         gamma_l = get_work2(&work, WORK_DXX, WORK_DERIVS_V, derivs_layers);
                    }

                    build_txr(n_quad, n_stokes, n_derivs, r_p, t_p, tpr, tmr, r_p_l, t_p_l, tpr_l, tmr_l, derivs_layers, work);

                    build_gamma(n_quad_v, n_derivs, tpr, tmr, gamma, tpr_l, tmr_l, gamma_l, derivs_layers, work);

                    build_source_vectors_solar_classic_1n(n_quad, n_stokes, n_derivs, qx_v, F_0, omega, omega_l, as_0, as_0_l, P_q0_mm, P_q0_pm, tpr, tmr, gamma, Fs_p, Fs_m, P_q0_mm_l, P_q0_pm_l, tpr_l, tmr_l, gamma_l, Fs_p_l, Fs_m_l, derivs_layers, derivs_beam, save_tree, work);
               }
               else
                    build_source_vectors_solar_classic_2n(n_quad_v, n_derivs, qx_v, F_0, omega, omega_l, as_0, as_0_l, P_q0_mm, P_q0_pm, r_p, t_p, Fs_p, Fs_m, P_q0_mm_l, P_q0_pm_l, r_p_l, t_p_l, Fs_p_l, Fs_m_l, derivs_layers, derivs_beam, work);

               atran = exp(-ltau * as_0);

               if (flags_or(derivs_beam, n_derivs)) {
                    for (i = 0; i < n_derivs; ++i) {
                         if (! derivs_beam[i])
                              continue;

                         atran_l[i] = -(ltau_l[i] * as_0 + ltau * as_0_l[i]) * atran;
                    }
               }

               build_global_source_solar(n_quad_v, n_derivs, atran, atran_l, R_p, T_p, Fs_p, Fs_m, S_p, S_m, R_p_l, T_p_l, Fs_p_l, Fs_m_l, S_p_l, S_m_l, derivs_layers, derivs_beam, work, NULL, NULL, NULL, NULL);

          }

          if (! thermal) {
/*
               dvec_zero(Sl_p, n_quad_v);
               dvec_zero(Sl_m, n_quad_v);
*/
          }
          else {
               build_source_vectors_thermal2(n_quad, n_stokes, n_derivs, qx_v, planck0, planck1, planck0_l, planck1_l, omega, omega_l, ltau, ltau_l, r_p, t_p, r_m, t_m, Ft0_p, Ft0_m, Ft1_p, Ft1_m, r_p_l, t_p_l, r_m_l, t_m_l, Ft0_p_l, Ft0_m_l, Ft1_p_l, Ft1_m_l, derivs_layers, derivs_thermal, work);

               build_global_source_thermal  (n_quad_v, n_derivs, R_p, T_p, R_m, T_m, Ft0_p, Ft0_m, Ft1_p, Ft1_m, Sl_p, Sl_m, R_p_l, T_p_l, R_m_l, T_m_l, Ft0_p_l, Ft0_m_l, Ft1_p_l, Ft1_m_l, Sl_p_l, Sl_m_l, derivs_layers, derivs_thermal, work);
          }
     }
}



/*******************************************************************************
 *
 ******************************************************************************/
#ifdef INCLUDE_DEV_SOURCE
#include "../src2/xrtm_doub_rts2.c"
#endif
