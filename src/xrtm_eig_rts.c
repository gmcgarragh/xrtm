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
#include "xrtm_eig_rts.h"
#include "xrtm_eig_rts_a.h"
#include "xrtm_eig_util.h"
#include "xrtm_matrix.h"
#include "xrtm_save_tree.h"
#include "xrtm_scatter.h"
#include "xrtm_source.h"
#include "xrtm_utility.h"


/*******************************************************************************
 *
 ******************************************************************************/
/*
#define REAL
*/
#define TYPE		double
#define TYPE_PREFIX	d
#define TYPE_POSTFIX	d
#define WORK_XX		WORK_DX
#define WORK_XXX	WORK_DXX
#define XABS		fabs
#define XEXP		exp
#define XREAL(x)	(x)

#define CALC_GLOBAL_R_AND_T				calc_global_r_and_t
#define CALC_GLOBAL_R_AND_T_TL_WITH_AD			calc_global_r_and_t_tl_with_ad
#define SAVE_TREE_STRING				"calc_global_r_and_t"
#define FORWARD_SAVE_CALC_GLOBAL_R_AND_T_DATA		forward_save_calc_global_r_and_t_data
#define FORWARD_SAVE_CALC_GLOBAL_R_AND_T_ALLOC		forward_save_calc_global_r_and_t_alloc

#include "type_set.h"

#include "xrtm_eig_rts_x.c"

#include "type_unset.h"
/*
#undef REAL
*/
#undef TYPE
#undef TYPE_PREFIX
#undef TYPE_POSTFIX
#undef WORK_XX
#undef WORK_XXX
#undef XABS
#undef XEXP
#undef XREAL

#undef CALC_GLOBAL_R_AND_T
#undef CALC_GLOBAL_R_AND_T_TL_WITH_AD
#undef SAVE_TREE_STRING
#undef FORWARD_SAVE_CALC_GLOBAL_R_AND_T_DATA
#undef FORWARD_SAVE_CALC_GLOBAL_R_AND_T_ALLOC



/*******************************************************************************
 *
 ******************************************************************************/
#define TYPE		dcomplex
#define TYPE_PREFIX	z
#define TYPE_POSTFIX	dc
#define WORK_XX		WORK_ZX
#define WORK_XXX	WORK_ZXX
#define XABS		cabs
#define XEXP		cexp
#define XREAL(x)	(creal(x))

#define CALC_GLOBAL_R_AND_T				calc_global_r_and_t2
#define CALC_GLOBAL_R_AND_T_TL_WITH_AD			calc_global_r_and_t_tl_with_ad2
#define SAVE_TREE_STRING				"calc_global_r_and_t2"
#define FORWARD_SAVE_CALC_GLOBAL_R_AND_T_DATA		forward_save_calc_global_r_and_t_data2
#define FORWARD_SAVE_CALC_GLOBAL_R_AND_T_ALLOC		forward_save_calc_global_r_and_t_alloc2

#include "type_set.h"

#include "xrtm_eig_rts_x.c"

#include "type_unset.h"

#undef TYPE
#undef TYPE_PREFIX
#undef TYPE_POSTFIX
#undef WORK_XX
#undef WORK_XXX
#undef XABS
#undef XEXP
#undef XREAL

#undef CALC_GLOBAL_R_AND_T
#undef CALC_GLOBAL_R_AND_T_TL_WITH_AD
#undef SAVE_TREE_STRING
#undef FORWARD_SAVE_CALC_GLOBAL_R_AND_T_DATA
#undef FORWARD_SAVE_CALC_GLOBAL_R_AND_T_ALLOC



/*******************************************************************************
 *
 ******************************************************************************/
void rtm_eig_rts(int n_quad, int n_stokes, int n_derivs, double F_0,
                 double *qx_v, double *qw_v,
                 double planck0, double planck1, double *planck0_l, double *planck1_l,
                 double omega, double *omega_l, double ltau, double *ltau_l,
                 double as_0, double *as_0_l, double atran, double *atran_l,
                 double  *P_x0_p, double  *P_x0_m,
                 double  **r_p, double  **t_p, double  **r_m, double  **t_m,
                 double  **R_p, double  **T_p, double  **R_m, double  **T_m,
                 double  *S_p, double  *S_m, double  *Sl_p, double  *Sl_m,
                 double **P_x0_p_l, double **P_x0_m_l,
                 double ***r_p_l, double ***t_p_l, double ***r_m_l, double ***t_m_l,
                 double ***R_p_l, double ***T_p_l, double ***R_m_l, double ***T_m_l,
                 double **S_p_l, double **S_m_l, double **Sl_p_l, double **Sl_m_l,
                 int greens, int symmetric, int solar, int thermal, int vector,
                 int eigen_solver_real, int eigen_solver_complex,
                 uchar *derivs_layers, uchar *derivs_beam, uchar *derivs_thermal,
                 save_tree_data save_tree, work_data work) {

     int i;

     int n_quad_v;

     double **tmr;
     double **tpr;
     double **gamma;
/*
     double **aux = NULL;
*/
     double ***tmr_l;
     double ***tpr_l;
     double ***gamma_l;
/*
     double ***aul = NULL;
*/
     double *Fs_p;
     double *Fs_m;

     double *Fs1_p;
     double *Fs1_m;

     double *Ft0_p;
     double *Ft0_m;
     double *Ft1_p;
     double *Ft1_m;

     double **Fs_p_l;
     double **Fs_m_l;

     double **Fs1_p_l;
     double **Fs1_m_l;

     double **Ft0_p_l;
     double **Ft0_m_l;
     double **Ft1_p_l;
     double **Ft1_m_l;

     double *nu;
     double **nu_l;

     double **X_p;
     double **X_m;

     double ***X_p_l;
     double ***X_m_l;

     dcomplex *nu_c;
     dcomplex **nu_l_c;

     dcomplex **X_p_c;
     dcomplex **X_m_c;

     dcomplex ***X_p_l_c;
     dcomplex ***X_m_l_c;

     forward_save_rtm_eig_rts_data *save;


     n_quad_v = n_quad * n_stokes;


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (vector)
          symmetric = 0;


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (save_tree.t) {
          save_tree_encode_s(&save_tree, "rtm_eig_rts");

          if (save_tree_retrieve_data(&save_tree, forward_save_rtm_eig_rts_data, &save))
               forward_save_rtm_eig_rts_alloc(save, n_quad_v, vector);
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     tmr   = get_work1(&work, WORK_DXX);
     tpr   = get_work1(&work, WORK_DXX);
     gamma = get_work1(&work, WORK_DXX);

     if (solar) {
          Fs_p = get_work1(&work, WORK_DX);
          Fs_m = get_work1(&work, WORK_DX);

          if (greens) {
               Fs1_p = get_work1(&work, WORK_DX);
               Fs1_m = get_work1(&work, WORK_DX);
          }
     }

     if (thermal) {
          Ft0_p = get_work1(&work, WORK_DX);
          Ft0_m = get_work1(&work, WORK_DX);
          Ft1_p = get_work1(&work, WORK_DX);
          Ft1_m = get_work1(&work, WORK_DX);

          if (greens) {

          }
     }

     if (solar) {
          if (flags_or(derivs_beam, n_derivs)) {
               Fs_p_l = get_work2(&work, WORK_DX, WORK_DERIVS_V, derivs_beam);
               Fs_m_l = get_work2(&work, WORK_DX, WORK_DERIVS_V, derivs_beam);

               if (greens) {
                    Fs1_p_l = get_work2(&work, WORK_DX, WORK_DERIVS_V, derivs_beam);
                    Fs1_m_l = get_work2(&work, WORK_DX, WORK_DERIVS_V, derivs_beam);
               }
          }
     }

     if (thermal) {
          if (flags_or(derivs_thermal, n_derivs)) {
               Ft0_p_l = get_work2(&work, WORK_DX, WORK_DERIVS_V, derivs_thermal);
               Ft0_m_l = get_work2(&work, WORK_DX, WORK_DERIVS_V, derivs_thermal);
               Ft1_p_l = get_work2(&work, WORK_DX, WORK_DERIVS_V, derivs_thermal);
               Ft1_m_l = get_work2(&work, WORK_DX, WORK_DERIVS_V, derivs_thermal);

               if (greens) {

               }
          }
     }

     if (flags_or(derivs_layers, n_derivs)) {
          tmr_l   = get_work2(&work, WORK_DXX, WORK_DERIVS_V, derivs_layers);
          tpr_l   = get_work2(&work, WORK_DXX, WORK_DERIVS_V, derivs_layers);
          gamma_l = get_work2(&work, WORK_DXX, WORK_DERIVS_V, derivs_layers);
     }

     if (! vector) {
          nu  = get_work1(&work, WORK_DX);

          X_p = get_work1(&work, WORK_DXX);
          X_m = get_work1(&work, WORK_DXX);

          if (flags_or(derivs_layers, n_derivs)) {
               nu_l  = get_work2(&work, WORK_DX, WORK_DERIVS_V, derivs_layers);

               X_p_l = get_work2(&work, WORK_DXX, WORK_DERIVS_V, derivs_layers);
               X_m_l = get_work2(&work, WORK_DXX, WORK_DERIVS_V, derivs_layers);
          }
     }
     else {
          nu_c  = get_work1(&work, WORK_ZX);

          X_p_c = get_work1(&work, WORK_ZXX);
          X_m_c = get_work1(&work, WORK_ZXX);

          if (flags_or(derivs_layers, n_derivs)) {
               nu_l_c  = get_work2(&work, WORK_ZX, WORK_DERIVS_V, derivs_layers);

               X_p_l_c = get_work2(&work, WORK_ZXX, WORK_DERIVS_V, derivs_layers);
               X_m_l_c = get_work2(&work, WORK_ZXX, WORK_DERIVS_V, derivs_layers);
          }
     }

     if (symmetric) {
/*
          aux = get_work1(&work, WORK_DXX);

          if (flags_or(derivs_layers, n_derivs))
               aul = get_work2(&work, WORK_DXX, WORK_DERIVS_V, derivs_layers);
*/
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     build_txr(n_quad, n_stokes, n_derivs, r_p, t_p, tpr, tmr, r_p_l, t_p_l, tpr_l, tmr_l, derivs_layers, work);

     build_gamma(n_quad_v, n_derivs, tpr, tmr, gamma, tpr_l, tmr_l, gamma_l, derivs_layers, work);


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (! vector) {
/*
if (! symmetric)
*/
          eig_2n_gen_real(n_quad_v, n_derivs, tpr, tmr, gamma, nu, X_p, X_m, tpr_l, tmr_l, gamma_l, nu_l, X_p_l, X_m_l, eigen_solver_real, derivs_layers, save_tree, work);
/*
else
          eig_2n_sym_real(n_quad_v, n_derivs, tpr, tmr, gamma, nu, X_p, X_m, tpr_l, tmr_l, gamma_l, nu_l, X_p_l, X_m_l, aux, aul, derivs, work);
*/
          calc_global_r_and_t(n_quad_v, n_derivs, ltau, ltau_l, nu, X_p, X_m, R_p, T_p, nu_l, X_p_l, X_m_l, R_p_l, T_p_l, symmetric, derivs_layers, save_tree, work);
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     else {
          eig_2n_gen_complex(n_quad_v, n_derivs, tpr, tmr, gamma, nu_c, X_p_c, X_m_c, tpr_l, tmr_l, gamma_l, nu_l_c, X_p_l_c, X_m_l_c, eigen_solver_complex, derivs_layers, save_tree, work);

          calc_global_r_and_t2(n_quad_v, n_derivs, ltau, ltau_l, nu_c, X_p_c, X_m_c, R_p, T_p, nu_l_c, X_p_l_c, X_m_l_c, R_p_l, T_p_l, symmetric, derivs_layers, save_tree, work);

          dmat_mul_D_A(n_quad, n_stokes, R_p, R_p);
          phase_matrix_symmetry2(n_quad, n_stokes, T_p, R_p, T_m, R_m, 1.);

          if (flags_or(derivs_layers, n_derivs)) {
               for (i = 0; i < n_derivs; ++i) {
                    if (! derivs_layers[i])
                         continue;

                    dmat_mul_D_A(n_quad, n_stokes, R_p_l[i], R_p_l[i]);
                    phase_matrix_symmetry2(n_quad, n_stokes, T_p_l[i], R_p_l[i], T_m_l[i], R_m_l[i], 1.);
               }
          }
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (solar) {
          if (! CLASSICAL_PARTICULAR_SOLUTION_USE_2D) {
               if (! greens)
                    build_source_vectors_solar_classic_1n(n_quad, n_stokes, n_derivs, qx_v, F_0, omega, omega_l, as_0, as_0_l, P_x0_p, P_x0_m, tpr, tmr, gamma, Fs_p, Fs_m, P_x0_p_l, P_x0_m_l, tpr_l, tmr_l, gamma_l, Fs_p_l, Fs_m_l, derivs_layers, derivs_beam, save_tree, work);
               else {
if (! vector)
                    build_source_vectors_solar_green_s_1n (n_quad, n_stokes, n_derivs, qx_v, qw_v, F_0, omega, omega_l, ltau, ltau_l, as_0, as_0_l, atran, atran_l, P_x0_p, P_x0_m, nu, X_p, X_m, Fs_p, Fs_m, Fs1_p, Fs1_m, P_x0_p_l, P_x0_m_l, nu_l, X_p_l, X_m_l, Fs_p_l, Fs_m_l, Fs1_p_l, Fs1_m_l, derivs_layers, derivs_beam, save_tree, work);
else
                    build_source_vectors_solar_green_s_1n2(n_quad, n_stokes, n_derivs, qx_v, qw_v, F_0, omega, omega_l, ltau, ltau_l, as_0, as_0_l, atran, atran_l, P_x0_p, P_x0_m, nu_c, X_p_c, X_m_c, Fs_p, Fs_m, Fs1_p, Fs1_m, P_x0_p_l, P_x0_m_l, nu_l_c, X_p_l_c, X_m_l_c, Fs_p_l, Fs_m_l, Fs1_p_l, Fs1_m_l, derivs_layers, derivs_beam, save_tree, work);
               }
          }
          else {
               if (! vector)
                    build_source_vectors_solar_classic_2n (n_quad_v, n_derivs, qx_v, F_0, omega, omega_l, as_0, as_0_l, P_x0_p, P_x0_m, r_p, t_p, Fs_p, Fs_m, P_x0_p_l, P_x0_m_l, r_p_l, t_p_l, Fs_p_l, Fs_m_l, derivs_layers, derivs_beam, work);
               else
                    build_source_vectors_solar_classic_2n2(n_quad_v, n_derivs, qx_v, F_0, omega, omega_l, as_0, as_0_l, P_x0_p, P_x0_m, r_p, t_p, r_m, t_m, Fs_p, Fs_m, P_x0_p_l, P_x0_m_l, r_p_l, t_p_l, r_m_l, t_m_l, Fs_p_l, Fs_m_l, derivs_layers, derivs_beam, work);
          }

          if (! vector) {
               if (! greens)
                    build_global_source_solar (n_quad_v, n_derivs, atran, atran_l, R_p, T_p, Fs_p, Fs_m, S_p, S_m, R_p_l, T_p_l, Fs_p_l, Fs_m_l, S_p_l, S_m_l, derivs_layers, derivs_beam, work, Fs1_p, Fs1_m, Fs1_p_l, Fs1_m_l);
               else
                    build_global_source_solar (n_quad_v, n_derivs, -1.,   atran_l, R_p, T_p, Fs_p, Fs_m, S_p, S_m, R_p_l, T_p_l, Fs_p_l, Fs_m_l, S_p_l, S_m_l, derivs_layers, derivs_beam, work, Fs1_p, Fs1_m, Fs1_p_l, Fs1_m_l);
          }
          else {
               if (! greens)
                    build_global_source_solar2(n_quad_v, n_derivs, atran, atran_l, R_p, T_p, R_m, T_m, Fs_p, Fs_m, S_p, S_m, R_p_l, T_p_l, R_m_l, T_m_l, Fs_p_l, Fs_m_l, S_p_l, S_m_l, derivs_layers, derivs_beam, work, Fs1_p, Fs1_m, Fs1_p_l, Fs1_m_l);
               else
                    build_global_source_solar2(n_quad_v, n_derivs, -1.,   atran_l, R_p, T_p, R_m, T_m, Fs_p, Fs_m, S_p, S_m, R_p_l, T_p_l, R_m_l, T_m_l, Fs_p_l, Fs_m_l, S_p_l, S_m_l, derivs_layers, derivs_beam, work, Fs1_p, Fs1_m, Fs1_p_l, Fs1_m_l);
          }
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (thermal) {
          build_source_vectors_thermal2(n_quad, n_stokes, n_derivs, qx_v, planck0, planck1, planck0_l, planck1_l, omega, omega_l, ltau, ltau_l, r_p, t_p, r_m, t_m, Ft0_p, Ft0_m, Ft1_p, Ft1_m, r_p_l, t_p_l, r_m_l, t_m_l, Ft0_p_l, Ft0_m_l, Ft1_p_l, Ft1_m_l, derivs_layers, derivs_thermal, work);

          if (! vector)
               build_global_source_thermal(n_quad_v, n_derivs, R_p, T_p, R_p, T_p, Ft0_p, Ft0_m, Ft1_p, Ft1_m, Sl_p, Sl_m, R_p_l, T_p_l, R_p_l, T_p_l, Ft0_p_l, Ft0_m_l, Ft1_p_l, Ft1_m_l, Sl_p_l, Sl_m_l, derivs_layers, derivs_thermal, work);
          else
               build_global_source_thermal(n_quad_v, n_derivs, R_p, T_p, R_m, T_m, Ft0_p, Ft0_m, Ft1_p, Ft1_m, Sl_p, Sl_m, R_p_l, T_p_l, R_m_l, T_m_l, Ft0_p_l, Ft0_m_l, Ft1_p_l, Ft1_m_l, Sl_p_l, Sl_m_l, derivs_layers, derivs_thermal, work);
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (save_tree.t) {
          dvec_copy(save->Fs_p,  Fs_p,  n_quad_v);
          dvec_copy(save->Fs_m,  Fs_m,  n_quad_v);

          dmat_copy(save->tpr,   tpr,   n_quad_v, n_quad_v);
          dmat_copy(save->tmr,   tmr,   n_quad_v, n_quad_v);
          dmat_copy(save->gamma, gamma, n_quad_v, n_quad_v);

          if (! vector) {
               dvec_copy(save->nu,    nu,    n_quad_v);
               dmat_copy(save->X_p,   X_p,   n_quad_v, n_quad_v);
               dmat_copy(save->X_m,   X_m,   n_quad_v, n_quad_v);
          }
          else {
               zvec_copy(save->nu_c,  nu_c,  n_quad_v);
               zmat_copy(save->X_p_c, X_p_c, n_quad_v, n_quad_v);
               zmat_copy(save->X_m_c, X_m_c, n_quad_v, n_quad_v);
          }
     }
#ifdef USE_AD_FOR_TL_CALC_RTM_EIG_RTS
     rtm_eig_rts_tl_with_ad(n_quad, n_stokes, n_derivs, F_0,
                            qx_v, qw_v, planck0, planck1,
                            omega, omega_l, ltau, ltau_l,
                            as_0, as_0_l, atran, atran_l,
                            P_x0_p, P_x0_m,
                            r_p, t_p, r_m, t_m,
                            R_p, T_p, R_m, T_m,
                            S_p, S_m,
                            P_x0_p_l, P_x0_m_l,
                            r_p_l, t_p_l, r_m_l, t_m_l,
                            R_p_l, T_p_l, R_m_l, T_m_l,
                            S_p_l, S_m_l,
                            symmetric, thermal, vector,
                            eigen_solver_real, eigen_solver_complex,
                            derivs_layers, derivs_beam, save_tree, work);
#endif
}
