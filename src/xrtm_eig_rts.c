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
                 double *qx_v, double *qw_v, double planck0, double planck1,
                 double omega, double *omega_l, double ltau, double *ltau_l,
                 double as_0, double *as_0_l, double atran, double *atran_l,
                 double  *P_0p, double  *P_0m,
                 double  **r_p, double  **t_p, double  **r_m, double  **t_m,
                 double  **R_p, double  **T_p, double  **R_m, double  **T_m,
                 double  *S_p, double  *S_m,
                 double **P_0p_l, double **P_0m_l,
                 double ***r_p_l, double ***t_p_l, double ***r_m_l, double ***t_m_l,
                 double ***R_p_l, double ***T_p_l, double ***R_m_l, double ***T_m_l,
                 double **S_p_l, double **S_m_l,
                 int symmetric, int thermal, int vector,
                 int eigen_solver_real, int eigen_solver_complex,
                 uchar *derivs_h, uchar *derivs_p,
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
     double *F_p;
     double *F_m;

     double *F0_p;
     double *F0_m;
     double *F1_p;
     double *F1_m;

     double *St_p;
     double *St_m;

     double **F_p_l;
     double **F_m_l;

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
     F_p = get_work1(&work, WORK_DX);
     F_m = get_work1(&work, WORK_DX);

     if (thermal) {
          F0_p = get_work1(&work, WORK_DX);
          F0_m = get_work1(&work, WORK_DX);
          F1_p = get_work1(&work, WORK_DX);
          F1_m = get_work1(&work, WORK_DX);

          St_p = get_work1(&work, WORK_DX);
          St_m = get_work1(&work, WORK_DX);
     }

     tmr   = get_work1(&work, WORK_DXX);
     tpr   = get_work1(&work, WORK_DXX);
     gamma = get_work1(&work, WORK_DXX);

     if (flags_or(derivs_p, n_derivs)) {
          F_p_l = get_work2(&work, WORK_DX, WORK_DERIVS_V, derivs_p);
          F_m_l = get_work2(&work, WORK_DX, WORK_DERIVS_V, derivs_p);
     }

     if (flags_or(derivs_h, n_derivs)) {
          tmr_l   = get_work2(&work, WORK_DXX, WORK_DERIVS_V, derivs_h);
          tpr_l   = get_work2(&work, WORK_DXX, WORK_DERIVS_V, derivs_h);
          gamma_l = get_work2(&work, WORK_DXX, WORK_DERIVS_V, derivs_h);
     }

     if (! vector) {
          nu  = get_work1(&work, WORK_DX);

          X_p = get_work1(&work, WORK_DXX);
          X_m = get_work1(&work, WORK_DXX);

          if (flags_or(derivs_h, n_derivs)) {
               nu_l  = get_work2(&work, WORK_DX, WORK_DERIVS_V, derivs_h);

               X_p_l = get_work2(&work, WORK_DXX, WORK_DERIVS_V, derivs_h);
               X_m_l = get_work2(&work, WORK_DXX, WORK_DERIVS_V, derivs_h);
          }
     }
     else {
          nu_c  = get_work1(&work, WORK_ZX);

          X_p_c = get_work1(&work, WORK_ZXX);
          X_m_c = get_work1(&work, WORK_ZXX);

          if (flags_or(derivs_h, n_derivs)) {
               nu_l_c  = get_work2(&work, WORK_ZX, WORK_DERIVS_V, derivs_h);

               X_p_l_c = get_work2(&work, WORK_ZXX, WORK_DERIVS_V, derivs_h);
               X_m_l_c = get_work2(&work, WORK_ZXX, WORK_DERIVS_V, derivs_h);
          }
     }

     if (symmetric) {
/*
          aux = get_work1(&work, WORK_DXX);

          if (flags_or(derivs_h, n_derivs))
               aul = get_work2(&work, WORK_DXX, WORK_DERIVS_V, derivs_h);
*/
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     build_txr(n_quad, n_stokes, n_derivs, r_p, t_p, tpr, tmr, r_p_l, t_p_l, tpr_l, tmr_l, derivs_h, work);

     build_gamma(n_quad_v, n_derivs, tpr, tmr, gamma, tpr_l, tmr_l, gamma_l, derivs_h, work);


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (! vector) {
/*
if (! symmetric)
*/
          eig_2n_gen_real(n_quad_v, n_derivs, tpr, tmr, gamma, nu, X_p, X_m, tpr_l, tmr_l, gamma_l, nu_l, X_p_l, X_m_l, eigen_solver_real, derivs_h, save_tree, work);
/*
else
          eig_2n_sym_real(n_quad_v, n_derivs, tpr, tmr, gamma, nu, X_p, X_m, tpr_l, tmr_l, gamma_l, nu_l, X_p_l, X_m_l, aux, aul, derivs, work);
*/
          calc_global_r_and_t(n_quad_v, n_derivs, ltau, ltau_l, nu, X_p, X_m, R_p, T_p, nu_l, X_p_l, X_m_l, R_p_l, T_p_l, symmetric, derivs_h, save_tree, work);
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     else {
          eig_2n_gen_complex(n_quad_v, n_derivs, tpr, tmr, gamma, nu_c, X_p_c, X_m_c, tpr_l, tmr_l, gamma_l, nu_l_c, X_p_l_c, X_m_l_c, eigen_solver_complex, derivs_h, save_tree, work);

          calc_global_r_and_t2(n_quad_v, n_derivs, ltau, ltau_l, nu_c, X_p_c, X_m_c, R_p, T_p, nu_l_c, X_p_l_c, X_m_l_c, R_p_l, T_p_l, symmetric, derivs_h, save_tree, work);

          dmat_mul_D_A(n_quad, n_stokes, R_p, R_p);
          phase_matrix_symmetry2(n_quad, n_stokes, T_p, R_p, T_m, R_m, 1.);

          if (flags_or(derivs_h, n_derivs)) {
               for (i = 0; i < n_derivs; ++i) {
                    if (! derivs_h[i])
                         continue;

                    dmat_mul_D_A(n_quad, n_stokes, R_p_l[i], R_p_l[i]);
                    phase_matrix_symmetry2(n_quad, n_stokes, T_p_l[i], R_p_l[i], T_m_l[i], R_m_l[i], 1.);
               }
          }
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
if (1)
     build_source_vectors_1n(n_quad, n_stokes, n_derivs, qx_v, F_0, omega, omega_l, as_0, as_0_l, P_0p, P_0m, tpr, tmr, gamma, F_p, F_m, P_0p_l, P_0m_l, tpr_l, tmr_l, gamma_l, F_p_l, F_m_l, derivs_h, derivs_p, save_tree, work);
else {
     if (! vector)
          build_source_vectors_2n(n_quad_v, n_derivs, qx_v, F_0, omega, omega_l, as_0, as_0_l, P_0p, P_0m, r_p, t_p, F_p, F_m, P_0p_l, P_0m_l, r_p_l, t_p_l, F_p_l, F_m_l, derivs_h, derivs_p, work);
     else
          build_source_vectors_2n2(n_quad_v, n_derivs, qx_v, F_0, omega, omega_l, as_0, as_0_l, P_0p, P_0m, r_p, t_p, r_m, t_m, F_p, F_m, P_0p_l, P_0m_l, r_p_l, t_p_l, r_m_l, t_m_l, F_p_l, F_m_l, derivs_h, derivs_p, work);
}
     if (! vector)
          build_global_source (n_quad_v, n_derivs, atran, atran_l, R_p, T_p, F_p, F_m, S_p, S_m, R_p_l, T_p_l, F_p_l, F_m_l, S_p_l, S_m_l, derivs_h, derivs_p, work);
     else
          build_global_source2(n_quad_v, n_derivs, atran, atran_l, R_p, T_p, R_m, T_m, F_p, F_m, S_p, S_m, R_p_l, T_p_l, R_m_l, T_m_l, F_p_l, F_m_l, S_p_l, S_m_l, derivs_h, derivs_p, work);


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
if (thermal) {
     build_source_vectors_thermal(n_quad, n_stokes, n_derivs, qx_v, F_0, planck0, planck1, omega, omega_l, ltau, ltau_l, as_0, as_0_l, P_0p, P_0m, r_p, t_p, r_m, t_m, F0_p, F0_m, F1_p, F1_m, P_0p_l, P_0m_l, r_p_l, t_p_l, r_m_l, t_m_l, NULL, NULL, NULL, NULL, thermal, derivs_h, derivs_p, work);

     if (! vector)
          build_global_source_thermal(n_quad_v, n_derivs, atran, atran_l, R_p, T_p, R_p, T_p, F0_p, F0_m, F1_p, F1_m, St_p, St_m, R_p_l, T_p_l, R_p_l, T_p_l, NULL, NULL, NULL, NULL, S_p_l, S_m_l, derivs_h, derivs_p, work);
     else
          build_global_source_thermal(n_quad_v, n_derivs, atran, atran_l, R_p, T_p, R_m, T_m, F0_p, F0_m, F1_p, F1_m, St_p, St_m, R_p_l, T_p_l, R_m_l, T_m_l, NULL, NULL, NULL, NULL, S_p_l, S_m_l, derivs_h, derivs_p, work);
}

     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (F_0 == 0. && ! thermal) {
          dvec_zero(S_p, n_quad_v);
          dvec_zero(S_m, n_quad_v);
     }
     else
     if (F_0 == 0. && thermal) {
          dvec_copy(S_p, St_p, n_quad_v);
          dvec_copy(S_m, St_m, n_quad_v);
     }
     else
     if (F_0 != 0. && thermal) {
          dvec_add (S_p, St_p, S_p, n_quad_v);
          dvec_add (S_m, St_m, S_m, n_quad_v);
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (save_tree.t) {
          dvec_copy(save->F_p,   F_p,   n_quad_v);
          dvec_copy(save->F_m,   F_m,   n_quad_v);

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
                            P_0p, P_0m,
                            r_p, t_p, r_m, t_m,
                            R_p, T_p, R_m, T_m,
                            S_p, S_m,
                            P_0p_l, P_0m_l,
                            r_p_l, t_p_l, r_m_l, t_m_l,
                            R_p_l, T_p_l, R_m_l, T_m_l,
                            S_p_l, S_m_l,
                            symmetric, thermal, vector,
                            eigen_solver_real, eigen_solver_complex,
                            derivs_h, derivs_p, save_tree, work);
#endif
}
