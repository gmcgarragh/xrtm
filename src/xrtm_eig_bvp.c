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
#include "xrtm_eig_bvp.h"
#include "xrtm_eig_bvp_a.h"
#include "xrtm_eig_util.h"
#include "xrtm_matrix.h"
#include "xrtm_save_tree.h"
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
#define WORK_XLAYERS	WORK_DLAYERS
#define XEXP		exp
#define XREAL(x)	(x)

#define SOLVE_BVP					solve_bvp
#define SOLVE_BVP_SAVE_TREE_STRING			"solve_bvp"
#define SOLVE_BVP_TL_WITH_AD				solve_bvp_tl_with_ad
#define FORWARD_SAVE_SOLVE_BVP_DATA			forward_save_solve_bvp_data
#define FORWARD_SAVE_SOLVE_BVP_ALLOC			forward_save_solve_bvp_alloc
#define CALC_RADIANCE_LEVELS				calc_radiance_levels
#define CALC_RADIANCE_LEVELS_SAVE_TREE_STRING		"calc_radiance_levels"
#define FORWARD_SAVE_CALC_RADIANCE_LEVELS_DATA		forward_save_calc_radiance_levels_data
#define FORWARD_SAVE_CALC_RADIANCE_LEVELS_ALLOC		forward_save_calc_radiance_levels_alloc
#define CALC_RADIANCE_TAUS				calc_radiance_taus

#include "type_set.h"

#include "xrtm_eig_bvp_x.c"

#include "type_unset.h"

#undef TYPE
#undef TYPE_PREFIX
#undef TYPE_POSTFIX
#undef WORK_XX	
#undef WORK_XXX
#undef WORK_XLAYERS
#undef XEXP
#undef XREAL

#undef SOLVE_BVP
#undef SOLVE_BVP_TL_WITH_AD
#undef SOLVE_BVP_SAVE_TREE_STRING
#undef FORWARD_SAVE_SOLVE_BVP_DATA
#undef FORWARD_SAVE_SOLVE_BVP_ALLOC
#undef CALC_RADIANCE_LEVELS
#undef CALC_RADIANCE_LEVELS_SAVE_TREE_STRING
#undef FORWARD_SAVE_CALC_RADIANCE_LEVELS_DATA
#undef FORWARD_SAVE_CALC_RADIANCE_LEVELS_ALLOC
#undef CALC_RADIANCE_TAUS



/*******************************************************************************
 *
 ******************************************************************************/
#define TYPE		dcomplex
#define TYPE_PREFIX	z
#define TYPE_POSTFIX	dc
#define WORK_XX		WORK_ZX
#define WORK_XXX	WORK_ZXX
#define WORK_XLAYERS	WORK_ZLAYERS
#define XEXP		cexp
#define XREAL(x)	(creal(x))

#define SOLVE_BVP					solve_bvp2
#define SOLVE_BVP_TL_WITH_AD				solve_bvp_tl_with_ad2
#define SOLVE_BVP_SAVE_TREE_STRING			"solve_bvp2"
#define FORWARD_SAVE_SOLVE_BVP_DATA			forward_save_solve_bvp_data2
#define FORWARD_SAVE_SOLVE_BVP_ALLOC			forward_save_solve_bvp_alloc2
#define CALC_RADIANCE_LEVELS				calc_radiance_levels2
#define CALC_RADIANCE_LEVELS_SAVE_TREE_STRING		"calc_radiance_levels2"
#define FORWARD_SAVE_CALC_RADIANCE_LEVELS_DATA		forward_save_calc_radiance_levels_data2
#define FORWARD_SAVE_CALC_RADIANCE_LEVELS_ALLOC		forward_save_calc_radiance_levels_alloc2
#define CALC_RADIANCE_TAUS				calc_radiance_taus2

#include "type_set.h"

#include "xrtm_eig_bvp_x.c"

#include "type_unset.h"

#undef TYPE
#undef TYPE_PREFIX
#undef TYPE_POSTFIX
#undef WORK_XX	
#undef WORK_XXX
#undef WORK_XLAYERS
#undef XEXP
#undef XREAL

#undef SOLVE_BVP
#undef SOLVE_BVP_TL_WITH_AD
#undef SOLVE_BVP_SAVE_TREE_STRING
#undef FORWARD_SAVE_SOLVE_BVP_DATA
#undef FORWARD_SAVE_SOLVE_BVP_ALLOC
#undef CALC_RADIANCE_LEVELS
#undef CALC_RADIANCE_LEVELS_SAVE_TREE_STRING
#undef FORWARD_SAVE_CALC_RADIANCE_LEVELS_DATA
#undef FORWARD_SAVE_CALC_RADIANCE_LEVELS_ALLOC
#undef CALC_RADIANCE_TAUS



/*******************************************************************************
 *
 ******************************************************************************/
void rtm_eig_bvp(int i_four, int n_quad, int n_stokes, int n_derivs, int n_layers,
       double qf, double *qx_v, double *qw_v, double F_0, double mu_0,
       int *ulevels, double *utaus, int n_ulevels,
       double *umus, int n_umus, double *planck,
       double *omega, double **omega_l, double *ltau, double **ltau_l,
       double **Rs_qq, double ***Rs_qq_l,
       double *Rs_u0, double **Rs_u0_l, double **Rs_uq, double ***Rs_uq_l, 
       double *btran, double **btran_l,
       double *as_0, double **as_0_l, double *atran, double **atran_l,
       double **P_q0_mm, double **P_q0_pm, double **P_u0_mm, double **P_u0_pm, 
       double ***P_uq_pp, double ***P_uq_mp, double ***P_uq_mm, double ***P_uq_pm, 
       double ***r_p, double ***t_p, double ***r_m, double ***t_m,
       double ***P_q0_mm_l, double ***P_q0_pm_l, double ***P_u0_mm_l, double ***P_u0_pm_l,
       double ****P_uq_pp_l, double ****P_uq_mp_l, double ****P_uq_mm_l, double ****P_uq_pm_l, 
       double ****r_p_l, double ****t_p_l, double ****r_m_l, double ****t_m_l,
       double *I1_m, double **I1_m_l, double *In_p, double **In_p_l,
       double **I_p, double **I_m, double ***I_p_l, double ***I_m_l,
       int sfi, int surface, int thermal, int upwelling,
       int downwelling, int utau_output, int vector,
       int eigen_solver_real, int eigen_solver_complex,
       uchar **derivs_h, uchar **derivs_p, save_tree_data save_tree, work_data work) {

     uchar *derivs_h2;
     uchar *derivs_p2;

     int i;
     int j;

     int n_quad_x;
     int n_quad_v;
     int n_quad_v_x;
     int n_umus_v;

     int i_quad2;
     int n_quad2;

     int n_comp;

     double *omega_l2;
     double *ltau_l2;
     double *btran_l2;
     double *as_0_l2;

     double **P_q0_mm_l2;
     double **P_q0_pm_l2;

     double ***r_p_l2;
     double ***t_p_l2;

     double ***r_m_l2;
     double ***t_m_l2;

     double **tmr;
     double **tpr;
     double **gamma;

     double ***tmr_l;
     double ***tpr_l;
     double ***gamma_l;

     double **F_p;
     double **F_m;

     double **F0_p;
     double **F0_m;
     double **F1_p;
     double **F1_m;

     double **F_p_l2;
     double **F_m_l2;

     double ***F_p_l;
     double ***F_m_l;

     double *B;
     double **B_l;

     double **nu_l2;

     double **nu;
     double ***nu_l;

     double ***X_p;
     double ***X_m;

     double ***X_p_l2;
     double ***X_m_l2;

     double ****X_p_l;
     double ****X_m_l;

     double *I_p2;
     double *I_m2;
     double **I_p_l2;
     double **I_m_l2;

     dcomplex *B_c;
     dcomplex **B_l_c;

     dcomplex **nu_l_c2;

     dcomplex **nu_c;
     dcomplex ***nu_l_c;

     dcomplex ***X_p_l_c2;
     dcomplex ***X_m_l_c2;

     dcomplex ***X_p_c;
     dcomplex ***X_m_c;

     dcomplex ****X_p_l_c;
     dcomplex ****X_m_l_c;

     forward_save_rtm_eig_bvp_data *save;

     work_data work2;


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     calc_quad_and_umus(n_quad, n_umus, n_stokes, &n_quad_x, NULL, &n_quad_v, &n_quad_v_x, NULL, &n_umus_v, sfi);

     n_comp = 2 * n_quad_v_x * n_layers;


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (save_tree.t) {
          save_tree_encode_s(&save_tree, "rtm_eig_bvp");

          if (save_tree_retrieve_data(&save_tree, forward_save_rtm_eig_bvp_data, &save))
               forward_save_rtm_eig_bvp_alloc(save, n_layers, n_quad_v_x, n_comp, thermal, vector);
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     tmr   = get_work1(&work, WORK_DXX);
     tpr   = get_work1(&work, WORK_DXX);
     gamma = get_work1(&work, WORK_DXX);

     F_p = get_work2(&work, WORK_DX, WORK_LAYERS_V, NULL);
     F_m = get_work2(&work, WORK_DX, WORK_LAYERS_V, NULL);

     if (thermal) {
          F0_p = get_work2(&work, WORK_DX, WORK_LAYERS_V, NULL);
          F0_m = get_work2(&work, WORK_DX, WORK_LAYERS_V, NULL);
          F1_p = get_work2(&work, WORK_DX, WORK_LAYERS_V, NULL);
          F1_m = get_work2(&work, WORK_DX, WORK_LAYERS_V, NULL);
     }

     if (flags_or2(derivs_p, n_layers, n_derivs)) {
          F_p_l = get_work3(&work, WORK_DX, WORK_BOTH_V, derivs_p);
          F_m_l = get_work3(&work, WORK_DX, WORK_BOTH_V, derivs_p);
     }

     if (! vector) {
          B = get_work_d1(&work, n_comp);
          if (n_derivs > 0)
               B_l = get_work_d2(&work, n_derivs, n_comp);

          nu  = get_work2(&work, WORK_DX, WORK_LAYERS_V, NULL);

          X_p = get_work2(&work, WORK_DXX, WORK_LAYERS_V, NULL);
          X_m = get_work2(&work, WORK_DXX, WORK_LAYERS_V, NULL);

          if (flags_or2(derivs_h, n_layers, n_derivs)) {
               nu_l  = get_work3(&work, WORK_DX, WORK_BOTH_V, derivs_h);

               X_p_l = get_work3(&work, WORK_DXX, WORK_BOTH_V, derivs_h);
               X_m_l = get_work3(&work, WORK_DXX, WORK_BOTH_V, derivs_h);
          }
     }
     else {
          B_c = get_work_dc1(&work, n_comp);
          if (n_derivs > 0)
               B_l_c = get_work_dc2(&work, n_derivs, n_comp);

          nu_c  = get_work2(&work, WORK_ZX, WORK_LAYERS_V, NULL);

          X_p_c = get_work2(&work, WORK_ZXX, WORK_LAYERS_V, NULL);
          X_m_c = get_work2(&work, WORK_ZXX, WORK_LAYERS_V, NULL);

          if (flags_or2(derivs_h, n_layers, n_derivs)) {
               nu_l_c  = get_work3(&work, WORK_ZX, WORK_BOTH_V, derivs_h);

               X_p_l_c = get_work3(&work, WORK_ZXX, WORK_BOTH_V, derivs_h);
               X_m_l_c = get_work3(&work, WORK_ZXX, WORK_BOTH_V, derivs_h);
          }
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     for (i = 0; i < n_layers; ++i) {
          work2 = work;

          if (save_tree.t)
               save_tree_recode_i(&save_tree, i, i == 0);


          if (flags_or2(derivs_h, n_layers, n_derivs)) {
               tmr_l   = get_work2(&work2, WORK_DXX, WORK_DERIVS_V, derivs_h[i]);
               tpr_l   = get_work2(&work2, WORK_DXX, WORK_DERIVS_V, derivs_h[i]);
               gamma_l = get_work2(&work2, WORK_DXX, WORK_DERIVS_V, derivs_h[i]);
          }

          if (n_derivs > 0) {
               derivs_h2 = derivs_h[i];
               derivs_p2 = derivs_p[i];
          }

          if (n_derivs > 0 && flags_or(derivs_h[i], n_derivs)) {
               omega_l2   = omega_l  [i];
               ltau_l2    = ltau_l   [i];
               P_q0_mm_l2 = P_q0_mm_l[i];
               P_q0_pm_l2 = P_q0_pm_l[i];
               r_p_l2     = r_p_l    [i];
               t_p_l2     = t_p_l    [i];
          }

          if (n_derivs > 0 && flags_or(derivs_p[i], n_derivs)) {
               btran_l2   = btran_l[i];
               as_0_l2    = as_0_l [i];
               F_p_l2     = F_p_l  [i];
               F_m_l2     = F_m_l  [i];
          }


          /*--------------------------------------------------------------------
           *
           *------------------------------------------------------------------*/
          build_txr(n_quad_x, n_stokes, n_derivs, r_p[i], t_p[i], tpr, tmr, r_p_l2, t_p_l2, tpr_l, tmr_l, derivs_h2, work2);

          build_gamma(n_quad_v_x, n_derivs, tpr, tmr, gamma, tpr_l, tmr_l, gamma_l, derivs_h2, work2);


          /*--------------------------------------------------------------------
           *
           *------------------------------------------------------------------*/
          if (! vector) {
               if (n_derivs > 0 && flags_or(derivs_h[i], n_derivs)) {
                    nu_l2      = nu_l [i];
                    X_p_l2     = X_p_l[i];
                    X_m_l2     = X_m_l[i];
               }

               eig_2n_gen_real(n_quad_v_x, n_derivs, tpr, tmr, gamma, nu[i], X_p[i], X_m[i], tpr_l, tmr_l, gamma_l, nu_l2, X_p_l2, X_m_l2, eigen_solver_real, derivs_h2, save_tree, work2);
if (1)
               build_source_vectors_1n(n_quad_x, n_stokes, n_derivs, qx_v, F_0, omega[i], omega_l2, as_0[i], as_0_l2, P_q0_mm[i], P_q0_pm[i], tpr, tmr, gamma, F_p[i], F_m[i], P_q0_mm_l2, P_q0_pm_l2, tpr_l, tmr_l, gamma_l, F_p_l2, F_m_l2, derivs_h2, derivs_p2, save_tree, work2);
else
               build_source_vectors_2n(n_quad_v_x, n_derivs, qx_v, F_0, omega[i], omega_l2, as_0[i], as_0_l2, P_q0_mm[i], P_q0_pm[i], r_p[i], t_p[i], F_p[i], F_m[i], P_q0_mm_l2, P_q0_pm_l2, r_p_l2, t_p_l2, F_p_l2, F_m_l2, derivs_h2, derivs_p2, work2);

               if (save_tree.t) {
                    copy_array1_d(save->F_p [i], F_p[i], n_quad_v_x);
                    copy_array1_d(save->F_m [i], F_m[i], n_quad_v_x);
               }

               scale_source_vectors(n_quad_v_x, n_derivs, btran[i], btran_l2, F_p[i], F_m[i], F_p[i], F_m[i], F_p_l2, F_m_l2, F_p_l2, F_m_l2, derivs_h2, derivs_p2, work2);

               if (save_tree.t) {
                    copy_array1_d(save->F_p2[i], F_p[i], n_quad_v_x);
                    copy_array1_d(save->F_m2[i], F_m[i], n_quad_v_x);
               }

if (thermal)
               build_source_vectors_thermal(n_quad_x, n_stokes, n_derivs, qx_v, F_0, planck[i], planck[i+1], omega[i], omega_l2, ltau[i], ltau_l2, as_0[i], as_0_l2, P_q0_mm[i], P_q0_pm[i], r_p[i], t_p[i], r_p[i], t_p[i], F0_p[i], F0_m[i], F1_p[i], F1_m[i], P_q0_mm_l2, P_q0_pm_l2, r_p_l2, t_p_l2, r_m_l2, t_m_l2, NULL, NULL, NULL, NULL, thermal, derivs_h2, derivs_p2, work2);
          }


          /*--------------------------------------------------------------------
           *
           *------------------------------------------------------------------*/
          else {
               if (n_derivs > 0 && flags_or(derivs_h[i], n_derivs)) {
                    r_m_l2     = r_m_l[i];
                    t_m_l2     = t_m_l[i];
               }

               if (n_derivs > 0 && flags_or(derivs_h[i], n_derivs)) {
                    nu_l_c2    = nu_l_c [i];
                    X_p_l_c2   = X_p_l_c[i];
                    X_m_l_c2   = X_m_l_c[i];
               }

               eig_2n_gen_complex(n_quad_v_x, n_derivs, tpr, tmr, gamma, nu_c[i], X_p_c[i], X_m_c[i], tpr_l, tmr_l, gamma_l, nu_l_c2, X_p_l_c2, X_m_l_c2, eigen_solver_complex, derivs_h2, save_tree, work2);
if (1) {
               build_source_vectors_1n(n_quad_x, n_stokes, n_derivs, qx_v, F_0, omega[i], omega_l2, as_0[i], as_0_l2, P_q0_mm[i], P_q0_pm[i], tpr, tmr, gamma, F_p[i], F_m[i], P_q0_mm_l2, P_q0_pm_l2, tpr_l, tmr_l, gamma_l, F_p_l2, F_m_l2, derivs_h2, derivs_p2, save_tree, work2);

               dm_v_mul_D_A(n_quad_x, n_stokes, F_m[i], F_m[i]);
               if (n_derivs > 0 && flags_or(derivs_p[i], n_derivs)) {
                    for (j = 0; j < n_derivs; ++j) {
                         if (! derivs_p[i][j])
                              continue;

                         dm_v_mul_D_A(n_quad_x, n_stokes, F_m_l[i][j], F_m_l[i][j]);
                    }
               }
}
else
               build_source_vectors_2n2(n_quad_v_x, n_derivs, qx_v, F_0, omega[i], omega_l2, as_0[i], as_0_l[i], P_q0_mm[i], P_q0_pm[i], r_p[i], t_p[i], r_m[i], t_m[i], F_p[i], F_m[i], P_q0_mm_l2, P_q0_pm_l2, r_p_l2, t_p_l2, r_m_l2, t_m_l2, F_p_l2, F_m_l2, derivs_h2, derivs_p2, work2);

               if (save_tree.t) {
                    copy_array1_d(save->F_p [i], F_p[i], n_quad_v_x);
                    copy_array1_d(save->F_m [i], F_m[i], n_quad_v_x);
               }

               scale_source_vectors(n_quad_v_x, n_derivs, btran[i], btran_l2, F_p[i], F_m[i], F_p[i], F_m[i], F_p_l2, F_m_l2, F_p_l2, F_m_l2, derivs_h2, derivs_p2, work2);
               if (save_tree.t) {
                    copy_array1_d(save->F_p2[i], F_p[i], n_quad_v_x);
                    copy_array1_d(save->F_m2[i], F_m[i], n_quad_v_x);
               }

if (thermal)
               build_source_vectors_thermal(n_quad_x, n_stokes, n_derivs, qx_v, F_0, planck[i], planck[i+1], omega[i], omega_l2, ltau[i], ltau_l2, as_0[i], as_0_l2, P_q0_mm[i], P_q0_pm[i], r_p[i], t_p[i], r_m[i], t_m[i], F0_p[i], F0_m[i], F1_p[i], F1_m[i], P_q0_mm_l2, P_q0_pm_l2, r_p_l2, t_p_l2, r_m_l2, t_m_l2, NULL, NULL, NULL, NULL, thermal, derivs_h2, derivs_p2, work2);
          }

          if (save_tree.t) {
               copy_array2_d(save->tpr  [i], tpr,   n_quad_v_x, n_quad_v_x);
               copy_array2_d(save->tmr  [i], tmr,   n_quad_v_x, n_quad_v_x);
               copy_array2_d(save->gamma[i], gamma, n_quad_v_x, n_quad_v_x);
          }
     }

     if (save_tree.t)
          save_tree_decode_i(&save_tree);


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (! vector) {
          solve_bvp (n_quad_x, n_stokes, n_derivs, n_layers, ltau, ltau_l, Rs_qq, Rs_qq_l, atran, atran_l, nu, X_p, X_m, F_p, F_m, F0_p, F0_m, F1_p, F1_m, nu_l, X_p_l, X_m_l, F_p_l, F_m_l, NULL, NULL, NULL, NULL, B, B_l, I1_m, I1_m_l, In_p, In_p_l, surface, thermal, derivs_h, derivs_p, save_tree, work);

if (! sfi) {
          if (! utau_output)
               calc_radiance_levels(n_quad_v_x, n_layers, n_derivs, n_ulevels, ulevels, ltau, ltau_l, atran, atran_l, nu, X_p, X_m, F_p, F_m, F0_p, F0_m, F1_p, F1_m, nu_l, X_p_l, X_m_l, F_p_l, F_m_l, NULL, NULL, NULL, NULL, B, B_l, I_p, I_m, I_p_l, I_m_l, thermal, derivs_h, derivs_p, save_tree, work);
          else
               calc_radiance_taus  (n_quad_v_x, n_layers, n_derivs, n_ulevels, ulevels, utaus, ltau, ltau_l, as_0, as_0_l, atran, atran_l, nu, X_p, X_m, F_p, F_m, nu_l, X_p_l, X_m_l, F_p_l, F_m_l, B, B_l, I_p, I_m, I_p_l, I_m_l, derivs_h, derivs_p, save_tree, work);
}
     }
     else {
          solve_bvp2(n_quad_x, n_stokes, n_derivs, n_layers, ltau, ltau_l, Rs_qq, Rs_qq_l, atran, atran_l, nu_c, X_p_c, X_m_c, F_p, F_m, F0_p, F0_m, F1_p, F1_m, nu_l_c, X_p_l_c, X_m_l_c, F_p_l, F_m_l, NULL, NULL, NULL, NULL, B_c, B_l_c, I1_m, I1_m_l, In_p, In_p_l, surface, thermal, derivs_h, derivs_p, save_tree, work);
if (! sfi) {
          if (! utau_output)
               calc_radiance_levels2(n_quad_v_x, n_layers, n_derivs, n_ulevels, ulevels, ltau, ltau_l, atran, atran_l, nu_c, X_p_c, X_m_c, F_p, F_m, F0_p, F0_m, F1_p, F1_m, nu_l_c, X_p_l_c, X_m_l_c, F_p_l, F_m_l, NULL, NULL, NULL, NULL, B_c, B_l_c, I_p, I_m, I_p_l, I_m_l, thermal, derivs_h, derivs_p, save_tree, work);
          else
               calc_radiance_taus2  (n_quad_v_x, n_layers, n_derivs, n_ulevels, ulevels, utaus, ltau, ltau_l, as_0, as_0_l, atran, atran_l, nu_c, X_p_c, X_m_c, F_p, F_m, nu_l_c, X_p_l_c, X_m_l_c, F_p_l, F_m_l, B_c, B_l_c, I_p, I_m, I_p_l, I_m_l, derivs_h, derivs_p, save_tree, work);
}
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (save_tree.t) {
          for (i = 0; i < n_layers; ++i) {
               if (thermal) {
                    copy_array1_d(save->F0_p[i], F0_p[i], n_quad_v_x);
                    copy_array1_d(save->F0_m[i], F0_m[i], n_quad_v_x);
                    copy_array1_d(save->F1_p[i], F1_p[i], n_quad_v_x);
                    copy_array1_d(save->F1_m[i], F1_m[i], n_quad_v_x);
               }

               if (! vector) {
                    copy_array1_d(save->nu[i],  nu[i],  n_quad_v_x);
                    copy_array2_d(save->X_p[i], X_p[i], n_quad_v_x, n_quad_v_x);
                    copy_array2_d(save->X_m[i], X_m[i], n_quad_v_x, n_quad_v_x);
               }
               else {
                    copy_array1_dc(save->nu_c[i],  nu_c[i],  n_quad_v_x);
                    copy_array2_dc(save->X_p_c[i], X_p_c[i], n_quad_v_x, n_quad_v_x);
                    copy_array2_dc(save->X_m_c[i], X_m_c[i], n_quad_v_x, n_quad_v_x);
               }
          }

          if (! vector)
               copy_array1_d (save->B,   B,   n_comp);
          else
               copy_array1_dc(save->B_c, B_c, n_comp);
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (sfi && n_umus > 0) {
          I_p2 = get_work1(&work, WORK_DX);
          I_m2 = get_work1(&work, WORK_DX);

          if (n_derivs > 0) {
               I_p_l2 = get_work2(&work, WORK_DX, WORK_DERIVS_V, NULL);
               I_m_l2 = get_work2(&work, WORK_DX, WORK_DERIVS_V, NULL);
          }

          if (! vector)
               calc_radiance_levels (n_quad_v_x, n_layers, n_derivs, 1, &n_layers, ltau, ltau_l, atran, atran_l, nu, X_p, X_m, F_p, F_m, F0_p, F0_m, F1_p, F1_m, nu_l, X_p_l, X_m_l, F_p_l, F_m_l, NULL, NULL, NULL, NULL, B, B_l, &I_p2, &I_m2, &I_p_l2, &I_m_l2, thermal, derivs_h, derivs_p, save_tree, work);
          else
               calc_radiance_levels2(n_quad_v_x, n_layers, n_derivs, 1, &n_layers, ltau, ltau_l, atran, atran_l, nu_c, X_p_c, X_m_c, F_p, F_m, F0_p, F0_m, F1_p, F1_m, nu_l_c, X_p_l_c, X_m_l_c, F_p_l, F_m_l, NULL, NULL, NULL, NULL, B_c, B_l_c, &I_p2, &I_m2, &I_p_l2, &I_m_l2, thermal, derivs_h, derivs_p, save_tree, work);
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (! vector) {
          if (sfi && n_umus > 0) {
               if (upwelling)
                    sfi_up (i_four, n_quad, n_stokes, n_derivs, n_layers, n_umus, qf, qx_v, qw_v, F_0, mu_0, n_ulevels, ulevels, utaus, umus, omega, omega_l, ltau, ltau_l, Rs_u0, Rs_u0_l, Rs_uq, Rs_uq_l, btran, btran_l, as_0, as_0_l, atran, atran_l, P_u0_pm, P_uq_pp, P_uq_mp, P_uq_mm, P_uq_pm, nu,   X_p,   X_m,   F_p, F_m, P_u0_pm_l, P_uq_pp_l, P_uq_mp_l, P_uq_mm_l, P_uq_pm_l, nu_l,   X_p_l,   X_m_l,   F_p_l, F_m_l, B,   B_l,   &I_m2, &I_m_l2, I_p, n_quad_v, I_p_l, surface, utau_output, derivs_h, derivs_p, work);
               else {
                    for (i = 0; i < n_ulevels; ++i) {
                         init_array1_d(I_p[i]+n_quad_v, n_umus_v, 0.);
                         for (j = 0; j < n_derivs; ++j)
                              init_array1_d(I_p_l[i][j]+n_quad_v, n_umus_v, 0.);
                    }
               }

               if (downwelling)
                    sfi_dn (i_four, n_quad, n_stokes, n_derivs, n_layers, n_umus, qf, qx_v, qw_v, F_0,       n_ulevels, ulevels, utaus, umus, omega, omega_l, ltau, ltau_l,                             btran, btran_l, as_0, as_0_l, atran, atran_l, P_u0_mm, P_uq_pp, P_uq_mp, P_uq_mm, P_uq_pm, nu,   X_p,   X_m,   F_p, F_m, P_u0_mm_l, P_uq_pp_l, P_uq_mp_l, P_uq_mm_l, P_uq_pm_l, nu_l,   X_p_l,   X_m_l,   F_p_l, F_m_l, B,   B_l,   &I_m2, &I_m_l2, I_m, n_quad_v, I_m_l,          utau_output, derivs_h, derivs_p, work);
               else {
                    for (i = 0; i < n_ulevels; ++i) {
                         init_array1_d(I_m[i]+n_quad_v, n_umus_v, 0.);
                         for (j = 0; j < n_derivs; ++j)
                              init_array1_d(I_m_l[i][j]+n_quad_v, n_umus_v, 0.);
                    }
               }
          }
     }
     else {
          if (sfi && n_umus > 0) {
               if (upwelling)
                    sfi_up2(i_four, n_quad, n_stokes, n_derivs, n_layers, n_umus, qf, qx_v, qw_v, F_0, mu_0, n_ulevels, ulevels, utaus, umus, omega, omega_l, ltau, ltau_l, Rs_u0, Rs_u0_l, Rs_uq, Rs_uq_l, btran, btran_l, as_0, as_0_l, atran, atran_l, P_u0_pm, P_uq_pp, P_uq_mp, P_uq_mm, P_uq_pm, nu_c, X_p_c, X_m_c, F_p, F_m, P_u0_pm_l, P_uq_pp_l, P_uq_mp_l, P_uq_mm_l, P_uq_pm_l, nu_l_c, X_p_l_c, X_m_l_c, F_p_l, F_m_l, B_c, B_l_c, &I_m2, &I_m_l2, I_p, n_quad_v, I_p_l, surface, utau_output, derivs_h, derivs_p, work);
               else {
                    for (i = 0; i < n_ulevels; ++i) {
                         init_array1_d(I_p[i]+n_quad_v, n_umus_v, 0.);
                         for (j = 0; j < n_derivs; ++j)
                              init_array1_d(I_p_l[i][j]+n_quad_v, n_umus_v, 0.);
                    }
               }

               if (downwelling)
                    sfi_dn2(i_four, n_quad, n_stokes, n_derivs, n_layers, n_umus, qf, qx_v, qw_v, F_0,       n_ulevels, ulevels, utaus, umus, omega, omega_l, ltau, ltau_l,                             btran, btran_l, as_0, as_0_l, atran, atran_l, P_u0_mm, P_uq_pp, P_uq_mp, P_uq_mm, P_uq_pm, nu_c, X_p_c, X_m_c, F_p, F_m, P_u0_mm_l, P_uq_pp_l, P_uq_mp_l, P_uq_mm_l, P_uq_pm_l, nu_l_c, X_p_l_c, X_m_l_c, F_p_l, F_m_l, B_c, B_l_c, &I_m2, &I_m_l2, I_m, n_quad_v, I_m_l,          utau_output, derivs_h, derivs_p, work);
               else {
                    for (i = 0; i < n_ulevels; ++i) {
                         init_array1_d(I_m[i]+n_quad_v, n_umus_v, 0.);
                         for (j = 0; j < n_derivs; ++j)
                              init_array1_d(I_m_l[i][j]+n_quad_v, n_umus_v, 0.);
                    }
               }
          }
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (vector) {
          i_quad2 = n_umus == 0 ? 0 : n_quad_v;
          n_quad2 = n_umus == 0 ? n_quad : n_umus;

          for (i = 0; i < n_ulevels; ++i) {
               dm_v_mul_D_A(n_quad2, n_stokes, I_m[i]+i_quad2, I_m[i]+i_quad2);

               for (j = 0; j < n_derivs; ++j) {
                    dm_v_mul_D_A(n_quad2, n_stokes, I_m_l[i][j]+i_quad2, I_m_l[i][j]+i_quad2);
               }
          }
     }

#ifdef USE_AD_FOR_TL_CALC_RTM_EIG_BVP
     rtm_eig_bvp_tl_with_ad(i_four, n_quad, n_stokes, n_derivs, n_layers,
                            qf, qx_v, qw_v, F_0, mu_0,
                            ulevels, utaus, n_ulevels,
                            umus, n_umus, planck,
                            omega, omega_l, ltau, ltau_l,
                            Rs_qq, Rs_qq_l,
                            Rs_u0, Rs_u0_l, Rs_uq, Rs_uq_l, 
                            btran, btran_l,
                            as_0, as_0_l, atran, atran_l,
                            P_q0_mm, P_q0_pm, P_u0_mm, P_u0_pm, 
                            P_uq_pp, P_uq_mp, P_uq_mm, P_uq_pm, 
                            r_p, t_p, r_m, t_m,
                            P_q0_mm_l, P_q0_pm_l, P_u0_mm_l, P_u0_pm_l,
                            P_uq_pp_l, P_uq_mp_l, P_uq_mm_l, P_uq_pm_l, 
                            r_p_l, t_p_l, r_m_l, t_m_l,
                            I1_m, I1_m_l, In_p, In_p_l,
                            I_p, I_m, I_p_l, I_m_l,
                            sfi, surface, thermal, upwelling,
                            downwelling, utau_output, vector,
                            eigen_solver_real, eigen_solver_complex,
                            derivs_h, derivs_p, save_tree, work);
#endif
}



/*******************************************************************************
 *
 ******************************************************************************/
static void sfi_layer(int i, double ptau, int i_four, int n_quad, int n_stokes, int n_derivs, int n_layers, int n_umus, double qf, double *qx_v, double *qw_v, double F_0, int n_ulevels, int *ulevels, double *utaus, double *umus, double *omega, double **omega_l, double *ltau, double **ltau_l, double *btran, double **btran_l, double *as_0, double **as_0_l, double *atran, double **atran_l, double **P_u0_pm, double ***P_uq_pp, double ***P_uq_mp, double **nu, double ***X_p, double ***X_m, double **F_p, double **F_m, double ***P_u0_pm_l, double ****P_uq_pp_l, double ****P_uq_mp_l, double ***nu_l, double ****X_p_l, double ****X_m_l, double ***F_p_l, double ***F_m_l, double *B, double **B_l, double *I_0, double **I_0_l, double *I_u, int offset, double **I_u_l, int utau_output, uchar **derivs_h, uchar **derivs_p, work_data work, int flag) {

     int ii;
     int j;
     int k;
     int l;

     int n_quad_v;

     double a;
     double b;

     double delfac;
     double solfac;

     double *q1;
     double *q2;
     double *q3;

     double *u1;
     double *u2;
     double *u3;
     double *u4;

     double *A;

     double *e1;
     double **e1_l;
     double e2;
     double e2_l;

     double *ptau_l;

     double *e3;
     double *e4;
     double *e5;

     double *Fu_0;
     double *Du_0;
     double *Gu_0;
     double  Eu_0;

     double **qq1;
     double **qq2;
     double **qq3;
     double **qq4;

     double **uq1;
     double **uq2;

     double **Xu_p;
     double **Xu_m;
     double **Yu_p;
     double **Yu_m;

     double **Du_p;
     double **Du_m;
     double   Eu_p;
     double   Eu_m;


     n_quad_v = n_quad * n_stokes;


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     q1   = get_work1(&work, WORK_DX);
     q2   = get_work1(&work, WORK_DX);
     q3   = get_work1(&work, WORK_DX);

     u1   = get_work1(&work, WORK_DU);
     u2   = get_work1(&work, WORK_DU);
     u3   = get_work1(&work, WORK_DU);
     u4   = get_work1(&work, WORK_DU);

     A    = get_work1(&work, WORK_DX);

     e1   = get_work1(&work, WORK_DU);

     Fu_0 = get_work1(&work, WORK_DU);
     Du_0 = get_work1(&work, WORK_DU);

     e3   = get_work1(&work, WORK_DX);
     e4   = get_work1(&work, WORK_DX);
     e5   = get_work1(&work, WORK_DX);

     qq1  = get_work1(&work, WORK_DXX);
     qq2  = get_work1(&work, WORK_DXX);
     qq3  = get_work1(&work, WORK_DXX);
     qq4  = get_work1(&work, WORK_DXX);

     uq1  = get_work1(&work, WORK_DUX);
     uq2  = get_work1(&work, WORK_DUX);

     Xu_p = get_work1(&work, WORK_DUX);
     Xu_m = get_work1(&work, WORK_DUX);

     Du_p = get_work1(&work, WORK_DUX);
     Du_m = get_work1(&work, WORK_DUX);

     if (n_derivs > 0) {
          ptau_l = get_work1(&work, WORK_DDERIVS);

          e1_l   = get_work2(&work, WORK_DU, WORK_DERIVS_V, NULL);

          Gu_0   = get_work1(&work, WORK_DU);

          Yu_p   = get_work1(&work, WORK_DUX);
          Yu_m   = get_work1(&work, WORK_DUX);
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     delfac = 2. - (i_four == 0 ? 1. : 0.);

     solfac = F_0 / (4. * PI);


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     dvec_scale(1. / (2. * delfac), qw_v, A, n_quad_v);


     for (j = 0; j < n_umus; ++j) {
          if (! flag)
               e1[j] = exp(-(ltau[i] - ptau) / umus[j]);
          else
               e1[j] = exp(-           ptau  / umus[j]);
     }

     for (j = 0; j < n_derivs; ++j) {
          ptau_l[j] = 0.;
          if (derivs_h[i][j])
              ptau_l[j] = ptau * ltau_l[i][j] / ltau[i];

          for (k = 0; k < n_umus; ++k) {
               e1_l[j][k] = 0.;
               if (derivs_h[i][j]) {
                    if (! flag)
                         e1_l[j][k] = -(ltau_l[i][j] - ptau_l[j]) / umus[k] * e1[k];
                    else
                         e1_l[j][k] = -                ptau_l[j]  / umus[k] * e1[k];
               }
          }
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     dm_v_diag_mul(A, F_p[i], q1, n_quad_v);
     dm_v_diag_mul(A, F_m[i], q2, n_quad_v);

     dm_v_mul(P_uq_pp[i], q1, n_umus, n_quad_v, u1);
     dm_v_mul(P_uq_mp[i], q2, n_umus, n_quad_v, u2);
     dvec_add(u1, u2, u1, n_umus);

     a = btran[i] * solfac;

     for (j = 0; j < n_umus; ++j) {
          u1[j] += a * P_u0_pm[i][j];
          Fu_0[j] = omega[i] * u1[j];
     }

     e2 = exp(-ptau * as_0[i]);

     for (j = 0; j < n_umus; ++j) {
          if (! flag)
               Du_0[j] = (e2 - atran[i] * e1[j]) / (1. + umus[j] * as_0[i]);
          else
               Du_0[j] = (e2 -            e1[j]) / (1. - umus[j] * as_0[i]);
     }

     for (j = 0; j < n_derivs; ++j) {
          e2_l = 0.;
          if (derivs_h[i][j])
               e2_l = (-ptau_l[j] * as_0[i] - ptau * as_0_l[i][j]) * e2;

          dvec_zero(Gu_0, n_umus);

          if (derivs_h[i][j]) {
               dm_v_mul(P_uq_pp_l[i][j], q1, n_umus, n_quad_v, u2);
               dm_v_mul(P_uq_mp_l[i][j], q2, n_umus, n_quad_v, u3);
               dvec_add(u2, u3, u2, n_umus);

               for (k = 0; k < n_umus; ++k)
                    Gu_0[k] += omega_l[i][j] * u1[k] + omega[i] * (a * P_u0_pm_l[i][j][k] + u2[k]);
          }

          if (derivs_p[i][j]) {
               dm_v_diag_mul(A, F_p_l[i][j], q3, n_quad_v);
               dm_v_mul(P_uq_pp[i], q3, n_umus, n_quad_v, u3);

               dm_v_diag_mul(A, F_m_l[i][j], q3, n_quad_v);
               dm_v_mul(P_uq_mp[i], q3, n_umus, n_quad_v, u4);
               dvec_add(u3, u4, u3, n_umus);

               b = btran_l[i][j] * solfac;

               for (k = 0; k < n_umus; ++k)
                    Gu_0[k] +=                         omega[i] * (b * P_u0_pm[i]   [k] + u3[k]);
          }

          for (k = 0; k < n_umus; ++k) {
               I_u_l[j][offset+k] = I_0_l[j][k] * e1[k] + Gu_0[k] * Du_0[k];

               if (derivs_h[i][j]) {
                    if (! flag)
                         Eu_0 = ((e2_l - atran_l[i][j] * e1     [k] - atran[i] * e1_l[j][k]) - Du_0[k] * ( umus[k] * as_0_l[i][j])) / (1. + umus[k] * as_0[i]);
                    else
                         Eu_0 = ((e2_l -                 e1_l[j][k])                         - Du_0[k] * (-umus[k] * as_0_l[i][j])) / (1. - umus[k] * as_0[i]);

                    I_u_l[j][offset+k] += I_0[k] * e1_l[j][k] + Fu_0[k] * Eu_0;
               }
          }
     }

     for (j = 0; j < n_umus; ++j)
          I_u[offset+j] = I_0[j] * e1[j] + Fu_0[j] * Du_0[j];


     /*-------------------------------------------------------------------------
      *
      *------------------------------------------------------------------------*/
     ii = i * n_quad_v * 2;

     dmat_diag_mul(A, X_p[i], qq1, n_quad_v, n_quad_v);
     dmat_diag_mul(A, X_m[i], qq2, n_quad_v, n_quad_v);

     dmat_mul(P_uq_pp[i], qq1, n_umus, n_quad_v, n_quad_v, uq1);
     dmat_mul(P_uq_mp[i], qq2, n_umus, n_quad_v, n_quad_v, uq2);
     dmat_sub(uq1, uq2, Xu_p, n_umus, n_quad_v);

     dmat_mul(P_uq_pp[i], qq2, n_umus, n_quad_v, n_quad_v, uq1);
     dmat_mul(P_uq_mp[i], qq1, n_umus, n_quad_v, n_quad_v, uq2);
     dmat_sub(uq1, uq2, Xu_m, n_umus, n_quad_v);

     for (j = 0; j < n_quad_v; ++j) {
          e3[j] = exp(-nu[i][j] *  ptau);
          e4[j] = exp(-nu[i][j] *  ltau[i]);
          e5[j] = exp(-nu[i][j] * (ltau[i] - ptau));
     }

     for (j = 0; j < n_umus; ++j) {
          for (k = 0; k < n_quad_v; ++k) {
               if (! flag) {
                    Du_p[j][k] = (e3[k] - e4[k] * e1[j]) / (1. + umus[j] * nu[i][k]);
                    Du_m[j][k] = (e5[k] -         e1[j]) / (1. - umus[j] * nu[i][k]);
               }
               else {
                    Du_p[j][k] = (e3[k] -         e1[j]) / (1. - umus[j] * nu[i][k]);
                    Du_m[j][k] = (e5[k] - e4[k] * e1[j]) / (1. + umus[j] * nu[i][k]);
               }
          }
     }

     for (j = 0; j < n_derivs; ++j) {
          if (derivs_h[i][j]) {
               dmat_diag_mul(A, X_p_l[i][j], qq3, n_quad_v, n_quad_v);
               dmat_diag_mul(A, X_m_l[i][j], qq4, n_quad_v, n_quad_v);

               dmat_mul(P_uq_pp_l[i][j], qq1, n_umus, n_quad_v, n_quad_v, uq1);
               dmat_mul(P_uq_mp_l[i][j], qq2, n_umus, n_quad_v, n_quad_v, uq2);
               dmat_sub(uq1, uq2, Yu_p, n_umus, n_quad_v);

               dmat_mul(P_uq_pp[i], qq3, n_umus, n_quad_v, n_quad_v, uq1);
               dmat_add(Yu_p, uq1, Yu_p, n_umus, n_quad_v);
               dmat_mul(P_uq_mp[i], qq4, n_umus, n_quad_v, n_quad_v, uq1);
               dmat_sub(Yu_p, uq1, Yu_p, n_umus, n_quad_v);

               dmat_mul(P_uq_pp_l[i][j], qq2, n_umus, n_quad_v, n_quad_v, uq1);
               dmat_mul(P_uq_mp_l[i][j], qq1, n_umus, n_quad_v, n_quad_v, uq2);
               dmat_sub(uq1, uq2, Yu_m, n_umus, n_quad_v);

               dmat_mul(P_uq_pp[i], qq4, n_umus, n_quad_v, n_quad_v, uq1);
               dmat_add(Yu_m, uq1, Yu_m, n_umus, n_quad_v);
               dmat_mul(P_uq_mp[i], qq3, n_umus, n_quad_v, n_quad_v, uq1);
               dmat_sub(Yu_m, uq1, Yu_m, n_umus, n_quad_v);
          }

          for (k = 0; k < n_umus; ++k) {
               for (l = 0; l < n_quad_v; ++l) {
                    I_u_l[j][offset+k] += omega_l[i][j] * (B     [ii + l] * Xu_p[k][l] * Du_p[k][l] + B     [ii + l + n_quad_v] * Xu_m[k][l] * Du_m[k][l]);
                    I_u_l[j][offset+k] += omega  [i]    * (B_l[j][ii + l] * Xu_p[k][l] * Du_p[k][l] + B_l[j][ii + l + n_quad_v] * Xu_m[k][l] * Du_m[k][l]);

                    if (derivs_h[i][j]) {
                         if (! flag) {
                              Eu_p = (((-nu_l[i][j][l] *            ptau  - nu[i][l] *                 ptau_l[j])  * e3[l] - ((-nu_l[i][j][l] * ltau[i] - nu[i][l] * ltau_l[i][j]) * e4[l] * e1[k] + e4[l] * e1_l[j][k])) + Du_p[k][l] * -umus[k] * nu_l[i][j][l]) / (1. + umus[k] * nu[i][l]);
                              Eu_m = (((-nu_l[i][j][l] * (ltau[i] - ptau) - nu[i][l] * (ltau_l[i][j] - ptau_l[j])) * e5[l] -                                                                                 e1_l[j][k])  + Du_m[k][l] *  umus[k] * nu_l[i][j][l]) / (1. - umus[k] * nu[i][l]);
                         }
                         else {
                              Eu_p = (((-nu_l[i][j][l] *            ptau  - nu[i][l] *                 ptau_l[j])  * e3[l] -                                                                                 e1_l[j][k])  + Du_p[k][l] *  umus[k] * nu_l[i][j][l]) / (1. - umus[k] * nu[i][l]);
                              Eu_m = (((-nu_l[i][j][l] * (ltau[i] - ptau) - nu[i][l] * (ltau_l[i][j] - ptau_l[j])) * e5[l] - ((-nu_l[i][j][l] * ltau[i] - nu[i][l] * ltau_l[i][j]) * e4[l] * e1[k] + e4[l] * e1_l[j][k])) + Du_m[k][l] * -umus[k] * nu_l[i][j][l]) / (1. + umus[k] * nu[i][l]);
                         }

                         I_u_l[j][offset+k] += omega[i] * (B[ii + l           ] * Yu_p[k][l] * Du_p[k][l] + B[ii + l           ] * Xu_p[k][l] * Eu_p +
                                                           B[ii + l + n_quad_v] * Yu_m[k][l] * Du_m[k][l] + B[ii + l + n_quad_v] * Xu_m[k][l] * Eu_m);
                    }
               }
          }
     }

     for (j = 0; j < n_umus; ++j) {
          for (k = 0; k < n_quad_v; ++k) {
               I_u[offset+j] += omega[i] * (B[ii + k] * Xu_p[j][k] * Du_p[j][k] + B[ii + k + n_quad_v] * Xu_m[j][k] * Du_m[j][k]);
          }
     }
}



/*******************************************************************************
 *
 ******************************************************************************/
static void sfi(int i_four, int n_quad, int n_stokes, int n_derivs, int n_layers, int n_umus, double qf, double *qx_v, double *qw_v, double F_0, int n_ulevels, int *ulevels, double *utaus, double *umus, double *omega, double **omega_l, double *ltau, double **ltau_l, double *btran, double **btran_l, double *as_0, double **as_0_l, double *atran, double **atran_l, double **P_u0_pm, double ***P_uq_pp, double ***P_uq_mp, double **nu, double ***X_p, double ***X_m, double **F_p, double **F_m, double ***P_u0_pm_l, double ****P_uq_pp_l, double ****P_uq_mp_l, double ***nu_l, double ****X_p_l, double ****X_m_l, double ***F_p_l, double ***F_m_l, double *B, double **B_l, double *I_0, double **I_0_l, double **I_u, int offset, double ***I_u_l, int utau_output, uchar **derivs_h, uchar **derivs_p, work_data work, int flag) {

     int i;
     int i1;
     int i2;
     int j;
     int k;

     int inc;

     int i_out_level;

     double ptau;

     double **I_u_l2;


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (! flag) {
          i1          = n_layers - 1;
          i2          = ulevels[0] - 1;
          inc         = -1;
          i_out_level = n_ulevels - 1;
     }
     else {
          i1          = 0;
          if (! utau_output)
               i2 = ulevels[n_ulevels - 1];
          else
               i2 = ulevels[n_ulevels - 1] + 1;
          inc         =  1;
          i_out_level = 0;
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if ((! flag && ! utau_output && i1 + 1 == ulevels[i_out_level]) || (flag && ! utau_output && i1 == ulevels[i_out_level]) || (! flag && utau_output && i1 + 1 == ulevels[i_out_level] && utaus[i_out_level] == 0.) || (flag && utau_output && i1 == ulevels[i_out_level] && utaus[i_out_level] == 0.)) {
          for (i = 0; i < n_umus; ++i) {
               I_u[i_out_level][offset + i] = I_0[i];
          }
          for (i = 0; i < n_derivs; ++i) {
               for (j = 0; j < n_umus; ++j) {
                    I_u_l[i_out_level][i][offset + j] = I_0_l[i][j];
               }
          }

          i_out_level += inc;
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     for (i = i1; i != i2; i += inc) {

          if ((! flag && utau_output && i == ulevels[i_out_level]) || (flag && utau_output && i == ulevels[i_out_level])) {
               while (i_out_level >= 0 && i_out_level < n_ulevels && i == ulevels[i_out_level]) {
                    ptau = utaus[i_out_level];

                    if (n_derivs > 0)
                         I_u_l2 = I_u_l[i_out_level];

                    sfi_layer(i, ptau, i_four, n_quad, n_stokes, n_derivs, n_layers, n_umus, qf, qx_v, qw_v, F_0, n_ulevels, ulevels, utaus, umus, omega, omega_l, ltau, ltau_l, btran, btran_l, as_0, as_0_l, atran, atran_l, P_u0_pm, P_uq_pp, P_uq_mp, nu, X_p, X_m, F_p, F_m, P_u0_pm_l, P_uq_pp_l, P_uq_mp_l, nu_l, X_p_l, X_m_l, F_p_l, F_m_l, B, B_l, I_0, I_0_l, I_u[i_out_level], offset, I_u_l2, utau_output, derivs_h, derivs_p, work, flag);

                    i_out_level += inc;
               }
          }

          if (! flag)
               ptau = 0.;
          else
               ptau = ltau[i];

          sfi_layer(i, ptau, i_four, n_quad, n_stokes, n_derivs, n_layers, n_umus, qf, qx_v, qw_v, F_0, n_ulevels, ulevels, utaus, umus, omega, omega_l, ltau, ltau_l, btran, btran_l, as_0, as_0_l, atran, atran_l, P_u0_pm, P_uq_pp, P_uq_mp, nu, X_p, X_m, F_p, F_m, P_u0_pm_l, P_uq_pp_l, P_uq_mp_l, nu_l, X_p_l, X_m_l, F_p_l, F_m_l, B, B_l, I_0, I_0_l, I_0, 0, I_0_l, utau_output, derivs_h, derivs_p, work, flag);

          if ((! flag && ! utau_output && i == ulevels[i_out_level]) || (flag && ! utau_output && i + 1 == ulevels[i_out_level])) {
               for (j = 0; j < n_umus; ++j) {
                    I_u[i_out_level][offset + j] = I_0[j];
               }
               for (j = 0; j < n_derivs; ++j) {
                    for (k = 0; k < n_umus; ++k) {
                         I_u_l[i_out_level][j][offset + k] = I_0_l[j][k];
                    }
               }

               i_out_level += inc;
          }
     }
}



/*******************************************************************************
 *
 ******************************************************************************/
void sfi_up(int i_four, int n_quad, int n_stokes, int n_derivs, int n_layers, int n_umus, double qf, double *qx_v, double *qw_v, double F_0, double mu_0, int n_ulevels, int *ulevels, double *utaus, double *umus, double *omega, double **omega_l, double *ltau, double **ltau_l, double *Rs_u0, double **Rs_u0_l, double **Rs_uq, double ***Rs_uq_l, double *btran, double **btran_l, double *as_0, double **as_0_l, double *atran, double **atran_l, double **P_u0_pm, double ***P_uq_pp, double ***P_uq_mp, double ***P_uq_mm, double ***P_uq_pm, double **nu, double ***X_p, double ***X_m, double **F_p, double **F_m, double ***P_u0_pm_l, double ****P_uq_pp_l, double ****P_uq_mp_l, double ****P_uq_mm_l, double ****P_uq_pm_l, double ***nu_l, double ****X_p_l, double ****X_m_l, double ***F_p_l, double ***F_m_l, double *B, double **B_l, double **I_m, double ***K_m, double **I_u, int offset, double ***I_u_l, int surface, int utau_output, uchar **derivs_h, uchar **derivs_p, work_data work) {

     int i;
     int j;

     int n_quad_v;

     double a;
     double b;
     double c;

     double *u1;

     double *I_0;
     double **I_0_l;


     n_quad_v = n_quad * n_stokes;


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     u1  = get_work1(&work, WORK_DU);

     I_0 = get_work1(&work, WORK_DU);

     if (n_derivs > 0)
          I_0_l = get_work2(&work, WORK_DU, WORK_DERIVS_V, NULL);


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (surface) {
          dm_v_mul(Rs_uq, I_m[0], n_umus, n_quad_v, I_0);

          a = qf * mu_0 * F_0 / PI;

          b = a * btran[n_layers];

          dvec_scale(b, Rs_u0, u1, n_umus);
          dvec_add(I_0, u1, I_0, n_umus);

          for (i = 0; i < n_derivs; ++i) {
               init_array1_d(I_0_l[i], n_umus, 0.);

               if (derivs_h[n_layers][i]) {
                    dm_v_mul(Rs_uq_l[i], I_m[0], n_umus, n_quad_v, u1);
                    dvec_add(I_0_l[i], u1, I_0_l[i], n_umus);

                    for (j = 0; j < n_umus; ++j)
                         I_0_l[i][j] += b * Rs_u0_l[i][j];
/*
                    dvec_scale(b, Rs_u0_l[i], u1, n_umus);
                    dvec_add(I_0_l[i], u1, I_0_l[i], n_umus);
*/
               }
               if (derivs_p[n_layers][i]) {
                    dm_v_mul(Rs_uq, K_m[0][i], n_umus, n_quad_v, u1);
                    dvec_add(I_0_l[i], u1, I_0_l[i], n_umus);

                    c = a * btran_l[n_layers][i];

                    for (j = 0; j < n_umus; ++j)
                         I_0_l[i][j] += c * Rs_u0   [j];
/*
                    dvec_scale(c, Rs_u0,    u1, n_umus);
                    dvec_add(I_u_l[i], u1, I_u_l[i], n_umus);
*/
               }
          }
     }
     else {
          for (i = 0; i < n_umus; ++i)
               I_0[i] = 0.;

          for (i = 0; i < n_derivs; ++i) {
               for (j = 0; j < n_umus; ++j) {
                    I_0_l[i][j] = 0.;
               }
          }
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     sfi(i_four, n_quad, n_stokes, n_derivs, n_layers, n_umus, qf, qx_v, qw_v, F_0, n_ulevels, ulevels, utaus, umus, omega, omega_l, ltau, ltau_l, btran, btran_l, as_0, as_0_l, atran, atran_l, P_u0_pm, P_uq_pp, P_uq_mp, nu, X_p, X_m, F_p, F_m, P_u0_pm_l, P_uq_pp_l, P_uq_mp_l, nu_l, X_p_l, X_m_l, F_p_l, F_m_l, B, B_l, I_0, I_0_l, I_u, offset, I_u_l, utau_output, derivs_h, derivs_p, work, 0);
}



/*******************************************************************************
 *
 ******************************************************************************/
void sfi_dn(int i_four, int n_quad, int n_stokes, int n_derivs, int n_layers, int n_umus, double qf, double *qx_v, double *qw_v, double F_0, int n_ulevels, int *ulevels, double *utaus, double *umus, double *omega, double **omega_l, double *ltau, double **ltau_l, double *btran, double **btran_l, double *as_0, double **as_0_l, double *atran, double **atran_l, double **P_u0_mm, double ***P_uq_pp, double ***P_uq_mp, double ***P_uq_mm, double ***P_uq_pm, double **nu, double ***X_p, double ***X_m, double **F_p, double **F_m, double ***P_u0_mm_l, double ****P_uq_pp_l, double ****P_uq_mp_l, double ****P_uq_mm_l, double ****P_uq_pm_l, double ***nu_l, double ****X_p_l, double ****X_m_l, double ***F_p_l, double ***F_m_l, double *B, double **B_l, double **I_m, double ***K_m, double **I_u, int offset, double ***I_u_l, int utau_output, uchar **derivs_h, uchar **derivs_p, work_data work) {

     int i;
     int j;

     double *I_0;
     double **I_0_l;


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     I_0 = get_work1(&work, WORK_DU);

     if (n_derivs > 0)
          I_0_l = get_work2(&work, WORK_DU, WORK_DERIVS_V, NULL);


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     for (i = 0; i < n_umus; ++i)
          I_0[i] = 0.;

     for (i = 0; i < n_derivs; ++i) {
          for (j = 0; j < n_umus; ++j) {
               I_0_l[i][j] = 0.;
          }
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     sfi(i_four, n_quad, n_stokes, n_derivs, n_layers, n_umus, qf, qx_v, qw_v, F_0, n_ulevels, ulevels, utaus, umus, omega, omega_l, ltau, ltau_l, btran, btran_l, as_0, as_0_l, atran, atran_l, P_u0_mm, P_uq_mp, P_uq_pp, nu, X_p, X_m, F_p, F_m, P_u0_mm_l, P_uq_mp_l, P_uq_pp_l, nu_l, X_p_l, X_m_l, F_p_l, F_m_l, B, B_l, I_0, I_0_l, I_u, offset, I_u_l, utau_output, derivs_h, derivs_p, work, 1);
}



/*******************************************************************************
 *
 ******************************************************************************/
static void sfi_layer2(int i, double ptau, int i_four, int n_quad, int n_stokes, int n_derivs, int n_layers, int n_umus, double qf, double *qx_v, double *qw_v, double F_0, int n_ulevels, int *ulevels, double *utaus, double *umus, double *omega, double **omega_l, double *ltau, double **ltau_l, double *btran, double **btran_l, double *as_0, double **as_0_l, double *atran, double **atran_l, double **P_u0_pm, double ***P_uq_pp, double ***P_uq_mp, dcomplex **nu, dcomplex ***X_p, dcomplex ***X_m, double **F_p, double **F_m, double ***P_u0_pm_l, double ****P_uq_pp_l, double ****P_uq_mp_l, dcomplex ***nu_l, dcomplex ****X_p_l, dcomplex ****X_m_l, double ***F_p_l, double ***F_m_l, dcomplex *B, dcomplex **B_l, double *I_0, double **I_0_l, double *I_u, int offset, double **I_u_l, int utau_output, uchar **derivs_h, uchar **derivs_p, work_data work, int flag) {

     int ii;
     int j;
     int k;
     int l;

     int n_quad_v;
     int n_umus_v;

     double a;
     double b;

     double delfac;
     double solfac;

     double *dq1;
     double *dq2;
     double *dq3;

     double *du1;
     double *du2;
     double *du3;
     double *du4;

     double *A;

     double *e1;
     double **e1_l;
     double e2;
     double e2_l;

     double *ptau_l;

     double **dqq1;

     double **duq1;
     double **duq2;

     dcomplex *e3;
     dcomplex *e4;
     dcomplex *e5;

     double *Fu_0;
     double *Du_0;
     double *Gu_0;
     double  Eu_0;

     dcomplex **zqq1;
     dcomplex **zqq2;
     dcomplex **zqq3;
     dcomplex **zqq4;

     dcomplex **zuq1;
     dcomplex **zuq2;

     dcomplex **Xu_p;
     dcomplex **Xu_m;
     dcomplex **Yu_p;
     dcomplex **Yu_m;

     dcomplex **Du_p;
     dcomplex **Du_m;
     dcomplex   Eu_p;
     dcomplex   Eu_m;


     n_quad_v = n_quad * n_stokes;

     n_umus_v = n_umus  * n_stokes;


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     dq1  = get_work1(&work, WORK_DX);
     dq2  = get_work1(&work, WORK_DX);
     dq3  = get_work1(&work, WORK_DX);

     du1  = get_work1(&work, WORK_DU);
     du2  = get_work1(&work, WORK_DU);
     du3  = get_work1(&work, WORK_DU);
     du4  = get_work1(&work, WORK_DU);

     A    = get_work1(&work, WORK_DX);

     e1   = get_work1(&work, WORK_DU);

     Fu_0 = get_work1(&work, WORK_DU);
     Du_0 = get_work1(&work, WORK_DU);

     e3   = get_work1(&work, WORK_ZX);
     e4   = get_work1(&work, WORK_ZX);
     e5   = get_work1(&work, WORK_ZX);

     dqq1 = get_work1(&work, WORK_DXX);

     duq1 = get_work1(&work, WORK_DUX);
     duq2 = get_work1(&work, WORK_DUX);

     zqq1 = get_work1(&work, WORK_ZXX);
     zqq2 = get_work1(&work, WORK_ZXX);
     zqq3 = get_work1(&work, WORK_ZXX);
     zqq4 = get_work1(&work, WORK_ZXX);

     zuq1 = get_work1(&work, WORK_ZUX);
     zuq2 = get_work1(&work, WORK_ZUX);

     Xu_p = get_work1(&work, WORK_ZUX);
     Xu_m = get_work1(&work, WORK_ZUX);

     Du_p = get_work1(&work, WORK_ZUX);
     Du_m = get_work1(&work, WORK_ZUX);

     if (n_derivs > 0) {
          ptau_l = get_work1(&work, WORK_DDERIVS);

          e1_l   = get_work2(&work, WORK_DU, WORK_DERIVS_V, NULL);

          Gu_0   = get_work1(&work, WORK_DU);

          Yu_p   = get_work1(&work, WORK_ZUX);
          Yu_m   = get_work1(&work, WORK_ZUX);
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     delfac = 2. - (i_four == 0 ? 1. : 0.);

     solfac = F_0 / (4. * PI);


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     dvec_scale(1. / (2. * delfac), qw_v, A, n_quad_v);


     for (j = 0; j < n_umus_v; ++j) {
          if (! flag)
               e1[j] = exp(-(ltau[i] - ptau) / umus[j]);
          else
               e1[j] = exp(-           ptau  / umus[j]);
     }

     for (j = 0; j < n_derivs; ++j) {
          ptau_l[j] = 0.;
          if (derivs_h[i][j])
              ptau_l[j] = ptau * ltau_l[i][j] / ltau[i];

          for (k = 0; k < n_umus_v; ++k) {
               e1_l[j][k] = 0.;
               if (derivs_h[i][j]) {
                    if (! flag)
                         e1_l[j][k] = -(ltau_l[i][j] - ptau_l[j]) / umus[k] * e1[k];
                    else
                         e1_l[j][k] = -                ptau_l[j]  / umus[k] * e1[k];
               }
          }
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     dm_v_diag_mul(A, F_p[i], dq1, n_quad_v);
     dm_v_diag_mul(A, F_m[i], dq2, n_quad_v);

     dm_v_mul(P_uq_pp[i], dq1, n_umus_v, n_quad_v, du1);
     dm_v_mul(P_uq_mp[i], dq2, n_umus_v, n_quad_v, du2);
     dvec_add(du1, du2, du1, n_umus_v);

     a = btran[i] * solfac;

     for (j = 0; j < n_umus_v; ++j) {
          du1[j] += a * P_u0_pm[i][j];
          Fu_0[j] = omega[i] * du1[j];
     }

     e2 = exp(-ptau * as_0[i]);

     for (j = 0; j < n_umus_v; ++j) {
          if (! flag)
               Du_0[j] = (e2 - atran[i] * e1[j]) / (1. + umus[j] * as_0[i]);
          else
               Du_0[j] = (e2 -            e1[j]) / (1. - umus[j] * as_0[i]);
     }

     for (j = 0; j < n_derivs; ++j) {
          e2_l = 0.;
          if (derivs_h[i][j])
               e2_l = (-ptau_l[j] * as_0[i] - ptau * as_0_l[i][j]) * e2;

          dvec_zero(Gu_0, n_umus_v);

          if (derivs_h[i][j]) {
               dm_v_mul(P_uq_pp_l[i][j], dq1, n_umus_v, n_quad_v, du2);
               dm_v_mul(P_uq_mp_l[i][j], dq2, n_umus_v, n_quad_v, du3);
               dvec_add(du2, du3, du2, n_umus_v);

               for (k = 0; k < n_umus_v; ++k)
                    Gu_0[k] += omega_l[i][j] * du1[k] + omega[i] * (a * P_u0_pm_l[i][j][k] + du2[k]);
          }

          if (derivs_p[i][j]) {
               dm_v_diag_mul(A, F_p_l[i][j], dq3, n_quad_v);
               dm_v_mul(P_uq_pp[i], dq3, n_umus_v, n_quad_v, du3);

               dm_v_diag_mul(A, F_m_l[i][j], dq3, n_quad_v);
               dm_v_mul(P_uq_mp[i], dq3, n_umus_v, n_quad_v, du4);
               dvec_add(du3, du4, du3, n_umus_v);

               b = btran_l[i][j] * solfac;

               for (k = 0; k < n_umus_v; ++k)
                    Gu_0[k] +=                         omega[i] * (b * P_u0_pm[i]   [k] + du3[k]);
          }

          for (k = 0; k < n_umus_v; ++k) {
               I_u_l[j][offset+k] = I_0_l[j][k] * e1[k] + Gu_0[k] * Du_0[k];

               if (derivs_h[i][j]) {
                    if (! flag)
                         Eu_0 = ((e2_l - atran_l[i][j] * e1     [k] - atran[i] * e1_l[j][k]) - Du_0[k] * ( umus[k] * as_0_l[i][j])) / (1. + umus[k] * as_0[i]);
                    else
                         Eu_0 = ((e2_l -                 e1_l[j][k])                         - Du_0[k] * (-umus[k] * as_0_l[i][j])) / (1. - umus[k] * as_0[i]);

                    I_u_l[j][offset+k] += I_0[k] * e1_l[j][k] + Fu_0[k] * Eu_0;
               }
          }
     }

     for (j = 0; j < n_umus_v; ++j)
          I_u[offset+j] = I_0[j] * e1[j] + Fu_0[j] * Du_0[j];


     /*--------------------------------------------------------------------
      *
      *-------------------------------------------------------------------*/
     ii = i * n_quad_v * 2;

     dzmat_diag_mul(A, X_p[i], zqq1, n_quad_v, n_quad_v);
     dzmat_diag_mul(A, X_m[i], zqq2, n_quad_v, n_quad_v);

     dzmat_mul2(P_uq_pp[i], zqq1, n_umus_v, n_quad_v, n_quad_v, zuq1, dqq1, duq1, duq2);
     dzmat_mul2(P_uq_mp[i], zqq2, n_umus_v, n_quad_v, n_quad_v, zuq2, dqq1, duq1, duq2);
     zmat_sub(zuq1, zuq2, Xu_p, n_umus_v, n_quad_v);

     dzmat_mul2(P_uq_pp[i], zqq2, n_umus_v, n_quad_v, n_quad_v, zuq1, dqq1, duq1, duq2);
     dzmat_mul2(P_uq_mp[i], zqq1, n_umus_v, n_quad_v, n_quad_v, zuq2, dqq1, duq1, duq2);
     zmat_sub(zuq1, zuq2, Xu_m, n_umus_v, n_quad_v);

     for (j = 0; j < n_quad_v; ++j) {
          e3[j] = cexp(-nu[i][j] *  ptau);
          e4[j] = cexp(-nu[i][j] *  ltau[i]);
          e5[j] = cexp(-nu[i][j] * (ltau[i] - ptau));
     }

     for (j = 0; j < n_umus_v; ++j) {
          for (k = 0; k < n_quad_v; ++k) {
               if (! flag) {
                    Du_p[j][k] = (e3[k] - e4[k] * e1[j]) / (1. + umus[j] * nu[i][k]);
                    Du_m[j][k] = (e5[k] -         e1[j]) / (1. - umus[j] * nu[i][k]);
               }
               else {
                    Du_p[j][k] = (e3[k] -         e1[j]) / (1. - umus[j] * nu[i][k]);
                    Du_m[j][k] = (e5[k] - e4[k] * e1[j]) / (1. + umus[j] * nu[i][k]);
               }
          }
     }

     for (j = 0; j < n_derivs; ++j) {
          if (derivs_h[i][j]) {
               dzmat_diag_mul(A, X_p_l[i][j], zqq3, n_quad_v, n_quad_v);
               dzmat_diag_mul(A, X_m_l[i][j], zqq4, n_quad_v, n_quad_v);

               dzmat_mul2(P_uq_pp_l[i][j], zqq1, n_umus_v, n_quad_v, n_quad_v, zuq1, dqq1, duq1, duq2);
               dzmat_mul2(P_uq_mp_l[i][j], zqq2, n_umus_v, n_quad_v, n_quad_v, zuq2, dqq1, duq1, duq2);
               zmat_sub(zuq1, zuq2, Yu_p, n_umus_v, n_quad_v);

               dzmat_mul2(P_uq_pp[i], zqq3, n_umus_v, n_quad_v, n_quad_v, zuq1, dqq1, duq1, duq2);
               zmat_add(Yu_p, zuq1, Yu_p, n_umus_v, n_quad_v);
               dzmat_mul2(P_uq_mp[i], zqq4, n_umus_v, n_quad_v, n_quad_v, zuq1, dqq1, duq1, duq2);
               zmat_sub(Yu_p, zuq1, Yu_p, n_umus_v, n_quad_v);

               dzmat_mul2(P_uq_pp_l[i][j], zqq2, n_umus_v, n_quad_v, n_quad_v, zuq1, dqq1, duq1, duq2);
               dzmat_mul2(P_uq_mp_l[i][j], zqq1, n_umus_v, n_quad_v, n_quad_v, zuq2, dqq1, duq1, duq2);
               zmat_sub(zuq1, zuq2, Yu_m, n_umus_v, n_quad_v);

               dzmat_mul2(P_uq_pp[i], zqq4, n_umus_v, n_quad_v, n_quad_v, zuq1, dqq1, duq1, duq2);
               zmat_add(Yu_m, zuq1, Yu_m, n_umus_v, n_quad_v);
               dzmat_mul2(P_uq_mp[i], zqq3, n_umus_v, n_quad_v, n_quad_v, zuq1, dqq1, duq1, duq2);
               zmat_sub(Yu_m, zuq1, Yu_m, n_umus_v, n_quad_v);
          }

          for (k = 0; k < n_umus_v; ++k) {
               for (l = 0; l < n_quad_v; ++l) {
                    I_u_l[j][offset+k] += omega_l[i][j] * creal(B     [ii + l] * Xu_p[k][l] * Du_p[k][l] + B     [ii + l + n_quad_v] * Xu_m[k][l] * Du_m[k][l]);
                    I_u_l[j][offset+k] += omega  [i]    * creal(B_l[j][ii + l] * Xu_p[k][l] * Du_p[k][l] + B_l[j][ii + l + n_quad_v] * Xu_m[k][l] * Du_m[k][l]);

                    if (derivs_h[i][j]) {
                         if (! flag) {
                              Eu_p = (((-nu_l[i][j][l] *            ptau  - nu[i][l] *                 ptau_l[j])  * e3[l] - ((-nu_l[i][j][l] * ltau[i] - nu[i][l] * ltau_l[i][j]) * e4[l] * e1[k] + e4[l] * e1_l[j][k])) + Du_p[k][l] * -umus[k] * nu_l[i][j][l]) / (1. + umus[k] * nu[i][l]);
                              Eu_m = (((-nu_l[i][j][l] * (ltau[i] - ptau) - nu[i][l] * (ltau_l[i][j] - ptau_l[j])) * e5[l] -                                                                                 e1_l[j][k])  + Du_m[k][l] *  umus[k] * nu_l[i][j][l]) / (1. - umus[k] * nu[i][l]);
                         }
                         else {
                              Eu_p = (((-nu_l[i][j][l] *            ptau  - nu[i][l] *                 ptau_l[j])  * e3[l] -                                                                                 e1_l[j][k])  + Du_p[k][l] *  umus[k] * nu_l[i][j][l]) / (1. - umus[k] * nu[i][l]);
                              Eu_m = (((-nu_l[i][j][l] * (ltau[i] - ptau) - nu[i][l] * (ltau_l[i][j] - ptau_l[j])) * e5[l] - ((-nu_l[i][j][l] * ltau[i] - nu[i][l] * ltau_l[i][j]) * e4[l] * e1[k] + e4[l] * e1_l[j][k])) + Du_m[k][l] * -umus[k] * nu_l[i][j][l]) / (1. + umus[k] * nu[i][l]);
                         }

                         I_u_l[j][offset+k] += omega[i] * creal(B[ii + l           ] * Yu_p[k][l] * Du_p[k][l] + B[ii + l           ] * Xu_p[k][l] * Eu_p +
                                                                B[ii + l + n_quad_v] * Yu_m[k][l] * Du_m[k][l] + B[ii + l + n_quad_v] * Xu_m[k][l] * Eu_m);
                    }
               }
          }
     }

     for (j = 0; j < n_umus_v; ++j) {
          for (k = 0; k < n_quad_v; ++k) {
               I_u[offset+j] += omega[i] * creal(B[ii + k] * Xu_p[j][k] * Du_p[j][k] + B[ii + k + n_quad_v] * Xu_m[j][k] * Du_m[j][k]);
          }
     }
}



/*******************************************************************************
 *
 ******************************************************************************/
static void sfi2(int i_four, int n_quad, int n_stokes, int n_derivs, int n_layers, int n_umus, double qf, double *qx_v, double *qw_v, double F_0, int n_ulevels, int *ulevels, double *utaus, double *umus, double *omega, double **omega_l, double *ltau, double **ltau_l, double *btran, double **btran_l, double *as_0, double **as_0_l, double *atran, double **atran_l, double **P_u0_pm, double ***P_uq_pp, double ***P_uq_mp, dcomplex **nu, dcomplex ***X_p, dcomplex ***X_m, double **F_p, double **F_m, double ***P_u0_pm_l, double ****P_uq_pp_l, double ****P_uq_mp_l, dcomplex ***nu_l, dcomplex ****X_p_l, dcomplex ****X_m_l, double ***F_p_l, double ***F_m_l, dcomplex *B, dcomplex **B_l, double *I_0, double **I_0_l, double **I_u, int offset, double ***I_u_l, int utau_output, uchar **derivs_h, uchar **derivs_p, work_data work, int flag) {

     int i;
     int i1;
     int i2;
     int j;
     int k;

     int n_umus_v;

     int inc;

     int i_out_level;

     double ptau;

     double **I_u_l2;


     n_umus_v = n_umus  * n_stokes;


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (! flag) {
          i1          = n_layers - 1;
          i2          = ulevels[0] - 1;
          inc         = -1;
          i_out_level = n_ulevels - 1;
     }
     else {
          i1          = 0;
          if (! utau_output)
               i2 = ulevels[n_ulevels - 1];
          else
               i2 = ulevels[n_ulevels - 1] + 1;
          inc         =  1;
          i_out_level = 0;
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if ((! flag && ! utau_output && i1 + 1 == ulevels[i_out_level]) || (flag && ! utau_output && i1 == ulevels[i_out_level]) || (! flag && utau_output && i1 + 1 == ulevels[i_out_level] && utaus[i_out_level] == 0.) || (flag && utau_output && i1 == ulevels[i_out_level] && utaus[i_out_level] == 0.)) {
          for (i = 0; i < n_umus_v; ++i) {
               I_u[i_out_level][offset + i] = I_0[i];
          }
          for (i = 0; i < n_derivs; ++i) {
               for (j = 0; j < n_umus_v; ++j) {
                    I_u_l[i_out_level][i][offset + j] = I_0_l[i][j];
               }
          }

          i_out_level += inc;
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     for (i = i1; i != i2; i += inc) {
          if ((! flag && utau_output && i == ulevels[i_out_level]) || (flag && utau_output && i == ulevels[i_out_level])) {
               while (i_out_level >= 0 && i_out_level < n_ulevels && i == ulevels[i_out_level]) {
                    ptau = utaus[i_out_level];

                    if (n_derivs > 0)
                         I_u_l2 = I_u_l[i_out_level];

                    sfi_layer2(i, ptau, i_four, n_quad, n_stokes, n_derivs, n_layers, n_umus, qf, qx_v, qw_v, F_0, n_ulevels, ulevels, utaus, umus, omega, omega_l, ltau, ltau_l, btran, btran_l, as_0, as_0_l, atran, atran_l, P_u0_pm, P_uq_pp, P_uq_mp, nu, X_p, X_m, F_p, F_m, P_u0_pm_l, P_uq_pp_l, P_uq_mp_l, nu_l, X_p_l, X_m_l, F_p_l, F_m_l, B, B_l, I_0, I_0_l, I_u[i_out_level], offset, I_u_l2, utau_output, derivs_h, derivs_p, work, flag);

                    i_out_level += inc;
               }
          }

          if (! flag)
               ptau = 0.;
          else
               ptau = ltau[i];

          sfi_layer2(i, ptau, i_four, n_quad, n_stokes, n_derivs, n_layers, n_umus, qf, qx_v, qw_v, F_0, n_ulevels, ulevels, utaus, umus, omega, omega_l, ltau, ltau_l, btran, btran_l, as_0, as_0_l, atran, atran_l, P_u0_pm, P_uq_pp, P_uq_mp, nu, X_p, X_m, F_p, F_m, P_u0_pm_l, P_uq_pp_l, P_uq_mp_l, nu_l, X_p_l, X_m_l, F_p_l, F_m_l, B, B_l, I_0, I_0_l, I_0, 0, I_0_l, utau_output, derivs_h, derivs_p, work, flag);

          if ((! flag && ! utau_output && i == ulevels[i_out_level]) || (flag && ! utau_output && i + 1 == ulevels[i_out_level])) {
               for (j = 0; j < n_umus_v; ++j) {
                    I_u[i_out_level][offset + j] = I_0[j];
               }
               for (j = 0; j < n_derivs; ++j) {
                    for (k = 0; k < n_umus_v; ++k) {
                         I_u_l[i_out_level][j][offset + k] = I_0_l[j][k];
                    }
               }

               i_out_level += inc;
          }
     }
}



/*******************************************************************************
 *
 ******************************************************************************/
void sfi_up2(int i_four, int n_quad, int n_stokes, int n_derivs, int n_layers, int n_umus, double qf, double *qx_v, double *qw_v, double F_0, double mu_0, int n_ulevels, int *ulevels, double *utaus, double *umus, double *omega, double **omega_l, double *ltau, double **ltau_l, double *Rs_u0, double **Rs_u0_l, double **Rs_uq, double ***Rs_uq_l, double *btran, double **btran_l, double *as_0, double **as_0_l, double *atran, double **atran_l, double **P_u0_pm, double ***P_uq_pp, double ***P_uq_mp, double ***P_uq_mm, double ***P_uq_pm, dcomplex **nu, dcomplex ***X_p, dcomplex ***X_m, double **F_p, double **F_m, double ***P_u0_pm_l, double ****P_uq_pp_l, double ****P_uq_mp_l, double ****P_uq_mm_l, double ****P_uq_pm_l, dcomplex ***nu_l, dcomplex ****X_p_l, dcomplex ****X_m_l, double ***F_p_l, double ***F_m_l, dcomplex *B, dcomplex **B_l, double **I_m, double ***K_m, double **I_u, int offset, double ***I_u_l, int surface, int utau_output, uchar **derivs_h, uchar **derivs_p, work_data work) {

     int i;
     int j;

     int n_quad_v;
     int n_umus_v;

     double a;
     double b;
     double c;

     double *u1;

     double *I_0;
     double **I_0_l;


     n_quad_v = n_quad * n_stokes;

     n_umus_v = n_umus  * n_stokes;


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     u1 = get_work1(&work, WORK_DU);

     I_0 = get_work1(&work, WORK_DU);

     if (n_derivs > 0)
          I_0_l = get_work2(&work, WORK_DU, WORK_DERIVS_V, NULL);


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (surface) {
          dm_v_mul(Rs_uq, I_m[0], n_umus_v, n_quad_v, I_0);

          a = qf * mu_0 * F_0 / PI;

          b = a * btran[n_layers];

          dvec_scale(b, Rs_u0, u1, n_umus_v);
          dvec_add(I_0, u1, I_0, n_umus_v);

          for (i = 0; i < n_derivs; ++i) {
               init_array1_d(I_0_l[i], n_umus_v, 0.);

               if (derivs_h[n_layers][i]) {
                    dm_v_mul(Rs_uq_l[i], I_m[0], n_umus_v, n_quad_v, u1);
                    dvec_add(I_0_l[i], u1, I_0_l[i], n_umus_v);

                    for (j = 0; j < n_umus_v; ++j)
                         I_0_l[i][j] += b * Rs_u0_l[i][j];
/*
                    dvec_scale(b, Rs_u0_l[i], u1, n_umus_v);
                    dvec_add(I_0_l[i], u1, I_0_l[i], n_umus_v);
*/
               }
               if (derivs_p[n_layers][i]) {
                    dm_v_mul(Rs_uq, K_m[0][i], n_umus_v, n_quad_v, u1);
                    dvec_add(I_0_l[i], u1, I_0_l[i], n_umus_v);

                    c = a * btran_l[n_layers][i];

                    for (j = 0; j < n_umus_v; ++j)
                         I_0_l[i][j] += c * Rs_u0   [j];
/*
                    dvec_scale(c, Rs_u0,    u1, n_umus_v);
                    dvec_add(I_0_l[i], u1, I_0_l[i], n_umus_v);
*/
               }
          }
     }
     else {
          for (i = 0; i < n_umus_v; ++i)
               I_0[i] = 0.;

          for (i = 0; i < n_derivs; ++i) {
               for (j = 0; j < n_umus_v; ++j) {
                    I_0_l[i][j] = 0.;
               }
          }
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     for (i = 0; i < n_layers; ++i) {
           dmat_mul_D_A2(n_umus, n_stokes, n_quad, n_stokes, P_uq_mp[i], P_uq_mp[i]);

           for (j = 0; j < n_derivs; ++j) {
                if (derivs_h[i][j]) {
                     dmat_mul_D_A2(n_umus, n_stokes, n_quad, n_stokes, P_uq_mp_l[i][j], P_uq_mp_l[i][j]);
                }
           }
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     sfi2(i_four, n_quad, n_stokes, n_derivs, n_layers, n_umus, qf, qx_v, qw_v, F_0, n_ulevels, ulevels, utaus, umus, omega, omega_l, ltau, ltau_l, btran, btran_l, as_0, as_0_l, atran, atran_l, P_u0_pm, P_uq_pp, P_uq_mp, nu, X_p, X_m, F_p, F_m, P_u0_pm_l, P_uq_pp_l, P_uq_mp_l, nu_l, X_p_l, X_m_l, F_p_l, F_m_l, B, B_l, I_0, I_0_l, I_u, offset, I_u_l, utau_output, derivs_h, derivs_p, work, 0);


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     for (i = 0; i < n_layers; ++i) {
           dmat_mul_D_A2(n_umus, n_stokes, n_quad, n_stokes, P_uq_mp[i], P_uq_mp[i]);

           for (j = 0; j < n_derivs; ++j) {
                if (derivs_h[i][j]) {
                     dmat_mul_D_A2(n_umus, n_stokes, n_quad, n_stokes, P_uq_mp_l[i][j], P_uq_mp_l[i][j]);
                }
           }
     }
}



/*******************************************************************************
 *
 ******************************************************************************/
void sfi_dn2(int i_four, int n_quad, int n_stokes, int n_derivs, int n_layers, int n_umus, double qf, double *qx_v, double *qw_v, double F_0, int n_ulevels, int *ulevels, double *utaus, double *umus, double *omega, double **omega_l, double *ltau, double **ltau_l, double *btran, double **btran_l, double *as_0, double **as_0_l, double *atran, double **atran_l, double **P_u0_mm, double ***P_uq_pp, double ***P_uq_mp, double ***P_uq_mm, double ***P_uq_pm, dcomplex **nu, dcomplex ***X_p, dcomplex ***X_m, double **F_p, double **F_m, double ***P_u0_mm_l, double ****P_uq_pp_l, double ****P_uq_mp_l, double ****P_uq_mm_l, double ****P_uq_pm_l, dcomplex ***nu_l, dcomplex ****X_p_l, dcomplex ****X_m_l, double ***F_p_l, double ***F_m_l, dcomplex *B, dcomplex **B_l, double **I_m, double ***K_m, double **I_u, int offset, double ***I_u_l, int utau_output, uchar **derivs_h, uchar **derivs_p, work_data work) {

     int i;
     int j;

     int n_umus_v;

     double *I_0;
     double **I_0_l;


     n_umus_v = n_umus  * n_stokes;


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     I_0 = get_work1(&work, WORK_DU);

     if (n_derivs > 0)
          I_0_l = get_work2(&work, WORK_DU, WORK_DERIVS_V, NULL);


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     for (i = 0; i < n_umus_v; ++i)
          I_0[i] = 0.;

     for (i = 0; i < n_derivs; ++i) {
          for (j = 0; j < n_umus_v; ++j) {
               I_0_l[i][j] = 0.;
          }
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     for (i = 0; i < n_layers; ++i) {
           dm_v_mul_D_A (n_umus, n_stokes,                   P_u0_mm[i], P_u0_mm[i]);
           dmat_mul_D_A2(n_umus, n_stokes, n_quad, n_stokes, P_uq_mp[i], P_uq_mp[i]);

           for (j = 0; j < n_derivs; ++j) {
                if (derivs_h[i][j]) {
                     dm_v_mul_D_A (n_umus, n_stokes,                   P_u0_mm_l[i][j], P_u0_mm_l[i][j]);
                     dmat_mul_D_A2(n_umus, n_stokes, n_quad, n_stokes, P_uq_mp_l[i][j], P_uq_mp_l[i][j]);
                }
           }
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     sfi2(i_four, n_quad, n_stokes, n_derivs, n_layers, n_umus, qf, qx_v, qw_v, F_0, n_ulevels, ulevels, utaus, umus, omega, omega_l, ltau, ltau_l, btran, btran_l, as_0, as_0_l, atran, atran_l, P_u0_mm, P_uq_mp, P_uq_pp, nu, X_p, X_m, F_p, F_m, P_u0_mm_l, P_uq_mp_l, P_uq_pp_l, nu_l, X_p_l, X_m_l, F_p_l, F_m_l, B, B_l, I_0, I_0_l, I_u, offset, I_u_l, utau_output, derivs_h, derivs_p, work, 1);


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     for (i = 0; i < n_layers; ++i) {
           dm_v_mul_D_A (n_umus, n_stokes,                   P_u0_mm[i], P_u0_mm[i]);
           dmat_mul_D_A2(n_umus, n_stokes, n_quad, n_stokes, P_uq_mp[i], P_uq_mp[i]);

           for (j = 0; j < n_derivs; ++j) {
                if (derivs_h[i][j]) {
                     dm_v_mul_D_A (n_umus, n_stokes,                   P_u0_mm_l[i][j], P_u0_mm_l[i][j]);
                     dmat_mul_D_A2(n_umus, n_stokes, n_quad, n_stokes, P_uq_mp_l[i][j], P_uq_mp_l[i][j]);
                }
           }
     }
}
