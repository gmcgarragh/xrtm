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
#include "xrtm_eig_bvp.h"
#include "xrtm_eig_bvp_a.h"
#include "xrtm_eig_util.h"
#include "xrtm_matrix.h"
#include "xrtm_save_tree.h"
#include "xrtm_sfi.h"
#include "xrtm_source.h"
#include "xrtm_utility.h"


#define LEGACY_MODE 1


/*******************************************************************************
 *
 ******************************************************************************/
#define REAL

#define TYPE		double
#define TYPE_PREFIX	d
#define TYPE_POSTFIX	d
#define WORK_XX		WORK_DX
#define WORK_XXX	WORK_DXX
#define WORK_XLAYERS	WORK_DLAYERS
#define XEXP		exp
#define XREAL(x)	(x)

#define SOLVE_BVP				solve_bvp
#define SOLVE_BVP_SAVE_TREE_STRING		"solve_bvp"
#define SOLVE_BVP_TL_WITH_AD			solve_bvp_tl_with_ad
#define FORWARD_SAVE_SOLVE_BVP_DATA		forward_save_solve_bvp_data
#define FORWARD_SAVE_SOLVE_BVP_ALLOC		forward_save_solve_bvp_alloc
#define CALC_RADIANCE_LEVELS			calc_radiance_levels
#define CALC_RADIANCE_LEVELS_SAVE_TREE_STRING	"calc_radiance_levels"
#define FORWARD_SAVE_CALC_RADIANCE_LEVELS_DATA	forward_save_calc_radiance_levels_data
#define FORWARD_SAVE_CALC_RADIANCE_LEVELS_ALLOC	forward_save_calc_radiance_levels_alloc
#define CALC_RADIANCE_TAUS			calc_radiance_taus

#include "type_set.h"

#include "xrtm_eig_bvp_x.c"

#include "type_unset.h"

#undef REAL

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

#define SOLVE_BVP				solve_bvp2
#define SOLVE_BVP_TL_WITH_AD			solve_bvp_tl_with_ad2
#define SOLVE_BVP_SAVE_TREE_STRING		"solve_bvp2"
#define FORWARD_SAVE_SOLVE_BVP_DATA		forward_save_solve_bvp_data2
#define FORWARD_SAVE_SOLVE_BVP_ALLOC		forward_save_solve_bvp_alloc2
#define CALC_RADIANCE_LEVELS			calc_radiance_levels2
#define CALC_RADIANCE_LEVELS_SAVE_TREE_STRING	"calc_radiance_levels2"
#define FORWARD_SAVE_CALC_RADIANCE_LEVELS_DATA	forward_save_calc_radiance_levels_data2
#define FORWARD_SAVE_CALC_RADIANCE_LEVELS_ALLOC	forward_save_calc_radiance_levels_alloc2
#define CALC_RADIANCE_TAUS			calc_radiance_taus2

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
void rtm_eig_bvp(int i_four,
                 int n_quad, int n_stokes, int n_derivs, int n_layers,
                 double qf, double *qx_v, double *qw_v, double F_0, double mu_0,
                 int n_ulevels, int *ulevels, double *utaus, int n_umus, double *umus,
                 double *planck, double **planck_l,
                 double *omega, double **omega_l, double *ltau, double **ltau_l,
                 double surface_b, double *surface_b_l,
                 double *btran, double **btran_l,
                 double *as_0, double **as_0_l, double *atran, double **atran_l,
                 double **P_q0_mm, double **P_q0_pm, double **P_u0_mm, double **P_u0_pm,
                 double ***P_uq_pp, double ***P_uq_mp, double ***P_uq_mm, double ***P_uq_pm,
                 double ***r_p, double ***t_p, double ***r_m, double ***t_m,
                 double ***P_q0_mm_l, double ***P_q0_pm_l, double ***P_u0_mm_l, double ***P_u0_pm_l,
                 double ****P_uq_pp_l, double ****P_uq_mp_l, double ****P_uq_mm_l, double ****P_uq_pm_l,
                 double ****r_p_l, double ****t_p_l, double ****r_m_l, double ****t_m_l,
                 double **Rs_qq, double ***Rs_qq_l,
                 double *Rs_u0, double **Rs_u0_l, double **Rs_uq, double ***Rs_uq_l,
                 double *I1_m, double **I1_m_l, double *In_p, double **In_p_l,
                 double F_iso_top, double *F_iso_top_l, double F_iso_bot, double *F_iso_bot_l,
                 double **I_p, double **I_m, double ***I_p_l, double ***I_m_l,
                 int add_single_scattering, int greens, int sfi, int solar, int thermal,
                 int surface, int upwelling, int downwelling, int utau_output, int vector,
                 int eigen_solver_real, int eigen_solver_complex,
                 derivs_data *derivs, save_tree_data save_tree, work_data work) {

     uchar *derivs_layers2;
     uchar *derivs_beam2;
     uchar *derivs_thermal2;

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
     double *atran_l2;

     double *planck0_l2;
     double *planck1_l2;

     double  **p_d_tau;
     double ***p_d_tau_l;

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

     double **Fs_p;
     double **Fs_m;

     double **Fs1_p;
     double **Fs1_m;

     double **Fs_p_l2;
     double **Fs_m_l2;
     double **Fs1_p_l2;
     double **Fs1_m_l2;

     double ***Fs_p_l;
     double ***Fs_m_l;

     double ***Fs1_p_l;
     double ***Fs1_m_l;

     double ***At_p;
     double ***At_m;

     double ***At_p_l2;
     double ***At_m_l2;

     double ****At_p_l;
     double ****At_m_l;

     double **Ft0_p;
     double **Ft0_m;
     double **Ft1_p;
     double **Ft1_m;

     double ***Ft0_p_l;
     double ***Ft0_m_l;
     double ***Ft1_p_l;
     double ***Ft1_m_l;

     double **Ft0_p_l2;
     double **Ft0_m_l2;
     double **Ft1_p_l2;
     double **Ft1_m_l2;

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
double *atran2;
double **atran2_l;
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

     if (solar) {
          Fs_p = get_work2(&work, WORK_DX, WORK_LAYERS_V, NULL);
          Fs_m = get_work2(&work, WORK_DX, WORK_LAYERS_V, NULL);

          if (! greens) {
               Fs1_p = Fs_p;
               Fs1_m = Fs_m;
          }
          else {
               Fs1_p = get_work2(&work, WORK_DX, WORK_LAYERS_V, NULL);
               Fs1_m = get_work2(&work, WORK_DX, WORK_LAYERS_V, NULL);
          }
     }

     if (thermal) {
          Ft0_p = get_work2(&work, WORK_DX, WORK_LAYERS_V, NULL);
          Ft0_m = get_work2(&work, WORK_DX, WORK_LAYERS_V, NULL);
          Ft1_p = get_work2(&work, WORK_DX, WORK_LAYERS_V, NULL);
          Ft1_m = get_work2(&work, WORK_DX, WORK_LAYERS_V, NULL);

          if (utau_output || sfi) {
               At_p = alloc_array3_d(n_layers, 2, n_quad_v_x);
               At_m = alloc_array3_d(n_layers, 2, n_quad_v_x);
          }

          p_d_tau = get_work_d2(&work, n_layers, 2);
     }

     if (solar) {
          if (flags_or2(derivs->beam, n_layers, n_derivs)) {
               Fs_p_l = get_work3(&work, WORK_DX, WORK_BOTH_V, derivs->beam);
               Fs_m_l = get_work3(&work, WORK_DX, WORK_BOTH_V, derivs->beam);

               if (! greens) {
                    Fs1_p_l = Fs_p_l;
                    Fs1_m_l = Fs_m_l;
               }
               else {
                    Fs1_p_l = get_work3(&work, WORK_DX, WORK_BOTH_V, derivs->beam);
                    Fs1_m_l = get_work3(&work, WORK_DX, WORK_BOTH_V, derivs->beam);
               }
          }
     }

     if (thermal) {
          if (flags_or2(derivs->thermal, n_layers, n_derivs)) {
               Ft0_p_l = get_work3(&work, WORK_DX, WORK_BOTH_V, derivs->thermal);
               Ft0_m_l = get_work3(&work, WORK_DX, WORK_BOTH_V, derivs->thermal);
               Ft1_p_l = get_work3(&work, WORK_DX, WORK_BOTH_V, derivs->thermal);
               Ft1_m_l = get_work3(&work, WORK_DX, WORK_BOTH_V, derivs->thermal);

               if (utau_output || sfi) {
                    At_p_l = alloc_array4_d(n_layers, n_derivs, 2, n_quad_v_x);
                    At_m_l = alloc_array4_d(n_layers, n_derivs, 2, n_quad_v_x);
               }

               p_d_tau_l = get_work_d3(&work, n_layers, n_derivs, 2);
          }
     }

     if (! vector) {
          B = get_work_d1(&work, n_comp);
          if (n_derivs > 0)
               B_l = get_work_d2(&work, n_derivs, n_comp);

          nu  = get_work2(&work, WORK_DX, WORK_LAYERS_V, NULL);

          X_p = get_work2(&work, WORK_DXX, WORK_LAYERS_V, NULL);
          X_m = get_work2(&work, WORK_DXX, WORK_LAYERS_V, NULL);

          if (flags_or2(derivs->layers, n_layers, n_derivs)) {
               nu_l  = get_work3(&work, WORK_DX, WORK_BOTH_V, derivs->layers);

               X_p_l = get_work3(&work, WORK_DXX, WORK_BOTH_V, derivs->layers);
               X_m_l = get_work3(&work, WORK_DXX, WORK_BOTH_V, derivs->layers);
          }
     }
     else {
          B_c = get_work_dc1(&work, n_comp);
          if (n_derivs > 0)
               B_l_c = get_work_dc2(&work, n_derivs, n_comp);

          nu_c  = get_work2(&work, WORK_ZX, WORK_LAYERS_V, NULL);

          X_p_c = get_work2(&work, WORK_ZXX, WORK_LAYERS_V, NULL);
          X_m_c = get_work2(&work, WORK_ZXX, WORK_LAYERS_V, NULL);

          if (flags_or2(derivs->layers, n_layers, n_derivs)) {
               nu_l_c  = get_work3(&work, WORK_ZX, WORK_BOTH_V, derivs->layers);

               X_p_l_c = get_work3(&work, WORK_ZXX, WORK_BOTH_V, derivs->layers);
               X_m_l_c = get_work3(&work, WORK_ZXX, WORK_BOTH_V, derivs->layers);
          }
     }
if (! greens) {
     atran2   = atran;
     atran2_l = atran_l;
}
else {
     atran2 = alloc_array1_d(n_layers);
     init_array1_d(atran2, n_layers, 1.);

     if (n_derivs > 0) {
          atran2_l = alloc_array2_d(n_layers, n_derivs);
          init_array2_d(atran2_l, n_layers, n_derivs, 0.);
     }
}
     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     for (i = 0; i < n_layers; ++i) {
          work2 = work;

          if (save_tree.t)
               save_tree_recode_i(&save_tree, i, i == 0);


          if (flags_or2(derivs->layers, n_layers, n_derivs)) {
               tmr_l   = get_work2(&work2, WORK_DXX, WORK_DERIVS_V, derivs->layers[i]);
               tpr_l   = get_work2(&work2, WORK_DXX, WORK_DERIVS_V, derivs->layers[i]);
               gamma_l = get_work2(&work2, WORK_DXX, WORK_DERIVS_V, derivs->layers[i]);
          }

          if (n_derivs > 0) {
               derivs_layers2 = derivs->layers[i];
               if (solar)
                    derivs_beam2 = derivs->beam[i];
               if (thermal)
                    derivs_thermal2 = derivs->thermal[i];
          }

          if (n_derivs > 0 && flags_or(derivs->layers[i], n_derivs)) {
               omega_l2 = omega_l[i];
               ltau_l2  = ltau_l [i];
               r_p_l2   = r_p_l  [i];
               t_p_l2   = t_p_l  [i];
          }

          if (solar) {
               if (n_derivs > 0 && flags_or(derivs->beam[i], n_derivs)) {
                    btran_l2  = btran_l[i];
                    as_0_l2   = as_0_l [i];
                    atran_l2  = atran_l[i];
                    Fs_p_l2   = Fs_p_l [i];
                    Fs_m_l2   = Fs_m_l [i];
                    Fs1_p_l2  = Fs1_p_l[i];
                    Fs1_m_l2  = Fs1_m_l[i];
               }

               if (n_derivs > 0 && flags_or(derivs->layers[i], n_derivs)) {
                    P_q0_mm_l2 = P_q0_mm_l[i];
                    P_q0_pm_l2 = P_q0_pm_l[i];
               }
          }

          if (thermal) {
               if (n_derivs > 0 && flags_or(derivs->thermal[i], n_derivs)) {
                    planck0_l2 = planck_l[i];
                    planck1_l2 = planck_l[i+1];
               }

               if (n_derivs > 0 && (flags_or(derivs->layers[i], n_derivs) || flags_or(derivs->thermal[i], n_derivs))) {
                    Ft0_p_l2 = Ft0_p_l[i];
                    Ft0_m_l2 = Ft0_m_l[i];
                    Ft1_p_l2 = Ft1_p_l[i];
                    Ft1_m_l2 = Ft1_m_l[i];
                    if (utau_output || sfi) {
                         At_p_l2 = At_p_l[i];
                         At_m_l2 = At_m_l[i];
                    }
               }
          }


          /*--------------------------------------------------------------------
           *
           *------------------------------------------------------------------*/
          build_txr(n_quad_x, n_stokes, n_derivs, r_p[i], t_p[i], tpr, tmr, r_p_l2, t_p_l2, tpr_l, tmr_l, derivs_layers2, work2);

          build_gamma(n_quad_v_x, n_derivs, tpr, tmr, gamma, tpr_l, tmr_l, gamma_l, derivs_layers2, work2);


          /*--------------------------------------------------------------------
           *
           *------------------------------------------------------------------*/
          if (! vector) {
               if (n_derivs > 0 && flags_or(derivs->layers[i], n_derivs)) {
                    nu_l2  = nu_l [i];
                    X_p_l2 = X_p_l[i];
                    X_m_l2 = X_m_l[i];
               }

               eig_2n_gen_real(n_quad_v_x, n_derivs, tpr, tmr, gamma, nu[i], X_p[i], X_m[i], tpr_l, tmr_l, gamma_l, nu_l2, X_p_l2, X_m_l2, eigen_solver_real, derivs_layers2, save_tree, work2);

               if (solar) {
                    if (! CLASSICAL_PARTICULAR_SOLUTION_USE_2D) {

                         if (! greens)
                              build_source_vectors_solar_classic_1n(n_quad_x, n_stokes, n_derivs, qx_v, F_0, omega[i], omega_l2, as_0[i], as_0_l2, P_q0_mm[i], P_q0_pm[i], tpr, tmr, gamma, Fs_p[i], Fs_m[i], P_q0_mm_l2, P_q0_pm_l2, tpr_l, tmr_l, gamma_l, Fs_p_l2, Fs_m_l2, derivs_layers2, derivs_beam2, save_tree, work2);
                         else
                              build_source_vectors_solar_green_s_1n(n_quad_x, n_stokes, n_derivs, qx_v, qw_v, F_0, omega[i], omega_l2, ltau[i], ltau_l2, as_0[i], as_0_l2, atran[i], atran_l2, P_q0_mm[i], P_q0_pm[i], nu[i], X_p[i], X_m[i], Fs_p[i], Fs_m[i], Fs1_p[i], Fs1_m[i], P_q0_mm_l2, P_q0_pm_l2, nu_l2, X_p_l2, X_m_l2, Fs_p_l2, Fs_m_l2, Fs1_p_l2, Fs1_m_l2, derivs_layers2, derivs_beam2, save_tree, work2);
                    }
                    else
                         build_source_vectors_solar_classic_2n(n_quad_v_x, n_derivs, qx_v, F_0, omega[i], omega_l2, as_0[i], as_0_l2, P_q0_mm[i], P_q0_pm[i], r_p[i], t_p[i], Fs_p[i], Fs_m[i], P_q0_mm_l2, P_q0_pm_l2, r_p_l2, t_p_l2, Fs_p_l2, Fs_m_l2, derivs_layers2, derivs_beam2, work2);

                    if (save_tree.t) {
                         copy_array1_d(save->Fs_p [i], Fs_p[i], n_quad_v_x);
                         copy_array1_d(save->Fs_m [i], Fs_m[i], n_quad_v_x);
                    }

                    scale_source_vectors_solar(n_quad_v_x, n_derivs, btran[i], btran_l2, Fs_p[i], Fs_m[i], Fs_p[i], Fs_m[i], Fs_p_l2, Fs_m_l2, Fs_p_l2, Fs_m_l2, derivs_layers2, derivs_beam2, work2);

                    if (save_tree.t) {
                         copy_array1_d(save->Fs_p2[i], Fs_p[i], n_quad_v_x);
                         copy_array1_d(save->Fs_m2[i], Fs_m[i], n_quad_v_x);
                    }
               }

               if (thermal) {
                    if (! utau_output && ! sfi)
                         build_source_vectors_thermal2    (n_quad_x, n_stokes, n_derivs, qx_v, planck[i], planck[i+1], planck0_l2, planck1_l2, omega[i], omega_l2, ltau[i], ltau_l2, r_p[i], t_p[i], r_p[i], t_p[i], Ft0_p[i], Ft0_m[i], Ft1_p[i], Ft1_m[i], r_p_l2, t_p_l2, r_p_l2, t_p_l2, Ft0_p_l2, Ft0_m_l2, Ft1_p_l2, Ft1_m_l2, derivs_layers2, derivs_thermal2, work2);
                    else {
                         build_thermal_source_coefficients(n_quad_x, n_stokes, n_derivs, qx_v, planck[i], planck[i+1], planck0_l2, planck1_l2, omega[i], omega_l2, ltau[i], ltau_l2, r_p[i], t_p[i], r_p[i], t_p[i], At_p[i], At_m[i],                       r_p_l2, t_p_l2, r_p_l2, t_p_l2, At_p_l2, At_m_l2,                       derivs_layers2, derivs_thermal2, work);

                         build_source_vectors_thermal     (n_quad_x, n_stokes, n_derivs, ltau[i], ltau_l2, At_p[i], At_m[i], Ft0_p[i], Ft0_m[i], Ft1_p[i], Ft1_m[i], At_p_l2, At_m_l2, Ft0_p_l2, Ft0_m_l2, Ft1_p_l2, Ft1_m_l2, derivs_layers2, derivs_thermal2, work);
                    }
               }
          }


          /*--------------------------------------------------------------------
           *
           *------------------------------------------------------------------*/
          else {
               if (n_derivs > 0 && flags_or(derivs->layers[i], n_derivs)) {
                    r_m_l2 = r_m_l[i];
                    t_m_l2 = t_m_l[i];
               }

               if (n_derivs > 0 && flags_or(derivs->layers[i], n_derivs)) {
                    nu_l_c2  = nu_l_c [i];
                    X_p_l_c2 = X_p_l_c[i];
                    X_m_l_c2 = X_m_l_c[i];
               }

               eig_2n_gen_complex(n_quad_v_x, n_derivs, tpr, tmr, gamma, nu_c[i], X_p_c[i], X_m_c[i], tpr_l, tmr_l, gamma_l, nu_l_c2, X_p_l_c2, X_m_l_c2, eigen_solver_complex, derivs_layers2, save_tree, work2);

               if (solar) {
                    if (! CLASSICAL_PARTICULAR_SOLUTION_USE_2D) {
                         if (! greens)
                              build_source_vectors_solar_classic_1n(n_quad_x, n_stokes, n_derivs, qx_v, F_0, omega[i], omega_l2, as_0[i], as_0_l2, P_q0_mm[i], P_q0_pm[i], tpr, tmr, gamma, Fs_p[i], Fs_m[i], P_q0_mm_l2, P_q0_pm_l2, tpr_l, tmr_l, gamma_l, Fs_p_l2, Fs_m_l2, derivs_layers2, derivs_beam2, save_tree, work2);
                         else
                              build_source_vectors_solar_green_s_1n2(n_quad_x, n_stokes, n_derivs, qx_v, qw_v, F_0, omega[i], omega_l2, ltau[i], ltau_l2, as_0[i], as_0_l2, atran[i], atran_l2, P_q0_mm[i], P_q0_pm[i], nu_c[i], X_p_c[i], X_m_c[i], Fs_p[i], Fs_m[i], Fs1_p[i], Fs1_m[i], P_q0_mm_l2, P_q0_pm_l2, nu_l_c2, X_p_l_c2, X_m_l_c2, Fs_p_l2, Fs_m_l2, Fs1_p_l2, Fs1_m_l2, derivs_layers2, derivs_beam2, save_tree, work2);

                         dm_v_mul_D_A(n_quad_x, n_stokes, Fs_m[i], Fs_m[i]);
                         if (n_derivs > 0 && flags_or(derivs->beam[i], n_derivs)) {
                              for (j = 0; j < n_derivs; ++j) {
                                   if (! derivs->beam[i][j])
                                        continue;

                                   dm_v_mul_D_A(n_quad_x, n_stokes, Fs_m_l[i][j], Fs_m_l[i][j]);
                              }
                         }
                    }
                    else
                         build_source_vectors_solar_classic_2n2(n_quad_v_x, n_derivs, qx_v, F_0, omega[i], omega_l2, as_0[i], as_0_l[i], P_q0_mm[i], P_q0_pm[i], r_p[i], t_p[i], r_m[i], t_m[i], Fs_p[i], Fs_m[i], P_q0_mm_l2, P_q0_pm_l2, r_p_l2, t_p_l2, r_m_l2, t_m_l2, Fs_p_l2, Fs_m_l2, derivs_layers2, derivs_beam2, work2);

                    if (save_tree.t) {
                         copy_array1_d(save->Fs_p [i], Fs_p[i], n_quad_v_x);
                         copy_array1_d(save->Fs_m [i], Fs_m[i], n_quad_v_x);
                    }

                    scale_source_vectors_solar(n_quad_v_x, n_derivs, btran[i], btran_l2, Fs_p[i], Fs_m[i], Fs_p[i], Fs_m[i], Fs_p_l2, Fs_m_l2, Fs_p_l2, Fs_m_l2, derivs_layers2, derivs_beam2, work2);
                    if (save_tree.t) {
                         copy_array1_d(save->Fs_p2[i], Fs_p[i], n_quad_v_x);
                         copy_array1_d(save->Fs_m2[i], Fs_m[i], n_quad_v_x);
                    }
               }

               if (thermal) {
                    if (! utau_output && ! sfi)
                         build_source_vectors_thermal2    (n_quad_x, n_stokes, n_derivs, qx_v, planck[i], planck[i+1], planck0_l2, planck1_l2, omega[i], omega_l2, ltau[i], ltau_l2, r_p[i], t_p[i], r_m[i], t_m[i], Ft0_p[i], Ft0_m[i], Ft1_p[i], Ft1_m[i], r_p_l2, t_p_l2, r_m_l2, t_m_l2, Ft0_p_l2, Ft0_m_l2, Ft1_p_l2, Ft1_m_l2, derivs_layers2, derivs_thermal2, work2);
                    else {
                         build_thermal_source_coefficients(n_quad_x, n_stokes, n_derivs, qx_v, planck[i], planck[i+1], planck0_l2, planck1_l2, omega[i], omega_l2, ltau[i], ltau_l2, r_p[i], t_p[i], r_m[i], t_m[i], At_p[i], At_m[i],                       r_p_l2, t_p_l2, r_m_l2, t_m_l2, At_p_l2, At_m_l2,                       derivs_layers2, derivs_thermal2, work);

                         build_source_vectors_thermal     (n_quad_x, n_stokes, n_derivs, ltau[i], ltau_l2, At_p[i], At_m[i], Ft0_p[i], Ft0_m[i], Ft1_p[i], Ft1_m[i], At_p_l2, At_m_l2, Ft0_p_l2, Ft0_m_l2, Ft1_p_l2, Ft1_m_l2, derivs_layers2, derivs_thermal2, work);
                    }
               }
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
     if (thermal) {
          for (i = 0; i < n_layers; ++i) {
               p_d_tau[i][0] =  planck[i];
               p_d_tau[i][1] = (planck[i+1] - planck[i]) / ltau[i];
               if (n_derivs > 0 && flags_or(derivs->thermal[i], n_derivs)) {
                    for (j = 0; j < n_derivs; ++j) {
                         if (! derivs->thermal[i][j])
                              continue;

                         p_d_tau_l[i][j][0] =  planck_l[i][j];
                         p_d_tau_l[i][j][1] = ((planck_l[i+1][j] - planck_l[i][j]) - p_d_tau[i][1] * ltau_l[i][j]) / ltau[i];
                    }
               }
          }
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (! vector) {
          solve_bvp (n_quad_x, n_stokes, n_derivs, n_layers,
                     ltau, ltau_l,
                     atran2, atran2_l,
                     nu, X_p, X_m,
                     Fs_p, Fs_m,
                     Fs1_p, Fs1_m,
                     Ft0_p, Ft0_m, Ft1_p, Ft1_m,
                     nu_l, X_p_l, X_m_l,
                     Fs_p_l, Fs_m_l,
                     Fs1_p_l, Fs1_m_l,
                     Ft0_p_l, Ft0_m_l, Ft1_p_l, Ft1_m_l,
                     Rs_qq, Rs_qq_l,
                     B, B_l,
                     I1_m, I1_m_l, In_p, In_p_l,
                     surface, solar, thermal,
                     derivs->layers, derivs->beam, derivs->thermal,
                     save_tree, work);

          if (! sfi) {
               if (! utau_output)
                    calc_radiance_levels(n_quad_v_x, n_layers, n_derivs,
                                         n_ulevels, ulevels,
                                         ltau, ltau_l,
                                         atran2, atran2_l,
                                         nu, X_p, X_m,
                                         Fs_p, Fs_m, Fs1_p, Fs1_m,
                                         Ft0_p, Ft0_m, Ft1_p, Ft1_m,
                                         nu_l, X_p_l, X_m_l,
                                         Fs_p_l, Fs_m_l, Fs1_p_l, Fs1_m_l,
                                         Ft0_p_l, Ft0_m_l, Ft1_p_l, Ft1_m_l,
                                         B, B_l,
                                         I_p, I_m, I_p_l, I_m_l,
                                         solar, thermal,
                                         derivs->layers, derivs->beam, derivs->thermal,
                                         save_tree, work);
               else
                    calc_radiance_taus  (n_quad_v_x, n_layers, n_derivs,
                                         n_ulevels, ulevels, utaus,
                                         ltau, ltau_l,
                                         as_0, as_0_l, atran2, atran2_l,
                                         nu, X_p, X_m,
                                         Fs_p, Fs_m,
                                         nu_l, X_p_l, X_m_l,
                                         Fs_p_l, Fs_m_l,
                                         B, B_l,
                                         I_p, I_m, I_p_l, I_m_l,
                                         solar, thermal,
                                         derivs->layers, derivs->beam, derivs->thermal,
                                         save_tree, work);
          }
     }
     else {
          solve_bvp2(n_quad_x, n_stokes, n_derivs, n_layers,
                     ltau, ltau_l,
                     atran2, atran2_l,
                     nu_c, X_p_c, X_m_c,
                     Fs_p, Fs_m,
                     Fs1_p, Fs1_m,
                     Ft0_p, Ft0_m, Ft1_p, Ft1_m,
                     nu_l_c, X_p_l_c, X_m_l_c,
                     Fs_p_l, Fs_m_l,
                     Fs1_p_l, Fs1_m_l,
                     Ft0_p_l, Ft0_m_l, Ft1_p_l, Ft1_m_l,
                     Rs_qq, Rs_qq_l,
                     B_c, B_l_c,
                     I1_m, I1_m_l, In_p, In_p_l,
                     surface, solar, thermal,
                     derivs->layers, derivs->beam, derivs->thermal,
                     save_tree, work);

          if (! sfi) {
               if (! utau_output)
                    calc_radiance_levels2(n_quad_v_x, n_layers, n_derivs,
                                          n_ulevels, ulevels,
                                          ltau, ltau_l,
                                          atran2, atran2_l,
                                          nu_c, X_p_c, X_m_c,
                                          Fs_p, Fs_m, Fs1_p, Fs1_m,
                                          Ft0_p, Ft0_m, Ft1_p, Ft1_m,
                                          nu_l_c, X_p_l_c, X_m_l_c,
                                          Fs_p_l, Fs_m_l, Fs1_p_l, Fs1_m_l,
                                          Ft0_p_l, Ft0_m_l, Ft1_p_l, Ft1_m_l,
                                          B_c, B_l_c,
                                          I_p, I_m, I_p_l, I_m_l,
                                          solar, thermal,
                                          derivs->layers, derivs->beam, derivs->thermal,
                                          save_tree, work);
               else
                    calc_radiance_taus2  (n_quad_v_x, n_layers, n_derivs,
                                          n_ulevels, ulevels, utaus,
                                          ltau, ltau_l,
                                          as_0, as_0_l, atran2, atran2_l,
                                          nu_c, X_p_c, X_m_c,
                                          Fs_p, Fs_m,
                                          nu_l_c, X_p_l_c, X_m_l_c,
                                          Fs_p_l, Fs_m_l,
                                          B_c, B_l_c,
                                          I_p, I_m, I_p_l, I_m_l,
                                          solar, thermal,
                                          derivs->layers, derivs->beam, derivs->thermal,
                                          save_tree, work);
               }
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (save_tree.t) {
          for (i = 0; i < n_layers; ++i) {
               if (thermal) {
                    copy_array1_d(save->Ft0_p[i], Ft0_p[i], n_quad_v_x);
                    copy_array1_d(save->Ft0_m[i], Ft0_m[i], n_quad_v_x);
                    copy_array1_d(save->Ft1_p[i], Ft1_p[i], n_quad_v_x);
                    copy_array1_d(save->Ft1_m[i], Ft1_m[i], n_quad_v_x);
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
               calc_radiance_levels (n_quad_v_x, n_layers, n_derivs,
                                     1, &n_layers,
                                     ltau, ltau_l,
                                     atran2, atran2_l,
                                     nu, X_p, X_m,
                                     Fs_p, Fs_m, Fs1_p, Fs1_m,
                                     Ft0_p, Ft0_m, Ft1_p, Ft1_m,
                                     nu_l, X_p_l, X_m_l,
                                     Fs_p_l, Fs_m_l, Fs1_p_l, Fs1_m_l,
                                     Ft0_p_l, Ft0_m_l, Ft1_p_l, Ft1_m_l,
                                     B, B_l,
                                     &I_p2, &I_m2, &I_p_l2, &I_m_l2,
                                     solar, thermal,
                                     derivs->layers, derivs->beam, derivs->thermal,
                                     save_tree, work);
          else
               calc_radiance_levels2(n_quad_v_x, n_layers, n_derivs,
                                     1, &n_layers,
                                     ltau, ltau_l,
                                     atran2, atran2_l,
                                     nu_c, X_p_c, X_m_c,
                                     Fs_p, Fs_m, Fs1_p, Fs1_m,
                                     Ft0_p, Ft0_m, Ft1_p, Ft1_m,
                                     nu_l_c, X_p_l_c, X_m_l_c,
                                     Fs_p_l, Fs_m_l, Fs1_p_l, Fs1_m_l,
                                     Ft0_p_l, Ft0_m_l, Ft1_p_l, Ft1_m_l,
                                     B_c, B_l_c,
                                     &I_p2, &I_m2, &I_p_l2, &I_m_l2,
                                     solar, thermal,
                                     derivs->layers, derivs->beam, derivs->thermal,
                                     save_tree, work);
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (! vector) {
          if (sfi && n_umus > 0) {
               if (upwelling)
                    sfi_up (i_four,
                            n_quad, n_stokes, n_derivs, n_layers,
                            qf, qx_v, qw_v, F_0, mu_0,
                            n_ulevels, ulevels, utaus, n_umus, umus,
                            omega, omega_l, ltau, ltau_l,
                            surface_b, surface_b_l,
                            p_d_tau, p_d_tau_l,
                            btran, btran_l,
                            as_0, as_0_l, atran, atran_l,
                            P_q0_mm, P_q0_pm,
                            P_u0_pm,
                            P_uq_pp, P_uq_mp, P_uq_mm, P_uq_pm,
                            nu, X_p, X_m,
                            Fs_p, Fs_m,
                            At_p, At_m,
                            Ft0_p, Ft0_m, Ft1_p, Ft1_m,
                            P_q0_mm_l, P_q0_pm_l,
                            P_u0_pm_l,
                            P_uq_pp_l, P_uq_mp_l, P_uq_mm_l, P_uq_pm_l,
                            nu_l, X_p_l, X_m_l,
                            Fs_p_l, Fs_m_l,
                            At_p_l, At_m_l,
                            Rs_u0, Rs_u0_l, Rs_uq, Rs_uq_l,
                            B, B_l,
                            F_iso_bot, F_iso_bot_l,
                            &I_m2, &I_m_l2,
                            I_p, I_p_l, n_quad_v,
                            add_single_scattering, greens, surface, solar, thermal, utau_output,
                            derivs, work);
               else {
                    for (i = 0; i < n_ulevels; ++i) {
                         init_array1_d(I_p[i]+n_quad_v, n_umus_v, 0.);
                         for (j = 0; j < n_derivs; ++j)
                              init_array1_d(I_p_l[i][j]+n_quad_v, n_umus_v, 0.);
                    }
               }

               if (downwelling)
                    sfi_dn (i_four,
                            n_quad, n_stokes, n_derivs, n_layers,
                            qf, qx_v, qw_v, F_0,
                            n_ulevels, ulevels, utaus, n_umus, umus,
                            omega, omega_l, ltau, ltau_l,
                            p_d_tau, p_d_tau_l,
                            btran, btran_l,
                            as_0, as_0_l, atran, atran_l,
                            P_q0_mm, P_q0_pm,
                            P_u0_mm,
                            P_uq_pp, P_uq_mp, P_uq_mm, P_uq_pm,
                            nu, X_p, X_m,
                            Fs_p, Fs_m,
                            At_p, At_m,
                            Ft0_p, Ft0_m, Ft1_p, Ft1_m,
                            P_q0_mm_l, P_q0_pm_l,
                            P_u0_mm_l,
                            P_uq_pp_l, P_uq_mp_l, P_uq_mm_l, P_uq_pm_l,
                            nu_l, X_p_l, X_m_l,
                            Fs_p_l, Fs_m_l,
                            At_p_l, At_m_l,
                            B, B_l,
                            F_iso_top, F_iso_top_l,
                            &I_m2, &I_m_l2,
                            I_m, I_m_l, n_quad_v,
                            add_single_scattering, greens, solar, thermal, utau_output,
                            derivs, work);
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
                    sfi_up2(i_four,
                            n_quad, n_stokes, n_derivs, n_layers,
                            qf, qx_v, qw_v, F_0, mu_0,
                            n_ulevels, ulevels, utaus, n_umus, umus,
                            omega, omega_l, ltau, ltau_l,
                            surface_b, surface_b_l,
                            p_d_tau, p_d_tau_l,
                            btran, btran_l,
                            as_0, as_0_l, atran, atran_l,
                            P_q0_mm, P_q0_pm,
                            P_u0_pm,
                            P_uq_pp, P_uq_mp, P_uq_mm, P_uq_pm,
                            nu_c, X_p_c, X_m_c,
                            Fs_p, Fs_m,
                            At_p, At_m,
                            Ft0_p, Ft0_m, Ft1_p, Ft1_m,
                            P_q0_mm_l, P_q0_pm_l,
                            P_u0_pm_l,
                            P_uq_pp_l, P_uq_mp_l, P_uq_mm_l, P_uq_pm_l,
                            nu_l_c, X_p_l_c, X_m_l_c,
                            Fs_p_l, Fs_m_l,
                            At_p_l, At_m_l,
                            Rs_u0, Rs_u0_l, Rs_uq, Rs_uq_l,
                            B_c, B_l_c,
                            F_iso_bot, F_iso_bot_l,
                            &I_m2, &I_m_l2,
                            I_p, I_p_l, n_quad_v,
                            add_single_scattering, greens, surface, solar, thermal, utau_output,
                            derivs, work);
               else {
                    for (i = 0; i < n_ulevels; ++i) {
                         init_array1_d(I_p[i]+n_quad_v, n_umus_v, 0.);
                         for (j = 0; j < n_derivs; ++j)
                              init_array1_d(I_p_l[i][j]+n_quad_v, n_umus_v, 0.);
                    }
               }

               if (downwelling)
                    sfi_dn2(i_four,
                            n_quad, n_stokes, n_derivs, n_layers,
                            qf, qx_v, qw_v, F_0,
                            n_ulevels, ulevels, utaus, n_umus, umus,
                            omega, omega_l, ltau, ltau_l,
                            p_d_tau, p_d_tau_l,
                            btran, btran_l,
                            as_0, as_0_l, atran, atran_l,
                            P_q0_mm, P_q0_pm,
                            P_u0_mm,
                            P_uq_pp, P_uq_mp, P_uq_mm, P_uq_pm,
                            nu_c, X_p_c, X_m_c,
                            Fs_p, Fs_m,
                            At_p, At_m,
                            Ft0_p, Ft0_m, Ft1_p, Ft1_m,
                            P_q0_mm_l, P_q0_pm_l,
                            P_u0_mm_l,
                            P_uq_pp_l, P_uq_mp_l, P_uq_mm_l, P_uq_pm_l,
                            nu_l_c, X_p_l_c, X_m_l_c,
                            Fs_p_l, Fs_m_l,
                            At_p_l, At_m_l,
                            B_c, B_l_c,
                            F_iso_top, F_iso_top_l,
                            &I_m2, &I_m_l2,
                            I_m, I_m_l, n_quad_v,
                            add_single_scattering, greens, solar, thermal, utau_output,
                            derivs, work);
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
                            derivs_layers, derivs_beam, save_tree, work);
#endif
     if (thermal) {
          if (utau_output || sfi) {
               free_array3_d(At_p);
               free_array3_d(At_m);

               if (flags_or2(derivs->thermal, n_layers, n_derivs)) {
                    free_array4_d(At_p_l);
                    free_array4_d(At_m_l);
               }
          }
     }
     if (greens)
          free_array1_d(atran2);
}
