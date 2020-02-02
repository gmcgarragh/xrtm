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
#include "xrtm_eig_bvp_a.h"
#include "xrtm_eig_util_a.h"
#include "xrtm_matrix.h"
#include "xrtm_save_tree.h"
#include "xrtm_source_a.h"
#include "xrtm_utility.h"
#include "xrtm_utility_a.h"


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

#define SOLVE_BVP_A					solve_bvp_a
#define SOLVE_BVP_TL_WITH_AD				solve_bvp_tl_with_ad
#define SOLVE_BVP_SAVE_TREE_STRING			"solve_bvp"
#define FORWARD_SAVE_SOLVE_BVP_DATA			forward_save_solve_bvp_data
#define FORWARD_SAVE_SOLVE_BVP_ALLOC			forward_save_solve_bvp_alloc
#define FORWARD_SAVE_SOLVE_BVP_FREE			forward_save_solve_bvp_free
#define CALC_RADIANCE_LEVELS_A				calc_radiance_levels_a
#define CALC_RADIANCE_LEVELS_SAVE_TREE_STRING		"calc_radiance_levels"
#define FORWARD_SAVE_CALC_RADIANCE_LEVELS_DATA		forward_save_calc_radiance_levels_data
#define FORWARD_SAVE_CALC_RADIANCE_LEVELS_ALLOC		forward_save_calc_radiance_levels_alloc
#define FORWARD_SAVE_CALC_RADIANCE_LEVELS_FREE		forward_save_calc_radiance_levels_free
#define CALC_RADIANCE_TAUS_A				calc_radiance_taus_a

#include "type_set.h"

#include "xrtm_eig_bvp_x_a.c"

#include "type_unset.h"

#undef TYPE
#undef TYPE_PREFIX
#undef TYPE_POSTFIX
#undef WORK_XX
#undef WORK_XXX
#undef WORK_XLAYERS
#undef XEXP
#undef XREAL

#undef SOLVE_BVP_A
#undef SOLVE_BVP_TL_WITH_AD
#undef SOLVE_BVP_SAVE_TREE_STRING
#undef FORWARD_SAVE_SOLVE_BVP_DATA
#undef FORWARD_SAVE_SOLVE_BVP_ALLOC
#undef FORWARD_SAVE_SOLVE_BVP_FREE
#undef CALC_RADIANCE_LEVELS_A
#undef CALC_RADIANCE_LEVELS_SAVE_TREE_STRING
#undef FORWARD_SAVE_CALC_RADIANCE_LEVELS_DATA
#undef FORWARD_SAVE_CALC_RADIANCE_LEVELS_ALLOC
#undef FORWARD_SAVE_CALC_RADIANCE_LEVELS_FREE
#undef CALC_RADIANCE_TAUS_A


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

#define SOLVE_BVP_A					solve_bvp_a2
#define SOLVE_BVP_TL_WITH_AD				solve_bvp_tl_with_ad2
#define SOLVE_BVP_SAVE_TREE_STRING			"solve_bvp2"
#define FORWARD_SAVE_SOLVE_BVP_DATA			forward_save_solve_bvp_data2
#define FORWARD_SAVE_SOLVE_BVP_ALLOC			forward_save_solve_bvp_alloc2
#define FORWARD_SAVE_SOLVE_BVP_FREE			forward_save_solve_bvp_free2
#define CALC_RADIANCE_LEVELS_A				calc_radiance_levels_a2
#define CALC_RADIANCE_LEVELS_SAVE_TREE_STRING		"calc_radiance_levels2"
#define FORWARD_SAVE_CALC_RADIANCE_LEVELS_DATA		forward_save_calc_radiance_levels_data2
#define FORWARD_SAVE_CALC_RADIANCE_LEVELS_ALLOC		forward_save_calc_radiance_levels_alloc2
#define FORWARD_SAVE_CALC_RADIANCE_LEVELS_FREE		forward_save_calc_radiance_levels_free2
#define CALC_RADIANCE_TAUS_A				calc_radiance_taus_a2

#include "type_set.h"

#include "xrtm_eig_bvp_x_a.c"

#include "type_unset.h"

#undef TYPE
#undef TYPE_PREFIX
#undef TYPE_POSTFIX
#undef WORK_XX
#undef WORK_XXX
#undef WORK_XLAYERS
#undef XEXP
#undef XREAL

#undef SOLVE_BVP_A
#undef SOLVE_BVP_TL_WITH_AD
#undef SOLVE_BVP_SAVE_TREE_STRING
#undef FORWARD_SAVE_SOLVE_BVP_DATA
#undef FORWARD_SAVE_SOLVE_BVP_ALLOC
#undef FORWARD_SAVE_SOLVE_BVP_FREE
#undef CALC_RADIANCE_LEVELS_A
#undef CALC_RADIANCE_LEVELS_SAVE_TREE_STRING
#undef FORWARD_SAVE_CALC_RADIANCE_LEVELS_DATA
#undef FORWARD_SAVE_CALC_RADIANCE_LEVELS_ALLOC
#undef FORWARD_SAVE_CALC_RADIANCE_LEVELS_FREE
#undef CALC_RADIANCE_TAUS_A


/*******************************************************************************
 *
 ******************************************************************************/
static void forward_save_rtm_eig_bvp_free(forward_save_rtm_eig_bvp_data *d, int n_layers, int n_quad_v, int n_comp, int thermal, int vector);

int forward_save_rtm_eig_bvp_alloc(forward_save_rtm_eig_bvp_data *d, int n_layers, int n_quad_v, int n_comp, int thermal, int vector) {

     d->free = (void (*)(void *)) forward_save_rtm_eig_bvp_free;

     d->thermal = thermal;
     d->vector  = vector;

     d->tpr   = alloc_array3_d(n_layers, n_quad_v, n_quad_v);
     d->tmr   = alloc_array3_d(n_layers, n_quad_v, n_quad_v);
     d->gamma = alloc_array3_d(n_layers, n_quad_v, n_quad_v);

     d->Fs_p  = alloc_array2_d(n_layers, n_quad_v);
     d->Fs_m  = alloc_array2_d(n_layers, n_quad_v);
     d->Fs_p2 = alloc_array2_d(n_layers, n_quad_v);
     d->Fs_m2 = alloc_array2_d(n_layers, n_quad_v);

     if (thermal) {
          d->Ft0_p = alloc_array2_d(n_layers, n_quad_v);
          d->Ft0_m = alloc_array2_d(n_layers, n_quad_v);
          d->Ft1_p = alloc_array2_d(n_layers, n_quad_v);
          d->Ft1_m = alloc_array2_d(n_layers, n_quad_v);
     }

     if (! vector) {
          d->nu    = alloc_array2_d(n_layers, n_quad_v);
          d->X_p   = alloc_array3_d(n_layers, n_quad_v, n_quad_v);
          d->X_m   = alloc_array3_d(n_layers, n_quad_v, n_quad_v);
          d->B     = alloc_array1_d(n_comp);
     }
     else {
          d->nu_c  = alloc_array2_dc(n_layers, n_quad_v);
          d->X_p_c = alloc_array3_dc(n_layers, n_quad_v, n_quad_v);
          d->X_m_c = alloc_array3_dc(n_layers, n_quad_v, n_quad_v);
          d->B_c   = alloc_array1_dc(n_comp);
     }

     return 0;
}



static void forward_save_rtm_eig_bvp_free(forward_save_rtm_eig_bvp_data *d, int n_layers, int n_quad_v, int n_comp, int thermal, int vector) {

     free_array3_d(d->tpr);
     free_array3_d(d->tmr);
     free_array3_d(d->gamma);

     free_array2_d(d->Fs_p);
     free_array2_d(d->Fs_m);
     free_array2_d(d->Fs_p2);
     free_array2_d(d->Fs_m2);

     if (d->thermal) {
          free_array2_d(d->Ft0_p);
          free_array2_d(d->Ft0_m);
          free_array2_d(d->Ft1_p);
          free_array2_d(d->Ft1_m);
     }

     if (! d->vector) {
          free_array2_d(d->nu);
          free_array3_d(d->X_p);
          free_array3_d(d->X_m);
          free_array1_d(d->B);
     }
     else {
          free_array2_dc(d->nu_c);
          free_array3_dc(d->X_p_c);
          free_array3_dc(d->X_m_c);
          free_array1_dc(d->B_c);
     }
}



/*******************************************************************************
 *
 ******************************************************************************/
void rtm_eig_bvp_a(int i_four, int n_quad, int n_stokes, int n_layers,
       double qf, double *qx_v, double *qw_v, double F_0, double mu_0,
       int *ulevels, double *utaus, int n_ulevels,
       double *umus, int n_umus, double *planck,
       double *omega, double *omega_a, double *ltau, double *ltau_a,
       double **Rs_qq, double **Rs_qq_a,
       double *Rs_u0, double *Rs_u0_a, double **Rs_uq, double **Rs_uq_a,
       double *btran, double *btran_a,
       double *as_0, double *as_0_a, double *atran, double *atran_a,
       double **P_q0_mm, double **P_q0_pm, double **P_u0_mm, double **P_u0_pm,
       double ***P_uq_pp, double ***P_uq_mp, double ***P_uq_mm, double ***P_uq_pm,
       double ***r_p, double ***t_p, double ***r_m, double ***t_m,
       double **P_q0_mm_a, double **P_q0_pm_a, double **P_u0_mm_a, double **P_u0_pm_a,
       double ***P_uq_pp_a, double ***P_uq_mp_a, double ***P_uq_mm_a, double ***P_uq_pm_a,
       double ***r_p_a, double ***t_p_a, double ***r_m_a, double ***t_m_a,
       double *I1_m, double *I1_m_a, double *In_p, double *In_p_a,
       double **I_p, double **I_m, double **I_p_a, double **I_m_a,
       int sfi, int surface, int thermal, int upwelling,
       int downwelling, int utau_output, int vector,
       uchar *derivs_layers, uchar *derivs_beam, save_tree_data save_tree, work_data work) {

     int i;

     int n_quad_x;
     int n_quad_v;
     int n_quad_v_x;
     int n_umus_v;

     int n_comp;

     double **tmr_a;
     double **tpr_a;
     double **gamma_a;

     double **Fs_p_a;
     double **Fs_m_a;
     double  *Fs_p_a2;
     double  *Fs_m_a2;

     double *B_a;

     double **nu_a;

     double ***X_p_a;
     double ***X_m_a;

     dcomplex *B_a_c;

     dcomplex **nu_a_c;

     dcomplex ***X_p_a_c;
     dcomplex ***X_m_a_c;

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
     save_tree_encode_s(&save_tree, "rtm_eig_bvp");

     save_tree_retrieve_data(&save_tree, forward_save_rtm_eig_bvp_data, &save);


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     tmr_a   = get_work1(&work, WORK_DXX);
     tpr_a   = get_work1(&work, WORK_DXX);
     gamma_a = get_work1(&work, WORK_DXX);

     Fs_p_a   = get_work2(&work, WORK_DX, WORK_LAYERS_V, NULL);
     Fs_m_a   = get_work2(&work, WORK_DX, WORK_LAYERS_V, NULL);
     Fs_p_a2  = get_work1(&work, WORK_DX);
     Fs_m_a2  = get_work1(&work, WORK_DX);

     if (! vector) {
          B_a   = get_work_d1(&work, n_comp);

          nu_a  = get_work2(&work, WORK_DX, WORK_LAYERS_V, NULL);

          X_p_a = get_work2(&work, WORK_DXX, WORK_LAYERS_V, NULL);
          X_m_a = get_work2(&work, WORK_DXX, WORK_LAYERS_V, NULL);
     }
     else {
          B_a_c   = get_work_dc1(&work, n_comp);

          nu_a_c  = get_work2(&work, WORK_ZX, WORK_LAYERS_V, NULL);

          X_p_a_c = get_work2(&work, WORK_ZXX, WORK_LAYERS_V, NULL);
          X_m_a_c = get_work2(&work, WORK_ZXX, WORK_LAYERS_V, NULL);
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     for (i = 0; i < n_layers; ++i) {
          dvec_zero(Fs_p_a [i], n_quad_v_x);
          dvec_zero(Fs_m_a [i], n_quad_v_x);

          if (! vector) {
               dvec_zero(nu_a[i], n_quad_v_x);
               dmat_zero(X_p_a[i], n_quad_v_x, n_quad_v_x);
               dmat_zero(X_m_a[i], n_quad_v_x, n_quad_v_x);
          }
          else {
               zvec_zero(nu_a_c[i], n_quad_v_x);
               zmat_zero(X_p_a_c[i], n_quad_v_x, n_quad_v_x);
               zmat_zero(X_m_a_c[i], n_quad_v_x, n_quad_v_x);
          }
     }

     if (! vector)
          dvec_zero(B_a,   n_comp);
     else
          zvec_zero(B_a_c, n_comp);


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (vector) {
          for (i = 0; i < n_ulevels; ++i)
               dm_v_mul_D_A(n_quad + n_umus, n_stokes, I_m_a[i], I_m_a[i]);
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (! vector) {
          if (! utau_output)
               calc_radiance_levels_a(n_quad_v_x, n_layers, n_ulevels, ulevels, ltau, ltau_a, atran, atran_a, save->nu, save->X_p, save->X_m, save->Fs_p2, save->Fs_m2, save->Ft0_p, save->Ft0_m, save->Ft1_p, save->Ft1_m, nu_a, X_p_a, X_m_a, Fs_p_a, Fs_m_a, NULL, NULL, NULL, NULL, save->B, B_a, I_p, I_m, I_p_a, I_m_a, derivs_layers, derivs_beam, thermal, save_tree, work);
          else
               calc_radiance_taus_a  (n_quad_v_x, n_layers, n_ulevels, ulevels, utaus, ltau, ltau_a, as_0, as_0_a, atran, atran_a, save->nu, save->X_p, save->X_m, save->Fs_p2, save->Fs_m2, nu_a, X_p_a, X_m_a, Fs_p_a, Fs_m_a, save->B, B_a, I_p, I_m, I_p_a, I_m_a, derivs_layers, derivs_beam, save_tree, work);

          solve_bvp_a(n_quad_x, n_stokes, n_layers, ltau, ltau_a, Rs_qq, Rs_qq_a, atran, atran_a, save->nu, save->X_p, save->X_m, save->Fs_p2, save->Fs_m2, save->Ft0_p, save->Ft0_m, save->Ft1_p, save->Ft1_m, nu_a, X_p_a, X_m_a, Fs_p_a, Fs_m_a, NULL, NULL, NULL, NULL, save->B, B_a, I1_m, I1_m_a, In_p, In_p_a, surface, thermal, derivs_layers, derivs_beam, save_tree, work);
     }
     else {
          if (! utau_output)
               calc_radiance_levels_a2(n_quad_v_x, n_layers, n_ulevels, ulevels, ltau, ltau_a, atran, atran_a, save->nu_c, save->X_p_c, save->X_m_c, save->Fs_p2, save->Fs_m2, save->Ft0_p, save->Ft0_m, save->Ft1_p, save->Ft1_m, nu_a_c, X_p_a_c, X_m_a_c, Fs_p_a, Fs_m_a, NULL, NULL, NULL, NULL, save->B_c, B_a_c, I_p, I_m, I_p_a, I_m_a, derivs_layers, derivs_beam, thermal, save_tree, work);
          else
               calc_radiance_taus_a2  (n_quad_v_x, n_layers, n_ulevels, ulevels, utaus, ltau, ltau_a, as_0, as_0_a, atran, atran_a, save->nu_c, save->X_p_c, save->X_m_c, save->Fs_p2, save->Fs_m2, nu_a_c, X_p_a_c, X_m_a_c, Fs_p_a, Fs_m_a, save->B_c, B_a_c, I_p, I_m, I_p_a, I_m_a, derivs_layers, derivs_beam, save_tree, work);

          solve_bvp_a2(n_quad_x, n_stokes, n_layers, ltau, ltau_a, Rs_qq, Rs_qq_a, atran, atran_a, save->nu_c, save->X_p_c, save->X_m_c, save->Fs_p2, save->Fs_m2, save->Ft0_p, save->Ft0_m, save->Ft1_p, save->Ft1_m, nu_a_c, X_p_a_c, X_m_a_c, Fs_p_a, Fs_m_a, NULL, NULL, NULL, NULL, save->B_c, B_a_c, I1_m, I1_m_a, In_p, In_p_a, surface, thermal, derivs_layers, derivs_beam, save_tree, work);
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     for (i = 0; i < n_layers; ++i) {
          work2 = work;

          save_tree_recode_i(&save_tree, i, i == 0);

if (derivs_beam[i]) {
          dvec_zero(Fs_p_a2, n_quad_v_x);
          dvec_zero(Fs_m_a2, n_quad_v_x);
          scale_source_vectors_a(n_quad_v_x, btran[i], &btran_a[i], save->Fs_p[i], save->Fs_m[i], NULL, NULL, Fs_p_a2, Fs_m_a2, Fs_p_a[i], Fs_m_a[i], work2);
          dvec_copy(Fs_p_a[i], Fs_p_a2, n_quad_v_x);
          dvec_copy(Fs_m_a[i], Fs_m_a2, n_quad_v_x);
}
if (derivs_layers[i]) {
          dmat_zero(tmr_a, n_quad_v_x, n_quad_v_x);
          dmat_zero(tpr_a, n_quad_v_x, n_quad_v_x);
          dmat_zero(gamma_a, n_quad_v_x, n_quad_v_x);
}
if (derivs_beam[i]) {
          if (vector)
               dm_v_mul_D_A(n_quad_x, n_stokes, Fs_m_a[i], Fs_m_a[i]);

          build_source_vectors_1n_a(n_quad_x, n_stokes, qx_v, F_0, omega[i], &omega_a[i], as_0[i], &as_0_a[i], P_q0_mm[i], P_q0_pm[i], save->tpr[i], save->tmr[i], save->gamma[i], NULL, NULL, P_q0_mm_a[i], P_q0_pm_a[i], tpr_a, tmr_a, gamma_a, Fs_p_a[i], Fs_m_a[i], save_tree, work2);
}
if (derivs_layers[i]) {
          /*--------------------------------------------------------------------
           *
           *------------------------------------------------------------------*/
          if (! vector)
               eig_2n_gen_real_a(n_quad_v_x, save->tpr[i], save->tmr[i], save->gamma[i], save->nu[i], save->X_p[i], save->X_m[i], tpr_a, tmr_a, gamma_a, nu_a[i], X_p_a[i], X_m_a[i], save_tree, work2);
          else
               eig_2n_gen_complex_a(n_quad_v_x, save->tpr[i], save->tmr[i], save->gamma[i], save->nu_c[i], save->X_p_c[i], save->X_m_c[i], tpr_a, tmr_a, gamma_a, nu_a_c[i], X_p_a_c[i], X_m_a_c[i], save_tree, work2);


          /*--------------------------------------------------------------------
           *
           *------------------------------------------------------------------*/
          build_gamma_a(n_quad_v_x, save->tpr[i], save->tmr[i], save->gamma[i], tpr_a, tmr_a, gamma_a, work2);

          build_txr_a(n_quad_x, n_stokes, r_p[i], t_p[i], save->tpr[i], save->tmr[i], r_p_a[i], t_p_a[i], tpr_a, tmr_a, work2);
}
     }

     save_tree_decode_i(&save_tree);


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     init_array2_d(I_p_a, n_ulevels, n_quad_v_x, 0.);
     init_array2_d(I_m_a, n_ulevels, n_quad_v_x, 0.);
}



/*******************************************************************************
 *
 ******************************************************************************/
#ifdef INCLUDE_DEV_SOURCE
#include "../src2/xrtm_eig_bvp_a2.c"
#endif
