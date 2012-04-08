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
#include "xrtm_eig_rts_a.h"
#include "xrtm_eig_util_a.h"
#include "xrtm_matrix.h"
#include "xrtm_save_tree.h"
#include "xrtm_scatter_a.h"
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
#define XABS		fabs
#define XEXP		exp
#define XREAL(x)	(x)

#define CALC_GLOBAL_R_AND_T_A				calc_global_r_and_t_a
#define CALC_GLOBAL_R_AND_T_TL_WITH_AD			calc_global_r_and_t_tl_with_ad
#define SAVE_TREE_STRING				"calc_global_r_and_t"
#define FORWARD_SAVE_CALC_GLOBAL_R_AND_T_DATA		forward_save_calc_global_r_and_t_data
#define FORWARD_SAVE_CALC_GLOBAL_R_AND_T_ALLOC		forward_save_calc_global_r_and_t_alloc
#define FORWARD_SAVE_CALC_GLOBAL_R_AND_T_FREE		forward_save_calc_global_r_and_t_free

#include "type_set.h"

#include "xrtm_eig_rts_x_a.c"

#include "type_unset.h"

#undef TYPE
#undef TYPE_PREFIX
#undef TYPE_POSTFIX
#undef WORK_XX	
#undef WORK_XXX
#undef XABS
#undef XEXP
#undef XREAL

#undef CALC_GLOBAL_R_AND_T_A
#undef CALC_GLOBAL_R_AND_T_TL_WITH_AD
#undef SAVE_TREE_STRING
#undef FORWARD_SAVE_CALC_GLOBAL_R_AND_T_DATA
#undef FORWARD_SAVE_CALC_GLOBAL_R_AND_T_ALLOC
#undef FORWARD_SAVE_CALC_GLOBAL_R_AND_T_FREE



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

#define CALC_GLOBAL_R_AND_T_A				calc_global_r_and_t_a2
#define CALC_GLOBAL_R_AND_T_TL_WITH_AD			calc_global_r_and_t_tl_with_ad2
#define SAVE_TREE_STRING				"calc_global_r_and_t2"
#define FORWARD_SAVE_CALC_GLOBAL_R_AND_T_DATA		forward_save_calc_global_r_and_t_data2
#define FORWARD_SAVE_CALC_GLOBAL_R_AND_T_ALLOC		forward_save_calc_global_r_and_t_alloc2
#define FORWARD_SAVE_CALC_GLOBAL_R_AND_T_FREE		forward_save_calc_global_r_and_t_free2

#include "type_set.h"

#include "xrtm_eig_rts_x_a.c"

#include "type_unset.h"

#undef TYPE
#undef TYPE_PREFIX
#undef TYPE_POSTFIX
#undef WORK_XX	
#undef WORK_XXX
#undef XABS
#undef XEXP
#undef XREAL

#undef CALC_GLOBAL_R_AND_T_A
#undef CALC_GLOBAL_R_AND_T_TL_WITH_AD
#undef SAVE_TREE_STRING
#undef FORWARD_SAVE_CALC_GLOBAL_R_AND_T_DATA
#undef FORWARD_SAVE_CALC_GLOBAL_R_AND_T_ALLOC
#undef FORWARD_SAVE_CALC_GLOBAL_R_AND_T_FREE



/*******************************************************************************
 *
 ******************************************************************************/
static void forward_save_rtm_eig_rts_free(forward_save_rtm_eig_rts_data *d);

int forward_save_rtm_eig_rts_alloc(forward_save_rtm_eig_rts_data *d, int n_quad_v, int flag) {

     d->free = (void (*)(void *)) forward_save_rtm_eig_rts_free;

     d->flag  = flag;

     d->F_p   = alloc_array1_d(n_quad_v);
     d->F_m   = alloc_array1_d(n_quad_v);

     d->tpr   = alloc_array2_d(n_quad_v, n_quad_v);
     d->tmr   = alloc_array2_d(n_quad_v, n_quad_v);
     d->gamma = alloc_array2_d(n_quad_v, n_quad_v);

     if (! flag) {
          d->nu    = alloc_array1_d(n_quad_v);
          d->X_p   = alloc_array2_d(n_quad_v, n_quad_v);
          d->X_m   = alloc_array2_d(n_quad_v, n_quad_v);
     }
     else {
          d->nu_c  = alloc_array1_dc(n_quad_v);
          d->X_p_c = alloc_array2_dc(n_quad_v, n_quad_v);
          d->X_m_c = alloc_array2_dc(n_quad_v, n_quad_v);
     }

     return 0;
}



static void forward_save_rtm_eig_rts_free(forward_save_rtm_eig_rts_data *d) {

     free_array1_d(d->F_p);
     free_array1_d(d->F_m);

     free_array2_d(d->tpr);
     free_array2_d(d->tmr);
     free_array2_d(d->gamma);

     if (! d->flag) {
          free_array1_d(d->nu);
          free_array2_d(d->X_p);
          free_array2_d(d->X_m);
     }
     else {
          free_array1_dc(d->nu_c);
          free_array2_dc(d->X_p_c);
          free_array2_dc(d->X_m_c);
     }
}



/*******************************************************************************
 *
 ******************************************************************************/
void rtm_eig_rts_a(int n_quad, int n_stokes, double F_0,
                   double *qx_v, double *qw_v, double planck0, double planck1,
                   double omega, double *omega_a, double ltau, double *ltau_a,
                   double as_0, double *as_0_a, double atran, double *atran_a,
                   double *P_0p, double *P_0m,
                   double **r_p, double **t_p, double **r_m, double **t_m,
                   double **R_p, double **T_p, double **R_m, double **T_m,
                   double *S_p, double *S_m,
                   double *P_0p_a, double *P_0m_a,
                   double **r_p_a, double **t_p_a, double **r_m_a, double **t_m_a,
                   double **R_p_a, double **T_p_a, double **R_m_a, double **T_m_a,
                   double *S_p_a, double *S_m_a,
                   int symmetric, int thermal, int vector,
                   int eigen_solver_real, int eigen_solver_complex,
                   uchar derivs_h, uchar derivs_p,
                   save_tree_data save_tree, work_data work) {

     int n_quad_v;

     double **tmr;
     double **tpr;
     double **gamma;

     double **tmr_a;
     double **tpr_a;
     double **gamma_a;

     double *F_p;
     double *F_m;

     double *F_p_a;
     double *F_m_a;

     double *nu;
     double *nu_a;

     double **X_p;
     double **X_m;

     double **X_p_a;
     double **X_m_a;

     dcomplex *nu_c;
     dcomplex *nu_a_c;

     dcomplex **X_p_c;
     dcomplex **X_m_c;

     dcomplex **X_p_a_c;
     dcomplex **X_m_a_c;

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
     save_tree_encode_s(&save_tree, "rtm_eig_rts");

     save_tree_retrieve_data(&save_tree, forward_save_rtm_eig_rts_data, &save);


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     F_p_a   = get_work1(&work, WORK_DX);
     F_m_a   = get_work1(&work, WORK_DX);

     tmr_a   = get_work1(&work, WORK_DXX);
     tpr_a   = get_work1(&work, WORK_DXX);
     gamma_a = get_work1(&work, WORK_DXX);

     if (! vector) {
          nu_a    = get_work1(&work, WORK_DX);
          X_p_a   = get_work1(&work, WORK_DXX);
          X_m_a   = get_work1(&work, WORK_DXX);
     }
     else {
          nu_a_c  = get_work1(&work, WORK_ZX);
          X_p_a_c = get_work1(&work, WORK_ZXX);
          X_m_a_c = get_work1(&work, WORK_ZXX);
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     F_p   = save->F_p;
     F_m   = save->F_m;

     tpr   = save->tpr;
     tmr   = save->tmr;
     gamma = save->gamma;

     if (! vector) {
          nu    = save->nu;
          X_p   = save->X_p;
          X_m   = save->X_m;
     }
     else {
          nu_c  = save->nu_c;
          X_p_c = save->X_p_c;
          X_m_c = save->X_m_c;
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     dvec_zero(F_p_a, n_quad_v);
     dvec_zero(F_m_a, n_quad_v);

     dmat_zero(tmr_a, n_quad_v, n_quad_v);
     dmat_zero(tpr_a, n_quad_v, n_quad_v);
     dmat_zero(gamma_a, n_quad_v, n_quad_v);

     if (! vector) {
          dvec_zero(nu_a, n_quad_v);
          dmat_zero(X_p_a, n_quad_v, n_quad_v);
          dmat_zero(X_m_a, n_quad_v, n_quad_v);
     }
     else {
          zvec_zero(nu_a_c, n_quad_v);
          zmat_zero(X_p_a_c, n_quad_v, n_quad_v);
          zmat_zero(X_m_a_c, n_quad_v, n_quad_v);
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
if (derivs_p) {
     if (! vector)
          build_global_source_a(n_quad_v, atran, atran_a, R_p, T_p, F_p, F_m, S_p, S_m, R_p_a, T_p_a, F_p_a, F_m_a, S_p_a, S_m_a, work);
     else
          build_global_source_a2(n_quad_v, atran, atran_a, R_p, T_p, R_m, T_m, F_p, F_m, S_p, S_m, R_p_a, T_p_a, R_m_a, T_m_a, F_p_a, F_m_a, S_p_a, S_m_a, work);

     build_source_vectors_1n_a(n_quad, n_stokes, qx_v, F_0, omega, omega_a, as_0, as_0_a, P_0p, P_0m, tpr, tmr, gamma, F_p, F_m, P_0p_a, P_0m_a, tpr_a, tmr_a, gamma_a, F_p_a, F_m_a, save_tree, work);
}
if (derivs_h) {
     if (! vector) {
          calc_global_r_and_t_a(n_quad_v, ltau, ltau_a, nu, X_p, X_m, R_p, T_p, nu_a, X_p_a, X_m_a, R_p_a, T_p_a, symmetric, save_tree, work);

          eig_2n_gen_real_a(n_quad_v, tpr, tmr, gamma, nu, X_p, X_m, tpr_a, tmr_a, gamma_a, nu_a, X_p_a, X_m_a, save_tree, work);
     }
     else {
          phase_matrix_symmetry_a3(n_quad, n_stokes, n_quad, n_stokes, T_p_a, R_p_a, T_m_a, R_m_a, 1.);
          dmat_mul_D_A(n_quad, n_stokes, R_p_a, R_p_a);

          calc_global_r_and_t_a2(n_quad_v, ltau, ltau_a, nu_c, X_p_c, X_m_c, NULL, NULL, nu_a_c, X_p_a_c, X_m_a_c, R_p_a, T_p_a, symmetric, save_tree, work);

          eig_2n_gen_complex_a(n_quad_v, tpr, tmr, gamma, nu_c, X_p_c, X_m_c, tpr_a, tmr_a, gamma_a, nu_a_c, X_p_a_c, X_m_a_c, save_tree, work);
     }

     build_gamma_a(n_quad_v, tpr, tmr, gamma, tpr_a, tmr_a, gamma_a, work);

     build_txr_a(n_quad, n_stokes, r_p, t_p, tpr, tmr, r_p_a, t_p_a, tpr_a, tmr_a, work);
}

     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     dmat_zero(R_p_a, n_quad_v, n_quad_v);
     dmat_zero(T_p_a, n_quad_v, n_quad_v);
if (vector) {
     dmat_zero(R_m_a, n_quad_v, n_quad_v);
     dmat_zero(T_m_a, n_quad_v, n_quad_v);
}
     dvec_zero(S_p_a, n_quad_v);
     dvec_zero(S_m_a, n_quad_v);
}



/*******************************************************************************
 *
 ******************************************************************************/
#ifdef INCLUDE_DEV_SOURCE
#include "../src2/xrtm_eig_rts_a2.c"
#endif
