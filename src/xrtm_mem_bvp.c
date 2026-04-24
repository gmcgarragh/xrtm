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
#include "xrtm_eig_util.h"
#include "xrtm_matrix.h"
#include "xrtm_mem_bvp.h"
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
#define XEXP		exp
#define XREAL(x)	(x)

#define SOLVE_BVP		solve_bvp
#define CALC_RADIANCE_LEVELS	calc_radiance_levels
#define CALC_RADIANCE_TAUS	calc_radiance_taus

#include "type_set.h"

#include "xrtm_mem_bvp_x.c"

#include "type_unset.h"

#undef TYPE
#undef TYPE_PREFIX
#undef TYPE_POSTFIX
#undef WORK_XX
#undef WORK_XXX
#undef XEXP
#undef XREAL

#undef SOLVE_BVP
#undef CALC_RADIANCE_LEVELS
#undef CALC_RADIANCE_TAUS



/*******************************************************************************
 *
 ******************************************************************************/
#define TYPE		dcomplex
#define TYPE_PREFIX	z
#define TYPE_POSTFIX	dc
#define WORK_XX		WORK_ZX
#define WORK_XXX	WORK_ZXX
#define XEXP		cexp
#define XREAL(x)	(creal(x))

#define SOLVE_BVP		solve_bvp2
#define CALC_RADIANCE_LEVELS	calc_radiance_levels2
#define CALC_RADIANCE_TAUS	calc_radiance_taus2

#include "type_set.h"

#include "xrtm_mem_bvp_x.c"

#include "type_unset.h"

#undef TYPE
#undef TYPE_PREFIX
#undef TYPE_POSTFIX
#undef WORK_XX
#undef WORK_XXX
#undef XEXP
#undef XREAL

#undef SOLVE_BVP
#undef CALC_RADIANCE_LEVELS
#undef CALC_RADIANCE_TAUS



/*******************************************************************************
 *
 ******************************************************************************/
void rtm_mem_bvp(int i_four,
                 int n_quad, int n_stokes, int n_derivs, int n_layers,
                 double qf, double *qx_v, double *qw_v, double F_0, double mu_0,
                 int n_ulevels, int *ulevels, double *utaus,
                 int n_umus, double *umus,
                 double *omega, double **omega_l, double *ltau, double **ltau_l,
                 double *btau, double **btau_l, double *btran, double **btran_l,
                 double *as_0, double **as_0_l, double *atran, double **atran_l,
                 double **P_q0_mm, double **P_q0_pm, double **P_u0_mm, double **P_u0_pm,
                 double ***P_uq_pp, double ***P_uq_mp, double ***P_uq_mm, double ***P_uq_pm,
                 double ***r_p, double ***t_p, double ***r_m, double ***t_m,
                 double ***P_q0_mm_l, double ***P_q0_pm_l, double ***P_u0_mm_l, double ***P_u0_pm_l,
                 double ****P_uq_pp_l, double ****P_uq_mp_l, double ****P_uq_mm_l, double ****P_uq_pm_l,
                 double ****r_p_l, double ****t_p_l, double ****r_m_l, double ****t_m_l,
                 double **Rs_qq, double ***Rs_qq_l,
                 double *Rs_u0, double **Us_u0, double **Rs_uq, double ***Us_uq,
                 double *I1_m, double **I1_m_l, double *In_p, double **In_p_l,
                 double **I_p, double **I_m, double ***I_p_l, double ***I_m_l,
                 int sfi, int surface, int utau_output, int vector,
                 int eigen_solver_real, int eigen_solver_complex,
                 derivs_data *derivs, save_tree_data save_tree, work_data work) {

     uchar *derivs_layers2;
     uchar *derivs_beam2;

     int i;
     int j;

     int n_quad_x;
     int n_quad_v;
     int n_quad_v_x;
     int n_umus_v;

     int i_quad2;
     int n_quad2;

     int n_comp;
/*
     double *ltau_l2;
*/
     double *omega_l2;
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

     double **Fs_p;
     double **Fs_m;

     double **Fs_p_l2;
     double **Fs_m_l2;

     double ***Fs_p_l;
     double ***Fs_m_l;

     double *B;
     double **B_l;

     double **nu_l2;

     double **nu;
     double ***nu_l;

     double ***X_p;
     double ***X_m;

     double ***X_i_11;
     double ***X_i_12;
     double **sigma_p;
     double **sigma_m;

     double ***X_p_l2;
     double ***X_m_l2;

     double ****X_p_l;
     double ****X_m_l;

     double ****X_i_11_l;
     double ****X_i_12_l;
     double ***sigma_p_l;
     double ***sigma_m_l;

     dcomplex *B_c;
     dcomplex **B_l_c;

     dcomplex **nu_l_c2;

     dcomplex **nu_c;
     dcomplex ***nu_l_c;

     dcomplex ***X_p_c;
     dcomplex ***X_m_c;

     dcomplex ***X_i_11_c;
     dcomplex ***X_i_12_c;
     dcomplex **sigma_p_c;
     dcomplex **sigma_m_c;

     dcomplex ***X_p_l_c2;
     dcomplex ***X_m_l_c2;

     dcomplex ****X_p_l_c;
     dcomplex ****X_m_l_c;

     dcomplex ****X_i_11_l_c;
     dcomplex ****X_i_12_l_c;
     dcomplex ***sigma_p_l_c;
     dcomplex ***sigma_m_l_c;

     work_data work2;


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     calc_quad_and_umus(n_quad, n_umus, n_stokes, &n_quad_x, NULL, &n_quad_v, &n_quad_v_x, NULL, &n_umus_v, sfi);

     n_comp = 2 * n_quad_v_x * (n_layers + 1);


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     Fs_p  = get_work2(&work, WORK_DX, WORK_LAYERS_V, NULL);
     Fs_m  = get_work2(&work, WORK_DX, WORK_LAYERS_V, NULL);

     tmr   = get_work1(&work, WORK_DXX);
     tpr   = get_work1(&work, WORK_DXX);
     gamma = get_work1(&work, WORK_DXX);

     if (flags_or2(derivs->beam, n_layers, n_derivs)) {
          Fs_p_l = get_work3(&work, WORK_DX, WORK_BOTH_V, derivs->beam);
          Fs_m_l = get_work3(&work, WORK_DX, WORK_BOTH_V, derivs->beam);
     }

     if (! vector) {
          B = get_work_d1(&work, n_comp);
          if (n_derivs > 0)
               B_l = get_work_d2(&work, n_derivs, n_comp);

          nu  = get_work2(&work, WORK_DX, WORK_LAYERS_V, NULL);

          X_p = get_work2(&work, WORK_DXX, WORK_LAYERS_V, NULL);
          X_m = get_work2(&work, WORK_DXX, WORK_LAYERS_V, NULL);

          X_i_11  = get_work2(&work, WORK_DXX, WORK_LAYERS_V, NULL);
          X_i_12  = get_work2(&work, WORK_DXX, WORK_LAYERS_V, NULL);

          sigma_p = get_work2(&work, WORK_DX,  WORK_LAYERS_V, NULL);
          sigma_m = get_work2(&work, WORK_DX,  WORK_LAYERS_V, NULL);

          if (flags_or2(derivs->layers, n_layers, n_derivs)) {
               nu_l  = get_work3(&work, WORK_DX, WORK_BOTH_V, derivs->layers);

               X_p_l = get_work3(&work, WORK_DXX, WORK_BOTH_V, derivs->layers);
               X_m_l = get_work3(&work, WORK_DXX, WORK_BOTH_V, derivs->layers);

               X_i_11_l  = get_work3(&work, WORK_DXX, WORK_BOTH_V, derivs->layers);
               X_i_12_l  = get_work3(&work, WORK_DXX, WORK_BOTH_V, derivs->layers);

               sigma_p_l = get_work3(&work, WORK_DX,  WORK_BOTH_V, derivs->beam);
               sigma_m_l = get_work3(&work, WORK_DX,  WORK_BOTH_V, derivs->beam);
          }
     }
     else {
          B_c = get_work_dc1(&work, n_comp);
          if (n_derivs > 0)
               B_l_c = get_work_dc2(&work, n_derivs, n_comp);

          nu_c  = get_work2(&work, WORK_ZX, WORK_LAYERS_V, NULL);

          X_p_c = get_work2(&work, WORK_ZXX, WORK_LAYERS_V, NULL);
          X_m_c = get_work2(&work, WORK_ZXX, WORK_LAYERS_V, NULL);

          X_i_11_c  = get_work2(&work, WORK_ZXX, WORK_LAYERS_V, NULL);
          X_i_12_c  = get_work2(&work, WORK_ZXX, WORK_LAYERS_V, NULL);

          sigma_p_c = get_work2(&work, WORK_ZX,  WORK_LAYERS_V, NULL);
          sigma_m_c = get_work2(&work, WORK_ZX,  WORK_LAYERS_V, NULL);

          if (flags_or2(derivs->layers, n_layers, n_derivs)) {
               nu_l_c  = get_work3(&work, WORK_ZX, WORK_BOTH_V, derivs->layers);

               X_p_l_c = get_work3(&work, WORK_ZXX, WORK_BOTH_V, derivs->layers);
               X_m_l_c = get_work3(&work, WORK_ZXX, WORK_BOTH_V, derivs->layers);

               X_i_11_l_c  = get_work3(&work, WORK_ZXX, WORK_BOTH_V, derivs->layers);
               X_i_12_l_c  = get_work3(&work, WORK_ZXX, WORK_BOTH_V, derivs->layers);

               sigma_p_l_c = get_work3(&work, WORK_ZX,  WORK_BOTH_V, derivs->beam);
               sigma_m_l_c = get_work3(&work, WORK_ZX,  WORK_BOTH_V, derivs->beam);
          }
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     for (i = 0; i < n_layers; ++i) {
          work2 = work;

          if (flags_or2(derivs->layers, n_layers, n_derivs)) {
               tmr_l   = get_work2(&work2, WORK_DXX, WORK_DERIVS_V, derivs->layers[i]);
               tpr_l   = get_work2(&work2, WORK_DXX, WORK_DERIVS_V, derivs->layers[i]);
               gamma_l = get_work2(&work2, WORK_DXX, WORK_DERIVS_V, derivs->layers[i]);
          }

          if (n_derivs > 0) {
               derivs_layers2 = derivs->layers[i];
               derivs_beam2   = derivs->beam[i];
          }

          if (n_derivs > 0 && flags_or(derivs->layers[i], n_derivs)) {
               omega_l2   = omega_l  [i];
/*
               ltau_l2    = ltau_l   [i];
*/
               P_q0_mm_l2 = P_q0_mm_l[i];
               P_q0_pm_l2 = P_q0_pm_l[i];
               r_p_l2     = r_p_l    [i];
               t_p_l2     = t_p_l    [i];
          }

          if (n_derivs > 0 && flags_or(derivs->beam[i], n_derivs)) {
               btran_l2   = btran_l  [i];
               as_0_l2    = as_0_l   [i];
               Fs_p_l2     = Fs_p_l    [i];
               Fs_m_l2     = Fs_m_l    [i];
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
                    nu_l2      = nu_l [i];
                    X_p_l2     = X_p_l[i];
                    X_m_l2     = X_m_l[i];
               }

               eig_2n_gen_real(n_quad_v_x, n_derivs, tpr, tmr, gamma, nu[i], X_p[i], X_m[i], tpr_l, tmr_l, gamma_l, nu_l2, X_p_l2, X_m_l2, eigen_solver_real, derivs_layers2, save_tree, work2);
if (1)
               build_source_vectors_solar_classic_1n(n_quad_x, n_stokes, n_derivs, qx_v, F_0, omega[i], omega_l2, as_0[i], as_0_l2, P_q0_mm[i], P_q0_pm[i], tpr, tmr, gamma, Fs_p[i], Fs_m[i], P_q0_mm_l2, P_q0_pm_l2, tpr_l, tmr_l, gamma_l, Fs_p_l2, Fs_m_l2, derivs_layers2, derivs_beam2, save_tree, work2);
else
               build_source_vectors_solar_classic_2n(n_quad_v_x, n_derivs, qx_v, F_0, omega[i], omega_l2, as_0[i], as_0_l2, P_q0_mm[i], P_q0_pm[i], r_p[i], t_p[i], Fs_p[i], Fs_m[i], P_q0_mm_l2, P_q0_pm_l2, r_p_l2, t_p_l2, Fs_p_l2, Fs_m_l2, derivs_layers2, derivs_beam2, work2);

               scale_source_vectors_solar(n_quad_v_x, n_derivs, btran[i], btran_l2, Fs_p[i], Fs_m[i], Fs_p[i], Fs_m[i], Fs_p_l2, Fs_m_l2, Fs_p_l2, Fs_m_l2, derivs_layers2, derivs_beam2, work);
          }


          /*--------------------------------------------------------------------
           *
           *------------------------------------------------------------------*/
          else {
               if (n_derivs > 0 && flags_or(derivs->layers[i], n_derivs)) {
                    r_m_l2     = r_m_l    [i];
                    t_m_l2     = t_m_l    [i];
               }

               if (n_derivs > 0 && flags_or(derivs->layers[i], n_derivs)) {
                    nu_l_c2    = nu_l_c [i];
                    X_p_l_c2   = X_p_l_c[i];
                    X_m_l_c2   = X_m_l_c[i];
               }

               eig_2n_gen_complex(n_quad_v_x, n_derivs, tpr, tmr, gamma, nu_c[i], X_p_c[i], X_m_c[i], tpr_l, tmr_l, gamma_l, nu_l_c2, X_p_l_c2, X_m_l_c2, eigen_solver_complex, derivs_layers2, save_tree, work2);
if (1) {
               build_source_vectors_solar_classic_1n(n_quad_x, n_stokes, n_derivs, qx_v, F_0, omega[i], omega_l2, as_0[i], as_0_l2, P_q0_mm[i], P_q0_pm[i], tpr, tmr, gamma, Fs_p[i], Fs_m[i], P_q0_mm_l2, P_q0_pm_l2, tpr_l, tmr_l, gamma_l, Fs_p_l2, Fs_m_l2, derivs_layers2, derivs_beam2, save_tree, work2);

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

               scale_source_vectors_solar(n_quad_v_x, n_derivs, btran[i], btran_l2, Fs_p[i], Fs_m[i], Fs_p[i], Fs_m[i], Fs_p_l2, Fs_m_l2, Fs_p_l2, Fs_m_l2, derivs_layers2, derivs_beam2, work);
          }
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (! vector) {
          solve_bvp (n_quad_x, n_stokes, n_derivs, n_layers,
                     qf, qx_v, qw_v, F_0,
                     n_ulevels, ulevels,
                     omega, omega_l, ltau, ltau_l,
                     btau, btau_l, btran, btran_l,
                     atran, atran_l,
                     P_q0_mm, P_q0_pm,
                     nu, X_p, X_m,
                     Fs_p, Fs_m,
                     P_q0_mm_l, P_q0_pm_l,
                     nu_l, X_p_l, X_m_l,
                     Fs_p_l, Fs_m_l,
                     Rs_qq, Rs_qq_l,
                     B, B_l,
                     I1_m, I1_m_l, In_p, In_p_l,
                     surface,
                     derivs->layers, derivs->beam, work,
                     X_i_11, X_i_12, sigma_p, sigma_m,
                     X_i_11_l, X_i_12_l, sigma_p_l, sigma_m_l);

          if (! utau_output)
               calc_radiance_levels (n_quad_v_x, n_layers, n_derivs,
                                     n_ulevels, ulevels, B, B_l,
                                     I_p, I_m, I_p_l, I_m_l, derivs->beam);
          else
               calc_radiance_taus   (n_quad_v_x, n_layers, n_derivs, n_ulevels, ulevels, utaus, ltau, ltau_l, btau, btau_l, nu,   X_p,   X_m,   X_i_11,   X_i_12,   sigma_p,   sigma_m,   nu_l,   X_p_l,   X_m_l,   X_i_11_l,   X_i_12_l,   sigma_p_l,   sigma_m_l,   B,   B_l,   I_p, I_m, I_p_l, I_m_l, derivs->layers, derivs->beam, work);
     }
     else {
          solve_bvp2(n_quad_x, n_stokes, n_derivs, n_layers,
                     qf, qx_v, qw_v, F_0,
                     n_ulevels, ulevels,
                     omega, omega_l, ltau, ltau_l,
                     btau, btau_l, btran, btran_l,
                     atran, atran_l,
                     P_q0_mm, P_q0_pm,
                     nu_c, X_p_c, X_m_c,
                     Fs_p, Fs_m,
                     P_q0_mm_l, P_q0_pm_l,
                     nu_l_c, X_p_l_c, X_m_l_c,
                     Fs_p_l, Fs_m_l,
                     Rs_qq, Rs_qq_l,
                     B_c, B_l_c,
                     I1_m, I1_m_l, In_p, In_p_l,
                     surface,
                     derivs->layers, derivs->beam, work,
                     X_i_11_c, X_i_12_c, sigma_p_c, sigma_m_c,
                     X_i_11_l_c, X_i_12_l_c, sigma_p_l_c, sigma_m_l_c);

          if (! utau_output)
               calc_radiance_levels2(n_quad_v_x, n_layers, n_derivs,
                                     n_ulevels, ulevels, B_c, B_l_c,
                                     I_p, I_m, I_p_l, I_m_l, derivs->beam);
          else
               calc_radiance_taus2  (n_quad_v_x, n_layers, n_derivs, n_ulevels, ulevels, utaus, ltau, ltau_l, btau, btau_l, nu_c, X_p_c, X_m_c, X_i_11_c, X_i_12_c, sigma_p_c, sigma_m_c, nu_l_c, X_p_l_c, X_m_l_c, X_i_11_l_c, X_i_12_l_c, sigma_p_l_c, sigma_m_l_c, B_c, B_l_c, I_p, I_m, I_p_l, I_m_l, derivs->layers, derivs->beam, work);
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
}
