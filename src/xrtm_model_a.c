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

#include <rtutil_math.h>
#include <rtutil_support.h>

#include "xrtm.h"
#include "xrtm_adding.h"
#include "xrtm_adding_a.h"
#include "xrtm_brdf.h"
#include "xrtm_brdf_a.h"
#include "xrtm_eig_bvp_a.h"
#include "xrtm_eig_rts_a.h"
#include "xrtm_external.h"
#include "xrtm_interface.h"
#include "xrtm_model.h"
#include "xrtm_model_a.h"
#include "xrtm_radiance_a.h"
#include "xrtm_save_tree.h"
#include "xrtm_single_a.h"
#include "xrtm_scatter.h"
#include "xrtm_scatter_a.h"
#include "xrtm_source_a.h"
#include "xrtm_stacks.h"
#include "xrtm_support.h"
#include "xrtm_utility.h"
#include "xrtm_utility_a.h"


/*******************************************************************************
 *
 ******************************************************************************/
static void get_delta_m_f_a(xrtm_data *d, int i_layer, double a, double f, double *f_a) {

     if (d->solvers & XRTM_SOLVERS_USE_G) {

     }
     else {
          if (d->n_coef_layer[i_layer] <= d->n_coef2) {

          }
          else {
               d->coef0_a[i_layer][0][d->n_coef2] += *f_a / a;
          }

          *f_a = 0.;
     }
}



/*******************************************************************************
 *
 ******************************************************************************/
static void opt_props_init_a(xrtm_data *d, int i_layer, int n_coef, double *omega_a, double *ltau_a, double ***coef_a, work_data work) {

     int j;
     int k;

     omega_a[i_layer] = 0.;
     ltau_a [i_layer] = 0.;

     for (j = 0; j < d->n_elem; ++j) {
          for (k = 0; k < MIN(n_coef, d->n_coef_layer[i_layer]); ++k) {
               coef_a[i_layer][j][k] = 0.;
          }
     }

     if (d->options & XRTM_OPTION_DELTA_M)
          coef_a[i_layer][0][d->n_coef2] = 0.;
}



static void opt_props_init_up_a(xrtm_data *d, int i_layer, int n_coef, double *omega_a, double *ltau_a, double ***coef_a, work_data work) {

     int i;

     for (i = 0; i <= i_layer; ++i)
          opt_props_init_a(d, i, n_coef, omega_a, ltau_a, coef_a, work);
}


#ifdef INCLUDE_DEV_SOURCE
static void init_opt_props0_a(xrtm_data *d, int i_layer, int n_coef, work_data work) {

     opt_props_init_a(d, i_layer, n_coef, d->omega0_a, d->ltau0_a, d->coef0_a, work);
}



static void init_opt_props0_up_a(xrtm_data *d, int i_layer, int n_coef, work_data work) {

     opt_props_init_up_a(d, i_layer, n_coef, d->omega0_a, d->ltau0_a, d->coef0_a, work);
}
#endif

/*
static
*/
void init_opt_props0_all_a(xrtm_data *d, int n_coef, work_data work) {

     opt_props_init_up_a(d, d->n_layers - 1, n_coef, d->omega0_a, d->ltau0_a, d->coef0_a, work);
}


/*
static void init_opt_props_a(xrtm_data *d, int i_layer, int n_coef, work_data work) {

     opt_props_init_a(d, i_layer, n_coef, d->omega_a, d->ltau_a, d->coef_a, work);
}
*/

#ifdef INCLUDE_DEV_SOURCE
static void init_opt_props_up_a(xrtm_data *d, int i_layer, int n_coef, work_data work) {

     opt_props_init_up_a(d, i_layer, n_coef, d->omega_a, d->ltau_a, d->coef_a, work);
}
#endif
/*
static
*/

void init_opt_props_all_a(xrtm_data *d, int n_coef, work_data work) {

     opt_props_init_up_a(d, d->n_layers - 1, n_coef, d->omega_a, d->ltau_a, d->coef_a, work);
}



static void opt_props_update_a(xrtm_data *d, int i_layer, work_data work) {

     int i;

     int flag;

     int coef_type;

     double a;

     double f;
     double f_a;


     if (! (d->options & XRTM_OPTION_DELTA_M))
          return;


     i = i_layer;

     a = 2. * d->n_coef2 + 1.;

     get_delta_m_f(d, i, 0, 0, a, &f, NULL);
     f_a = 0.;

     coef_type = xrtm_phase_type_to_coef_type(d);

     flag = d->n_stokes != 1 ? 1 : 0;

     delta_m_ltau_a(f, &f_a, d->omega0[i], &d->omega0_a[i], d->ltau0[i], &d->ltau0_a[i], &d->ltau_a[i]);
     delta_m_omega_a(f, &f_a, d->omega0[i], &d->omega0_a[i], &d->omega_a[i]);
     delta_m_coef_a(d->n_coef2, f, &f_a, d->coef0[i], d->coef0_a[i], d->coef_a[i], flag, coef_type);

     get_delta_m_f_a(d, i, a, f, &f_a);
/*
     init_opt_props_a(d, i_layer, d->n_coef2, work);
*/
}



static void update_opt_props_a(xrtm_data *d, int i_layer, work_data work) {

     int i;

     for (i = 0; i <= i_layer; ++i)
          opt_props_update_a(d, i, work);
}



static void update_opt_props_all_a(xrtm_data *d, work_data work) {

     int i;

     for (i = 0; i < d->n_layers; ++i)
          opt_props_update_a(d, i, work);
}



/*******************************************************************************
 *
 ******************************************************************************/
static void beam_params_init_a(xrtm_data *d, int i_level, int i_layer, double *btau_a, double *btran_a, double *as_0_a, double *atran_a, work_data work) {

     int i;

     int i_level2;

     if (i_level == d->n_layers)
          i_level2 = i_level;
     else
          i_level2 = i_level + 1;

     for (i = 0; i <= i_level2; ++i) {
          btau_a [i] = 0.;
          btran_a[i] = 0.;
     }

     for (i = 0; i <= i_layer; ++i) {
          as_0_a [i] = 0.;
          atran_a[i] = 0.;
     }
}



static void init_beam_params_up_a(xrtm_data *d, int i_level, work_data work) {

     int i_layer = i_level < d->n_layers ? i_level : d->n_layers - 1;

     beam_params_init_a(d, i_level, i_layer, d->btau_a, d->btran_a, d->as_0_a, d->atran_a, work);
}

/*
static
*/
void init_beam_params_all_a(xrtm_data *d, work_data work) {

     beam_params_init_a(d, d->n_layers, d->n_layers - 1, d->btau_a, d->btran_a, d->as_0_a, d->atran_a, work);
}



static void beam_params_update_a(xrtm_data *d, int i_level, double *ltau, double *ltau_a, double *btau, double *btau_a, double *btran, double *btran_a, double *as_0, double *as_0_a, double *atran, double *atran_a, work_data work) {

     int i;
     int j;

     int i_level2;

     double t;

     double **chapman;

     if (i_level == d->n_layers)
          i_level2 = i_level;
     else
          i_level2 = i_level + 1;

     if (d->options & XRTM_OPTION_PSA)
          get_chapman(d, &chapman);

     for (i = i_level2; i >= 1; --i) {
          t = -atran_a[i-1] * atran[i-1];

          ltau_a[i-1] += t * as_0[i-1];

          as_0_a[i-1] += t * ltau[i-1];

          t = as_0_a[i-1] / ltau[i-1];

          btau_a[i] += t;

          btau_a[i-1] += -t;

          ltau_a[i-1] += -t * as_0[i-1];

          btau_a[i] += -btran_a[i] * btran[i];

          if (! (d->options & XRTM_OPTION_PSA)) {
               btau_a[i-1] += btau_a[i];
               ltau_a[i-1] += btau_a[i] / d->mu_0;
          }
          else {
               for (j = i - 1; j >= 0; --j)
                    ltau_a[j] += btau_a[i] * d->chapman[i-1][j];
          }
     }

     init_beam_params_up_a(d, i_level, work);
}



static void update_beam_params_a(xrtm_data *d, int i_level, work_data work) {

     beam_params_update_a(d, i_level, d->ltau, d->ltau_a, d->btau, d->btau_a, d->btran, d->btran_a, d->as_0, d->as_0_a, d->atran, d->atran_a, work);
}



static void update_beam_params_all_a(xrtm_data *d, work_data work) {

     beam_params_update_a(d, d->n_layers, d->ltau, d->ltau_a, d->btau, d->btau_a, d->btran, d->btran_a, d->as_0, d->as_0_a, d->atran, d->atran_a, work);
}



/*******************************************************************************
 *
 ******************************************************************************/
static void phase_func_get_a(xrtm_data *d, int n_coef, int i_layer, double mu, double phi, void *polys, double ***coef_a, double *P_a, work_data work) {

     double c;

     if (d->options & XRTM_OPTION_PHASE_SCALAR)
          phase_func_a(n_coef, (double *) polys, coef_a[i_layer][0], *P_a);
     else
     if (d->options & XRTM_OPTION_PHASE_MATRIX_GC) {
          c = scat_angle(d->mu_0, d->phi_0, mu, phi);

          scat_vector_rotate_hovenier_a(d->n_stokes, d->mu_0, mu, c, phi - d->phi_0, P_a, P_a);
          build_scat_vector_gc_a(n_coef, d->n_stokes, (double **) polys, coef_a[i_layer], P_a);
     }
     else
     if (d->options & XRTM_OPTION_PHASE_MATRIX_LC) {

     }
#ifdef DEBUG
     else {
         fprintf(stderr, "ERROR: get_phase_func(): end of if / else if\n");
         exit(1);
     }
#endif
}



/*******************************************************************************
 *
 ******************************************************************************/
static void get_phase_func_trun_a(xrtm_data *d, int i_layer, double mu, double phi, void *polys, double *P_trun_a, work_data work) {

     int n_coef = MIN(d->n_coef_layer[i_layer], d->n_coef2);

     phase_func_get_a(d, n_coef, i_layer, mu, phi, polys, d->coef_a,  P_trun_a, work);
}



/*******************************************************************************
 *
 ******************************************************************************/
static int get_phase_func_trun_all_a(xrtm_data *d, double mu, double phi, void *polys, double **P_trun_a, work_data work) {

     int i;

     for (i = 0; i < d->n_layers; ++i)
          get_phase_func_trun_a(d, i, mu, phi, polys, P_trun_a[i], work);

     return 0;
}



/*******************************************************************************
 *
 ******************************************************************************/
static void get_phase_func_full_a(xrtm_data *d, int i_layer, double mu, double phi, void *polys, double *P_full_a, work_data work) {

     int n_coef =       d->n_coef_layer[i_layer];

     phase_func_get_a(d,  n_coef, i_layer, mu, phi, polys, d->coef0_a, P_full_a, work);
}



/*******************************************************************************
 *
 ******************************************************************************/
static int get_phase_func_full_all_a(xrtm_data *d, double mu, double phi, void *polys, double **P_full_a, work_data work) {

     int i;

     for (i = 0; i < d->n_layers; ++i)
          get_phase_func_full_a(d, i, mu, phi, polys, P_full_a[i], work);

     return 0;
}



/*******************************************************************************
 *
 ******************************************************************************/
static int phase_vecs_get_a(xrtm_data *d, int i_four, int i_layer, int n_mus1, double *mus1, double mu2, double **Y1, double *Y2, double *alpha1, double ***P_pp_a0, double ***P_pm_a0, double *P_pp_a, double *P_pm_a, int flag, int flag2, int flag3, int work_type, work_data *work) {

     int i;

     int n_coef;

     int n_mus_v1;

     work_data work2;


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (d->options & XRTM_OPTION_SAVE_PHASE_MATS) {

     }
     else if (flag) {

     }


     work2 = *work;


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (flag3 && d->misc_input.use_symmetric_form) {
          n_mus_v1 = n_mus1 * d->n_stokes;

          if (d->options & XRTM_OPTION_PHASE_SCALAR ||
              d->options & XRTM_OPTION_PHASE_MATRIX_GC)
               vec_sim_trans2(n_mus_v1, P_pp_a, P_pm_a, alpha1);
          else {
               for (i = 0; i < d->n_four; ++i)
                    vec_sim_trans2(n_mus_v1, P_pp_a0[i][i_layer], P_pm_a0[i][i_layer], alpha1);
          }
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     n_coef = MIN(d->n_coef_layer[i_layer], d->n_coef2);

     if (d->options & XRTM_OPTION_PHASE_SCALAR)
          build_phase_vecs_scalar_a(i_four, n_coef, n_mus1, Y1, Y2, d->n_coef2, d->coef_a[i_layer][0], P_pm_a, P_pp_a);
     else
     if (d->options & XRTM_OPTION_PHASE_MATRIX_GC)
          build_phase_vecs_vector_gc_a(i_four, n_coef, n_mus1, d->n_stokes, mus1, mu2, d->coef_a[i_layer], P_pm_a, P_pp_a, work2);
/*
#ifdef INCLUDE_DEV_SOURCE
     else
     if (d->options & XRTM_OPTION_PHASE_MATRIX_LC)
          build_phase_vecs_vector_lc_a(d->n_four2, n_coef, n_mus1, d->n_stokes, mus1, mu2, d->coef_a[i_layer], P_pm_a, P_pp_a, work2, d->n_layers * d->n_quad_v_x);
#endif
*/
#ifdef DEBUG
     else {
         fprintf(stderr, "ERROR: phase_vecs_get_a(): end of if / else if\n");
         exit(1);
     }
#endif

     return 0;
}



/*******************************************************************************
 *
 ******************************************************************************/
static int phase_vecs_get_all_a(xrtm_data *d, int i_four, int n_mus1, double *mus1, double mu2, double **Y1, double *Y2, double *alpha1, double **P_pp_a, double **P_pm_a, int flag, int flag2, int flag3, enum work_type work_type, work_data *work) {

     int i;

     if (d->options & XRTM_OPTION_SAVE_PHASE_MATS) {

     }
     else {
          if (flag) {

          }

          for (i = 0; i < d->n_layers; ++i) {
               if (phase_vecs_get_a(d, i_four, i, n_mus1, mus1, mu2, Y1, Y2, alpha1, NULL, NULL, P_pp_a[i], P_pm_a[i], 0, flag2, flag3, work_type, work)) {
                    fprintf(stderr, "ERROR: get_phase_vecs_a(), i_four = %d, i_layer = %d\n", i_four, i);
                    return -1;
               }
          }
     }

     return 0;
}



/*******************************************************************************
 *
 ******************************************************************************/
static int get_phase_vecs_q0_a(xrtm_data *d, int i_four, int i_layer,
                               double *P_q0_mm_a, double *P_q0_pm_a,
                               int flag, work_data *work) {

     double *Y_0;

     double **Y_p;

     if (d->options & XRTM_OPTION_PHASE_SCALAR) {
          get_Y_p(d, i_four, &Y_p);
          get_Y_0(d, i_four, &Y_0);
     }

     if (phase_vecs_get_a(d, i_four, i_layer, d->n_quad_x, d->qx, -d->mu_0, Y_p, Y_0, d->alpha1, NULL, NULL, P_q0_mm_a, P_q0_pm_a, flag, 1, 1, WORK_DX, work)) {
          fprintf(stderr, "ERROR: phase_vecs_get_a(), i_four = %d, i_layer = %d\n", i_four, i_layer);
          return -1;
     }

     return 0;
}



static int get_phase_vecs_q0_all_a(xrtm_data *d, int i_four,
                                   double **P_q0_mm_a, double **P_q0_pm_a,
                                   int flag, work_data *work) {

     double *Y_0;

     double **Y_p;

     if (d->options & XRTM_OPTION_PHASE_SCALAR) {
          get_Y_p(d, i_four, &Y_p);
          get_Y_0(d, i_four, &Y_0);
     }

     if (phase_vecs_get_all_a(d, i_four, d->n_quad_x, d->qx, -d->mu_0, Y_p, Y_0, d->alpha1, P_q0_mm_a, P_q0_pm_a, flag, 1, 1, WORK_DX, work)) {
          fprintf(stderr, "ERROR: phase_vecs_get_all_a(), i_four = %d\n", i_four);
          return -1;
     }

     return 0;
}



/*******************************************************************************
 *
 ******************************************************************************/
/*
static int get_phase_vecs_u0_a(xrtm_data *d, int i_four, int i_layer,
                               double *P_u0_mm_a, double *P_u0_pm_a,
                               int flag, work_data *work) {

     double *Y_0;

     double **Y_u;

     if (d->options & XRTM_OPTION_PHASE_SCALAR) {
          get_Y_u(d, i_four, &Y_u);
          get_Y_0(d, i_four, &Y_0);
     }

     if (phase_vecs_get_a(d, i_four, i_layer, d->n_umus, d->umus, -d->mu_0, Y_u, Y_0, d->alpha1 + d->n_quad_v, NULL, NULL, P_u0_mm_a, P_u0_pm_a, flag, 0, 1, WORK_DU, work)) {
          fprintf(stderr, "ERROR: phase_vecs_get_a(), i_four = %d, i_layer = %d\n", i_four, i_layer);
          return -1;
     }

     return 0;
}
*/


static int get_phase_vecs_u0_all_a(xrtm_data *d, int i_four,
                                   double **P_u0_mm_a, double **P_u0_pm_a,
                                   int flag, work_data *work) {

     double *Y_0;

     double **Y_u;

     if (d->options & XRTM_OPTION_PHASE_SCALAR) {
          get_Y_u(d, i_four, &Y_u);
          get_Y_0(d, i_four, &Y_0);
     }

     if (phase_vecs_get_all_a(d, i_four, d->n_umus, d->umus, -d->mu_0, Y_u, Y_0, d->alpha1 + d->n_quad_v, P_u0_mm_a, P_u0_pm_a, flag, 0, 1, WORK_DU, work)) {
          fprintf(stderr, "ERROR: phase_vecs_get_all_a(), i_four = %d\n", i_four);
          return -1;
     }

     return 0;
}



/*******************************************************************************
 *
 ******************************************************************************/
static int phase_mats_get_a(xrtm_data *d, int i_four, int i_layer, int n_mus1, int n_mus2, double *mus1, double *mus2, double **Y1, double **Y2, double *alpha1, double *alpha2, double ****P_pp_a0, double ****P_mp_a0, double ****P_mm_a0, double ****P_pm_a0, double **P_pp_a, double **P_mp_a, double **P_mm_a, double **P_pm_a, int flag, int flag2, int flag3, int work_type, work_data *work) {

     int i;

     int n_coef;

     int n_mus_v1;
     int n_mus_v2;

     work_data work2;


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (d->options & XRTM_OPTION_SAVE_PHASE_MATS) {

     }
     else if (flag) {

     }


     work2 = *work;


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
    if (flag3 && d->misc_input.use_symmetric_form) {
          n_mus_v1 = n_mus1 * d->n_stokes;
          n_mus_v2 = n_mus2 * d->n_stokes;

          if (d->options & XRTM_OPTION_PHASE_SCALAR ||
              d->options & XRTM_OPTION_PHASE_MATRIX_GC) {
               mat_sim_trans2(n_mus_v1, n_mus_v2, P_pp_a, P_mp_a, alpha1, alpha2);

               if (d->options & XRTM_OPTION_VECTOR)
                    mat_sim_trans2(n_mus_v1, n_mus_v2, P_mm_a, P_pm_a, alpha1, alpha2);
          }
          else {
               for (i = 0; i < d->n_four; ++i) {
                    mat_sim_trans2(n_mus_v1, n_mus_v2, P_pp_a0[i][i_layer], P_mp_a0[i][i_layer], alpha1, alpha2);

                    if (d->options & XRTM_OPTION_VECTOR)
                         mat_sim_trans2(n_mus_v1, n_mus_v2, P_mm_a0[i][i_layer], P_pm_a0[i][i_layer], alpha1, alpha2);
               }
          }
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     n_coef = MIN(d->n_coef_layer[i_layer], d->n_coef2);

     if (d->options & XRTM_OPTION_PHASE_SCALAR)
          build_phase_mats_scalar_a(i_four, n_coef, n_mus1, n_mus2, Y1, Y2, d->n_coef2, d->coef_a[i_layer][0], P_pp_a, P_mp_a);
     else
     if (d->options & XRTM_OPTION_PHASE_MATRIX_GC)
          build_phase_mats_vector_gc_a(i_four, n_coef, n_mus1, n_mus2, d->n_stokes, mus1, mus2, d->coef_a[i_layer], P_pp_a, P_mp_a, P_mm_a, P_pm_a, d->options & XRTM_OPTION_VECTOR, work2);
/*
#ifdef INCLUDE_DEV_SOURCE
     else
     if (d->options & XRTM_OPTION_PHASE_MATRIX_LC)
          build_phase_mats_vector_lc_a(d->n_four2, n_coef, n_mus1, n_mus2, d->n_stokes, mus1, mus2, d->coef_a[i_layer], P_pp_a, P_mp_a, P_mm_a, P_pm_a, d->options & XRTM_OPTION_VECTOR, work2, d->n_layers * d->n_quad_v_x * d->n_quad_v_x);
#endif
*/
#ifdef DEBUG
     else {
         fprintf(stderr, "ERROR: phase_mats_get(): end of if / else if\n");
         exit(1);
     }
#endif

     return 0;
}



/*******************************************************************************
 *
 ******************************************************************************/
static int phase_mats_get_all_a(xrtm_data *d, int i_four, int n_mus1, int n_mus2, double *mus1, double *mus2, double **Y1, double **Y2, double *alpha1, double *alpha2, double ***P_pp0_a, double ***P_mp_a0, double ***P_mm_a0, double ***P_pm0_a, double ***P_pp_a, double ***P_mp_a, double ***P_mm_a, double ***P_pm_a, int flag, int flag2, int flag3, enum work_type work_type, work_data *work) {

     int i;

     if (d->options & XRTM_OPTION_SAVE_PHASE_MATS) {

     }
     else {
          if (flag) {

          }

          for (i = 0; i < d->n_layers; ++i) {
               if (phase_mats_get_a(d, i_four, i, n_mus1, n_mus2, mus1, mus2, Y1, Y2, alpha1, alpha2, NULL, NULL, NULL, NULL, P_pp_a[i], P_mp_a[i], P_mm_a[i], P_pm_a[i], 0, flag2, flag3, work_type, work)) {
                    fprintf(stderr, "ERROR: phase_mats_get_a(), i_four = %d, i_layer = %d\n", i_four, i);
                    return -1;
               }
          }
     }

     return 0;
}



/*******************************************************************************
 *
 ******************************************************************************/
static int get_phase_mats_qq_a(xrtm_data *d, int i_four, int i_layer,
                               double **P_qq_pp_a, double **P_qq_mp_a,
                               double **P_qq_mm_a, double **P_qq_pm_a,
                               int flag, work_data *work) {

     double **Y_p;

     if (d->options & XRTM_OPTION_PHASE_SCALAR)
          get_Y_p(d, i_four, &Y_p);

     if (phase_mats_get_a(d, i_four, i_layer, d->n_quad_x, d->n_quad_x, d->qx, d->qx, Y_p, Y_p, d->alpha1, d->alpha2, NULL, NULL, NULL, NULL, P_qq_pp_a, P_qq_mp_a, P_qq_mm_a, P_qq_pm_a, flag, 1, 1, WORK_DXX, work)) {
          fprintf(stderr, "ERROR: phase_mats_get_a(), i_four = %d, i_layer = %d,\n", i_four, i_layer);
          return -1;
     }

     return 0;
}



static int get_phase_mats_qq_all_a(xrtm_data *d, int i_four,
                                   double ***P_qq_pp_a, double ***P_qq_mp_a,
                                   double ***P_qq_mm_a, double ***P_qq_pm_a,
                                   int flag, work_data *work) {

     double **Y_p;

     if (d->options & XRTM_OPTION_PHASE_SCALAR)
          get_Y_p(d, i_four, &Y_p);

     if (phase_mats_get_all_a(d, i_four, d->n_quad_x, d->n_quad_x, d->qx, d->qx, Y_p, Y_p, d->alpha1, d->alpha2, NULL, NULL, NULL, NULL, P_qq_pp_a, P_qq_mp_a, P_qq_mm_a, P_qq_pm_a, flag, 1, 1, WORK_DXX, work)) {
          fprintf(stderr, "ERROR: phase_mats_get_all(), i_four = %d\n", i_four);
          return -1;
     }

     return 0;
}



/*******************************************************************************
 *
 ******************************************************************************/
/*
static int get_phase_mats_uq_a(xrtm_data *d, int i_four, int i_layer,
                               double **P_uq_pp_a, double **P_uq_mp_a,
                               double **P_uq_mm_a, double **P_uq_pm_a,
                               int flag, work_data *work) {

     double **Y_p;
     double **Y_u;

     if (d->options & XRTM_OPTION_PHASE_SCALAR) {
          get_Y_p(d, i_four, &Y_p);
          get_Y_u(d, i_four, &Y_u);
     }

     if (phase_mats_get_a(d, i_four, i_layer, d->n_umus, d->n_quad_x, d->umus, d->qx, Y_u, Y_p, d->alpha1 + d->n_quad_v, d->alpha2, NULL, NULL, NULL, NULL, P_uq_pp_a, P_uq_mp_a, P_uq_mm_a, P_uq_pm_a, flag, 0, 1, WORK_DUX, work)) {
          fprintf(stderr, "ERROR: phase_mats_get_a(), i_four = %d, i_layer = %d\n", i_four, i_layer);
          return -1;
     }

     return 0;
}
*/


static int get_phase_mats_uq_all_a(xrtm_data *d, int i_four,
                                   double ***P_uq_pp_a, double ***P_uq_mp_a,
                                   double ***P_uq_mm_a, double ***P_uq_pm_a,
                                   int flag, work_data *work) {

     double **Y_p;
     double **Y_u;

     if (d->options & XRTM_OPTION_PHASE_SCALAR) {
          get_Y_p(d, i_four, &Y_p);
          get_Y_u(d, i_four, &Y_u);
     }

     if (phase_mats_get_all_a(d, i_four, d->n_umus, d->n_quad_x, d->umus, d->qx, Y_u, Y_p, d->alpha1 + d->n_quad_v, d->alpha2, NULL, NULL, NULL, NULL, P_uq_pp_a, P_uq_mp_a, P_uq_mm_a, P_uq_pm_a, flag, 0, 1, WORK_DUX, work)) {
          fprintf(stderr, "ERROR: phase_mats_get_all_a(), i_four = %d\n", i_four);
          return -1;
     }

     return 0;
}



/*******************************************************************************
 *
 ******************************************************************************/
static void forward_save_get_local_r_t_u_w_free(forward_save_get_local_r_t_u_w_data *d);

int forward_save_get_local_r_t_u_w_alloc(forward_save_get_local_r_t_u_w_data *d, int options, int n_quad_v_x) {

     d->free = (void (*)(void *)) forward_save_get_local_r_t_u_w_free;

     d->options = options;

     d->P_qq_pp = alloc_array2_d(n_quad_v_x, n_quad_v_x);
     d->P_qq_mp = alloc_array2_d(n_quad_v_x, n_quad_v_x);

     if (options & XRTM_OPTION_VECTOR) {
          d->P_qq_mm = alloc_array2_d(n_quad_v_x, n_quad_v_x);
          d->P_qq_pm = alloc_array2_d(n_quad_v_x, n_quad_v_x);
     }

     return 0;
}



static void forward_save_get_local_r_t_u_w_free(forward_save_get_local_r_t_u_w_data *d) {

     free_array2_d(d->P_qq_pp);
     free_array2_d(d->P_qq_mp);

     if (d->options & XRTM_OPTION_VECTOR) {
          free_array2_d(d->P_qq_mm);
          free_array2_d(d->P_qq_pm);
     }
}



/*******************************************************************************
 *
 ******************************************************************************/
static int get_local_r_t_u_w_a
     (xrtm_data *d, int i_four, int i_layer,
      double **r_p, double **t_p, double **r_m, double **t_m,
      double **r_p_a, double **t_p_a, double **r_m_a, double **t_m_a,
      int flag, save_tree_data save_tree, work_data *work) {

     double **P_qq_pp_a;
     double **P_qq_mp_a;
     double **P_qq_mm_a;
     double **P_qq_pm_a;

     forward_save_get_local_r_t_u_w_data *save;

     work_data work2;


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     save_tree_encode_s(&save_tree, "get_local_r_t_u_w");

     save_tree_retrieve_data(&save_tree, forward_save_get_local_r_t_u_w_data, &save);


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (d->options & XRTM_OPTION_SAVE_LOCAL_R_T) {

     }
     else if (flag) {

     }


     work2 = *work;


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     P_qq_pp_a = get_work1(&work2, WORK_DXX);
     P_qq_mp_a = get_work1(&work2, WORK_DXX);

     if (d->options & XRTM_OPTION_VECTOR) {
          P_qq_mm_a = get_work1(&work2, WORK_DXX);
          P_qq_pm_a = get_work1(&work2, WORK_DXX);
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     dmat_zero(P_qq_pp_a, d->n_quad_v_x, d->n_quad_v_x);
     dmat_zero(P_qq_mp_a, d->n_quad_v_x, d->n_quad_v_x);

     if (d->options & XRTM_OPTION_VECTOR) {
          dmat_zero(P_qq_mm_a, d->n_quad_v_x, d->n_quad_v_x);
          dmat_zero(P_qq_pm_a, d->n_quad_v_x, d->n_quad_v_x);
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     build_local_r_and_t_a(i_four, d->n_quad_v_x, d->qx_v, d->qw_v, d->omega[i_layer], &d->omega_a[i_layer], save->P_qq_pp, save->P_qq_mp, r_p, t_p, P_qq_pp_a, P_qq_mp_a, r_p_a, t_p_a, work2);

     if (d->options & XRTM_OPTION_VECTOR)
          build_local_r_and_t_a(i_four, d->n_quad_v_x, d->qx_v, d->qw_v, d->omega[i_layer], &d->omega_a[i_layer], save->P_qq_mm, save->P_qq_pm, r_m, t_m, P_qq_mm_a, P_qq_pm_a, r_m_a, t_m_a, work2);


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (get_phase_mats_qq_a(d, i_four, i_layer, P_qq_pp_a, P_qq_mp_a, P_qq_mm_a, P_qq_pm_a, 1, &work2)) {
          fprintf(stderr, "ERROR: get_phase_mats_qq()\n");
          return -1;
     }


     return 0;
}



/*******************************************************************************
 *
 ******************************************************************************/
static int get_local_r_t_u_w_all_a
     (xrtm_data *d, int i_four,
      double ***r_p, double ***t_p, double ***r_m, double ***t_m,
      double ***r_p_a, double ***t_p_a, double ***r_m_a, double ***t_m_a,
      int flag, save_tree_data save_tree, work_data *work) {

     int i;

     double **r_m_a2;
     double **t_m_a2;


     save_tree_encode_s(&save_tree, "get_local_r_t_u_w_all");


     if (d->options & XRTM_OPTION_SAVE_LOCAL_R_T) {

     }
     else {
          if (flag) {

          }

          for (i = 0; i < d->n_layers; ++i) {
               save_tree_recode_i(&save_tree, i, i == 0);

               if (d->derivs.layers_union[i]) {
                    if (d->options & XRTM_OPTION_VECTOR) {
                         r_m_a2 = r_m_a[i];
                         t_m_a2 = t_m_a[i];
                    }

                    if (get_local_r_t_u_w_a(d, i_four, i, NULL, NULL, NULL, NULL, r_p_a[i], t_p_a[i], r_m_a2, t_m_a2, 0, save_tree, work)) {
                         fprintf(stderr, "ERROR: get_local_r_t_u_w_a()\n");
                         return -1;
                    }
               }
          }

          save_tree_decode_i(&save_tree);
     }

     return 0;
}



/*******************************************************************************
 *
 ******************************************************************************/
static void call_solver_a
     (xrtm_data *d, int solver, int i_four,
      double planck0, double planck1,
      double omega, double *omega_a, double ltau, double *ltau_a,
      double as_0, double *as_0_a, double atran, double *atran_a,
      double *P_q0_mm, double *P_q0_pm,
      double  **r_p, double  **t_p, double  **r_m, double  **t_m,
      double  **R_p, double  **T_p, double  **R_m, double  **T_m,
      double  *S_p,  double  *S_m,
      double *P_q0_mm_a, double *P_q0_pm_a,
      double **r_p_a, double **t_p_a, double **r_m_a, double **t_m_a,
      double **R_p_a, double **T_p_a, double **R_m_a, double **T_m_a,
      double *S_p_a,  double *S_m_a, uchar derivs_layers, uchar derivs_beam,
      save_tree_data save_tree, work_data work) {

     int thermal;
     int vector;

     int symmetric;
     int eigen_solver_real;
     int eigen_solver_complex;
/*
     int check_condition;
*/
     double **R_x_a;
     double **T_x_a;

     thermal = d->options & XRTM_OPTION_SOURCE_THERMAL && i_four == 0;
     vector  = d->options & XRTM_OPTION_VECTOR;

     symmetric            = d->misc_input.use_symmetric_form;
     eigen_solver_real    = d->misc_input.eigen_solver_gen_real;
     eigen_solver_complex = d->misc_input.eigen_solver_gen_complex;
/*
     check_condition      = d->misc_input.use_pade_check_condition;
*/
     if (! (d->options & XRTM_OPTION_VECTOR)) {
          dmat_add(R_p_a, R_m_a, R_p_a, d->n_quad_v_x, d->n_quad_v_x);
          dmat_add(T_p_a, T_m_a, T_p_a, d->n_quad_v_x, d->n_quad_v_x);

          R_x_a = NULL;
          T_x_a = NULL;
     }
     else {
          R_x_a = R_m_a;
          T_x_a = T_m_a;
     }
/*
     if (solver & XRTM_SOLVER_DOUB_ADD) {
          rtm_doub_rts_a (d->n_quad_x, d->n_stokes, d->F_0, d->qx_v, d->qw_v, planck0, planck1, omega, omega_a, ltau, ltau_a, as_0, as_0_a, P_q0_mm, P_q0_pm, r_p, t_p, r_m, t_m, R_p, T_p, R_m, T_m, S_p, S_m, P_q0_mm_a, P_q0_pm_a, r_p_a, t_p_a, r_m_a, t_m_a, R_p_a, T_p_a, R_x_a, T_x_a, S_p_a, S_m_a, d->doub_d_tau, symmetric, thermal, vector, work);
     }
     else
*/
     if (solver & XRTM_SOLVER_EIG_ADD) {
          rtm_eig_rts_a  (d->n_quad_x, d->n_stokes, d->F_0, d->qx_v, d->qw_v, planck0, planck1, omega, omega_a, ltau, ltau_a, as_0, as_0_a, atran, atran_a, P_q0_mm, P_q0_pm, r_p, t_p, r_m, t_m, R_p, T_p, R_m, T_m, S_p, S_m, P_q0_mm_a, P_q0_pm_a, r_p_a, t_p_a, r_m_a, t_m_a, R_p_a, T_p_a, R_x_a, T_x_a, S_p_a, S_m_a, symmetric, thermal, vector, eigen_solver_real, eigen_solver_complex, derivs_layers, derivs_beam, save_tree, work);
     }
/*
     else
     if (solver & XRTM_SOLVER_PADE_ADD) {
          rtm_pade_rts_a2(d->n_quad_x, d->n_stokes, d->F_0, d->qx_v, d->qw_v, d->umus, d->n_umus, planck0, planck1, omega, omega_a, ltau, ltau_a, as_0, as_0_a, atran, atran_a, P_q0_mm, P_q0_pm, r_p, t_p, r_m, t_m, R_p, T_p, R_m, T_m, S_p, S_m, P_q0_mm_a, P_q0_pm_a, r_p_a, t_p_a, r_m_a, t_m_a, R_p_a, T_p_a, R_x_a, T_x_a, S_p_a, S_m_a, d->pade_s, d->pade_r, check_condition, symmetric, thermal, vector, &d->misc_output.pade_condition, work);
     }
*/
#ifdef DEBUG
     else {
         fprintf(stderr, "ERROR: call_solver(): end of if / else if\n");
         exit(1);
     }
#endif

     if (! (d->options & XRTM_OPTION_VECTOR)) {
          dmat_copy(r_m_a, r_p_a, d->n_quad_v_x, d->n_quad_v_x);
          dmat_copy(t_m_a, t_p_a, d->n_quad_v_x, d->n_quad_v_x);
     }
}



/*******************************************************************************
 *
 ******************************************************************************/
static void forward_save_get_layer_R_T_S_U_W_V_free(forward_save_get_layer_R_T_S_U_W_V_data *d);

int forward_save_get_layer_R_T_S_U_W_V_alloc(forward_save_get_layer_R_T_S_U_W_V_data *d, int options, int n_quad_v_x) {

     d->free = (void (*)(void *)) forward_save_get_layer_R_T_S_U_W_V_free;

     d->options = options;

     d->P_q0_mm = alloc_array1_d(n_quad_v_x);
     d->P_q0_pm = alloc_array1_d(n_quad_v_x);

     d->r_p     = alloc_array2_d(n_quad_v_x, n_quad_v_x);
     d->t_p     = alloc_array2_d(n_quad_v_x, n_quad_v_x);

     if (options & XRTM_OPTION_VECTOR) {
          d->r_m = alloc_array2_d(n_quad_v_x, n_quad_v_x);
          d->t_m = alloc_array2_d(n_quad_v_x, n_quad_v_x);
     }

     return 0;
}



static void forward_save_get_layer_R_T_S_U_W_V_free(forward_save_get_layer_R_T_S_U_W_V_data *d) {

     free_array1_d(d->P_q0_mm);
     free_array1_d(d->P_q0_pm);

     free_array2_d(d->r_p);
     free_array2_d(d->t_p);

     if (d->options & XRTM_OPTION_VECTOR) {
          free_array2_d(d->r_m);
          free_array2_d(d->t_m);
     }
}



/*******************************************************************************
 *
 ******************************************************************************/
static int get_layer_R_T_S_U_W_V_a
     (xrtm_data *d, int i_four, int i_layer, int solver,
      double **R_p, double **T_p, double **R_m, double **T_m,
      double *S_p, double *S_m,
      double **R_p_a, double **T_p_a, double **R_m_a, double **T_m_a,
      double *S_p_a, double *S_m_a, int flag, save_tree_data save_tree, work_data *work) {

     double b0;
     double b1;

     double *P_q0_mm_a;
     double *P_q0_pm_a;

     double **r_p_a;
     double **t_p_a;
     double **r_m_a;
     double **t_m_a;

     forward_save_get_layer_R_T_S_U_W_V_data *save;

     work_data work2;


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     save_tree_encode_s(&save_tree, "get_layer_R_T_S_U_W_V");

     save_tree_retrieve_data(&save_tree, forward_save_get_layer_R_T_S_U_W_V_data, &save);


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (d->options & XRTM_OPTION_SAVE_LAYER_R_T_S) {

     }
     else if (flag) {

     }


     work2 = *work;


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     P_q0_mm_a = get_work1(&work2, WORK_DX);
     P_q0_pm_a = get_work1(&work2, WORK_DX);

     r_p_a = get_work1(&work2, WORK_DXX);
     t_p_a = get_work1(&work2, WORK_DXX);
/*
     if (d->options & XRTM_OPTION_VECTOR) {
*/
          r_m_a = get_work1(&work2, WORK_DXX);
          t_m_a = get_work1(&work2, WORK_DXX);
/*
     }
*/

     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     dvec_zero(P_q0_mm_a, d->n_quad_v_x);
     dvec_zero(P_q0_pm_a, d->n_quad_v_x);

     dmat_zero(r_p_a, d->n_quad_v_x, d->n_quad_v_x);
     dmat_zero(t_p_a, d->n_quad_v_x, d->n_quad_v_x);
/*
     if (d->options & XRTM_OPTION_VECTOR) {
*/
          dmat_zero(r_m_a, d->n_quad_v_x, d->n_quad_v_x);
          dmat_zero(t_m_a, d->n_quad_v_x, d->n_quad_v_x);
/*
     }
*/

     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
/*
     if (d->omega[i_layer] == 0. || i_four >= d->n_coef_layer[i_layer])
*/

     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
/*
     else {
*/
          if (d->options & XRTM_OPTION_SOURCE_THERMAL) {
               b0 = d->levels_b[i_layer];
               b1 = d->levels_b[i_layer + 1];
          }

          call_solver_a(d, solver, i_four, b0, b1, d->omega[i_layer], &d->omega_a[i_layer], d->ltau[i_layer], &d->ltau_a[i_layer], d->as_0[i_layer], &d->as_0_a[i_layer], d->atran[i_layer], &d->atran_a[i_layer], save->P_q0_mm, save->P_q0_pm, save->r_p, save->t_p, save->r_m, save->t_m, R_p, T_p, R_m, T_m, S_p, S_m, P_q0_mm_a, P_q0_pm_a, r_p_a, t_p_a, r_m_a, t_m_a, R_p_a, T_p_a, R_m_a, T_m_a, S_p_a, S_m_a, d->derivs.layers_union[i_layer], d->derivs.beam_union[i_layer], save_tree, work2);
if (d->derivs.layers_union[i_layer]) {
          if (get_local_r_t_u_w_a(d, i_four, i_layer, save->r_p, save->t_p, save->r_m, save->t_m, r_p_a, t_p_a, r_m_a, t_m_a, 1, save_tree, &work2)) {
               fprintf(stderr, "ERROR: get_local_r_t_u_w()\n");
               return -1;
          }

          if (get_phase_vecs_q0_a(d, i_four, i_layer, P_q0_mm_a, P_q0_pm_a, 1, &work2)) {
               fprintf(stderr, "ERROR: get_phase_vecs_q0()\n");
               return -1;
          }
}
          update_beam_params_a(d, i_layer, work2);
/*
     }
*/
     update_opt_props_a(d, i_layer, work2);


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
/*
     if (i_four == 0) {
          if (check_R_and_T_norm_a(d->n_quad, d->n_stokes, R_p_a, T_p_a) != 0) {
               fprintf(stderr, "ERROR: check_R_and_T_norm_a()\n");
               return -1;
          }
     }
*/

     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     layer_zero_all(R_m_a, T_m_a, S_m_a, R_p_a, T_p_a, S_p_a, d->n_quad_v_x);


     return 0;
}



/*******************************************************************************
 *
 ******************************************************************************/
static void forward_save_get_total_TOA_R_S_U_V_free(forward_save_get_total_TOA_R_S_U_V_data *d);

int forward_save_get_total_TOA_R_S_U_V_alloc(forward_save_get_total_TOA_R_S_U_V_data *d, int n_quad_v_x, int n_layers) {

     d->free = (void (*)(void *)) forward_save_get_total_TOA_R_S_U_V_free;

     d->R1_m = alloc_array3_d(n_layers, n_quad_v_x, n_quad_v_x);
     d->S1_p = alloc_array2_d(n_layers, n_quad_v_x);

     d->R2_p = alloc_array3_d(n_layers, n_quad_v_x, n_quad_v_x);
     d->T2_p = alloc_array3_d(n_layers, n_quad_v_x, n_quad_v_x);
     d->R2_m = alloc_array3_d(n_layers, n_quad_v_x, n_quad_v_x);
     d->T2_m = alloc_array3_d(n_layers, n_quad_v_x, n_quad_v_x);
     d->S2_p = alloc_array2_d(n_layers, n_quad_v_x);
     d->S2_m = alloc_array2_d(n_layers, n_quad_v_x);

     d->S3_p = alloc_array2_d(n_layers, n_quad_v_x);
     d->S3_m = alloc_array2_d(n_layers, n_quad_v_x);

     return 0;
}



static void forward_save_get_total_TOA_R_S_U_V_free(forward_save_get_total_TOA_R_S_U_V_data *d) {

     free_array3_d(d->R1_m);
     free_array2_d(d->S1_p);

     free_array3_d(d->R2_p);
     free_array3_d(d->T2_p);
     free_array3_d(d->R2_m);
     free_array3_d(d->T2_m);
     free_array2_d(d->S2_p);
     free_array2_d(d->S2_m);

     free_array2_d(d->S3_p);
     free_array2_d(d->S3_m);
}



/*******************************************************************************
 *
 ******************************************************************************/
static int get_total_TOA_R_S_U_V_a
     (xrtm_data *d, int i_four, int solver,
      double **Rs_qq, double *S_s, double **Rs_qq_a, double *S_s_a,
      double **R_m, double *S_p, double **R_m_a, double *S_p_a,
      int surface, int flag, save_tree_data save_tree, work_data *work) {

     uchar *derivs_adding_up_union;

     int i;

     double trans_a;

     double *S1_p_a;

     double *S2_p_a;
     double *S2_m_a;

     double *S3_p_a;
     double *S3_m_a;

     double *S4_p_a;

     double **R1_m_a;

     double **R2_p_a;
     double **T2_p_a;
     double **R2_m_a;
     double **T2_m_a;

     double **R3_m_a;

     forward_save_get_total_TOA_R_S_U_V_data *save;

     work_data work2;
     work_data work3;


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     save_tree_encode_s(&save_tree, "get_total_TOA_R_S_U_V");

     save_tree_retrieve_data(&save_tree, forward_save_get_total_TOA_R_S_U_V_data, &save);


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (flag) {

     }


     work2 = *work;


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     R1_m_a = R_m_a;
     S1_p_a = S_p_a;

     R2_p_a = get_work1(&work2, WORK_DXX);
     T2_p_a = get_work1(&work2, WORK_DXX);
     R2_m_a = get_work1(&work2, WORK_DXX);
     T2_m_a = get_work1(&work2, WORK_DXX);
     S2_p_a = get_work1(&work2, WORK_DX);
     S2_m_a = get_work1(&work2, WORK_DX);

     S3_p_a = get_work1(&work2, WORK_DX);
     S3_m_a = get_work1(&work2, WORK_DX);

     R3_m_a = get_work1(&work2, WORK_DXX);
     S4_p_a = get_work1(&work2, WORK_DX);


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (i_four == 0)
          derivs_adding_up_union = d->derivs.adding_up_union;
     else
          derivs_adding_up_union = d->derivs.adding_up_m_union;


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     dmat_zero(R2_p_a, d->n_quad_v_x, d->n_quad_v_x);
     dmat_zero(T2_p_a, d->n_quad_v_x, d->n_quad_v_x);
     dmat_zero(R2_m_a, d->n_quad_v_x, d->n_quad_v_x);
     dmat_zero(T2_m_a, d->n_quad_v_x, d->n_quad_v_x);
     dvec_zero(S2_p_a, d->n_quad_v_x);
     dvec_zero(S2_m_a, d->n_quad_v_x);

     dvec_zero(S3_p_a, d->n_quad_v_x);
     dvec_zero(S3_m_a, d->n_quad_v_x);

     dmat_zero(R3_m_a, d->n_quad_v_x, d->n_quad_v_x);
     dvec_zero(S4_p_a, d->n_quad_v_x);


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     trans_a = 0.;

     for (i = 0; i < d->n_layers; ++i) {
          work3 = work2;

          save_tree_recode_i(&save_tree, i, i == 0);

          if (! surface && i == d->n_layers - 1)
               layer_copy_ref(R2_m_a, S3_p_a, R1_m_a, S1_p_a, d->n_quad_v_x);
          else {
               layer_add_ref_a(save->R2_m[i], save->T2_m[i], save->S3_m[i], save->R2_p[i], save->T2_p[i], save->S3_p[i], save->R1_m[i], save->S1_p[i], NULL, NULL, R2_m_a, T2_m_a, S3_m_a, R2_p_a, T2_p_a, S3_p_a, R3_m_a, S4_p_a, R1_m_a, S1_p_a, d->n_quad_v_x, 1., &trans_a, derivs_adding_up_union[i], 0, 1, save_tree, work3);

               layer_copy_ref(R1_m_a, S1_p_a, R3_m_a, S4_p_a, d->n_quad_v_x);
               layer_zero_ref(R3_m_a, S4_p_a, d->n_quad_v_x);
          }

          scale_source_vectors_a(d->n_quad_v_x, d->btran[i], &d->btran_a[i], save->S2_p[i], save->S2_m[i], save->S3_p[i], save->S3_m[i], S2_p_a, S2_m_a, S3_p_a, S3_m_a, work3);

          update_beam_params_a(d, i, work3);

          if (get_layer_R_T_S_U_W_V_a(d, i_four, i, solver, save->R2_p[i], save->T2_p[i], save->R2_m[i], save->T2_m[i], save->S2_p[i], save->S2_m[i], R2_p_a, T2_p_a, R2_m_a, T2_m_a, S2_p_a, S2_m_a, 1, save_tree, &work3)) {
               fprintf(stderr, "ERROR: get_layer_R_T_S_U_W_V()\n");
               return -1;
          }
     }

     save_tree_decode_i(&save_tree);

     if (surface) {
          dmat_add(Rs_qq_a, R1_m_a, Rs_qq_a, d->n_quad_v_x, d->n_quad_v_x);
          dvec_add(S_s_a, S1_p_a, S_s_a, d->n_quad_v_x);
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     layer_zero_ref(R_m_a, S_p_a, d->n_quad_v_x);


     return 0;
}



/*******************************************************************************
 *
 ******************************************************************************/
static void forward_save_get_total_R_T_S_U_W_V_free(forward_save_get_total_R_T_S_U_W_V_data *d);

int forward_save_get_total_R_T_S_U_W_V_alloc(forward_save_get_total_R_T_S_U_W_V_data *d, int n_quad_v_x, int n_layers) {

     d->free = (void (*)(void *)) forward_save_get_total_R_T_S_U_W_V_free;

     d->R1_p = alloc_array3_d(n_layers, n_quad_v_x, n_quad_v_x);
     d->T1_p = alloc_array3_d(n_layers, n_quad_v_x, n_quad_v_x);
     d->R1_m = alloc_array3_d(n_layers, n_quad_v_x, n_quad_v_x);
     d->T1_m = alloc_array3_d(n_layers, n_quad_v_x, n_quad_v_x);
     d->S1_p = alloc_array2_d(n_layers, n_quad_v_x);
     d->S1_m = alloc_array2_d(n_layers, n_quad_v_x);

     d->R2_p = alloc_array3_d(n_layers, n_quad_v_x, n_quad_v_x);
     d->T2_p = alloc_array3_d(n_layers, n_quad_v_x, n_quad_v_x);
     d->R2_m = alloc_array3_d(n_layers, n_quad_v_x, n_quad_v_x);
     d->T2_m = alloc_array3_d(n_layers, n_quad_v_x, n_quad_v_x);
     d->S2_p = alloc_array2_d(n_layers, n_quad_v_x);
     d->S2_m = alloc_array2_d(n_layers, n_quad_v_x);

     d->S3_p = alloc_array2_d(n_layers, n_quad_v_x);
     d->S3_m = alloc_array2_d(n_layers, n_quad_v_x);

     return 0;
}



static void forward_save_get_total_R_T_S_U_W_V_free(forward_save_get_total_R_T_S_U_W_V_data *d) {

     free_array3_d(d->R1_p);
     free_array3_d(d->T1_p);
     free_array3_d(d->R1_m);
     free_array3_d(d->T1_m);
     free_array2_d(d->S1_p);
     free_array2_d(d->S1_m);

     free_array3_d(d->R2_p);
     free_array3_d(d->T2_p);
     free_array3_d(d->R2_m);
     free_array3_d(d->T2_m);
     free_array2_d(d->S2_p);
     free_array2_d(d->S2_m);

     free_array2_d(d->S3_p);
     free_array2_d(d->S3_m);
}



/*******************************************************************************
 *
 ******************************************************************************/
static int get_total_R_T_S_U_W_V_a
     (xrtm_data *d, int i_four, int solver,
      double **R_p, double **T_p, double **R_m, double **T_m,
      double *S_p, double *S_m,
      double **R_p_a, double **T_p_a, double **R_m_a, double **T_m_a,
      double *S_p_a, double *S_m_a,
      int flag, save_tree_data save_tree, work_data *work) {

     int i;

     double trans_a;

     double *S1_p_a;
     double *S1_m_a;

     double *S2_p_a;
     double *S2_m_a;

     double *S3_p_a;
     double *S3_m_a;

     double *S4_p_a;
     double *S4_m_a;

     double **R1_p_a;
     double **T1_p_a;
     double **R1_m_a;
     double **T1_m_a;

     double **R2_p_a;
     double **T2_p_a;
     double **R2_m_a;
     double **T2_m_a;

     double **R3_p_a;
     double **T3_p_a;
     double **R3_m_a;
     double **T3_m_a;

     forward_save_get_total_R_T_S_U_W_V_data *save;

     work_data work2;
     work_data work3;


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     save_tree_encode_s(&save_tree, "get_total_R_T_S_U_W_V");

     save_tree_retrieve_data(&save_tree, forward_save_get_total_R_T_S_U_W_V_data, &save);


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (d->options & XRTM_OPTION_SAVE_TOTAL_R_T_S) {

     }
     else if (flag) {

     }


     work2 = *work;


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     R1_p_a = R_p_a;     T1_p_a = T_p_a;
     R1_m_a = R_m_a;     T1_m_a = T_m_a;
     S1_p_a = S_p_a;     S1_m_a = S_m_a;

     R2_p_a = get_work1(&work2, WORK_DXX);
     T2_p_a = get_work1(&work2, WORK_DXX);
     R2_m_a = get_work1(&work2, WORK_DXX);
     T2_m_a = get_work1(&work2, WORK_DXX);
     S2_p_a = get_work1(&work2, WORK_DX);
     S2_m_a = get_work1(&work2, WORK_DX);

     S3_p_a = get_work1(&work2, WORK_DX);
     S3_m_a = get_work1(&work2, WORK_DX);

     R3_p_a = get_work1(&work2, WORK_DXX);
     T3_p_a = get_work1(&work2, WORK_DXX);
     R3_m_a = get_work1(&work2, WORK_DXX);
     T3_m_a = get_work1(&work2, WORK_DXX);
     S4_p_a = get_work1(&work2, WORK_DX);
     S4_m_a = get_work1(&work2, WORK_DX);


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     dmat_zero(R2_p_a, d->n_quad_v_x, d->n_quad_v_x);
     dmat_zero(T2_p_a, d->n_quad_v_x, d->n_quad_v_x);
     dmat_zero(R2_m_a, d->n_quad_v_x, d->n_quad_v_x);
     dmat_zero(T2_m_a, d->n_quad_v_x, d->n_quad_v_x);
     dvec_zero(S2_p_a, d->n_quad_v_x);
     dvec_zero(S2_m_a, d->n_quad_v_x);

     dmat_zero(R3_p_a, d->n_quad_v_x, d->n_quad_v_x);
     dmat_zero(T3_p_a, d->n_quad_v_x, d->n_quad_v_x);
     dmat_zero(R3_m_a, d->n_quad_v_x, d->n_quad_v_x);
     dmat_zero(T3_m_a, d->n_quad_v_x, d->n_quad_v_x);

     dvec_zero(S3_p_a, d->n_quad_v_x);
     dvec_zero(S3_m_a, d->n_quad_v_x);

     dvec_zero(S4_p_a, d->n_quad_v_x);
     dvec_zero(S4_m_a, d->n_quad_v_x);


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     trans_a = 0.;

     for (i = d->n_layers - 1; i >= 0; --i) {
          work3 = work2;

          save_tree_recode_i(&save_tree, i, i == d->n_layers - 1);

          if (i == 0)
               layer_copy_all(R2_m_a, T2_m_a, S3_m_a, R2_p_a, T2_p_a, S3_p_a,
                              R1_m_a, T1_m_a, S1_m_a, R1_p_a, T1_p_a, S1_p_a,
                              d->n_quad_v_x);
          else {
               layer_add_all_a(save->R1_m[i], save->T1_m[i], save->S1_m[i], save->R1_p[i], save->T1_p[i], save->S1_p[i], save->R2_m[i], save->T2_m[i], save->S3_m[i], save->R2_p[i], save->T2_p[i], save->S3_p[i], NULL, NULL, NULL, NULL, NULL, NULL, R3_m_a, T3_m_a, S4_m_a, R3_p_a, T3_p_a, S4_p_a, R2_m_a, T2_m_a, S3_m_a, R2_p_a, T2_p_a, S3_p_a, R1_m_a, T1_m_a, S1_m_a, R1_p_a, T1_p_a, S1_p_a, d->n_quad_v_x, 1., &trans_a, d->derivs.adding_down_union[i], save_tree, work3);

               layer_copy_all(R1_m_a, T1_m_a, S1_m_a, R1_p_a, T1_p_a, S1_p_a, R3_m_a, T3_m_a, S4_m_a, R3_p_a, T3_p_a, S4_p_a, d->n_quad_v_x);
               layer_zero_all(R3_m_a, T3_m_a, S4_m_a, R3_p_a, T3_p_a, S4_p_a, d->n_quad_v_x);
          }

          scale_source_vectors_a(d->n_quad_v_x, d->btran[i], &d->btran_a[i], save->S2_p[i], save->S2_m[i], save->S3_p[i], save->S3_m[i], S2_p_a, S2_m_a, S3_p_a, S3_m_a, work3);

          update_beam_params_a(d, i, work3);

          if (get_layer_R_T_S_U_W_V_a(d, i_four, i, solver, save->R2_p[i], save->T2_p[i], save->R2_m[i], save->T2_m[i], save->S2_p[i], save->S2_m[i], R2_p_a, T2_p_a, R2_m_a, T2_m_a, S2_p_a, S2_m_a, 1, save_tree, &work3)) {
               fprintf(stderr, "ERROR: get_layer_R_T_S_U_W_V_a()\n");
               return -1;
          }
     }

     save_tree_decode_i(&save_tree);


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     layer_zero_all(R_m_a, T_m_a, S_m_a, R_p_a, T_p_a, S_p_a, d->n_quad_v_x);


     return 0;
}



/*******************************************************************************
 *
 ******************************************************************************/
static void brdf_vecs_get_a(xrtm_data *d, int i_four, int i_offset, int n_mus, int j_offset, double **kernel_f, double **kernel_f_a, double *alpha1, double *R, double *R_a, work_data work) {

     int n_mus_v;

     if (d->misc_input.use_symmetric_form) {
          n_mus_v = n_mus * d->n_stokes;

          vec_sim_trans(n_mus_v, R_a, alpha1);
     }

     init_array2_d(kernel_f_a, d->n_kernel_quad * 2 + 1, n_mus, 0.);

     brdf_build_ref_vec_a(i_four, i_offset, n_mus, j_offset, d->n_stokes, d->qf, d->qx_v, d->qw_v, d->n_kernels, d->n_kernel_quad, d->kernels, d->kernel_qx, d->kernel_qw, kernel_f_a, R_a, work);
/*
     if (i_four == 0)
*/
          brdf_build_kernel_vecs_a(i_offset, n_mus, j_offset, d->n_stokes, d->n_kernels, d->n_kernel_quad, d->kernels, d->kernel_ampfac, d->kernel_ampfac_a, d->kernel_params, d->kernel_params_a, d->kernel_qx, d->kernel_qw, &d->kernel_aux, kernel_f, kernel_f_a, work);
}



/*******************************************************************************
 *
 ******************************************************************************/
static void get_brdf_vecs_q0_a(xrtm_data *d, int i_four, double *Rs_q0, double *Rs_q0_a, work_data work) {

     brdf_vecs_get_a(d, i_four, 0, d->n_quad_x, d->n_quad + d->n_umus, d->kernel_f_q0, d->kernel_f_q0_a, d->alpha1, Rs_q0, Rs_q0_a, work);
}



/*******************************************************************************
 *
 ******************************************************************************/
static void get_brdf_vecs_u0_a(xrtm_data *d, int i_four, double *Rs_u0, double *Rs_u0_a, work_data work) {

     brdf_vecs_get_a(d, i_four, d->n_quad, d->n_umus, d->n_quad + d->n_umus, d->kernel_f_u0, d->kernel_f_u0_a, d->alpha1 + d->n_quad_v, Rs_u0, Rs_u0_a, work);
}



/*******************************************************************************
 *
 ******************************************************************************/
static void brdf_mats_get_a(xrtm_data *d, int i_four, int i_offset, int n_mus1, int j_offset, int n_mus2, double ***kernel_f, double ***kernel_f_a, double *alpha1, double *alpha2, double **R, double **R_a, work_data work) {

     int n_mus_v1;
     int n_mus_v2;

     if (d->misc_input.use_symmetric_form) {
          n_mus_v1 = n_mus1 * d->n_stokes;
          n_mus_v2 = n_mus2 * d->n_stokes;

          mat_sim_trans(n_mus_v1, n_mus_v2, R_a, alpha1, alpha2);
     }

     init_array3_d(kernel_f_a, d->n_kernel_quad * 2 + 1, n_mus1, n_mus2, 0.);

     brdf_build_ref_mat_a(i_four, i_offset, n_mus1, j_offset, n_mus2, d->n_stokes, d->qf, d->qx_v, d->qw_v, d->n_kernels, d->n_kernel_quad, d->kernels, d->kernel_qx, d->kernel_qw, kernel_f_a, R_a, work);
/*
     if (i_four == 0)
*/
          brdf_build_kernel_mats_a(i_offset, n_mus1, j_offset, n_mus2, d->n_stokes, d->n_kernels, d->n_kernel_quad, d->kernels, d->kernel_ampfac, d->kernel_ampfac_a, d->kernel_params, d->kernel_params_a, d->kernel_qx, d->kernel_qw, &d->kernel_aux, kernel_f, kernel_f_a, work);
}



/*******************************************************************************
 *
 ******************************************************************************/
static void get_brdf_mats_qq_a(xrtm_data *d, int i_four, double **Rs_qq, double **Rs_qq_a, work_data work) {

     brdf_mats_get_a(d, i_four, 0, d->n_quad_x, 0, d->n_quad_x, d->kernel_f_qq, d->kernel_f_qq_a, d->alpha1, d->alpha2, Rs_qq, Rs_qq_a, work);
}



/*******************************************************************************
 *
 ******************************************************************************/
static void get_brdf_mats_uq_a(xrtm_data *d, int i_four, double **Rs_uq, double **Rs_uq_a, work_data work) {

     brdf_mats_get_a(d, i_four, d->n_quad, d->n_umus, 0, d->n_quad_x, d->kernel_f_uq, d->kernel_f_uq_a, d->alpha1 + d->n_quad_v, d->alpha2, Rs_uq, Rs_uq_a, work);
}



/*******************************************************************************
 *
 ******************************************************************************/
static void diff_bound_input_init_a(xrtm_data *d, int i_four, double **In_p_a, double **I1_m_a, work_data work) {

     init_array1_d(In_p_a[i_four], d->n_quad_v + d->n_umus_v, 0.);
     init_array1_d(I1_m_a[i_four], d->n_quad_v + d->n_umus_v, 0.);
}


#ifdef INCLUDE_DEV_SOURCE
static void init_diff_bound_input_a(xrtm_data *d, int i_four, work_data work) {

     diff_bound_input_init_a(d, i_four, d->In_p_a, d->I1_m_a, work);
}
#endif
/*
static
*/
void init_diff_bound_input_all_a(xrtm_data *d, int n_four, work_data work) {

     int i;

     for (i = 0; i < n_four; ++i)
          diff_bound_input_init_a(d, i, d->In_p_a, d->I1_m_a, work);
}



/*******************************************************************************
 *
 ******************************************************************************/
static void forward_save_update_diff_bound_input_free(forward_save_update_diff_bound_input_data *d);

int forward_save_update_diff_bound_input_alloc(forward_save_update_diff_bound_input_data *d, int n_quad_v, int n_umus_v) {

     d->free = (void (*)(void *)) forward_save_update_diff_bound_input_free;

     d->Rs_q0 = alloc_array1_d(n_quad_v + n_umus_v);
     d->Rs_u0 = alloc_array1_d(n_quad_v + n_umus_v);

     return 0;
}



static void forward_save_update_diff_bound_input_free(forward_save_update_diff_bound_input_data *d) {

     free_array1_d(d->Rs_q0);
     free_array1_d(d->Rs_u0);
}


/*******************************************************************************
 *
 ******************************************************************************/
static void diff_bound_input_update_a(xrtm_data *d, int i_four, double btran, double *btran_a, double **In_p, double **I1_m, double **In_p_a, double **I1_m_a, save_tree_data save_tree, work_data work) {

     int i;
     int ii;

     double a;

     double *Rs_q0_a;
     double *Rs_u0_a;

     forward_save_update_diff_bound_input_data *save;


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (d->F_0 != 0. && solar_surface_for_i_four(d, i_four)) {
          save_tree_encode_s(&save_tree, "update_diff_bound_input");
          save_tree_retrieve_data(&save_tree, forward_save_update_diff_bound_input_data, &save);

          Rs_q0_a = get_work1(&work, WORK_DX);

          if (d->options & XRTM_OPTION_SFI && d->n_umus > 0)
               Rs_u0_a = get_work1(&work, WORK_DU);

          a = d->qf * d->mu_0 * d->F_0 / PI;

          dvec_zero(Rs_q0_a, d->n_quad_v_x);

          if (d->options & XRTM_OPTION_SFI && d->n_umus > 0)
               dvec_zero(Rs_u0_a, d->n_umus_v);

          for (i = 0; i < d->n_quad_v_x; ++i) {
              *btran_a    += In_p_a[i_four][i] * a * save->Rs_q0[i];
               Rs_q0_a[i] += In_p_a[i_four][i] * a * btran;
          }

          for (ii = 0; i < d->n_quad_v + d->n_umus_v; ++i, ++ii) {
              *btran_a     += In_p_a[i_four][i] * a * save->Rs_u0[ii];
               Rs_u0_a[ii] += In_p_a[i_four][i] * a * btran;
          }

          get_brdf_vecs_q0_a(d, i_four, save->Rs_q0, Rs_q0_a, work);

          if (d->options & XRTM_OPTION_SFI && d->n_umus > 0)
               get_brdf_vecs_u0_a(d, i_four, save->Rs_u0, Rs_u0_a, work);
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (d->options & XRTM_OPTION_SOURCE_THERMAL && thermal_surface_for_i_four(d, i_four)) {

     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     init_array1_d(In_p_a[i_four], d->n_quad_v + d->n_umus_v, 0.);
     init_array1_d(I1_m_a[i_four], d->n_quad_v + d->n_umus_v, 0.);
}



static void update_diff_bound_input_a(xrtm_data *d, int i_four, save_tree_data save_tree, work_data work) {

     diff_bound_input_update_a(d, i_four, d->btran[d->n_layers], &d->btran_a[d->n_layers], d->In_p, d->I1_m, d->In_p_a, d->I1_m_a, save_tree, work);

     update_beam_params_a(d, d->n_layers, work);
}



/*******************************************************************************
 *
 ******************************************************************************/
static void forward_save_get_single_free(forward_save_get_single_data *d);

int forward_save_get_single_alloc(forward_save_get_single_data *d, int options, int n_coef, int n_quad_v, int n_umus_v, int n_mus2, int n_phis, int n_layers, int n_stokes) {

     d->free = (void (*)(void *)) forward_save_get_single_free;

     d->options = options;

     if (options & XRTM_OPTION_PHASE_SCALAR) {
          d->polys_up = (void ***) alloc_array3_d(n_mus2, n_phis, n_coef);
          d->polys_dn = (void ***) alloc_array3_d(n_mus2, n_phis, n_coef);
     }
     else
     if (options & XRTM_OPTION_PHASE_MATRIX_GC) {
          d->polys_up = (void ***) alloc_array4_d(n_mus2, n_phis, 4, n_coef);
          d->polys_dn = (void ***) alloc_array4_d(n_mus2, n_phis, 4, n_coef);
     }
     else
     if (options & XRTM_OPTION_PHASE_MATRIX_LC) {
          d->polys_up = (void ***) alloc_array3_d(n_mus2, n_phis, n_coef);
          d->polys_dn = (void ***) alloc_array3_d(n_mus2, n_phis, n_coef);
     }

     d->P_full_up = alloc_array4_d(n_mus2, n_phis, n_layers, n_stokes);
     d->P_full_dn = alloc_array4_d(n_mus2, n_phis, n_layers, n_stokes);

     return 0;
}



static void forward_save_get_single_free(forward_save_get_single_data *d) {

     if (d->options & XRTM_OPTION_PHASE_SCALAR) {
          free_array3_d((double ***) d->polys_up);
          free_array3_d((double ***) d->polys_dn);
     }
     else
     if (d->options & XRTM_OPTION_PHASE_MATRIX_GC) {
          free_array4_d((double ****) d->polys_up);
          free_array4_d((double ****) d->polys_dn);
     }
     else
     if (d->options & XRTM_OPTION_PHASE_MATRIX_LC) {
          free_array3_d((double ***) d->polys_up);
          free_array3_d((double ***) d->polys_dn);
     }

     free_array4_d(d->P_full_up);
     free_array4_d(d->P_full_dn);
}



/*******************************************************************************
 *
 ******************************************************************************/
static void get_single_a(xrtm_data *d, int n_phis, double **phis, double ****I_p, double ****I_m, double ****I_p_a, double ****I_m_a, save_tree_data save_tree, work_data work) {

     int i;
     int ii;
     int j;
     int k;
/*
     int n_quad_x;
*/
     int i_mus2_v;
     int n_mus2;

     double a;

     double f;
     double f_a;

     double *mus2;

     double *omega2;
     double *omega2_a;

     double **P_full_a;

     double **I_a;

     forward_save_get_single_data *save;


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
/*
     n_quad_x = d->n_quad + d->n_umus;
*/
     if (d->n_umus == 0) {
          mus2     = d->qx;
          i_mus2_v = 0;
          n_mus2   = d->n_quad;
     }
     else {
          mus2     = d->umus;
          i_mus2_v = d->n_quad_v;
          n_mus2   = d->n_umus;
     }
/*
     else {
          mus2     = d->qx;
          i_mus2_v = 0;
          n_mus2   = d->n_quad + d->n_umus;
     }
*/

     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     save_tree_encode_s(&save_tree, "get_single");

     save_tree_retrieve_data(&save_tree, forward_save_get_single_data, &save);


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     omega2   = get_work1(&work, WORK_DLAYERS);

     omega2_a = get_work1(&work, WORK_DLAYERS);

     P_full_a = get_work2(&work, WORK_DS, WORK_LAYERS_V, NULL);

     I_a      = get_work_d2(&work, d->n_ulevels, d->n_stokes);


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
/*
     update_utaus(d, work);
*/

     if (! (d->options & XRTM_OPTION_N_T_TMS))
           omega2 = d->omega0;
     else {
          a = 2. * d->n_coef2 + 1.;

          for (i = 0; i < d->n_layers; ++i) {
               get_delta_m_f(d, i, 0, 0, a, &f, NULL);

               n_t_tms_scaling(0, f, NULL, d->omega0[i], &omega2[i], NULL, NULL);
          }
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     init_array1_d(omega2_a, d->n_layers, 0.);


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     ii = i_mus2_v;

     for (i = 0; i < n_mus2; ++i) {
          for (j = 0; j < n_phis; ++j) {
               save_tree_recode_i_j(&save_tree, i, j, i == 0 && j == 0);


               for (k = 0; k < d->n_layers; ++k)
                    init_array1_d(P_full_a[k], d->n_stokes, 0.);

               for (k = 0; k < d->n_ulevels; ++k)
                    dvec_copy(I_a[k], I_p_a[k][i][j], d->n_stokes);

               single_scattered_radiance_up_a(d->n_stokes, d->n_layers, d->F_0, d->n_ulevels, d->ulevels, d->utaus2, &mus2[i], 1, omega2, omega2_a, d->ltau, d->ltau_a, d->btran, d->btran_a, d->as_0, d->as_0_a, d->atran, d->atran_a, save->P_full_up[i][j], P_full_a, d->In_p[0]+ii, d->In_p_a[0]+ii, NULL, I_a, d->options & XRTM_OPTION_OUTPUT_AT_TAUS, save_tree, work);

               get_phase_func_full_all_a(d, -mus2[i], phis[i][j], save->polys_up[i][j], P_full_a, work);


               for (k = 0; k < d->n_layers; ++k)
                    init_array1_d(P_full_a[k], d->n_stokes, 0.);

               for (k = 0; k < d->n_ulevels; ++k)
                    dvec_copy(I_a[k], I_m_a[k][i][j], d->n_stokes);

               single_scattered_radiance_dn_a(d->n_stokes, d->n_layers, d->F_0, d->n_ulevels, d->ulevels, d->utaus2, &mus2[i], 1, omega2, omega2_a, d->ltau, d->ltau_a, d->btran, d->btran_a, d->as_0, d->as_0_a, d->atran, d->atran_a, save->P_full_dn[i][j], P_full_a, d->I1_m[0]+ii, d->I1_m_a[0]+ii, NULL, I_a, d->options & XRTM_OPTION_OUTPUT_AT_TAUS, save_tree, work);

               get_phase_func_full_all_a(d,  mus2[i], phis[i][j], save->polys_dn[i][j], P_full_a, work);
          }

          ii += d->n_stokes;
     }

     save_tree_decode_i_j(&save_tree);


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (! (d->options & XRTM_OPTION_N_T_TMS)) {
           for (i = 0; i < d->n_layers; ++i)
                d->omega0_a[i] += omega2_a[i];
     }
     else {
          a = 2. * d->n_coef2 + 1.;

          f_a = 0.;

          for (i = 0; i < d->n_layers; ++i) {
               get_delta_m_f(d, i, 0, 0, a, &f, NULL);

               n_t_tms_scaling_a(f, &f_a, d->omega0[i], omega2[i], &d->omega0_a[i], omega2_a[i]);

               get_delta_m_f_a(d, i, a, f, &f_a);
          }
     }

     update_diff_bound_input_a(d, 0, save_tree, work);

     update_beam_params_all_a(d, work);

     update_opt_props_all_a(d, work);


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     init_array4_d(I_p_a, d->n_ulevels, n_mus2, n_phis, d->n_stokes, 0.);
     init_array4_d(I_m_a, d->n_ulevels, n_mus2, n_phis, d->n_stokes, 0.);
}



/*******************************************************************************
 *
 ******************************************************************************/
static void forward_save_fourier_get_add_bottom_up_free(forward_save_fourier_get_add_bottom_up_data *d);

int forward_save_fourier_get_add_bottom_up_alloc(forward_save_fourier_get_add_bottom_up_data *d, int n_quad_v_x) {

     d->free = (void (*)(void *)) forward_save_fourier_get_add_bottom_up_free;

     d->R_m   = alloc_array2_d(n_quad_v_x, n_quad_v_x);
     d->S_p   = alloc_array1_d(n_quad_v_x);

     d->Rs_qq = alloc_array2_d(n_quad_v_x, n_quad_v_x);

     return 0;
}



static void forward_save_fourier_get_add_bottom_up_free(forward_save_fourier_get_add_bottom_up_data *d) {

     free_array2_d(d->R_m);
     free_array1_d(d->S_p);

     free_array2_d(d->Rs_qq);
}


/*******************************************************************************
 *
 ******************************************************************************/
static int fourier_get_add_bottom_up_a(xrtm_data *d, int i_four, int solver, double **I_p, double **I_m, double **I_p_a, double **I_m_a, save_tree_data save_tree, work_data work) {

     int surface;

     double **Rs_qq_a;

     double *S_p_a;

     double **R_m_a;

     forward_save_fourier_get_add_bottom_up_data *save;


     surface = solar_surface_for_i_four(d, i_four);


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     save_tree_encode_s(&save_tree, "fourier_get_add_bottom_up");

     save_tree_retrieve_data(&save_tree, forward_save_fourier_get_add_bottom_up_data, &save);


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     R_m_a = get_work1(&work, WORK_DXX);
     S_p_a = get_work1(&work, WORK_DX);

     if (surface)
          Rs_qq_a = get_work1(&work, WORK_DXX);


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     dmat_zero(R_m_a, d->n_quad_v_x, d->n_quad_v_x);
     dvec_zero(S_p_a, d->n_quad_v_x);

     if (surface)
          dmat_zero(Rs_qq_a, d->n_quad_v_x, d->n_quad_v_x);


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     dvec_add(d->I1_m_a[i_four], I_m_a[0], d->I1_m_a[i_four], d->n_quad_v_x);

     radiance_toa_ref_a(d->n_quad_v_x, save->R_m, save->S_p, R_m_a, S_p_a, d->I1_m[i_four], d->I1_m_a[i_four], I_p[0], I_p_a[0], work);


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if ( ! (d->options & XRTM_OPTION_CALC_DERIVS) || ! (d->options & XRTM_OPTION_STACK_REUSE_ADDING)) {
          if (get_total_TOA_R_S_U_V_a(d, i_four, solver, save->Rs_qq, d->In_p[i_four], Rs_qq_a, d->In_p_a[i_four], save->R_m, save->S_p, R_m_a, S_p_a, surface, 1, save_tree, &work)) {
               fprintf(stderr, "ERROR: get_total_TOA_R_S_U_V_a()\n");
               return -1;
          }
     }
     else {

     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (surface)
          get_brdf_mats_qq_a(d, i_four, NULL, Rs_qq_a, work);


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     update_diff_bound_input_a(d, i_four, save_tree, work);


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     init_array2_d(I_p_a, d->n_ulevels, d->n_quad_v_x, 0.);
     init_array2_d(I_m_a, d->n_ulevels, d->n_quad_v_x, 0.);


     return 0;
}



/*******************************************************************************
 *
 ******************************************************************************/
static void forward_save_fourier_get_add_both_ways_free(forward_save_fourier_get_add_both_ways_data *d);

int forward_save_fourier_get_add_both_ways_alloc(forward_save_fourier_get_add_both_ways_data *d, int n_quad_v_x) {

     d->free = (void (*)(void *)) forward_save_fourier_get_add_both_ways_free;

     d->R_p   = alloc_array2_d(n_quad_v_x, n_quad_v_x);
     d->T_p   = alloc_array2_d(n_quad_v_x, n_quad_v_x);
     d->R_m   = alloc_array2_d(n_quad_v_x, n_quad_v_x);
     d->T_m   = alloc_array2_d(n_quad_v_x, n_quad_v_x);
     d->S_p   = alloc_array1_d(n_quad_v_x);
     d->S_m   = alloc_array1_d(n_quad_v_x);

     d->I_p1  = alloc_array1_d(n_quad_v_x);
     d->I_m1  = alloc_array1_d(n_quad_v_x);

     d->Rs_qq = alloc_array2_d(n_quad_v_x, n_quad_v_x);

     return 0;
}



static void forward_save_fourier_get_add_both_ways_free(forward_save_fourier_get_add_both_ways_data *d) {

     free_array2_d(d->R_p);
     free_array2_d(d->T_p);
     free_array2_d(d->R_m);
     free_array2_d(d->T_m);
     free_array1_d(d->S_p);
     free_array1_d(d->S_m);

     free_array1_d(d->I_p1);
     free_array1_d(d->I_m1);

     free_array2_d(d->Rs_qq);
}



/*******************************************************************************
 *
 ******************************************************************************/
static int fourier_get_add_both_ways_a(xrtm_data *d, int i_four, int solver, double **I_p, double **I_m, double **I_p_a, double **I_m_a, save_tree_data save_tree, work_data work) {

     int i;

     int surface;

     double **Rs_qq_a;

     double *S_p_a;
     double *S_m_a;

     double **R_p_a;
     double **T_p_a;
     double **R_m_a;
     double **T_m_a;

     double *I_p1_a;
     double *I_m1_a;

     forward_save_fourier_get_add_both_ways_data *save;


     surface = solar_surface_for_i_four(d, i_four);


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     save_tree_encode_s(&save_tree, "fourier_get_add_both_ways");

     save_tree_retrieve_data(&save_tree, forward_save_fourier_get_add_both_ways_data, &save);


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     R_p_a = get_work1(&work, WORK_DXX);
     T_p_a = get_work1(&work, WORK_DXX);
     R_m_a = get_work1(&work, WORK_DXX);
     T_m_a = get_work1(&work, WORK_DXX);
     S_p_a = get_work1(&work, WORK_DX);
     S_m_a = get_work1(&work, WORK_DX);

     if (surface)
          Rs_qq_a = get_work1(&work, WORK_DXX);


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     dmat_zero(R_p_a, d->n_quad_v_x, d->n_quad_v_x);
     dmat_zero(T_p_a, d->n_quad_v_x, d->n_quad_v_x);
     dmat_zero(R_m_a, d->n_quad_v_x, d->n_quad_v_x);
     dmat_zero(T_m_a, d->n_quad_v_x, d->n_quad_v_x);
     dvec_zero(S_p_a, d->n_quad_v_x);
     dvec_zero(S_m_a, d->n_quad_v_x);

     if (surface)
          dmat_zero(Rs_qq_a, d->n_quad_v_x, d->n_quad_v_x);


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (! surface) {
          for (i = 0; i < d->n_ulevels; ++i) {
               if (d->ulevels[i] == 0) {
                    dvec_add(d->I1_m_a[i_four], I_m_a[i], d->I1_m_a[i_four], d->n_quad_v_x);

                    radiance_slab_a(d->n_quad_v_x, save->R_m, save->T_p, save->S_p, R_m_a, T_p_a, S_p_a, d->I1_m[i_four], d->In_p[i_four], d->I1_m_a[i_four], d->In_p_a[i_four], I_p[i], I_p_a[i], work);
               }
               else
               if (d->ulevels[i] == d->n_layers)
                    radiance_slab_a(d->n_quad_v_x, save->R_p, save->T_m, save->S_m, R_p_a, T_m_a, S_m_a, d->I1_m[i_four], d->In_p[i_four], d->I1_m_a[i_four], d->In_p_a[i_four], I_m[i], I_m_a[i], work);
#ifdef DEBUG
               else {
                   fprintf(stderr, "ERROR: fourier_get_add_both_ways(): end of if / else if\n");
                   exit(1);
               }
#endif
          }
     }
     else {
          if (d->ulevels[0] == d->n_layers) {
               I_p1_a = I_p_a[0];
               I_m1_a = I_m_a[0];
          }
          else
          if (d->n_ulevels > 1 &&
              d->ulevels[1] == d->n_layers) {
               I_p1_a = I_p_a[1];
               I_m1_a = I_m_a[1];
          }
          else {
               I_p1_a = get_work_d1(&work, d->n_quad_v_x);
               init_array1_d(I_p1_a, d->n_quad_v_x, 0.);
               I_m1_a = get_work_d1(&work, d->n_quad_v_x);
               init_array1_d(I_m1_a, d->n_quad_v_x, 0.);
          }

          if (d->ulevels[0] == 0)
               radiance_toa_all_a(d->n_quad_v_x, save->R_m, save->T_p, save->S_p, R_m_a, T_p_a, S_p_a, d->I1_m[i_four], save->I_p1, d->I1_m_a[i_four], I_p1_a, I_p[0], I_m[0], I_p_a[0], I_m_a[0], work);

          radiance_boa_all_a(d->n_quad_v_x, save->R_p, save->T_m, save->Rs_qq, save->S_m, R_p_a, T_m_a, Rs_qq_a, S_m_a, d->In_p[i_four], d->I1_m[i_four], d->In_p_a[i_four], d->I1_m_a[i_four], save->I_p1, save->I_m1, I_p1_a, I_m1_a, save_tree, work);
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (surface)
          get_brdf_mats_qq_a(d, i_four, NULL, Rs_qq_a, work);


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     update_diff_bound_input_a(d, i_four, save_tree, work);


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if ( ! (d->options & XRTM_OPTION_CALC_DERIVS) ||
          ! (d->options & XRTM_OPTION_STACK_REUSE_ADDING)) {
          if (get_total_R_T_S_U_W_V_a(d, i_four, solver, save->R_p, save->T_p, save->R_m, save->T_m, save->S_p, save->S_m, R_p_a, T_p_a, R_m_a, T_m_a, S_p_a, S_m_a, 1, save_tree, &work)) {
               fprintf(stderr, "ERROR: get_total_R_T_S_U_W_V_a()\n");
               return -1;
          }
     }
     else {

     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     init_array2_d(I_p_a, d->n_ulevels, d->n_quad_v_x, 0.);
     init_array2_d(I_m_a, d->n_ulevels, d->n_quad_v_x, 0.);


     return 0;
}



/*******************************************************************************
 *
 ******************************************************************************/
static void forward_save_fourier_get_bvp_free(forward_save_fourier_get_bvp_data *d);

int forward_save_fourier_get_bvp_alloc(forward_save_fourier_get_bvp_data *d, int options, int solver, int surface, int n_layers, int n_quad_v_x, int n_umus_v) {

     d->free = (void (*)(void *)) forward_save_fourier_get_bvp_free;

     d->options  = options;
     d->solver   = solver;
     d->surface  = surface;
     d->n_umus_v = n_umus_v;

     d->P_q0_mm = alloc_array2_d(n_layers, n_quad_v_x);
     d->P_q0_pm = alloc_array2_d(n_layers, n_quad_v_x);

     if (solver & XRTM_SOLVER_SOS || (solver & XRTM_SOLVER_TWO_OS && ! (options & XRTM_OPTION_SFI))) {
          d->P_qq_pp = alloc_array3_d(n_layers, n_quad_v_x, n_quad_v_x);
          d->P_qq_mp = alloc_array3_d(n_layers, n_quad_v_x, n_quad_v_x);

          if (options & XRTM_OPTION_VECTOR) {
               d->P_qq_mm = alloc_array3_d(n_layers, n_quad_v_x, n_quad_v_x);
               d->P_qq_pm = alloc_array3_d(n_layers, n_quad_v_x, n_quad_v_x);
          }
     }

     if (options & XRTM_OPTION_SFI && n_umus_v > 0) {
          d->P_u0_mm = alloc_array2_d(n_layers, n_umus_v);
          d->P_u0_pm = alloc_array2_d(n_layers, n_umus_v);

          d->P_uq_pp = alloc_array3_d(n_layers, n_umus_v, n_quad_v_x);
          d->P_uq_mp = alloc_array3_d(n_layers, n_umus_v, n_quad_v_x);

          if (options & XRTM_OPTION_VECTOR) {
               d->P_uq_mm = alloc_array3_d(n_layers, n_umus_v, n_quad_v_x);
               d->P_uq_pm = alloc_array3_d(n_layers, n_umus_v, n_quad_v_x);
          }
     }

     d->r_p = alloc_array3_d(n_layers, n_umus_v, n_quad_v_x);
     d->t_p = alloc_array3_d(n_layers, n_umus_v, n_quad_v_x);

     if (options & XRTM_OPTION_VECTOR) {
          d->r_m = alloc_array3_d(n_layers, n_umus_v, n_quad_v_x);
          d->t_m = alloc_array3_d(n_layers, n_umus_v, n_quad_v_x);
     }

     if (surface) {
          d->Rs_q0 = alloc_array1_d(n_quad_v_x);
          d->Rs_qq = alloc_array2_d(n_quad_v_x, n_quad_v_x);

          if (n_umus_v > 0) {
               d->Rs_u0 = alloc_array1_d(n_umus_v);
               d->Rs_uq = alloc_array2_d(n_umus_v, n_quad_v_x);
          }
     }

     return 0;
}



static void forward_save_fourier_get_bvp_free(forward_save_fourier_get_bvp_data *d) {

     free_array2_d(d->P_q0_mm);
     free_array2_d(d->P_q0_pm);

     if (d->solver & XRTM_SOLVER_SOS || (d->solver & XRTM_SOLVER_TWO_OS && ! (d->options & XRTM_OPTION_SFI))) {
          free_array3_d(d->P_qq_pp);
          free_array3_d(d->P_qq_mp);

          if (d->options & XRTM_OPTION_VECTOR) {
               free_array3_d(d->P_qq_mm);
               free_array3_d(d->P_qq_pm);
          }
     }

     if (d->options & XRTM_OPTION_SFI && d->n_umus_v > 0) {
          free_array2_d(d->P_u0_mm);
          free_array2_d(d->P_u0_pm);

          free_array3_d(d->P_uq_pp);
          free_array3_d(d->P_uq_mp);

          if (d->options & XRTM_OPTION_VECTOR) {
               free_array3_d(d->P_uq_mm);
               free_array3_d(d->P_uq_pm);
          }
     }

     free_array3_d(d->r_p);
     free_array3_d(d->t_p);

     if (d->options & XRTM_OPTION_VECTOR) {
          free_array3_d(d->r_m);
          free_array3_d(d->t_m);
     }

     if (d->surface) {
          free_array1_d(d->Rs_q0);
          free_array2_d(d->Rs_qq);

          if (d->n_umus_v > 0) {
               free_array1_d(d->Rs_u0);
               free_array2_d(d->Rs_uq);
          }
     }
}



/*******************************************************************************
 *
 ******************************************************************************/
static int fourier_get_bvp_a(xrtm_data *d, int i_four, int solver, double **I_p, double **I_m, double **I_p_a, double **I_m_a, save_tree_data save_tree, work_data work) {

     int i;

     int surface;
     int thermal;

     double **P_q0_mm_a;
     double **P_q0_pm_a;

     double ***P_qq_pp_a;
     double ***P_qq_mp_a;
     double ***P_qq_mm_a;
     double ***P_qq_pm_a;

     double **P_u0_mm_a;
     double **P_u0_pm_a;

     double ***P_uq_pp_a;
     double ***P_uq_mp_a;
     double ***P_uq_mm_a;
     double ***P_uq_pm_a;

     double ***r_p_a;
     double ***t_p_a;
     double ***r_m_a;
     double ***t_m_a;

     double *Rs_q0_a;
     double *Rs_u0_a;
     double **Rs_qq_a;
     double **Rs_uq_a;

     forward_save_fourier_get_bvp_data *save;


     surface = solar_surface_for_i_four(d, i_four);

     thermal = d->options & XRTM_OPTION_SOURCE_THERMAL && i_four == 0;


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     save_tree_encode_s(&save_tree, "fourier_get_bvp");

     save_tree_retrieve_data(&save_tree, forward_save_fourier_get_bvp_data, &save);


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     P_q0_mm_a = get_work2(&work, WORK_DX, WORK_LAYERS_V, NULL);
     P_q0_pm_a = get_work2(&work, WORK_DX, WORK_LAYERS_V, NULL);

     if (solver & XRTM_SOLVER_SOS || (solver & XRTM_SOLVER_TWO_OS && ! (d->options & XRTM_OPTION_SFI))) {
          P_qq_pp_a = get_work2(&work, WORK_DXX, WORK_LAYERS_V, NULL);
          P_qq_mp_a = get_work2(&work, WORK_DXX, WORK_LAYERS_V, NULL);

          if (d->options & XRTM_OPTION_VECTOR) {
               P_qq_mm_a = get_work2(&work, WORK_DXX, WORK_LAYERS_V, NULL);
               P_qq_pm_a = get_work2(&work, WORK_DXX, WORK_LAYERS_V, NULL);
          }
     }

     if (d->options & XRTM_OPTION_SFI && d->n_umus > 0) {
          P_u0_mm_a = get_work2(&work, WORK_DU, WORK_LAYERS_V, NULL);
          P_u0_pm_a = get_work2(&work, WORK_DU, WORK_LAYERS_V, NULL);

          P_uq_pp_a = get_work2(&work, WORK_DUX, WORK_LAYERS_V, NULL);
          P_uq_mp_a = get_work2(&work, WORK_DUX, WORK_LAYERS_V, NULL);

          if (d->options & XRTM_OPTION_VECTOR) {
               P_uq_mm_a = get_work2(&work, WORK_DUX, WORK_LAYERS_V, NULL);
               P_uq_pm_a = get_work2(&work, WORK_DUX, WORK_LAYERS_V, NULL);
          }
     }

     r_p_a = get_work2(&work, WORK_DXX, WORK_LAYERS_V, NULL);
     t_p_a = get_work2(&work, WORK_DXX, WORK_LAYERS_V, NULL);

     if (d->options & XRTM_OPTION_VECTOR) {
          r_m_a = get_work2(&work, WORK_DXX, WORK_LAYERS_V, NULL);
          t_m_a = get_work2(&work, WORK_DXX, WORK_LAYERS_V, NULL);
     }

     if (surface) {
          Rs_q0_a = get_work1(&work, WORK_DX);
          Rs_qq_a = get_work1(&work, WORK_DXX);

          if (d->n_umus > 0) {
               Rs_u0_a = get_work1(&work, WORK_DU);
               Rs_uq_a = get_work1(&work, WORK_DUX);
          }
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     for (i = 0; i < d->n_layers; ++i) {
          init_array1_d(P_q0_mm_a[i], d->n_quad_v_x, 0.);
          init_array1_d(P_q0_pm_a[i], d->n_quad_v_x, 0.);

          if (solver & XRTM_SOLVER_SOS || (solver & XRTM_SOLVER_TWO_OS && ! (d->options & XRTM_OPTION_SFI))) {
               init_array2_d(P_qq_pp_a[i], d->n_quad_v_x, d->n_quad_v_x, 0.);
               init_array2_d(P_qq_mp_a[i], d->n_quad_v_x, d->n_quad_v_x, 0.);

               if (d->options & XRTM_OPTION_VECTOR) {
                    init_array2_d(P_qq_mm_a[i], d->n_quad_v_x, d->n_quad_v_x, 0.);
                    init_array2_d(P_qq_pm_a[i], d->n_quad_v_x, d->n_quad_v_x, 0.);
               }
          }

          if (d->options & XRTM_OPTION_SFI && d->n_umus_v > 0) {
               init_array1_d(P_u0_mm_a[i], d->n_umus_v, 0.);
               init_array1_d(P_u0_pm_a[i], d->n_umus_v, 0.);

               init_array2_d(P_uq_pp_a[i], d->n_umus_v, d->n_quad_v_x, 0.);
               init_array2_d(P_uq_mp_a[i], d->n_umus_v, d->n_quad_v_x, 0.);

               if (d->options & XRTM_OPTION_VECTOR) {
                    init_array2_d(P_uq_mm_a[i], d->n_umus_v, d->n_quad_v_x, 0.);
                    init_array2_d(P_uq_pm_a[i], d->n_umus_v, d->n_quad_v_x, 0.);
               }
          }

          init_array2_d(r_p_a[i], d->n_quad_v_x, d->n_quad_v_x, 0.);
          init_array2_d(t_p_a[i], d->n_quad_v_x, d->n_quad_v_x, 0.);

          if (d->options & XRTM_OPTION_VECTOR) {
               init_array2_d(r_m_a[i], d->n_quad_v_x, d->n_quad_v_x, 0.);
               init_array2_d(t_m_a[i], d->n_quad_v_x, d->n_quad_v_x, 0.);
          }
      }

      if (surface) {
           init_array1_d(Rs_q0_a, d->n_quad_v_x, 0.);
           init_array2_d(Rs_qq_a, d->n_quad_v_x, d->n_quad_v_x, 0.);

           if (d->n_umus_v > 0) {
                init_array1_d(Rs_u0_a, d->n_umus_v, 0.);
                init_array2_d(Rs_uq_a, d->n_umus_v, d->n_quad_v_x, 0.);
           }
      }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (solver & XRTM_SOLVER_EIG_BVP)
          rtm_eig_bvp_a(i_four, d->n_quad,   d->n_stokes, d->n_layers, d->qf, d->qx_v, d->qw_v, d->F_0, d->mu_0, d->ulevels, d->utaus2, d->n_ulevels, d->umus_v, d->n_umus, d->levels_b, d->omega, d->omega_a, d->ltau, d->ltau_a, save->Rs_qq, Rs_qq_a, save->Rs_u0, Rs_u0_a, save->Rs_uq, Rs_uq_a,                     d->btran, d->btran_a, d->as_0, d->as_0_a, d->atran, d->atran_a, save->P_q0_mm, save->P_q0_pm, save->P_u0_mm, save->P_u0_pm, save->P_uq_pp, save->P_uq_mp, save->P_uq_mm, save->P_uq_pm, save->r_p, save->t_p, save->r_m, save->t_m, P_q0_mm_a, P_q0_pm_a, P_u0_mm_a, P_u0_pm_a, P_uq_pp_a, P_uq_mp_a, P_uq_mm_a, P_uq_pm_a, r_p_a, t_p_a, r_m_a, t_m_a, d->I1_m[i_four], d->I1_m_a[i_four], d->In_p[i_four], d->In_p_a[i_four], I_p, I_m, I_p_a, I_m_a, d->options & XRTM_OPTION_SFI, surface, thermal, d->options & XRTM_OPTION_UPWELLING_OUTPUT, d->options & XRTM_OPTION_DOWNWELLING_OUTPUT, d->options & XRTM_OPTION_OUTPUT_AT_TAUS, d->options & XRTM_OPTION_VECTOR, d->derivs.layers_union, d->derivs.beam_union, save_tree, work);
/*
     else
     if (solver & XRTM_SOLVER_MEM_BVP)
          rtm_mem_bvp_a(i_four, d->n_quad,   d->n_stokes, d->n_layers, d->qf, d->qx_v, d->qw_v, d->F_0, d->mu_0, d->ulevels, d->utaus2, d->n_ulevels, d->umus_v, d->n_umus,              d->omega, d->omega_a, d->ltau, d->ltau_a, Rs_qq, Rs_qq_a, Rs_u0, Rs_u0_a, Rs_uq, Rs_uq_a, d->btau, d->btau_a, d->btran, d->btran_a, d->as_0, d->as_0_a, d->atran, d->atran_a, P_q0_mm, P_q0_pm, P_u0_mm, P_u0_pm, P_uq_pp, P_uq_mp, P_uq_mm, P_uq_pm, r_p, t_p, r_m, t_m, P_q0_mm_a, P_q0_pm_a, P_u0_mm_a, P_u0_pm_a, P_uq_pp_a, P_uq_mp_a, P_uq_mm_a, P_uq_pm_a, r_p_a, t_p_a, r_m_a, t_m_a, d->I1_m[i_four], I1_m_a, d->In_p[i_four], In_p_a, I_p, I_m, I_p_a, I_m_a, d->options & XRTM_OPTION_SFI, surface,                                                                                                  d->options & XRTM_OPTION_OUTPUT_AT_TAUS, d->options & XRTM_OPTION_VECTOR, d->misc_input.eigen_solver_gen_real, d->misc_input.eigen_solver_gen_complex, d->derivs.layers_union, d->derivs.beam_union, save_tree, work);
     else
     if (solver & XRTM_SOLVER_SOS)
          rtm_sos_a    (i_four, d->n_quad_x, d->n_stokes, d->n_layers, d->qf, d->qx_v, d->qw_v, d->F_0, d->mu_0, d->ulevels, d->utaus2, d->n_ulevels, d->omega, d->omega_a, d->ltau, d->ltau_a, Rs_qq,               Rs_qq_a, d->btran, d->btran_a, d->as_0, d->as_0_a, d->atran, d->atran_a, P_q0_mm, P_q0_pm, P_qq_pp, P_qq_mp, P_qq_mm, P_qq_pm, P_q0_mm_a, P_q0_pm_a, P_qq_pp_a, P_qq_mp_a, P_qq_mm_a, P_qq_pm_a, d->I1_m[i_four], I1_m_a, d->In_p[i_four], In_p_a, I_p, I_m, I_p_a, I_m_a, d->sos_max_os, d->sos_max_tau, d->sos_tol, d->options & XRTM_OPTION_SFI, surface, d->options & XRTM_OPTION_OUTPUT_AT_TAUS, d->options & XRTM_OPTION_VECTOR, d->derivs.layers_union, d->derivs.beam_union, work);
     else
     if (solver & XRTM_SOLVER_TWO_OS) {
          if (! (d->options & XRTM_OPTION_SFI))
               rtm_two_os_a(i_four, d->n_quad,   d->n_stokes, d->n_layers, d->qf, d->qx_v, d->qw_v, d->F_0, d->mu_0, d->ulevels, d->utaus2, d->n_ulevels, d->umus_v, d->n_umus, d->omega, d->omega_a, d->ltau, d->ltau_a, Rs_q0, Rs_q0_a, Rs_qq, Rs_qq_a, d->btran, d->btran_a, d->as_0, d->as_0_a, d->atran, d->atran_a, P_q0_mm, P_q0_pm, P_qq_pp, P_qq_mp, P_qq_mm, P_qq_pm, P_q0_mm_a, P_q0_pm_a, P_qq_pp_a, P_qq_mp_a, P_qq_mm_a, P_qq_pm_a,                         I_p, I_m, I_p_a, I_m_a,                                            d->options & XRTM_OPTION_SFI, surface, d->options & XRTM_OPTION_OUTPUT_AT_TAUS, d->options & XRTM_OPTION_VECTOR, 0, d->derivs.layers_union, d->derivs.beam_union, work);
          else
               rtm_two_os_a(i_four, d->n_quad_x, d->n_stokes, d->n_layers, d->qf, d->qx_v, d->qw_v, d->F_0, d->mu_0, d->ulevels, d->utaus2, d->n_ulevels, d->umus_v, d->n_umus, d->omega, d->omega_a, d->ltau, d->ltau_a, Rs_q0, Rs_q0_a, Rs_uq, Rs_uq_a, d->btran, d->btran_a, d->as_0, d->as_0_a, d->atran, d->atran_a, P_q0_mm, P_q0_pm, P_uq_pp, P_uq_mp, P_uq_mm, P_uq_pm, P_q0_mm_a, P_q0_pm_a, P_uq_pp_a, P_uq_mp_a, P_uq_mm_a, P_uq_pm_a,                         I_p, I_m, I_p_a, I_m_a,                                            d->options & XRTM_OPTION_SFI, surface, d->options & XRTM_OPTION_OUTPUT_AT_TAUS, d->options & XRTM_OPTION_VECTOR, 1, d->derivs.layers_union, d->derivs.beam_union, work);
     }
*/
#ifdef DEBUG
     else {
         fprintf(stderr, "ERROR: fourier_get_bvp(): end of if / else if\n");
         exit(1);
     }
#endif

     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
/*
     update_utaus(d, work);
*/
     if (surface) {
          get_brdf_vecs_q0_a(d, i_four, NULL, Rs_q0_a, work);
          get_brdf_mats_qq_a(d, i_four, NULL, Rs_qq_a, work);
     }

     if (surface && d->options & XRTM_OPTION_SFI && d->n_umus > 0) {
          get_brdf_vecs_u0_a(d, i_four, NULL, Rs_u0_a, work);
          get_brdf_mats_uq_a(d, i_four, NULL, Rs_uq_a, work);
     }

     update_diff_bound_input_a(d, i_four, save_tree, work);

     if (solver & XRTM_SOLVER_EIG_BVP || solver & XRTM_SOLVER_MEM_BVP) {
          if (get_local_r_t_u_w_all_a(d, i_four, NULL, NULL, NULL, NULL, r_p_a, t_p_a, r_m_a, t_m_a, 1, save_tree, &work)) {
               fprintf(stderr, "ERROR: get_local_r_t_u_w_all_a()\n");
               return -1;
          }
     }

     if (get_phase_vecs_q0_all_a(d, i_four, P_q0_mm_a, P_q0_pm_a, 1, &work)) {
          fprintf(stderr, "ERROR: get_phase_vecs_q0_all_a()\n");
          return -1;
     }

     if (solver & XRTM_SOLVER_SOS || (solver & XRTM_SOLVER_TWO_OS && ! (d->options & XRTM_OPTION_SFI))) {
          if (get_phase_mats_qq_all_a(d, i_four, P_qq_pp_a, P_qq_mp_a, P_qq_mm_a, P_qq_pm_a, 1, &work)) {
               fprintf(stderr, "ERROR: get_phase_mats_qq_all_a()\n");
               return -1;
          }
     }

     if (d->options & XRTM_OPTION_SFI && d->n_umus > 0) {
          if (get_phase_vecs_u0_all_a(d, i_four, P_u0_mm_a, P_u0_pm_a, 1, &work)) {
               fprintf(stderr, "ERROR: get_phase_vecs_u0_all_a()\n");
               return -1;
          }
     }

     if (d->options & XRTM_OPTION_SFI && d->n_umus > 0) {
          if (get_phase_mats_uq_all_a(d, i_four, P_uq_pp_a, P_uq_mp_a, P_uq_mm_a, P_uq_pm_a, 1, &work)) {
               fprintf(stderr, "ERROR: get_phase_mats_uq_all_a()\n");
               return -1;
          }
     }

     update_beam_params_all_a(d, work);

     update_opt_props_all_a(d, work);


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     init_array2_d(I_p_a, d->n_ulevels, d->n_quad_v_x, 0.);
     init_array2_d(I_m_a, d->n_ulevels, d->n_quad_v_x, 0.);


     return 0;
}



/*******************************************************************************
 *
 ******************************************************************************/
static int get_fourier_a(xrtm_data *d, int i_four, int solver, double **I_p, double **I_m, double **I_p_a, double **I_m_a, save_tree_data save_tree, work_data work) {

     int i;


     /*--------------------------------------------------------------------
      *
      *------------------------------------------------------------------*/
     if (d->misc_input.use_symmetric_form) {
          for (i = 0; i < d->n_ulevels; ++i)
               vec_sim_trans2(d->n_quad_v_d, I_p_a[i], I_m_a[i], d->alpha2);
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (solver & XRTM_SOLVERS_ADDING) {
          if (d->options & XRTM_OPTION_TOP_DOWN_ADDING) {
/*
               if (fourier_get_add_top_down_a(d, i_four, solver, I_p, I_m, I_p_a, I_m_a, save_tree, work)) {
                    fprintf(stderr, "ERROR: fourier_get_add_top_down_a()\n");
                    return -1;
               }
*/
          }
          else
          if (d->options & XRTM_OPTION_BOTTOM_UP_ADDING) {
               if (fourier_get_add_bottom_up_a(d, i_four, solver, I_p, I_m, I_p_a, I_m_a, save_tree, work)) {
                    fprintf(stderr, "ERROR: fourier_get_add_bottom_up_a()\n");
                    return -1;
               }
          }
          else {
               if (fourier_get_add_both_ways_a(d, i_four, solver, I_p, I_m, I_p_a, I_m_a, save_tree, work)) {
                    fprintf(stderr, "ERROR: fourier_get_add_both_ways_a()\n");
                    return -1;
               }
          }
     }
     else if (solver & XRTM_SOLVERS_BVP) {
          if (fourier_get_bvp_a(d, i_four, solver, I_p, I_m, I_p_a, I_m_a, save_tree, work)) {
               fprintf(stderr, "ERROR: fourier_get_bvp_a()\n");
               return -1;
          }
     }
#ifdef DEBUG
     else {
         fprintf(stderr, "ERROR: get_fourier(): end of if / else if\n");
         exit(1);
     }
#endif

     return 0;
}



/*******************************************************************************
 *
 ******************************************************************************/
static void forward_save_apply_corrections_radiance_free(forward_save_apply_corrections_radiance_data *d);

int forward_save_apply_corrections_radiance_alloc(forward_save_apply_corrections_radiance_data *d, int options, int n_coef, int n_quad_v_x, int n_mus2, int n_phis, int n_layers, int n_stokes) {

     d->free = (void (*)(void *)) forward_save_apply_corrections_radiance_free;

     d->options = options;

     if (options & XRTM_OPTION_PHASE_SCALAR) {
          d->polys_up = (void ***) alloc_array3_d(n_mus2, n_phis, n_coef);
          d->polys_dn = (void ***) alloc_array3_d(n_mus2, n_phis, n_coef);
     }
     else
     if (options & XRTM_OPTION_PHASE_MATRIX_GC) {
          d->polys_up = (void ***) alloc_array4_d(n_mus2, n_phis, 4, n_coef);
          d->polys_dn = (void ***) alloc_array4_d(n_mus2, n_phis, 4, n_coef);
     }
     else
     if (options & XRTM_OPTION_PHASE_MATRIX_LC) {
          d->polys_up = (void ***) alloc_array3_d(n_mus2, n_phis, n_coef);
          d->polys_dn = (void ***) alloc_array3_d(n_mus2, n_phis, n_coef);
     }

     d->P_trun_up = alloc_array4_d(n_mus2, n_phis, n_layers, n_stokes);
     d->P_full_up = alloc_array4_d(n_mus2, n_phis, n_layers, n_stokes);

     d->P_trun_dn = alloc_array4_d(n_mus2, n_phis, n_layers, n_stokes);
     d->P_full_dn = alloc_array4_d(n_mus2, n_phis, n_layers, n_stokes);

     return 0;
}



static void forward_save_apply_corrections_radiance_free(forward_save_apply_corrections_radiance_data *d) {

     if (d->options & XRTM_OPTION_PHASE_SCALAR) {
          free_array3_d((double ***) d->polys_up);
          free_array3_d((double ***) d->polys_dn);
     }
     else
     if (d->options & XRTM_OPTION_PHASE_MATRIX_GC) {
          free_array4_d((double ****) d->polys_up);
          free_array4_d((double ****) d->polys_dn);
     }
     else
     if (d->options & XRTM_OPTION_PHASE_MATRIX_LC) {
          free_array3_d((double ***) d->polys_up);
          free_array3_d((double ***) d->polys_dn);
     }

     free_array4_d(d->P_trun_up);
     free_array4_d(d->P_full_up);

     free_array4_d(d->P_trun_dn);
     free_array4_d(d->P_full_dn);
}



/*******************************************************************************
 *
 ******************************************************************************/
static int apply_corrections_radiance_a(xrtm_data *d, int n_phis, double **phis, double ****I_p, double ****I_m, double ****I_p_a, double ****I_m_a, save_tree_data save_tree, work_data work) {

     int i;
     int ii;
     int j;
     int k;

     int i_mus2_v;
     int n_mus2;

     double a;

     double f;
     double f_a;

     double *mus2;

     double *omega_tms;
     double *omega_tms_a;

     double **P_trun_a;
     double **P_full_a;

     double **I_a;

     forward_save_apply_corrections_radiance_data *save;


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (d->options & XRTM_OPTION_N_T_TMS) {

          if (d->n_umus == 0) {
               mus2     = d->qx;
               i_mus2_v = 0;
               n_mus2   = d->n_quad;
          }
          else {
               mus2     = d->umus;
               i_mus2_v = d->n_quad_v;
               n_mus2   = d->n_umus;
          }
/*
          else {
               mus2     = d->qx;
               i_mus2_v = 0;
               n_mus2   = d->n_quad + d->n_umus;
          }
*/

          /*--------------------------------------------------------------------
           *
           *------------------------------------------------------------------*/
          save_tree_encode_s(&save_tree, "apply_corrections_radiance");

          save_tree_retrieve_data(&save_tree, forward_save_apply_corrections_radiance_data, &save);


          /*--------------------------------------------------------------------
           *
           *------------------------------------------------------------------*/
          omega_tms   = get_work1(&work, WORK_DLAYERS);

          omega_tms_a = get_work1(&work, WORK_DLAYERS);

          P_trun_a = get_work2(&work, WORK_DS, WORK_LAYERS_V, NULL);
          P_full_a = get_work2(&work, WORK_DS, WORK_LAYERS_V, NULL);

          I_a      = get_work_d2(&work, d->n_ulevels, d->n_stokes);


          /*--------------------------------------------------------------------
           *
           *------------------------------------------------------------------*/
/*
          update_utaus(d, work);
*/

          a = 2. * d->n_coef2 + 1.;

          for (i = 0; i < d->n_layers; ++i) {
               get_delta_m_f(d, i, 0, 0, a, &f, NULL);

               n_t_tms_scaling(0, f, NULL, d->omega0[i], &omega_tms[i], NULL, NULL);
          }


          /*--------------------------------------------------------------------
           *
           *------------------------------------------------------------------*/
          init_array1_d(omega_tms_a, d->n_layers, 0.);


          /*--------------------------------------------------------------------
           *
           *------------------------------------------------------------------*/
          ii = i_mus2_v;

          for (i = 0; i < n_mus2; ++i) {
               for (j = 0; j < n_phis; ++j) {
                    save_tree_recode_i_j(&save_tree, i, j, i == 0 && j == 0);


                    if (d->options & XRTM_OPTION_UPWELLING_OUTPUT) {
                         for (k = 0; k < d->n_ulevels; ++k)
                              dvec_copy(I_a[k], I_p_a[k][i][j], d->n_stokes);

                         if (d->options & XRTM_OPTION_FOUR_CONV_OLD || d->fourier_tol == 0.) {
                              for (k = 0; k < d->n_layers; ++k) {
                                   init_array1_d(P_trun_a[k], d->n_stokes, 0.);
                                   init_array1_d(P_full_a[k], d->n_stokes, 0.);
                              }

                              n_t_tms_correction_up_a(d->n_stokes, d->n_layers, d->F_0, d->n_ulevels, d->ulevels, d->utaus2, &mus2[i], 1, d->omega, d->omega_a, omega_tms, omega_tms_a, d->ltau, d->ltau_a, d->btran, d->btran_a, d->as_0, d->as_0_a, d->atran, d->atran_a, save->P_trun_up[i][j], save->P_full_up[i][j], P_trun_a, P_full_a, d->In_p[0]+ii, d->In_p_a[0]+ii, NULL, I_a, d->options & XRTM_OPTION_OUTPUT_AT_TAUS, save_tree, work);

                              get_phase_func_trun_all_a(d, -mus2[i], phis[i][j], save->polys_up[i][j], P_trun_a, work);
                              get_phase_func_full_all_a(d, -mus2[i], phis[i][j], save->polys_up[i][j], P_full_a, work);
                         }
                         else {
                              for (k = 0; k < d->n_layers; ++k)
                                   init_array1_d(P_full_a[k], d->n_stokes, 0.);

                              single_scattered_radiance_up_a(d->n_stokes, d->n_layers, d->F_0, d->n_ulevels, d->ulevels, d->utaus2, &mus2[i], 1, omega_tms, omega_tms_a, d->ltau, d->ltau_a, d->btran, d->btran_a, d->as_0, d->as_0_a, d->atran, d->atran_a, save->P_full_up[i][j], P_full_a, d->In_p[0]+ii, d->In_p_a[0]+ii, NULL, I_a, d->options & XRTM_OPTION_OUTPUT_AT_TAUS, save_tree, work);

                              get_phase_func_full_all_a(d, -mus2[i], phis[i][j], save->polys_up[i][j], P_full_a, work);
                         }
                    }


                    if (d->options & XRTM_OPTION_DOWNWELLING_OUTPUT) {
                         for (k = 0; k < d->n_ulevels; ++k)
                              dvec_copy(I_a[k], I_m_a[k][i][j], d->n_stokes);

                         if (d->options & XRTM_OPTION_FOUR_CONV_OLD || d->fourier_tol == 0.) {
                              for (k = 0; k < d->n_layers; ++k) {
                                   init_array1_d(P_trun_a[k], d->n_stokes, 0.);
                                   init_array1_d(P_full_a[k], d->n_stokes, 0.);
                              }

                              n_t_tms_correction_dn_a(d->n_stokes, d->n_layers, d->F_0, d->n_ulevels, d->ulevels, d->utaus2, &mus2[i], 1, d->omega, d->omega_a, omega_tms, omega_tms_a, d->ltau, d->ltau_a, d->btran, d->btran_a, d->as_0, d->as_0_a, d->atran, d->atran_a, save->P_trun_dn[i][j], save->P_full_dn[i][j], P_trun_a, P_full_a, d->I1_m[0]+ii, d->I1_m_a[0]+ii, NULL, I_a, d->options & XRTM_OPTION_OUTPUT_AT_TAUS, save_tree, work);

                              get_phase_func_trun_all_a(d,  mus2[i], phis[i][j], save->polys_dn[i][j], P_trun_a, work);
                              get_phase_func_full_all_a(d,  mus2[i], phis[i][j], save->polys_dn[i][j], P_full_a, work);
                         }
                         else {
                              for (k = 0; k < d->n_layers; ++k)
                                   init_array1_d(P_full_a[k], d->n_stokes, 0.);

                              single_scattered_radiance_dn_a(d->n_stokes, d->n_layers, d->F_0, d->n_ulevels, d->ulevels, d->utaus2, &mus2[i], 1, omega_tms, omega_tms_a, d->ltau, d->ltau_a, d->btran, d->btran_a, d->as_0, d->as_0_a, d->atran, d->atran_a, save->P_full_dn[i][j], P_full_a, d->I1_m[0]+ii, d->I1_m_a[0]+ii, NULL, I_a, d->options & XRTM_OPTION_OUTPUT_AT_TAUS, save_tree, work);

                              get_phase_func_full_all_a(d,  mus2[i], phis[i][j], save->polys_dn[i][j], P_full_a, work);
                         }
                    }
               }

               ii += d->n_stokes;
          }

          save_tree_decode_i_j(&save_tree);


          /*--------------------------------------------------------------------
           *
           *------------------------------------------------------------------*/
          a = 2. * d->n_coef2 + 1.;

          f_a = 0.;

          for (i = 0; i < d->n_layers; ++i) {
               get_delta_m_f(d, i, 0, 0, a, &f, NULL);

               n_t_tms_scaling_a(f, &f_a, d->omega0[i], omega_tms[i], &d->omega0_a[i], omega_tms_a[i]);

               get_delta_m_f_a(d, i, a, f, &f_a);
          }

          update_diff_bound_input_a(d, 0, save_tree, work);

          update_beam_params_all_a(d, work);

          update_opt_props_all_a(d, work);
     }


     return 0;
}



/*******************************************************************************
 *
 ******************************************************************************/
static void forward_save_apply_corrections_fourier_free(forward_save_apply_corrections_fourier_data *d);

int forward_save_apply_corrections_fourier_alloc(forward_save_apply_corrections_fourier_data *d, int n_quad_v, int n_quad_v_x, int n_umus_v, int n_layers) {

     d->free = (void (*)(void *)) forward_save_apply_corrections_fourier_free;

     d->In_p      = alloc_array1_d(n_quad_v + n_umus_v);
     d->I1_m      = alloc_array1_d(n_quad_v + n_umus_v);

     d->P_q0_mm   = alloc_array2_d(n_layers, n_quad_v_x);
     d->P_q0_pm   = alloc_array2_d(n_layers, n_quad_v_x);

     d->P_u0_mm   = alloc_array2_d(n_layers, n_umus_v);
     d->P_u0_pm   = alloc_array2_d(n_layers, n_umus_v);

     return 0;
}



static void forward_save_apply_corrections_fourier_free(forward_save_apply_corrections_fourier_data *d) {

     free_array1_d(d->In_p);
     free_array1_d(d->I1_m);

     free_array2_d(d->P_q0_mm);
     free_array2_d(d->P_q0_pm);

     free_array2_d(d->P_u0_mm);
     free_array2_d(d->P_u0_pm);
}


/*******************************************************************************
 *
 ******************************************************************************/
static int apply_corrections_fourier_a(xrtm_data *d, int i_four, double **I_p, double **I_m, double **I_p_a, double **I_m_a, save_tree_data save_tree, work_data work) {

     int i;
     int j;

     double **I_a;

     double **P_q0_mm_a;
     double **P_q0_pm_a;

     double **P_u0_mm_a;
     double **P_u0_pm_a;

     forward_save_apply_corrections_fourier_data *save;


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (d->options & XRTM_OPTION_N_T_TMS && d->options & XRTM_OPTION_FOUR_CONV_NEW && d->fourier_tol > 0.) {
          save_tree_encode_s(&save_tree, "apply_corrections_fourier");

          save_tree_retrieve_data(&save_tree, forward_save_apply_corrections_fourier_data, &save);


          /*--------------------------------------------------------------------
           *
           *------------------------------------------------------------------*/
          I_a = get_work_d2(&work, d->n_ulevels, d->n_quad_v_x);

          P_q0_mm_a = get_work2(&work, WORK_DX, WORK_LAYERS_V, NULL);
          P_q0_pm_a = get_work2(&work, WORK_DX, WORK_LAYERS_V, NULL);

          if (d->options & XRTM_OPTION_SFI && d->n_umus > 0) {
              P_u0_mm_a = get_work2(&work, WORK_DU, WORK_LAYERS_V, NULL);
              P_u0_pm_a = get_work2(&work, WORK_DU, WORK_LAYERS_V, NULL);
          }


          /*--------------------------------------------------------------------
           *
           *------------------------------------------------------------------*/
          for (i = 0; i < d->n_layers; ++i) {
               dvec_zero(P_q0_mm_a[i], d->n_quad_v_x);
               dvec_zero(P_q0_pm_a[i], d->n_quad_v_x);
          }

          if (d->options & XRTM_OPTION_SFI && d->n_umus > 0) {
               for (i = 0; i < d->n_layers; ++i) {
                    dvec_zero(P_u0_mm_a[i], d->n_umus_v);
                    dvec_zero(P_u0_pm_a[i], d->n_umus_v);
               }
          }


          /*--------------------------------------------------------------------
           *
           *------------------------------------------------------------------*/

          if (d->options & XRTM_OPTION_UPWELLING_OUTPUT) {
               for (j = 0; j < d->n_ulevels; ++j)
                    dvec_scale(-1., I_p_a[j], I_a[j], d->n_quad_v_x);

               single_scattered_radiance_up_a(d->n_stokes, d->n_layers, d->F_0, d->n_ulevels, d->ulevels, d->utaus2, d->qx, d->n_quad_x, d->omega, d->omega_a, d->ltau, d->ltau_a, d->btran, d->btran_a, d->as_0, d->as_0_a, d->atran, d->atran_a, save->P_q0_pm, P_q0_pm_a, save->In_p, d->In_p_a[i_four], NULL, I_a, d->options & XRTM_OPTION_OUTPUT_AT_TAUS, save_tree, work);

               if (d->options & XRTM_OPTION_SFI && d->n_umus > 0) {

               }
          }

          if (d->options & XRTM_OPTION_DOWNWELLING_OUTPUT) {
               for (j = 0; j < d->n_ulevels; ++j)
                    dvec_scale(-1., I_m_a[j], I_a[j], d->n_quad_v_x);

               single_scattered_radiance_dn_a(d->n_stokes, d->n_layers, d->F_0, d->n_ulevels, d->ulevels, d->utaus2, d->qx, d->n_quad_x, d->omega, d->omega_a, d->ltau, d->ltau_a, d->btran, d->btran_a, d->as_0, d->as_0_a, d->atran, d->atran_a, save->P_q0_mm, P_q0_mm_a, save->I1_m, d->I1_m_a[i_four], NULL, I_a, d->options & XRTM_OPTION_OUTPUT_AT_TAUS, save_tree, work);

               if (d->options & XRTM_OPTION_SFI && d->n_umus > 0) {

               }
          }


          /*--------------------------------------------------------------------
           *
           *------------------------------------------------------------------*/
          update_diff_bound_input_a(d, i_four, save_tree, work);

          if (d->options & XRTM_OPTION_SFI && d->n_umus > 0) {
               if (get_phase_vecs_u0_all_a(d, i_four, P_u0_mm_a, P_u0_pm_a, 1, &work)) {
                    fprintf(stderr, "ERROR: get_phase_vecs_u0_all_a()\n");
                    return -1;
               }
          }

          if (get_phase_vecs_q0_all_a(d, i_four, P_q0_mm_a, P_q0_pm_a, 1, &work)) {
               fprintf(stderr, "ERROR: get_phase_vecs_q0_all_a()\n");
               return -1;
          }

          update_beam_params_all_a(d, work);

          update_opt_props_all_a(d, work);
     }


     return 0;
}



/*******************************************************************************
 *
 ******************************************************************************/
static void forward_save_get_solution_internal_free(forward_save_get_solution_internal_data *d);

int forward_save_get_solution_internal_alloc(forward_save_get_solution_internal_data *d, int n_four, int n_ulevels, int n_quad_v_x) {

     d->free = (void (*)(void *)) forward_save_get_solution_internal_free;

     d->i_p = alloc_array3_d(n_four, n_ulevels, n_quad_v_x);
     d->i_m = alloc_array3_d(n_four, n_ulevels, n_quad_v_x);

     return 0;
}



static void forward_save_get_solution_internal_free(forward_save_get_solution_internal_data *d) {

     free_array3_d(d->i_p);
     free_array3_d(d->i_m);
}



/*******************************************************************************
 *
 ******************************************************************************/
static void azimuth_modes_a(int n_quad, int n_stokes,
                            double a, double b, double *v1, double *v2) {

     int i;
     int ii;

     ii = 0;
     if (n_stokes == 1) {
          for (i = 0; i < n_quad; ++i) {
               v1[ii] += v2[ii] * a; ii++;
          }
     }
     else
     if (n_stokes == 2) {
          for (i = 0; i < n_quad; ++i) {
               v1[ii] += v2[ii] * a; ii++;
               v1[ii] += v2[ii] * a; ii++;
          }
     }
     else
     if (n_stokes == 3) {
          for (i = 0; i < n_quad; ++i) {
               v1[ii] += v2[ii] * a; ii++;
               v1[ii] += v2[ii] * a; ii++;
               v1[ii] += v2[ii] * b; ii++;
          }
     }
     else {
          for (i = 0; i < n_quad; ++i) {
               v1[ii] += v2[ii] * a; ii++;
               v1[ii] += v2[ii] * a; ii++;
               v1[ii] += v2[ii] * b; ii++;
               v1[ii] += v2[ii] * b; ii++;
          }
     }
}



/*******************************************************************************
 *
 ******************************************************************************/
/*
static
*/
int get_solution_internal_a(xrtm_data *d, enum xrtm_solver_mask solver, int solutions, int n_phis, double **phis, double ****I_p, double ****I_m, double ****I_p_a, double ****I_m_a, double *mean_p, double *mean_m, double *mean_p_a, double *mean_m_a, double *flux_p, double *flux_m, double *flux_p_a, double *flux_m_a, double *flux_div, double *flux_div_a, save_tree_data save_tree, work_data work) {

     int i;
     int j;
     int k;
     int kk;
     int l;
/*
     int n_quad_x;
*/
     int n_quad_v_x;

     int i_mus2;
     int n_mus2;

     double a;
     double b;
     double c;

     double *v1;
     double *v2;

     double **i_p_a;
     double **i_m_a;

     forward_save_get_solution_internal_data *save;


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     save_tree_encode_s(&save_tree, "get_solution_internal");

     save_tree_retrieve_data(&save_tree, forward_save_get_solution_internal_data, &save);


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (solver & XRTM_SOLVER_SINGLE) {
          get_single_a(d, n_phis, phis, I_p, I_m, I_p_a, I_m_a, save_tree, work);
          return 0;
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
/*
     n_quad_x   = d->n_quad   + d->n_umus;
*/
     n_quad_v_x = d->n_quad_v + d->n_umus_v;

     if (d->n_umus == 0) {
          i_mus2   = 0;
          n_mus2   = d->n_quad;
     }
     else {
          i_mus2   = d->n_quad;
          n_mus2   = d->n_umus;
     }
/*
     else {
          i_mus2   = 0;
          n_mus2   = d->n_quad + d->n_umus;
     }
*/


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     v1 = get_work1(&work, WORK_DD);
     v2 = get_work1(&work, WORK_DD);

     i_p_a = get_work_d2(&work, d->n_ulevels, n_quad_v_x);
     i_m_a = get_work_d2(&work, d->n_ulevels, n_quad_v_x);


     /*--------------------------------------------------------------------
      *
      *------------------------------------------------------------------*/
     init_array2_d(i_p_a, d->n_ulevels, n_quad_v_x, 0.);
     init_array2_d(i_m_a, d->n_ulevels, n_quad_v_x, 0.);


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     for (i = 0; i < save->count; ++i) {

          save_tree_recode_i(&save_tree, i, i == 0);


          if (solutions & XRTM_OUTPUT_RADIANCE) {
               for (j = 0; j < d->n_ulevels; ++j) {
                    kk = i_mus2 * d->n_stokes;
                    for (k = 0; k < n_mus2; ++k) {
                         for (l = 0; l < n_phis; ++l) {
                              a = i * (phis[k][l] - d->phi_0)*D2R;

                               b = cos(a);
                               if (d->n_stokes > 2)
                                    c = sin(a);

                              if (d->options & XRTM_OPTION_UPWELLING_OUTPUT) {
                                   dvec_copy(v1, I_p_a[j][k][l], d->n_stokes);
                                   azimuth_modes_a(1, d->n_stokes, b, c, &i_p_a[j][kk], v1);
                              }

                              if (d->options & XRTM_OPTION_DOWNWELLING_OUTPUT) {
                                   dvec_copy(v2, I_m_a[j][k][l], d->n_stokes);
                                   azimuth_modes_a(1, d->n_stokes, b, c, &i_m_a[j][kk], v2);
                              }
                         }

                         kk += d->n_stokes;
                    }
               }
          }


          /*--------------------------------------------------------------------
           *
           *------------------------------------------------------------------*/
/*
          if (solutions & XRTM_OUTPUT_RADIANCE_MEAN && i == 0) {
               for (j = 0; j < d->n_ulevels; ++j) {
                    jj = j == 0 ? 0 : d->n_layers;

                    mean_p[j] = radiance_to_mean(d->n_quad_v, d->qx_v, d->qw_v, d->F_0, d->mu_0, d->btran[jj], i_p[j]),
                    mean_m[j] = radiance_to_mean(d->n_quad_v, d->qx_v, d->qw_v, d->F_0, d->mu_0, d->btran[jj], i_m[j]) + d->F_0 * d->btran[jj] / (4. * PI);
               }
          }

          if (solutions & XRTM_OUTPUT_FLUX          && i == 0) {
               for (j = 0; j < d->n_ulevels; ++j) {
                    jj = j == 0 ? 0 : d->n_layers;

                    flux_p[j] = radiance_to_flux(d->n_quad_v, d->qx_v, d->qw_v, d->F_0, d->mu_0, d->btran[jj], i_p[j]),
                    flux_m[j] = radiance_to_flux(d->n_quad_v, d->qx_v, d->qw_v, d->F_0, d->mu_0, d->btran[jj], i_m[j]) + d->F_0 * d->mu_0 * d->btran[jj];
               }
          }
*/

          /*--------------------------------------------------------------------
           *
           *------------------------------------------------------------------*/
          if (apply_corrections_fourier_a(d, i, save->i_p[i], save->i_m[i], i_p_a, i_m_a, save_tree, work)) {
               fprintf(stderr, "ERROR: apply_corrections_fourier_a()\n");
               return -1;
          }


          /*--------------------------------------------------------------------
           *
           *------------------------------------------------------------------*/
          if (get_fourier_a(d, i, solver, save->i_p[i], save->i_m[i], i_p_a, i_m_a, save_tree, work)) {
               fprintf(stderr, "ERROR: get_fourier()\n");
               return -1;
          }
     }

     save_tree_decode_i(&save_tree);


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (apply_corrections_radiance_a(d, n_phis, phis, I_p, I_m, I_p_a, I_m_a, save_tree, work)) {
          fprintf(stderr, "ERROR: apply_corrections_radiance_a()\n");
          return -1;
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     init_array4_d(I_p_a, d->n_ulevels, n_mus2, n_phis, d->n_stokes, 0.);
     init_array4_d(I_m_a, d->n_ulevels, n_mus2, n_phis, d->n_stokes, 0.);


     return 0;
}



/*******************************************************************************
 *
 ******************************************************************************/
static int get_solution_a(xrtm_data *d, enum xrtm_solver_mask solver, int solutions, int n_phis, double **phis, double ****I_p, double ****I_m, double ****I_p_a, double ****I_m_a, double *mean_p, double *mean_m, double *mean_p_a, double *mean_m_a, double *flux_p, double *flux_m, double *flux_p_a, double *flux_m_a, double *flux_div, double *flux_div_a, save_tree_data save_tree, work_data work) {

     if (solver & XRTM_SOLVERS_INTERNAL) {
          if (get_solution_internal_a(d, solver, solutions, n_phis, phis, I_p, I_m, I_p_a, I_m_a, mean_p, mean_m, mean_p_a, mean_m_a, flux_p, flux_m, flux_p_a, flux_m_a, flux_div, flux_div_a, save_tree, work)) {
               fprintf(stderr, "ERROR: get_solution_internal_a()\n");
               return -1;
          }
     }
#ifdef INCLUDE_DEV_SOURCE
     else {

     }
#endif
     return 0;
}



/*******************************************************************************
 *
 ******************************************************************************/
int xrtm_solution_a(xrtm_data *d, enum xrtm_solver_mask solver, int solutions, int n_phis, double **phis, double ****I_p, double ****I_m, double ****I_p_a, double ****I_m_a, double *mean_p, double *mean_m, double *mean_p_a, double *mean_m_a, double *flux_p, double *flux_m, double *flux_p_a, double *flux_m_a, double *flux_div, double *flux_div_a) {

     if (! (solver & XRTM_SOLVERS_ALL)) {
          fprintf(stderr, "ERROR: invalid value for solver\n");
          return -1;
     }

     if (! (solver & d->solvers)) {
          fprintf(stderr, "ERROR: model not initialized for solver %s\n",
                  xrtm_solver_mask_to_name(solver));
          return -1;
     }

     if (! (solutions & XRTM_SOLUTION_ALL)) {
          fprintf(stderr, "ERROR: invalid value for solution\n");
          return -1;
     }

     if (n_phis <= 0) {
          fprintf(stderr, "ERROR: invalid value for n_phis: %d, must be > zero\n", n_phis);
          return -1;
     }
#ifdef INCLUDE_DEV_SOURCE
/*
     if (check_dev_solvers_solution(d, n_phis)) {
          fprintf(stderr, "ERROR: check_dev_solvers_solution()\n");
          return -1;
     }
*/
#endif
/*
     SOLUTION_PREPROCESS();
*/
     if (get_solution_a(d, solver, solutions, n_phis, phis, I_p, I_m, I_p_a, I_m_a, mean_p, mean_m, mean_p_a, mean_m_a, flux_p, flux_m, flux_p_a, flux_m_a, flux_div, flux_div_a, d->save_tree, d->work)) {
          fprintf(stderr, "ERROR: get_solution_a()\n");
          return -1;
     }

     return 0;
}



/*******************************************************************************
 *
 ******************************************************************************/
#ifdef INCLUDE_DEV_SOURCE
#include "../src2/xrtm_model_a2.c"
#endif
