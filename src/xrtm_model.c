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
#include "xrtm_brdf.h"
#include "xrtm_derivs.h"
#include "xrtm_doub_rts.h"
#include "xrtm_eig_bvp.h"
#include "xrtm_eig_rts.h"
#include "xrtm_external.h"
#include "xrtm_four_stream.h"
#include "xrtm_interface.h"
#include "xrtm_mem_bvp.h"
#include "xrtm_model.h"
#include "xrtm_model_a.h"
#include "xrtm_pade_rts.h"
#include "xrtm_radiance.h"
#include "xrtm_save_tree.h"
#include "xrtm_single.h"
#include "xrtm_scatter.h"
#include "xrtm_scatter_a.h"
#include "xrtm_six_stream.h"
#include "xrtm_sos.h"
#include "xrtm_source.h"
#include "xrtm_stacks.h"
#include "xrtm_support.h"
#include "xrtm_two_os.h"
#include "xrtm_two_stream.h"
#include "xrtm_utility.h"

/*
#define USE_NEW_SFI_FOUR_CONV
*/

/*******************************************************************************
 *
 ******************************************************************************/
static int resolve_deps0(xrtm_data *d,
                         int i, int *dep_flags, int deps) {

     if (dep_flags[i] &  deps) {
         dep_flags[i] &= (deps ^ DEP_MASK_MASK);
         return 1;
     }

     return 0;
}



static int resolve_deps1(xrtm_data *d,
                         int i, int j, int **dep_flags, int deps) {

     if (dep_flags[i][j] &  deps) {
         dep_flags[i][j] &= (deps ^ DEP_MASK_MASK);
         return 1;
     }

     return 0;
}



static int resolve_deps2(xrtm_data *d,
                         int i_four, int i_layer, int **dep_flags, int deps) {

     int i;

     if (dep_flags[i_four][i_layer] &  deps) {
         for (i = 0; i < d->n_four; ++i)
              dep_flags[i][i_layer] &= (deps ^ DEP_MASK_MASK);
         return 1;
     }

     return 0;
}



/*
static int resolve_deps3(xrtm_data *d,
                         int i_four, int i_layer, int **dep_flags, int deps) {

     int i;

     if (dep_flags[i_four][i_layer] &  deps) {
         for (i = 0; i < d->n_layers; ++i)
              dep_flags[i_four][i] &= (deps ^ DEP_MASK_MASK);
         return 1;
     }

     return 0;
}
*/


/*******************************************************************************
 *
 ******************************************************************************/
int xrtm_phase_type_to_coef_type(xrtm_data *d) {

     if (d->options & XRTM_OPTION_PHASE_SCALAR)
          return SCAT_COEF_GC;
     else
     if (d->options & XRTM_OPTION_PHASE_MATRIX_GC)
          return SCAT_COEF_GC;
     else
#ifdef DEBUG
     if (d->options & XRTM_OPTION_PHASE_MATRIX_LC)
#endif
          return SCAT_COEF_LC;
#ifdef DEBUG
     else {
         fprintf(stderr, "ERROR: xrtm_phase_type_to_coef_type(): end of if / else if\n");
         exit(1);
     }
#endif
}



/*******************************************************************************
 *
 ******************************************************************************/
int solar_surface_for_i_four(xrtm_data *d, int i_four) {

     return d->n_kernels > 0 && (i_four == 0 || brdf_needs_fourier(d->kernels, d->n_kernels));
}



int thermal_surface_for_i_four(xrtm_data *d, int i_four) {

     return d->n_kernels > 0 &&  i_four == 0;
}



/*******************************************************************************
 *
 ******************************************************************************/
/*
static
*/
void get_delta_m_f(xrtm_data *d, int i_layer, int i_deriv, int n_derivs, double a, double *f, double *f_l) {

     int i;
     int ii;

     if (d->solvers & XRTM_SOLVERS_USE_G) {
          *f = d->g0[i_layer] * d->g0[i_layer];

           for (i = i_deriv, ii = 0; i < n_derivs; ++i, ++ii)
               f_l[ii] = 2. * d->g0[i_layer] * d->g0_l[i][i_layer];
     }
     else {
          if (d->n_coef_layer[i_layer] <= d->n_coef2) {
               *f = 0.;

               for (i = i_deriv, ii = 0; i < n_derivs; ++i, ++ii)
                    f_l[ii] = 0.;
          }
          else {
               *f = d->coef0[i_layer][0][d->n_coef2] / a;

               for (i = i_deriv, ii = 0; i < n_derivs; ++i, ++ii)
                    f_l[ii] = d->coef0_l[i_layer][i][0][d->n_coef2] / a;
          }
     }
}



/*******************************************************************************
 *
 ******************************************************************************/
static int utaus_deps(xrtm_data *d) {

     return DEP_MASK_UTAUS | DEP_MASK_LTAU | DEP_MASK_COEF;
}



static void update_utaus(xrtm_data *d, work_data work) {

     int i;
     int ii;

     double a;

     double f;

     if (! (d->options & XRTM_OPTION_OUTPUT_AT_TAUS))
          return;

     if (d->dep_flags_utaus & DEP_MASK_LTAU) {
          d->ttau[0] = 0;
          for (i  = 0; i < d->n_layers; ++i)
               d->ttau[i+1] = d->ttau[i] + d->ltau0[i];
     }

     if (d->dep_flags_utaus & (DEP_MASK_UTAUS | DEP_MASK_LTAU)) {
          ii = 0;
          for (i = 0; i < d->n_ulevels - 1; ++i) {
               for ( ; d->utaus[i] >= d->ttau[ii]; ++ii) ;

               d->ulevels[i] = ii - 1;
               d->utaus20[i] = d->utaus[i] - d->ttau[ii - 1];
          }

          if (d->utaus[d->n_ulevels - 1] != d->ttau[d->n_layers])
               for ( ; d->utaus[i] >= d->ttau[ii]; ++ii) ;

          d->ulevels[i] = ii - 1;
          d->utaus20[i] = d->utaus[i] - d->ttau[ii - 1];
     }

     if (d->options & XRTM_OPTION_DELTA_M &&
         d->dep_flags_utaus & (DEP_MASK_UTAUS | DEP_MASK_LTAU | DEP_MASK_COEF)) {
          a = 2. * d->n_coef2 + 1.;

          for (i = 0; i < d->n_ulevels; ++i) {
               ii = d->ulevels[i];

               get_delta_m_f(d, ii, 0, 0, a, &f, NULL);

               delta_m_ltau(d->n_derivs, f, NULL, d->omega0[ii], NULL, d->utaus20[i], &d->utaus2[i], NULL, NULL, 0);
          }
     }


     d->dep_flags_utaus &= (utaus_deps(d) ^ DEP_MASK_MASK);
}



/*******************************************************************************
 *
 ******************************************************************************/
static int opt_props_deps(xrtm_data *d) {

     if (d->solvers & XRTM_SOLVERS_USE_G)
          return DEP_MASK_G | DEP_MASK_G_L | DEP_MASK_OMEGA | DEP_MASK_OMEGA_L |
                 DEP_MASK_LTAU | DEP_MASK_LTAU_L | utaus_deps(d);
     else
          return DEP_MASK_COEF | DEP_MASK_COEF_L | DEP_MASK_OMEGA | DEP_MASK_OMEGA_L |
                 DEP_MASK_LTAU | DEP_MASK_LTAU_L | utaus_deps(d);
}



static void update_opt_props(xrtm_data *d, int i_layer, work_data work) {

     int i;

     int flag;

     int deps;

     int coef_type;

     double a;

     double f;
     double *f_l;

     if (! (d->options & XRTM_OPTION_DELTA_M))
          return;


     deps = opt_props_deps(d);

     if (! (d->dep_flags_opt_props[i_layer] & deps))
          return;


     i = i_layer;

     a = 2. * d->n_coef2 + 1.;

     if (! (d->options & XRTM_OPTION_CALC_DERIVS))
          get_delta_m_f(d, i, 0, 0, a, &f, NULL);
     else {
          f_l = get_work1(&work, WORK_DDERIVS);

          get_delta_m_f(d, i, 0, d->n_derivs, a, &f, f_l);
     }

     coef_type = xrtm_phase_type_to_coef_type(d);

     flag = d->n_stokes != 1 ? 1 : 0;

     if (! (d->options & XRTM_OPTION_CALC_DERIVS)) {
          if (d->solvers & XRTM_SOLVERS_USE_G &&
              d->dep_flags_opt_props[i_layer] & DEP_MASK_G)
               delta_m_g(d->n_derivs, f, NULL, d->g0[i], &d->g[i], NULL, NULL, 0);
          else
          if (d->solvers & XRTM_SOLVERS_USE_COEF &&
              d->dep_flags_opt_props[i_layer] & DEP_MASK_COEF)
               delta_m_coef(d->n_coef2, d->n_derivs, f, NULL, d->coef0[i], d->coef[i], NULL, NULL, flag, coef_type, 0);

          if (d->dep_flags_opt_props[i_layer] & DEP_MASK_OMEGA)
               delta_m_omega(d->n_derivs, f, NULL, d->omega0[i], &d->omega[i], NULL, NULL, 0);

          if (d->dep_flags_opt_props[i_layer] & DEP_MASK_LTAU)
               delta_m_ltau(d->n_derivs, f, NULL, d->omega0[i], NULL, d->ltau0[i], &d->ltau[i], NULL, NULL, 0);
     }
     else {
          if (d->solvers & XRTM_SOLVERS_USE_G &&
              d->dep_flags_opt_props[i_layer] & DEP_MASK_G)
               delta_m_g(d->n_derivs, f, f_l, d->g0[i], &d->g[i], d->g0_l[i], d->g_l[i], 1);
          else
          if (d->solvers & XRTM_SOLVERS_USE_G &&
              d->dep_flags_opt_props[i_layer] & DEP_MASK_G_L)
               delta_m_g_l(d->n_derivs, f, f_l, d->g0[i], d->g0_l[i], d->g_l[i]);
          else
          if (d->solvers & XRTM_SOLVERS_USE_COEF &&
              d->dep_flags_opt_props[i_layer] & DEP_MASK_COEF)
               delta_m_coef(d->n_coef2, d->n_derivs, f, f_l, d->coef0[i], d->coef[i], d->coef0_l[i], d->coef_l[i], flag, coef_type, 1);
          else
          if (d->solvers & XRTM_SOLVERS_USE_COEF &&
              d->dep_flags_opt_props[i_layer] & DEP_MASK_COEF_L)
               delta_m_coef_l(d->n_coef2, d->n_derivs, f, f_l, d->coef0[i], d->coef0_l[i], d->coef_l[i], flag, coef_type);

          if (d->dep_flags_opt_props[i_layer] & DEP_MASK_OMEGA)
               delta_m_omega(d->n_derivs, f, f_l, d->omega0[i], &d->omega[i], d->omega0_l[i], d->omega_l[i], 1);
          else
          if (d->dep_flags_opt_props[i_layer] & DEP_MASK_OMEGA_L)
               delta_m_omega_l(d->n_derivs, f, f_l, d->omega0[i], d->omega0_l[i], d->omega_l[i]);

          if (d->dep_flags_opt_props[i_layer] & DEP_MASK_LTAU)
               delta_m_ltau(d->n_derivs, f, f_l, d->omega0[i], d->omega0_l[i], d->ltau0[i], &d->ltau[i], d->ltau0_l[i], d->ltau_l[i], 1);
          else
          if (d->dep_flags_opt_props[i_layer] & DEP_MASK_LTAU_L)
               delta_m_ltau_l(d->n_derivs, f, f_l, d->omega0[i], d->omega0_l[i], d->ltau0[i], d->ltau0_l[i], d->ltau_l[i]);
     }

     d->dep_flags_opt_props[i_layer] &= (deps ^ DEP_MASK_MASK);
#ifdef USE_AD_FOR_TL_CALC_UPDATE_OPT_PROPS
          update_opt_props_tl_with_ad(d, i_layer, work);
#endif
}



static void update_opt_props_all(xrtm_data *d, work_data work) {

     int i;

     for (i = 0; i < d->n_layers; ++i)
          update_opt_props(d, i, work);
}



/*******************************************************************************
 *
 ******************************************************************************/
static int chapman_deps(xrtm_data *d) {

     return DEP_MASK_MU_0 | DEP_MASK_PLANET_R | DEP_MASK_LEVELS_Z;
}


/*
static
*/
void get_chapman(xrtm_data *d, double ***chapman) {

     int deps;

     deps = chapman_deps(d);

     if (d->dep_flags_chapman & deps) {
          d->dep_flags_chapman &= (deps ^ DEP_MASK_MASK);

          chapman_functions(d->n_layers, d->mu_0, d->planet_r, d->levels_z, d->chapman);
     }

     *chapman = d->chapman;
}



/*******************************************************************************
 *
 ******************************************************************************/
static void beam_params_update(xrtm_data *d, int i_level, int *dep_flags, int deps, double *ltau, double **ltau_l, double *btau, double **btau_l, double *btran, double **btran_l, double *as_0, double **as_0_l, double *atran, double **atran_l, work_data work) {

     int i;
     int j;
     int k;

     int i_level2;

     double **chapman;
if (d->options & XRTM_OPTION_SOURCE_SOLAR) {
     if (i_level == d->n_layers)
          i_level2 = i_level;
     else
          i_level2 = i_level + 1;

     for (i = 1; i <= i_level2; ++i) {
          if (dep_flags[i-1] & deps)
               break;
     }

     if (d->options & XRTM_OPTION_PSA)
          get_chapman(d, &chapman);

     for (     ; i <= i_level2; ++i) {
          dep_flags[i-1] &= (deps ^ DEP_MASK_MASK);

          update_opt_props(d, i-1, work);

          if (! (d->options & XRTM_OPTION_PSA))
               btau[i] = btau[i-1] + ltau[i-1] / d->mu_0;
          else {
               btau[i] = 0;
               for (j = 0; j < i; ++j)
                    btau[i] += ltau[j] * d->chapman[i-1][j];
          }

          btran[i] = exp(-btau[i]);

          as_0[i-1] = (btau[i] - btau[i-1]) / ltau[i-1];

          atran[i-1] = exp(-ltau[i-1] * as_0[i-1]);

          if (d->options & XRTM_OPTION_CALC_DERIVS) {
               for (j = 0; j < d->n_derivs; ++j) {
                    if (! (d->options & XRTM_OPTION_PSA))
                         btau_l[i][j] = btau_l[i-1][j] +
                                        ltau_l[i-1][j] / d->mu_0;
                    else {
                         btau_l[i][j] = 0;
                         for (k = 0; k < i; ++k)
                               btau_l[i][j] +=
                               ltau_l[k][j] * d->chapman[i-1][k];
                    }

                    btran_l[i][j] = -btau_l[i][j] * btran[i];

                    as_0_l[i-1][j] = (btau_l[i][j] - btau_l[i-1][j] -
                                      as_0[i-1] * ltau_l[i-1][j]) / ltau[i-1];

                    atran_l[i-1][j] = -(ltau_l[i-1][j] * as_0[i-1] +
                                        ltau[i-1] * as_0_l[i-1][j]) * atran[i-1];
               }
          }
     }
}
}



/*******************************************************************************
 *
 ******************************************************************************/
static int beam_params_deps(xrtm_data *d) {

     return opt_props_deps(d) | chapman_deps(d) |
            DEP_MASK_MU_0 | DEP_MASK_LTAU | DEP_MASK_LTAU_L;
}



static void update_beam_params(xrtm_data *d, int i_level, work_data work) {

     beam_params_update(d, i_level, d->dep_flags_beam_params, beam_params_deps(d), d->ltau, d->ltau_l, d->btau, d->btau_l, d->btran, d->btran_l, d->as_0, d->as_0_l, d->atran, d->atran_l, work);
}



static void update_beam_params_all(xrtm_data *d, work_data work) {

     beam_params_update(d, d->n_layers, d->dep_flags_beam_params, beam_params_deps(d), d->ltau, d->ltau_l, d->btau, d->btau_l, d->btran, d->btran_l, d->as_0, d->as_0_l, d->atran, d->atran_l, work);
}



/*******************************************************************************
 *
 ******************************************************************************/
/*
static int beam_params0_deps(xrtm_data *d) {

     return chapman_deps(d) |
            DEP_MASK_MU_0 | DEP_MASK_LTAU | DEP_MASK_LTAU_L;
}



static void update_beam_params0(xrtm_data *d, int i_level, work_data work) {

     if (! (d->options & XRTM_OPTION_DELTA_M))
          update_beam_params(d, i_level, btau, btau_l, btran, btran_l, as_0, as_0_l, atran, atran_l, work);

     beam_params_update(d, i_level, d->dep_flags_beam_params0, beam_params0_deps(d), d->ltau0, d->ltau0_l, d->btau0, d->btau0_l, d->btran0, d->btran0_l, d->as_00, d->as_00_l, d->atran0, d->atran0_l, work);
}



static void update_beam_params0_all(xrtm_data *d, work_data work) {

     if (! (d->options & XRTM_OPTION_DELTA_M))
          update_beam_params_all(d, btau, btau_l, btran, btran_l, as_0, as_0_l, atran, atran_l, work);

     beam_params_update(d, d->n_layers, d->dep_flags_beam_params0, beam_params0_deps(d), d->ltau0, d->ltau0_l, d->btau0, d->btau0_l, d->btran0, d->btran0_l, d->as_00, d->as_00_l, d->atran0, d->atran0_l, work);
}
*/


/*******************************************************************************
 *
 ******************************************************************************/
static void Y_x_get(xrtm_data *d, int i_four, int *dep_flags, int deps, int n_mus, double *mus, int *order_x, double *c_x, double *d_x, void **Y_x0, void **Y_x, int flag) {

     int i;

     int flag2;

     flag2 = resolve_deps0(d, i_four, dep_flags, deps);

     if (! (d->options & XRTM_OPTION_SAVE_LEG_POLYS)) {
          *Y_x = Y_x0[0];

          flag2 = flag2 || i_four != *order_x;
     }
     else
          *Y_x = Y_x0[i_four];

     if (! flag2)
          return;

     if (i_four == *order_x) {
          if (! flag)
               leg_poly_assoc (i_four, d->n_coef2, *mus, c_x, d_x, (double *) *Y_x, 0);
          else
               leg_poly_assoc2(i_four, n_mus, d->n_coef2, mus, c_x, d_x, (double **) *Y_x, 0);
     }
     else {
          if (i_four < *order_x)
               *order_x = -1;

          for (i = *order_x + 1; i <= i_four; ++i) {
               if (! flag)
                    leg_poly_assoc (i, d->n_coef2, *mus, c_x, d_x, (double *) *Y_x, 1);
               else
                    leg_poly_assoc2(i, n_mus, d->n_coef2, mus, c_x, d_x, (double **) *Y_x, 1);
          }

          *order_x = i_four;
     }
}



/*******************************************************************************
 *
 ******************************************************************************/
static int Y_p_deps(xrtm_data *d) {

     int deps = DEP_MASK_INIT;

     if (! (d->options & XRTM_OPTION_SFI))
          deps |= DEP_MASK_UMUS;

     return deps;
}


/*
static
*/
void get_Y_p(xrtm_data *d, int i_four, double ***Y_p) {

     Y_x_get(d, i_four, d->dep_flags_Y_p, Y_p_deps(d), d->n_quad_x, d->qx, &d->order_p, &d->c_p, &d->d_p, (void **) d->Y_p, (void **) Y_p, 1);
}



/*******************************************************************************
 *
 ******************************************************************************/
static int Y_0_deps(xrtm_data *d) {

     int deps = DEP_MASK_INIT | DEP_MASK_MU_0;

     if (! (d->options & XRTM_OPTION_SFI))
          deps |= DEP_MASK_UMUS;

     return deps;
}


/*
static
*/
void get_Y_0(xrtm_data *d, int i_four, double **Y_0) {

     double a = -d->mu_0;

     Y_x_get(d, i_four, d->dep_flags_Y_0, Y_0_deps(d), 1, &a, &d->order_0, &d->c_0, &d->d_0, (void **) d->Y_0, (void **) Y_0, 0);
}



/*******************************************************************************
 *
 ******************************************************************************/
static int Y_u_deps(xrtm_data *d) {

     return DEP_MASK_INIT | DEP_MASK_UMUS;
}


/*
static
*/
void get_Y_u(xrtm_data *d, int i_four, double ***Y_u) {

     Y_x_get(d, i_four, d->dep_flags_Y_u, Y_u_deps(d), d->n_umus, d->umus, &d->order_u, &d->c_u, &d->d_u, (void **) d->Y_u, (void **) Y_u, 1);
}



/*******************************************************************************
 *
 ******************************************************************************/
static void phase_polys_get(xrtm_data *d, int n_coef, double mu1, double phi1, double mu2, double phi2, void *polys) {

     double c;

     double **polys2;

     c = scat_angle(mu1, phi1, mu2, phi2);

     if (d->options & XRTM_OPTION_PHASE_SCALAR)
          leg_poly(n_coef, c, (double *) polys);
     else
     if (d->options & XRTM_OPTION_PHASE_MATRIX_GC) {
          polys2 = (double **) polys;
          gen_spher_funcs_prt(n_coef, c, polys2[0], polys2[1], polys2[2], polys2[3]);
     }
     else
     if (d->options & XRTM_OPTION_PHASE_MATRIX_LC)
          leg_poly(n_coef, c, (double *) polys);
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
static void phase_func_get(xrtm_data *d, int n_coef, int i_layer, double mu, double phi, void *polys, double ***coef, double ****coef_l, double *P, double **P_l, work_data work) {

     int i;

     double c;

     if (d->options & XRTM_OPTION_PHASE_SCALAR) {
          phase_func(n_coef, (double *) polys, coef[i_layer][0], P);

          if (d->options & XRTM_OPTION_CALC_DERIVS && flags_or(d->derivs.layers[i_layer], d->n_derivs)) {
               for (i = 0; i < d->n_derivs; ++i) {
                    if (! d->derivs.layers[i_layer][i])
                         continue;

                    phase_func(n_coef, (double *) polys, coef_l[i_layer][i][0], P_l[i]);
               }
          }
     }
     else
     if (d->options & XRTM_OPTION_PHASE_MATRIX_GC) {
          c = scat_angle(d->mu_0, d->phi_0, mu, phi);

          build_scat_vector_gc(n_coef, d->n_stokes, (double **) polys, coef[i_layer], P);
          scat_vector_rotate_hovenier(d->n_stokes, d->mu_0, mu, c, phi - d->phi_0, P, P);

          if (d->options & XRTM_OPTION_CALC_DERIVS && flags_or(d->derivs.layers[i_layer], d->n_derivs)) {
               for (i = 0; i < d->n_derivs; ++i) {
                    if (! d->derivs.layers[i_layer][i])
                         continue;

                    build_scat_vector_gc(n_coef, d->n_stokes, (double **) polys, coef_l[i_layer][i], P_l[i]);
                    scat_vector_rotate_hovenier(d->n_stokes, d->mu_0, mu, c, phi - d->phi_0, P_l[i], P_l[i]);
               }
          }
     }
     else
     if (d->options & XRTM_OPTION_PHASE_MATRIX_LC) {
          c = scat_angle(d->mu_0, d->phi_0, mu, phi);

          build_scat_vector_lc(n_coef, d->n_stokes, (double *) polys, coef[i_layer], P);
          scat_vector_rotate_hovenier(d->n_stokes, d->mu_0, mu, c, phi - d->phi_0, P, P);

          if (d->options & XRTM_OPTION_CALC_DERIVS && flags_or(d->derivs.layers[i_layer], d->n_derivs)) {
               for (i = 0; i < d->n_derivs; ++i) {
                    if (! d->derivs.layers[i_layer][i])
                         continue;

                    build_scat_vector_lc(n_coef, d->n_stokes, (double *) polys, coef_l[i_layer][i], P_l[i]);
                    scat_vector_rotate_hovenier(d->n_stokes, d->mu_0, mu, c, phi - d->phi_0, P_l[i], P_l[i]);
               }
          }
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
static void get_phase_func_trun(xrtm_data *d, int i_layer, double mu, double phi, void *polys, double *P_trun, double **P_trun_l, work_data work) {

     int n_coef = MIN(d->n_coef_layer[i_layer], d->n_coef2);

     phase_func_get(d, n_coef, i_layer, mu, phi, polys, d->coef,  d->coef_l,  P_trun, P_trun_l, work);
}



/*******************************************************************************
 *
 ******************************************************************************/
static int get_phase_func_trun_all(xrtm_data *d, double mu, double phi, void *polys, double **P_trun, double ***P_trun_l, work_data work) {

     int i;

     double **P_trun_l2;

     for (i = 0; i < d->n_layers; ++i) {
          if (d->options & XRTM_OPTION_CALC_DERIVS && flags_or(d->derivs.layers[i], d->n_derivs))
               P_trun_l2 = P_trun_l[i];

          get_phase_func_trun(d, i, mu, phi, polys, P_trun[i], P_trun_l2, work);
     }

     return 0;
}



/*******************************************************************************
 *
 ******************************************************************************/
static void get_phase_func_full(xrtm_data *d, int i_layer, double mu, double phi, void *polys, double *P_full, double **P_full_l, work_data work) {

     int n_coef =       d->n_coef_layer[i_layer];

     phase_func_get(d,  n_coef, i_layer, mu, phi, polys, d->coef0, d->coef0_l, P_full, P_full_l, work);
}



/*******************************************************************************
 *
 ******************************************************************************/
static int get_phase_func_full_all(xrtm_data *d, double mu, double phi, void *polys, double **P_full, double ***P_full_l, work_data work) {

     int i;

     double **P_full_l2;

     for (i = 0; i < d->n_layers; ++i) {
          if (d->options & XRTM_OPTION_CALC_DERIVS && flags_or(d->derivs.layers[i], d->n_derivs))
               P_full_l2 = P_full_l[i];

          get_phase_func_full(d, i, mu, phi, polys, P_full[i], P_full_l2, work);
     }

     return 0;
}



/*******************************************************************************
 *
 ******************************************************************************/
static int phase_vecs_get(xrtm_data *d, int i_four, int i_layer, int n_mus1, double *mus1, double mu2, double **Y1, double *Y2, double *alpha1, int **dep_flags_phase_vecs_xx, int (*deps_func)(xrtm_data *), double ***P_pp0, double ***P_pm0, double ****P_pp0_l, double ****P_pm0_l, double **P_pp, double **P_pm, double ***P_pp_l, double ***P_pm_l, int flag, int flag2, int flag3, enum work_type work_type, work_data *work) {

     int i;
     int j;

     int flag4;

     int deps;

     int n_coef;

     int n_mus_v1;

     work_data work2;


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (d->options & XRTM_OPTION_SAVE_PHASE_MATS) {
          deps = deps_func(d);

          if (! (d->options & XRTM_OPTION_PHASE_MATRIX_LC))
               flag4 = resolve_deps1(d, i_four, i_layer,
                                     dep_flags_phase_vecs_xx, deps);
          else
               flag4 = resolve_deps2(d, i_four, i_layer,
                                     dep_flags_phase_vecs_xx, deps);

          *P_pp = P_pp0[i_four][i_layer];
          *P_pm = P_pm0[i_four][i_layer];

          if (d->options & XRTM_OPTION_CALC_DERIVS && flags_or(d->derivs.layers[i_layer], d->n_derivs)) {
               *P_pp_l = P_pp0_l[i_four][i_layer];
               *P_pm_l = P_pm0_l[i_four][i_layer];
          }

          if (! flag4)
               return 0;
     }
     else if (flag) {
#ifndef __cplusplus
          *P_pp = get_work1(work, work_type);
          *P_pm = get_work1(work, work_type);

          if (d->options & XRTM_OPTION_CALC_DERIVS && flags_or(d->derivs.layers[i_layer], d->n_derivs)) {
               *P_pp_l = get_work2(work, work_type, WORK_DERIVS_V, d->derivs.layers[i_layer]);
               *P_pm_l = get_work2(work, work_type, WORK_DERIVS_V, d->derivs.layers[i_layer]);
          }
#else
          int n_mus_v1 = n_mus1 * d->n_stokes;
          *P_pp = get_work_d1(work, n_mus_v1);
          *P_pm = get_work_d1(work, n_mus_v1);

          if (d->options & XRTM_OPTION_CALC_DERIVS && flags_or(d->derivs.layers[i_layer], d->n_derivs)) {
               *P_pp_l = get_work_d2(work, d->n_derivs, n_mus_v1);
               *P_pm_l = get_work_d2(work, d->n_derivs, n_mus_v1);
          }
#endif
     }


     work2 = *work;


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     n_coef = MIN(d->n_coef_layer[i_layer], d->n_coef2);

     if (d->options & XRTM_OPTION_PHASE_SCALAR) {
          build_phase_vecs_scalar(i_four, n_coef, n_mus1, Y1, Y2, d->n_coef2, d->coef[i_layer][0], *P_pm, *P_pp);

          if (d->options & XRTM_OPTION_CALC_DERIVS && flags_or(d->derivs.layers[i_layer], d->n_derivs)) {
               for (i = 0; i < d->n_derivs; ++i) {
                    if (! d->derivs.layers[i_layer][i])
                         continue;

                    build_phase_vecs_scalar(i_four, n_coef, n_mus1, Y1, Y2, d->n_coef2, d->coef_l[i_layer][i][0], (*P_pm_l)[i], (*P_pp_l)[i]);
               }
          }
     }
     else
     if (d->options & XRTM_OPTION_PHASE_MATRIX_GC) {
          build_phase_vecs_vector_gc(i_four, n_coef, n_mus1, d->n_stokes, mus1, mu2, d->coef[i_layer], *P_pm, *P_pp, work2);
          if (d->options & XRTM_OPTION_CALC_DERIVS && flags_or(d->derivs.layers[i_layer], d->n_derivs)) {
               for (i = 0; i < d->n_derivs; ++i) {
                    if (! d->derivs.layers[i_layer][i])
                         continue;
#ifndef USE_AD_FOR_TL_CALC_BUILD_PHASE_VECS_VECTOR_GC
                    build_phase_vecs_vector_gc(i_four, n_coef, n_mus1, d->n_stokes, mus1, mu2, d->coef_l[i_layer][i], (*P_pm_l)[i], (*P_pp_l)[i], work2);
#else
                    build_phase_vecs_vector_gc_tl_with_ad(i_four, n_coef, n_mus1, d->n_stokes, mus1, mu2, d->coef_l[i_layer][i], (*P_pm_l)[i], (*P_pp_l)[i], work2);
#endif
               }
          }
     }
#ifdef INCLUDE_DEV_SOURCE
     else
     if (d->options & XRTM_OPTION_PHASE_MATRIX_LC) {
          int n_mus_v1 = n_mus1 * d->n_stokes;

          build_phase_vecs_vector_lc (d->n_four2, n_coef, n_mus1, d->n_stokes, mus1, mu2, d->coef[i_layer], P_pm, P_pp, work2, d->n_layers * n_mus_v1);
          if (d->options & XRTM_OPTION_CALC_DERIVS && flags_or(d->derivs.layers[i_layer], d->n_derivs)) {
               for (i = 0; i < d->n_derivs; ++i) {
                    if (! d->derivs.layers[i_layer][i])
                         continue;

                    build_phase_vecs_vector_lc(d->n_four2, n_coef, n_mus1, d->n_stokes, mus1, mu2, d->coef_l[i_layer][i], *P_pm_l+i, *P_pp_l+i, work2, d->n_layers * d->n_derivs * n_mus_v1);
               }
          }
     }
#endif
#ifdef DEBUG
     else {
         fprintf(stderr, "ERROR: phase_vecs_get(): end of if / else if\n");
         exit(1);
     }
#endif


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (flag2 && i_four == 0 && d->omega[i_layer] != 0.) {
          if (check_phase_vecs_norm(n_mus1, d->n_stokes, d->qx, d->qw, *P_pp, *P_pm) != 0) {
               fprintf(stderr, "ERROR: check_phase_vecs_norm(), i_four = %d, i_layer = %d\n", i_four, i_layer);
               return -1;
          }

          if (d->options & XRTM_OPTION_CALC_DERIVS && flags_or(d->derivs.layers[i_layer], d->n_derivs)) {
               for (i = 0; i < d->n_derivs; ++i) {
                    if (! d->derivs.layers[i_layer][i])
                         continue;

                    if (check_phase_vecs_norm_l(n_mus1, d->n_stokes, d->qx, d->qw, (*P_pp_l)[i], (*P_pm_l)[i]) != 0) {
                         fprintf(stderr, "ERROR: check_phase_vecs_norm_l(), i_four = %d, i_layer = %d\n", i_four, i_layer);
                         return -1;
                    }
               }
          }
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (flag3 && d->misc_input.use_symmetric_form) {
          n_mus_v1 = n_mus1 * d->n_stokes;

          if (d->options & XRTM_OPTION_PHASE_SCALAR ||
              d->options & XRTM_OPTION_PHASE_MATRIX_GC) {
               vec_sim_trans2(n_mus_v1, *P_pp, *P_pm, alpha1);

               if (d->options & XRTM_OPTION_CALC_DERIVS && flags_or(d->derivs.layers[i_layer], d->n_derivs)) {
                    for (i = 0; i < d->n_derivs; ++i) {
                         if (! d->derivs.layers[i_layer][i])
                              continue;

                         vec_sim_trans2(n_mus_v1, (*P_pp_l)[i], (*P_pm_l)[i], alpha1);
                    }
               }
          }
          else {
               for (i = 0; i < d->n_four; ++i) {
                    vec_sim_trans2(n_mus_v1, P_pp0[i][i_layer], P_pm0[i][i_layer], alpha1);
                    if (d->options & XRTM_OPTION_CALC_DERIVS && flags_or(d->derivs.layers[i_layer], d->n_derivs)) {
                         for (j = 0; j < d->n_derivs; ++j) {
                              if (! d->derivs.layers[i_layer][j])
                                   continue;

                              vec_sim_trans2(n_mus_v1, P_pp0_l[i][i_layer][j], P_pm0_l[i][i_layer][j], alpha1);
                         }
                    }
               }
          }
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
/*
     if (d->options & XRTM_OPTION_VECTOR) {
          if (d->options & XRTM_OPTION_PHASE_SCALAR ||
              d->options & XRTM_OPTION_PHASE_MATRIX_GC) {
               dm_v_mul_D_A_2(n_mus1, d->n_stokes, *P_pp, *P_pp);

               if (d->options & XRTM_OPTION_CALC_DERIVS && flags_or(d->derivs.layers[i_layer], d->n_derivs)) {
                    for (i = 0; i < d->n_derivs; ++i) {
                         if (! d->derivs.layers[i_layer][i])
                              continue;

                         dm_v_mul_D_A_2(n_mus1, d->n_stokes, (*P_pp_l)[i], (*P_pp_l)[i]);
                    }
               }
          }
          else {
               for (i = 0; i < d->n_four; ++i) {
                    dm_v_mul_D_A_2(n_mus1, d->n_stokes, P_pp0[i][i_layer], P_pp0[i][i_layer]);

                    if (d->options & XRTM_OPTION_CALC_DERIVS && flags_or(d->derivs.layers[i_layer], d->n_derivs)) {
                         for (j = 0; j < d->n_derivs; ++j) {
                              if (! d->derivs.layers[i_layer][j])
                                   continue;

                              dm_v_mul_D_A_2(n_mus1, d->n_stokes, P_pp0_l[i][i_layer][j], P_pp0_l[i][i_layer][j]);
                         }
                    }
               }
          }
     }
*/

     return 0;
}



/*******************************************************************************
 *
 ******************************************************************************/
static int phase_vecs_get_all(xrtm_data *d, int i_four, int n_mus1, double *mus1, double mu2, double **Y1, double *Y2, double *alpha1, int **dep_flags_phase_vecs_xx, int (*deps_func)(xrtm_data *), double ***P_pp0, double ***P_pm0, double ****P_pp0_l, double ****P_pm0_l, double ***P_pp, double ***P_pm, double ****P_pp_l, double ****P_pm_l, int flag, int flag2, int flag3, enum work_type work_type, work_data *work) {

     int i;

     double *P_pp2;
     double *P_pm2;
     double **P_pp_l2;
     double **P_pm_l2;

     if (d->options & XRTM_OPTION_SAVE_PHASE_MATS) {
          *P_pp = P_pp0[i_four];
          *P_pm = P_pm0[i_four];

          if (d->options & XRTM_OPTION_CALC_DERIVS) {
               *P_pp_l = P_pp0_l[i_four];
               *P_pm_l = P_pm0_l[i_four];
          }

          for (i = 0; i < d->n_layers; ++i) {
               if (phase_vecs_get(d, i_four, i, n_mus1, mus1, mu2, Y1, Y2, alpha1, dep_flags_phase_vecs_xx, deps_func, P_pp0, P_pm0, P_pp0_l, P_pm0_l, &P_pp2, &P_pm2, &P_pp_l2, &P_pm_l2, 0, flag2, flag3, work_type, work)) {
                    fprintf(stderr, "ERROR: get_phase_vecs(), i_four = %d, i_layer = %d\n", i_four, i);
                    return -1;
               }
          }
     }
     else {
          if (flag) {
#ifndef __cplusplus
               *P_pp = get_work2(work, work_type, WORK_LAYERS_V, NULL);
               *P_pm = get_work2(work, work_type, WORK_LAYERS_V, NULL);

               if (d->options & XRTM_OPTION_CALC_DERIVS) {
                    *P_pp_l = get_work3(work, work_type, WORK_BOTH_V, d->derivs.layers);
                    *P_pm_l = get_work3(work, work_type, WORK_BOTH_V, d->derivs.layers);
               }
#else
               int n_mus_v1 = n_mus1 * d->n_stokes;

               *P_pp = get_work_d2(work, d->n_layers, n_mus_v1);
               *P_pm = get_work_d2(work, d->n_layers, n_mus_v1);

               if (d->options & XRTM_OPTION_CALC_DERIVS) {
                    *P_pp_l = get_work_d3(work, d->n_layers, d->n_derivs, n_mus_v1);
                    *P_pm_l = get_work_d3(work, d->n_layers, d->n_derivs, n_mus_v1);
               }
#endif
          }

          for (i = 0; i < d->n_layers; ++i) {
               if (d->options & XRTM_OPTION_CALC_DERIVS && flags_or(d->derivs.layers[i], d->n_derivs)) {
                    P_pp_l2 = (*P_pp_l)[i];
                    P_pm_l2 = (*P_pm_l)[i];
               }

               if (phase_vecs_get(d, i_four, i, n_mus1, mus1, mu2, Y1, Y2, alpha1, dep_flags_phase_vecs_xx, deps_func, NULL, NULL, NULL, NULL, *P_pp+i, *P_pm+i, &P_pp_l2, &P_pm_l2, 0, flag2, flag3, work_type, work)) {
                    fprintf(stderr, "ERROR: get_phase_vecs(), i_four = %d, i_layer = %d\n", i_four, i);
                    return -1;
               }
          }
     }

     return 0;
}



/*******************************************************************************
 *
 ******************************************************************************/
static int phase_vecs_deps_q0(xrtm_data *d) {

     int deps = DEP_MASK_MU_0 | DEP_MASK_COEF | DEP_MASK_COEF_L;

     if (! (d->options & XRTM_OPTION_SFI))
          deps |= DEP_MASK_UMUS;

     return deps;
}



static int get_phase_vecs_q0(xrtm_data *d, int i_four, int i_layer,
                             double **P_q0_mm, double **P_q0_pm,
                             double ***P_q0_mm_l, double ***P_q0_pm_l,
                             int flag, work_data *work) {

     double *Y_0;

     double **Y_p;

     if (d->options & XRTM_OPTION_PHASE_SCALAR) {
          get_Y_p(d, i_four, &Y_p);
          get_Y_0(d, i_four, &Y_0);
     }

     if (phase_vecs_get(d, i_four, i_layer, d->n_quad_x, d->qx, -d->mu_0, Y_p, Y_0, d->alpha1, d->dep_flags_phase_vecs_q0, phase_vecs_deps_q0, d->P_q0_mm, d->P_q0_pm, d->P_q0_mm_l, d->P_q0_pm_l, P_q0_mm, P_q0_pm, P_q0_mm_l, P_q0_pm_l, flag, 1, 1, WORK_DX, work)) {
          fprintf(stderr, "ERROR: phase_vecs_get(), i_four = %d, i_layer = %d\n", i_four, i_layer);
          return -1;
     }

     return 0;
}



static int get_phase_vecs_q0_all(xrtm_data *d, int i_four,
                                 double ***P_q0_mm, double ***P_q0_pm,
                                 double ****P_q0_mm_l, double ****P_q0_pm_l,
                                 int flag, work_data *work) {

     double *Y_0;

     double **Y_p;

     if (d->options & XRTM_OPTION_PHASE_SCALAR) {
          get_Y_p(d, i_four, &Y_p);
          get_Y_0(d, i_four, &Y_0);
     }

     if (phase_vecs_get_all(d, i_four, d->n_quad_x, d->qx, -d->mu_0, Y_p, Y_0, d->alpha1, d->dep_flags_phase_vecs_q0, phase_vecs_deps_q0, d->P_q0_mm, d->P_q0_pm, d->P_q0_mm_l, d->P_q0_pm_l, P_q0_mm, P_q0_pm, P_q0_mm_l, P_q0_pm_l, flag, 1, 1, WORK_DX, work)) {
          fprintf(stderr, "ERROR: phase_vecs_get_all(), i_four = %d\n", i_four);
          return -1;
     }

     return 0;
}



/*******************************************************************************
 *
 ******************************************************************************/
static int phase_vecs_deps_u0(xrtm_data *d) {

     return DEP_MASK_MU_0 | DEP_MASK_UMUS | DEP_MASK_COEF | DEP_MASK_COEF_L;
}


/*
static int get_phase_vecs_u0(xrtm_data *d, int i_four, int i_layer,
                             double **P_u0_mm, double **P_u0_pm,
                             double ***P_u0_mm_l, double ***P_u0_pm_l,
                             int flag, work_data *work) {

     double *Y_0;

     double **Y_u;

     if (d->options & XRTM_OPTION_PHASE_SCALAR) {
          get_Y_u(d, i_four, &Y_u);
          get_Y_0(d, i_four, &Y_0);
     }

     if (phase_vecs_get(d, i_four, i_layer, d->n_umus, d->umus, -d->mu_0, Y_u, Y_0, d->alpha1 + d->n_quad_v, d->dep_flags_phase_vecs_u0, phase_vecs_deps_u0, d->P_u0_mm, d->P_u0_pm, d->P_u0_mm_l, d->P_u0_pm_l, P_u0_mm, P_u0_mp, P_u0_mm_l, Q_u0_mp, flag, 0, 1, WORK_DU, work)) {
          fprintf(stderr, "ERROR: phase_vecs_get(), i_four = %d, i_layer = %d\n", i_four, i_layer);
          return -1;
     }

     return 0;
}
*/


static int get_phase_vecs_u0_all(xrtm_data *d, int i_four,
                                 double ***P_u0_mm, double ***P_u0_pm,
                                 double ****P_u0_mm_l, double ****P_u0_pm_l,
                                 int flag, work_data *work) {

     double *Y_0;

     double **Y_u;

     if (d->options & XRTM_OPTION_PHASE_SCALAR) {
          get_Y_u(d, i_four, &Y_u);
          get_Y_0(d, i_four, &Y_0);
     }

     if (phase_vecs_get_all(d, i_four, d->n_umus, d->umus, -d->mu_0, Y_u, Y_0, d->alpha1 + d->n_quad_v, d->dep_flags_phase_vecs_u0, phase_vecs_deps_u0, d->P_u0_mm, d->P_u0_pm, d->P_u0_mm_l, d->P_u0_pm_l, P_u0_mm, P_u0_pm, P_u0_mm_l, P_u0_pm_l, flag, 0, 1, WORK_DU, work)) {
          fprintf(stderr, "ERROR: phase_vecs_get_all(), i_four = %d\n", i_four);
          return -1;
     }

     return 0;
}



/*******************************************************************************
 *
 ******************************************************************************/
static int phase_mats_get(xrtm_data *d, int i_four, int i_layer, int n_mus1, int n_mus2, double *mus1, double *mus2, double **Y1, double **Y2, double *alpha1, double *alpha2, int **dep_flags_phase_mats_xx, int (*deps_func)(xrtm_data *), double ****P_pp0, double ****P_mp0, double ****P_mm0, double ****P_pm0, double *****P_pp0_l, double *****P_mp0_l, double *****P_mm0_l, double *****P_pm0_l, double ***P_pp, double ***P_mp, double ***P_mm, double ***P_pm, double ****P_pp_l, double ****P_mp_l, double ****P_mm_l, double ****P_pm_l, int flag, int flag2, int flag3, enum work_type work_type, work_data *work) {

     int i;
     int j;

     int flag4;

     int deps;

     int n_coef;

     int n_mus_v1;
     int n_mus_v2;

     work_data work2;


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (d->options & XRTM_OPTION_SAVE_PHASE_MATS) {
          deps = deps_func(d);

          if (! (d->options & XRTM_OPTION_PHASE_MATRIX_LC))
               flag4 = resolve_deps1(d, i_four, i_layer,
                                     dep_flags_phase_mats_xx, deps);
          else
               flag4 = resolve_deps2(d, i_four, i_layer,
                                     dep_flags_phase_mats_xx, deps);

          *P_pp = P_pp0[i_four][i_layer];
          *P_mp = P_mp0[i_four][i_layer];

          if (d->options & XRTM_OPTION_VECTOR) {
               *P_mm = P_mm0[i_four][i_layer];
               *P_pm = P_pm0[i_four][i_layer];
          }

          if (d->options & XRTM_OPTION_CALC_DERIVS && flags_or(d->derivs.layers[i_layer], d->n_derivs)) {
               *P_pp_l = P_pp0_l[i_four][i_layer];
               *P_mp_l = P_mp0_l[i_four][i_layer];

               if (d->options & XRTM_OPTION_VECTOR) {
                    *P_mm_l = P_mm0_l[i_four][i_layer];
                    *P_pm_l = P_pm0_l[i_four][i_layer];
               }
          }

          if (! flag4)
               return 0;
     }
     else if (flag) {
#ifndef __cplusplus
          *P_pp = get_work1(work, work_type);
          *P_mp = get_work1(work, work_type);

          if (d->options & XRTM_OPTION_VECTOR) {
               *P_mm = get_work1(work, work_type);
               *P_pm = get_work1(work, work_type);
          }

          if (d->options & XRTM_OPTION_CALC_DERIVS && flags_or(d->derivs.layers[i_layer], d->n_derivs)) {
               *P_pp_l = get_work2(work, work_type, WORK_DERIVS_V, d->derivs.layers[i_layer]);
               *P_mp_l = get_work2(work, work_type, WORK_DERIVS_V, d->derivs.layers[i_layer]);

               if (d->options & XRTM_OPTION_VECTOR) {
                    *P_mm_l = get_work2(work, work_type, WORK_DERIVS_V, d->derivs.layers[i_layer]);
                    *P_pm_l = get_work2(work, work_type, WORK_DERIVS_V, d->derivs.layers[i_layer]);
               }
          }
#else
          int n_mus_v1 = n_mus1 * d->n_stokes;
          int n_mus_v2 = n_mus2 * d->n_stokes;

          *P_pp = get_work_d2(work, n_mus_v1, n_mus_v2);
          *P_mp = get_work_d2(work, n_mus_v1, n_mus_v2);

          if (d->options & XRTM_OPTION_VECTOR) {
               *P_mm = get_work_d2(work, n_mus_v1, n_mus_v2);
               *P_pm = get_work_d2(work, n_mus_v1, n_mus_v2);
          }

          if (d->options & XRTM_OPTION_CALC_DERIVS && flags_or(d->derivs.layers[i_layer], d->n_derivs)) {
               *P_pp_l = get_work_d3(work, d->n_derivs, n_mus_v1, n_mus_v2);
               *P_mp_l = get_work_d3(work, d->n_derivs, n_mus_v1, n_mus_v2);

               if (d->options & XRTM_OPTION_VECTOR) {
                    *P_mm_l = get_work_d3(work, d->n_derivs, n_mus_v1, n_mus_v2);
                    *P_pm_l = get_work_d3(work, d->n_derivs, n_mus_v1, n_mus_v2);
               }
          }
#endif
     }


     work2 = *work;


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     n_coef = MIN(d->n_coef_layer[i_layer], d->n_coef2);

     if (d->options & XRTM_OPTION_PHASE_SCALAR) {
          build_phase_mats_scalar(i_four, n_coef, n_mus1, n_mus2, Y1, Y2, d->n_coef2, d->coef[i_layer][0], *P_pp, *P_mp);

          if (d->options & XRTM_OPTION_CALC_DERIVS && flags_or(d->derivs.layers[i_layer], d->n_derivs)) {
               for (i = 0; i < d->n_derivs; ++i) {
                    if (! d->derivs.layers[i_layer][i])
                         continue;
#ifndef USE_AD_FOR_TL_CALC_BUILD_PHASE_MATS_SCALAR
                    build_phase_mats_scalar(i_four, n_coef, n_mus1, n_mus2, Y1, Y2, d->n_coef2, d->coef_l[i_layer][i][0], (*P_pp_l)[i], (*P_mp_l)[i]);
#else
                    build_phase_mats_scalar_tl_with_ad(i_four, n_coef, n_mus1, n_mus2, Y1, Y2, d->n_coef2, d->coef_l[i_layer][i][0], (*P_pp_l)[i], (*P_mp_l)[i], work2);
#endif
               }
          }
     }
     else
     if (d->options & XRTM_OPTION_PHASE_MATRIX_GC) {
          build_phase_mats_vector_gc(i_four, n_coef, n_mus1, n_mus2, d->n_stokes, mus1, mus2, d->coef[i_layer], *P_pp, *P_mp, *P_mm, *P_pm, d->options & XRTM_OPTION_VECTOR, work2);
          if (d->options & XRTM_OPTION_CALC_DERIVS && flags_or(d->derivs.layers[i_layer], d->n_derivs)) {
               for (i = 0; i < d->n_derivs; ++i) {
                    if (! d->derivs.layers[i_layer][i])
                         continue;
#ifndef USE_AD_FOR_TL_CALC_BUILD_PHASE_MATS_VECTOR_GC
                    build_phase_mats_vector_gc(i_four, n_coef, n_mus1, n_mus2, d->n_stokes, mus1, mus2, d->coef_l[i_layer][i], (*P_pp_l)[i], (*P_mp_l)[i], (*P_mm_l)[i], (*P_pm_l)[i], d->options & XRTM_OPTION_VECTOR, work2);
#else
                    build_phase_mats_vector_gc_tl_with_ad(i_four, n_coef, n_mus1, n_mus2, d->n_stokes, mus1, mus2, d->coef_l[i_layer][i], (*P_pp_l)[i], (*P_mp_l)[i], (*P_mm_l)[i], (*P_pm_l)[i], d->options & XRTM_OPTION_VECTOR, work2);
#endif
               }
          }
     }
#ifdef INCLUDE_DEV_SOURCE
     else
     if (d->options & XRTM_OPTION_PHASE_MATRIX_LC) {
          int n_mus_v1 = n_mus1 * d->n_stokes;

          build_phase_mats_vector_lc(d->n_four2, n_coef, n_mus1, n_mus2, d->n_stokes, mus1, mus2, d->coef[i_layer], P_pp, P_mp, P_mm, P_pm, d->options & XRTM_OPTION_VECTOR, work2, d->n_layers * n_mus_v1 * d->n_quad_v_x);
          if (d->options & XRTM_OPTION_CALC_DERIVS && flags_or(d->derivs.layers[i_layer], d->n_derivs)) {
               for (i = 0; i < d->n_derivs; ++i) {
                    if (! d->derivs.layers[i_layer][i])
                         continue;

                    build_phase_mats_vector_lc(d->n_four2, n_coef, n_mus1, n_mus2, d->n_stokes, mus1, mus2, d->coef_l[i_layer][i], *P_pp_l+i, *P_mp_l+i, *P_mm_l+i, *P_pm_l+i, d->options & XRTM_OPTION_VECTOR, work2, d->n_layers * d->n_derivs * n_mus_v1 * d->n_quad_v_x);
               }
          }
     }
#endif
#ifdef DEBUG
     else {
         fprintf(stderr, "ERROR: phase_mats_get(): end of if / else if\n");
         exit(1);
     }
#endif

     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (flag2 && i_four == 0 && d->omega[i_layer] != 0.) {
          if (check_phase_mats_norm(n_mus1, n_mus2, d->n_stokes, d->qx, d->qw, *P_pp, *P_mp) != 0) {
               fprintf(stderr, "ERROR: check_phase_mats_norm(), i_four = %d, i_layer = %d\n", i_four, i_layer);
               return -1;
          }

          if (d->options & XRTM_OPTION_CALC_DERIVS && flags_or(d->derivs.layers[i_layer], d->n_derivs)) {
               for (i = 0; i < d->n_derivs; ++i) {
                    if (! d->derivs.layers[i_layer][i])
                         continue;

                    if (check_phase_mats_norm_l(n_mus1, n_mus2, d->n_stokes, d->qx, d->qw, (*P_pp_l)[i], (*P_mp_l)[i]) != 0) {
                         fprintf(stderr, "ERROR: check_phase_mats_norm_l(), i_four = %d, i_layer = %d\n", i_four, i_layer);
                         return -1;
                    }
               }
          }
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (flag3 && d->misc_input.use_symmetric_form) {
          n_mus_v1 = n_mus1 * d->n_stokes;
          n_mus_v2 = n_mus2 * d->n_stokes;

          if (d->options & XRTM_OPTION_PHASE_SCALAR ||
              d->options & XRTM_OPTION_PHASE_MATRIX_GC) {
               mat_sim_trans2(n_mus_v1, n_mus_v2, *P_pp, *P_mp, alpha1, alpha2);

               if (d->options & XRTM_OPTION_VECTOR)
                    mat_sim_trans2(n_mus_v1, n_mus_v2, *P_mm, *P_pm, alpha1, alpha2);

               if (d->options & XRTM_OPTION_CALC_DERIVS && flags_or(d->derivs.layers[i_layer], d->n_derivs)) {
                    for (i = 0; i < d->n_derivs; ++i) {
                         if (! d->derivs.layers[i_layer][i])
                              continue;

                         mat_sim_trans2(n_mus_v1, n_mus_v2, (*P_pp_l)[i], (*P_mp_l)[i], alpha1, alpha2);

                         if (d->options & XRTM_OPTION_VECTOR)
                              mat_sim_trans2(n_mus_v1, n_mus_v2, (*P_mm_l)[i], (*P_pm_l)[i], alpha1, alpha2);
                    }
               }
          }
          else {
               for (i = 0; i < d->n_four; ++i) {
                    mat_sim_trans2(n_mus_v1, n_mus_v2, P_pp0[i][i_layer], P_mp0[i][i_layer], alpha1, alpha2);

                    if (d->options & XRTM_OPTION_VECTOR)
                         mat_sim_trans2(n_mus_v1, n_mus_v2, P_mm0[i][i_layer], P_pm0[i][i_layer], alpha1, alpha2);

                    if (d->options & XRTM_OPTION_CALC_DERIVS && flags_or(d->derivs.layers[i_layer], d->n_derivs)) {
                         for (j = 0; j < d->n_derivs; ++j) {
                              if (! d->derivs.layers[i_layer][j])
                                   continue;

                              mat_sim_trans2(n_mus_v1, n_mus_v2, P_pp0_l[i][i_layer][j], P_mp0_l[i][i_layer][j], alpha1, alpha2);

                              if (d->options & XRTM_OPTION_VECTOR)
                                   mat_sim_trans2(n_mus_v1, n_mus_v2, P_mm0_l[i][i_layer][j], P_pm0_l[i][i_layer][j], alpha1, alpha2);
                         }
                    }
               }
          }
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
/*
     if (d->options & XRTM_OPTION_VECTOR) {
          if (d->options & XRTM_OPTION_PHASE_SCALAR ||
              d->options & XRTM_OPTION_PHASE_MATRIX_GC) {
               dmat_mul_D_A2_2(n_mus1, d->n_stokes, n_mus2, d->n_stokes, *P_mp, *P_mp);

               if (d->options & XRTM_OPTION_CALC_DERIVS && flags_or(d->derivs.layers[i_layer], d->n_derivs)) {
                    for (i = 0; i < d->n_derivs; ++i) {
                         if (! d->derivs.layers[i_layer][i])
                              continue;

                         dmat_mul_D_A2_2(n_mus1, d->n_stokes, n_mus2, d->n_stokes, (*P_mp_l)[i], (*P_mp_l)[i]);
                    }
               }
          }
          else {
               for (i = 0; i < d->n_four; ++i) {
                    dmat_mul_D_A2_2(n_mus1, d->n_stokes, n_mus2, d->n_stokes, P_mp0[i][i_layer], P_mp0[i][i_layer]);

                    if (d->options & XRTM_OPTION_CALC_DERIVS && flags_or(d->derivs.layers[i_layer], d->n_derivs)) {
                         for (j = 0; j < d->n_derivs; ++j) {
                              if (! d->derivs.layers[i_layer][j])
                                   continue;

                              dmat_mul_D_A2_2(n_mus1, d->n_stokes, n_mus2, d->n_stokes, P_mp0_l[i][i_layer][j], P_mp0_l[i][i_layer][j]);
                         }
                    }
               }
          }
     }
*/

     return 0;
}



/*******************************************************************************
 *
 ******************************************************************************/
static int phase_mats_get_all(xrtm_data *d, int i_four, int n_mus1, int n_mus2, double *mus1, double *mus2, double **Y1, double **Y2, double *alpha1, double *alpha2, int **dep_flags_phase_mats_xx, int (*deps_func)(xrtm_data *), double ****P_pp0, double ****P_mp0, double ****P_mm0, double ****P_pm0, double *****P_pp_l0, double *****P_mp0_l, double *****P_mm0_l, double *****P_pm_l0, double ****P_pp, double ****P_mp, double ****P_mm, double ****P_pm, double *****P_pp_l, double *****P_mp_l, double *****P_mm_l, double *****P_pm_l, int flag, int flag2, int flag3, enum work_type work_type, work_data *work) {

     int i;

     double **P_pp2;
     double **P_mp2;
     double **P_mm2;
     double **P_pm2;

     double ***P_pp_l2;
     double ***P_mp_l2;
     double ***P_mm_l2;
     double ***P_pm_l2;

     if (d->options & XRTM_OPTION_SAVE_PHASE_MATS) {
          *P_pp = P_pp0[i_four];
          *P_mp = P_mp0[i_four];

          if (d->options & XRTM_OPTION_VECTOR) {
               *P_mm = P_mm0[i_four];
               *P_pm = P_pm0[i_four];
          }

          if (d->options & XRTM_OPTION_CALC_DERIVS) {
               *P_pp_l = P_pp_l0[i_four];
               *P_mp_l = P_mp0_l[i_four];

               if (d->options & XRTM_OPTION_VECTOR) {
                    *P_mm_l = P_mm0_l[i_four];
                    *P_pm_l = P_pm_l0[i_four];
               }
          }

          for (i = 0; i < d->n_layers; ++i) {
               if (phase_mats_get(d, i_four, i, n_mus1, n_mus2, mus1, mus2, Y1, Y2, alpha1, alpha2, dep_flags_phase_mats_xx, deps_func, P_pp0, P_mp0, P_mm0, P_pm0, P_pp_l0, P_mp0_l, P_mm0_l, P_pm_l0, &P_pp2, &P_mp2, &P_mm2, &P_pm2, &P_pp_l2, &P_mp_l2, &P_mm_l2, &P_pm_l2, 0, flag2, flag3, work_type, work)) {
                    fprintf(stderr, "ERROR: phase_mats_get(), i_four = %d, i_layer = %d\n", i_four, i);
                    return -1;
               }
          }
     }
     else {
          if (flag) {
#ifndef __cplusplus
               *P_pp = get_work2(work, work_type, WORK_LAYERS_V, NULL);
               *P_mp = get_work2(work, work_type, WORK_LAYERS_V, NULL);

               if (d->options & XRTM_OPTION_VECTOR) {
                    *P_mm = get_work2(work, work_type, WORK_LAYERS_V, NULL);
                    *P_pm = get_work2(work, work_type, WORK_LAYERS_V, NULL);
               }

               if (d->options & XRTM_OPTION_CALC_DERIVS) {
                    *P_pp_l = get_work3(work, work_type, WORK_BOTH_V, d->derivs.layers);
                    *P_mp_l = get_work3(work, work_type, WORK_BOTH_V, d->derivs.layers);

                    if (d->options & XRTM_OPTION_VECTOR) {
                         *P_mm_l = get_work3(work, work_type, WORK_BOTH_V, d->derivs.layers);
                         *P_pm_l = get_work3(work, work_type, WORK_BOTH_V, d->derivs.layers);
                    }
               }
#else
               int n_mus_v1 = n_mus1 * d->n_stokes;
               int n_mus_v2 = n_mus2 * d->n_stokes;

               *P_pp = get_work_d3(work, d->n_layers, n_mus_v1, n_mus_v2);
               *P_mp = get_work_d3(work, d->n_layers, n_mus_v1, n_mus_v2);

               if (d->options & XRTM_OPTION_VECTOR) {
                    *P_mm = get_work_d3(work, d->n_layers, n_mus_v1, n_mus_v2);
                    *P_pm = get_work_d3(work, d->n_layers, n_mus_v1, n_mus_v2);
               }

               if (d->options & XRTM_OPTION_CALC_DERIVS) {
                    *P_pp_l = get_work_d4(work, d->n_layers, d->n_derivs, n_mus_v1, n_mus_v2);
                    *P_mp_l = get_work_d4(work, d->n_layers, d->n_derivs, n_mus_v1, n_mus_v2);

                    if (d->options & XRTM_OPTION_VECTOR) {
                         *P_mm_l = get_work_d4(work, d->n_layers, d->n_derivs, n_mus_v1, n_mus_v2);
                         *P_pm_l = get_work_d4(work, d->n_layers, d->n_derivs, n_mus_v1, n_mus_v2);
                    }
               }
#endif
          }

          for (i = 0; i < d->n_layers; ++i) {
               if (d->options & XRTM_OPTION_VECTOR) {
                    P_mm2 = (*P_mm)[i];
                    P_pm2 = (*P_pm)[i];
               }

               if (d->options & XRTM_OPTION_CALC_DERIVS && flags_or(d->derivs.layers[i], d->n_derivs)) {
                    P_pp_l2 = (*P_pp_l)[i];
                    P_mp_l2 = (*P_mp_l)[i];

                    if (d->options & XRTM_OPTION_VECTOR) {
                         P_mm_l2 = (*P_mm_l)[i];
                         P_pm_l2 = (*P_pm_l)[i];
                    }
               }

               if (phase_mats_get(d, i_four, i, n_mus1, n_mus2, mus1, mus2, Y1, Y2, alpha1, alpha2, dep_flags_phase_mats_xx, deps_func, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, *P_pp+i, *P_mp+i, &P_mm2, &P_pm2, &P_pp_l2, &P_mp_l2, &P_mm_l2, &P_pm_l2, 0, flag2, flag3, work_type, work)) {
                    fprintf(stderr, "ERROR: phase_mats_get(), i_four = %d, i_layer = %d\n", i_four, i);
                    return -1;
               }
          }
     }

     return 0;
}



/*******************************************************************************
 *
 ******************************************************************************/
static int phase_mats_deps_qq(xrtm_data *d) {

     int deps = DEP_MASK_COEF | DEP_MASK_COEF_L;

     if (! (d->options & XRTM_OPTION_SFI))
          deps |= DEP_MASK_UMUS;

     return deps;
}



static int get_phase_mats_qq(xrtm_data *d, int i_four, int i_layer,
                             double ***P_qq_pp, double ***P_qq_mp,
                             double ***P_qq_mm, double ***P_qq_pm,
                             double ****P_qq_pp_l, double ****P_qq_mp_l,
                             double ****P_qq_mm_l, double ****P_qq_pm_l,
                             int flag, work_data *work) {

     double **Y_p;

     if (d->options & XRTM_OPTION_PHASE_SCALAR)
          get_Y_p(d, i_four, &Y_p);

     if (phase_mats_get(d, i_four, i_layer, d->n_quad_x, d->n_quad_x, d->qx, d->qx, Y_p, Y_p, d->alpha1, d->alpha2, d->dep_flags_phase_mats_qq, phase_mats_deps_qq, d->P_qq_pp, d->P_qq_mp, d->P_qq_mm, d->P_qq_pm, d->P_qq_pp_l, d->P_qq_mp_l, d->P_qq_mm_l, d->P_qq_pm_l, P_qq_pp, P_qq_mp, P_qq_mm, P_qq_pm, P_qq_pp_l, P_qq_mp_l, P_qq_mm_l, P_qq_pm_l, flag, 1, 1, WORK_DXX, work)) {
          fprintf(stderr, "ERROR: phase_mats_get(), i_four = %d, i_layer = %d,\n", i_four, i_layer);
          return -1;
     }

     return 0;
}



static int get_phase_mats_qq_all(xrtm_data *d, int i_four,
                                 double ****P_qq_pp, double ****P_qq_mp,
                                 double ****P_qq_mm, double ****P_qq_pm,
                                 double *****P_qq_pp_l, double *****P_qq_mp_l,
                                 double *****P_qq_mm_l, double *****P_qq_pm_l,
                                 int flag, work_data *work) {

     double **Y_p;

     if (d->options & XRTM_OPTION_PHASE_SCALAR)
          get_Y_p(d, i_four, &Y_p);

     if (phase_mats_get_all(d, i_four, d->n_quad_x, d->n_quad_x, d->qx, d->qx, Y_p, Y_p, d->alpha1, d->alpha2, d->dep_flags_phase_mats_qq, phase_mats_deps_qq, d->P_qq_pp, d->P_qq_mp, d->P_qq_mm, d->P_qq_pm, d->P_qq_pp_l, d->P_qq_mp_l, d->P_qq_mm_l, d->P_qq_pm_l, P_qq_pp, P_qq_mp, P_qq_mm, P_qq_pm, P_qq_pp_l, P_qq_mp_l, P_qq_mm_l, P_qq_pm_l, flag, 1, 1, WORK_DXX, work)) {
          fprintf(stderr, "ERROR: phase_mats_get_all(), i_four = %d\n", i_four);
          return -1;
     }

     return 0;
}



/*******************************************************************************
 *
 ******************************************************************************/
static int phase_mats_deps_uq(xrtm_data *d) {

     return DEP_MASK_UMUS | DEP_MASK_COEF | DEP_MASK_COEF_L;
}


/*
static int get_phase_mats_uq(xrtm_data *d, int i_four, int i_layer,
                             double ***P_uq_pp, double ***P_uq_mp,
                             double ***P_uq_mm, double ***P_uq_pm,
                             double ****P_uq_pp_l, double ****P_uq_mp_l,
                             double ****P_uq_mm_l, double ****P_uq_pm_l,
                             int flag, work_data *work) {

     double **Y_p;
     double **Y_u;

     if (d->options & XRTM_OPTION_PHASE_SCALAR) {
          get_Y_p(d, i_four, &Y_p);
          get_Y_u(d, i_four, &Y_u);
     }

     if (phase_mats_get(d, i_four, i_layer, d->n_umus, d->n_quad_x, d->umus, d->qx, Y_u, Y_p, d->alpha1 + d->n_quad_v, d->alpha2, d->dep_flags_phase_mats_uq, phase_mats_deps_uq, d->P_uq_pp, d->P_uq_mp, d->P_uq_mm, d->P_uq_pm, d->P_uq_pp_l, d->P_uq_mp_l, d->P_uq_mm_l, d->P_uq_pm_l, P_uq_pp, P_uq_mp, P_uq_mm, P_uq_pm, P_uq_pp_l, P_uq_mp_l, P_uq_mm_l, P_uq_pm_l, flag, 0, 1, WORK_DUX, work)) {
          fprintf(stderr, "ERROR: phase_mats_get(), i_four = %d, i_layer = %d\n", i_four, i_layer);
          return -1;
     }

     return 0;
}
*/


static int get_phase_mats_uq_all(xrtm_data *d, int i_four,
                                 double ****P_uq_pp, double ****P_uq_mp,
                                 double ****P_uq_mm, double ****P_uq_pm,
                                 double *****P_uq_pp_l, double *****P_uq_mp_l,
                                 double *****P_uq_mm_l, double *****P_uq_pm_l,
                                 int flag, work_data *work) {

     double **Y_p;
     double **Y_u;

     if (d->options & XRTM_OPTION_PHASE_SCALAR) {
          get_Y_p(d, i_four, &Y_p);
          get_Y_u(d, i_four, &Y_u);
     }

     if (phase_mats_get_all(d, i_four, d->n_umus, d->n_quad_x, d->umus, d->qx, Y_u, Y_p, d->alpha1 + d->n_quad_v, d->alpha2, d->dep_flags_phase_mats_uq, phase_mats_deps_uq, d->P_uq_pp, d->P_uq_mp, d->P_uq_mm, d->P_uq_pm, d->P_uq_pp_l, d->P_uq_mp_l, d->P_uq_mm_l, d->P_uq_pm_l, P_uq_pp, P_uq_mp, P_uq_mm, P_uq_pm, P_uq_pp_l, P_uq_mp_l, P_uq_mm_l, P_uq_pm_l, flag, 0, 1, WORK_DUX, work)) {
          fprintf(stderr, "ERROR: phase_mats_get_all(), i_four = %d\n", i_four);
          return -1;
     }

     return 0;
}



/*******************************************************************************
 *
 ******************************************************************************/
static int local_r_t_u_w_deps(xrtm_data *d) {

     return phase_mats_deps_qq(d) | DEP_MASK_OMEGA | DEP_MASK_OMEGA_L;
}


/*
static
*/
int get_local_r_t_u_w
     (xrtm_data *d, int i_four, int i_layer,
      double ***r_p, double ***t_p, double ***r_m, double ***t_m,
      double ****r_p_l, double ****t_p_l, double ****r_m_l, double ****t_m_l,
      int flag, save_tree_data save_tree, work_data *work) {

     uchar *derivs_layers;

     int flag2;

     int deps;

     double *omega_l;

     double **P_qq_pp;
     double **P_qq_mp;
     double **P_qq_mm;
     double **P_qq_pm;

     double ***P_qq_pp_l;
     double ***P_qq_mp_l;
     double ***P_qq_mm_l;
     double ***P_qq_pm_l;

     forward_save_get_local_r_t_u_w_data *save;

     work_data work2;


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (d->options & XRTM_OPTION_REVERSE_DERIVS) {
          save_tree_encode_s(&save_tree, "get_local_r_t_u_w");

          if (save_tree_retrieve_data(&save_tree, forward_save_get_local_r_t_u_w_data, &save))
               forward_save_get_local_r_t_u_w_alloc(save, d->options, d->n_quad_v_x);
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (d->options & XRTM_OPTION_SAVE_LOCAL_R_T) {
          deps = local_r_t_u_w_deps(d);

          flag2 = resolve_deps1(d, i_four, i_layer,
                                d->dep_flags_local_r_t_u_w, deps);

          *r_p = d->r_p[i_four][i_layer];
          *t_p = d->t_p[i_four][i_layer];

          if (d->options & XRTM_OPTION_VECTOR) {
               *r_m = d->r_m[i_four][i_layer];
               *t_m = d->t_m[i_four][i_layer];
          }

          if (d->options & XRTM_OPTION_CALC_DERIVS && flags_or(d->derivs.layers[i_layer], d->n_derivs)) {
               *r_p_l = d->r_p_l[i_four][i_layer];
               *t_p_l = d->t_p_l[i_four][i_layer];

               if (d->options & XRTM_OPTION_VECTOR) {
                    *r_m_l = d->r_m_l[i_four][i_layer];
                    *t_m_l = d->t_m_l[i_four][i_layer];
               }
          }

          if (! flag2)
               return 0;
     }
     else if (flag) {
          *r_p = get_work1(work, WORK_DXX);
          *t_p = get_work1(work, WORK_DXX);

          if (d->options & XRTM_OPTION_VECTOR) {
               *r_m = get_work1(work, WORK_DXX);
               *t_m = get_work1(work, WORK_DXX);
          }

          if (d->options & XRTM_OPTION_CALC_DERIVS && flags_or(d->derivs.layers[i_layer], d->n_derivs)) {
               *r_p_l = get_work2(work, WORK_DXX, WORK_DERIVS_V, d->derivs.layers[i_layer]);
               *t_p_l = get_work2(work, WORK_DXX, WORK_DERIVS_V, d->derivs.layers[i_layer]);

               if (d->options & XRTM_OPTION_VECTOR) {
                    *r_m_l = get_work2(work, WORK_DXX, WORK_DERIVS_V, d->derivs.layers[i_layer]);
                    *t_m_l = get_work2(work, WORK_DXX, WORK_DERIVS_V, d->derivs.layers[i_layer]);
               }
          }
     }


     work2 = *work;


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     update_opt_props(d, i_layer, work2);


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (get_phase_mats_qq(d, i_four, i_layer, &P_qq_pp, &P_qq_mp, &P_qq_mm, &P_qq_pm, &P_qq_pp_l, &P_qq_mp_l, &P_qq_mm_l, &P_qq_pm_l, 1, &work2)) {
          fprintf(stderr, "ERROR: get_phase_mats_qq()\n");
          return -1;
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (d->options & XRTM_OPTION_CALC_DERIVS)
          derivs_layers = d->derivs.layers[i_layer];

     if (d->options & XRTM_OPTION_CALC_DERIVS && flags_or(d->derivs.layers[i_layer], d->n_derivs))
          omega_l  = d->omega_l [i_layer];

     build_local_r_and_t(i_four, d->n_quad_v_x, d->n_derivs, d->qx_v, d->qw_v, d->omega[i_layer], omega_l, P_qq_pp, P_qq_mp, *r_p, *t_p, P_qq_pp_l, P_qq_mp_l, *r_p_l, *t_p_l, derivs_layers, work2);

     if (d->options & XRTM_OPTION_VECTOR)
          build_local_r_and_t(i_four, d->n_quad_v_x, d->n_derivs, d->qx_v, d->qw_v, d->omega[i_layer], omega_l, P_qq_mm, P_qq_pm, *r_m, *t_m, P_qq_mm_l, P_qq_pm_l, *r_m_l, *t_m_l, derivs_layers, work2);


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (d->options & XRTM_OPTION_REVERSE_DERIVS) {
          dmat_copy(save->P_qq_pp, P_qq_pp, d->n_quad_v_x, d->n_quad_v_x);
          dmat_copy(save->P_qq_mp, P_qq_mp, d->n_quad_v_x, d->n_quad_v_x);

          if (d->options & XRTM_OPTION_VECTOR) {
               dmat_copy(save->P_qq_mm, P_qq_mm, d->n_quad_v_x, d->n_quad_v_x);
               dmat_copy(save->P_qq_pm, P_qq_pm, d->n_quad_v_x, d->n_quad_v_x);
          }
     }

#ifdef USE_AD_FOR_TL_CALC_GET_LOCAL_R_T_U_W
     get_local_r_t_u_w_tl_with_ad(d, i_four, i_layer, *r_p, *t_p, *r_m, *t_m, *r_p_l, *t_p_l, *r_m_l, *t_m_l, flag, save_tree, work2);
#endif

     return 0;
}



/*******************************************************************************
 *
 ******************************************************************************/
static int get_local_r_t_u_w_all
     (xrtm_data *d, int i_four,
      double ****r_p, double ****t_p, double ****r_m, double ****t_m,
      double *****r_p_l, double *****t_p_l, double *****r_m_l, double *****t_m_l,
      int flag, save_tree_data save_tree, work_data *work) {

     int i;

     double **r_p2;
     double **t_p2;
     double **r_m2;
     double **t_m2;

     double ***r_p_l2;
     double ***t_p_l2;
     double ***r_m_l2;
     double ***t_m_l2;


     if (d->options & XRTM_OPTION_REVERSE_DERIVS)
          save_tree_encode_s(&save_tree, "get_local_r_t_u_w_all");


     if (d->options & XRTM_OPTION_SAVE_LOCAL_R_T) {
          *r_p = d->r_p[i_four];
          *t_p = d->t_p[i_four];

          if (d->options & XRTM_OPTION_VECTOR) {
               *r_m = d->r_m[i_four];
               *t_m = d->t_m[i_four];
          }

          if (d->options & XRTM_OPTION_CALC_DERIVS) {
               *r_p_l = d->r_p_l[i_four];
               *t_p_l = d->t_p_l[i_four];

               if (d->options & XRTM_OPTION_VECTOR) {
                    *r_m_l = d->r_m_l[i_four];
                    *t_m_l = d->t_m_l[i_four];
               }
          }

          for (i = 0; i < d->n_layers; ++i) {
               if (get_local_r_t_u_w(d, i_four, i, &r_p2, &t_p2, &r_m2, &t_m2, &r_p_l2, &t_p_l2, &r_m_l2, &t_m_l2, 0, save_tree, work)) {
                    fprintf(stderr, "ERROR: get_local_r_t_u_w()\n");
                    return -1;
               }
          }
     }
     else {
          if (flag) {
               *r_p = get_work2(work, WORK_DXX, WORK_LAYERS_V, NULL);
               *t_p = get_work2(work, WORK_DXX, WORK_LAYERS_V, NULL);

               if (d->options & XRTM_OPTION_VECTOR) {
                    *r_m = get_work2(work, WORK_DXX, WORK_LAYERS_V, NULL);
                    *t_m = get_work2(work, WORK_DXX, WORK_LAYERS_V, NULL);
               }

               if (d->options & XRTM_OPTION_CALC_DERIVS) {
                    *r_p_l = get_work3(work, WORK_DXX, WORK_BOTH_V, d->derivs.layers);
                    *t_p_l = get_work3(work, WORK_DXX, WORK_BOTH_V, d->derivs.layers);

                    if (d->options & XRTM_OPTION_VECTOR) {
                         *r_m_l = get_work3(work, WORK_DXX, WORK_BOTH_V, d->derivs.layers);
                         *t_m_l = get_work3(work, WORK_DXX, WORK_BOTH_V, d->derivs.layers);
                    }
               }
          }

          for (i = 0; i < d->n_layers; ++i) {
               if (d->options & XRTM_OPTION_REVERSE_DERIVS)
                   save_tree_recode_i(&save_tree, i, i == 0);

               if (d->options & XRTM_OPTION_VECTOR) {
                    r_m2 = (*r_m)[i];
                    t_m2 = (*t_m)[i];
               }

               if (d->options & XRTM_OPTION_CALC_DERIVS && flags_or(d->derivs.layers[i], d->n_derivs)) {
                    r_p_l2 = (*r_p_l)[i];
                    t_p_l2 = (*t_p_l)[i];

                    if (d->options & XRTM_OPTION_VECTOR) {
                         r_m_l2 = (*r_m_l)[i];
                         t_m_l2 = (*t_m_l)[i];
                    }
               }

               if (get_local_r_t_u_w(d, i_four, i, *r_p+i, *t_p+i, &r_m2, &t_m2, &r_p_l2, &t_p_l2, &r_m_l2, &t_m_l2, 0, save_tree, work)) {
                    fprintf(stderr, "ERROR: get_local_r_t_u_w()\n");
                    return -1;
               }
          }

          if (d->options & XRTM_OPTION_REVERSE_DERIVS)
               save_tree_decode_i(&save_tree);
     }

     return 0;
}



/*******************************************************************************
 *
 ******************************************************************************/
static void call_solver
     (xrtm_data *d, int solver, int i_four,
      double planck0, double planck1, double *planck0_l, double *planck1_l,
      double omega, double *omega_l, double ltau, double *ltau_l,
      double as_0, double *as_0_l, double atran, double *atran_l,
      double *P_q0_mm, double *P_q0_pm,
      double  **r_p, double  **t_p, double  **r_m, double  **t_m,
      double  **R_p, double  **T_p, double  **R_m, double  **T_m,
      double  *S_p,  double  *S_m, double  *Sl_p,  double  *Sl_m,
      double **P_q0_mm_l, double **P_q0_pm_l,
      double ***r_p_l, double ***t_p_l, double ***r_m_l, double ***t_m_l,
      double ***R_p_l, double ***T_p_l, double ***R_m_l, double ***T_m_l,
      double **S_p_l,  double **S_m_l, double **Sl_p_l,  double **Sl_m_l,
      uchar *derivs_layers, uchar *derivs_beam, uchar *derivs_thermal,
      uchar *derivs_sources, save_tree_data save_tree, work_data work) {

     int i;

     int greens;
     int solar;
     int thermal;
     int vector;

     int symmetric;
     int eigen_solver_real;
     int eigen_solver_complex;
     int check_condition;

     double **r_x;
     double **t_x;
     double ***r_x_l;
     double ***t_x_l;

     greens  = d->options & XRTM_OPTION_PART_SOL_GREENS;
     solar   = d->options & XRTM_OPTION_SOURCE_SOLAR;
     thermal = d->options & XRTM_OPTION_SOURCE_THERMAL && i_four == 0;
     vector  = d->options & XRTM_OPTION_VECTOR;

     symmetric            = d->misc_input.use_symmetric_form;
     eigen_solver_real    = d->misc_input.eigen_solver_gen_real;
     eigen_solver_complex = d->misc_input.eigen_solver_gen_complex;
     check_condition      = d->misc_input.use_pade_check_condition;

     if (! (d->options & XRTM_OPTION_VECTOR)) {
          r_x   = r_p;
          t_x   = t_p;
          r_x_l = r_p_l;
          t_x_l = t_p_l;
     }
     else {
          r_x   = r_m;
          t_x   = t_m;
          r_x_l = r_m_l;
          t_x_l = t_m_l;
     }

     if (solver & XRTM_SOLVER_DOUB_ADD) {
          rtm_doub_rts (d->n_quad_x, d->n_stokes, d->n_derivs, d->F_0, d->qx_v, d->qw_v, planck0, planck1, planck0_l, planck1_l, omega, omega_l, ltau, ltau_l, as_0, as_0_l, P_q0_mm, P_q0_pm, r_p, t_p, r_x, t_x, R_p, T_p, R_m, T_m, S_p, S_m, Sl_p, Sl_m, P_q0_mm_l, P_q0_pm_l, r_p_l, t_p_l, r_x_l, t_x_l, R_p_l, T_p_l, R_m_l, T_m_l, S_p_l, S_m_l, Sl_p_l, Sl_m_l, d->doub_d_tau, symmetric, solar, thermal, vector, derivs_layers, derivs_beam, derivs_thermal, save_tree, work);
     }
     else
     if (solver & XRTM_SOLVER_EIG_ADD) {
          rtm_eig_rts  (d->n_quad_x, d->n_stokes, d->n_derivs, d->F_0, d->qx_v, d->qw_v, planck0, planck1, planck0_l, planck1_l, omega, omega_l, ltau, ltau_l, as_0, as_0_l, atran, atran_l, P_q0_mm, P_q0_pm, r_p, t_p, r_x, t_x, R_p, T_p, R_m, T_m, S_p, S_m, Sl_p, Sl_m, P_q0_mm_l, P_q0_pm_l, r_p_l, t_p_l, r_x_l, t_x_l, R_p_l, T_p_l, R_m_l, T_m_l, S_p_l, S_m_l, Sl_p_l, Sl_m_l, greens, symmetric, solar, thermal, vector, eigen_solver_real, eigen_solver_complex, derivs_layers, derivs_beam, derivs_thermal, save_tree, work);
     }
     else
     if (solver & XRTM_SOLVER_PADE_ADD) {
          rtm_pade_rts2(d->n_quad_x, d->n_stokes, d->n_derivs, d->F_0, d->qx_v, d->qw_v, d->umus, d->n_umus, planck0, planck1, planck0_l, planck1_l, omega, omega_l, ltau, ltau_l, as_0, as_0_l, atran, atran_l, P_q0_mm, P_q0_pm, r_p, t_p, r_x, t_x, R_p, T_p, R_m, T_m, S_p, S_m, Sl_p, Sl_m, P_q0_mm_l, P_q0_pm_l, r_p_l, t_p_l, r_x_l, t_x_l, R_p_l, T_p_l, R_m_l, T_m_l, S_p_l, S_m_l, Sl_p_l, Sl_m_l, d->pade_s, d->pade_r, check_condition, symmetric, solar, thermal, vector, 0,                                 &d->misc_output.pade_condition, derivs_layers, derivs_beam, derivs_thermal, derivs_sources, save_tree, work);
/*
          rtm_pade_rts2(d->n_quad_x, d->n_stokes, d->n_derivs, d->F_0, d->qx_v, d->qw_v, d->umus, d->n_umus, planck0, planck1, planck0_l, planck1_l, omega, omega_l, ltau, ltau_l, as_0, as_0_l, atran, atran_l, P_q0_mm, P_q0_pm, r_p, t_p, r_x, t_x, R_p, T_p, R_m, T_m, S_p, S_m, Sl_p, Sl_m, P_q0_mm_l, P_q0_pm_l, r_p_l, t_p_l, r_x_l, t_x_l, R_p_l, T_p_l, R_m_l, T_m_l, S_p_l, S_m_l, Sl_p_l, Sl_m_l, d->pade_s, d->pade_r, check_condition, symmetric, solar, thermal, vector, d->misc_input.use_pade_gamma_init, &d->misc_output.pade_condition, derivs_layers, derivs_beam, derivs_thermal, derivs_sources, save_tree, work);
*/
     }
#ifdef DEBUG
     else {
         fprintf(stderr, "ERROR: call_solver(): end of if / else if\n");
         exit(1);
     }
#endif

     if (! (d->options & XRTM_OPTION_VECTOR)) {
          dmat_copy(R_m, R_p, d->n_quad_v_x, d->n_quad_v_x);
          dmat_copy(T_m, T_p, d->n_quad_v_x, d->n_quad_v_x);

          if (d->options & XRTM_OPTION_CALC_DERIVS && flags_or(derivs_layers, d->n_derivs)) {
               for (i = 0; i < d->n_derivs; ++i) {
                    if (! derivs_layers[i])
                         continue;

                    dmat_copy(R_m_l[i], R_p_l[i], d->n_quad_v_x, d->n_quad_v_x);
                    dmat_copy(T_m_l[i], T_p_l[i], d->n_quad_v_x, d->n_quad_v_x);
               }
          }
     }
}



/*******************************************************************************
 *
 ******************************************************************************/
static int layer_R_T_S_U_W_V_deps(xrtm_data *d) {

     int deps = opt_props_deps(d) | beam_params_deps(d) | phase_vecs_deps_q0(d) |
                local_r_t_u_w_deps(d);

     if (d->solvers & XRTM_SOLVER_PADE_ADD)
          deps |= DEP_MASK_PADE_PARAMS;

     return deps;
}


static int layer_R_T_S_U_W_V_resolve(xrtm_data *d, int i_four, int i_layer) {

     return resolve_deps1(d, i_four, i_layer, d->dep_flags_layer_R_T_S_U_W_V, layer_R_T_S_U_W_V_deps(d));
}



static int get_layer_R_T_S_U_W_V
     (xrtm_data *d, int i_four, int i_layer, int solver,
      double ***R_p, double ***T_p, double ***R_m, double ***T_m,
      double **S_p, double **S_m, double **Sl_p, double **Sl_m,
      double ****R_p_l, double ****T_p_l, double ****R_m_l, double ****T_m_l,
      double ***S_p_l, double ***S_m_l, double ***Sl_p_l, double ***Sl_m_l,
      int flag, save_tree_data save_tree, work_data *work) {

     uchar *derivs_layers;
     uchar *derivs_beam;
     uchar *derivs_thermal;
     uchar *derivs_sources;

     int i;

     int flag2;

     int deps;

     double b0;
     double b1;

     double *b0_l;
     double *b1_l;

     double *ltau_l;
     double *omega_l;
     double *as_0_l;
     double *atran_l;

     double *P_q0_mm;
     double *P_q0_pm;

     double **P_q0_mm_l;
     double **P_q0_pm_l;

     double **r_p;
     double **t_p;
     double **r_m;
     double **t_m;

     double ***r_p_l;
     double ***t_p_l;
     double ***r_m_l;
     double ***t_m_l;

     forward_save_get_layer_R_T_S_U_W_V_data *save;

     work_data work2;


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (d->options & XRTM_OPTION_REVERSE_DERIVS) {
          save_tree_encode_s(&save_tree, "get_layer_R_T_S_U_W_V");

          if (save_tree_retrieve_data(&save_tree, forward_save_get_layer_R_T_S_U_W_V_data, &save))
               forward_save_get_layer_R_T_S_U_W_V_alloc(save, d->options, d->n_quad_v_x);
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (d->options & XRTM_OPTION_SAVE_LAYER_R_T_S) {
          deps = layer_R_T_S_U_W_V_deps(d);

          flag2 = resolve_deps1(d, i_four, i_layer,
                                d->dep_flags_layer_R_T_S_U_W_V, deps);

          *R_p  = d->Rl_p[i_four][i_layer];
          *T_p  = d->Tl_p[i_four][i_layer];
          *R_m  = d->Rl_m[i_four][i_layer];
          *T_m  = d->Tl_m[i_four][i_layer];

          *S_p  = d->Sl_p[i_four][i_layer];
          *S_m  = d->Sl_m[i_four][i_layer];

          *Sl_p = d->Sll_p[i_four][i_layer];
          *Sl_m = d->Sll_m[i_four][i_layer];

          if (d->options & XRTM_OPTION_CALC_DERIVS && flags_or(d->derivs.layers[i_layer], d->n_derivs)) {
               *R_p_l  = d->Rl_p_l[i_four][i_layer];
               *T_p_l  = d->Tl_p_l[i_four][i_layer];
               *R_m_l  = d->Rl_m_l[i_four][i_layer];
               *T_m_l  = d->Tl_m_l[i_four][i_layer];
          }
          if (d->options & XRTM_OPTION_SOURCE_SOLAR && d->options & XRTM_OPTION_CALC_DERIVS && flags_or(d->derivs.beam[i_layer], d->n_derivs)) {
               *S_p_l  = d->Sl_p_l[i_four][i_layer];
               *S_m_l  = d->Sl_m_l[i_four][i_layer];
          }

          if (d->options & XRTM_OPTION_SOURCE_THERMAL && d->options & XRTM_OPTION_CALC_DERIVS && flags_or(d->derivs.thermal[i_layer], d->n_derivs)) {
               *Sl_p_l = d->Sll_p_l[i_four][i_layer];
               *Sl_m_l = d->Sll_m_l[i_four][i_layer];
          }

          if (! flag2)
               return 0;
     }
     else if (flag) {
          *R_p  = get_work1(work, WORK_DXX);
          *T_p  = get_work1(work, WORK_DXX);
          *R_m  = get_work1(work, WORK_DXX);
          *T_m  = get_work1(work, WORK_DXX);

          *S_p  = get_work1(work, WORK_DX);
          *S_m  = get_work1(work, WORK_DX);

          *Sl_p = get_work1(work, WORK_DX);
          *Sl_m = get_work1(work, WORK_DX);

          if (d->options & XRTM_OPTION_CALC_DERIVS && flags_or(d->derivs.layers[i_layer], d->n_derivs)) {
               *R_p_l  = get_work2(work, WORK_DXX, WORK_DERIVS_V, d->derivs.layers[i_layer]);
               *T_p_l  = get_work2(work, WORK_DXX, WORK_DERIVS_V, d->derivs.layers[i_layer]);
               *R_m_l  = get_work2(work, WORK_DXX, WORK_DERIVS_V, d->derivs.layers[i_layer]);
               *T_m_l  = get_work2(work, WORK_DXX, WORK_DERIVS_V, d->derivs.layers[i_layer]);
          }
          if (d->options & XRTM_OPTION_SOURCE_SOLAR && d->options & XRTM_OPTION_CALC_DERIVS && flags_or(d->derivs.beam[i_layer], d->n_derivs)) {
               *S_p_l  = get_work2(work, WORK_DX, WORK_DERIVS_V, d->derivs.beam[i_layer]);
               *S_m_l  = get_work2(work, WORK_DX, WORK_DERIVS_V, d->derivs.beam[i_layer]);
          }

          if (d->options & XRTM_OPTION_SOURCE_THERMAL && d->options & XRTM_OPTION_CALC_DERIVS && flags_or(d->derivs.thermal[i_layer], d->n_derivs)) {
               *Sl_p_l = get_work2(work, WORK_DX, WORK_DERIVS_V, d->derivs.thermal[i_layer]);
               *Sl_m_l = get_work2(work, WORK_DX, WORK_DERIVS_V, d->derivs.thermal[i_layer]);
          }
     }


     work2 = *work;


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (d->options & XRTM_OPTION_CALC_DERIVS) {
          derivs_layers  = d->derivs.layers[i_layer];
          if (d->options & XRTM_OPTION_SOURCE_SOLAR)
               derivs_beam    = d->derivs.beam[i_layer];
          if (d->options & XRTM_OPTION_SOURCE_THERMAL)
               derivs_thermal = d->derivs.thermal[i_layer];
          derivs_sources = d->derivs.sources[i_layer];
     }

     if (d->options & XRTM_OPTION_CALC_DERIVS && flags_or(d->derivs.layers[i_layer], d->n_derivs)) {
          ltau_l  = d->ltau_l [i_layer];
          omega_l = d->omega_l[i_layer];
     }

     if (d->options & XRTM_OPTION_SOURCE_SOLAR && d->options & XRTM_OPTION_CALC_DERIVS && flags_or(d->derivs.beam[i_layer], d->n_derivs)) {
          as_0_l  = d->as_0_l [i_layer];
          atran_l = d->atran_l[i_layer];
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     update_opt_props(d, i_layer, work2);


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (d->omega[i_layer] == 0. || i_four >= d->n_coef_layer[i_layer])
          no_scatter_r_t_s(d->n_quad_v_x, d->n_derivs, d->qx_v, d->qw_v, d->ltau[i_layer], ltau_l, *R_p, *T_p, *R_m, *T_m, *S_p, *S_m, *Sl_p, *Sl_m, *R_p_l, *T_p_l, *R_m_l, *T_m_l, *S_p_l, *S_m_l, *Sl_p_l, *Sl_m_l, d->options & XRTM_OPTION_SOURCE_SOLAR, d->options & XRTM_OPTION_SOURCE_THERMAL && i_four == 0, derivs_layers, derivs_beam, derivs_thermal);


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     else {
          update_beam_params(d, i_layer, work2);
if (d->options & XRTM_OPTION_SOURCE_SOLAR) {
          if (get_phase_vecs_q0(d, i_four, i_layer, &P_q0_mm, &P_q0_pm, &P_q0_mm_l, &P_q0_pm_l, 1, &work2)) {
               fprintf(stderr, "ERROR: get_phase_vecs_q0()\n");
               return -1;
          }
}
          if (get_local_r_t_u_w(d, i_four, i_layer, &r_p, &t_p, &r_m, &t_m, &r_p_l, &t_p_l, &r_m_l, &t_m_l, 1, save_tree, &work2)) {
               fprintf(stderr, "ERROR: get_local_r_t_u_w()\n");
               return -1;
          }

          if (d->options & XRTM_OPTION_SOURCE_THERMAL) {
               b0 = d->levels_b[i_layer];
               b1 = d->levels_b[i_layer+1];

               if (d->options & XRTM_OPTION_CALC_DERIVS) {
                    b0_l = d->levels_b_l[i_layer];
                    b1_l = d->levels_b_l[i_layer+1];
               }
          }
          call_solver(d, solver, i_four, b0, b1, b0_l, b1_l, d->omega[i_layer], omega_l, d->ltau [i_layer], ltau_l, d->as_0[i_layer], as_0_l, d->atran[i_layer], atran_l, P_q0_mm, P_q0_pm, r_p, t_p, r_m, t_m, *R_p, *T_p, *R_m, *T_m, *S_p, *S_m, *Sl_p, *Sl_m, P_q0_mm_l, P_q0_pm_l, r_p_l, t_p_l, r_m_l, t_m_l, *R_p_l, *T_p_l, *R_m_l, *T_m_l, *S_p_l, *S_m_l, *Sl_p_l, *Sl_m_l, derivs_layers, derivs_beam, derivs_thermal, derivs_sources, save_tree, work2);


          if (d->options & XRTM_OPTION_REVERSE_DERIVS) {
               dvec_copy(save->P_q0_mm, P_q0_mm, d->n_quad_v_x);
               dvec_copy(save->P_q0_pm, P_q0_pm, d->n_quad_v_x);

               dmat_copy(save->r_p, r_p, d->n_quad_v_x, d->n_quad_v_x);
               dmat_copy(save->t_p, t_p, d->n_quad_v_x, d->n_quad_v_x);

               if (d->options & XRTM_OPTION_VECTOR) {
                    dmat_copy(save->r_m, r_m, d->n_quad_v_x, d->n_quad_v_x);
                    dmat_copy(save->t_m, t_m, d->n_quad_v_x, d->n_quad_v_x);
               }
          }
     }
#ifdef USE_AD_FOR_TL_CALC_GET_LAYER_R_T_S_U_W_V
     get_layer_R_T_S_U_W_V_tl_with_ad(d, i_four, i_layer, solver, *R_p, *T_p, *R_m, *T_m, *S_p, *S_m, *R_p_l, *T_p_l, *R_m_l, *T_m_l, *S_p_l, *S_m_l, flag, save_tree, work2);
#endif
     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
if (0) {
     if (i_four == 0) {
          if (check_R_and_T_norm(d->n_quad, d->n_stokes, *R_p, *T_p) != 0) {
               fprintf(stderr, "ERROR: check_R_and_T_norm()\n");
               return -1;
          }

          for (i = 0; i < d->n_derivs; ++i) {
               if (d->derivs.layers[i_layer][i]) {
                    if (check_R_and_T_norm_l(d->n_quad, d->n_stokes, (*R_p_l)[i], (*T_p_l)[i]) != 0) {
                         fprintf(stderr, "ERROR: check_R_and_T_norm_l()\n");
                         return -1;
                    }
               }
          }
     }
}

     return 0;
}



/*******************************************************************************
 *
 ******************************************************************************/
static void source_vectors_scale_add
     (xrtm_data *d, int i_four, int i_layer, int i_deriv, int n_derivs,
      double  *S_p1, double  *S_m1, double  *Sl_p1, double  *Sl_m1,
      double  *S_p2, double  *S_m2,
      double **S_p_l1, double **S_m_l1, double **Sl_p_l1, double **Sl_m_l1,
      double **S_p_l2, double **S_m_l2, work_data work) {

     uchar *derivs_layers;
     uchar *derivs_beam;
     uchar *derivs_thermal;
     uchar *derivs_sources;

     double *btran_l;

     update_beam_params(d, i_layer, work);

     if (d->options & XRTM_OPTION_CALC_DERIVS) {
          btran_l        = d->btran_l       [i_layer]+i_deriv;
          derivs_layers  = d->derivs.layers [i_layer]+i_deriv;
          if (d->options & XRTM_OPTION_SOURCE_SOLAR)
               derivs_beam    = d->derivs.beam[i_layer]+i_deriv;
          if (d->options & XRTM_OPTION_SOURCE_THERMAL)
               derivs_thermal = d->derivs.thermal[i_layer]+i_deriv;
          derivs_sources = d->derivs.sources[i_layer]+i_deriv;
     }

     scale_add_source_vectors(d->n_quad_v_x, n_derivs, d->btran[i_layer], btran_l, S_p1, S_m1, Sl_p1, Sl_m1, S_p2, S_m2, S_p_l1, S_m_l1, Sl_p_l1, Sl_m_l1, S_p_l2, S_m_l2, d->options & XRTM_OPTION_SOURCE_SOLAR, d->options & XRTM_OPTION_SOURCE_THERMAL && i_four == 0, derivs_layers, derivs_beam, derivs_thermal, derivs_sources, work);
}



/*******************************************************************************
 *
 ******************************************************************************/
static int build_layer_stacks(xrtm_data *d, int i_four, int solver, save_tree_data save_tree, work_data work) {

     uchar *derivs_stacks;

     int i;

     double *v1;
     double *v2;
     double *v3;
     double *v4;

     double *S_p1;
     double *S_m1;
     double **S_p_l1;
     double **S_m_l1;

     double *S_p2;
     double *S_m2;
     double **S_p_l2;
     double **S_m_l2;

     stack_data *s;

     stack_data *s0;
     stack_data *s1;
     stack_data *s2;

     work_data work2;


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     v1 = get_work1(&work, WORK_DX);
     v2 = get_work1(&work, WORK_DX);
     v3 = get_work1(&work, WORK_DX);
     v4 = get_work1(&work, WORK_DX);


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (! (derivs_stacks = flags_alloc(d->n_derivs))) {
          fprintf(stderr, "ERROR: flags_alloc()\n");
          return -1;
     }

     for (i = 0; i < d->n_stacks; ++i) {

          s = &d->stack_chain[i];

          if (i < d->n_layers) {
               if (layer_R_T_S_U_W_V_resolve(d, i_four, i) ||
                   d->misc_input.use_rebuild_stacks) {

                    if (get_layer_R_T_S_U_W_V(d, i_four, i, solver,
                                              &s->R_p[i_four], &s->T_p[i_four], &s->R_m[i_four], &s->T_m[i_four],
                                              &s->S_p[i_four], &s->S_m[i_four], &s->Sl_p[i_four], &s->Sl_m[i_four],
                                              &s->R_p_l[i_four], &s->T_p_l[i_four], &s->R_m_l[i_four], &s->T_m_l[i_four],
                                              &s->S_p_l[i_four], &s->S_m_l[i_four], &s->Sl_p_l[i_four], &s->Sl_m_l[i_four],
                                              0, save_tree, &work)) {

                         fprintf(stderr, "ERROR: get_layer_R_T_S_U_W_V()\n");
                         return -1;
                    }

                    s->touched = 1;
               }
               else
                    s->touched = 0;
          }
          else {
               s1 = s->p1;
               s2 = s->p2;
               s0 = s;

               s->touched = s1->touched || s2->touched || s->new_[i_four] || d->misc_input.use_rebuild_stacks;

               derivs_merge_s(d->n_derivs, d->derivs.stacks[s->i1], d->derivs.stacks[s->i2], derivs_stacks);

               work2 = work;

               if (s1->i1 == s1->i2) {
                    S_p1 = v1;
                    S_m1 = v2;

                    if (d->options & XRTM_OPTION_CALC_DERIVS && flags_or(d->derivs.sources[s1->i1], d->n_derivs)) {
                         S_p_l1 = get_work2(&work2, WORK_DX, WORK_DERIVS_V, d->derivs.sources[s1->i1]);
                         S_m_l1 = get_work2(&work2, WORK_DX, WORK_DERIVS_V, d->derivs.sources[s1->i1]);
                    }

                    source_vectors_scale_add(d, i_four, s1->i1, 0, d->n_derivs, s1->S_p[i_four], s1->S_m[i_four], s1->Sl_p[i_four], s1->Sl_m[i_four], S_p1, S_m1, s1->S_p_l[i_four], s1->S_m_l[i_four], s1->Sl_p_l[i_four], s1->Sl_m_l[i_four], S_p_l1, S_m_l1, work2);
               }
               else {
                    S_p1 = s1->S_p[i_four];     S_p_l1 = s1->S_p_l[i_four];
                    S_m1 = s1->S_m[i_four];     S_m_l1 = s1->S_m_l[i_four];
               }

               if (s2->i1 == s2->i2) {
                    S_p2 = v3;
                    S_m2 = v4;

                    if (d->options & XRTM_OPTION_CALC_DERIVS && flags_or(d->derivs.sources[s2->i1], d->n_derivs)) {
                         S_p_l2 = get_work2(&work2, WORK_DX, WORK_DERIVS_V, d->derivs.sources[s2->i1]);
                         S_m_l2 = get_work2(&work2, WORK_DX, WORK_DERIVS_V, d->derivs.sources[s2->i1]);
                    }

                    source_vectors_scale_add(d, i_four, s2->i1, 0, d->n_derivs, s2->S_p[i_four], s2->S_m[i_four], s2->Sl_p[i_four], s2->Sl_m[i_four], S_p2, S_m2, s2->S_p_l[i_four], s2->S_m_l[i_four], s2->Sl_p_l[i_four], s2->Sl_m_l[i_four], S_p_l2, S_m_l2, work2);
               }
               else {
                    S_p2 = s2->S_p[i_four];     S_p_l2 = s2->S_p_l[i_four];
                    S_m2 = s2->S_m[i_four];     S_m_l2 = s2->S_m_l[i_four];
               }

               if (! s->new_[i_four] && ! s->touched)
                    source_add_all2(s1->R_m[i_four], s1->T_m[i_four], S_m1, s1->R_p[i_four], s1->T_p[i_four], S_p1, s2->R_m[i_four], s2->T_m[i_four], S_m2, s2->R_p[i_four], s2->T_p[i_four], S_p2, s0->R_m[i_four], s0->T_m[i_four], s0->S_m[i_four], s0->R_p[i_four], s0->T_p[i_four], s0->S_p[i_four], NULL, NULL, S_m_l1, NULL, NULL, S_p_l1, NULL, NULL, S_m_l2, NULL, NULL, S_p_l2, NULL, NULL, s0->S_m_l[i_four], NULL, NULL, s0->S_p_l[i_four], d->n_quad_v_x, d->n_derivs, 1., NULL, derivs_stacks, 1, s->touched, work2, &s->d[i_four]);
               else
                    layer_add_all2 (s1->R_m[i_four], s1->T_m[i_four], S_m1, s1->R_p[i_four], s1->T_p[i_four], S_p1, s2->R_m[i_four], s2->T_m[i_four], S_m2, s2->R_p[i_four], s2->T_p[i_four], S_p2, s0->R_m[i_four], s0->T_m[i_four], s0->S_m[i_four], s0->R_p[i_four], s0->T_p[i_four], s0->S_p[i_four], NULL, NULL, S_m_l1, NULL, NULL, S_p_l1, NULL, NULL, S_m_l2, NULL, NULL, S_p_l2, NULL, NULL, s0->S_m_l[i_four], NULL, NULL, s0->S_p_l[i_four], d->n_quad_v_x, d->n_derivs, 1., NULL, derivs_stacks, 1, s->touched, work2, &s->d[i_four]);

               s->new_[i_four] = 0;
          }
     }


     flags_free(derivs_stacks);


     return 0;
}



/*******************************************************************************
 *
 ******************************************************************************/
static int total_R_T_S_U_W_V_deps(xrtm_data *d) {

     return layer_R_T_S_U_W_V_deps(d);
}



/*******************************************************************************
 *
 ******************************************************************************/
static int get_total_BOA_R_S_U_V
     (xrtm_data *d, int i_four, int solver,
      double  *S_t, double  **S_t_l,
      double ***R_p, double **S_m, double ****R_p_l, double ***S_m_l,
      int surface, int flag, save_tree_data save_tree, work_data *work) {

     uchar *derivs_layers;
     uchar *derivs_thermal;
     uchar *derivs_beam;
     uchar *derivs_sources;
     uchar *derivs_adding_down;

     int i;

     int deps;

     int flag2;

     double *trans_l;

     double *btran_l;

     double *S1_m;

     double *S2_p;
     double *S2_m;
     double *Sl2_p;
     double *Sl2_m;

     double *S3_p;
     double *S3_m;

     double **R1_p;

     double **R2_p;
     double **T2_p;
     double **R2_m;
     double **T2_m;

     double **S1_m_l;

     double **S2_p_l;
     double **S2_m_l;
     double **Sl2_p_l;
     double **Sl2_m_l;

     double **S3_p_l;
     double **S3_m_l;

     double ***R1_p_l;

     double ***R2_p_l;
     double ***T2_p_l;
     double ***R2_m_l;
     double ***T2_m_l;

     work_data work2;
     work_data work3;


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (d->options & XRTM_OPTION_SAVE_TOTAL_R_T_S) {
          deps = total_R_T_S_U_W_V_deps(d);

          flag2 = resolve_deps0(d, i_four, d->dep_flags_total_R_T_S_U_W_V, deps);

          *R_p = d->Ra_p[i_four];
          *S_m = d->Sa_m[i_four];

          if (d->options & XRTM_OPTION_CALC_DERIVS) {
               *R_p_l = d->Ra_p_l[i_four];
               *S_m_l = d->Sa_m_l[i_four];
          }

          if (! flag2)
               return 0;
     }
     else if (flag) {
          *R_p = get_work1(work, WORK_DXX);

          *S_m = get_work1(work, WORK_DX);

          if (d->options & XRTM_OPTION_CALC_DERIVS && flags_or(d->derivs.total[0], d->n_derivs)) {
               *R_p_l = get_work2(work, WORK_DXX, WORK_DERIVS_V, d->derivs.total[0]);

               *S_m_l = get_work2(work, WORK_DX,  WORK_DERIVS_V, d->derivs.total[0]);
          }
     }


     work2 = *work;


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     R1_p = *R_p;
     S1_m = *S_m;

     if (d->options & XRTM_OPTION_CALC_DERIVS) {
          R1_p_l = *R_p_l;
          S1_m_l = *S_m_l;
     }

     S3_p = get_work1(&work2, WORK_DX);
     S3_m = get_work1(&work2, WORK_DX);


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (d->options & XRTM_OPTION_CALC_DERIVS) {
          trans_l = get_work1(&work2, WORK_DDERIVS);

          for (i = 0; i < d->n_derivs; ++i)
               trans_l[i] = 0.;
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
/*
     dvec_copy(S1_m, S_t, d->n_quad_v_x);

     if (d->options & XRTM_OPTION_CALC_DERIVS && flags_or(d->derivs.sources[d->n_layers], d->n_derivs)) {
          for (i = 0; i < d->n_derivs; ++i) {
               if (d->derivs.sources[d->n_layers][i])
                    dvec_copy(S1_m_l[i], S_t_l[i], d->n_quad_v_x);
          }
     }
*/
     for (i = 0; i < d->n_layers; ++i) {
          work3 = work2;

          if (d->options & XRTM_OPTION_CALC_DERIVS && flags_or(d->derivs.sources[i], d->n_derivs)) {
               S3_p_l = get_work2(&work3, WORK_DX, WORK_DERIVS_V, d->derivs.sources[i]);
               S3_m_l = get_work2(&work3, WORK_DX, WORK_DERIVS_V, d->derivs.sources[i]);
          }

          if (get_layer_R_T_S_U_W_V(d, i_four, i, solver, &R2_p, &T2_p, &R2_m, &T2_m, &S2_p, &S2_m, &Sl2_p, &Sl2_m, &R2_p_l, &T2_p_l, &R2_m_l, &T2_m_l, &S2_p_l, &S2_m_l, &Sl2_p_l, &Sl2_m_l, 1, save_tree, &work3)) {
               fprintf(stderr, "ERROR: get_layer_R_T_S_U_W_V()\n");
               return -1;
          }

          if (d->options & XRTM_OPTION_CALC_DERIVS) {
               derivs_layers      = d->derivs.layers[i];
               if (d->options & XRTM_OPTION_SOURCE_SOLAR)
                    derivs_beam    = d->derivs.beam[i];
               if (d->options & XRTM_OPTION_SOURCE_THERMAL)
                    derivs_thermal = d->derivs.thermal[i];
               derivs_sources     = d->derivs.sources[i];
               derivs_adding_down = d->derivs.adding_down[i];
          }

          update_beam_params(d, i, work3);

          if (d->options & XRTM_OPTION_CALC_DERIVS)
               btran_l = d->btran_l[i];

          scale_add_source_vectors(d->n_quad_v_x, d->n_derivs, d->btran[i], btran_l, S2_p, S2_m, Sl2_p, Sl2_m, S3_p, S3_m, S2_p_l, S2_m_l, Sl2_p_l, Sl2_m_l, S3_p_l, S3_m_l, d->options & XRTM_OPTION_SOURCE_SOLAR, d->options & XRTM_OPTION_SOURCE_THERMAL && i_four == 0, derivs_layers, derivs_beam, derivs_thermal, derivs_sources, work3);

          if (i == 0)
               layer_copy_ref_l(R1_p, S1_m, R2_p, S3_m, R1_p_l, S1_m_l, R2_p_l, S3_m_l, d->n_quad_v_x, d->n_derivs, derivs_layers, derivs_sources);
          else
               layer_add_ref   (R2_p, T2_p, S3_p, R2_m, T2_m, S3_m, R1_p, S1_m, R1_p, S1_m, R2_p_l, T2_p_l, S3_p_l, R2_m_l, T2_m_l, S3_m_l, R1_p_l, S1_m_l, R1_p_l, S1_m_l, d->n_quad_v_x, d->n_derivs, 1., trans_l, derivs_adding_down, 1, 1, save_tree, work3);
     }


     return 0;
}



/*******************************************************************************
 *
 ******************************************************************************/
static int get_total_TOA_R_S_U_V
     (xrtm_data *d, int i_four, int solver,
      double  **Rs_qq, double  *S_s, double ***Rs_qq_l, double  **S_s_l,
      double ***R_m, double **S_p, double ****R_m_l, double ***S_p_l,
      int surface, int flag, save_tree_data save_tree, work_data *work) {

     uchar **derivs_layers;
     uchar **derivs_thermal;
     uchar **derivs_beam;
     uchar **derivs_sources;
     uchar **derivs_total;
     uchar **derivs_adding_up;

     uchar  *derivs_layers2;
     uchar  *derivs_thermal2;
     uchar  *derivs_beam2;
     uchar  *derivs_sources2;
     uchar  *derivs_adding_up2;

     int i;

     double *trans_l;

     double *btran_l;

     double *S1_p;

     double *S2_p;
     double *S2_m;
     double *Sl2_p;
     double *Sl2_m;

     double *S3_p;
     double *S3_m;

     double **R1_m;

     double **R2_p;
     double **T2_p;
     double **R2_m;
     double **T2_m;

     double **S1_p_l;

     double **S2_p_l;
     double **S2_m_l;
     double **Sl2_p_l;
     double **Sl2_m_l;

     double **S3_p_l;
     double **S3_m_l;

     double ***R1_m_l;

     double ***R2_p_l;
     double ***T2_p_l;
     double ***R2_m_l;
     double ***T2_m_l;

     forward_save_get_total_TOA_R_S_U_V_data *save;

     work_data work2;
     work_data work3;


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (d->options & XRTM_OPTION_REVERSE_DERIVS) {
          save_tree_encode_s(&save_tree, "get_total_TOA_R_S_U_V");

          if (save_tree_retrieve_data(&save_tree, forward_save_get_total_TOA_R_S_U_V_data, &save))
               forward_save_get_total_TOA_R_S_U_V_alloc(save, d->n_quad_v_x, d->n_layers);
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (flag) {
          *R_m = get_work1(work, WORK_DXX);

          *S_p = get_work1(work, WORK_DX);

          if (d->options & XRTM_OPTION_CALC_DERIVS && flags_or(d->derivs.total[0], d->n_derivs)) {
               *R_m_l = get_work2(work, WORK_DXX, WORK_DERIVS_V, d->derivs.total[0]);

               *S_p_l = get_work2(work, WORK_DX,  WORK_DERIVS_V, d->derivs.total[0]);
          }
     }


     work2 = *work;


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     R1_m = *R_m;
     S1_p = *S_p;

     if (d->options & XRTM_OPTION_CALC_DERIVS) {
          R1_m_l = *R_m_l;
          S1_p_l = *S_p_l;
     }

     S3_p = get_work1(&work2, WORK_DX);
     S3_m = get_work1(&work2, WORK_DX);


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (i_four == 0) {
          derivs_layers    = d->derivs.layers;
          derivs_beam      = d->derivs.beam;
          derivs_thermal   = d->derivs.thermal;
          derivs_sources   = d->derivs.sources;
          derivs_total     = d->derivs.total;
          derivs_adding_up = d->derivs.adding_up;
     }
     else {
          derivs_layers    = d->derivs.layers_m;
          derivs_beam      = d->derivs.beam_m;
          derivs_thermal   = d->derivs.thermal_m;
          derivs_sources   = d->derivs.sources_m;
          derivs_total     = d->derivs.total_m;
          derivs_adding_up = d->derivs.adding_up_m;
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (d->options & XRTM_OPTION_CALC_DERIVS) {
          trans_l = get_work1(&work2, WORK_DDERIVS);

          for (i = 0; i < d->n_derivs; ++i)
               trans_l[i] = 0.;
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (surface) {
          dmat_copy(R1_m, Rs_qq, d->n_quad_v_x, d->n_quad_v_x);
          dvec_copy(S1_p, S_s, d->n_quad_v_x);

          if (d->options & XRTM_OPTION_CALC_DERIVS && flags_or(derivs_total[0], d->n_derivs)) {
               for (i = 0; i < d->n_derivs; ++i) {
                    if (derivs_layers[d->n_layers][i])
                         dmat_copy(R1_m_l[i], Rs_qq_l[i], d->n_quad_v_x, d->n_quad_v_x);

                    if (derivs_sources[d->n_layers][i])
                         dvec_copy(S1_p_l[i], S_s_l[i], d->n_quad_v_x);
               }
          }
     }

     for (i = d->n_layers - 1; i >= 0; --i) {
          work3 = work2;

          if (d->options & XRTM_OPTION_REVERSE_DERIVS)
               save_tree_recode_i(&save_tree, i, i == d->n_layers - 1);

          if (d->options & XRTM_OPTION_CALC_DERIVS && flags_or(d->derivs.sources[i], d->n_derivs)) {
               S3_p_l = get_work2(&work3, WORK_DX, WORK_DERIVS_V, d->derivs.sources[i]);
               S3_m_l = get_work2(&work3, WORK_DX, WORK_DERIVS_V, d->derivs.sources[i]);
          }

          if (get_layer_R_T_S_U_W_V(d, i_four, i, solver, &R2_p, &T2_p, &R2_m, &T2_m, &S2_p, &S2_m, &Sl2_p, &Sl2_m, &R2_p_l, &T2_p_l, &R2_m_l, &T2_m_l, &S2_p_l, &S2_m_l, &Sl2_p_l, &Sl2_m_l, 1, save_tree, &work3)) {
               fprintf(stderr, "ERROR: get_layer_R_T_S_U_W_V()\n");
               return -1;
          }

          if (d->options & XRTM_OPTION_CALC_DERIVS) {
               derivs_layers2    = derivs_layers[i];
               if (d->options & XRTM_OPTION_SOURCE_SOLAR)
                    derivs_beam2 = derivs_beam[i];
               if (d->options & XRTM_OPTION_SOURCE_THERMAL)
                    derivs_thermal2   = derivs_thermal[i];
               derivs_sources2   = derivs_sources[i];
               derivs_adding_up2 = derivs_adding_up[i];
          }

          update_beam_params(d, i, work3);

          if (d->options & XRTM_OPTION_CALC_DERIVS)
               btran_l = d->btran_l[i];

          scale_add_source_vectors(d->n_quad_v_x, d->n_derivs, d->btran[i], btran_l, S2_p, S2_m, Sl2_p, Sl2_m, S3_p, S3_m, S2_p_l, S2_m_l, Sl2_p_l, Sl2_m_l, S3_p_l, S3_m_l, d->options & XRTM_OPTION_SOURCE_SOLAR, d->options & XRTM_OPTION_SOURCE_THERMAL && i_four == 0, derivs_layers2, derivs_beam2, derivs_thermal2, derivs_sources2, work3);

          if (d->options & XRTM_OPTION_REVERSE_DERIVS) {
               layer_copy_ref(save->R1_m[i], save->S1_p[i], R1_m, S1_p, d->n_quad_v_x);

               layer_copy_all(save->R2_m[i], save->T2_m[i], save->S2_m[i], save->R2_p[i], save->T2_p[i], save->S2_p[i], R2_m, T2_m, S2_m, R2_p, T2_p, S2_p, d->n_quad_v_x);

               dvec_copy(save->S3_p[i], S3_p, d->n_quad_v_x);
               dvec_copy(save->S3_m[i], S3_m, d->n_quad_v_x);
          }

          if (! surface && i == d->n_layers - 1)
               layer_copy_ref_l(R1_m, S1_p, R2_m, S3_p, R1_m_l, S1_p_l, R2_m_l, S3_p_l, d->n_quad_v_x, d->n_derivs, derivs_layers2, derivs_sources2);
          else
               layer_add_ref   (R2_m, T2_m, S3_m, R2_p, T2_p, S3_p, R1_m, S1_p, R1_m, S1_p, R2_m_l, T2_m_l, S3_m_l, R2_p_l, T2_p_l, S3_p_l, R1_m_l, S1_p_l, R1_m_l, S1_p_l, d->n_quad_v_x, d->n_derivs, 1., trans_l, derivs_adding_up2, 0, 1, save_tree, work3);
     }

     if (d->options & XRTM_OPTION_REVERSE_DERIVS)
          save_tree_decode_i(&save_tree);

#ifdef USE_AD_FOR_TL_CALC_GET_TOTAL_TOA_R_S_U_V
     get_total_TOA_R_S_U_V_tl_with_ad(d, i_four, solver, Rs_qq, S_s, Rs_qq_l, S_s_l, *R_m, *S_p, *R_m_l, *S_p_l, surface, flag, save_tree, work2);
#endif

     return 0;
}



/*******************************************************************************
 *
 ******************************************************************************/
static int get_total_TOA_R_S_U_V_2
     (xrtm_data *d, int i_four, int solver,
      double  **Rs_qq, double  *S_s, double  ***Rs_qq_l, double  **S_s_l,
      double ***R_m, double **S_p, double ****R_m_l, double ***S_p_l,
      int surface, int flag, save_tree_data save_tree, work_data *work) {

     uchar **derivs_layers;
     uchar **derivs_sources;

     uchar d_h[2];
     uchar d_d[2];
     uchar d_a[2];

     int i;
     int i1;
     int i2;
     int ii;
     int j;

     int flag2;
     int flag3;

     int n_c;
     int n_a;
/*
     int i_stack;
*/
     double *trans_l;

     double *v1;
     double *v2;
     double *v3;
     double *v4;

     double *S_p1;
     double *S_m1;

     double **S_p_l1;
     double **S_m_l1;

     double ***R_p_l1;
     double ***T_p_l1;
     double ***R_m_l1;
     double ***T_m_l1;


     layer_add_aux_data *layer_add_aux;

     stack_data *s;

     stack_data *s_c[2];
     stack_data *s_a[2];

     work_data work2;


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (flag) {
          *R_m = get_work1(work, WORK_DXX);

          *S_p = get_work1(work, WORK_DX);

          if (d->options & XRTM_OPTION_CALC_DERIVS && flags_or(d->derivs.total[0], d->n_derivs)) {
               *R_m_l = get_work2(work, WORK_DXX, WORK_DERIVS_V, d->derivs.total[0]);

               *S_p_l = get_work2(work, WORK_DX,  WORK_DERIVS_V, d->derivs.total[0]);
          }
     }


     work2 = *work;


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (i_four == 0) {
          derivs_layers  = d->derivs.layers;
          derivs_sources = d->derivs.sources;
     }
     else {
          derivs_layers  = d->derivs.layers_m;
          derivs_sources = d->derivs.sources_m;
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (build_layer_stacks(d, i_four, solver, save_tree, work2) < 0) {
          fprintf(stderr, "ERROR: build_layer_stacks()\n");
          return -1;
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     v1 = get_work1(&work2, WORK_DX);
     v2 = get_work1(&work2, WORK_DX);
     v3 = get_work1(&work2, WORK_DX);
     v4 = get_work1(&work2, WORK_DX);


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (d->options & XRTM_OPTION_CALC_DERIVS) {
          trans_l = get_work1(&work2, WORK_DDERIVS);

          for (i = 0; i < d->n_derivs; ++i)
               trans_l[i] = 0.;
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     for (j = 0; j < d->n_derivs; ++j) {

          if (surface) {
               dmat_copy(*R_m, Rs_qq, d->n_quad_v_x, d->n_quad_v_x);
               dvec_copy(*S_p, S_s, d->n_quad_v_x);

               i = d->n_layers;
               if (derivs_layers[i][j])
                    dmat_copy((*R_m_l)[j], Rs_qq_l[j], d->n_quad_v_x, d->n_quad_v_x);

               if (derivs_sources[i][j])
                    dvec_copy((*S_p_l)[j], S_s_l[j], d->n_quad_v_x);
          }

          i1    = d->n_layers - 1;
          flag2 = 0;

          for (i = d->n_layers - 1; i >= 0; --i) {

               if (derivs_layers[i][j]    && ! flag2 && i == d->n_layers - 1 && ! surface) {
                    /* copy linearized layer */
                    n_c    = 1;
                    d_h[0] = derivs_layers[i][j];
                    d_d[0] = derivs_sources[i][j];
                    s_c[0] = &d->stack_chain[i];
                    n_a    = 0;
                    flag2  = 1;
                    i1     = i - 1;
               }
               else
               if (derivs_layers[i][j]    &&  ! flag2 && i == d->n_layers - 1 &&   surface) {
                    /* add linearized layer */
                    n_c    = 0;
                    n_a    = 1;
                    d_a[0] = (uchar) (derivs_layers[d->n_layers][j] ? ADDING_B_B : ADDING_B_S);
                    s_a[0] = &d->stack_chain[i];
                    flag2  = 1;
                    i1     = i - 1;
               }
               else
               if (derivs_layers[i][j]    && ! flag2 && i != d->n_layers - 1 && ! surface) {
                    /* copy unlinearized stack, add linearized layer */
                    i2     = i + 1;
                    n_c    = 1;
                    d_h[0] = 0;
                    d_d[0] = 1;
                    s_c[0] = d->stack_grid[i2][i1];
                    n_a    = 1;
                    d_a[0] = ADDING_B_S;
                    s_a[0] = &d->stack_chain[i];
                    flag2  = 1;
                    i1     = i - 1;
               }
               else
               if (derivs_layers[i][j]    && ! flag2 && i != d->n_layers - 1 &&   surface) {
                    /* add unlinearized stack, add linearized layer */
                    i2     = i + 1;
                    n_c    = 0;
                    n_a    = 2;
                    d_a[0] = (uchar) (derivs_layers[d->n_layers][j] ? ADDING_S_B : ADDING_S_S);
                    s_a[0] = d->stack_grid[i2][i1];
                    d_a[1] = (uchar) (derivs_layers[d->n_layers][j] ? ADDING_B_B : ADDING_B_S);
                    s_a[1] = &d->stack_chain[i];
                    flag2  = 1;
                    i1     = i - 1;
               }
               else
               if (derivs_layers[i][j]    &&   flag2 && i != d->n_layers - 1 && i == i1) {
                    /* add linearized layer */
                    n_c    = 0;
                    n_a    = 1;
                    d_a[0] = ADDING_B_B;
                    s_a[0] = &d->stack_chain[i];
                    i1     = i - 1;
               }
               else
               if (derivs_layers[i][j]    &&   flag2 && i != d->n_layers - 1 && i != i1) {
                    /* add  unlinearized stack, add linearized layer */
                    i2     = i + 1;
                    n_c    = 0;
                    n_a    = 2;
                    d_a[0] = ADDING_S_B;
                    s_a[0] = d->stack_grid[i2][i1];
                    d_a[1] = ADDING_B_B;
                    s_a[1] = &d->stack_chain[i];
                    i1     = i - 1;
               }
               else
               if (i == 0 && ! flag2 && ! surface) {
                    /* copy unlinearized stack */
                    i2     = i;
                    n_c    = 1;
                    d_h[0] = 0;
                    d_d[0] = 0;
                    s_c[0] = d->stack_grid[i2][i1];
                    n_a    = 0;
                    flag2  = 1;
                    i1     = i - 1;
               }
               else
               if (i == 0 && ! flag2 &&   surface) {
                    /* add unlinearized stack */
                    i2     = i;
                    n_c    = 0;
                    n_a    = 1;
                    d_a[0] = (uchar) (derivs_layers[d->n_layers][j] ? ADDING_U_B : ADDING_U_U);
                    s_a[0] = d->stack_grid[i2][i1];
                    flag2  = 1;
                    i1     = i - 1;
               }
               else
               if (i == 0 &&   flag2) {
                    /* add  unlinearized stack */
                    i2     = i;
                    n_c    = 0;
                    n_a    = 1;
                    d_a[0] = ADDING_U_B;
                    s_a[0] = d->stack_grid[i2][i1];
                    i1     = i - 1;
               }
               else {
                    n_c  = 0;
                    n_a  = 0;
               }

               for (ii = 0; ii < n_c; ++ii) {
                    s = s_c[ii];

                    if (s->i1 == s->i2) {
                         S_p1 = v1;
                         S_m1 = v2;

                         if (d->derivs.sources[s->i1][j]) {
                             S_p_l1 = &v3;
                             S_m_l1 = &v4;
                         }

                         source_vectors_scale_add(d, i_four, s->i1, j, 1, s->S_p[i_four], s->S_m[i_four], s->Sl_p[i_four], s->Sl_m[i_four], S_p1, S_m1, s->S_p_l[i_four]+j, s->S_m_l[i_four]+j, s->Sl_p_l[i_four]+j, s->Sl_m_l[i_four]+j, S_p_l1, S_m_l1, work2);
                    }
                    else {
                         S_p1 = s->S_p[i_four];     S_p_l1 = s->S_p_l[i_four]+j;
                         S_m1 = s->S_m[i_four];     S_m_l1 = s->S_m_l[i_four]+j;
                    }

                    if (s->i1 == s->i2) {
                         if (d->derivs.layers[s->i1][j]) {
                              R_m_l1 = s->R_m_l[i_four]+j;
                         }
                    }

                    layer_copy_ref_l(*R_m, *S_p, s->R_m[i_four], S_p1, (*R_m_l)+j, (*S_p_l)+j, R_m_l1, S_p_l1, d->n_quad_v_x, 1, &d_h[ii], &d_d[ii]);
               }

               for (ii = 0; ii < n_a; ++ii) {
                    s = s_a[ii];

                    if (s->i1 == s->i2) {
                         S_p1 = v1;
                         S_m1 = v2;

                         if (d->derivs.sources[s->i1][j]) {
                             S_p_l1 = &v3;
                             S_m_l1 = &v4;
                         }

                         source_vectors_scale_add(d, i_four, s->i1, j, 1, s->S_p[i_four], s->S_m[i_four], s->Sl_p[i_four], s->Sl_m[i_four], S_p1, S_m1, s->S_p_l[i_four]+j, s->S_m_l[i_four]+j, s->Sl_p_l[i_four]+j, s->Sl_m_l[i_four]+j, S_p_l1, S_m_l1, work2);
                    }
                    else {
                         S_p1 = s->S_p[i_four];     S_p_l1 = s->S_p_l[i_four]+j;
                         S_m1 = s->S_m[i_four];     S_m_l1 = s->S_m_l[i_four]+j;
                    }

                    flag3 = ! (j != d->n_derivs - 1 && s_a[ii]->i2 == 0);
/*
                    if ((i_stack = stack_chain_find2(d->n_stacks, d->stack_chain, s->i1, s->i2, s->i2 + 1, d->n_layers - 1)) >= 0)
                         layer_add_aux = &(d->stack_chain[i_stack].d[i_four]);
                    else
*/
                         layer_add_aux = NULL;

                    if (s->i1 == s->i2) {
                         if (d->derivs.layers[s->i1][j]) {
                              R_p_l1 = s->R_p_l[i_four]+j;
                              T_p_l1 = s->T_p_l[i_four]+j;
                              R_m_l1 = s->R_m_l[i_four]+j;
                              T_m_l1 = s->T_m_l[i_four]+j;
                         }
                    }

                    layer_add_ref2(s->R_m[i_four], s->T_m[i_four], S_m1, s->R_p[i_four], s->T_p[i_four], S_p1, *R_m, *S_p, *R_m, *S_p, R_m_l1, T_m_l1, S_m_l1, R_p_l1, T_p_l1, S_p_l1, (*R_m_l)+j, (*S_p_l)+j, (*R_m_l)+j, (*S_p_l)+j, d->n_quad_v_x, 1, 1., trans_l, &d_a[ii], 0, flag3, layer_add_aux == NULL, work2, layer_add_aux);
               }
          }
     }


     return 0;
}



/*******************************************************************************
 *
 ******************************************************************************/
static int get_total_R_T_S_U_W_V
     (xrtm_data *d, int i_four, int solver,
      double ***R_p, double ***T_p, double ***R_m, double ***T_m,
      double **S_p, double **S_m,
      double ****R_p_l, double ****T_p_l, double ****R_m_l, double ****T_m_l,
      double ***S_p_l, double ***S_m_l,
      int flag, save_tree_data save_tree, work_data *work) {

     uchar *derivs_layers;
     uchar *derivs_thermal;
     uchar *derivs_beam;
     uchar *derivs_sources;
     uchar *derivs_adding_down;

     int i;

     int deps;

     int flag2;

     double *trans_l;

     double *btran_l;

     double *S1_p;
     double *S1_m;

     double *S2_p;
     double *S2_m;
     double *Sl2_p;
     double *Sl2_m;

     double *S3_p;
     double *S3_m;

     double **R1_p;
     double **T1_p;
     double **R1_m;
     double **T1_m;

     double **R2_p;
     double **T2_p;
     double **R2_m;
     double **T2_m;

     double **S1_p_l;
     double **S1_m_l;

     double **S2_p_l;
     double **S2_m_l;
     double **Sl2_p_l;
     double **Sl2_m_l;

     double **S3_p_l;
     double **S3_m_l;

     double ***R1_p_l;
     double ***T1_p_l;
     double ***R1_m_l;
     double ***T1_m_l;

     double ***R2_p_l;
     double ***T2_p_l;
     double ***R2_m_l;
     double ***T2_m_l;

     forward_save_get_total_R_T_S_U_W_V_data *save;

     work_data work2;
     work_data work3;


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (d->options & XRTM_OPTION_REVERSE_DERIVS) {
          save_tree_encode_s(&save_tree, "get_total_R_T_S_U_W_V");

          if (save_tree_retrieve_data(&save_tree, forward_save_get_total_R_T_S_U_W_V_data, &save))
               forward_save_get_total_R_T_S_U_W_V_alloc(save, d->n_quad_v_x, d->n_layers);
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (d->options & XRTM_OPTION_SAVE_TOTAL_R_T_S) {
          deps = total_R_T_S_U_W_V_deps(d);

          flag2 = resolve_deps0(d, i_four, d->dep_flags_total_R_T_S_U_W_V, deps);

          *R_p = d->Ra_p[i_four];     *T_p = d->Ta_p[i_four];
          *R_m = d->Ra_m[i_four];     *T_m = d->Ta_m[i_four];
          *S_p = d->Sa_p[i_four];     *S_m = d->Sa_m[i_four];

          if (d->options & XRTM_OPTION_CALC_DERIVS) {
               *R_p_l = d->Ra_p_l[i_four];     *T_p_l = d->Ta_p_l[i_four];
               *R_m_l = d->Ra_m_l[i_four];     *T_m_l = d->Ta_m_l[i_four];
               *S_p_l = d->Sa_p_l[i_four];     *S_m_l = d->Sa_m_l[i_four];
          }

          if (! flag2)
               return 0;
     }
     else if (flag) {
          *R_p = get_work1(work, WORK_DXX);
          *T_p = get_work1(work, WORK_DXX);
          *R_m = get_work1(work, WORK_DXX);
          *T_m = get_work1(work, WORK_DXX);

          *S_p = get_work1(work, WORK_DX);
          *S_m = get_work1(work, WORK_DX);

          if (d->options & XRTM_OPTION_CALC_DERIVS && flags_or(d->derivs.total[0], d->n_derivs)) {
               *R_p_l = get_work2(work, WORK_DXX, WORK_DERIVS_V, d->derivs.total[0]);
               *T_p_l = get_work2(work, WORK_DXX, WORK_DERIVS_V, d->derivs.total[0]);
               *R_m_l = get_work2(work, WORK_DXX, WORK_DERIVS_V, d->derivs.total[0]);
               *T_m_l = get_work2(work, WORK_DXX, WORK_DERIVS_V, d->derivs.total[0]);

               *S_p_l = get_work2(work, WORK_DX,  WORK_DERIVS_V, d->derivs.total[0]);
               *S_m_l = get_work2(work, WORK_DX,  WORK_DERIVS_V, d->derivs.total[0]);
          }
     }


     work2 = *work;


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     R1_p = *R_p;     T1_p = *T_p;
     R1_m = *R_m;     T1_m = *T_m;
     S1_p = *S_p;     S1_m = *S_m;

     if (d->options & XRTM_OPTION_CALC_DERIVS) {
          R1_p_l = *R_p_l;     T1_p_l = *T_p_l;
          R1_m_l = *R_m_l;     T1_m_l = *T_m_l;
          S1_p_l = *S_p_l;     S1_m_l = *S_m_l;
     }

     S3_p = get_work1(&work2, WORK_DX);
     S3_m = get_work1(&work2, WORK_DX);


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (d->options & XRTM_OPTION_CALC_DERIVS) {
          trans_l = get_work1(&work2, WORK_DDERIVS);

          for (i = 0; i < d->n_derivs; ++i)
               trans_l[i] = 0.;
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     for (i = 0; i < d->n_layers; ++i) {
          work3 = work2;

          if (d->options & XRTM_OPTION_REVERSE_DERIVS)
               save_tree_recode_i(&save_tree, i, i == 0);

          if (d->options & XRTM_OPTION_CALC_DERIVS && flags_or(d->derivs.sources[i], d->n_derivs)) {
               S3_p_l = get_work2(&work3, WORK_DX, WORK_DERIVS_V, d->derivs.sources[i]);
               S3_m_l = get_work2(&work3, WORK_DX, WORK_DERIVS_V, d->derivs.sources[i]);
          }

          if (get_layer_R_T_S_U_W_V(d, i_four, i, solver, &R2_p, &T2_p, &R2_m, &T2_m, &S2_p, &S2_m, &Sl2_p, &Sl2_m, &R2_p_l, &T2_p_l, &R2_m_l, &T2_m_l, &S2_p_l, &S2_m_l, &Sl2_p_l, &Sl2_m_l, 1, save_tree, &work3)) {
               fprintf(stderr, "ERROR: get_layer_R_T_S_U_W_V()\n");
               return -1;
          }

          if (d->options & XRTM_OPTION_CALC_DERIVS) {
               derivs_layers      = d->derivs.layers[i];
               if (d->options & XRTM_OPTION_SOURCE_SOLAR)
                    derivs_beam    = d->derivs.beam[i];
               if (d->options & XRTM_OPTION_SOURCE_THERMAL)
                    derivs_thermal = d->derivs.thermal[i];
               derivs_sources     = d->derivs.sources[i];
               derivs_adding_down = d->derivs.adding_down[i];
        }

          update_beam_params(d, i, work3);

          if (d->options & XRTM_OPTION_CALC_DERIVS)
               btran_l = d->btran_l[i];

          scale_add_source_vectors(d->n_quad_v_x, d->n_derivs, d->btran[i], btran_l, S2_p, S2_m, Sl2_p, Sl2_m, S3_p, S3_m, S2_p_l, S2_m_l, Sl2_p_l, Sl2_m_l, S3_p_l, S3_m_l, d->options & XRTM_OPTION_SOURCE_SOLAR, d->options & XRTM_OPTION_SOURCE_THERMAL && i_four == 0, derivs_layers, derivs_beam, derivs_thermal, derivs_sources, work3);

          if (d->options & XRTM_OPTION_REVERSE_DERIVS) {
               layer_copy_all(save->R1_m[i], save->T1_m[i], save->S1_m[i], save->R1_p[i], save->T1_p[i], save->S1_p[i], R1_m, T1_m, S1_m, R1_p, T1_p, S1_p, d->n_quad_v_x);

               layer_copy_all(save->R2_m[i], save->T2_m[i], save->S2_m[i], save->R2_p[i], save->T2_p[i], save->S2_p[i], R2_m, T2_m, S2_m, R2_p, T2_p, S2_p, d->n_quad_v_x);

               dvec_copy(save->S3_p[i], S3_p, d->n_quad_v_x);
               dvec_copy(save->S3_m[i], S3_m, d->n_quad_v_x);
          }

          if (i == 0)
               layer_copy_all_l(R1_m, T1_m, S1_m, R1_p, T1_p, S1_p,
                                R2_m, T2_m, S3_m, R2_p, T2_p, S3_p,
                                R1_m_l, T1_m_l, S1_m_l, R1_p_l, T1_p_l, S1_p_l,
                                R2_m_l, T2_m_l, S3_m_l, R2_p_l, T2_p_l, S3_p_l,
                                d->n_quad_v_x, d->n_derivs, derivs_adding_down);
          else
               layer_add_all   (R1_m, T1_m, S1_m, R1_p, T1_p, S1_p,
                                R2_m, T2_m, S3_m, R2_p, T2_p, S3_p,
                                R1_m, T1_m, S1_m, R1_p, T1_p, S1_p,
                                R1_m_l, T1_m_l, S1_m_l, R1_p_l, T1_p_l, S1_p_l,
                                R2_m_l, T2_m_l, S3_m_l, R2_p_l, T2_p_l, S3_p_l,
                                R1_m_l, T1_m_l, S1_m_l, R1_p_l, T1_p_l, S1_p_l,
                                d->n_quad_v_x, d->n_derivs, 1., trans_l,
                                derivs_adding_down, 1, save_tree, work3);
     }

     if (d->options & XRTM_OPTION_REVERSE_DERIVS)
          save_tree_decode_i(&save_tree);

#ifdef USE_AD_FOR_TL_CALC_GET_TOTAL_R_T_S_U_W_V
     get_total_R_T_S_U_W_V_tl_with_ad(d, i_four, solver, *R_p, *T_p, *R_m, *T_m, *S_p, *S_m, *R_p_l, *T_p_l, *R_m_l, *T_m_l, *S_p_l, *S_m_l, flag, save_tree, work2);
#endif

     return 0;
}



/*******************************************************************************
 *
 ******************************************************************************/
static int get_total_R_T_S_U_W_V_2
     (xrtm_data *d, int i_four, int solver,
      double  ***R_p, double  ***T_p, double  ***R_m, double  ***T_m,
      double  **S_p, double  **S_m,
      double ****R_p_l, double ****T_p_l, double ****R_m_l, double ****T_m_l,
      double ***S_p_l, double ***S_m_l,
      int flag, save_tree_data save_tree, work_data *work) {

     uchar d_c[2];
     uchar d_a[2];

     int i;
     int i1;
     int i2;
     int ii;
     int j;

     int deps;

     int flag2;
     int flag3;

     int n_c;
     int n_a;

     int i_stack;

     double *trans_l;

     double *v1;
     double *v2;
     double *v3;
     double *v4;

     double *S_p1;
     double *S_m1;

     double **S_p_l1;
     double **S_m_l1;

     double ***R_p_l1;
     double ***T_p_l1;
     double ***R_m_l1;
     double ***T_m_l1;

     layer_add_aux_data *layer_add_aux;

     stack_data *s;

     stack_data *s_c[2];
     stack_data *s_a[2];

     work_data work2;


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (d->options & XRTM_OPTION_SAVE_TOTAL_R_T_S) {
          deps = total_R_T_S_U_W_V_deps(d);

          flag2 = resolve_deps0(d, i_four, d->dep_flags_total_R_T_S_U_W_V, deps);

          *R_p = d->Ra_p[i_four];     *T_p = d->Ta_p[i_four];
          *R_m = d->Ra_m[i_four];     *T_m = d->Ta_m[i_four];
          *S_p = d->Sa_p[i_four];     *S_m = d->Sa_m[i_four];

          if (d->options & XRTM_OPTION_CALC_DERIVS) {
               *R_p_l = d->Ra_p_l[i_four];     *T_p_l = d->Ta_p_l[i_four];
               *R_m_l = d->Ra_m_l[i_four];     *T_m_l = d->Ta_m_l[i_four];
               *S_p_l = d->Sa_p_l[i_four];     *S_m_l = d->Sa_m_l[i_four];
          }

          if (! flag2)
               return 0;
     }
     else if (flag) {
          *R_p = get_work1(work, WORK_DXX);
          *T_p = get_work1(work, WORK_DXX);
          *R_m = get_work1(work, WORK_DXX);
          *T_m = get_work1(work, WORK_DXX);

          *S_p = get_work1(work, WORK_DX);
          *S_m = get_work1(work, WORK_DX);

          if (d->options & XRTM_OPTION_CALC_DERIVS && flags_or(d->derivs.total[0], d->n_derivs)) {
               *R_p_l = get_work2(work, WORK_DXX, WORK_DERIVS_V, d->derivs.total[0]);
               *T_p_l = get_work2(work, WORK_DXX, WORK_DERIVS_V, d->derivs.total[0]);
               *R_m_l = get_work2(work, WORK_DXX, WORK_DERIVS_V, d->derivs.total[0]);
               *T_m_l = get_work2(work, WORK_DXX, WORK_DERIVS_V, d->derivs.total[0]);

               *S_p_l = get_work2(work, WORK_DX,  WORK_DERIVS_V, d->derivs.total[0]);
               *S_m_l = get_work2(work, WORK_DX,  WORK_DERIVS_V, d->derivs.total[0]);
          }
     }


     work2 = *work;


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (build_layer_stacks(d, i_four, solver, save_tree, work2)) {
          fprintf(stderr, "ERROR: build_layer_stacks()\n");
          return -1;
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     v1 = get_work1(&work2, WORK_DX);
     v2 = get_work1(&work2, WORK_DX);
     v3 = get_work1(&work2, WORK_DX);
     v4 = get_work1(&work2, WORK_DX);


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (d->options & XRTM_OPTION_CALC_DERIVS) {
          trans_l = get_work1(&work2, WORK_DDERIVS);

          for (i = 0; i < d->n_derivs; ++i)
               trans_l[i] = 0.;
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     for (j = 0; j < d->n_derivs; ++j) {
          i1    = 0;
          flag2 = 0;

          for (i = 0; i < d->n_layers; ++i) {

               if (d->derivs.layers[i][j]    && ! flag2 && i == 0) {
                    /* copy linearized layer */
                    n_c    = 1;
                    d_c[0] = ADDING_U_B;
                    s_c[0] = &d->stack_chain[i];
                    n_a    = 0;
                    flag2  = 1;
                    i1     = i + 1;
               }
               else
               if (d->derivs.layers[i][j]    && ! flag2 && i != 0) {
                    /* copy unlinearized stack, add linearized layer */
                    i2     = i - 1;
                    n_c    = 1;
                    d_c[0] = 0;
                    s_c[0] = d->stack_grid[i1][i2];
                    n_a    = 1;
                    d_a[0] = ADDING_U_B;
                    s_a[0] = &d->stack_chain[i];
                    flag2  = 1;
                    i1     = i + 1;
               }
               else
               if (d->derivs.layers[i][j]    &&   flag2 && i != 0 && i == i1) {
                    /* add linearized layer */
                    n_c    = 0;
                    n_a    = 1;
                    d_a[0] = ADDING_B_B;
                    s_a[0] = &d->stack_chain[i];
                    i1     = i + 1;
               }
               else
               if (d->derivs.layers[i][j]    &&   flag2 && i != 0 && i != i1) {
                    /* add  unlinearized stack, add linearized layer */
                    i2     = i - 1;
                    n_c    = 0;
                    n_a    = 2;
                    d_a[0] = ADDING_B_S;
                    s_a[0] = d->stack_grid[i1][i2];
                    d_a[1] = ADDING_B_B;
                    s_a[1] = &d->stack_chain[i];
                    i1     = i + 1;
               }
               else
               if (i == d->n_layers - 1 && ! flag2) {
                    /* copy unlinearized stack */
                    i2     = i;
                    n_c    = 1;
                    d_c[0] = 0;
                    s_c[0] = d->stack_grid[i1][i2];
                    n_a    = 0;
                    flag2  = 1;
                    i1     = i + 1;
               }
               else
               if (i == d->n_layers - 1 &&   flag2) {
                    /* add  unlinearized stack */
                    i2     = i;
                    n_c    = 0;
                    n_a    = 1;
                    d_a[0] = ADDING_B_S;
                    s_a[0] = d->stack_grid[i1][i2];
                    i1     = i + 1;
               }
               else {
                    n_c  = 0;
                    n_a  = 0;
               }

               for (ii = 0; ii < n_c; ++ii) {
                    s = s_c[ii];

                    if (s->i1 == s->i2) {
                         S_p1 = v1;
                         S_m1 = v2;

                         if (d->derivs.sources[s->i1][j]) {
                             S_p_l1 = &v3;
                             S_m_l1 = &v4;
                         }

                         source_vectors_scale_add(d, i_four, s->i1, j, 1, s->S_p[i_four], s->S_m[i_four], s->Sl_p[i_four], s->Sl_m[i_four], S_p1, S_m1, s->S_p_l[i_four]+j, s->S_m_l[i_four]+j, s->Sl_p_l[i_four]+j, s->Sl_m_l[i_four]+j, S_p_l1, S_m_l1, work2);
                    }
                    else {
                         S_p1 = s->S_p[i_four];     S_p_l1 = s->S_p_l[i_four]+j;
                         S_m1 = s->S_m[i_four];     S_m_l1 = s->S_m_l[i_four]+j;
                    }

                    if (s->i1 == s->i2) {
                         if (d->derivs.layers[s->i1][j]) {
                              R_m_l1 = s->R_m_l[i_four]+j;
                              T_m_l1 = s->T_m_l[i_four]+j;
                              R_p_l1 = s->R_p_l[i_four]+j;
                              T_p_l1 = s->T_p_l[i_four]+j;
                         }
                    }

                    layer_copy_all_l(*R_m, *T_m, *S_m, *R_p, *T_p, *S_p, s->R_m[i_four], s->T_m[i_four], S_m1, s->R_p[i_four], s->T_p[i_four], S_p1, (*R_m_l)+j, (*T_m_l)+j, (*S_m_l)+j, (*R_p_l)+j, (*T_p_l)+j, (*S_p_l)+j, R_m_l1, T_m_l1, S_m_l1, R_p_l1, T_p_l1, S_p_l1, d->n_quad_v_x, 1, &d_c[ii]);
               }

               for (ii = 0; ii < n_a; ++ii) {
                    s = s_a[ii];

                    if (s->i1 == s->i2) {
                         S_p1 = v1;
                         S_m1 = v2;

                         if (d->derivs.sources[s->i1][j]) {
                             S_p_l1 = &v3;
                             S_m_l1 = &v4;
                         }

                         source_vectors_scale_add(d, i_four, s->i1, j, 1, s->S_p[i_four], s->S_m[i_four], s->Sl_p[i_four], s->Sl_m[i_four], S_p1, S_m1, s->S_p_l[i_four]+j, s->S_m_l[i_four]+j, s->Sl_p_l[i_four]+j, s->Sl_m_l[i_four]+j, S_p_l1, S_m_l1, work2);
                    }
                    else {
                         S_p1 = s->S_p[i_four];     S_p_l1 = s->S_p_l[i_four]+j;
                         S_m1 = s->S_m[i_four];     S_m_l1 = s->S_m_l[i_four]+j;
                    }

                    flag3 = ! (j != d->n_derivs - 1 && s_a[ii]->i2 == d->n_layers - 1);

                    if ((i_stack = stack_chain_find2(d->n_stacks, d->stack_chain, 0, s->i1 - 1, s->i1, s->i2)) >= 0)
                         layer_add_aux = &(d->stack_chain[i_stack].d[i_four]);
                    else
                         layer_add_aux = NULL;

                    if (s->i1 == s->i2) {
                         if (d->derivs.layers[s->i1][j]) {
                              R_m_l1 = s->R_m_l[i_four]+j;
                              T_m_l1 = s->T_m_l[i_four]+j;
                              R_p_l1 = s->R_p_l[i_four]+j;
                              T_p_l1 = s->T_p_l[i_four]+j;
                         }
                    }

                    layer_add_all2(*R_m, *T_m, *S_m, *R_p, *T_p, *S_p, s->R_m[i_four], s->T_m[i_four], S_m1, s->R_p[i_four], s->T_p[i_four], S_p1, *R_m, *T_m, *S_m, *R_p, *T_p, *S_p, (*R_m_l)+j, (*T_m_l)+j, (*S_m_l)+j, (*R_p_l)+j, (*T_p_l)+j, (*S_p_l)+j, R_m_l1, T_m_l1, S_m_l1, R_p_l1, T_p_l1, S_p_l1, (*R_m_l)+j, (*T_m_l)+j, (*S_m_l)+j, (*R_p_l)+j, (*T_p_l)+j, (*S_p_l)+j, d->n_quad_v_x, 1, 1., trans_l, &d_a[ii], flag3, layer_add_aux == NULL, work2, layer_add_aux);
               }
          }
     }


     return 0;
}



/*******************************************************************************
 *
 ******************************************************************************/
static void brdf_vecs_get(xrtm_data *d, int i_four, int i_offset, int n_mus, int j_offset, enum xrtm_brdf_geometry brdf_geom, double **kernel_f, double ***kernel_f_l, double *alpha1, double *R, double **R_l, work_data work) {

     int i;

     int n_mus_v;

     uchar *derivs_layers;

     if (i_four == 0) {
          if (d->options & XRTM_OPTION_CALC_DERIVS)
               derivs_layers = d->derivs.layers[d->n_layers];
          brdf_build_kernel_vecs(i_offset, n_mus, j_offset, d->n_stokes, d->n_derivs, d->n_kernels, d->n_kernel_quad, d->kernels, d->kernel_ampfac, d->kernel_ampfac_l, d->kernel_params, d->kernel_params_l, d->kernel_func, d->kernel_qx, d->kernel_qw, &d->kernel_aux, brdf_geom, kernel_f, kernel_f_l, derivs_layers, work);
     }

     brdf_build_ref_vec(i_four, i_offset, n_mus, j_offset, d->n_stokes, d->qf, d->qx_v, d->qw_v, d->n_kernels, d->n_kernel_quad, d->kernels, d->kernel_qx, d->kernel_qw, kernel_f, R, work);

     if (d->options & XRTM_OPTION_CALC_DERIVS && flags_or(d->derivs.layers[d->n_layers], d->n_derivs)) {
          for (i = 0; i < d->n_derivs; ++i) {
               if (! d->derivs.layers[d->n_layers][i])
                    continue;

               brdf_build_ref_vec(i_four, i_offset, n_mus, j_offset, d->n_stokes, d->qf, d->qx_v, d->qw_v, d->n_kernels, d->n_kernel_quad, d->kernels, d->kernel_qx, d->kernel_qw, kernel_f_l[i], R_l[i], work);
          }
     }

     if (d->misc_input.use_symmetric_form) {
          n_mus_v = n_mus * d->n_stokes;

          vec_sim_trans(n_mus_v, R, alpha1);

          if (d->options & XRTM_OPTION_CALC_DERIVS && flags_or(d->derivs.layers[d->n_layers], d->n_derivs)) {
               for (i = 0; i < d->n_derivs; ++i) {
                    if (! d->derivs.layers[d->n_layers][i])
                         continue;

                    vec_sim_trans(n_mus_v, R_l[i], alpha1);
               }
          }
     }
}



/*******************************************************************************
 *
 ******************************************************************************/
static void get_brdf_vecs_q0(xrtm_data *d, int i_four, double *Rs_q0, double **Rs_q0_l, work_data work) {

     brdf_vecs_get(d, i_four, 0, d->n_quad_x, d->n_quad + d->n_umus, XRTM_BRDF_GEOMETRY_Q0, d->kernel_f_q0, d->kernel_f_q0_l, d->alpha1, Rs_q0, Rs_q0_l, work);
}



/*******************************************************************************
 *
 ******************************************************************************/
static void get_brdf_vecs_u0(xrtm_data *d, int i_four, double *Rs_u0, double **Rs_u0_l, work_data work) {

     brdf_vecs_get(d, i_four, d->n_quad, d->n_umus, d->n_quad + d->n_umus, XRTM_BRDF_GEOMETRY_U0, d->kernel_f_u0, d->kernel_f_u0_l, d->alpha1 + d->n_quad_v, Rs_u0, Rs_u0_l, work);
}



/*******************************************************************************
 *
 ******************************************************************************/
static void brdf_mats_get(xrtm_data *d, int i_four, int i_offset, int n_mus1, int j_offset, int n_mus2, enum xrtm_brdf_geometry brdf_geom, double ***kernel_f, double ****kernel_f_l, double *alpha1, double *alpha2, double **R, double ***R_l, work_data work) {

     int i;

     int n_mus_v1;
     int n_mus_v2;

     uchar *derivs_layers;

     if (i_four == 0) {
          if (d->options & XRTM_OPTION_CALC_DERIVS)
               derivs_layers = d->derivs.layers[d->n_layers];

          brdf_build_kernel_mats(i_offset, n_mus1, j_offset, n_mus2, d->n_stokes, d->n_derivs, d->n_kernels, d->n_kernel_quad, d->kernels, d->kernel_ampfac, d->kernel_ampfac_l, d->kernel_params, d->kernel_params_l, d->kernel_func, d->kernel_qx, d->kernel_qw, &d->kernel_aux, brdf_geom, kernel_f, kernel_f_l, derivs_layers, work);
     }

     brdf_build_ref_mat(i_four, i_offset, n_mus1, j_offset, n_mus2, d->n_stokes, d->qf, d->qx_v, d->qw_v, d->n_kernels, d->n_kernel_quad, d->kernels, d->kernel_qx, d->kernel_qw, kernel_f, R, work);

     if (d->options & XRTM_OPTION_CALC_DERIVS && flags_or(d->derivs.layers[d->n_layers], d->n_derivs)) {
          for (i = 0; i < d->n_derivs; ++i) {
               if (! d->derivs.layers[d->n_layers][i])
                    continue;

               brdf_build_ref_mat(i_four, i_offset, n_mus1, j_offset, n_mus2, d->n_stokes, d->qf, d->qx_v, d->qw_v, d->n_kernels, d->n_kernel_quad, d->kernels, d->kernel_qx, d->kernel_qw, kernel_f_l[i], R_l[i], work);
          }
     }

     if (d->misc_input.use_symmetric_form) {
          n_mus_v1 = n_mus1 * d->n_stokes;
          n_mus_v2 = n_mus2 * d->n_stokes;

          mat_sim_trans(n_mus_v1, n_mus_v2, R, alpha1, alpha2);

          if (d->options & XRTM_OPTION_CALC_DERIVS && flags_or(d->derivs.layers[d->n_layers], d->n_derivs)) {
               for (i = 0; i < d->n_derivs; ++i) {
                    if (! d->derivs.layers[d->n_layers][i])
                         continue;

                    mat_sim_trans(n_mus_v1, n_mus_v2, R_l[i], alpha1, alpha2);
               }
          }
     }
}



/*******************************************************************************
 *
 ******************************************************************************/
static void get_brdf_mats_qq(xrtm_data *d, int i_four, double **Rs_qq, double ***Rs_qq_l, work_data work) {

     brdf_mats_get(d, i_four, 0, d->n_quad_x, 0, d->n_quad_x, XRTM_BRDF_GEOMETRY_QQ, d->kernel_f_qq, d->kernel_f_qq_l, d->alpha1, d->alpha2, Rs_qq, Rs_qq_l, work);
}



/*******************************************************************************
 *
 ******************************************************************************/
static void get_brdf_mats_uq(xrtm_data *d, int i_four, double **Rs_uq, double ***Rs_uq_l, work_data work) {

     brdf_mats_get(d, i_four, d->n_quad, d->n_umus, 0, d->n_quad_x, XRTM_BRDF_GEOMETRY_UQ, d->kernel_f_uq, d->kernel_f_uq_l, d->alpha1 + d->n_quad_v, d->alpha2, Rs_uq, Rs_uq_l, work);
}



/*******************************************************************************
 *
 ******************************************************************************/
static void diff_bound_input_update(xrtm_data *d, int i_four, int *dep_flags, int deps, double btran, double *btran_l, save_tree_data save_tree, work_data work) {

     char temp[sizeof("update_diff_bound_input") + 6];

     int i;
     int ii;
     int iii;
     int j;
     int jj;
     int jjj;
     int k;
     int kk;

     int flag;

     int solar_surface;
     int thermal_surface;

     double a;
     double b;
     double c;
     double e;

     double *Rs_q0;
     double **Rs_q0_l;

     double **Rs_qq;
     double ***Rs_qq_l;

     double *Rs_u0;
     double **Rs_u0_l;
     double **Rs_uq;
     double ***Rs_uq_l;

     forward_save_update_diff_bound_input_data *save;


     solar_surface   = solar_surface_for_i_four(d, i_four);
     thermal_surface = thermal_surface_for_i_four(d, i_four);


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (d->options & XRTM_OPTION_REVERSE_DERIVS && (solar_surface || thermal_surface)) {
          snprintf(temp, sizeof("update_diff_bound_input") + 6, "%s_%04d", "update_diff_bound_input", i_four);

          save_tree_encode_s(&save_tree, "update_diff_bound_input");

          if (save_tree_retrieve_proxy(&save_tree, temp, forward_save_update_diff_bound_input_data, &save))
               forward_save_update_diff_bound_input_alloc(save, d->n_quad_v, d->n_umus_v);
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     flag = resolve_deps0(d, i_four, dep_flags, deps);

     if (flag) {
          init_array1_d(d->In_p[i_four], d->n_quad_v + d->n_umus_v, 0.);
          init_array1_d(d->I1_m[i_four], d->n_quad_v + d->n_umus_v, 0.);

          if (d->options & XRTM_OPTION_CALC_DERIVS) {
               init_array2_d(d->In_p_l[i_four], d->n_derivs, d->n_quad_v + d->n_umus_v, 0.);
               init_array2_d(d->I1_m_l[i_four], d->n_derivs, d->n_quad_v + d->n_umus_v, 0.);
          }


          /*--------------------------------------------------------------------
           *
           *------------------------------------------------------------------*/
if (i_four == 0) {
          ii = 0;
          for (i = 0; i < d->n_quad + d->n_umus; ++i) {
               d->I1_m[i_four][ii] += d->F_iso_top;
               d->In_p[i_four][ii] += d->F_iso_bot;
               ii += d->n_stokes;
          }

          for (i = 0; i < d->n_derivs; ++i) {
               jj = 0;
               for (j = 0; j < d->n_quad + d->n_umus; ++j) {
                    d->I1_m_l[i_four][i][jj] += d->F_iso_top_l[i];
                    d->In_p_l[i_four][i][jj] += d->F_iso_bot_l[i];
                    jj += d->n_stokes;
               }
          }
}

          /*--------------------------------------------------------------------
           *
           *------------------------------------------------------------------*/
          if (d->options & XRTM_OPTION_SOURCE_SOLAR && solar_surface) {
               Rs_q0 = get_work1(&work, WORK_DX);

               if (d->options & XRTM_OPTION_CALC_DERIVS)
                    Rs_q0_l = get_work2(&work, WORK_DX, WORK_DERIVS_V, NULL);

               get_brdf_vecs_q0(d, i_four, Rs_q0, Rs_q0_l, work);

               if (d->options & XRTM_OPTION_SFI && d->n_umus > 0) {
                    Rs_u0 = get_work1(&work, WORK_DU);

                    if (d->options & XRTM_OPTION_CALC_DERIVS)
                         Rs_u0_l = get_work2(&work, WORK_DU, WORK_DERIVS_V, NULL);

                    get_brdf_vecs_u0(d, i_four, Rs_u0, Rs_u0_l, work);
               }

               a = d->qf * d->mu_0 * d->F_0 / PI;

               b = a * btran;

               for (i = 0; i < d->n_quad_v_x; ++i)
                    d->In_p[i_four][i] += b * Rs_q0[i];

               ii = 0;
               for ( ; i < d->n_quad_v + d->n_umus_v; ++i)
                    d->In_p[i_four][i] += b * Rs_u0[ii++];

               for (i = 0; i < d->n_derivs; ++i) {
                    c = a * btran_l[i];

                    for (j = 0; j < d->n_quad_v_x; ++j) {

                         d->In_p_l[i_four][i][j] += c * Rs_q0[j];

                         if (d->derivs.layers[d->n_layers][i])
                              d->In_p_l[i_four][i][j] += b * Rs_q0_l[i][j];
                    }

                    jj = 0;
                    for (; j < d->n_quad_v + d->n_umus_v; ++j) {

                         d->In_p_l[i_four][i][j] += c * Rs_u0[jj];

                         if (d->derivs.layers[d->n_layers][i])
                              d->In_p_l[i_four][i][j] += b * Rs_u0_l[i][jj];
                         jj++;
                    }
               }

               if (d->options & XRTM_OPTION_REVERSE_DERIVS) {
                    dvec_copy(save->Rs_q0, Rs_q0, d->n_quad_v_x);

                    if (d->options & XRTM_OPTION_SFI && d->n_umus > 0)
                         dvec_copy(save->Rs_u0, Rs_u0, d->n_umus_v);
               }
          }


          /*--------------------------------------------------------------------
           *
           *------------------------------------------------------------------*/
          if (d->options & XRTM_OPTION_SOURCE_THERMAL && thermal_surface) {
               Rs_qq = get_work1(&work, WORK_DXX);

               if (d->options & XRTM_OPTION_CALC_DERIVS)
                    Rs_qq_l = get_work2(&work, WORK_DXX, WORK_DERIVS_V, NULL);

               get_brdf_mats_qq(d, i_four, Rs_qq, Rs_qq_l, work);

               if (d->options & XRTM_OPTION_SFI && d->n_umus > 0) {
                    Rs_uq = get_work1(&work, WORK_DUX);

                    if (d->options & XRTM_OPTION_CALC_DERIVS)
                         Rs_uq_l = get_work2(&work, WORK_DUX, WORK_DERIVS_V, NULL);

                    get_brdf_mats_uq(d, i_four, Rs_uq, Rs_uq_l, work);
               }

               ii = 0;
               for (i = 0; i < d->n_quad_x; ++i) {
                    a  = 0.;
                    jj = 0;
                    for (j = 0; j < d->n_quad_x; ++j) {
                         a  += Rs_qq[ii][jj];
                         jj += d->n_stokes;
                    }
                    d->In_p[i_four][ii] += d->surface_b * (1. - a);

                    ii += d->n_stokes;
               }

               iii = ii; ii = 0;
               for ( ; i < d->n_quad + d->n_umus; ++i) {
                    b  = 0.;
                    jj = 0;
                    for (j = 0; j < d->n_quad_x; ++j) {
                         b  += Rs_uq[ii][jj];
                         jj += d->n_stokes;
                    }
                    d->In_p[i_four][iii] += d->surface_b * (1. - b);

                    ii  += d->n_stokes;
                    iii += d->n_stokes;
               }

               for (i = 0; i < d->n_derivs; ++i) {
/*
                    if (d->derivs.layers[d->n_layers][i] || d->surface_b_l[i] > 0.) {

                    if (d->derivs.layers[d->n_layers][i]) {
*/
                         c = d->surface_b_l[i] * (1. - a);

                         jj = 0;
                         for (j = 0; j < d->n_quad_x; ++j) {
                              e  = 0.;
                              kk = 0;
if (d->derivs.layers[d->n_layers][i]) {
                              for (k = 0; k < d->n_quad_x; ++k) {
                                   e  += Rs_qq_l[i][jj][kk];
                                   kk += d->n_stokes;
                              }
}
                              d->In_p_l[i_four][i][jj] += c + d->surface_b * -e;

                              jj += d->n_stokes;
                         }

                         c = d->surface_b_l[i] * (1. - b);

                         jjj = jj; jj = 0;
                         for ( ; j < d->n_quad + d->n_umus; ++j) {
                              e  = 0.;
                              kk = 0;
if (d->derivs.layers[d->n_layers][i]) {
                              for (k = 0; k < d->n_quad_x; ++k) {
                                   e += Rs_uq_l[i][jj][kk];
                                   kk += d->n_stokes;
                              }
}
                              d->In_p_l[i_four][i][jjj] += c + d->surface_b * -e;

                              jj  += d->n_stokes;
                              jjj += d->n_stokes;
                         }
/*
                   }
*/
               }
          }
     }
}



static int diff_bound_input_deps(xrtm_data *d) {

     return opt_props_deps(d) | beam_params_deps(d) |
            DEP_MASK_F_ISO_TOP | DEP_MASK_F_ISO_BOT | DEP_MASK_F_0 | DEP_MASK_MU_0 | DEP_MASK_BRDF | DEP_MASK_BRDF_L;
}



static void update_diff_bound_input(xrtm_data *d, int i_four, save_tree_data save_tree, work_data work) {

     double *btran_l;

     update_beam_params(d, d->n_layers, work);

     if (d->options & XRTM_OPTION_CALC_DERIVS)
          btran_l = d->btran_l[d->n_layers];

     diff_bound_input_update(d, i_four, d->dep_flags_diff_bound_input, diff_bound_input_deps(d), d->btran[d->n_layers], btran_l, save_tree, work);
#ifdef USE_AD_FOR_TL_CALC_UPDATE_DIFF_BOUND_INPUT
     update_diff_bound_input_tl_with_ad(d, i_four, save_tree, work);
#endif
}


/*
static int diff_bound_input0_deps(xrtm_data *d) {

     return beam_params_deps(d) |
            DEP_MASK_F_0 | DEP_MASK_MU_0 | DEP_MASK_BRDF | DEP_MASK_BRDF_L;
}



static void update_diff_bound_input0(xrtm_data *d, int i_four, work_data work) {

     double *btran0_l;

     if (! (d->options & XRTM_OPTION_DELTA_M))
          update_diff_bound_input(d, i_four, work);

     update_beam_params0(d, d->n_layers, work);

     if (d->options & XRTM_OPTION_CALC_DERIVS)
          btran0_l = d->btran0_l[d->n_layers];

     diff_bound_input_update(d, i_four, d->dep_flags_diff_bound_input0, diff_bound_input0_deps(d), d->btran0[d->n_layers], btran0_l, work);
}
*/


/*******************************************************************************
 *
 ******************************************************************************/
static void copy_polys(xrtm_data *d, void *polys2, void *polys1) {

     if (d->options & XRTM_OPTION_PHASE_SCALAR)
          dvec_copy((double *) polys2, (double *) polys1, d->n_coef);
     else
     if (d->options & XRTM_OPTION_PHASE_MATRIX_GC)
          dmat_copy((double **) polys2, (double **) polys1, 4, d->n_coef);
     else
     if (d->options & XRTM_OPTION_PHASE_MATRIX_LC)
          dvec_copy((double *) polys2, (double *) polys1, d->n_coef);
}



/*******************************************************************************
 *
 ******************************************************************************/
static void get_single(xrtm_data *d, int n_phis, double **phis, double ****I_p, double ****I_m, double *****I_p_l, double *****I_m_l, save_tree_data save_tree, work_data work) {

     int i;
     int ii;
     int j;
     int k;
     int l;
/*
     int n_quad_x;
*/
     int i_mus2_v;
     int n_mus2;

     double a;

     double *mus2;

     double f;
     double *f_l;

     double *omega2;
     double **omega2_l;

     double **In_p_l;
     double **I1_m_l;

     double **P_full;
     double ***P_full_l;

     double **I;
     double ***I_l;

     void *polys;

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
     if (d->options & XRTM_OPTION_REVERSE_DERIVS) {
          save_tree_encode_s(&save_tree, "get_single");

          if (save_tree_retrieve_data(&save_tree, forward_save_get_single_data, &save))
               forward_save_get_single_alloc(save, d->options, d->n_coef, d->n_quad_v, d->n_umus_v, n_mus2, n_phis, d->n_layers, d->n_stokes);
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (d->options & XRTM_OPTION_N_T_TMS) {
          if (d->options & XRTM_OPTION_CALC_DERIVS)
               f_l = get_work1(&work, WORK_DDERIVS);

          omega2 = get_work1(&work, WORK_DLAYERS);

          if (d->options & XRTM_OPTION_CALC_DERIVS)
               omega2_l = get_work1(&work, WORK_DBOTH);
     }

     if (d->options & XRTM_OPTION_PHASE_SCALAR)
          polys = get_work_d1(&work, d->n_coef);
     else
     if (d->options & XRTM_OPTION_PHASE_MATRIX_GC)
          polys = get_work_d2(&work, 4, d->n_coef);
     else
     if (d->options & XRTM_OPTION_PHASE_MATRIX_LC)
          polys = get_work_d1(&work, d->n_coef);

     P_full = get_work2(&work, WORK_DS, WORK_LAYERS_V, NULL);

     I      = get_work_d2(&work, d->n_ulevels, d->n_stokes);

     if (d->options & XRTM_OPTION_CALC_DERIVS) {
          P_full_l = get_work3(&work, WORK_DS, WORK_BOTH_V, d->derivs.layers);

          In_p_l   = get_work2(&work, WORK_DD, WORK_DERIVS_V, NULL);
          I1_m_l   = get_work2(&work, WORK_DD, WORK_DERIVS_V, NULL);

          I_l      = get_work_d3(&work, d->n_ulevels, d->n_derivs, d->n_stokes);
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     update_utaus(d, work);

     update_opt_props_all(d, work);

     update_beam_params_all(d, work);

     update_diff_bound_input(d, 0, save_tree, work);

     if (! (d->options & XRTM_OPTION_N_T_TMS)) {
           omega2   = d->omega0;
           omega2_l = d->omega0_l;
     }
     else {
          a = 2. * d->n_coef2 + 1.;

          for (i = 0; i < d->n_layers; ++i) {
               get_delta_m_f(d, i, 0, d->n_derivs, a, &f, f_l);

               if (! (d->options & XRTM_OPTION_CALC_DERIVS))
                    n_t_tms_scaling(d->n_derivs, f, f_l, d->omega0[i], &omega2[i], NULL, NULL);
               else
                    n_t_tms_scaling(d->n_derivs, f, f_l, d->omega0[i], &omega2[i], d->omega0_l[i], omega2_l[i]);
          }
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
/*
     init_array4_d(I_p, d->n_ulevels, n_phis, n_mus2, d->n_stokes, 0.);
*/
     init_array4_d(I_m, d->n_ulevels, n_phis, n_mus2, d->n_stokes, 0.);

     if (d->options & XRTM_OPTION_CALC_DERIVS) {
          init_array5_d(I_p_l, d->n_ulevels, d->n_derivs, n_phis, n_mus2, d->n_stokes, 0.);
/*
          init_array5_d(I_m_l, d->n_ulevels, d->n_derivs, n_phis, n_mus2, d->n_stokes, 0.);
*/
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     ii = i_mus2_v;

     for (i = 0; i < n_mus2; ++i) {
          for (j = 0; j < n_phis; ++j) {

          if (d->options & XRTM_OPTION_REVERSE_DERIVS)
               save_tree_recode_i_j(&save_tree, i, j, i == 0 && j == 0);


               phase_polys_get(d, d->n_coef, d->mu_0, d->phi_0, -mus2[i], phis[i][j], polys);

               get_phase_func_full_all(d, -mus2[i], phis[i][j], polys, P_full, P_full_l, work);

               if (d->options & XRTM_OPTION_REVERSE_DERIVS) {
                    copy_polys(d, save->polys_up[i][j], polys);

                    for (k = 0; k < d->n_layers; ++k)
                         dvec_copy(save->P_full_up[i][j][k], P_full[k], d->n_stokes);
               }

               if (d->options & XRTM_OPTION_CALC_DERIVS) {
                    for (k = 0; k < d->n_derivs; ++k) {
                         for (l = 0; l < d->n_stokes; ++l) {
                              In_p_l[k][l] = d->In_p_l[0][k][ii + l];
                         }
                    }
               }

               single_scattered_radiance_up(d->n_stokes, d->n_derivs, d->n_layers,
                                            d->F_0,
                                            d->n_ulevels, d->ulevels, d->utaus2,
                                            1, &mus2[i],
                                            omega2, omega2_l,
                                            d->ltau, d->ltau_l,
                                            P_full, P_full_l,
                                            d->btran, d->btran_l,
                                            d->as_0, d->as_0_l,
                                            d->atran, d->atran_l,
                                            d->In_p[0]+ii, In_p_l,
                                            I, I_l,
                                            d->options & XRTM_OPTION_OUTPUT_AT_TAUS,
                                            d->derivs.layers, d->derivs.beam,
                                            save_tree, work);
               for (k = 0; k < d->n_ulevels; ++k) {
                    dvec_copy(I_p[k][i][j], I[k], d->n_stokes);
                    for (l = 0; l < d->n_derivs; ++l) {
                         dvec_copy(I_p_l[k][l][i][j], I_l[k][l], d->n_stokes);
                    }
               }


               phase_polys_get(d, d->n_coef, d->mu_0, d->phi_0,  mus2[i], phis[i][j], polys);

               get_phase_func_full_all(d,  mus2[i], phis[i][j], polys, P_full, P_full_l, work);

               if (d->options & XRTM_OPTION_REVERSE_DERIVS) {
                    copy_polys(d, save->polys_dn[i][j], polys);

                    for (k = 0; k < d->n_layers; ++k)
                         dvec_copy(save->P_full_dn[i][j][k], P_full[k], d->n_stokes);
               }

               if (d->options & XRTM_OPTION_CALC_DERIVS) {
                    for (k = 0; k < d->n_derivs; ++k) {
                         for (l = 0; l < d->n_stokes; ++l) {
                              I1_m_l[k][l] = d->I1_m_l[0][k][ii + l];
                         }
                    }
               }

               single_scattered_radiance_dn(d->n_stokes, d->n_derivs, d->n_layers,
                                            d->F_0,
                                            d->n_ulevels, d->ulevels, d->utaus2,
                                            1, &mus2[i],
                                            omega2, omega2_l,
                                            d->ltau, d->ltau_l,
                                            P_full, P_full_l,
                                            d->btran, d->btran_l,
                                            d->as_0, d->as_0_l,
                                            d->atran, d->atran_l,
                                            d->I1_m[0]+ii, I1_m_l,
                                            I, I_l,
                                            d->options & XRTM_OPTION_OUTPUT_AT_TAUS,
                                            d->derivs.layers, d->derivs.beam,
                                            save_tree, work);
               for (k = 0; k < d->n_ulevels; ++k) {
                    dvec_copy(I_m[k][i][j], I[k], d->n_stokes);
                    for (l = 0; l < d->n_derivs; ++l) {
                         dvec_copy(I_m_l[k][l][i][j], I_l[k][l], d->n_stokes);
                    }
               }
          }

          ii += d->n_stokes;
     }

     if (d->options & XRTM_OPTION_REVERSE_DERIVS)
          save_tree_decode_i_j(&save_tree);

#ifdef USE_AD_FOR_TL_CALC_GET_SINGLE
     get_single_tl_with_ad(d, n_phis, phis, I_p, I_m, I_p_l, I_m_l, save_tree, work);
#endif
}



/*******************************************************************************
 *
 ******************************************************************************/
static int fourier_get_add_top_down(xrtm_data *d, int i_four, int solver, double **I_p, double **I_m, double ***I_p_l, double ***I_m_l, save_tree_data save_tree, work_data work) {

     uchar *derivs_adding_down;

     int i;

     int surface;

     double **Rs_qq;
     double ***Rs_qq_l;

     double **In_p_l;
     double **I1_m_l;

     double **I_p_l1;
     double **I_m_l1;

     double *S_m;

     double **S_m_l;

     double **R_p;

     double ***R_p_l;


     surface = solar_surface_for_i_four(d, i_four);


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (surface) {
          Rs_qq = get_work1(&work, WORK_DXX);

          if (d->options & XRTM_OPTION_CALC_DERIVS && flags_or(d->derivs.layers[d->n_layers], d->n_derivs))
               Rs_qq_l = get_work2(&work, WORK_DXX, WORK_DERIVS_V, d->derivs.layers[d->n_layers]);
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     update_diff_bound_input(d, i_four, save_tree, work);

     if (surface)
          get_brdf_mats_qq(d, i_four, Rs_qq, Rs_qq_l, work);


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (d->options & XRTM_OPTION_CALC_DERIVS) {
          In_p_l = d->In_p_l[i_four];
          I1_m_l = d->I1_m_l[i_four];
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if ( ! (d->options & XRTM_OPTION_CALC_DERIVS) || ! (d->options & XRTM_OPTION_STACK_REUSE_ADDING)) {
          if (get_total_BOA_R_S_U_V  (d, i_four, solver, d->I1_m[i_four], I1_m_l, &R_p, &S_m, &R_p_l, &S_m_l, surface, 1, save_tree, &work)) {
               fprintf(stderr, "ERROR: get_total_BOA_R_S_U_V()\n");
               return -1;
          }
     }
     else {
/*
          if (get_total_BOA_R_S_U_V_2(d, i_four, solver, d->I1_m[i_four], I1_m_l, &R_p, &S_m, &R_p_l, &S_m_l, surface, 1, &work)) {
               fprintf(stderr, "ERROR: get_total_BOA_R_S_U_V_2()\n");
               return -1;
          }
*/
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (! surface) {
          init_array1_d(I_p[0], d->n_quad_v_x, 0.);
          copy_array1_d(I_m[0], S_m, d->n_quad_v_x);

          if (d->options & XRTM_OPTION_CALC_DERIVS) {
               for (i = 0; i < d->n_derivs; ++i) {
                    init_array1_d(I_p_l[0][i], d->n_quad_v_x, 0.);
                    if (! (d->derivs.adding_down[d->n_layers][i] & DERIV_TOP_S))
                         init_array1_d(I_m_l[0][i], d->n_quad_v_x, 0.);
                    else
                         copy_array1_d(I_m_l[0][i], S_m_l[i], d->n_quad_v_x);
               }
          }
     }
     else {
          if (d->options & XRTM_OPTION_CALC_DERIVS) {
               I_p_l1 = I_p_l[0]; I_m_l1 = I_m_l[0];

               derivs_adding_down = d->derivs.adding_down[d->n_layers];
          }

          radiance_boa_ref(d->n_quad_v_x, d->n_derivs, R_p, S_m, R_p_l, S_m_l, Rs_qq, Rs_qq_l, d->In_p[i_four], In_p_l, I_p[0], I_p_l1, I_m[0], I_m_l1, derivs_adding_down, work);
     }


     return 0;
}



/*******************************************************************************
 *
 ******************************************************************************/
static int fourier_get_add_bottom_up(xrtm_data *d, int i_four, int solver, double **I_p, double **I_m, double ***I_p_l, double ***I_m_l, save_tree_data save_tree, work_data work) {

     uchar  *derivs_adding_up;

     int i;

     int surface;

     double **Rs_qq;
     double ***Rs_qq_l;

     double **In_p_l;
     double **I1_m_l;

     double **I_p_l0;

     double *S_p;

     double **S_p_l;

     double **R_m;

     double ***R_m_l;

     forward_save_fourier_get_add_bottom_up_data *save;


     surface = solar_surface_for_i_four(d, i_four);


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (d->options & XRTM_OPTION_REVERSE_DERIVS) {
          save_tree_encode_s(&save_tree, "fourier_get_add_bottom_up");

          if (save_tree_retrieve_data(&save_tree, forward_save_fourier_get_add_bottom_up_data, &save))
               forward_save_fourier_get_add_bottom_up_alloc(save, d->n_quad_v_x);
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (surface) {
          Rs_qq = get_work1(&work, WORK_DXX);

          if (d->options & XRTM_OPTION_CALC_DERIVS && flags_or(d->derivs.layers[d->n_layers], d->n_derivs))
               Rs_qq_l = get_work2(&work, WORK_DXX, WORK_DERIVS_V, d->derivs.layers[d->n_layers]);
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     update_diff_bound_input(d, i_four, save_tree, work);

     if (surface)
          get_brdf_mats_qq(d, i_four, Rs_qq, Rs_qq_l, work);


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (d->options & XRTM_OPTION_CALC_DERIVS) {
          In_p_l = d->In_p_l[i_four];
          I1_m_l = d->I1_m_l[i_four];
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if ( ! (d->options & XRTM_OPTION_CALC_DERIVS) || ! (d->options & XRTM_OPTION_STACK_REUSE_ADDING)) {
          if (get_total_TOA_R_S_U_V  (d, i_four, solver, Rs_qq, d->In_p[i_four], Rs_qq_l, In_p_l, &R_m, &S_p, &R_m_l, &S_p_l, surface, 1, save_tree, &work)) {
               fprintf(stderr, "ERROR: get_total_TOA_R_S_U_V()\n");
               return -1;
          }
     }
     else {
          if (get_total_TOA_R_S_U_V_2(d, i_four, solver, Rs_qq, d->In_p[i_four], Rs_qq_l, In_p_l, &R_m, &S_p, &R_m_l, &S_p_l, surface, 1, save_tree, &work)) {
               fprintf(stderr, "ERROR: get_total_TOA_R_S_U_V_2()\n");
               return -1;
          }
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (d->options & XRTM_OPTION_CALC_DERIVS) {
          I_p_l0 = I_p_l[0];

          if (i_four == 0)
               derivs_adding_up = d->derivs.adding_up[0];
          else
               derivs_adding_up = d->derivs.adding_up_m[0];
     }

     radiance_toa_ref(d->n_quad_v_x, d->n_derivs, R_m, S_p, R_m_l, S_p_l, d->I1_m[i_four], I1_m_l, I_p[0], I_p_l0, derivs_adding_up, work);

     dvec_copy(I_m[0], d->I1_m[i_four], d->n_quad_v_x);

     if (d->options & XRTM_OPTION_CALC_DERIVS) {
          for (i = 0; i < d->n_derivs; ++i) {
               dvec_copy(I_m_l[0][i], I1_m_l[i], d->n_quad_v_x);
          }
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (d->options & XRTM_OPTION_REVERSE_DERIVS) {
          dmat_copy(save->R_m, R_m, d->n_quad_v_x, d->n_quad_v_x);
          dvec_copy(save->S_p, S_p, d->n_quad_v_x);

          if (surface)
               dmat_copy(save->Rs_qq, Rs_qq, d->n_quad_v_x, d->n_quad_v_x);
     }


     return 0;
}



/*******************************************************************************
 *
 ******************************************************************************/
static int fourier_get_add_both_ways(xrtm_data *d, int i_four, int solver, double **I_p, double **I_m, double ***I_p_l, double ***I_m_l, save_tree_data save_tree, work_data work) {

     uchar *derivs_adding_down;

     int i;
     int j;

     int surface;

     double **Rs_qq;
     double ***Rs_qq_l;

     double **In_p_l;
     double **I1_m_l;

     double *I_p1;
     double *I_m1;

     double **I_p_l0;
     double **I_m_l0;
     double **I_p_l1;
     double **I_m_l1;

     double *S_p;
     double *S_m;

     double **S_p_l;
     double **S_m_l;

     double **R_p;
     double **T_p;
     double **R_m;
     double **T_m;

     double ***R_p_l;
     double ***T_p_l;
     double ***R_m_l;
     double ***T_m_l;

     forward_save_fourier_get_add_both_ways_data *save;


     surface = solar_surface_for_i_four(d, i_four);


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (d->options & XRTM_OPTION_REVERSE_DERIVS) {
          save_tree_encode_s(&save_tree, "fourier_get_add_both_ways");

          if (save_tree_retrieve_data(&save_tree, forward_save_fourier_get_add_both_ways_data, &save))
               forward_save_fourier_get_add_both_ways_alloc(save, d->n_quad_v_x);
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (surface) {
          Rs_qq = get_work1(&work, WORK_DXX);

          if (d->options & XRTM_OPTION_CALC_DERIVS && flags_or(d->derivs.layers[d->n_layers], d->n_derivs))
               Rs_qq_l = get_work2(&work, WORK_DXX, WORK_DERIVS_V, d->derivs.layers[d->n_layers]);
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if ( ! (d->options & XRTM_OPTION_CALC_DERIVS) ||
          ! (d->options & XRTM_OPTION_STACK_REUSE_ADDING)) {
          if (get_total_R_T_S_U_W_V  (d, i_four, solver, &R_p, &T_p, &R_m, &T_m, &S_p, &S_m, &R_p_l, &T_p_l, &R_m_l, &T_m_l, &S_p_l, &S_m_l, 1, save_tree, &work)) {
               fprintf(stderr, "ERROR: get_total_R_T_S_U_W_V()\n");
               return -1;
          }
     }
     else {
          if (get_total_R_T_S_U_W_V_2(d, i_four, solver, &R_p, &T_p, &R_m, &T_m, &S_p, &S_m, &R_p_l, &T_p_l, &R_m_l, &T_m_l, &S_p_l, &S_m_l, 1, save_tree, &work)) {
               fprintf(stderr, "ERROR: get_total_R_T_S_U_W_V_2()\n");
               return -1;
          }
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     update_diff_bound_input(d, i_four, save_tree, work);

     if (surface)
          get_brdf_mats_qq(d, i_four, Rs_qq, Rs_qq_l, work);


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (d->options & XRTM_OPTION_CALC_DERIVS) {
          In_p_l = d->In_p_l[i_four];
          I1_m_l = d->I1_m_l[i_four];
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (d->options & XRTM_OPTION_CALC_DERIVS)
          derivs_adding_down = d->derivs.adding_down[d->n_layers];

     if (! surface) {
          for (i = 0; i < d->n_ulevels; ++i) {
               if (d->ulevels[i] == 0) {
                    if (d->options & XRTM_OPTION_CALC_DERIVS)
                         I_p_l0 = I_p_l[i];

                    radiance_slab(d->n_quad_v_x, d->n_derivs, R_m, T_p, S_p, R_m_l, T_p_l, S_p_l, d->I1_m[i_four], d->In_p[i_four], I1_m_l, In_p_l, I_p[i], I_p_l0, derivs_adding_down, work);

                    dvec_copy(I_m[i], d->I1_m[i_four], d->n_quad_v_x);

                    if (d->options & XRTM_OPTION_CALC_DERIVS) {
                         for (j = 0; j < d->n_derivs; ++j) {
                              dvec_copy(I_m_l[i][j], d->I1_m_l[i_four][j], d->n_quad_v_x);
                         }
                    }
               }
               else
               if (d->ulevels[i] == d->n_layers) {
                    if (d->options & XRTM_OPTION_CALC_DERIVS)
                         I_m_l1 = I_m_l[i];

                    radiance_slab(d->n_quad_v_x, d->n_derivs, R_p, T_m, S_m, R_p_l, T_m_l, S_m_l, d->In_p[i_four], d->I1_m[i_four], In_p_l, I1_m_l, I_m[i], I_m_l1, derivs_adding_down, work);

                    dvec_copy(I_p[i], d->In_p[i_four], d->n_quad_v_x);

                    if (d->options & XRTM_OPTION_CALC_DERIVS) {
                         for (j = 0; j < d->n_derivs; ++j) {
                              dvec_copy(I_p_l[i][j], d->In_p_l[i_four][j], d->n_quad_v_x);
                         }
                    }
               }
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
               I_p1 = I_p[0];
               I_m1 = I_m[0];
               if (d->options & XRTM_OPTION_CALC_DERIVS) {
                    I_p_l1 = I_p_l[0];
                    I_m_l1 = I_m_l[0];
               }
          }
          else
          if (d->n_ulevels > 1 &&
              d->ulevels[1] == d->n_layers) {
               I_p1 = I_p[1];
               I_m1 = I_m[1];
               if (d->options & XRTM_OPTION_CALC_DERIVS) {
                    I_p_l1 = I_p_l[1];
                    I_m_l1 = I_m_l[1];
               }
          }
          else {
               I_p1 = get_work_d1(&work, d->n_quad_v_x);
               I_m1 = get_work_d1(&work, d->n_quad_v_x);
               if (d->options & XRTM_OPTION_CALC_DERIVS) {
                    I_p_l1 = get_work_d2(&work, d->n_derivs, d->n_quad_v_x);
                    I_m_l1 = get_work_d2(&work, d->n_derivs, d->n_quad_v_x);
               }
          }

          radiance_boa_all(d->n_quad_v_x, d->n_derivs, R_p, T_m, Rs_qq, S_m, R_p_l, T_m_l, Rs_qq_l, S_m_l, d->In_p[i_four], d->I1_m[i_four], In_p_l, I1_m_l, I_p1, I_m1, I_p_l1, I_m_l1, derivs_adding_down, save_tree, work);

          if (d->ulevels[0] == 0) {
               if (d->options & XRTM_OPTION_CALC_DERIVS) {
                    I_p_l0 = I_p_l[0];
                    I_m_l0 = I_m_l[0];
               }

               radiance_toa_all(d->n_quad_v_x, d->n_derivs, R_m, T_p, S_p, R_m_l, T_p_l, S_p_l, d->I1_m[i_four], I_p1, I1_m_l, I_p_l1, I_p[0], I_m[0], I_p_l0, I_m_l0, derivs_adding_down, work);
          }
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (d->options & XRTM_OPTION_REVERSE_DERIVS) {
          dmat_copy(save->R_p, R_p, d->n_quad_v_x, d->n_quad_v_x);
          dmat_copy(save->T_p, T_p, d->n_quad_v_x, d->n_quad_v_x);
          dmat_copy(save->R_m, R_m, d->n_quad_v_x, d->n_quad_v_x);
          dmat_copy(save->T_m, T_m, d->n_quad_v_x, d->n_quad_v_x);

          dvec_copy(save->S_p, S_p, d->n_quad_v_x);
          dvec_copy(save->S_m, S_m, d->n_quad_v_x);

          if (surface) {
               dmat_copy(save->Rs_qq, Rs_qq, d->n_quad_v_x, d->n_quad_v_x);

               dvec_copy(save->I_p1,  I_p1, d->n_quad_v_x);
               dvec_copy(save->I_m1,  I_m1, d->n_quad_v_x);
          }
     }

#ifdef USE_AD_FOR_TL_CALC_FOURIER_GET_ADD_BOTH_WAYS_A
     fourier_get_add_both_ways_tl_with_ad(d, i_four, solver, I_p, I_m, I_p_l, I_m_l, save_tree, work);
#endif

     return 0;
}



/*******************************************************************************
 *
 ******************************************************************************/
static int fourier_get_bvp(xrtm_data *d, int i_four, int solver, double **I_p, double **I_m, double ***I_p_l, double ***I_m_l, save_tree_data save_tree, work_data work) {

     int i;

     int surface;
     int solar;
     int thermal;

     int add_single_scattering;

     double **In_p_l;
     double **I1_m_l;

     double **P_q0_mm;
     double **P_q0_pm;

     double ***P_q0_mm_l;
     double ***P_q0_pm_l;

     double **P_u0_mm;
     double **P_u0_pm;

     double ***P_u0_mm_l;
     double ***P_u0_pm_l;

     double ***P_qq_pp;
     double ***P_qq_mp;
     double ***P_qq_mm;
     double ***P_qq_pm;

     double ****P_qq_pp_l;
     double ****P_qq_mp_l;
     double ****P_qq_mm_l;
     double ****P_qq_pm_l;

     double ***P_uq_pp;
     double ***P_uq_mp;
     double ***P_uq_mm;
     double ***P_uq_pm;

     double ****P_uq_pp_l;
     double ****P_uq_mp_l;
     double ****P_uq_mm_l;
     double ****P_uq_pm_l;

     double ***r_p;
     double ***t_p;
     double ***r_m;
     double ***t_m;

     double ****r_p_l;
     double ****t_p_l;
     double ****r_m_l;
     double ****t_m_l;

     double *Rs_q0;
     double **Rs_q0_l;

     double *Rs_u0;
     double **Rs_u0_l;

     double **Rs_qq;
     double ***Rs_qq_l;

     double **Rs_uq;
     double ***Rs_uq_l;


     forward_save_fourier_get_bvp_data *save;


     surface = solar_surface_for_i_four(d, i_four);

     thermal = d->options & XRTM_OPTION_SOURCE_THERMAL && i_four == 0;
     solar   = d->options & XRTM_OPTION_SOURCE_SOLAR;


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (d->options & XRTM_OPTION_REVERSE_DERIVS) {
          save_tree_encode_s(&save_tree, "fourier_get_bvp");

          if (save_tree_retrieve_data(&save_tree, forward_save_fourier_get_bvp_data, &save))
               forward_save_fourier_get_bvp_alloc(save, d->options, solver, surface, d->n_layers, d->n_quad_v_x, d->n_umus_v);
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (surface) {
          Rs_q0 = get_work1(&work, WORK_DX);
          Rs_qq = get_work1(&work, WORK_DXX);

          if (d->n_umus > 0) {
               Rs_u0 = get_work1(&work, WORK_DU);
               Rs_uq = get_work1(&work, WORK_DUX);
          }

          if (d->options & XRTM_OPTION_CALC_DERIVS) {
               Rs_q0_l = get_work2(&work, WORK_DX,  WORK_DERIVS_V, d->derivs.layers[d->n_layers]);
               Rs_qq_l = get_work2(&work, WORK_DXX, WORK_DERIVS_V, d->derivs.layers[d->n_layers]);

               if (d->n_umus > 0) {
                    Rs_u0_l = get_work2(&work, WORK_DU,  WORK_DERIVS_V, d->derivs.layers[d->n_layers]);
                    Rs_uq_l = get_work2(&work, WORK_DUX, WORK_DERIVS_V, d->derivs.layers[d->n_layers]);
               }
          }
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     update_utaus(d, work);

     update_opt_props_all(d, work);

     update_beam_params_all(d, work);
if (! (solver & XRTM_SOLVER_TWO_STREAM || solver & XRTM_SOLVER_FOUR_STREAM || solver & XRTM_SOLVER_SIX_STREAM)) {
if (d->options & XRTM_OPTION_SOURCE_SOLAR) {
     if (get_phase_vecs_q0_all(d, i_four, &P_q0_mm, &P_q0_pm, &P_q0_mm_l, &P_q0_pm_l, 1, &work)) {
          fprintf(stderr, "ERROR: get_phase_vecs_q0_all()\n");
          return -1;
     }
}
     if (solver & XRTM_SOLVER_SOS || (solver & XRTM_SOLVER_TWO_OS && ! (d->options & XRTM_OPTION_SFI))) {
          if (get_phase_mats_qq_all(d, i_four, &P_qq_pp, &P_qq_mp, &P_qq_mm, &P_qq_pm, &P_qq_pp_l, &P_qq_mp_l, &P_qq_mm_l, &P_qq_pm_l, 1, &work)) {
               fprintf(stderr, "ERROR: get_phase_mats_qq_all()\n");
               return -1;
          }
     }
if (d->options & XRTM_OPTION_SOURCE_SOLAR) {
     if (d->options & XRTM_OPTION_SFI && d->n_umus > 0) {
          if (get_phase_vecs_u0_all(d, i_four, &P_u0_mm, &P_u0_pm, &P_u0_mm_l, &P_u0_pm_l, 1, &work)) {
               fprintf(stderr, "ERROR: get_phase_vecs_u0_all()\n");
               return -1;
          }
     }
}
     if (d->options & XRTM_OPTION_SFI && d->n_umus > 0) {
          if (get_phase_mats_uq_all(d, i_four, &P_uq_pp, &P_uq_mp, &P_uq_mm, &P_uq_pm, &P_uq_pp_l, &P_uq_mp_l, &P_uq_mm_l, &P_uq_pm_l, 1, &work)) {
               fprintf(stderr, "ERROR: get_phase_mats_uq_all()\n");
               return -1;
          }
     }

     if (solver & XRTM_SOLVER_EIG_BVP || solver & XRTM_SOLVER_MEM_BVP) {
          if (get_local_r_t_u_w_all(d, i_four, &r_p, &t_p, &r_m, &t_m, &r_p_l, &t_p_l, &r_m_l, &t_m_l, 1, save_tree, &work)) {
               fprintf(stderr, "ERROR: get_local_r_t_u_w_all()\n");
               return -1;
          }
     }
}
if (! (d->solvers & XRTM_SOLVER_TWO_STREAM)) {
     update_diff_bound_input(d, i_four, save_tree, work);

     if (d->options & XRTM_OPTION_CALC_DERIVS) {
          In_p_l = d->In_p_l[i_four];
          I1_m_l = d->I1_m_l[i_four];
     }
}
     if (surface) {
          get_brdf_vecs_q0(d, i_four, Rs_q0, Rs_q0_l, work);
          get_brdf_mats_qq(d, i_four, Rs_qq, Rs_qq_l, work);
     }

     if (surface && d->options & XRTM_OPTION_SFI && d->n_umus > 0) {
          get_brdf_vecs_u0(d, i_four, Rs_u0, Rs_u0_l, work);
          get_brdf_mats_uq(d, i_four, Rs_uq, Rs_uq_l, work);
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     add_single_scattering = 1;

     if (d->misc_input.do_not_add_sfi_ss)
          add_single_scattering = 0;

     if (solver & XRTM_SOLVER_EIG_BVP) {
#ifndef USE_NEW_SFI_FOUR_CONV
          int add_single_scattering = 1;

          if (d->misc_input.do_not_add_sfi_ss)
               add_single_scattering = 0;
#else
          int add_single_scattering = (! (d->options & XRTM_OPTION_N_T_TMS)) ||
                                      (   d->options & XRTM_OPTION_N_T_TMS &&   (d->options & XRTM_OPTION_FOUR_CONV_OLD || d->fourier_tol == 0.)) ||
                                      (   d->options & XRTM_OPTION_N_T_TMS && ! (d->options & XRTM_OPTION_FOUR_CONV_OLD || d->fourier_tol == 0.) && ! (d->options & XRTM_OPTION_SFI));

          if (d->misc_input.do_not_add_sfi_ss)
               add_single_scattering = 0;
#endif
          rtm_eig_bvp(i_four,
                      d->n_quad,   d->n_stokes, d->n_derivs, d->n_layers,
                      d->qf, d->qx_v, d->qw_v, d->F_0, d->mu_0,
                      d->n_ulevels, d->ulevels, d->utaus2, d->n_umus, d->umus_v,
                      d->levels_b, d->levels_b_l,
                      d->omega, d->omega_l, d->ltau, d->ltau_l,
                      d->surface_b, d->surface_b_l,
                      d->btran, d->btran_l,
                      d->as_0, d->as_0_l, d->atran, d->atran_l,
                      P_q0_mm, P_q0_pm, P_u0_mm, P_u0_pm,
                      P_uq_pp, P_uq_mp, P_uq_mm, P_uq_pm,
                      r_p, t_p, r_m, t_m,
                      P_q0_mm_l, P_q0_pm_l, P_u0_mm_l, P_u0_pm_l,
                      P_uq_pp_l, P_uq_mp_l, P_uq_mm_l, P_uq_pm_l,
                      r_p_l, t_p_l, r_m_l, t_m_l,
                      Rs_qq, Rs_qq_l,
                      Rs_u0, Rs_u0_l, Rs_uq, Rs_uq_l,
                      d->I1_m[i_four], I1_m_l, d->In_p[i_four], In_p_l,
                      d->F_iso_top, d->F_iso_top_l, d->F_iso_bot, d->F_iso_bot_l,
                      I_p, I_m, I_p_l, I_m_l,
                      add_single_scattering,
                      d->options & XRTM_OPTION_PART_SOL_GREENS,
                      d->options & XRTM_OPTION_SFI,
                      solar, thermal, surface,
                      d->options & XRTM_OPTION_UPWELLING_OUTPUT, d->options & XRTM_OPTION_DOWNWELLING_OUTPUT,
                      d->options & XRTM_OPTION_OUTPUT_AT_TAUS, d->options & XRTM_OPTION_VECTOR,
                      d->misc_input.eigen_solver_gen_real, d->misc_input.eigen_solver_gen_complex,
                      &d->derivs, save_tree, work);
     }
     else
     if (solver & XRTM_SOLVER_MEM_BVP)
          rtm_mem_bvp(i_four,
                      d->n_quad, d->n_stokes, d->n_derivs, d->n_layers,
                      d->qf, d->qx_v, d->qw_v, d->F_0, d->mu_0,
                      d->n_ulevels, d->ulevels, d->utaus2, d->n_umus, d->umus_v,
                      d->omega, d->omega_l, d->ltau, d->ltau_l,
                      d->btau, d->btau_l, d->btran, d->btran_l,
                      d->as_0, d->as_0_l, d->atran, d->atran_l,
                      P_q0_mm, P_q0_pm, P_u0_mm, P_u0_pm,
                      P_uq_pp, P_uq_mp, P_uq_mm, P_uq_pm,
                      r_p, t_p, r_m, t_m,
                      P_q0_mm_l, P_q0_pm_l, P_u0_mm_l, P_u0_pm_l,
                      P_uq_pp_l, P_uq_mp_l, P_uq_mm_l, P_uq_pm_l,
                      r_p_l, t_p_l, r_m_l, t_m_l,
                      Rs_qq, Rs_qq_l,
                      Rs_u0, Rs_u0_l, Rs_uq, Rs_uq_l,
                      d->I1_m[i_four], I1_m_l, d->In_p[i_four], In_p_l,
                      I_p, I_m, I_p_l, I_m_l,
                      d->options & XRTM_OPTION_SFI, surface,
                      d->options & XRTM_OPTION_OUTPUT_AT_TAUS, d->options & XRTM_OPTION_VECTOR,
                      d->misc_input.eigen_solver_gen_real, d->misc_input.eigen_solver_gen_complex,
                      &d->derivs, save_tree, work);
     else
     if (solver & XRTM_SOLVER_SOS)
          rtm_sos    (i_four,
                      d->n_quad_x, d->n_stokes, d->n_derivs, d->n_layers,
                      d->qf, d->qx_v, d->qw_v, d->F_0, d->mu_0,
                      d->n_ulevels, d->ulevels, d->utaus2,
                      d->omega, d->omega_l, d->ltau, d->ltau_l,
                      d->btran, d->btran_l,
                      d->as_0, d->as_0_l, d->atran, d->atran_l,
                      P_q0_mm, P_q0_pm,
                      P_qq_pp, P_qq_mp, P_qq_mm, P_qq_pm,
                      P_q0_mm_l, P_q0_pm_l,
                      P_qq_pp_l, P_qq_mp_l, P_qq_mm_l, P_qq_pm_l,
                      Rs_qq, Rs_qq_l,
                      d->I1_m[i_four], I1_m_l, d->In_p[i_four], In_p_l,
                      I_p, I_m, I_p_l, I_m_l,
                      d->sos_max_os, d->sos_max_tau, d->sos_tol,
                      d->options & XRTM_OPTION_SFI, surface, d->options & XRTM_OPTION_OUTPUT_AT_TAUS, d->options & XRTM_OPTION_VECTOR,
                      &d->derivs, work);
     else
     if (solver & XRTM_SOLVER_TWO_OS) {
          if (! (d->options & XRTM_OPTION_SFI))
               rtm_two_os(i_four,
                          d->n_quad, d->n_stokes, d->n_derivs, d->n_layers,
                          d->qf, d->qx_v, d->qw_v, d->F_0, d->mu_0,
                          d->n_ulevels, d->ulevels, d->utaus2, d->n_umus, d->umus_v,
                          d->omega, d->omega_l,
                          d->ltau, d->ltau_l,
                          Rs_q0, Rs_q0_l, Rs_qq, Rs_qq_l,
                          d->btran, d->btran_l,
                          d->as_0, d->as_0_l,
                          d->atran, d->atran_l,
                          P_q0_mm, P_q0_pm,
                          P_qq_pp, P_qq_mp, P_qq_mm, P_qq_pm,
                          P_q0_mm_l, P_q0_pm_l,
                          P_qq_pp_l, P_qq_mp_l, P_qq_mm_l, P_qq_pm_l,
                          I_p, I_m, I_p_l, I_m_l,
                          d->options & XRTM_OPTION_SFI, surface, d->options & XRTM_OPTION_OUTPUT_AT_TAUS, d->options & XRTM_OPTION_VECTOR, 0,
                          &d->derivs, work);
          else
               rtm_two_os(i_four,
                          d->n_quad_x, d->n_stokes, d->n_derivs, d->n_layers,
                          d->qf, d->qx_v, d->qw_v, d->F_0, d->mu_0,
                          d->n_ulevels, d->ulevels, d->utaus2,
                          d->n_umus, d->umus_v,
                          d->omega, d->omega_l, d->ltau, d->ltau_l,
                          Rs_q0, Rs_q0_l, Rs_uq, Rs_uq_l,
                          d->btran, d->btran_l,
                          d->as_0, d->as_0_l, d->atran, d->atran_l,
                          P_q0_mm, P_q0_pm,
                          P_uq_pp, P_uq_mp, P_uq_mm, P_uq_pm,
                          P_q0_mm_l, P_q0_pm_l,
                          P_uq_pp_l, P_uq_mp_l, P_uq_mm_l, P_uq_pm_l,
                          I_p, I_m, I_p_l, I_m_l,
                          d->options & XRTM_OPTION_SFI, surface, d->options & XRTM_OPTION_OUTPUT_AT_TAUS, d->options & XRTM_OPTION_VECTOR, 1,
                          &d->derivs, work);
     }
     else
     if (solver & XRTM_SOLVER_TWO_STREAM) {
#ifdef USE_NEW_SFI_FOUR_CONV
          int add_single_scattering = (! (d->options & XRTM_OPTION_N_T_TMS)) ||
                                      (   d->options & XRTM_OPTION_N_T_TMS &&   (d->options & XRTM_OPTION_FOUR_CONV_OLD || d->fourier_tol == 0.)) ||
                                      (   d->options & XRTM_OPTION_N_T_TMS && ! (d->options & XRTM_OPTION_FOUR_CONV_OLD || d->fourier_tol == 0.) && ! (d->options & XRTM_OPTION_SFI));

          if (d->misc_input.do_not_add_sfi_ss)
               add_single_scattering = 0;
#endif
          rtm_two_stream(i_four,
                         d->n_derivs, d->n_layers,
                         d->qx[0], d->F_0, d->mu_0,
                         d->n_ulevels, d->ulevels, d->utaus2,
                         d->n_umus, d->umus,
                         d->coef, d->coef_l, d->omega, d->omega_l, d->ltau, d->ltau_l,
                         d->btran, d->btran_l,
                         d->as_0, d->as_0_l, d->atran, d->atran_l,
                         Rs_q0, Rs_q0_l, Rs_qq, Rs_qq_l,
                         Rs_u0, Rs_u0_l, Rs_uq, Rs_uq_l,
                         I_p, I_m, I_p_l, I_m_l,
                         add_single_scattering,
                         solar, thermal, surface,
                         d->options & XRTM_OPTION_UPWELLING_OUTPUT, d->options & XRTM_OPTION_DOWNWELLING_OUTPUT,
                         d->options & XRTM_OPTION_OUTPUT_AT_TAUS,
                         &d->derivs,
                         save_tree, work);
     }
     else
     if (solver & XRTM_SOLVER_FOUR_STREAM)
          rtm_four_stream(i_four,
                          d->n_derivs, d->n_layers,
                          d->qx, d->qw, d->F_0, d->mu_0,
                          d->n_ulevels, d->ulevels, d->utaus2,
                          d->n_umus, d->umus,
                          d->coef, d->coef_l, d->omega, d->omega_l, d->ltau, d->ltau_l,
                          d->btran, d->btran_l,
                          d->as_0, d->as_0_l, d->atran, d->atran_l,
                          Rs_q0, Rs_q0_l, Rs_qq, Rs_qq_l,
                          Rs_u0, Rs_u0_l, Rs_uq, Rs_uq_l,
                          d->I1_m[i_four], I1_m_l, d->In_p[i_four], In_p_l,
                          I_p, I_m, I_p_l, I_m_l,
                          add_single_scattering,
                          solar, thermal, surface,
                          d->options & XRTM_OPTION_UPWELLING_OUTPUT, d->options & XRTM_OPTION_DOWNWELLING_OUTPUT,
                          d->options & XRTM_OPTION_OUTPUT_AT_TAUS,
                          &d->derivs, save_tree, work);
     else
     if (solver & XRTM_SOLVER_SIX_STREAM)
          rtm_six_stream (i_four,
                          d->n_derivs, d->n_layers,
                          d->qx, d->qw, d->F_0, d->mu_0,
                          d->n_ulevels, d->ulevels, d->utaus2,
                          d->n_umus, d->umus,
                          d->coef, d->coef_l,
                          d->omega, d->omega_l, d->ltau, d->ltau_l,
                          d->btran, d->btran_l,
                          d->as_0, d->as_0_l, d->atran, d->atran_l,
                          Rs_q0, Rs_q0_l, Rs_qq, Rs_qq_l,
                          Rs_u0, Rs_u0_l, Rs_uq, Rs_uq_l,
                          d->I1_m[i_four], I1_m_l, d->In_p[i_four], In_p_l,
                          I_p, I_m, I_p_l, I_m_l,
                          add_single_scattering,
                          solar, thermal, surface,
                          d->options & XRTM_OPTION_UPWELLING_OUTPUT, d->options & XRTM_OPTION_DOWNWELLING_OUTPUT,
                          d->options & XRTM_OPTION_OUTPUT_AT_TAUS,
                          &d->derivs, save_tree, work);
#ifdef DEBUG
     else {
         fprintf(stderr, "ERROR: fourier_get_bvp(): end of if / else if\n");
         exit(1);
     }
#endif

     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (d->options & XRTM_OPTION_REVERSE_DERIVS) {
          for (i = 0; i < d->n_layers; ++i) {
               copy_array1_d(save->P_q0_mm[i], P_q0_mm[i], d->n_quad_v_x);
               copy_array1_d(save->P_q0_pm[i], P_q0_pm[i], d->n_quad_v_x);

               if (solver & XRTM_SOLVER_SOS || (solver & XRTM_SOLVER_TWO_OS && ! (d->options & XRTM_OPTION_SFI))) {
                    copy_array2_d(save->P_qq_pp[i], P_qq_pp[i], d->n_quad_v_x, d->n_quad_v_x);
                    copy_array2_d(save->P_qq_mp[i], P_qq_mp[i], d->n_quad_v_x, d->n_quad_v_x);

                    if (d->options & XRTM_OPTION_VECTOR) {
                         copy_array2_d(save->P_qq_mm[i], P_qq_mm[i], d->n_quad_v_x, d->n_quad_v_x);
                         copy_array2_d(save->P_qq_pm[i], P_qq_pm[i], d->n_quad_v_x, d->n_quad_v_x);
                    }
               }

               if (d->options & XRTM_OPTION_SFI && d->n_umus_v > 0) {
                    copy_array1_d(save->P_u0_mm[i], P_u0_mm[i], d->n_umus_v);
                    copy_array1_d(save->P_u0_pm[i], P_u0_pm[i], d->n_umus_v);

                    copy_array2_d(save->P_uq_pp[i], P_uq_pp[i], d->n_umus_v, d->n_quad_v_x);
                    copy_array2_d(save->P_uq_mp[i], P_uq_mp[i], d->n_umus_v, d->n_quad_v_x);

                    if (d->options & XRTM_OPTION_VECTOR) {
                         copy_array2_d(save->P_uq_mm[i], P_uq_mm[i], d->n_umus_v, d->n_quad_v_x);
                         copy_array2_d(save->P_uq_pm[i], P_uq_pm[i], d->n_umus_v, d->n_quad_v_x);
                    }
               }

               copy_array2_d(save->r_p[i], r_p[i], d->n_umus_v, d->n_quad_v_x);
               copy_array2_d(save->t_p[i], t_p[i], d->n_umus_v, d->n_quad_v_x);

               if (d->options & XRTM_OPTION_VECTOR) {
                    copy_array2_d(save->r_m[i], r_m[i], d->n_umus_v, d->n_quad_v_x);
                    copy_array2_d(save->t_m[i], t_m[i], d->n_umus_v, d->n_quad_v_x);
               }

          }

          if (surface) {
               copy_array1_d(save->Rs_q0, Rs_q0, d->n_quad_v_x);
               copy_array2_d(save->Rs_qq, Rs_qq, d->n_quad_v_x, d->n_quad_v_x);

               if (d->n_umus_v > 0) {
                    copy_array1_d(save->Rs_u0, Rs_u0, d->n_umus_v);
                    copy_array2_d(save->Rs_uq, Rs_uq, d->n_umus_v, d->n_quad_v_x);
               }
          }
     }

     return 0;
}



/*******************************************************************************
 *
 ******************************************************************************/
static int get_fourier(xrtm_data *d, int i_four, int solver, double **I_p, double **I_m, double ***I_p_l, double ***I_m_l, save_tree_data save_tree, work_data work) {

     int i;
     int j;


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (solver & XRTM_SOLVERS_ADDING) {
          if (d->options & XRTM_OPTION_TOP_DOWN_ADDING) {
               if (fourier_get_add_top_down(d, i_four, solver, I_p, I_m, I_p_l, I_m_l, save_tree, work)) {
                    fprintf(stderr, "ERROR: fourier_get_add_top_down()\n");
                    return -1;
               }
          }
          else
          if (d->options & XRTM_OPTION_BOTTOM_UP_ADDING) {
               if (fourier_get_add_bottom_up(d, i_four, solver, I_p, I_m, I_p_l, I_m_l, save_tree, work)) {
                    fprintf(stderr, "ERROR: fourier_get_add_bottom_up()\n");
                    return -1;
               }
          }
          else {
               if (fourier_get_add_both_ways(d, i_four, solver, I_p, I_m, I_p_l, I_m_l, save_tree, work)) {
                    fprintf(stderr, "ERROR: fourier_get_add_both_ways()\n");
                    return -1;
               }
          }
     }
     else if (solver & XRTM_SOLVERS_BVP) {
          if (fourier_get_bvp(d, i_four, solver, I_p, I_m, I_p_l, I_m_l, save_tree, work)) {
               fprintf(stderr, "ERROR: fourier_get_bvp()\n");
               return -1;
          }
     }
#ifdef DEBUG
     else {
         fprintf(stderr, "ERROR: get_fourier(): end of if / else if\n");
         exit(1);
     }
#endif

     /*--------------------------------------------------------------------
      *
      *------------------------------------------------------------------*/
/*
     if (d->options & XRTM_OPTION_VECTOR) {
          for (i = 0; i < d->n_ulevels; ++i) {
               dm_v_mul_D_A_2(d->n_quad_d, d->n_stokes, I_m[i], I_m[i]);

               if (d->options & XRTM_OPTION_CALC_DERIVS) {
                    for (j = 0; j < d->n_derivs; ++j) {
                         dm_v_mul_D_A_2(d->n_quad_d, d->n_stokes, I_m_l[i][j], I_m_l[i][j]);
                    }
               }
          }
     }
*/

     /*--------------------------------------------------------------------
      *
      *------------------------------------------------------------------*/
     if (d->misc_input.use_symmetric_form) {
          for (i = 0; i < d->n_ulevels; ++i) {
               vec_sim_trans2(d->n_quad_v_d, I_p[i], I_m[i], d->alpha2);

               if (d->options & XRTM_OPTION_CALC_DERIVS) {
                    for (j = 0; j < d->n_derivs; ++j) {
                         vec_sim_trans2(d->n_quad_v_d, I_p_l[i][j], I_m_l[i][j], d->alpha2);
                    }
               }
          }
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
if (0) {
     if (i_four == 0) {
          update_beam_params_all(d, work);

          if (check_conserve_energy(d->n_quad, d->n_stokes,
                                    d->qx, d->qw, d->F_0, d->mu_0,
                                    d->btran[d->n_layers], I_p[0], I_m[1]) > 0) {
               fprintf(stderr, "ERROR: check_conserve_energy()\n");
               return -1;
          }

          if (d->options & XRTM_OPTION_CALC_DERIVS) {
               for (i = 0; i < d->n_derivs; ++i) {
                    if (check_conserve_energy_l(d->n_quad, d->n_stokes,
                                                d->qx, d->qw, d->F_0, d->mu_0,
                                                d->btran_l[d->n_layers][i], I_p_l[0][i], I_m_l[1][i]) > 0) {
                         fprintf(stderr, "ERROR: check_conserve_energy_l(): i_derivs = %d\n", i);
                         return -1;
                    }
               }
          }
     }
}

     return 0;
}



/*******************************************************************************
 *
 ******************************************************************************/
static int apply_corrections_radiance(xrtm_data *d, int n_phis, double **phis, double ****I_p, double ****I_m, double *****I_p_l, double *****I_m_l, save_tree_data save_tree, work_data work) {

     int i;
     int ii;
     int j;
     int k;
     int l;

     int i_mus2_v;
     int n_mus2;

     double a;

     double *mus2;

     double f;
     double *f_l = NULL;

     double *omega_tms;
     double **omega_tms_l;

     double **In_p_l;
     double **I1_m_l;

     double **I;
     double ***I_l;

     double **P_trun;
     double **P_full;
     double ***P_trun_l;
     double ***P_full_l;

     void *polys;

     forward_save_apply_corrections_radiance_data *save;


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (d->options & XRTM_OPTION_N_T_TMS && d->options & XRTM_OPTION_SOURCE_SOLAR) {

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
          if (d->options & XRTM_OPTION_REVERSE_DERIVS) {
               save_tree_encode_s(&save_tree, "apply_corrections_radiance");

               if (save_tree_retrieve_data(&save_tree, forward_save_apply_corrections_radiance_data, &save))
                    forward_save_apply_corrections_radiance_alloc(save, d->options, d->n_coef, d->n_quad_v + d->n_umus_v, n_mus2, n_phis, d->n_layers, d->n_stokes);
          }


          /*--------------------------------------------------------------------
           *
           *------------------------------------------------------------------*/
          if (d->options & XRTM_OPTION_CALC_DERIVS)
               f_l = get_work1(&work, WORK_DDERIVS);

          omega_tms = get_work1(&work, WORK_DLAYERS);

          if (d->options & XRTM_OPTION_CALC_DERIVS)
               omega_tms_l = get_work1(&work, WORK_DBOTH);

          if (d->options & XRTM_OPTION_PHASE_SCALAR)
               polys = get_work_d1(&work, d->n_coef);
          else
          if (d->options & XRTM_OPTION_PHASE_MATRIX_GC)
               polys = get_work_d2(&work, 4, d->n_coef);
          else
          if (d->options & XRTM_OPTION_PHASE_MATRIX_LC)
               polys = get_work_d1(&work, d->n_coef);
#ifdef DEBUG
          else {
              fprintf(stderr, "ERROR: get_phase_func(): end of if / else if\n");
              exit(1);
          }
#endif
          P_trun = get_work2(&work, WORK_DS, WORK_LAYERS_V, NULL);
          P_full = get_work2(&work, WORK_DS, WORK_LAYERS_V, NULL);

          I      = get_work_d2(&work, d->n_ulevels, d->n_stokes);

          if (d->options & XRTM_OPTION_CALC_DERIVS) {
               P_trun_l = get_work3(&work, WORK_DS, WORK_BOTH_V, d->derivs.layers);
               P_full_l = get_work3(&work, WORK_DS, WORK_BOTH_V, d->derivs.layers);

               In_p_l   = get_work2(&work, WORK_DS, WORK_DERIVS_V, NULL);
               I1_m_l   = get_work2(&work, WORK_DS, WORK_DERIVS_V, NULL);

               I_l      = get_work_d3(&work, d->n_ulevels, d->n_derivs, d->n_stokes);
          }


          /*--------------------------------------------------------------------
           *
           *------------------------------------------------------------------*/
          update_utaus(d, work);

          update_opt_props_all(d, work);

          update_beam_params_all(d, work);

          update_diff_bound_input(d, 0, save_tree, work);


          a = 2. * d->n_coef2 + 1.;

          for (i = 0; i < d->n_layers; ++i) {
               get_delta_m_f(d, i, 0, d->n_derivs, a, &f, f_l);

               if (! (d->options & XRTM_OPTION_CALC_DERIVS))
                    n_t_tms_scaling(d->n_derivs, f, f_l, d->omega0[i], &omega_tms[i], NULL, NULL);
               else
                    n_t_tms_scaling(d->n_derivs, f, f_l, d->omega0[i], &omega_tms[i], d->omega0_l[i], omega_tms_l[i]);
          }


          /*--------------------------------------------------------------------
           *
           *------------------------------------------------------------------*/
          ii = i_mus2_v;

          for (i = 0; i < n_mus2; ++i) {
               for (j = 0; j < n_phis; ++j) {

               if (d->options & XRTM_OPTION_REVERSE_DERIVS)
                    save_tree_recode_i_j(&save_tree, i, j, i == 0 && j == 0);


                    if (d->options & XRTM_OPTION_UPWELLING_OUTPUT) {
                         phase_polys_get(d, d->n_coef, d->mu_0, d->phi_0, -mus2[i], phis[i][j], polys);

                         if (d->options & XRTM_OPTION_REVERSE_DERIVS)
                              copy_polys(d, save->polys_up[i][j], polys);

                         if (d->options & XRTM_OPTION_CALC_DERIVS) {
                              for (k = 0; k < d->n_derivs; ++k) {
                                   for (l = 0; l < d->n_stokes; ++l) {
                                        In_p_l[k][l] = d->In_p_l[0][k][ii + l];
                                   }
                              }
                         }

                         if (d->options & XRTM_OPTION_FOUR_CONV_OLD || d->fourier_tol == 0.) {
                              get_phase_func_trun_all(d, -mus2[i], phis[i][j], polys, P_trun, P_trun_l, work);
                              get_phase_func_full_all(d, -mus2[i], phis[i][j], polys, P_full, P_full_l, work);

                              if (d->options & XRTM_OPTION_REVERSE_DERIVS) {
                                   for (k = 0; k < d->n_layers; ++k) {
                                        dvec_copy(save->P_trun_up[i][j][k], P_trun[k], d->n_stokes);
                                        dvec_copy(save->P_full_up[i][j][k], P_full[k], d->n_stokes);
                                   }
                              }

                              n_t_tms_correction_up(d->n_stokes, d->n_derivs, d->n_layers,
                                                    d->F_0,
                                                    d->n_ulevels, d->ulevels, d->utaus2,
                                                    1, &mus2[i],
                                                    d->omega, d->omega_l,
                                                    omega_tms, omega_tms_l,
                                                    d->ltau, d->ltau_l,
                                                    P_trun, P_full,
                                                    P_trun_l, P_full_l,
                                                    d->btran, d->btran_l,
                                                    d->as_0, d->as_0_l,
                                                    d->atran, d->atran_l,
                                                    d->In_p[0]+ii, In_p_l,
                                                    I, I_l,
                                                    d->options & XRTM_OPTION_OUTPUT_AT_TAUS,
                                                    d->derivs.layers, d->derivs.beam,
                                                    save_tree, work);
                         }
                         else {
                              get_phase_func_full_all(d, -mus2[i], phis[i][j], polys, P_full, P_full_l, work);

                              if (d->options & XRTM_OPTION_REVERSE_DERIVS) {
                                   for (k = 0; k < d->n_layers; ++k)
                                        dvec_copy(save->P_full_up[i][j][k], P_full[k], d->n_stokes);
                              }

                              single_scattered_radiance_up(d->n_stokes, d->n_derivs, d->n_layers,
                                                           d->F_0,
                                                           d->n_ulevels, d->ulevels, d->utaus2,
                                                           1, &mus2[i],
                                                           omega_tms, omega_tms_l,
                                                           d->ltau, d->ltau_l,
                                                           P_full, P_full_l,
                                                           d->btran, d->btran_l,
                                                           d->as_0, d->as_0_l,
                                                           d->atran, d->atran_l,
                                                           d->In_p[0]+ii, In_p_l,
                                                           I, I_l,
                                                           d->options & XRTM_OPTION_OUTPUT_AT_TAUS,
                                                           d->derivs.layers, d->derivs.beam,
                                                           save_tree, work);
                         }

                         for (k = 0; k < d->n_ulevels; ++k) {
                              dvec_add(I_p[k][i][j], I[k], I_p[k][i][j], d->n_stokes);
                              if (d->options & XRTM_OPTION_CALC_DERIVS) {
                                   for (l = 0; l < d->n_derivs; ++l) {
                                        dvec_add(I_p_l[k][l][i][j], I_l[k][l], I_p_l[k][l][i][j], d->n_stokes);
                                   }
                              }
                         }
                    }


                    if (d->options & XRTM_OPTION_DOWNWELLING_OUTPUT) {
                         phase_polys_get(d, d->n_coef, d->mu_0, d->phi_0,  mus2[i], phis[i][j], polys);
                         if (d->options & XRTM_OPTION_REVERSE_DERIVS)
                              copy_polys(d, save->polys_dn[i][j], polys);

                         if (d->options & XRTM_OPTION_CALC_DERIVS) {
                              for (k = 0; k < d->n_derivs; ++k) {
                                   for (l = 0; l < d->n_stokes; ++l) {
                                        I1_m_l[k][l] = d->I1_m_l[0][k][ii + l];
                                   }
                              }
                         }

                         if (d->options & XRTM_OPTION_FOUR_CONV_OLD || d->fourier_tol == 0.) {
                              get_phase_func_trun_all(d,  mus2[i], phis[i][j], polys, P_trun, P_trun_l, work);
                              get_phase_func_full_all(d,  mus2[i], phis[i][j], polys, P_full, P_full_l, work);

                              if (d->options & XRTM_OPTION_REVERSE_DERIVS) {
                                   for (k = 0; k < d->n_layers; ++k) {
                                        dvec_copy(save->P_trun_dn[i][j][k], P_trun[k], d->n_stokes);
                                        dvec_copy(save->P_full_dn[i][j][k], P_full[k], d->n_stokes);
                                   }
                              }

                              n_t_tms_correction_dn(d->n_stokes, d->n_derivs, d->n_layers,
                                                    d->F_0,
                                                    d->n_ulevels, d->ulevels, d->utaus2,
                                                    1., &mus2[i],
                                                    d->omega, d->omega_l,
                                                    omega_tms, omega_tms_l,
                                                    d->ltau, d->ltau_l,
                                                    P_trun, P_full,
                                                    P_trun_l, P_full_l,
                                                    d->btran, d->btran_l,
                                                    d->as_0, d->as_0_l,
                                                    d->atran, d->atran_l,
                                                    d->I1_m[0]+ii, I1_m_l,
                                                    I, I_l,
                                                    d->options & XRTM_OPTION_OUTPUT_AT_TAUS,
                                                    d->derivs.layers, d->derivs.beam,
                                                    save_tree, work);
                         }
                         else {
                              get_phase_func_full_all(d,  mus2[i], phis[i][j], polys, P_full, P_full_l, work);

                              if (d->options & XRTM_OPTION_REVERSE_DERIVS) {
                                   for (k = 0; k < d->n_layers; ++k)
                                        dvec_copy(save->P_full_dn[i][j][k], P_full[k], d->n_stokes);
                              }

                              single_scattered_radiance_dn(d->n_stokes, d->n_derivs, d->n_layers,
                                                           d->F_0,
                                                           d->n_ulevels, d->ulevels, d->utaus2,
                                                           1, &mus2[i],
                                                           omega_tms, omega_tms_l,
                                                           d->ltau, d->ltau_l,
                                                           P_full, P_full_l,
                                                           d->btran, d->btran_l,
                                                           d->as_0, d->as_0_l,
                                                           d->atran, d->atran_l,
                                                           d->I1_m[0]+ii, I1_m_l,
                                                           I, I_l,
                                                           d->options & XRTM_OPTION_OUTPUT_AT_TAUS,
                                                           d->derivs.layers, d->derivs.beam,
                                                           save_tree, work);
                         }

                         for (k = 0; k < d->n_ulevels; ++k) {
                              dvec_add(I_m[k][i][j], I[k], I_m[k][i][j], d->n_stokes);
                              if (d->options & XRTM_OPTION_CALC_DERIVS) {
                                   for (l = 0; l < d->n_derivs; ++l) {
                                        dvec_add(I_m_l[k][l][i][j], I_l[k][l], I_m_l[k][l][i][j], d->n_stokes);
                                   }
                              }
                         }
                    }
               }

               ii += d->n_stokes;
          }
     }


     return 0;
}



/*******************************************************************************
 *
 ******************************************************************************/
static int apply_corrections_fourier(xrtm_data *d, int i_four, double **I_p, double **I_m, double ***I_p_l, double ***I_m_l, save_tree_data save_tree, work_data work) {

     int i;
     int k;

     double **m1;

     double **I;
     double ***I_l;

     double **P_q0_mm;
     double **P_q0_pm;

     double ***P_q0_mm_l;
     double ***P_q0_pm_l;

     double **P_u0_mm;
     double **P_u0_pm;

     double ***P_u0_mm_l;
     double ***P_u0_pm_l;

     double **In_p_l;
     double **I1_m_l;

     forward_save_apply_corrections_fourier_data *save;


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
#ifndef USE_NEW_SFI_FOUR_CONV
     if (d->options & XRTM_OPTION_N_T_TMS && d->options & XRTM_OPTION_SOURCE_SOLAR && (d->options & XRTM_OPTION_FOUR_CONV_NEW && d->fourier_tol > 0.)                                    ) {
#else
     if (d->options & XRTM_OPTION_N_T_TMS && d->options & XRTM_OPTION_SOURCE_SOLAR && (d->options & XRTM_OPTION_FOUR_CONV_NEW && d->fourier_tol > 0.) && ! (d->options & XRTM_OPTION_SFI)) {
#endif
          if (d->options & XRTM_OPTION_REVERSE_DERIVS) {
               save_tree_encode_s(&save_tree, "apply_corrections_fourier");

               if (save_tree_retrieve_data(&save_tree, forward_save_apply_corrections_fourier_data, &save))
                    forward_save_apply_corrections_fourier_alloc(save, d->n_quad_v, d->n_quad_v_x, d->n_umus_v, d->n_layers);
          }


          /*--------------------------------------------------------------------
           *
           *------------------------------------------------------------------*/
          m1 = work_get2(&work, WORK_DU, WORK_DERIVS_V, NULL);

          I  = get_work_d2(&work, d->n_ulevels, d->n_quad_v + d->n_umus_v);

          if (d->options & XRTM_OPTION_CALC_DERIVS)
               I_l = get_work_d3(&work, d->n_ulevels, d->n_derivs, d->n_quad_v + d->n_umus_v);


          /*--------------------------------------------------------------------
           *
           *------------------------------------------------------------------*/
          update_opt_props_all(d, work);

          update_beam_params_all(d, work);

          if (get_phase_vecs_q0_all(d, i_four, &P_q0_mm, &P_q0_pm, &P_q0_mm_l, &P_q0_pm_l, 1, &work)) {
               fprintf(stderr, "ERROR: get_phase_vecs_q0_all()\n");
               return -1;
          }

          if (d->options & XRTM_OPTION_SFI && d->n_umus > 0) {
               if (get_phase_vecs_u0_all(d, i_four, &P_u0_mm, &P_u0_pm, &P_u0_mm_l, &P_u0_pm_l, 1, &work)) {
                    fprintf(stderr, "ERROR: get_phase_vecs_u0_all()\n");
                    return -1;
               }
          }

          update_diff_bound_input(d, i_four, save_tree, work);

          if (d->options & XRTM_OPTION_CALC_DERIVS) {
               In_p_l = d->In_p_l[i_four];
               I1_m_l = d->I1_m_l[i_four];
          }


          /*--------------------------------------------------------------------
           *
           *------------------------------------------------------------------*/
          if (d->options & XRTM_OPTION_REVERSE_DERIVS) {
               for (k = 0; k < d->n_layers; ++k) {
                    dvec_copy(save->P_q0_mm[k], P_q0_mm[k], d->n_quad_v_x);
                    dvec_copy(save->P_q0_pm[k], P_q0_pm[k], d->n_quad_v_x);

               }
               if (d->options & XRTM_OPTION_SFI && d->n_umus > 0) {
                    for (k = 0; k < d->n_layers; ++k) {
                         dvec_copy(save->P_u0_mm[k], P_u0_mm[k], d->n_umus_v);
                         dvec_copy(save->P_u0_pm[k], P_u0_pm[k], d->n_umus_v);
                    }
               }

               dvec_copy(save->In_p, d->In_p[i_four], d->n_quad_v + d->n_umus_v);
               dvec_copy(save->I1_m, d->I1_m[i_four], d->n_quad_v + d->n_umus_v);
          }


          /*--------------------------------------------------------------------
           *
           *------------------------------------------------------------------*/
          if (d->options & XRTM_OPTION_UPWELLING_OUTPUT) {
               single_scattered_radiance_up(d->n_stokes, d->n_derivs, d->n_layers,
                                            d->F_0,
                                            d->n_ulevels, d->ulevels, d->utaus2,
                                            d->n_quad_x, d->qx,
                                            d->omega, d->omega_l,
                                            d->ltau, d->ltau_l,
                                            P_q0_pm, P_q0_pm_l,
                                            d->btran, d->btran_l,
                                            d->as_0, d->as_0_l,
                                            d->atran, d->atran_l,
                                            d->In_p[i_four], In_p_l,
                                            I, I_l,
                                            d->options & XRTM_OPTION_OUTPUT_AT_TAUS,
                                            d->derivs.layers, d->derivs.beam,
                                            save_tree, work);

               for (k = 0; k < d->n_ulevels; ++k) {
                    dvec_sub(I_p[k], I[k], I_p[k], d->n_quad_v_x);
                    if (d->options & XRTM_OPTION_CALC_DERIVS) {
                         for (i = 0; i < d->n_derivs; ++i) {
                              dvec_sub(I_p_l[k][i], I_l[k][i], I_p_l[k][i], d->n_quad_v_x);
                         }
                    }
               }

               if (d->options & XRTM_OPTION_SFI && d->n_umus > 0) {
                    if (d->options & XRTM_OPTION_CALC_DERIVS) {
                         for (k = 0; k < d->n_derivs; ++k)
                              dvec_copy(m1[k], In_p_l[k] + d->n_quad_v, d->n_umus_v);
                    }

                    single_scattered_radiance_up(d->n_stokes, d->n_derivs, d->n_layers,
                                                 d->F_0,
                                                 d->n_ulevels, d->ulevels, d->utaus2,
                                                 d->n_umus, d->umus,
                                                 d->omega, d->omega_l,
                                                 d->ltau, d->ltau_l,
                                                 P_u0_pm, P_u0_pm_l,
                                                 d->btran, d->btran_l,
                                                 d->as_0, d->as_0_l,
                                                 d->atran, d->atran_l,
                                                 d->In_p[i_four] + d->n_quad_v, m1,
                                                 I, I_l, d->options & XRTM_OPTION_OUTPUT_AT_TAUS,
                                                 d->derivs.layers, d->derivs.beam,
                                                 save_tree, work);

                    for (k = 0; k < d->n_ulevels; ++k) {
                         dvec_sub(I_p[k]+d->n_quad_v, I[k], I_p[k]+d->n_quad_v, d->n_umus_v);
                         if (d->options & XRTM_OPTION_CALC_DERIVS) {
                              for (i = 0; i < d->n_derivs; ++i) {
                                   dvec_sub(I_p_l[k][i]+d->n_quad_v, I_l[k][i], I_p_l[k][i]+d->n_quad_v, d->n_umus_v);
                              }
                         }
                    }
               }
          }

          if (d->options & XRTM_OPTION_DOWNWELLING_OUTPUT) {
               single_scattered_radiance_dn(d->n_stokes, d->n_derivs, d->n_layers,
                                            d->F_0, d->n_ulevels, d->ulevels, d->utaus2,
                                            d->n_quad_x, d->qx,
                                            d->omega, d->omega_l,
                                            d->ltau, d->ltau_l,
                                            P_q0_mm, P_q0_mm_l,
                                            d->btran, d->btran_l,
                                            d->as_0, d->as_0_l,
                                            d->atran, d->atran_l,
                                            d->I1_m[i_four], I1_m_l,
                                            I, I_l,
                                            d->options & XRTM_OPTION_OUTPUT_AT_TAUS,
                                            d->derivs.layers, d->derivs.beam,
                                            save_tree, work);

               for (k = 0; k < d->n_ulevels; ++k) {
                    dvec_sub(I_m[k], I[k], I_m[k], d->n_quad_v_x);
                    if (d->options & XRTM_OPTION_CALC_DERIVS) {
                         for (i = 0; i < d->n_derivs; ++i) {
                              dvec_sub(I_m_l[k][i], I_l[k][i], I_m_l[k][i], d->n_quad_v_x);
                         }
                    }
               }

               if (d->options & XRTM_OPTION_SFI && d->n_umus > 0) {
                    if (d->options & XRTM_OPTION_CALC_DERIVS) {
                         for (k = 0; k < d->n_derivs; ++k)
                              dvec_copy(m1[k], I1_m_l[k] + d->n_quad_v, d->n_umus_v);
                    }

                    single_scattered_radiance_dn(d->n_stokes, d->n_derivs, d->n_layers,
                                                 d->F_0,
                                                 d->n_ulevels, d->ulevels, d->utaus2,
                                                 d->n_umus, d->umus,
                                                 d->omega, d->omega_l,
                                                 d->ltau, d->ltau_l,
                                                 P_u0_mm, P_u0_mm_l,
                                                 d->btran, d->btran_l,
                                                 d->as_0, d->as_0_l,
                                                 d->atran, d->atran_l,
                                                 d->I1_m[i_four] + d->n_quad_v, m1,
                                                 I, I_l,
                                                 d->options & XRTM_OPTION_OUTPUT_AT_TAUS,
                                                 d->derivs.layers, d->derivs.beam,
                                                 save_tree, work);

                    for (k = 0; k < d->n_ulevels; ++k) {
                         dvec_sub(I_m[k]+d->n_quad_v, I[k], I_m[k]+d->n_quad_v, d->n_umus_v);
                         if (d->options & XRTM_OPTION_CALC_DERIVS) {
                              for (i = 0; i < d->n_derivs; ++i) {
                                   dvec_sub(I_m_l[k][i]+d->n_quad_v, I_l[k][i], I_m_l[k][i]+d->n_quad_v, d->n_umus_v);
                              }
                         }
                    }
               }
          }
     }


     return 0;
}



/*******************************************************************************
 *
 ******************************************************************************/
static void azimuth_modes(int n_quad, int n_stokes,
                          double a, double b, double *v1, double *v2) {

     int i;
     int ii;

     ii = 0;
     if (n_stokes == 1) {
          for (i = 0; i < n_quad; ++i) {
               v2[ii] = v1[ii] * a; ii++;
          }
     }
     else
     if (n_stokes == 2) {
          for (i = 0; i < n_quad; ++i) {
               v2[ii] = v1[ii] * a; ii++;
               v2[ii] = v1[ii] * a; ii++;
          }
     }
     else
     if (n_stokes == 3) {
          for (i = 0; i < n_quad; ++i) {
               v2[ii] = v1[ii] * a; ii++;
               v2[ii] = v1[ii] * a; ii++;
               v2[ii] = v1[ii] * b; ii++;
          }
     }
     else {
          for (i = 0; i < n_quad; ++i) {
               v2[ii] = v1[ii] * a; ii++;
               v2[ii] = v1[ii] * a; ii++;
               v2[ii] = v1[ii] * b; ii++;
               v2[ii] = v1[ii] * b; ii++;
          }
     }
}



static double max_radiance_ratio(double I1, double I2) {

     double a;

     a = 0.;

     if (I1 != 0.)
          a = MAX(a, fabs(I2 / I1));

     return a;
}



static int get_solution_internal(xrtm_data *d, enum xrtm_solver_mask solver, int solutions, int n_phis, double **phis, double ****I_p, double ****I_m, double *****I_p_l, double *****I_m_l, double *mean_p, double *mean_m, double **mean_p_l, double **mean_m_l, double *flux_p, double *flux_m, double **flux_p_l, double **flux_m_l, double *flux_div, double **flux_div_l, save_tree_data save_tree, work_data work) {

     int i;
     int j;
     int jj;
     int k;
     int kk;
     int l;
     int m;

     int count;

     int n_four;
/*
     int n_quad_x;
*/
     int n_quad_v_x;

     int i_mus2;
     int n_mus2;

     double a;
     double b;
     double c;
     double e;

     double *v1;
     double *v2;

     double **i_p;
     double **i_m;
     double ***i_p_l;
     double ***i_m_l;

     forward_save_get_solution_internal_data *save;


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
/*
     n_quad_x   = d->n_quad   + d->n_umus;
*/
     n_quad_v_x = d->n_quad_v + d->n_umus_v;

     if (d->n_umus == 0) {
          i_mus2  = 0;
          n_mus2  = d->n_quad;
     }
     else {
          i_mus2  = d->n_quad;
          n_mus2  = d->n_umus;
     }
/*
     else {
          i_mus2  = 0;
          n_mus2  = n_quad_x;
     }
*/

     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (solutions & XRTM_OUTPUT_RADIANCE)
          n_four = d->n_four2;
     else
          n_four = 1;


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (d->options & XRTM_OPTION_REVERSE_DERIVS) {
          save_tree_encode_s(&save_tree, "get_solution_internal");

          if (save_tree_retrieve_data(&save_tree, forward_save_get_solution_internal_data, &save))
               forward_save_get_solution_internal_alloc(save, n_four, d->n_ulevels, n_quad_v_x);
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (solver & XRTM_SOLVER_SINGLE) {
          get_single(d, n_phis, phis, I_p, I_m, I_p_l, I_m_l, save_tree, work);
          return 0;
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     v1 = get_work1(&work, WORK_DD);
     v2 = get_work1(&work, WORK_DD);

     i_p = get_work_d2(&work, d->n_ulevels, n_quad_v_x);
     i_m = get_work_d2(&work, d->n_ulevels, n_quad_v_x);

     if (d->options & XRTM_OPTION_CALC_DERIVS) {
          i_p_l = get_work_d3(&work, d->n_ulevels, d->n_derivs, n_quad_v_x);
          i_m_l = get_work_d3(&work, d->n_ulevels, d->n_derivs, n_quad_v_x);
     }

     if (solutions & XRTM_OUTPUT_RADIANCE) {
          init_array4_d(I_p, d->n_ulevels, n_mus2, n_phis, d->n_stokes, 0.);
          init_array4_d(I_m, d->n_ulevels, n_mus2, n_phis, d->n_stokes, 0.);

          if (d->options & XRTM_OPTION_CALC_DERIVS) {
               init_array5_d(I_p_l, d->n_ulevels, d->n_derivs, n_mus2, n_phis, d->n_stokes, 0.);
               init_array5_d(I_m_l, d->n_ulevels, d->n_derivs, n_mus2, n_phis, d->n_stokes, 0.);
          }
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (solutions & XRTM_OUTPUT_RADIANCE) {
          if (apply_corrections_radiance(d, n_phis, phis, I_p, I_m, I_p_l, I_m_l, save_tree, work)) {
               fprintf(stderr, "ERROR: apply_corrections_radiance()\n");
               return -1;
          }
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     count = 0;

     for (i = 0; i < n_four; ++i) {
          if (d->options & XRTM_OPTION_REVERSE_DERIVS)
               save_tree_recode_i(&save_tree, i, i == 0);


          if (get_fourier(d, i, solver, i_p, i_m, i_p_l, i_m_l, save_tree, work)) {
               fprintf(stderr, "ERROR: get_fourier()\n");
               return -1;
          }


          /*--------------------------------------------------------------------
           *
           *------------------------------------------------------------------*/
          if (d->options & XRTM_OPTION_REVERSE_DERIVS) {
               dmat_copy(save->i_p[i], i_p, d->n_ulevels, n_quad_v_x);
               dmat_copy(save->i_m[i], i_m, d->n_ulevels, n_quad_v_x);
          }


          /*--------------------------------------------------------------------
           *
           *------------------------------------------------------------------*/
          if (solutions & XRTM_OUTPUT_RADIANCE_MEAN   && i == 0) {
               for (j = 0; j < d->n_ulevels; ++j) {
                    jj = j == 0 ? 0 : d->n_layers;

                    if (d->options & XRTM_OPTION_UPWELLING_OUTPUT)
                         mean_p[j] = radiance_to_mean(d->n_quad_v, d->qx_v, d->qw_v, d->F_0, d->mu_0, d->btran[jj], i_p[j]);
                    if (d->options & XRTM_OPTION_DOWNWELLING_OUTPUT)
                         mean_m[j] = radiance_to_mean(d->n_quad_v, d->qx_v, d->qw_v, d->F_0, d->mu_0, d->btran[jj], i_m[j]) + d->F_0 * d->btran[jj] / (4. * PI);
               }
          }


          /*--------------------------------------------------------------------
           *
           *------------------------------------------------------------------*/
          if (solutions & XRTM_OUTPUT_FLUX            && i == 0) {
               for (j = 0; j < d->n_ulevels; ++j) {
                    jj = j == 0 ? 0 : d->n_layers;

                    if (d->options & XRTM_OPTION_UPWELLING_OUTPUT)
                         flux_p[j] = radiance_to_flux(d->n_quad_v, d->qx_v, d->qw_v, d->F_0, d->mu_0, d->btran[jj], i_p[j]);
                    if (d->options & XRTM_OPTION_DOWNWELLING_OUTPUT)
                         flux_m[j] = radiance_to_flux(d->n_quad_v, d->qx_v, d->qw_v, d->F_0, d->mu_0, d->btran[jj], i_m[j]) + d->F_0 * d->mu_0 * d->btran[jj];
               }
          }


          /*--------------------------------------------------------------------
           *
           *------------------------------------------------------------------*/
          if (solutions & XRTM_OUTPUT_FLUX_DIVERGENCE && i == 0) {
               for (j = 0; j < d->n_ulevels; ++j) {
                    jj = j == 0 ? 0 : d->n_layers;

               }
          }


          /*--------------------------------------------------------------------
           *
           *------------------------------------------------------------------*/
          if (apply_corrections_fourier(d, i, i_p, i_m, i_p_l, i_m_l, save_tree, work)) {
               fprintf(stderr, "ERROR: apply_corrections_fourier()\n");
               return -1;
          }


          /*--------------------------------------------------------------------
           *
           *------------------------------------------------------------------*/
          if (solutions & XRTM_OUTPUT_RADIANCE) {
               e = 0.;
               for (j = 0; j < d->n_ulevels; ++j) {
                    kk = i_mus2 * d->n_stokes;
                    for (k = 0; k < n_mus2; ++k) {
                         for (l = 0; l < n_phis; ++l) {
                              a = i * (phis[k][l] - d->phi_0)*D2R;

                              b = cos(a);
                              if (d->n_stokes > 2)
                                   c = sin(a);

                              if (d->options & XRTM_OPTION_UPWELLING_OUTPUT) {
                                   azimuth_modes(1, d->n_stokes, b, c, &i_p[j][kk], v1);
                                   dvec_add(I_p[j][k][l], v1, I_p[j][k][l], d->n_stokes);
                              }

                              if (d->options & XRTM_OPTION_DOWNWELLING_OUTPUT) {
                                   azimuth_modes(1, d->n_stokes, b, c, &i_m[j][kk], v2);
                                   dvec_add(I_m[j][k][l], v2, I_m[j][k][l], d->n_stokes);
                              }

                              if (d->fourier_tol > 0. && i > 0) {
                                   if (d->options & XRTM_OPTION_UPWELLING_OUTPUT)
                                        e = MAX(e, max_radiance_ratio(I_p[j][k][l][0], v1[0]));

                                   if (d->options & XRTM_OPTION_DOWNWELLING_OUTPUT)
                                        e = MAX(e, max_radiance_ratio(I_m[j][k][l][0], v2[0]));
                              }

                              if (d->options & XRTM_OPTION_CALC_DERIVS) {
                                   for (m = 0; m < d->n_derivs; ++m) {
                                        if (d->options & XRTM_OPTION_UPWELLING_OUTPUT) {
                                             azimuth_modes(1, d->n_stokes, b, c, &i_p_l[j][m][kk], v1);
                                             dvec_add(I_p_l[j][m][k][l], v1, I_p_l[j][m][k][l], d->n_stokes);
                                        }

                                        if (d->options & XRTM_OPTION_DOWNWELLING_OUTPUT) {
                                             azimuth_modes(1, d->n_stokes, b, c, &i_m_l[j][m][kk], v2);
                                             dvec_add(I_m_l[j][m][k][l], v2, I_m_l[j][m][k][l], d->n_stokes);
                                        }
                                   }
                              }
                         }

                         kk += d->n_stokes;
                    }
               }

               if (i > 0 && e < d->fourier_tol)
                    count++;
               else
                    count = 0;

               if (count >= 2) {
                    i++; break;
               }
          }
     }

     if (d->options & XRTM_OPTION_REVERSE_DERIVS)
          save_tree_decode_i(&save_tree);

     d->misc_output.fourier_count = i;

     if (d->options & XRTM_OPTION_REVERSE_DERIVS)
          save->count = i;

#ifdef USE_AD_FOR_TL_CALC_GET_SOLUTION_INTERNAL
     get_solution_internal_tl_with_ad(d, solver, solutions, n_phis, phis, I_p, I_m, I_p_l, I_m_l, mean_p, mean_m, mean_p_l, mean_m_l, flux_p, flux_m, flux_p_l, flux_m_l, flux_div, flux_div_l, save_tree, work);
#endif

     return 0;
}



/*******************************************************************************
 *
 ******************************************************************************/
int get_solution(xrtm_data *d, enum xrtm_solver_mask solver, int solutions, int n_phis, double **phis, double ****I_p, double ****I_m, double *****I_p_l, double *****I_m_l, double *mean_p, double *mean_m, double **mean_p_l, double **mean_m_l, double *flux_p, double *flux_m, double **flux_p_l, double **flux_m_l, double *flux_div, double **flux_div_l, save_tree_data save_tree, work_data work) {

     if (solver & XRTM_SOLVERS_INTERNAL) {
          if (get_solution_internal(d, solver, solutions, n_phis, phis, I_p, I_m, I_p_l, I_m_l, mean_p, mean_m, mean_p_l, mean_m_l, flux_p, flux_m, flux_p_l, flux_m_l, flux_div, flux_div_l, save_tree, work)) {
               fprintf(stderr, "ERROR: get_radiance_internal()\n");
               return -1;
          }
     }
#ifdef INCLUDE_DEV_SOURCE
     else {
          if (get_solution_external(d, solver, solutions, n_phis, phis, I_p, I_m, I_p_l, I_m_l, mean_p, mean_m, mean_p_l, mean_m_l, flux_p, flux_m, flux_p_l, flux_m_l, flux_div, flux_div_l, work)) {
               fprintf(stderr, "ERROR: get_radiance_external()\n");
               return -1;
          }
     }
#endif
     return 0;
}



/*******************************************************************************
 *
 ******************************************************************************/
#ifdef INCLUDE_DEV_SOURCE
#include "../src2/xrtm_model2.c"
#endif
