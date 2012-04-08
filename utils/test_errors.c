/******************************************************************************%
**
**    Copyright (C) 2007-2012 Greg McGarragh <gregm@atmos.colostate.edu>
**
**    This source code is licensed under the GNU General Public License (GPL),
**    Version 3.  See the file COPYING for more details.
**
*******************************************************************************/

#include <gutil.h>

#include <xrtm_interface.h>
#include <xrtm_support.h>

#include "test.h"
#include "test_macros.h"


/*******************************************************************************
 *
 ******************************************************************************/
#define FAIL_LOCATION() do {								\
     fprintf (stderr, "FAIL: see error test in \"%s\" at line %d\n", __FILE__, __LINE__);\
     return -1;										\
} while (0)


#define TEST_ERROR(x) do {								\
     if (! x)										\
          FAIL_LOCATION();								\
} while (0)

#define TEST_ERROR2(x, y) do {								\
     if (x > y)										\
          FAIL_LOCATION();								\
} while (0)



/*******************************************************************************
 *
 ******************************************************************************/
int test_errors(test_data *t) {

     int options;
     int solvers;

     int n_coef_layer;

     enum xrtm_kernel_type kernels;

     int s_pade;
     int r_pade;

     double thetas;

     double g;
     double omega;
     double ltau;

     xrtm_data d;


     /**************************************************************************
      *
      *************************************************************************/
/*
     XRTM_CREATE(&d, options, solvers, 32, 16, 1, 0, 1, 0, 32, NULL, 0, 1);
*/

     /**************************************************************************
      *
      *************************************************************************/
     options = 0;
     solvers = XRTM_SOLVER_DOUB_ADD;

     /* d->n_coef < 0 */
     TEST_ERROR(xrtm_create(&d, options, solvers, -1, 16, 1, 0, 1, 0, 32, NULL, 2, 1));

     /* d->n_quad <= 0 */
     TEST_ERROR(xrtm_create(&d, options, solvers, 32, 0, 1, 0, 1, 0, 32, NULL, 2, 1));

     /* d->n_stokes <= 0 || d->n_stokes > 4 */
     TEST_ERROR(xrtm_create(&d, options, solvers, 32, 16, 0, 0, 1, 0, 32, NULL, 2, 1));
     TEST_ERROR(xrtm_create(&d, options, solvers, 32, 16, 5, 0, 1, 0, 32, NULL, 2, 1));

     /* d->n_derivs < 0 */
     TEST_ERROR(xrtm_create(&d, options, solvers, 32, 16, 1, -1, 1, 0, 32, NULL, 2, 1));

     /* d->n_layers <= 0 */
     TEST_ERROR(xrtm_create(&d, options, solvers, 32, 16, 1, 0, -1, 0, 32, NULL, 2, 1));

     /* d->n_kernels < 0 */
     TEST_ERROR(xrtm_create(&d, options, solvers, 32, 16, 1, 0, 1, -1, 32, NULL, 2, 1));

     /* d->n_kernels > 0 && d->n_kernel_quad <= 0 */
     TEST_ERROR(xrtm_create(&d, options, solvers, 32, 16, 1, 0, 1, 1, 0, NULL, 2, 1));

     /* d->n_out_levels <= 0 */
     TEST_ERROR(xrtm_create(&d, options, solvers, 32, 16, 1, 0, 1, 0, 32, NULL, 0, 1));

     /* d->n_mus < 0 */
     TEST_ERROR(xrtm_create(&d, options, solvers, 32, 16, 1, 0, 1, 0, 32, NULL, 2, -1));


     /**************************************************************************
      *
      *************************************************************************/
     solvers = XRTM_SOLVER_DOUB_ADD;

     /* ! XRTM_OPTION_CALC_DERIVS && d->n_derivs > 0 */
     options = 0;
     TEST_ERROR(xrtm_create(&d, options, solvers, 32, 16, 1, 1, 1, 0, 32, NULL, 2, 1));

     /* XRTM_OPTION_CALC_DERIVS && d->n_derivs == 0 */
     options = XRTM_OPTION_CALC_DERIVS;
     TEST_ERROR(xrtm_create(&d, options, solvers, 32, 16, 1, 0, 1, 0, 32, NULL, 2, 1));


     /* ! XRTM_OPTION_N_T_TMS && XRTM_OPTION_FOUR_CONV_NEW */
/*
     options = XRTM_OPTION_FOUR_CONV_NEW;
     TEST_ERROR(xrtm_create(&d, options, solvers, 32, 16, 1, 0, 1, 0, 32, NULL, 2, 1));
*/
     /* XRTM_OPTION_N_T_TMS && XRTM_OPTION_NO_AZIMUTHAL */
     options = XRTM_OPTION_N_T_TMS | XRTM_OPTION_NO_AZIMUTHAL;
     TEST_ERROR(xrtm_create(&d, options, solvers, 33, 16, 1, 0, 1, 0, 32, NULL, 2, 1));


     /* XRTM_OPTION_FOUR_CONV_OLD && XRTM_OPTION_FOUR_CONV_NEW */
     options = XRTM_OPTION_FOUR_CONV_OLD | XRTM_OPTION_FOUR_CONV_NEW;
     TEST_ERROR(xrtm_create(&d, options, solvers, 32, 16, 1, 0, 1, 0, 32, NULL, 2, 1));


     /* XRTM_OPTION_OUTPUT_AT_LEVELS && XRTM_OPTION_OUTPUT_AT_TAUS */
     options = XRTM_OPTION_OUTPUT_AT_LEVELS | XRTM_OPTION_OUTPUT_AT_TAUS;
     TEST_ERROR(xrtm_create(&d, options, solvers, 32, 16, 1, 0, 1, 0, 32, NULL, 2, 1));


     /* XRTM_OPTION_PHASE_SCALAR && XRTM_OPTION_PHASE_MATRIX_GC */
     options = XRTM_OPTION_PHASE_SCALAR    | XRTM_OPTION_PHASE_MATRIX_GC;
     TEST_ERROR(xrtm_create(&d, options, solvers, 32, 16, 1, 0, 1, 0, 32, NULL, 2, 1));

     /* XRTM_OPTION_PHASE_SCALAR && XRTM_OPTION_PHASE_MATRIX_LC */
     options = XRTM_OPTION_PHASE_SCALAR    | XRTM_OPTION_PHASE_MATRIX_LC;
     TEST_ERROR(xrtm_create(&d, options, solvers, 32, 16, 1, 0, 1, 0, 32, NULL, 2, 1));

     /* XRTM_OPTION_PHASE_MATRIX_GC && XRTM_OPTION_PHASE_MATRIX_LC */
     options = XRTM_OPTION_PHASE_MATRIX_GC | XRTM_OPTION_PHASE_MATRIX_LC;
     TEST_ERROR(xrtm_create(&d, options, solvers, 32, 16, 1, 0, 1, 0, 32, NULL, 2, 1));


     /* XRTM_OPTION_QUAD_NORM_GAUS_LEG && XRTM_OPTION_QUAD_DOUB_GAUS_LEG */
     options = XRTM_OPTION_QUAD_NORM_GAUS_LEG | XRTM_OPTION_QUAD_DOUB_GAUS_LEG;
     TEST_ERROR(xrtm_create(&d, options, solvers, 32, 16, 1, 0, 1, 0, 32, NULL, 2, 1));

     /* XRTM_OPTION_QUAD_NORM_GAUS_LEG && XRTM_OPTION_QUAD_LOBATTO */
     options = XRTM_OPTION_QUAD_NORM_GAUS_LEG | XRTM_OPTION_QUAD_LOBATTO;
     TEST_ERROR(xrtm_create(&d, options, solvers, 32, 16, 1, 0, 1, 0, 32, NULL, 2, 1));

     /* XRTM_OPTION_QUAD_DOUB_GAUS_LEG && XRTM_OPTION_QUAD_LOBATTO */
     options = XRTM_OPTION_QUAD_DOUB_GAUS_LEG | XRTM_OPTION_QUAD_LOBATTO;
     TEST_ERROR(xrtm_create(&d, options, solvers, 32, 16, 1, 0, 1, 0, 32, NULL, 2, 1));


     /* XRTM_OPTION_SAVE_LAYER_R_T_S && XRTM_OPTION_STACK_REUSE_ADDING */
     options = XRTM_OPTION_SAVE_LAYER_R_T_S | XRTM_OPTION_STACK_REUSE_ADDING;
     TEST_ERROR(xrtm_create(&d, options, solvers, 32, 16, 1, 0, 1, 0, 32, NULL, 2, 1));

     /* XRTM_OPTION_SAVE_LAYER_R_T_S && ! XRTM_OPTION_CALC_DERIVS */
     options = XRTM_OPTION_STACK_REUSE_ADDING;
     TEST_ERROR(xrtm_create(&d, options, solvers, 32, 16, 1, 0, 1, 0, 32, NULL, 2, 1));


     /* XRTM_OPTION_TOP_DOWN_ADDING && n_out_levels != 1 */
     options = XRTM_OPTION_TOP_DOWN_ADDING;
     TEST_ERROR(xrtm_create(&d, options, solvers, 32, 16, 1, 0, 1, 0, 32, NULL, 2, 1));

     /* XRTM_OPTION_BOTTOM_UP_ADDING && n_out_levels != 1 */
     options = XRTM_OPTION_BOTTOM_UP_ADDING;
     TEST_ERROR(xrtm_create(&d, options, solvers, 32, 16, 1, 0, 1, 0, 32, NULL, 2, 1));


     /* ! XRTM_OPTION_VECTOR | n_stokes > 1 */
/*
     options = 0;
     TEST_ERROR(xrtm_create(&d, options, solvers, 32, 16, 4, 0, 1, 0, 32, NULL, 2, 1));
*/
     /* ! XRTM_OPTION_VECTOR && (XRTM_OPTION_PHASE_MATRIX_GC | XRTM_OPTION_PHASE_MATRIX_LC) */
     options = XRTM_OPTION_PHASE_MATRIX_GC;
     options = XRTM_OPTION_PHASE_MATRIX_LC;
     TEST_ERROR(xrtm_create(&d, options, solvers, 32, 16, 1, 0, 1, 0, 32, NULL, 2, 1));

     /* XRTM_OPTION_VECTOR & ! (XRTM_OPTION_PHASE_MATRIX_GC | XRTM_OPTION_PHASE_MATRIX_LC) */
/*
     options = XRTM_OPTION_VECTOR;
     TEST_ERROR(xrtm_create(&d, options, solvers, 32, 16, 4, 0, 1, 0, 32, NULL, 2, 1));
*/

     /**************************************************************************
      *
      *************************************************************************/
#ifdef INCLUDE_DEV_SOURCE
     if (test_errors_dev_solvers(t)) {
          eprintf("ERROR: test_errors_dev_solvers()\n");
          return -1;
     }
#endif

     /**************************************************************************
      *
      *************************************************************************/
     /* (d->options & XRTM_OPTION_DELTA_M) && d->n_coef < 2 * d->n_quad + 1 */
     options = XRTM_OPTION_DELTA_M;
     solvers = 0;
     TEST_ERROR(xrtm_create(&d, options, solvers, 2 * 16 + 1 - 1, 16, 1, 0, 1, 0, 32, NULL, 2, 1));


     /**************************************************************************
      *
      *************************************************************************/
     options = XRTM_OPTION_SOURCE_THERMAL;
     solvers = XRTM_SOLVER_DOUB_ADD | XRTM_SOLVER_PADE_ADD;
     kernels = XRTM_KERNEL_LI_SPARSE;
     XRTM_CREATE(&d, options, solvers, 32, 16, 1, 0, 1, 1, 32, &kernels, 2, 1);


     TEST_ERROR2(xrtm_get_doub_d_tau(&d), -1.);
     TEST_ERROR(xrtm_set_doub_d_tau(&d, -1.e9));
     TEST_ERROR(xrtm_set_doub_d_tau(&d, 0.));


     TEST_ERROR2(xrtm_get_pade_params(&d, &s_pade, &r_pade), -1.);
/*
     TEST_ERROR(xrtm_set_pade_params(&d, -1, 1));
     TEST_ERROR(xrtm_set_pade_params(&d, 1, 0));
*/

     TEST_ERROR2(xrtm_get_lambda(&d), -1.);
     TEST_ERROR(xrtm_set_lambda(&d, 0.));


     TEST_ERROR2(xrtm_get_F_0(&d), -1.);
     TEST_ERROR(xrtm_set_F_0(&d, -1.));


     TEST_ERROR2(xrtm_get_theta_0(&d), -1.);
     TEST_ERROR(xrtm_set_theta_0(&d, -1.));
     TEST_ERROR(xrtm_set_theta_0(&d, 90.));

     TEST_ERROR2(xrtm_get_phi_0(&d), -1.);
     TEST_ERROR(xrtm_set_phi_0(&d, -1.));


     xrtm_destroy(&d);


     /**************************************************************************
      *
      *************************************************************************/
     options = XRTM_OPTION_OUTPUT_AT_TAUS;
     solvers = XRTM_SOLVER_EIG_BVP;
     XRTM_CREATE(&d, options, solvers, 32, 16, 1, 0, 1, 0, 32, NULL, 2, 1);

     TEST_ERROR(xrtm_set_out_levels(&d, NULL));

     xrtm_destroy(&d);


     options = XRTM_OPTION_OUTPUT_AT_LEVELS;
     XRTM_CREATE(&d, options, solvers, 32, 16, 1, 0, 1, 0, 32, NULL, 2, 1);

     TEST_ERROR(xrtm_get_out_levels(&d, NULL));

     {
          int out_levels[] = {-1, 1};
          TEST_ERROR(xrtm_set_out_levels(&d, out_levels));
     } {
          int out_levels[] = {0, 2};
          TEST_ERROR(xrtm_set_out_levels(&d, out_levels));
     }
     xrtm_destroy(&d);


     /**************************************************************************
      *
      *************************************************************************/
     options = XRTM_OPTION_OUTPUT_AT_LEVELS;
     solvers = XRTM_SOLVER_EIG_BVP;
     XRTM_CREATE(&d, options, solvers, 32, 16, 1, 0, 1, 0, 32, NULL, 2, 1);

     TEST_ERROR(xrtm_set_out_taus(&d, NULL));

     xrtm_destroy(&d);


     options = XRTM_OPTION_OUTPUT_AT_TAUS;
     XRTM_CREATE(&d, options, solvers, 32, 16, 1, 0, 1, 0, 32, NULL, 2, 1);

     TEST_ERROR(xrtm_get_out_taus(&d, NULL));
/*
     {
          double out_taus[] = {-1., 1.};
          TEST_ERROR(xrtm_set_out_taus(&d, out_taus));
     } {
          double out_taus[] = {0., 2.};
          TEST_ERROR(xrtm_set_out_taus(&d, out_taus));
     }
*/
     xrtm_destroy(&d);


     /**************************************************************************
      *
      *************************************************************************/
     options = 0;
     solvers = XRTM_SOLVER_DOUB_ADD;
     XRTM_CREATE(&d, options, solvers, 32, 16, 1, 0, 1, 0, 32, NULL, 2, 1);

     TEST_ERROR(xrtm_get_out_thetas(&d, NULL));

     thetas = -1.;
     TEST_ERROR(xrtm_set_out_thetas(&d, &thetas));

     thetas = 90.;
     TEST_ERROR(xrtm_set_out_thetas(&d, &thetas));

     xrtm_destroy(&d);


     /**************************************************************************
      *
      *************************************************************************/
     options = 0;
     solvers = XRTM_SOLVER_DOUB_ADD;
     XRTM_CREATE(&d, options, solvers, 32, 16, 1, 0, 1, 0, 32, NULL, 2, 1);

     TEST_ERROR(xrtm_set_planet_r(&d, 0.));
     TEST_ERROR(xrtm_set_levels_z(&d, NULL));

     xrtm_destroy(&d);


     options = XRTM_OPTION_PSA;
     XRTM_CREATE(&d, options, solvers, 32, 16, 1, 0, 1, 0, 32, NULL, 2, 1);

     TEST_ERROR2(xrtm_get_planet_r(&d), -1.);
     TEST_ERROR(xrtm_set_planet_r(&d, 0.));

     TEST_ERROR(xrtm_get_levels_z(&d, NULL));
     {
          double levels_z[] = {0., 0.};
          TEST_ERROR(xrtm_set_levels_z(&d, levels_z));
     } {
          double levels_z[] = {0., 1.};
          TEST_ERROR(xrtm_set_levels_z(&d, levels_z));
     }
     xrtm_destroy(&d);


     /**************************************************************************
      *
      *************************************************************************/
     options = 0;
     solvers = XRTM_SOLVER_DOUB_ADD;
     XRTM_CREATE(&d, options, solvers, 32, 16, 1, 0, 1, 0, 32, NULL, 2, 1);

     TEST_ERROR(xrtm_set_levels_b(&d, NULL));

     xrtm_destroy(&d);


     options = XRTM_OPTION_SOURCE_THERMAL;
     XRTM_CREATE(&d, options, solvers, 32, 16, 1, 0, 1, 0, 32, NULL, 2, 1);

     TEST_ERROR(xrtm_get_levels_b(&d, NULL));
     {
          double levels_b[] = {273., -1.};
          TEST_ERROR(xrtm_set_levels_b(&d, levels_b));
     }
     xrtm_destroy(&d);


     /**************************************************************************
      *
      *************************************************************************/
     options = XRTM_OPTION_QUAD_NORM_GAUS_LEG;
     solvers = XRTM_SOLVER_DOUB_ADD;
     XRTM_CREATE(&d, options, solvers, 2, 1, 1, 0, 1, 0, 32, NULL, 2, 1);

     TEST_ERROR(xrtm_set_g_1(&d, 0, 0.));
     TEST_ERROR(xrtm_set_g_n(&d, &g));

     xrtm_destroy(&d);


     solvers = XRTM_SOLVERS_USE_G;
     XRTM_CREATE(&d, options, solvers, 2, 1, 1, 0, 1, 0, 32, NULL, 2, 1);

     TEST_ERROR2(xrtm_get_g(&d, 0), -1.);

     TEST_ERROR(xrtm_set_g_1(&d, -1, 0.));
     TEST_ERROR(xrtm_set_g_1(&d, 1, 0.));
     TEST_ERROR(xrtm_set_g_1(&d, 0, -1.1));
     TEST_ERROR(xrtm_set_g_1(&d, 0, 1.1));

     g = -1.1;
     TEST_ERROR(xrtm_set_g_n(&d, &g));

     g = 1.1;
     TEST_ERROR(xrtm_set_g_n(&d, &g));

     TEST_ERROR(xrtm_set_g_l_11(&d, 0, 0, 0.));
     TEST_ERROR(xrtm_set_g_l_n1(&d, 0, NULL));
     TEST_ERROR(xrtm_set_g_l_1n(&d, 0, NULL));
     TEST_ERROR(xrtm_set_g_l_nn(&d, NULL));

     xrtm_destroy(&d);

/*
     options = XRTM_OPTION_CALC_DERIVS;
     solvers = XRTM_SOLVERS_USE_G;
     XRTM_CREATE(&d, options, solvers, 2, 1, 1, 1, 1, 0, 32, NULL, 2, 1);

     TEST_ERROR(xrtm_set_g_l_11(&d, -1, 0, 0.));
     TEST_ERROR(xrtm_set_g_l_11(&d, 1, 0, 0.));
     TEST_ERROR(xrtm_set_g_l_11(&d, 0, -1, 0.));
     TEST_ERROR(xrtm_set_g_l_11(&d, 0, 1, 0.));
     TEST_ERROR(xrtm_set_g_l_n1(&d, -1, NULL));
     TEST_ERROR(xrtm_set_g_l_n1(&d, 1, NULL));
     TEST_ERROR(xrtm_set_g_l_1n(&d, -1, NULL));
     TEST_ERROR(xrtm_set_g_l_1n(&d, 1, NULL));

     xrtm_destroy(&d);
*/

     /**************************************************************************
      *
      *************************************************************************/
     options = XRTM_OPTION_QUAD_NORM_GAUS_LEG;
     solvers = XRTM_SOLVERS_USE_G;
     XRTM_CREATE(&d, options, solvers, 2, 1, 1, 0, 1, 0, 32, NULL, 2, 1);

     TEST_ERROR(xrtm_set_coef_1(&d, 0, 0, NULL));
     TEST_ERROR(xrtm_set_coef_n(&d, 0, NULL));

     xrtm_destroy(&d);


     options = 0;
     solvers = XRTM_SOLVER_DOUB_ADD;
     XRTM_CREATE(&d, options, solvers, 4, 2, 1, 0, 1, 0, 32, NULL, 2, 1);

     TEST_ERROR(xrtm_get_n_coef(&d, 0));

     TEST_ERROR(xrtm_set_coef_1(&d, -1, 0, NULL));
     TEST_ERROR(xrtm_set_coef_1(&d, 1, 0, NULL));

     TEST_ERROR(xrtm_set_coef_1(&d, 0, 5, NULL));
     {
          double chi[] = {1.1, 0., 0., 0.};
          double *chi2  = chi;
          TEST_ERROR(xrtm_set_coef_1(&d, 0, 4, (double **) &chi2));
     }

     n_coef_layer = -1;
     TEST_ERROR(xrtm_set_coef_n(&d, &n_coef_layer, NULL));
     n_coef_layer = 5;
     TEST_ERROR(xrtm_set_coef_n(&d, &n_coef_layer, NULL));
     {
          n_coef_layer = 4;
          double chi[] = {1.1, 0., 0., 0.};
          double *chi2 = chi;
          double **chi3 = &chi2;
          TEST_ERROR(xrtm_set_coef_n(&d, &n_coef_layer, &chi3));
     }

     TEST_ERROR(xrtm_set_coef_l_11(&d, 0, 0, NULL));
     TEST_ERROR(xrtm_set_coef_l_n1(&d, 0, NULL));
     TEST_ERROR(xrtm_set_coef_l_1n(&d, 0, NULL));
     TEST_ERROR(xrtm_set_coef_l_nn(&d, NULL));

     xrtm_destroy(&d);


     options = XRTM_OPTION_CALC_DERIVS;
     solvers = XRTM_SOLVER_DOUB_ADD;
     XRTM_CREATE(&d, options, solvers, 4, 2, 1, 1, 1, 0, 32, NULL, 2, 1);

     TEST_ERROR(xrtm_set_coef_l_11(&d, -1, 0, NULL));
     TEST_ERROR(xrtm_set_coef_l_11(&d, 1, 0, NULL));
     TEST_ERROR(xrtm_set_coef_l_11(&d, 0, -1, NULL));
     TEST_ERROR(xrtm_set_coef_l_11(&d, 0, 1, NULL));
     TEST_ERROR(xrtm_set_coef_l_n1(&d, -1, NULL));
     TEST_ERROR(xrtm_set_coef_l_n1(&d, 1, NULL));
     TEST_ERROR(xrtm_set_coef_l_1n(&d, -1, NULL));
     TEST_ERROR(xrtm_set_coef_l_1n(&d, 1, NULL));

     xrtm_destroy(&d);


     /**************************************************************************
      *
      *************************************************************************/
     options = 0;
     solvers = XRTM_SOLVER_DOUB_ADD;
     XRTM_CREATE(&d, options, solvers, 32, 16, 1, 0, 1, 0, 32, NULL, 2, 1);

     TEST_ERROR2(xrtm_get_omega(&d, 0), -1.);

     TEST_ERROR(xrtm_set_omega_1(&d, -1, 1.));
     TEST_ERROR(xrtm_set_omega_1(&d, 1, 1.));
     TEST_ERROR(xrtm_set_omega_1(&d, 0, -0.1));
     TEST_ERROR(xrtm_set_omega_1(&d, 0, 1.1));

     omega = -0.1;
     TEST_ERROR(xrtm_set_omega_n(&d, &omega));

     omega = 1.1;
     TEST_ERROR(xrtm_set_omega_n(&d, &omega));

     TEST_ERROR(xrtm_set_omega_l_11(&d, 0, 0, 0.));
     TEST_ERROR(xrtm_set_omega_l_n1(&d, 0, NULL));
     TEST_ERROR(xrtm_set_omega_l_1n(&d, 0, NULL));
     TEST_ERROR(xrtm_set_omega_l_nn(&d, NULL));

     xrtm_destroy(&d);


     options = XRTM_OPTION_CALC_DERIVS;
     solvers = XRTM_SOLVER_DOUB_ADD;
     XRTM_CREATE(&d, options, solvers, 32, 16, 1, 1, 1, 0, 32, NULL, 2, 1);

     TEST_ERROR(xrtm_set_omega_l_11(&d, -1, 0, 0.));
     TEST_ERROR(xrtm_set_omega_l_11(&d, 1, 0, 0.));
     TEST_ERROR(xrtm_set_omega_l_11(&d, 0, -1, 0.));
     TEST_ERROR(xrtm_set_omega_l_11(&d, 0, 1, 0.));
     TEST_ERROR(xrtm_set_omega_l_n1(&d, -1, NULL));
     TEST_ERROR(xrtm_set_omega_l_n1(&d, 1, NULL));
     TEST_ERROR(xrtm_set_omega_l_1n(&d, -1, NULL));
     TEST_ERROR(xrtm_set_omega_l_1n(&d, 1, NULL));

     xrtm_destroy(&d);


     /**************************************************************************
      *
      *************************************************************************/
     options = 0;
     solvers = XRTM_SOLVER_DOUB_ADD;
     XRTM_CREATE(&d, options, solvers, 32, 16, 1, 0, 1, 0, 32, NULL, 2, 1);

     TEST_ERROR2(xrtm_get_ltau(&d, 0), -1.);

     TEST_ERROR(xrtm_set_ltau_1(&d, -1, 1.));
     TEST_ERROR(xrtm_set_ltau_1(&d, 1, 1.));
     TEST_ERROR(xrtm_set_ltau_1(&d, 0, -1.));

     ltau = -1.;
     TEST_ERROR(xrtm_set_ltau_n(&d, &ltau));

     TEST_ERROR(xrtm_set_ltau_l_11(&d, 0, 0, 0.));
     TEST_ERROR(xrtm_set_ltau_l_n1(&d, 0, NULL));
     TEST_ERROR(xrtm_set_ltau_l_1n(&d, 0, NULL));
     TEST_ERROR(xrtm_set_ltau_l_nn(&d, NULL));

     xrtm_destroy(&d);


     options = XRTM_OPTION_CALC_DERIVS;
     solvers = XRTM_SOLVER_DOUB_ADD;
     XRTM_CREATE(&d, options, solvers, 32, 16, 1, 1, 1, 0, 32, NULL, 2, 1);

     TEST_ERROR(xrtm_set_ltau_l_11(&d, -1, 0, 0.));
     TEST_ERROR(xrtm_set_ltau_l_11(&d, 1, 0, 0.));
     TEST_ERROR(xrtm_set_ltau_l_11(&d, 0, -1, 0.));
     TEST_ERROR(xrtm_set_ltau_l_11(&d, 0, 1, 0.));
     TEST_ERROR(xrtm_set_ltau_l_n1(&d, -1, NULL));
     TEST_ERROR(xrtm_set_ltau_l_n1(&d, 1, NULL));
     TEST_ERROR(xrtm_set_ltau_l_1n(&d, -1, NULL));
     TEST_ERROR(xrtm_set_ltau_l_1n(&d, 1, NULL));

     xrtm_destroy(&d);


     /**************************************************************************
      *
      *************************************************************************/
     options = 0;
     solvers = XRTM_SOLVER_DOUB_ADD;
     XRTM_CREATE(&d, options, solvers, 32, 16, 1, 0, 1, 0, 32, NULL, 2, 1);

     TEST_ERROR(xrtm_set_surface_b(&d, 0.));

     xrtm_destroy(&d);


     options = XRTM_OPTION_SOURCE_THERMAL;
     XRTM_CREATE(&d, options, solvers, 32, 16, 1, 0, 1, 0, 32, NULL, 2, 1);

     TEST_ERROR2(xrtm_get_surface_b(&d), -1.);

     TEST_ERROR(xrtm_set_surface_b(&d, -1.));

     xrtm_destroy(&d);


     /**************************************************************************
      *
      *************************************************************************/
     options = 0;
     solvers = XRTM_SOLVER_DOUB_ADD;
     XRTM_CREATE(&d, options, solvers, 32, 16, 1, 0, 1, 0, 32, NULL, 2, 1);

     TEST_ERROR(xrtm_set_kernel_ampfac(&d, 0, 1.));
     TEST_ERROR(xrtm_set_kernel_params_1(&d, 0, 0, 0.));
     TEST_ERROR(xrtm_set_kernel_params_n(&d, 0, NULL));

     xrtm_destroy(&d);


     kernels = XRTM_KERNEL_LI_SPARSE;
     XRTM_CREATE(&d, options, solvers, 32, 16, 1, 0, 1, 1, 32, &kernels, 2, 1);

     TEST_ERROR2(xrtm_get_kernel_ampfac(&d, 0), -1.);
     TEST_ERROR2(xrtm_get_kernel_params(&d, 0, 0), -1.);

     TEST_ERROR(xrtm_set_kernel_ampfac(&d, -1, 1.));
     TEST_ERROR(xrtm_set_kernel_ampfac(&d, 1, 1.));
     TEST_ERROR(xrtm_set_kernel_ampfac(&d, 0, -0.1));
     TEST_ERROR(xrtm_set_kernel_ampfac(&d, 0, 1.1));

     TEST_ERROR(xrtm_set_kernel_params_1(&d, -1, 0, 0.));
     TEST_ERROR(xrtm_set_kernel_params_1(&d, 1, 0, 0.));

     TEST_ERROR(xrtm_set_kernel_params_n(&d, -1, NULL));
     TEST_ERROR(xrtm_set_kernel_params_n(&d, 1, NULL));

     xrtm_destroy(&d);


     /**************************************************************************
      *
      *************************************************************************/
     options = 0;
     solvers = XRTM_SOLVER_DOUB_ADD;
     XRTM_CREATE(&d, options, solvers, 32, 16, 1, 0, 1, 0, 32, NULL, 2, 1);
/*
     TEST_ERROR(xrtm_local_r_t_u_w(&d, -1, 0, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL));
     TEST_ERROR(xrtm_local_r_t_u_w(&d, 32, 0, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL));
     TEST_ERROR(xrtm_local_r_t_u_w(&d, 0, -1, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL));
     TEST_ERROR(xrtm_local_r_t_u_w(&d, 0, 1, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL));
     TEST_ERROR(xrtm_local_r_t_u_w(&d, 0, 0, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL));
*/
     TEST_ERROR(xrtm_solution(&d, XRTM_SOLVER_EIG_ADD, XRTM_OUTPUT_RADIANCE, 1, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL));
     TEST_ERROR(xrtm_solution(&d, XRTM_SOLVER_DOUB_ADD, XRTM_OUTPUT_RADIANCE, 1, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL));

     TEST_ERROR(xrtm_radiance(&d, XRTM_SOLVER_EIG_ADD, 1, NULL, NULL, NULL, NULL, NULL));
     TEST_ERROR(xrtm_radiance(&d, XRTM_SOLVER_DOUB_ADD, 1, NULL, NULL, NULL, NULL, NULL));

     TEST_ERROR(xrtm_mean_radiance(&d, XRTM_SOLVER_EIG_ADD, NULL, NULL, NULL, NULL));
     TEST_ERROR(xrtm_mean_radiance(&d, XRTM_SOLVER_DOUB_ADD, NULL, NULL, NULL, NULL));

     TEST_ERROR(xrtm_flux(&d, XRTM_SOLVER_EIG_ADD, NULL, NULL, NULL, NULL));
     TEST_ERROR(xrtm_flux(&d, XRTM_SOLVER_DOUB_ADD, NULL, NULL, NULL, NULL));

     TEST_ERROR(xrtm_flux_divergence(&d, XRTM_SOLVER_EIG_ADD, NULL, NULL));
     TEST_ERROR(xrtm_flux_divergence(&d, XRTM_SOLVER_DOUB_ADD, NULL, NULL));

     xrtm_destroy(&d);


     return 0;
}



/*******************************************************************************
 *
 ******************************************************************************/
#ifdef INCLUDE_DEV_SOURCE
#include "../utils2/test_errors2.c"
#endif
