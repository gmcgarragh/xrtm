/*******************************************************************************
**
**    Copyright (C) 2007-2020 Greg McGarragh <greg.mcgarragh@colostate.edu>
**
**    This source code is licensed under the GNU General Public License (GPL),
**    Version 3.  See the file COPYING for more details.
**
*******************************************************************************/

#include <gutil.h>

#include <xrtm_interface.h>
#include <xrtm_support.h>

#include "test.h"
#include "test_core.h"
#include "test_loop.h"
#include "test_macros.h"


/*******************************************************************************
 *
 ******************************************************************************/
static int set_xrtm_defaults(test_xrtm_data *test_xrtm, xrtm_input_data *xrtm_input) {

     if (test_xrtm->solvers & XRTM_SOLVERS_DOUBLING) {
          xrtm_input->set_flag_doub_d_tau = 1;
          xrtm_input->doub_d_tau = DEFAULT_DOUB_D_TAU;
     }

     if (test_xrtm->solvers & XRTM_SOLVER_PADE_ADD) {
          xrtm_input->set_flag_pade_params = 1;
          xrtm_input->pade_s = -1;
          xrtm_input->pade_r = -1;
     }

     if (test_xrtm->solvers & XRTM_SOLVER_SOS) {
          xrtm_input->set_flag_sos_params = 1;
          xrtm_input->sos_max_os  = DEFAULT_SOS_PARAMS_MAX_OS;
          xrtm_input->sos_max_tau = DEFAULT_SOS_PARAMS_MAX_TAU;
          xrtm_input->sos_tol     = DEFAULT_SOS_PARAMS_SOS_TOL;
     }

     xrtm_input->set_flag_fourier_tol = 1;
     xrtm_input->fourier_tol = DEFAULT_FOURIER_TOL;

     if (test_xrtm->options & XRTM_OPTION_PSA) {
          xrtm_input->set_flag_planet_r = 1;
          xrtm_input->planet_r = DEFAULT_PLANET_R;
     }

     xrtm_input->set_flag_F_iso_top = 1;
     xrtm_input->F_iso_top = DEFAULT_F_ISO_TOP;
     xrtm_input->set_flag_F_iso_bot = 1;
     xrtm_input->F_iso_bot = DEFAULT_F_ISO_BOT;

     if (test_xrtm->options & XRTM_OPTION_SOURCE_SOLAR) {
          xrtm_input->set_flag_phi_0 = 1;
          xrtm_input->phi_0 = DEFAULT_PHI_0;
     }

     if (test_xrtm->options & XRTM_OPTION_SOURCE_THERMAL) {
          xrtm_input->set_flag_lambda = 1;
          xrtm_input->lambda = 1.;
     }

     return 0;
}



/*******************************************************************************
 *
 ******************************************************************************/
static int setup_bounds_test(test_data *td,
                             int n_solvers, long *solver_list,
                             int (*set_xrtm_defaults)(test_xrtm_data *test_xrtm, xrtm_input_data *xrtm_input),
                             int n_layers, int n_derivs,
                             int n_out_levels, int *out_levels,
                             int n_out_taus, double *out_taus,
                             int n_out_thetas, double *out_thetas,
                             int n_phis, double *phis,
                             double F_0,
                             int n_levels_z_bounds, double *levels_z,
                             int n_theta_0_bounds, double *theta_0,
                             int n_levels_b_bounds, double *levels_b,
                             int n_surface_b_bounds, double *surface_b,
                             int n_coefs_bounds, int *n_coefs, double ***coefs,
                             int n_omega_bounds, double *omega,
                             int n_ltau_bounds, double *ltau,
                             int n_albedo_bounds, double *albedo,
                             int n_deriv_bounds,
                             double *levels_b_l,
                             double ***coefs_l, double *omega_l, double *ltau_l,
                             double *albedo_l,
                             double *tol, double *tol_l,
                             int n_ignore_index_solver,
                             int *ignore_list_index_solver,
                             int *ignore_mask_index_solver,
                             int n_ignore_deriv,
                             int *ignore_list_deriv,
                             int *ignore_mask_deriv,
                             int (*ignore_mask_solver)(double theta_0, double *theta)) {

     int i;
     int j;
     int k;
     int l;
     int m;
     int n;
     int o;
     int p;
     int q;

     int n_kernels = 1;

     test_xrtm_data *test_xrtm;

     xrtm_input_data *xrtm_input;


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     test_xrtm = list_last_elem(&td->test_xrtm_list);


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     test_xrtm->n_solvers                = n_solvers;
     test_xrtm->solver_list              = dup_array1_l(solver_list, n_solvers);

     test_xrtm->n_phi_bounds             = n_phis;
     test_xrtm->phi                      = dup_array1_d(phis, n_phis);

     test_xrtm->tol                      = dup_array1_d(tol,   n_solvers);
     if (n_deriv_bounds > 0)
          test_xrtm->tol_l               = dup_array1_d(tol_l, n_solvers);

     test_xrtm->n_ignore_index_solver    = n_ignore_index_solver;
     test_xrtm->ignore_list_index_solver = dup_array1_i(ignore_list_index_solver, n_ignore_index_solver);
     test_xrtm->ignore_mask_index_solver = dup_array1_i(ignore_mask_index_solver, n_ignore_index_solver);

     test_xrtm->n_ignore_deriv           = n_ignore_deriv;
     test_xrtm->ignore_list_deriv        = dup_array1_i(ignore_list_deriv,        n_ignore_deriv);
     test_xrtm->ignore_mask_deriv        = dup_array1_i(ignore_mask_deriv,        n_ignore_deriv);

     test_xrtm->ignore_mask_solver       = ignore_mask_solver;


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     for (i = 0; i < n_theta_0_bounds; ++i) {
          for (j = 0; j < n_levels_b_bounds; ++j) {
               for (k = 0; k < n_coefs_bounds; ++k) {
                    for (l = 0; l < n_omega_bounds; ++l) {
                         for (m = 0; m < n_ltau_bounds; ++m) {
                              for (n = 0; n < n_albedo_bounds; ++n) {

                                   xrtm_input = malloc(sizeof(xrtm_input_data));
                                   xrtm_input_init(xrtm_input);
                                   list_append(&test_xrtm->xrtm_input_list, xrtm_input, 0);

                                   if (set_xrtm_defaults(test_xrtm, xrtm_input)) {
                                        fprintf(stderr, "ERROR: set_xrtm_defaults()\n");
                                        return -1;
                                   }

                                   if (test_xrtm->options & XRTM_OPTION_PSA) {
                                        xrtm_input->set_flag_levels_z = 1;
                                        xrtm_input->levels_z = dup_array1_d(&levels_z[0 * (n_layers + 1)], n_layers + 1);
                                   }

                                   if (test_xrtm->options & XRTM_OPTION_OUTPUT_AT_LEVELS) {
                                        xrtm_input->set_flag_out_levels = 1;
                                        xrtm_input->out_levels = dup_array1_i(out_levels, n_out_levels);
                                   }
                                   else {
                                        xrtm_input->set_flag_out_taus   = 1;
                                        xrtm_input->out_taus   = dup_array1_d(out_taus,   n_out_taus);
                                   }

                                   if (n_out_thetas > 0) {
                                        xrtm_input->set_flag_out_thetas = 1;
                                        xrtm_input->out_thetas = dup_array1_d(out_thetas, n_out_thetas);
                                   }

                                   xrtm_input->set_flag_F_0 = 1;
                                   xrtm_input->F_0 = F_0;

                                   xrtm_input->set_flag_theta_0 = 1;
                                   xrtm_input->theta_0 = theta_0[i];

                                   xrtm_input->set_flag_levels_b = 1;
                                   xrtm_input->levels_b = dup_array1_d(&levels_b[j * (n_layers + 1)], n_layers + 1);

                                   xrtm_input->set_flag_surface_b = 1;
                                   xrtm_input->surface_b = surface_b[0];

                                   xrtm_input->set_flag_coef = 1;
                                   xrtm_input->n_coef = dup_array1_i(&n_coefs[k * n_layers], n_layers);
                                   xrtm_input->coef = (double ***) dup_array1(&coefs[k * n_layers], n_layers, sizeof(double **));

                                   xrtm_input->set_flag_omega = 1;
                                   xrtm_input->omega = dup_array1_d(&omega[l * n_layers], n_layers);

                                   xrtm_input->set_flag_ltau = 1;
                                   xrtm_input->ltau = dup_array1_d(&ltau  [m * n_layers], n_layers);

                                   xrtm_input->set_flag_ampfac = 1;
                                   xrtm_input->ampfac = dup_array1_d(&albedo[n * n_kernels], n_kernels);

                                   for (o = 0; o < n_deriv_bounds; ++o) {
                                        xrtm_input->set_flag_levels_b_l = 1;
                                        xrtm_input->levels_b_l = alloc_array2_d(n_layers + 1, n_derivs);
                                        xrtm_input->set_flag_coef_l     = 1;
                                        xrtm_input->coef_l     = (double ****) alloc_array2(n_layers, n_derivs, sizeof(double **));
                                        xrtm_input->set_flag_omega_l    = 1;
                                        xrtm_input->omega_l    = alloc_array2_d(n_layers, n_derivs);
                                        xrtm_input->set_flag_ltau_l     = 1;
                                        xrtm_input->ltau_l     = alloc_array2_d(n_layers, n_derivs);

                                        for (p = 0; p < n_layers + 1; ++p)
                                             copy_array1_d(xrtm_input->levels_b_l[p], &levels_b_l[o * (n_layers + 1) * n_derivs + p * n_derivs], n_derivs);

                                        for (p = 0; p < n_layers; ++p) {
                                             for (q = 0; q < n_derivs; ++q) {
                                                  if (! coefs_l[k * n_deriv_bounds * n_layers * n_derivs + o * n_layers * n_derivs + p * n_derivs + q])
                                                       xrtm_input->coef_l[p][q] = NULL;
                                                  else
                                                       xrtm_input->coef_l[p][q] = coefs_l[k * n_deriv_bounds * n_layers * n_derivs + o * n_layers * n_derivs + p * n_derivs + q];
                                             }

                                             copy_array1_d(xrtm_input->omega_l[p], &omega_l[o * n_layers * n_derivs + p * n_derivs], n_derivs);

                                             copy_array1_d(xrtm_input->ltau_l [p], &ltau_l [o * n_layers * n_derivs + p * n_derivs], n_derivs);
                                        }

                                        xrtm_input->set_flag_ampfac_l = 1;
                                        xrtm_input->ampfac_l = alloc_array2_d(n_kernels, n_derivs);
                                        copy_array1_d(xrtm_input->ampfac_l[0], &albedo_l[o * n_derivs], n_derivs);
                                   }
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
#define N_LAYERS_SMALL_STACK			3
#define N_LAYERS_LARGE_STACK			31

#define N_DERIVS_ONE_LAYER_BOUND_DERIVS		17

#define N_DERIVS_SMALL_STACK_BOUND_DERIVS	9


int test_core(test_data *td) {

     uint i;
     uint ii;
     uint j;

     uint count;

     uint mask;

     int n_solvers;
     int n_solvers2;
     int solvers_mask;
     long *solvers_list;

     int n_derivs;
     int n_layers;
     int n_kernels;
     int n_out_levels;
     int n_out_taus;

     int n_derivs2;

     int n_out_levels2;

     int *out_levels;

     double a;

     double *out_taus;
     double *out_thetas;
     double *out_phis;

     int n_levels_z_bounds;
     int n_theta_0_bounds;
     int n_levels_b_bounds;
     int n_surface_b_bounds;
     int n_coefs_bounds;
     int n_omega_bounds;
     int n_ltau_bounds;
     int n_albedo_bounds;

     int n_deriv_bounds = 1;

     double F_0;

     double *levels_z;
     double *theta_0;
     double *levels_b;
     double *surface_b;
     int *n_coefs;
     double ***coefs;
     double *omega;
     double *ltau;
     double *albedo;

     double *levels_b_l;
     double ***coefs_l;
     double *omega_l;
     double *ltau_l;
     double *albedo_l;

     double *tol;
     double *tol_l;

     enum xrtm_kernel_type kernels[N_XRTM_KERNEL_TYPES];

     test_xrtm_data *test_xrtm;


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     double   simple_one_layer_levels_z [][1 + 1] = {{11.,
                                                      9.}};
     double   simple_one_layer_theta_0  []        = { 35. };
     double   simple_one_layer_levels_0 [][1 + 1] = {{0.,
                                                      0.}};
     double   simple_one_layer_surface_0[]        = { 0. };
     double   simple_one_layer_levels_b [][1 + 1] = {{.1,
                                                      .2}};
     double   simple_one_layer_surface_b[]        = { .3 };
     int      simple_one_layer_n_coefs  [][1]     = {{td->n_gc[2]}};
     double **simple_one_layer_coefs    [][1]     = {{td->  gc[2]}};
     double   simple_one_layer_omega    [][1]     = {{.7}};
     double   simple_one_layer_ltau     [][1]     = {{1.}};
     double   simple_one_layer_albedo   []        = { .3 };


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     double   simple_small_stack_levels_z [][N_LAYERS_SMALL_STACK + 1] = {{20.,
                                                                           11.,
                                                                           9.,
                                                                           0.}};
     double   simple_small_stack_theta_0  []                           = { 35. };
     double   simple_small_stack_levels_0 [][N_LAYERS_SMALL_STACK + 1] = {{0.,
                                                                           0.,
                                                                           0.,
                                                                           0.}};
     double   simple_small_stack_surface_0[]                           = { 0. };
     double   simple_small_stack_levels_b [][N_LAYERS_SMALL_STACK + 1] = {{.05,
                                                                           .1,
                                                                           .15,
                                                                           .2}};

     double   simple_small_stack_surface_b[]                           = { .25};
/*
     double   simple_small_stack_surface_b[]                           = { .3 };
*/
     int      simple_small_stack_n_coefs  [][N_LAYERS_SMALL_STACK]     = {{td->n_gc[0],
                                                                           td->n_gc[2],
                                                                           td->n_gc[4]}};
     double **simple_small_stack_coefs    [][N_LAYERS_SMALL_STACK]     = {{td->  gc[0],
                                                                           td->  gc[2],
                                                                           td->  gc[4]}};
     double   simple_small_stack_omega    [][N_LAYERS_SMALL_STACK]     = {{.7,
                                                                           .8,
                                                                           .9}};
     double   simple_small_stack_ltau     [][N_LAYERS_SMALL_STACK]     = {{.01,
                                                                           .1,
                                                                           1.}};
     double   simple_small_stack_albedo   []                           = { .3 };


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     double   bounds_one_layer_levels_z [][1 + 1] = {{11.,
                                                      9.}};
     double   bounds_one_layer_theta_0  []        = { 0., 35., 89. };
     double   bounds_one_layer_levels_0 [][1 + 1] = {{0.,
                                                      0.}};
     double   bounds_one_layer_surface_0[]        = { 0. };
     double   bounds_one_layer_levels_b [][1 + 1] = {{.1,
                                                      .2}};
     double   bounds_one_layer_surface_b[]        = { .3 };
     int      bounds_one_layer_n_coefs  [][1]     = {{td->n_gc[0]}, {td->n_gc[2]}, {td->n_gc[4]}};
     double **bounds_one_layer_coefs    [][1]     = {{td->  gc[0]}, {td->  gc[2]}, {td->  gc[4]}};
     double   bounds_one_layer_omega    [][1]     = {{.0001}, {.7}, {.9999}, {1.}};
     double   bounds_one_layer_ltau     [][1]     = {{.001}, {1.}, {100.}};
     double   bounds_one_layer_albedo   []        = { 0., .3, 1. };


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     double   bounds_small_stack_levels_z [][N_LAYERS_SMALL_STACK + 1] = {{20.,
                                                                           11.,
                                                                           9.,
                                                                           0.}};
     double   bounds_small_stack_theta_0  []                           = { 0., 35., 89. };
     double   bounds_small_stack_levels_0 [][N_LAYERS_SMALL_STACK + 1] = {{0.,
                                                                           0.,
                                                                           0.,
                                                                           0.}};
     double   bounds_small_stack_surface_0[]                           =  { 0., 0., 0. };
     double   bounds_small_stack_levels_b [][N_LAYERS_SMALL_STACK + 1] = {{.05,
                                                                           .1,
                                                                           .15,
                                                                           .2}};
/*
     double   bounds_small_stack_surface_b[]                           = { .3 };
*/
     double   bounds_small_stack_surface_b[]                           = { 0., .3, 1. };
     int      bounds_small_stack_n_coefs  [][N_LAYERS_SMALL_STACK]     = {{td->n_gc[3],
                                                                           td->n_gc[0],
                                                                           td->n_gc[3]},
		                                                          {td->n_gc[3],
                                                                           td->n_gc[2],
                                                                           td->n_gc[3]},
		                                                          {td->n_gc[3],
                                                                           td->n_gc[4],
                                                                           td->n_gc[3]}};
     double **bounds_small_stack_coefs    [][N_LAYERS_SMALL_STACK]     = {{td->  gc[3],
                                                                           td->  gc[0],
                                                                           td->  gc[3]},
		                                                          {td->  gc[3],
                                                                           td->  gc[2],
                                                                           td->  gc[3]},
		                                                          {td->  gc[3],
                                                                           td->  gc[4],
                                                                           td->  gc[3]}};
     double   bounds_small_stack_omega    [][N_LAYERS_SMALL_STACK]     = {{.9,
                                                                           .0001,
                                                                           .9},
		                                                          {.9,
                                                                           .7,
                                                                           .9},
		                                                          {.9,
                                                                           .9999,
                                                                           .9},
		                                                          {.9,
                                                                           1.,
                                                                           .9}};
     double   bounds_small_stack_ltau     [][N_LAYERS_SMALL_STACK]     = {{.01,
                                                                           .001,
                                                                           .01},
		                                                          {.01,
                                                                           1.,
                                                                           .01},
		                                                          {.01,
                                                                           100.,
                                                                           .01}};
     double   bounds_small_stack_albedo   []                           =  { 0., .3, 1. };


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     double   derivs_one_layer_one_derivs_levels_b_l [][1 + 1][1] = {{{1.},
                                                                      {1.}}};
/*
     double   derivs_one_layer_one_derivs_surface_b_l[]       [1] = { {1.} };
*/
     double **derivs_one_layer_one_derivs_coefs_l    [][1    ][1] = {{{td->gc_l[2]}}};
     double   derivs_one_layer_one_derivs_omega_l    [][1    ][1] = {{{1.}}};
     double   derivs_one_layer_one_derivs_ltau_l     [][1    ][1] = {{{1.}}};
     double   derivs_one_layer_one_derivs_albedo_l   []       [1] = { {1.} };


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/

     double   derivs_one_layer_bound_derivs_levels_b_l[][1 + 1][N_DERIVS_ONE_LAYER_BOUND_DERIVS] = {{{0.,   2.,   0.,          4.,          0.,   6.,          0.,          8.,          0.,   10.,  0.,          12.,         0.,   14.,  0.,          16.,         0.  },
                                                                                                     {0.,   2.,   0.,          4.,          0.,   6.,          0.,          8.,          0.,   10.,  0.,          12.,         0.,   14.,  0.,          16.,         0.  }}};
     double **derivs_one_layer_bound_derivs_coefs_l   [][1    ][N_DERIVS_ONE_LAYER_BOUND_DERIVS] = {{{NULL, NULL, td->gc_l[2], td->gc_l[2], NULL, NULL,        td->gc_l[2], td->gc_l[2], NULL, NULL, td->gc_l[2], td->gc_l[2], NULL, NULL, td->gc_l[2], td->gc_l[2], NULL}}};
     double   derivs_one_layer_bound_derivs_omega_l   [][1    ][N_DERIVS_ONE_LAYER_BOUND_DERIVS] = {{{0.,   0.,   0.,          0.,          5.,   6.,          7.,          8.,          0.,   0.,   0.,          0.,          13.,  14.,  15.,         16.,         0.  }}};
     double   derivs_one_layer_bound_derivs_ltau_l    [][1    ][N_DERIVS_ONE_LAYER_BOUND_DERIVS] = {{{0.,   0.,   0.,          0.,          0.,   0.,          0.,          0.,          9.,   10.,  11.,         12.,         13.,  14.,  15.,         16.,         0.  }}};
     double   derivs_one_layer_bound_derivs_albedo_l  []       [N_DERIVS_ONE_LAYER_BOUND_DERIVS] = { {0.,   2.,   3.,          4.,          5.,   6.,          7.,          8.,          9.,   10.,  11.,         12.,         13.,  14.,  15.,         0.,          17. } };
/*
     double   derivs_one_layer_bound_derivs_levels_b_l [][1 + 1][N_DERIVS_ONE_LAYER_BOUND_DERIVS] = {{{0.,   1.,   0.,   3.,   0.,          5.,          0.,          7.,          0.,   9.,   0.,   11.,  0.,          13.,         0.,          15.,         0.,   17.,  0.  },
                                                                                                      {0.,   1.,   0.,   3.,   0.,          5.,          0.,          7.,          0.,   9.,   0.,   11.,  0.,          13.,         0.,          15.,         0.,   17.,  0.  }}};
     double   derivs_one_layer_bound_derivs_surface_b_l[]       [N_DERIVS_ONE_LAYER_BOUND_DERIVS] = { {0.,   0.,   2.,   3.,   0.,          0.,          6.,          7.,          0.,   0.,   10.,  11.,  0.,          0.,          14.,         15.,         0.,   0.,   18. } };
     double **derivs_one_layer_bound_derivs_coefs_l    [][1    ][N_DERIVS_ONE_LAYER_BOUND_DERIVS] = {{{NULL, NULL, NULL, NULL, td->gc_l[2], td->gc_l[2], td->gc_l[2], td->gc_l[2], NULL, NULL, NULL, NULL, td->gc_l[2], td->gc_l[2], td->gc_l[2], td->gc_l[2], NULL, NULL, NULL}}};
     double   derivs_one_layer_bound_derivs_omega_l    [][1    ][N_DERIVS_ONE_LAYER_BOUND_DERIVS] = {{{0.,   0.,   0.,   0.,   4.,          0.,          0.,          0.,          8.,   9.,   10.,  11.,  12.,         13.,         14.,         15.,         0.,   0.,   0.  }}};
     double   derivs_one_layer_bound_derivs_ltau_l     [][1    ][N_DERIVS_ONE_LAYER_BOUND_DERIVS] = {{{0.,   0.,   0.,   0.,   0.,          0.,          0.,          0.,          0.,   0.,   0.,   0.,   0.,          0.,          0.,          0.,          16.,  17.,  18. }}};
     double   derivs_one_layer_bound_derivs_albedo_l   []       [N_DERIVS_ONE_LAYER_BOUND_DERIVS] = { {0.,   0.,   0.,   0.,   0.,          0.,          0.,          0.,          0.,   0.,   0.,   0.,   0.,          0.,          0.,          0.,          0.,   0.,   0.  } };
*/

     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     double   derivs_small_stack_one_derivs_levels_b_l[][N_LAYERS_SMALL_STACK + 1][1] = {{{0.},
                                                                                          {1.},
                                                                                          {1.},
                                                                                          {0.}}};
     double **derivs_small_stack_one_derivs_coefs_l   [][N_LAYERS_SMALL_STACK    ][1] = {{{td->gc_l[0]},
                                                                                          {td->gc_l[2]},
                                                                                          {td->gc_l[0]}}};
     double   derivs_small_stack_one_derivs_omega_l   [][N_LAYERS_SMALL_STACK    ][1] = {{{0.},
                                                                                          {1.},
                                                                                          {0.}}};
     double   derivs_small_stack_one_derivs_ltau_l    [][N_LAYERS_SMALL_STACK    ][1] = {{{0.},
                                                                                          {1.},
                                                                                          {0.}}};
     double   derivs_small_stack_one_derivs_albedo_l  []                          [1] = { {1.} };


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     double   derivs_small_stack_bound_derivs_levels_b_l[][N_LAYERS_SMALL_STACK + 1][N_DERIVS_SMALL_STACK_BOUND_DERIVS] = {{{0.,          2.,          0.,          4.,          0.,           6.,          7.,           8.,          0.},
		                                                                                                            {0.,          2.,          3.,          4.,          0.,           6.,          7.,           8.,          0.},
		                                                                                                            {0.,          0.,          3.,          4.,          5.,           6.,          7.,           8.,          0.},
                                                                                                                            {0.,          0.,          0.,          0.,          5.,           6.,          7.,           8.,          0.}}};
     double **derivs_small_stack_bound_derivs_coefs_l   [][N_LAYERS_SMALL_STACK    ][N_DERIVS_SMALL_STACK_BOUND_DERIVS] = {{{NULL,        td->gc_l[0], NULL,        td->gc_l[0], NULL,         td->gc_l[0], NULL,         td->gc_l[0], NULL},
		                                                                                                            {NULL,        NULL,        td->gc_l[2], td->gc_l[2], NULL,         NULL,        td->gc_l[2],  td->gc_l[2], NULL},
		                                                                                                            {NULL,        NULL,        NULL,        NULL,        td->gc_l[4],  td->gc_l[4], td->gc_l[4],  td->gc_l[4], NULL}}};
     double   derivs_small_stack_bound_derivs_omega_l   [][N_LAYERS_SMALL_STACK    ][N_DERIVS_SMALL_STACK_BOUND_DERIVS] = {{{0.,          2.,          0.,          4.,          0.,           6.,          0.,           8.,          0.},
		                                                                                                            {0.,          0.,          3.,          4.,          0.,           0.,          7.,           8.,          0.},
		                                                                                                            {0.,          0.,          0.,          0.,          5.,           6.,          7.,           8.,          0.}}};
     double   derivs_small_stack_bound_derivs_ltau_l    [][N_LAYERS_SMALL_STACK    ][N_DERIVS_SMALL_STACK_BOUND_DERIVS] = {{{0.,          2.,          3.,          4.,          0.,           6.,          0.,           8.,          0.},
		                                                                                                            {0.,          0.,          0.,          4.,          0.,           0.,          7.,           8.,          0.},
		                                                                                                            {0.,          0.,          0.,          0.,          5.,           6.,          7.,           8.,          0.}}};
     double   derivs_small_stack_bound_derivs_albedo_l  []                          [N_DERIVS_SMALL_STACK_BOUND_DERIVS] =  {{0.,          2.,          3.,          4.,          5.,           6.,          7.,           0.,          10.}};


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     n_solvers = xrtm_solver_n();

     solvers_list = alloc_array1_l(n_solvers);


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (td->stack_type == TEST_STACK_ONE_LAYER)
          n_layers = 1;
     else
     if (td->stack_type == TEST_STACK_SMALL_STACK)
          n_layers = N_LAYERS_SMALL_STACK;
     else
     if (td->stack_type == TEST_STACK_LARGE_STACK)
          n_layers = N_LAYERS_LARGE_STACK;

     n_kernels = 1;

     kernels[0] = XRTM_KERNEL_LAMBERTIAN;

     if (td->derivs_type == TEST_DERIVS_NO_DERIVS)
          n_derivs = 0;
     else
     if (td->derivs_type == TEST_DERIVS_ONE_DERIVS)
          n_derivs = 1;
     else
     if (td->derivs_type == TEST_DERIVS_BOUND_DERIVS) {
          if (td->stack_type == TEST_STACK_ONE_LAYER)
               n_derivs = N_DERIVS_ONE_LAYER_BOUND_DERIVS;
          else
          if (td->stack_type == TEST_STACK_SMALL_STACK)
               n_derivs = N_DERIVS_SMALL_STACK_BOUND_DERIVS;
/*
          else
          if (td->stack_type == TEST_STACK_LARGE_STACK)
               n_derivs = N_DERIVS_LARGE_STACK_BOUND_DERIVS;
*/
     }

     if (td->output_at_levels_type == TEST_OUTPUT_AT_LEVELS_TOP)
          n_out_levels = 1;
     else
     if (td->output_at_levels_type == TEST_OUTPUT_AT_LEVELS_BOTTOM)
          n_out_levels = 1;
     else
     if (td->output_at_levels_type == TEST_OUTPUT_AT_LEVELS_TOP_AND_BOTTOM)
          n_out_levels = 2;
     else
     if (td->output_at_levels_type == TEST_OUTPUT_AT_LEVELS_ALL_LEVELS)
          n_out_levels = n_layers + 1;

     if (td->output_at_taus_type   == TEST_OUTPUT_AT_TAUS_TWO_END_MIDDLES) {
          if (n_layers == 1)
               n_out_taus = 1;
          else
               n_out_taus = 2;
     }
     else
     if (td->output_at_taus_type   == TEST_OUTPUT_AT_TAUS_BOUNDS_AND_MIDDLES)
          n_out_taus   = 2 + n_layers;


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (td->core_zero_solar_source)
          F_0 = 0.;
     else
          F_0 = DEFAULT_F_0;

     if (td->bound_type == TEST_BOUND_SIMPLE) {
          if (td->stack_type == TEST_STACK_ONE_LAYER) {
               n_levels_z_bounds  = LENGTH(simple_one_layer_levels_z);
               n_theta_0_bounds   = LENGTH(simple_one_layer_theta_0);
               n_levels_b_bounds  = LENGTH(simple_one_layer_levels_b);
               n_surface_b_bounds = LENGTH(simple_one_layer_surface_b);
               n_coefs_bounds     = LENGTH(simple_one_layer_n_coefs);
               n_omega_bounds     = LENGTH(simple_one_layer_omega);
               n_ltau_bounds      = LENGTH(simple_one_layer_ltau);
               n_albedo_bounds    = LENGTH(simple_one_layer_albedo);

               levels_z  = *simple_one_layer_levels_z;
               theta_0   =  simple_one_layer_theta_0;
               if (td->core_zero_thermal_source) {
                    levels_b  = *simple_one_layer_levels_0;
                    surface_b =  simple_one_layer_surface_0;
               }
               else {
                    levels_b  = *simple_one_layer_levels_b;
                    surface_b =  simple_one_layer_surface_b;
               }
               n_coefs   = *simple_one_layer_n_coefs;
               coefs     = *simple_one_layer_coefs;
               omega     = *simple_one_layer_omega;
               ltau      = *simple_one_layer_ltau;
               albedo    =  simple_one_layer_albedo;
          }
          else
          if (td->stack_type == TEST_STACK_SMALL_STACK) {
               n_levels_z_bounds  = LENGTH(simple_small_stack_levels_z);
               n_theta_0_bounds   = LENGTH(simple_small_stack_theta_0);
               n_levels_b_bounds  = LENGTH(simple_small_stack_levels_b);
               n_surface_b_bounds = LENGTH(simple_small_stack_surface_b);
               n_coefs_bounds     = LENGTH(simple_small_stack_n_coefs);
               n_omega_bounds     = LENGTH(simple_small_stack_omega);
               n_ltau_bounds      = LENGTH(simple_small_stack_ltau);
               n_albedo_bounds    = LENGTH(simple_small_stack_albedo);

               levels_z  = *simple_small_stack_levels_z;
               theta_0   =  simple_small_stack_theta_0;
               if (td->core_zero_thermal_source) {
                    levels_b  = *simple_small_stack_levels_0;
                    surface_b =  simple_small_stack_surface_0;
               }
               else {
                    levels_b  = *simple_small_stack_levels_b;
                    surface_b =  simple_small_stack_surface_b;
               }
               n_coefs   = *simple_small_stack_n_coefs;
               coefs     = *simple_small_stack_coefs;
               omega     = *simple_small_stack_omega;
               ltau      = *simple_small_stack_ltau;
               albedo    =  simple_small_stack_albedo;
          }
          else
          if (td->stack_type == TEST_STACK_LARGE_STACK)
               ;
     }
     else
     if (td->bound_type == TEST_BOUND_BOUNDS) {
          if (td->stack_type == TEST_STACK_ONE_LAYER) {
               n_levels_z_bounds  = LENGTH(bounds_one_layer_levels_z);
               n_theta_0_bounds   = LENGTH(bounds_one_layer_theta_0);
               n_levels_b_bounds  = LENGTH(bounds_one_layer_levels_b);
               n_surface_b_bounds = LENGTH(bounds_one_layer_surface_b);
               n_coefs_bounds     = LENGTH(bounds_one_layer_n_coefs);
               n_omega_bounds     = LENGTH(bounds_one_layer_omega);
               n_ltau_bounds      = LENGTH(bounds_one_layer_ltau);
               n_albedo_bounds    = LENGTH(bounds_one_layer_albedo);

               levels_z  = *bounds_one_layer_levels_z;
               theta_0   =  bounds_one_layer_theta_0;
               if (td->core_zero_thermal_source) {
                    levels_b  = *bounds_one_layer_levels_0;
                    surface_b =  bounds_one_layer_surface_0;
               }
               else {
                    levels_b  = *bounds_one_layer_levels_b;
                    surface_b =  bounds_one_layer_surface_b;
               }
               n_coefs   = *bounds_one_layer_n_coefs;
               coefs     = *bounds_one_layer_coefs;
               omega     = *bounds_one_layer_omega;
               ltau      = *bounds_one_layer_ltau;
               albedo    =  bounds_one_layer_albedo;
          }
          else
          if (td->stack_type == TEST_STACK_SMALL_STACK) {
               n_levels_z_bounds  = LENGTH(bounds_small_stack_levels_z);
               n_theta_0_bounds   = LENGTH(bounds_small_stack_theta_0);
               n_levels_b_bounds  = LENGTH(bounds_small_stack_levels_b);
               n_surface_b_bounds = LENGTH(bounds_small_stack_surface_b);
               n_coefs_bounds     = LENGTH(bounds_small_stack_n_coefs);
               n_omega_bounds     = LENGTH(bounds_small_stack_omega);
               n_ltau_bounds      = LENGTH(bounds_small_stack_ltau);
               n_albedo_bounds    = LENGTH(bounds_small_stack_albedo);

               levels_z  = *bounds_small_stack_levels_z;
               theta_0   =  bounds_small_stack_theta_0;
               if (td->core_zero_thermal_source) {
                    levels_b  = *bounds_small_stack_levels_0;
                    surface_b =  bounds_small_stack_surface_0;
               }
               else {
                    levels_b  = *bounds_small_stack_levels_b;
                    surface_b =  bounds_small_stack_surface_b;
               }
               n_coefs   = *bounds_small_stack_n_coefs;
               coefs     = *bounds_small_stack_coefs;
               omega     = *bounds_small_stack_omega;
               ltau      = *bounds_small_stack_ltau;
               albedo    =  bounds_small_stack_albedo;
          }
          else
          if (td->stack_type == TEST_STACK_LARGE_STACK)
               ;
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (td->derivs_type == TEST_DERIVS_NO_DERIVS)
          n_deriv_bounds = 0;
     else
     if (td->derivs_type == TEST_DERIVS_ONE_DERIVS) {
          n_deriv_bounds = 1;

          if (td->stack_type == TEST_STACK_ONE_LAYER) {
               levels_b_l = **derivs_one_layer_one_derivs_levels_b_l;
               coefs_l    = **derivs_one_layer_one_derivs_coefs_l;
               omega_l    = **derivs_one_layer_one_derivs_omega_l;
               ltau_l     = **derivs_one_layer_one_derivs_ltau_l;
               albedo_l   =  *derivs_one_layer_one_derivs_albedo_l;
          }
          else
          if (td->stack_type == TEST_STACK_SMALL_STACK) {
               levels_b_l = **derivs_small_stack_one_derivs_levels_b_l;
               coefs_l    = **derivs_small_stack_one_derivs_coefs_l;
               omega_l    = **derivs_small_stack_one_derivs_omega_l;
               ltau_l     = **derivs_small_stack_one_derivs_ltau_l;
               albedo_l   =  *derivs_small_stack_one_derivs_albedo_l;
          }
          else
          if (td->stack_type == TEST_STACK_LARGE_STACK)
               ;
     }
     else
     if (td->derivs_type == TEST_DERIVS_BOUND_DERIVS) {
          if (td->stack_type == TEST_STACK_ONE_LAYER) {
               levels_b_l = **derivs_one_layer_bound_derivs_levels_b_l;
               coefs_l    = **derivs_one_layer_bound_derivs_coefs_l;
               omega_l    = **derivs_one_layer_bound_derivs_omega_l;
               ltau_l     = **derivs_one_layer_bound_derivs_ltau_l;
               albedo_l   =  *derivs_one_layer_bound_derivs_albedo_l;
          }
          else
          if (td->stack_type == TEST_STACK_SMALL_STACK) {
               levels_b_l = **derivs_small_stack_bound_derivs_levels_b_l;
               coefs_l    = **derivs_small_stack_bound_derivs_coefs_l;
               omega_l    = **derivs_small_stack_bound_derivs_omega_l;
               ltau_l     = **derivs_small_stack_bound_derivs_ltau_l;
               albedo_l   =  *derivs_small_stack_bound_derivs_albedo_l;
          }
          else
          if (td->stack_type == TEST_STACK_LARGE_STACK)
               ;
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     out_levels = alloc_array1_i(n_out_levels);

     if (td->output_at_levels_type == TEST_OUTPUT_AT_LEVELS_TOP) {
          out_levels[0] = 0;
     }
     else
     if (td->output_at_levels_type == TEST_OUTPUT_AT_LEVELS_BOTTOM) {
          out_levels[0] = n_layers;
     }
     else
     if (td->output_at_levels_type == TEST_OUTPUT_AT_LEVELS_TOP_AND_BOTTOM) {
          out_levels[0] = 0;
          out_levels[1] = n_layers;
     }
     else
     if (td->output_at_levels_type == TEST_OUTPUT_AT_LEVELS_ALL_LEVELS) {
          for (i = 0; i < n_layers + 1; ++i)
               out_levels[i] = i;
     }

     out_taus = alloc_array1_d(n_out_taus);

     if (td->output_at_taus_type   == TEST_OUTPUT_AT_TAUS_TWO_END_MIDDLES) {
          ii = 0;
          out_taus[ii++] = ltau[0] / 2.;
          if (n_layers > 1) {
               a = 0.;
               for (i = 0; i < n_layers; ++i)
                    a += ltau[i];
               out_taus[ii++] = a - ltau[n_layers - 1] / 2.;
          }
     }
     else
     if (td->output_at_taus_type   == TEST_OUTPUT_AT_TAUS_BOUNDS_AND_MIDDLES) {
          ii = 0;
          out_taus[ii++] = 0.;
          a = 0.;
          for (i = 0; i < n_layers; ++i) {
               out_taus[ii++] = a + ltau[i] / 2.;
               a += ltau[i];
          }
          out_taus[ii++] = a;
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     out_thetas = alloc_array1_d(td->n_out_thetas);

     if (td->n_out_thetas == 1)
          out_thetas[0] = 45.;
     else {
          a = 89. / (td->n_out_thetas - 1);

          out_thetas[0] = 0.;

          for (i = 1; i < td->n_out_thetas; ++i)
               out_thetas[i] = out_thetas[i - 1] + a;
     }


     out_phis = alloc_array1_d(td->n_out_phis);

     if (td->n_out_phis == 1)
          out_phis[0] = 45.;
     else {
          a = 180. / (td->n_out_phis - 1);

          out_phis[0] = 0.;

          for (i = 1; i < td->n_out_phis; ++i)
               out_phis[i] = out_phis[i - 1] + a;
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     tol = alloc_array1_d(n_solvers);
     init_array1_d(tol, n_solvers, DEF_TOL);
     if (n_derivs > 0) {
          tol_l = alloc_array1_d(n_solvers);
          init_array1_d(tol_l, n_solvers, DEF_TOL_L);
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     count = 0;

     for (i = 0; i < MAX_XRTM_OPTION_MASK; ++i) {
          if (xrtm_check_options(i, 0) != 0)
               continue;

          if (i & XRTM_OPTION_REVERSE_DERIVS)
               continue;

          if (i & XRTM_OPTION_PART_SOL_GREENS)
               continue;

          if (i & XRTM_OPTION_PHASE_MATRIX_LC)
               continue;

          if (! (i & XRTM_OPTION_SOURCE_SOLAR) && (i & XRTM_OPTION_STACK_REUSE_ADDING))
               continue;

          if (td->exact_options > 0 && i != td->exact_options)
               continue;
          else {
               if ((i & td->required_options) != td->required_options)
                    continue;

               if ((i & td->unwanted_options) != 0)
                    continue;
          }


          /*--------------------------------------------------------------------
           *
           *------------------------------------------------------------------*/
          if (i & XRTM_OPTION_CALC_DERIVS &&
              td->derivs_type == TEST_DERIVS_NO_DERIVS)
               continue;

          if (i & XRTM_OPTION_TOP_DOWN_ADDING &&
              td->output_at_levels_type != TEST_OUTPUT_AT_LEVELS_BOTTOM)
               continue;
          if (i & XRTM_OPTION_BOTTOM_UP_ADDING &&
              td->output_at_levels_type != TEST_OUTPUT_AT_LEVELS_TOP)
               continue;


          /*--------------------------------------------------------------------
           *
           *------------------------------------------------------------------*/
          n_derivs2 = 0;
          if (i & XRTM_OPTION_CALC_DERIVS)
               n_derivs2 = n_derivs;

          if (i & XRTM_OPTION_OUTPUT_AT_LEVELS)
               n_out_levels2 = n_out_levels;
          else
          if (i & XRTM_OPTION_OUTPUT_AT_TAUS)
               n_out_levels2 = n_out_taus;


          if (xrtm_check_counts(i, td->max_coef, td->n_quad, td->n_stokes, n_derivs2, n_layers, 1, n_kernels, 2 * td->n_quad, kernels, n_out_levels2, td->n_out_thetas, 0) != 0)
               continue;


          /*--------------------------------------------------------------------
           *
           *------------------------------------------------------------------*/
          n_solvers2 = 0;

          solvers_mask = 0;
          for (j = 0; j < n_solvers; ++j) {
               mask = xrtm_solver_index_to_mask(j);

               if (mask & XRTM_SOLVER_TWO_STREAM)
                    continue;

               if (mask & XRTM_SOLVER_FOUR_STREAM)
                    continue;

               if (mask & XRTM_SOLVER_SIX_STREAM)
                    continue;

               if (mask & XRTM_SOLVERS_EXT_EXTERNAL && ! td->core_include_dev_solvers)
                    continue;

               if (i & XRTM_OPTION_CALC_DERIVS) {
#ifdef INCLUDE_DEV_SOURCE
                    if (mask & XRTM_SOLVER_LIDORT)
                         continue;

                    if (mask & XRTM_SOLVER_VLIDORT)
                         continue;
#endif
               }
#ifdef INCLUDE_DEV_SOURCE
               if (mask & XRTM_SOLVER_POLRAD)
                    continue;

               if (mask & XRTM_SOLVER_RADIANT)
                    continue;

               if (mask & XRTM_SOLVER_SOI)
                    continue;
#endif

               if (xrtm_check_solvers(i, mask, td->max_coef, td->n_quad, td->n_stokes, n_derivs2, n_layers, 1, n_kernels, 2 * td->n_quad, kernels, n_out_levels2, td->n_out_thetas, 0) == 0) {
                    solvers_mask |= mask;
                    n_solvers2++;
               }
          }

          if (td->exact_solvers != 0) {
               solvers_mask &= td->exact_solvers;

               if (solvers_mask != td->exact_solvers)
                    continue;
          }
          else {
               if ((solvers_mask & td->required_solvers) != td->required_solvers)
                    continue;

               solvers_mask &= ~(solvers_mask & td->unwanted_solvers);
          }

          if (solvers_mask == 0)
               continue;

          n_solvers2 = xrtm_solver_mask_to_value_list(solvers_mask, solvers_list, n_solvers2);


          /*--------------------------------------------------------------------
           *
           *------------------------------------------------------------------*/
          test_xrtm = malloc(sizeof(test_xrtm_data));
          test_xrtm_init(test_xrtm);

          list_append(&td->test_xrtm_list, test_xrtm, 0);


          /*--------------------------------------------------------------------
           *
           *------------------------------------------------------------------*/
          test_xrtm_fill(test_xrtm, i, solvers_mask, td->max_coef, td->n_quad, td->n_stokes, n_derivs2, n_layers, n_kernels, 2 * td->n_quad, kernels, n_out_levels2, td->n_out_thetas);


          /*--------------------------------------------------------------------
           *
           *------------------------------------------------------------------*/
          HANDLE_RETURN(setup_bounds_test(td,
                                          n_solvers2, solvers_list,
                                          set_xrtm_defaults,
                                          n_layers, n_derivs2,
                                          n_out_levels, out_levels,
                                          n_out_taus, out_taus,
                                          td->n_out_thetas, out_thetas,
                                          td->n_out_phis, out_phis,
                                          F_0,
                                          n_levels_z_bounds, levels_z,
                                          n_theta_0_bounds, theta_0,
                                          n_levels_b_bounds, levels_b,
                                          n_surface_b_bounds, surface_b,
                                          n_coefs_bounds, n_coefs, coefs,
                                          n_omega_bounds, omega,
                                          n_ltau_bounds, ltau,
                                          n_albedo_bounds, albedo,
                                          n_deriv_bounds,
                                          levels_b_l,
                                          coefs_l, omega_l, ltau_l,
                                          albedo_l,
                                          tol, tol_l,
                                          0,
                                          NULL,
                                          NULL,
                                          0,
                                          NULL,
                                          NULL,
                                          NULL), "bounds_test()");

          count++;
     }
/*
     printf("%d\n", count);
*/

     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     free_array1_i((int *) solvers_list);

     free_array1_i(out_levels);
     free_array1_d(out_taus);

     free_array1_d(out_thetas);
     free_array1_d(out_phis);

     free_array1_d(tol);
     if (n_derivs > 0)
          free_array1_d(tol_l);


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
/*
     if (test_loop_outer(td)) {
          fprintf(stderr, "ERROR: test_loop_outer()\n");
          return -1;
     }
*/
     if (test_loop_outer2(td)) {
          fprintf(stderr, "ERROR: test_loop_outer2()\n");
          return -1;
     }

     return 0;
}
