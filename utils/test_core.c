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
#include "test_util.h"


/*******************************************************************************
 *
 ******************************************************************************/
#define N_TEST_SIMPLE_DERIVS1		2
#define N_TEST_SIMPLE_DERIVS2		1

#define N_TEST_BOUNDS_DERIVS1		9
#define N_TEST_BOUNDS_DERIVS2		7

#define N_TEST_LAYERS			3

#define N_TEST_SIMPLE_LAYERS_DERIVS1	2
#define N_TEST_SIMPLE_LAYERS_DERIVS2	1

#define N_TEST_BOUNDS_LAYERS_DERIVS1	5
#define N_TEST_BOUNDS_LAYERS_DERIVS2	11

#define TEST_PLANET_R			6371.0071823

static       double TEST_LEVELS_Z       [1             + 1] = {     20., 10.};
static       double TEST_LEVELS_Z_LAYERS[N_TEST_LAYERS + 1] = {30., 20., 10., 0.};
/*
static const double TEST_LEVELS_Z       [1             + 1] = {     20., 10.};
static const double TEST_LEVELS_Z_LAYERS[N_TEST_LAYERS + 1] = {30., 20., 10., 0.};
*/


/*******************************************************************************
 *
 ******************************************************************************/
typedef int (*test_func     )(xrtm_data *gd, test_data *td, int n_solvers, enum xrtm_solver_mask *solvers, double *tol, double *tol_l);

typedef struct {
     int n_derivs1;
     int n_derivs2;
     int n_layers;
     int n_layers_derivs1;
     int n_layers_derivs2;

     test_func func;
     test_func func_derivs1;
     test_func func_derivs2;
     test_func func_layers;
     test_func func_layers_derivs1;
     test_func func_layers_derivs2;
} test_group_data;


typedef int (*test_func_vary)(xrtm_data *gd, test_data *td, int n_solvers, enum xrtm_solver_mask *solvers, double *tol, double *tol_l, int n_x);

typedef struct {
     int n_derivs1;
     int n_derivs2;
     int n_layers;
     int n_layers_derivs1;
     int n_layers_derivs2;

     test_func_vary func_vary_layers;
     test_func_vary func_vary_layers_derivs1;
     test_func_vary func_vary_layers_derivs2;

     test_func_vary func_vary_derivs1;
     test_func_vary func_vary_derivs2;
     test_func_vary func_vary_derivs_layers1;
     test_func_vary func_vary_derivs_layers2;

     test_func_vary func_vary_thetas;
     test_func_vary func_vary_thetas_derivs1;
     test_func_vary func_vary_thetas_layers;
     test_func_vary func_vary_thetas_layers_derivs1;
} test_group_vary_data;



/*******************************************************************************
 *
 ******************************************************************************/
static int set_xrtm_defaults(xrtm_data *d, int options, int solvers_mask) {

     int out_levels[2];

     if (solvers_mask & XRTM_SOLVERS_DOUBLING)
          XRTM_SET_DOUB_D_TAU(d, DEFAULT_DOUB_D_TAU);

     if (solvers_mask & XRTM_SOLVER_PADE_ADD)
          XRTM_SET_PADE_PARAMS(d, -1, -1);

     if (solvers_mask & XRTM_SOLVER_SOS)
          XRTM_SET_SOS_PARAMS(d, DEFAULT_SOS_PARAMS_MAX_OS, DEFAULT_SOS_PARAMS_MAX_TAU, DEFAULT_SOS_PARAMS_SOS_TOL);

     XRTM_SET_FOURIER_TOL(d, DEFAULT_FOURIER_TOL);

     XRTM_SET_F_0(d, 1.);

     XRTM_SET_PHI_0(d, 0.);

     out_levels[0] = 0;
     out_levels[1] = xrtm_get_n_layers(d);

     xrtm_set_out_levels(d, out_levels);

     if (options & XRTM_OPTION_PSA) {
          XRTM_SET_PLANET_R(d, TEST_PLANET_R);
          XRTM_SET_LEVELS_Z(d, TEST_LEVELS_Z);
     }

     return 0;
}



static int set_xrtm_defaults_layers(xrtm_data *d, int options, int solvers_mask) {

     int out_levels[2];

     if (solvers_mask & XRTM_SOLVERS_DOUBLING)
          XRTM_SET_DOUB_D_TAU(d, DEFAULT_DOUB_D_TAU);

     if (solvers_mask & XRTM_SOLVER_PADE_ADD)
          XRTM_SET_PADE_PARAMS(d, -1, -1);

     if (solvers_mask & XRTM_SOLVER_SOS)
          XRTM_SET_SOS_PARAMS(d, DEFAULT_SOS_PARAMS_MAX_OS, DEFAULT_SOS_PARAMS_MAX_TAU, DEFAULT_SOS_PARAMS_SOS_TOL);

     XRTM_SET_FOURIER_TOL(d, DEFAULT_FOURIER_TOL);

     XRTM_SET_F_0(d, 1.);

     XRTM_SET_PHI_0(d, 0.);

     out_levels[0] = 0;
     out_levels[1] = xrtm_get_n_layers(d);

     xrtm_set_out_levels(d, out_levels);

     if (options & XRTM_OPTION_PSA) {
          XRTM_SET_PLANET_R(d, TEST_PLANET_R);
          XRTM_SET_LEVELS_Z(d, TEST_LEVELS_Z_LAYERS);
     }

     return 0;
}



/*******************************************************************************
 *
 ******************************************************************************/
static int initiate_test(test_data *td, int n_solvers, enum xrtm_solver_mask *solvers_array, int options, int solvers_mask, int max_coef, int n_quad, int n_stokes, int n_derivs, int n_layers, int n_kernels, int n_kernel_quad, enum xrtm_kernel_type *kernels, int n_out_thetas, test_func func, const char *name, double *tol, double *tol_l) {

     int i;

     xrtm_data *xrtm;

     xrtm = (xrtm_data *) malloc(td->n_threads * sizeof(xrtm_data));

     for (i = 0; i < td->n_threads; ++i) {
          XRTM_CREATE(&xrtm[i], options, solvers_mask, max_coef, n_quad, n_stokes, n_derivs, n_layers, n_kernels, n_kernel_quad, kernels, 2, n_out_thetas);

          if (set_xrtm_defaults(&xrtm[i], options, solvers_mask)) {
               eprintf("ERROR: set_xrtm_defaults()\n");
               return -1;
          }
     }

     HANDLE_RETURN(func(xrtm, td, n_solvers, solvers_array, tol, tol_l), name);

     for (i = 0; i < td->n_threads; ++i)
          XRTM_DESTROY(&xrtm[i]);

     free(xrtm);

     return 0;
}



static int initiate_test_layers(test_data *td, int n_solvers, enum xrtm_solver_mask *solvers_array, int options, int solvers_mask, int max_coef, int n_quad, int n_stokes, int n_derivs, int n_layers, int n_kernels, int n_kernel_quad, enum xrtm_kernel_type *kernels, int n_out_thetas, test_func func, const char *name, double *tol, double *tol_l) {

     int i;

     xrtm_data *xrtm;

     xrtm = (xrtm_data *) malloc(td->n_threads * sizeof(xrtm_data));

     for (i = 0; i < td->n_threads; ++i) {
          XRTM_CREATE(&xrtm[i], options, solvers_mask, max_coef, n_quad, n_stokes, n_derivs, n_layers, n_kernels, n_kernel_quad, kernels, 2, n_out_thetas);

          if (set_xrtm_defaults_layers(&xrtm[i], options, solvers_mask)) {
               eprintf("ERROR: set_xrtm_defaults_layers()\n");
               return -1;
          }
     }

     HANDLE_RETURN(func(xrtm, td, n_solvers, solvers_array, tol, tol_l), name);

     for (i = 0; i < td->n_threads; ++i)
          XRTM_DESTROY(&xrtm[i]);

     free(xrtm);

     return 0;
}



static int initiate_test_vary(test_data *td, int n_solvers, enum xrtm_solver_mask *solvers_array, int options, int solvers_mask, int max_coef, int n_quad, int n_stokes, int n_derivs, int n_layers, int n_kernels, int n_kernel_quad, enum xrtm_kernel_type *kernels, int n_out_thetas, test_func_vary func, const char *name, double *tol, double *tol_l, int n_vary) {

     int i;

     xrtm_data *xrtm;

     xrtm = (xrtm_data *) malloc(td->n_threads * sizeof(xrtm_data));

     for (i = 0; i < td->n_threads; ++i) {
          XRTM_CREATE(&xrtm[i], options, solvers_mask, max_coef, n_quad, n_stokes, n_derivs, n_layers, n_kernels, n_kernel_quad, kernels, 2, n_out_thetas);

          if (set_xrtm_defaults(&xrtm[i], options, solvers_mask)) {
               eprintf("ERROR: set_xrtm_defaults()\n");
               return -1;
          }
     }

     HANDLE_RETURN(func(xrtm, td, n_solvers, solvers_array, tol, tol_l, n_vary), name);

     for (i = 0; i < td->n_threads; ++i)
          XRTM_DESTROY(&xrtm[i]);

     free(xrtm);

     return 0;
}



static int initiate_test_vary_layers(test_data *td, int n_solvers, enum xrtm_solver_mask *solvers_array, int options, int solvers_mask, int max_coef, int n_quad, int n_stokes, int n_derivs, int n_layers, int n_kernels, int n_kernel_quad, enum xrtm_kernel_type *kernels, int n_out_thetas, test_func_vary func, const char *name, double *tol, double *tol_l, int n_vary) {

     int i;

     xrtm_data *xrtm;

     xrtm = (xrtm_data *) malloc(td->n_threads * sizeof(xrtm_data));

     for (i = 0; i < td->n_threads; ++i) {
          XRTM_CREATE(&xrtm[i], options, solvers_mask, max_coef, n_quad, n_stokes, n_derivs, n_layers, n_kernels, n_kernel_quad, kernels, 2, n_out_thetas);

          if (set_xrtm_defaults_layers(&xrtm[i], options, solvers_mask)) {
               eprintf("ERROR: set_xrtm_defaults_layers()\n");
               return -1;
          }
     }

     HANDLE_RETURN(func(xrtm, td, n_solvers, solvers_array, tol, tol_l, n_vary), name);

     for (i = 0; i < td->n_threads; ++i)
          XRTM_DESTROY(&xrtm[i]);

     free(xrtm);

     return 0;
}



/*******************************************************************************
 *
 ******************************************************************************/
static int ignore_mask_solver(double theta_0, double *theta) {

     return 0;
}


/*******************************************************************************
 *
 ******************************************************************************/
static int test_simple(xrtm_data *gd, test_data *td, int n_solvers, enum xrtm_solver_mask *solvers, double *tol, double *tol_l) {

     double theta_0[]    = {35.};

     double theta  [][1] = {{15.}};

     int n_chi     [][1] = {{td->n_gc[3]}};
     double **chi  [][1] = {{td->  gc[3]}};

     double omega  [][1] = {{.7}};

     double ltau   [][1] = {{1.}};

     double albedo []    = {.3};

     double phi    []    = {45.};

     int ignore_list_index_solver[] = {-1};
     int ignore_mask_index_solver[] = {-1};

     HANDLE_RETURN(bounds_test(gd, td, n_solvers, solvers, 1, 1, LENGTH(theta_0), theta_0, LENGTH(theta), *theta, LENGTH(n_chi), *n_chi, *chi, LENGTH(omega), *omega, LENGTH(ltau), *ltau, LENGTH(albedo), albedo, LENGTH(phi), phi, tol, tol_l, LENGTH(ignore_list_index_solver) - 1, ignore_list_index_solver, ignore_mask_index_solver, ignore_mask_solver), "bounds_test()");

     return 0;
}



static int test_bounds(xrtm_data *gd, test_data *td, int n_solvers, enum xrtm_solver_mask *solvers, double *tol, double *tol_l) {

     double theta_0[]    = {0., 45., 70., 88.};

     double theta  [][1] = {{0.}, {45.}, {70.}, {88.}};

     int n_chi     [][1] = {{td->n_gc[0]}, {td->n_gc[2]}, {td->n_gc[4]}};
     double **chi  [][1] = {{td->  gc[0]}, {td->  gc[2]}, {td->  gc[4]}};

     double omega  [][1] = {{.0001}, {.7}, {.9999}, {1.}};

     double ltau   [][1] = {{.001}, {1.}, {100.}};

     double albedo []    = {0., .2, 1.};

     double phi    []    = {0., 45., 90., 135., 180.};

     int ignore_list_index_solver[] = {-1};
     int ignore_mask_index_solver[] = {-1};

     HANDLE_RETURN(bounds_test(gd, td, n_solvers, solvers, 1, 1, LENGTH(theta_0), theta_0, LENGTH(theta), *theta, LENGTH(n_chi), *n_chi, *chi, LENGTH(omega), *omega, LENGTH(ltau), *ltau, LENGTH(albedo), albedo, LENGTH(phi), phi, tol, tol_l, LENGTH(ignore_list_index_solver) - 1, ignore_list_index_solver, ignore_mask_index_solver, ignore_mask_solver), "bounds_test()");

     return 0;
}



/*******************************************************************************
 *
 ******************************************************************************/
static int test_simple_derivs1(xrtm_data *gd, test_data *td, int n_solvers, enum xrtm_solver_mask *solvers, double *tol, double *tol_l) {

     double theta_0[]    = {35.};

     double theta  [][1] = {{15.}};

     int n_chi     [][1] = {{td->n_gc[2]}};
     double **chi  [][1] = {{td->  gc[2]}};

     double omega  [][1] = {{.7}};

     double ltau   [][1] = {{1.}};

     double albedo []    = {.3};

     double phi    []    = {45.};

     double **chi_l[][1][1][N_TEST_SIMPLE_DERIVS1] = {{{{td->gc_l[2], NULL}}}};

     double omega_l  [ ][1][N_TEST_SIMPLE_DERIVS1] = {{{1., 0.}}};

     double ltau_l   [ ][1][N_TEST_SIMPLE_DERIVS1] = {{{1., 0.}}};

     double albedo_l [ ]   [N_TEST_SIMPLE_DERIVS1] = { {0., 2.} };

     int ignore_list_index_solver[] = {-1};
     int ignore_mask_index_solver[] = {-1};

     int ignore_list_deriv       [] = {-1};
     int ignore_mask_deriv       [] = {-1};

     HANDLE_RETURN(bounds_test_derivs(gd, td, n_solvers, solvers, 1, N_TEST_SIMPLE_DERIVS1, 1, LENGTH(theta_0), theta_0, LENGTH(theta), *theta, LENGTH(n_chi), *n_chi, *chi, LENGTH(omega), *omega, LENGTH(ltau), *ltau, LENGTH(albedo), albedo, LENGTH(phi), phi, LENGTH(chi_l[0]), ***chi_l, **omega_l, **ltau_l, *albedo_l, tol, tol_l, LENGTH(ignore_list_index_solver) - 1, ignore_list_index_solver, ignore_mask_index_solver, LENGTH(ignore_list_deriv) - 1, ignore_list_deriv, ignore_mask_deriv, ignore_mask_solver), "bounds_test_derivs()");

     return 0;
}



static int test_simple_derivs2(xrtm_data *gd, test_data *td, int n_solvers, enum xrtm_solver_mask *solvers, double *tol, double *tol_l) {

     double theta_0[]    = {35.};

     double theta  [][1] = {{15.}};

     int n_chi     [][1] = {{td->n_gc[2]}};
     double **chi  [][1] = {{td->  gc[2]}};

     double omega  [][1] = {{.7}};

     double ltau   [][1] = {{1.}};

     double albedo []    = {.3};

     double phi    []    = {45.};

     double **chi_l[][1][1][N_TEST_SIMPLE_DERIVS2] = {{{{td->gc_l[2]}}}};

     double omega_l  [ ][1][N_TEST_SIMPLE_DERIVS2] = {{{1.}}};

     double ltau_l   [ ][1][N_TEST_SIMPLE_DERIVS2] = {{{1.}}};

     double albedo_l [ ]   [N_TEST_SIMPLE_DERIVS2] = { {1.} };

     int ignore_list_index_solver[] = {-1};
     int ignore_mask_index_solver[] = {-1};

     int ignore_list_deriv       [] = {-1};
     int ignore_mask_deriv       [] = {-1};

     HANDLE_RETURN(bounds_test_derivs(gd, td, n_solvers, solvers, 1, N_TEST_SIMPLE_DERIVS2, 1, LENGTH(theta_0), theta_0, LENGTH(theta), *theta, LENGTH(n_chi), *n_chi, *chi, LENGTH(omega), *omega, LENGTH(ltau), *ltau, LENGTH(albedo), albedo, LENGTH(phi), phi, LENGTH(chi_l[0]), ***chi_l, **omega_l, **ltau_l, *albedo_l, tol, tol_l, LENGTH(ignore_list_index_solver) - 1, ignore_list_index_solver, ignore_mask_index_solver, LENGTH(ignore_list_deriv) - 1, ignore_list_deriv, ignore_mask_deriv, ignore_mask_solver), "bounds_test_derivs()");

     return 0;
}



static int test_bounds_derivs1(xrtm_data *gd, test_data *td, int n_solvers, enum xrtm_solver_mask *solvers, double *tol, double *tol_l) {

     double theta_0[ ]    = {0., 45., 70., 88.};

     double theta  [ ][1] = {{0.}, {45.}, {70.}, {88.}};

     int n_chi     [ ][1] = {{td->n_gc[0]}, {td->n_gc[2]}, {td->n_gc[4]}};
     double **chi  [ ][1] = {{td->  gc[0]}, {td->  gc[2]}, {td->  gc[4]}};

     double omega  [ ][1] = {{.0001}, {.7}, {.9999}, {1.}};

     double ltau   [ ][1] = {{.001}, {1.}, {100.}};

     double albedo [ ]    = {0., .2, 1.};

     double phi    [ ]    = {0., 45., 90., 135., 180.};

     double **chi_l[][1][1][N_TEST_BOUNDS_DERIVS1] = {{{{NULL, td->gc_l[0], NULL, td->gc_l[0], NULL, td->gc_l[0], NULL, td->gc_l[0], NULL}}},
                                                      {{{NULL, td->gc_l[2], NULL, td->gc_l[2], NULL, td->gc_l[2], NULL, td->gc_l[2], NULL}}},
                                                      {{{NULL, td->gc_l[4], NULL, td->gc_l[4], NULL, td->gc_l[4], NULL, td->gc_l[4], NULL}}}};

     double omega_l  [ ][1][N_TEST_BOUNDS_DERIVS1] = {{{0., 0., 2., 3., 0., 0., 6., 7., 0.}}};

     double ltau_l   [ ][1][N_TEST_BOUNDS_DERIVS1] = {{{0., 0., 0., 0., 4., 5., 6., 7., 0.}}};

     double albedo_l [ ]   [N_TEST_BOUNDS_DERIVS1] = { {0., 0., 0., 0., 0., 0., 0., 0., 8.} };

     int ignore_list_index_solver[] = {-1};
     int ignore_mask_index_solver[] = {-1};

     int ignore_list_deriv       [] = {-1};
     int ignore_mask_deriv       [] = {-1};

     HANDLE_RETURN(bounds_test_derivs(gd, td, n_solvers, solvers, 1, N_TEST_BOUNDS_DERIVS1, 1, LENGTH(theta_0), theta_0, LENGTH(theta), *theta, LENGTH(n_chi), *n_chi, *chi, LENGTH(omega), *omega, LENGTH(ltau), *ltau, LENGTH(albedo), albedo, LENGTH(phi), phi, LENGTH(chi_l[0]), ***chi_l, **omega_l, **ltau_l, *albedo_l, tol, tol_l, LENGTH(ignore_list_index_solver) - 1, ignore_list_index_solver, ignore_mask_index_solver, LENGTH(ignore_list_deriv) - 1, ignore_list_deriv, ignore_mask_deriv, ignore_mask_solver), "bounds_test_derivs()");

     return 0;
}



static int test_bounds_derivs2(xrtm_data *gd, test_data *td, int n_solvers, enum xrtm_solver_mask *solvers, double *tol, double *tol_l) {

     double theta_0[ ]    = {0., 45., 70., 88.};

     double theta  [ ][1] = {{0.}, {45.}, {70.}, {88.}};

     int n_chi     [ ][1] = {{td->n_gc[0]}, {td->n_gc[2]}, {td->n_gc[4]}};
     double **chi  [ ][1] = {{td->  gc[0]}, {td->  gc[2]}, {td->  gc[4]}};

     double omega  [ ][1] = {{.0001}, {.7}, {.9999}, {1.}};

     double ltau   [ ][1] = {{.001}, {1.}, {100.}};

     double albedo [ ]    = {0., .2, 1.};

     double phi    [ ]    = {0., 45., 90., 135., 180.};

     double **chi_l[][1][1][N_TEST_BOUNDS_DERIVS2] = {{{{td->gc_l[0],  NULL, td->gc_l[0],  NULL, td->gc_l[0],  NULL, td->gc_l[0]}}},
                                                      {{{td->gc_l[2],  NULL, td->gc_l[2],  NULL, td->gc_l[2],  NULL, td->gc_l[2]}}},
                                                      {{{td->gc_l[4],  NULL, td->gc_l[4],  NULL, td->gc_l[4],  NULL, td->gc_l[4]}}}};

     double omega_l  [ ][1][N_TEST_BOUNDS_DERIVS2] = {{{0., 10., 11.,  0., 0.,  14., 15.}}};

     double ltau_l   [ ][1][N_TEST_BOUNDS_DERIVS2] = {{{0.,  0.,  0., 12., 13., 14., 15.}}};

     double albedo_l [ ]   [N_TEST_BOUNDS_DERIVS2] = { {9., 10., 11., 12., 13., 14., 15.} };

     int ignore_list_index_solver[] = {-1};
     int ignore_mask_index_solver[] = {-1};

     int ignore_list_deriv       [] = {-1};
     int ignore_mask_deriv       [] = {-1};

     HANDLE_RETURN(bounds_test_derivs(gd, td, n_solvers, solvers, 1, N_TEST_BOUNDS_DERIVS2, 1, LENGTH(theta_0), theta_0, LENGTH(theta), *theta, LENGTH(n_chi), *n_chi, *chi, LENGTH(omega), *omega, LENGTH(ltau), *ltau, LENGTH(albedo), albedo, LENGTH(phi), phi, LENGTH(chi_l[0]), ***chi_l, **omega_l, **ltau_l, *albedo_l, tol, tol_l, LENGTH(ignore_list_index_solver) - 1, ignore_list_index_solver, ignore_mask_index_solver, LENGTH(ignore_list_deriv) - 1, ignore_list_deriv, ignore_mask_deriv, ignore_mask_solver), "bounds_test_derivs()");

     return 0;
}



/*******************************************************************************
 *
 ******************************************************************************/
static int test_simple_layers(xrtm_data *gd, test_data *td, int n_solvers, enum xrtm_solver_mask *solvers, double *tol, double *tol_l) {

     double theta_0[]                = {35.};

     double theta  [][1]             = {{15.}};

     int n_chi     [][N_TEST_LAYERS] = {{td->n_gc[0], td->n_gc[2], td->n_gc[4]}};

     double **chi  [][N_TEST_LAYERS] = {{td->  gc[0], td->  gc[2], td->  gc[4]}};

     double omega  [][N_TEST_LAYERS] = {{.7,  .8, .9}};

     double ltau   [][N_TEST_LAYERS] = {{.01,  .1, 1.}};

     double albedo []                = {.3};

     double phi    []                = {45.};

     int ignore_list_index_solver[] = {-1};
     int ignore_mask_index_solver[] = {-1};

     HANDLE_RETURN(bounds_test(gd, td, n_solvers, solvers, N_TEST_LAYERS, 1, LENGTH(theta_0), theta_0, LENGTH(theta), *theta, LENGTH(n_chi), *n_chi, *chi, LENGTH(omega), *omega, LENGTH(ltau), *ltau, LENGTH(albedo), albedo, LENGTH(phi), phi, tol, tol_l, LENGTH(ignore_list_index_solver) - 1, ignore_list_index_solver, ignore_mask_index_solver, ignore_mask_solver), "bounds_test()");

     return 0;
}



static int test_bounds_layers(xrtm_data *gd, test_data *td, int n_solvers, enum xrtm_solver_mask *solvers, double *tol, double *tol_l) {

     double theta_0[]                = {0., 45., 70., 88.};

     double theta  [][1]             = {{0.}, {45.}, {70.}, {88.}};

     int n_chi     [][N_TEST_LAYERS] = {{td->n_gc[3], td->n_gc[0], td->n_gc[3]},
                                        {td->n_gc[3], td->n_gc[2], td->n_gc[3]},
                                        {td->n_gc[3], td->n_gc[4], td->n_gc[3]}};

     double **chi  [][N_TEST_LAYERS] = {{td->  gc[3], td->  gc[0], td->  gc[3]},
                                        {td->  gc[3], td->  gc[2], td->  gc[3]},
                                        {td->  gc[3], td->  gc[4], td->  gc[3]}};

     double omega  [][N_TEST_LAYERS] = {{.9,  .0001, .9},
                                        {.9,  .7,    .9},
                                        {.9,  .9999, .9},
                                        {.9, 1.,     .9}};

     double ltau   [][N_TEST_LAYERS] = {{.01,    .001, .01},
                                        {.01,   1.,    .01},
                                        {.01, 100.,    .01}};

     double albedo []                = {0., .2, 1.};

     double phi    []                = {0., 45., 90., 135., 180.};

     int ignore_list_index_solver[] = {-1};
     int ignore_mask_index_solver[] = {-1};

     HANDLE_RETURN(bounds_test(gd, td, n_solvers, solvers, N_TEST_LAYERS, 1, LENGTH(theta_0), theta_0, LENGTH(theta), *theta, LENGTH(n_chi), *n_chi, *chi, LENGTH(omega), *omega, LENGTH(ltau), *ltau, LENGTH(albedo), albedo, LENGTH(phi), phi, tol, tol_l, LENGTH(ignore_list_index_solver) - 1, ignore_list_index_solver, ignore_mask_index_solver, ignore_mask_solver), "bounds_test()");

     return 0;
}



/*******************************************************************************
 *
 ******************************************************************************/
static int test_simple_layers_derivs1(xrtm_data *gd, test_data *td, int n_solvers, enum xrtm_solver_mask *solvers, double *tol, double *tol_l) {

     double theta_0[]                = {35.};

     double theta  [][1]             = {{15.}};

     int n_chi     [][N_TEST_LAYERS] = {{td->n_gc[0], td->n_gc[2], td->n_gc[4]}};

     double **chi  [][N_TEST_LAYERS] = {{td->  gc[0], td->  gc[2], td->  gc[4]}};

     double omega  [][N_TEST_LAYERS] = {{.7,  .8, .9}};

     double ltau   [][N_TEST_LAYERS] = {{.01,  .1, 1.}};

     double albedo []                = {.3};

     double phi    []                = {45.};

     double **chi_l[1][1][N_TEST_LAYERS][N_TEST_SIMPLE_LAYERS_DERIVS1] = {{{{NULL,        NULL},
                                                                            {td->gc_l[2], NULL},
                                                                            {NULL,        NULL}}}};

     double omega_l   [ ][N_TEST_LAYERS][N_TEST_SIMPLE_LAYERS_DERIVS1] = {{{0., 0.},
                                                                           {1., 0.},
                                                                           {0., 0.}}};

     double ltau_l    [ ][N_TEST_LAYERS][N_TEST_SIMPLE_LAYERS_DERIVS1] = {{{0., 0.},
                                                                           {1., 0.},
                                                                           {0., 0.}}};

     double albedo_l  [ ]               [N_TEST_SIMPLE_LAYERS_DERIVS1] =  {{0., 1.}};

     int ignore_list_index_solver[] = {-1};
     int ignore_mask_index_solver[] = {-1};

     int ignore_list_deriv       [] = {-1};
     int ignore_mask_deriv       [] = {-1};

     HANDLE_RETURN(bounds_test_derivs(gd, td, n_solvers, solvers, N_TEST_LAYERS, N_TEST_SIMPLE_LAYERS_DERIVS1, 1, LENGTH(theta_0), theta_0, LENGTH(theta), *theta, LENGTH(n_chi), *n_chi, *chi, LENGTH(omega), *omega, LENGTH(ltau), *ltau, LENGTH(albedo), albedo, LENGTH(phi), phi, LENGTH(chi_l[0]), ***chi_l, **omega_l, **ltau_l, *albedo_l, tol, tol_l, LENGTH(ignore_list_index_solver) - 1, ignore_list_index_solver, ignore_mask_index_solver, LENGTH(ignore_list_deriv) - 1, ignore_list_deriv, ignore_mask_deriv, ignore_mask_solver), "bounds_test_derivs()");

     return 0;
}



static int test_simple_layers_derivs2(xrtm_data *gd, test_data *td, int n_solvers, enum xrtm_solver_mask *solvers, double *tol, double *tol_l) {

     double theta_0[]                = {35.};

     double theta  [][1]             = {{15.}};

     int n_chi     [][N_TEST_LAYERS] = {{td->n_gc[0], td->n_gc[2], td->n_gc[4]}};

     double **chi  [][N_TEST_LAYERS] = {{td->  gc[0], td->  gc[2], td->  gc[4]}};

     double omega  [][N_TEST_LAYERS] = {{.7,  .8, .9}};

     double ltau   [][N_TEST_LAYERS] = {{.01,  .1, 1.}};

     double albedo []                = {.3};

     double phi    []                = {45.};

     double **chi_l[1][1][N_TEST_LAYERS][N_TEST_SIMPLE_LAYERS_DERIVS2] = {{{{NULL       },
                                                                            {td->gc_l[2]},
                                                                            {NULL       }}}};

     double omega_l   [ ][N_TEST_LAYERS][N_TEST_SIMPLE_LAYERS_DERIVS2] = {{{0.},
                                                                           {1.},
                                                                           {0.}}};

     double ltau_l    [ ][N_TEST_LAYERS][N_TEST_SIMPLE_LAYERS_DERIVS2] = {{{0.},
                                                                           {1.},
                                                                           {0.}}};

     double albedo_l  [ ]               [N_TEST_SIMPLE_LAYERS_DERIVS2] =  {{1.}};

     int ignore_list_index_solver[] = {-1};
     int ignore_mask_index_solver[] = {-1};

     int ignore_list_deriv       [] = {-1};
     int ignore_mask_deriv       [] = {-1};

     HANDLE_RETURN(bounds_test_derivs(gd, td, n_solvers, solvers, N_TEST_LAYERS, N_TEST_SIMPLE_LAYERS_DERIVS2, 1, LENGTH(theta_0), theta_0, LENGTH(theta), *theta, LENGTH(n_chi), *n_chi, *chi, LENGTH(omega), *omega, LENGTH(ltau), *ltau, LENGTH(albedo), albedo, LENGTH(phi), phi, LENGTH(chi_l[0]), ***chi_l, **omega_l, **ltau_l, *albedo_l, tol, tol_l, LENGTH(ignore_list_index_solver) - 1, ignore_list_index_solver, ignore_mask_index_solver, LENGTH(ignore_list_deriv) - 1, ignore_list_deriv, ignore_mask_deriv, ignore_mask_solver), "bounds_test_derivs()");

     return 0;
}



static int test_bounds_layers_derivs1(xrtm_data *gd, test_data *td, int n_solvers, enum xrtm_solver_mask *solvers, double *tol, double *tol_l) {

     double theta_0[]                = {0., 45., 70., 88.};

     double theta  [][1]             = {{0.}, {45.}, {70.}, {88.}};

     int n_chi     [][N_TEST_LAYERS] = {{td->n_gc[3], td->n_gc[0], td->n_gc[3]},
                                        {td->n_gc[3], td->n_gc[2], td->n_gc[3]},
                                        {td->n_gc[3], td->n_gc[4], td->n_gc[3]}};

     double **chi  [][N_TEST_LAYERS] = {{td->gc  [3], td->gc  [0], td->gc  [3]},
                                        {td->gc  [3], td->gc  [2], td->gc  [3]},
                                        {td->gc  [3], td->gc  [4], td->gc  [3]}};

     double omega  [][N_TEST_LAYERS] = {{.9,  .0001, .9},
                                        {.9,  .7,    .9},
                                        {.9,  .9999, .9},
                                        {.9, 1.,     .9}};

     double ltau   [][N_TEST_LAYERS] = {{.01,    .001, .01},
                                        {.01,   1.,    .01},
                                        {.01, 100.,    .01}};

     double albedo []                = {0., .2, 1.};

     double phi    []                = {0., 45., 90., 135., 180.};

     double **chi_l[3][1][N_TEST_LAYERS][N_TEST_BOUNDS_LAYERS_DERIVS1] = {{{{NULL, td->gc_l[3], NULL,         NULL,         NULL},
                                                                            {NULL, NULL,        td->gc_l[0],  NULL,         NULL},
                                                                            {NULL, NULL,        NULL,         td->gc_l[3],  NULL}}},
                                                                          {{{NULL, td->gc_l[3], NULL,         NULL,         NULL},
                                                                            {NULL, NULL,        td->gc_l[2],  NULL,         NULL},
                                                                            {NULL, NULL,        NULL,         td->gc_l[3],  NULL}}},
                                                                          {{{NULL, td->gc_l[3], NULL,         NULL,         NULL},
                                                                            {NULL, NULL,        td->gc_l[4],  NULL,         NULL},
                                                                            {NULL, NULL,        NULL,         td->gc_l[3],  NULL}}}};

     double omega_l   [ ][N_TEST_LAYERS][N_TEST_BOUNDS_LAYERS_DERIVS1] = {{{0., 1., 0., 0., 0.},
                                                                           {0., 0., 2., 0., 0.},
                                                                           {0., 0., 0., 4., 0.}}};

     double ltau_l    [ ][N_TEST_LAYERS][N_TEST_BOUNDS_LAYERS_DERIVS1] = {{{0., 1., 0., 0., 0.},
                                                                           {0., 0., 2., 0., 0.},
                                                                           {0., 0., 0., 4., 0.}}};

     double albedo_l  [ ]               [N_TEST_BOUNDS_LAYERS_DERIVS1] =  {{0., 0., 0., 0., 8.}};

     int ignore_list_index_solver[] = {-1};
     int ignore_mask_index_solver[] = {-1};

     int ignore_list_deriv       [] = {-1};
     int ignore_mask_deriv       [] = {-1};

     HANDLE_RETURN(bounds_test_derivs(gd, td, n_solvers, solvers, N_TEST_LAYERS, N_TEST_BOUNDS_LAYERS_DERIVS1, 1, LENGTH(theta_0), theta_0, LENGTH(theta), *theta, LENGTH(n_chi), *n_chi, *chi, LENGTH(omega), *omega, LENGTH(ltau), *ltau, LENGTH(albedo), albedo, LENGTH(phi), phi, LENGTH(chi_l[0]), ***chi_l, **omega_l, **ltau_l, *albedo_l, tol, tol_l, LENGTH(ignore_list_index_solver) - 1, ignore_list_index_solver, ignore_mask_index_solver, LENGTH(ignore_list_deriv) - 1, ignore_list_deriv, ignore_mask_deriv, ignore_mask_solver), "bounds_test_derivs()");

     return 0;
}



static int test_bounds_layers_derivs2(xrtm_data *gd, test_data *td, int n_solvers, enum xrtm_solver_mask *solvers, double *tol, double *tol_l) {

     double theta_0[]                = {0., 45., 70., 88.};

     double theta  [][1]             = {{0.}, {45.}, {70.}, {88.}};

     int n_chi     [][N_TEST_LAYERS] = {{td->n_gc[3], td->n_gc[0], td->n_gc[3]},
                                        {td->n_gc[3], td->n_gc[2], td->n_gc[3]},
                                        {td->n_gc[3], td->n_gc[4], td->n_gc[3]}};

     double **chi  [][N_TEST_LAYERS] = {{td->gc  [3], td->gc  [0], td->gc  [3]},
                                        {td->gc  [3], td->gc  [2], td->gc  [3]},
                                        {td->gc  [3], td->gc  [4], td->gc  [3]}};

     double omega  [][N_TEST_LAYERS] = {{.9,  .0001, .9},
                                        {.9,  .7,    .9},
                                        {.9,  .9999, .9},
                                        {.9, 1.,     .9}};

     double ltau   [][N_TEST_LAYERS] = {{.01,    .001, .01},
                                        {.01,   1.,    .01},
                                        {.01, 100.,    .01}};

     double albedo []                = {0., .2, 1.};

     double phi    []                = {0., 45., 90., 135., 180.};

     double **chi_l[3][1][N_TEST_LAYERS][N_TEST_BOUNDS_LAYERS_DERIVS2] = {{{{td->gc_l[3], td->gc_l[3], NULL,        td->gc_l[3], td->gc_l[3],  NULL,        td->gc_l[3],  NULL,        td->gc_l[3], NULL,        td->gc_l[3]},
                                                                            {td->gc_l[0], NULL,        td->gc_l[0], td->gc_l[0], NULL,         td->gc_l[0], td->gc_l[0],  NULL,        NULL,        td->gc_l[0], td->gc_l[0]},
                                                                            {NULL,        td->gc_l[3], td->gc_l[3], td->gc_l[3], NULL,         NULL,        NULL,         td->gc_l[3], td->gc_l[3], td->gc_l[3], td->gc_l[3]}}},
                                                                          {{{td->gc_l[3], td->gc_l[3], NULL,        td->gc_l[3], td->gc_l[3],  NULL,        td->gc_l[3],  NULL,        td->gc_l[3], NULL,        td->gc_l[3]},
                                                                            {td->gc_l[2], NULL,        td->gc_l[2], td->gc_l[2], NULL,         td->gc_l[2], td->gc_l[2],  NULL,        NULL,        td->gc_l[2], td->gc_l[2]},
                                                                            {NULL,        td->gc_l[3], td->gc_l[3], td->gc_l[3], NULL,         NULL,        NULL,         td->gc_l[3], td->gc_l[3], td->gc_l[3], td->gc_l[3]}}},
                                                                          {{{td->gc_l[3], td->gc_l[3], NULL,        td->gc_l[3], td->gc_l[3],  NULL,        td->gc_l[3],  NULL,        td->gc_l[3], NULL,        td->gc_l[3]},
                                                                            {td->gc_l[4], NULL,        td->gc_l[4], td->gc_l[4], NULL,         td->gc_l[4], td->gc_l[4],  NULL,        NULL,        td->gc_l[4], td->gc_l[4]},
                                                                            {NULL,        td->gc_l[3], td->gc_l[3], td->gc_l[3], NULL,         NULL,        NULL,         td->gc_l[3], td->gc_l[3], td->gc_l[3], td->gc_l[3]}}}};

     double omega_l   [ ][N_TEST_LAYERS][N_TEST_BOUNDS_LAYERS_DERIVS2] = {{{3., 5., 0., 7., 9.,  0., 11.,  0., 13.,  0., 15.},
                                                                           {3., 0., 6., 7., 0., 10., 11.,  0.,  0., 14., 15.},
                                                                           {0., 5., 6., 7., 0.,  0.,  0., 12., 13., 14., 15.}}};

     double ltau_l    [ ][N_TEST_LAYERS][N_TEST_BOUNDS_LAYERS_DERIVS2] = {{{3., 5., 0., 7., 9.,  0., 11.,  0., 13.,  0., 15.},
                                                                           {3., 0., 6., 7., 0., 10., 11.,  0.,  0., 14., 15.},
                                                                           {0., 5., 6., 7., 0.,  0.,  0., 12., 13., 14., 15.}}};

     double albedo_l  [ ]               [N_TEST_BOUNDS_LAYERS_DERIVS2] =  {{0., 0., 0., 0., 9., 10., 11., 12., 13., 14., 15.}};

     int ignore_list_index_solver[] = {-1};
     int ignore_mask_index_solver[] = {-1};

     int ignore_list_deriv       [] = {-1};
     int ignore_mask_deriv       [] = {-1};

     HANDLE_RETURN(bounds_test_derivs(gd, td, n_solvers, solvers, N_TEST_LAYERS, N_TEST_BOUNDS_LAYERS_DERIVS2, 1, LENGTH(theta_0), theta_0, LENGTH(theta), *theta, LENGTH(n_chi), *n_chi, *chi, LENGTH(omega), *omega, LENGTH(ltau), *ltau, LENGTH(albedo), albedo, LENGTH(phi), phi, LENGTH(chi_l[0]), ***chi_l, **omega_l, **ltau_l, *albedo_l, tol, tol_l, LENGTH(ignore_list_index_solver) - 1, ignore_list_index_solver, ignore_mask_index_solver, LENGTH(ignore_list_deriv) - 1, ignore_list_deriv, ignore_mask_deriv, ignore_mask_solver), "bounds_test_derivs()");

     return 0;
}



/*******************************************************************************
 *
 ******************************************************************************/
static int test_group(test_data *td, test_group_data *tcd, int options_mask, int solvers_mask1, int n_quad, int n_stokes) {

     int options;

     int n_solvers;

     int solvers_mask2;

     enum xrtm_kernel_type kernels;

     enum xrtm_solver_mask solvers_array2[N_XRTM_SOLVERS];

     double tol  [N_XRTM_SOLVERS];
     double tol_l[N_XRTM_SOLVERS];


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     options = options_mask | XRTM_OPTION_DELTA_M;

     if ((n_solvers = make_solvers(solvers_mask1, &solvers_mask2, solvers_array2, tol,
                                    GROUP_SOLVERS_DEV
                                    XRTM_SOLVER_DOUB_ADD, DEF_TOL,
                                    XRTM_SOLVER_EIG_ADD,  DEF_TOL,
                                    XRTM_SOLVER_EIG_BVP,  DEF_TOL,
                                    XRTM_SOLVER_MEM_BVP,  DEF_TOL,
                                    XRTM_SOLVER_PADE_ADD, 1.e-3,
                                    XRTM_SOLVER_SOS,      DEF_TOL, 0)) < 0) {
          eprintf("ERROR: make_solvers()\n");
          return -1;
     }

     kernels = XRTM_KERNEL_LAMBERTIAN;

     if (initiate_test(td, n_solvers, solvers_array2, options, solvers_mask2, td->max_coef, n_quad, n_stokes, 0, 1, 1, 32, &kernels, 1, tcd->func, "test_group_data->func()", tol, NULL)) {
          eprintf("ERROR: initiate_test()\n");
          return -1;
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     options = options_mask | XRTM_OPTION_CALC_DERIVS | XRTM_OPTION_DELTA_M;

     if ((n_solvers = make_solvers2(solvers_mask1, &solvers_mask2, solvers_array2, tol, tol_l,
                                    GROUP_SOLVERS_DEV_DERIVS1
                                    XRTM_SOLVER_DOUB_ADD, DEF_TOL, DEF_TOL,
                                    XRTM_SOLVER_EIG_ADD,  DEF_TOL, DEF_TOL,
                                    XRTM_SOLVER_EIG_BVP,  DEF_TOL, DEF_TOL,
                                    XRTM_SOLVER_MEM_BVP,  DEF_TOL, DEF_TOL,
                                    XRTM_SOLVER_PADE_ADD, 1.e-3, 1.e-2,
                                    XRTM_SOLVER_SOS,      DEF_TOL, DEF_TOL, 0)) < 0) {
          eprintf("ERROR: make_solvers2()\n");
          return -1;
     }

     kernels = XRTM_KERNEL_LAMBERTIAN;

     if (initiate_test(td, n_solvers, solvers_array2, options, solvers_mask2, td->max_coef, n_quad, n_stokes, tcd->n_derivs1, 1, 1, 32, &kernels, 1, tcd->func_derivs1, "test_group_data->func_derivs1()", tol, tol_l)) {
          eprintf("ERROR: initiate_test()\n");
          return -1;
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     options = options_mask | XRTM_OPTION_CALC_DERIVS | XRTM_OPTION_DELTA_M;

     if ((n_solvers = make_solvers2(solvers_mask1, &solvers_mask2, solvers_array2, tol, tol_l,
                                    XRTM_SOLVER_DOUB_ADD, DEF_TOL, DEF_TOL,
                                    XRTM_SOLVER_EIG_ADD,  DEF_TOL, DEF_TOL,
                                    XRTM_SOLVER_EIG_BVP,  DEF_TOL, DEF_TOL,
                                    XRTM_SOLVER_MEM_BVP,  DEF_TOL, DEF_TOL,
                                    XRTM_SOLVER_PADE_ADD, 1.e-3, 1.e-2,
                                    XRTM_SOLVER_SOS,      DEF_TOL, DEF_TOL, 0)) < 0) {
          eprintf("ERROR: make_solvers2()\n");
          return -1;
     }

     kernels = XRTM_KERNEL_LAMBERTIAN;

     if (initiate_test(td, n_solvers, solvers_array2, options, solvers_mask2, td->max_coef, n_quad, n_stokes, tcd->n_derivs2, 1, 1, 32, &kernels, 1, tcd->func_derivs2, "test_group_data->func_derivs2()", tol, tol_l)) {
          eprintf("ERROR: initiate_test()\n");
          return -1;
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     options = options_mask | XRTM_OPTION_DELTA_M;

     if ((n_solvers = make_solvers(solvers_mask1, &solvers_mask2, solvers_array2, tol,
                                    GROUP_SOLVERS_DEV
                                    XRTM_SOLVER_DOUB_ADD, DEF_TOL,
                                    XRTM_SOLVER_EIG_ADD,  DEF_TOL,
                                    XRTM_SOLVER_EIG_BVP,  DEF_TOL,
                                    XRTM_SOLVER_MEM_BVP,  DEF_TOL,
                                    XRTM_SOLVER_PADE_ADD, 1.e-3,
                                    XRTM_SOLVER_SOS,      DEF_TOL, 0)) < 0) {
          eprintf("ERROR: make_solvers()\n");
          return -1;
     }

     kernels = XRTM_KERNEL_LAMBERTIAN;

     if (initiate_test_layers(td, n_solvers, solvers_array2, options, solvers_mask2, td->max_coef, n_quad, n_stokes, 0, tcd->n_layers, 1, 32, &kernels, 1, tcd->func_layers, "test_group_data->func_layers()", tol, NULL)) {
          eprintf("ERROR: initiate_test_layers()\n");
          return -1;
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     options = options_mask | XRTM_OPTION_CALC_DERIVS | XRTM_OPTION_DELTA_M;

     if ((n_solvers = make_solvers2(solvers_mask1, &solvers_mask2, solvers_array2, tol, tol_l,
                                    GROUP_SOLVERS_DEV_DERIVS1
                                    XRTM_SOLVER_DOUB_ADD, DEF_TOL, DEF_TOL,
                                    XRTM_SOLVER_EIG_ADD,  DEF_TOL, DEF_TOL,
                                    XRTM_SOLVER_EIG_BVP,  DEF_TOL, DEF_TOL,
                                    XRTM_SOLVER_MEM_BVP,  DEF_TOL, DEF_TOL,
                                    XRTM_SOLVER_PADE_ADD, 1.e-3, 1.e-2,
                                    XRTM_SOLVER_SOS,      DEF_TOL, DEF_TOL, 0)) < 0) {
          eprintf("ERROR: make_solvers2()\n");
          return -1;
     }

     kernels = XRTM_KERNEL_LAMBERTIAN;

     if (initiate_test_layers(td, n_solvers, solvers_array2, options, solvers_mask2, td->max_coef, n_quad, n_stokes, tcd->n_layers_derivs1, tcd->n_layers, 1, 32, &kernels, 1, tcd->func_layers_derivs1, "test_group_data->func_layers_derivs1()", tol, tol_l)) {
          eprintf("ERROR: initiate_test_layers()\n");
          return -1;
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     options = options_mask | XRTM_OPTION_CALC_DERIVS | XRTM_OPTION_DELTA_M;

     if ((n_solvers = make_solvers2(solvers_mask1, &solvers_mask2, solvers_array2, tol, tol_l,
                                    XRTM_SOLVER_DOUB_ADD, DEF_TOL, DEF_TOL,
                                    XRTM_SOLVER_EIG_ADD,  DEF_TOL, DEF_TOL,
                                    XRTM_SOLVER_EIG_BVP,  DEF_TOL, DEF_TOL,
                                    XRTM_SOLVER_MEM_BVP,  DEF_TOL, DEF_TOL,
                                    XRTM_SOLVER_PADE_ADD, 1.e-3, 1.e-2,
                                    XRTM_SOLVER_SOS,      DEF_TOL, DEF_TOL, 0)) < 0) {
          eprintf("ERROR: make_solvers2()\n");
          return -1;
     }

     kernels = XRTM_KERNEL_LAMBERTIAN;

     if (initiate_test_layers(td, n_solvers, solvers_array2, options, solvers_mask2, td->max_coef, n_quad, n_stokes, tcd->n_layers_derivs2, tcd->n_layers, 1, 32, &kernels, 1, tcd->func_layers_derivs2, "test_group_data->func_layers_derivs2()", tol, tol_l)) {
          eprintf("ERROR: initiate_test_layers()\n");
          return -1;
     }


     return 0;
}



/*******************************************************************************
 *
 ******************************************************************************/
static int test_group2(test_data *td, test_group_data *tcd, int option_mask, int solvers_mask1, int n_quad, int n_quad_stokes) {

     int solvers_mask2;

     if (test_group(td, tcd, option_mask, solvers_mask1, n_quad, 1)) {
          eprintf("ERROR: test_group()\n");
          return -1;
     }

     solvers_mask2 = solvers_mask1 & XRTM_SOLVERS_VECTOR;

     if (test_group(td, tcd, option_mask | XRTM_OPTION_PHASE_MATRIX_GC | XRTM_OPTION_VECTOR, solvers_mask2, n_quad_stokes, 4)) {
          eprintf("ERROR: test_group()\n");
          return -1;
     }

     return 0;
}


/*******************************************************************************
 *
 ******************************************************************************/
static int test_group_init_simple(test_group_data *d) {

     d->n_derivs1           = N_TEST_SIMPLE_DERIVS1;
     d->n_derivs2           = N_TEST_SIMPLE_DERIVS2;
     d->n_layers            = N_TEST_LAYERS;
     d->n_layers_derivs1    = N_TEST_SIMPLE_LAYERS_DERIVS1;
     d->n_layers_derivs2    = N_TEST_SIMPLE_LAYERS_DERIVS2;

     d->func                = test_simple;
     d->func_derivs1        = test_simple_derivs1;
     d->func_derivs2        = test_simple_derivs2;
     d->func_layers         = test_simple_layers;
     d->func_layers_derivs1 = test_simple_layers_derivs1;
     d->func_layers_derivs2 = test_simple_layers_derivs2;

     return 0;
}



static int test_group_init_bounds(test_group_data *d) {

     d->n_derivs1           = N_TEST_BOUNDS_DERIVS1;
     d->n_derivs2           = N_TEST_BOUNDS_DERIVS2;
     d->n_layers            = N_TEST_LAYERS;
     d->n_layers_derivs1    = N_TEST_BOUNDS_LAYERS_DERIVS1;
     d->n_layers_derivs2    = N_TEST_BOUNDS_LAYERS_DERIVS2;

     d->func                = test_bounds;
     d->func_derivs1        = test_bounds_derivs1;
     d->func_derivs2        = test_bounds_derivs2;
     d->func_layers         = test_bounds_layers;
     d->func_layers_derivs1 = test_bounds_layers_derivs1;
     d->func_layers_derivs2 = test_bounds_layers_derivs2;

     return 0;
}



/*******************************************************************************
 *
 ******************************************************************************/
static int test_simple_vary_derivs1(xrtm_data *gd, test_data *td, int n_solvers, enum xrtm_solver_mask *solvers, double *tol, double *tol_l, int n_derivs) {

     int i;
     int n;

     double theta_0[]    = {35.};

     double theta  [][1] = {{15.}};

     int n_chi     [][1] = {{td->n_gc[2]}};
     double **chi  [][1] = {{td->  gc[2]}};

     double omega  [][1] = {{.7}};

     double ltau   [][1] = {{1.}};

     double albedo []    = {.3};

     double phi    []    = {45.};

     double ******chi_l = (double ******) alloc_array4(1, 1, 1, n_derivs, sizeof(double **));

     double ***omega_l  = alloc_array3_d(1, 1, n_derivs);

     double ***ltau_l   = alloc_array3_d(1, 1, n_derivs);

     double **albedo_l  = alloc_array2_d(1, n_derivs);

     int ignore_list_index_solver[] = {-1};
     int ignore_mask_index_solver[] = {-1};

     int ignore_list_deriv       [] = {-1};
     int ignore_mask_deriv       [] = {-1};

     n = (n_derivs - 1) / 2 + 1;

     for (i = 0; i < n; ++i) {
          chi_l  [0][0][0][i] = td->gc_l[2];
          omega_l   [0][0][i] = i;
          ltau_l    [0][0][i] = i;

          albedo_l  [0]   [i] = 0.;
     }

     for ( ; i < n_derivs; ++i) {
          chi_l  [0][0][0][i] = NULL;
          omega_l   [0][0][i] = 0.;
          ltau_l    [0][0][i] = 0.;

          albedo_l  [0]   [i] = i;
     }

     HANDLE_RETURN(bounds_test_derivs(gd, td, n_solvers, solvers, 1, n_derivs, 1, LENGTH(theta_0), theta_0, LENGTH(theta), *theta, LENGTH(n_chi), *n_chi, *chi, LENGTH(omega), *omega, LENGTH(ltau), *ltau, LENGTH(albedo), albedo, LENGTH(phi), phi, LENGTH(chi_l[0]), ***chi_l, **omega_l, **ltau_l, *albedo_l, tol, tol_l, LENGTH(ignore_list_index_solver) - 1, ignore_list_index_solver, ignore_mask_index_solver, LENGTH(ignore_list_deriv) - 1, ignore_list_deriv, ignore_mask_deriv, ignore_mask_solver), "bounds_test_derivs()");

     free_array4((void ****) chi_l);
     free_array3_d(omega_l);
     free_array3_d(ltau_l);
     free_array2_d(albedo_l);

     return 0;
}



static int test_simple_vary_derivs2(xrtm_data *gd, test_data *td, int n_solvers, enum xrtm_solver_mask *solvers, double *tol, double *tol_l, int n_derivs) {

     int i;

     double theta_0[]    = {35.};

     double theta  [][1] = {{15.}};

     int n_chi     [][1] = {{td->n_gc[2]}};
     double **chi  [][1] = {{td->  gc[2]}};

     double omega  [][1] = {{.7}};

     double ltau   [][1] = {{1.}};

     double albedo []    = {.3};

     double phi    []    = {45.};

     double ******chi_l = (double ******) alloc_array4(1, 1, 1, n_derivs, sizeof(double **));

     double ***omega_l  = alloc_array3_d(1, 1, n_derivs);

     double ***ltau_l   = alloc_array3_d(1, 1, n_derivs);

     double **albedo_l  = alloc_array2_d(1, n_derivs);

     int ignore_list_index_solver[] = {-1};
     int ignore_mask_index_solver[] = {-1};

     int ignore_list_deriv       [] = {-1};
     int ignore_mask_deriv       [] = {-1};

     for (i = 0; i < n_derivs; ++i) {
          chi_l  [0][0][0][i] = td->gc_l[2];
          omega_l   [0][0][i] = i;
          ltau_l    [0][0][i] = i;

          albedo_l  [0]   [i] = i;
     }

     HANDLE_RETURN(bounds_test_derivs(gd, td, n_solvers, solvers, 1, n_derivs, 1, LENGTH(theta_0), theta_0, LENGTH(theta), *theta, LENGTH(n_chi), *n_chi, *chi, LENGTH(omega), *omega, LENGTH(ltau), *ltau, LENGTH(albedo), albedo, LENGTH(phi), phi, LENGTH(chi_l[0]), ***chi_l, **omega_l, **ltau_l, *albedo_l, tol, tol_l, LENGTH(ignore_list_index_solver) - 1, ignore_list_index_solver, ignore_mask_index_solver, LENGTH(ignore_list_deriv) - 1, ignore_list_deriv, ignore_mask_deriv, ignore_mask_solver), "bounds_test_derivs()");

     free_array4((void ****) chi_l);
     free_array3_d(omega_l);
     free_array3_d(ltau_l);
     free_array2_d(albedo_l);

     return 0;
}



/*******************************************************************************
 *
 ******************************************************************************/
static int test_simple_vary_derivs_layers1(xrtm_data *gd, test_data *td, int n_solvers, enum xrtm_solver_mask *solvers, double *tol, double *tol_l, int n_derivs) {

     int i;
     int j;
     int n;

     double theta_0[]                = {35.};

     double theta  [][1]             = {{15.}};

     int n_chi     [][N_TEST_LAYERS] = {{td->n_gc[0], td->n_gc[2], td->n_gc[4]}};

     double **chi  [][N_TEST_LAYERS] = {{td->  gc[0], td->  gc[2], td->  gc[4]}};

     double omega  [][N_TEST_LAYERS] = {{.7,  .8, .9}};

     double ltau   [][N_TEST_LAYERS] = {{.01,  .1, 1.}};

     double albedo []                = {.3};

     double phi    []                = {45.};

     double ******chi_l = (double ******) alloc_array4(1, 1, N_TEST_LAYERS, n_derivs, sizeof(double **));

     double ***omega_l  = alloc_array3_d(1, N_TEST_LAYERS, n_derivs);

     double ***ltau_l   = alloc_array3_d(1, N_TEST_LAYERS, n_derivs);

     double **albedo_l  = alloc_array2_d(1, n_derivs);

     int ignore_list_index_solver[] = {-1};
     int ignore_mask_index_solver[] = {-1};

     int ignore_list_deriv       [] = {-1};
     int ignore_mask_deriv       [] = {-1};

     for (i = 0; i < N_TEST_LAYERS; ++i) {
          for (j = 0; j < n_derivs; ++j) {
               chi_l  [0][0][i][j] = NULL;
               omega_l   [0][i][j] = 0.;
               ltau_l    [0][i][j] = 0.;

               albedo_l  [0]   [j] = 0.;
          }
     }

     i = (N_TEST_LAYERS - 1) / 2 + 1;

     n = (n_derivs - 1) / 2 + 1;

     for (j = 0; j < n; ++j) {
          chi_l  [0][0][i][j] = td->gc_l[2];
          omega_l   [0][i][j] = i;
          ltau_l    [0][i][j] = i;

          albedo_l  [0]   [j] = 0.;
     }

     for ( ; j < n_derivs; ++j) {
          chi_l  [0][0][i][j] = NULL;
          omega_l   [0][i][j] = 0.;
          ltau_l    [0][i][j] = 0.;

          albedo_l  [0]   [j] = i;
     }

     HANDLE_RETURN(bounds_test_derivs(gd, td, n_solvers, solvers, N_TEST_LAYERS, n_derivs, 1, LENGTH(theta_0), theta_0, LENGTH(theta), *theta, LENGTH(n_chi), *n_chi, *chi, LENGTH(omega), *omega, LENGTH(ltau), *ltau, LENGTH(albedo), albedo, LENGTH(phi), phi, LENGTH(chi_l[0]), ***chi_l, **omega_l, **ltau_l, *albedo_l, tol, tol_l, LENGTH(ignore_list_index_solver) - 1, ignore_list_index_solver, ignore_mask_index_solver, LENGTH(ignore_list_deriv) - 1, ignore_list_deriv, ignore_mask_deriv, ignore_mask_solver), "bounds_test_derivs()");

     free_array4((void ****) chi_l);
     free_array3_d(omega_l);
     free_array3_d(ltau_l);
     free_array2_d(albedo_l);

     return 0;
}



static int test_simple_vary_derivs_layers2(xrtm_data *gd, test_data *td, int n_solvers, enum xrtm_solver_mask *solvers, double *tol, double *tol_l, int n_derivs) {

     int i;
     int j;

     double theta_0[]                = {35.};

     double theta  [][1]             = {{15.}};

     int n_chi     [][N_TEST_LAYERS] = {{td->n_gc[0], td->n_gc[2], td->n_gc[4]}};

     double **chi  [][N_TEST_LAYERS] = {{td->  gc[0], td->  gc[2], td->  gc[4]}};

     double omega  [][N_TEST_LAYERS] = {{.7,  .8, .9}};

     double ltau   [][N_TEST_LAYERS] = {{.01,  .1, 1.}};

     double albedo []                = {.3};

     double phi    []                = {45.};

     double ******chi_l = (double ******) alloc_array4(1, 1, N_TEST_LAYERS, n_derivs, sizeof(double **));

     double ***omega_l  = alloc_array3_d(1, N_TEST_LAYERS, n_derivs);

     double ***ltau_l   = alloc_array3_d(1, N_TEST_LAYERS, n_derivs);

     double **albedo_l  = alloc_array2_d(1, n_derivs);

     int ignore_list_index_solver[] = {-1};
     int ignore_mask_index_solver[] = {-1};

     int ignore_list_deriv       [] = {-1};
     int ignore_mask_deriv       [] = {-1};

     for (i = 0; i < N_TEST_LAYERS; ++i) {
          for (j = 0; j < n_derivs; ++j) {
               chi_l  [0][0][i][j] = NULL;
               omega_l   [0][i][j] = 0.;
               ltau_l    [0][i][j] = 0.;

               albedo_l  [0]   [j] = 0.;
          }
     }

     i = (N_TEST_LAYERS - 1) / 2 + 1;

     for (j = 0; j < n_derivs; ++j) {
          chi_l  [0][0][i][j] = td->gc_l[2];
          omega_l   [0][i][j] = i;
          ltau_l    [0][i][j] = i;

          albedo_l  [0]   [j] = i;
     }

     HANDLE_RETURN(bounds_test_derivs(gd, td, n_solvers, solvers, N_TEST_LAYERS, n_derivs, 1, LENGTH(theta_0), theta_0, LENGTH(theta), *theta, LENGTH(n_chi), *n_chi, *chi, LENGTH(omega), *omega, LENGTH(ltau), *ltau, LENGTH(albedo), albedo, LENGTH(phi), phi, LENGTH(chi_l[0]), ***chi_l, **omega_l, **ltau_l, *albedo_l, tol, tol_l, LENGTH(ignore_list_index_solver) - 1, ignore_list_index_solver, ignore_mask_index_solver, LENGTH(ignore_list_deriv) - 1, ignore_list_deriv, ignore_mask_deriv, ignore_mask_solver), "bounds_test_derivs()");

     free_array4((void ****) chi_l);
     free_array3_d(omega_l);
     free_array3_d(ltau_l);
     free_array2_d(albedo_l);

     return 0;
}



/*******************************************************************************
 *
 ******************************************************************************/
static int test_group_vary_derivs(test_data *td, test_group_vary_data *tcd, int options_mask, int solvers_mask1, int n_quad, int n_stokes, int n_derivs) {

     int options;

     int n_solvers;

     int solvers_mask2;

     enum xrtm_kernel_type kernels;

     enum xrtm_solver_mask solvers_array2[N_XRTM_SOLVERS];

     double tol  [N_XRTM_SOLVERS];
     double tol_l[N_XRTM_SOLVERS];


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     options = options_mask | XRTM_OPTION_CALC_DERIVS | XRTM_OPTION_DELTA_M;

     if ((n_solvers = make_solvers2(solvers_mask1, &solvers_mask2, solvers_array2, tol, tol_l,
                                    GROUP_SOLVERS_DEV_DERIVS1
                                    XRTM_SOLVER_DOUB_ADD, DEF_TOL, DEF_TOL,
                                    XRTM_SOLVER_EIG_ADD,  DEF_TOL, DEF_TOL,
                                    XRTM_SOLVER_EIG_BVP,  DEF_TOL, DEF_TOL,
                                    XRTM_SOLVER_MEM_BVP,  DEF_TOL, DEF_TOL,
                                    XRTM_SOLVER_PADE_ADD, 1.e-3, 1.e-2,
                                    XRTM_SOLVER_SOS,      DEF_TOL, DEF_TOL, 0)) < 0) {
          eprintf("ERROR: make_solvers2()\n");
          return -1;
     }

     kernels = XRTM_KERNEL_LAMBERTIAN;

     if (initiate_test_vary(td, n_solvers, solvers_array2, options, solvers_mask2, td->max_coef, n_quad, n_stokes, n_derivs, 1, 1, 32, &kernels, 1, tcd->func_vary_derivs1, "test_group_data->func_vary_derivs1()", tol, tol_l, n_derivs)) {
          eprintf("ERROR: initiate_test_vary()\n");
          return -1;
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     options = options_mask | XRTM_OPTION_CALC_DERIVS | XRTM_OPTION_DELTA_M;

     if ((n_solvers = make_solvers2(solvers_mask1, &solvers_mask2, solvers_array2, tol, tol_l,
                                    XRTM_SOLVER_DOUB_ADD, DEF_TOL, DEF_TOL,
                                    XRTM_SOLVER_EIG_ADD,  DEF_TOL, DEF_TOL,
                                    XRTM_SOLVER_EIG_BVP,  DEF_TOL, DEF_TOL,
                                    XRTM_SOLVER_MEM_BVP,  DEF_TOL, DEF_TOL,
                                    XRTM_SOLVER_PADE_ADD, 1.e-3, 1.e-2,
                                    XRTM_SOLVER_SOS,      DEF_TOL, DEF_TOL, 0)) < 0) {
          eprintf("ERROR: make_solvers2()\n");
          return -1;
     }

     kernels = XRTM_KERNEL_LAMBERTIAN;

     if (initiate_test_vary(td, n_solvers, solvers_array2, options, solvers_mask2, td->max_coef, n_quad, n_stokes, n_derivs, 1, 1, 32, &kernels, 1, tcd->func_vary_derivs2, "test_group_data->func_vary_derivs2()", tol, tol_l, n_derivs)) {
          eprintf("ERROR: initiate_test_vary()\n");
          return -1;
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     options = options_mask | XRTM_OPTION_CALC_DERIVS | XRTM_OPTION_DELTA_M;

     if ((n_solvers = make_solvers2(solvers_mask1, &solvers_mask2, solvers_array2, tol, tol_l,
                                    GROUP_SOLVERS_DEV_DERIVS1
                                    XRTM_SOLVER_DOUB_ADD, DEF_TOL, DEF_TOL,
                                    XRTM_SOLVER_EIG_ADD,  DEF_TOL, DEF_TOL,
                                    XRTM_SOLVER_EIG_BVP,  DEF_TOL, DEF_TOL,
                                    XRTM_SOLVER_MEM_BVP,  DEF_TOL, DEF_TOL,
                                    XRTM_SOLVER_PADE_ADD, 1.e-3, 1.e-2,
                                    XRTM_SOLVER_SOS,      DEF_TOL, DEF_TOL, 0)) < 0) {
          eprintf("ERROR: make_solvers2()\n");
          return -1;
     }

     kernels = XRTM_KERNEL_LAMBERTIAN;

     if (initiate_test_vary_layers(td, n_solvers, solvers_array2, options, solvers_mask2, td->max_coef, n_quad, n_stokes, n_derivs, tcd->n_layers, 1, 32, &kernels, 1, tcd->func_vary_derivs_layers1, "test_group_data->func_vary_derivs_layers1()", tol, tol_l, n_derivs)) {
          eprintf("ERROR: initiate_test_vary_layers(()\n");
          return -1;
     }



     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     options = options_mask | XRTM_OPTION_CALC_DERIVS | XRTM_OPTION_DELTA_M;

     if ((n_solvers = make_solvers2(solvers_mask1, &solvers_mask2, solvers_array2, tol, tol_l,
                                    XRTM_SOLVER_DOUB_ADD, DEF_TOL, DEF_TOL,
                                    XRTM_SOLVER_EIG_ADD,  DEF_TOL, DEF_TOL,
                                    XRTM_SOLVER_EIG_BVP,  DEF_TOL, DEF_TOL,
                                    XRTM_SOLVER_MEM_BVP,  DEF_TOL, DEF_TOL,
                                    XRTM_SOLVER_PADE_ADD, 1.e-3, 1.e-2,
                                    XRTM_SOLVER_SOS,      DEF_TOL, DEF_TOL, 0)) < 0) {
          eprintf("ERROR: make_solvers2()\n");
          return -1;
     }

     kernels = XRTM_KERNEL_LAMBERTIAN;

     if (initiate_test_vary_layers(td, n_solvers, solvers_array2, options, solvers_mask2, td->max_coef, n_quad, n_stokes, n_derivs, tcd->n_layers, 1, 32, &kernels, 1, tcd->func_vary_derivs_layers2, "test_group_data->func_vary_derivs_layers2()", tol, tol_l, n_derivs)) {
          eprintf("ERROR: initiate_test_vary_layers()\n");
          return -1;
     }


     return 0;
}



/*******************************************************************************
 *
 ******************************************************************************/
static int test_group_vary_derivs2(test_data *td, test_group_vary_data *tcd, int option_mask, int solvers_mask1, int n_quad, int n_quad_stokes, int n_derivs) {

     int solvers_mask2;

     if (test_group_vary_derivs(td, tcd, option_mask, solvers_mask1, n_quad, 1, n_derivs)) {
          eprintf("ERROR: test_group_vary_derivs()\n");
          return -1;
     }

     solvers_mask2 = solvers_mask1 & XRTM_SOLVERS_VECTOR;

     if (test_group_vary_derivs(td, tcd, option_mask | XRTM_OPTION_PHASE_MATRIX_GC | XRTM_OPTION_VECTOR, solvers_mask2, n_quad_stokes, 4, n_derivs)) {
          eprintf("ERROR: test_group_vary_derivs()\n");
          return -1;
     }

     return 0;
}



/*******************************************************************************
 *
 ******************************************************************************/
static int test_simple_vary_layers(xrtm_data *gd, test_data *td, int n_solvers, enum xrtm_solver_mask *solvers, double *tol, double *tol_l, int n_layers) {

     int i;

     double theta_0[ ]    = {35.};

     double theta  [ ][1] = {{15.}};

     int **n_chi    = alloc_array2_i(1, n_layers);

     double ****chi = (double ****) alloc_array2(1, n_layers, sizeof(double **));

     double **omega = alloc_array2_d(1, n_layers);

     double **ltau  = alloc_array2_d(1, n_layers);

     double albedo [ ] = {.3};

     double phi    [ ] = {45.};

     int ignore_list_index_solver[] = {-1};
     int ignore_mask_index_solver[] = {-1};

     for (i = 0; i < n_layers; ++i) {
          n_chi[0][i] = td->n_gc[2];
          chi  [0][i] = td->  gc[2];
          omega[0][i] = .7;
          ltau [0][i] = .1;
     }

     HANDLE_RETURN(bounds_test(gd, td, n_solvers, solvers, n_layers, 1, LENGTH(theta_0), theta_0, LENGTH(theta), *theta, LENGTH(n_chi), *n_chi, *chi, LENGTH(omega), *omega, LENGTH(ltau), *ltau, LENGTH(albedo), albedo, LENGTH(phi), phi, tol, tol_l, LENGTH(ignore_list_index_solver) - 1, ignore_list_index_solver, ignore_mask_index_solver, ignore_mask_solver), "bounds_test()");

     free_array2_i(n_chi);
     free_array2((void **) chi);
     free_array2_d(omega);
     free_array2_d(ltau);

     return 0;
}



/*******************************************************************************
 *
 ******************************************************************************/
static int test_simple_vary_layers_derivs1(xrtm_data *gd, test_data *td, int n_solvers, enum xrtm_solver_mask *solvers, double *tol, double *tol_l, int n_layers) {

     int i;
     int j;

     double theta_0[ ]    = {35.};

     double theta  [ ][1] = {{15.}};

     int **n_chi    = alloc_array2_i(1, n_layers);

     double ****chi = (double ****) alloc_array2(1, n_layers, sizeof(double **));

     double **omega = alloc_array2_d(1, n_layers);

     double **ltau  = alloc_array2_d(1, n_layers);

     double albedo [ ] = {.3};

     double phi    [ ] = {45.};

     double ******chi_l = (double ******) alloc_array4(1, 1, n_layers, N_TEST_SIMPLE_LAYERS_DERIVS1, sizeof(double **));

     double ***omega_l  = alloc_array3_d(1, n_layers, N_TEST_SIMPLE_LAYERS_DERIVS1);

     double ***ltau_l   = alloc_array3_d(1, n_layers, N_TEST_SIMPLE_LAYERS_DERIVS1);

     double **albedo_l  = alloc_array2_d(1, N_TEST_SIMPLE_LAYERS_DERIVS1);

     int ignore_list_index_solver[] = {-1};
     int ignore_mask_index_solver[] = {-1};

     int ignore_list_deriv       [] = {-1};
     int ignore_mask_deriv       [] = {-1};

     for (i = 0; i < n_layers; ++i) {
          n_chi[0][i] = td->n_gc[2];
          chi  [0][i] = td->  gc[2];
          omega[0][i] = .7;
          ltau [0][i] = .1;

          for (j = 0; j < N_TEST_SIMPLE_LAYERS_DERIVS1; ++j) {
               chi_l  [0][0][i][j] = NULL;
               omega_l   [0][i][j] = 0.;
               ltau_l    [0][i][j] = 0.;
          }
     }

     for (j = 0; j < N_TEST_SIMPLE_LAYERS_DERIVS1; ++j)
          albedo_l[0][j] = 0.;

     i = n_layers / 2;

     chi_l  [0][0][i][0] = td->gc_l[2];
     omega_l   [0][i][0] = 1.;
     ltau_l    [0][i][0] = 1.;

     albedo_l[0][1] = 1.;

     HANDLE_RETURN(bounds_test_derivs(gd, td, n_solvers, solvers, n_layers, N_TEST_SIMPLE_LAYERS_DERIVS1, 1, LENGTH(theta_0), theta_0, LENGTH(theta), *theta, LENGTH(n_chi), *n_chi, *chi, LENGTH(omega), *omega, LENGTH(ltau), *ltau, LENGTH(albedo), albedo, LENGTH(phi), phi, LENGTH(chi_l[0]), ***chi_l, **omega_l, **ltau_l, *albedo_l, tol, tol_l, LENGTH(ignore_list_index_solver) - 1, ignore_list_index_solver, ignore_mask_index_solver, LENGTH(ignore_list_deriv) - 1, ignore_list_deriv, ignore_mask_deriv, ignore_mask_solver), "bounds_test_derivs()");

     free_array2_i(n_chi);
     free_array2((void **) chi);
     free_array2_d(omega);
     free_array2_d(ltau);

     free_array4((void ****) chi_l);
     free_array3_d(omega_l);
     free_array3_d(ltau_l);
     free_array2_d(albedo_l);

     return 0;
}



static int test_simple_vary_layers_derivs2(xrtm_data *gd, test_data *td, int n_solvers, enum xrtm_solver_mask *solvers, double *tol, double *tol_l, int n_layers) {

     int i;
     int j;

     double theta_0[ ]    = {35.};

     double theta  [ ][1] = {{15.}};

     int **n_chi    = alloc_array2_i(1, n_layers);

     double ****chi = (double ****) alloc_array2(1, n_layers, sizeof(double **));

     double **omega = alloc_array2_d(1, n_layers);

     double **ltau  = alloc_array2_d(1, n_layers);

     double albedo [ ] = {.3};

     double phi    [ ] = {45.};

     double ******chi_l = (double ******) alloc_array4(1, 1, n_layers, N_TEST_SIMPLE_LAYERS_DERIVS2, sizeof(double **));

     double ***omega_l  = alloc_array3_d(1, n_layers, N_TEST_SIMPLE_LAYERS_DERIVS2);

     double ***ltau_l   = alloc_array3_d(1, n_layers, N_TEST_SIMPLE_LAYERS_DERIVS2);

     double **albedo_l  = alloc_array2_d(1, N_TEST_SIMPLE_LAYERS_DERIVS2);

     int ignore_list_index_solver[] = {-1};
     int ignore_mask_index_solver[] = {-1};

     int ignore_list_deriv       [] = {-1};
     int ignore_mask_deriv       [] = {-1};

     for (i = 0; i < n_layers; ++i) {
          n_chi[0][i] = td->n_gc[2];
          chi  [0][i] = td->  gc[2];
          omega[0][i] = .7;
          ltau [0][i] = .1;

          for (j = 0; j < N_TEST_SIMPLE_LAYERS_DERIVS2; ++j) {
               chi_l  [0][0][i][j] = NULL;
               omega_l   [0][i][j] = 0.;
               ltau_l    [0][i][j] = 0.;
          }
     }

     for (j = 0; j < N_TEST_SIMPLE_LAYERS_DERIVS2; ++j)
          albedo_l[0][j] = 0.;

     i = n_layers / 2;

     chi_l  [0][0][i][0] = td->gc_l[2];
     omega_l   [0][i][0] = 1.;
     ltau_l    [0][i][0] = 1.;

     albedo_l[0][0] = 1.;

     HANDLE_RETURN(bounds_test_derivs(gd, td, n_solvers, solvers, n_layers, N_TEST_SIMPLE_LAYERS_DERIVS2, 1, LENGTH(theta_0), theta_0, LENGTH(theta), *theta, LENGTH(n_chi), *n_chi, *chi, LENGTH(omega), *omega, LENGTH(ltau), *ltau, LENGTH(albedo), albedo, LENGTH(phi), phi, LENGTH(chi_l[0]), ***chi_l, **omega_l, **ltau_l, *albedo_l, tol, tol_l, LENGTH(ignore_list_index_solver) - 1, ignore_list_index_solver, ignore_mask_index_solver, LENGTH(ignore_list_deriv) - 1, ignore_list_deriv, ignore_mask_deriv, ignore_mask_solver), "bounds_test_derivs()");

     free_array2_i(n_chi);
     free_array2((void **) chi);
     free_array2_d(omega);
     free_array2_d(ltau);

     free_array4((void ****) chi_l);
     free_array3_d(omega_l);
     free_array3_d(ltau_l);
     free_array2_d(albedo_l);

     return 0;
}



/*******************************************************************************
 *
 ******************************************************************************/
static int test_group_vary_layers(test_data *td, test_group_vary_data *tcd, int options_mask, int solvers_mask1, int n_quad, int n_stokes, int n_layers) {

     int options;

     int n_solvers;

     int solvers_mask2;

     enum xrtm_kernel_type kernels;

     enum xrtm_solver_mask solvers_array2[N_XRTM_SOLVERS];

     double tol  [N_XRTM_SOLVERS];
     double tol_l[N_XRTM_SOLVERS];


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     options = options_mask | XRTM_OPTION_DELTA_M;

     if ((n_solvers = make_solvers(solvers_mask1, &solvers_mask2, solvers_array2, tol,
                                    GROUP_SOLVERS_DEV
                                    XRTM_SOLVER_DOUB_ADD, DEF_TOL,
                                    XRTM_SOLVER_EIG_ADD,  DEF_TOL,
                                    XRTM_SOLVER_EIG_BVP,  DEF_TOL,
                                    XRTM_SOLVER_MEM_BVP,  DEF_TOL,
                                    XRTM_SOLVER_PADE_ADD, 1.e-3,
                                    XRTM_SOLVER_SOS,      DEF_TOL, 0)) < 0) {
          eprintf("ERROR: make_solvers()\n");
          return -1;
     }

     kernels = XRTM_KERNEL_LAMBERTIAN;

     if (initiate_test_vary(td, n_solvers, solvers_array2, options, solvers_mask2, td->max_coef, n_quad, n_stokes, 0, n_layers, 1, 32, &kernels, 1, tcd->func_vary_layers, "test_group_data->func_vary_layers()", tol, NULL, n_layers)) {
          eprintf("ERROR: initiate_test_vary()\n");
          return -1;
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     options = options_mask | XRTM_OPTION_CALC_DERIVS | XRTM_OPTION_DELTA_M;

     if ((n_solvers = make_solvers2(solvers_mask1, &solvers_mask2, solvers_array2, tol, tol_l,
                                    GROUP_SOLVERS_DEV_DERIVS1
                                    XRTM_SOLVER_DOUB_ADD, DEF_TOL, DEF_TOL,
                                    XRTM_SOLVER_EIG_ADD,  DEF_TOL, DEF_TOL,
                                    XRTM_SOLVER_EIG_BVP,  DEF_TOL, DEF_TOL,
                                    XRTM_SOLVER_MEM_BVP,  DEF_TOL, DEF_TOL,
                                    XRTM_SOLVER_PADE_ADD, 1.e-3, 1.e-2,
                                    XRTM_SOLVER_SOS,      DEF_TOL, DEF_TOL, 0)) < 0) {
          eprintf("ERROR: make_solvers2()\n");
          return -1;
     }

     kernels = XRTM_KERNEL_LAMBERTIAN;

     if (initiate_test_vary(td, n_solvers, solvers_array2, options, solvers_mask2, td->max_coef, n_quad, n_stokes, tcd->n_derivs1, n_layers, 1, 32, &kernels, 1, tcd->func_vary_layers_derivs1, "test_group_data->func_vary_layers_derivs1()", tol, tol_l, n_layers)) {
          eprintf("ERROR: initiate_test_vary()\n");
          return -1;
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     options = options_mask | XRTM_OPTION_CALC_DERIVS | XRTM_OPTION_DELTA_M;

     if ((n_solvers = make_solvers2(solvers_mask1, &solvers_mask2, solvers_array2, tol, tol_l,
                                    XRTM_SOLVER_DOUB_ADD, DEF_TOL, DEF_TOL,
                                    XRTM_SOLVER_EIG_ADD,  DEF_TOL, DEF_TOL,
                                    XRTM_SOLVER_EIG_BVP,  DEF_TOL, DEF_TOL,
                                    XRTM_SOLVER_MEM_BVP,  DEF_TOL, DEF_TOL,
                                    XRTM_SOLVER_PADE_ADD, 1.e-3, 1.e-2,
                                    XRTM_SOLVER_SOS,      DEF_TOL, DEF_TOL, 0)) < 0) {
          eprintf("ERROR: make_solvers2()\n");
          return -1;
     }

     kernels = XRTM_KERNEL_LAMBERTIAN;

     if (initiate_test_vary(td, n_solvers, solvers_array2, options, solvers_mask2, td->max_coef, n_quad, n_stokes, tcd->n_derivs2, n_layers, 1, 32, &kernels, 1, tcd->func_vary_layers_derivs2, "test_group_data->func_vary_layers_derivs2()", tol, tol_l, n_layers)) {
          eprintf("ERROR: initiate_test_vary()\n");
          return -1;
     }


     return 0;
}



/*******************************************************************************
 *
 ******************************************************************************/
static int test_group_vary_layers2(test_data *td, test_group_vary_data *tcd, int option_mask, int solvers_mask1, int n_quad, int n_quad_stokes, int n_layers) {

     int solvers_mask2;

     if (test_group_vary_layers(td, tcd, option_mask, solvers_mask1, n_quad, 1, n_layers)) {
          eprintf("ERROR: test_group_vary_layers()\n");
          return -1;
     }

     solvers_mask2 = solvers_mask1 & XRTM_SOLVERS_VECTOR;

     if (test_group_vary_layers(td, tcd, option_mask | XRTM_OPTION_PHASE_MATRIX_GC | XRTM_OPTION_VECTOR, solvers_mask2, n_quad_stokes, 4, n_layers)) {
          eprintf("ERROR: test_group_vary_layers()\n");
          return -1;
     }

     return 0;
}



/*******************************************************************************
 *
 ******************************************************************************/
void init_thetas(int n_thetas, double *thetas) {

     int i;

     double a;
     double d_theta;

     d_theta = (0. - 90.) / n_thetas;

     a = 90.;
     for (i = 0; i < n_thetas; ++i) {
          thetas[i] = a + d_theta / 2.;
          a += d_theta;
     }
}



/*******************************************************************************
 *
 ******************************************************************************/
static int test_simple_vary_thetas(xrtm_data *gd, test_data *td, int n_solvers, enum xrtm_solver_mask *solvers, double *tol, double *tol_l, int n_thetas) {

     double theta_0[ ]    = {35.};

     double **theta = alloc_array2_d(1, n_thetas);

     int n_chi     [ ][1] = {{td->n_gc[3]}};
     double **chi  [ ][1] = {{td->  gc[3]}};

     double omega  [ ][1] = {{.7}};

     double ltau   [ ][1] = {{1.}};

     double albedo [ ]    = {.3};

     double phi    [ ]    = {45.};

     int ignore_list_index_solver[] = {-1};
     int ignore_mask_index_solver[] = {-1};

     init_thetas(n_thetas, theta[0]);
     
     HANDLE_RETURN(bounds_test(gd, td, n_solvers, solvers, 1, n_thetas, LENGTH(theta_0), theta_0, LENGTH(theta), *theta, LENGTH(n_chi), *n_chi, *chi, LENGTH(omega), *omega, LENGTH(ltau), *ltau, LENGTH(albedo), albedo, LENGTH(phi), phi, tol, tol_l, LENGTH(ignore_list_index_solver) - 1, ignore_list_index_solver, ignore_mask_index_solver, ignore_mask_solver), "bounds_test_derivs()");

     free_array2_d(theta);

     return 0;
}



/*******************************************************************************
 *
 ******************************************************************************/
static int test_simple_vary_thetas_derivs1(xrtm_data *gd, test_data *td, int n_solvers, enum xrtm_solver_mask *solvers, double *tol, double *tol_l, int n_thetas) {

     double theta_0[ ]    = {35.};

     double **theta = alloc_array2_d(1, n_thetas);

     int n_chi     [ ][1] = {{td->n_gc[2]}};
     double **chi  [ ][1] = {{td->  gc[2]}};

     double omega  [ ][1] = {{.7}};

     double ltau   [ ][1] = {{1.}};

     double albedo [ ]    = {.3};

     double phi    [ ]    = {45.};

     double **chi_l[][1][1][N_TEST_SIMPLE_DERIVS1] = {{{{td->gc_l[2], NULL}}}};

     double omega_l  [ ][1][N_TEST_SIMPLE_DERIVS1] = {{{1., 0.}}};

     double ltau_l   [ ][1][N_TEST_SIMPLE_DERIVS1] = {{{1., 0.}}};

     double albedo_l [ ]   [N_TEST_SIMPLE_DERIVS1] = { {0., 2.} };

     int ignore_list_index_solver[] = {-1};
     int ignore_mask_index_solver[] = {-1};

     int ignore_list_deriv       [] = {-1};
     int ignore_mask_deriv       [] = {-1};

     init_thetas(n_thetas, theta[0]);

     HANDLE_RETURN(bounds_test_derivs(gd, td, n_solvers, solvers, 1, N_TEST_SIMPLE_DERIVS1, n_thetas, LENGTH(theta_0), theta_0, LENGTH(theta), *theta, LENGTH(n_chi), *n_chi, *chi, LENGTH(omega), *omega, LENGTH(ltau), *ltau, LENGTH(albedo), albedo, LENGTH(phi), phi, LENGTH(chi_l[0]), ***chi_l, **omega_l, **ltau_l, *albedo_l, tol, tol_l, LENGTH(ignore_list_index_solver) - 1, ignore_list_index_solver, ignore_mask_index_solver, LENGTH(ignore_list_deriv) - 1, ignore_list_deriv, ignore_mask_deriv, ignore_mask_solver), "bounds_test_derivs()");

     free_array2_d(theta);

     return 0;
}



/*******************************************************************************
 *
 ******************************************************************************/
static int test_simple_vary_thetas_layers(xrtm_data *gd, test_data *td, int n_solvers, enum xrtm_solver_mask *solvers, double *tol, double *tol_l, int n_thetas) {

     double theta_0[ ]                = {35.};

     double **theta = alloc_array2_d(1, n_thetas);

     int n_chi     [ ][N_TEST_LAYERS] = {{td->n_gc[0], td->n_gc[2], td->n_gc[4]}};

     double **chi  [ ][N_TEST_LAYERS] = {{td->  gc[0], td->  gc[2], td->  gc[4]}};

     double omega  [ ][N_TEST_LAYERS] = {{.7,  .8, .9}};

     double ltau   [ ][N_TEST_LAYERS] = {{.01,  .1, 1.}};

     double albedo [ ]                = {.3};

     double phi    [ ]                = {45.};

     int ignore_list_index_solver[] = {-1};
     int ignore_mask_index_solver[] = {-1};

     init_thetas(n_thetas, theta[0]);

     HANDLE_RETURN(bounds_test(gd, td, n_solvers, solvers, N_TEST_LAYERS, n_thetas, LENGTH(theta_0), theta_0, LENGTH(theta), *theta, LENGTH(n_chi), *n_chi, *chi, LENGTH(omega), *omega, LENGTH(ltau), *ltau, LENGTH(albedo), albedo, LENGTH(phi), phi, tol, tol_l, LENGTH(ignore_list_index_solver) - 1, ignore_list_index_solver, ignore_mask_index_solver, ignore_mask_solver), "bounds_test()");

     free_array2_d(theta);

     return 0;
}



/*******************************************************************************
 *
 ******************************************************************************/
static int test_simple_vary_thetas_layers_derivs1(xrtm_data *gd, test_data *td, int n_solvers, enum xrtm_solver_mask *solvers, double *tol, double *tol_l, int n_thetas) {

     double theta_0[ ]                = {35.};

     double **theta = alloc_array2_d(1, n_thetas);

     int n_chi     [ ][N_TEST_LAYERS] = {{td->n_gc[0], td->n_gc[2], td->n_gc[4]}};

     double **chi  [ ][N_TEST_LAYERS] = {{td->  gc[0], td->  gc[2], td->  gc[4]}};

     double omega  [ ][N_TEST_LAYERS] = {{.7,  .8, .9}};

     double ltau   [ ][N_TEST_LAYERS] = {{.01,  .1, 1.}};

     double albedo [ ]                = {.3};

     double phi    [ ]                = {45.};

     double **chi_l[1][1][N_TEST_LAYERS][N_TEST_SIMPLE_LAYERS_DERIVS1] = {{{{NULL,        NULL},
                                                                            {td->gc_l[2], NULL},
                                                                            {NULL,        NULL}}}};

     double omega_l   [ ][N_TEST_LAYERS][N_TEST_SIMPLE_LAYERS_DERIVS1] = {{{0., 0.},
                                                                           {1., 0.},
                                                                           {0., 0.}}};

     double ltau_l    [ ][N_TEST_LAYERS][N_TEST_SIMPLE_LAYERS_DERIVS1] = {{{0., 0.},
                                                                           {1., 0.},
                                                                           {0., 0.}}};

     double albedo_l  [ ]               [N_TEST_SIMPLE_LAYERS_DERIVS1] =  {{0., 1.}};

     int ignore_list_index_solver[] = {-1};
     int ignore_mask_index_solver[] = {-1};

     int ignore_list_deriv       [] = {-1};
     int ignore_mask_deriv       [] = {-1};

     init_thetas(n_thetas, theta[0]);

     HANDLE_RETURN(bounds_test_derivs(gd, td, n_solvers, solvers, N_TEST_LAYERS, N_TEST_SIMPLE_LAYERS_DERIVS1, n_thetas, LENGTH(theta_0), theta_0, LENGTH(theta), *theta, LENGTH(n_chi), *n_chi, *chi, LENGTH(omega), *omega, LENGTH(ltau), *ltau, LENGTH(albedo), albedo, LENGTH(phi), phi, LENGTH(chi_l[0]), ***chi_l, **omega_l, **ltau_l, *albedo_l, tol, tol_l, LENGTH(ignore_list_index_solver) - 1, ignore_list_index_solver, ignore_mask_index_solver, LENGTH(ignore_list_deriv) - 1, ignore_list_deriv, ignore_mask_deriv, ignore_mask_solver), "bounds_test_derivs()");

     free_array2_d(theta);

     return 0;
}



/*******************************************************************************
 *
 ******************************************************************************/
static int test_group_vary_thetas(test_data *td, test_group_vary_data *tcd, int options_mask, int solvers_mask1, int n_quad, int n_stokes, int n_thetas) {

     int options;

     int n_solvers;

     int solvers_mask2;

     enum xrtm_kernel_type kernels;

     enum xrtm_solver_mask solvers_array2[N_XRTM_SOLVERS];

     double tol  [N_XRTM_SOLVERS];
     double tol_l[N_XRTM_SOLVERS];


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     options = options_mask | XRTM_OPTION_DELTA_M;

     if ((n_solvers = make_solvers(solvers_mask1, &solvers_mask2, solvers_array2, tol,
                                    GROUP_SOLVERS_DEV
                                    XRTM_SOLVER_DOUB_ADD, DEF_TOL,
                                    XRTM_SOLVER_EIG_ADD,  DEF_TOL,
                                    XRTM_SOLVER_EIG_BVP,  DEF_TOL,
                                    XRTM_SOLVER_MEM_BVP,  DEF_TOL,
                                    XRTM_SOLVER_PADE_ADD, 1.e-3,
                                    XRTM_SOLVER_SOS,      DEF_TOL, 0)) < 0) {
          eprintf("ERROR: make_solvers()\n");
          return -1;
     }

     kernels = XRTM_KERNEL_LAMBERTIAN;

     if (initiate_test_vary(td, n_solvers, solvers_array2, options, solvers_mask2, td->max_coef, n_quad, n_stokes, 0, 1, 1, 32, &kernels, n_thetas, tcd->func_vary_thetas, "test_group_data->func_vary_thetas()", tol, NULL, n_thetas)) {
          eprintf("ERROR: initiate_test_vary()\n");
          return -1;
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     options = options_mask | XRTM_OPTION_CALC_DERIVS | XRTM_OPTION_DELTA_M;

     if ((n_solvers = make_solvers2(solvers_mask1, &solvers_mask2, solvers_array2, tol, tol_l,
                                    GROUP_SOLVERS_DEV_DERIVS1
                                    XRTM_SOLVER_DOUB_ADD, DEF_TOL, DEF_TOL,
                                    XRTM_SOLVER_EIG_ADD,  DEF_TOL, DEF_TOL,
                                    XRTM_SOLVER_EIG_BVP,  DEF_TOL, DEF_TOL,
                                    XRTM_SOLVER_MEM_BVP,  DEF_TOL, DEF_TOL,
                                    XRTM_SOLVER_PADE_ADD, 1.e-3, 1.e-2,
                                    XRTM_SOLVER_SOS,      DEF_TOL, DEF_TOL, 0)) < 0) {
          eprintf("ERROR: make_solvers2()\n");
          return -1;
     }

     kernels = XRTM_KERNEL_LAMBERTIAN;

     if (initiate_test_vary(td, n_solvers, solvers_array2, options, solvers_mask2, td->max_coef, n_quad, n_stokes, tcd->n_derivs1, 1, 1, 32, &kernels, n_thetas, tcd->func_vary_thetas_derivs1, "test_group_data->func_vary_thetas_derivs1()", tol, tol_l, n_thetas)) {
          eprintf("ERROR: initiate_test_vary()\n");
          return -1;
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     options = options_mask | XRTM_OPTION_DELTA_M;

     if ((n_solvers = make_solvers(solvers_mask1, &solvers_mask2, solvers_array2, tol,
                                    GROUP_SOLVERS_DEV
                                    XRTM_SOLVER_DOUB_ADD, DEF_TOL,
                                    XRTM_SOLVER_EIG_ADD,  DEF_TOL,
                                    XRTM_SOLVER_EIG_BVP,  DEF_TOL,
                                    XRTM_SOLVER_MEM_BVP,  DEF_TOL,
                                    XRTM_SOLVER_PADE_ADD, 1.e-3,
                                    XRTM_SOLVER_SOS,      DEF_TOL, 0)) < 0) {
          eprintf("ERROR: make_solvers()\n");
          return -1;
     }

     kernels = XRTM_KERNEL_LAMBERTIAN;

     if (initiate_test_vary(td, n_solvers, solvers_array2, options, solvers_mask2, td->max_coef, n_quad, n_stokes, 0, tcd->n_layers, 1, 32, &kernels, n_thetas, tcd->func_vary_thetas_layers, "test_group_data->func_vary_thetas_layers()", tol, NULL, n_thetas)) {
          eprintf("ERROR: initiate_test_vary()\n");
          return -1;
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     options = options_mask | XRTM_OPTION_CALC_DERIVS | XRTM_OPTION_DELTA_M;

     if ((n_solvers = make_solvers2(solvers_mask1, &solvers_mask2, solvers_array2, tol, tol_l,
                                    GROUP_SOLVERS_DEV_DERIVS12
                                    XRTM_SOLVER_DOUB_ADD, DEF_TOL, DEF_TOL,
                                    XRTM_SOLVER_EIG_ADD,  DEF_TOL, DEF_TOL,
                                    XRTM_SOLVER_EIG_BVP,  DEF_TOL, DEF_TOL,
                                    XRTM_SOLVER_MEM_BVP,  DEF_TOL, DEF_TOL,
                                    XRTM_SOLVER_PADE_ADD, 1.e-3, 1.e-2,
                                    XRTM_SOLVER_SOS,      DEF_TOL, DEF_TOL, 0)) < 0) {
          eprintf("ERROR: make_solvers2()\n");
          return -1;
     }

     kernels = XRTM_KERNEL_LAMBERTIAN;

     if (initiate_test_vary(td, n_solvers, solvers_array2, options, solvers_mask2, td->max_coef, n_quad, n_stokes, tcd->n_derivs1, tcd->n_layers, 1, 32, &kernels, n_thetas, tcd->func_vary_thetas_layers_derivs1, "test_group_data->func_vary_thetas_layers_derivs1()", tol, tol_l, n_thetas)) {
          eprintf("ERROR: initiate_test_vary()\n");
          return -1;
     }


     return 0;
}



/*******************************************************************************
 *
 ******************************************************************************/
static int test_group_vary_thetas2(test_data *td, test_group_vary_data *tcd, int option_mask, int solvers_mask1, int n_quad, int n_quad_stokes, int n_thetas) {

     int solvers_mask2;

     if (test_group_vary_thetas(td, tcd, option_mask, solvers_mask1, n_quad, 1, n_thetas)) {
          eprintf("ERROR: test_group_vary_thetas()\n");
          return -1;
     }

     solvers_mask2 = solvers_mask1 & XRTM_SOLVERS_VECTOR;

     if (test_group_vary_thetas(td, tcd, option_mask | XRTM_OPTION_PHASE_MATRIX_GC | XRTM_OPTION_VECTOR, solvers_mask2, n_quad_stokes, 4, n_thetas)) {
          eprintf("ERROR: test_group_vary_thetas()\n");
          return -1;
     }

     return 0;
}



/*******************************************************************************
 *
 ******************************************************************************/
static int test_group_vary_init_simple(test_group_vary_data *d) {

     d->n_derivs1                       = N_TEST_SIMPLE_DERIVS1;
     d->n_derivs2                       = N_TEST_SIMPLE_DERIVS2;
     d->n_layers                        = N_TEST_LAYERS;
     d->n_layers_derivs1                = N_TEST_SIMPLE_LAYERS_DERIVS1;
     d->n_layers_derivs2                = N_TEST_SIMPLE_LAYERS_DERIVS2;

     d->func_vary_layers                = test_simple_vary_layers;
     d->func_vary_layers_derivs1        = test_simple_vary_layers_derivs1;
     d->func_vary_layers_derivs2        = test_simple_vary_layers_derivs2;

     d->func_vary_derivs1               = test_simple_vary_derivs1;
     d->func_vary_derivs2               = test_simple_vary_derivs2;
     d->func_vary_derivs_layers1        = test_simple_vary_derivs_layers1;
     d->func_vary_derivs_layers2        = test_simple_vary_derivs_layers2;

     d->func_vary_thetas                = test_simple_vary_thetas;
     d->func_vary_thetas_derivs1        = test_simple_vary_thetas_derivs1;
     d->func_vary_thetas_layers         = test_simple_vary_thetas_layers;
     d->func_vary_thetas_layers_derivs1 = test_simple_vary_thetas_layers_derivs1;

     return 0;
}



static int test_group_vary_init_bounds(test_group_vary_data *d) {

     d->n_derivs1                = N_TEST_BOUNDS_DERIVS1;
     d->n_derivs2                = N_TEST_BOUNDS_DERIVS2;
     d->n_layers                 = N_TEST_LAYERS;
     d->n_layers_derivs1         = N_TEST_BOUNDS_LAYERS_DERIVS1;
     d->n_layers_derivs2         = N_TEST_BOUNDS_LAYERS_DERIVS2;

     return 0;
}



/*******************************************************************************
 *
 ******************************************************************************/
int test_core(test_data *td) {

     int i;
     int n;

     int mask;

     int max_pow;

     int solvers_mask;
     int solvers_mask2;

     test_group_data test_group_simple;
     test_group_data test_group_bounds;

     test_group_vary_data test_group_vary_simple;
     test_group_vary_data test_group_vary_bounds;


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     test_group_init_simple(&test_group_simple);
     test_group_init_bounds(&test_group_bounds);

     test_group_vary_init_simple(&test_group_vary_simple);
     test_group_vary_init_bounds(&test_group_vary_bounds);


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     printf("Range over n_quad: index = %d\n", td->index);

     solvers_mask = XRTM_SOLVER_DOUB_ADD | XRTM_SOLVER_EIG_ADD | XRTM_SOLVER_EIG_BVP | XRTM_SOLVER_MEM_BVP | XRTM_SOLVER_PADE_ADD | XRTM_SOLVER_SINGLE | XRTM_SOLVER_SOS | XRTM_SOLVER_TWO_OS | DEV_RANGE_OVER_N_QUAD;

     max_pow = MAX_N_QUAD_POWER_SCALER;
     if (td->quick_run)
          max_pow = 1;

     for (i = 0; i <= max_pow; ++i) {
          n = 1 << i;

          if (test_group(td, &test_group_simple, 0, solvers_mask, n, 1)) {
               eprintf("ERROR: test_group()\n");
               return -1;
          }
     }

     max_pow = MAX_N_QUAD_POWER_VECTOR;
     if (td->quick_run)
          max_pow = 1;

     for (i = 0; i <= max_pow; ++i) {
          n = 1 << i;

          if (test_group(td, &test_group_simple, XRTM_OPTION_PHASE_MATRIX_GC | XRTM_OPTION_VECTOR, solvers_mask, n, 4)) {
               eprintf("ERROR: test_group()\n");
               return -1;
          }
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     printf("Range over n_stokes: index = %d\n", td->index);

     solvers_mask = XRTM_SOLVER_DOUB_ADD | XRTM_SOLVER_EIG_ADD | XRTM_SOLVER_EIG_BVP | XRTM_SOLVER_MEM_BVP | XRTM_SOLVER_PADE_ADD | XRTM_SOLVER_SINGLE | XRTM_SOLVER_SOS | XRTM_SOLVER_TWO_OS | DEV_RANGE_OVER_N_STOKES;

     for (i = 1; i <= 4; ++i) {
          if (test_group(td, &test_group_simple, XRTM_OPTION_PHASE_MATRIX_GC | XRTM_OPTION_VECTOR, solvers_mask, DEFAULT_N_QUAD_VECTOR, i)) {
               eprintf("ERROR: test_group()\n");
               return -1;
          }
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     printf("Range over n_derivs: index = %d\n", td->index);

     solvers_mask = XRTM_SOLVER_DOUB_ADD | XRTM_SOLVER_EIG_ADD | XRTM_SOLVER_EIG_BVP | XRTM_SOLVER_MEM_BVP | XRTM_SOLVER_PADE_ADD | XRTM_SOLVER_SINGLE | XRTM_SOLVER_SOS | XRTM_SOLVER_TWO_OS | DEV_RANGE_OVER_N_DERIVS;

     max_pow = MAX_N_DERIV_POWER;
     if (td->quick_run)
          max_pow = 1;

     for (i = 0; i <= max_pow; ++i) {
          n = 1 << i;

          if (test_group_vary_derivs2(td, &test_group_vary_simple, 0, solvers_mask, DEFAULT_N_QUAD_SCALER, DEFAULT_N_QUAD_VECTOR, n)) {
               eprintf("ERROR: test_group_vary_derivs2()\n");
               return -1;
          }
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     printf("Range over n_layers: index = %d\n", td->index);

     solvers_mask = XRTM_SOLVER_DOUB_ADD | XRTM_SOLVER_EIG_ADD | XRTM_SOLVER_EIG_BVP | XRTM_SOLVER_MEM_BVP | XRTM_SOLVER_PADE_ADD | XRTM_SOLVER_SINGLE | XRTM_SOLVER_SOS | XRTM_SOLVER_TWO_OS | DEV_RANGE_OVER_N_LAYERS;

     max_pow = MAX_N_LAYER_POWER;
     if (td->quick_run)
          max_pow = 1;

     for (i = 0; i <= max_pow; ++i) {
          n = 1 << i;

          if (test_group_vary_layers2(td, &test_group_vary_simple, 0, solvers_mask, DEFAULT_N_QUAD_SCALER, DEFAULT_N_QUAD_VECTOR, n)) {
               eprintf("ERROR: test_group_vary_layers2()\n");
               return -1;
          }
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     printf("Range over n_out_levels: index = %d\n", td->index);


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     printf("Range over n_out_thetas: index = %d\n", td->index);

     solvers_mask = XRTM_SOLVER_DOUB_ADD | XRTM_SOLVER_EIG_ADD | XRTM_SOLVER_EIG_BVP | XRTM_SOLVER_MEM_BVP | XRTM_SOLVER_PADE_ADD | XRTM_SOLVER_SINGLE | XRTM_SOLVER_SOS | XRTM_SOLVER_TWO_OS | DEV_RANGE_OVER_N_OUT_THETAS;

     max_pow = MAX_N_THETA_POWER;
     if (td->quick_run)
          max_pow = 1;

     for (i = 0; i <= max_pow; ++i) {
          n = 1 << i;

          if (test_group_vary_thetas2(td, &test_group_vary_simple, 0, solvers_mask, DEFAULT_N_QUAD_SCALER, DEFAULT_N_QUAD_VECTOR, n)) {
               eprintf("ERROR: test_group_vary_thetas2()\n");
               return -1;
          }
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     printf("Range over optical property inputs: index = %d\n", td->index);

     if (test_group2(td, &test_group_simple, 0, XRTM_SOLVER_DOUB_ADD | XRTM_SOLVER_EIG_ADD | XRTM_SOLVER_EIG_BVP | XRTM_SOLVER_MEM_BVP | XRTM_SOLVER_PADE_ADD | XRTM_SOLVER_SINGLE | XRTM_SOLVER_TWO_OS | XRTM_SOLVER_SOS | DEV_RANGE_OVER_OPTICAL_PROPERTY_INPUTS, DEFAULT_N_QUAD_SCALER, DEFAULT_N_QUAD_VECTOR)) {
          eprintf("ERROR: test_group2()\n");
          return -1;
     }
/*
     if (test_group2(td, &test_group_bounds, 0, XRTM_SOLVER_DOUB_ADD | XRTM_SOLVER_EIG_ADD | XRTM_SOLVER_EIG_BVP | XRTM_SOLVER_MEM_BVP | XRTM_SOLVER_PADE_ADD | XRTM_SOLVER_SINGLE | XRTM_SOLVER_TWO_OS | XRTM_SOLVER_SOS,                                          DEFAULT_N_QUAD_SCALER, DEFAULT_N_QUAD_VECTOR)) {
          eprintf("ERROR: test_group2()\n");
          return -1;
     }

     if (test_group2(td, &test_group_bounds, 0, XRTM_SOLVER_DOUB_ADD | XRTM_SOLVER_EIG_ADD | XRTM_SOLVER_EIG_BVP | XRTM_SOLVER_MEM_BVP | XRTM_SOLVER_PADE_ADD | XRTM_SOLVER_SINGLE | XRTM_SOLVER_TWO_OS | XRTM_SOLVER_SOS | DEV_RANGE_OVER_OPTICAL_PROPERTY_INPUTS, DEFAULT_N_QUAD_SCALER, DEFAULT_N_QUAD_VECTOR)) {
          eprintf("ERROR: test_group2()\n");
          return -1;
     }
*/

     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
/*
     printf("Range over eigen sovers: index = %d\n", td->index);

     misc_input_init(&misc_input);

     misc_input.eigen_solver_gen_real = EIGEN_SOLVER_GEN_REAL_ASYMTX;

     if (test_group(td, &test_group_simple, 0, XRTM_SOLVER_EIG_ADD | XRTM_SOLVER_EIG_BVP | XRTM_SOLVER_MEM_BVP | DEV_RANGE_OVER_EIGEN_SOVERS, DEFAULT_N_QUAD_SCALER, 1)) {
          eprintf("ERROR: test_group()\n");
          return -1;
     }

     misc_input.eigen_solver_gen_real    = EIGEN_SOLVER_GEN_REAL_EISPACK;
     misc_input.eigen_solver_gen_complex = EIGEN_SOLVER_GEN_COMPLEX_EISPACK;

     if (test_group2(td, &test_group_simple, 0, XRTM_SOLVER_EIG_ADD | XRTM_SOLVER_EIG_BVP | XRTM_SOLVER_MEM_BVP | DEV_RANGE_OVER_EIGEN_SOVERS, DEFAULT_N_QUAD_SCALER, DEFAULT_N_QUAD_VECTOR)) {
          eprintf("ERROR: test_group2()\n");
          return -1;
     }

     misc_input.eigen_solver_gen_real    = EIGEN_SOLVER_GEN_REAL_LAPACK;
     misc_input.eigen_solver_gen_complex = EIGEN_SOLVER_GEN_COMPLEX_LAPACK;

     if (test_group2(td, &test_group_simple, 0, XRTM_SOLVER_EIG_ADD | XRTM_SOLVER_EIG_BVP | XRTM_SOLVER_MEM_BVP | DEV_RANGE_OVER_EIGEN_SOVERS, DEFAULT_N_QUAD_SCALER, DEFAULT_N_QUAD_VECTOR)) {
          eprintf("ERROR: test_group2()\n");
          return -1;
     }
*/

     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     printf("Range over options: index = %d\n", td->index);

     for (i = 0; i < XRTM_OPTION_VECTOR; ++i) {
          if (i & XRTM_OPTION_CALC_DERIVS)
               continue;
/*
          if (i & XRTM_OPTION_FORWARD_DERIVS)
               continue;
*/
          if (i & XRTM_OPTION_REVERSE_DERIVS)
               continue;

          if (i & XRTM_OPTION_OUTPUT_AT_TAUS)
               continue;

          if (i & XRTM_OPTION_PHASE_MATRIX_LC)
               continue;

          if (i & XRTM_OPTION_QUAD_LOBATTO)
               continue;

          if (i & XRTM_OPTION_SAVE_LEG_POLYS)
               continue;
          if (i & XRTM_OPTION_SAVE_PHASE_MATS)
               continue;
          if (i & XRTM_OPTION_SAVE_LOCAL_R_T)
               continue;
          if (i & XRTM_OPTION_SAVE_LAYER_R_T_S)
               continue;
          if (i & XRTM_OPTION_SAVE_TOTAL_R_T_S)
               continue;

          if (i & XRTM_OPTION_SOURCE_THERMAL)
               continue;

          if (i & XRTM_OPTION_STACK_REUSE_ADDING)
               continue;

          if (i & XRTM_OPTION_TOP_DOWN_ADDING)
               continue;
          if (i & XRTM_OPTION_BOTTOM_UP_ADDING)
               continue;

          if (i & XRTM_OPTION_UPWELLING_OUTPUT)
               continue;
          if (i & XRTM_OPTION_DOWNWELLING_OUTPUT)
               continue;

          if (i & XRTM_OPTION_VECTOR)
               continue;

          if (! (i & XRTM_OPTION_DELTA_M) && (i & XRTM_OPTION_N_T_TMS))
               continue;
/*
          if (! (i & XRTM_OPTION_N_T_TMS) && (i & XRTM_OPTION_FOUR_CONV_NEW))
               continue;
*/
          if ((i & XRTM_OPTION_N_T_TMS) && (i & XRTM_OPTION_NO_AZIMUTHAL))
               continue;
          mask = i & (XRTM_OPTION_FOUR_CONV_OLD | XRTM_OPTION_FOUR_CONV_NEW);
          if (mask == 0 || mask == (XRTM_OPTION_FOUR_CONV_OLD | XRTM_OPTION_FOUR_CONV_NEW))
               continue;
          mask = i & (XRTM_OPTION_OUTPUT_AT_LEVELS | XRTM_OPTION_OUTPUT_AT_TAUS);
          if (mask == 0 || mask == (XRTM_OPTION_OUTPUT_AT_LEVELS | XRTM_OPTION_OUTPUT_AT_TAUS))
              continue;
          mask = i & (XRTM_OPTION_QUAD_NORM_GAUS_LEG | XRTM_OPTION_QUAD_DOUB_GAUS_LEG | XRTM_OPTION_QUAD_LOBATTO);
          if (mask == 0 || mask == (XRTM_OPTION_QUAD_NORM_GAUS_LEG | XRTM_OPTION_QUAD_DOUB_GAUS_LEG) ||
                           mask == (XRTM_OPTION_QUAD_NORM_GAUS_LEG | XRTM_OPTION_QUAD_LOBATTO      ) ||
                           mask == (XRTM_OPTION_QUAD_DOUB_GAUS_LEG | XRTM_OPTION_QUAD_LOBATTO      ))
               continue;
          mask = i & (XRTM_OPTION_PHASE_SCALAR | XRTM_OPTION_PHASE_MATRIX_GC | XRTM_OPTION_PHASE_MATRIX_LC);
          if (mask == 0 || mask == (XRTM_OPTION_PHASE_SCALAR    | XRTM_OPTION_PHASE_MATRIX_GC)
                        || mask == (XRTM_OPTION_PHASE_SCALAR    | XRTM_OPTION_PHASE_MATRIX_LC)
                        || mask == (XRTM_OPTION_PHASE_MATRIX_GC | XRTM_OPTION_PHASE_MATRIX_LC))
              continue;

          solvers_mask = 0xffffffff;

          if (i & XRTM_OPTION_NO_AZIMUTHAL)
               solvers_mask &= 0xffffffff ^ (XRTM_SOLVER_SINGLE);
          if (i & XRTM_OPTION_SFI)
               solvers_mask &= 0xffffffff ^ (XRTM_SOLVER_DOUB_ADD | XRTM_SOLVER_EIG_ADD | XRTM_SOLVER_MEM_BVP | XRTM_SOLVER_PADE_ADD | XRTM_SOLVER_SINGLE | XRTM_SOLVER_TWO_OS | XRTM_SOLVER_SOS);

          DEV_OPT_OUT_RANGE_OVER_OPTIONS(solvers_mask);

          if (i & XRTM_OPTION_PHASE_SCALAR) {
               if (test_group(td, &test_group_simple, i, solvers_mask, DEFAULT_N_QUAD_SCALER, 1)) {
                    eprintf("ERROR: test_group()\n");
                    return -1;
               }
          }
          else {
               solvers_mask2 = solvers_mask & XRTM_SOLVERS_VECTOR;

               if (test_group(td, &test_group_simple, i, solvers_mask2, DEFAULT_N_QUAD_VECTOR, 4)) {
                    eprintf("ERROR: test_group()\n");
                    return -1;
               }
          }
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     printf("Range over saving options: index = %d\n", td->index);

     for (i = 0; i < XRTM_OPTION_VECTOR; ++i) {
          if (i & (0xffffffff ^ (XRTM_OPTION_SAVE_LEG_POLYS | XRTM_OPTION_SAVE_PHASE_MATS | XRTM_OPTION_SAVE_LOCAL_R_T | XRTM_OPTION_SAVE_LAYER_R_T_S | XRTM_OPTION_SAVE_TOTAL_R_T_S)))
               continue;

          solvers_mask = XRTM_SOLVER_DOUB_ADD | XRTM_SOLVER_EIG_ADD | XRTM_SOLVER_EIG_BVP | XRTM_SOLVER_MEM_BVP | XRTM_SOLVER_PADE_ADD | XRTM_SOLVER_SINGLE | XRTM_SOLVER_SOS | XRTM_SOLVER_TWO_OS | DEV_RANGE_OVER_SAVING_OPTIONS;

          if (test_group2(td, &test_group_simple, 0, 0xffffffff, DEFAULT_N_QUAD_SCALER, DEFAULT_N_QUAD_VECTOR)) {
               eprintf("ERROR: test_group2()\n");
               return -1;
          }
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     printf("Layer saving over multiple calls: index = %d\n", td->index);


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     printf("Stack layer saving configurations and over multiple calls: index = %d\n", td->index);


     return 0;
}
