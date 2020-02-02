/*******************************************************************************
**
**    Copyright (C) 2007-2020 Greg McGarragh <greg.mcgarragh@colostate.edu>
**
**    This source code is licensed under the GNU General Public License (GPL),
**    Version 3.  See the file COPYING for more details.
**
*******************************************************************************/

#ifndef TEST_XRTM_H
#define TEST_XRTM_H

#include <xrtm_interface.h>

#include "test.h"

#ifdef __cplusplus
extern "C" {
#endif


typedef struct {
     char *name;
     struct list_data *prev;
     struct list_data *next;

     int set_flag_doub_d_tau;
     double doub_d_tau;

     int set_flag_pade_params;
     int pade_s;
     int pade_r;

     int set_flag_sos_params;
     int sos_max_os;
     double sos_max_tau;
     double sos_tol;

     int set_flag_fourier_tol;
     double fourier_tol;

     int set_flag_lambda;
     double lambda;

     int set_flag_planet_r;
     double planet_r;

     int set_flag_levels_z;
     double *levels_z;

     int set_flag_out_levels;
     int *out_levels;

     int set_flag_out_taus;
     double *out_taus;

     int set_flag_out_thetas;
     double *out_thetas;

     int set_flag_F_iso_top;
     double F_iso_top;

     int set_flag_F_iso_bot;
     double F_iso_bot;

     int set_flag_F_0;
     double F_0;

     int set_flag_theta_0;
     double theta_0;

     int set_flag_phi_0;
     double phi_0;

     int set_flag_levels_b;
     double *levels_b;
     int set_flag_levels_b_l;
     double **levels_b_l;

     int set_flag_surface_b;
     double surface_b;

     int set_flag_g;
     double *g;			/* n_layers		*/
     int set_flag_g_l;
     double **g_l;		/* n_layers, n_derivs	*/

     int set_flag_coef;
     int *n_coef;		/* n_layers		*/
     double ***coef;		/* n_layers		*/
     int set_flag_coef_l;
     double ****coef_l;		/* n_layers, n_derivs	*/

     int set_flag_omega;
     double *omega;		/* n_layers		*/
     int set_flag_omega_l;
     double **omega_l;		/* n_layers, n_derivs	*/

     int set_flag_ltau;
     double *ltau;		/* n_layers		*/
     int set_flag_ltau_l;
     double **ltau_l;		/* n_layers, n_derivs	*/

     int set_flag_ampfac;
     double *ampfac;		/* n_kernels				*/
     int set_flag_param;
     double **param;		/* n_kernels, n_params			*/
     int set_flag_ampfac_l;
     double **ampfac_l;		/* n_kernels, n_derivs			*/
     int set_flag_param_l;
     double ***param_l;		/* n_kernels, n_derivs, n_params	*/
} xrtm_input_data;


typedef struct {
     char *name;
     struct list_data *prev;
     struct list_data *next;

     int n_solvers;
     long *solver_list;

     int n_phi_bounds;
     double *phi;

     double *tol;
     double *tol_l;

     int n_ignore_index_solver;
     int *ignore_list_index_solver;
     int *ignore_mask_index_solver;

     int n_ignore_deriv;
     int *ignore_list_deriv;
     int *ignore_mask_deriv;

     int (*ignore_mask_solver)(double theta_0, double *theta);

     int options;
     enum xrtm_solver_mask solvers;
     int max_coef;
     int n_quad;
     int n_stokes;
     int n_derivs;
     int n_layers;
     int n_kernels;
     int n_kernel_quad;
     enum xrtm_kernel_type *kernels;	/* n_kernels	*/
     int n_out_levels;
     int n_out_thetas;

     xrtm_input_data xrtm_input_list;

     test_cmd_data test_cmd_list;
     test_result_data test_result_list;
} test_xrtm_data;


#include "prototypes/test_xrtm_p.h"


#ifdef __cplusplus
}
#endif

#endif /* TEST_XRTM_H */
