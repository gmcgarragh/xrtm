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

#include "test.h"
#include "test_xrtm.h"


/*******************************************************************************
 *
 ******************************************************************************/
int xrtm_input_init(xrtm_input_data *d) {

     d->set_flag_doub_d_tau  = 0;
     d->set_flag_pade_params = 0;
     d->set_flag_sos_params  = 0;
     d->set_flag_fourier_tol = 0;
     d->set_flag_lambda      = 0;
     d->set_flag_planet_r    = 0;
     d->set_flag_levels_z    = 0;
     d->levels_z             = NULL;
     d->set_flag_out_levels  = 0;
     d->out_levels           = NULL;
     d->set_flag_out_taus    = 0;
     d->out_taus             = NULL;
     d->set_flag_out_thetas  = 0;
     d->out_thetas           = NULL;
     d->set_flag_F_iso_top   = 0;
     d->set_flag_F_iso_bot   = 0;
     d->set_flag_F_0         = 0;
     d->set_flag_theta_0     = 0;
     d->set_flag_phi_0       = 0;
     d->set_flag_levels_b    = 0;
     d->levels_b             = NULL;
     d->set_flag_surface_b   = 0;
     d->levels_b_l           = NULL;
     d->set_flag_g           = 0;
     d->g                    = NULL;
     d->set_flag_g_l         = 0;
     d->g_l                  = NULL;
     d->set_flag_coef        = 0;
     d->n_coef               = NULL;
     d->coef                 = NULL;
     d->set_flag_coef_l      = 0;
     d->coef_l               = NULL;
     d->set_flag_omega       = 0;
     d->omega                = NULL;
     d->set_flag_omega_l     = 0;
     d->omega_l              = NULL;
     d->set_flag_ltau        = 0;
     d->ltau                 = NULL;
     d->set_flag_ltau_l      = 0;
     d->ltau_l               = NULL;
     d->set_flag_ampfac      = 0;
     d->ampfac               = NULL;
     d->set_flag_param       = 0;
     d->param                = NULL;
     d->set_flag_ampfac_l    = 0;
     d->ampfac_l             = NULL;
     d->set_flag_param_l     = 0;
     d->param_l              = NULL;

     return 0;
}



void xrtm_input_free(xrtm_input_data *d) {

     free_array1_d(d->levels_z);
     free_array1_i(d->out_levels);
     free_array1_d(d->out_taus);
     free_array1_d(d->out_thetas);
     free_array1_d(d->levels_b);
     free_array2_d(d->levels_b_l);
     free_array1_d(d->g);
     free_array2_d(d->g_l);
     free_array1_i(d->n_coef);
     free_array1((void * ) d->coef);
     free_array2((void **) d->coef_l);
     free_array1_d(d->omega);
     free_array2_d(d->omega_l);
     free_array1_d(d->ltau);
     free_array2_d(d->ltau_l);
     free_array1_d(d->ampfac);
     free_array2_d(d->param);
     free_array2_d(d->ampfac_l);
     free_array3_d(d->param_l);
}



/*******************************************************************************
 *
 ******************************************************************************/
int test_xrtm_init(test_xrtm_data *d) {

     d->solver_list              = NULL;

     d->phi                      = NULL;

     d->tol                      = NULL;
     d->tol_l                    = NULL;

     d->ignore_list_index_solver = NULL;
     d->ignore_mask_index_solver = NULL;

     d->ignore_list_deriv        = NULL;
     d->ignore_mask_deriv        = NULL;

     list_init(&d->xrtm_input_list);
     list_init(&d->test_cmd_list);
     list_init(&d->test_result_list);

     return 0;
}



void test_xrtm_free(test_xrtm_data *d) {

     xrtm_input_data *xrtm_input;
     test_cmd_data *test_cmd;
     test_result_data *test_result;

     free_array1_l(d->solver_list);

     free_array1_d(d->phi);

     free_array1_d(d->tol);
     free_array1_d(d->tol_l);

     free_array1_i(d->ignore_list_index_solver);
     free_array1_i(d->ignore_mask_index_solver);

     free_array1_i(d->ignore_list_deriv);
     free_array1_i(d->ignore_mask_deriv);

     list_for_each(&d->xrtm_input_list, xrtm_input)
          xrtm_input_free(xrtm_input);
     list_free(&d->xrtm_input_list);

     list_for_each(&d->test_cmd_list, test_cmd)
          test_cmd_free(test_cmd);
     list_free(&d->test_cmd_list);

     list_for_each(&d->test_result_list, test_result)
          test_result_free(test_result);
     list_free(&d->test_result_list);

     free_array1_i((int *) d->kernels);
}



int test_xrtm_fill(test_xrtm_data *d, int options, int solvers, int max_coef, int n_quad, int n_stokes, int n_derivs, int n_layers, int n_kernels, int n_kernel_quad, enum xrtm_kernel_type *kernels, int n_out_levels, int n_out_thetas) {

     d->options       = options;
     d->solvers       = solvers;
     d->max_coef      = max_coef;
     d->n_quad        = n_quad;
     d->n_stokes      = n_stokes;
     d->n_derivs      = n_derivs;
     d->n_layers      = n_layers;
     d->n_kernels     = n_kernels;
     d->n_kernel_quad = n_kernel_quad;
     d->kernels       = (enum xrtm_kernel_type *) dup_array1_i((int *) kernels, n_kernels);
     d->n_out_levels  = n_out_levels;
     d->n_out_thetas  = n_out_thetas;

     return 0;
}
