/*******************************************************************************
**
**    Copyright (C) 2007-2020 Greg McGarragh <greg.mcgarragh@colostate.edu>
**
**    This source code is licensed under the GNU General Public License (GPL),
**    Version 3.  See the file COPYING for more details.
**
*******************************************************************************/

#include <gutil.h>
#include <glist.h>

#include <xrtm_support.h>

#include "test.h"
#include "test_execute.h"
#include "test_loop.h"
#include "test_macros.h"
#include "test_util.h"

#include <omp.h>


/*******************************************************************************
 *
 ******************************************************************************/
static void set_coefs_file(misc_data *md, test_data *td, int i_layer, double **coefs) {

     int i;

     for (i = 0; i < MAX_PFS; ++i) {
          if (coefs == td->gc[i])
               break;
     }

     md->coef_files[i_layer] = strdup(td->gf[i]);
}



static void set_coefs_file_derivs(misc_data *md, test_data *td, int i_layer, int n_derivs, double ***coefs_l) {

     int i;
     int j;

     for (i = 0; i < n_derivs; ++i) {
          for (j = 0; j < MAX_PFS; ++j) {
               if (coefs_l[i] == td->gc_l[j])
                    break;
          }

          if (j == MAX_PFS)
               md->coef_files_l[i_layer][i] = strdup("none");
          else
               md->coef_files_l[i_layer][i] = strdup(td->gf_l[j]);
     }
}



static void set_coefs_file_layers(misc_data *md, test_data *td, int n_layers, double ***coefs) {

     int i;

     for (i = 0; i < n_layers; ++i)
          set_coefs_file(md, td, i, coefs[i]);

}



/*******************************************************************************
 *
 ******************************************************************************/
static int test_loop_inner(test_data *td, int index, test_xrtm_data *test_xrtm) {

     int p;
     int q;

     int index2;
/*
     int count;
*/
     int ignore_mask_solver2;

     misc_data misc;

     xrtm_data xrtm;

     xrtm_input_data *xrtm_input;

     XRTM_CREATE(&xrtm, test_xrtm->options, test_xrtm->solvers, test_xrtm->max_coef, test_xrtm->n_quad, test_xrtm->n_stokes, test_xrtm->n_derivs, test_xrtm->n_layers, 1, test_xrtm->n_kernels, test_xrtm->n_kernel_quad, test_xrtm->kernels, test_xrtm->n_out_levels, test_xrtm->n_out_thetas);

     misc_init(&xrtm, &misc, -1);
/*
     count = list_count(&test_xrtm->xrtm_input_list);
*/
     index2 = 0;

     list_for_each(&test_xrtm->xrtm_input_list, xrtm_input) {
/*
          if (i % 1 == 0)
               printf("     %d %d\n", count, i);
*/
          if (test_xrtm->solvers & XRTM_SOLVERS_DOUBLING &&
               xrtm_input->set_flag_doub_d_tau) {
               XRTM_SET_DOUB_D_TAU(&xrtm, xrtm_input->doub_d_tau);
          }

          if (test_xrtm->solvers & XRTM_SOLVER_PADE_ADD &&
               xrtm_input->set_flag_pade_params) {
               XRTM_SET_PADE_PARAMS(&xrtm, xrtm_input->pade_s, xrtm_input->pade_r);
          }

          if (test_xrtm->solvers & XRTM_SOLVER_SOS &&
               xrtm_input->set_flag_sos_params) {
               XRTM_SET_SOS_PARAMS(&xrtm, xrtm_input->sos_max_os, xrtm_input->sos_max_tau, xrtm_input->sos_tol);
          }

          if (xrtm_input->set_flag_fourier_tol)
               XRTM_SET_FOURIER_TOL(&xrtm, xrtm_input->fourier_tol);

          if (test_xrtm->options & XRTM_OPTION_PSA) {
               if (xrtm_input->set_flag_planet_r)
                    XRTM_SET_PLANET_R(&xrtm, xrtm_input->planet_r);
               if (xrtm_input->set_flag_levels_z)
                    XRTM_SET_LEVELS_Z(&xrtm, xrtm_input->levels_z);
          }

          if (xrtm_input->set_flag_out_levels)
               XRTM_SET_OUT_LEVELS(&xrtm, xrtm_input->out_levels);

          if (xrtm_input->set_flag_out_taus)
               XRTM_SET_OUT_TAUS(&xrtm, xrtm_input->out_taus);

          if (xrtm_input->set_flag_out_thetas)
               XRTM_SET_OUT_THETAS(&xrtm, xrtm_input->out_thetas);

          if (xrtm_input->set_flag_F_iso_top)
               XRTM_SET_F_ISO_TOP(&xrtm, xrtm_input->F_iso_top);

          if (xrtm_input->set_flag_F_iso_bot)
               XRTM_SET_F_ISO_BOT(&xrtm, xrtm_input->F_iso_bot);

          if (test_xrtm->options & XRTM_OPTION_SOURCE_SOLAR) {
               if (xrtm_input->set_flag_F_0)
                    XRTM_SET_F_0(&xrtm, xrtm_input->F_0);

               if (xrtm_input->set_flag_theta_0)
                    XRTM_SET_THETA_0(&xrtm, xrtm_input->theta_0);

               if (xrtm_input->set_flag_phi_0)
                    XRTM_SET_PHI_0  (&xrtm, xrtm_input->phi_0);
          }

          if (test_xrtm->options & XRTM_OPTION_SOURCE_THERMAL) {
               if (xrtm_input->set_flag_lambda)
                    XRTM_SET_LAMBDA(&xrtm, xrtm_input->lambda);

               if (xrtm_input->set_flag_levels_b)
                    XRTM_SET_LEVELS_B(&xrtm, xrtm_input->levels_b);

               if (xrtm_input->set_flag_surface_b)
                    XRTM_SET_SURFACE_B(&xrtm, xrtm_input->surface_b);

          }

          if (xrtm_input->set_flag_coef) {
               XRTM_SET_COEF_N(&xrtm, xrtm_input->n_coef, xrtm_input->coef);
               set_coefs_file_layers(&misc, td, test_xrtm->n_layers, xrtm_input->coef);
          }

          if (xrtm_input->set_flag_omega)
               XRTM_SET_OMEGA_N(&xrtm, xrtm_input->omega);

          if (xrtm_input->set_flag_ltau)
               XRTM_SET_LTAU_N(&xrtm, xrtm_input->ltau);

          if (xrtm_input->set_flag_ampfac)
               XRTM_SET_KERNEL_AMPFAC(&xrtm, 0, xrtm_input->ampfac[0]);

          if (test_xrtm->n_derivs > 0) {
               for (p = 0; p < test_xrtm->n_layers; ++p) {
                    for (q = 0; q < test_xrtm->n_derivs; ++q) {
                         if (xrtm_input->set_flag_coef &&
                             xrtm_input->coef_l[p][q])
                              XRTM_SET_COEF_L_11(&xrtm, p, q, xrtm_input->coef_l[p][q]);
                    }

                    set_coefs_file_derivs(&misc, td, p, test_xrtm->n_derivs, xrtm_input->coef_l[p]);

                    if (xrtm_input->set_flag_omega)
                         XRTM_SET_OMEGA_L_1N(&xrtm, p, xrtm_input->omega_l[p]);
                    if (xrtm_input->set_flag_ltau)
                         XRTM_SET_LTAU_L_1N (&xrtm, p, xrtm_input->ltau_l [p]);
               }

               if (xrtm_input->set_flag_ampfac_l)
                    XRTM_SET_KERNEL_AMPFAC_L_N(&xrtm, 0, xrtm_input->ampfac_l[0]);

               if (xrtm_update_varied_layers(&xrtm)) {
                    fprintf(stderr, "ERROR: xrtm_update_varied_layers()\n");
                    return -1;
               }
          }

          ignore_mask_solver2 = 0;
          if (test_xrtm->ignore_mask_solver)
               ignore_mask_solver2 = test_xrtm->ignore_mask_solver(xrtm_input->theta_0, xrtm_input->out_thetas);

          HANDLE_RETURN(test_execute(&xrtm, &misc, td, index, index2, &test_xrtm->test_cmd_list, &test_xrtm->test_result_list, test_xrtm->n_solvers, test_xrtm->solver_list, test_xrtm->n_phi_bounds, test_xrtm->phi, test_xrtm->tol, test_xrtm->tol_l, test_xrtm->n_ignore_index_solver, test_xrtm->ignore_list_index_solver, test_xrtm->ignore_mask_index_solver, ignore_mask_solver2, ignore_mask_solver2, test_xrtm->n_ignore_deriv, test_xrtm->ignore_list_deriv, test_xrtm->ignore_mask_deriv), "test_execute()");

          index2++;
     }

     misc_free(&misc);

     XRTM_DESTROY(&xrtm);

     return 0;
}



/*******************************************************************************
 *
 ******************************************************************************/
int test_loop_outer(test_data *d) {

     int i;
/*
     int count;
*/
     test_xrtm_data *test_xrtm;
/*
     count = list_count(&d->test_xrtm_list);
*/
     i = 0;
     list_for_each(&d->test_xrtm_list, test_xrtm) {
/*
          if (i % 1 == 0)
               printf("%d %d\n", count, i);
*/
          if (test_loop_inner(d, i, test_xrtm)) {
               fprintf(stderr, "ERROR: test_loop_inner()\n");
               return -1;
          }

          ++i;
     }

     return 0;
}



int test_loop_outer2(test_data *td) {

     int i;

     int count;

     test_xrtm_data *test_xrtm;

     test_xrtm_data **test_xrtm_array;

     count = list_count(&td->test_xrtm_list);

     test_xrtm_array = malloc(count * sizeof(test_xrtm_data *));

     i = 0;
     list_for_each(&td->test_xrtm_list, test_xrtm) {
          test_xrtm_array[i] = test_xrtm;
          ++i;
     }
#pragma omp parallel
{
#pragma omp for
     for (i = 0; i < count; ++i) {
/*
          if (i % 1 == 0)
               printf("%d %d\n", count, i);
*/
          if (test_loop_inner(td, i, test_xrtm_array[i])) {
               fprintf(stderr, "ERROR: test_loop_inner()\n");
               exit(1);
          }
     }
}

     free(test_xrtm_array);

     return 0;
}
