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

#include <omp.h>

#include "input_util.h"

#include "test.h"
#include "test_result.h"
#include "test_util.h"
#include "test_write.h"


/*******************************************************************************
 *
 ******************************************************************************/
/*
static int stokes_element_zero(int i, double theta_0, double theta, double phi) {

     if (! (fabs(phi -   0.) < DBL_EPSILON) && ! (fabs(phi -   45.) < DBL_EPSILON) && ! (fabs(phi -   90.) < DBL_EPSILON) && ! (fabs(phi -   135.) < DBL_EPSILON) && ! (fabs(phi -   180.) < DBL_EPSILON)) {
          if (theta_0 < DBL_EPSILON && theta < DBL_EPSILON) {
               if (i == 0 || i == 1 || i == 2) return 0;
          }
          if (theta_0 > DBL_EPSILON && theta < DBL_EPSILON) {
               if (i == 0 || i == 1 || i == 2) return 0;
          }
          if (theta_0 < DBL_EPSILON && theta > DBL_EPSILON) {
               if (i == 0 || i == 1) return 0;
          }
          if (theta_0 > DBL_EPSILON && theta > DBL_EPSILON) {
               if (i == 0 || i == 1 || i == 2 || i == 3) return 0;
          }

          if (theta_0 < DBL_EPSILON && theta < DBL_EPSILON) {
               if (i == 0 || i == 1 || i == 2) return 0;
          }
          if (theta_0 > DBL_EPSILON && theta < DBL_EPSILON) {
               if (i == 0 || i == 1 || i == 2) return 0;
          }
          if (theta_0 < DBL_EPSILON && theta > DBL_EPSILON) {
               if (i == 0 || i == 1) return 0;
          }
          if (theta_0 > DBL_EPSILON && theta > DBL_EPSILON) {
               if (i == 0 || i == 1 || i == 2 || i == 3) return 0;
          }
     }
     else {
          if (theta_0 < DBL_EPSILON && theta < DBL_EPSILON && fabs(phi -   0.) < DBL_EPSILON) {
               if (i == 0 || i == 1) return 0;
          }
          if (theta_0 > DBL_EPSILON && theta < DBL_EPSILON && fabs(phi -   0.) < DBL_EPSILON) {
               if (i == 0 || i == 1) return 0;
          }
          if (theta_0 < DBL_EPSILON && theta > DBL_EPSILON && fabs(phi -   0.) < DBL_EPSILON) {
               if (i == 0 || i == 1) return 0;
          }
          if (theta_0 > DBL_EPSILON && theta > DBL_EPSILON && fabs(phi -   0.) < DBL_EPSILON) {
               if (i == 0 || i == 1) return 0;
          }

          if (theta_0 < DBL_EPSILON && theta < DBL_EPSILON && fabs(phi -  45.) < DBL_EPSILON) {
               if (i == 0 || i == 2) return 0;
          }
          if (theta_0 > DBL_EPSILON && theta < DBL_EPSILON && fabs(phi -  45.) < DBL_EPSILON) {
               if (i == 0 || i == 2) return 0;
          }
          if (theta_0 < DBL_EPSILON && theta > DBL_EPSILON && fabs(phi -  45.) < DBL_EPSILON) {
               if (i == 0 || i == 1) return 0;
          }
          if (theta_0 > DBL_EPSILON && theta > DBL_EPSILON && fabs(phi -  45.) < DBL_EPSILON) {
               if (i == 0 || i == 1 || i == 2 || i == 3) return 0;
          }

          if (theta_0 < DBL_EPSILON && theta < DBL_EPSILON && fabs(phi -  90.) < DBL_EPSILON) {
               if (i == 0 || i == 1) return 0;
          }
          if (theta_0 > DBL_EPSILON && theta < DBL_EPSILON && fabs(phi -  90.) < DBL_EPSILON) {
               if (i == 0 || i == 1) return 0;
          }
          if (theta_0 < DBL_EPSILON && theta > DBL_EPSILON && fabs(phi -  90.) < DBL_EPSILON) {
               if (i == 0 || i == 1) return 0;
          }
          if (theta_0 > DBL_EPSILON && theta > DBL_EPSILON && fabs(phi -  90.) < DBL_EPSILON) {
               if (i == 0 || i == 1 || i == 2 || i == 3) return 0;
          }

          if (theta_0 < DBL_EPSILON && theta < DBL_EPSILON && fabs(phi - 135.) < DBL_EPSILON) {
               if (i == 0 || i == 2) return 0;
          }
          if (theta_0 > DBL_EPSILON && theta < DBL_EPSILON && fabs(phi - 135.) < DBL_EPSILON) {
               if (i == 0 || i == 2) return 0;
          }
          if (theta_0 < DBL_EPSILON && theta > DBL_EPSILON && fabs(phi - 135.) < DBL_EPSILON) {
               if (i == 0 || i == 1) return 0;
          }
          if (theta_0 > DBL_EPSILON && theta > DBL_EPSILON && fabs(phi - 135.) < DBL_EPSILON) {
               if (i == 0 || i == 1 || i == 2 || i == 3) return 0;
          }

          if (theta_0 < DBL_EPSILON && theta < DBL_EPSILON && fabs(phi - 180.) < DBL_EPSILON) {
               if (i == 0 || i == 1) return 0;
          }
          if (theta_0 > DBL_EPSILON && theta < DBL_EPSILON && fabs(phi - 180.) < DBL_EPSILON) {
               if (i == 0 || i == 1) return 0;
          }
          if (theta_0 < DBL_EPSILON && theta > DBL_EPSILON && fabs(phi - 180.) < DBL_EPSILON) {
               if (i == 0 || i == 1) return 0;
          }
          if (theta_0 > DBL_EPSILON && theta > DBL_EPSILON && fabs(phi - 180.) < DBL_EPSILON) {
               if (i == 0 || i == 1) return 0;
          }
     }

     return 1;
}
*/


/*******************************************************************************
 *
 ******************************************************************************/
static int cmp_p_i2(int *a, int *b) {

     if (*a < *b)
          return -1;
     if (*a > *b)
          return 1;

     return 0;
}



/*******************************************************************************
 *
 ******************************************************************************/
int test_execute(xrtm_data *gd, misc_data *md, test_data *td, int index, int index2, test_cmd_data *test_cmd_list, test_result_data *test_result_list, int n_solvers, long *solvers, int n_out_phis, double *out_phis, double *tol, double *tol_l, int n_ignore_index_solver, int *ignore_list_index_solver, int *ignore_mask_index_solver, int ignore_solver_mask_ref, int ignore_solver_mask_tran, int n_ignore_deriv, int *ignore_list_deriv, int *ignore_mask_deriv) {

     char temp[1024];
     char temp2[MAX_COMMAND_LENGTH];

     int i;

     int l;

     int flag;

     int options;

     int n_quad;
     int n_stokes;
     int n_derivs;
     int n_out_levels;
     int n_out_thetas;

     int n_quad_x;

     int n_mus;

     int n_solvers_all;

     int i_out_level;
     int i_deriv;
     int i_mu;
     int i_mu_offset;
     int i_out_phi;
     int i_stokes;
     int i_solver;
     int i_solver2;

     int fail_mask;

     int *i_ptr_index = NULL;
     int *i_ptr_deriv = NULL;
/*
     double theta_0;
     double phi_0;
*/
     double out_theta;

     double d1;
     double d2;

     double *qx;
     double *qw;

     double **out_phis2;

     double *****I_p;
     double *****I_m;

     double *****I_p_l2;
     double *****I_m_l2;

     double ******I_p_l;
     double ******I_m_l;

     test_cmd_data *test_cmd;

     test_result_data *test_result;


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (td->core_echo_xrtm_cmd) {
          if (get_callxrtm_cmd(gd, md, n_out_phis, out_phis, temp2, MAX_COMMAND_LENGTH)) {
               fprintf(stderr, "ERROR: write_callxrtm_cmd()\n");
               return -1;
          }
          fprintf(stderr, "index = %d, index2 = %d, ../xrtm/utils/%s\n", index, index2, temp2);
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (td->core_dont_execute)
          return 0;


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     options      = xrtm_get_options(gd);

     n_quad       = xrtm_get_n_quad(gd);
     n_stokes     = xrtm_get_n_stokes(gd);
     n_derivs     = xrtm_get_n_derivs(gd);
     n_out_levels = xrtm_get_n_out_levels(gd);
     n_out_thetas = xrtm_get_n_out_thetas(gd);

     n_quad_x  = n_quad + n_out_thetas;

     if (n_out_thetas == 0)
          n_mus = n_quad;
     else
          n_mus = n_out_thetas;


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     qx = alloc_array1_d(n_quad_x);
     if (xrtm_qx(gd, qx)) {
          fprintf(stderr, "ERROR: xrtm_qx()\n");
          return -1;
     }

     qw = alloc_array1_d(n_quad_x);
     if (xrtm_qw(gd, qw)) {
          fprintf(stderr, "ERROR: xrtm_qw()\n");
          return -1;
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     out_phis2 = (double **) alloc_array1(n_mus, sizeof(double *));

     for (i = 0; i < n_mus; ++i)
          out_phis2[i] = out_phis;


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     I_p = alloc_array5_d(n_solvers, n_out_levels, n_mus, n_out_phis, n_stokes);
     I_m = alloc_array5_d(n_solvers, n_out_levels, n_mus, n_out_phis, n_stokes);

     if (options & XRTM_OPTION_CALC_DERIVS) {
          I_p_l = alloc_array6_d(n_solvers, n_out_levels, n_derivs, n_mus, n_out_phis, n_stokes);
          I_m_l = alloc_array6_d(n_solvers, n_out_levels, n_derivs, n_mus, n_out_phis, n_stokes);
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     for (i = 0; i < n_solvers; ++i) {
          if (options & XRTM_OPTION_CALC_DERIVS) {
               I_p_l2 = I_p_l[i];
               I_m_l2 = I_m_l[i];
          }
#ifndef INCLUDE_DEV_SOURCE
          if ( ! omp_in_parallel()) {
#else
          if ( ! omp_in_parallel() ||
              (! (solvers[i] & (XRTM_SOLVER_DISORT |
                                XRTM_SOLVER_LIDORT |
                                XRTM_SOLVER_POLRAD |
                                XRTM_SOLVER_RADIANT |
                                XRTM_SOLVER_RADTRAN3 |
                                XRTM_SOLVER_VLIDORT)))) {
#endif
               if (xrtm_radiance(gd, solvers[i], n_out_phis, out_phis2, I_p[i], I_m[i], I_p_l2, I_m_l2)) {
                    if (get_callxrtm_cmd(gd, md, n_out_phis, out_phis, temp2, MAX_COMMAND_LENGTH)) {
                         fprintf(stderr, "ERROR: write_callxrtm_cmd()\n");
                         return -1;
                    }
                    fprintf(stderr, "ERROR: index = %d, index2 = %d, callxrtm = %s\n", index, index2, temp2);
                    fprintf(stderr, "ERROR: xrtm_radiance1()\n");
                    return -1;
               }
          }
          else {
               flag = 0;
#pragma omp critical
               if (xrtm_radiance(gd, solvers[i], n_out_phis, out_phis2, I_p[i], I_m[i], I_p_l2, I_m_l2)) {
                    flag = 1;
                    if (get_callxrtm_cmd(gd, md, n_out_phis, out_phis, temp2, MAX_COMMAND_LENGTH)) {
                         fprintf(stderr, "ERROR: write_callxrtm_cmd()\n");
                    }
                    else {
                         fprintf(stderr, "ERROR: index = %d, index2 = %d, callxrtm = %s\n", index, index2, temp2);
                         fprintf(stderr, "ERROR: xrtm_radiance2()\n");
                    }
               }

               if (flag) return -1;
          }
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (td->core_write_xrtm_cmd) {
          if (get_callxrtm_cmd(gd, md, n_out_phis, out_phis, temp2, MAX_COMMAND_LENGTH)) {
               fprintf(stderr, "ERROR: write_callxrtm_cmd()\n");
               return -1;
          }

          test_cmd = malloc(sizeof(test_cmd_data));

          test_cmd->index  = index;
          test_cmd->index2 = index2;

          test_cmd->cmd    = strdup(temp2);

          list_append(test_cmd_list, test_cmd, 0);
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     n_solvers_all = xrtm_solver_n();
/*
     theta_0 = xrtm_get_theta_0(gd);
     phi_0   = xrtm_get_phi_0(gd);
*/
     if (n_ignore_index_solver > 0)
          i_ptr_index = (int *) bsearch(&index, ignore_list_index_solver, n_ignore_index_solver, sizeof(int), (int (*)(const void *, const void *)) cmp_p_i2);

     flag = 0;

     fail_mask = 0;


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (n_out_thetas == 0)
          i_mu_offset = 0;
     else
          i_mu_offset = n_quad;

     for (i_out_level = 0; i_out_level < n_out_levels; ++i_out_level) {
          for (i_mu = 0; i_mu < n_mus; ++i_mu) {
               out_theta = acos(qx[i_mu_offset + i_mu]) * R2D;

               for (i_out_phi = 0; i_out_phi < n_out_phis; ++i_out_phi) {

                    for (i_stokes = 0; i_stokes < n_stokes; ++i_stokes) {
                         test_result = malloc(sizeof(test_result_data));
                         test_result_init(test_result, n_solvers);
                         list_append(test_result_list, test_result, 0);

                         save_result_header(test_result, index, index2, i_out_level, -1, out_theta, out_phis[i_out_phi], i_stokes);

                         for (i_solver = 0; i_solver < n_solvers_all; ++i_solver) {
                              for (i_solver2 = 0; i_solver2 < n_solvers; ++i_solver2) {
                                   if (xrtm_solver_index_to_mask(i_solver) == solvers[i_solver2])
                                        break;
                              }

                              if (i_solver2 == n_solvers)
                                   continue;

                              save_result_solver_mask(test_result, i_solver2, solvers[i_solver2]);
/*
                              if (stokes_element_zero(i_stokes, theta_0, out_theta, out_phis[i_out_phi] - phi_0)) {
                                   I_p[i_solver2][i_out_level][i_mu][i_out_phi][i_stokes] = 0.;
                                   I_m[i_solver2][i_out_level][i_mu][i_out_phi][i_stokes] = 0.;
                              }
*/
                              if (solvers[i_solver] == XRTM_SOLVER_EIG_BVP &&
                                  solvers[i_solver] == XRTM_SOLVER_MEM_BVP &&
                                  i_out_level == n_out_levels - 1 &&
                                  fabs(I_p[i_solver2][i_out_level][i_mu][i_out_phi][i_stokes]) < 1.e-15)
                                   I_p[i_solver2][i_out_level][i_mu][i_out_phi][i_stokes] = 0.;

                              if (solvers[i_solver] == XRTM_SOLVER_EIG_BVP &&
                                  solvers[i_solver] == XRTM_SOLVER_MEM_BVP &&
                                  i_out_level == 0 &&
                                  fabs(I_m[i_solver2][i_out_level][i_mu][i_out_phi][i_stokes]) < 1.e-15)
                                   I_m[i_solver2][i_out_level][i_mu][i_out_phi][i_stokes] = 0.;

                              save_result_solver_values(test_result, i_solver2, I_p[i_solver2][i_out_level][i_mu][i_out_phi][i_stokes], I_m[i_solver2][i_out_level][i_mu][i_out_phi][i_stokes]);

                              if (td->core_check_diffs) {
                                   d1 = 0.;
                                   if (fabs(I_p[0][i_out_level][i_mu][i_out_phi][i_stokes]) > DBL_EPSILON)
                                        d1 = (I_p[0][i_out_level][i_mu][i_out_phi][i_stokes] - I_p[i_solver2][i_out_level][i_mu][i_out_phi][i_stokes]) / I_p[0][i_out_level][i_mu][i_out_phi][i_stokes];
                                   d2 = 0.;
                                   if (fabs(I_m[0][i_out_level][i_mu][i_out_phi][i_stokes]) > DBL_EPSILON)
                                        d2 = (I_m[0][i_out_level][i_mu][i_out_phi][i_stokes] - I_m[i_solver2][i_out_level][i_mu][i_out_phi][i_stokes]) / I_m[0][i_out_level][i_mu][i_out_phi][i_stokes];

                                   if (! i_ptr_index || ! (ignore_mask_index_solver[i_ptr_index - ignore_list_index_solver] & solvers[i_solver2])) {
                                        if ((! (ignore_solver_mask_ref & solvers[i_solver2]) && fabs(d1) > tol[i_solver2]) || (! (ignore_solver_mask_tran & solvers[i_solver2]) && fabs(d2) > tol[i_solver2])) {
                                             fail_mask |= solvers[i_solver2];
                                             flag = 1;
                                        }
                                   }

                                   save_result_solver_deltas(test_result, i_solver2, d1, d2);
                              }
                         }
                    }
               }
          }


          /*--------------------------------------------------------------------
           *
           *------------------------------------------------------------------*/
          for (i_deriv = 0; i_deriv < n_derivs; ++i_deriv) {
               for (i_mu = 0; i_mu < n_mus; ++i_mu) {
                    out_theta = acos(qx[i_mu_offset + i_mu]) * R2D;

                    for (i_out_phi = 0; i_out_phi < n_out_phis; ++i_out_phi) {

                         for (i_stokes = 0; i_stokes < n_stokes; ++i_stokes) {
                              test_result = malloc(sizeof(test_result_data));
                              test_result_init(test_result, n_solvers);
                              list_append(test_result_list, test_result, 0);

                              save_result_header(test_result, index, index2, i_out_level, i_deriv, out_theta, out_phis[i_out_phi], i_stokes);

                              for (i_solver = 0; i_solver < n_solvers_all; ++i_solver) {
                                   for (i_solver2 = 0; i_solver2 < n_solvers; ++i_solver2) {
                                        if (xrtm_solver_index_to_mask(i_solver) == solvers[i_solver2])
                                             break;
                                   }

                                   if (i_solver2 == n_solvers)
                                        continue;

                                   save_result_solver_mask(test_result, i_solver2, solvers[i_solver2]);
/*
                                   if (stokes_element_zero(i_stokes, theta_0, out_theta, out_phis[i_out_phi] - phi_0)) {
                                        I_p_l[i_solver2][i_out_level][i_deriv][i_mu][i_out_phi][i_stokes] = 0.;
                                        I_m_l[i_solver2][i_out_level][i_deriv][i_mu][i_out_phi][i_stokes] = 0.;
                                   }
*/
                                   save_result_solver_values(test_result, i_solver2, I_p_l[i_solver2][i_out_level][i_deriv][i_mu][i_out_phi][i_stokes], I_m_l[i_solver2][i_out_level][i_deriv][i_mu][i_out_phi][i_stokes]);

                                   if (td->core_check_diffs) {
                                        d1 = 0.;
                                        if (fabs(I_p_l[0][i_out_level][i_deriv][i_mu][i_out_phi][i_stokes]) > DBL_EPSILON)
                                             d1 = (I_p_l[0][i_out_level][i_deriv][i_mu][i_out_phi][i_stokes] - I_p_l[i_solver2][i_out_level][i_deriv][i_mu][i_out_phi][i_stokes]) / I_p_l[0][i_out_level][i_deriv][i_mu][i_out_phi][i_stokes];
                                        d2 = 0.;
                                        if (fabs(I_m_l[0][i_out_level][i_deriv][i_mu][i_out_phi][i_stokes]) > DBL_EPSILON)
                                             d2 = (I_m_l[0][i_out_level][i_deriv][i_mu][i_out_phi][i_stokes] - I_m_l[i_solver2][i_out_level][i_deriv][i_mu][i_out_phi][i_stokes]) / I_m_l[0][i_out_level][i_deriv][i_mu][i_out_phi][i_stokes];

                                        if (n_ignore_deriv > 0)
                                             i_ptr_deriv = (int *) bsearch(&l, ignore_list_deriv, n_ignore_deriv, sizeof(int), (int (*)(const void *, const void *)) cmp_p_i2);

                                        if (! i_ptr_index || ! (ignore_mask_index_solver[i_ptr_index - ignore_list_index_solver] & solvers[i_solver2])) {
                                             if ((! (ignore_solver_mask_ref & solvers[i_solver2]) && fabs(d1) > tol_l[i_solver2]) || (! (ignore_solver_mask_tran & solvers[i_solver2]) && fabs(d2) > tol_l[i_solver2])) {
                                                  if (! i_ptr_deriv || ! (ignore_mask_deriv[i_ptr_deriv - ignore_list_deriv] & solvers[i_solver2])) {
                                                       fail_mask |= solvers[i_solver2];
                                                       flag = 1;
                                                  }
                                             }
                                        }

                                        save_result_solver_deltas(test_result, i_solver2, d1, d2);
                                   }
                              }
                         }
                    }
               }
          }
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (flag) {
          if (td->core_on_failed_diff_write_xrtm_cmd || td->core_on_failed_diff_stop_with_error)
               fprintf(stderr, "failed: index = %06d, phi = %e, solvers = %s\n", index, out_phis[i_out_phi], xrtm_solver_mask_to_name_list(fail_mask, temp, 1024));

          if (td->core_on_failed_diff_write_xrtm_cmd) {
               if (get_callxrtm_cmd(gd, md, n_out_phis, out_phis, temp2, MAX_COMMAND_LENGTH)) {
                    fprintf(stderr, "ERROR: write_callxrtm_cmd()\n");
                    return -1;
               }
               fprintf(stderr, "%s\n", temp2);

               if (xrtm_fwrite_input_fn(gd, NULL, NULL, "failed_input.txt", '\n', 1, 0, 0, 1)) {
                    fprintf(stderr, "ERROR: xrtm_fwrite_input_fn(): %s\n", "test_input.txt");
                    return -1;
               }
          }

          if (td->core_on_failed_diff_stop_with_error)
               return 1;
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     free_array1_d(qx);
     free_array1_d(qw);

     free_array1(out_phis2);

     free_array5_d(I_p);
     free_array5_d(I_m);

     if (options & XRTM_OPTION_CALC_DERIVS) {
          free_array6_d(I_p_l);
          free_array6_d(I_m_l);
     }


     return 0;
}
