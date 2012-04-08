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

#include "input_util.h"

#include "test.h"

#ifdef HAVE_PTHREADS_LIBRARY
#include <pthread.h>
#endif


/*******************************************************************************
 *
 ******************************************************************************/
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


/*
typedef struct {
     
} results_data;



static int write_results(FILE *fp, results_data *d) {

     return 0;
}
*/


static int write_callxrtm_cmd(FILE *fp, xrtm_data *gd, misc_data *md, double phi) {

     fprintf(fp, "callxrtm -input_string \"");

     if (xrtm_fwrite_input_fp(gd, NULL, md, fp, ' ', 0, 1, 1, 0)) {
          eprintf("ERROR: xrtm_fwrite_input(): %s\n", "test_input.txt");
          return -1;
     }

     fprintf(fp, "\" -out_phis %e\n", phi);

     return 0;
}



int test_execute(xrtm_data *gd, misc_data *md, test_data *td, int index, int n_solvers, enum xrtm_solver_mask *solvers, int n_phis, double *phis, double *tol, double *tol_l, int n_ignore_index_solver, int *ignore_list_index_solver, int *ignore_mask_index_solver, int ignore_solver_mask_ref, int ignore_solver_mask_tran, int n_ignore_deriv, int *ignore_list_deriv, int *ignore_mask_deriv, void *mutex) {

     char temp[1024];

     int i;
     int ii;
     int j;
     int k;
     int l;

     int flag;

     int options;

     int n_quad;
     int n_stokes;
     int n_derivs;
     int n_out_levels;
     int n_out_thetas;

     int n_quad_x;

     int n_mus2;

     int i_phis;

     int fail_mask;

     int *i_ptr_index = NULL;
     int *i_ptr_deriv = NULL;

     double theta_0;
     double phi_0;

     double d1;
     double d2;

     double **phis2;

     double *qx;
     double *qw;

     double *****I_p;
     double *****I_m;

     double *****K_p2;
     double *****K_m2;

     double ******K_p;
     double ******K_m;


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
          n_mus2 = n_quad;
     else
          n_mus2 = n_out_thetas;


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     qx = alloc_array1_d(n_quad_x);
     if (xrtm_qx(gd, qx)) {
          eprintf("ERROR: xrtm_qx()\n");
          return -1;
     }

     qw = alloc_array1_d(n_quad_x);
     if (xrtm_qw(gd, qw)) {
          eprintf("ERROR: xrtm_qw()\n");
          return -1;
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     phis2 = (double **) alloc_array1(n_mus2, sizeof(double *));

     for (i = 0; i < n_mus2; ++i)
          phis2[i] = phis;


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     I_p = alloc_array5_d(n_solvers, n_out_levels, n_mus2, n_phis, n_stokes);
     I_m = alloc_array5_d(n_solvers, n_out_levels, n_mus2, n_phis, n_stokes);

     if (options & XRTM_OPTION_CALC_DERIVS) {
          K_p = alloc_array6_d(n_solvers, n_out_levels, n_derivs, n_mus2, n_phis, n_stokes);
          K_m = alloc_array6_d(n_solvers, n_out_levels, n_derivs, n_mus2, n_phis, n_stokes);
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     for (i = 0; i < n_solvers; ++i) {
          if (options & XRTM_OPTION_CALC_DERIVS) {
               K_p2 = K_p[i];
               K_m2 = K_m[i];
          }

          if (xrtm_radiance(gd, solvers[i], n_phis, phis2, I_p[i], I_m[i], K_p2, K_m2)) {
               eprintf("ERROR: xrtm_radiance()\n");
               return -1;
          }
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     theta_0 = xrtm_get_theta_0(gd);
     phi_0   = xrtm_get_phi_0(gd);
#ifdef HAVE_PTHREADS_LIBRARY
     if (mutex)
          pthread_mutex_lock((pthread_mutex_t *) mutex);
#endif
     for (i_phis = 0; i_phis < n_phis; ++i_phis) {
          flag = 0;

          if (n_out_thetas == 0)
               ii = n_quad - 1;
          else
               ii = 0;

          fail_mask = 0;

          if (n_ignore_index_solver > 0)
               i_ptr_index = (int *) bsearch(&index, ignore_list_index_solver, n_ignore_index_solver, sizeof(int), (int (*)(const void *, const void *)) cmp_p_i2);

          for (i = 0; i < n_mus2; ++i) {
               for (j = 0; j < n_stokes; ++j) {

                    if (td->write_results)
                         fprintf(td->fp_results, "%06d, theta_0 = %9f, theta = %9f, phi = %10f, i_stokes = %d: ", index, theta_0, acos(qx[n_quad + ii])*R2D, phis[i_phis], j);

                    for (k = 0; k < n_solvers; ++k) {
                         if (stokes_element_zero(j, theta_0, acos(qx[n_quad + ii])*R2D, phis[i_phis] - phi_0)) {
                              I_p[k][0][ii][i_phis][j] = 0.;
                              I_m[k][1][ii][i_phis][j] = 0.;
                              for (l = 0; l < n_derivs; ++l) {
                                   K_p[k][0][l][ii][i_phis][j] = 0.;
                                   K_m[k][1][l][ii][i_phis][j] = 0.;
                              }
                         }

                         if (td->write_results) {
                              fprintf(td->fp_results, "   %s -> ", xrtm_solver_name2(solvers[k]));

                              fprintf(td->fp_results, " % e % e",
                                      I_p[k][0][ii][i_phis][j], I_m[k][1][ii][i_phis][j]);
                              for (l = 0; l < n_derivs; ++l) {
                                   fprintf(td->fp_results, " % e % e",
                                           K_p[k][0][l][ii][i_phis][j], K_m[k][1][l][ii][i_phis][j]);
                              }
                         }

                         if (td->check_diffs) {
                              if (td->write_results)
                                   fprintf(td->fp_results, "   %s delta -> ", xrtm_solver_name2(solvers[k]));

                              d1 = 0.;
                              if (fabs(I_p[0][0][ii][i_phis][j]) > DBL_EPSILON)
                                   d1 = (I_p[0][0][ii][i_phis][j] - I_p[k][0][ii][i_phis][j]) / I_p[0][0][ii][i_phis][j];
                              d2 = 0.;
                              if (fabs(I_m[0][1][ii][i_phis][j]) > DBL_EPSILON)
                                   d2 = (I_m[0][1][ii][i_phis][j] - I_m[k][1][ii][i_phis][j]) / I_m[0][1][ii][i_phis][j];

                              if (! i_ptr_index || ! (ignore_mask_index_solver[i_ptr_index - ignore_list_index_solver] & solvers[k])) {
                                   if ((! (ignore_solver_mask_ref & solvers[k]) && fabs(d1) > tol[k]) || (! (ignore_solver_mask_tran & solvers[k]) && fabs(d2) > tol[k])) {
                                        fail_mask |= solvers[k];
                                        flag = 1;
                                   }
                              }

                              if (td->write_results)
                                   fprintf(td->fp_results, " % e % e", d1, d2);

                              for (l = 0; l < n_derivs; ++l) {
                                   d1 = 0.;
                                   if (fabs(K_p[0][0][l][ii][i_phis][j]) > DBL_EPSILON)
                                        d1 = (K_p[0][0][l][ii][i_phis][j] - K_p[k][0][l][ii][i_phis][j]) / K_p[0][0][l][ii][i_phis][j];
                                   d2 = 0.;
                                   if (fabs(K_m[0][1][l][ii][i_phis][j]) > DBL_EPSILON)
                                        d2 = (K_m[0][1][l][ii][i_phis][j] - K_m[k][1][l][ii][i_phis][j]) / K_m[0][1][l][ii][i_phis][j];

                                   if (n_ignore_deriv > 0)
                                        i_ptr_deriv = (int *) bsearch(&l, ignore_list_deriv, n_ignore_deriv, sizeof(int), (int (*)(const void *, const void *)) cmp_p_i2);

                                   if (! i_ptr_index || ! (ignore_mask_index_solver[i_ptr_index - ignore_list_index_solver] & solvers[k])) {
                                        if ((! (ignore_solver_mask_ref & solvers[k]) && fabs(d1) > tol_l[k]) || (! (ignore_solver_mask_tran & solvers[k]) && fabs(d2) > tol_l[k])) {
                                             if (! i_ptr_deriv || ! (ignore_mask_deriv[i_ptr_deriv - ignore_list_deriv] & solvers[k])) {
                                                  fail_mask |= solvers[k];
                                                  flag = 1;
                                             }
                                        }
                                   }

                                   if (td->write_results)
                                        fprintf(td->fp_results, " % e % e", d1, d2);
                              }
                         }
                    }

                    if (td->write_results)
                         fprintf(td->fp_results, "\n");
               }

               if (n_out_thetas == 0)
                    ii--;
               else
                    ii++;
          }


          /*--------------------------------------------------------------------
           *
           *------------------------------------------------------------------*/
          fprintf(td->fp_xrtm_cmd, "%06d: ", index);

          if (write_callxrtm_cmd(td->fp_xrtm_cmd, gd, md, phis[i_phis])) {
               eprintf("ERROR: write_callxrtm_cmd()\n");
               return -1;
          }


          /*--------------------------------------------------------------------
           *
           *------------------------------------------------------------------*/
          if (flag) {
               if (td->on_failed_diff_write_xrtm_cmd || td->on_failed_diff_stop)
                    fprintf(stderr, "failed: index = %06d, theta_0 = %e, phi = %e, solvers = %s\n", index, theta_0, phis[i_phis], xrtm_solver_list(fail_mask, temp));

               if (td->on_failed_diff_write_xrtm_cmd) {
                    if (write_callxrtm_cmd(stderr, gd, md, phis[i_phis])) {
                         eprintf("ERROR: write_callxrtm_cmd()\n");
                         return -1;
                    }
/*
                    if (xrtm_fwrite_input_fn(gd, NULL, NULL, "failed_input.txt", '\n', 1, 0, 0, 1)) {
                         eprintf("ERROR: xrtm_fwrite_input_fn(): %s\n", "test_input.txt");
                         exit(1);
                    }
*/
               }

               if (td->on_failed_diff_stop)
                    return 1;
          }


          index++;
     }
#ifdef HAVE_PTHREADS_LIBRARY
     if (mutex)
          pthread_mutex_unlock((pthread_mutex_t *) mutex);
#endif

     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     free_array1_d(qx);
     free_array1_d(qw);

     free_array1(phis2);

     free_array5_d(I_p);
     free_array5_d(I_m);

     if (options & XRTM_OPTION_CALC_DERIVS) {
          free_array6_d(K_p);
          free_array6_d(K_m);
     }


     return 0;
}
