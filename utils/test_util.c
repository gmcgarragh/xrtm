/******************************************************************************%
**
**    Copyright (C) 2007-2012 Greg McGarragh <gregm@atmos.colostate.edu>
**
**    This source code is licensed under the GNU General Public License (GPL),
**    Version 3.  See the file COPYING for more details.
**
*******************************************************************************/

#include <gutil.h>

#include <xrtm_support.h>

#include "test.h"
#include "test_util.h"
#include "test_macros.h"

#ifdef HAVE_PTHREADS_LIBRARY
#include <pthread.h>
#endif


/*******************************************************************************
 *
 ******************************************************************************/
int make_solvers(int solvers_mask1, int *solvers_mask2, enum xrtm_solver_mask *solvers_array2, double *tolerance, ...) {

     int i;
     int ii;

     enum xrtm_solver_mask solver;

     double a;

     va_list ap;

     va_start(ap, tolerance);

    *solvers_mask2 = 0;

     ii = 0;
     for (i = 0; (solver = (enum xrtm_solver_mask) va_arg(ap, int)) != 0; ++i) {
          a = va_arg(ap, double);
          if (solvers_mask1 & solver) {
              *solvers_mask2 |= solvers_mask1 & solver;
               if (solvers_array2)
                    solvers_array2[ii] = solver;
               tolerance[ii] = a;
               ii++;
          }
     }

     va_end(ap);

     return ii;
}



int make_solvers2(int solvers_mask1, int *solvers_mask2, enum xrtm_solver_mask *solvers_array2, double *tolerance_ref, double *tolerance_tran, ...) {

     int i;
     int ii;

     enum xrtm_solver_mask solver;

     double a;
     double b;

     va_list ap;

     va_start(ap, tolerance_tran);

    *solvers_mask2 = 0;

     ii = 0;
     for (i = 0; (solver = (enum xrtm_solver_mask) va_arg(ap, int)) != 0; ++i) {
          a = va_arg(ap, double);
          b = va_arg(ap, double);
          if (solvers_mask1 & solver) {
              *solvers_mask2 |= solvers_mask1 & solver;
               if (solvers_array2)
                    solvers_array2[ii] = solver;
               tolerance_ref [ii] = a;
               tolerance_tran[ii] = b;
               ii++;
          }
     }

     va_end(ap);

     return ii;
}



/*******************************************************************************
 *
 ******************************************************************************/
int make_tolerance(double *tol, ...) {

     int i;

     va_list ap;

     va_start(ap, tol);

     for (i = 0; (tol[i] = va_arg(ap, double)) != 0.; ++i) ;

     va_end(ap);

     return i;
}


/*******************************************************************************
 *
 ******************************************************************************/
static void set_chi_file(misc_data *md, test_data *td, int i_layer, double **chi) {

     int i;

     for (i = 0; i < MAX_PFS; ++i) {
          if (chi == td->gc[i])
               break;
     }

     md->coef_files[i_layer] = strdup(td->gf[i]);
}



static void set_chi_file_derivs(misc_data *md, test_data *td, int i_layer, int n_derivs, double ***chi_l) {

     int i;
     int j;

     for (i = 0; i < n_derivs; ++i) {
          for (j = 0; j < MAX_PFS; ++j) {
               if (chi_l[i] == td->gc_l[j])
                    break;
          }

          if (j == MAX_PFS)
               md->coef_files_l[i_layer][i] = strdup("none");
          else
               md->coef_files_l[i_layer][i] = strdup(td->gf_l[j]);
     }
}



static void set_chi_file_layers(misc_data *md, test_data *td, int n_layers, double ***chi) {

     int i;

     for (i = 0; i < n_layers; ++i)
          set_chi_file(md, td, i, chi[i]);

}



#ifdef HAVE_PTHREADS_LIBRARY

/*******************************************************************************
 *
 ******************************************************************************/
typedef struct {
     double theta_0;
     double *theta;
     int *n_chi;
     double ***chi;
     void *chi_l;
     double *omega;
     void *omega_l;
     double *ltau;
     void *ltau_l;
     double albedo;
     void *albedo_l;
     int ignore_mask_solver;
} test_case_data;


typedef struct {
     xrtm_data *gd;
     test_data *td;
     misc_data *md;

     int n_solvers;
     enum xrtm_solver_mask *solvers;

     int n_layers;
     int n_derivs;

     int *i_case;
     int  n_cases;
     test_case_data *tcd;

     int n_phi;
     double *phi;

     double *tol;
     double *tol_l;

     int n_ignore_index_solver;
     int *ignore_list_index_solver;
     int *ignore_mask_index_solver;

     pthread_mutex_t *mutex;
     pthread_mutex_t *mutex2;
} bounds_test_thread_data;


/*******************************************************************************
 *
 ******************************************************************************/
static int bounds_test_thread_loop(void *v) {

     int i;
     int j;
     int k;

     int flag;

     int index;

     bounds_test_thread_data *d;

     misc_data misc;

     d = (bounds_test_thread_data *) v;

     misc_init(d->gd, &misc, -1);

     flag  = 0;

     while (1) {
          pthread_mutex_lock(&d->mutex[0]);

          if (*d->i_case >= d->n_cases)
               flag = 1;
          else
               i = (*d->i_case)++;

          index = d->td->index;

          if (! flag)
               d->td->index += d->n_phi;

          pthread_mutex_unlock(&d->mutex[0]);

          if (flag)
               break;
/*
          printf("%d %d\n", d->n_cases, i);
*/
          XRTM_SET_THETA_0(d->gd, d->tcd[i].theta_0);

          XRTM_SET_OUT_THETAS (d->gd, d->tcd[i].theta);

          XRTM_SET_COEF_N (d->gd, d->tcd[i].n_chi, d->tcd[i].chi);
          set_chi_file_layers(&misc, d->td, d->n_layers, d->tcd[i].chi);

          XRTM_SET_OMEGA_N(d->gd, d->tcd[i].omega);

          XRTM_SET_LTAU_N (d->gd, d->tcd[i].ltau);

          XRTM_SET_KERNEL_AMPFAC(d->gd, 0, d->tcd[i].albedo);

          if (d->n_derivs > 0) {
               for (j = 0; j < d->n_layers; ++j) {
                    for (k = 0; k < d->n_derivs; ++k) {
                         if (((double ****) d->tcd[i].chi_l)[j][k])
                              XRTM_SET_COEF_L_11(d->gd, j, k, ((double ****) d->tcd[i].chi_l)[j][k]);
                    }

                    set_chi_file_derivs(&misc, d->td, j, d->n_derivs, ((double ****) d->tcd[i].chi_l)[j]);

                    XRTM_SET_OMEGA_L_1N(d->gd, j, ((double **) d->tcd[i].omega_l)[j]);
                    XRTM_SET_LTAU_L_1N (d->gd, j, ((double **) d->tcd[i].ltau_l )[j]);
               }

               XRTM_SET_KERNEL_AMPFAC_L_N(d->gd, 0, (double *) d->tcd[i].albedo_l);

               xrtm_update_varied_layers(d->gd);
          }

          if (test_execute(d->gd, &misc, d->td, index, d->n_solvers, d->solvers, d->n_phi, d->phi, d->tol, d->tol_l, d->n_ignore_index_solver, d->ignore_list_index_solver, d->ignore_mask_index_solver, d->tcd[i].ignore_mask_solver, d->tcd[i].ignore_mask_solver, 0, NULL, 0, d->mutex2)) {
               eprintf("ERROR: test_execute()\n");
               return 0;
/*
               pthread_exit(NULL);
*/
          }
     }

     misc_free(&misc);

     return -1;
/*
     pthread_exit(v);
*/
}



/*******************************************************************************
 *
 ******************************************************************************/
static int bounds_test_threads(xrtm_data *gd, test_data *td, int n_solvers, enum xrtm_solver_mask *solvers, int n_layers, int n_derivs, int n_cases, test_case_data *tcd, int n_phi, double *phi, double *tol, double *tol_l, int n_ignore_index_solver, int *ignore_list_index_solver, int *ignore_mask_index_solver, int n_threads) {

     int i;

     int i_case = 0;

     void *t_result;

     pthread_t *threads;

     pthread_mutex_t mutex;
     pthread_mutex_t mutex2;

     bounds_test_thread_data *d;


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     threads = (pthread_t *) malloc(n_threads * sizeof(pthread_t));
     if (threads == NULL) {
          eprintf("ERROR: memory allocation failed\n");
          exit(1);
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     pthread_mutex_init(&mutex, NULL);
     pthread_mutex_init(&mutex2, NULL);


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     d = (bounds_test_thread_data *) malloc(td->n_threads * sizeof(bounds_test_thread_data));

     for (i = 0; i < n_threads; ++i) {
          d[i].gd                       = &gd[i];
          d[i].td                       = td;
          d[i].n_solvers                = n_solvers;
          d[i].solvers                  = solvers;
          d[i].n_layers                 = n_layers;
          d[i].n_derivs                 = n_derivs;
          d[i].i_case                   = &i_case;
          d[i].n_cases                  = n_cases;
          d[i].tcd                      = tcd;
          d[i].n_phi                    = n_phi;
          d[i].phi                      = phi;
          d[i].tol                      = tol;
          d[i].tol_l                    = tol_l;
          d[i].n_ignore_index_solver    = n_ignore_index_solver;
          d[i].ignore_list_index_solver = ignore_list_index_solver;
          d[i].ignore_mask_index_solver = ignore_mask_index_solver;
          d[i].mutex                    = &mutex;
          d[i].mutex2                   = &mutex2;
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     for (i = 0; i < n_threads; ++i)
          pthread_create(&threads[i], NULL, (void *(*)(void *)) bounds_test_thread_loop, &d[i]);

     for (i = 0; i < n_threads; ++i) {
          pthread_join  ( threads[i], &t_result);
          if (t_result == NULL) {
               eprintf("ERROR: pthread_join()\n");
               return -1;
          }
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     free(threads);

     pthread_mutex_destroy(&mutex);
     pthread_mutex_destroy(&mutex2);

     free(d);


     return 0;
}

#endif

/*******************************************************************************
 *
 ******************************************************************************/
int bounds_test(xrtm_data *gd, test_data *td,
                int n_solvers, enum xrtm_solver_mask *solvers,
                int n_layers, int n_thetas,
                int n_theta_0_bounds, double *theta_0,
                int n_theta_bounds, double *theta,
                int n_chi_bounds, int *n_chi, double ***chi,
                int n_omega_bounds, double *omega,
                int n_ltau_bounds, double *ltau,
                int n_albedo_bounds, double *albedo,
                int n_phi_bounds, double *phi,
                double *tol, double *tol_l,
                int n_ignore_index_solver,
                int *ignore_list_index_solver,
                int *ignore_mask_index_solver,
                int (*ignore_mask_solver)(double theta_0, double *theta)) {

     int i;
     int j;
     int k;
     int l;
     int m;
     int n;

     int ignore_mask_solver2;

     misc_data misc;
#ifdef HAVE_PTHREADS_LIBRARY
     int n_cases;

     int i_case;

     test_case_data *test_cases;

     n_cases = n_theta_0_bounds * n_theta_bounds * n_chi_bounds * n_omega_bounds * n_ltau_bounds * n_albedo_bounds;

    if (td->n_threads == 1 || n_cases == 1) {
#endif
          misc_init(gd, &misc, -1);

          for (i = 0; i < n_theta_0_bounds; ++i) {
            XRTM_SET_THETA_0(gd, theta_0[i]);

            for (j = 0; j < n_theta_bounds; ++j) {
              XRTM_SET_OUT_THETAS(gd, &theta[j * n_thetas]);

              for (k = 0; k < n_chi_bounds; ++k) {
                XRTM_SET_COEF_N(gd, &n_chi[k * n_layers], &chi[k * n_layers]);

                set_chi_file_layers(&misc, td, n_layers, &chi[k * n_layers]);

                for (l = 0; l < n_omega_bounds; ++l) {
                  XRTM_SET_OMEGA_N(gd, &omega[l * n_layers]);

                  for (m = 0; m < n_ltau_bounds; ++m) {
                    XRTM_SET_LTAU_N(gd, &ltau[m * n_layers]);

                    for (n = 0; n < n_albedo_bounds; ++n) {
                      XRTM_SET_KERNEL_AMPFAC(gd, 0, albedo[n]);
/*
                      printf("%d %d,  %d %d,  %d %d,  %d %d,  %d %d,  %d %d\n", n_theta_0_bounds, i, n_theta_bounds, j, n_chi_bounds, k, n_omega_bounds, l, n_ltau_bounds, m, n_albedo_bounds, n); 
*/
                      ignore_mask_solver2 = ignore_mask_solver(theta_0[i], &theta[j * n_thetas]);

                      HANDLE_RETURN(test_execute(gd, &misc, td, td->index, n_solvers, solvers, n_phi_bounds, phi, tol, tol_l, n_ignore_index_solver, ignore_list_index_solver, ignore_mask_index_solver, ignore_mask_solver2, ignore_mask_solver2, 0, NULL, 0, NULL), "test_execute()");
                      td->index += n_phi_bounds;

                      if (td->quick_run)
                           goto L1;
                    }
                  }
                }
              }
            }
          }

L1:  misc_free(&misc);
#ifdef HAVE_PTHREADS_LIBRARY
    }
    else {
          test_cases = (test_case_data *) malloc(n_cases * sizeof(test_case_data));

          i_case = 0;
          for (i = 0; i < n_theta_0_bounds; ++i) {
            for (j = 0; j < n_theta_bounds; ++j) {
              for (k = 0; k < n_chi_bounds; ++k) {
                for (l = 0; l < n_omega_bounds; ++l) {
                  for (m = 0; m < n_ltau_bounds; ++m) {
                    for (n = 0; n < n_albedo_bounds; ++n) {
                      test_cases[i_case].theta_0 = theta_0[i];
                      test_cases[i_case].theta   = &theta  [j * n_thetas];
                      test_cases[i_case].n_chi   = &n_chi  [k * n_layers];
                      test_cases[i_case].chi     = &chi    [k * n_layers];
                      test_cases[i_case].omega   = &omega  [l * n_layers];
                      test_cases[i_case].ltau    = &ltau   [m * n_layers];
                      test_cases[i_case].albedo  = albedo [n];

                      test_cases[i_case].ignore_mask_solver = ignore_mask_solver(theta_0[i], &theta[j * n_thetas]);

                      i_case++;
                    }
                  }
                }
              }
            }
          }

          if (bounds_test_threads(gd, td, n_solvers, solvers, n_layers, 0, n_cases, test_cases, n_phi_bounds, phi, tol, tol_l, n_ignore_index_solver, ignore_list_index_solver, ignore_mask_index_solver, td->n_threads)) {
               eprintf("ERROR: bounds_test_threads()\n");
               return -1;
          }

          free(test_cases);
     }
#endif
     return 0;
}



/*******************************************************************************
 *
 ******************************************************************************/
int bounds_test_derivs(xrtm_data *gd, test_data *td,
                       int n_solvers, enum xrtm_solver_mask *solvers,
                       int n_layers, int n_derivs, int n_thetas,
                       int n_theta_0_bounds, double *theta_0,
                       int n_theta_bounds, double *theta,
                       int n_chi_bounds, int *n_chi, double ***chi,
                       int n_omega_bounds, double *omega,
                       int n_ltau_bounds, double *ltau,
                       int n_albedo_bounds, double *albedo,
                       int n_phi_bounds, double *phi,
                       int n_deriv_bounds,
                       double ***chi_l, double *omega_l, double *ltau_l,
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

     int ignore_mask_solver2;

     misc_data misc;
#ifdef HAVE_PTHREADS_LIBRARY
     int n_cases;

     int i_case;

     test_case_data *test_cases;

     n_cases = n_theta_0_bounds * n_theta_bounds * n_chi_bounds * n_omega_bounds * n_ltau_bounds * n_albedo_bounds;

     if (td->n_threads == 1 || n_cases == 1) {
#endif
          misc_init(gd, &misc, -1);

          for (i = 0; i < n_theta_0_bounds; ++i) {
            XRTM_SET_THETA_0(gd, theta_0[i]);

            for (j = 0; j < n_theta_bounds; ++j) {
              XRTM_SET_OUT_THETAS(gd, &theta[j * n_thetas]);

              for (k = 0; k < n_chi_bounds; ++k) {
                XRTM_SET_COEF_N(gd, &n_chi[k * n_layers], &chi[k * n_layers]);

                set_chi_file_layers(&misc, td, n_layers, &chi[k * n_layers]);

                for (l = 0; l < n_omega_bounds; ++l) {
                  XRTM_SET_OMEGA_N(gd, &omega[l * n_layers]);

                  for (m = 0; m < n_ltau_bounds; ++m) {
                    XRTM_SET_LTAU_N(gd, &ltau[m * n_layers]);

                    for (n = 0; n < n_albedo_bounds; ++n) {
                      XRTM_SET_KERNEL_AMPFAC(gd, 0, albedo[n]);

                      for (o = 0; o < n_deriv_bounds; ++o) {
/*
                        printf("%d %d,  %d %d,  %d %d,  %d %d,  %d %d,  %d %d,  %d %d\n", n_theta_0_bounds, i, n_theta_bounds, j, n_chi_bounds, k, n_omega_bounds, l, n_ltau_bounds, m, n_albedo_bounds, n, n_deriv_bounds, o); 
*/
                        for (p = 0; p < n_layers; ++p) {
                             for (q = 0; q < n_derivs; ++q) {
                                  if (chi_l[k * n_deriv_bounds * n_layers * n_derivs + o * n_layers * n_derivs + p * n_derivs + q])
                                       XRTM_SET_COEF_L_11(gd, p, q, chi_l[k * n_deriv_bounds * n_layers * n_derivs + o * n_layers * n_derivs + p * n_derivs + q]);
                             }

                             set_chi_file_derivs(&misc, td, p, n_derivs, &chi_l[k * n_deriv_bounds * n_layers * n_derivs + o * n_layers * n_derivs + p * n_derivs]);

                             XRTM_SET_OMEGA_L_1N(gd, p, &omega_l[o * n_layers * n_derivs + p * n_derivs]);
                             XRTM_SET_LTAU_L_1N (gd, p, &ltau_l [o * n_layers * n_derivs + p * n_derivs]);
                        }

                        XRTM_SET_KERNEL_AMPFAC_L_N(gd, 0, &albedo_l[o * n_derivs]);

                        xrtm_update_varied_layers(gd);

                        ignore_mask_solver2 = ignore_mask_solver(theta_0[i], &theta[j * n_thetas]);

                        HANDLE_RETURN(test_execute(gd, &misc, td, td->index, n_solvers, solvers, n_phi_bounds, phi, tol, tol_l, n_ignore_index_solver, ignore_list_index_solver, ignore_mask_index_solver, ignore_mask_solver2, ignore_mask_solver2, n_ignore_deriv, ignore_list_deriv, ignore_mask_deriv, NULL), "test_execute()");
                        td->index += n_phi_bounds;

                        if (td->quick_run)
                             goto L1;
                      }
                    }
                  }
                }
              }
            }
          }

L1:  misc_free(&misc);
#ifdef HAVE_PTHREADS_LIBRARY
     }
     else {
          test_cases = (test_case_data *) malloc(n_cases * sizeof(test_case_data));

          i_case = 0;
          for (i = 0; i < n_theta_0_bounds; ++i) {
            for (j = 0; j < n_theta_bounds; ++j) {
              for (k = 0; k < n_chi_bounds; ++k) {
                for (l = 0; l < n_omega_bounds; ++l) {
                  for (m = 0; m < n_ltau_bounds; ++m) {
                    for (n = 0; n < n_albedo_bounds; ++n) {
                      for (o = 0; o < n_deriv_bounds; ++o) {
                           test_cases[i_case].theta_0 = theta_0[i];
                           test_cases[i_case].theta   = &theta [j * n_thetas];
                           test_cases[i_case].n_chi   = &n_chi [k * n_layers];
                           test_cases[i_case].chi     = &chi   [k * n_layers];
                           test_cases[i_case].omega   = &omega [l * n_layers];
                           test_cases[i_case].ltau    = &ltau  [m * n_layers];
                           test_cases[i_case].albedo  = albedo [n];

                           test_cases[i_case].chi_l     = &chi_l   [k][o];
                           test_cases[i_case].omega_l   = &omega_l [o * n_layers * n_derivs];
                           test_cases[i_case].ltau_l    = &ltau_l  [o * n_layers * n_derivs];
                           test_cases[i_case].albedo_l  = &albedo_l[o * n_derivs];
                      }

                      test_cases[i_case].ignore_mask_solver = ignore_mask_solver(theta_0[i], &theta[j * n_thetas]);

                      i_case++;
                    }
                  }
                }
              }
            }
          }

          if (bounds_test_threads(gd, td, n_solvers, solvers, n_layers, n_derivs, n_cases, test_cases, n_phi_bounds, phi, tol, tol_l, n_ignore_index_solver, ignore_list_index_solver, ignore_mask_index_solver, td->n_threads)) {
               eprintf("ERROR: bounds_test_threads()\n");
               return -1;
          }

          free(test_cases);
     }
#endif
     return 0;
}
