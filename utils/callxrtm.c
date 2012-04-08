/******************************************************************************%
**
**    Copyright (C) 2007-2012 Greg McGarragh <gregm@atmos.colostate.edu>
**
**    This source code is licensed under the GNU General Public License (GPL),
**    Version 3.  See the file COPYING for more details.
**
*******************************************************************************/

#include <gutil.h>

#include <rtutil_support.h>

#include <xrtm_fd.h>
#include <xrtm_interface.h>
#include <xrtm_support.h>

#include "input_util.h"


#ifdef NAME
#undef NAME
#endif
#define NAME    "callxrtm"


typedef struct {
     int derivs;
     int help;
     int input_file;
     int input_string;
     int radiance;
     int radiance_mean;
     int flux;
     int flux_divergence;
     int output_results;
     int out_phis;
     int timing;
     int version;
} options_data;


void usage();
void version();


int main(int argc, char *argv[]) {

     char *input;
     char *input2;
     char *input3;

     const char *defaults  = "-n_quad 16 -n_stokes 1 -n_derivs 0 -n_layers 1 -n_kernels 0 -n_kernel_quad 16 -n_out_levels 1 -n_out_thetas 0";

     const char *defaults2 = "-doub_d_tau 1.e-9 -pade_params 4,5 -sos_params 128,.005,1.e-9 -fourier_tol 0. -F_0 1. -theta_0 0. -phi_0 0.";

     int i;
     int ii;
     int j;
     int k;
     int l;
     int m;
     int n;
     int n1;
     int n2;
     int n21;
     int n22;

     int mask;
     int mask2;

     int options;
     int solver;
     int solutions;

     int n_quad;
     int n_stokes;
     int n_derivs;
     int n_out_levels;
     int n_out_thetas;

     int i_out_thetas2;
     int n_out_thetas2;

     int n_quad_d;

     int n_out_phis;

     int i_timing_loop;
     int n_timing_loop;

     int output_up;
     int output_down;

     int *out_levels;

     double dzero = 0.;

     double *out_phis;
     double **out_phis2;

     double *qx;
     double *qw;

     double *out_taus;

     double ****I_p;
     double ****I_m;

     double *****K_p;
     double *****K_m;

     double *mean_p;
     double *mean_m;
     double **mean_p_l;
     double **mean_m_l;

     double *flux_p;
     double *flux_m;
     double **flux_p_l;
     double **flux_m_l;

     double *flux_div;
     double **flux_div_l;

     clock_t time1;
     clock_t time2;

     xrtm_data d;
     xrtm_fd_data fd;
     misc_data md;

     options_data o;

#if PLATFORM == WIN32_MSVC
     _set_output_format(_TWO_DIGIT_EXPONENT);
#endif

     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     o.derivs          = 0;
     o.help            = 0;
     o.input_file      = 0;
     o.input_string    = 0;
     o.output_results  = 1;
     o.out_phis        = 0;
     o.radiance        = 1;
     o.radiance_mean   = 0;
     o.flux            = 0;
     o.flux_divergence = 0;
     o.timing          = 0;
     o.version         = 0;

     solver        = 0;

     n_out_phis    = 1;
     out_phis      = &dzero;

     n_timing_loop = 1;

     n = 0;
     for (i = 1; i < argc; ++i) {
          if (argv[i][0] == '-') {
               if (strcmp(argv[i], "-help") == 0) {
                    usage();
                    exit(0);
               }
               else if (strcmp(argv[i], "-input_file") == 0) {
                    check_arg_count(i, argc, 1, argv[i]);
                    o.input_file   = 1;
                    o.input_string = 0;
                    input = argv[++i];
               }
               else if (strcmp(argv[i], "-input_string") == 0) {
                    check_arg_count(i, argc, 1, argv[i]);
                    o.input_file   = 0;
                    o.input_string = 1;
                    input = argv[++i];
               }
               else if (strcmp(argv[i], "-output_results") == 0)
                    o.output_results = 1;
               else if (strcmp(argv[i], "-no_output_results") == 0)
                    o.output_results = 0;
               else if (strcmp(argv[i], "-out_phis") == 0) {
                    check_arg_count(i, argc, 1, argv[i]);
                    o.out_phis = 1;
                    n_out_phis = get_value_list_d(argv[i], argv[i+1], &out_phis);
                    i++;
               }
               else if (strcmp(argv[i], "-radiance") == 0)
                    o.radiance = 1;
               else if (strcmp(argv[i], "-no_radiance") == 0)
                    o.radiance = 0;
               else if (strcmp(argv[i], "-radiance_mean") == 0)
                    o.radiance_mean = 1;
               else if (strcmp(argv[i], "-flux") == 0)
                    o.flux = 1;
               else if (strcmp(argv[i], "-flux_divergence") == 0)
                    o.flux_divergence = 1;
               else if (strcmp(argv[i], "-solver") == 0) {
                    check_arg_count(i, argc, 1, argv[i]);
                    if ((solver = xrtm_solver_mask2(argv[++i])) == 0) {
                         printf("ERROR: xrtm_solver_mask2()\n");
                         exit(1);
                    }
               }
               else if (strcmp(argv[i], "-timing") == 0) {
                    ii = i;
                    check_arg_count(i, argc, 1, argv[ii]);
                    o.timing = 1;
                    n_timing_loop = strtoi_errmsg_exit(argv[++i], argv[ii]);
               }
               else if (strcmp(argv[i], "-version") == 0) {
                    version();
                    exit(0);
               }
               else {
                    printf("Invalid option: %s, use -help for more information\n", argv[i]);
                    exit(1);
               }
          }
          else {
/*
               if (n == 0)

               else {
*/
                    printf("Too many arguments, use -help for more information\n");
                    exit(1);
/*
               }
*/
               ++n;
          }
     }

     if (n < 0) {
          printf("Not enough arguments, use -help for more information\n");
          exit(1);
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (o.input_file) {
          if (xrtm_fread_input(&d, &fd, &md, input, 1, 0, 0, 1)) {
               eprintf("ERROR: xrtm_fread_input(): %s\n", input);
               exit(1);
          }
     }
     else
     if (o.input_string) {
          input2 = strstr(input, "::");
          if (! input2) {
               eprintf("ERROR: must have :: in input string\n");
               exit(1);
          }

          *input2 = '\0';
           input2 += 2;

          n1  = strlen(defaults);
          n2  = strlen(input);
          n21 = strlen(defaults2);
          n22 = strlen(input2);

          input3 = (char *) malloc(n1 + 1 + n2 + 1 + n21 + 1 + n22 + 1);

          ii = 0;
          for (i = 0; i < n1; ++i)
               input3[ii++] = defaults[i];

          input3[ii++] = ' ';

          for (i = 0; i < n2; ++i)
               input3[ii++] = input   [i];

          for (i = 0; i < n21; ++i)
               input3[ii++] = defaults2[i];

          input3[ii++] = ' ';

          for (i = 0; i < n22; ++i)
               input3[ii++] = input2   [i];

          input3[ii++] = '\0';

          if (xrtm_sread_input(&d, &fd, &md, input3, 0, 1, 1, 0)) {
               eprintf("ERROR: xrtm_sread_input()\n");
               exit(1);
          }

          free(input3);
     }
     else {
          if (xrtm_fread_input(&d, &fd, &md, "-", 1, 0, 0, 1)) {
               eprintf("ERROR: xrtm_fread_input()\n");
               exit(1);
          }
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     options      = xrtm_get_options(&d);

     o.derivs     = options & XRTM_OPTION_CALC_DERIVS || md.fd_method >= 0;

     if (solver == 0) {
          mask = xrtm_get_solvers(&d);

          n = xrtm_solvers_n();

          for (i = 0; i < n; ++i) {
               mask2 = xrtm_solver_mask(i);
               if (mask & mask2) {
                    solver = mask2;
                    break;
               }
          }
     }

     n_quad       = xrtm_get_n_quad(&d);
     n_stokes     = xrtm_get_n_stokes(&d);
     n_derivs     = 0;
     if (options & XRTM_OPTION_CALC_DERIVS)
          n_derivs     = xrtm_get_n_derivs(&d);
     else
     if (md.fd_method >= 0)
          n_derivs     = xrtm_fd_get_n_derivs(&fd);
     n_out_levels = xrtm_get_n_out_levels(&d);
     n_out_thetas = xrtm_get_n_out_thetas(&d);
     if (n_out_thetas == 0)
          n_out_thetas2 = n_quad;
     else
          n_out_thetas2 = n_out_thetas;

     n_quad_d     = n_quad   + n_out_thetas;

     if (o.radiance) {
          qx = alloc_array1_d(n_quad_d);
          xrtm_qx(&d, qx);

          qw = alloc_array1_d(n_quad_d);
          xrtm_qw(&d, qw);
     }

     if (options & XRTM_OPTION_OUTPUT_AT_LEVELS) {
          out_levels = alloc_array1_i(n_out_levels);
          xrtm_get_out_levels(&d, out_levels);
     }
     else {
          out_taus   = alloc_array1_d(n_out_levels);
          xrtm_get_out_taus  (&d, out_taus);
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (options & XRTM_OPTION_CALC_DERIVS) {
          if (xrtm_update_varied_layers(&d)) {
               eprintf("ERROR: xrtm_update_varied_layers()\n");
               exit(1);
          }
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     out_phis2 = (double **) alloc_array1(n_out_thetas2, sizeof(double *));

     for (i = 0; i < n_out_thetas2; ++i)
               out_phis2[i] = out_phis;


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (o.radiance) {
          I_p = alloc_array4_d(n_out_levels, n_out_thetas2, n_out_phis, n_stokes);
          I_m = alloc_array4_d(n_out_levels, n_out_thetas2, n_out_phis, n_stokes);

          if (o.derivs) {
               K_p = alloc_array5_d(n_out_levels, n_derivs, n_out_thetas2, n_out_phis, n_stokes);
               K_m = alloc_array5_d(n_out_levels, n_derivs, n_out_thetas2, n_out_phis, n_stokes);
          }
     }

     if (o.radiance_mean) {
          mean_p = alloc_array1_d(n_out_levels);
          mean_m = alloc_array1_d(n_out_levels);

          if (o.derivs) {
               mean_p_l = alloc_array2_d(n_out_levels, n_derivs);
               mean_m_l = alloc_array2_d(n_out_levels, n_derivs);
          }
     }

     if (o.flux) {
          flux_p = alloc_array1_d(n_out_levels);
          flux_m = alloc_array1_d(n_out_levels);

          if (o.derivs) {
               flux_p_l = alloc_array2_d(n_out_levels, n_derivs);
               flux_m_l = alloc_array2_d(n_out_levels, n_derivs);
          }
     }

     if (o.flux_divergence) {
          flux_div = alloc_array1_d(n_out_levels);

          if (o.derivs)
               flux_div_l = alloc_array2_d(n_out_levels, n_derivs);
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     solutions = 0;

     if (o.radiance)
          solutions |= XRTM_OUTPUT_RADIANCE;
     if (o.radiance_mean)
          solutions |= XRTM_OUTPUT_RADIANCE_MEAN;
     if (o.flux)
          solutions |= XRTM_OUTPUT_FLUX;
     if (o.flux_divergence)
          solutions |= XRTM_OUTPUT_FLUX_DIVERGENCE;

     if (md.fd_method < 0) {
/*
          if (1) {
*/
          if (! (options & XRTM_OPTION_REVERSE_DERIVS)) {
               if (xrtm_solution(&d, (enum xrtm_solver_mask) solver, solutions, n_out_phis, out_phis2, I_p, I_m, K_p, K_m, mean_p, mean_m, mean_p_l, mean_m_l, flux_p, flux_m, flux_p_l, flux_m_l, flux_div, flux_div_l)) {
                    eprintf("ERROR: xrtm_solution()\n");
                    exit(1);
               }
          }
          else {
               if (xrtm_solution_2(&d, (enum xrtm_solver_mask) solver, solutions, n_out_phis, out_phis2, I_p, I_m, K_p, K_m, mean_p, mean_m, mean_p_l, mean_m_l, flux_p, flux_m, flux_p_l, flux_m_l, flux_div, flux_div_l)) {
                    eprintf("ERROR: xrtm_solution_2()\n");
                    exit(1);
               }
          }
     }
     else {
          if (xrtm_fd_solution(&fd, (enum xrtm_solver_mask) solver, solutions, FD_METHOD_CENTRAL, n_out_phis, out_phis2, I_p, I_m, K_p, K_m, mean_p, mean_m, mean_p_l, mean_m_l, flux_p, flux_m, flux_p_l, flux_m_l, flux_div, flux_div_l)) {
               eprintf("ERROR: xrtm_fd_solution()\n");
               exit(1);
          }
     }

     if (o.timing) {
          time1 = clock();

          for (i_timing_loop = 0; i_timing_loop < n_timing_loop; ++i_timing_loop) {
               if (md.fd_method < 0) {
                    if (! (options & XRTM_OPTION_REVERSE_DERIVS)) {
                         if (xrtm_solution(&d, (enum xrtm_solver_mask) solver, solutions, n_out_phis, out_phis2, I_p, I_m, K_p, K_m, mean_p, mean_m, mean_p_l, mean_m_l, flux_p, flux_m, flux_p_l, flux_m_l, flux_div, flux_div_l)) {
                              eprintf("ERROR: xrtm_solution()\n");
                              exit(1);
                         }
                    }
                    else {
                         if (xrtm_solution_2(&d, (enum xrtm_solver_mask) solver, solutions, n_out_phis, out_phis2, I_p, I_m, K_p, K_m, mean_p, mean_m, mean_p_l, mean_m_l, flux_p, flux_m, flux_p_l, flux_m_l, flux_div, flux_div_l)) {
                              eprintf("ERROR: xrtm_solution_2()\n");
                              exit(1);
                         }
                    }
               }
               else {
                    if (xrtm_fd_solution(&fd, (enum xrtm_solver_mask) solver, solutions, FD_METHOD_CENTRAL, n_out_phis, out_phis2, I_p, I_m, K_p, K_m, mean_p, mean_m, mean_p_l, mean_m_l, flux_p, flux_m, flux_p_l, flux_m_l, flux_div, flux_div_l)) {
                         eprintf("ERROR: xrtm_fd_solution()\n");
                         exit(1);
                    }
               }
          }

          time2 = clock();
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (o.output_results) {
          output_up   = ! (d.options & (XRTM_OPTION_UPWELLING_OUTPUT | XRTM_OPTION_DOWNWELLING_OUTPUT)) || d.options & XRTM_OPTION_UPWELLING_OUTPUT;
          output_down = ! (d.options & (XRTM_OPTION_UPWELLING_OUTPUT | XRTM_OPTION_DOWNWELLING_OUTPUT)) || d.options & XRTM_OPTION_DOWNWELLING_OUTPUT;

          if (o.radiance) {
               if (n_out_thetas == 0) {
                    n_out_thetas2 = n_quad;
                    i_out_thetas2 = 0;
               }
               else {
                    n_out_thetas2 = n_out_thetas;
                    i_out_thetas2 = n_quad;
               }

               if (options & XRTM_OPTION_OUTPUT_AT_LEVELS)
                    printf("level         ");
               else
                    printf("tau           ");
               printf("  theta         phi           i_stokes  ");
               if (output_up)
                    printf("I_up          ");
               if (output_down)
                    printf("I_dn          ");
               for (i = 0; i < n_derivs; ++i) {
                    if (i < 10) {
                         if (output_up)
                              printf("L(I_up)[%d]    ", i);
                         if (output_down)
                              printf("L(I_dn)[%d]    ", i);
                    }
                    else {
                         if (output_up)
                              printf("L(I_up)[%d]   ", i);
                         if (output_down)
                              printf("L(I_dn)[%d]   ", i);
                    }
               }
               printf("\n");

               for (i = 0; i < n_out_levels; ++i) {
                    for (j = 0; j < n_out_thetas2; ++j) {
                         for (k = 0; k < n_out_phis; ++k) {
                              for (l = 0; l < n_stokes; ++l) {
                                   if (options & XRTM_OPTION_OUTPUT_AT_LEVELS)
                                        printf("%14d", out_levels[i]);
                                   else
                                        printf("%14e", out_taus  [i]); 
                                   printf("%14e%14e%10d", acos(qx[i_out_thetas2 + j])*R2D, out_phis[k], l);
                                   if (output_up)
                                        printf("%14e", I_p[i][j][k][l]);
                                   if (output_down)
                                        printf("%14e", I_m[i][j][k][l]);
                                   if (o.derivs) {
                                        for (m = 0; m < n_derivs; ++m) {
                                             if (output_up)
                                                  printf("%14e", K_p[i][m][j][k][l]);
                                             if (output_down)
                                                  printf("%14e", K_m[i][m][j][k][l]);
                                        }
                                   }
                                   printf("\n");
                              }
                         }
                    }
               }

               printf("\n");
          }


          /*--------------------------------------------------------------------
           *
           *------------------------------------------------------------------*/
          if (o.radiance_mean) {
               if (options & XRTM_OPTION_OUTPUT_AT_LEVELS)
                    printf("level         ");
               else
                    printf("tau           ");
               if (output_up)
                    printf("  mean_up      ");
               if (output_down)
                    printf("  mean_dn      ");
               for (i = 0; i < n_derivs; ++i) {
                    if (output_up)
                         printf("L(mean_up)[%d] ", i);
                    if (output_down)
                         printf("L(mean_dn)[%d] ", i);
               }
               printf("\n");

               for (i = 0; i < n_out_levels; ++i) {
                    if (options & XRTM_OPTION_OUTPUT_AT_LEVELS)
                         printf("%14d", out_levels[i]);
                    else
                         printf("%14e", out_taus  [i]); 
                    if (output_up)
                         printf("%14e", mean_p[i]);
                    if (output_down)
                         printf("%14e", mean_m[i]);
                    if (o.derivs) {
                         for (m = 0; m < n_derivs; ++m) {
                              if (output_up)
                                   printf("%14e", 0.);
/*
                                   printf("%14e", mean_p_l[i][m]);
*/
                              if (output_down)
                                   printf("%14e", 0.);
/*
                                   printf("%14e", mean_m_l[i][m]);
*/
                         }
                    }
                    printf("\n");
               }

               printf("\n");
          }


          /*--------------------------------------------------------------------
           *
           *------------------------------------------------------------------*/
          if (o.flux) {
               if (options & XRTM_OPTION_OUTPUT_AT_LEVELS)
                    printf("level         ");
               else
                    printf("tau           ");
               if (output_up)
                    printf("  flux_up      ");
               if (output_down)
                    printf("  flux_dn      ");
               for (i = 0; i < n_derivs; ++i) {
                    if (output_up)
                         printf("L(flux_up)[%d] ", i);
                    if (output_down)
                         printf("L(flux_dn)[%d] ", i);
               }
               printf("\n");

               for (i = 0; i < n_out_levels; ++i) {
                    if (options & XRTM_OPTION_OUTPUT_AT_LEVELS)
                         printf("%14d", out_levels[i]);
                    else
                         printf("%14e", out_taus  [i]); 
                    if (output_up)
                         printf("%14e", flux_p[i]);
                    if (output_down)
                         printf("%14e", flux_m[i]);
                    if (o.derivs) {
                         for (m = 0; m < n_derivs; ++m) {
                              if (output_up)
                                   printf("%14e", 0.);
/*
                                   printf("%14e", flux_p_l[i][m]);
*/
                              if (output_down)
                                   printf("%14e", 0.);
/*
                                   printf("%14e", flux_m_l[i][m]);
*/
                         }
                    }
                    printf("\n");
               }

               printf("\n");
          }

          if (o.flux_divergence) {

          }

          if (o.timing)
               printf("\n");
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (o.timing)
          printf("%.2f\n", (double) (time2 - time1) / (double) CLOCKS_PER_SEC);


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     xrtm_destroy(&d);

     if (md.fd_method >= 0)
          xrtm_fd_destroy(&fd);

     misc_free(&md);

     if (o.out_phis)
          free_array1_d(out_phis);

     if (o.radiance) {
          free_array1_d(qx);
          free_array1_d(qw);
     }

     if (options & XRTM_OPTION_OUTPUT_AT_LEVELS)
          free_array1_i(out_levels);
     else
          free_array1_d(out_taus);

     free_array1(out_phis2);

     if (o.radiance) {
          free_array4_d(I_p);
          free_array4_d(I_m);

          if (o.derivs) {
               free_array5_d(K_p);
               free_array5_d(K_m);
          }
     }

     if (o.radiance_mean) {
          free_array1_d(mean_p);
          free_array1_d(mean_m);

          if (o.derivs) {
               free_array2_d(mean_p_l);
               free_array2_d(mean_m_l);
          }
     }

     if (o.flux) {
          free_array1_d(flux_p);
          free_array1_d(flux_m);

          if (o.derivs) {
               free_array2_d(flux_p_l);
               free_array2_d(flux_m_l);
          }
     }

     if (o.flux_divergence) {
          free_array1_d(flux_div);

          if (o.derivs)
               free_array2_d(flux_div_l);
     }


     exit(0);
}



void usage() {

     printf("%s <options>\n", NAME);
}



void version() {

     printf("%s version %s\n", NAME, VERSION);
}
