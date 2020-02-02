/*******************************************************************************
**
**    Copyright (C) 2007-2020 Greg McGarragh <greg.mcgarragh@colostate.edu>
**
**    This source code is licensed under the GNU General Public License (GPL),
**    Version 3.  See the file COPYING for more details.
**
*******************************************************************************/

#include <gutil.h>
#include <rtutil_scat_io.h>

#include <xrtm_interface.h>

#include <sys/stat.h>

#include "input_util.h"

#include "test.h"
#include "test_core.h"
#include "test_macros.h"
#include "test_util.h"
#include "test_write.h"

#include <omp.h>


#ifdef NAME
#undef NAME
#endif
#define NAME   "testxrtm"

#define N_TEMP 1024


typedef struct {
     int core;
     int help;
     int version;
} options_data;


void usage();
void version();


int main(int argc, char *argv[]) {

     char temp[N_TEMP];

     const char *prefix;

     int i;
     int ii;

     int n;

     int flag;

     options_data options;

     test_data t;

     test_xrtm_data *test_xrtm;
#if PLATFORM == WIN32_MSVC
     _set_output_format(_TWO_DIGIT_EXPONENT);
#endif

     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     options.core    = 0;
     options.help    = 0;
     options.version = 0;

     t.core_check_diffs                    = 1;
     t.core_dont_execute                   = 0;
     t.core_echo_xrtm_cmd                  = 0;
     t.core_fill_blanks                    = 1;
     t.core_include_dev_solvers            = 1;
     t.core_on_failed_diff_write_xrtm_cmd  = 0;
     t.core_on_failed_diff_stop_with_error = 0;
     t.core_write_results_bin              = 1;
     t.core_write_results_text             = 1;
     t.core_write_xrtm_cmd                 = 1;
     t.core_zero_solar_source              = 0;
     t.core_zero_thermal_source            = 0;

     t.n_threads                           = 8;

     t.exact_options                       = 0;
     t.required_options                    = 0;
     t.unwanted_options                    = 0;

     t.exact_solvers                       = 0;
     t.required_solvers                    = 0;
     t.unwanted_solvers                    = 0;

     t.n_quad                              = 3;
     t.n_stokes                            = 1;
     t.n_out_thetas                        = 1;
     t.n_out_phis                          = 1;

     n = 0;
     for (i = 1; i < argc; ++i) {
          if (argv[i][0] == '-') {
               if (strcmp(argv[i], "-core") == 0) {
                    options.core = 1;
                    if ((int) (t.bound_type = test_bound_name_to_value(argv[++i])) < 0) {
                         fprintf(stderr, "ERROR: Invalid test bound type: %s\n",  argv[i]);
                         exit(1);
                    }
                    if ((int) (t.stack_type = test_stack_name_to_value(argv[++i])) < 0) {
                         fprintf(stderr, "ERROR: Invalid test stack type: %s\n",  argv[i]);
                         exit(1);
                    }
                    if ((int) (t.derivs_type = test_derivs_name_to_value(argv[++i])) < 0) {
                         fprintf(stderr, "ERROR: Invalid test derivs type: %s\n", argv[i]);
                         exit(1);
                    }
                    if ((int) (t.output_at_levels_type = test_output_at_levels_name_to_value(argv[++i])) < 0) {
                         fprintf(stderr, "ERROR: Invalid test output_at_levels type: %s\n", argv[i]);
                         exit(1);
                    }
                    if ((int) (t.output_at_taus_type   = test_output_at_taus_name_to_value  (argv[++i])) < 0) {
                         fprintf(stderr, "ERROR: Invalid test output_at_taus type: %s\n",   argv[i]);
                         exit(1);
                    }
               }
               else if (strcmp(argv[i], "-no-core") == 0)
                    options.core = 0;
               else if (strcmp(argv[i], "-core_check_diffs") == 0)
                    t.core_check_diffs = 1;
               else if (strcmp(argv[i], "-no-core_check_diffs") == 0)
                    t.core_check_diffs = 0;
               else if (strcmp(argv[i], "-core_dont_execute") == 0)
                    t.core_dont_execute = 1;
               else if (strcmp(argv[i], "-no-core_dont_execute") == 0)
                    t.core_dont_execute = 0;
               else if (strcmp(argv[i], "-core_echo_xrtm_cmd") == 0)
                    t.core_echo_xrtm_cmd = 1;
               else if (strcmp(argv[i], "-no-core_echo_xrtm_cmd") == 0)
                    t.core_echo_xrtm_cmd = 0;
               else if (strcmp(argv[i], "-core_fill_blanks") == 0)
                    t.core_fill_blanks = 1;
               else if (strcmp(argv[i], "-no-core_fill_blanks") == 0)
                    t.core_fill_blanks = 0;
               else if (strcmp(argv[i], "-core_include_dev_solvers") == 0)
                    t.core_include_dev_solvers = 1;
               else if (strcmp(argv[i], "-no-core_include_dev_solvers") == 0)
                    t.core_include_dev_solvers = 0;
               else if (strcmp(argv[i], "-core_on_failed_diff_write_xrtm_cmd") == 0)
                    t.core_on_failed_diff_write_xrtm_cmd = 1;
               else if (strcmp(argv[i], "-no-core_on_failed_diff_write_xrtm_cmd") == 0)
                    t.core_on_failed_diff_write_xrtm_cmd = 0;
               else if (strcmp(argv[i], "-core_on_failed_diff_stop_with_error") == 0)
                    t.core_on_failed_diff_stop_with_error = 1;
               else if (strcmp(argv[i], "-core_write_results_bin") == 0)
                    t.core_write_results_bin = 1;
               else if (strcmp(argv[i], "-no-core_write_results_bin") == 0)
                    t.core_write_results_bin = 0;
               else if (strcmp(argv[i], "-core_write_results_text") == 0)
                    t.core_write_results_text = 1;
               else if (strcmp(argv[i], "-no-core_write_results_text") == 0)
                    t.core_write_results_text = 0;
               else if (strcmp(argv[i], "-core_write_xrtm_cmd") == 0)
                    t.core_write_xrtm_cmd = 1;
               else if (strcmp(argv[i], "-no-core_write_xrtm_cmd") == 0)
                    t.core_write_xrtm_cmd = 0;
               else if (strcmp(argv[i], "-core_zero_solar_source") == 0)
                    t.core_zero_solar_source = 1;
               else if (strcmp(argv[i], "-core_zero_thermal_source") == 0)
                    t.core_zero_thermal_source = 1;

               else if (strcmp(argv[i], "-help") == 0) {
                    usage();
                    exit(0);
               }
               else if (strcmp(argv[i], "-n_quad") == 0) {
                    ii = i;
                    check_arg_count(i, argc, 1, argv[ii]);
                    t.n_quad = strtoi_errmsg_exit(argv[++i], argv[ii]);
               }
               else if (strcmp(argv[i], "-n_stokes") == 0) {
                    ii = i;
                    check_arg_count(i, argc, 1, argv[ii]);
                    t.n_stokes = strtoi_errmsg_exit(argv[++i], argv[ii]);
               }
               else if (strcmp(argv[i], "-n_out_thetas") == 0) {
                    ii = i;
                    check_arg_count(i, argc, 1, argv[ii]);
                    t.n_out_thetas = strtoi_errmsg_exit(argv[++i], argv[ii]);
               }
               else if (strcmp(argv[i], "-n_out_phis") == 0) {
                    ii = i;
                    check_arg_count(i, argc, 1, argv[ii]);
                    t.n_out_phis = strtoi_errmsg_exit(argv[++i], argv[ii]);
               }
               else if (strcmp(argv[i], "-n_threads") == 0) {
                    ii = i;
                    check_arg_count(i, argc, 1, argv[ii]);
                    t.n_threads = strtoi_errmsg_exit(argv[++i], argv[ii]);
               }
               else if (strcmp(argv[i], "-exact_options") == 0) {
                    ii = i;
                    check_arg_count(i, argc, 1, argv[ii]);
                    if ((t.exact_options |= xrtm_option_name_list_to_mask(argv[++i])) < 0) {
                         fprintf(stderr, "ERROR: xrtm_option_name_list_to_mask()\n");
                         exit(1);
                    }
               }
               else if (strcmp(argv[i], "-required_options") == 0) {
                    ii = i;
                    check_arg_count(i, argc, 1, argv[ii]);
                    if ((t.required_options |= xrtm_option_name_list_to_mask(argv[++i])) < 0) {
                         fprintf(stderr, "ERROR: xrtm_option_name_list_to_mask()\n");
                         exit(1);
                    }
               }
               else if (strcmp(argv[i], "-unwanted_options") == 0) {
                    ii = i;
                    check_arg_count(i, argc, 1, argv[ii]);
                    if ((t.unwanted_options |= xrtm_option_name_list_to_mask(argv[++i])) < 0) {
                         fprintf(stderr, "ERROR: xrtm_option_name_list_to_mask()\n");
                         exit(1);
                    }
               }
               else if (strcmp(argv[i], "-exact_solvers") == 0) {
                    ii = i;
                    check_arg_count(i, argc, 1, argv[ii]);
                    if ((t.exact_solvers |= xrtm_solver_name_list_to_mask(argv[++i])) < 0) {
                         fprintf(stderr, "ERROR: xrtm_solver_name_list_to_mask()\n");
                         exit(1);
                    }
               }
               else if (strcmp(argv[i], "-required_solvers") == 0) {
                    ii = i;
                    check_arg_count(i, argc, 1, argv[ii]);
                    if ((t.required_solvers |= xrtm_solver_name_list_to_mask(argv[++i])) < 0) {
                         fprintf(stderr, "ERROR: xrtm_solver_name_list_to_mask()\n");
                         exit(1);
                    }
               }
               else if (strcmp(argv[i], "-unwanted_solvers") == 0) {
                    ii = i;
                    check_arg_count(i, argc, 1, argv[ii]);
                    if ((t.unwanted_solvers |= xrtm_solver_name_list_to_mask(argv[++i])) < 0) {
                         fprintf(stderr, "ERROR: xrtm_solver_name_list_to_mask()\n");
                         exit(1);
                    }
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
               if (n == 0)
                    prefix = argv[i];
               else {
                    printf("Too many arguments, use -help for more information\n");
                    exit(1);
               }

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
     t.index = 0;

     list_init(&t.test_xrtm_list);

     omp_set_num_threads(t.n_threads);


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     t.max_coef = 0;
     for (i = 0; i < MAX_PFS; ++i) {
          snprintf(temp, N_TEMP, "gc6:pfs/phase_test%d.gc", i + 1);

          t.gf[i] = strdup(temp);

          if ((t.n_gc[i] = load_scat_coefs(temp, -1, &t.gc[i], &flag)) < 0) {
               fprintf(stderr, "ERROR: load_scat_coefs()\n");
               exit(1);
          }

          if (t.n_gc[i] > t.max_coef)
               t.max_coef = t.n_gc[i];

          snprintf(temp, N_TEMP, "gc6:pfs/phase_test%d_l.gc", i + 1);

          t.gf_l[i] = strdup(temp);

          if ((t.n_gc[i] = load_scat_coefs(temp, -1, &t.gc_l[i], &flag)) < 0) {
               fprintf(stderr, "ERROR: load_scat_coefs()\n");
               exit(1);
          }
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     list_init(&t.test_xrtm_list);


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (options.core) {
          if (t.core_write_results_bin) {
               snprintf(temp, N_TEMP, "%stestxrtm_results_bin.out", prefix);
               if ((t.fp_core_results_bin = fopen(temp, "w")) == NULL) {
                    fprintf(stderr, "ERROR: Error opening file for writing: %s ... %s\n",
                            temp, strerror(errno));
                    exit(1);
               }

               snprintf(temp, N_TEMP, "%stestxrtm_results_bin_w_deltas.out", prefix);
               if ((t.fp_core_results_bin_w_deltas = fopen(temp, "w")) == NULL) {
                    fprintf(stderr, "ERROR: Error opening file for writing: %s ... %s\n",
                            temp, strerror(errno));
                    exit(1);
               }
          }

          if (t.core_write_results_text) {
               snprintf(temp, N_TEMP, "%stestxrtm_results_text.out", prefix);
               if ((t.fp_core_results_text = fopen(temp, "w")) == NULL) {
                    fprintf(stderr, "ERROR: Error opening file for writing: %s ... %s\n",
                            temp, strerror(errno));
                    exit(1);
               }

               snprintf(temp, N_TEMP, "%stestxrtm_results_text_w_deltas.out", prefix);
               if ((t.fp_core_results_text_w_deltas = fopen(temp, "w")) == NULL) {
                    fprintf(stderr, "ERROR: Error opening file for writing: %s ... %s\n",
                            temp, strerror(errno));
                    exit(1);
               }
          }

          if (t.core_write_xrtm_cmd) {
               snprintf(temp, N_TEMP, "%stestxrtm_xrtm_cmd.out", prefix);
               if ((t.fp_core_xrtm_cmd = fopen(temp, "w")) == NULL) {
                    fprintf(stderr, "ERROR: Error opening file for writing: %s ... %s\n",
                            temp, strerror(errno));
                    exit(1);
               }
          }

          HANDLE_RETURN(test_core(&t), "test_core()");

          if (write_cmds(&t)) {
               fprintf(stderr, "ERROR: write_cmds()\n");
               exit(1);
          }

          if (write_results(&t)) {
               fprintf(stderr, "ERROR: write_results()\n");
               exit(1);
          }

          if (t.core_write_results_bin) {
               fclose(t.fp_core_results_bin);
               fclose(t.fp_core_results_bin_w_deltas);
          }

          if (t.core_write_results_text) {
               fclose(t.fp_core_results_text);
               fclose(t.fp_core_results_text_w_deltas);
          }

          if (t.core_write_xrtm_cmd)
               fclose(t.fp_core_xrtm_cmd);
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     for (i = 0; i < MAX_PFS; ++i) {
          free(t.gf[i]);
          free_array2_d(t.gc[i]);

          free(t.gf_l[i]);
          free_array2_d(t.gc_l[i]);
     }

     list_for_each(&t.test_xrtm_list, test_xrtm)
          test_xrtm_free(test_xrtm);
     list_free(&t.test_xrtm_list);


     exit(0);
}


void usage() {

     printf("%s <flags ...>\n", NAME);
     printf("\n");

     printf("See the XRTM documentation for details on running %s.\n", NAME);
}



void version() {

     printf("%s, version %s, build SHA-1: %s, build date: %s\n",
            NAME, VERSION, build_sha_1, build_date);
}
