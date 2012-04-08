/******************************************************************************%
**
**    Copyright (C) 2007-2012 Greg McGarragh <gregm@atmos.colostate.edu>
**
**    This source code is licensed under the GNU General Public License (GPL),
**    Version 3.  See the file COPYING for more details.
**
*******************************************************************************/

#include <gutil.h>
#include <rtutil_scat_io.h>

#include <xrtm_interface.h>

#include "input_util.h"

#include "test.h"
#include "test_macros.h"
#ifdef HAVE_PTHREADS_LIBRARY
#include <unistd.h>
#endif

#ifdef NAME
#undef NAME
#endif
#define NAME    "testxrtm"


typedef struct {
     int check_diffs;
     int on_failed_diff_write_xrtm_cmd;
     int on_failed_diff_stop;
     int help;
     int quick_run;
     int regression_check;
     int on_regression_write_xrtm_cmd;
     int on_regression_stop;
     int test_core;
     int test_errors;
     int write_results;
     int write_xrtm_cmd;
     int version;
} options_data;


void usage();
void version();


int main(int argc, char *argv[]) {

     char temp[1024];

     const char *filename;

     const char *cmd;

     int i;
     int n;

     int flag;

     int n_threads = 1;
/*
     FILE *fp;
*/
     options_data o;

     test_data t;

#if PLATFORM == WIN32_MSVC
     _set_output_format(_TWO_DIGIT_EXPONENT);
#endif

     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     o.check_diffs                   = 1;
     o.on_failed_diff_write_xrtm_cmd = 0;
     o.on_failed_diff_stop           = 0;
     o.help                          = 0;
     o.quick_run                     = 0;
     o.regression_check              = 1;
     o.on_regression_write_xrtm_cmd  = 0;
     o.on_regression_stop            = 0;
     o.test_core                     = 0;
     o.test_errors                   = 0;
     o.write_results                 = 1;
     o.write_xrtm_cmd                = 1;
     o.version                       = 0;
#ifdef HAVE_PTHREADS_LIBRARY
     n_threads = sysconf(_SC_NPROCESSORS_ONLN);
#endif
     n = 0;
     for (i = 1; i < argc; ++i) {
          if (argv[i][0] == '-') {
               if (strcmp(argv[i], "-check_diffs") == 0)
                    o.check_diffs = 1;
               if (strcmp(argv[i], "-on_failed_diff_write_xrtm_cmd") == 0)
                    o.on_failed_diff_write_xrtm_cmd = 1;
               if (strcmp(argv[i], "-on_failed_diff_stop") == 0)
                    o.on_failed_diff_stop = 1;
               else if (strcmp(argv[i], "-help") == 0) {
                    usage();
                    exit(0);
               }
               else if (strcmp(argv[i], "-quick_run") == 0)
                    o.quick_run = 1;
               else if (strcmp(argv[i], "-regression_check") == 0)
                    o.regression_check = 1;
               else if (strcmp(argv[i], "-on_regression_write_xrtm_cmd") == 0)
                    o.on_regression_write_xrtm_cmd = 1;
               else if (strcmp(argv[i], "-on_regression_stop") == 0)
                    o.on_regression_stop = 1;
               else if (strcmp(argv[i], "-test_core") == 0)
                    o.test_core = 1;
               else if (strcmp(argv[i], "-test_errors") == 0)
                    o.test_errors = 1;
               else if (strcmp(argv[i], "-write_results") == 0)
                    o.write_results = 1;
               else if (strcmp(argv[i], "-write_xrtm_cmd") == 0)
                    o.write_xrtm_cmd = 1;
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
     t.index                         = 0;

     t.check_diffs                   = o.check_diffs;
     t.on_failed_diff_write_xrtm_cmd = o.on_failed_diff_write_xrtm_cmd;
     t.on_failed_diff_stop           = o.on_failed_diff_stop;
     t.quick_run                     = o.quick_run;
     t.on_regression_write_xrtm_cmd  = o.on_regression_write_xrtm_cmd;
     t.on_regression_stop            = o.on_regression_stop;
     t.write_results                 = o.write_results;
     t.write_xrtm_cmd                = o.write_xrtm_cmd;

     t.n_threads                     = n_threads;


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     t.max_coef = 0;
     for (i = 0; i < MAX_PFS; ++i) {
          sprintf(temp, "gc6:pfs/phase_test%d.gc", i + 1);

          t.gf[i] = strdup(temp);

          if ((t.n_gc[i] = load_scat_coefs(temp, -1, &t.gc[i], &flag)) < 0) {
               eprintf("ERROR: load_scat_coefs()\n");
               exit(1);
          }

          if (t.n_gc[i] > t.max_coef)
               t.max_coef = t.n_gc[i];

          sprintf(temp, "gc6:pfs/phase_test%d_l.gc", i + 1);

          t.gf_l[i] = strdup(temp);

          if ((t.n_gc[i] = load_scat_coefs(temp, -1, &t.gc_l[i], &flag)) < 0) {
               eprintf("ERROR: load_scat_coefs()\n");
               exit(1);
          }
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (o.test_core) {
          if (t.write_results) {
               filename = "testxrtm_results.out";
               if ((t.fp_results = fopen(filename, "w")) == NULL) {
                    printf("ERROR: Error opening file for writing: %s ... %s\n",
                           filename, strerror(errno));
                    exit(1);
               }
          }

          if (t.write_xrtm_cmd) {
               filename = "testxrtm_xrtm_cmd.out";
               if ((t.fp_xrtm_cmd = fopen(filename, "w")) == NULL) {
                    printf("ERROR: Error opening file for writing: %s ... %s\n",
                           filename, strerror(errno));
                    exit(1);
               }
          }

          HANDLE_RETURN(test_core(&t), "test_core()");

          if (t.write_results)
               fclose(t.fp_results);

          if (t.write_xrtm_cmd)
               fclose(t.fp_xrtm_cmd);

          if (n_threads > 1) {
               cmd = "sort -o testxrtm_results.out testxrtm_results.out";
               if (system(cmd)) {
                   printf("ERROR: system(): '%s'\n", cmd);
                   exit(1);
               }
          }
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (o.test_errors) {
          if (t.write_results) {
               filename = "testxrtm_errors.out";

               if (freopen(filename, "w", stderr) == NULL) {
                    printf("ERROR: Error opening file for writing: %s ... %s\n",
                           filename, strerror(errno));
                    exit(1);
               }
/*
               if ((t.fp_results = fopen(filename, "w")) == NULL) {
                    printf("ERROR: Error opening file for writing: %s ... %s\n",
                           filename, strerror(errno));
                    exit(1);
               }

               fp = stderr;
               stderr = t.fp_results;
*/
          }

          HANDLE_RETURN(test_errors(&t), "test_errors()");

          if (t.write_results) {
/*
               stderr = fp;
               fclose(t.fp_results);
*/
          }
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


     exit(0);
}


void usage() {

     printf("%s <options>\n", NAME);
}



void version() {

     printf("%s version %s\n", NAME, VERSION);
}
