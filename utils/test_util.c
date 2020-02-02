/*******************************************************************************
**
**    Copyright (C) 2007-2020 Greg McGarragh <greg.mcgarragh@colostate.edu>
**
**    This source code is licensed under the GNU General Public License (GPL),
**    Version 3.  See the file COPYING for more details.
**
*******************************************************************************/

#include <gutil.h>
#include <gindex_name_value.h>

#include <xrtm_interface.h>

#include "test.h"
#include "test_util.h"


/*******************************************************************************
 *
 ******************************************************************************/
static const char *test_bound_names[] = {
     "simple",
     "bounds"
};

static long test_bound_types[] = {
     TEST_BOUND_SIMPLE,
     TEST_BOUND_BOUNDS
};

GINDEX_NAME_VALUE_TEMPLATE(test_bound, "test bound", N_TEST_BOUNDS)


static const char *test_stack_names[] = {
     "one_layer",
     "short_stack",
     "large_stack"
};

static long test_stack_types[] = {
     TEST_STACK_ONE_LAYER,
     TEST_STACK_SMALL_STACK,
     TEST_STACK_LARGE_STACK
};

GINDEX_NAME_VALUE_TEMPLATE(test_stack, "test stack", N_TEST_STACKS)


static const char *test_derivs_names[] = {
     "no_derivs",
     "one_derivs",
     "bound_derivs"
};

static long test_derivs_types[] = {
     TEST_DERIVS_NO_DERIVS,
     TEST_DERIVS_ONE_DERIVS,
     TEST_DERIVS_BOUND_DERIVS
};

GINDEX_NAME_VALUE_TEMPLATE(test_derivs, "test derivs", N_TEST_DERIVS)


static const char *test_output_at_levels_names[] = {
     "top",
     "bottom",
     "top_and_bottom",
     "all_levels"
};

static long test_output_at_levels_types[] = {
     TEST_OUTPUT_AT_LEVELS_TOP,
     TEST_OUTPUT_AT_LEVELS_BOTTOM,
     TEST_OUTPUT_AT_LEVELS_TOP_AND_BOTTOM,
     TEST_OUTPUT_AT_LEVELS_ALL_LEVELS
};

GINDEX_NAME_VALUE_TEMPLATE(test_output_at_levels, "output at levels", N_TEST_OUTPUT_AT_LEVELS)


static const char *test_output_at_taus_names[] = {
     "two_end_middles",
     "bounds_and_middles"
};

static long test_output_at_taus_types[] = {
     TEST_OUTPUT_AT_TAUS_TWO_END_MIDDLES,
     TEST_OUTPUT_AT_TAUS_BOUNDS_AND_MIDDLES
};

GINDEX_NAME_VALUE_TEMPLATE(test_output_at_taus, "output at taus", N_TEST_OUTPUT_AT_TAUS)


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
int get_callxrtm_cmd(xrtm_data *gd, misc_data *md, int n_phis, double *phis, char *s, size_t size) {

     int i;

     FILE *fp;

     if ((fp = fmemopen(s, size, "w")) == NULL) {
          fprintf(stderr, "ERROR: Error opening file for writing: %s\n", strerror(errno));
          return -1;
     }

     fprintf(fp, "callxrtm -input_string \"");

     if (xrtm_fwrite_input_fp(gd, NULL, md, fp, ' ', 0, 1, 1, 0)) {
          fprintf(stderr, "ERROR: xrtm_fwrite_input(): %s\n", "test_input.txt");
          return -1;
     }

     fprintf(fp, "\" -out_phis ");

     for (i = 0; i < n_phis; ++i) {
          fprintf(fp, "%e", phis[i]);
          if (i + 1 < n_phis)
               fprintf(fp, ", ");
     }

     fclose(fp);

     return 0;
}
