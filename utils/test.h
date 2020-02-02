/*******************************************************************************
**
**    Copyright (C) 2007-2020 Greg McGarragh <greg.mcgarragh@colostate.edu>
**
**    This source code is licensed under the GNU General Public License (GPL),
**    Version 3.  See the file COPYING for more details.
**
*******************************************************************************/

#ifndef TEST_H
#define TEST_H

#include <xrtm_interface.h>

#include "input_util.h"

#include "test_result.h"
#include "test_xrtm.h"

#ifdef __cplusplus
extern "C" {
#endif


#define MAX_PFS				5


#define DEFAULT_DOUB_D_TAU		1.e-8

#define DEFAULT_SOS_PARAMS_MAX_OS	128
#define DEFAULT_SOS_PARAMS_MAX_TAU	.05
#define DEFAULT_SOS_PARAMS_SOS_TOL	1.e-9

#define DEFAULT_FOURIER_TOL		1.e-9

#define DEFAULT_PLANET_R		6371.0071823

#define DEFAULT_F_0			1.
#define DEFAULT_PHI_0			0.

#define DEFAULT_F_ISO_TOP		0.
#define DEFAULT_F_ISO_BOT		0.


#define DEF_TOL				1.e-4
#define DEF_TOL_L			1.e-3


#define MAX_COMMAND_LENGTH		8192


enum test_bound_type {
     TEST_BOUND_SIMPLE,
     TEST_BOUND_BOUNDS,

     N_TEST_BOUNDS
};


enum test_stack_type {
     TEST_STACK_ONE_LAYER,
     TEST_STACK_SMALL_STACK,
     TEST_STACK_LARGE_STACK,

     N_TEST_STACKS
};


enum test_derivs_type {
     TEST_DERIVS_NO_DERIVS,
     TEST_DERIVS_ONE_DERIVS,
     TEST_DERIVS_BOUND_DERIVS,

     N_TEST_DERIVS
};


enum test_output_at_levels_type {
     TEST_OUTPUT_AT_LEVELS_TOP,
     TEST_OUTPUT_AT_LEVELS_BOTTOM,
     TEST_OUTPUT_AT_LEVELS_TOP_AND_BOTTOM,
     TEST_OUTPUT_AT_LEVELS_ALL_LEVELS,

     N_TEST_OUTPUT_AT_LEVELS
};


enum test_output_at_taus_type {
     TEST_OUTPUT_AT_TAUS_TWO_END_MIDDLES,
     TEST_OUTPUT_AT_TAUS_BOUNDS_AND_MIDDLES,

     N_TEST_OUTPUT_AT_TAUS
};


typedef struct {
     int index;

     int core_check_diffs;
     int core_dont_execute;
     int core_echo_xrtm_cmd;
     int core_fill_blanks;
     int core_include_dev_solvers;
     int core_on_failed_diff_write_xrtm_cmd;
     int core_on_failed_diff_stop_with_error;
     int core_write_results_bin;
     int core_write_results_text;
     int core_write_xrtm_cmd;
     int core_zero_solar_source;
     int core_zero_thermal_source;

     int n_threads;

     int exact_options;
     int required_options;
     int unwanted_options;

     int exact_solvers;
     int required_solvers;
     int unwanted_solvers;

     int max_coef;

     char *gf[MAX_PFS];
     char *gf_l[MAX_PFS];
     int n_gc[MAX_PFS];
     double **gc[MAX_PFS];
     double **gc_l[MAX_PFS];

     char *lf[MAX_PFS];
     char *lf_l[MAX_PFS];
     int n_lc[MAX_PFS];
     double **lc[MAX_PFS];
     double **lc_l[MAX_PFS];

     int n_quad;
     int n_stokes;
     int n_out_thetas;
     int n_out_phis;

     enum test_bound_type bound_type;
     enum test_stack_type stack_type;

     enum test_derivs_type derivs_type;

     enum test_output_at_levels_type output_at_levels_type;
     enum test_output_at_taus_type   output_at_taus_type;

     FILE *fp_core_results_bin;
     FILE *fp_core_results_bin_w_deltas;
     FILE *fp_core_results_text;
     FILE *fp_core_results_text_w_deltas;
     FILE *fp_core_xrtm_cmd;

     test_xrtm_data test_xrtm_list;
} test_data;


#ifdef INCLUDE_DEV_SOURCE
#include "../utils2/test2.h"
#endif


#ifndef INCLUDE_DEV_SOURCE
#define GROUP_SOLVERS_DEV
#define GROUP_SOLVERS_DEV_DERIVS1
#define GROUP_SOLVERS_DEV_DERIVS12

#define DEV_RANGE_OVER_N_QUAD			(0)
#define DEV_RANGE_OVER_N_STOKES			(0)
#define DEV_RANGE_OVER_N_DERIVS			(0)
#define DEV_RANGE_OVER_N_LAYERS			(0)
#define DEV_RANGE_OVER_N_OUT_THETAS		(0)
#define DEV_RANGE_OVER_OPTICAL_PROPERTY_INPUTS	(0)
#define DEV_RANGE_OVER_EIGEN_SOVERS		(0)
#define DEV_OPT_OUT_RANGE_OVER_OPTIONS(MASK)
#define DEV_RANGE_OVER_SAVING_OPTIONS		(0)
#endif


#ifdef __cplusplus
}
#endif

#endif /* TEST_H */
