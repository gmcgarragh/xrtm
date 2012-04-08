/******************************************************************************%
**
**    Copyright (C) 2007-2012 Greg McGarragh <gregm@atmos.colostate.edu>
**
**    This source code is licensed under the GNU General Public License (GPL),
**    Version 3.  See the file COPYING for more details.
**
*******************************************************************************/

#ifndef TEST_H
#define TEST_H

#include <xrtm_interface.h>

#include "input_util.h"

#ifdef __cplusplus
extern "C" {
#endif


#include "../config.h"


#define MAX_PFS				5

#define MAX_N_QUAD_POWER_SCALER		5
#define MAX_N_QUAD_POWER_VECTOR		3
/*
#define MAX_N_QUAD_POWER_SCALER		8
#define MAX_N_QUAD_POWER_VECTOR		6
*/
#define MAX_N_DERIV_POWER		8

#define MAX_N_LAYER_POWER		8

#define MAX_N_LEVEL_POWER		8

#define MAX_N_THETA_POWER		5
/*
#define MAX_N_THETA_POWER		8
*/

#define DEFAULT_N_QUAD_SCALER		8
#define DEFAULT_N_QUAD_VECTOR		4

#define DEFAULT_DOUB_D_TAU		1.e-8

#define DEFAULT_SOS_PARAMS_MAX_OS	128

#define DEFAULT_SOS_PARAMS_MAX_TAU	.05
/*
#define DEFAULT_SOS_PARAMS_MAX_TAU	.005
*/
#define DEFAULT_SOS_PARAMS_SOS_TOL	1.e-9

#define DEFAULT_FOURIER_TOL		1.e-9

#define DEF_TOL				1.e-4
#define DEF_TOL_L			1.e-3


typedef struct {
     int index;

     int check_diffs;
     int on_failed_diff_write_xrtm_cmd;
     int on_failed_diff_stop;
     int quick_run;
     int on_regression_write_xrtm_cmd;
     int on_regression_stop;
     int write_results;
     int write_xrtm_cmd;

     int n_threads;

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

     FILE *fp_results;
     FILE *fp_xrtm_cmd;
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


int test_core(test_data *t);

int test_errors(test_data *t);

int test_execute(xrtm_data *gd, misc_data *md, test_data *td, int index, int n_solvers, enum xrtm_solver_mask *solvers, int n_phis, double *phis, double *tol, double *tol_l, int n_ignore_index_solver, int *ignore_index_solver_list, int *ignore_index_solver_mask, int ignore_solver_mask_ref, int ignore_solver_mask_tran, int n_ignore_deriv, int *ignore_deriv_list, int *ignore_deriv_mask, void *mutex);


#ifdef __cplusplus
}
#endif

#endif /* TEST_H */
