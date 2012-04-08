/******************************************************************************%
**
**    Copyright (C) 2007-2012 Greg McGarragh <gregm@atmos.colostate.edu>
**
**    This source code is licensed under the GNU General Public License (GPL),
**    Version 3.  See the file COPYING for more details.
**
*******************************************************************************/

#ifndef TEST_UTIL_H
#define TEST_UTIL_H

#include <xrtm_model.h>

#include "test.h"

#ifdef __cplusplus
extern "C" {
#endif


int make_solvers(int solvers_mask1, int *solvers_mask2, enum xrtm_solver_mask *solvers_array2, double *tolerance, ...);
int make_solvers2(int solvers_mask1, int *solvers_mask2, enum xrtm_solver_mask *solvers_array2, double *tolerance_ref, double *tolerance_tran, ...);

int make_tolerance(double *tol, ...);

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
                int (*ignore_mask_solver)(double theta_0, double *theta));
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
                       double ***chi_l, double *omega_l,
                       double *ltau_l, double *albedo_l,
                       double *tol, double *tol_l,
                       int n_ignore_index_solver,
                       int *ignore_list_index_solver,
                       int *ignore_mask_index_solver,
                       int n_ignore_deriv,
                       int *ignore_list_deriv,
                       int *ignore_mask_deriv,
                       int (*ignore_mask_solver)(double theta_0, double *theta));


#ifdef __cplusplus
}
#endif

#endif /* TEST_UTIL_H */
