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

#include "test_result.h"


/*******************************************************************************
 *
 ******************************************************************************/
int test_result_init(test_result_data *d, int n_solvers) {

     d->n_solvers = n_solvers;

     d->solvers = malloc(n_solvers * sizeof(test_result_solver_data));

     return 0;
}



void test_result_free(test_result_data *d) {

     free(d->solvers);
}



/*******************************************************************************
 *
 ******************************************************************************/
int test_cmd_init(test_cmd_data *d, int n_solvers) {

     return 0;
}



void test_cmd_free(test_cmd_data *d) {

     free(d->cmd);
}



/*******************************************************************************
 *
 ******************************************************************************/
int save_result_header(test_result_data *d, int index, int index2, int i_out_level, int i_deriv, double out_theta, double out_phi, int i_stokes) {

     d->index       = index;
     d->index2      = index2;
     d->i_out_level = i_out_level;
     d->i_deriv     = i_deriv;
     d->out_theta   = out_theta;
     d->out_phi     = out_phi;
     d->i_stokes    = i_stokes;

     return 0;
}



int save_result_solver_mask(test_result_data *d, int i_solver, enum xrtm_solver_mask mask) {

     d->solvers[i_solver].mask = mask;

     return 0;
}


int save_result_solver_values(test_result_data *d, int i_solver, double x_p, double x_m) {

     d->solvers[i_solver].x_p = x_p;
     d->solvers[i_solver].x_m = x_m;

     return 0;
}


int save_result_solver_deltas(test_result_data *d, int i_solver, double x_p, double x_m) {

     d->solvers[i_solver].x_p_delta = x_p;
     d->solvers[i_solver].x_m_delta = x_m;

     return 0;
}
