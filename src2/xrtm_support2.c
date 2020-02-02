/*******************************************************************************
**
**    Copyright (C) 2007-2020 Greg McGarragh <greg.mcgarragh@colostate.edu>
**
**    This source code is licensed under the GNU General Public License (GPL),
**    Version 3.  See the file COPYING for more details.
**
*******************************************************************************/


/*******************************************************************************
 *
 ******************************************************************************/
void xrtm_misc_input_init_dev(misc_input_data *d) {

     d->use_lidort_quad_output = USE_LIDORT_QUAD_OUTPUT;

     d->use_pade_gamma_init    = USE_PADE_GAMMA_INIT;

     d->eigen_solver_sym_real  = EIGEN_SOLVER_SYM_REAL;
}



/*******************************************************************************
 *
 ******************************************************************************/
static const char *eigen_solver_sym_real_names[] = {
     "eispack",
     "lapack"
};

int eigen_solver_sym_real_code(char *name) {
     return name_to_index(name, eigen_solver_sym_real_names,
          N_EIGEN_SOLVER_SYM_REAL_TYPES, "symmetric real eigen solver types");
}

const char *eigen_solver_sym_real_name(enum eigen_solver_sym_real_type code) {
     return index_to_name(code, eigen_solver_sym_real_names,
          N_EIGEN_SOLVER_SYM_REAL_TYPES, "symmetric real eigen solver types");
}
