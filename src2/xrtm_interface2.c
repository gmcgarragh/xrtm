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
/*
          void radtran3_lobatto_quadrature_(int *, double *, double *);

          radtran3_lobatto_quadrature_(&d->n_quad, d->qx, d->qw);
*/


/*******************************************************************************
 *
 ******************************************************************************/
int check_dev_solvers_create(int options, int solvers, int max_coef, int n_quad, int n_stokes, int n_derivs, int n_layers, int n_theta_0s, int n_kernels, int n_kernel_quad, enum xrtm_kernel_type *kernels, int n_out_levels, int n_out_thetas, int report) {

     int i;

     if (n_quad > 1) {
         if (check_solvers(solvers, 0, report, "option", "n_quad > 1", XRTM_SOLVER_TWOSTR, 0))
              return XRTM_INT_ERROR;
     }

     if (n_stokes == 1) {
         if (check_solvers(solvers, 0, report, "option", "n_stokes == 1", XRTM_SOLVER_LRAD, 0))
              return XRTM_INT_ERROR;
     }

     if (n_stokes > 3) {
         if (check_solvers(solvers, 0, report, "option", "n_stokes > 3", XRTM_SOLVER_LRAD, 0))
              return XRTM_INT_ERROR;
     }

     if (n_theta_0s != 1) {
          if (check_solvers(solvers, 0, report, "option", "n_theta_0s != 1", XRTM_SOLVER_DISORT, XRTM_SOLVER_LIDORT, XRTM_SOLVER_LRAD, XRTM_SOLVER_POLRAD, XRTM_SOLVER_RADIANT, XRTM_SOLVER_RADTRAN3, XRTM_SOLVER_SOI, XRTM_SOLVER_TWOSTR, XRTM_SOLVER_VLIDORT, 0))
              return XRTM_INT_ERROR;
     }
/*
     if (n_out_levels != 1) {
          if (check_solvers(solvers, 0, report, "option", "n_out_levels != 1", XRTM_SOLVER_LRAD, XRTM_SOLVER_SOI, 0))
              return XRTM_INT_ERROR;
     }
*/
     if (n_out_levels != 2) {
          if (check_solvers(solvers, 0, report, "option", "n_out_levels != 2", XRTM_SOLVER_LRAD, XRTM_SOLVER_SOI, 0))
              return XRTM_INT_ERROR;
     }

     if (n_out_levels != 2 && ! (options & XRTM_OPTION_OUTPUT_AT_TAUS)) {
          if (check_solvers(solvers, 0, report, "option", "n_out_levels != 2 without output_at_taus", XRTM_SOLVER_RADIANT, 0))
              return XRTM_INT_ERROR;
     }

     if (n_out_thetas == 0) {
         if (check_solvers(solvers, 0, report, "option", "n_out_thetas == 0", XRTM_SOLVER_LRAD, XRTM_SOLVER_SOI, 0))
              return XRTM_INT_ERROR;
     }

     if (n_out_thetas > 1) {
         if (check_solvers(solvers, 0, report, "option", "n_out_thetas > 1", XRTM_SOLVER_LRAD, XRTM_SOLVER_RADIANT, XRTM_SOLVER_SOI, 0))
              return XRTM_INT_ERROR;
     }

     if (options & XRTM_OPTION_CALC_DERIVS) {
          if (check_solvers(solvers, 0, report, "option", xrtm_option_mask_to_name(XRTM_OPTION_CALC_DERIVS), XRTM_SOLVER_DISORT, XRTM_SOLVER_RADTRAN3, XRTM_SOLVER_SOI, XRTM_SOLVER_TWOSTR, 0))
               return XRTM_INT_ERROR;
     }

     if (options & XRTM_OPTION_DELTA_M && ! (options & XRTM_OPTION_N_T_TMS)) {
          if (check_solvers(solvers, 0, report, "option", "delta_m without n_t_tms", XRTM_SOLVER_LRAD, 0))
               return XRTM_INT_ERROR;
     }

     if (options & XRTM_OPTION_N_T_TMS) {
          if (check_solvers(solvers, 0, report, "option", xrtm_option_mask_to_name(XRTM_OPTION_N_T_TMS), XRTM_SOLVER_RADTRAN3, XRTM_SOLVER_TWOSTR, 0))
               return XRTM_INT_ERROR;
     }

     if (options & XRTM_OPTION_NO_AZIMUTHAL) {
          if (check_solvers(solvers, 0, report, "option", xrtm_option_mask_to_name(XRTM_OPTION_NO_AZIMUTHAL), XRTM_SOLVER_DISORT, XRTM_SOLVER_LRAD, XRTM_SOLVER_TWOSTR, XRTM_SOLVER_VLIDORT, 0))
               return XRTM_INT_ERROR;
     }

     if (options & XRTM_OPTION_OUTPUT_AT_TAUS) {
          if (check_solvers(solvers, 0, report, "option", xrtm_option_mask_to_name(XRTM_OPTION_OUTPUT_AT_TAUS), XRTM_SOLVER_LRAD, XRTM_SOLVER_RADTRAN3, XRTM_SOLVER_SOI, 0))
               return XRTM_INT_ERROR;
     }

     if (options & XRTM_OPTION_QUAD_NORM_GAUS_LEG) {
          if (check_solvers(solvers, 0, report, "option", xrtm_option_mask_to_name(XRTM_OPTION_QUAD_NORM_GAUS_LEG), XRTM_SOLVER_LRAD, XRTM_SOLVER_DISORT, XRTM_SOLVER_LIDORT, XRTM_SOLVER_SOI, XRTM_SOLVER_VLIDORT, 0))
               return XRTM_INT_ERROR;
     }

     if (options & XRTM_OPTION_QUAD_DOUB_GAUS_LEG) {
          if (check_solvers(solvers, 0, report, "option", xrtm_option_mask_to_name(XRTM_OPTION_QUAD_DOUB_GAUS_LEG), XRTM_SOLVER_TWOSTR, 0))
               return XRTM_INT_ERROR;
     }

     if (options & XRTM_OPTION_QUAD_LOBATTO) {
          if (check_solvers(solvers, 0, report, "option", xrtm_option_mask_to_name(XRTM_OPTION_QUAD_LOBATTO), XRTM_SOLVER_LRAD, XRTM_SOLVER_DISORT, XRTM_SOLVER_LIDORT, XRTM_SOLVER_SOI, XRTM_SOLVER_TWOSTR, XRTM_SOLVER_VLIDORT, 0))
               return XRTM_INT_ERROR;
     }

     if (options & XRTM_OPTION_PSA) {
          if (check_solvers(solvers, 0, report, "option", xrtm_option_mask_to_name(XRTM_OPTION_PSA), XRTM_SOLVER_DISORT, XRTM_SOLVER_RADTRAN3, 0))
               return XRTM_INT_ERROR;
     }

     if (options & XRTM_OPTION_SFI) {
          if (check_solvers(solvers, 0, report, "option", xrtm_option_mask_to_name(XRTM_OPTION_SFI), XRTM_SOLVER_LRAD, XRTM_SOLVER_RADIANT, XRTM_SOLVER_RADTRAN3, XRTM_SOLVER_SOI, 0))
               return XRTM_INT_ERROR;
     }

     if (! (options & XRTM_OPTION_SOURCE_SOLAR)) {
          if (check_solvers(solvers, 1, report, "option", xrtm_option_mask_to_name(XRTM_OPTION_SOURCE_SOLAR), XRTM_SOLVER_LIDORT, XRTM_SOLVER_VLIDORT, 0))
               return XRTM_INT_ERROR;
     }

     if (options & XRTM_OPTION_SOURCE_THERMAL) {
          if (check_solvers(solvers, 0, report, "option", xrtm_option_mask_to_name(XRTM_OPTION_SOURCE_THERMAL), XRTM_SOLVER_LIDORT, XRTM_SOLVER_LRAD, XRTM_SOLVER_SOI, XRTM_SOLVER_VLIDORT, 0))
               return XRTM_INT_ERROR;
     }

     if (options & XRTM_OPTION_TOP_DOWN_ADDING) {
          if (check_solvers(solvers, 0, report, "option", xrtm_option_mask_to_name(XRTM_OPTION_TOP_DOWN_ADDING), XRTM_SOLVER_DISORT, XRTM_SOLVER_LIDORT, XRTM_SOLVER_LRAD, XRTM_SOLVER_RADIANT, XRTM_SOLVER_RADTRAN3, XRTM_SOLVER_SOI, XRTM_SOLVER_TWOSTR, XRTM_SOLVER_VLIDORT, 0))
               return XRTM_INT_ERROR;
     }

     if (options & XRTM_OPTION_BOTTOM_UP_ADDING) {
          if (check_solvers(solvers, 0, report, "option", xrtm_option_mask_to_name(XRTM_OPTION_BOTTOM_UP_ADDING), XRTM_SOLVER_DISORT, XRTM_SOLVER_LIDORT, XRTM_SOLVER_LRAD, XRTM_SOLVER_RADIANT, XRTM_SOLVER_RADTRAN3, XRTM_SOLVER_SOI, XRTM_SOLVER_TWOSTR, XRTM_SOLVER_VLIDORT, 0))
               return XRTM_INT_ERROR;
     }

     if (options & XRTM_OPTION_VECTOR) {
          if (check_solvers(solvers, 0, report, "option", xrtm_option_mask_to_name(XRTM_OPTION_VECTOR), XRTM_SOLVER_DISORT, XRTM_SOLVER_LIDORT, XRTM_SOLVER_RADIANT, XRTM_SOLVER_SOI, XRTM_SOLVER_TWOSTR, 0))
               return XRTM_INT_ERROR;
     }

     for (i = 0; i < n_kernels; ++i) {
          if ((solvers & XRTM_SOLVER_DISORT) && (kernels[i] != XRTM_KERNEL_LAMBERTIAN)) {
               fprintf(stderr, "ERROR: solver \"disort\" only support a lambertian kernel\n");
               return XRTM_INT_ERROR;
          }

          if ((solvers & XRTM_SOLVER_RADTRAN3) && (kernels[i] != XRTM_KERNEL_LAMBERTIAN)) {
               fprintf(stderr, "ERROR: solver \"radtran3\" only support a lambertian kernel\n");
               return XRTM_INT_ERROR;
          }
     }

     return 0;
}



/*******************************************************************************
 *
 ******************************************************************************/
int check_dev_solvers_set_ulevels(xrtm_data *d, int *ulevels) {

     if (ulevels[0] != 0) {
         if (check_solvers(d->solvers, 0, 1, "option", "ulevels[0] != 0", XRTM_SOLVER_LRAD, XRTM_SOLVER_SOI, 0))
              return XRTM_INT_ERROR;
     }

     return 0;
}



/*******************************************************************************
 *
 ******************************************************************************/
int check_dev_solvers_solution(xrtm_data *d, int solver, int solutions, int n_phis, double **phis) {

     int i;
     int j;

     if (n_phis > 1) {
          if (check_solvers(solver, 0, 1, "option", "n_phis > 1", XRTM_SOLVER_LRAD, XRTM_SOLVER_RADIANT, XRTM_SOLVER_SOI, 0))
               return XRTM_INT_ERROR;
     }

     if (solver & XRTM_SOLVERS_EXT_EXTERNAL) {
          for (i = 1; i < d->n_umus; ++i) {
               for (j = 0; j < n_phis; ++j) {
                    if (phis[i][j] != phis[0][j]) {
                         fprintf(stderr, "ERROR: phis[i][j] must equal phis[0][j] for external solvers\n");
                         return XRTM_INT_ERROR;
                    }
               }
          }
     }

     if (solutions & XRTM_OUTPUT_RADIANCE) {
          if (check_solvers(solver, 0, 1, "output option", xrtm_output_mask_to_name(XRTM_OUTPUT_RADIANCE), XRTM_SOLVER_TWOSTR, 0))
               return XRTM_INT_ERROR;
     }

     if (solutions & XRTM_OUTPUT_RADIANCE_MEAN) {
          if (check_solvers(solver, 0, 1, "output option", xrtm_output_mask_to_name(XRTM_OUTPUT_RADIANCE_MEAN), XRTM_SOLVER_LRAD, XRTM_SOLVER_RADIANT, XRTM_SOLVER_RADTRAN3, 0))
               return XRTM_INT_ERROR;
     }

     if (solutions & XRTM_OUTPUT_FLUX) {
          if (check_solvers(solver, 0, 1, "output option", xrtm_output_mask_to_name(XRTM_OUTPUT_FLUX), XRTM_SOLVER_SOI, XRTM_SOLVER_LRAD, 0))
               return XRTM_INT_ERROR;
     }

     if (solutions & XRTM_OUTPUT_FLUX_DIVERGENCE) {
          if (check_solvers(solver, 0, 1, "output option", xrtm_output_mask_to_name(XRTM_OUTPUT_FLUX_DIVERGENCE), XRTM_SOLVER_DISORT, XRTM_SOLVER_LIDORT, XRTM_SOLVER_LRAD, XRTM_SOLVER_RADIANT, XRTM_SOLVER_RADTRAN3, XRTM_SOLVER_VLIDORT, 0))
               return XRTM_INT_ERROR;
     }

     return 0;
}
