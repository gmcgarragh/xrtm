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
int get_solution_external(xrtm_data *d, int solver, int solutions, int n_phis, double **phis, double ****I_p, double ****I_m, double *****I_p_l, double *****I_m_l, double *mean_p, double *mean_m, double **mean_p_l, double **mean_m_l, double *flux_p, double *flux_m, double **flux_p_l, double **flux_m_l, double *flux_div, double **flux_div_l, work_data work) {

     double a;

     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (solver & XRTM_SOLVER_DISORT) {
          if (d->n_kernels == 0)
               a = 0.;
          else
               a = d->kernel_ampfac[0];

          if (d->options & XRTM_OPTION_CALC_DERIVS) {
               init_array5_d(I_p_l, d->n_layers + 1, d->n_derivs, n_phis, d->n_quad_x, d->n_stokes, 0.);
               init_array5_d(I_m_l, d->n_layers + 1, d->n_derivs, n_phis, d->n_quad_x, d->n_stokes, 0.);
          }

          if (call_disort(d->n_coef, d->n_quad, d->n_layers, d->qx, d->lambda, d->F_0, d->mu_0, d->phi_0, d->ulevels, d->utaus, d->n_ulevels, d->umus, d->n_umus, phis[0], n_phis, d->F_iso_top, d->levels_b, d->surface_b, a, d->n_coef_layer, d->coef0, d->omega0, d->ltau0, I_p, I_m, mean_p, mean_m, flux_p, flux_m, d->fourier_tol, d->options & XRTM_OPTION_DELTA_M, d->options & XRTM_OPTION_N_T_TMS, solutions & XRTM_OUTPUT_RADIANCE, solutions & XRTM_OUTPUT_RADIANCE_MEAN, solutions & XRTM_OUTPUT_FLUX, d->options & XRTM_OPTION_SOURCE_SOLAR, d->options & XRTM_OPTION_SOURCE_THERMAL, d->options & XRTM_OPTION_OUTPUT_AT_TAUS)) {
               fprintf(stderr, "ERROR: call_disort()\n");
               return -1;
          }
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     else if (solver & XRTM_SOLVER_LIDORT) {
          if (call_lidort(d->n_four2, d->n_elem, d->n_coef, d->n_quad, d->n_derivs, d->n_layers, d->qx, d->F_0, d->mu_0, d->phi_0, d->ulevels, d->utaus, d->n_ulevels, d->umus, d->n_umus, phis[0], n_phis, d->planet_r, d->levels_z, d->n_kernels, d->n_kernel_quad, d->kernels, d->kernel_ampfac, d->kernel_ampfac_l, d->kernel_params, d->kernel_params_l, d->n_coef_layer, d->coef0, d->coef0_l, d->omega0, d->omega0_l, d->ltau0, d->ltau0_l, I_p, I_m, I_p_l, I_m_l, mean_p, mean_m, mean_p_l, mean_m_l, flux_p, flux_m, flux_p_l, flux_m_l, d->fourier_tol, d->options & XRTM_OPTION_DELTA_M, d->options & XRTM_OPTION_N_T_TMS, d->options & XRTM_OPTION_PSA, d->misc_input.use_lidort_quad_output, solutions & XRTM_OUTPUT_RADIANCE, solutions & XRTM_OUTPUT_RADIANCE_MEAN, solutions & XRTM_OUTPUT_FLUX, d->options & XRTM_OPTION_OUTPUT_AT_TAUS, d->derivs.layers)) {
               fprintf(stderr, "ERROR: call_lidort()\n");
               return -1;
          }
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     else if (solver & XRTM_SOLVER_LRAD) {
          update_opt_props_all(d, work);

          if (call_lrad(d->n_four2, d->n_elem, d->n_coef, d->n_quad, d->n_stokes, d->n_derivs, d->n_layers, d->qx, d->lambda, d->F_0, d->mu_0, d->phi_0, d->n_ulevels, d->umus, d->n_umus, phis[0], n_phis, d->planet_r, d->levels_z, d->levels_b, d->surface_b, d->n_kernels, d->n_kernel_quad, d->kernels, d->kernel_ampfac, d->kernel_ampfac_l, d->kernel_params, d->kernel_params_l, d->n_coef_layer, d->coef0, d->coef0_l, d->coef, d->coef_l, d->omega0, d->omega0_l, d->omega, d->omega_l, d->ltau, d->ltau_l, I_p, I_m, I_p_l, I_m_l, d->fourier_tol, d->options & XRTM_OPTION_DELTA_M, d->options & XRTM_OPTION_N_T_TMS, d->options & XRTM_OPTION_PSA, d->derivs.layers)) {
               fprintf(stderr, "ERROR: call_lrad()\n");
               return -1;
          }
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     else if (solver & XRTM_SOLVER_POLRAD) {
          if (d->n_kernels == 0)
               a = 0.;
          else
               a = d->kernel_ampfac[0];

          if (call_polrad(d->n_four2, d->n_coef, d->n_quad, d->n_stokes, d->n_layers, d->qx, d->F_0, d->mu_0, d->phi_0, d->n_ulevels, d->umus, d->n_umus, phis[0], n_phis, a, d->n_coef_layer, d->coef0, d->omega0, d->ltau0, I_p, I_m, d->fourier_tol, d->options & XRTM_OPTION_DELTA_M)) {
               fprintf(stderr, "ERROR: call_polrad()\n");
               return -1;
          }
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     else if (solver & XRTM_SOLVER_RADIANT) {
          if (call_radiant(d->n_four2, d->n_elem, d->n_coef, d->n_quad, d->n_derivs, d->n_layers, d->qx, d->lambda, d->F_0, d->mu_0, d->phi_0, d->utaus, d->n_ulevels, d->umus, d->n_umus, phis[0], n_phis, d->F_iso_top, d->planet_r, d->levels_z, d->levels_b, d->surface_b, d->n_kernels, d->n_kernel_quad, d->kernels, d->kernel_ampfac, d->kernel_ampfac_l, d->kernel_params, d->kernel_params_l, d->n_coef_layer, d->coef0, d->coef0_l, d->omega0, d->omega0_l, d->ltau0, d->ltau0_l, I_p, I_m, I_p_l, I_m_l, flux_p, flux_m, d->fourier_tol, d->quad_type, d->options & XRTM_OPTION_DELTA_M, d->options & XRTM_OPTION_N_T_TMS, d->options & XRTM_OPTION_PSA, solutions & XRTM_OUTPUT_RADIANCE, solutions & XRTM_OUTPUT_FLUX, d->options & XRTM_OPTION_SOURCE_THERMAL, d->options & XRTM_OPTION_OUTPUT_AT_TAUS, d->derivs.layers)) {
               fprintf(stderr, "ERROR: call_radiant()\n");
               return -1;
          }
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     else if (solver & XRTM_SOLVER_RADTRAN3) {
          if (d->n_kernels == 0)
               a = 0.;
          else
               a = d->kernel_ampfac[0];

          if (d->options & XRTM_OPTION_CALC_DERIVS) {
               init_array5_d(I_p_l, d->n_layers + 1, d->n_derivs, n_phis, d->n_quad_x, d->n_stokes, 0.);
               init_array5_d(I_m_l, d->n_layers + 1, d->n_derivs, n_phis, d->n_quad_x, d->n_stokes, 0.);
          }

          if (call_radtran3(d->n_four2, d->n_elem, d->n_coef2, d->n_quad, d->n_stokes, d->n_layers, d->qx, d->lambda, d->F_0, d->mu_0, d->phi_0, d->ulevels, d->n_ulevels, d->umus, d->n_umus, phis[0], n_phis, 0., d->levels_b, d->surface_b, a, d->n_coef_layer, d->coef0, d->omega0, d->ltau0, I_p, I_m, flux_p, flux_m, d->doub_d_tau, d->quad_type, d->options & XRTM_OPTION_DELTA_M, solutions & XRTM_OUTPUT_RADIANCE, solutions & XRTM_OUTPUT_FLUX, d->options & XRTM_OPTION_SOURCE_SOLAR, d->options & XRTM_OPTION_SOURCE_THERMAL)) {
               fprintf(stderr, "ERROR: call_radtran3()\n");
               return -1;
          }
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     else if (solver & XRTM_SOLVER_SOI) {
          if (call_soi(d->n_four2, d->n_elem, d->n_coef, d->n_quad, d->n_layers, d->qx, d->F_0, d->mu_0, d->phi_0, d->umus, d->n_umus, phis[0], n_phis, d->planet_r, d->levels_z, d->n_kernels, d->n_kernel_quad, d->kernels, d->kernel_ampfac, d->kernel_params, d->n_coef_layer, d->coef0, d->omega0, d->ltau0, I_p, I_m, d->doub_d_tau, d->fourier_tol, d->options & XRTM_OPTION_DELTA_M, d->options & XRTM_OPTION_N_T_TMS, d->options & XRTM_OPTION_PSA)) {
               fprintf(stderr, "ERROR: call_soi()\n");
               return -1;
          }
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     else if (solver & XRTM_SOLVER_TWOSTR) {
          if (d->n_kernels == 0)
               a = 0.;
          else
               a = d->kernel_ampfac[0];

          if (d->options & XRTM_OPTION_CALC_DERIVS) {
               init_array5_d(I_p_l, d->n_layers + 1, d->n_derivs, n_phis, d->n_quad_x, d->n_stokes, 0.);
               init_array5_d(I_m_l, d->n_layers + 1, d->n_derivs, n_phis, d->n_quad_x, d->n_stokes, 0.);
          }

          if (call_twostr(d->n_layers, d->lambda, d->F_0, d->mu_0, d->phi_0, d->ulevels, d->utaus, d->n_ulevels, d->F_iso_top, d->planet_r, d->levels_z, d->levels_b, d->surface_b, a, d->g0, d->omega0, d->ltau0, mean_p, mean_m, flux_p, flux_m, flux_div, d->options & XRTM_OPTION_DELTA_M, solutions & XRTM_OUTPUT_RADIANCE_MEAN, solutions & XRTM_OUTPUT_FLUX, solutions & XRTM_OUTPUT_FLUX_DIVERGENCE, d->options & XRTM_OPTION_PSA, d->options & XRTM_OPTION_SOURCE_THERMAL, d->options & XRTM_OPTION_OUTPUT_AT_TAUS)) {
               fprintf(stderr, "ERROR: call_twostr()\n");
               return -1;
          }
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     else if (solver & XRTM_SOLVER_VLIDORT) {
          if (call_vlidort(d->n_four2, d->n_elem, d->n_coef, d->n_quad, d->n_stokes, d->n_derivs, d->n_layers, d->qx, d->lambda, d->F_0, d->mu_0, d->phi_0, d->ulevels, d->utaus, d->n_ulevels, d->umus, d->n_umus, phis[0], n_phis, d->F_iso_top, d->planet_r, d->levels_z, d->levels_b, d->surface_b, d->n_kernels, d->n_kernel_quad, d->kernels, d->kernel_ampfac, d->kernel_ampfac_l, d->kernel_params, d->kernel_params_l, d->n_coef_layer, d->coef0, d->coef0_l, d->omega0, d->omega0_l, d->ltau0, d->ltau0_l, I_p, I_m, I_p_l, I_m_l, mean_p, mean_m, mean_p_l, mean_m_l, flux_p, flux_m, flux_p_l, flux_m_l, d->fourier_tol, d->options & XRTM_OPTION_DELTA_M, d->options & XRTM_OPTION_N_T_TMS, d->options & XRTM_OPTION_PSA, d->misc_input.use_lidort_quad_output, solutions & XRTM_OUTPUT_RADIANCE, solutions & XRTM_OUTPUT_RADIANCE_MEAN, solutions & XRTM_OUTPUT_FLUX, d->options & XRTM_OPTION_SOURCE_THERMAL, d->options & XRTM_OPTION_OUTPUT_AT_TAUS, d->derivs.layers)) {
               fprintf(stderr, "ERROR: call_vlidort()\n");
               return -1;
          }
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
#ifdef DEBUG
     else {
         fprintf(stderr, "ERROR: get_radiance_external(): end of if / else if\n");
         exit(1);
     }
#endif

     return 0;
}
