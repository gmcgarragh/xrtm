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
#include <gmath_matrix.h>

#include <rtutil_math.h>
#include <rtutil_support.h>

#include <setjmp.h>

#include "xrtm.h"
#include "xrtm_brdf.h"
#include "xrtm_derivs.h"
#include "xrtm_interface.h"
#include "xrtm_model.h"
#include "xrtm_model_a.h"
#include "xrtm_support.h"
#include "xrtm_utility.h"
#include "version.h"


/*******************************************************************************
 *
 ******************************************************************************/
static const char *xrtm_option_names[] = {
     "calc_derivs",
/*
     "forward_derivs",
*/
     "reverse_derivs",
     "delta_m",
     "n_t_tms",
     "four_conv_old",
     "four_conv_new",
     "no_azimuthal",
     "output_at_levels",
     "output_at_taus",
     "part_sol_classical",
     "part_sol_greens",
     "phase_scalar",
     "phase_matrix_gc",
     "phase_matrix_lc",
     "psa",
     "quad_norm_gaus_leg",
     "quad_doub_gaus_leg",
     "quad_lobatto",
     "save_leg_polys",
     "save_phase_mats",
     "save_local_r_t",
     "save_layer_r_t_s",
     "save_total_r_t_s",
     "sfi",
     "source_solar",
     "source_thermal",
     "stack_reuse_adding",
     "top_down_adding",
     "bottom_up_adding",
     "upwelling_output",
     "downwelling_output",
     "vector"
};


static long xrtm_option_masks[] = {
     XRTM_OPTION_CALC_DERIVS,
/*
     XRTM_OPTION_FORWARD_DERIVS,
*/
     XRTM_OPTION_REVERSE_DERIVS,
     XRTM_OPTION_DELTA_M,
     XRTM_OPTION_N_T_TMS,
     XRTM_OPTION_FOUR_CONV_OLD,
     XRTM_OPTION_FOUR_CONV_NEW,
     XRTM_OPTION_NO_AZIMUTHAL,
     XRTM_OPTION_OUTPUT_AT_LEVELS,
     XRTM_OPTION_OUTPUT_AT_TAUS,
     XRTM_OPTION_PART_SOL_CLASSICAL,
     XRTM_OPTION_PART_SOL_GREENS,
     XRTM_OPTION_PHASE_SCALAR,
     XRTM_OPTION_PHASE_MATRIX_GC,
     XRTM_OPTION_PHASE_MATRIX_LC,
     XRTM_OPTION_PSA,
     XRTM_OPTION_QUAD_NORM_GAUS_LEG,
     XRTM_OPTION_QUAD_DOUB_GAUS_LEG,
     XRTM_OPTION_QUAD_LOBATTO,
     XRTM_OPTION_SAVE_LEG_POLYS,
     XRTM_OPTION_SAVE_PHASE_MATS,
     XRTM_OPTION_SAVE_LOCAL_R_T,
     XRTM_OPTION_SAVE_LAYER_R_T_S,
     XRTM_OPTION_SAVE_TOTAL_R_T_S,
     XRTM_OPTION_SFI,
     XRTM_OPTION_SOURCE_SOLAR,
     XRTM_OPTION_SOURCE_THERMAL,
     XRTM_OPTION_STACK_REUSE_ADDING,
     XRTM_OPTION_TOP_DOWN_ADDING,
     XRTM_OPTION_BOTTOM_UP_ADDING,
     XRTM_OPTION_UPWELLING_OUTPUT,
     XRTM_OPTION_DOWNWELLING_OUTPUT,
     XRTM_OPTION_VECTOR
};


GINDEX_NAME_MASK_TEMPLATE(xrtm_option, "xrtm option", N_XRTM_OPTIONS)


static const char *xrtm_solver_names[] = {
     "doub_add",
     "eig_add",
     "eig_bvp",
     "four_stream",
     "mem_bvp",
     "pade_add",
     "single",
     "six_stream",
     "sos",
     "two_os",
     "two_stream",
#ifdef INCLUDE_DEV_SOURCE
     XRTM_DEV_SOLVER_NAMES
#endif
};


static long xrtm_solver_masks[] = {
     XRTM_SOLVER_DOUB_ADD,
     XRTM_SOLVER_EIG_ADD,
     XRTM_SOLVER_EIG_BVP,
     XRTM_SOLVER_FOUR_STREAM,
     XRTM_SOLVER_MEM_BVP,
     XRTM_SOLVER_PADE_ADD,
     XRTM_SOLVER_SINGLE,
     XRTM_SOLVER_SIX_STREAM,
     XRTM_SOLVER_SOS,
     XRTM_SOLVER_TWO_OS,
     XRTM_SOLVER_TWO_STREAM,
#ifdef INCLUDE_DEV_SOURCE
     XRTM_DEV_SOLVER_MASKS
#endif
};


GINDEX_NAME_MASK_TEMPLATE(xrtm_solver, "xrtm solver", N_XRTM_SOLVERS)


static const char *xrtm_input_names[] = {
     "",
     "doub_d_tau",
     "pade_params",
     "sos_params",
     "fourier_tol",
     "lambda",
     "planet_r",
     "levels_z",
     "out_levels",
     "out_taus",
     "out_thetas",
     "F_iso_top",
     "F_iso_top_l",
     "F_iso_bot",
     "F_iso_bot_l",
     "F_0",
     "theta_0",
     "phi_0",
     "levels_b",
     "levels_b_l",
     "surface_b",
     "surface_b_l",
     "g",
     "g_l",
     "coef",
     "coef_l",
     "omega",
     "omega_l",
     "ltau",
     "ltau_l",
     "kernel_ampfac",
     "kernel_ampfac_l",
     "kernel_params",
     "kernel_params_l"
     "kernel_funcs"
};


enum xrtm_input_type {
     INPUT_DOUB_D_TAU = 1,
     INPUT_PADE_PARAMS,
     INPUT_SOS_PARAMS,
     INPUT_FOURIER_TOL,
     INPUT_LAMBDA,
     INPUT_PLANET_R,
     INPUT_LEVELS_Z,
     INPUT_OUT_LEVELS,
     INPUT_OUT_TAUS,
     INPUT_OUT_THETAS,
     INPUT_F_ISO_TOP,
     INPUT_F_ISO_TOP_L,
     INPUT_F_ISO_BOT,
     INPUT_F_ISO_BOT_L,
     INPUT_F_0,
     INPUT_THETA_0,
     INPUT_PHI_0,
     INPUT_LEVELS_B,
     INPUT_LEVELS_B_L,
     INPUT_SURFACE_B,
     INPUT_SURFACE_B_L,
     INPUT_G,
     INPUT_G_L,
     INPUT_COEF,
     INPUT_COEF_L,
     INPUT_OMEGA,
     INPUT_OMEGA_L,
     INPUT_LTAU,
     INPUT_LTAU_L,
     INPUT_KERNEL_AMPFAC,
     INPUT_KERNEL_AMPFAC_L,
     INPUT_KERNEL_PARAMS,
     INPUT_KERNEL_PARAMS_L,
     INPUT_KERNEL_FUNCS,

     N_INPUTS
};

static long xrtm_input_types[] = {
     INPUT_DOUB_D_TAU,
     INPUT_PADE_PARAMS,
     INPUT_SOS_PARAMS,
     INPUT_FOURIER_TOL,
     INPUT_LAMBDA,
     INPUT_PLANET_R,
     INPUT_LEVELS_Z,
     INPUT_OUT_LEVELS,
     INPUT_OUT_TAUS,
     INPUT_OUT_THETAS,
     INPUT_F_ISO_TOP,
     INPUT_F_ISO_TOP_L,
     INPUT_F_ISO_BOT,
     INPUT_F_ISO_BOT_L,
     INPUT_F_0,
     INPUT_THETA_0,
     INPUT_PHI_0,
     INPUT_LEVELS_B,
     INPUT_LEVELS_B_L,
     INPUT_SURFACE_B,
     INPUT_SURFACE_B_L,
     INPUT_G,
     INPUT_G_L,
     INPUT_COEF,
     INPUT_COEF_L,
     INPUT_OMEGA,
     INPUT_OMEGA_L,
     INPUT_LTAU,
     INPUT_LTAU_L,
     INPUT_KERNEL_AMPFAC,
     INPUT_KERNEL_AMPFAC_L,
     INPUT_KERNEL_PARAMS,
     INPUT_KERNEL_PARAMS_L
};


GINDEX_NAME_VALUE_TEMPLATE_STATIC(xrtm_input, "xrtm input", N_INPUTS)



static const char *xrtm_output_names[] = {
     "radiance",
     "radiance_mean",
     "flux",
     "flux_divergence"
};

static long xrtm_output_masks[] = {
     XRTM_OUTPUT_RADIANCE,
     XRTM_OUTPUT_RADIANCE_MEAN,
     XRTM_OUTPUT_FLUX,
     XRTM_OUTPUT_FLUX_DIVERGENCE
};


GINDEX_NAME_MASK_TEMPLATE(xrtm_output, "xrtm output", N_XRTM_OUTPUTS)


/*******************************************************************************
 *
 ******************************************************************************/
void xrtm_misc_input_init(misc_input_data *d) {

     d->do_not_add_sfi_ss          = DO_NOT_ADD_SFI_SS;

     d->use_pade_check_condition   = USE_PADE_CHECK_CONDITION;

     d->use_rebuild_stacks         = USE_REBUILD_STACKS;

     d->use_symmetric_form         = USE_SYMMETRIC_FORM;

     d->eigen_solver_gen_real      = EIGEN_SOLVER_GEN_REAL;
     d->eigen_solver_gen_complex   = EIGEN_SOLVER_GEN_COMPLEX;

     d->threshold_sym_mul_block    = THRESHOLD_SYM_MUL_BLOCK;

     d->threshold_mu_0_singlarity  = THRESHOLD_MU_0_SINGLARITY;
     d->threshold_omega_singlarity = THRESHOLD_OMEGA_SINGLARITY;
#ifdef INCLUDE_DEV_SOURCE
     xrtm_misc_input_init_dev(d);
#endif
}



/*******************************************************************************
 *
 ******************************************************************************/
int check_solvers(int solvers, int not_, int report, const char *type, const char *option, ...) {

     int i;

     int solver;

     va_list ap;

     va_start(ap, option);

     for (i = 0; (solver = va_arg(ap, int)) != 0; ++i) {
          if (solvers & solver) {
               if (report) {
                    if (not_)
                         fprintf(stderr, "ERROR: solver \"%s\" requires %s \"%s\"\n", xrtm_solver_mask_to_name((enum xrtm_solver_mask) solver), type, option);
                    else
                         fprintf(stderr, "ERROR: solver \"%s\" does not support %s \"%s\"\n", xrtm_solver_mask_to_name((enum xrtm_solver_mask) solver), type, option);
               }
               return -1;
          }
     }

     va_end(ap);

     return 0;
}



/*******************************************************************************
 *
 ******************************************************************************/
static int xrtm_quad_to_quad_type(xrtm_data *d) {

     if (d->options & XRTM_OPTION_QUAD_NORM_GAUS_LEG)
          return QUAD_NORM_GAUS_LEG;
     else
     if (d->options & XRTM_OPTION_QUAD_DOUB_GAUS_LEG)
          return QUAD_DOUB_GAUS_LEG;
     else
#ifdef DEBUG
     if (d->options & XRTM_OPTION_QUAD_LOBATTO)
#endif
          return QUAD_LOBATTO;
#ifdef DEBUG
     else {
         fprintf(stderr, "ERROR: xrtm_quad_type_to_quad_type(): end of if / else if\n");
         exit(1);
     }
#endif
}



/*******************************************************************************
 *
 ******************************************************************************/
static void init_deps(xrtm_data *d, int dep, int i_layer, int n_layers) {

     int i;
     int j;

     d->dep_flags_utaus   = dep;
     d->dep_flags_chapman = dep;

     for (i = 0; i < d->n_four; ++i) {
          d->dep_flags_Y_p[i]               = dep;
          d->dep_flags_Y_0[i]               = dep;
          d->dep_flags_Y_u[i]               = dep;
/*
          d->dep_flags_diff_bound_input0[i] = dep;
*/
          d->dep_flags_diff_bound_input [i] = dep;
     }

     for (i = i_layer; i < n_layers; ++i) {
          d->dep_flags_opt_props   [i]      = dep;
/*
          d->dep_flags_beam_params0[i]      = dep;
*/
          d->dep_flags_beam_params [i]      = dep;
     }

     for (i = 0; i < d->n_four; ++i) {
          d->dep_flags_total_R_T_S_U_W_V[i] = dep;

          for (j = i_layer; j < n_layers; ++j) {
               d->dep_flags_phase_mats_qq    [i][j] = dep;
               d->dep_flags_phase_mats_uq    [i][j] = dep;
               d->dep_flags_phase_vecs_q0    [i][j] = dep;
               d->dep_flags_phase_vecs_u0    [i][j] = dep;
               d->dep_flags_local_r_t_u_w    [i][j] = dep;
               d->dep_flags_layer_R_T_S_U_W_V[i][j] = dep;
          }
     }
}



static void set_deps(xrtm_data *d, int dep, int i_layer, int n_layers) {

     int i;
     int j;

     d->dep_flags_utaus   |= dep;
     d->dep_flags_chapman |= dep;

     for (i = 0; i < d->n_four; ++i) {
          d->dep_flags_Y_p[i]               |= dep;
          d->dep_flags_Y_0[i]               |= dep;
          d->dep_flags_Y_u[i]               |= dep;
/*
          d->dep_flags_diff_bound_input0[i] |= dep;
*/
          d->dep_flags_diff_bound_input [i] |= dep;
     }

     for (i = i_layer; i < n_layers; ++i) {
          d->dep_flags_opt_props   [i]      |= dep;
/*
          d->dep_flags_beam_params0[i]      |= dep;
*/
          d->dep_flags_beam_params [i]      |= dep;
     }

     for (i = 0; i < d->n_four; ++i) {
          d->dep_flags_total_R_T_S_U_W_V[i] |= dep;

          for (j = i_layer; j < n_layers; ++j) {
               d->dep_flags_phase_mats_qq    [i][j] |= dep;
               d->dep_flags_phase_mats_uq    [i][j] |= dep;
               d->dep_flags_phase_vecs_q0    [i][j] |= dep;
               d->dep_flags_phase_vecs_u0    [i][j] |= dep;
               d->dep_flags_local_r_t_u_w    [i][j] |= dep;
               d->dep_flags_layer_R_T_S_U_W_V[i][j] |= dep;
          }
     }
}



/*******************************************************************************
 *
 ******************************************************************************/
const char *xrtm_get_build_sha_1() {

     return build_sha_1;
}


const char *xrtm_get_build_date() {

     return build_date;
}



/*******************************************************************************
 *
 ******************************************************************************/
int xrtm_check_options(int options, int report) {
/*
     if (! (options & XRTM_OPTION_N_T_TMS) && (options & XRTM_OPTION_FOUR_CONV_NEW)) {
          if (report) fprintf(stderr, "ERROR: option \"%s\" required for option \"%s\"\n", xrtm_option_mask_to_name(XRTM_OPTION_N_T_TMS), xrtm_option_mask_to_name(XRTM_OPTION_FOUR_CONV_NEW));
          return XRTM_INT_ERROR;
     }
*/
     if ((options & XRTM_OPTION_N_T_TMS) && (options & XRTM_OPTION_NO_AZIMUTHAL)) {
          if (report) fprintf(stderr, "ERROR: options \"%s\" and \"%s\" cannot be used together\n", xrtm_option_mask_to_name(XRTM_OPTION_N_T_TMS), xrtm_option_mask_to_name(XRTM_OPTION_NO_AZIMUTHAL));
          return XRTM_INT_ERROR;
     }

     if (! (options & (XRTM_OPTION_FOUR_CONV_OLD | XRTM_OPTION_FOUR_CONV_NEW))) {
          if (report) fprintf(stderr, "ERROR: must specify a convergence type as an option\n");
          return XRTM_INT_ERROR;
     }

     if ((options & XRTM_OPTION_FOUR_CONV_OLD) && (options & XRTM_OPTION_FOUR_CONV_NEW)) {
          if (report) fprintf(stderr, "ERROR: options \"%s\" and \"%s\" cannot be used together\n", xrtm_option_mask_to_name(XRTM_OPTION_FOUR_CONV_OLD), xrtm_option_mask_to_name(XRTM_OPTION_FOUR_CONV_NEW));
          return XRTM_INT_ERROR;
     }

     if (! (options & (XRTM_OPTION_OUTPUT_AT_LEVELS | XRTM_OPTION_OUTPUT_AT_TAUS))) {
          if (report) fprintf(stderr, "ERROR: must specify an output level type (\"%s\" or \"%s\") as an option\n", xrtm_option_mask_to_name(XRTM_OPTION_OUTPUT_AT_LEVELS), xrtm_option_mask_to_name(XRTM_OPTION_OUTPUT_AT_TAUS));
          return XRTM_INT_ERROR;
     }

     if ((options & XRTM_OPTION_OUTPUT_AT_LEVELS) && (options & XRTM_OPTION_OUTPUT_AT_TAUS)) {
          if (report) fprintf(stderr, "ERROR: options \"%s\" and \"%s\" cannot be used together\n", xrtm_option_mask_to_name(XRTM_OPTION_OUTPUT_AT_LEVELS), xrtm_option_mask_to_name(XRTM_OPTION_OUTPUT_AT_TAUS));
          return XRTM_INT_ERROR;
     }

     if (! (options & (XRTM_OPTION_PART_SOL_CLASSICAL | XRTM_OPTION_PART_SOL_GREENS))) {
          if (report) fprintf(stderr, "ERROR: must specify a particular solution type (\"%s\" or \"%s\") as an option\n", xrtm_option_mask_to_name(XRTM_OPTION_PART_SOL_CLASSICAL), xrtm_option_mask_to_name(XRTM_OPTION_PART_SOL_GREENS));
          return XRTM_INT_ERROR;
     }

     if ((options & XRTM_OPTION_PART_SOL_CLASSICAL) && (options & XRTM_OPTION_PART_SOL_GREENS)) {
          if (report) fprintf(stderr, "ERROR: options \"%s\" and \"%s\" cannot be used together\n", xrtm_option_mask_to_name(XRTM_OPTION_PART_SOL_CLASSICAL), xrtm_option_mask_to_name(XRTM_OPTION_PART_SOL_GREENS));
          return XRTM_INT_ERROR;
     }

     if (! (options & (XRTM_OPTION_PHASE_SCALAR | XRTM_OPTION_PHASE_MATRIX_GC | XRTM_OPTION_PHASE_MATRIX_LC))) {
          if (report) fprintf(stderr, "ERROR: must specify a phase function type as an option\n");
          return XRTM_INT_ERROR;
     }

     if (((options & XRTM_OPTION_PHASE_SCALAR)    && (options & XRTM_OPTION_PHASE_MATRIX_GC)) ||
         ((options & XRTM_OPTION_PHASE_SCALAR)    && (options & XRTM_OPTION_PHASE_MATRIX_LC)) ||
         ((options & XRTM_OPTION_PHASE_MATRIX_GC) && (options & XRTM_OPTION_PHASE_MATRIX_LC))) {
          if (report) fprintf(stderr, "ERROR: options \"%s\", \"%s\", and \"%s\" cannot be used together\n", xrtm_option_mask_to_name(XRTM_OPTION_PHASE_SCALAR), xrtm_option_mask_to_name(XRTM_OPTION_PHASE_MATRIX_GC), xrtm_option_mask_to_name(XRTM_OPTION_PHASE_MATRIX_LC));
          return XRTM_INT_ERROR;
     }

     if (! (options & (XRTM_OPTION_QUAD_NORM_GAUS_LEG | XRTM_OPTION_QUAD_DOUB_GAUS_LEG | XRTM_OPTION_QUAD_LOBATTO))) {
          if (report) fprintf(stderr, "ERROR: must specify the quadrature type as an option\n");
          return XRTM_INT_ERROR;
     }

     if (((options & XRTM_OPTION_QUAD_NORM_GAUS_LEG) && (options & XRTM_OPTION_QUAD_DOUB_GAUS_LEG)) ||
         ((options & XRTM_OPTION_QUAD_NORM_GAUS_LEG) && (options & XRTM_OPTION_QUAD_LOBATTO)) ||
         ((options & XRTM_OPTION_QUAD_DOUB_GAUS_LEG) && (options & XRTM_OPTION_QUAD_LOBATTO))) {
          if (report) fprintf(stderr, "ERROR: options \"%s\", \"%s\", and \"%s\" cannot be used together\n", xrtm_option_mask_to_name(XRTM_OPTION_QUAD_NORM_GAUS_LEG), xrtm_option_mask_to_name(XRTM_OPTION_QUAD_DOUB_GAUS_LEG), xrtm_option_mask_to_name(XRTM_OPTION_QUAD_LOBATTO));
          return XRTM_INT_ERROR;
     }

     if (options & XRTM_OPTION_SAVE_LAYER_R_T_S && options & XRTM_OPTION_STACK_REUSE_ADDING) {
          if (report) fprintf(stderr, "ERROR: cannot use option \"%s\" with option \"%s\"\n", xrtm_option_mask_to_name(XRTM_OPTION_SAVE_LAYER_R_T_S), xrtm_option_mask_to_name(XRTM_OPTION_STACK_REUSE_ADDING));
          return XRTM_INT_ERROR;
     }

     if (options & XRTM_OPTION_STACK_REUSE_ADDING && ! (options & XRTM_OPTION_CALC_DERIVS)) {
          if (report) fprintf(stderr, "ERROR: option \"%s\" can only be used with \"%s\"\n", xrtm_option_mask_to_name(XRTM_OPTION_STACK_REUSE_ADDING), xrtm_option_mask_to_name(XRTM_OPTION_CALC_DERIVS));
          return XRTM_INT_ERROR;
     }

     if (options & XRTM_OPTION_STACK_REUSE_ADDING && (options & XRTM_OPTION_TOP_DOWN_ADDING)) {
          if (report) fprintf(stderr, "ERROR: cannot use option \"%s\" with option \"%s\"\n", xrtm_option_mask_to_name(XRTM_OPTION_STACK_REUSE_ADDING), xrtm_option_mask_to_name(XRTM_OPTION_TOP_DOWN_ADDING));
          return XRTM_INT_ERROR;
     }

     if ((options & XRTM_OPTION_TOP_DOWN_ADDING) && (options & XRTM_OPTION_BOTTOM_UP_ADDING)) {
          if (report) fprintf(stderr, "ERROR: options \"%s\" and \"%s\" cannot be used together\n", xrtm_option_mask_to_name(XRTM_OPTION_TOP_DOWN_ADDING), xrtm_option_mask_to_name(XRTM_OPTION_BOTTOM_UP_ADDING));
          return XRTM_INT_ERROR;
     }

     if (! (options & (XRTM_OPTION_UPWELLING_OUTPUT | XRTM_OPTION_DOWNWELLING_OUTPUT))) {
          if (report) fprintf(stderr, "ERROR: must specify an output direction (\"%s\" or \"%s\") as an option\n", xrtm_option_mask_to_name(XRTM_OPTION_UPWELLING_OUTPUT), xrtm_option_mask_to_name(XRTM_OPTION_DOWNWELLING_OUTPUT));
          return XRTM_INT_ERROR;
     }

     if (! (options & XRTM_OPTION_VECTOR) &&
         (options & (XRTM_OPTION_PHASE_MATRIX_GC | XRTM_OPTION_PHASE_MATRIX_LC))) {
          if (report) fprintf(stderr, "ERROR: options \"%s\" or \"%s\" require option \"%s\"\n", xrtm_option_mask_to_name(XRTM_OPTION_PHASE_MATRIX_GC), xrtm_option_mask_to_name(XRTM_OPTION_PHASE_MATRIX_LC), xrtm_option_mask_to_name(XRTM_OPTION_VECTOR));
          return XRTM_INT_ERROR;
     }

     if (options & XRTM_OPTION_VECTOR &&
         ! (options & (XRTM_OPTION_PHASE_MATRIX_GC | XRTM_OPTION_PHASE_MATRIX_LC))) {
          if (report) fprintf(stderr, "ERROR: option \"%s\" requires option \"%s\" or \"%s\"\n", xrtm_option_mask_to_name(XRTM_OPTION_VECTOR), xrtm_option_mask_to_name(XRTM_OPTION_PHASE_MATRIX_GC), xrtm_option_mask_to_name(XRTM_OPTION_PHASE_MATRIX_LC));
          return XRTM_INT_ERROR;
     }

     return 0;
}



/*******************************************************************************
 *
 ******************************************************************************/
int xrtm_check_counts(int options, int max_coef, int n_quad, int n_stokes, int n_derivs, int n_layers, int n_theta_0s, int n_kernels, int n_kernel_quad, enum xrtm_kernel_type *kernels, int n_out_levels, int n_out_thetas, int report) {

     if (! (options & XRTM_OPTION_CALC_DERIVS) && n_derivs > 0) {
          if (report) fprintf(stderr, "ERROR: option \"%s\" required for n_derivs > 0\n", xrtm_option_mask_to_name(XRTM_OPTION_CALC_DERIVS));
          return XRTM_INT_ERROR;
     }

     if ((options & XRTM_OPTION_CALC_DERIVS) && n_derivs == 0) {
          if (report) fprintf(stderr, "ERROR: n_derivs > 0 required for option \"%s\"\n", xrtm_option_mask_to_name(XRTM_OPTION_CALC_DERIVS));
          return XRTM_INT_ERROR;
     }

     if (options & XRTM_OPTION_DELTA_M && max_coef < 2 * n_quad + 1) {
          if (report) fprintf(stderr, "ERROR: invalid value for max_coef: %d, must be >= (2 * n_quad + 1) = %d for \"%s\"\n", max_coef, 2 * n_quad + 1, xrtm_option_mask_to_name(XRTM_OPTION_DELTA_M));
          return XRTM_INT_ERROR;
     }

     if (options & XRTM_OPTION_SFI && n_out_thetas == 0) {
          if (report) fprintf(stderr, "ERROR: n_out_thetas cannot be zero with option \"%s\"\n", xrtm_option_mask_to_name(XRTM_OPTION_SFI));
          return XRTM_INT_ERROR;
     }

     if (options & XRTM_OPTION_TOP_DOWN_ADDING && n_out_levels != 1) {
          if (report) fprintf(stderr, "ERROR: option \"%s\" only works for n_out_levels == 1\n", xrtm_option_mask_to_name(XRTM_OPTION_TOP_DOWN_ADDING));
          return XRTM_INT_ERROR;
     }

     if (options & XRTM_OPTION_BOTTOM_UP_ADDING && n_out_levels != 1) {
          if (report) fprintf(stderr, "ERROR: option \"%s\" only works for n_out_levels == 1\n", xrtm_option_mask_to_name(XRTM_OPTION_BOTTOM_UP_ADDING));
          return XRTM_INT_ERROR;
     }

     if (! (options & XRTM_OPTION_VECTOR) && n_stokes > 1) {
          if (report) fprintf(stderr, "ERROR: option \"%s\" required for n_stokes > 1\n", xrtm_option_mask_to_name(XRTM_OPTION_VECTOR));
          return XRTM_INT_ERROR;
     }

     return 0;
}



/*******************************************************************************
 *
 ******************************************************************************/
int xrtm_check_solvers(int options, int solvers, int max_coef, int n_quad, int n_stokes, int n_derivs, int n_layers, int n_theta_0s, int n_kernels, int n_kernel_quad, enum xrtm_kernel_type *kernels, int n_out_levels, int n_out_thetas, int report) {

     if (n_quad != 1) {
          if (check_solvers(solvers, 0, report, "option", "n_quad != 1", XRTM_SOLVER_TWO_STREAM, 0))
              return XRTM_INT_ERROR;
     }

     if (n_quad != 2) {
          if (check_solvers(solvers, 0, report, "option", "n_quad != 2", XRTM_SOLVER_FOUR_STREAM, 0))
              return XRTM_INT_ERROR;
     }

     if (n_quad != 3) {
          if (check_solvers(solvers, 0, report, "option", "n_quad != 3", XRTM_SOLVER_SIX_STREAM, 0))
              return XRTM_INT_ERROR;
     }

     if (n_theta_0s != 1) {
          if (check_solvers(solvers, 0, report, "option", "n_theta_0s != 1", XRTM_SOLVER_EIG_ADD, XRTM_SOLVER_EIG_BVP, XRTM_SOLVER_FOUR_STREAM, XRTM_SOLVER_MEM_BVP, XRTM_SOLVER_PADE_ADD, XRTM_SOLVER_SINGLE, XRTM_SOLVER_SIX_STREAM, XRTM_SOLVER_SOS, XRTM_SOLVER_TWO_OS, XRTM_SOLVER_TWO_STREAM, 0))
              return XRTM_INT_ERROR;
     }

     if (n_out_levels > 2) {
          if (check_solvers(solvers, 0, report, "option", "n_out_levels > 2", XRTM_SOLVER_DOUB_ADD, XRTM_SOLVER_EIG_ADD, XRTM_SOLVER_PADE_ADD, XRTM_SOLVER_SOS, XRTM_SOLVER_TWO_OS, 0))
              return XRTM_INT_ERROR;
     }

     if (n_out_thetas > 0 && ! (options & XRTM_OPTION_SFI)) {
          if (check_solvers(solvers, 0, report, "option", "n_out_thetas > 0 without sfi turned on", XRTM_SOLVER_TWO_STREAM, XRTM_SOLVER_FOUR_STREAM, XRTM_SOLVER_SIX_STREAM, 0))
              return XRTM_INT_ERROR;
     }

     if (options & XRTM_OPTION_DELTA_M && ! (options & XRTM_OPTION_N_T_TMS)) {
          if (check_solvers(solvers, 0, report, "option", "delta_m without n_t_tms", XRTM_SOLVER_SINGLE, 0))
               return XRTM_INT_ERROR;
     }

     if (options & XRTM_OPTION_N_T_TMS) {
          if (check_solvers(solvers, 0, report, "option", xrtm_option_mask_to_name(XRTM_OPTION_N_T_TMS), XRTM_SOLVER_TWO_OS, 0))
               return XRTM_INT_ERROR;
     }

     if (options & XRTM_OPTION_NO_AZIMUTHAL) {
          if (check_solvers(solvers, 0, report, "option", xrtm_option_mask_to_name(XRTM_OPTION_NO_AZIMUTHAL), XRTM_SOLVER_SINGLE, 0))
               return XRTM_INT_ERROR;
     }

     if (options & XRTM_OPTION_OUTPUT_AT_TAUS) {
          if (check_solvers(solvers, 0, report, "option", xrtm_option_mask_to_name(XRTM_OPTION_OUTPUT_AT_TAUS), XRTM_SOLVER_DOUB_ADD, XRTM_SOLVER_EIG_ADD, XRTM_SOLVER_PADE_ADD, XRTM_SOLVER_SOS, 0))
               return XRTM_INT_ERROR;
     }

     if (options & XRTM_OPTION_QUAD_LOBATTO) {
          if (check_solvers(solvers, 0, report, "option", xrtm_option_mask_to_name(XRTM_OPTION_QUAD_LOBATTO), XRTM_SOLVER_DOUB_ADD, XRTM_SOLVER_EIG_ADD, XRTM_SOLVER_EIG_BVP, XRTM_SOLVER_FOUR_STREAM, XRTM_SOLVER_MEM_BVP, XRTM_SOLVER_PADE_ADD, XRTM_SOLVER_SINGLE, XRTM_SOLVER_SIX_STREAM, XRTM_SOLVER_SOS, XRTM_SOLVER_TWO_OS, XRTM_SOLVER_TWO_STREAM, 0))
               return -1;
     }

     if (options & XRTM_OPTION_SFI) {
          if (check_solvers(solvers, 0, report, "option", xrtm_option_mask_to_name(XRTM_OPTION_SFI), XRTM_SOLVER_DOUB_ADD, XRTM_SOLVER_EIG_ADD, XRTM_SOLVER_MEM_BVP, XRTM_SOLVER_PADE_ADD,                     XRTM_SOLVER_SOS, 0))
/*
          if (check_solvers(solvers, 0, report, "option", xrtm_option_mask_to_name(XRTM_OPTION_SFI), XRTM_SOLVER_DOUB_ADD, XRTM_SOLVER_EIG_ADD, XRTM_SOLVER_MEM_BVP, XRTM_SOLVER_PADE_ADD, XRTM_SOLVER_SINGLE, XRTM_SOLVER_SOS, 0))
*/
               return XRTM_INT_ERROR;
     }

     if (! (options & XRTM_OPTION_SOURCE_SOLAR)) {
          if (check_solvers(solvers, 0, report, "option", xrtm_option_mask_to_name(XRTM_OPTION_SOURCE_SOLAR), XRTM_SOLVER_MEM_BVP, XRTM_SOLVER_SINGLE, XRTM_SOLVER_SOS, XRTM_SOLVER_TWO_OS, 0))
               return XRTM_INT_ERROR;
     }

     if (options & XRTM_OPTION_SOURCE_THERMAL) {
          if (check_solvers(solvers, 0, report, "option", xrtm_option_mask_to_name(XRTM_OPTION_SOURCE_THERMAL), XRTM_SOLVER_FOUR_STREAM, XRTM_SOLVER_MEM_BVP, XRTM_SOLVER_SINGLE, XRTM_SOLVER_SIX_STREAM, XRTM_SOLVER_SOS, XRTM_SOLVER_TWO_OS, XRTM_SOLVER_TWO_STREAM, 0))
               return XRTM_INT_ERROR;
     }

     if (options & XRTM_OPTION_STACK_REUSE_ADDING) {
          if (check_solvers(solvers, 0, report, "option", xrtm_option_mask_to_name(XRTM_OPTION_STACK_REUSE_ADDING), XRTM_SOLVER_EIG_BVP, XRTM_SOLVER_MEM_BVP, XRTM_SOLVER_SINGLE, XRTM_SOLVER_SOS, XRTM_SOLVER_TWO_OS, 0))
               return XRTM_INT_ERROR;
     }

     if (options & XRTM_OPTION_TOP_DOWN_ADDING) {
          if (check_solvers(solvers, 0, report, "option", xrtm_option_mask_to_name(XRTM_OPTION_TOP_DOWN_ADDING), XRTM_SOLVER_EIG_BVP, XRTM_SOLVER_MEM_BVP, XRTM_SOLVER_SINGLE, XRTM_SOLVER_SOS, XRTM_SOLVER_TWO_OS, 0))
               return XRTM_INT_ERROR;
     }

     if (options & XRTM_OPTION_BOTTOM_UP_ADDING) {
          if (check_solvers(solvers, 0, report, "option", xrtm_option_mask_to_name(XRTM_OPTION_BOTTOM_UP_ADDING), XRTM_SOLVER_EIG_BVP, XRTM_SOLVER_MEM_BVP, XRTM_SOLVER_SINGLE, XRTM_SOLVER_SOS, XRTM_SOLVER_TWO_OS, 0))
               return XRTM_INT_ERROR;
     }
#ifdef INCLUDE_DEV_SOURCE
     if (check_dev_solvers_create(options, solvers, max_coef, n_quad, n_stokes, n_derivs, n_layers, n_theta_0s, n_kernels, n_kernel_quad, kernels, n_out_levels, n_out_thetas, report)) {
          if (report) fprintf(stderr, "ERROR: check_dev_solvers_create()\n");
          return XRTM_INT_ERROR;
     }
#endif
     return 0;
}



/*******************************************************************************
 *
 ******************************************************************************/
static void update_n_four2(xrtm_data *d, int flag1, int flag2) {

     if (d->options & XRTM_OPTION_NO_AZIMUTHAL ||
        (! (d->options & XRTM_OPTION_SOURCE_SOLAR) && (d->options & XRTM_OPTION_SOURCE_THERMAL)) ||
        (d->n_stokes == 1 && ((flag1 && d->mu_0_0 == 1.) || (flag2 && d->n_umus == 1 && d->umus[0] == 1.))))
/*
     if (d->options & XRTM_OPTION_NO_AZIMUTHAL || (! (d->options & XRTM_OPTION_SOURCE_SOLAR) && (d->options & XRTM_OPTION_SOURCE_THERMAL)) || (d->n_stokes == 1 && ((flag1 && d->mu_0_0 == 1.) || (flag2 && d->n_umus == 1 && d->umus[0] == 1.))) || ((flag1 && d->mu_0_0 == 1.) && (flag2 && d->n_umus == 1 && d->umus[*] != 1.)))
*/
          d->n_four2 = 1;
     else
          d->n_four2 = d->n_four;
}



/*******************************************************************************
 *
 ******************************************************************************/
int xrtm_create(xrtm_data *d, int options, int solvers, int max_coef, int n_quad, int n_stokes, int n_derivs, int n_layers, int n_theta_0s, int n_kernels, int n_kernel_quad, enum xrtm_kernel_type *kernels, int n_out_levels, int n_out_thetas) {

     misc_input_data misc_input;

     xrtm_misc_input_init(&misc_input);

     return xrtm_create_misc(d, options, solvers, max_coef, n_quad, n_stokes, n_derivs, n_layers, n_theta_0s, n_kernels, n_kernel_quad, kernels, n_out_levels, n_out_thetas, &misc_input);
}



xrtm_data *xrtm_create2(int options, int solvers, int max_coef, int n_quad, int n_stokes, int n_derivs, int n_layers, int n_theta_0s, int n_kernels, int n_kernel_quad, enum xrtm_kernel_type *kernels, int n_out_levels, int n_out_thetas) {

     xrtm_data *d;

     d = (xrtm_data *) malloc(sizeof(xrtm_data));

     if (xrtm_create(d, options, solvers, max_coef, n_quad, n_stokes, n_derivs, n_layers, n_theta_0s, n_kernels, n_kernel_quad, kernels, n_out_levels, n_out_thetas) == XRTM_INT_ERROR) {
          fprintf(stderr, "ERROR: xrtm_create()\n");
          return NULL;
     }

    return d;
}



int xrtm_create_misc(xrtm_data *d, int options, int solvers, int max_coef, int n_quad, int n_stokes, int n_derivs, int n_layers, int n_theta_0s, int n_kernels, int n_kernel_quad, enum xrtm_kernel_type *kernels, int n_out_levels, int n_out_thetas, misc_input_data *misc_input) {

     int i;
     int ii;
     int j;
     int jj;

     int two_n_quad;

     int r;

     jmp_buf env;


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     memset(d, 0xffffffff, sizeof(xrtm_data));


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     d->options        = options;
     d->solvers        = solvers;
     d->n_coef         = max_coef;
     d->n_quad         = n_quad;
     d->n_stokes       = n_stokes;
     d->n_derivs       = n_derivs;
     d->n_layers       = n_layers;
     d->n_theta_0s     = n_theta_0s;
     d->n_kernels      = n_kernels;
     d->n_kernel_quad  = n_kernel_quad;
     d->n_ulevels      = n_out_levels;
     d->n_out_levels   = n_out_levels;
     d->n_umus         = n_out_thetas;
     d->n_out_thetas   = n_out_thetas;
     d->n_out_thetas_2 = d->n_out_thetas == 0 ? d->n_quad : d->n_out_thetas;

     d->misc_input    = *misc_input;


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (d->n_coef < 0) {
          fprintf(stderr, "ERROR: invalid value for max_coef: %d, must be > zero\n",
                  d->n_coef);
          return XRTM_INT_ERROR;
     }

     if (d->n_quad <= 0) {
          fprintf(stderr, "ERROR: invalid value for n_quad: %d, must be > zero\n",
                  d->n_quad);
          return XRTM_INT_ERROR;
     }

     if (d->n_stokes <= 0 || d->n_stokes > 4) {
          fprintf(stderr, "ERROR: invalid value for n_stokes: %d, must be > zero\n",
                  d->n_stokes);
          return XRTM_INT_ERROR;
     }

     if (d->n_derivs < 0) {
          fprintf(stderr, "ERROR: invalid value for n_derivs: %d, must be >= zero\n",
                  d->n_stokes);
          return XRTM_INT_ERROR;
     }

     if (d->n_layers <= 0) {
          fprintf(stderr, "ERROR: invalid value for n_layers: %d, must be > zero\n",
                  d->n_layers);
          return XRTM_INT_ERROR;
     }

     if (d->n_theta_0s <= 0) {
          fprintf(stderr, "ERROR: invalid value for n_theta_0s: %d, must be > zero\n",
                  d->n_theta_0s);
          return XRTM_INT_ERROR;
     }

     if (d->n_kernels < 0) {
          fprintf(stderr, "ERROR: invalid value for n_kernels: %d, must be > zero\n",
                  d->n_kernels);
          return XRTM_INT_ERROR;
     }

     if (d->n_kernels > 0 && d->n_kernel_quad <= 0) {
          fprintf(stderr, "ERROR: invalid value for n_kernel_quad: %d, must be > zero\n",
                  d->n_kernel_quad);
          return XRTM_INT_ERROR;
     }

     for (i = 0; i < d->n_kernels; ++i) {
          if (! kernel_is_valid(kernels[i])) {
               fprintf(stderr, "ERROR: i_kernel %d is invalid\n", i);
               return XRTM_INT_ERROR;
          }
     }

     if (d->n_ulevels <= 0) {
          fprintf(stderr, "ERROR: invalid value for n_out_levels: %d, must be > zero\n",
                  d->n_ulevels);
          return XRTM_INT_ERROR;
     }

     if (d->n_umus < 0) {
          fprintf(stderr, "ERROR: invalid value for n_out_thetas: %d, must be >= zero\n",
                  d->n_umus);
          return XRTM_INT_ERROR;
     }


     /*-------------------------------------------------------------------------
      * Implicit configuration settings.
      *-----------------------------------------------------------------------*/
     if (d->n_stokes > 1)
          d->options |= XRTM_OPTION_VECTOR;

     if (d->options & XRTM_OPTION_N_T_TMS)
          d->options |= XRTM_OPTION_DELTA_M;

     if (! (d->options & (XRTM_OPTION_FOUR_CONV_OLD | XRTM_OPTION_FOUR_CONV_NEW)))
          d->options |= XRTM_OPTION_FOUR_CONV_NEW;

     if (! (d->options & (XRTM_OPTION_PART_SOL_CLASSICAL | XRTM_OPTION_PART_SOL_GREENS)))
          d->options |= XRTM_OPTION_PART_SOL_CLASSICAL;

     if (! (d->options & (XRTM_OPTION_PHASE_SCALAR | XRTM_OPTION_PHASE_MATRIX_GC | XRTM_OPTION_PHASE_MATRIX_LC))) {
          if (! (d->options & XRTM_OPTION_VECTOR))
               d->options |= XRTM_OPTION_PHASE_SCALAR;
          else
               d->options |= XRTM_OPTION_PHASE_MATRIX_GC;
     }

     if (! (d->options & (XRTM_OPTION_QUAD_NORM_GAUS_LEG | XRTM_OPTION_QUAD_DOUB_GAUS_LEG | XRTM_OPTION_QUAD_LOBATTO)))
          d->options |= XRTM_OPTION_QUAD_DOUB_GAUS_LEG;

     if (! (d->options & (XRTM_OPTION_UPWELLING_OUTPUT | XRTM_OPTION_DOWNWELLING_OUTPUT)))
          d->options |= XRTM_OPTION_UPWELLING_OUTPUT | XRTM_OPTION_DOWNWELLING_OUTPUT;


     /*-------------------------------------------------------------------------
      * Check if options are valid.
      *-----------------------------------------------------------------------*/
     if (xrtm_check_options(d->options, 1)) {
          fprintf(stderr, "ERROR: xrtm_check_options()\n");
          return XRTM_INT_ERROR;
     }


     /*-------------------------------------------------------------------------
      * Check for count incompatibilities.
      *-----------------------------------------------------------------------*/
    if (xrtm_check_counts(d->options, d->n_coef, d->n_quad, d->n_stokes, d->n_derivs, d->n_layers, d->n_theta_0s, d->n_kernels, d->n_kernel_quad, kernels, d->n_ulevels, d->n_umus, 1)) {
          fprintf(stderr, "ERROR: xrtm_check_counts()\n");
          return XRTM_INT_ERROR;
     }


     /*-------------------------------------------------------------------------
      * Check for solver incompatibilities.
      *-----------------------------------------------------------------------*/
    if (xrtm_check_solvers(d->options, d->solvers, d->n_coef, d->n_quad, d->n_stokes, d->n_derivs, d->n_layers, d->n_theta_0s, d->n_kernels, d->n_kernel_quad, kernels, d->n_ulevels, d->n_umus, 1)) {
          fprintf(stderr, "ERROR: xrtm_check_solvers()\n");
          return XRTM_INT_ERROR;
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
#if defined USE_AD_FOR_TL_CALC_UPDATE_OPT_PROPS || \
    defined USE_AD_FOR_TL_CALC_BUILD_PHASE_MATS_SCALAR || \
    defined USE_AD_FOR_TL_CALC_BUILD_PHASE_VECS_VECTOR_GC || \
    defined USE_AD_FOR_TL_CALC_BUILD_PHASE_MATS_VECTOR_GC || \
    defined USE_AD_FOR_TL_CALC_BUILD_LOCAL_R_AND_T || \
    defined USE_AD_FOR_TL_CALC_BUILD_TXR || \
    defined USE_AD_FOR_TL_CALC_BUILD_GAMMA || \
    defined USE_AD_FOR_TL_CALC_EIG_1N_GEN_REAL || \
    defined USE_AD_FOR_TL_CALC_EIG_1N_GEN_COMPLEX || \
    defined USE_AD_FOR_TL_CALC_EIG_1N_TO_2N_REAL || \
    defined USE_AD_FOR_TL_CALC_EIG_1N_TO_2N_COMPLEX || \
    defined USE_AD_FOR_TL_CALC_CALC_SOLVE_BVP || \
    defined USE_AD_FOR_TL_CALC_CALC_GLOBAL_R_AND_T || \
    defined USE_AD_FOR_TL_CALC_GAMMA_SOURCE_VECTORS || \
    defined USE_AD_FOR_TL_CALC_BUILD_GLOBAL_SOURCE || \
    defined USE_AD_FOR_TL_CALC_RTM_EIG_RTS || \
    defined USE_AD_FOR_TL_LAYER_ADD_ALL || \
    defined USE_AD_FOR_TL_SSR_UP_LAYER || \
    defined USE_AD_FOR_TL_SINGLE_SCATTERED_RADIANCE_UP || \
    defined USE_AD_FOR_TL_SSR_DN_LAYER || \
    defined USE_AD_FOR_TL_SINGLE_SCATTERED_RADIANCE_DN || \
    defined USE_AD_FOR_TL_CALC_RADIANCE_SLAB || \
    defined USE_AD_FOR_TL_CALC_GET_LOCAL_R_T_U_W || \
    defined USE_AD_FOR_TL_CALC_GET_LAYER_R_T_S_U_W_V || \
    defined USE_AD_FOR_TL_CALC_GET_TOTAL_TOA_R_S_U_V || \
    defined USE_AD_FOR_TL_CALC_GET_TOTAL_R_T_S_U_W_V || \
    defined USE_AD_FOR_TL_CALC_GET_DIFF_BOUND_INPUT || \
    defined USE_AD_FOR_TL_CALC_GET_SINGLE || \
    defined USE_AD_FOR_TL_CALC_FOURIER_GET_ADD_BOTH_WAYS_A || \
    defined USE_AD_FOR_TL_CALC_GET_SOLUTION_INTERNAL
     if (! (d->options & XRTM_OPTION_CALC_DERIVS && d->options & XRTM_OPTION_REVERSE_DERIVS)) {
          fprintf(stderr, "ERROR: option \"%s\" can only be used with \"%s\"\n", xrtm_option_mask_to_name(XRTM_OPTION_CALC_DERIVS), xrtm_option_mask_to_name(XRTM_OPTION_REVERSE_DERIVS));
          return XRTM_INT_ERROR;
     }
#endif

     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (! (d->options & XRTM_OPTION_CALC_DERIVS))
          d->n_derivs = 0;

     if (   d->options & XRTM_OPTION_PHASE_MATRIX_LC)
          d->options |= XRTM_OPTION_SAVE_PHASE_MATS;


     /*-------------------------------------------------------------------------
      * Value implicitly set to defaults.
      *-----------------------------------------------------------------------*/
     d->phi_0 = 0.;


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     two_n_quad = 2 * d->n_quad;

     if (d->n_coef == 0) {
          d->n_four  = 1;
          d->n_coef2 = 0;
     }
     else
     if (d->n_coef <= two_n_quad) {
          d->n_four  = d->n_coef;
          d->n_coef2 = d->n_coef;
     }
     else {
          d->n_four  = two_n_quad;
          d->n_coef2 = two_n_quad;
     }

     if (d->options & XRTM_OPTION_NO_AZIMUTHAL)
          d->n_four  = 1;

     if (d->n_stokes == 1)
          d->n_elem = 1;
     else
          d->n_elem = 6;

     update_n_four2(d, 0, 0);


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if ((r = setjmp(env)) == 0) {


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     d->set_flags_g     = ALLOC_ARRAY1_UC(d->n_layers);
     d->set_flags_coef  = ALLOC_ARRAY1_UC(d->n_layers);
     d->set_flags_omega = ALLOC_ARRAY1_UC(d->n_layers);
     d->set_flags_ltau  = ALLOC_ARRAY1_UC(d->n_layers);

     if (d->n_kernels > 0) {
          d->set_flags_kernel_ampfac = ALLOC_ARRAY1_UC(d->n_kernels);
          d->set_flags_kernel_params = ALLOC_ARRAY2_UC(d->n_kernels, MAX_KERNEL_PARAMS);
          d->set_flags_kernel_func   = ALLOC_ARRAY1_UC(d->n_kernels);
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     d->dep_flags_Y_p               = ALLOC_ARRAY1_I(d->n_four);
     d->dep_flags_Y_0               = ALLOC_ARRAY1_I(d->n_four);
     d->dep_flags_Y_u               = ALLOC_ARRAY1_I(d->n_four);

     d->dep_flags_opt_props         = ALLOC_ARRAY1_I(d->n_layers);
/*
     d->dep_flags_beam_params0      = ALLOC_ARRAY1_I(d->n_layers);
*/
     d->dep_flags_beam_params       = ALLOC_ARRAY1_I(d->n_layers);
/*
     d->dep_flags_diff_bound_input0 = ALLOC_ARRAY1_I(d->n_four);
*/
     d->dep_flags_diff_bound_input  = ALLOC_ARRAY1_I(d->n_four);

     d->dep_flags_total_R_T_S_U_W_V = ALLOC_ARRAY1_I(d->n_four);

     d->dep_flags_phase_mats_qq     = ALLOC_ARRAY2_I(d->n_four, d->n_layers);
     d->dep_flags_phase_mats_uq     = ALLOC_ARRAY2_I(d->n_four, d->n_layers);
     d->dep_flags_phase_vecs_q0     = ALLOC_ARRAY2_I(d->n_four, d->n_layers);
     d->dep_flags_phase_vecs_u0     = ALLOC_ARRAY2_I(d->n_four, d->n_layers);
     d->dep_flags_local_r_t_u_w     = ALLOC_ARRAY2_I(d->n_four, d->n_layers);
     d->dep_flags_layer_R_T_S_U_W_V = ALLOC_ARRAY2_I(d->n_four, d->n_layers);


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (! (d->options & XRTM_OPTION_SFI))
          d->n_quad_x = d->n_quad + d->n_umus;
     else
          d->n_quad_x = d->n_quad;

     d->n_quad_d   = d->n_quad + d->n_umus;

     d->n_quad_v   = d->n_quad   * d->n_stokes;
     d->n_quad_v_x = d->n_quad_x * d->n_stokes;
     d->n_quad_v_d = d->n_quad_d * d->n_stokes;

     d->n_umus_v   = d->n_umus   * d->n_stokes;


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     d->quad_type = xrtm_quad_to_quad_type(d);

     d->qx = ALLOC_ARRAY1_D(d->n_quad_d);
     d->qw = ALLOC_ARRAY1_D(d->n_quad_d);

     if (d->options & XRTM_OPTION_QUAD_NORM_GAUS_LEG)
          gauss_leg_quad (d->n_quad, d->qx, d->qw);
     else
     if (d->options & XRTM_OPTION_QUAD_DOUB_GAUS_LEG)
          gauss_leg_quad2(d->n_quad, d->qx, d->qw);
/*
          gauss_leg_quadx(d->n_quad, 0., 1., d->qx, d->qw);
*/
     else
     if (d->options & XRTM_OPTION_QUAD_LOBATTO) {
          init_array1_d(d->qx, d->n_quad, -DBL_MAX);
          init_array1_d(d->qw, d->n_quad, -DBL_MAX);
     }
#ifdef DEBUG
     else {
         fprintf(stderr, "ERROR: xrtm_create(): end of if / else if\n");
         exit(1);
     }
#endif
     d->umus = d->qx + d->n_quad;

     d->qf = 0.;
     for (i = 0; i < d->n_quad; ++i)
          d->qf += d->qx[i] * d->qw[i];
     d->qf = .5 / d->qf;


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (d->solvers & XRTM_SOLVERS_INTERNAL) {
          d->qx_v = ALLOC_ARRAY1_D(d->n_quad_v_d);
          d->qw_v = ALLOC_ARRAY1_D(d->n_quad_v_d);

          ii = 0;
          for (i = 0; i < d->n_quad; ++i) {
               jj = ii;
               for (j = 0; j < d->n_stokes; ++j) {
                    d->qx_v[jj] = d->qx[i];
                    d->qw_v[jj] = d->qw[i];
                    jj++;
               }
               ii += d->n_stokes;
          }

          d->umus_v = d->qx_v + d->n_quad_v;
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (d->solvers & XRTM_SOLVERS_INTERNAL && d->misc_input.use_symmetric_form) {
          d->alpha1 = ALLOC_ARRAY1_D(d->n_quad_v_d);
          d->alpha2 = ALLOC_ARRAY1_D(d->n_quad_v_d);

          build_sim_vectors(d->n_quad_v, 0, d->qx_v, d->qw_v, d->alpha1, d->alpha2);
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (d->solvers & XRTM_SOLVERS_USE_COEF) {
          if (! (d->options & XRTM_OPTION_SAVE_LEG_POLYS)) {
               d->Y_p = ALLOC_ARRAY3_D(1, d->n_quad_x, d->n_coef2);
               d->Y_0 = ALLOC_ARRAY2_D(1, d->n_coef2);

               if (d->n_umus > 0)
                    d->Y_u = ALLOC_ARRAY3_D(1, d->n_umus, d->n_coef2);
          }
          else {
               d->Y_p = ALLOC_ARRAY3_D(d->n_four, d->n_quad_x, d->n_coef2);
               d->Y_0 = ALLOC_ARRAY2_D(d->n_four, d->n_coef2);

               if (d->n_umus > 0)
                    d->Y_u = ALLOC_ARRAY3_D(d->n_four, d->n_umus, d->n_coef2);
          }
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (d->solvers & XRTM_SOLVERS_INTERNAL && d->n_kernels > 0) {
          d->kernel_qx = ALLOC_ARRAY1_D(d->n_kernel_quad);
          d->kernel_qw = ALLOC_ARRAY1_D(d->n_kernel_quad);

          gauss_leg_quad2(d->n_kernel_quad, d->kernel_qx, d->kernel_qw);

          for (i = 0; i < d->n_kernel_quad; ++i)
               d->kernel_qx[i] *= PI;

          if (brdf_aux_alloc(&d->kernel_aux, d->n_quad + d->n_umus, d->n_kernel_quad)) {
               fprintf(stderr, "ERROR: brdf_aux_alloc()\n");
               return XRTM_INT_ERROR;
          }
          brdf_aux_calc_base(&d->kernel_aux, d->n_quad, d->qx, d->n_kernel_quad, d->kernel_qx);
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (d->options & XRTM_OPTION_PSA)
          d->levels_z = ALLOC_ARRAY1_D(d->n_layers + 1);

     if (d->options & XRTM_OPTION_OUTPUT_AT_LEVELS) {
          d->ulevels = ALLOC_ARRAY1_I(d->n_ulevels);
          if (d->n_ulevels > 0)
               d->ulevels[0] = 0;
          if (d->n_ulevels > 1)
               d->ulevels[1] = d->n_layers;
     }
     else {
          d->ulevels = ALLOC_ARRAY1_I(d->n_ulevels);
          d->utaus   = ALLOC_ARRAY1_D(d->n_ulevels);
          d->utaus2  = ALLOC_ARRAY1_D(d->n_ulevels);

          if (! (d->options & XRTM_OPTION_DELTA_M))
               d->utaus20 = d->utaus2;
          else
               d->utaus20 = ALLOC_ARRAY1_D(d->n_ulevels);
     }

     if (d->options & XRTM_OPTION_SOURCE_THERMAL)
          d->levels_b = ALLOC_ARRAY1_D(d->n_layers + 1);


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     d->ltau = ALLOC_ARRAY1_D(d->n_layers);

     if (d->options & XRTM_OPTION_OUTPUT_AT_TAUS)
          d->ttau = ALLOC_ARRAY1_D(d->n_layers + 1);

     d->omega = ALLOC_ARRAY1_D(d->n_layers);

     if (d->solvers & XRTM_SOLVERS_USE_G)
          d->g = ALLOC_ARRAY1_D(d->n_layers);
     if (d->solvers & XRTM_SOLVERS_USE_COEF) {
          d->n_coef_layer  = ALLOC_ARRAY1_I(d->n_layers);
          d->n_coef_layer2 = ALLOC_ARRAY1_I(d->n_layers);

          d->coef = ALLOC_ARRAY3_D(d->n_layers, d->n_elem, d->n_coef);
     }

     if (! (d->options & XRTM_OPTION_DELTA_M)) {
          d->g0     = d->g;
          d->coef0  = d->coef;
          d->ltau0  = d->ltau;
          d->omega0 = d->omega;
     }
     else {
          d->ltau0  = ALLOC_ARRAY1_D(d->n_layers);
          d->omega0 = ALLOC_ARRAY1_D(d->n_layers);

          if (d->solvers & XRTM_SOLVERS_USE_G)
               d->g0    = ALLOC_ARRAY1_D(d->n_layers);
          if (d->solvers & XRTM_SOLVERS_USE_COEF)
               d->coef0 = ALLOC_ARRAY3_D(d->n_layers, d->n_elem, d->n_coef);
     }

     if (d->n_kernels > 0) {
          d->kernels = (enum xrtm_kernel_type *) ALLOC_ARRAY1_I(d->n_kernels);
          for (i = 0; i < d->n_kernels; ++i)
               d->kernels[i] = kernels[i];

          d->kernel_ampfac = ALLOC_ARRAY1_D(d->n_kernels);
          d->kernel_params = ALLOC_ARRAY2_D(d->n_kernels, MAX_KERNEL_PARAMS);
          d->kernel_func   = ALLOC_ARRAY1  (d->n_kernels, brdf_kernel_func_data);
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (d->options & XRTM_OPTION_CALC_DERIVS) {
          d->F_iso_top_l = ALLOC_ARRAY1_D(d->n_derivs);
          d->F_iso_bot_l = ALLOC_ARRAY1_D(d->n_derivs);

          if (d->options & XRTM_OPTION_SOURCE_THERMAL) {
               d->levels_b_l  = ALLOC_ARRAY2_D(d->n_layers + 1, d->n_derivs);
               d->surface_b_l = ALLOC_ARRAY1_D(d->n_derivs);
          }

          d->ltau_l  = ALLOC_ARRAY2_D(d->n_layers, d->n_derivs);
          d->omega_l = ALLOC_ARRAY2_D(d->n_layers, d->n_derivs);

          if (d->solvers & XRTM_SOLVERS_USE_G)
               d->g_l    = ALLOC_ARRAY2_D(d->n_layers, d->n_derivs);
          if (d->solvers & XRTM_SOLVERS_USE_COEF)
               d->coef_l = ALLOC_ARRAY4_D(d->n_layers, d->n_derivs, d->n_elem, d->n_coef);

          if (! (d->options & XRTM_OPTION_DELTA_M)) {
               d->g0_l     = d->g_l;
               d->coef0_l  = d->coef_l;
               d->ltau0_l  = d->ltau_l;
               d->omega0_l = d->omega_l;
          }
          else {
               d->ltau0_l  = ALLOC_ARRAY2_D(d->n_layers, d->n_derivs);
               d->omega0_l = ALLOC_ARRAY2_D(d->n_layers, d->n_derivs);

               if (d->solvers & XRTM_SOLVERS_USE_G)
                    d->g0_l    = ALLOC_ARRAY2_D(d->n_layers, d->n_derivs);
               if (d->solvers & XRTM_SOLVERS_USE_COEF)
                    d->coef0_l = ALLOC_ARRAY4_D(d->n_layers, d->n_derivs, d->n_elem, d->n_coef);
          }
     }

     if (d->options & XRTM_OPTION_CALC_DERIVS && d->n_kernels > 0) {
          d->kernel_ampfac_l = ALLOC_ARRAY2_D(d->n_kernels, d->n_derivs);

          d->kernel_params_l = ALLOC_ARRAY3_D(d->n_kernels, d->n_derivs, MAX_KERNEL_PARAMS);
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (d->options & XRTM_OPTION_REVERSE_DERIVS) {
          d->ltau_a  = ALLOC_ARRAY1_D(d->n_layers);
          d->omega_a = ALLOC_ARRAY1_D(d->n_layers);

          if (d->solvers & XRTM_SOLVERS_USE_G)
               d->g_a    = ALLOC_ARRAY1_D(d->n_layers);
          if (d->solvers & XRTM_SOLVERS_USE_COEF)
               d->coef_a = ALLOC_ARRAY3_D(d->n_layers, d->n_elem, d->n_coef);

          if (! (d->options & XRTM_OPTION_DELTA_M)) {
               d->g0_a     = d->g_a;
               d->coef0_a  = d->coef_a;
               d->ltau0_a  = d->ltau_a;
               d->omega0_a = d->omega_a;
          }
          else {
               d->ltau0_a  = ALLOC_ARRAY1_D(d->n_layers);
               d->omega0_a = ALLOC_ARRAY1_D(d->n_layers);

               if (d->solvers & XRTM_SOLVERS_USE_G)
                    d->g0_a    = ALLOC_ARRAY1_D(d->n_layers);
               if (d->solvers & XRTM_SOLVERS_USE_COEF)
                    d->coef0_a = ALLOC_ARRAY3_D(d->n_layers, d->n_elem, d->n_coef);
          }
     }

     if (d->options & XRTM_OPTION_REVERSE_DERIVS && d->n_kernels > 0) {
          d->kernel_ampfac_a = ALLOC_ARRAY1_D(d->n_kernels);

          d->kernel_params_a = ALLOC_ARRAY2_D(d->n_kernels, MAX_KERNEL_PARAMS);
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (d->options & XRTM_OPTION_CALC_DERIVS) {
          d->derivs.layers        = flags_alloc2(d->n_layers + 1, d->n_derivs);
          d->derivs.total         = flags_alloc2(              1, d->n_derivs);
          d->derivs.stacks        = flags_alloc2(d->n_layers + 1, d->n_derivs);
          d->derivs.sources       = flags_alloc2(d->n_layers + 1, d->n_derivs);

          d->derivs.adding_up     = flags_alloc2(d->n_layers + 1, d->n_derivs);
          d->derivs.adding_down   = flags_alloc2(d->n_layers + 1, d->n_derivs);

          d->derivs.beam          = (uchar **) 1;
          if (d->options & XRTM_OPTION_SOURCE_SOLAR)
               d->derivs.beam        = flags_alloc2(d->n_layers + 1, d->n_derivs);

          d->derivs.thermal       = (uchar **) 1;
          if (d->options & XRTM_OPTION_SOURCE_THERMAL)
               d->derivs.thermal     = flags_alloc2(d->n_layers + 1, d->n_derivs);

          d->derivs.layers_m      = flags_alloc2(d->n_layers + 1, d->n_derivs);
          d->derivs.total_m       = flags_alloc2(              1, d->n_derivs);
          d->derivs.stacks_m      = flags_alloc2(d->n_layers + 1, d->n_derivs);
          d->derivs.sources_m     = flags_alloc2(d->n_layers + 1, d->n_derivs);

          d->derivs.adding_up_m   = flags_alloc2(d->n_layers + 1, d->n_derivs);
          d->derivs.adding_down_m = flags_alloc2(d->n_layers + 1, d->n_derivs);

          d->derivs.beam_m        = (uchar **) 1;
          if (d->options & XRTM_OPTION_SOURCE_SOLAR)
               d->derivs.beam_m      = flags_alloc2(d->n_layers + 1, d->n_derivs);

          d->derivs.thermal_m     = (uchar **) 1;
          if (d->options & XRTM_OPTION_SOURCE_THERMAL)
               d->derivs.thermal_m   = flags_alloc2(d->n_layers + 1, d->n_derivs);

          if (0) {
               fprintf(stderr, "ERROR: flags_alloc2()\n");
               return XRTM_INT_ERROR;
          }
     }

     if (d->options & XRTM_OPTION_REVERSE_DERIVS) {
          d->derivs.layers_union        = ALLOC_ARRAY1_UC(d->n_layers + 1);
          d->derivs.total_union         = ALLOC_ARRAY1_UC(d->n_layers + 1);
          d->derivs.stacks_union        = ALLOC_ARRAY1_UC(d->n_layers + 1);
          d->derivs.sources_union       = ALLOC_ARRAY1_UC(d->n_layers + 1);

          if (d->options & XRTM_OPTION_SOURCE_SOLAR)
               d->derivs.beam_union        = ALLOC_ARRAY1_UC(d->n_layers + 1);

          if (d->options & XRTM_OPTION_SOURCE_THERMAL)
               d->derivs.thermal_union     = ALLOC_ARRAY1_UC(d->n_layers + 1);

          d->derivs.adding_up_union     = ALLOC_ARRAY1_UC(d->n_layers + 1);
          d->derivs.adding_down_union   = ALLOC_ARRAY1_UC(d->n_layers + 1);

          d->derivs.layers_m_union      = ALLOC_ARRAY1_UC(d->n_layers + 1);
          d->derivs.total_m_union       = ALLOC_ARRAY1_UC(d->n_layers + 1);
          d->derivs.stacks_m_union      = ALLOC_ARRAY1_UC(d->n_layers + 1);
          d->derivs.sources_m_union     = ALLOC_ARRAY1_UC(d->n_layers + 1);

          if (d->options & XRTM_OPTION_SOURCE_SOLAR)
               d->derivs.beam_m_union      = ALLOC_ARRAY1_UC(d->n_layers + 1);

          if (d->options & XRTM_OPTION_SOURCE_THERMAL)
               d->derivs.thermal_m_union   = ALLOC_ARRAY1_UC(d->n_layers + 1);

          d->derivs.adding_up_m_union   = ALLOC_ARRAY1_UC(d->n_layers + 1);
          d->derivs.adding_down_m_union = ALLOC_ARRAY1_UC(d->n_layers + 1);
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (d->solvers & XRTM_SOLVERS_INTERNAL) {
          d->btau  = ALLOC_ARRAY1_D(d->n_layers + 1);
          d->btran = ALLOC_ARRAY1_D(d->n_layers + 1);
          d->as_0  = ALLOC_ARRAY1_D(d->n_layers);
          d->atran = ALLOC_ARRAY1_D(d->n_layers);
/*
          if (! (d->options & XRTM_OPTION_DELTA_M)) {
               d->btau0  = d->btau;
               d->btran0 = d->btran;
               d->as_00  = d->as_0;
               d->atran0 = d->atran;
          }
          else {
               d->btau0  = ALLOC_ARRAY1_D(d->n_layers + 1);
               d->btran0 = ALLOC_ARRAY1_D(d->n_layers + 1);
               d->as_00  = ALLOC_ARRAY1_D(d->n_layers);
               d->atran0 = ALLOC_ARRAY1_D(d->n_layers);
          }
*/
          d->chapman = ALLOC_ARRAY2_D(d->n_layers, d->n_layers);

          d->In_p    = ALLOC_ARRAY2_D(d->n_four, d->n_quad_v + d->n_umus_v);
          d->I1_m    = ALLOC_ARRAY2_D(d->n_four, d->n_quad_v + d->n_umus_v);

          if (d->n_kernels > 0) {
               d->kernel_f_qq = ALLOC_ARRAY3_D(d->n_kernel_quad * 2 + 1, d->n_quad_x, d->n_quad_x);
               d->kernel_f_q0 = ALLOC_ARRAY2_D(d->n_kernel_quad * 2 + 1, d->n_quad_x);

               if (d->options & XRTM_OPTION_SFI && d->n_umus > 0) {
                    d->kernel_f_uq = ALLOC_ARRAY3_D(d->n_kernel_quad * 2 + 1, d->n_umus, d->n_quad_x);
                    d->kernel_f_u0 = ALLOC_ARRAY2_D(d->n_kernel_quad * 2 + 1, d->n_umus);
               }
          }
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (d->solvers & XRTM_SOLVERS_INTERNAL) {
          if (d->options & XRTM_OPTION_CALC_DERIVS) {
               d->btau_l  = ALLOC_ARRAY2_D(d->n_layers + 1, d->n_derivs);
               d->btran_l = ALLOC_ARRAY2_D(d->n_layers + 1, d->n_derivs);
               d->as_0_l  = ALLOC_ARRAY2_D(d->n_layers, d->n_derivs);
               d->atran_l = ALLOC_ARRAY2_D(d->n_layers, d->n_derivs);

               d->In_p_l  = ALLOC_ARRAY3_D(d->n_four, d->n_derivs, d->n_quad_v + d->n_umus_v);
               d->I1_m_l  = ALLOC_ARRAY3_D(d->n_four, d->n_derivs, d->n_quad_v + d->n_umus_v);
/*
               if (! (d->options & XRTM_OPTION_DELTA_M)) {
                    d->btau0_l  = d->btau_l;
                    d->btran0_l = d->btran_l;
                    d->as_00_l  = d->as_0_l;
                    d->atran0_l = d->atran_l;
               }
               else {
                    d->btau0_l  = ALLOC_ARRAY2_D(d->n_layers + 1, d->n_derivs);
                    d->btran0_l = ALLOC_ARRAY2_D(d->n_layers + 1, d->n_derivs);
                    d->as_00_l  = ALLOC_ARRAY2_D(d->n_layers, d->n_derivs);
                    d->atran0_l = ALLOC_ARRAY2_D(d->n_layers, d->n_derivs);
               }
*/
               if (d->n_kernels > 0) {
                    d->kernel_f_qq_l = ALLOC_ARRAY4_D(d->n_derivs, d->n_kernel_quad * 2 + 1, d->n_quad_x, d->n_quad_x);
                    d->kernel_f_q0_l = ALLOC_ARRAY3_D(d->n_derivs, d->n_kernel_quad * 2 + 1, d->n_quad_x);

                    if (d->options & XRTM_OPTION_SFI && d->n_umus > 0) {
                         d->kernel_f_uq_l = ALLOC_ARRAY4_D(d->n_derivs, d->n_kernel_quad * 2 + 1, d->n_umus, d->n_quad_x);
                         d->kernel_f_u0_l = ALLOC_ARRAY3_D(d->n_derivs, d->n_kernel_quad * 2 + 1, d->n_umus);
                    }
               }
          }
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (d->solvers & XRTM_SOLVERS_INTERNAL) {
          if (d->options & XRTM_OPTION_REVERSE_DERIVS) {
               d->btau_a  = ALLOC_ARRAY1_D(d->n_layers + 1);
               d->btran_a = ALLOC_ARRAY1_D(d->n_layers + 1);
               d->as_0_a  = ALLOC_ARRAY1_D(d->n_layers);
               d->atran_a = ALLOC_ARRAY1_D(d->n_layers);

               d->In_p_a  = ALLOC_ARRAY2_D(d->n_four, d->n_quad_v + d->n_umus_v);
               d->I1_m_a  = ALLOC_ARRAY2_D(d->n_four, d->n_quad_v + d->n_umus_v);
/*
               if (! (d->options & XRTM_OPTION_DELTA_M)) {
                    d->btau0_a  = d->btau_a;
                    d->btran0_a = d->btran_a;
                    d->as_00_a  = d->as_0_a;
                    d->atran0_a = d->atran_a;
               }
               else {
                    d->btau0_a  = ALLOC_ARRAY1_D(d->n_layers + 1);
                    d->btran0_a = ALLOC_ARRAY1_D(d->n_layers + 1);
                    d->as_00_a  = ALLOC_ARRAY1_D(d->n_layers);
                    d->atran0_a = ALLOC_ARRAY1_D(d->n_layers);
               }
*/
               if (d->n_kernels > 0) {
                    d->kernel_f_qq_a = ALLOC_ARRAY3_D(d->n_kernel_quad * 2 + 1, d->n_quad_x, d->n_quad_x);
                    d->kernel_f_q0_a = ALLOC_ARRAY2_D(d->n_kernel_quad * 2 + 1, d->n_quad_x);

                    if (d->options & XRTM_OPTION_SFI && d->n_umus > 0) {
                         d->kernel_f_uq_a = ALLOC_ARRAY3_D(d->n_kernel_quad * 2 + 1, d->n_umus, d->n_quad_x);
                         d->kernel_f_u0_a = ALLOC_ARRAY2_D(d->n_kernel_quad * 2 + 1, d->n_umus);
                    }
               }
          }
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (d->solvers & XRTM_SOLVERS_INTERNAL) {
          if (d->options & XRTM_OPTION_SAVE_PHASE_MATS) {
               d->P_qq_pp = ALLOC_ARRAY4_D(d->n_four, d->n_layers, d->n_quad_v_x, d->n_quad_v_x);
               d->P_qq_mp = ALLOC_ARRAY4_D(d->n_four, d->n_layers, d->n_quad_v_x, d->n_quad_v_x);

               d->P_q0_mm = ALLOC_ARRAY3_D(d->n_four, d->n_layers, d->n_quad_v_x);
               d->P_q0_pm = ALLOC_ARRAY3_D(d->n_four, d->n_layers, d->n_quad_v_x);

               if (d->options & XRTM_OPTION_SFI && d->n_umus > 0) {
                    d->P_uq_pp = ALLOC_ARRAY4_D(d->n_four, d->n_layers, d->n_umus_v, d->n_quad_v_x);
                    d->P_uq_mp = ALLOC_ARRAY4_D(d->n_four, d->n_layers, d->n_umus_v, d->n_quad_v_x);

                    d->P_u0_mm = ALLOC_ARRAY3_D(d->n_four, d->n_layers, d->n_umus_v);
                    d->P_u0_pm = ALLOC_ARRAY3_D(d->n_four, d->n_layers, d->n_umus_v);
               }
          }

          if (d->options & XRTM_OPTION_SAVE_LOCAL_R_T) {
               d->r_p = ALLOC_ARRAY4_D(d->n_four, d->n_layers, d->n_quad_v_x, d->n_quad_v_x);
               d->t_p = ALLOC_ARRAY4_D(d->n_four, d->n_layers, d->n_quad_v_x, d->n_quad_v_x);
          }

          if (d->options & XRTM_OPTION_SAVE_PHASE_MATS && d->options & XRTM_OPTION_VECTOR) {
               d->P_qq_mm = ALLOC_ARRAY4_D(d->n_four, d->n_layers, d->n_quad_v_x, d->n_quad_v_x);
               d->P_qq_pm = ALLOC_ARRAY4_D(d->n_four, d->n_layers, d->n_quad_v_x, d->n_quad_v_x);

               if (d->options & XRTM_OPTION_SFI && d->n_umus > 0) {
                    d->P_uq_mm = ALLOC_ARRAY4_D(d->n_four, d->n_layers, d->n_umus_v, d->n_quad_v_x);
                    d->P_uq_pm = ALLOC_ARRAY4_D(d->n_four, d->n_layers, d->n_umus_v, d->n_quad_v_x);
               }
          }

          if (d->options & XRTM_OPTION_SAVE_LOCAL_R_T && d->options & XRTM_OPTION_VECTOR) {
               d->r_m = ALLOC_ARRAY4_D(d->n_four, d->n_layers, d->n_quad_v_x, d->n_quad_v_x);
               d->t_m = ALLOC_ARRAY4_D(d->n_four, d->n_layers, d->n_quad_v_x, d->n_quad_v_x);
          }

          if (d->options & XRTM_OPTION_CALC_DERIVS && d->options & XRTM_OPTION_SAVE_PHASE_MATS) {
               d->P_qq_pp_l = ALLOC_ARRAY5_D(d->n_four, d->n_layers, d->n_derivs, d->n_quad_v_x, d->n_quad_v_x);
               d->P_qq_mp_l = ALLOC_ARRAY5_D(d->n_four, d->n_layers, d->n_derivs, d->n_quad_v_x, d->n_quad_v_x);

               d->P_q0_mm_l = ALLOC_ARRAY4_D(d->n_four, d->n_layers, d->n_derivs, d->n_quad_v_x);
               d->P_q0_pm_l = ALLOC_ARRAY4_D(d->n_four, d->n_layers, d->n_derivs, d->n_quad_v_x);

               if (d->options & XRTM_OPTION_SFI && d->n_umus > 0) {
                    d->P_uq_pp_l = ALLOC_ARRAY5_D(d->n_four, d->n_layers, d->n_derivs, d->n_umus_v, d->n_quad_v_x);
                    d->P_uq_mp_l = ALLOC_ARRAY5_D(d->n_four, d->n_layers, d->n_derivs, d->n_umus_v, d->n_quad_v_x);

                    d->P_u0_mm_l = ALLOC_ARRAY4_D(d->n_four, d->n_layers, d->n_derivs, d->n_umus_v);
                    d->P_u0_pm_l = ALLOC_ARRAY4_D(d->n_four, d->n_layers, d->n_derivs, d->n_umus_v);
               }
          }

          if (d->options & XRTM_OPTION_CALC_DERIVS && d->options & XRTM_OPTION_SAVE_LOCAL_R_T) {
               d->r_p_l = ALLOC_ARRAY5_D(d->n_four, d->n_layers, d->n_derivs, d->n_quad_v_x, d->n_quad_v_x);
               d->t_p_l = ALLOC_ARRAY5_D(d->n_four, d->n_layers, d->n_derivs, d->n_quad_v_x, d->n_quad_v_x);
          }

          if (d->options & XRTM_OPTION_CALC_DERIVS && d->options & XRTM_OPTION_SAVE_PHASE_MATS &&
              d->options & XRTM_OPTION_VECTOR) {
               d->P_qq_mm_l = ALLOC_ARRAY5_D(d->n_four, d->n_layers, d->n_derivs, d->n_quad_v_x, d->n_quad_v_x);
               d->P_qq_pm_l = ALLOC_ARRAY5_D(d->n_four, d->n_layers, d->n_derivs, d->n_quad_v_x, d->n_quad_v_x);

               if (d->options & XRTM_OPTION_SFI && d->n_umus > 0) {
                    d->P_uq_mm_l = ALLOC_ARRAY5_D(d->n_four, d->n_layers, d->n_derivs, d->n_umus_v, d->n_quad_v_x);
                    d->P_uq_pm_l = ALLOC_ARRAY5_D(d->n_four, d->n_layers, d->n_derivs, d->n_umus_v, d->n_quad_v_x);
               }
          }

          if (d->options & XRTM_OPTION_CALC_DERIVS && d->options & XRTM_OPTION_SAVE_LOCAL_R_T &&
              d->options & XRTM_OPTION_VECTOR) {
               d->r_m_l = ALLOC_ARRAY5_D(d->n_four, d->n_layers, d->n_derivs, d->n_quad_v_x, d->n_quad_v_x);
               d->t_m_l = ALLOC_ARRAY5_D(d->n_four, d->n_layers, d->n_derivs, d->n_quad_v_x, d->n_quad_v_x);
          }
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (d->solvers & (XRTM_SOLVERS_INTERNAL & XRTM_SOLVERS_ADDING)) {
          if (d->options & XRTM_OPTION_STACK_REUSE_ADDING) {
               d->stack_chain = NULL;

               d->stack_grid  = (stack_data ***) ALLOC_ARRAY2(d->n_layers, d->n_layers, stack_data *);
          }

          if (d->options & XRTM_OPTION_SAVE_LAYER_R_T_S) {
               d->Rl_p  = ALLOC_ARRAY4_D(d->n_four, d->n_layers, d->n_quad_v_x, d->n_quad_v_x);
               d->Tl_p  = ALLOC_ARRAY4_D(d->n_four, d->n_layers, d->n_quad_v_x, d->n_quad_v_x);
               d->Rl_m  = ALLOC_ARRAY4_D(d->n_four, d->n_layers, d->n_quad_v_x, d->n_quad_v_x);
               d->Tl_m  = ALLOC_ARRAY4_D(d->n_four, d->n_layers, d->n_quad_v_x, d->n_quad_v_x);

               d->Sl_p  = ALLOC_ARRAY3_D(d->n_four, d->n_layers, d->n_quad_v_x);
               d->Sl_m  = ALLOC_ARRAY3_D(d->n_four, d->n_layers, d->n_quad_v_x);

               d->Sll_p = ALLOC_ARRAY3_D(d->n_four, d->n_layers, d->n_quad_v_x);
               d->Sll_m = ALLOC_ARRAY3_D(d->n_four, d->n_layers, d->n_quad_v_x);
          }

          if (d->options & XRTM_OPTION_SAVE_TOTAL_R_T_S) {
               d->Ra_p  = ALLOC_ARRAY3_D(d->n_four, d->n_quad_v_x, d->n_quad_v_x);
               d->Ta_p  = ALLOC_ARRAY3_D(d->n_four, d->n_quad_v_x, d->n_quad_v_x);
               d->Ra_m  = ALLOC_ARRAY3_D(d->n_four, d->n_quad_v_x, d->n_quad_v_x);
               d->Ta_m  = ALLOC_ARRAY3_D(d->n_four, d->n_quad_v_x, d->n_quad_v_x);

               d->Sa_p  = ALLOC_ARRAY2_D(d->n_four, d->n_quad_v_x);
               d->Sa_m  = ALLOC_ARRAY2_D(d->n_four, d->n_quad_v_x);

               d->Sla_p = ALLOC_ARRAY2_D(d->n_four, d->n_quad_v_x);
               d->Sla_m = ALLOC_ARRAY2_D(d->n_four, d->n_quad_v_x);
          }

          if (d->options & XRTM_OPTION_SAVE_LAYER_R_T_S && d->options & XRTM_OPTION_CALC_DERIVS) {
               d->Rl_p_l  = ALLOC_ARRAY5_D(d->n_four, d->n_layers, d->n_derivs, d->n_quad_v_x, d->n_quad_v_x);
               d->Tl_p_l  = ALLOC_ARRAY5_D(d->n_four, d->n_layers, d->n_derivs, d->n_quad_v_x, d->n_quad_v_x);
               d->Rl_m_l  = ALLOC_ARRAY5_D(d->n_four, d->n_layers, d->n_derivs, d->n_quad_v_x, d->n_quad_v_x);
               d->Tl_m_l  = ALLOC_ARRAY5_D(d->n_four, d->n_layers, d->n_derivs, d->n_quad_v_x, d->n_quad_v_x);

               d->Sl_p_l  = ALLOC_ARRAY4_D(d->n_four, d->n_layers, d->n_derivs, d->n_quad_v_x);
               d->Sl_m_l  = ALLOC_ARRAY4_D(d->n_four, d->n_layers, d->n_derivs, d->n_quad_v_x);

               d->Sll_p_l = ALLOC_ARRAY4_D(d->n_four, d->n_layers, d->n_derivs, d->n_quad_v_x);
               d->Sll_m_l = ALLOC_ARRAY4_D(d->n_four, d->n_layers, d->n_derivs, d->n_quad_v_x);
          }

          if (d->options & XRTM_OPTION_SAVE_TOTAL_R_T_S && d->options & XRTM_OPTION_CALC_DERIVS) {
               d->Ra_p_l  = ALLOC_ARRAY4_D(d->n_four, d->n_derivs, d->n_quad_v_x, d->n_quad_v_x);
               d->Ta_p_l  = ALLOC_ARRAY4_D(d->n_four, d->n_derivs, d->n_quad_v_x, d->n_quad_v_x);
               d->Ra_m_l  = ALLOC_ARRAY4_D(d->n_four, d->n_derivs, d->n_quad_v_x, d->n_quad_v_x);
               d->Ta_m_l  = ALLOC_ARRAY4_D(d->n_four, d->n_derivs, d->n_quad_v_x, d->n_quad_v_x);

               d->Sa_p_l  = ALLOC_ARRAY3_D(d->n_four, d->n_derivs, d->n_quad_v_x);
               d->Sa_m_l  = ALLOC_ARRAY3_D(d->n_four, d->n_derivs, d->n_quad_v_x);

               d->Sla_p_l = ALLOC_ARRAY3_D(d->n_four, d->n_derivs, d->n_quad_v_x);
               d->Sla_m_l = ALLOC_ARRAY3_D(d->n_four, d->n_derivs, d->n_quad_v_x);
          }
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     }
     else if (r < 0)
          return XRTM_INT_ERROR;


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
#ifndef INCLUDE_DEV_SOURCE
     if (d->solvers & (XRTM_SOLVERS_INTERNAL)) {
#else
     if (d->solvers & (XRTM_SOLVERS_INTERNAL | XRTM_SOLVERS_EXT_REQUIRES_WORK)) {
#endif
          if (work_init(&d->work, d->n_quad, d->n_quad_x, d->n_stokes, d->n_derivs, d->n_layers, d->n_umus)) {
               fprintf(stderr, "ERROR: work_init()\n");
               return XRTM_INT_ERROR;
          }
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     d->save_tree.t = NULL;
     if (d->options & XRTM_OPTION_REVERSE_DERIVS)
          save_tree_init(&d->save_tree);


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     d->inputs_initialized = 1;
     d->varied_layers_updated = 1;


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     d->set_flags_doub_d_tau  = 0;
     d->set_flags_pade_params = 0;
     d->set_flags_sos_params  = 0;
     d->set_flags_fourier_tol = 0;
     d->set_flags_lambda      = 0;
     d->set_flags_planet_r    = 0;
     d->set_flags_levels_z    = 0;
     d->set_flags_ulevels     = 0;
     d->set_flags_utaus       = 0;
     d->set_flags_umus        = 0;
     d->set_flags_F_iso_top   = 0;
     d->set_flags_F_iso_bot   = 0;
     d->set_flags_F_0         = 0;
     d->set_flags_mu_0        = 0;
     d->set_flags_phi_0       = 0;
     d->set_flags_levels_b    = 0;
     d->set_flags_surface_b   = 0;

     init_array1_uc(d->set_flags_g,     d->n_layers, 0);
     init_array1_uc(d->set_flags_coef,  d->n_layers, 0);
     init_array1_uc(d->set_flags_omega, d->n_layers, 0);
     init_array1_uc(d->set_flags_ltau,  d->n_layers, 0);

     if (d->n_kernels > 0) {
          init_array1_uc(d->set_flags_kernel_ampfac, d->n_kernels, 0);
          init_array2_uc(d->set_flags_kernel_params, d->n_kernels, MAX_KERNEL_PARAMS, 0);
          init_array1_uc(d->set_flags_kernel_func ,  d->n_kernels, 0);
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     init_deps(d, DEP_MASK_MASK, 0, d->n_layers);


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     d->order_p = -1;
     d->order_0 = -1;
     d->order_u = -1;


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (d->solvers & XRTM_SOLVERS_USE_COEF) {
          init_array1_i(d->n_coef_layer,  d->n_layers, d->n_coef);
          init_array1_i(d->n_coef_layer2, d->n_layers, d->n_coef2);
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (d->options & XRTM_OPTION_CALC_DERIVS) {
          init_array1_d(d->F_iso_top_l, d->n_derivs, 0.);
          init_array1_d(d->F_iso_bot_l, d->n_derivs, 0.);

          if (d->options & XRTM_OPTION_SOURCE_THERMAL) {
               init_array2_d(d->levels_b_l, d->n_layers + 1, d->n_derivs, 0.);
               init_array1_d(d->surface_b_l, d->n_derivs, 0.);

          }
          init_array2_d(d->ltau_l, d->n_layers, d->n_derivs, 0.);
          init_array2_d(d->omega_l, d->n_layers, d->n_derivs, 0.);

          if (d->solvers & XRTM_SOLVERS_USE_COEF)
               init_array4_d(d->coef_l, d->n_layers, d->n_derivs, d->n_elem, d->n_coef, 0.);

          if (d->options & XRTM_OPTION_DELTA_M) {
               init_array2_d(d->ltau0_l, d->n_layers, d->n_derivs, 0.);
               init_array2_d(d->omega0_l, d->n_layers, d->n_derivs, 0.);

               if (d->solvers & XRTM_SOLVERS_USE_COEF)
                    init_array4_d(d->coef0_l, d->n_layers, d->n_derivs, d->n_elem, d->n_coef, 0.);
          }

          if (d->n_kernels > 0) {
               init_array2_d(d->kernel_ampfac_l, d->n_kernels, d->n_derivs, 0.);

               init_array3_d(d->kernel_params_l, d->n_kernels, d->n_derivs, MAX_KERNEL_PARAMS, 0.);
          }
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (d->solvers & XRTM_SOLVERS_INTERNAL) {
          d->btau [0] = 0.;
          d->btran[0] = 1.;
/*
          d->btau0[0] = 0.;
          d->btran0[0] = 1.;
*/
          for (i = 0; i < d->n_derivs; ++i) {
               d->btau_l [0][i] = 0.;
               d->btran_l[0][i] = 0.;
/*
               d->btau0_l[0][i] = 0.;
               d->btran0_l[0][i] = 0.;
*/
          }
     }


     return 0;
}



/*******************************************************************************
 *
 ******************************************************************************/
int xrtm_destroy(xrtm_data *d) {

     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (d->solvers & (XRTM_SOLVERS_INTERNAL & XRTM_SOLVERS_ADDING)) {
          if (d->options & XRTM_OPTION_STACK_REUSE_ADDING) {
              stack_chain_free(d->n_four, d->n_quad_v_x, d->n_derivs, d->n_layers, d->n_stacks, d->stack_chain, d->derivs.layers, d->derivs.stacks, d->derivs.adding_down);
          }
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (d->inputs_initialized) {
          free_array1_uc(d->set_flags_g);
          free_array1_uc(d->set_flags_coef);
          free_array1_uc(d->set_flags_omega);
          free_array1_uc(d->set_flags_ltau);

          if (d->n_kernels > 0) {
               free_array1_uc(d->set_flags_kernel_ampfac);
               free_array2_uc(d->set_flags_kernel_params);
               free_array1_uc(d->set_flags_kernel_func);
          }
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     free_array1_i(d->dep_flags_Y_p);
     free_array1_i(d->dep_flags_Y_0);
     free_array1_i(d->dep_flags_Y_u);

     free_array1_i(d->dep_flags_opt_props);
/*
     free_array1_i(d->dep_flags_beam_params0);
*/
     free_array1_i(d->dep_flags_beam_params);
/*
     free_array1_i(d->dep_flags_diff_bound_input0);
*/
     free_array1_i(d->dep_flags_diff_bound_input);

     free_array1_i(d->dep_flags_total_R_T_S_U_W_V);

     free_array2_i(d->dep_flags_phase_mats_qq);
     free_array2_i(d->dep_flags_phase_mats_uq);
     free_array2_i(d->dep_flags_phase_vecs_q0);
     free_array2_i(d->dep_flags_phase_vecs_u0);
     free_array2_i(d->dep_flags_local_r_t_u_w);
     free_array2_i(d->dep_flags_layer_R_T_S_U_W_V);


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     free_array1_d(d->qx);
     free_array1_d(d->qw);


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (d->solvers & XRTM_SOLVERS_INTERNAL) {
          free_array1_d(d->qx_v);
          free_array1_d(d->qw_v);
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (d->solvers & XRTM_SOLVERS_INTERNAL && d->misc_input.use_symmetric_form) {
          free_array1_d(d->alpha1);
          free_array1_d(d->alpha2);
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (d->solvers & XRTM_SOLVERS_USE_COEF) {
          free_array3_d(d->Y_p);
          free_array2_d(d->Y_0);
          if (d->n_umus > 0)
               free_array3_d(d->Y_u);
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (d->solvers & XRTM_SOLVERS_INTERNAL && d->n_kernels > 0) {
          free_array1_d(d->kernel_qx);
          free_array1_d(d->kernel_qw);

          brdf_aux_free(&d->kernel_aux);
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (d->options & XRTM_OPTION_PSA)
          free_array1_d(d->levels_z);

     if (d->options & XRTM_OPTION_OUTPUT_AT_LEVELS)
          free_array1_i(d->ulevels);
     else {
          free_array1_i(d->ulevels);
          free_array1_d(d->utaus);
          free_array1_d(d->utaus2);

          if (d->options & XRTM_OPTION_DELTA_M)
               free_array1_d(d->utaus20);
     }

     if (d->options & XRTM_OPTION_SOURCE_THERMAL)
          free_array1_d(d->levels_b);

     free_array1_d(d->ltau);

     if (d->options & XRTM_OPTION_OUTPUT_AT_TAUS)
          free_array1_d(d->ttau);

     free_array1_d(d->omega);

     if (d->solvers & XRTM_SOLVERS_USE_G)
          free_array1_d(d->g);
     if (d->solvers & XRTM_SOLVERS_USE_COEF) {
          free_array1_i(d->n_coef_layer);
          free_array1_i(d->n_coef_layer2);

          free_array3_d(d->coef);
     }

     if (d->options & XRTM_OPTION_DELTA_M) {
          free_array1_d(d->ltau0);
          free_array1_d(d->omega0);

          if (d->solvers & XRTM_SOLVERS_USE_G)
               free_array1_d(d->g0);
          if (d->solvers & XRTM_SOLVERS_USE_COEF)
               free_array3_d(d->coef0);
     }

     if (d->n_kernels > 0) {
          free_array1_i((int *) d->kernels);

          free_array1_d(d->kernel_ampfac);
          free_array2_d(d->kernel_params);
          free_array1  (d->kernel_func);
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (d->options & XRTM_OPTION_CALC_DERIVS) {
          free_array1_d(d->F_iso_top_l);
          free_array1_d(d->F_iso_bot_l);

          if (d->options & XRTM_OPTION_SOURCE_THERMAL) {
               free_array2_d(d->levels_b_l);
               free_array1_d(d->surface_b_l);
          }

          free_array2_d(d->ltau_l);
          free_array2_d(d->omega_l);

          if (d->solvers & XRTM_SOLVERS_USE_G)
               free_array2_d(d->g_l);
          if (d->solvers & XRTM_SOLVERS_USE_COEF)
               free_array4_d(d->coef_l);

          if (d->options & XRTM_OPTION_DELTA_M) {
               free_array2_d(d->ltau0_l);
               free_array2_d(d->omega0_l);

               if (d->solvers & XRTM_SOLVERS_USE_G)
                    free_array2_d(d->g0_l);
               if (d->solvers & XRTM_SOLVERS_USE_COEF)
                    free_array4_d(d->coef0_l);
          }
     }

     if (d->options & XRTM_OPTION_CALC_DERIVS && d->n_kernels > 0) {
          free_array2_d(d->kernel_ampfac_l);

          free_array3_d(d->kernel_params_l);
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (d->options & XRTM_OPTION_REVERSE_DERIVS) {
          free_array1_d(d->ltau_a);
          free_array1_d(d->omega_a);

          if (d->solvers & XRTM_SOLVERS_USE_G)
               free_array1_d(d->g_a);
          if (d->solvers & XRTM_SOLVERS_USE_COEF)
               free_array3_d(d->coef_a);

          if (d->options & XRTM_OPTION_DELTA_M) {
               free_array1_d(d->ltau0_a);
               free_array1_d(d->omega0_a);

               if (d->solvers & XRTM_SOLVERS_USE_G)
                    free_array1_d(d->g0_a);
               if (d->solvers & XRTM_SOLVERS_USE_COEF)
                    free_array3_d(d->coef0_a);
          }
     }

     if (d->options & XRTM_OPTION_REVERSE_DERIVS && d->n_kernels > 0) {
          free_array1_d(d->kernel_ampfac_a);

          free_array2_d(d->kernel_params_a);
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (d->options & XRTM_OPTION_CALC_DERIVS) {
          flags_free2(d->derivs.layers);
          flags_free2(d->derivs.total);
          flags_free2(d->derivs.stacks);
          flags_free2(d->derivs.sources);
          flags_free2(d->derivs.adding_down);
          flags_free2(d->derivs.adding_up);

          if (d->options & XRTM_OPTION_SOURCE_SOLAR)
               flags_free2(d->derivs.beam);

          if (d->options & XRTM_OPTION_SOURCE_THERMAL)
               flags_free2(d->derivs.thermal);

          flags_free2(d->derivs.layers_m);
          flags_free2(d->derivs.total_m);
          flags_free2(d->derivs.stacks_m);
          flags_free2(d->derivs.sources_m);
          flags_free2(d->derivs.adding_down_m);
          flags_free2(d->derivs.adding_up_m);

          if (d->options & XRTM_OPTION_SOURCE_SOLAR)
               flags_free2(d->derivs.beam_m);

          if (d->options & XRTM_OPTION_SOURCE_THERMAL)
               flags_free2(d->derivs.thermal_m);
     }

     if (d->options & XRTM_OPTION_REVERSE_DERIVS) {
          free_array1_uc(d->derivs.layers_union);
          free_array1_uc(d->derivs.total_union);
          free_array1_uc(d->derivs.stacks_union);
          free_array1_uc(d->derivs.sources_union);

          if (d->options & XRTM_OPTION_SOURCE_SOLAR)
               free_array1_uc(d->derivs.beam_union);

          if (d->options & XRTM_OPTION_SOURCE_THERMAL)
               free_array1_uc(d->derivs.thermal_union);

          free_array1_uc(d->derivs.adding_up_union);
          free_array1_uc(d->derivs.adding_down_union);

          free_array1_uc(d->derivs.layers_m_union);
          free_array1_uc(d->derivs.total_m_union);
          free_array1_uc(d->derivs.stacks_m_union);
          free_array1_uc(d->derivs.sources_m_union);

          if (d->options & XRTM_OPTION_SOURCE_SOLAR)
               free_array1_uc(d->derivs.beam_m_union);

          if (d->options & XRTM_OPTION_SOURCE_THERMAL)
               free_array1_uc(d->derivs.thermal_m_union);

          free_array1_uc(d->derivs.adding_up_m_union);
          free_array1_uc(d->derivs.adding_down_m_union);
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (d->solvers & XRTM_SOLVERS_INTERNAL) {
          free_array1_d(d->btau);
          free_array1_d(d->btran);
          free_array1_d(d->as_0);
          free_array1_d(d->atran);
/*
          if (d->options & XRTM_OPTION_DELTA_M) {
               free_array1_d(d->btau0);
               free_array1_d(d->btran0);
               free_array1_d(d->as_00);
               free_array1_d(d->atran0);
          }
*/
          free_array2_d(d->chapman);

          free_array2_d(d->In_p);
          free_array2_d(d->I1_m);

          if (d->n_kernels > 0) {
               free_array3_d(d->kernel_f_qq);
               free_array2_d(d->kernel_f_q0);

               if (d->options & XRTM_OPTION_SFI && d->n_umus > 0) {
                    free_array3_d(d->kernel_f_uq);
                    free_array2_d(d->kernel_f_u0);
               }
          }
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (d->solvers & XRTM_SOLVERS_INTERNAL) {
          if (d->options & XRTM_OPTION_CALC_DERIVS) {
               free_array2_d(d->btau_l);
               free_array2_d(d->btran_l);
               free_array2_d(d->as_0_l);
               free_array2_d(d->atran_l);
/*
               if (d->options & XRTM_OPTION_DELTA_M) {
                    free_array2_d(d->btau0_l);
                    free_array2_d(d->btran0_l);
                    free_array2_d(d->as_00_l);
                    free_array2_d(d->atran0_l);
               }
*/
               free_array3_d(d->In_p_l);
               free_array3_d(d->I1_m_l);

               if (d->n_kernels > 0) {
                    free_array4_d(d->kernel_f_qq_l);
                    free_array3_d (d->kernel_f_q0_l);

                    if (d->options & XRTM_OPTION_SFI && d->n_umus > 0) {
                         free_array4_d(d->kernel_f_uq_l);
                         free_array3_d (d->kernel_f_u0_l);
                    }
               }
          }
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (d->solvers & XRTM_SOLVERS_INTERNAL) {
          if (d->options & XRTM_OPTION_REVERSE_DERIVS) {
               free_array1_d(d->btau_a);
               free_array1_d(d->btran_a);
               free_array1_d(d->as_0_a);
               free_array1_d(d->atran_a);
/*
               if (d->options & XRTM_OPTION_DELTA_M) {
                    free_array1_d(d->btau0_a);
                    free_array1_d(d->btran0_a);
                    free_array1_d(d->as_00_a);
                    free_array1_d(d->atran0_a);
               }
*/
               free_array2_d(d->In_p_a);
               free_array2_d(d->I1_m_a);

               if (d->n_kernels > 0) {
                    free_array3_d(d->kernel_f_qq_a);
                    free_array2_d (d->kernel_f_q0_a);

                    if (d->options & XRTM_OPTION_SFI && d->n_umus > 0) {
                         free_array3_d(d->kernel_f_uq_a);
                         free_array2_d (d->kernel_f_u0_a);
                    }
               }
          }
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (d->solvers & XRTM_SOLVERS_INTERNAL) {
          if (d->options & XRTM_OPTION_SAVE_PHASE_MATS) {
               free_array4_d(d->P_qq_pp);
               free_array4_d(d->P_qq_mp);

               free_array3_d (d->P_q0_mm);
               free_array3_d (d->P_q0_pm);

               if (d->options & XRTM_OPTION_SFI && d->n_umus > 0) {
                    free_array4_d(d->P_uq_pp);
                    free_array4_d(d->P_uq_mp);

                    free_array3_d (d->P_u0_mm);
                    free_array3_d (d->P_u0_pm);
               }
          }

          if (d->options & XRTM_OPTION_SAVE_LOCAL_R_T) {
               free_array4_d(d->r_p);
               free_array4_d(d->t_p);
          }

          if (d->options & XRTM_OPTION_SAVE_PHASE_MATS && d->options & XRTM_OPTION_VECTOR) {
               free_array4_d(d->P_qq_mm);
               free_array4_d(d->P_qq_pm);

               if (d->options & XRTM_OPTION_SFI && d->n_umus > 0) {
                    free_array4_d(d->P_uq_mm);
                    free_array4_d(d->P_uq_pm);
               }
          }

          if (d->options & XRTM_OPTION_SAVE_LOCAL_R_T  && d->options & XRTM_OPTION_VECTOR) {
               free_array4_d(d->r_m);
               free_array4_d(d->t_m);
          }

          if (d->options & XRTM_OPTION_CALC_DERIVS && d->options & XRTM_OPTION_SAVE_PHASE_MATS) {
               free_array5_d(d->P_qq_pp_l);
               free_array5_d(d->P_qq_mp_l);

               free_array4_d(d->P_q0_mm_l);
               free_array4_d(d->P_q0_pm_l);

               if (d->options & XRTM_OPTION_SFI && d->n_umus > 0) {
                    free_array5_d(d->P_uq_pp_l);
                    free_array5_d(d->P_uq_mp_l);

                    free_array4_d(d->P_u0_mm_l);
                    free_array4_d(d->P_u0_pm_l);
               }
          }

          if (d->options & XRTM_OPTION_CALC_DERIVS && d->options & XRTM_OPTION_SAVE_LOCAL_R_T) {
               free_array5_d(d->r_p_l);
               free_array5_d(d->t_p_l);
          }

          if (d->options & XRTM_OPTION_CALC_DERIVS && d->options & XRTM_OPTION_SAVE_PHASE_MATS &&
              d->options & XRTM_OPTION_VECTOR) {
               free_array5_d(d->P_qq_mm_l);
               free_array5_d(d->P_qq_pm_l);

               if (d->options & XRTM_OPTION_SFI && d->n_umus > 0) {
                    free_array5_d(d->P_uq_mm_l);
                    free_array5_d(d->P_uq_pm_l);
               }
          }

          if (d->options & XRTM_OPTION_CALC_DERIVS && d->options & XRTM_OPTION_SAVE_LOCAL_R_T &&
              d->options & XRTM_OPTION_VECTOR) {
               free_array5_d(d->r_m_l);
               free_array5_d(d->t_m_l);
          }
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (d->solvers & (XRTM_SOLVERS_INTERNAL & XRTM_SOLVERS_ADDING)) {
          if (d->options & XRTM_OPTION_STACK_REUSE_ADDING) {
               free_array1((void *)  d->stack_chain);
               free_array2((void **) d->stack_grid);
          }

          if (d->options & XRTM_OPTION_SAVE_LAYER_R_T_S) {
               free_array4_d(d->Rl_p);
               free_array4_d(d->Tl_p);
               free_array4_d(d->Rl_m);
               free_array4_d(d->Tl_m);

               free_array3_d (d->Sl_p);
               free_array3_d (d->Sl_m);

               free_array3_d (d->Sll_p);
               free_array3_d (d->Sll_m);
          }

          if (d->options & XRTM_OPTION_SAVE_TOTAL_R_T_S) {
               free_array3_d(d->Ra_p);
               free_array3_d(d->Ta_p);
               free_array3_d(d->Ra_m);
               free_array3_d(d->Ta_m);

               free_array2_d(d->Sa_p);
               free_array2_d(d->Sa_m);

               free_array2_d(d->Sla_p);
               free_array2_d(d->Sla_m);
          }

          if (d->options & XRTM_OPTION_SAVE_LAYER_R_T_S && d->options & XRTM_OPTION_CALC_DERIVS) {
               free_array5_d(d->Rl_p_l);
               free_array5_d(d->Tl_p_l);
               free_array5_d(d->Rl_m_l);
               free_array5_d(d->Tl_m_l);

               free_array4_d(d->Sl_p_l);
               free_array4_d(d->Sl_m_l);

               free_array4_d(d->Sll_p_l);
               free_array4_d(d->Sll_m_l);
          }

          if (d->options & XRTM_OPTION_SAVE_TOTAL_R_T_S && d->options & XRTM_OPTION_CALC_DERIVS) {
               free_array4_d(d->Ra_p_l);
               free_array4_d(d->Ta_p_l);
               free_array4_d(d->Ra_m_l);
               free_array4_d(d->Ta_m_l);

               free_array3_d (d->Sa_p_l);
               free_array3_d (d->Sa_m_l);

               free_array3_d (d->Sla_p_l);
               free_array3_d (d->Sla_m_l);
          }
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (d->solvers & XRTM_SOLVERS_INTERNAL)
          work_free(&d->work);


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (d->options & XRTM_OPTION_REVERSE_DERIVS)
          save_tree_free(&d->save_tree);


     return 0;
}



int xrtm_destroy2(xrtm_data *d) {

     if (xrtm_destroy(d) == XRTM_INT_ERROR) {
          fprintf(stderr, "ERROR: xrtm_destroy()\n");
          return XRTM_INT_ERROR;
     }

     free(d);

     return 0;
}



/*******************************************************************************
 *
 ******************************************************************************/
static int check_inputs_initialized(xrtm_data *d) {

     int i;
     int j;
     int n;

     if (d->solvers & (XRTM_SOLVERS_DOUBLING)) {
          if (! d->set_flags_doub_d_tau) {
               fprintf(stderr, "ERROR: %s not set\n", xrtm_input_index_to_name(INPUT_DOUB_D_TAU));
               return INPUT_DOUB_D_TAU;
          }
     }

     if (d->solvers & XRTM_SOLVER_PADE_ADD) {
          if (! d->set_flags_pade_params) {
               fprintf(stderr, "ERROR: %s not set\n", xrtm_input_index_to_name(INPUT_PADE_PARAMS));
               return INPUT_PADE_PARAMS;
          }
     }

     if (d->solvers & XRTM_SOLVER_SOS) {
          if (! d->set_flags_sos_params) {
               fprintf(stderr, "ERROR: %s not set\n", xrtm_input_index_to_name(INPUT_SOS_PARAMS));
               return INPUT_SOS_PARAMS;
          }
     }

     if (! d->set_flags_fourier_tol) {
          fprintf(stderr, "ERROR: %s not set\n", xrtm_input_index_to_name(INPUT_FOURIER_TOL));
          return INPUT_FOURIER_TOL;
     }
/*
     if (d->options & XRTM_OPTION_SOURCE_THERMAL) {
          if (! d->set_flags_lambda) {
               fprintf(stderr, "ERROR: %s not set\n", xrtm_input_index_to_name(INPUT_LAMBDA));
               return INPUT_LAMBDA;
          }
     }
*/
     if (d->options & XRTM_OPTION_PSA) {
          if (! d->set_flags_planet_r) {
               fprintf(stderr, "ERROR: %s not set\n", xrtm_input_index_to_name(INPUT_PLANET_R));
               return INPUT_PLANET_R;
          }
          if (! d->set_flags_levels_z) {
               fprintf(stderr, "ERROR: %s not set\n", xrtm_input_index_to_name(INPUT_LEVELS_Z));
               return INPUT_LEVELS_Z;
          }
     }

     if (d->options & XRTM_OPTION_OUTPUT_AT_LEVELS && ! d->set_flags_ulevels) {
          fprintf(stderr, "ERROR: %s not set\n", xrtm_input_index_to_name(INPUT_OUT_LEVELS));
          return INPUT_OUT_LEVELS;
     }
     else
     if (d->options & XRTM_OPTION_OUTPUT_AT_TAUS && ! d->set_flags_utaus) {
          fprintf(stderr, "ERROR: %s not set, or if %s has been reset then %s must be reset\n", xrtm_input_index_to_name(INPUT_OUT_TAUS), xrtm_input_index_to_name(INPUT_LTAU), xrtm_input_index_to_name(INPUT_OUT_TAUS));
          return INPUT_OUT_TAUS;
     }

     if (d->n_umus > 0 && ! d->set_flags_umus) {
          fprintf(stderr, "ERROR: %s not set\n", xrtm_input_index_to_name(INPUT_OUT_THETAS));
          return INPUT_OUT_THETAS;
     }

     if (! d->set_flags_F_iso_top) {
          fprintf(stderr, "ERROR: %s not set\n", xrtm_input_index_to_name(INPUT_F_ISO_TOP));
          return INPUT_F_ISO_TOP;
     }

     if (! d->set_flags_F_iso_bot) {
          fprintf(stderr, "ERROR: %s not set\n", xrtm_input_index_to_name(INPUT_F_ISO_BOT));
          return INPUT_F_ISO_BOT;
     }

     if (d->options & XRTM_OPTION_SOURCE_SOLAR) {
          if (! d->set_flags_F_0) {
               fprintf(stderr, "ERROR: %s not set\n", xrtm_input_index_to_name(INPUT_F_0));
               return INPUT_F_0;
          }

          if (! d->set_flags_mu_0) {
               fprintf(stderr, "ERROR: %s not set\n", xrtm_input_index_to_name(INPUT_THETA_0));
               return INPUT_THETA_0;
          }

          if (! d->set_flags_phi_0) {
               fprintf(stderr, "ERROR: %s not set\n", xrtm_input_index_to_name(INPUT_PHI_0));
               return INPUT_PHI_0;
          }
     }

     if (d->options & XRTM_OPTION_SOURCE_THERMAL) {
          if (! d->set_flags_levels_b) {
               fprintf(stderr, "ERROR: %s not set\n", xrtm_input_index_to_name(INPUT_LEVELS_B));
               return INPUT_LEVELS_B;
          }

          if (! d->set_flags_surface_b) {
               fprintf(stderr, "ERROR: %s not set\n", xrtm_input_index_to_name(INPUT_SURFACE_B));
               return INPUT_SURFACE_B;
          }
     }

     for (i = 0; i < d->n_layers; ++i) {
          if (d->solvers & XRTM_SOLVERS_USE_G) {
               if (! d->set_flags_g[i]) {
                    fprintf(stderr, "ERROR: %s not set for layer %d\n", xrtm_input_index_to_name(INPUT_G), i);
                    return INPUT_G;
               }
          }
          if (d->solvers & XRTM_SOLVERS_USE_COEF) {
               if (! d->set_flags_coef[i]) {
                    fprintf(stderr, "ERROR: %s not set for layer %d\n", xrtm_input_index_to_name(INPUT_COEF), i);
                    return INPUT_COEF;
               }
          }
          if (! d->set_flags_omega[i]) {
               fprintf(stderr, "ERROR: %s not set for layer %d\n", xrtm_input_index_to_name(INPUT_OMEGA), i);
               return INPUT_OMEGA;
          }
          if (! d->set_flags_ltau[i]) {
               fprintf(stderr, "ERROR: %s not set for layer %d\n", xrtm_input_index_to_name(INPUT_LTAU ), i);
               return INPUT_LTAU;
          }
     }

     for (i = 0; i < d->n_kernels; ++i) {
          if (! d->set_flags_kernel_ampfac[i]) {
               fprintf(stderr, "ERROR: %s not set for kernel %d\n", xrtm_input_index_to_name(INPUT_KERNEL_AMPFAC), i);
               return INPUT_KERNEL_AMPFAC;
          }

          n = kernel_n_params(d->kernels[i]);

          for (j = 0; j < n; ++j) {
               if (! d->set_flags_kernel_params[i][j]) {
                    fprintf(stderr, "ERROR: %s not set for kernel %d, param %d\n", xrtm_input_index_to_name(INPUT_KERNEL_PARAMS), i, j);
                    return INPUT_KERNEL_PARAMS;
               }
          }

          if (d->kernels[i] == XRTM_KERNEL_USER_DEFINED) {
               if (! d->set_flags_kernel_func[i]) {
                    fprintf(stderr, "ERROR: %s not set for kernel %d\n", xrtm_input_index_to_name(INPUT_KERNEL_FUNCS), i);
                    return INPUT_KERNEL_FUNCS;
               }
          }
     }

     free_array1_uc(d->set_flags_g);
     free_array1_uc(d->set_flags_coef);
     free_array1_uc(d->set_flags_omega);
     free_array1_uc(d->set_flags_ltau);

     if (d->n_kernels > 0) {
          free_array1_uc(d->set_flags_kernel_ampfac);
          free_array2_uc(d->set_flags_kernel_params);
     }

     d->inputs_initialized = 0;

     return 0;
}



static int check_inputs_post_checks(xrtm_data *d) {

     int i;

     double a;

     if (d->options & XRTM_OPTION_OUTPUT_AT_TAUS &&
         d->dep_flags_utaus & (DEP_MASK_UTAUS | DEP_MASK_LTAU)) {

          a = 0;
          for (i = 0; i < d->n_layers; ++i)
               a += d->ltau0[i];

          if (d->utaus[0] < 0. || d->utaus[0] > a) {
               fprintf(stderr, "ERROR: invalid value for utaus[%d]: %e, must be >= zero and <= total optical depth of all layers\n", 0, d->utaus[0]);
               return -1;
          }

          for (i = 1; i < d->n_ulevels; ++i) {
               if (d->utaus[i] <= d->utaus[i-1] || d->utaus[i] > a) {
                    fprintf(stderr, "ERROR: invalid value for utaus[%d]: %e, must be > utaus[%d] and <= total optical depth of all layers\n", i, d->utaus[i], i - 1);
                    return -1;
               }
          }
     }

     return 0;
}



static int solution_preprocess(xrtm_data *d) {

     if (d->inputs_initialized) {
          if (check_inputs_initialized(d)) {
               fprintf(stderr, "ERROR: check_inputs_initialized()\n");
               return -1;
          }
     }

     if (check_inputs_post_checks(d)) {
          fprintf(stderr, "ERROR: check_inputs_post_checks()\n");
          return -1;
     }

     if (d->options & XRTM_OPTION_CALC_DERIVS) {
          if (d->varied_layers_updated) {
               fprintf(stderr, "ERROR: xrtm_update_varied_layers() has not been called\n");
               return -1;
          }
     }

     return 0;
}



/*******************************************************************************
 *
 ******************************************************************************/
static int mu_0_singularity_adjust(xrtm_data *d) {

     int i;
     int n;

     int count;

     double a;

     d->mu_0 = d->mu_0_0;

     n = d->set_flags_umus ? d->n_quad + d->n_umus : d->n_quad;

     count = 0;
     for (i = 0; i < n; ++i) {
          a = fabs(singularity_adjust_auto(d->qx[i], d->mu_0_0,
                   d->misc_input.threshold_mu_0_singlarity));
          if (a != d->mu_0_0) {
               d->mu_0 = a;
               count++;
          }
          if (count > 1) {
               fprintf(stderr, "ERROR: mu_0 singularity adjustment conflict\n");
               return -1;
          }
     }

     if (count > 0)
          set_deps(d, DEP_MASK_MU_0, 0, d->n_layers);

     return 0;
}



/*******************************************************************************
 *
 ******************************************************************************/
int xrtm_get_options(xrtm_data *d) {

     return d->options;
}



/*******************************************************************************
 *
 ******************************************************************************/
int xrtm_get_solvers(xrtm_data *d) {

     return d->solvers;
}



/*******************************************************************************
 *
 ******************************************************************************/
int xrtm_get_max_coef(xrtm_data *d) {

     return d->n_coef;
}



/*******************************************************************************
 *
 ******************************************************************************/
int xrtm_get_n_quad(xrtm_data *d) {

     return d->n_quad;
}



/*******************************************************************************
 *
 ******************************************************************************/
int xrtm_get_n_stokes(xrtm_data *d) {

     return d->n_stokes;
}



/*******************************************************************************
 *
 ******************************************************************************/
int xrtm_get_n_elem(xrtm_data *d) {

     return d->n_elem;
}



/*******************************************************************************
 *
 ******************************************************************************/
int xrtm_get_n_derivs(xrtm_data *d) {

     return d->n_derivs;
}



/*******************************************************************************
 *
 ******************************************************************************/
int xrtm_get_n_layers(xrtm_data *d) {

     return d->n_layers;
}



/*******************************************************************************
 *
 ******************************************************************************/
int xrtm_get_n_kernels(xrtm_data *d) {

     return d->n_kernels;
}



/*******************************************************************************
 *
 ******************************************************************************/
int xrtm_get_n_kernel_quad(xrtm_data *d) {

     return d->n_kernel_quad;
}



/*******************************************************************************
 *
 ******************************************************************************/
enum xrtm_kernel_type xrtm_get_kernel(xrtm_data *d, int i_kernel) {

     if (i_kernel > d->n_kernels - 1) {
          fprintf(stderr, "ERROR: invalid kernel index: %d\n", i_kernel);
          return XRTM_INT_ERROR;
     }

     return d->kernels[i_kernel];
}



/*******************************************************************************
 *
 ******************************************************************************/
int xrtm_get_n_out_levels(xrtm_data *d) {

     return d->n_ulevels;
}



/*******************************************************************************
 *
 ******************************************************************************/
int xrtm_get_n_out_thetas(xrtm_data *d) {

     return d->n_umus;
}



int xrtm_get_n_out_thetas_2(xrtm_data *d) {

     return d->n_out_thetas_2;
}



/*******************************************************************************
 *
 ******************************************************************************/
int xrtm_set_doub_d_tau(xrtm_data *d, double doub_d_tau) {

     if (doub_d_tau <= 0.) {
          fprintf(stderr, "ERROR: invalid value for d_tau: %e, "
                  "must be > zero\n", doub_d_tau);
          return XRTM_INT_ERROR;
     }

     d->doub_d_tau = doub_d_tau;

     d->set_flags_doub_d_tau = 1;

     set_deps(d, DEP_MASK_DOUB_D_TAU, 0, d->n_layers);

     return 0;
}



/*******************************************************************************
 *
 ******************************************************************************/
double xrtm_get_doub_d_tau(xrtm_data *d) {

     if (! (d->solvers & XRTM_SOLVERS_DOUBLING)) {
          fprintf(stderr, "ERROR: model not initialized for doubling solvers\n");
          return XRTM_DBL_ERROR;
     }

     GET_INPUTS_CHECK_SET_FLAGS(d->set_flags_doub_d_tau, xrtm_input_index_to_name(INPUT_DOUB_D_TAU), XRTM_DBL_ERROR);

     return d->doub_d_tau;
}



/*******************************************************************************
 *
 ******************************************************************************/
int xrtm_set_pade_params(xrtm_data *d, int pade_s, int pade_r) {
/*
     if (pade_s < 0) {
          fprintf(stderr, "ERROR: invalid value for s: %d, "
                  "must be >= zero\n", pade_s);
          return XRTM_INT_ERROR;
     }

     if (pade_r <= 0) {
          fprintf(stderr, "ERROR: invalid value for r: %d, "
                  "must be > zero\n", pade_r);
          return XRTM_INT_ERROR;
     }
*/
     d->pade_s = pade_s;
     d->pade_r = pade_r;

     d->set_flags_pade_params = 1;

     set_deps(d, DEP_MASK_PADE_PARAMS, 0, d->n_layers);

     return 0;
}



/*******************************************************************************
 *
 ******************************************************************************/
int xrtm_get_pade_params(xrtm_data *d, int *pade_s, int *pade_r) {

     if (! (d->solvers & XRTM_SOLVER_PADE_ADD)) {
          fprintf(stderr, "ERROR: model not initialized for solver \"pade_add\"\n");
          return XRTM_INT_ERROR;
     }

     GET_INPUTS_CHECK_SET_FLAGS(d->set_flags_pade_params, xrtm_input_index_to_name(INPUT_PADE_PARAMS), XRTM_INT_ERROR);

     *pade_s = d->pade_s;
     *pade_r = d->pade_r;

     return 0;
}



/*******************************************************************************
 *
 ******************************************************************************/
int xrtm_set_sos_params(xrtm_data *d, int sos_max_os, double sos_max_tau, double sos_tol) {

     if (sos_max_os <= 0) {
          fprintf(stderr, "ERROR: invalid value for max_os: %d, "
                  "must be > zero\n", sos_max_os);
          return XRTM_INT_ERROR;
     }

     if (sos_max_tau <= 0.) {
          fprintf(stderr, "ERROR: invalid value for max_tau: %e, "
                  "must be > zero\n", sos_max_tau);
          return XRTM_INT_ERROR;
     }

     if (sos_tol < 0.) {
          fprintf(stderr, "ERROR: invalid value for sos_tol: %e, "
                  "must be >= zero\n", sos_tol);
          return XRTM_INT_ERROR;
     }

     d->sos_max_os  = sos_max_os;
     d->sos_max_tau = sos_max_tau;
     d->sos_tol     = sos_tol;

     d->set_flags_sos_params = 1;
/*
     set_deps(d, DEP_MASK_SOS_PARAMS, 0, d->n_layers);
*/
     return 0;
}



int xrtm_get_sos_params(xrtm_data *d, int *sos_max_os, double *sos_max_tau, double *sos_tol) {

     if (! (d->solvers & XRTM_SOLVER_SOS)) {
          fprintf(stderr, "ERROR: model not initialized for solver \"sos\"\n");
          return XRTM_INT_ERROR;
     }

     GET_INPUTS_CHECK_SET_FLAGS(d->set_flags_sos_params, xrtm_input_index_to_name(INPUT_SOS_PARAMS), XRTM_INT_ERROR);

     *sos_max_os  = d->sos_max_os;
     *sos_max_tau = d->sos_max_tau;
     *sos_tol     = d->sos_tol;

     return 0;
}



/*******************************************************************************
 *
 ******************************************************************************/
int xrtm_set_fourier_tol(xrtm_data *d, double fourier_tol) {

     if (fourier_tol < 0.) {
          fprintf(stderr, "ERROR: invalid value for fourier_tol: %e, "
                  "must be >= zero\n", fourier_tol);
          return XRTM_INT_ERROR;
     }

     d->fourier_tol = fourier_tol;

     d->set_flags_fourier_tol = 1;

     return 0;
}



double xrtm_get_fourier_tol(xrtm_data *d) {

     GET_INPUTS_CHECK_SET_FLAGS(d->set_flags_fourier_tol, xrtm_input_index_to_name(INPUT_FOURIER_TOL), XRTM_DBL_ERROR);

     return d->fourier_tol;
}



/*******************************************************************************
 *
 ******************************************************************************/
int xrtm_set_lambda(xrtm_data *d, double lambda) {

     if (! (d->options & XRTM_OPTION_SOURCE_THERMAL)) {
          fprintf(stderr, "ERROR: model not initialized for thermal sources\n");
          return XRTM_INT_ERROR;
     }

     if (lambda <= 0.) {
          fprintf(stderr, "ERROR: invalid value for lambda: %e, "
                  "must be > zero\n", lambda);
          return XRTM_INT_ERROR;
     }

     d->lambda = lambda;

     d->set_flags_lambda = 1;

     set_deps(d, DEP_MASK_LAMBDA, 0, d->n_layers);

     return 0;
}



double xrtm_get_lambda(xrtm_data *d) {

     if (! (d->options & XRTM_OPTION_SOURCE_THERMAL)) {
          fprintf(stderr, "ERROR: model not initialized for thermal sources\n");
          return XRTM_DBL_ERROR;
     }

     GET_INPUTS_CHECK_SET_FLAGS(d->set_flags_lambda, xrtm_input_index_to_name(INPUT_LAMBDA), XRTM_DBL_ERROR);

     return d->lambda;
}


/*******************************************************************************
 *
 ******************************************************************************/
int xrtm_set_planet_r(xrtm_data *d, double planet_r) {

     if (! (d->options & XRTM_OPTION_PSA)) {
          fprintf(stderr, "ERROR: model not initialized for psa\n");
          return XRTM_INT_ERROR;
     }

     if (planet_r <= 0.) {
          fprintf(stderr, "ERROR: invalid value for planet_r: %e, "
                  "must be > zero\n", planet_r);
          return XRTM_INT_ERROR;
     }

     d->planet_r = planet_r;

     d->set_flags_planet_r = 1;

     set_deps(d, DEP_MASK_PLANET_R, 0, d->n_layers);

     return 0;
}



double xrtm_get_planet_r(xrtm_data *d) {

     if (! (d->options & XRTM_OPTION_PSA)) {
          fprintf(stderr, "ERROR: model not initialized for psa\n");
          return XRTM_DBL_ERROR;
     }

     GET_INPUTS_CHECK_SET_FLAGS(d->set_flags_planet_r, xrtm_input_index_to_name(INPUT_PLANET_R), XRTM_DBL_ERROR);

     return d->planet_r;
}



/*******************************************************************************
 *
 ******************************************************************************/
int xrtm_set_levels_z(xrtm_data *d, double *levels_z) {

     int i;

     if (! (d->options & XRTM_OPTION_PSA)) {
          fprintf(stderr, "ERROR: model not initialized for psa\n");
          return XRTM_INT_ERROR;
     }
/*
     for (i = 0; i < d->n_layers + 1; ++i) {
          if (levels_z[i] + d->planet_r <= 0.) {
               fprintf(stderr, "ERROR: invalid value for levels_z[%d]: %e, "
                       "must be > -planet_r\n", i, levels_z[i]);
               return XRTM_INT_ERROR;
          }
     }
*/
     for (i = 1; i < d->n_layers + 1; ++i) {
          if (levels_z[i] > levels_z[i-1]) {
               fprintf(stderr, "ERROR: levels_z[%d]: %e, is not < levels_z[%d]: %e\n",
                       i, levels_z[i], i-1, levels_z[i-1]);
               return XRTM_INT_ERROR;
          }
     }

     for (i = 0; i < d->n_layers + 1; ++i)
          d->levels_z[i] = levels_z[i];

     d->set_flags_levels_z = 1;

     set_deps(d, DEP_MASK_LEVELS_Z, 0, d->n_layers);

     return 0;
}



int xrtm_get_levels_z(xrtm_data *d, double *levels_z) {

     int i;

     if (! (d->options & XRTM_OPTION_PSA)) {
          fprintf(stderr, "ERROR: model not initialized for psa\n");
          return XRTM_INT_ERROR;
     }

     GET_INPUTS_CHECK_SET_FLAGS(d->set_flags_levels_z, xrtm_input_index_to_name(INPUT_LEVELS_Z), XRTM_INT_ERROR);

     for (i = 0; i < d->n_layers + 1; ++i)
          levels_z[i] = d->levels_z[i];

     return 0;
}



/*******************************************************************************
 *
 ******************************************************************************/
int xrtm_set_out_levels(xrtm_data *d, int *out_levels) {

     int i;

     if (! (d->options & XRTM_OPTION_OUTPUT_AT_LEVELS)) {
          fprintf(stderr, "ERROR: model not initialized for output at levels\n");
          return XRTM_INT_ERROR;
     }

     if (out_levels[0] < 0 || out_levels[0] > d->n_layers) {
          fprintf(stderr, "ERROR: invalid value for out_levels[%d]: %d, must be >= zero and <= n_layers\n", 0, out_levels[0]);
          return XRTM_INT_ERROR;
     }

     for (i = 1; i < d->n_ulevels; ++i) {
          if (out_levels[i] <= out_levels[i-1] || out_levels[i] > d->n_layers) {
               fprintf(stderr, "ERROR: invalid value for out_levels[%d]: %d, must be > out_levels[%d] and <= n_layers\n", i, out_levels[i], i - 1);
               return XRTM_INT_ERROR;
          }
     }

     for (i = 0; i < MIN(2, d->n_ulevels); ++i) {
          if (out_levels[i] != 0 && out_levels[i] != d->n_layers) {
               if (check_solvers(d->solvers, 0, 1, "option", "out_levels[i] != 0 || out_levels[i] != d->n_layers", XRTM_SOLVER_DOUB_ADD, XRTM_SOLVER_EIG_ADD, XRTM_SOLVER_PADE_ADD, XRTM_SOLVER_SOS, XRTM_SOLVER_TWO_OS, 0))
                    return XRTM_INT_ERROR;
          }
     }

     if (d->options & XRTM_OPTION_TOP_DOWN_ADDING) {
          if (out_levels[0] != d->n_layers) {
               fprintf(stderr, "ERROR: out_levels[0] must be equal to n_layers for option %s\n", xrtm_option_mask_to_name(XRTM_OPTION_TOP_DOWN_ADDING));
               return XRTM_INT_ERROR;
          }
     }

     if (d->options & XRTM_OPTION_BOTTOM_UP_ADDING) {
          if (out_levels[0] != 0) {
               fprintf(stderr, "ERROR: out_levels[0] must be equal to 0        for option %s\n", xrtm_option_mask_to_name(XRTM_OPTION_BOTTOM_UP_ADDING));
               return XRTM_INT_ERROR;
          }
     }
#ifdef INCLUDE_DEV_SOURCE
     if (check_dev_solvers_set_ulevels(d, out_levels)) {
          fprintf(stderr, "ERROR: check_dev_solvers_set_ulevels()\n");
          return XRTM_INT_ERROR;
     }
#endif
     for (i = 0; i < d->n_ulevels; ++i)
          d->ulevels[i] = out_levels[i];

     d->set_flags_ulevels = 1;

     set_deps(d, DEP_MASK_ULEVELS, 0, d->n_layers);

     return 0;
}



int xrtm_get_out_levels(xrtm_data *d, int *out_levels) {

     int i;

     if (! (d->options & XRTM_OPTION_OUTPUT_AT_LEVELS)) {
          fprintf(stderr, "ERROR: model not initialized for ulevels output\n");
          return XRTM_INT_ERROR;
     }

     GET_INPUTS_CHECK_SET_FLAGS(d->set_flags_ulevels, xrtm_input_index_to_name(INPUT_OUT_LEVELS), XRTM_INT_ERROR);

     for (i = 0; i < d->n_ulevels; ++i)
          out_levels[i] = d->ulevels[i];

     return 0;
}



/*******************************************************************************
 *
 ******************************************************************************/
int xrtm_set_out_taus(xrtm_data *d, double *out_taus) {

     int i;

     if (! (d->options & XRTM_OPTION_OUTPUT_AT_TAUS)) {
          fprintf(stderr, "ERROR: model not initialized for output at taus\n");
          return XRTM_INT_ERROR;
     }

     for (i = 0; i < d->n_ulevels; ++i) {
          if (out_taus[i] < 0.) {
               fprintf(stderr, "ERROR: invalid value for out_taus[%d]: %e, must be >= zero\n", i, out_taus[i]);
               return XRTM_INT_ERROR;
          }

          d->utaus[i] = out_taus[i];
     }

     d->set_flags_utaus = 1;

     set_deps(d, DEP_MASK_UTAUS, 0, d->n_layers);

     return 0;
}



int xrtm_get_out_taus(xrtm_data *d, double *out_taus) {

     int i;

     if (! (d->options & XRTM_OPTION_OUTPUT_AT_TAUS)) {
          fprintf(stderr, "ERROR: model not initialized for utaus output\n");
          return XRTM_INT_ERROR;
     }

     GET_INPUTS_CHECK_SET_FLAGS(d->set_flags_utaus, xrtm_input_index_to_name(INPUT_OUT_TAUS), XRTM_INT_ERROR);

     for (i = 0; i < d->n_ulevels; ++i)
          out_taus[i] = d->utaus[i];

     return 0;
}



/*******************************************************************************
 *
 ******************************************************************************/
int xrtm_set_out_thetas(xrtm_data *d, double *out_thetas) {

     int i;
     int ii;
     int j;
     int jj;

     if (d->n_umus == 0) {
          fprintf(stderr, "ERROR: model not initialized for user defined output angles (n_out_thetas == 0)\n");
          return XRTM_INT_ERROR;
     }

     for (i = 0; i < d->n_umus; ++i) {
          if (out_thetas[i] < 0. || out_thetas[i] >= 90.) {
               fprintf(stderr, "ERROR: invalid value for out_thetas[%d]: %e, "
                       "must be >= zero and < 90.0\n", i, out_thetas[i]);
               return XRTM_INT_ERROR;
          }
     }

     ii = d->n_quad;
     for (i = 0; i < d->n_umus; ++i) {
          d->qx[ii] = cos(out_thetas[i]*D2R);
          d->qw[ii] = 0.;
          ii++;
     }

     if (d->solvers & XRTM_SOLVERS_INTERNAL) {
          ii = d->n_quad_v;
          for (i = 0; i < d->n_umus; ++i) {
               jj = ii;
               for (j = 0; j < d->n_stokes; ++j) {
                    d->qx_v[jj] = cos(out_thetas[i]*D2R);
                    d->qw_v[jj] = 0.;
                    jj++;
               }
               ii += d->n_stokes;
          }
     }

     d->set_flags_umus = 1;

     update_n_four2(d, d->set_flags_mu_0, 1);

     if (d->set_flags_mu_0) {
          if (mu_0_singularity_adjust(d)) {
               fprintf(stderr, "ERROR: mu_0_singularity_adjust()\n");
               return XRTM_INT_ERROR;
          }
     }

     if (d->solvers & XRTM_SOLVERS_INTERNAL && d->misc_input.use_symmetric_form)
          build_sim_vectors(0, d->n_umus_v, d->qx_v + d->n_quad_v, d->qw_v + d->n_quad_v, d->alpha1 + d->n_quad_v, d->alpha2 + d->n_quad_v);

     if (d->solvers & XRTM_SOLVERS_INTERNAL && d->n_kernels > 0) {
          brdf_aux_calc_part(&d->kernel_aux, d->n_quad, d->n_quad + d->n_umus, d->qx, d->n_kernel_quad, d->kernel_qx);

          if (d->set_flags_mu_0)
               brdf_aux_calc_mu_0(&d->kernel_aux, d->n_quad + d->n_umus, d->mu_0, d->n_kernel_quad, d->kernel_qx);
     }

     set_deps(d, DEP_MASK_UMUS, 0, d->n_layers);

     return 0;
}



int xrtm_get_out_thetas(xrtm_data *d, double *out_thetas) {

     int i;
     int ii;

     if (d->n_umus == 0) {
          fprintf(stderr, "ERROR: model not initialized for user defined output angles (n_out_thetas == 0)\n");
          return XRTM_INT_ERROR;
     }

     GET_INPUTS_CHECK_SET_FLAGS(d->set_flags_umus, xrtm_input_index_to_name(INPUT_OUT_THETAS), XRTM_INT_ERROR);

     ii = d->n_quad;
     for (i = 0; i < d->n_umus; ++i) {
          out_thetas[i] = acos(d->qx[ii])*R2D;
          ii++;
     }

     return 0;
}



/*******************************************************************************
 *
 ******************************************************************************/
int xrtm_set_F_iso_top(xrtm_data *d, double F_iso_top) {

     if (F_iso_top < 0.) {
          fprintf(stderr, "ERROR: invalid value for F_iso_top: %e, "
                  "must be > zero\n", F_iso_top);
          return XRTM_INT_ERROR;
     }

     d->F_iso_top = F_iso_top;

     d->set_flags_F_iso_top = 1;

     set_deps(d, DEP_MASK_F_ISO_TOP, 0, d->n_layers);

     return 0;
}



double xrtm_get_F_iso_top(xrtm_data *d) {

     GET_INPUTS_CHECK_SET_FLAGS(d->set_flags_F_iso_top, xrtm_input_index_to_name(INPUT_F_ISO_TOP), XRTM_DBL_ERROR);

     return d->F_iso_top;
}



int xrtm_set_F_iso_top_l_1(xrtm_data *d, int i_deriv, double F_iso_top_l) {

     CHECK_INIT_FOR_DERIVS();

     CHECK_INDEX_RANGE(i_deriv, d->n_derivs, "i_deriv", "n_derivs", XRTM_INT_ERROR);

     d->F_iso_top_l[i_deriv] = F_iso_top_l;

     set_deps(d, DEP_MASK_F_ISO_TOP_L, 0, d->n_layers);

     return 0;
}



int xrtm_set_F_iso_top_l_n(xrtm_data *d, double *F_iso_top_l) {

     int i;

     CHECK_INIT_FOR_DERIVS();

     for (i = 0; i < d->n_derivs; ++i)
          d->F_iso_top_l[i] = F_iso_top_l[i];

     set_deps(d, DEP_MASK_F_ISO_TOP_L, 0, d->n_layers);

     return 0;
}



double xrtm_get_F_iso_top_l(xrtm_data *d, int i_deriv) {

     CHECK_INDEX_RANGE(i_deriv, d->n_derivs, "i_deriv", "n_derivs", XRTM_DBL_ERROR);

     return d->F_iso_top_l[i_deriv];
}



/*******************************************************************************
 *
 ******************************************************************************/
int xrtm_set_F_iso_bot(xrtm_data *d, double F_iso_bot) {

     if (F_iso_bot < 0.) {
          fprintf(stderr, "ERROR: invalid value for F_iso_bot: %e, "
                  "must be > zero\n", F_iso_bot);
          return XRTM_INT_ERROR;
     }

     d->F_iso_bot = F_iso_bot;

     d->set_flags_F_iso_bot = 1;

     set_deps(d, DEP_MASK_F_ISO_BOT, 0, d->n_layers);

     return 0;
}



double xrtm_get_F_iso_bot(xrtm_data *d) {

     GET_INPUTS_CHECK_SET_FLAGS(d->set_flags_F_iso_bot, xrtm_input_index_to_name(INPUT_F_ISO_BOT), XRTM_DBL_ERROR);

     return d->F_iso_bot;
}



int xrtm_set_F_iso_bot_l_1(xrtm_data *d, int i_deriv, double F_iso_bot_l) {

     CHECK_INIT_FOR_DERIVS();

     CHECK_INDEX_RANGE(i_deriv, d->n_derivs, "i_deriv", "n_derivs", XRTM_INT_ERROR);

     d->F_iso_bot_l[i_deriv] = F_iso_bot_l;

     set_deps(d, DEP_MASK_F_ISO_BOT_L, 0, d->n_layers);

     return 0;
}



int xrtm_set_F_iso_bot_l_n(xrtm_data *d, double *F_iso_bot_l) {

     int i;

     CHECK_INIT_FOR_DERIVS();

     for (i = 0; i < d->n_derivs; ++i)
          d->F_iso_bot_l[i] = F_iso_bot_l[i];

     set_deps(d, DEP_MASK_F_ISO_BOT_L, 0, d->n_layers);

     return 0;
}



double xrtm_get_F_iso_bot_l(xrtm_data *d, int i_deriv) {

     CHECK_INDEX_RANGE(i_deriv, d->n_derivs, "i_deriv", "n_derivs", XRTM_DBL_ERROR);

     return d->F_iso_bot_l[i_deriv];
}



/*******************************************************************************
 *
 ******************************************************************************/
int xrtm_set_F_0(xrtm_data *d, double F_0) {

     if (! (d->options & XRTM_OPTION_SOURCE_SOLAR)) {
          fprintf(stderr, "ERROR: model not initialized for solar sources\n");
          return XRTM_INT_ERROR;
     }

     if (F_0 < 0.) {
          fprintf(stderr, "ERROR: invalid value for F_0: %e, "
                  "must be >= zero\n", F_0);
          return XRTM_INT_ERROR;
     }

     d->F_0 = F_0;

     d->set_flags_F_0 = 1;

     set_deps(d, DEP_MASK_F_0, 0, d->n_layers);

     return 0;
}



double xrtm_get_F_0(xrtm_data *d) {

     if (! (d->options & XRTM_OPTION_SOURCE_SOLAR)) {
          fprintf(stderr, "ERROR: model not initialized for solar sources\n");
          return XRTM_DBL_ERROR;
     }

     GET_INPUTS_CHECK_SET_FLAGS(d->set_flags_F_0, xrtm_input_index_to_name(INPUT_F_0), XRTM_DBL_ERROR);

     return d->F_0;
}



/*******************************************************************************
 *
 ******************************************************************************/
int xrtm_set_theta_0(xrtm_data *d, double theta_0) {

     if (! (d->options & XRTM_OPTION_SOURCE_SOLAR)) {
          fprintf(stderr, "ERROR: model not initialized for solar sources\n");
          return XRTM_INT_ERROR;
     }

     if (theta_0 < 0. || theta_0 >= 90.) {
          fprintf(stderr, "ERROR: invalid value for theta_0: %e, "
                  "must be >= zero and < 90.0\n", theta_0);
          return XRTM_INT_ERROR;
     }

     d->mu_0_0 = cos(theta_0*D2R);

     update_n_four2(d, 1, d->set_flags_umus);

     if (mu_0_singularity_adjust(d)) {
          fprintf(stderr, "ERROR: mu_0_singularity_adjust()\n");
          return XRTM_INT_ERROR;
     }

     if (d->solvers & XRTM_SOLVERS_INTERNAL && d->n_kernels > 0)
          brdf_aux_calc_mu_0(&d->kernel_aux, d->n_quad + d->n_umus, d->mu_0, d->n_kernel_quad, d->kernel_qx);

     d->set_flags_mu_0 = 1;

     set_deps(d, DEP_MASK_MU_0, 0, d->n_layers);

     return 0;
}



double xrtm_get_theta_0(xrtm_data *d) {

     if (! (d->options & XRTM_OPTION_SOURCE_SOLAR)) {
          fprintf(stderr, "ERROR: model not initialized for solar sources\n");
          return XRTM_INT_ERROR;
     }

     GET_INPUTS_CHECK_SET_FLAGS(d->set_flags_mu_0, xrtm_input_index_to_name(INPUT_THETA_0), XRTM_DBL_ERROR);

     return acos(d->mu_0_0)*R2D;
}



/*******************************************************************************
 *
 ******************************************************************************/
int xrtm_set_phi_0(xrtm_data *d, double phi_0) {
/*
     if (! (d->options & XRTM_OPTION_SOURCE_SOLAR)) {
          fprintf(stderr, "ERROR: model not initialized for solar sources\n");
          return XRTM_INT_ERROR;
     }
*/
     d->phi_0 = phi_0;

     set_deps(d, DEP_MASK_PHI_0, 0, d->n_layers);

     d->set_flags_phi_0 = 1;

     return 0;
}



double xrtm_get_phi_0(xrtm_data *d) {
/*
     if (! (d->options & XRTM_OPTION_SOURCE_SOLAR)) {
          fprintf(stderr, "ERROR: model not initialized for solar sources\n");
          return XRTM_INT_ERROR;
     }
*/
     GET_INPUTS_CHECK_SET_FLAGS(d->set_flags_phi_0, xrtm_input_index_to_name(INPUT_PHI_0), XRTM_DBL_ERROR);

     return d->phi_0;
}



/*******************************************************************************
 *
 ******************************************************************************/
int xrtm_set_levels_b(xrtm_data *d, double *levels_b) {

     int i;

     if (! (d->options & XRTM_OPTION_SOURCE_THERMAL)) {
          fprintf(stderr, "ERROR: model not initialized for thermal sources\n");
          return XRTM_INT_ERROR;
     }

     for (i = 0; i < d->n_layers + 1; ++i) {
          if (levels_b[i] < 0.) {
               fprintf(stderr, "ERROR: invalid value for levels_b[%d]: %e, "
                       "must be >= zero\n", i, levels_b[i]);
               return XRTM_INT_ERROR;
          }
     }

     for (i = 0; i < d->n_layers + 1; ++i)
          d->levels_b[i] = levels_b[i];

     d->set_flags_levels_b = 1;

     set_deps(d, DEP_MASK_LEVELS_B, 0, d->n_layers);

     return 0;
}



int xrtm_get_levels_b(xrtm_data *d, double *levels_b) {

     int i;

     if (! (d->options & XRTM_OPTION_SOURCE_THERMAL)) {
          fprintf(stderr, "ERROR: model not initialized for thermal sources\n");
          return XRTM_INT_ERROR;
     }

     GET_INPUTS_CHECK_SET_FLAGS(d->set_flags_levels_b, xrtm_input_index_to_name(INPUT_LEVELS_B), XRTM_INT_ERROR);

     for (i = 0; i < d->n_layers + 1; ++i)
          levels_b[i] = d->levels_b[i];

     return 0;
}



int xrtm_set_levels_b_l_1(xrtm_data *d, int i_deriv, double *levels_b_l) {

     int i;

     if (! (d->options & XRTM_OPTION_SOURCE_THERMAL)) {
          fprintf(stderr, "ERROR: model not initialized for thermal sources\n");
          return XRTM_INT_ERROR;
     }

     CHECK_INIT_FOR_DERIVS();

     CHECK_INDEX_RANGE(i_deriv, d->n_derivs, "i_deriv", "n_derivs", XRTM_INT_ERROR);

     for (i = 0; i < d->n_layers + 1; ++i)
          d->levels_b_l[i][i_deriv] = levels_b_l[i];

     set_deps(d, DEP_MASK_LEVELS_B_L, 0, d->n_layers);

     return 0;
}



int xrtm_set_levels_b_l_n(xrtm_data *d, double **levels_b_l) {

     int i;
     int j;

     if (! (d->options & XRTM_OPTION_SOURCE_THERMAL)) {
          fprintf(stderr, "ERROR: model not initialized for thermal sources\n");
          return XRTM_INT_ERROR;
     }

     CHECK_INIT_FOR_DERIVS();

     for (i = 0; i < d->n_layers + 1; ++i) {
          for (j = 0; j < d->n_derivs; ++j)
               d->levels_b_l[i][j] = levels_b_l[i][j];
     }

     set_deps(d, DEP_MASK_LEVELS_B_L, 0, d->n_layers);

     return 0;
}



double xrtm_get_levels_b_l(xrtm_data *d, int i_deriv, double *levels_b_l) {

     int i;

     if (! (d->options & XRTM_OPTION_SOURCE_THERMAL)) {
          fprintf(stderr, "ERROR: model not initialized for thermal sources\n");
          return XRTM_DBL_ERROR;
     }

     CHECK_INDEX_RANGE(i_deriv, d->n_derivs, "i_deriv", "n_derivs", XRTM_DBL_ERROR);

     for (i = 0; i < d->n_layers + 1; ++i)
          levels_b_l[i] = d->levels_b_l[i][i_deriv];

     return 0;
}



/*******************************************************************************
 *
 ******************************************************************************/
int xrtm_set_surface_b(xrtm_data *d, double surface_b) {

     if (! (d->options & XRTM_OPTION_SOURCE_THERMAL)) {
          fprintf(stderr, "ERROR: model not initialized for thermal sources\n");
          return XRTM_INT_ERROR;
     }

     if (surface_b < 0.) {
          fprintf(stderr, "ERROR: invalid value for surface_b: %e, "
                  "must be > zero\n", surface_b);
          return XRTM_INT_ERROR;
     }

     d->surface_b = surface_b;

     d->set_flags_surface_b = 1;

     set_deps(d, DEP_MASK_SURFACE_B, 0, d->n_layers);

     return 0;
}



double xrtm_get_surface_b(xrtm_data *d) {

     if (! (d->options & XRTM_OPTION_SOURCE_THERMAL)) {
          fprintf(stderr, "ERROR: model not initialized for thermal sources\n");
          return XRTM_INT_ERROR;
     }

     GET_INPUTS_CHECK_SET_FLAGS(d->set_flags_surface_b, xrtm_input_index_to_name(INPUT_SURFACE_B), XRTM_DBL_ERROR);

     return d->surface_b;
}



int xrtm_set_surface_b_l_1(xrtm_data *d, int i_deriv, double surface_b_l) {

     if (! (d->options & XRTM_OPTION_SOURCE_THERMAL)) {
          fprintf(stderr, "ERROR: model not initialized for thermal sources\n");
          return XRTM_INT_ERROR;
     }

     CHECK_INIT_FOR_DERIVS();

     CHECK_INDEX_RANGE(i_deriv, d->n_derivs, "i_deriv", "n_derivs", XRTM_INT_ERROR);

     d->surface_b_l[i_deriv] = surface_b_l;

     set_deps(d, DEP_MASK_SURFACE_B_L, 0, d->n_layers);

     return 0;
}



int xrtm_set_surface_b_l_n(xrtm_data *d, double *surface_b_l) {

     int i;

     if (! (d->options & XRTM_OPTION_SOURCE_THERMAL)) {
          fprintf(stderr, "ERROR: model not initialized for thermal sources\n");
          return XRTM_INT_ERROR;
     }

     CHECK_INIT_FOR_DERIVS();

     for (i = 0; i < d->n_derivs; ++i)
          d->surface_b_l[i] = surface_b_l[i];

     set_deps(d, DEP_MASK_SURFACE_B_L, 0, d->n_layers);

     return 0;
}



double xrtm_get_surface_b_l(xrtm_data *d, int i_deriv) {

     if (! (d->options & XRTM_OPTION_SOURCE_THERMAL)) {
          fprintf(stderr, "ERROR: model not initialized for thermal sources\n");
          return XRTM_INT_ERROR;
     }

     CHECK_INDEX_RANGE(i_deriv, d->n_derivs, "i_deriv", "n_derivs", XRTM_DBL_ERROR);

     return d->surface_b_l[i_deriv];

     return 0;
}



/*******************************************************************************
 *
 ******************************************************************************/
static int xrtm_set_g_x(xrtm_data *d, double *g,
                        int i_layer, int n_layers) {

     int i;

     int deps;

     if (! (d->solvers & XRTM_SOLVERS_USE_G)) {
          fprintf(stderr, "ERROR: no two-stream solvers have been specified\n");
          return -1;
     }

     CHECK_INDEX_RANGE(i_layer, d->n_layers, "i_layer", "n_layers", XRTM_INT_ERROR);

     for (i = i_layer; i < n_layers; ++i) {
          if (g[i - i_layer] < -1. || g[i - i_layer] > 1.) {
               fprintf(stderr, "ERROR: invalid value for g[%d]: %e, must be >= -one and <= one\n", i, g[i - i_layer]);
               return -1;
          }
     }

     if (xrtm_set_layer_x(d, d->g0, g,
                          i_layer, n_layers, d->set_flags_g)) {
          fprintf(stderr, "ERROR: xrtm_set_layer_x()\n");
          return -1;
     }

     deps = DEP_MASK_G;

     set_deps(d, deps, i_layer, n_layers);

     return 0;
}



int xrtm_set_g_1(xrtm_data *d, int i_layer, double g) {

     if (xrtm_set_g_x(d, &g, i_layer, i_layer + 1)) {
          fprintf(stderr, "ERROR: xrtm_set_g_x()\n");
          return XRTM_INT_ERROR;
     }

     return 0;
}



int xrtm_set_g_n(xrtm_data *d, double *g) {

     if (xrtm_set_g_x(d, g, 0, d->n_layers)) {
          fprintf(stderr, "ERROR: xrtm_set_g_x()\n");
          return XRTM_INT_ERROR;
     }

     return 0;
}



double xrtm_get_g(xrtm_data *d, int i_layer) {

     if (! (d->solvers & XRTM_SOLVERS_USE_G)) {
          fprintf(stderr, "ERROR: no two-stream solvers have been specified\n");
          return XRTM_DBL_ERROR;
     }

     GET_INPUTS_CHECK_SET_FLAGS(d->set_flags_g[i_layer], xrtm_input_index_to_name(INPUT_G), XRTM_DBL_ERROR);

     CHECK_INDEX_RANGE(i_layer, d->n_layers, "i_layer", "n_layers", XRTM_DBL_ERROR);

     return d->g0[i_layer];
}



static int xrtm_set_g_x_l(xrtm_data *d, void *g_l, int i_layer,
                          int n_layers, int i_deriv, int n_derivs, int type) {

     int deps;

     CHECK_INIT_FOR_DERIVS();

     if (! (d->solvers & XRTM_SOLVERS_USE_G)) {
          fprintf(stderr, "ERROR: no two-stream solvers have been specified\n");
          return -1;
     }

     CHECK_INDEX_RANGE(i_layer, d->n_layers, "i_layer", "n_layers", XRTM_INT_ERROR);
     CHECK_INDEX_RANGE(i_deriv, d->n_derivs, "i_deriv", "n_derivs", XRTM_INT_ERROR);

     if (xrtm_set_layer_x_l(d, d->g0_l, g_l,
                            i_layer, n_layers, i_deriv, n_derivs, type)) {
          fprintf(stderr, "ERROR: xrtm_set_layer_x_l()\n");
          return -1;
     }

     deps = DEP_MASK_G_L;

     set_deps(d, deps, i_layer, n_layers);

     return 0;
}



int xrtm_set_g_l_11(xrtm_data *d, int i_layer, int i_deriv, double g_l) {

     if (xrtm_set_g_x_l(d, &g_l, i_layer, i_layer + 1, i_deriv, i_deriv + 1, 0)) {
          fprintf(stderr, "ERROR: xrtm_set_g_x_l()\n");
          return XRTM_INT_ERROR;
     }

     return 0;
}



int xrtm_set_g_l_n1(xrtm_data *d, int i_deriv, double *g_l) {

     if (xrtm_set_g_x_l(d, g_l, 0, d->n_layers, i_deriv, i_deriv + 1, 1)) {
          fprintf(stderr, "ERROR: xrtm_set_g_x_l()\n");
          return XRTM_INT_ERROR;
     }

     return 0;
}



int xrtm_set_g_l_1n(xrtm_data *d, int i_layer, double *g_l) {

     if (xrtm_set_g_x_l(d, g_l, i_layer, i_layer + 1, 0, d->n_derivs, 2)) {
          fprintf(stderr, "ERROR: xrtm_set_g_x_l()\n");
          return XRTM_INT_ERROR;
     }

     return 0;
}



int xrtm_set_g_l_nn(xrtm_data *d, double **g_l) {

     if (xrtm_set_g_x_l(d, g_l, 0, d->n_layers, 0, d->n_derivs, 3)) {
          fprintf(stderr, "ERROR: xrtm_set_g_x_l()\n");
          return XRTM_INT_ERROR;
     }

     return 0;
}



double xrtm_get_g_l(xrtm_data *d, int i_layer, int i_deriv) {

     if (! (d->solvers & XRTM_SOLVERS_USE_G)) {
          fprintf(stderr, "ERROR: no two-stream solvers have been specified\n");
          return XRTM_DBL_ERROR;
     }

     CHECK_INDEX_RANGE(i_layer, d->n_layers, "i_layer", "n_layers", XRTM_DBL_ERROR);
     CHECK_INDEX_RANGE(i_deriv, d->n_derivs, "i_deriv", "n_derivs", XRTM_DBL_ERROR);

     return d->g0_l[i_layer][i_deriv];
}



/*******************************************************************************
 *
 ******************************************************************************/
static int xrtm_set_coef_x(xrtm_data *d, double ***coef,
                           int i_layer, int n_layers, int *n_coef_layer) {

     int i;
     int ii;
     int j;
     int k;

     int deps;

     if (! (d->solvers & XRTM_SOLVERS_USE_COEF)) {
          fprintf(stderr, "ERROR: no multi-stream solvers have been specified\n");
          return -1;
     }

     CHECK_INDEX_RANGE(i_layer, d->n_layers, "i_layer", "n_layers", XRTM_INT_ERROR);

     if (n_coef_layer) {
          for (i = i_layer, ii = 0; i < n_layers; ++i, ++ii) {
               if (n_coef_layer[ii] < 0) {
                    fprintf(stderr, "ERROR: n_coef for layer %d must be >= 0\n", i_layer);
                    return -1;
               }

               if (n_coef_layer[ii] > d->n_coef) {
                    fprintf(stderr, "ERROR: n_coef for layer %d > max # of coefs for which model was initiaized (%d > %d)\n", i_layer, n_coef_layer[ii], d->n_coef);
                    return -1;
               }
          }
     }

     for (i = i_layer, ii = 0; i < n_layers; ++i, ++ii) {
          if (! n_coef_layer) {
               d->n_coef_layer [i] = d->n_coef;
               d->n_coef_layer2[i] = d->n_coef2;
          }
          else {
               d->n_coef_layer [i] = n_coef_layer[ii];
               d->n_coef_layer2[i] = MIN(n_coef_layer[ii], d->n_coef2);
          }

          if (d->n_coef_layer[i] > 1) {
               if (fabs(1. - coef[ii][0][0]) > DBL_EPSILON) {
                    fprintf(stderr, "ERROR: the first phase function coefficient for layer %d is not = one: %.16e\n", i, coef[ii][0][0]);
                    return -1;
               }
          }

          for (j = 0; j < d->n_elem; ++j) {
               for (k = 0; k < d->n_coef_layer[i]; ++k) {
                    d->coef0[i][j][k] = coef[ii][j][k];
               }
          }

          if (d->inputs_initialized)
               d->set_flags_coef[i] = 1;
     }

     deps = DEP_MASK_COEF;

     if (d->solvers & XRTM_SOLVERS_INTERNAL && d->options & XRTM_OPTION_DELTA_M) {
          deps |= DEP_MASK_OMEGA;
          deps |= DEP_MASK_LTAU;

          if (d->options & XRTM_OPTION_CALC_DERIVS) {
               deps |= DEP_MASK_COEF_L;
               deps |= DEP_MASK_OMEGA_L;
               deps |= DEP_MASK_LTAU_L;
          }
     }

     set_deps(d, deps, i_layer, n_layers);

     return 0;
}



int xrtm_set_coef_1(xrtm_data *d, int i_layer, int n_coef_layer, double **coef) {

     int *p;

     if (n_coef_layer < 0)
          p = NULL;
     else
          p = &n_coef_layer;

     if (xrtm_set_coef_x(d, &coef, i_layer, i_layer + 1, p)) {
          fprintf(stderr, "ERROR: xrtm_set_coef_x()\n");
          return XRTM_INT_ERROR;
     }

     return 0;
}



int xrtm_set_coef_n(xrtm_data *d, int *n_coef_layer, double ***coef) {

     if (xrtm_set_coef_x(d, coef, 0, d->n_layers, n_coef_layer)) {
          fprintf(stderr, "ERROR: xrtm_set_coef_x()\n");
          return XRTM_INT_ERROR;
     }

     return 0;
}



int xrtm_get_n_coef(xrtm_data *d, int i_layer) {

     if (! (d->solvers & XRTM_SOLVERS_USE_COEF)) {
          fprintf(stderr, "ERROR: no multi-stream solvers have been specified\n");
          return XRTM_INT_ERROR;
     }

     GET_INPUTS_CHECK_SET_FLAGS(d->set_flags_coef[i_layer], "n_coef", XRTM_INT_ERROR);

     CHECK_INDEX_RANGE(i_layer, d->n_layers, "i_layer", "n_layers", XRTM_INT_ERROR);

     return d->n_coef_layer[i_layer];
}



double xrtm_get_coef(xrtm_data *d, int i_layer, int i_elem, int i_coef) {

     if (! (d->solvers & XRTM_SOLVERS_USE_COEF)) {
          fprintf(stderr, "ERROR: no multi-stream solvers have been specified\n");
          return XRTM_DBL_ERROR;
     }

     GET_INPUTS_CHECK_SET_FLAGS(d->set_flags_coef[i_layer], xrtm_input_index_to_name(INPUT_COEF), XRTM_DBL_ERROR);

     CHECK_INDEX_RANGE(i_layer, d->n_layers,              "i_layer", "n_layers",              XRTM_DBL_ERROR);
     CHECK_INDEX_RANGE(i_elem,  d->n_elem,                "i_elem",  "n_elem",                XRTM_DBL_ERROR);
     CHECK_INDEX_RANGE(i_coef,  d->n_coef_layer[i_layer], "i_coef",  "n_coef_layer[i_layer]", XRTM_DBL_ERROR);

     return d->coef0[i_layer][i_elem][i_coef];
}



static int xrtm_set_coef_x_l(xrtm_data *d, void *coef_l, int i_layer, int n_layers,
                             int i_deriv, int n_derivs, int type) {

     int i;
     int ii;
     int j;
     int jj;
     int k;
     int l;

     int deps;

     double a;

     CHECK_INIT_FOR_DERIVS();

     if (! (d->solvers & XRTM_SOLVERS_USE_COEF)) {
          fprintf(stderr, "ERROR: no multi-stream solvers have been specified\n");
          return -1;
     }

     CHECK_INDEX_RANGE(i_layer, d->n_layers, "i_layer", "n_layers", XRTM_INT_ERROR);
     CHECK_INDEX_RANGE(i_deriv, d->n_derivs, "i_deriv", "n_derivs", XRTM_INT_ERROR);

     for (i = i_layer, ii = 0; i < n_layers; ++i, ++ii) {
          for (j = i_deriv, jj = 0; j < n_derivs; ++j, ++jj) {
               for (k = 0; k < d->n_elem; ++k) {
                    for (l = 0; l < d->n_coef_layer[i]; ++l) {
                         switch(type) {
                              case 0:
                                   a = ((double **)   coef_l)[k][l];
                                   break;
                              case 1:
                                   a = ((double ***)  coef_l)[ii][k][l];
                                   break;
                              case 2:
                                   a = ((double ***)  coef_l)[jj][k][l];
                                   break;
                              case 3:
                                   a = ((double ****) coef_l)[ii][jj][k][l];
                                   break;
                              default:
#ifdef DEBUG
                                   fprintf(stderr, "ERROR: xrtm_set_coef_x_l(): end of switch\n");
                                   exit(1);
#endif
                                   break;
                         }

                         d->coef0_l[i][j][k][l] = a;
                    }
               }
          }
     }

     deps = DEP_MASK_COEF_L;

     if (d->solvers & XRTM_SOLVERS_INTERNAL && d->options & XRTM_OPTION_DELTA_M) {
          deps |= DEP_MASK_OMEGA_L;
          deps |= DEP_MASK_LTAU_L;
     }

     set_deps(d, deps, i_layer, n_layers);

     return 0;
}



int xrtm_set_coef_l_11(xrtm_data *d, int i_layer, int i_deriv, double **coef_l) {

     if (xrtm_set_coef_x_l(d, coef_l, i_layer, i_layer + 1, i_deriv, i_deriv + 1, 0)) {
          fprintf(stderr, "ERROR: xrtm_set_coef_x_l()\n");
          return XRTM_INT_ERROR;
     }

     return 0;
}



int xrtm_set_coef_l_n1(xrtm_data *d, int i_deriv, double ***coef_l) {

     if (xrtm_set_coef_x_l(d, coef_l, 0, d->n_layers, i_deriv, i_deriv + 1, 1)) {
          fprintf(stderr, "ERROR: xrtm_set_coef_x_l()\n");
          return XRTM_INT_ERROR;
     }

     return 0;
}



int xrtm_set_coef_l_1n(xrtm_data *d, int i_layer, double ***coef_l) {

     if (xrtm_set_coef_x_l(d, coef_l, i_layer, i_layer + 1, 0, d->n_derivs, 2)) {
          fprintf(stderr, "ERROR: xrtm_set_coef_x_l()\n");
          return XRTM_INT_ERROR;
     }

     return 0;
}



int xrtm_set_coef_l_nn(xrtm_data *d, double ****coef_l) {

     if (xrtm_set_coef_x_l(d, coef_l, 0, d->n_layers, 0, d->n_derivs, 3)) {
          fprintf(stderr, "ERROR: xrtm_set_coef_x_l()\n");
          return XRTM_INT_ERROR;
     }

     return 0;
}



double xrtm_get_coef_l(xrtm_data *d, int i_layer, int i_deriv, int i_elem, int i_coef) {

     if (! (d->solvers & XRTM_SOLVERS_USE_COEF)) {
          fprintf(stderr, "ERROR: no multi-stream solvers have been specified\n");
          return XRTM_DBL_ERROR;
     }

     CHECK_INDEX_RANGE(i_layer, d->n_layers,              "i_layer", "n_layers",              XRTM_DBL_ERROR);
     CHECK_INDEX_RANGE(i_deriv, d->n_derivs,              "i_deriv", "n_derivs",              XRTM_DBL_ERROR);
     CHECK_INDEX_RANGE(i_elem,  d->n_elem,                "i_elem",  "n_elem",                XRTM_DBL_ERROR);
     CHECK_INDEX_RANGE(i_coef,  d->n_coef_layer[i_layer], "i_coef",  "n_coef_layer[i_layer]", XRTM_DBL_ERROR);

     return d->coef0_l[i_layer][i_deriv][i_elem][i_coef];
}



/*******************************************************************************
 *
 ******************************************************************************/
static int xrtm_set_omega_x(xrtm_data *d, double *omega,
                            int i_layer, int n_layers) {

     int i;

     int deps;

     CHECK_INDEX_RANGE(i_layer, d->n_layers, "i_layer", "n_layers", XRTM_INT_ERROR);

     for (i = i_layer; i < n_layers; ++i) {
          if (omega[i - i_layer] < 0. || omega[i - i_layer] > 1.) {
               fprintf(stderr, "ERROR: invalid value for omega[%d]: %e, must be >= zero and <= one\n", i, omega[i - i_layer]);
               return -1;
          }
     }

     if (xrtm_set_layer_x(d, d->omega0, omega,
                          i_layer, n_layers, d->set_flags_omega)) {
          fprintf(stderr, "ERROR: xrtm_set_layer_x()\n");
          return -1;
     }

     for (i = i_layer; i < n_layers; ++i)
          d->omega0[i] = singularity_adjust_down(1., d->omega0[i],
                              d->misc_input.threshold_omega_singlarity);

     deps = DEP_MASK_OMEGA;

     if (d->solvers & XRTM_SOLVERS_INTERNAL && d->options & XRTM_OPTION_DELTA_M) {
          deps |= DEP_MASK_LTAU;

          if (d->options & XRTM_OPTION_CALC_DERIVS) {
               deps |= DEP_MASK_OMEGA_L;
               deps |= DEP_MASK_LTAU_L;
          }
     }

     set_deps(d, deps, i_layer, n_layers);

     return 0;
}



int xrtm_set_omega_1(xrtm_data *d, int i_layer, double omega) {

     if (xrtm_set_omega_x(d, &omega, i_layer, i_layer + 1)) {
          fprintf(stderr, "ERROR: xrtm_set_omega_x()\n");
          return XRTM_INT_ERROR;
     }

     return 0;
}



int xrtm_set_omega_n(xrtm_data *d, double *omega) {

     if (xrtm_set_omega_x(d, omega, 0, d->n_layers)) {
          fprintf(stderr, "ERROR: xrtm_set_omega_x()\n");
          return XRTM_INT_ERROR;
     }

     return 0;
}



double xrtm_get_omega(xrtm_data *d, int i_layer) {

     GET_INPUTS_CHECK_SET_FLAGS(d->set_flags_omega[i_layer], xrtm_input_index_to_name(INPUT_OMEGA), XRTM_DBL_ERROR);

     CHECK_INDEX_RANGE(i_layer, d->n_layers, "i_layer", "n_layers", XRTM_DBL_ERROR);

     return d->omega0[i_layer];
}



static int xrtm_set_omega_x_l(xrtm_data *d, void *omega_l, int i_layer,
                              int n_layers, int i_deriv, int n_derivs, int type) {

     int deps;

     CHECK_INIT_FOR_DERIVS();

     CHECK_INDEX_RANGE(i_layer, d->n_layers, "i_layer", "n_layers", XRTM_INT_ERROR);
     CHECK_INDEX_RANGE(i_deriv, d->n_derivs, "i_deriv", "n_derivs", XRTM_INT_ERROR);

     if (xrtm_set_layer_x_l(d, d->omega0_l, omega_l,
                            i_layer, n_layers, i_deriv, n_derivs, type)) {
          fprintf(stderr, "ERROR: xrtm_set_layer_x_l()\n");
          return -1;
     }

     deps = DEP_MASK_OMEGA_L;

     if (d->options & XRTM_OPTION_DELTA_M)
          deps |= DEP_MASK_LTAU_L;

     set_deps(d, deps, i_layer, n_layers);

     return 0;
}



int xrtm_set_omega_l_11(xrtm_data *d, int i_layer, int i_deriv, double omega_l) {

     if (xrtm_set_omega_x_l(d, &omega_l, i_layer, i_layer + 1, i_deriv, i_deriv + 1, 0)) {
          fprintf(stderr, "ERROR: xrtm_set_omega_x_l()\n");
          return XRTM_INT_ERROR;
     }

     return 0;
}



int xrtm_set_omega_l_n1(xrtm_data *d, int i_deriv, double *omega_l) {

     if (xrtm_set_omega_x_l(d, omega_l, 0, d->n_layers, i_deriv, i_deriv + 1, 1)) {
          fprintf(stderr, "ERROR: xrtm_set_omega_x_l()\n");
          return XRTM_INT_ERROR;
     }

     return 0;
}



int xrtm_set_omega_l_1n(xrtm_data *d, int i_layer, double *omega_l) {

     if (xrtm_set_omega_x_l(d, omega_l, i_layer, i_layer + 1, 0, d->n_derivs, 2)) {
          fprintf(stderr, "ERROR: xrtm_set_omega_x_l()\n");
          return XRTM_INT_ERROR;
     }

     return 0;
}



int xrtm_set_omega_l_nn(xrtm_data *d, double **omega_l) {

     if (xrtm_set_omega_x_l(d, omega_l, 0, d->n_layers, 0, d->n_derivs, 3)) {
          fprintf(stderr, "ERROR: xrtm_set_omega_x_l()\n");
          return XRTM_INT_ERROR;
     }

     return 0;
}



double xrtm_get_omega_l(xrtm_data *d, int i_layer, int i_deriv) {

     CHECK_INDEX_RANGE(i_layer, d->n_layers, "i_layer", "n_layers", XRTM_DBL_ERROR);
     CHECK_INDEX_RANGE(i_deriv, d->n_derivs, "i_deriv", "n_derivs", XRTM_DBL_ERROR);

     return d->omega0_l[i_layer][i_deriv];
}



/*******************************************************************************
 *
 ******************************************************************************/
static int xrtm_set_ltau_x(xrtm_data *d, double *ltau,
                           int i_layer, int n_layers) {

     int i;

     CHECK_INDEX_RANGE(i_layer, d->n_layers, "i_layer", "n_layers", XRTM_INT_ERROR);

     for (i = i_layer; i < n_layers; ++i) {
          if (ltau[i - i_layer] <= 0.) {
               fprintf(stderr, "ERROR: invalid value for ltau[%d]: %e, must be > zero\n", i, ltau[i - i_layer]);
               return -1;
          }
     }

     if (xrtm_set_layer_x(d, d->ltau0, ltau,
                          i_layer, n_layers, d->set_flags_ltau)) {
          fprintf(stderr, "ERROR: xrtm_set_layer_x()\n");
          return -1;
     }

     set_deps(d, DEP_MASK_LTAU, i_layer, n_layers);

     return 0;
}



int xrtm_set_ltau_1(xrtm_data *d, int i_layer, double ltau) {

     if (xrtm_set_ltau_x(d, &ltau, i_layer, i_layer + 1)) {
          fprintf(stderr, "ERROR: xrtm_set_ltau_x()\n");
          return XRTM_INT_ERROR;
     }

     return 0;
}



int xrtm_set_ltau_n(xrtm_data *d, double *ltau) {

     if (xrtm_set_ltau_x(d, ltau, 0, d->n_layers)) {
          fprintf(stderr, "ERROR: xrtm_set_ltau_x()\n");
          return XRTM_INT_ERROR;
     }

     return 0;
}



double xrtm_get_ltau(xrtm_data *d, int i_layer) {

     GET_INPUTS_CHECK_SET_FLAGS(d->set_flags_ltau[i_layer], xrtm_input_index_to_name(INPUT_LTAU), XRTM_DBL_ERROR);

     CHECK_INDEX_RANGE(i_layer, d->n_layers, "i_layer", "n_layers", XRTM_DBL_ERROR);

     return d->ltau0[i_layer];
}



static int xrtm_set_ltau_x_l(xrtm_data *d, void *ltau_l, int i_layer,
                             int n_layers, int i_deriv, int n_derivs, int type) {

     CHECK_INIT_FOR_DERIVS();

     CHECK_INDEX_RANGE(i_layer, d->n_layers, "i_layer", "n_layers", XRTM_INT_ERROR);
     CHECK_INDEX_RANGE(i_deriv, d->n_derivs, "i_deriv", "n_derivs", XRTM_INT_ERROR);

     if (xrtm_set_layer_x_l(d, d->ltau0_l, ltau_l,
                            i_layer, n_layers, i_deriv, n_derivs, type)) {
          fprintf(stderr, "ERROR: xrtm_set_layer_x_l()\n");
          return -1;
     }

     set_deps(d, DEP_MASK_LTAU_L, i_layer, n_layers);

     return 0;
}



int xrtm_set_ltau_l_11(xrtm_data *d, int i_layer, int i_deriv, double ltau_l) {

     if (xrtm_set_ltau_x_l(d, &ltau_l, i_layer, i_layer + 1, i_deriv, i_deriv + 1, 0)) {
          fprintf(stderr, "ERROR: xrtm_set_ltau_x_l()\n");
          return XRTM_INT_ERROR;
     }

     return 0;
}



int xrtm_set_ltau_l_n1(xrtm_data *d, int i_deriv, double *ltau_l) {

     if (xrtm_set_ltau_x_l(d, ltau_l, 0, d->n_layers, i_deriv, i_deriv + 1, 1)) {
          fprintf(stderr, "ERROR: xrtm_set_ltau_x_l()\n");
          return XRTM_INT_ERROR;
     }

     return 0;
}



int xrtm_set_ltau_l_1n(xrtm_data *d, int i_layer, double *ltau_l) {

     if (xrtm_set_ltau_x_l(d, ltau_l, i_layer, i_layer + 1, 0, d->n_derivs, 2)) {
          fprintf(stderr, "ERROR: xrtm_set_ltau_x_l()\n");
          return XRTM_INT_ERROR;
     }

     return 0;
}



int xrtm_set_ltau_l_nn(xrtm_data *d, double **ltau_l) {

     if (xrtm_set_ltau_x_l(d, ltau_l, 0, d->n_layers, 0, d->n_derivs, 3)) {
          fprintf(stderr, "ERROR: xrtm_set_ltau_x_l()\n");
          return XRTM_INT_ERROR;
     }

     return 0;
}



double xrtm_get_ltau_l(xrtm_data *d, int i_layer, int i_deriv) {

     CHECK_INDEX_RANGE(i_layer, d->n_layers, "i_layer", "n_layers", XRTM_DBL_ERROR);
     CHECK_INDEX_RANGE(i_deriv, d->n_derivs, "i_deriv", "n_derivs", XRTM_DBL_ERROR);

     return d->ltau0_l[i_layer][i_deriv];
}



/*******************************************************************************
 *
 ******************************************************************************/
int xrtm_set_kernel_ampfac(xrtm_data *d, int i_kernel, double ampfac) {

     CHECK_N_KERNELS_NOT_ZERO(d, XRTM_INT_ERROR);

     CHECK_INDEX_RANGE(i_kernel, d->n_kernels, "i_kernel", "n_kernels", XRTM_INT_ERROR);

     if (ampfac < 0. || ampfac > 1.) {
          fprintf(stderr, "ERROR: invalid value for ampfac: %e, "
                   "must be >= zero and <= one\n", ampfac);
          return XRTM_INT_ERROR;
     }

     d->kernel_ampfac[i_kernel] = ampfac;

     if (d->inputs_initialized)
          d->set_flags_kernel_ampfac[i_kernel] = 1;

     set_deps(d, DEP_MASK_BRDF, 0, d->n_layers);

     return 0;
}



double xrtm_get_kernel_ampfac(xrtm_data *d, int i_kernel) {

     CHECK_N_KERNELS_NOT_ZERO(d, XRTM_DBL_ERROR);

     GET_INPUTS_CHECK_SET_FLAGS(d->set_flags_kernel_ampfac[i_kernel], xrtm_input_index_to_name(INPUT_KERNEL_AMPFAC),  XRTM_DBL_ERROR);

     CHECK_INDEX_RANGE(i_kernel, d->n_kernels, "i_kernel", "n_kernels", XRTM_DBL_ERROR);

     return d->kernel_ampfac[i_kernel];
}



int xrtm_set_kernel_params_1(xrtm_data *d, int i_kernel, int i_param, double param) {

     int n_params;

     CHECK_N_KERNELS_NOT_ZERO(d, XRTM_INT_ERROR);

     CHECK_INDEX_RANGE(i_kernel, d->n_kernels, "i_kernel", "n_kernels", XRTM_INT_ERROR);

     n_params = kernel_n_params((enum xrtm_kernel_type) d->kernels[i_kernel]);

     CHECK_INDEX_RANGE(i_param,     n_params,  "i_param",  "n_params", XRTM_INT_ERROR);

     d->kernel_params[i_kernel][i_param] = param;

     if (d->inputs_initialized)
          d->set_flags_kernel_params[i_kernel][i_param] = 1;

     set_deps(d, DEP_MASK_BRDF, 0, d->n_layers);

     return 0;
}



int xrtm_set_kernel_params_n(xrtm_data *d, int i_kernel, double *params) {

     int i;

     int n_params;

     CHECK_N_KERNELS_NOT_ZERO(d, XRTM_INT_ERROR);

     CHECK_INDEX_RANGE(i_kernel, d->n_kernels, "i_kernel", "n_kernels", XRTM_INT_ERROR);

     n_params = kernel_n_params((enum xrtm_kernel_type) d->kernels[i_kernel]);

     for (i = 0; i < n_params; ++i) {
          d->kernel_params[i_kernel][i] = params[i];

          if (d->inputs_initialized)
               d->set_flags_kernel_params[i_kernel][i] = 1;
     }

     set_deps(d, DEP_MASK_BRDF, 0, d->n_layers);

     return 0;
}



double xrtm_get_kernel_params(xrtm_data *d, int i_kernel, int i_param) {

     int n_params;

     CHECK_N_KERNELS_NOT_ZERO(d, XRTM_DBL_ERROR);

     GET_INPUTS_CHECK_SET_FLAGS(d->set_flags_kernel_params[i_kernel][i_param], xrtm_input_index_to_name(INPUT_KERNEL_PARAMS), XRTM_DBL_ERROR);

     CHECK_INDEX_RANGE(i_kernel, d->n_kernels, "i_kernel", "n_kernels", XRTM_DBL_ERROR);

     n_params = kernel_n_params((enum xrtm_kernel_type) d->kernels[i_kernel]);

     CHECK_INDEX_RANGE(i_param,     n_params,  "i_param",  "n_params",  XRTM_DBL_ERROR);

     return d->kernel_params[i_kernel][i_param];
}



int xrtm_set_kernel_ampfac_l_1(xrtm_data *d, int i_kernel, int i_deriv, double ampfac_l) {

     CHECK_N_KERNELS_NOT_ZERO(d, XRTM_INT_ERROR);

     CHECK_INIT_FOR_DERIVS();

     CHECK_INDEX_RANGE(i_kernel, d->n_kernels, "i_kernel", "n_kernels", XRTM_INT_ERROR);
     CHECK_INDEX_RANGE(i_deriv,  d->n_derivs,  "i_deriv",  "n_derivs",  XRTM_INT_ERROR);

     d->kernel_ampfac_l[i_kernel][i_deriv] = ampfac_l;

     set_deps(d, DEP_MASK_BRDF_L, 0, d->n_layers);

     return 0;
}



int xrtm_set_kernel_ampfac_l_n(xrtm_data *d, int i_kernel, double *ampfac_l) {

     int i;

     CHECK_N_KERNELS_NOT_ZERO(d, XRTM_INT_ERROR);

     CHECK_INIT_FOR_DERIVS();

     CHECK_INDEX_RANGE(i_kernel, d->n_kernels, "i_kernel", "n_kernels", XRTM_INT_ERROR);

     for (i = 0; i < d->n_derivs; ++i)
          d->kernel_ampfac_l[i_kernel][i] = ampfac_l[i];

     set_deps(d, DEP_MASK_BRDF_L, 0, d->n_layers);

     return 0;
}



double xrtm_get_kernel_ampfac_l(xrtm_data *d, int i_kernel, int i_deriv) {

     CHECK_N_KERNELS_NOT_ZERO(d, XRTM_DBL_ERROR);

     CHECK_INDEX_RANGE(i_kernel, d->n_kernels, "i_kernel", "n_kernels", XRTM_DBL_ERROR);
     CHECK_INDEX_RANGE(i_deriv,  d->n_derivs,  "i_deriv",  "n_derivs",  XRTM_DBL_ERROR);

     return d->kernel_ampfac_l[i_kernel][i_deriv];
}



int xrtm_set_kernel_params_l_11(xrtm_data *d, int i_kernel, int i_deriv, int i_param, double param_l) {

     int n_params;

     CHECK_N_KERNELS_NOT_ZERO(d, XRTM_INT_ERROR);

     CHECK_INIT_FOR_DERIVS();

     CHECK_INDEX_RANGE(i_kernel, d->n_kernels, "i_kernel", "n_kernels", XRTM_INT_ERROR);

     n_params = kernel_n_params((enum xrtm_kernel_type) d->kernels[i_kernel]);

     CHECK_INDEX_RANGE(i_deriv,  d->n_derivs,  "i_deriv",  "n_derivs",  XRTM_INT_ERROR);
     CHECK_INDEX_RANGE(i_param,     n_params,  "i_param",  "n_params",  XRTM_INT_ERROR);

     d->kernel_params_l[i_kernel][i_deriv][i_param] = param_l;

     set_deps(d, DEP_MASK_BRDF_L, 0, d->n_layers);

     return 0;
}



int xrtm_set_kernel_params_l_1n(xrtm_data *d, int i_kernel, int i_deriv, double *params_l) {

     int i;

     int n_params;

     CHECK_N_KERNELS_NOT_ZERO(d, XRTM_INT_ERROR);

     CHECK_INIT_FOR_DERIVS();

     CHECK_INDEX_RANGE(i_kernel, d->n_kernels, "i_kernel", "n_kernels", XRTM_INT_ERROR);
     CHECK_INDEX_RANGE(i_deriv,  d->n_derivs,  "i_deriv",  "n_derivs",  XRTM_INT_ERROR);

     n_params = kernel_n_params((enum xrtm_kernel_type) d->kernels[i_kernel]);

     for (i = 0; i < n_params; ++i)
          d->kernel_params_l[i_kernel][i_deriv][i] = params_l[i];

     set_deps(d, DEP_MASK_BRDF_L, 0, d->n_layers);

     return 0;
}



int xrtm_set_kernel_params_l_n1(xrtm_data *d, int i_kernel, int i_param, double *params_l) {

     int i;

     int n_params;

     CHECK_N_KERNELS_NOT_ZERO(d, XRTM_INT_ERROR);

     CHECK_INIT_FOR_DERIVS();

     CHECK_INDEX_RANGE(i_kernel, d->n_kernels, "i_kernel", "n_kernels", XRTM_INT_ERROR);

     n_params = kernel_n_params((enum xrtm_kernel_type) d->kernels[i_kernel]);

     CHECK_INDEX_RANGE(i_param,     n_params,  "i_param",  "n_params",  XRTM_INT_ERROR);

     for (i = 0; i < d->n_derivs; ++i)
          d->kernel_params_l[i_kernel][i][i_param] = params_l[i];

     set_deps(d, DEP_MASK_BRDF_L, 0, d->n_layers);

     return 0;
}



int xrtm_set_kernel_params_l_nn(xrtm_data *d, int i_kernel, double **params_l) {

     int i;
     int j;

     int n_params;

     CHECK_N_KERNELS_NOT_ZERO(d, XRTM_INT_ERROR);

     CHECK_INIT_FOR_DERIVS();

     CHECK_INDEX_RANGE(i_kernel, d->n_kernels, "i_kernel", "n_kernels", XRTM_INT_ERROR);

     n_params = kernel_n_params((enum xrtm_kernel_type) d->kernels[i_kernel]);

     for (i = 0; i < d->n_derivs; ++i) {
          for (j = 0; j < n_params; ++j) {
               d->kernel_params_l[i_kernel][i][j] = params_l[i][j];
          }
     }

     set_deps(d, DEP_MASK_BRDF_L, 0, d->n_layers);

     return 0;
}



double xrtm_get_kernel_params_l(xrtm_data *d, int i_kernel, int i_deriv, int i_param) {

     int n_params;

     CHECK_N_KERNELS_NOT_ZERO(d, XRTM_DBL_ERROR);

     CHECK_INDEX_RANGE(i_kernel, d->n_kernels, "i_kernel", "n_kernels", XRTM_DBL_ERROR);

     n_params = kernel_n_params((enum xrtm_kernel_type) d->kernels[i_kernel]);

     CHECK_INDEX_RANGE(i_deriv,  d->n_derivs,  "i_deriv",  "n_derivs",  XRTM_DBL_ERROR);
     CHECK_INDEX_RANGE(i_param,     n_params,  "i_param",  "n_params",  XRTM_DBL_ERROR);

     return d->kernel_params_l[i_kernel][i_deriv][i_param];
}



int xrtm_set_kernel_func(xrtm_data *d, int i_kernel, void (*aux)(brdf_aux_data *aux, void *aux2, int i, int j, int n_derivs, uchar *derivs, double *p, double **p_l),
                         double (*kernel)(brdf_aux_data *aux, void *aux2, enum xrtm_brdf_geometry brdf_geom, int i, int j, int k, double phi, int n_derivs, uchar *derivs, double *p, double **p_l, double *f_l)) {

     CHECK_N_KERNELS_NOT_ZERO(d, XRTM_INT_ERROR);

     CHECK_INDEX_RANGE(i_kernel, d->n_kernels, "i_kernel", "n_kernels", XRTM_INT_ERROR);

     if (d->kernels[i_kernel] != XRTM_KERNEL_USER_DEFINED) {
          fprintf(stderr, "ERROR: kernel %d is not a user defined kernel\n", i_kernel);
          return XRTM_INT_ERROR;
     }

     d->kernel_func[i_kernel].aux    = aux;
     d->kernel_func[i_kernel].kernel = kernel;

     if (d->inputs_initialized)
          d->set_flags_kernel_func[i_kernel] = 1;

     set_deps(d, DEP_MASK_BRDF, 0, d->n_layers);

     return 0;
}



void *xrtm_get_kernel_func(xrtm_data *d, int i_kernel) {

     return NULL;
}



/*******************************************************************************
 *
 ******************************************************************************/
void xrtm_clear_saved(xrtm_data *d) {

     init_deps(d, DEP_MASK_MASK, 0, d->n_layers);
}



/*******************************************************************************
 *
 ******************************************************************************/
static void build_derivs_layers(xrtm_data *d, uchar **derivs_layers, int flag) {

     int i;
     int j;
     int k;
     int l;

     int n_params;

     for (i = 0; i < d->n_layers; ++i) {
          for (j = 0; j < d->n_derivs; ++j) {
               derivs_layers[i][j] = 0;

               if (d->ltau0_l[i][j] != 0.) {
                    derivs_layers[i][j] = 1;
                    goto L1;
               }
               if (d->omega0_l[i][j] != 0.) {
                    derivs_layers[i][j] = 1;
                    goto L1;
               }

               for (k = 0; k < d->n_elem; ++k) {
                    for (l = 0; l < d->n_coef_layer[i]; ++l) {
                         if (d->coef0_l[i][j][k][l] != 0.) {
                              derivs_layers[i][j] = 1;
                              goto L1;
                         }
                    }
               }

               derivs_layers[i][j] = 0;

L1:            continue;
          }
     }

     for (j = 0; j < d->n_derivs; ++j) {
          derivs_layers[i][j] = 0;

          if (flag || brdf_needs_fourier(d->kernels, d->n_kernels)) {
               for (k = 0; k < d->n_kernels; ++k) {
                    if (d->kernel_ampfac_l[k][j] != 0.) {
                         derivs_layers[i][j] = 1;
                         goto L2;
                    }

                    n_params = kernel_n_params((enum xrtm_kernel_type) d->kernels[k]);

                    for (l = 0; l < n_params; ++l) {
                         if (d->kernel_params_l[k][j][l] != 0.) {
                              derivs_layers[i][j] = 1;
                              goto L2;
                         }
                    }
               }
          }

L2:       continue;
     }

     flags_fill_meta2(derivs_layers, d->n_layers + 1, d->n_derivs);
}



static void build_derivs_total(xrtm_data *d, uchar **derivs_layers, uchar **derivs_sources, uchar **derivs_total) {

     int i;
     int j;

     for (j = 0; j < d->n_derivs; ++j)
          derivs_total[0][j] = 0;

     for (i = 0; i < d->n_layers + 1; ++i) {
          for (j = 0; j < d->n_derivs; ++j) {
               if (derivs_layers[i][j] || derivs_sources[i][j])
                    derivs_total[0][j] = 1;
          }
     }

     flags_fill_meta2(derivs_total, 1, d->n_derivs);
}



static void build_derivs_beam(xrtm_data *d, uchar **derivs_layers, uchar **derivs_beam) {

     int i;
     int j;

     int flag;

     for (j = 0; j < d->n_derivs; ++j) {
          flag = 0;

          for (i = 0; i <  d->n_layers + 1; ++i) {
               if (derivs_layers[i][j])
                    flag = 1;

               derivs_beam[i][j] = 0;
               if (flag)
                    derivs_beam[i][j] = 1;
          }
     }

     flags_fill_meta2(derivs_beam, d->n_layers + 1, d->n_derivs);
}



static void build_derivs_thermal(xrtm_data *d, uchar **derivs_layers, uchar **derivs_thermal) {

     int i;
     int j;

     for (j = 0; j < d->n_derivs; ++j) {
          for (i = 0; i < d->n_layers; ++i) {
               derivs_thermal[i][j] = 0;

               if (derivs_layers[i][j] ||
                   d->levels_b_l[i][j] != 0. || d->levels_b_l[i + 1][j] != 0.)
                    derivs_thermal[i][j] = 1;
          }

          derivs_thermal[i][j] = 0;
          if (derivs_layers[i][j] || d->surface_b_l[j] != 0.)
               derivs_thermal[i][j] = 1;
     }


     flags_fill_meta2(derivs_thermal, d->n_layers + 1, d->n_derivs);
}



static void build_derivs_sources(xrtm_data *d, uchar **derivs_layers, uchar **derivs_beam, uchar **derivs_thermal, uchar **derivs_sources) {

     int i;
     int j;

     for (j = 0; j < d->n_derivs; ++j) {
          for (i = 0; i < d->n_layers + 1; ++i) {
               derivs_sources[i][j] = 0;

               if (d->options & XRTM_OPTION_SOURCE_SOLAR)
                    derivs_sources[i][j] |= derivs_beam[i][j];
               if (d->options & XRTM_OPTION_SOURCE_THERMAL)
                    derivs_sources[i][j] |= derivs_thermal[i][j];
          }
     }

     flags_fill_meta2(derivs_sources, d->n_layers + 1, d->n_derivs);
}



static void build_derivs_adding_down(xrtm_data *d, uchar **derivs_layers, uchar **derivs_sources, uchar **derivs_adding_down) {

     int i;
     int j;

     int flag;
     int flag2;

     for (j = 0; j < d->n_derivs; ++j) {
          flag = 0;
          flag2 = 0;

          derivs_adding_down[0][j] = 255;

          if (! derivs_layers[0][j] && ! derivs_sources[0][j])
               derivs_adding_down[0][j] = ADDING_U_U;
          if (! derivs_layers[0][j] &&   derivs_sources[0][j])
               derivs_adding_down[0][j] = ADDING_U_S;
          if (  derivs_layers[0][j] && ! derivs_sources[0][j])
               derivs_adding_down[0][j] = ADDING_U_L;
          if (  derivs_layers[0][j] &&   derivs_sources[0][j])
               derivs_adding_down[0][j] = ADDING_U_B;

          if (derivs_layers[0][j])
               flag = 1;
          if (derivs_sources[0][j])
               flag2 = 1;

          for (i = 1; i < d->n_layers + 1; ++i) {
               derivs_adding_down[i][j] = 255;

               if (! flag && ! flag2 && ! derivs_layers[i][j] && ! derivs_sources[i][j])
                    derivs_adding_down[i][j] = ADDING_U_U;
               else
               if (! flag && ! flag2 && ! derivs_layers[i][j] &&   derivs_sources[i][j])
                    derivs_adding_down[i][j] = ADDING_U_S;
               else
               if (! flag && ! flag2 &&   derivs_layers[i][j] && ! derivs_sources[i][j])
                    derivs_adding_down[i][j] = ADDING_U_L;
               else
               if (! flag && ! flag2 &&   derivs_layers[i][j] &&   derivs_sources[i][j])
                    derivs_adding_down[i][j] = ADDING_U_B;
               else


               if (! flag &&   flag2 && ! derivs_layers[i][j] && ! derivs_sources[i][j])
                    derivs_adding_down[i][j] = ADDING_S_U;
               else
               if (! flag &&   flag2 && ! derivs_layers[i][j] &&   derivs_sources[i][j])
                    derivs_adding_down[i][j] = ADDING_S_S;
               else
               if (! flag &&   flag2 &&   derivs_layers[i][j] && ! derivs_sources[i][j])
                    derivs_adding_down[i][j] = ADDING_S_L;
               else
               if (! flag &&   flag2 &&   derivs_layers[i][j] &&   derivs_sources[i][j])
                    derivs_adding_down[i][j] = ADDING_S_B;
               else


               if (  flag && ! flag2 && ! derivs_layers[i][j] && ! derivs_sources[i][j])
                    derivs_adding_down[i][j] = ADDING_L_U;
               else
               if (  flag && ! flag2 && ! derivs_layers[i][j] &&   derivs_sources[i][j])
                    derivs_adding_down[i][j] = ADDING_L_S;
               else
               if (  flag && ! flag2 &&   derivs_layers[i][j] && ! derivs_sources[i][j])
                    derivs_adding_down[i][j] = ADDING_L_L;
               else
               if (  flag && ! flag2 &&   derivs_layers[i][j] &&   derivs_sources[i][j])
                    derivs_adding_down[i][j] = ADDING_L_B;
               else


               if (  flag &&   flag2 && ! derivs_layers[i][j] && ! derivs_sources[i][j])
                    derivs_adding_down[i][j] = ADDING_B_U;
               else
               if (  flag &&   flag2 && ! derivs_layers[i][j] &&   derivs_sources[i][j])
                    derivs_adding_down[i][j] = ADDING_B_S;
               else
               if (  flag &&   flag2 &&   derivs_layers[i][j] && ! derivs_sources[i][j])
                    derivs_adding_down[i][j] = ADDING_B_L;
               else
               if (  flag &&   flag2 &&   derivs_layers[i][j])
                    derivs_adding_down[i][j] = ADDING_B_B;
#ifdef DEBUG
               else {
                    fprintf(stderr, "ERROR: invalid downward layer adding combination\n");
                    exit(1);
               }
#endif
               if (derivs_layers[i][j])
                    flag = 1;
               if (derivs_sources[i][j])
                    flag2 = 1;
          }
     }

     flags_fill_meta2(derivs_adding_down, d->n_layers + 1, d->n_derivs);
}



static void build_derivs_adding_up(xrtm_data *d, uchar **derivs_layers, uchar **derivs_sources, uchar **derivs_adding_up) {

     int i;
     int j;

     int flag;
     int flag2;

     for (j = 0; j < d->n_derivs; ++j) {
          flag  = 0;
          flag2 = 0;

          derivs_adding_up[d->n_layers][j] = 255;

          if (! derivs_layers[d->n_layers][j] && ! derivs_sources[d->n_layers][j])
               derivs_adding_up[d->n_layers][j] = ADDING_U_U;
          if (! derivs_layers[d->n_layers][j] &&   derivs_sources[d->n_layers][j])
               derivs_adding_up[d->n_layers][j] = ADDING_S_U;
          if (  derivs_layers[d->n_layers][j] && ! derivs_sources[d->n_layers][j])
               derivs_adding_up[d->n_layers][j] = ADDING_L_U;
          if (  derivs_layers[d->n_layers][j] &&   derivs_sources[d->n_layers][j])
               derivs_adding_up[d->n_layers][j] = ADDING_B_U;

          if (derivs_layers[d->n_layers][j])
               flag = 1;

          if (derivs_sources[d->n_layers][j])
               flag2 = 1;

          for (i = d->n_layers - 1; i >= 0; --i) {
               derivs_adding_up[i][j] = 255;

               if (! derivs_layers[i][j] && ! derivs_sources[i][j] && ! flag && ! flag2)
                    derivs_adding_up[i][j] = ADDING_U_U;
               else
               if (! derivs_layers[i][j] &&   derivs_sources[i][j] && ! flag && ! flag2)
                    derivs_adding_up[i][j] = ADDING_S_U;
               else
               if (  derivs_layers[i][j] && ! derivs_sources[i][j] && ! flag && ! flag2)
                    derivs_adding_up[i][j] = ADDING_L_U;
               else
               if (  derivs_layers[i][j] &&   derivs_sources[i][j] && ! flag && ! flag2)
                    derivs_adding_up[i][j] = ADDING_B_U;
               else


               if (! derivs_layers[i][j] && ! derivs_sources[i][j] && ! flag &&   flag2)
                    derivs_adding_up[i][j] = ADDING_U_S;
               else
               if (! derivs_layers[i][j] &&   derivs_sources[i][j] && ! flag &&   flag2)
                    derivs_adding_up[i][j] = ADDING_S_S;
               else
               if (  derivs_layers[i][j] && ! derivs_sources[i][j] && ! flag &&   flag2)
                    derivs_adding_up[i][j] = ADDING_L_S;
               else
               if (  derivs_layers[i][j] &&   derivs_sources[i][j] && ! flag &&   flag2)
                    derivs_adding_up[i][j] = ADDING_B_S;
               else


               if (! derivs_layers[i][j] && ! derivs_sources[i][j] &&   flag && ! flag2)
                    derivs_adding_up[i][j] = ADDING_U_L;
               else
               if (! derivs_layers[i][j] &&   derivs_sources[i][j] &&   flag && ! flag2)
                    derivs_adding_up[i][j] = ADDING_S_L;
               else
               if (  derivs_layers[i][j] && ! derivs_sources[i][j] &&   flag && ! flag2)
                    derivs_adding_up[i][j] = ADDING_L_L;
               else
               if (  derivs_layers[i][j] &&   derivs_sources[i][j] &&   flag && ! flag2)
                    derivs_adding_up[i][j] = ADDING_B_L;
               else


               if (! derivs_layers[i][j] && ! derivs_sources[i][j] &&   flag &&   flag2)
                    derivs_adding_up[i][j] = ADDING_U_B;
               else
               if (! derivs_layers[i][j] &&   derivs_sources[i][j] &&   flag &&   flag2)
                    derivs_adding_up[i][j] = ADDING_S_B;
               else
               if (  derivs_layers[i][j] && ! derivs_sources[i][j] &&   flag &&   flag2)
                    derivs_adding_up[i][j] = ADDING_L_B;
               else
               if (  derivs_layers[i][j] &&   derivs_sources[i][j] &&   flag &&   flag2)
                    derivs_adding_up[i][j] = ADDING_B_B;
#ifdef DEBUG
               else {
                    fprintf(stderr, "ERROR: invalid upward layer adding combination\n");
                    exit(1);
               }
#endif
               if (derivs_layers[i][j])
                    flag = 1;

               if (derivs_sources[i][j])
                    flag2 = 1;
          }
     }

     flags_fill_meta2(derivs_adding_up, d->n_layers + 1, d->n_derivs);
}



static void build_derivs_stacks(xrtm_data *d, uchar **derivs_layers, uchar **derivs_stacks) {

     int i;
     int j;

     int flag;

     for (j = 0; j < d->n_derivs; ++j) {
          flag = 0;

          derivs_stacks[0][j] = ADDING_U_U;
          if (derivs_layers[0][j])
               derivs_stacks[0][j] = ADDING_U_S;

          if (derivs_layers[0][j])
               flag = 1;

          for (i = 1; i <  d->n_layers + 1; ++i) {
               if (! derivs_layers[i-1][j] && ! flag)
                    derivs_stacks[i][j] = ADDING_U_U;
               else
               if (  derivs_layers[i-1][j] && ! flag)
                    derivs_stacks[i][j] = ADDING_U_S;
               else
               if (! derivs_layers[i-1][j] &&  flag)
                    derivs_stacks[i][j] = ADDING_S_S;
               else
               if (  derivs_layers[i-1][j] &&  flag)
                    derivs_stacks[i][j] = ADDING_S_S;
#ifdef DEBUG
               else {
                    fprintf(stderr, "ERROR: invalid layer adding combination\n");
                    exit(1);
               }
#endif
               if (derivs_layers[i][j])
                    flag = 1;

          }
     }

     flags_fill_meta2(derivs_stacks, d->n_layers + 1, d->n_derivs);
}



/*******************************************************************************
 *
 ******************************************************************************/
int xrtm_update_varied_layers(xrtm_data *d) {

     int i;
     int j;

     int n_stacks;

     size_t size;

     stack_data *stack_chain;

     if (! (d->options & XRTM_OPTION_CALC_DERIVS)) {
          fprintf(stderr, "ERROR: model not initialized for linearized values\n");
          return XRTM_INT_ERROR;
     }

     if (d->inputs_initialized) {
          if (check_inputs_initialized(d)) {
               fprintf(stderr, "ERROR: check_inputs_initialized()\n");
               return XRTM_INT_ERROR;
          }
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     build_derivs_layers(d, d->derivs.layers, 1);
     build_derivs_stacks(d, d->derivs.layers, d->derivs.stacks);

     if (d->options & XRTM_OPTION_SOURCE_SOLAR)
          build_derivs_beam(d, d->derivs.layers, d->derivs.beam);

     if (d->options & XRTM_OPTION_SOURCE_THERMAL)
          build_derivs_thermal(d, d->derivs.layers, d->derivs.thermal);

     build_derivs_sources(d, d->derivs.layers, d->derivs.beam, d->derivs.thermal, d->derivs.sources);

     build_derivs_total(d, d->derivs.layers, d->derivs.sources, d->derivs.total);

     build_derivs_adding_down(d, d->derivs.layers, d->derivs.sources, d->derivs.adding_down);
     build_derivs_adding_up  (d, d->derivs.layers, d->derivs.sources, d->derivs.adding_up);


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     build_derivs_layers(d, d->derivs.layers_m, 0);
     build_derivs_stacks(d, d->derivs.layers_m, d->derivs.stacks_m);

     if (d->options & XRTM_OPTION_SOURCE_SOLAR)
          build_derivs_beam(d, d->derivs.layers_m, d->derivs.beam_m);

     if (d->options & XRTM_OPTION_SOURCE_THERMAL)
          build_derivs_thermal(d, d->derivs.layers_m, d->derivs.thermal_m);

     build_derivs_sources(d, d->derivs.layers_m, d->derivs.beam_m, d->derivs.thermal_m, d->derivs.sources_m);

     build_derivs_total(d, d->derivs.layers_m, d->derivs.sources_m, d->derivs.total_m);

     build_derivs_adding_down(d, d->derivs.layers_m, d->derivs.sources_m, d->derivs.adding_down_m);
     build_derivs_adding_up  (d, d->derivs.layers_m, d->derivs.sources_m, d->derivs.adding_up_m);


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (d->options & XRTM_OPTION_REVERSE_DERIVS) {
          derivs_union_logical_or2(d->n_layers + 1, d->n_derivs, d->derivs.layers,  d->derivs.layers_union);
/*
          derivs_union_logical_or2(d->n_layers + 1, d->n_derivs, d->derivs.total,   d->derivs.total_union);
*/
          derivs_union_logical_or2(d->n_layers + 1, d->n_derivs, d->derivs.stacks,  d->derivs.stacks_union);
          derivs_union_logical_or2(d->n_layers + 1, d->n_derivs, d->derivs.sources, d->derivs.sources_union);
          if (d->options & XRTM_OPTION_SOURCE_SOLAR)
               derivs_union_logical_or2(d->n_layers + 1, d->n_derivs, d->derivs.beam, d->derivs.beam_union);

          if (d->options & XRTM_OPTION_SOURCE_THERMAL)
               derivs_union_logical_or2(d->n_layers + 1, d->n_derivs, d->derivs.thermal, d->derivs.thermal_union);

          derivs_union_bitwise_or2(d->n_layers + 1, d->n_derivs, d->derivs.adding_up,   d->derivs.adding_up_union);
          derivs_union_bitwise_or2(d->n_layers + 1, d->n_derivs, d->derivs.adding_down, d->derivs.adding_down_union);


          /*--------------------------------------------------------------------
           *
           *------------------------------------------------------------------*/
          derivs_union_logical_or2(d->n_layers + 1, d->n_derivs, d->derivs.layers,      d->derivs.layers_m_union);
/*
          derivs_union_logical_or2(d->n_layers + 1, d->n_derivs, d->derivs.total,       d->derivs.total_m_union);
*/
          derivs_union_logical_or2(d->n_layers + 1, d->n_derivs, d->derivs.stacks,      d->derivs.stacks_m_union);
          derivs_union_logical_or2(d->n_layers + 1, d->n_derivs, d->derivs.sources,     d->derivs.sources_m_union);
          if (d->options & XRTM_OPTION_SOURCE_SOLAR)
               derivs_union_logical_or2(d->n_layers + 1, d->n_derivs, d->derivs.beam, d->derivs.beam_union);

          if (d->options & XRTM_OPTION_SOURCE_THERMAL)
               derivs_union_logical_or2(d->n_layers + 1, d->n_derivs, d->derivs.thermal, d->derivs.thermal_m_union);

          derivs_union_bitwise_or2(d->n_layers + 1, d->n_derivs, d->derivs.adding_up_m,   d->derivs.adding_up_m_union);
          derivs_union_bitwise_or2(d->n_layers + 1, d->n_derivs, d->derivs.adding_down_m, d->derivs.adding_down_m_union);
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (d->solvers & (XRTM_SOLVERS_INTERNAL & XRTM_SOLVERS_ADDING) && d->options & XRTM_OPTION_STACK_REUSE_ADDING) {
          size = d->n_layers * d->n_layers;
          if (! (stack_chain = (stack_data *) alloc_array1(size, sizeof(stack_data)))) {
               fprintf(stderr, "ERROR: alloc_array1(%ld, %ld)\n", size, sizeof(stack_data));
               return XRTM_INT_ERROR;
          }

          n_stacks = build_stack_chain(d->n_layers, d->n_derivs, d->derivs.layers, stack_chain);

          if (stack_chain_alloc(d->n_four, d->n_quad_v_x, d->n_derivs, d->n_layers, d->n_stacks, d->stack_chain, n_stacks, stack_chain, d->derivs.layers, d->derivs.stacks, d->derivs.adding_down)) {
               fprintf(stderr, "ERROR: stack_chain_alloc()\n");
               return XRTM_INT_ERROR;
          }

          stack_chain_free(d->n_four, d->n_quad_v_x, d->n_derivs, d->n_layers, d->n_stacks, d->stack_chain, d->derivs.layers, d->derivs.stacks, d->derivs.adding_down);

          free_array1((void * ) d->stack_chain);

          d->n_stacks    = n_stacks;
          d->stack_chain = stack_chain;

          for (i = 0; i < d->n_layers; ++i) {
               for (j = 0; j < d->n_layers; ++j) {
                   d->stack_grid[i][j] = NULL;
               }
          }

          for (i = 0; i < d->n_stacks; ++i)
               d->stack_grid[d->stack_chain[i].i1][d->stack_chain[i].i2] = &(d->stack_chain[i]);
     }

     d->varied_layers_updated = 0;

     return 0;
}



/*******************************************************************************
 *
 ******************************************************************************/
int xrtm_qx(xrtm_data *d, double *qx) {

     int i;
/*
     SOLUTION_PREPROCESS();
*/
     for (i = 0; i < d->n_quad + d->n_umus; ++i)
          qx[i] = d->qx[i];

     return 0;
}



int xrtm_qw(xrtm_data *d, double *qw) {

     int i;
/*
     SOLUTION_PREPROCESS();
*/
     for (i = 0; i < d->n_quad + d->n_umus; ++i)
          qw[i] = d->qw[i];

     return 0;
}



int xrtm_kernel_qx(xrtm_data *d, double *kernel_qx) {

     int i;
/*
     SOLUTION_PREPROCESS();
*/
     for (i = 0; i < d->n_kernel_quad; ++i)
          kernel_qx[i] = d->kernel_qx[i];

     return 0;
}



int xrtm_kernel_qw(xrtm_data *d, double *kernel_qw) {

     int i;
/*
     SOLUTION_PREPROCESS();
*/
     for (i = 0; i < d->n_kernel_quad; ++i)
          kernel_qw[i] = d->kernel_qw[i];

     return 0;
}



/*******************************************************************************
 *
 ******************************************************************************/
#include "xrtm_pade_rts.h"

int xrtm_pade_params(xrtm_data *d, int i_four, int i_layer, int *pade_s, int *pade_r) {

     int flag = 0;

     double **r_p;
     double **t_p;
     double **r_m;
     double **t_m;

     double ***r_p_l;
     double ***t_p_l;
     double ***r_m_l;
     double ***t_m_l;

     SOLUTION_PREPROCESS();

     CHECK_INDEX_RANGE(i_four,  d->n_four,   "i_four",  "n_four",   XRTM_INT_ERROR);
     CHECK_INDEX_RANGE(i_layer, d->n_layers, "i_layer", "n_layers", XRTM_INT_ERROR);

     if (get_local_r_t_u_w(d, i_four, i_layer, &r_p, &t_p, &r_m, &t_m, &r_p_l, &t_p_l, &r_m_l, &t_m_l, 1, d->save_tree, &d->work)) {
          fprintf(stderr, "ERROR: get_local_r_t_u_w()\n");
          return XRTM_INT_ERROR;
     }

     if (d->options & XRTM_OPTION_CALC_DERIVS)
          flag = flags_or(d->derivs.layers, d->n_derivs) | flags_or(d->derivs.sources, d->n_derivs);

     pade_get_s_and_r(r_p, t_p, d->n_quad_v, d->ltau[i_layer], d->n_umus, d->umus, pade_s, pade_r, flag);

     return 0;
}



/*******************************************************************************
 *
 ******************************************************************************/
int xrtm_local_r_t_u_w
     (xrtm_data *d, int i_four, int i_layer,
      double  **r_p_, double  **t_p_, double  **r_m_, double  **t_m_,
      double ***r_p_l_, double ***t_p_l_, double ***r_m_l_, double ***t_m_l_) {

     int i;

     double **r_p;
     double **t_p;
     double **r_m;
     double **t_m;

     double ***r_p_l;
     double ***t_p_l;
     double ***r_m_l;
     double ***t_m_l;

     SOLUTION_PREPROCESS();

     CHECK_INDEX_RANGE(i_four,  d->n_four,   "i_four",  "n_four",   XRTM_INT_ERROR);
     CHECK_INDEX_RANGE(i_layer, d->n_layers, "i_layer", "n_layers", XRTM_INT_ERROR);

     if (get_local_r_t_u_w(d, i_four, i_layer, &r_p, &t_p, &r_m, &t_m, &r_p_l, &t_p_l, &r_m_l, &t_m_l, 1, d->save_tree, &d->work)) {
          fprintf(stderr, "ERROR: get_local_r_t_u_w()\n");
          return XRTM_INT_ERROR;
     }

     dmat_copy(r_p_, r_p, d->n_quad_v_x, d->n_quad_v_x);
     dmat_copy(t_p_, t_p, d->n_quad_v_x, d->n_quad_v_x);

     if (d->options & XRTM_OPTION_VECTOR) {
          dmat_copy(r_m_, r_m, d->n_quad_v_x, d->n_quad_v_x);
          dmat_copy(t_m_, t_m, d->n_quad_v_x, d->n_quad_v_x);
     }

     if (d->options & XRTM_OPTION_CALC_DERIVS && flags_or(d->derivs.layers[i_layer], d->n_derivs)) {
          for (i = 0; i < d->n_derivs; ++i) {
               if (! d->derivs.layers[i_layer][i])
                    continue;

               dmat_copy(r_p_l_[i], r_p_l[i], d->n_quad_v_x, d->n_quad_v_x);
               dmat_copy(t_p_l_[i], t_p_l[i], d->n_quad_v_x, d->n_quad_v_x);

               if (d->options & XRTM_OPTION_VECTOR) {
                    dmat_copy(r_m_l_[i], r_m_l[i], d->n_quad_v_x, d->n_quad_v_x);
                    dmat_copy(t_m_l_[i], t_m_l[i], d->n_quad_v_x, d->n_quad_v_x);
               }
          }
     }

     return 0;
}



/*******************************************************************************
 *
 ******************************************************************************/
int xrtm_solution(xrtm_data *d, enum xrtm_solver_mask solver, int solutions, int n_out_phis, double **out_phis, double ****I_p, double ****I_m, double *****I_p_l, double *****I_m_l, double *mean_p, double *mean_m, double **mean_p_l, double **mean_m_l, double *flux_p, double *flux_m, double **flux_p_l, double **flux_m_l, double *flux_div, double **flux_div_l) {

     if (! (solver & XRTM_SOLVERS_ALL)) {
          fprintf(stderr, "ERROR: invalid value for solver: %d\n", solver);
          return XRTM_INT_ERROR;
     }

     if (! (solver & d->solvers)) {
          fprintf(stderr, "ERROR: model not initialized for solver %s\n",
                  xrtm_solver_mask_to_name(solver));
          return XRTM_INT_ERROR;
     }

     if (! (solutions & XRTM_SOLUTION_ALL)) {
          fprintf(stderr, "ERROR: invalid value for solution\n");
          return XRTM_INT_ERROR;
     }

     if (solutions & XRTM_OUTPUT_RADIANCE && n_out_phis <= 0) {
          fprintf(stderr, "ERROR: invalid value for n_out_phis: %d, must be > zero for radiance output\n", n_out_phis);
          return XRTM_INT_ERROR;
     }

     if (solutions & XRTM_OUTPUT_RADIANCE_MEAN) {
          if (check_solvers(solver, 0, 1, "output option", xrtm_output_mask_to_name(XRTM_OUTPUT_RADIANCE_MEAN), XRTM_SOLVER_SINGLE, 0))
               return XRTM_INT_ERROR;
     }

     if (solutions & XRTM_OUTPUT_FLUX) {
          if (check_solvers(solver, 0, 1, "output option", xrtm_output_mask_to_name(XRTM_OUTPUT_FLUX), XRTM_SOLVER_SINGLE, 0))
               return XRTM_INT_ERROR;
     }

     if (solutions & XRTM_OUTPUT_FLUX_DIVERGENCE) {
          if (check_solvers(solver, 0, 1, "output option", xrtm_output_mask_to_name(XRTM_OUTPUT_FLUX_DIVERGENCE), XRTM_SOLVER_SINGLE, 0))
               return XRTM_INT_ERROR;
     }

#ifdef INCLUDE_DEV_SOURCE
     if (check_dev_solvers_solution(d, solver, solutions, n_out_phis, out_phis)) {
          fprintf(stderr, "ERROR: check_dev_solvers_solution()\n");
          return XRTM_INT_ERROR;
     }
#endif
     SOLUTION_PREPROCESS();

     if (get_solution(d, solver, solutions, n_out_phis, out_phis, I_p, I_m, I_p_l, I_m_l, mean_p, mean_m, mean_p_l, mean_m_l, flux_p, flux_m, flux_p_l, flux_m_l, flux_div, flux_div_l, d->save_tree, d->work)) {
          fprintf(stderr, "ERROR: get_solution()\n");
          return XRTM_INT_ERROR;
     }

     return 0;
}



/*******************************************************************************
 *
 ******************************************************************************/
int xrtm_radiance(xrtm_data *d, enum xrtm_solver_mask solver, int n_out_phis, double **out_phis, double ****I_p, double ****I_m, double *****I_p_l, double *****I_m_l) {

     if (xrtm_solution(d, solver, XRTM_OUTPUT_RADIANCE, n_out_phis, out_phis, I_p, I_m, I_p_l, I_m_l, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL)) {
          fprintf(stderr, "ERROR: xrtm_solution()\n");
          return XRTM_INT_ERROR;
     }

     return 0;
}



/*******************************************************************************
 *
 ******************************************************************************/
int xrtm_mean_radiance(xrtm_data *d, enum xrtm_solver_mask solver, double *mean_p, double *mean_m, double **mean_p_l, double **mean_m_l) {

     if (xrtm_solution(d, solver, XRTM_OUTPUT_RADIANCE_MEAN, 0, NULL, NULL, NULL, NULL, NULL, mean_p, mean_m, mean_p_l, mean_m_l, NULL, NULL, NULL, NULL, NULL, NULL)) {
          fprintf(stderr, "ERROR: xrtm_solution()\n");
          return XRTM_INT_ERROR;
     }

     return 0;
}



/*******************************************************************************
 *
 ******************************************************************************/
int xrtm_flux(xrtm_data *d, enum xrtm_solver_mask solver, double *flux_p, double *flux_m, double **flux_p_l, double **flux_m_l) {

     if (xrtm_solution(d, solver, XRTM_OUTPUT_FLUX, 0, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, flux_p, flux_m, flux_p_l, flux_m_l, NULL, NULL)) {
          fprintf(stderr, "ERROR: xrtm_solution()\n");
          return XRTM_INT_ERROR;
     }

     return 0;
}



/*******************************************************************************
 *
 ******************************************************************************/
int xrtm_flux_divergence(xrtm_data *d, enum xrtm_solver_mask solver, double *flux_div, double **flux_div_l) {

     fprintf(stderr, "ERROR: flux_divergence not yet supported as an output\n");
     return XRTM_INT_ERROR;


     if (xrtm_solution(d, solver, XRTM_OUTPUT_FLUX_DIVERGENCE, 0, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, flux_div, flux_div_l)) {
          fprintf(stderr, "ERROR: xrtm_solution()\n");
          return XRTM_INT_ERROR;
     }

     return 0;
}



/*******************************************************************************
 *
 ******************************************************************************/
int xrtm_solution_2(xrtm_data *d, enum xrtm_solver_mask solver, int solutions, int n_out_phis, double **out_phis, double ****I_p, double ****I_m, double *****I_p_l, double *****I_m_l, double *mean_p, double *mean_m, double **mean_p_l, double **mean_m_l, double *flux_p, double *flux_m, double **flux_p_l, double **flux_m_l, double *flux_div, double **flux_div_l) {

     int i;
     int ii;
     int j;
     int k;
     int l;
     int m;

     int n_mus2;

     int options;

     int n_derivs;

     int n_dirs;

     int i_input;
     int n_inputs;

     int ii_input;
     int jj_input;
     int kk_input;
     int ll_input;

     int input[2];

     double ****I_p_a;
     double ****I_m_a;

     double *****a;

     work_data work2;


     work2 = d->work;


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (solver & XRTM_SOLVERS_EXT_EXTERNAL) {
          if (get_solution(d, solver, solutions, n_out_phis, out_phis, I_p, I_m, I_p_l, I_m_l, mean_p, mean_m, mean_p_l, mean_m_l, flux_p, flux_m, flux_p_l, flux_m_l, flux_div, flux_div_l, d->save_tree, work2)) {
               fprintf(stderr, "ERROR: get_solution()\n");
               return -1;
          }

          return 0;
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (! (d->options & XRTM_OPTION_CALC_DERIVS)) {
          fprintf(stderr, "ERROR: xrtm_solution() requires option \"%s\"\n", xrtm_option_mask_to_name(XRTM_OPTION_CALC_DERIVS));
          return -1;
     }

     if (! (d->options & XRTM_OPTION_REVERSE_DERIVS)) {
          fprintf(stderr, "ERROR: xrtm_solution() requires option \"%s\"\n", xrtm_option_mask_to_name(XRTM_OPTION_REVERSE_DERIVS));
          return -1;
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (d->n_umus == 0)
          n_mus2   = d->n_quad;
     else
          n_mus2   = d->n_umus;
/*
     else
          n_mus2   = d->n_quad + d->n_umus;
*/

     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     options  = d->options;
     if (d->options & XRTM_OPTION_CALC_DERIVS)
          d->options ^= XRTM_OPTION_CALC_DERIVS;

     n_derivs = d->n_derivs;
     d->n_derivs = 0;

     if (get_solution(d, solver, solutions, n_out_phis, out_phis, I_p, I_m, NULL, NULL, mean_p, mean_m, NULL, NULL, flux_p, flux_m, NULL, NULL, flux_div, NULL, d->save_tree, work2)) {
          fprintf(stderr, "ERROR: get_solution2()\n");
          return -1;
     }

     d->options  = options;
     d->n_derivs = n_derivs;


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     init_array5_d(I_p_l, d->n_ulevels, d->n_derivs, n_mus2, n_out_phis, d->n_stokes, 0.);
     init_array5_d(I_m_l, d->n_ulevels, d->n_derivs, n_mus2, n_out_phis, d->n_stokes, 0.);


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     I_p_a = get_work_d4(&work2, d->n_ulevels, n_mus2, n_out_phis, d->n_stokes);
     I_m_a = get_work_d4(&work2, d->n_ulevels, n_mus2, n_out_phis, d->n_stokes);


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     init_array4_d(I_p_a, d->n_ulevels, n_mus2, n_out_phis, d->n_stokes, 0.);
     init_array4_d(I_m_a, d->n_ulevels, n_mus2, n_out_phis, d->n_stokes, 0.);


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (d->options & XRTM_OPTION_UPWELLING_OUTPUT &&
         d->options & XRTM_OPTION_DOWNWELLING_OUTPUT) {
          n_dirs   = 2;
          input[0] = 0;
          input[1] = 1;
     }
     else
     if (d->options & XRTM_OPTION_UPWELLING_OUTPUT) {
          n_dirs   = 1;
          input[0] = 0;
     }
     else
     if (d->options & XRTM_OPTION_DOWNWELLING_OUTPUT) {
          n_dirs   = 1;
          input[0] = 1;
     }

     n_inputs = n_dirs * d->n_ulevels * n_mus2 * n_out_phis * d->n_stokes;

     for (i = 0; i < n_inputs; ++i) {
          if (i < d->n_ulevels * n_mus2 * n_out_phis * d->n_stokes) {
               i_input = input[0];
               ii = i;
               ii_input = ii / (n_mus2 * n_out_phis * d->n_stokes);
               jj_input = ii % (n_mus2 * n_out_phis * d->n_stokes) / (n_out_phis * d->n_stokes);
               kk_input = ii % (n_mus2 * n_out_phis * d->n_stokes) % (n_out_phis * d->n_stokes) / d->n_stokes;
               ll_input = ii % (n_mus2 * n_out_phis * d->n_stokes) % (n_out_phis * d->n_stokes) % d->n_stokes;
          }
          else {
               i_input = input[1];
               ii = i - d->n_ulevels * n_mus2 * n_out_phis * d->n_stokes;
               ii_input = ii / (n_mus2 * n_out_phis * d->n_stokes);
               jj_input = ii % (n_mus2 * n_out_phis * d->n_stokes) / (n_out_phis * d->n_stokes);
               kk_input = ii % (n_mus2 * n_out_phis * d->n_stokes) % (n_out_phis * d->n_stokes) / d->n_stokes;
               ll_input = ii % (n_mus2 * n_out_phis * d->n_stokes) % (n_out_phis * d->n_stokes) % d->n_stokes;
          }


          /*--------------------------------------------------------------------
           *
           *------------------------------------------------------------------*/
          init_opt_props0_all_a(d, d->n_coef, work2);

          if (d->options & XRTM_OPTION_DELTA_M)
               init_opt_props_all_a(d, d->n_coef2, work2);

          init_beam_params_all_a(d, work2);

          init_diff_bound_input_all_a(d, d->n_four, work2);

          for (j = 0; j < d->n_kernels; ++j)
               d->kernel_ampfac_a[j] = 0.;


          if (i_input == 0)
               I_p_a[ii_input][jj_input][kk_input][ll_input] = 1.;
          else
               I_m_a[ii_input][jj_input][kk_input][ll_input] = 1.;

          get_solution_internal_a(d, solver, solutions, n_out_phis, out_phis, I_p, I_m, I_p_a, I_m_a, mean_p, mean_m, NULL, NULL, flux_p, flux_m, NULL, NULL, NULL, NULL, d->save_tree, work2);


          /*--------------------------------------------------------------------
           *
           *------------------------------------------------------------------*/
          if (i_input == 0)
               a = I_p_l;
          else
               a = I_m_l;

          for (j = 0; j < d->n_derivs; ++j)
               a[ii_input][j][jj_input][kk_input][ll_input] = 0.;

          for (j = 0; j < d->n_kernels; ++j) {
               for (k = 0; k < d->n_derivs; ++k) {
                    a[ii_input][k][jj_input][kk_input][ll_input] += d->kernel_ampfac_l[j][k] * d->kernel_ampfac_a[j];
               }
          }

          for (j = 0; j < d->n_layers; ++j) {
               for (k = 0; k < d->n_derivs; ++k) {
                    a[ii_input][k][jj_input][kk_input][ll_input] += d->ltau0_l [j][k] * d->ltau0_a [j];

                    if (d->derivs.layers[j][k]) {
                         a[ii_input][k][jj_input][kk_input][ll_input] += d->omega0_l[j][k] * d->omega0_a[j];

                         for (l = 0; l < d->n_elem; ++l) {
                              for (m = 0; m < d->n_coef_layer[j]; ++m) {
                                   a[ii_input][k][jj_input][kk_input][ll_input] += d->coef0_l[j][k][l][m] * d->coef0_a[j][l][m];
                              }
                         }
                    }
               }
          }
     }


     return 0;
}



/*******************************************************************************
 *
 ******************************************************************************/
int xrtm_radiance_2(xrtm_data *d, enum xrtm_solver_mask solver, int n_out_phis, double **out_phis, double ****I_p, double ****I_m, double *****I_p_l, double *****I_m_l) {

     if (xrtm_solution_2(d, solver, XRTM_OUTPUT_RADIANCE, n_out_phis, out_phis, I_p, I_m, I_p_l, I_m_l, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL)) {
          fprintf(stderr, "ERROR: xrtm_solution_2()\n");
          return XRTM_INT_ERROR;
     }

     return 0;
}



/*******************************************************************************
 *
 ******************************************************************************/
#ifdef INCLUDE_DEV_SOURCE
#include "../src2/xrtm_interface2.c"
#endif
