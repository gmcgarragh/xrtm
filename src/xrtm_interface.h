/******************************************************************************%
**
**    Copyright (C) 2007-2012 Greg McGarragh <gregm@atmos.colostate.edu>
**
**    This source code is licensed under the GNU General Public License (GPL),
**    Version 3.  See the file COPYING for more details.
**
*******************************************************************************/

#ifndef XRTM_MODEL_H
#define XRTM_MODEL_H

#include "xrtm.h"
#include "xrtm_brdf.h"
#include "xrtm_save_tree.h"
#include "xrtm_stacks.h"
#include "xrtm_work.h"

#ifdef __cplusplus
extern "C" {
#endif


#ifdef INCLUDE_DEV_SOURCE
#include "../src2/xrtm_interface2.h"
#endif


#define XRTM_INT_ERROR	-INT_MAX
#define XRTM_DBL_ERROR	-DBL_MAX


#define N_XRTM_OPTIONS 30
/*
#define N_XRTM_OPTIONS 31
*/

enum xrtm_option_mask {
     XRTM_OPTION_CALC_DERIVS        = (1<<0),
/*
     XRTM_OPTION_FORWARD_DERIVS     = (1<<0),
*/
     XRTM_OPTION_REVERSE_DERIVS     = (1<<1),
     XRTM_OPTION_DELTA_M            = (1<<2),
     XRTM_OPTION_N_T_TMS            = (1<<3),
     XRTM_OPTION_FOUR_CONV_OLD      = (1<<4),
     XRTM_OPTION_FOUR_CONV_NEW      = (1<<5),
     XRTM_OPTION_NO_AZIMUTHAL       = (1<<6),
     XRTM_OPTION_OUTPUT_AT_LEVELS   = (1<<7),
     XRTM_OPTION_OUTPUT_AT_TAUS     = (1<<8),
     XRTM_OPTION_PHASE_SCALAR       = (1<<9),
     XRTM_OPTION_PHASE_MATRIX_GC    = (1<<10),
     XRTM_OPTION_PHASE_MATRIX_LC    = (1<<11),
     XRTM_OPTION_PSA                = (1<<12),
     XRTM_OPTION_QUAD_NORM_GAUS_LEG = (1<<13),
     XRTM_OPTION_QUAD_DOUB_GAUS_LEG = (1<<14),
     XRTM_OPTION_QUAD_LOBATTO       = (1<<15),
     XRTM_OPTION_SAVE_LEG_POLYS     = (1<<16),
     XRTM_OPTION_SAVE_PHASE_MATS    = (1<<17),
     XRTM_OPTION_SAVE_LOCAL_R_T     = (1<<18),
     XRTM_OPTION_SAVE_LAYER_R_T_S   = (1<<19),
     XRTM_OPTION_SAVE_TOTAL_R_T_S   = (1<<20),
     XRTM_OPTION_SFI                = (1<<21),
     XRTM_OPTION_SOURCE_SOLAR       = (1<<22),
     XRTM_OPTION_SOURCE_THERMAL     = (1<<23),
     XRTM_OPTION_STACK_REUSE_ADDING = (1<<24),
     XRTM_OPTION_TOP_DOWN_ADDING    = (1<<25),
     XRTM_OPTION_BOTTOM_UP_ADDING   = (1<<26),
     XRTM_OPTION_UPWELLING_OUTPUT   = (1<<27),
     XRTM_OPTION_DOWNWELLING_OUTPUT = (1<<28),
     XRTM_OPTION_VECTOR             = (1<<29)
};


#ifndef INCLUDE_DEV_SOURCE
#define N_XRTM_SOLVERS (8)
#else
#define N_XRTM_SOLVERS (8 + N_XRTM_DEV_SOLVERS)
#endif

enum xrtm_solver_mask {
     XRTM_SOLVER_DOUB_ADD    = (1<<0),
     XRTM_SOLVER_EIG_ADD     = (1<<1),
     XRTM_SOLVER_EIG_BVP     = (1<<2),
     XRTM_SOLVER_MEM_BVP     = (1<<3),
     XRTM_SOLVER_PADE_ADD    = (1<<4),
     XRTM_SOLVER_SINGLE      = (1<<5),
     XRTM_SOLVER_SOS         = (1<<6),
     XRTM_SOLVER_TWO_OS      = (1<<7)
#ifdef INCLUDE_DEV_SOURCE
,    XRTM_DEV_SOLVER_ENUMS
#endif
};


#define N_XRTM_SOLUTIONS 4

enum xrtm_output_mask {
     XRTM_OUTPUT_RADIANCE        = (1<<0),
     XRTM_OUTPUT_RADIANCE_MEAN   = (1<<1),
     XRTM_OUTPUT_FLUX            = (1<<2),
     XRTM_OUTPUT_FLUX_DIVERGENCE = (1<<3)
};


typedef struct {
     int use_pade_check_condition;

     int use_rebuild_stacks;

     int use_symmetric_form;

     int eigen_solver_gen_real;
     int eigen_solver_gen_complex;

     int threshold_sym_mul_block;

     double threshold_mu_0_singlarity;
     double threshold_omega_singlarity;
#ifdef INCLUDE_DEV_SOURCE
     XRTM_DEV_MISC_INPUT
#endif
} misc_input_data;


typedef struct {
     int fourier_count;

     double pade_condition;
} misc_output_data;


typedef struct {
     uchar set_flags_doub_d_tau;
     uchar set_flags_pade_params;
     uchar set_flags_sos_params;
     uchar set_flags_fourier_tol;
     uchar set_flags_lambda;
     uchar set_flags_F_0;
     uchar set_flags_mu_0;
     uchar set_flags_phi_0;
     uchar set_flags_ulevels;
     uchar set_flags_utaus;
     uchar set_flags_umus;
     uchar set_flags_top_b;
     uchar set_flags_planet_r;
     uchar set_flags_levels_z;
     uchar set_flags_levels_b;
     uchar set_flags_surface_b;
     uchar *set_flags_g;
     uchar *set_flags_coef;
     uchar *set_flags_omega;
     uchar *set_flags_ltau;
     uchar *set_flags_kernel_ampfac;
     uchar **set_flags_kernel_params;

     uchar **derivs_h;
     uchar **derivs_s;
     uchar **derivs_u;
     uchar **derivs_d;

     uchar **derivs_hm;
     uchar **derivs_sm;
     uchar **derivs_um;
     uchar **derivs_dm;

     uchar *derivs_h_union;
     uchar *derivs_u_union;
     uchar *derivs_d_union;

     uchar *derivs_dm_union;
     uchar *derivs_um_union;

     int initial_inputs;
     int initial_derivs;

     int solvers;
     int options;

     int n_four;
     int n_four2;

     int n_coef;
     int n_coef2;

     int n_elem;

     int n_quad;
     int n_quad_x;
     int n_quad_d;
     int n_quad_v;
     int n_quad_v_x;
     int n_quad_v_d;

     int n_stokes;

     int n_derivs;

     int n_layers;

     int n_kernels;
     int n_kernel_quad;

     int n_ulevels;

     int n_umus;
     int n_umus_v;

     int order_p;
     int order_0;
     int order_u;

     int quad_type;

     int pade_s;
     int pade_r;

     int sos_max_os;

     int n_stacks;

     enum xrtm_kernel_type *kernels;

     int *ulevels;

     int *n_coef_layer;
     int *n_coef_layer2;

     int dep_flags_utaus;
     int dep_flags_chapman;

     int *dep_flags_Y_p;
     int *dep_flags_Y_0;
     int *dep_flags_Y_u;
     int *dep_flags_opt_props;
     int *dep_flags_beam_params0;
     int *dep_flags_beam_params;
     int *dep_flags_diff_bound_input0;
     int *dep_flags_diff_bound_input;
     int *dep_flags_total_R_T_S_U_W_V;

     int **dep_flags_phase_mats_qq;
     int **dep_flags_phase_mats_uq;
     int **dep_flags_phase_vecs_q0;
     int **dep_flags_phase_vecs_u0;
     int **dep_flags_local_r_t_u_w;
     int **dep_flags_layer_R_T_S_U_W_V;

     double qf;

     double c_p;
     double d_p;

     double c_u;
     double d_u;

     double c_0;
     double d_0;

     double fourier_tol;

     double lambda;

     double F_0;
     double mu_0;
     double phi_0;

     double mu_0_0;

     double top_b;
     double *top_b_l;

     double planet_r;

     double doub_d_tau;

     double sos_max_tau;
     double sos_tol;

     double *qx;
     double *qw;

     double *qx_v;
     double *qw_v;

     double *kernel_qx;
     double *kernel_qw;

     double **Y_0;
     double ***Y_p;
     double ***Y_u;

     double *alpha1;
     double *alpha2;

     double *utaus;
     double *utaus20;
     double *utaus2;

     double *umus;
     double *umus_v;

     double *levels_z;
     double *levels_b;
     double **levels_b_l;

     double *g;
     double **g_l;
     double *g_a;

     double ***coef;
     double ****coef_l;
     double ***coef_a;

     double *omega;
     double **omega_l;
     double *omega_a;

     double *ltau;
     double **ltau_l;
     double *ltau_a;

     double *g0;
     double **g0_l;
     double *g0_a;

     double ***coef0;
     double ****coef0_l;
     double ***coef0_a;

     double *omega0;
     double **omega0_l;
     double *omega0_a;

     double *ltau0;
     double **ltau0_l;
     double *ltau0_a;

     double *ttau;

     double surface_b;
     double *surface_b_l;

     double *kernel_ampfac;
     double **kernel_ampfac_l;
     double *kernel_ampfac_a;

     double **kernel_params;
     double ***kernel_params_l;
     double **kernel_params_a;

     double **chapman;

     double *btau;
     double **btau_l;
     double *btau_a;

     double *btran;
     double **btran_l;
     double *btran_a;

     double *as_0;
     double **as_0_l;
     double *as_0_a;

     double *atran;
     double **atran_l;
     double *atran_a;
/*
     double *btau0;
     double **btau0_l;
     double *btau0_a;

     double *btran0;
     double **btran0_l;
     double *btran0_a;

     double *as_00;
     double **as_00_l;
     double *as_00_a;

     double *atran0;
     double **atran0_l;
     double *atran0_a;
*/
     double **I1_m;
     double **In_p;
     double ***I1_m_l;
     double ***In_p_l;
     double **I1_m_a;
     double **In_p_a;

     double ***P_q0_mm;
     double ***P_q0_pm;

     double ****P_q0_mm_l;
     double ****P_q0_pm_l;

     double ****P_qq_pp;
     double ****P_qq_mp;
     double ****P_qq_mm;
     double ****P_qq_pm;

     double *****P_qq_pp_l;
     double *****P_qq_mp_l;
     double *****P_qq_mm_l;
     double *****P_qq_pm_l;

     double ***P_u0_mm;
     double ***P_u0_pm;

     double ****P_u0_mm_l;
     double ****P_u0_pm_l;

     double ****P_uq_pp;
     double ****P_uq_mp;
     double ****P_uq_mm;
     double ****P_uq_pm;

     double *****P_uq_pp_l;
     double *****P_uq_mp_l;
     double *****P_uq_mm_l;
     double *****P_uq_pm_l;

     double ****r_p;
     double ****t_p;
     double ****r_m;
     double ****t_m;

     double *****r_p_l;
     double *****t_p_l;
     double *****r_m_l;
     double *****t_m_l;

     double ***Sl_p;
     double ***Sl_m;

     double ****Rl_p;
     double ****Tl_p;
     double ****Rl_m;
     double ****Tl_m;

     double ****Sl_p_l;
     double ****Sl_m_l;

     double *****Rl_p_l;
     double *****Tl_p_l;
     double *****Rl_m_l;
     double *****Tl_m_l;

     double **Sa_p;
     double **Sa_m;

     double ***Ra_p;
     double ***Ta_p;
     double ***Ra_m;
     double ***Ta_m;

     double ***Sa_p_l;
     double ***Sa_m_l;

     double ****Ra_p_l;
     double ****Ta_p_l;
     double ****Ra_m_l;
     double ****Ta_m_l;

     double **kernel_f_q0;
     double ***kernel_f_q0_l;
     double **kernel_f_q0_a;

     double ***kernel_f_qq;
     double ****kernel_f_qq_l;
     double ***kernel_f_qq_a;

     double **kernel_f_u0;
     double ***kernel_f_u0_l;
     double **kernel_f_u0_a;

     double ***kernel_f_uq;
     double ****kernel_f_uq_l;
     double ***kernel_f_uq_a;

     brdf_aux_data kernel_aux;

     misc_input_data misc_input;
     misc_output_data misc_output;

     save_tree_data save_tree;

     stack_data *stack_chain;
     stack_data ***stack_grid;

     work_data work;
} xrtm_data;


#include "prototypes/xrtm_interface_p.h"


#ifdef __cplusplus
}
#endif

#endif /* XRTM_MODEL_H */
