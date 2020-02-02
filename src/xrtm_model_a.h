/*******************************************************************************
**
**    Copyright (C) 2007-2020 Greg McGarragh <greg.mcgarragh@colostate.edu>
**
**    This source code is licensed under the GNU General Public License (GPL),
**    Version 3.  See the file COPYING for more details.
**
*******************************************************************************/

#ifndef XRTM_MODEL_A_H
#define XRTM_MODEL_A_H

#include "xrtm_interface.h"
#include "xrtm_save_tree.h"
#include "xrtm_work.h"

#ifdef __cplusplus
extern "C" {
#endif


typedef struct {
     void (*free)(void *);

     int options;

     double **P_qq_pp;
     double **P_qq_mp;
     double **P_qq_mm;
     double **P_qq_pm;
} forward_save_get_local_r_t_u_w_data;


typedef struct {
     void (*free)(void *);

     int options;

     double *P_q0_mm;
     double *P_q0_pm;

     double **r_p;
     double **t_p;
     double **r_m;
     double **t_m;
} forward_save_get_layer_R_T_S_U_W_V_data;


typedef struct {
     void (*free)(void *);

     double ***R1_m;
     double  **S1_p;

     double ***R2_p;
     double ***T2_p;
     double ***R2_m;
     double ***T2_m;
     double  **S2_p;
     double  **S2_m;

     double  **S3_p;
     double  **S3_m;
} forward_save_get_total_TOA_R_S_U_V_data;


typedef struct {
     void (*free)(void *);

     double ***R1_p;
     double ***T1_p;
     double ***R1_m;
     double ***T1_m;
     double  **S1_p;
     double  **S1_m;

     double ***R2_p;
     double ***T2_p;
     double ***R2_m;
     double ***T2_m;
     double  **S2_p;
     double  **S2_m;

     double  **S3_p;
     double  **S3_m;
} forward_save_get_total_R_T_S_U_W_V_data;


typedef struct {
     void (*free)(void *);

     double *Rs_q0;
     double *Rs_u0;
} forward_save_update_diff_bound_input_data;


typedef struct {
     void (*free)(void *);

     int options;

     void ***polys_up;
     void ***polys_dn;

     double ****P_full_up;
     double ****P_full_dn;
} forward_save_get_single_data;


typedef struct {
     void (*free)(void *);

     double **R_m;
     double  *S_p;

     double **Rs_qq;
} forward_save_fourier_get_add_bottom_up_data;


typedef struct {
     void (*free)(void *);

     double **R_p;
     double **T_p;
     double **R_m;
     double **T_m;
     double  *S_p;
     double  *S_m;

     double  *I_p1;
     double  *I_m1;

     double **Rs_qq;
} forward_save_fourier_get_add_both_ways_data;


typedef struct {
     void (*free)(void *);

     int options;
     int solver;
     int surface;
     int n_umus_v;

     double **P_q0_mm;
     double **P_q0_pm;

     double **P_u0_mm;
     double **P_u0_pm;

     double ***P_qq_pp;
     double ***P_qq_mp;
     double ***P_qq_mm;
     double ***P_qq_pm;

     double ***P_uq_pp;
     double ***P_uq_mp;
     double ***P_uq_mm;
     double ***P_uq_pm;

     double ***r_p;
     double ***t_p;
     double ***r_m;
     double ***t_m;

     double *Rs_q0;
     double *Rs_u0;
     double **Rs_qq;
     double **Rs_uq;
} forward_save_fourier_get_bvp_data;


typedef struct {
     void (*free)(void *);

     int options;

     void ***polys_up;
     void ***polys_dn;

     double ****P_trun_up;
     double ****P_full_up;
     double ****P_trun_dn;
     double ****P_full_dn;
} forward_save_apply_corrections_radiance_data;


typedef struct {
     void (*free)(void *);

     double  *In_p;
     double  *I1_m;

     double **P_q0_mm;
     double **P_q0_pm;

     double **P_u0_mm;
     double **P_u0_pm;
} forward_save_apply_corrections_fourier_data;


typedef struct {
     void (*free)(void *);

     int count;

     double ***i_p;
     double ***i_m;
} forward_save_get_solution_internal_data;


#include "prototypes/xrtm_model_a_p.h"


#ifdef __cplusplus
}
#endif

#endif /* XRTMMODEL_A_H */
