/******************************************************************************%
**
**    Copyright (C) 2007-2012 Greg McGarragh <gregm@atmos.colostate.edu>
**
**    This source code is licensed under the GNU General Public License (GPL),
**    Version 3.  See the file COPYING for more details.
**
*******************************************************************************/

#include <gutil.h>
#include <gmath_matrix.h>

#include "xrtm.h"
#include "xrtm_save_tree.h"
#include "xrtm_source_a.h"
#include "xrtm_utility.h"


/*******************************************************************************
 *
 ******************************************************************************/
static void forward_save_build_source_vectors_1n_free(forward_save_build_source_vectors_1n_data *d);

int forward_save_build_source_vectors_1n_alloc(forward_save_build_source_vectors_1n_data *d, int n_quad_v) {

     d->free = (void (*)(void *)) forward_save_build_source_vectors_1n_free;

     d->ip = alloc_array1_i(n_quad_v);

     d->d  = alloc_array1_d(n_quad_v);
     d->e  = alloc_array1_d(n_quad_v);
     d->p  = alloc_array1_d(n_quad_v);
     d->h  = alloc_array1_d(n_quad_v);
     d->f  = alloc_array2_d(n_quad_v, n_quad_v);

     return 0;
}



static void forward_save_build_source_vectors_1n_free(forward_save_build_source_vectors_1n_data *d) {

     free_array1_i(d->ip);

     free_array1_d(d->d);
     free_array1_d(d->e);
     free_array1_d(d->p);
     free_array1_d(d->h);
     free_array2_d(d->f);
}



/*******************************************************************************
 *
 ******************************************************************************/
void build_source_vectors_1n_a(int n_quad, int n_stokes, double *qx_v,
                               double F_0, double omega, double *omega_a,
                               double as_0, double *as_0_a,
                               double  *P_q0_mm, double  *P_q0_pm,
                               double  **tpr, double  **tmr, double  **gamma,
                               double  *F_p, double  *F_m,
                               double *P_q0_mm_a, double *P_q0_pm_a,
                               double **tpr_a, double **tmr_a, double **gamma_a,
                               double *F_p_a, double *F_m_a,
                               save_tree_data save_tree, work_data work) {

     int i;

     int n_quad_v;

     double *v1;

     double *a_a;
     double *b_a;
     double *c_a;
     double *d_a;
     double *e_a;
     double *p_a;
     double *h_a;

     double **w1;

     forward_save_build_source_vectors_1n_data *save;


     n_quad_v = n_quad * n_stokes;


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     save_tree_encode_s(&save_tree, "build_source_vectors_1n");

     save_tree_retrieve_data(&save_tree, forward_save_build_source_vectors_1n_data, &save);


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     v1  = get_work1(&work, WORK_DX);

     a_a = get_work1(&work, WORK_DX);
     b_a = get_work1(&work, WORK_DX);
     c_a = get_work1(&work, WORK_DX);
     d_a = get_work1(&work, WORK_DX);
     e_a = get_work1(&work, WORK_DX);
     p_a = get_work1(&work, WORK_DX);
     h_a = get_work1(&work, WORK_DX);

     w1  = get_work1(&work, WORK_DXX);


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     dm_v_mul_D_A(n_quad, n_stokes, F_m_a, F_m_a);

     dvec_add(F_p_a, F_m_a, F_p_a, n_quad_v);

     dvec_scale(-1., F_m_a, p_a, n_quad_v);

     *as_0_a += -dvec_dot(save->h, F_p_a, n_quad_v) / (2. * as_0 * as_0);

     dvec_scale(1. / (2. * as_0), F_p_a, h_a, n_quad_v);

     dvec_scale(1. / 2., F_p_a, v1, n_quad_v);
     dvec_add(p_a, v1, p_a, n_quad_v);

     dmat_vxvtmx(h_a, save->p, 1., tpr_a, 1., n_quad_v, n_quad_v);

     dmat_gxvxmx(1, tpr, h_a, 1., p_a, 1., n_quad_v, n_quad_v);

     dvec_copy(e_a, h_a, n_quad_v);

     dvec_copy(v1, p_a,  n_quad_v);
     dmat_getrs2('t', save->f, &v1, n_quad_v, 1, save->ip);

     dmat_vxvtmx(v1, save->p, -1., w1, 0., n_quad_v, n_quad_v);

     dmat_add(gamma_a, w1, gamma_a, n_quad_v, n_quad_v);

     for (i = 0; i < n_quad_v; ++i)
          *as_0_a += -2. * as_0 * w1[i][i];

     dmat_vxvtmx(v1, save->e, -1., tmr_a, 1., n_quad_v, n_quad_v);

     dmat_gxvxmx(1, tmr, v1, -1., e_a, 1., n_quad_v, n_quad_v);

     *as_0_a -= dvec_dot(v1, save->d, n_quad_v);

     dvec_scale(-as_0, v1, d_a, n_quad_v);

     dvec_copy(b_a, e_a, n_quad_v);
     dvec_scale(-1., e_a, c_a, n_quad_v);

     dvec_add(b_a, d_a, b_a, n_quad_v);
     dvec_add(c_a, d_a, c_a, n_quad_v);

     dm_v_mul_D_A(n_quad, n_stokes, c_a, c_a);

     for (i = 0; i < n_quad_v; ++i)
          a_a[i]  = c_a[i] * P_q0_mm[i];

     for (i = 0; i < n_quad_v; ++i)
          P_q0_mm_a[i] += F_0 * omega / (4. * PI) / qx_v[i] * c_a[i];
     
     for (i = 0; i < n_quad_v; ++i)
          a_a[i] += b_a[i] * P_q0_pm[i];

     for (i = 0; i < n_quad_v; ++i)
          P_q0_pm_a[i] += F_0 * omega / (4. * PI) / qx_v[i] * b_a[i];

     for (i = 0; i < n_quad_v; ++i)
          *omega_a += F_0 / (4. * PI) / qx_v[i] * a_a[i];

     dvec_zero(F_p_a, n_quad_v);
     dvec_zero(F_m_a, n_quad_v);
}



/*******************************************************************************
 *
 ******************************************************************************/
void build_global_source_a(int n_quad_v, double atran, double *atran_a,
                      double **R_p, double **T_p,
                      double *F_p, double *F_m, double *S_p, double *S_m,
                      double **R_p_a, double **T_p_a,
                      double *F_p_a, double *F_m_a, double *S_p_a, double *S_m_a,
                      work_data work) {

     build_global_source_a2(n_quad_v, atran, atran_a,
                            R_p, T_p, R_p, T_p,
                            F_p, F_m, S_p, S_m,
                            R_p_a, T_p_a, R_p_a, T_p_a,
                            F_p_a, F_m_a, S_p_a, S_m_a,
                            work);

}



void build_global_source_a2(int n_quad_v, double atran, double *atran_a,
                       double **R_p, double **T_p, double **R_m, double **T_m,
                       double *F_p, double *F_m, double *S_p, double *S_m,
                       double **R_p_a, double **T_p_a, double **R_m_a, double **T_m_a,
                       double *F_p_a, double *F_m_a, double *S_p_a, double *S_m_a,
                       work_data work) {

     double *v1;
     double *v2;

     v1 = get_work1(&work, WORK_DX);
     v2 = get_work1(&work, WORK_DX);

     dvec_scale(atran, S_m_a, v1, n_quad_v);
     dvec_add(F_m_a, v1, F_m_a, n_quad_v);

     *atran_a += dvec_dot(F_m, S_m_a, n_quad_v);

     dmat_vxvtmx(S_m_a, F_m, -1., T_m_a, 1., n_quad_v, n_quad_v);

     dmat_gxvxmx(1, T_m, S_m_a, -1., F_m_a, 1., n_quad_v, n_quad_v);

     dmat_vxvtmx(S_m_a, F_p, -atran, R_p_a, 1., n_quad_v, n_quad_v);

     dmat_gxvxmx(1, R_p, S_m_a, -1., v1, 0., n_quad_v, n_quad_v);

     dvec_add(F_p_a, S_p_a, F_p_a, n_quad_v);

     dmat_vxvtmx(S_p_a, F_p, -atran, T_p_a, 1., n_quad_v, n_quad_v);

     dmat_gxvxmx(1, T_p, S_p_a, -1., v1, 1., n_quad_v, n_quad_v);

     dmat_vxvtmx(S_p_a, F_m, -1., R_m_a, 1., n_quad_v, n_quad_v);

     dmat_gxvxmx(1, R_m, S_p_a, -1., F_m_a, 1., n_quad_v, n_quad_v);

     dvec_scale(atran, v1, v2, n_quad_v);
     dvec_add(F_p_a, v2, F_p_a, n_quad_v);

     *atran_a += dvec_dot(F_p, v1, n_quad_v);

     dvec_zero(S_p_a, n_quad_v);
     dvec_zero(S_m_a, n_quad_v);
}



/*******************************************************************************
 *
 ******************************************************************************/
void scale_source_vectors_a(int n_quad_v,
                            double btran, double *btran_a,
                            double *F_p1, double *F_m1,
                            double *F_p2, double *F_m2,
                            double *F_p_a1, double *F_m_a1,
                            double *F_p_a2, double *F_m_a2, work_data work) {

     double *v1;

     v1 = get_work1(&work, WORK_DX);

     *btran_a += dvec_dot(F_p1, F_p_a2, n_quad_v);

     dvec_scale(btran, F_p_a2, v1, n_quad_v);
     dvec_add(F_p_a1, v1, F_p_a1, n_quad_v);

     *btran_a += dvec_dot(F_m1, F_m_a2, n_quad_v);

     dvec_scale(btran, F_m_a2, v1, n_quad_v);
     dvec_add(F_m_a1, v1, F_m_a1, n_quad_v);

     dvec_zero(F_p_a2, n_quad_v);
     dvec_zero(F_m_a2, n_quad_v);
}



/*******************************************************************************
 *
 ******************************************************************************/
#ifdef INCLUDE_DEV_SOURCE
#include "../src2/xrtm_source_a2.c"
#endif
