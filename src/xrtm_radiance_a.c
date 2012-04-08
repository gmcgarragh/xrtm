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
#include "xrtm_adding.h"
#include "xrtm_radiance_a.h"
#include "xrtm_utility.h"


/*******************************************************************************
 *
 ******************************************************************************/
void radiance_slab_a(int n_quad,
                     double **R_m, double **T_p, double *S_p,
                     double **R_m_a, double **T_p_a, double *S_p_a,
                     double *I1_m, double *I2_p, double *I1_m_a,
                     double *I2_p_a, double *I1_p, double *I1_p_a,
                     work_data work) {
/*
     double *v1;

     double **w1;

     v1 = get_work1(&work, WORK_DX);

     w1 = get_work1(&work, WORK_DXX);
*/
     dmat_vxvtmx(I1_p_a, I1_m, 1., R_m_a, 1., n_quad, n_quad);
/*
     dmat_trans(R_m, w1, n_quad, n_quad);
     dm_v_mul(w1, I1_p_a, n_quad, n_quad, v1);
     dvec_add(I1_m_a, v1, I1_m_a, n_quad);
*/
     dmat_gxvxmx(1, R_m, I1_p_a, 1., I1_m_a, 1., n_quad, n_quad);

     dmat_vxvtmx(I1_p_a, I2_p, 1., T_p_a, 1., n_quad, n_quad);
/*
     dmat_trans(T_p, w1, n_quad, n_quad);
     dm_v_mul(w1, I1_p_a, n_quad, n_quad, v1);
     dvec_add(I2_p_a, v1, I2_p_a, n_quad);
*/
     dmat_gxvxmx(1, T_p, I1_p_a, 1., I2_p_a, 1., n_quad, n_quad);

     dvec_add(S_p_a, I1_p_a, S_p_a, n_quad);

     dvec_zero(I1_p_a, n_quad);
}



/*******************************************************************************
 *
 ******************************************************************************/
void radiance_toa_ref_a(int n_quad,
                        double **R_m, double *S_p, double **R_m_a, double *S_p_a,
                        double *I1_m, double *I1_m_a, double *I1_p, double *I1_p_a,
                        work_data work) {
/*
     int i;
     int j;

     double *v1;

     double **w1;

     v1 = get_work1(&work, WORK_DX);

     w1 = get_work1(&work, WORK_DXX);

     for (i = 0; i < n_quad; ++i) {
          for (j = 0; j < n_quad; ++j) {
               R_m_a[i][j] += I1_p_a[i] * I1_m[j];
          }
     }
*/
     dmat_vxvtmx(I1_p_a, I1_m, 1., R_m_a, 1., n_quad, n_quad);

/*
     dmat_trans(R_m, w1, n_quad, n_quad);
     dm_v_mul(w1, I1_p_a, n_quad, n_quad, v1);
     dvec_add(I1_m_a, v1, I1_m_a, n_quad);
*/
     dmat_gxvxmx(1, R_m, I1_p_a, 1., I1_m_a, 1., n_quad, n_quad);

     dvec_add(S_p_a, I1_p_a, S_p_a, n_quad);
}



/*******************************************************************************
 *
 ******************************************************************************/
void radiance_toa_all_a(int n_quad,
                        double **R_m, double **T_p, double *S_p,
                        double **R_m_a, double **T_p_a, double *S_p_a,
                        double *I1_m, double *I2_p, double *I1_m_a, double *I2_p_a,
                        double *I1_p, double *I1m_, double *I1_p_a, double *I1m_a_,
                        work_data work) {

     dvec_add(I1_m_a, I1m_a_, I1_m_a, n_quad);

     radiance_slab_a(n_quad, R_m, T_p, S_p, R_m_a, T_p_a, S_p_a,
                     I1_m, I2_p, I1_m_a, I2_p_a, I1_p, I1_p_a, work);
}



/*******************************************************************************
 *
 ******************************************************************************/
static void forward_save_radiance_boa_all_free(forward_save_radiance_boa_all_data *d);

int forward_save_radiance_boa_all_alloc(forward_save_radiance_boa_all_data *d, int n_quad) {

     d->free = (void (*)(void *)) forward_save_radiance_boa_all_free;

     d->ip = alloc_array1_i(n_quad);

     d->P  = alloc_array2_d(n_quad, n_quad);

     return 0;
}



static void forward_save_radiance_boa_all_free(forward_save_radiance_boa_all_data *d) {

     free_array1_i(d->ip);

     free_array2_d(d->P);
}



/*******************************************************************************
 *
 ******************************************************************************/
void radiance_boa_all_a(int n_quad,
        double **R12_p, double **T12_m, double **R23_m, double *S12_m,
        double **R12_p_a, double **T12_m_a, double **R23_m_a, double *S12_m_a,
        double *I3_p, double *I1_m, double *I3_p_a, double *I1_m_a,
        double *I2_p, double *I2_m, double *I2_p_a, double *I2_m_a,
        save_tree_data save_tree, work_data work) {
/*
     int i;
     int j;

     double *v1;
*/
     double *t;
/*
     double **w1;
     double **w2;
*/
     double **Q_a;

     forward_save_radiance_boa_all_data *save;


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     save_tree_encode_s(&save_tree, "radiance_boa");

     save_tree_retrieve_data(&save_tree, forward_save_radiance_boa_all_data, &save);


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
/*
     v1  = get_work1(&work, WORK_DX);
*/
     t   = get_work1(&work, WORK_DX);
/*
     w1  = get_work1(&work, WORK_DXX);
     w2  = get_work1(&work, WORK_DXX);
*/
     Q_a = get_work1(&work, WORK_DXX);


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
/*
     for (i = 0; i < n_quad; ++i) {
          for (j = 0; j < n_quad; ++j) {
               R23_m_a[i][j] += I2_p_a[i] * I2_m[j];
          }
     }
*/
     dmat_vxvtmx(I2_p_a, I2_m, 1., R23_m_a, 1., n_quad, n_quad);
/*
     dmat_trans(R23_m, w1, n_quad, n_quad);
     dm_v_mul(w1, I2_p_a, n_quad, n_quad, v1);
     dvec_add(I2_m_a, v1, I2_m_a, n_quad);
*/
     dmat_gxvxmx(1, R23_m, I2_p_a, 1., I2_m_a, 1., n_quad, n_quad);

     dvec_add(I3_p_a, I2_p_a, I3_p_a, n_quad);

     dvec_copy(t, I2_m_a, n_quad);
     dmat_getrs2('t', save->P, &t, n_quad, 1, save->ip);
/*
     for (i = 0; i < n_quad; ++i) {
          for (j = 0; j < n_quad; ++j) {
               Q_a[i][j] = t[i] * I2_m[j];
          }
     }
*/
     dmat_vxvtmx(t, I2_m, 1., Q_a, 0., n_quad, n_quad);
/*
     for (i = 0; i < n_quad; ++i) {
          for (j = 0; j < n_quad; ++j) {
               R12_p_a[i][j] += t[i] * I3_p[j];
          }
     }
*/
     dmat_vxvtmx(t, I3_p, 1., R12_p_a, 1., n_quad, n_quad);
/*
     dmat_trans(R12_p, w1, n_quad, n_quad);
     dm_v_mul(w1, t, n_quad, n_quad, v1);
     dvec_add(I3_p_a, v1, I3_p_a, n_quad);
*/
     dmat_gxvxmx(1, R12_p, t, 1., I3_p_a, 1., n_quad, n_quad);
/*
     for (i = 0; i < n_quad; ++i) {
          for (j = 0; j < n_quad; ++j) {
               T12_m_a[i][j] += t[i] * I1_m[j];
          }
     }
*/
     dmat_vxvtmx(t, I1_m, 1., T12_m_a, 1., n_quad, n_quad);
/*
     dmat_trans(T12_m, w1, n_quad, n_quad);
     dm_v_mul(w1, t, n_quad, n_quad, v1);
     dvec_add(I1_m_a, v1, I1_m_a, n_quad);
*/
     dmat_gxvxmx(1, T12_m, t, 1., I1_m_a, 1., n_quad, n_quad);

     dvec_add(S12_m_a, t, S12_m_a, n_quad);
/*
     dmat_trans(R23_m, w1, n_quad, n_quad);
     dmat_mul(Q_a, w1, n_quad, n_quad, n_quad, w2);
     dmat_add(R12_p_a, w2, R12_p_a, n_quad, n_quad);
*/
     dmat_gxgxmx(0, Q_a, 1, R23_m, 1., R12_p_a, 1., n_quad, n_quad, n_quad);
/*
     dmat_trans(R12_p, w1, n_quad, n_quad);
     dmat_mul(w1, Q_a, n_quad, n_quad, n_quad, w2);
     dmat_add(R23_m_a, w2, R23_m_a, n_quad, n_quad);
*/
     dmat_gxgxmx(1, R12_p, 0, Q_a, 1., R23_m_a, 1., n_quad, n_quad, n_quad);
}



/*******************************************************************************
 *
 ******************************************************************************/
#ifdef INCLUDE_DEV_SOURCE
#include "../src2/xrtm_radiance_a2.c"
#endif
