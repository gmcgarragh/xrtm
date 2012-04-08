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
#include "xrtm_radiance.h"
#include "xrtm_radiance_a.h"
#include "xrtm_utility.h"


/*******************************************************************************
 *
 ******************************************************************************/
void radiance_slab(int n_quad, int n_derivs,
                   double **R_m, double **T_p, double *S_p,
                   double ***R_m_l, double ***T_p_l, double **S_p_l,
                   double *I1_m, double *I2_p, double **I1_m_l,
                   double **I2_p_l, double *I1_p, double **I1_p_l,
                   uchar *derivs, work_data work) {

     int i;

     double *v1;

     v1 = get_work1(&work, WORK_DX);

     dm_v_mul(R_m, I1_m, n_quad, n_quad, I1_p);
     dm_v_mul(T_p, I2_p, n_quad, n_quad, v1);
     dvec_add(I1_p, v1, I1_p, n_quad);
     dvec_add(I1_p, S_p, I1_p, n_quad);

     for (i = 0; i < n_derivs; ++i) {
          if (! derivs[i]) {
               dm_v_mul(R_m, I1_m_l[i], n_quad, n_quad, I1_p_l[i]);
               dm_v_mul(T_p, I2_p_l[i], n_quad, n_quad, v1);
               dvec_add(I1_p_l[i], v1, I1_p_l[i], n_quad);
          }
          else {
               dm_v_mul(R_m_l[i], I1_m, n_quad, n_quad, I1_p_l[i]);
               dm_v_mul(R_m, I1_m_l[i], n_quad, n_quad, v1);
               dvec_add(I1_p_l[i], v1, I1_p_l[i], n_quad);
               dm_v_mul(T_p_l[i], I2_p, n_quad, n_quad, v1);
               dvec_add(I1_p_l[i], v1, I1_p_l[i], n_quad);
               dm_v_mul(T_p, I2_p_l[i], n_quad, n_quad, v1);
               dvec_add(I1_p_l[i], v1, I1_p_l[i], n_quad);

               dvec_add(I1_p_l[i], S_p_l[i], I1_p_l[i], n_quad);
          }
     }
#ifdef USE_AD_FOR_TL_CALC_RADIANCE_SLAB
     radiance_slab_tl_with_ad(n_quad, n_derivs,
                              R_m, T_p, S_p,
                              R_m_l, T_p_l, S_p_l,
                              I1_m, I2_p, I1_m_l,
                              I2_p_l, I1_p, I1_p_l,
                              derivs, work);
#endif
}



/*******************************************************************************
 *
 ******************************************************************************/
void radiance_toa_ref(int n_quad, int n_derivs,
                      double **R_m, double *S_p, double ***R_m_l, double **S_p_l,
                      double *I1_m, double **I1_m_l, double *I1_p, double **I1_p_l,
                      uchar *derivs, work_data work) {

     int i;

     double *v1;

     v1 = get_work1(&work, WORK_DX);

     dm_v_mul(R_m, I1_m, n_quad, n_quad, I1_p);
     dvec_add(I1_p, S_p, I1_p, n_quad);

     for (i = 0; i < n_derivs; ++i) {
          if (! derivs[i])
               dm_v_mul(R_m, I1_m_l[i], n_quad, n_quad, I1_p_l[i]);
          else {
               dm_v_mul(R_m_l[i], I1_m, n_quad, n_quad, I1_p_l[i]);
               dm_v_mul(R_m, I1_m_l[i], n_quad, n_quad, v1);
               dvec_add(I1_p_l[i], v1, I1_p_l[i], n_quad);

               dvec_add(I1_p_l[i], S_p_l[i], I1_p_l[i], n_quad);
          }
     }
}



/*******************************************************************************
 *
 ******************************************************************************/
void radiance_boa_ref(int n_quad, int n_derivs,
        double **R12_p, double *S12_m, double ***R12_p_l, double **S12_m_l,
        double **R23_m, double ***R23_m_l, double *I3_p, double **I3_p_l,
        double *I2_p, double **I2_p_l, double *I2_m, double **I2_m_l,
        uchar *derivs, work_data work) {

     int i;

     int *i1;

     double *v1;
     double *v2;

     double **w1;
     double **w2;
     double **w3;

     i1 = get_work1(&work, WORK_IX);

     v1 = get_work1(&work, WORK_DX);
     v2 = get_work1(&work, WORK_DX);

     w1 = get_work1(&work, WORK_DXX);

     dmat_mul(R12_p, R23_m, n_quad, n_quad, n_quad, w1);
     dmat_i_sub(w1, w1, n_quad);
     dmat_getrf(w1, n_quad, n_quad, i1);

     dm_v_mul(R12_p, I3_p, n_quad, n_quad, I2_m);
     dvec_add(I2_m, S12_m, I2_m, n_quad);

     dmat_getrs(w1, &I2_m, n_quad, 1, i1);

     dm_v_mul(R23_m, I2_m, n_quad, n_quad, I2_p);
     dvec_add(I2_p, I3_p, I2_p, n_quad);
/*
     if (flags_or(derivs, n_derivs)) {
*/
          w2 = get_work1(&work, WORK_DXX);
          w3 = get_work1(&work, WORK_DXX);

          for (i = 0; i < n_derivs; ++i) {
               if (derivs[i] == ADDING_U_U) {
                    dm_v_mul(R12_p, I3_p_l[i], n_quad, n_quad, I2_m_l[i]);

                    dmat_getrs(w1, &I2_m_l[i], n_quad, 1, i1);

                    dm_v_mul(R23_m, I2_m_l[i], n_quad, n_quad, I2_p_l[i]);

                    dvec_add(I2_p_l[i], I3_p_l[i], I2_p_l[i], n_quad);
               }
               else
               if (derivs[i] == ADDING_U_L) {
                    dmat_mul(R12_p, R23_m_l[i], n_quad, n_quad, n_quad, w2);

                    dm_v_mul(w2, I2_m, n_quad, n_quad, v1);

                    dm_v_mul(R12_p, I3_p_l[i], n_quad, n_quad, v2);
                    dvec_add(v1, v2, I2_m_l[i], n_quad);

                    dmat_getrs(w1, &I2_m_l[i], n_quad, 1, i1);

                    dm_v_mul(R23_m_l[i], I2_m, n_quad, n_quad, v1);
                    dm_v_mul(R23_m, I2_m_l[i], n_quad, n_quad, v2);
                    dvec_add(v1, v2, I2_p_l[i], n_quad);

                    dvec_add(I2_p_l[i], I3_p_l[i], I2_p_l[i], n_quad);
               }
               else
               if (derivs[i] == ADDING_L_P) {
                    dmat_mul(R12_p_l[i], R23_m, n_quad, n_quad, n_quad, w2);

                    dm_v_mul(w2, I2_m, n_quad, n_quad, v1);

                    dm_v_mul(R12_p_l[i], I3_p, n_quad, n_quad, v2);
                    dvec_add(v1, v2, v1, n_quad);
                    dm_v_mul(R12_p, I3_p_l[i], n_quad, n_quad, v2);
                    dvec_add(v1, v2, v1, n_quad);
                    dvec_add(v1, S12_m_l[i], I2_m_l[i], n_quad);

                    dmat_getrs(w1, &I2_m_l[i], n_quad, 1, i1);

                    dm_v_mul(R23_m, I2_m_l[i], n_quad, n_quad, I2_p_l[i]);

                    dvec_add(I2_p_l[i], I3_p_l[i], I2_p_l[i], n_quad);
               }
               else
               if (derivs[i] == ADDING_L_L) {
                    dmat_mul(R12_p_l[i], R23_m, n_quad, n_quad, n_quad, w2);
                    dmat_mul(R12_p, R23_m_l[i], n_quad, n_quad, n_quad, w3);
                    dmat_add(w2, w3, w2,  n_quad, n_quad);

                    dm_v_mul(w2, I2_m, n_quad, n_quad, v1);

                    dm_v_mul(R12_p_l[i], I3_p, n_quad, n_quad, v2);
                    dvec_add(v1, v2, v1, n_quad);
                    dm_v_mul(R12_p, I3_p_l[i], n_quad, n_quad, v2);
                    dvec_add(v1, v2, v1, n_quad);
                    dvec_add(v1, S12_m_l[i], I2_m_l[i], n_quad);

                    dmat_getrs(w1, &I2_m_l[i], n_quad, 1, i1);

                    dm_v_mul(R23_m_l[i], I2_m, n_quad, n_quad, v1);
                    dm_v_mul(R23_m, I2_m_l[i], n_quad, n_quad, v2);
                    dvec_add(v1, v2, I2_p_l[i], n_quad);

                    dvec_add(I2_p_l[i], I3_p_l[i], I2_p_l[i], n_quad);
               }
          }
/*
     }
*/
}



/*******************************************************************************
 *
 ******************************************************************************/
void radiance_toa_all(int n_quad, int n_derivs,
                      double  **R_m, double  **T_p, double  *S_p,
                      double ***R_m_l, double ***T_p_l, double **S_p_l,
                      double *I1_m, double *I2_p, double **I1_m_l, double **I2_p_l,
                      double *I1_p, double *I1m_, double **I1_p_l, double **I1m_l_,
                      uchar *derivs, work_data work) {

     int i;

     radiance_slab(n_quad, n_derivs, R_m, T_p, S_p, R_m_l, T_p_l, S_p_l,
                   I1_m, I2_p, I1_m_l, I2_p_l, I1_p, I1_p_l, derivs, work);

     dvec_copy(I1m_, I1_m, n_quad);

     for (i = 0; i < n_derivs; ++i)
          dvec_copy(I1m_l_[i], I1_m_l[i], n_quad);
}



/*******************************************************************************
 *
 ******************************************************************************/
void radiance_boa_all(int n_quad, int n_derivs,
        double  **R12_p, double  **T12_m, double  **R23_m, double  *S12_m,
        double ***R12_p_l, double ***T12_m_l, double ***R23_m_l, double **S12_m_l,
        double *I3_p, double *I1_m, double **I3_p_l, double **I1_m_l,
        double *I2_p, double *I2_m, double **I2_p_l, double **I2_m_l,
        uchar *derivs, save_tree_data save_tree, work_data work) {

     int i;

     int *i1;

     double *v1;
     double *v2;

     double **w1;
     double **w2;
     double **w3;

     forward_save_radiance_boa_all_data *save;


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (save_tree.t) {
          save_tree_encode_s(&save_tree, "radiance_boa");

          if (save_tree_retrieve_data(&save_tree, forward_save_radiance_boa_all_data, &save))
               forward_save_radiance_boa_all_alloc(save, n_quad);
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     i1 = get_work1(&work, WORK_IX);

     v1 = get_work1(&work, WORK_DX);
     v2 = get_work1(&work, WORK_DX);

     w1 = get_work1(&work, WORK_DXX);


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     dmat_mul(R12_p, R23_m, n_quad, n_quad, n_quad, w1);
     dmat_i_sub(w1, w1, n_quad);
     dmat_getrf(w1, n_quad, n_quad, i1);

     if (save_tree.t) {
          copy_array1_i(save->ip, i1, n_quad);

          dmat_copy(save->P, w1, n_quad, n_quad);
     }

     dm_v_mul(R12_p, I3_p, n_quad, n_quad, v1);
     dm_v_mul(T12_m, I1_m, n_quad, n_quad, v2);
     dvec_add(v1, v2, I2_m, n_quad);
     dvec_add(I2_m, S12_m, I2_m, n_quad);

     dmat_getrs(w1, &I2_m, n_quad, 1, i1);

     dm_v_mul(R23_m, I2_m, n_quad, n_quad, I2_p);
     dvec_add(I2_p, I3_p, I2_p, n_quad);
/*
     if (flags_or(derivs, n_derivs)) {
*/
          w2 = get_work1(&work, WORK_DXX);
          w3 = get_work1(&work, WORK_DXX);

          for (i = 0; i < n_derivs; ++i) {
               if (derivs[i] == ADDING_U_U) {
                    dm_v_mul(R12_p, I3_p_l[i], n_quad, n_quad, v1);
                    dm_v_mul(T12_m, I1_m_l[i], n_quad, n_quad, v2);
                    dvec_add(v1, v2, I2_m_l[i], n_quad);

                    dmat_getrs(w1, &I2_m_l[i], n_quad, 1, i1);

                    dm_v_mul(R23_m, I2_m_l[i], n_quad, n_quad, I2_p_l[i]);

                    dvec_add(I2_p_l[i], I3_p_l[i], I2_p_l[i], n_quad);
               }
               else
               if (derivs[i] == ADDING_U_L) {
                    dmat_mul(R12_p, R23_m_l[i], n_quad, n_quad, n_quad, w2);

                    dm_v_mul(w2, I2_m, n_quad, n_quad, v1);

                    dm_v_mul(R12_p, I3_p_l[i], n_quad, n_quad, v2);
                    dvec_add(v1, v2, v1, n_quad);
                    dm_v_mul(T12_m, I1_m_l[i], n_quad, n_quad, v2);
                    dvec_add(v1, v2, I2_m_l[i], n_quad);

                    dmat_getrs(w1, &I2_m_l[i], n_quad, 1, i1);

                    dm_v_mul(R23_m_l[i], I2_m, n_quad, n_quad, v1);
                    dm_v_mul(R23_m, I2_m_l[i], n_quad, n_quad, v2);
                    dvec_add(v1, v2, I2_p_l[i], n_quad);

                    dvec_add(I2_p_l[i], I3_p_l[i], I2_p_l[i], n_quad);
               }
               else
               if (derivs[i] == ADDING_L_P) {
                    dmat_mul(R12_p_l[i], R23_m, n_quad, n_quad, n_quad, w2);

                    dm_v_mul(w2, I2_m, n_quad, n_quad, v1);

                    dm_v_mul(R12_p_l[i], I3_p, n_quad, n_quad, v2);
                    dvec_add(v1, v2, v1, n_quad);
                    dm_v_mul(R12_p, I3_p_l[i], n_quad, n_quad, v2);
                    dvec_add(v1, v2, v1, n_quad);
                    dm_v_mul(T12_m_l[i], I1_m, n_quad, n_quad, v2);
                    dvec_add(v1, v2, v1, n_quad);
                    dm_v_mul(T12_m, I1_m_l[i], n_quad, n_quad, v2);
                    dvec_add(v1, v2, v1, n_quad);
                    dvec_add(v1, S12_m_l[i], I2_m_l[i], n_quad);

                    dmat_getrs(w1, &I2_m_l[i], n_quad, 1, i1);

                    dm_v_mul(R23_m, I2_m_l[i], n_quad, n_quad, I2_p_l[i]);

                    dvec_add(I2_p_l[i], I3_p_l[i], I2_p_l[i], n_quad);
               }
               else
               if (derivs[i] == ADDING_L_L) {
                    dmat_mul(R12_p_l[i], R23_m, n_quad, n_quad, n_quad, w2);
                    dmat_mul(R12_p, R23_m_l[i], n_quad, n_quad, n_quad, w3);
                    dmat_add(w2, w3, w2,  n_quad, n_quad);

                    dm_v_mul(w2, I2_m, n_quad, n_quad, v1);

                    dm_v_mul(R12_p_l[i], I3_p, n_quad, n_quad, v2);
                    dvec_add(v1, v2, v1, n_quad);
                    dm_v_mul(R12_p, I3_p_l[i], n_quad, n_quad, v2);
                    dvec_add(v1, v2, v1, n_quad);
                    dm_v_mul(T12_m_l[i], I1_m, n_quad, n_quad, v2);
                    dvec_add(v1, v2, v1, n_quad);
                    dm_v_mul(T12_m, I1_m_l[i], n_quad, n_quad, v2);
                    dvec_add(v1, v2, v1, n_quad);
                    dvec_add(v1, S12_m_l[i], I2_m_l[i], n_quad);

                    dmat_getrs(w1, &I2_m_l[i], n_quad, 1, i1);

                    dm_v_mul(R23_m_l[i], I2_m, n_quad, n_quad, v1);
                    dm_v_mul(R23_m, I2_m_l[i], n_quad, n_quad, v2);
                    dvec_add(v1, v2, I2_p_l[i], n_quad);

                    dvec_add(I2_p_l[i], I3_p_l[i], I2_p_l[i], n_quad);
               }
          }
/*
     }
*/
}



/*******************************************************************************
 *
 ******************************************************************************/
double radiance_to_mean(int n_quad, double *qx_v, double *qw_v,
                        double F_0, double mu_0, double btran, double *I) {

     int i;

     double a = 0.;

     for (i = 0; i < n_quad; ++i)
          a += qw_v[i] * I[i];

     return a / 2.; // + F_0 * btran / (4. * PI);
}



/*******************************************************************************
 *
 ******************************************************************************/
double radiance_to_flux(int n_quad, double *qx_v, double *qw_v,
                        double F_0, double mu_0, double btran, double *I) {

     int i;

     double a = 0.;

     for (i = 0; i < n_quad; ++i)
          a += qx_v[i] * qw_v[i] * I[i];

     return a * 2. * PI; // + F_0 * mu_0 * btran;
}



/*******************************************************************************
 *
 ******************************************************************************/
#ifdef INCLUDE_DEV_SOURCE
#include "../src2/xrtm_radiance2.c"
#endif
