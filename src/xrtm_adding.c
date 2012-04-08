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
#include "xrtm_adding_a.h"
#include "xrtm_matrix.h"
#include "xrtm_save_tree.h"
#include "xrtm_scatter.h"
#include "xrtm_utility.h"


/*******************************************************************************
 *
 ******************************************************************************/
int layer_add_aux_alloc(layer_add_aux_data *d, int n) {

     d->P13 = alloc_array2_d(n, n);
     d->P31 = alloc_array2_d(n, n);
     d->A13 = alloc_array2_d(n, n);
     d->A31 = alloc_array2_d(n, n);
     d->B13 = alloc_array2_d(n, n);
     d->B31 = alloc_array2_d(n, n);

     if (! (d->P13 && d->P31 && d->A13 && d->A31 && d->B13 && d->B31)) {
          eprintf("ERROR: alloc_array2_d(%d, %d)\n", n, n);
          return -1;
     }

     d->C13 = alloc_array1_d(n);
     d->C31 = alloc_array1_d(n);

     if (! (d->C13 && d->C31)) {
          eprintf("ERROR: alloc_array1_d(%d)\n", n);
          return -1;
     }

     return 0;
}



void layer_add_aux_free(layer_add_aux_data *d) {

     free_array2_d(d->P13);
     free_array2_d(d->P31);
     free_array2_d(d->A13);
     free_array2_d(d->A31);
     free_array2_d(d->B13);
     free_array2_d(d->B31);
     free_array1_d(d->C13);
     free_array1_d(d->C31);
}



/*******************************************************************************
 *
 ******************************************************************************/
void layer_zero_ref(double **R_m, double *S_p,
                    int n) {

     layer_zero_ref_l(R_m, S_p,
                     NULL, NULL, n, 0, NULL, NULL);
}



/*******************************************************************************
 *
 ******************************************************************************/
void layer_zero_ref_l(double **R_m, double *S_p,
                      double ***R_m_l, double **S_p_l,
                      int n, int n_derivs, uchar *derivs_h, uchar *derivs_p) {

     int i;

     dmat_zero(R_m, n, n);
     dvec_zero(S_p, n);

     for (i = 0; i < n_derivs; ++i) {
          if (derivs_h[i])
               dmat_zero(R_m_l[i], n, n);
          if (derivs_p[i])
               dvec_zero(S_p_l[i], n);
     }
}



/*******************************************************************************
 *
 ******************************************************************************/
void layer_copy_ref(double **R_m2, double *S_p2,
                    double **R_m1, double *S_p1,
                    int n) {

     layer_copy_ref_l(R_m2, S_p2, R_m1, S_p1,
                     NULL, NULL, NULL, NULL, n, 0, NULL, NULL);
}



/*******************************************************************************
 *
 ******************************************************************************/
void layer_copy_ref_l(double **R_m2, double *S_p2,
                      double **R_m1, double *S_p1,
                      double ***R_m_l2, double **S_p_l2,
                      double ***R_m_l1, double **S_p_l1,
                      int n, int n_derivs, uchar *derivs_h, uchar *derivs_p) {

     int i;

     dmat_copy(R_m2, R_m1, n, n);
     dvec_copy(S_p2, S_p1, n);

     for (i = 0; i < n_derivs; ++i) {
          if (derivs_h[i])
               dmat_copy(R_m_l2[i], R_m_l1[i], n, n);
          if (derivs_p[i])
               dvec_copy(S_p_l2[i], S_p_l1[i], n);
     }
}



/*******************************************************************************
 * 
 ******************************************************************************/
static void add_ref_u_p(double **AXX, double **R23, double *S32_l, double *S31_l, int n, double atran, double *v1, double *v2) {

     /* S31_l = A31 * S32_l * t12  */
     dvec_scale(atran, S32_l, v1, n);
     dm_v_mul(AXX, v1, n, n, S31_l);
}



static void add_ref_p_u(double **AXX, double **R23, double *S12_l, double *S21_l, double *S31_l, int n, double atran, double *v1, double *v2) {

     /* S31_l = S21_l + A31 * R23 * S12_l  */
     dm_v_mul(R23, S12_l, n, n, v1);

     dm_v_mul(AXX, v1, n, n, v2);

     dvec_add(S21_l, v2, S31_l, n);
}



static void add_ref_p_p(double **AXX, double **R23, double *S12_l, double *S21_l, double *S32_l, double *S31_l, int n, double atran, double *v1, double *v2) {

     /* S31_l = S21_l + A31 * (S32_l * t12 + R23 * S12_l)  */
     dvec_scale(atran, S32_l, v1, n);
     dm_v_mul(R23, S12_l, n, n, v2);
     dvec_add(v1, v2, v1, n);

     dm_v_mul(AXX, v1, n, n, v2);

     dvec_add(S21_l, v2, S31_l, n);
}



static void add_ref_u_l(double **PXX, int *i1, double **AXX, double **BXX, double *CXX, double **T12, double *S12, double **R21, double **R23, double **R23_l, double *S32_l, double **R13_l, double *S31_l, int n, double atran, double *v1, double *v2, double **w1, double **w2, work_data work) {

     double **DXX;
     double **EXX;

     DXX = get_work1(&work, WORK_DXX);
     EXX = get_work1(&work, WORK_DXX);

     /* D31 = A31 * R23_l */
     dmat_mul(AXX, R23_l, n, n, n, DXX);

     /* E31 = D31 * R21 * P31 */
     dmat_mul(DXX, R21, n, n, n, EXX);
     dmat_getrs(PXX, EXX, n, n, i1);

     /* S31_l = E31 * C31 + A31 * (S32_l * t12 + R23_l * S12) */
     dvec_scale(atran, S32_l, v1, n);
     dm_v_mul(R23_l, S12, n, n, v2);
     dvec_add(v1, v2, v1, n);

     dm_v_mul(AXX, v1, n, n, v2);

     dm_v_mul(EXX, CXX, n, n, v1);
     dvec_add(v1, v2, S31_l, n);

     /* R13_l = E31 * B13 + D31 * T12 */
     dmat_mul(EXX, BXX, n, n, n, w1);
     dmat_mul(DXX, T12, n, n, n, w2);
     dmat_add(w1, w2, R13_l, n, n);
}



static void add_ref_l_u(double **PXX, int *i1, double **AXX, double **BXX, double *CXX, double **R23, double *S32, double **R12_l, double **T12_l, double *S12_l, double **R21_l, double **T21_l, double *S21_l, double **R13_l, double *S31_l, int n, double atran, double atran_l, double *v1, double *v2, double **w1, double **w2, work_data work) {

     double **DXX;
     double **EXX;

     DXX = get_work1(&work, WORK_DXX);
     EXX = get_work1(&work, WORK_DXX);

     /* D31 = A31 * R23 */
     dmat_mul(AXX, R23, n, n, n, DXX);

     /* E31 = (T21_l + D31 * R21_l) * P31 */
     dmat_mul(DXX, R21_l, n, n, n, EXX);
     dmat_add(T21_l, EXX, EXX, n, n);
     dmat_getrs(PXX, EXX, n, n, i1);

     /* S31_l = S21_l + E31 * C31 + A31 * (S32 * L(t12) + R23 * S12_l)  */
     dvec_scale(atran_l, S32, v1, n);
     dm_v_mul(R23, S12_l, n, n, v2);
     dvec_add(v1, v2, v1, n);

     dm_v_mul(AXX, v1, n, n, v2);

     dvec_add(S21_l, v2, v1, n);

     dm_v_mul(EXX, CXX, n, n, v2);
     dvec_add(v1, v2, S31_l, n);

     /* R13_l = R12_l + E31 * B13 + D31 * T12_l */
     dmat_mul(EXX, BXX, n, n, n, w1);
     dmat_mul(DXX, T12_l, n, n, n, w2);
     dmat_add(w1, w2, w1, n, n);
     dmat_add(R12_l, w1, R13_l, n, n);
}



static void add_ref_p_l(double **PXX, int *i1, double **AXX, double **BXX, double *CXX, double **T12, double *S12, double **R21, double **R23, double *S12_l, double *S21_l, double **R23_l, double *S32_l, double **R13_l, double *S31_l, int n, double atran, double *v1, double *v2, double **w1, double **w2, work_data work) {

     double **DXX;
     double **EXX;

     DXX = get_work1(&work, WORK_DXX);
     EXX = get_work1(&work, WORK_DXX);

     /* D31 = A31 * R23_l */
     dmat_mul(AXX, R23_l, n, n, n, DXX);

     /* E31 = D31 * R21 * P31 */
     dmat_mul(DXX, R21, n, n, n, EXX);
     dmat_getrs(PXX, EXX, n, n, i1);

     /* S31_l = S21_l + E31 * C31 + A31 * (S32_l * t12 + R23_l * S12 + R23 * S12_l)  */
     dvec_scale(atran, S32_l, v1, n);
     dm_v_mul(R23_l, S12, n, n, v2);
     dvec_add(v1, v2, v1, n);
     dm_v_mul(R23, S12_l, n, n, v2);
     dvec_add(v1, v2, v1, n);

     dm_v_mul(AXX, v1, n, n, v2);

     dvec_add(S21_l, v2, v1, n);

     dm_v_mul(EXX, CXX, n, n, v2);
     dvec_add(v1, v2, S31_l, n);

     /* R13_l = E31 * B13 + D31 * T12 */
     dmat_mul(EXX, BXX, n, n, n, w1);
     dmat_mul(DXX, T12, n, n, n, w2);
     dmat_add(w1, w2, R13_l, n, n);
}



static void add_ref_l_p(double **PXX, int *i1, double **AXX, double **BXX, double *CXX, double **R23, double *S32, double **R12_l, double **T12_l, double *S12_l, double **R21_l, double **T21_l, double *S21_l, double *S32_l, double **R13_l, double *S31_l, int n, double atran, double atran_l, double *v1, double *v2, double **w1, double **w2, work_data work) {

     double **DXX;
     double **EXX;

     DXX = get_work1(&work, WORK_DXX);
     EXX = get_work1(&work, WORK_DXX);

     /* D31 = A31 * R23 */
     dmat_mul(AXX, R23, n, n, n, DXX);

     /* E31 = (T21_l + D31 * R21_l) * P31 */
     dmat_mul(DXX, R21_l, n, n, n, EXX);
     dmat_add(T21_l, EXX, EXX, n, n);
     dmat_getrs(PXX, EXX, n, n, i1);

     /* S31_l = S21_l + E31 * C31 + A31 * (S32_l * t12 + S32 * L(t12) + R23 * S12_l) */
     dvec_scale(atran, S32_l, v1, n);
     dvec_scale(atran_l, S32, v2, n);
     dvec_add(v1, v2, v1, n);
     dm_v_mul(R23, S12_l, n, n, v2);
     dvec_add(v1, v2, v1, n);

     dm_v_mul(AXX, v1, n, n, v2);

     dvec_add(S21_l, v2, v1, n);

     dm_v_mul(EXX, CXX, n, n, v2);
     dvec_add(v1, v2, S31_l, n);

     /* R13_l = R12_l + E31 * B13 + D31 * T12_l */
     dmat_mul(EXX, BXX, n, n, n, w1);
     dmat_mul(DXX, T12_l, n, n, n, w2);
     dmat_add(w1, w2, w1, n, n);

     dmat_add(R12_l, w1, R13_l, n, n);
}



static void add_ref_l_l(double **PXX, int *i1, double **AXX, double **BXX, double *CXX, double **T12, double *S12, double **R21, double **R23, double *S32, double **R12_l, double **T12_l, double *S12_l, double **R21_l, double **T21_l, double *S21_l, double **R23_l, double *S32_l, double **R13_l, double *S31_l, int n, double atran, double atran_l, double *v1, double *v2, double **w1, double **w2, work_data work) {

     double **DXX;

     DXX = get_work1(&work, WORK_DXX);

     /* D31 = (T21_l + A31 * (R23_l * R21 + R23 * R21_l)) * P31 */
     dmat_mul(R23_l, R21, n, n, n, w1);
     dmat_mul(R23, R21_l, n, n, n, w2);
     dmat_add(w1, w2, w1, n, n);
     dmat_mul(AXX, w1, n, n, n, w2);

     dmat_add(T21_l, w2, DXX, n, n);

     dmat_getrs(PXX, DXX, n, n, i1);

     /* S31_l = S21_l + D31 * C31 + A31 * (S32_l * t12 + S32 * L(t12) + R23_l * S12 + R23 * S12_l)  */
     dvec_scale(atran, S32_l, v1, n);
     dvec_scale(atran_l, S32, v2, n);
     dvec_add(v1, v2, v1, n);
     dm_v_mul(R23_l, S12, n, n, v2);
     dvec_add(v1, v2, v1, n);
     dm_v_mul(R23, S12_l, n, n, v2);
     dvec_add(v1, v2, v1, n);

     dm_v_mul(AXX, v1, n, n, v2);

     dvec_add(S21_l, v2, v1, n);

     dm_v_mul(DXX, CXX, n, n, v2);
     dvec_add(v1, v2, S31_l, n);

     /* R13_l = R12_l + D31 * B13 + A31 * (R23_l * T12 + R23 * T12_l) */
     dmat_mul(R23_l, T12, n, n, n, w1);
     dmat_mul(R23, T12_l, n, n, n, w2);
     dmat_add(w1, w2, w1, n, n);
     dmat_mul(AXX, w1, n, n, n, w2);
     dmat_add(R12_l, w2, w1, n, n);
     dmat_mul(DXX, BXX, n, n, n, w2);
     dmat_add(w1, w2, R13_l, n, n);
}



/*******************************************************************************
 *
 ******************************************************************************/
void layer_add_ref(double **R12, double **T12, double *S12,
                   double **R21, double **T21, double *S21,
                   double **R23, double *S32,
                   double **R13, double *S31,
                   double ***R12_l, double ***T12_l, double **S12_l,
                   double ***R21_l, double ***T21_l, double **S21_l,
                   double ***R23_l, double **S32_l,
                   double ***R13_l, double **S31_l,
                   int n, int n_derivs, double atran, double *atran_l,
                   uchar *derivs, int flag, int flag2,
                   save_tree_data save_tree, work_data work) {

     int i;

     int *i1;

     double *v1;
     double *v2;

     double *CXX;

     double **w1;
     double **w2;

     double **PXX;
     double **AXX;
     double **BXX;

     forward_save_layer_add_ref_data *save;


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (save_tree.t) {
          save_tree_encode_s(&save_tree, "layer_add_ref");

          if (save_tree_retrieve_data(&save_tree, forward_save_layer_add_ref_data, &save))
               forward_save_layer_add_ref_alloc(save, n);
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     i1  = get_work1(&work, WORK_IX);

     CXX = get_work1(&work, WORK_DX);

     v1  = get_work1(&work, WORK_DX);
     v2  = get_work1(&work, WORK_DX);

     w1  = get_work1(&work, WORK_DXX);

     if (n_derivs > 0)
          w2  = get_work1(&work, WORK_DXX);

     PXX = get_work1(&work, WORK_DXX);
     AXX = get_work1(&work, WORK_DXX);
     BXX = get_work1(&work, WORK_DXX);


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     /* w1 = LU of (E - R23 * R21)^T */
     dmat_mul(R23, R21, n, n, n, PXX);
     dmat_i_sub(PXX, w1, n);
     dmat_trans(w1, PXX, n, n);
     dmat_getrf(PXX, n, n, i1);

     /* A31 = T21 * P31 */
     dmat_copy(AXX, T21, n, n);
     dmat_getrs(PXX, AXX, n, n, i1);

     /* B13 = R23 * T12 */
     dmat_mul(R23, T12, n, n, n, BXX);

     /* C31 = S32 * t12 + R23 * S12 */
     dvec_scale(atran, S32, v1, n);
     dm_v_mul(R23, S12, n, n, v2);
     dvec_add(v1, v2, CXX, n);

     for (i = 0; i < n_derivs; ++i) {
          if ((! flag && derivs[i] == ADDING_U_P) || (flag && derivs[i] == ADDING_P_U))
               add_ref_u_p(AXX, R23, S32_l[i], S31_l[i], n, atran, v1, v2);
          else
          if ((! flag && derivs[i] == ADDING_P_U) || (flag && derivs[i] == ADDING_U_P))
               add_ref_p_u(AXX, R23, S12_l[i], S21_l[i], S31_l[i], n, atran, v1, v2);
          else
          if (derivs[i] == ADDING_P_P)
               add_ref_p_p(AXX, R23, S12_l[i], S21_l[i], S32_l[i], S31_l[i], n, atran, v1, v2);
          else
          if ((! flag && derivs[i] == ADDING_U_L) || (flag && derivs[i] == ADDING_L_U))
               add_ref_u_l(PXX, i1, AXX, BXX, CXX, T12, S12, R21, R23, R23_l[i], S32_l[i], R13_l[i], S31_l[i], n, atran, v1, v2, w1, w2, work);
          else
          if ((! flag && derivs[i] == ADDING_L_U) || (flag && derivs[i] == ADDING_U_L))
               add_ref_l_u(PXX, i1, AXX, BXX, CXX, R23, S32, R12_l[i], T12_l[i], S12_l[i], R21_l[i], T21_l[i], S21_l[i], R13_l[i], S31_l[i], n, atran, atran_l[i], v1, v2, w1, w2, work);
          else
          if ((! flag && derivs[i] == ADDING_P_L) || (flag && derivs[i] == ADDING_L_P))
               add_ref_p_l(PXX, i1, AXX, BXX, CXX, T12, S12, R21, R23, S12_l[i], S21_l[i], R23_l[i], S32_l[i], R13_l[i], S31_l[i], n, atran, v1, v2, w1, w2, work);
          else
          if ((! flag && derivs[i] == ADDING_L_P) || (flag && derivs[i] == ADDING_P_L))
               add_ref_l_p(PXX, i1, AXX, BXX, CXX, R23, S32, R12_l[i], T12_l[i], S12_l[i], R21_l[i], T21_l[i], S21_l[i], S32_l[i], R13_l[i], S31_l[i], n, atran, atran_l[i], v1, v2, w1, w2, work);
          else
          if (derivs[i] == ADDING_L_L)
               add_ref_l_l(PXX, i1, AXX, BXX, CXX, T12, S12, R21, R23, S32, R12_l[i], T12_l[i], S12_l[i], R21_l[i], T21_l[i], S21_l[i], R23_l[i], S32_l[i], R13_l[i], S31_l[i], n, atran, atran_l[i], v1, v2, w1, w2, work);
     }
if (flag2) {
     /* S31 = S21       + A31 * C31 */
     dm_v_mul(AXX, CXX, n, n, v1);
     dvec_add(S21, v1, S31, n);

     /* R13 = R12 + A31 * B13 */
     dmat_mul(AXX, BXX, n, n, n, w1);
     dmat_add(R12, w1, R13, n, n);
}

     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (save_tree.t) {
          memcpy(save->i13, i1, n * sizeof(int));
          dvec_copy(save->C13, CXX, n);
          dmat_copy(save->P13, PXX, n, n);
          dmat_copy(save->A13, AXX, n, n);
          dmat_copy(save->B31, BXX, n, n);
     }
}



/*******************************************************************************
 *
 ******************************************************************************/
void layer_add_ref2(double **R12, double **T12, double *S12,
                    double **R21, double **T21, double *S21,
                    double **R23, double *S32,
                    double **R13, double *S31,
                    double ***R12_l, double ***T12_l, double **S12_l,
                    double ***R21_l, double ***T21_l, double **S21_l,
                    double ***R23_l, double **S32_l,
                    double ***R13_l, double **S31_l,
                    int n, int n_derivs, double atran, double *atran_l,
                    uchar *derivs, int flag, int flag2, int flag3,
                    work_data work, layer_add_aux_data *d) {

     int i;

     int *i1;

     double *v1;
     double *v2;

     double *CXX;

     double **w1;
     double **w2;

     double **PXX;
     double **AXX;
     double **BXX;

     i1  = get_work1(&work, WORK_IX);

     CXX = get_work1(&work, WORK_DX);

     v1  = get_work1(&work, WORK_DX);
     v2  = get_work1(&work, WORK_DX);

     w1  = get_work1(&work, WORK_DXX);

     if (n_derivs > 0)
          w2  = get_work1(&work, WORK_DXX);
if (! d) {
     PXX = get_work1(&work, WORK_DXX);
     AXX = get_work1(&work, WORK_DXX);
     BXX = get_work1(&work, WORK_DXX);
}

     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
if (d) {
     PXX = d->P31;
     AXX = d->A31;
     BXX = d->B13;
}
if (! d || flag3) {
     /* w1 = LU of (E - R23 * R21)^T */
     dmat_mul(R23, R21, n, n, n, PXX);
     dmat_i_sub(PXX, w1, n);
     dmat_trans(w1, PXX, n, n);
     dmat_getrf(PXX, n, n, i1);

     /* A31 = T21 * P31 */
     dmat_copy(AXX, T21, n, n);
     dmat_getrs(PXX, AXX, n, n, i1);

     /* B13 = R23 * T12 */
     dmat_mul(R23, T12, n, n, n, BXX);
}
     /* C31 = S32 * t12 + R23 * S12 */
     dvec_scale(atran, S32, v1, n);
     dm_v_mul(R23, S12, n, n, v2);
     dvec_add(v1, v2, CXX, n);

     for (i = 0; i < n_derivs; ++i) {
          if ((! flag && derivs[i] == ADDING_U_P) || (flag && derivs[i] == ADDING_P_U))
               add_ref_u_p(AXX, R23, S32_l[i], S31_l[i], n, atran, v1, v2);
          else
          if ((! flag && derivs[i] == ADDING_P_U) || (flag && derivs[i] == ADDING_U_P))
               add_ref_p_u(AXX, R23, S12_l[i], S21_l[i], S31_l[i], n, atran, v1, v2);
          else
          if (derivs[i] == ADDING_P_P)
               add_ref_p_p(AXX, R23, S12_l[i], S21_l[i], S32_l[i], S31_l[i], n, atran, v1, v2);
          else
          if ((! flag && derivs[i] == ADDING_U_L) || (flag && derivs[i] == ADDING_L_U))
               add_ref_u_l(PXX, i1, AXX, BXX, CXX, T12, S12, R21, R23, R23_l[i], S32_l[i], R13_l[i], S31_l[i], n, atran, v1, v2, w1, w2, work);
          else
          if ((! flag && derivs[i] == ADDING_L_U) || (flag && derivs[i] == ADDING_U_L))
               add_ref_l_u(PXX, i1, AXX, BXX, CXX, R23, S32, R12_l[i], T12_l[i], S12_l[i], R21_l[i], T21_l[i], S21_l[i], R13_l[i], S31_l[i], n, atran, atran_l[i], v1, v2, w1, w2, work);
          else
          if ((! flag && derivs[i] == ADDING_P_L) || (flag && derivs[i] == ADDING_L_P))
               add_ref_p_l(PXX, i1, AXX, BXX, CXX, T12, S12, R21, R23, S12_l[i], S21_l[i], R23_l[i], S32_l[i], R13_l[i], S31_l[i], n, atran, v1, v2, w1, w2, work);
          else
          if ((! flag && derivs[i] == ADDING_L_P) || (flag && derivs[i] == ADDING_P_L))
               add_ref_l_p(PXX, i1, AXX, BXX, CXX, R23, S32, R12_l[i], T12_l[i], S12_l[i], R21_l[i], T21_l[i], S21_l[i], S32_l[i], R13_l[i], S31_l[i], n, atran, atran_l[i], v1, v2, w1, w2, work);
          else
          if (derivs[i] == ADDING_L_L)
               add_ref_l_l(PXX, i1, AXX, BXX, CXX, T12, S12, R21, R23, S32, R12_l[i], T12_l[i], S12_l[i], R21_l[i], T21_l[i], S21_l[i], R23_l[i], S32_l[i], R13_l[i], S31_l[i], n, atran, atran_l[i], v1, v2, w1, w2, work);
     }
if (flag2) {
     /* S31 = S21       + A31 * C31 */
     dm_v_mul(AXX, CXX, n, n, v1);
     dvec_add(S21, v1, S31, n);

     /* R13 = R12 + A31 * B13 */
     dmat_mul(AXX, BXX, n, n, n, w1);
     dmat_add(R12, w1, R13, n, n);
}
}



/*******************************************************************************
 *
 ******************************************************************************/
void layer_zero_all(double **R_m, double **T_m, double *S_m,
                    double **R_p, double **T_p, double *S_p,
                    int n) {

     layer_zero_all_l(R_m, T_m, S_m, R_p, T_p, S_p,
                      NULL, NULL, NULL, NULL, NULL, NULL, n, 0, NULL);
}



/*******************************************************************************
 *
 ******************************************************************************/
void layer_zero_all_l(double **R_m, double **T_m, double *S_m,
                      double **R_p, double **T_p, double *S_p,
                      double ***R_m_l, double ***T_m_l, double **S_m_l,
                      double ***R_p_l, double ***T_p_l, double **S_p_l,
                      int n, int n_derivs, uchar *derivs) {

     int i;

     dmat_zero(R_m, n, n);
     dmat_zero(T_m, n, n);
     dvec_zero(S_m, n);
     dmat_zero(R_p, n, n);
     dmat_zero(T_p, n, n);
     dvec_zero(S_p, n);

     for (i = 0; i < n_derivs; ++i) {
          if (! derivs[i])
               continue;

          dmat_zero(R_m_l[i], n, n);
          dmat_zero(T_m_l[i], n, n);
          dvec_zero(S_m_l[i], n);
          dmat_zero(R_p_l[i], n, n);
          dmat_zero(T_p_l[i], n, n);
          dvec_zero(S_p_l[i], n);
     }
}



/*******************************************************************************
 *
 ******************************************************************************/
void layer_copy_all(double **R_m2, double **T_m2, double *S_m2,
                    double **R_p2, double **T_p2, double *S_p2,
                    double **R_m1, double **T_m1, double *S_m1,
                    double **R_p1, double **T_p1, double *S_p1,
                    int n) {

     layer_copy_all_l(R_m2, T_m2, S_m2, R_p2, T_p2, S_p2,
                      R_m1, T_m1, S_m1, R_p1, T_p1, S_p1,
                      NULL, NULL, NULL, NULL, NULL, NULL,
                      NULL, NULL, NULL, NULL, NULL, NULL, n, 0, NULL);
}



/*******************************************************************************
 *
 ******************************************************************************/
void layer_copy_all_l(double **R_m2, double **T_m2, double *S_m2,
                      double **R_p2, double **T_p2, double *S_p2,
                      double **R_m1, double **T_m1, double *S_m1,
                      double **R_p1, double **T_p1, double *S_p1,
                      double ***R_m_l2, double ***T_m_l2, double **S_m_l2,
                      double ***R_p_l2, double ***T_p_l2, double **S_p_l2,
                      double ***R_m_l1, double ***T_m_l1, double **S_m_l1,
                      double ***R_p_l1, double ***T_p_l1, double **S_p_l1,
                      int n, int n_derivs, uchar *derivs) {

     int i;

     dmat_copy(R_m2, R_m1, n, n);
     dmat_copy(T_m2, T_m1, n, n);
     dvec_copy(S_m2, S_m1, n);
     dmat_copy(R_p2, R_p1, n, n);
     dmat_copy(T_p2, T_p1, n, n);
     dvec_copy(S_p2, S_p1, n);

     for (i = 0; i < n_derivs; ++i) {
          if (! derivs[i])
               continue;

          dmat_copy(R_m_l2[i], R_m_l1[i], n, n);
          dmat_copy(T_m_l2[i], T_m_l1[i], n, n);
          dvec_copy(S_m_l2[i], S_m_l1[i], n);
          dmat_copy(R_p_l2[i], R_p_l1[i], n, n);
          dmat_copy(T_p_l2[i], T_p_l1[i], n, n);
          dvec_copy(S_p_l2[i], S_p_l1[i], n);
     }
}



/*******************************************************************************
 * 
 ******************************************************************************/
static void add_all_u_p_up(double **AXX, double **R23, double *S23_l, double *S32_l, double *S31_l, int n, double atran, double *v1, double *v2) {

     /* S31_l = A31 * S32_l * t12  */
     dvec_scale(atran, S32_l, v1, n);
     dm_v_mul(AXX, v1, n, n, S31_l);
}



static void add_all_p_u_up(double **AXX, double **R23, double *S12_l, double *S21_l, double *S31_l, int n, double atran, double *v1, double *v2) {

     /* S31_l = S21_l + A31 * R23 * S12_l  */
     dm_v_mul(R23, S12_l, n, n, v1);

     dm_v_mul(AXX, v1, n, n, v2);

     dvec_add(S21_l, v2, S31_l, n);
}



static void add_all_p_p_up(double **AXX, double **R23, double *S12_l, double *S21_l, double *S23_l, double *S32_l, double *S31_l, int n, double atran, double *v1, double *v2) {

     /* S31_l = S21_l + A31 * (S32_l * t12 + R23 * S12_l)  */
     dvec_scale(atran, S32_l, v1, n);
     dm_v_mul(R23, S12_l, n, n, v2);
     dvec_add(v1, v2, v1, n);

     dm_v_mul(AXX, v1, n, n, v2);

     dvec_add(S21_l, v2, S31_l, n);
}



static void add_all_u_l_up(double **PXX, int *i1, double **AXX, double **BXX, double *CXX, double **T12, double *S12, double **R21, double **T32, double **R23_l, double **T32_l, double *S32_l, double **R13_l, double **T31_l, double *S31_l, int n, double atran, double *v1, double *v2, double **w1, double **w2, work_data work) {

     double **DXX;
     double **EXX;

     DXX = get_work1(&work, WORK_DXX);
     EXX = get_work1(&work, WORK_DXX);

     /* D31 = A31 * R23_l */
     dmat_mul(AXX, R23_l, n, n, n, DXX);

     /* E31 = D31 * R21 * P31 */
     dmat_mul(DXX, R21, n, n, n, EXX);
     dmat_getrs(PXX, EXX, n, n, i1);

     /* S31_l = E31 * C31 + A31 * (S32_l * t12 + R23_l * S12) */
     dvec_scale(atran, S32_l, v1, n);
     dm_v_mul(R23_l, S12, n, n, v2);
     dvec_add(v1, v2, v1, n);

     dm_v_mul(AXX, v1, n, n, v2);

     dm_v_mul(EXX, CXX, n, n, v1);
     dvec_add(v1, v2, S31_l, n);

     /* R13_l = E31 * B13 + D31 * T12 */
     dmat_mul(EXX, BXX, n, n, n, w1);
     dmat_mul(DXX, T12, n, n, n, w2);
     dmat_add(w1, w2, R13_l, n, n);

     /* T31_l = E31 * T32 + A31 * T32_l */
     dmat_mul(EXX, T32, n, n, n, w1);
     dmat_mul(AXX, T32_l, n, n, n, w2);
     dmat_add(w1, w2, T31_l, n, n);
}



static void add_all_l_u_up(double **PXX, int *i1, double **AXX, double **BXX, double *CXX, double **DXX, double **R23, double **T32, double *S32, double **R12_l, double **T12_l, double *S12_l, double **R21_l, double **T21_l, double *S21_l, double **R13_l, double **T31_l, double *S31_l, int n, double atran, double atran_l, double *v1, double *v2, double **w1, double **w2, work_data work) {
/*
     double **DXX;
*/
     double **EXX;
/*
     DXX = get_work1(&work, WORK_DXX);
*/
     EXX = get_work1(&work, WORK_DXX);

     /* D31 = A31 * R23 */
/*
     dmat_mul(AXX, R23, n, n, n, DXX);
*/
     /* E31 = (T21_l + D31 * R21_l) * P31 */
     dmat_mul(DXX, R21_l, n, n, n, EXX);
     dmat_add(T21_l, EXX, EXX, n, n);
     dmat_getrs(PXX, EXX, n, n, i1);

     /* S31_l = S21_l + E31 * C31 + A31 * (S32 * L(t12) + R23 * S12_l)  */
     dvec_scale(atran_l, S32, v1, n);
     dm_v_mul(R23, S12_l, n, n, v2);
     dvec_add(v1, v2, v1, n);

     dm_v_mul(AXX, v1, n, n, v2);

     dvec_add(S21_l, v2, v1, n);

     dm_v_mul(EXX, CXX, n, n, v2);
     dvec_add(v1, v2, S31_l, n);

     /* R13_l = R12_l + E31 * B13 + D31 * T12_l */
     dmat_mul(EXX, BXX, n, n, n, w1);
     dmat_mul(DXX, T12_l, n, n, n, w2);
     dmat_add(w1, w2, w1, n, n);
     dmat_add(R12_l, w1, R13_l, n, n);

     /* T31_l = E31 * T32 */
     dmat_mul(EXX, T32, n, n, n, T31_l);
}



static void add_all_p_l_up(double **PXX, int *i1, double **AXX, double **BXX, double *CXX, double **T12, double *S12, double **R21, double **R23, double **T32, double *S12_l, double *S21_l, double **R23_l, double **T32_l, double *S32_l, double **R13_l, double **T31_l, double *S31_l, int n, double atran, double *v1, double *v2, double **w1, double **w2, work_data work) {

     double **DXX;
     double **EXX;

     DXX = get_work1(&work, WORK_DXX);
     EXX = get_work1(&work, WORK_DXX);

     /* D31 = A31 * R23_l */
     dmat_mul(AXX, R23_l, n, n, n, DXX);

     /* E31 = D31 * R21 * P31 */
     dmat_mul(DXX, R21, n, n, n, EXX);
     dmat_getrs(PXX, EXX, n, n, i1);

     /* S31_l = S21_l + E31 * C31 + A31 * (S32_l * t12 + R23_l * S12 + R23 * S12_l) */
     dvec_scale(atran, S32_l, v1, n);
     dm_v_mul(R23_l, S12, n, n, v2);
     dvec_add(v1, v2, v1, n);
     dm_v_mul(R23, S12_l, n, n, v2);
     dvec_add(v1, v2, v1, n);

     dm_v_mul(AXX, v1, n, n, v2);

     dvec_add(S21_l, v2, S31_l, n);

     dm_v_mul(EXX, CXX, n, n, v1);
     dvec_add(S31_l, v1, S31_l, n);

     /* R13_l = E31 * B13 + D31 * T12 */
     dmat_mul(EXX, BXX, n, n, n, w1);
     dmat_mul(DXX, T12, n, n, n, w2);
     dmat_add(w1, w2, R13_l, n, n);

     /* T31_l = E31 * T32 + A31 * T32_l */
     dmat_mul(EXX, T32, n, n, n, w1);
     dmat_mul(AXX, T32_l, n, n, n, w2);
     dmat_add(w1, w2, T31_l, n, n);
}



static void add_all_l_p_up(double **PXX, int *i1, double **AXX, double **BXX, double *CXX, double **DXX, double **R23, double **T32, double *S32, double **R12_l, double **T12_l, double *S12_l, double **R21_l, double **T21_l, double *S21_l, double *S32_l, double **R13_l, double **T31_l, double *S31_l, int n, double atran, double atran_l, double *v1, double *v2, double **w1, double **w2, work_data work) {
/*
     double **DXX;
*/
     double **EXX;
/*
     DXX = get_work1(&work, WORK_DXX);
*/
     EXX = get_work1(&work, WORK_DXX);

     /* D31 = A31 * R23 */
/*
     dmat_mul(AXX, R23, n, n, n, DXX);
*/
     /* E31 = (T21_l + D31 * R21_l) * P31 */
     dmat_mul(DXX, R21_l, n, n, n, EXX);
     dmat_add(T21_l, EXX, EXX, n, n);
     dmat_getrs(PXX, EXX, n, n, i1);

     /* S31_l = S21_l + E31 * C31 + A31 * (S32_l * t12 + S32 * L(t12) + R23 * S12_l) */
     dvec_scale(atran, S32_l, v1, n);
     dvec_scale(atran_l, S32, v2, n);
     dvec_add(v1, v2, v1, n);
     dm_v_mul(R23, S12_l, n, n, v2);
     dvec_add(v1, v2, v1, n);

     dm_v_mul(AXX, v1, n, n, v2);

     dvec_add(S21_l, v2, S31_l, n);

     dm_v_mul(EXX, CXX, n, n, v1);
     dvec_add(S31_l, v1, S31_l, n);

     /* R13_l = R12_l + E31 * B13 + D31 * T12_l */
     dmat_mul(EXX, BXX, n, n, n, w1);
     dmat_mul(DXX, T12_l, n, n, n, w2);
     dmat_add(w1, w2, w1, n, n);

     dmat_add(R12_l, w1, R13_l, n, n);

     /* T31_l = E31 * T32 */
     dmat_mul(EXX, T32, n, n, n, T31_l);
}



static void add_all_l_l_up(double **PXX, int *i1, double **AXX, double **BXX, double *CXX, double **T12, double *S12, double **R21, double **R23, double **T32, double *S32, double **R12_l, double **T12_l, double *S12_l, double **R21_l, double **T21_l, double *S21_l, double **R23_l, double **T32_l, double *S32_l, double **R13_l, double **T31_l, double *S31_l, int n, double atran, double atran_l, double *v1, double *v2, double **w1, double **w2, work_data work) {

     double **DXX;

     DXX = get_work1(&work, WORK_DXX);

     /* D31 = (T21_l + A31 * (R23_l * R21 + R23 * R21_l)) * P31 */
     dmat_mul(R23_l, R21, n, n, n, w1);
     dmat_mul(R23, R21_l, n, n, n, w2);
     dmat_add(w1, w2, w1, n, n);
     dmat_mul(AXX, w1, n, n, n, w2);

     dmat_add(T21_l, w2, DXX, n, n);

     dmat_getrs(PXX, DXX, n, n, i1);

     /* S31_l = S21_l + D31 * C31 + A31 * (S32_l * t12 + S32 * L(t12) + R23_l * S12 + R23 * S12_l)  */
     dvec_scale(atran, S32_l, v1, n);
     dvec_scale(atran_l, S32, v2, n);
     dvec_add(v1, v2, v1, n);
     dm_v_mul(R23_l, S12, n, n, v2);
     dvec_add(v1, v2, v1, n);
     dm_v_mul(R23, S12_l, n, n, v2);
     dvec_add(v1, v2, v1, n);

     dm_v_mul(AXX, v1, n, n, v2);

     dvec_add(S21_l, v2, v1, n);

     dm_v_mul(DXX, CXX, n, n, v2);
     dvec_add(v1, v2, S31_l, n);

     /* R13_l = R12_l + D31 * B13 + A31 * (R23_l * T12 + R23 * T12_l) */
     dmat_mul(R23_l, T12, n, n, n, w1);
     dmat_mul(R23, T12_l, n, n, n, w2);
     dmat_add(w1, w2, w1, n, n);
     dmat_mul(AXX, w1, n, n, n, w2);
     dmat_add(R12_l, w2, w1, n, n);
     dmat_mul(DXX, BXX, n, n, n, w2);
     dmat_add(w1, w2, R13_l, n, n);

     /* T31_l = D31 * T32 + A31 * T32_l */
     dmat_mul(DXX, T32, n, n, n, w1);
     dmat_mul(AXX, T32_l, n, n, n, w2);
     dmat_add(w1, w2, T31_l, n, n);
}



/*******************************************************************************
 *
 ******************************************************************************/
void source_add_all(double **R12, double **T12, double *S12,
                    double **R21, double **T21, double *S21,
                    double **R23, double **T23, double *S23,
                    double **R32, double **T32, double *S32,
                    double **R13, double **T13, double *S13,
                    double **R31, double **T31, double *S31,
                    double ***R12_l, double ***T12_l, double **S12_l,
                    double ***R21_l, double ***T21_l, double **S21_l,
                    double ***R23_l, double ***T23_l, double **S23_l,
                    double ***R32_l, double ***T32_l, double **S32_l,
                    double ***R13_l, double ***T13_l, double **S13_l,
                    double ***R31_l, double ***T31_l, double **S31_l,
                    int n, int n_derivs, double atran, double *atran_l,
                    uchar *derivs, int flag, work_data work) {

     int i;

     int *i1;

     double *v1;
     double *v2;

     double *CXX;

     double **w1;

     double **PXX;
     double **AXX;

     i1  = get_work1(&work, WORK_IX);
if (flag)
     CXX = get_work1(&work, WORK_DX);

     v1  = get_work1(&work, WORK_DX);
     v2  = get_work1(&work, WORK_DX);

     w1  = get_work1(&work, WORK_DXX);

     PXX = get_work1(&work, WORK_DXX);
     AXX = get_work1(&work, WORK_DXX);

     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     /* w1 = LU of (E - R23 * R21)^T */
     dmat_mul(R23, R21, n, n, n, PXX);
     dmat_i_sub(PXX, w1, n);
     dmat_trans(w1, PXX, n, n);
     dmat_getrf(PXX, n, n, i1);

     /* A31 = T21 * P31 */
     dmat_copy(AXX, T21, n, n);
     dmat_getrs(PXX, AXX, n, n, i1);

     /* C31 = S32 * t12 + R23 * S12 */
if (flag) {
     dvec_scale(atran, S32, v1, n);
     dm_v_mul(R23, S12, n, n, v2);
     dvec_add(v1, v2, CXX, n);
}
     for (i = 0; i < n_derivs; ++i) {
          if (derivs[i] == ADDING_U_P)
               add_all_u_p_up(AXX, R23, S23_l[i], S32_l[i], S31_l[i], n, atran, v1, v2);
          else
          if (derivs[i] == ADDING_P_P)
               add_all_p_p_up(AXX, R23, S12_l[i], S21_l[i], S23_l[i], S32_l[i], S31_l[i], n, atran, v1, v2);
     }
if (flag) {
     /* S31 = S21       + A31 * C31 */
     dm_v_mul(AXX, CXX, n, n, v1);
     dvec_add(S21, v1, S31, n);
}

     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     /* w1 = LU of (E - R21 * R23)^T */
     dmat_mul(R21, R23, n, n, n, PXX);
     dmat_i_sub(PXX, w1, n);
     dmat_trans(w1, PXX, n, n);
     dmat_getrf(PXX, n, n, i1);

     /* A13 = T23 * P13 */
     dmat_copy(AXX, T23, n, n);
     dmat_getrs(PXX, AXX, n, n, i1);

     /* C13 = S12 + R21 * S32 * t12 */
if (flag) {
     dm_v_mul(R21, S32, n, n, v1);
     dvec_scale(atran, v1, v1, n);
     dvec_add(S12, v1, CXX, n);
}
     for (i = 0; i < n_derivs; ++i) {
          if (derivs[i] == ADDING_U_P)
               add_all_p_u_up(AXX, R21, S32_l[i], S23_l[i], S13_l[i], n, atran, v1, v2);
          else
          if (derivs[i] == ADDING_P_P)
               add_all_p_p_up(AXX, R21, S32_l[i], S23_l[i], S21_l[i], S12_l[i], S13_l[i], n, atran, v1, v2);
     }
if (flag) {
     /* S13 = S23 * t12 + A13 * C13 */
     dm_v_mul(AXX, CXX, n, n, v2);
     dvec_scale(atran, S23, v1, n);
     dvec_add(v1, v2, S13, n);
}
}



/*******************************************************************************
 *
 ******************************************************************************/
void layer_add_all(double **R12, double **T12, double *S12,
                   double **R21, double **T21, double *S21,
                   double **R23, double **T23, double *S23,
                   double **R32, double **T32, double *S32,
                   double **R13, double **T13, double *S13,
                   double **R31, double **T31, double *S31,
                   double ***R12_l, double ***T12_l, double **S12_l,
                   double ***R21_l, double ***T21_l, double **S21_l,
                   double ***R23_l, double ***T23_l, double **S23_l,
                   double ***R32_l, double ***T32_l, double **S32_l,
                   double ***R13_l, double ***T13_l, double **S13_l,
                   double ***R31_l, double ***T31_l, double **S31_l,
                   int n, int n_derivs, double atran, double *atran_l,
                   uchar *derivs, int flag, save_tree_data save_tree,
                   work_data work) {

     int i;

     int *i1;

     double *v1;
     double *v2;

     double *CXX;

     double **w1;
     double **w2;

     double **PXX;
     double **AXX;
     double **BXX;
     double **DXX;
#ifdef USE_AD_FOR_TL_LAYER_ADD_ALL
     double **R12_2;
     double **T12_2;
     double  *S12_2;
     double **R21_2;
     double **T21_2;
     double  *S21_2;

     double ***R12_l_2;
     double ***T12_l_2;
     double  **S12_l_2;
     double ***R21_l_2;
     double ***T21_l_2;
     double  **S21_l_2;
#endif
     forward_save_layer_add_all_data *save;


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (save_tree.t) {
          save_tree_encode_s(&save_tree, "layer_add_all");

          if (save_tree_retrieve_data(&save_tree, forward_save_layer_add_all_data, &save))
               forward_save_layer_add_all_alloc(save, n);
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     i1  = get_work1(&work, WORK_IX);

     CXX = get_work1(&work, WORK_DX);

     v1  = get_work1(&work, WORK_DX);
     v2  = get_work1(&work, WORK_DX);

     w1  = get_work1(&work, WORK_DXX);

     PXX = get_work1(&work, WORK_DXX);
     AXX = get_work1(&work, WORK_DXX);
     BXX = get_work1(&work, WORK_DXX);

     if (n_derivs > 0) {
          w2  = get_work1(&work, WORK_DXX);

          DXX = get_work1(&work, WORK_DXX);
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
#ifdef USE_AD_FOR_TL_LAYER_ADD_ALL
     R12_2 = get_work1(&work, WORK_DXX);
     T12_2 = get_work1(&work, WORK_DXX);
     S12_2 = get_work1(&work, WORK_DX);
     R21_2 = get_work1(&work, WORK_DXX);
     T21_2 = get_work1(&work, WORK_DXX);
     S21_2 = get_work1(&work, WORK_DX);

     R12_l_2 = get_work2(&work, WORK_DXX, WORK_DERIVS_V, derivs);
     T12_l_2 = get_work2(&work, WORK_DXX, WORK_DERIVS_V, derivs);
     S12_l_2 = get_work2(&work, WORK_DX,  WORK_DERIVS_V, derivs);
     R21_l_2 = get_work2(&work, WORK_DXX, WORK_DERIVS_V, derivs);
     T21_l_2 = get_work2(&work, WORK_DXX, WORK_DERIVS_V, derivs);
     S21_l_2 = get_work2(&work, WORK_DX,  WORK_DERIVS_V, derivs);

     dmat_copy(R12_2, R12, n, n);
     dmat_copy(T12_2, T12, n, n);
     dvec_copy(S12_2, S12, n);
     dmat_copy(R21_2, R21, n, n);
     dmat_copy(T21_2, T21, n, n);
     dvec_copy(S21_2, S21, n);

     for (i = 0; i < n_derivs; ++i) {
          if (! derivs[i])
               continue;

          dmat_copy(R12_l_2[i], R12_l[i], n, n);
          dmat_copy(T12_l_2[i], T12_l[i], n, n);
          dvec_copy(S12_l_2[i], S12_l[i], n);
          dmat_copy(R21_l_2[i], R21_l[i], n, n);
          dmat_copy(T21_l_2[i], T21_l[i], n, n);
          dvec_copy(S21_l_2[i], S21_l[i], n);
     }
#endif

     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     /* w1 = LU of (E - R23 * R21)^T */
     dmat_mul(R23, R21, n, n, n, PXX);
     dmat_i_sub(PXX, w1, n);
     dmat_trans(w1, PXX, n, n);
     dmat_getrf(PXX, n, n, i1);

     /* A31 = T21 * P31 */
     dmat_copy(AXX, T21, n, n);
     dmat_getrs(PXX, AXX, n, n, i1);

     /* B13 = R23 * T12 */
     dmat_mul(R23, T12, n, n, n, BXX);

     /* C31 = S32 * t12 + R23 * S12 */
     dvec_scale(atran, S32, v1, n);
     dm_v_mul(R23, S12, n, n, v2);
     dvec_add(v1, v2, CXX, n);

     if (n_derivs > 0)
          dmat_mul(AXX, R23, n, n, n, DXX);

     for (i = 0; i < n_derivs; ++i) {
          if (derivs[i] == ADDING_U_P)
               add_all_u_p_up(AXX, R23, S23_l[i], S32_l[i], S31_l[i], n, atran, v1, v2);
          else
          if (derivs[i] == ADDING_P_P)
               add_all_p_p_up(AXX, R23, S12_l[i], S21_l[i], S23_l[i], S32_l[i], S31_l[i], n, atran, v1, v2);
          else
          if (derivs[i] == ADDING_U_L)
               add_all_u_l_up(PXX, i1, AXX, BXX, CXX, T12, S12, R21, T32, R23_l[i], T32_l[i], S32_l[i], R13_l[i], T31_l[i], S31_l[i], n, atran, v1, v2, w1, w2, work);
          else
          if (derivs[i] == ADDING_L_P)
               add_all_l_p_up(PXX, i1, AXX, BXX, CXX, DXX, R23, T32, S32, R12_l[i], T12_l[i], S12_l[i], R21_l[i], T21_l[i], S21_l[i], S32_l[i], R13_l[i], T31_l[i], S31_l[i], n, atran, atran_l[i], v1, v2, w1, w2, work);
          else
          if (derivs[i] == ADDING_L_L)
               add_all_l_l_up(PXX, i1, AXX, BXX, CXX, T12, S12, R21, R23, T32, S32, R12_l[i], T12_l[i], S12_l[i], R21_l[i], T21_l[i], S21_l[i], R23_l[i], T32_l[i], S32_l[i], R13_l[i], T31_l[i], S31_l[i], n, atran, atran_l[i], v1, v2, w1, w2, work);
     }
if (flag) {
     /* S31 = S21       + A31 * C31 */
     dm_v_mul(AXX, CXX, n, n, v1);
     dvec_add(S21, v1, S31, n);

     /* R13 = R12 + A31 * B13 */
     dmat_mul(AXX, BXX, n, n, n, w1);
     dmat_add(R12, w1, R13, n, n);

     /* T31 = A31 * T32 */
     dmat_mul(AXX, T32, n, n, n, T31);
}

     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (save_tree.t) {
          memcpy(save->i31, i1, n * sizeof(int));
          dvec_copy(save->C31, CXX, n);
          dmat_copy(save->P31, PXX, n, n);
          dmat_copy(save->A31, AXX, n, n);
          dmat_copy(save->B13, BXX, n, n);
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     /* w1 = LU of (E - R21 * R23)^T */
     dmat_mul(R21, R23, n, n, n, PXX);
     dmat_i_sub(PXX, w1, n);
     dmat_trans(w1, PXX, n, n);
     dmat_getrf(PXX, n, n, i1);

     /* A13 = T23 * P13 */
     dmat_copy(AXX, T23, n, n);
     dmat_getrs(PXX, AXX, n, n, i1);

     /* B31 = R21 * T32 */
     dmat_mul(R21, T32, n, n, n, BXX);

     /* C13 = S12 + R21 * S32 * t12 */
     dm_v_mul(R21, S32, n, n, v1);
     dvec_scale(atran, v1, v1, n);
     dvec_add(S12, v1, CXX, n);

     if (n_derivs > 0)
          dmat_mul(AXX, R21, n, n, n, DXX);

     for (i = 0; i < n_derivs; ++i) {
          if (derivs[i] == ADDING_U_P)
               add_all_p_u_up(AXX, R21, S32_l[i], S23_l[i], S13_l[i], n, atran, v1, v2);
          else
          if (derivs[i] == ADDING_P_P)
               add_all_p_p_up(AXX, R21, S32_l[i], S23_l[i], S21_l[i], S12_l[i], S13_l[i], n, atran, v1, v2);
          else
          if (derivs[i] == ADDING_U_L)
               add_all_l_u_up(PXX, i1, AXX, BXX, CXX, DXX, R21, T12, S12, R32_l[i], T32_l[i], S32_l[i], R23_l[i], T23_l[i], S23_l[i], R31_l[i], T13_l[i], S13_l[i], n, atran, atran_l[i], v1, v2, w1, w2, work);
          else
          if (derivs[i] == ADDING_L_P)
               add_all_p_l_up(PXX, i1, AXX, BXX, CXX, T32, S32, R23, R21, T12, S32_l[i], S23_l[i], R21_l[i], T12_l[i], S12_l[i], R31_l[i], T13_l[i], S13_l[i], n, atran, v1, v2, w1, w2, work);
          else
          if (derivs[i] == ADDING_L_L)
               add_all_l_l_up(PXX, i1, AXX, BXX, CXX, T32, S32, R23, R21, T12, S12, R32_l[i], T32_l[i], S32_l[i], R23_l[i], T23_l[i], S23_l[i], R21_l[i], T12_l[i], S12_l[i], R31_l[i], T13_l[i], S13_l[i], n, atran, atran_l[i], v1, v2, w1, w2, work);
     }
if (flag) {
     /* S13 = S23 * t12 + A13 * C13 */
     dm_v_mul(AXX, CXX, n, n, v2);
     dvec_scale(atran, S23, v1, n);
     dvec_add(v1, v2, S13, n);

     /* R31 = R32 + A13 * B31 */
     dmat_mul(AXX, BXX, n, n, n, w1);
     dmat_add(R32, w1, R31, n, n);

     /* T13 = A13 * T12 */
/*
     dmat_mul(AXX, T12, n, n, n, T13);
*/
     dmat_mul(AXX, T12, n, n, n, w1);
     dmat_copy(T13, w1, n, n);
}

     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (save_tree.t) {
          memcpy(save->i13, i1, n * sizeof(int));
          dvec_copy(save->C13, CXX, n);
          dmat_copy(save->P13, PXX, n, n);
          dmat_copy(save->A13, AXX, n, n);
          dmat_copy(save->B31, BXX, n, n);
     }
#ifdef USE_AD_FOR_TL_LAYER_ADD_ALL
     layer_add_all_tl_with_ad(R12_2, T12_2, S12_2, R21_2, T21_2, S21_2, R23, T23, S23, R32, T32, S32, R13, T13, S13, R31, T31, S31, R12_l_2, T12_l_2, S12_l_2, R21_l_2, T21_l_2, S21_l_2, R23_l, T23_l, S23_l, R32_l, T32_l, S32_l, R13_l, T13_l, S13_l, R31_l, T31_l, S31_l, n, n_derivs, atran, atran_l, derivs, flag, save_tree, work);
#endif
}



/*******************************************************************************
 *
 ******************************************************************************/
void source_add_all2(double **R12, double **T12, double *S12,
                   double **R21, double **T21, double *S21,
                   double **R23, double **T23, double *S23,
                   double **R32, double **T32, double *S32,
                   double **R13, double **T13, double *S13,
                   double **R31, double **T31, double *S31,
                   double ***R12_l, double ***T12_l, double **S12_l,
                   double ***R21_l, double ***T21_l, double **S21_l,
                   double ***R23_l, double ***T23_l, double **S23_l,
                   double ***R32_l, double ***T32_l, double **S32_l,
                   double ***R13_l, double ***T13_l, double **S13_l,
                   double ***R31_l, double ***T31_l, double **S31_l,
                   int n, int n_derivs, double atran, double *atran_l,
                   uchar *derivs, int flag, int flag2, work_data work,
                   layer_add_aux_data *d) {

     int i;

     int *i1;

     double *v1;
     double *v2;

     double *CXX;

     double **w1;

     double **PXX;
     double **AXX;

     i1  = get_work1(&work, WORK_IX);

     CXX = get_work1(&work, WORK_DX);

     v1  = get_work1(&work, WORK_DX);
     v2  = get_work1(&work, WORK_DX);

     w1  = get_work1(&work, WORK_DXX);
if (! d) {
     PXX = get_work1(&work, WORK_DXX);
     AXX = get_work1(&work, WORK_DXX);
}

     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
if (d) {
     PXX = d->P31;
     AXX = d->A31;
}
if (! d || flag2) {
     /* w1 = LU of (E - R23 * R21)^T */
     dmat_mul(R23, R21, n, n, n, PXX);
     dmat_i_sub(PXX, w1, n);
     dmat_trans(w1, PXX, n, n);
     dmat_getrf(PXX, n, n, i1);

     /* A31 = T21 * P31 */
     dmat_copy(AXX, T21, n, n);
     dmat_getrs(PXX, AXX, n, n, i1);
}
     /* C31 = S32 * t12 + R23 * S12 */
     dvec_scale(atran, S32, v1, n);
     dm_v_mul(R23, S12, n, n, v2);
     dvec_add(v1, v2, CXX, n);

     for (i = 0; i < n_derivs; ++i) {
          if (derivs[i] == ADDING_U_P)
               add_all_u_p_up(AXX, R23, S23_l[i], S32_l[i], S31_l[i], n, atran, v1, v2);
          else
          if (derivs[i] == ADDING_P_P)
               add_all_p_p_up(AXX, R23, S12_l[i], S21_l[i], S23_l[i], S32_l[i], S31_l[i], n, atran, v1, v2);
     }
if (flag) {
     /* S31 = S21       + A31 * C31 */
     dm_v_mul(AXX, CXX, n, n, v1);
     dvec_add(S21, v1, S31, n);
}

     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
if (d) {
     PXX = d->P13;
     AXX = d->A13;
}
if (! d || flag2) {
     /* w1 = LU of (E - R21 * R23)^T */
     dmat_mul(R21, R23, n, n, n, PXX);
     dmat_i_sub(PXX, w1, n);
     dmat_trans(w1, PXX, n, n);
     dmat_getrf(PXX, n, n, i1);

     /* A13 = T23 * P13 */
     dmat_copy(AXX, T23, n, n);
     dmat_getrs(PXX, AXX, n, n, i1);
}
     /* C13 = S12 + R21 * S32 * t12 */
     dm_v_mul(R21, S32, n, n, v1);
     dvec_scale(atran, v1, v1, n);
     dvec_add(S12, v1, CXX, n);

     for (i = 0; i < n_derivs; ++i) {
          if (derivs[i] == ADDING_U_P)
               add_all_p_u_up(AXX, R21, S32_l[i], S23_l[i], S13_l[i], n, atran, v1, v2);
          else
          if (derivs[i] == ADDING_P_P)
               add_all_p_p_up(AXX, R21, S32_l[i], S23_l[i], S21_l[i], S12_l[i], S13_l[i], n, atran, v1, v2);
     }
if (flag) {
     /* S13 = S23 * t12 + A13 * C13 */
     dm_v_mul(AXX, CXX, n, n, v2);
     dvec_scale(atran, S23, v1, n);
     dvec_add(v1, v2, S13, n);
}
}



/*******************************************************************************
 *
 ******************************************************************************/
void layer_add_all2(double **R12, double **T12, double *S12,
                   double **R21, double **T21, double *S21,
                   double **R23, double **T23, double *S23,
                   double **R32, double **T32, double *S32,
                   double **R13, double **T13, double *S13,
                   double **R31, double **T31, double *S31,
                   double ***R12_l, double ***T12_l, double **S12_l,
                   double ***R21_l, double ***T21_l, double **S21_l,
                   double ***R23_l, double ***T23_l, double **S23_l,
                   double ***R32_l, double ***T32_l, double **S32_l,
                   double ***R13_l, double ***T13_l, double **S13_l,
                   double ***R31_l, double ***T31_l, double **S31_l,
                   int n, int n_derivs, double atran, double *atran_l,
                   uchar *derivs, int flag, int flag2, work_data work,
                  layer_add_aux_data *d) {

     int i;

     int *i1;

     double *v1;
     double *v2;

     double *CXX;

     double **w1;
     double **w2;

     double **PXX;
     double **AXX;
     double **BXX;

     i1  = get_work1(&work, WORK_IX);

     CXX = get_work1(&work, WORK_DX);

     v1  = get_work1(&work, WORK_DX);
     v2  = get_work1(&work, WORK_DX);

     w1  = get_work1(&work, WORK_DXX);

     if (n_derivs > 0)
          w2  = get_work1(&work, WORK_DXX);
if (! d) {
     PXX = get_work1(&work, WORK_DXX);
     AXX = get_work1(&work, WORK_DXX);
     BXX = get_work1(&work, WORK_DXX);
}

     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
if (d) {
     PXX = d->P31;
     AXX = d->A31;
     BXX = d->B13;
}
if (! d || flag2) {
     /* w1 = LU of (E - R23 * R21)^T */
     dmat_mul(R23, R21, n, n, n, PXX);
     dmat_i_sub(PXX, w1, n);
     dmat_trans(w1, PXX, n, n);
     dmat_getrf(PXX, n, n, i1);

     /* A31 = T21 * P31 */
     dmat_copy(AXX, T21, n, n);
     dmat_getrs(PXX, AXX, n, n, i1);

     /* B13 = R23 * T12 */
     dmat_mul(R23, T12, n, n, n, BXX);
}
     /* C31 = S32 * t12 + R23 * S12 */
     dvec_scale(atran, S32, v1, n);
     dm_v_mul(R23, S12, n, n, v2);
     dvec_add(v1, v2, CXX, n);

     for (i = 0; i < n_derivs; ++i) {
          if (derivs[i] == ADDING_U_P)
               add_all_u_p_up(AXX, R23, S23_l[i], S32_l[i], S31_l[i], n, atran, v1, v2);
          else
          if (derivs[i] == ADDING_P_P)
               add_all_p_p_up(AXX, R23, S12_l[i], S21_l[i], S23_l[i], S32_l[i], S31_l[i], n, atran, v1, v2);
          else
          if (derivs[i] == ADDING_U_L)
               add_all_u_l_up(PXX, i1, AXX, BXX, CXX, T12, S12, R21, T32, R23_l[i], T32_l[i], S32_l[i], R13_l[i], T31_l[i], S31_l[i], n, atran, v1, v2, w1, w2, work);
          else
          if (derivs[i] == ADDING_L_P)
               add_all_l_p_up(PXX, i1, AXX, BXX, CXX, NULL, R23, T32, S32, R12_l[i], T12_l[i], S12_l[i], R21_l[i], T21_l[i], S21_l[i], S32_l[i], R13_l[i], T31_l[i], S31_l[i], n, atran, atran_l[i], v1, v2, w1, w2, work);
          else
          if (derivs[i] == ADDING_L_L)
               add_all_l_l_up(PXX, i1, AXX, BXX, CXX, T12, S12, R21, R23, T32, S32, R12_l[i], T12_l[i], S12_l[i], R21_l[i], T21_l[i], S21_l[i], R23_l[i], T32_l[i], S32_l[i], R13_l[i], T31_l[i], S31_l[i], n, atran, atran_l[i], v1, v2, w1, w2, work);
     }
if (flag) {
     /* S31 = S21       + A31 * C31 */
     dm_v_mul(AXX, CXX, n, n, v1);
     dvec_add(S21, v1, S31, n);

     /* R13 = R12 + A31 * B13 */
     dmat_mul(AXX, BXX, n, n, n, w1);
     dmat_add(R12, w1, R13, n, n);

     /* T31 = A31 * T32 */
     dmat_mul(AXX, T32, n, n, n, T31);
}

     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
if (d) {
     PXX = d->P13;
     AXX = d->A13;
     BXX = d->B31;
}
if (! d || flag2) {
     /* w1 = LU of (E - R21 * R23)^T */
     dmat_mul(R21, R23, n, n, n, PXX);
     dmat_i_sub(PXX, w1, n);
     dmat_trans(w1, PXX, n, n);
     dmat_getrf(PXX, n, n, i1);

     /* A13 = T23 * P13 */
     dmat_copy(AXX, T23, n, n);
     dmat_getrs(PXX, AXX, n, n, i1);

     /* B31 = R21 * T32 */
     dmat_mul(R21, T32, n, n, n, BXX);
}
     /* C13 = S12 + R21 * S32 * t12 */
     dm_v_mul(R21, S32, n, n, v1);
     dvec_scale(atran, v1, v1, n);
     dvec_add(S12, v1, CXX, n);

     for (i = 0; i < n_derivs; ++i) {
          if (derivs[i] == ADDING_U_P)
               add_all_p_u_up(AXX, R21, S32_l[i], S23_l[i], S13_l[i], n, atran, v1, v2);
          else
          if (derivs[i] == ADDING_P_P)
               add_all_p_p_up(AXX, R21, S32_l[i], S23_l[i], S21_l[i], S12_l[i], S13_l[i], n, atran, v1, v2);
          else
          if (derivs[i] == ADDING_U_L)
               add_all_l_u_up(PXX, i1, AXX, BXX, CXX, NULL, R21, T12, S12, R32_l[i], T32_l[i], S32_l[i], R23_l[i], T23_l[i], S23_l[i], R31_l[i], T13_l[i], S13_l[i], n, atran, atran_l[i], v1, v2, w1, w2, work);
          else
          if (derivs[i] == ADDING_L_P)
               add_all_p_l_up(PXX, i1, AXX, BXX, CXX, T32, S32, R23, R21, T12, S32_l[i], S23_l[i], R21_l[i], T12_l[i], S12_l[i], R31_l[i], T13_l[i], S13_l[i], n, atran, v1, v2, w1, w2, work);
          else
          if (derivs[i] == ADDING_L_L)
               add_all_l_l_up(PXX, i1, AXX, BXX, CXX, T32, S32, R23, R21, T12, S12, R32_l[i], T32_l[i], S32_l[i], R23_l[i], T23_l[i], S23_l[i], R21_l[i], T12_l[i], S12_l[i], R31_l[i], T13_l[i], S13_l[i], n, atran, atran_l[i], v1, v2, w1, w2, work);
     }
if (flag) {
     /* S13 = S23 * t12 + A13 * C13 */
     dm_v_mul(AXX, CXX, n, n, v2);
     dvec_scale(atran, S23, v1, n);
     dvec_add(v1, v2, S13, n);

     /* R31 = R32 + A13 * B31 */
     dmat_mul(AXX, BXX, n, n, n, w1);
     dmat_add(R32, w1, R31, n, n);

     /* T13 = A13 * T12 */
/*
     dmat_mul(AXX, T12, n, n, n, T13);
*/
     dmat_mul(AXX, T12, n, n, n, w1);
     dmat_copy(T13, w1, n, n);
}
}



/*******************************************************************************
 *
 ******************************************************************************/
#ifdef INCLUDE_DEV_SOURCE
#include "../src2/xrtm_adding2.c"
#endif
