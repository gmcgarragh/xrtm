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
static void add_ref_u_p_a(double **AXX, double **R23, double *S32_a, double *S31_a, int n, double atran, double *v1, double *v2, double **w1) {

     dmat_trans(AXX, w1, n, n);
     dm_v_mul(w1, S31_a, n, n, v1);
     dvec_scale(atran, v1, v1, n);
     dvec_add(S32_a, v1, S32_a, n);

     dvec_zero(S31_a, n);
}



static void add_ref_p_u_a(double **AXX, double **R23, double *S12_a, double *S21_a, double *S31_a, int n, double atran, double *v1, double *v2) {

     dvec_zero(S31_a, n);
}



static void add_ref_p_p_a(double **AXX, double **R23, double *S12_a, double *S21_a, double *S32_a, double *S31_a, int n, double atran, double *v1, double *v2, double **w1) {

     dvec_add(S21_a, S31_a, S21_a, n);

     dmat_trans(AXX, w1, n, n);
     dm_v_mul(w1, S31_a, n, n, v1);

     dvec_scale(atran, v1, v2, n);
     dvec_add(S32_a, v2, S32_a, n);

     dmat_trans(R23, w1, n, n);
     dm_v_mul(w1, v1, n, n, v2);
     dvec_add(S12_a, v2, S12_a, n);

     dvec_zero(S31_a, n);
}



static void add_ref_u_l_a(double **PXX, int *i1, double **AXX, double **BXX, double *CXX, double **T12, double *S12, double **R21, double **R23, double **R23_a, double *S32_a, double **R13_a, double *S31_a, int n, double atran, double *v1, double *v2, double **w1, double **w2, work_data work) {

     double **w3;

     double **DXX_a;
     double **EXX_a;

     w3    = get_work1(&work, WORK_DXX);

     DXX_a = get_work1(&work, WORK_DXX);
     EXX_a = get_work1(&work, WORK_DXX);

     dmat_trans(BXX, w1, n, n);
     dmat_mul(R13_a, w1, n, n, n, EXX_a);

     dmat_trans(T12, w1, n, n);
     dmat_mul(R13_a, w1, n, n, n, DXX_a);

     dmat_vxvtmx(S31_a, CXX, 1., EXX_a, 1., n, n);

     dmat_trans(AXX, w1, n, n);
     dm_v_mul(w1, S31_a, n, n, v1);

     dvec_scale(atran, v1, v2, n);
     dvec_add(S32_a, v2, S32_a, n);

     dmat_vxvtmx(v1, S12, 1., R23_a, 1., n, n);

     dmat_copy(w1, EXX_a, n, n);
     dmat_getrs2('t', PXX, w1, n, n, i1);
     dmat_trans(R21, w2, n, n);
     dmat_mul(w1, w2, n, n, n, w3);
     dmat_add(DXX_a, w3, DXX_a, n, n);

     dmat_trans(AXX, w1, n, n);
     dmat_mul(w1, DXX_a, n, n, n, w2);
     dmat_add(R23_a, w2, R23_a, n, n);

     dmat_zero(R13_a, n, n);
     dvec_zero(S31_a, n);
}



static void add_ref_l_u_a(double **PXX, int *i1, double **AXX, double **BXX, double *CXX, double **R23, double *S32, double **R12_a, double **T12_a, double *S12_a, double **R21_a, double **T21_a, double *S21_a, double **R13_a, double *S31_a, int n, double atran, double *atran_a, double *v1, double *v2, double **w1, double **w2, work_data work) {

     dmat_zero(R13_a, n, n);
     dvec_zero(S31_a, n);
}



static void add_ref_p_l_a(double **PXX, int *i1, double **AXX, double **BXX, double *CXX, double **T12, double *S12, double **R21, double **R23, double *S12_a, double *S21_a, double **R23_a, double *S32_a, double **R13_a, double *S31_a, int n, double atran, double *v1, double *v2, double **w1, double **w2, work_data work) {

     double **w3;

     double **DXX_a;
     double **EXX_a;

     w3    = get_work1(&work, WORK_DXX);

     DXX_a = get_work1(&work, WORK_DXX);
     EXX_a = get_work1(&work, WORK_DXX);

     dmat_trans(BXX, w1, n, n);
     dmat_mul(R13_a, w1, n, n, n, EXX_a);

     dmat_trans(T12, w1, n, n);
     dmat_mul(R13_a, w1, n, n, n, DXX_a);

     dvec_add(S21_a, S31_a, S21_a, n);

     dmat_vxvtmx(S31_a, CXX, 1., EXX_a, 1., n, n);

     dmat_trans(AXX, w1, n, n);
     dm_v_mul(w1, S31_a, n, n, v1);

     dvec_scale(atran, v1, v2, n);
     dvec_add(S32_a, v2, S32_a, n);

     dmat_vxvtmx(v1, S12, 1., R23_a, 1., n, n);

     dmat_trans(R23, w1, n, n);
     dm_v_mul(w1, v1, n, n, v2);
     dvec_add(S12_a, v2, S12_a, n);

     dmat_copy(w1, EXX_a, n, n);
     dmat_getrs2('t', PXX, w1, n, n, i1);
     dmat_trans(R21, w2, n, n);
     dmat_mul(w1, w2, n, n, n, w3);
     dmat_add(DXX_a, w3, DXX_a, n, n);

     dmat_trans(AXX, w1, n, n);
     dmat_mul(w1, DXX_a, n, n, n, w2);
     dmat_add(R23_a, w2, R23_a, n, n);

     dmat_zero(R13_a, n, n);
     dvec_zero(S31_a, n);
}



static void add_ref_l_p_a(double **PXX, int *i1, double **AXX, double **BXX, double *CXX, double **DXX, double **R23, double *S32, double **R12_a, double **T12_a, double *S12_a, double **R21_a, double **T21_a, double *S21_a, double *S32_a, double **R13_a, double *S31_a, int n, double atran, double *atran_a, double *v1, double *v2, double **w1, double **w2, work_data work) {
/*
     double **w3;
*/
     double **EXX_a;
/*
     w3    = get_work1(&work, WORK_DXX);
*/
     EXX_a = get_work1(&work, WORK_DXX);

     dmat_add(R12_a, R13_a, R12_a, n, n);
/*
     dmat_trans(BXX, w1, n, n);
     dmat_mul(R13_a, w1, n, n, n, EXX_a);
*/
     dmat_gxgxmx(0, R13_a, 1, BXX, 1., EXX_a, 0., n, n, n);
/*
     dmat_trans(DXX, w1, n, n);
     dmat_mul(w1, R13_a, n, n, n, w2);
     dmat_add(T12_a, w2, T12_a, n, n);
*/
     dmat_gxgxmx(1, DXX, 0, R13_a, 1., T12_a, 1., n, n, n);

     dvec_add(S21_a, S31_a, S21_a, n);

     dmat_vxvtmx(S31_a, CXX, 1., EXX_a, 1., n, n);
/*
     dmat_trans(AXX, w1, n, n);
     dm_v_mul(w1, S31_a, n, n, v1);
*/
     dmat_gxvxmx(1, AXX, S31_a, 1., v1, 0., n, n);

     dvec_scale(atran, v1, v2, n);
     dvec_add(S32_a, v2, S32_a, n);

     *atran_a += dvec_dot(S32, v1, n);
/*
     dmat_trans(R23, w1, n, n);
     dm_v_mul(w1, v1, n, n, v2);
     dvec_add(S12_a, v2, S12_a, n);
*/
     dmat_gxvxmx(1, R23, v1, 1., S12_a, 1., n, n);

     dmat_copy(w1, EXX_a, n, n);
     dmat_getrs2('t', PXX, w1, n, n, i1);

     dmat_add(T21_a, w1, T21_a, n, n);
/*
     dmat_trans(DXX, w2, n, n);
     dmat_mul(w2, w1, n, n, n, w3);
     dmat_add(R21_a, w3, R21_a, n, n);
*/
     dmat_gxgxmx(1, DXX, 0, w1, 1., R21_a, 1., n, n, n);

     dmat_zero(R13_a, n, n);
     dvec_zero(S31_a, n);
}



static void add_ref_l_l_a(double **PXX, int *i1, double **AXX, double **BXX, double *CXX, double **T12, double *S12, double **R21, double **R23, double *S32, double **R12_a, double **T12_a, double *S12_a, double **R21_a, double **T21_a, double *S21_a, double **R23_a, double *S32_a, double **R13_a, double *S31_a, int n, double atran, double *atran_a, double *v1, double *v2, double **w1, double **w2, work_data work) {
/*
     double **w3;
*/
     double **DXX_a;
/*
     w3    = get_work1(&work, WORK_DXX);
*/
     DXX_a = get_work1(&work, WORK_DXX);

     dmat_add(R12_a, R13_a, R12_a, n, n);
/*
     dmat_trans(BXX, w1, n, n);
     dmat_mul(R13_a, w1, n, n, n, DXX_a);
*/
     dmat_gxgxmx(0, R13_a, 1, BXX, 1., DXX_a, 0., n, n, n);
/*
     dmat_trans(AXX, w1, n, n);
     dmat_mul(w1, R13_a, n, n, n, w2);
*/
     dmat_gxgxmx(1, AXX, 0, R13_a, 1., w2, 0., n, n, n);
/*
     dmat_trans(T12, w1, n, n);
     dmat_mul(w2, w1, n, n, n, w3);
     dmat_add(R23_a, w3, R23_a, n, n);
*/
     dmat_gxgxmx(0, w2, 1, T12, 1., R23_a, 1., n, n, n);
/*
     dmat_trans(R23, w1, n, n);
     dmat_mul(w1, w2, n, n, n, w3);
     dmat_add(T12_a, w3, T12_a, n, n);
*/
     dmat_gxgxmx(1, R23, 0, w2, 1., T12_a, 1., n, n, n);

     dvec_add(S21_a, S31_a, S21_a, n);

     dmat_vxvtmx(S31_a, CXX, 1., DXX_a, 1., n, n);
/*
     dmat_trans(AXX, w1, n, n);
     dm_v_mul(w1, S31_a, n, n, v1);
*/
     dmat_gxvxmx(1, AXX, S31_a, 1., v1, 0., n, n);

     dvec_scale(atran, v1, v2, n);
     dvec_add(S32_a, v2, S32_a, n);

     *atran_a += dvec_dot(S32, v1, n);     

     dmat_vxvtmx(v1, S12, 1., R23_a, 1., n, n);
/*
     dmat_trans(R23, w1, n, n);
     dm_v_mul(w1, v1, n, n, v2);
     dvec_add(S12_a, v2, S12_a, n);
*/
     dmat_gxvxmx(1, R23, v1, 1., S12_a, 1., n, n);

     dmat_copy(w1, DXX_a, n, n);
     dmat_getrs2('t', PXX, w1, n, n, i1);

     dmat_add(T21_a, w1, T21_a, n, n);
/*
     dmat_trans(AXX, w2, n, n);
     dmat_mul(w2, w1, n, n, n, w3);
*/
     dmat_gxgxmx(1, AXX, 0, w1, 1., w2, 0., n, n, n);
/*
     dmat_trans(R21, w1, n, n);
     dmat_mul(w3, w1, n, n, n, w2);
     dmat_add(R23_a, w2, R23_a, n, n);
*/
     dmat_gxgxmx(0, w2, 1, R21, 1., R23_a, 1., n, n, n);
/*
     dmat_trans(R23, w1, n, n);
     dmat_mul(w1, w3, n, n, n, w2);
     dmat_add(R21_a, w2, R21_a, n, n);
*/
     dmat_gxgxmx(1, R23, 0, w2, 1., R21_a, 1., n, n, n);

     dmat_zero(R13_a, n, n);
     dvec_zero(S31_a, n);
}



/*******************************************************************************
 *
 ******************************************************************************/
static void forward_save_layer_add_ref_free(forward_save_layer_add_ref_data *d);

int forward_save_layer_add_ref_alloc(forward_save_layer_add_ref_data *d, int n) {

     d->free = (void (*)(void *)) forward_save_layer_add_ref_free;

     d->i13 = alloc_array1_i(n);
     d->C13 = alloc_array1_d(n);
     d->P13 = alloc_array2_d(n, n);
     d->A13 = alloc_array2_d(n, n);
     d->B31 = alloc_array2_d(n, n);

     return 0;
}



static void forward_save_layer_add_ref_free(forward_save_layer_add_ref_data *d) {

     free_array1_i(d->i13);
     free_array1_d(d->C13);
     free_array2_d(d->P13);
     free_array2_d(d->A13);
     free_array2_d(d->B31);
}



/*******************************************************************************
 *
 ******************************************************************************/
void layer_add_ref_a(double **R12, double **T12, double *S12,
                     double **R21, double **T21, double *S21,
                     double **R23, double *S32,
                     double **R13, double *S31,
                     double **R12_a, double **T12_a, double *S12_a,
                     double **R21_a, double **T21_a, double *S21_a,
                     double **R23_a, double *S32_a,
                     double **R13_a, double *S31_a,
                     int n, double atran, double *atran_a,
                     uchar derivs, int flag, int flag2,
                     save_tree_data save_tree, work_data work) {

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

     forward_save_layer_add_ref_data *save;


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     save_tree_encode_s(&save_tree, "layer_add_ref");

     save_tree_retrieve_data(&save_tree, forward_save_layer_add_ref_data, &save);


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     v1  = get_work1(&work, WORK_DX);
     v2  = get_work1(&work, WORK_DX);

     w1  = get_work1(&work, WORK_DXX);
     w2  = get_work1(&work, WORK_DXX);

     DXX = get_work1(&work, WORK_DXX);


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     i1  = save->i13;

     CXX = save->C13;

     PXX = save->P13;
     AXX = save->A13;
     BXX = save->B31;

     dmat_mul(AXX, R23, n, n, n, DXX);


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if ((! flag && derivs == ADDING_U_P) || (flag && derivs == ADDING_P_U))
          add_ref_u_p_a(AXX, R23, S32_a, S31_a, n, atran, v1, v2, w1);
     else
     if ((! flag && derivs == ADDING_P_U) || (flag && derivs == ADDING_U_P))
          add_ref_p_u_a(AXX, R23, S12_a, S21_a, S31_a, n, atran, v1, v2);
     else
     if (derivs == ADDING_P_P)
          add_ref_p_p_a(AXX, R23, S12_a, S21_a, S32_a, S31_a, n, atran, v1, v2, w1);
     else
     if ((! flag && derivs == ADDING_U_L) || (flag && derivs == ADDING_L_U))
          add_ref_u_l_a(PXX, i1, AXX, BXX, CXX, T12, S12, R21, R23, R23_a, S32_a, R13_a, S31_a, n, atran, v1, v2, w1, w2, work);
     else
     if ((! flag && derivs == ADDING_L_U) || (flag && derivs == ADDING_U_L))
          add_ref_l_u_a(PXX, i1, AXX, BXX, CXX, R23, S32, R12_a, T12_a, S12_a, R21_a, T21_a, S21_a, R13_a, S31_a, n, atran, atran_a, v1, v2, w1, w2, work);
     else
     if ((! flag && derivs == ADDING_P_L) || (flag && derivs == ADDING_L_P))
          add_ref_p_l_a(PXX, i1, AXX, BXX, CXX, T12, S12, R21, R23, S12_a, S21_a, R23_a, S32_a, R13_a, S31_a, n, atran, v1, v2, w1, w2, work);
     else
     if ((! flag && derivs == ADDING_L_P) || (flag && derivs == ADDING_P_L))
          add_ref_l_p_a(PXX, i1, AXX, BXX, CXX, DXX, R23, S32, R12_a, T12_a, S12_a, R21_a, T21_a, S21_a, S32_a, R13_a, S31_a, n, atran, atran_a, v1, v2, w1, w2, work);
     else
     if (derivs == ADDING_L_L)
          add_ref_l_l_a(PXX, i1, AXX, BXX, CXX, T12, S12, R21, R23, S32, R12_a, T12_a, S12_a, R21_a, T21_a, S21_a, R23_a, S32_a, R13_a, S31_a, n, atran, atran_a, v1, v2, w1, w2, work);
}



/*******************************************************************************
 * 
 ******************************************************************************/
/*
static void add_all_u_p_up_a(double **AXX, double **R23, double *S23_a, double *S32_a, double *S31_a, int n, double atran, double *v1, double *v2, double **w1) {

     dmat_trans(AXX, w1, n, n);
     dm_v_mul(w1, S31_a, n, n, v1);
     dvec_scale(atran, v1, v1, n);
     dvec_add(S32_a, v1, S32_a, n);

     dvec_zero(S31_a, n);
}
*/

/*
static void add_all_p_u_up_a(double **AXX, double **R23, double *S12_a, double *S21_a, double *S31_a, int n, double atran, double *v1, double *v2) {

     dvec_add(S21_a, S31_a, S21_a, n);

     dvec_zero(S31_a, n);
}
*/

/*
static void add_all_p_p_up_a(double **AXX, double **R23, double *S12_a, double *S21_a, double *S23_a, double *S32_a, double *S31_a, int n, double atran, double *v1, double *v2, double **w1) {

     dvec_add(S21_a, S31_a, S21_a, n);

     dmat_trans(AXX, w1, n, n);
     dm_v_mul(w1, S31_a, n, n, v1);

     dvec_scale(atran, v1, v2, n);
     dvec_add(S32_a, v2, S32_a, n);

     dmat_trans(R23, w1, n, n);
     dm_v_mul(w1, v1, n, n, v2);
     dvec_add(S12_a, v2, S12_a, n);

     dvec_zero(S31_a, n);
}
*/


static void add_all_u_l_up_a(double **PXX, int *i1, double **AXX, double **BXX, double *CXX, double **T12, double *S12, double **R21, double **T32, double **R23_a, double **T32_a, double *S32_a, double **R13_a, double **T31_a, double *S31_a, int n, double atran, double *v1, double *v2, double **w1, double **w2, work_data work) {
/*
     double **w3;
*/
     double **DXX_a;
     double **EXX_a;
/*
     w3    = get_work1(&work, WORK_DXX);
*/
     DXX_a = get_work1(&work, WORK_DXX);
     EXX_a = get_work1(&work, WORK_DXX);
/*
     dmat_trans(T32, w1, n, n);
     dmat_mul(T31_a, w1, n, n, n, EXX_a);
*/
     dmat_gxgxmx(0, T31_a, 1, T32, 1., EXX_a, 0., n, n, n);
/*
     dmat_trans(AXX, w1, n, n);
     dmat_mul(w1, T31_a, n, n, n, w2);
     dmat_add(T32_a, w2, T32_a, n, n);
*/
     dmat_gxgxmx(1, AXX, 0, T31_a, 1., T32_a, 1., n, n, n);
/*
     dmat_trans(BXX, w1, n, n);
     dmat_mul(R13_a, w1, n, n, n, w2);
     dmat_add(EXX_a, w2, EXX_a, n, n);
*/
     dmat_gxgxmx(0, R13_a, 1, BXX, 1., EXX_a, 1., n, n, n);
/*
     dmat_trans(T12, w1, n, n);
     dmat_mul(R13_a, w1, n, n, n, DXX_a);
*/
     dmat_gxgxmx(0, R13_a, 1, T12, 1., DXX_a, 0., n, n, n);

     dmat_vxvtmx(S31_a, CXX, 1., EXX_a, 1., n, n);
/*
     dmat_trans(AXX, w1, n, n);
     dm_v_mul(w1, S31_a, n, n, v1);
*/
     dmat_gxvxmx(1, AXX, S31_a, 1., v1, 0., n, n);

     dvec_scale(atran, v1, v2, n);
     dvec_add(S32_a, v2, S32_a, n);

     dmat_vxvtmx(v1, S12, 1., R23_a, 1., n, n);

     dmat_copy(w1, EXX_a, n, n);
     dmat_getrs2('t', PXX, w1, n, n, i1);
/*
     dmat_trans(R21, w2, n, n);
     dmat_mul(w1, w2, n, n, n, w3);
     dmat_add(DXX_a, w3, DXX_a, n, n);
*/
     dmat_gxgxmx(0, w1, 1, R21, 1., DXX_a, 1., n, n, n);
/*
     dmat_trans(AXX, w1, n, n);
     dmat_mul(w1, DXX_a, n, n, n, w2);
     dmat_add(R23_a, w2, R23_a, n, n);
*/
     dmat_gxgxmx(1, AXX, 0, DXX_a, 1., R23_a, 1., n, n, n);

     dmat_zero(R13_a, n, n);
     dmat_zero(T31_a, n, n);
     dvec_zero(S31_a, n);
}



static void add_all_l_u_up_a(double **PXX, int *i1, double **AXX, double **BXX, double *CXX, double **DXX, double **R23, double **T32, double *S32, double **R12_a, double **T12_a, double *S12_a, double **R21_a, double **T21_a, double *S21_a, double **R13_a, double **T31_a, double *S31_a, int n, double atran, double *atran_a, double *v1, double *v2, double **w1, double **w2, work_data work) {
/*
     double **w3;
*/
     double **EXX_a;
/*
     w3    = get_work1(&work, WORK_DXX);
*/
     EXX_a = get_work1(&work, WORK_DXX);
/*
     dmat_trans(T32, w1, n, n);
     dmat_mul(T31_a, w1, n, n, n, EXX_a);
*/
     dmat_gxgxmx(0, T31_a, 1, T32, 1., EXX_a, 0., n, n, n);

     dmat_add(R12_a, R13_a, R12_a, n, n);
/*
     dmat_trans(BXX, w1, n, n);
     dmat_mul(R13_a, w1, n, n, n, w2);
     dmat_add(EXX_a, w2, EXX_a, n, n);
*/
     dmat_gxgxmx(0, R13_a, 1, BXX, 1., EXX_a, 1., n, n, n);
/*
     dmat_trans(DXX, w1, n, n);
     dmat_mul(w1, R13_a, n, n, n, w2);
     dmat_add(T12_a, w2, T12_a, n, n);
*/
     dmat_gxgxmx(1, DXX, 0, R13_a, 1., T12_a, 1., n, n, n);

     dvec_add(S21_a, S31_a, S21_a, n);

     dmat_vxvtmx(S31_a, CXX, 1., EXX_a, 1., n, n);
/*
     dmat_trans(AXX, w1, n, n);
     dm_v_mul(w1, S31_a, n, n, v1);
*/
     dmat_gxvxmx(1, AXX, S31_a, 1., v1, 0., n, n);

     *atran_a += dvec_dot(S32, v1, n);     
/*
     dmat_trans(R23, w1, n, n);
     dm_v_mul(w1, v1, n, n, v2);
     dvec_add(S12_a, v2, S12_a, n);
*/
     dmat_gxvxmx(1, R23, v1, 1., S12_a, 1., n, n);

     dmat_copy(w1, EXX_a, n, n);
     dmat_getrs2('t', PXX, w1, n, n, i1);

     dmat_add(T21_a, w1, T21_a, n, n);
/*
     dmat_trans(DXX, w2, n, n);
     dmat_mul(w2, w1, n, n, n, w3);
     dmat_add(R21_a, w3, R21_a, n, n);
*/
     dmat_gxgxmx(1, DXX, 0, w1, 1., R21_a, 1., n, n, n);

     dmat_zero(R13_a, n, n);
     dmat_zero(T31_a, n, n);
     dvec_zero(S31_a, n);
}



static void add_all_p_l_up_a(double **PXX, int *i1, double **AXX, double **BXX, double *CXX, double **T12, double *S12, double **R21, double **R23, double **T32, double *S12_a, double *S21_a, double **R23_a, double **T32_a, double *S32_a, double **R13_a, double **T31_a, double *S31_a, int n, double atran, double *v1, double *v2, double **w1, double **w2, work_data work) {
/*
     double **w3;
*/
     double **DXX_a;
     double **EXX_a;
/*
     w3    = get_work1(&work, WORK_DXX);
*/
     DXX_a = get_work1(&work, WORK_DXX);
     EXX_a = get_work1(&work, WORK_DXX);
/*
     dmat_trans(T32, w1, n, n);
     dmat_mul(T31_a, w1, n, n, n, EXX_a);
*/
     dmat_gxgxmx(0, T31_a, 1, T32, 1., EXX_a, 0., n, n, n);
/*
     dmat_trans(AXX, w1, n, n);
     dmat_mul(w1, T31_a, n, n, n, w2);
     dmat_add(T32_a, w2, T32_a, n, n);
*/
     dmat_gxgxmx(1, AXX, 0, T31_a, 1., T32_a, 1., n, n, n);
/*
     dmat_trans(BXX, w1, n, n);
     dmat_mul(R13_a, w1, n, n, n, w2);
     dmat_add(EXX_a, w2, EXX_a, n, n);
*/
     dmat_gxgxmx(0, R13_a, 1, BXX, 1., EXX_a, 1., n, n, n);
/*
     dmat_trans(T12, w1, n, n);
     dmat_mul(R13_a, w1, n, n, n, DXX_a);
*/
     dmat_gxgxmx(0, R13_a, 1, T12, 1., DXX_a, 0., n, n, n);

     dvec_add(S21_a, S31_a, S21_a, n);

     dmat_vxvtmx(S31_a, CXX, 1., EXX_a, 1., n, n);
/*
     dmat_trans(AXX, w1, n, n);
     dm_v_mul(w1, S31_a, n, n, v1);
*/
     dmat_gxvxmx(1, AXX, S31_a, 1., v1, 0., n, n);

     dvec_scale(atran, v1, v2, n);
     dvec_add(S32_a, v2, S32_a, n);

     dmat_vxvtmx(v1, S12, 1., R23_a, 1., n, n);
/*
     dmat_trans(R23, w1, n, n);
     dm_v_mul(w1, v1, n, n, v2);
     dvec_add(S12_a, v2, S12_a, n);
*/
     dmat_gxvxmx(1, R23, v1, 1., S12_a, 1., n, n);

     dmat_copy(w1, EXX_a, n, n);
     dmat_getrs2('t', PXX, w1, n, n, i1);
/*
     dmat_trans(R21, w2, n, n);
     dmat_mul(w1, w2, n, n, n, w3);
     dmat_add(DXX_a, w3, DXX_a, n, n);
*/
     dmat_gxgxmx(0, w1, 1, R21, 1., DXX_a, 1., n, n, n);
/*
     dmat_trans(AXX, w1, n, n);
     dmat_mul(w1, DXX_a, n, n, n, w2);
     dmat_add(R23_a, w2, R23_a, n, n);
*/
     dmat_gxgxmx(1, AXX, 0, DXX_a, 1., R23_a, 1., n, n, n);

     dmat_zero(R13_a, n, n);
     dmat_zero(T31_a, n, n);
     dvec_zero(S31_a, n);
}



static void add_all_l_p_up_a(double **PXX, int *i1, double **AXX, double **BXX, double *CXX, double **DXX, double **R23, double **T32, double *S32, double **R12_a, double **T12_a, double *S12_a, double **R21_a, double **T21_a, double *S21_a, double *S32_a, double **R13_a, double **T31_a, double *S31_a, int n, double atran, double *atran_a, double *v1, double *v2, double **w1, double **w2, work_data work) {
/*
     double **w3;
*/
     double **EXX_a;
/*
     w3    = get_work1(&work, WORK_DXX);
*/
     EXX_a = get_work1(&work, WORK_DXX);
/*
     dmat_trans(T32, w1, n, n);
     dmat_mul(T31_a, w1, n, n, n, EXX_a);
*/
     dmat_gxgxmx(0, T31_a, 1, T32, 1., EXX_a, 0., n, n, n);

     dmat_add(R12_a, R13_a, R12_a, n, n);
/*
     dmat_trans(BXX, w1, n, n);
     dmat_mul(R13_a, w1, n, n, n, w2);
     dmat_add(EXX_a, w2, EXX_a, n, n);
*/
     dmat_gxgxmx(0, R13_a, 1, BXX, 1., EXX_a, 1., n, n, n);
/*
     dmat_trans(DXX, w1, n, n);
     dmat_mul(w1, R13_a, n, n, n, w2);
     dmat_add(T12_a, w2, T12_a, n, n);
*/
     dmat_gxgxmx(1, DXX, 0, R13_a, 1., T12_a, 1., n, n, n);

     dvec_add(S21_a, S31_a, S21_a, n);

     dmat_vxvtmx(S31_a, CXX, 1., EXX_a, 1., n, n);
/*
     dmat_trans(AXX, w1, n, n);
     dm_v_mul(w1, S31_a, n, n, v1);
*/
     dmat_gxvxmx(1, AXX, S31_a, 1., v1, 0., n, n);

     dvec_scale(atran, v1, v2, n);
     dvec_add(S32_a, v2, S32_a, n);

     *atran_a += dvec_dot(S32, v1, n);
/*
     dmat_trans(R23, w1, n, n);
     dm_v_mul(w1, v1, n, n, v2);
     dvec_add(S12_a, v2, S12_a, n);
*/
     dmat_gxvxmx(1, R23, v1, 1., S12_a, 1., n, n);

     dmat_copy(w1, EXX_a, n, n);
     dmat_getrs2('t', PXX, w1, n, n, i1);

     dmat_add(T21_a, w1, T21_a, n, n);
/*
     dmat_trans(DXX, w2, n, n);
     dmat_mul(w2, w1, n, n, n, w3);
     dmat_add(R21_a, w3, R21_a, n, n);
*/
     dmat_gxgxmx(1, DXX, 0, w1, 1., R21_a, 1., n, n, n);

     dmat_zero(R13_a, n, n);
     dmat_zero(T31_a, n, n);
     dvec_zero(S31_a, n);
}



static void add_all_l_l_up_a(double **PXX, int *i1, double **AXX, double **BXX, double *CXX, double **T12, double *S12, double **R21, double **R23, double **T32, double *S32, double **R12_a, double **T12_a, double *S12_a, double **R21_a, double **T21_a, double *S21_a, double **R23_a, double **T32_a, double *S32_a, double **R13_a, double **T31_a, double *S31_a, int n, double atran, double *atran_a, double *v1, double *v2, double **w1, double **w2, work_data work) {
/*
     double **w3;
*/
     double **DXX_a;
/*
     w3    = get_work1(&work, WORK_DXX);
*/
     DXX_a = get_work1(&work, WORK_DXX);
/*
     dmat_trans(T32, w1, n, n);
     dmat_mul(T31_a, w1, n, n, n, DXX_a);
*/
     dmat_gxgxmx(0, T31_a, 1, T32, 1., DXX_a, 0., n, n, n);
/*
     dmat_trans(AXX, w1, n, n);
     dmat_mul(w1, T31_a, n, n, n, w2);
     dmat_add(T32_a, w2, T32_a, n, n);
*/
     dmat_gxgxmx(1, AXX, 0, T31_a, 1., T32_a, 1., n, n, n);

     dmat_add(R12_a, R13_a, R12_a, n, n);
/*
     dmat_trans(BXX, w1, n, n);
     dmat_mul(R13_a, w1, n, n, n, w2);
     dmat_add(DXX_a, w2, DXX_a, n, n);
*/
     dmat_gxgxmx(0, R13_a, 1, BXX, 1., DXX_a, 1., n, n, n);
/*
     dmat_trans(AXX, w1, n, n);
     dmat_mul(w1, R13_a, n, n, n, w2);
*/
     dmat_gxgxmx(1, AXX, 0, R13_a, 1., w2, 0., n, n, n);
/*
     dmat_trans(T12, w1, n, n);
     dmat_mul(w2, w1, n, n, n, w3);
     dmat_add(R23_a, w3, R23_a, n, n);
*/
     dmat_gxgxmx(0, w2, 1, T12, 1., R23_a, 1., n, n, n);
/*
     dmat_trans(R23, w1, n, n);
     dmat_mul(w1, w2, n, n, n, w3);
     dmat_add(T12_a, w3, T12_a, n, n);
*/
     dmat_gxgxmx(1, R23, 0, w2, 1., T12_a, 1., n, n, n);

     dvec_add(S21_a, S31_a, S21_a, n);

     dmat_vxvtmx(S31_a, CXX, 1., DXX_a, 1., n, n);
/*
     dmat_trans(AXX, w1, n, n);
     dm_v_mul(w1, S31_a, n, n, v1);
*/
     dmat_gxvxmx(1, AXX, S31_a, 1., v1, 0., n, n);

     dvec_scale(atran, v1, v2, n);
     dvec_add(S32_a, v2, S32_a, n);

     dmat_gxvxmx(1, AXX, S31_a, 1., v1, 0., n, n);

     *atran_a += dvec_dot(S32, v1, n);     

     dmat_vxvtmx(v1, S12, 1., R23_a, 1., n, n);
/*
     dmat_trans(R23, w1, n, n);
     dm_v_mul(w1, v1, n, n, v2);
     dvec_add(S12_a, v2, S12_a, n);
*/
     dmat_gxvxmx(1, R23, v1, 1., S12_a, 1., n, n);

     dmat_copy(w1, DXX_a, n, n);
     dmat_getrs2('t', PXX, w1, n, n, i1);

     dmat_add(T21_a, w1, T21_a, n, n);
/*
     dmat_trans(AXX, w2, n, n);
     dmat_mul(w2, w1, n, n, n, w3);
*/
     dmat_gxgxmx(1, AXX, 0, w1, 1., w2, 0., n, n, n);
/*
     dmat_trans(R21, w1, n, n);
     dmat_mul(w3, w1, n, n, n, w2);
     dmat_add(R23_a, w2, R23_a, n, n);
*/
     dmat_gxgxmx(0, w2, 1, R21, 1., R23_a, 1., n, n, n);
/*
     dmat_trans(R23, w1, n, n);
     dmat_mul(w1, w3, n, n, n, w2);
     dmat_add(R21_a, w2, R21_a, n, n);
*/
     dmat_gxgxmx(1, R23, 0, w2, 1., R21_a, 1., n, n, n);

     dmat_zero(R13_a, n, n);
     dmat_zero(T31_a, n, n);
     dvec_zero(S31_a, n);
}



/*******************************************************************************
 *
 ******************************************************************************/
static void forward_save_layer_add_all_free(forward_save_layer_add_all_data *d);

int forward_save_layer_add_all_alloc(forward_save_layer_add_all_data *d, int n) {

     d->free = (void (*)(void *)) forward_save_layer_add_all_free;

     d->i31 = alloc_array1_i(n);
     d->C31 = alloc_array1_d(n);
     d->P31 = alloc_array2_d(n, n);
     d->A31 = alloc_array2_d(n, n);
     d->B13 = alloc_array2_d(n, n);

     d->i13 = alloc_array1_i(n);
     d->C13 = alloc_array1_d(n);
     d->P13 = alloc_array2_d(n, n);
     d->A13 = alloc_array2_d(n, n);
     d->B31 = alloc_array2_d(n, n);

     return 0;
}



static void forward_save_layer_add_all_free(forward_save_layer_add_all_data *d) {

     free_array1_i(d->i31);
     free_array1_d(d->C31);
     free_array2_d(d->P31);
     free_array2_d(d->A31);
     free_array2_d(d->B13);

     free_array1_i(d->i13);
     free_array1_d(d->C13);
     free_array2_d(d->P13);
     free_array2_d(d->A13);
     free_array2_d(d->B31);
}



/*******************************************************************************
 *
 ******************************************************************************/
void layer_add_all_a(double **R12, double **T12, double *S12,
                     double **R21, double **T21, double *S21,
                     double **R23, double **T23, double *S23,
                     double **R32, double **T32, double *S32,
                     double **R13, double **T13, double *S13,
                     double **R31, double **T31, double *S31,
                     double **R12_a, double **T12_a, double *S12_a,
                     double **R21_a, double **T21_a, double *S21_a,
                     double **R23_a, double **T23_a, double *S23_a,
                     double **R32_a, double **T32_a, double *S32_a,
                     double **R13_a, double **T13_a, double *S13_a,
                     double **R31_a, double **T31_a, double *S31_a,
                     int n, double atran, double *atran_a, uchar derivs,
                     save_tree_data save_tree, work_data work) {

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

     forward_save_layer_add_all_data *save;


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     save_tree_encode_s(&save_tree, "layer_add_all");

     save_tree_retrieve_data(&save_tree, forward_save_layer_add_all_data, &save);


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     v1  = get_work1(&work, WORK_DX);
     v2  = get_work1(&work, WORK_DX);

     w1  = get_work1(&work, WORK_DXX);
     w2  = get_work1(&work, WORK_DXX);

     DXX = get_work1(&work, WORK_DXX);


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     i1  = save->i31;

     CXX = save->C31;

     PXX = save->P31;
     AXX = save->A31;
     BXX = save->B13;

     dmat_mul(AXX, R23, n, n, n, DXX);


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
/*
     if (derivs == ADDING_U_P)
          add_all_u_p_up_a(AXX, R23, S23_a, S32_a, S31_a, n, atran, v1, v2, w1);
     else
     if (derivs == ADDING_P_P)
          add_all_p_p_up_a(AXX, R23, S12_a, S21_a, S23_a, S32_a, S31_a, n, atran, v1, v2, w1);
     else
*/
     if (derivs == ADDING_U_L)
          add_all_u_l_up_a(PXX, i1, AXX, BXX, CXX, T12, S12, R21, T32, R23_a, T32_a, S32_a, R13_a, T31_a, S31_a, n, atran, v1, v2, w1, w2, work);
     else
     if (derivs == ADDING_L_P)
          add_all_l_p_up_a(PXX, i1, AXX, BXX, CXX, DXX, R23, T32, S32, R12_a, T12_a, S12_a, R21_a, T21_a, S21_a, S32_a, R13_a, T31_a, S31_a, n, atran, atran_a, v1, v2, w1, w2, work);
     else
/*
     if (derivs == ADDING_L_L)
*/
          add_all_l_l_up_a(PXX, i1, AXX, BXX, CXX, T12, S12, R21, R23, T32, S32, R12_a, T12_a, S12_a, R21_a, T21_a, S21_a, R23_a, T32_a, S32_a, R13_a, T31_a, S31_a, n, atran, atran_a, v1, v2, w1, w2, work);


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     i1  = save->i13;

     CXX = save->C13;

     PXX = save->P13;
     AXX = save->A13;
     BXX = save->B31;

     dmat_mul(AXX, R21, n, n, n, DXX);


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
/*
     if (derivs == ADDING_U_P)
          add_all_p_u_up_a(AXX, R21, S32_a, S23_a, S13_a, n, atran, v1, v2);
     else
     if (derivs == ADDING_P_P)
          add_all_p_p_up_a(AXX, R21, S32_a, S23_a, S21_a, S12_a, S13_a, n, atran, v1, v2, w1);
     else
*/
     if (derivs == ADDING_U_L)
          add_all_l_u_up_a(PXX, i1, AXX, BXX, CXX, DXX, R21, T12, S12, R32_a, T32_a, S32_a, R23_a, T23_a, S23_a, R31_a, T13_a, S13_a, n, atran, atran_a, v1, v2, w1, w2, work);
     else
     if (derivs == ADDING_L_P)
          add_all_p_l_up_a(PXX, i1, AXX, BXX, CXX, T32, S32, R23, R21, T12, S32_a, S23_a, R21_a, T12_a, S12_a, R31_a, T13_a, S13_a, n, atran, v1, v2, w1, w2, work);
     else
/*
     if (derivs == ADDING_L_L)
*/
          add_all_l_l_up_a(PXX, i1, AXX, BXX, CXX, T32, S32, R23, R21, T12, S12, R32_a, T32_a, S32_a, R23_a, T23_a, S23_a, R21_a, T12_a, S12_a, R31_a, T13_a, S13_a, n, atran, atran_a, v1, v2, w1, w2, work);
}



/*******************************************************************************
 *
 ******************************************************************************/
#ifdef INCLUDE_DEV_SOURCE
#include "../src2/xrtm_adding_a2.c"
#endif
