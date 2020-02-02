/*******************************************************************************
**
**    Copyright (C) 2007-2020 Greg McGarragh <greg.mcgarragh@colostate.edu>
**
**    This source code is licensed under the GNU General Public License (GPL),
**    Version 3.  See the file COPYING for more details.
**
*******************************************************************************/

#include <gutil.h>
#include <gmath_matrix.h>

#include "xrtm.h"
#include "xrtm_adding.h"
#include "xrtm_doubling.h"
#include "xrtm_matrix.h"
#include "xrtm_save_tree.h"
#include "xrtm_scatter.h"
#include "xrtm_utility.h"


/*******************************************************************************
 *
 ******************************************************************************/
/*
#define CUDA
*/
#ifdef  CUDA
#include "../src2/layer_double_cuda.c"
#else
void layer_double(double **R, double **T,
                  double *Se_m, double *Se_p, double *Sl_m, double *Sl_p,
                  double *L_m, double *L_p,
                  int n, double atran, double lin_fac,
                  int solar, int thermal, int initialize, work_data work) {

     int *i1;

     double *v1;
     double *v2;
     double *v3;

     double **w1;
     double **w2;

     double **A;


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     i1 = get_work1(&work, WORK_IX);

     v1 = get_work1(&work, WORK_DX);
     v2 = get_work1(&work, WORK_DX);
     v3 = get_work1(&work, WORK_DX);

     w1 = get_work1(&work, WORK_DXX);
     w2 = get_work1(&work, WORK_DXX);

     A  = get_work1(&work, WORK_DXX);


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     /* LU of (E - R * R)^T */
     dmat_mul(R, R, n, n, n, w1);
     dmat_i_sub(w1, w1, n);
     dmat_trans(w1, w2, n, n);
     dmat_getrf(w2, n, n, i1);

     /* A = T * P */
     dmat_copy(A, T, n, n);
     dmat_getrs(w2, A, n, n, i1);


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (solar) {
          /* v1 = A * (R * Se_m + Se_p * t) */
          dm_v_mul(R, Se_m, n, n, v1);
          dvec_scale(atran, Se_p, v2, n);
          dvec_add(v1, v2, v2, n);
          dm_v_mul(A, v2, n, n, v1);

          /* v2 = A * (R * Se_p * t + Se_m) */
          dm_v_mul(R, Se_p, n, n, v2);
          dvec_scale(atran, v2, v2, n);
          dvec_add(v2, Se_m, v3, n);
          dm_v_mul(A, v3, n, n, v2);

          /* Se_p = v1 + Se_p */
          dvec_add(v1, Se_p, Se_p, n);

          /* Se_m = v2 + Se_m * t */
          dvec_scale(atran, Se_m, Se_m, n);
          dvec_add(v2, Se_m, Se_m, n);
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (thermal) {
          if (initialize) {
               dvec_copy(L_p, Sl_p, n);
               dvec_copy(L_m, Sl_m, n);
          }

          /* v1 = A * (R * Sl_m + Sl_p + L_p * f) */
          dm_v_mul(R, Sl_m, n, n, v1);
          dvec_add(v1, Sl_p, v1, n);
          dvec_scale(lin_fac, L_p, v2, n);
          dvec_add(v1, v2, v2, n);
          dm_v_mul(A, v2, n, n, v1);

          /* v2 = A * (R * (Sl_p + L_p * f) + Sl_m) */
          dvec_scale(lin_fac, L_p, v2, n);
          dvec_add(Sl_p, v2, v2, n);
          dm_v_mul(R, v2, n, n, v3);
          dvec_add(v3, Sl_m, v3, n);
          dm_v_mul(A, v3, n, n, v2);

          /* Sl_p = v1 + Sl_p */
          dvec_add(v1, Sl_p, Sl_p, n);

          /* Sl_m = v2 + Sl_m + L_m * f */
          dvec_add(v2, Sl_m, Sl_m, n);
          dvec_scale(lin_fac, L_m, v1, n);
          dvec_add(Sl_m, v1, Sl_m, n);


          /* v1 = A * (R * L_m + L_p) */
          dm_v_mul(R, L_m, n, n, v1);
          dvec_add(v1, L_p, v2, n);
          dm_v_mul(A, v2, n, n, v1);

          /* v2 = A * (R * L_p + L_m) */
          dm_v_mul(R, L_p, n, n, v2);
          dvec_add(v2, L_m, v3, n);
          dm_v_mul(A, v3, n, n, v2);

          /* L_p = v1 + L_p */
          dvec_add(v1, L_p, L_p, n);

          /* L_m = v2 + L_m */
          dvec_add(v2, L_m, L_m, n);
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     /* R = R + A * B */
     dmat_mul(A, R, n, n, n, w1);
/*
     dmat_mul(w1, T, n, n, n, w2);
     dmat_add(R, w2, R, n, n);
*/
     dmat_gxgxmx(0, w1, 0, T, 1., R, 1., n, n, n);

     /* T = A * T */
     dmat_mul(A, T, n, n, n, w1);
     dmat_copy(T, w1, n, n);
}
#endif


/*******************************************************************************
 *
 ******************************************************************************/
void layer_double_s(double **R, double **T,
                    double *Se_m, double *Se_p, double *Sl_m, double *Sl_p,
                    double *L_m, double *L_p,
                    int n, double atran, double lin_fac,
                    int solar, int thermal, int initialize, work_data work) {

     int *i1;

     double *v1;
     double *v2;
     double *v3;

     double **w1;
     double **w2;

     double **A;


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     i1 = get_work1(&work, WORK_IX);

     v1 = get_work1(&work, WORK_DX);
     v2 = get_work1(&work, WORK_DX);
     v3 = get_work1(&work, WORK_DX);

     w1 = get_work1(&work, WORK_DXX);
     w2 = get_work1(&work, WORK_DXX);

     A  = get_work1(&work, WORK_DXX);


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
if (1) {
     /* LU of (E - R * R)^T */
     dsym_mul(R, R, n, n, w1, 1);
     dmat_i_sub(w1, w2, n);
     dmat_trans(w2, w1, n, n);
     dmat_getrf(w1, n, n, i1);

     /* A = T * P */
     dmat_copy(A, T, n, n);
     dmat_getrs(w1, A, n, n, i1);
}
else {
     /* LU of (E - R * R)^T */
     dsym_mul(R, R, n, n, w1, 1);
     dmat_i_sub(w1, w1, n);
     dmat_potrf(w1, n);

     /* A = T * P */
     dmat_copy(A, T, n, n);
     dmat_potrs(w1, A, n, n);
}

     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (solar) {
          /* v1 = (R * Se_m + Se_p * t) */
          dm_v_mul(R, Se_m, n, n, v1);
          dvec_scale(atran, Se_p, v2, n);
          dvec_add(v1, v2, v2, n);
          dm_v_mul(A, v2, n, n, v1);

          /* v2 = (R * Se_p * t + Se_m) */
          dm_v_mul(R, Se_p, n, n, v2);
          dvec_scale(atran, v2, v2, n);
          dvec_add(v2, Se_m, v3, n);
          dm_v_mul(A, v3, n, n, v2);

          /* Se_p = v1 + Se_p */
          dvec_add(v1, Se_p, Se_p, n);

          /* Se_m = v2 + Se_m * t */
          dvec_scale(atran, Se_m, Se_m, n);
          dvec_add(v2, Se_m, Se_m, n);
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (thermal) {
          if (initialize) {
               dvec_copy(L_p, Sl_p, n);
               dvec_copy(L_m, Sl_m, n);
          }
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     /* R = R + A * B */
     dmat_mul(A, R, n, n, n, w1);
     dsym_mul(w1, T, n, n, w2, 1);
     dmat_add(R, w2, R, n, n);

     /* T = A * T */
     dsym_mul(A, T, n, n, w1, 1);
     dmat_copy(T, w1, n, n);
}



/*******************************************************************************
 *
 ******************************************************************************/
void layer_double_l(double **R, double **T,
                    double *Se_m, double *Se_p, double *Sl_m, double *Sl_p,
                    double *L_m, double *L_p,
                    double ***R_l, double ***T_l,
                    double **Se_m_l, double **Se_p_l, double **Sl_m_l, double **Sl_p_l,
                    double **L_m_l, double **L_p_l,
                    int n, int n_derivs,
                    double atran, double *atran_l, double lin_fac, double *lin_fac_l,
                    int solar, int thermal, int initialize,
                    uchar *derivs_layers, uchar *derivs_beam, uchar *derivs_thermal, work_data work) {

     int i;

     int *i1;

     double *v1;
     double *v2;
     double *v3;
     double *v4;
     double *v5;

     double *a;
     double *b;
     double *c;
     double *d;

     double *e;
     double *f;

     double **w1;
     double **w2;

     double **P;
     double **A;
     double **B;

     double **A_l;

     i1 = get_work1(&work, WORK_IX);

     v1 = get_work1(&work, WORK_DX);
     v2 = get_work1(&work, WORK_DX);
     v3 = get_work1(&work, WORK_DX);

     a  = get_work1(&work, WORK_DX);
     b  = get_work1(&work, WORK_DX);
     c  = get_work1(&work, WORK_DX);
     d  = get_work1(&work, WORK_DX);

     e  = get_work1(&work, WORK_DX);
     f  = get_work1(&work, WORK_DX);

     w1 = get_work1(&work, WORK_DXX);

     P  = get_work1(&work, WORK_DXX);
     A  = get_work1(&work, WORK_DXX);
     B  = get_work1(&work, WORK_DXX);

     if (n_derivs > 0) {
          v4  = get_work1(&work, WORK_DX);
          v5  = get_work1(&work, WORK_DX);

          w2  = get_work1(&work, WORK_DXX);

          A_l = get_work1(&work, WORK_DXX);
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     /* P = LU of (E - R * R)^T */
     dmat_mul(R, R, n, n, n, P);
     dmat_i_sub(P, w1, n);
     dmat_trans(w1, P, n, n);
     dmat_getrf(P, n, n, i1);

     /* A = T * P */
     dmat_copy(A, T, n, n);
     dmat_getrs(P, A, n, n, i1);

     /* B = R * T */
     dmat_mul(R, T, n, n, n, B);

     if (solar) {
          /* a = R * Se_m + Se_p * t */
          dm_v_mul(R, Se_m, n, n, v1);
          dvec_scale(atran, Se_p, v2, n);
          dvec_add(v1, v2, a, n);

          /* b = R * Se_p * t + Se_m */
          dm_v_mul(R, Se_p, n, n, v1);
          dvec_scale(atran, v1, v1, n);
          dvec_add(v1, Se_m, b, n);
     }

     if (thermal) {
          if (initialize) {
               dvec_copy(L_p, Sl_p, n);
               dvec_copy(L_m, Sl_m, n);
          }

          /* e = R * Sl_m + Sl_p + L_p * f */
          dm_v_mul(R, Sl_m, n, n, v1);
          dvec_add(v1, Sl_p, v1, n);
          dvec_scale(lin_fac, L_p, v2, n);
          dvec_add(v1, v2, e, n);

          /* f = R * (Sl_p + L_p * f) + Sl_m */
          dvec_scale(lin_fac, L_p, v1, n);
          dvec_add(Sl_p, v1, v1, n);
          dm_v_mul(R, v1, n, n, v2);
          dvec_add(v2, Sl_m, f, n);

          /* c = R * L_m + L_p */
          dm_v_mul(R, L_m, n, n, v1);
          dvec_add(v1, L_p, c, n);

          /* d = R * L_p + L_m */
          dm_v_mul(R, L_p, n, n, v1);
          dvec_add(v1, L_m, d, n);
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     for (i = 0; i < n_derivs; ++i) {

          if (derivs_layers[i]) {
               /* A_l = T_l * P + A * (R_l * R + R * R_l) * P */
               dmat_copy(A_l, T_l[i], n, n);
               dmat_getrs(P, A_l, n, n, i1);

               dmat_mul(R_l[i], R, n, n, n, w1);
               dmat_mul(R, R_l[i], n, n, n, w2);
               dmat_add(w1, w2, w1, n, n);
               dmat_getrs(P, w1, n, n, i1);
               dmat_mul(A, w1, n, n, n, w2);

               dmat_add(A_l, w2, A_l, n, n);


               /*---------------------------------------------------------------
                *
                *-------------------------------------------------------------*/
               if (solar && derivs_beam[i]) {
                    /* v1 = A_l * a + A * (R_l * Se_m + R * Se_m_l + Se_p_l * t + Se_p * L(t)) */
                    dm_v_mul(A_l, a, n, n, v1);

                    dm_v_mul(R_l[i], Se_m, n, n, v2);
                    dm_v_mul(R, Se_m_l[i], n, n, v3);
                    dvec_add(v2, v3, v2, n);
                    dvec_scale(atran, Se_p_l[i], v3, n);
                    dvec_add(v2, v3, v2, n);
                    dvec_scale(atran_l[i], Se_p, v3, n);
                    dvec_add(v2, v3, v2, n);

                    dm_v_mul(A, v2, n, n, v3);
                    dvec_add(v1, v3, v1, n);

                    /* v2 = A_l * b + A * [R_l * Se_p * t + R * (Se_p_l * t + Se_p * L(t)) + Se_m_l] */
                    dm_v_mul(A_l, b, n, n, v2);

                    dvec_scale(atran, Se_p_l[i], v3, n);
                    dvec_scale(atran_l[i], Se_p, v4, n);
                    dvec_add(v3, v4, v3, n);
                    dm_v_mul(R, v3, n, n, v4);

                    dm_v_mul(R_l[i], Se_p, n, n, v3);
                    dvec_scale(atran, v3, v3, n);
                    dvec_add(v3, v4, v3, n);

                    dvec_add(v3, Se_m_l[i], v3, n);

                    dm_v_mul(A, v3, n, n, v4);
                    dvec_add(v2, v4, v2, n);

                    /* Se_p_l = v1 + Se_p_l */
                    dvec_add(v1, Se_p_l[i], Se_p_l[i], n);

                    /* Se_m_l = v2 + Se_m_l * t + Se_m * L(t) */
                    dvec_scale(atran, Se_m_l[i], v3, n);
                    dvec_add(v2, v3, v2, n);
                    dvec_scale(atran_l[i], Se_m, v3, n);
                    dvec_add(v2, v3, Se_m_l[i], n);
               }


               /*---------------------------------------------------------------
                *
                *-------------------------------------------------------------*/
               if (thermal && derivs_thermal[i]) {
                    if (initialize) {
                         dvec_copy(L_p_l[i], Sl_p_l[i], n);
                         dvec_copy(L_m_l[i], Sl_m_l[i], n);
                    }

                    /* v1 = */
                    dm_v_mul(A_l, e, n, n, v1);

                    dm_v_mul(R_l[i], Sl_m, n, n, v2);
                    dm_v_mul(R, Sl_m_l[i], n, n, v3);
                    dvec_add(v2, v3, v2, n);
                    dvec_add(v2, Sl_p_l[i], v2, n);
                    dvec_scale(lin_fac, L_p_l[i], v3, n);
                    dvec_add(v2, v3, v2, n);
                    dvec_scale(lin_fac_l[i], L_p, v3, n);
                    dvec_add(v2, v3, v2, n);

                    dm_v_mul(A, v2, n, n, v3);
                    dvec_add(v1, v3, v1, n);

                    /* v2 = */
                    dm_v_mul(A_l, f, n, n, v2);

                    dvec_scale(lin_fac, L_p, v3, n);
                    dvec_add(Sl_p, v3, v3, n);
                    dm_v_mul(R_l[i], v3, n, n, v4);

                    dvec_scale(lin_fac, L_p_l[i], v3, n);
                    dvec_add(Sl_p_l[i], v3, v3, n);
                    dvec_scale(lin_fac_l[i], L_p, v5, n);
                    dvec_add(v3, v5, v3, n);
                    dm_v_mul(R, v3, n, n, v5);

                    dvec_add(v4, v5, v3, n);

                    dvec_add(v3, Sl_m_l[i], v3, n);

                    dm_v_mul(A, v3, n, n, v4);
                    dvec_add(v2, v4, v2, n);

                    /* Sl_p = */
                    dvec_add(v1, Sl_p_l[i], Sl_p_l[i], n);

                    /* Sl_m = */
                    dvec_add(v2, Sl_m_l[i], Sl_m_l[i], n);

                    dvec_scale(lin_fac, L_m_l[i], v2, n);
                    dvec_add(Sl_m_l[i], v2, Sl_m_l[i], n);

                    dvec_scale(lin_fac_l[i], L_m, v2, n);
                    dvec_add(Sl_m_l[i], v2, Sl_m_l[i], n);


                    /* v1 = */
                    dm_v_mul(A_l, c, n, n, v1);

                    dm_v_mul(R_l[i], L_m, n, n, v2);
                    dm_v_mul(R, L_m_l[i], n, n, v3);
                    dvec_add(v2, v3, v2, n);
                    dvec_add(v2, L_p_l[i], v2, n);

                    dm_v_mul(A, v2, n, n, v3);
                    dvec_add(v1, v3, v1, n);

                    /* v2 = */
                    dm_v_mul(A_l, d, n, n, v2);

                    dm_v_mul(R_l[i], L_p, n, n, v3);
                    dm_v_mul(R, L_p_l[i], n, n, v4);
                    dvec_add(v3, v4, v3, n);
                    dvec_add(v3, L_m_l[i], v3, n);

                    dm_v_mul(A, v3, n, n, v4);
                    dvec_add(v2, v4, v2, n);

                    /* L_p = */
                    dvec_add(v1, L_p_l[i], L_p_l[i], n);

                    /* L_m = */
                    dvec_add(v2, L_m_l[i], L_m_l[i], n);
               }


               /*---------------------------------------------------------------
                *
                *-------------------------------------------------------------*/

               /* R_l = R_l + A_l * B + A * (R_l * T + R * T_l) */
               dmat_mul(R_l[i], T, n, n, n, w1);
               dmat_mul(R, T_l[i], n, n, n, w2);
               dmat_add(w1, w2, w1, n, n);
               dmat_mul(A, w1, n, n, n, w2);
               dmat_add(R_l[i], w2, R_l[i], n, n);

               dmat_mul(A_l, B, n, n, n, w1);
               dmat_add(R_l[i], w1, R_l[i], n, n);

               /* T_l = A_l * T + A * T_l */
               dmat_mul(A_l, T, n, n, n, w1);

               dmat_mul(A, T_l[i], n, n, n, w2);

               dmat_add(w1, w2, T_l[i], n, n);
          }


          /*--------------------------------------------------------------------
           *
           *------------------------------------------------------------------*/
          else {
               if (solar && derivs_beam[i]) {

                    /* v1 = A_l * a + A * (R_l * Se_m + R * Se_m_l + Se_p_l * t + Se_p * L(t)) */
                    /* v1 =           A * (             R * Se_m_l + Se_p_l * t + Se_p * L(t)) */
                    dm_v_mul(R, Se_m_l[i], n, n, v1);
                    dvec_scale(atran, Se_p_l[i], v2, n);
                    dvec_add(v1, v2, v1, n);
                    dvec_scale(atran_l[i], Se_p, v2, n);
                    dvec_add(v1, v2, v2, n);

                    dm_v_mul(A, v2, n, n, v1);

                    /* v2 = A_l * b + A * [R_l * Se_p * t + R * (Se_p_l * t + Se_p * L(t)) + Se_m_l] */
                    /* v2 =           A * [                 R * (Se_p_l * t + Se_p * L(t)) + Se_m_l] */
                    dvec_scale(atran, Se_p_l[i], v2, n);
                    dvec_scale(atran_l[i], Se_p, v3, n);
                    dvec_add(v2, v3, v2, n);
                    dm_v_mul(R, v2, n, n, v3);

                    dvec_add(v3, Se_m_l[i], v3, n);

                    dm_v_mul(A, v3, n, n, v2);

                    /* Se_p_l = v1 + Se_p_l */
                    dvec_add(v1, Se_p_l[i], Se_p_l[i], n);

                    /* Se_m_l = v2 + Se_m_l * t + Se_m * L(t) */
                    dvec_scale(atran, Se_m_l[i], v3, n);
                    dvec_add(v2, v3, v2, n);
                    dvec_scale(atran_l[i], Se_m, v3, n);
                    dvec_add(v2, v3, Se_m_l[i], n);
               }


               /*---------------------------------------------------------------
                *
                *-------------------------------------------------------------*/
               if (thermal && derivs_thermal[i]) {
                    if (initialize) {
                         dvec_copy(L_p_l[i], Sl_p_l[i], n);
                         dvec_copy(L_m_l[i], Sl_m_l[i], n);
                    }

                    /* v1 = */
                    dm_v_mul(R, Sl_m_l[i], n, n, v2);
                    dvec_add(v2, Sl_p_l[i], v2, n);
                    dvec_scale(lin_fac, L_p_l[i], v3, n);
                    dvec_add(v2, v3, v2, n);
                    dvec_scale(lin_fac_l[i], L_p, v3, n);
                    dvec_add(v2, v3, v2, n);

                    dm_v_mul(A, v2, n, n, v1);

                    /* v2 = */
                    dvec_scale(lin_fac, L_p_l[i], v3, n);
                    dvec_add(Sl_p_l[i], v3, v3, n);
                    dvec_scale(lin_fac_l[i], L_p, v5, n);
                    dvec_add(v3, v5, v3, n);
                    dm_v_mul(R, v3, n, n, v5);

                    dvec_add(v5, Sl_m_l[i], v3, n);

                    dm_v_mul(A, v3, n, n, v2);

                    /* Sl_p = */
                    dvec_add(v1, Sl_p_l[i], Sl_p_l[i], n);

                    /* Sl_m = */
                    dvec_add(v2, Sl_m_l[i], Sl_m_l[i], n);

                    dvec_scale(lin_fac, L_m_l[i], v2, n);
                    dvec_add(Sl_m_l[i], v2, Sl_m_l[i], n);

                    dvec_scale(lin_fac_l[i], L_m, v2, n);
                    dvec_add(Sl_m_l[i], v2, Sl_m_l[i], n);


                    /* v1 = */
                    dm_v_mul(R, L_m_l[i], n, n, v2);
                    dvec_add(v2, L_p_l[i], v2, n);

                    dm_v_mul(A, v2, n, n, v1);

                    /* v2 = */
                    dm_v_mul(R, L_p_l[i], n, n, v3);
                    dvec_add(v3, L_m_l[i], v3, n);

                    dm_v_mul(A, v3, n, n, v2);

                    /* L_p = */
                    dvec_add(v1, L_p_l[i], L_p_l[i], n);

                    /* L_m = */
                    dvec_add(v2, L_m_l[i], L_m_l[i], n);
               }
          }
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (solar) {
          /* Se_p = A * a + Se_p */
          dm_v_mul(A, a, n, n, v1);
          dvec_add(v1, Se_p, Se_p, n);

          /* Se_m = A * b + Se_m * t */
          dm_v_mul(A, b, n, n, v1);
          dvec_scale(atran, Se_m, Se_m, n);
          dvec_add(v1, Se_m, Se_m, n);
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (thermal) {
          /* Sl_p = A * e + Se_p */
          dm_v_mul(A, e, n, n, v1);
          dvec_add(v1, Sl_p, Sl_p, n);

          /* Sl_m = A * f + Se_m + L_m * f */
          dm_v_mul(A, f, n, n, v1);
          dvec_add(v1, Sl_m, Sl_m, n);
          dvec_scale(lin_fac, L_m, v1, n);
          dvec_add(Sl_m, v1, Sl_m, n);


          /* L_p = A * c + L_p */
          dm_v_mul(A, c, n, n, v1);
          dvec_add(v1, L_p, L_p, n);

          /* L_m = A * d + L_m */
          dm_v_mul(A, d, n, n, v1);
          dvec_add(v1, L_m, L_m, n);
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     /* R = R + A * B */
     dmat_mul(A, B, n, n, n, w1);
     dmat_add(R, w1, R, n, n);

     /* T = A * T */
     dmat_mul(A, T, n, n, n, w1);
     dmat_copy(T, w1, n, n);
}



/*******************************************************************************
 *
 ******************************************************************************/
void layer_double_s_l(double **R, double **T,
                      double *Se_m, double *Se_p, double *Sl_m, double *Sl_p,
                      double *L_m, double *L_p,
                      double ***R_l, double ***T_l,
                      double **Se_m_l, double **Se_p_l, double **Sl_m_l, double **Sl_p_l,
                      double **L_m_l, double **L_p_l,
                      int n, int n_derivs,
                      double atran, double *atran_l, double lin_fac, double *lin_fac_l,
                      int solar, int thermal, int initialize,
                      uchar *derivs_layers, uchar *derivs_beam, uchar *derivs_thermal, work_data work) {

     int i;

     int *i1;

     double *v1;
     double *v2;
     double *v3;
     double *v4;
/*
     double *a;
     double *b;
     double *c;
     double *d;

     double *e;
     double *f;
*/
     double **w1;
     double **w2;
     double **w3;

     double **P;
     double **A;
     double **B;
     double **C;
     double **D;

     i1 = get_work1(&work, WORK_IX);

     v1 = get_work1(&work, WORK_DX);
     v2 = get_work1(&work, WORK_DX);
     v3 = get_work1(&work, WORK_DX);

     w1 = get_work1(&work, WORK_DXX);

     P  = get_work1(&work, WORK_DXX);
     A  = get_work1(&work, WORK_DXX);
     B  = get_work1(&work, WORK_DXX);

     if (n_derivs > 0) {
          v4 = get_work1(&work, WORK_DX);

          w2 = get_work1(&work, WORK_DXX);
          w3 = get_work1(&work, WORK_DXX);

          C  = get_work1(&work, WORK_DXX);
          D  = get_work1(&work, WORK_DXX);
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
if (1) {
     /* P = LU of (E - R * R)^T */
     dsym_mul(R, R, n, n, P, 1);
     dmat_i_sub(P, w1, n);
     dmat_trans(w1, P, n, n);
     dmat_getrf(P, n, n, i1);

     /* A = T * P */
     dmat_copy(A, T, n, n);
     dmat_getrs(P, A, n, n, i1);

     /* B = R * T */
     dmat_mul(R, T, n, n, n, B);
}
else {
     /* w1 = LU of (E - R * R)^T */
     dsym_mul(R, R, n, n, P, 1);
     dmat_i_sub(P, P, n);
     dmat_potrf(P, n);

     /* A = T * P */
     dmat_copy(A, T, n, n);
     dmat_potrs(P, A, n, n);

     /* B = R * T */
     dmat_mul(R, T, n, n, n, B);
}
     if (solar) {

     }

     if (thermal) {
          if (initialize) {
               dvec_copy(L_p, Sl_p, n);
               dvec_copy(L_m, Sl_m, n);
          }
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     for (i = 0; i < n_derivs; ++i) {
          if (derivs_layers[i]) {
if (1) {
               /* C = T_l * P */
               dmat_copy(C, T_l[i], n, n);
               dmat_getrs(P, C, n, n, i1);

               /* D = A * (R_l * R + R * R_l) * P */
               dsym_mta(R_l[i], R, n, n, w2);
               dmat_getrs(P, w2, n, n, i1);
               dmat_mul(A, w2, n, n, n, D);
}
else {
               /* C = T_l * P */
               dmat_copy(C, T_l[i], n, n);
               dmat_potrs(P, C, n, n);

               /* D = A * (R_l * R + R * R_l) * P */
               dsym_mta(R_l[i], R, n, n, w2);
               dmat_potrs(P, w2, n, n);
               dmat_mul(A, w2, n, n, n, D);
}

               /*---------------------------------------------------------------
                *
                *-------------------------------------------------------------*/
               if (solar && derivs_beam[i]) {
                    /* an = R * Se_m + Se_p * t */
                    dm_v_mul(R, Se_m, n, n, v1);
                    dvec_scale(atran, Se_p, v2, n);
                    dvec_add(v1, v2, v1, n);

                    /* v1 = C * a + D * a + A * (R_l * Se_m + R * Se_m_l + Se_p_l * t + Se_p * L(t)) */
                    dm_v_mul(C, v1, n, n, v2);
                    dm_v_mul(D, v1, n, n, v3);
                    dvec_add(v2, v3, v1, n);

                    dm_v_mul(R_l[i], Se_m, n, n, v2);
                    dm_v_mul(R, Se_m_l[i], n, n, v3);
                    dvec_add(v2, v3, v2, n);
                    dvec_scale(atran, Se_p_l[i], v3, n);
                    dvec_add(v2, v3, v2, n);
                    dvec_scale(atran_l[i], Se_p, v3, n);
                    dvec_add(v2, v3, v2, n);

                    dm_v_mul(A, v2, n, n, v3);
                    dvec_add(v1, v3, v1, n);

                    /* bn = R * Se_m + Se_p * t */
                    dm_v_mul(R, Se_p, n, n, v2);
                    dvec_scale(atran, v2, v2, n);
                    dvec_add(v2, Se_m, v2, n);

                    /* v2 = C * b + D * b + A * [R_l * Se_p * t + R * (Se_p_l * t + Se_p * L(t)) + Se_m_l] */
                    dm_v_mul(C, v2, n, n, v3);
                    dm_v_mul(D, v2, n, n, v4);
                    dvec_add(v3, v4, v2, n);

                    dvec_scale(atran, Se_p_l[i], v3, n);
                    dvec_scale(atran_l[i], Se_p, v4, n);
                    dvec_add(v3, v4, v3, n);
                    dm_v_mul(R, v3, n, n, v4);

                    dm_v_mul(R_l[i], Se_p, n, n, v3);
                    dvec_scale(atran, v3, v3, n);
                    dvec_add(v3, v4, v3, n);

                    dvec_add(v3, Se_m_l[i], v3, n);

                    dm_v_mul(A, v3, n, n, v4);
                    dvec_add(v2, v4, v2, n);

                    /* V(n + 1)p = v1 + Se_p_l */
                    dvec_add(v1, Se_p_l[i], Se_p_l[i], n);

                    /* V(n + 1)m = v2 + Se_m_l * t + Se_m * L(t) */
                    dvec_scale(atran, Se_m_l[i], v3, n);
                    dvec_add(v2, v3, v2, n);
                    dvec_scale(atran_l[i], Se_m, v3, n);
                    dvec_add(v2, v3, Se_m_l[i], n);
               }


               /*---------------------------------------------------------------
                *
                *-------------------------------------------------------------*/
               if (thermal && derivs_thermal[i]) {

               }


               /*---------------------------------------------------------------
                *
                *-------------------------------------------------------------*/

               /* R_l = R_l + (C + D) * B + A * (R_l * T + R * T_l) */
               dmat_add(C, D, w1, n, n);

               dmat_mul(R_l[i], T, n, n, n, w2);
               dmat_mul(R, T_l[i], n, n, n, w3);
               dmat_add(w2, w3, w2, n, n);

               dsym_mma(w1, B, A, w2, n, n, w3, 0);

               dsym_add(R_l[i], w3, R_l[i], n, 1);

               /* T_l = (C + D) * T + A * T_l */
               dsym_mul(w1, T, n, n, w2, 0);
               dsym_mul(A , T_l[i], n, n, w3, 0);

               dsym_add(w2, w3, T_l[i], n, 1);
          }


          /*--------------------------------------------------------------------
           *
           *------------------------------------------------------------------*/
          else {
               if (solar && derivs_beam[i]) {
                    /* v1 = C * a + D * a + A * (R_l * Se_m + R * Se_m_l + Se_p_l * t + Se_p * L(t)) */
                    /* v1 =                 A * (             R * Se_m_l + Se_p_l * t + Se_p * L(t)) */
                    dm_v_mul(R, Se_m_l[i], n, n, v1);
                    dvec_scale(atran, Se_p_l[i], v2, n);
                    dvec_add(v1, v2, v1, n);
                    dvec_scale(atran_l[i], Se_p, v2, n);
                    dvec_add(v1, v2, v2, n);

                    dm_v_mul(A, v2, n, n, v1);

                    /* v2 = C * b + D * b + A * [R_l * Se_p * t + R * (Se_p_l * t + Se_p * L(t)) + Se_m_l] */
                    /* v2 =                 A * [                 R * (Se_p_l * t + Se_p * L(t)) + Se_m_l] */
                    dvec_scale(atran, Se_p_l[i], v2, n);
                    dvec_scale(atran_l[i], Se_p, v3, n);
                    dvec_add(v2, v3, v2, n);
                    dm_v_mul(R, v2, n, n, v3);

                    dvec_add(v3, Se_m_l[i], v3, n);

                    dm_v_mul(A, v3, n, n, v2);

                    /* V(n + 1)p = v1 + Se_p_l */
                    dvec_add(v1, Se_p_l[i], Se_p_l[i], n);

                    /* V(n + 1)m = v2 + Se_m_l * t + Se_m * L(t) */
                    dvec_scale(atran, Se_m_l[i], v3, n);
                    dvec_add(v2, v3, v2, n);
                    dvec_scale(atran_l[i], Se_m, v3, n);
                    dvec_add(v2, v3, Se_m_l[i], n);
               }


               /*---------------------------------------------------------------
                *
                *-------------------------------------------------------------*/
               if (thermal && derivs_thermal[i]) {

               }
          }
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (solar) {
          /* v1 = A * (R * Se_m + Se_p * t) */
          dm_v_mul(R, Se_m, n, n, v1);
          dvec_scale(atran, Se_p, v2, n);
          dvec_add(v1, v2, v2, n);
          dm_v_mul(A, v2, n, n, v1);

          /* v2 = A * (R * Se_p * t + Se_m) */
          dm_v_mul(R, Se_p, n, n, v2);
          dvec_scale(atran, v2, v2, n);
          dvec_add(v2, Se_m, v3, n);
          dm_v_mul(A, v3, n, n, v2);

          /* Se_p = v1 + Se_p */
          dvec_add(v1, Se_p, Se_p, n);

          /* Se_m = v2 + Se_m * t */
          dvec_scale(atran, Se_m, Se_m, n);
          dvec_add(v2, Se_m, Se_m, n);
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (thermal) {

     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     /* R = R + A * B */
     dsym_mul(A, B, n, n, w1, 1);
     dmat_add(R, w1, R, n, n);

     /* T = A * T */
     dsym_mul(A, T, n, n, w1, 1);
     dmat_copy(T, w1, n, n);
}



/*******************************************************************************
 *
 ******************************************************************************/
#ifdef INCLUDE_DEV_SOURCE
#include "../src2/xrtm_doubling2.c"
#endif
