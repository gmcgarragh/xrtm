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
                  int n, double atran, double lin_fac,
                  int flag1, int flag2, int flag3, work_data work) {

     int *i1;

     double *v1;
     double *v2;
     double *v3;

     double *a_p;
     double *a_m;

     double **w1;
     double **w2;

     double **A;


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     i1  = get_work1(&work, WORK_IX);

     v1  = get_work1(&work, WORK_DX);
     v2  = get_work1(&work, WORK_DX);
     v3  = get_work1(&work, WORK_DX);

     a_p = get_work1(&work, WORK_DX);
     a_m = get_work1(&work, WORK_DX);

     w1  = get_work1(&work, WORK_DXX);
     w2  = get_work1(&work, WORK_DXX);

     A   = get_work1(&work, WORK_DXX);


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     /* LU of (E - Rn * Rn)^T */
     dmat_mul(R, R, n, n, n, w1);
     dmat_i_sub(w1, w1, n);
     dmat_trans(w1, w2, n, n);
     dmat_getrf(w2, n, n, i1);

     /* An = Tn * Pn */
     dmat_copy(A, T, n, n);
     dmat_getrs(w2, A, n, n, i1);


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (flag1) {
          /* v1 = An * (Rn * Snm + Snp * t) */
          dm_v_mul(R, Se_m, n, n, v1);
          dvec_scale(atran, Se_p, v2, n);
          dvec_add(v1, v2, v2, n);
          dm_v_mul(A, v2, n, n, v1);

          /* v2 = An * (Rn * Snp * t + Snm) */
          dm_v_mul(R, Se_p, n, n, v2);
          dvec_scale(atran, v2, v2, n);
          dvec_add(v2, Se_m, v3, n);
          dm_v_mul(A, v3, n, n, v2);

          /* Snp = v1 + Snp */
          dvec_add(v1, Se_p, Se_p, n);

          /* Snm = v2 + Snm * t */
          dvec_scale(atran, Se_m, Se_m, n);
          dvec_add(v2, Se_m, Se_m, n);
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (flag2) {
          if (flag3) {
               dvec_copy(a_p, Sl_p, n);
               dvec_copy(a_m, Sl_m, n);
          }

          /* v1 = An * (Rn * Snm + Snp + anp * f) */
          dm_v_mul(R, Sl_m, n, n, v1);
          dvec_add(v1, Sl_p, v1, n);
          dvec_scale(lin_fac, a_p, v2, n);
          dvec_add(v1, v2, v2, n);
          dm_v_mul(A, v2, n, n, v1);

          /* v2 = An * (Rn * (Snp + anp * f) + Snm) */
          dvec_scale(lin_fac, a_p, v2, n);
          dvec_add(Sl_p, v2, v2, n);
          dm_v_mul(R, v2, n, n, v3);
          dvec_add(v3, Sl_m, v3, n);
          dm_v_mul(A, v3, n, n, v2);

          /* Snp = v1 + Snp */
          dvec_add(v1, Sl_p, Sl_p, n);

          /* Snm = v2 + Snm + a_m * f */
          dvec_add(v2, Sl_m, Sl_m, n);
          dvec_scale(lin_fac, a_m, v1, n);
          dvec_add(Sl_m, v1, Sl_m, n);


          /* v1 = An * (Rn * anm + anp) */
          dm_v_mul(R, a_m, n, n, v1);
          dvec_add(v1, a_p, v2, n);
          dm_v_mul(A, v2, n, n, v1);

          /* v2 = An * (Rn * anp + anm) */
          dm_v_mul(R, a_p, n, n, v2);
          dvec_add(v2, a_m, v3, n);
          dm_v_mul(A, v3, n, n, v2);

          /* anp = v1 + anp */
          dvec_add(v1, a_p, a_p, n);

          /* anm = v2 + anm */
          dvec_add(v2, a_m, a_m, n);
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     /* R(n+1) = Rn + An * Bn */
     dmat_mul(A, R, n, n, n, w1);
     dmat_mul(w1, T, n, n, n, w2);
     dmat_add(R, w2, R, n, n);

     /* T(n+1) = An * Tn */
     dmat_mul(A, T, n, n, n, w1);
     dmat_copy(T, w1, n, n);
}
#endif


/*******************************************************************************
 *
 ******************************************************************************/
void layer_double_s(double **R, double **T, double *S_m, double *S_p,
                    int n, double atran, work_data work) {

     int *i1;

     double *v1;
     double *v2;
     double *v3;

     double **w1;
     double **w2;

     double **A;

     i1 = get_work1(&work, WORK_IX);

     v1 = get_work1(&work, WORK_DX);
     v2 = get_work1(&work, WORK_DX);
     v3 = get_work1(&work, WORK_DX);

     w1 = get_work1(&work, WORK_DXX);
     w2 = get_work1(&work, WORK_DXX);

     A = get_work1(&work, WORK_DXX);
if (1) {
     /* LU of (E - Rn * Rn)^T */
     dsym_mul(R, R, n, n, w1, 1);
     dmat_i_sub(w1, w2, n);
     dmat_trans(w2, w1, n, n);
     dmat_getrf(w1, n, n, i1);

     /* An = Tn * Pn */
     dmat_copy(A, T, n, n);
     dmat_getrs(w1, A, n, n, i1);
}
else {
     /* LU of (E - Rn * Rn)^T */
     dsym_mul(R, R, n, n, w1, 1);
     dmat_i_sub(w1, w1, n);
     dmat_potrf(w1, n);

     /* An = Tn * Pn */
     dmat_copy(A, T, n, n);
     dmat_potrs(w1, A, n, n);
}
if (S_p) {
     /* v1 = (Rn * Snm + Snp * t) */
     dm_v_mul(R, S_m, n, n, v1);
     dvec_scale(atran, S_p, v2, n);
     dvec_add(v1, v2, v2, n);
     dm_v_mul(A, v2, n, n, v1);
}
if (S_m) {
     /* v2 = (Rn * Snp * t + Snm) */
     dm_v_mul(R, S_p, n, n, v2);
     dvec_scale(atran, v2, v2, n);
     dvec_add(v2, S_m, v3, n);
     dm_v_mul(A, v3, n, n, v2);
}
if (S_p) {
     /* Snp = v1 + Snp */
     dvec_add(v1, S_p, S_p, n);
}
if (S_m) {
     /* Snm = v2 + Snm * t */
     dvec_scale(atran, S_m, S_m, n);
     dvec_add(v2, S_m, S_m, n);
}
     /* R(n+1) = Rn + An * Bn */
     dmat_mul(A, R, n, n, n, w1);
     dsym_mul(w1, T, n, n, w2, 1);
     dmat_add(R, w2, R, n, n);

     /* T(n+1) = An * Tn */
     dsym_mul(A, T, n, n, w1, 1);
     dmat_copy(T, w1, n, n);
}



/*******************************************************************************
 *
 ******************************************************************************/
void layer_double_l(double **R, double **T, double *S_m, double *S_p,
                    double ***R_l, double ***T_l, double **S_m_l, double **S_p_l,
                    int n, int n_derivs, double atran, double *atran_l,
                    uchar *derivs_h, uchar *derivs_p, work_data work) {

     int i;

     int *i1;

     double *v1;
     double *v2;
     double *v3;
     double *v4;

     double **w1;
     double **w2;

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
     w2 = get_work1(&work, WORK_DXX);

     P  = get_work1(&work, WORK_DXX);
     A  = get_work1(&work, WORK_DXX);
     B  = get_work1(&work, WORK_DXX);

     if (n_derivs > 0) {
          v4 = get_work1(&work, WORK_DX);

          C  = get_work1(&work, WORK_DXX);
          D  = get_work1(&work, WORK_DXX);
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     /* P = LU of (E - Rn * Rn)^T */
     dmat_mul(R, R, n, n, n, P);
     dmat_i_sub(P, w1, n);
     dmat_trans(w1, P, n, n);
     dmat_getrf(P, n, n, i1);

     /* An = Tn * Pn */
     dmat_copy(A, T, n, n);
     dmat_getrs(P, A, n, n, i1);

     /* Bn = Rn * Tn */
     dmat_mul(R, T, n, n, n, B);


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     for (i = 0; i < n_derivs; ++i) {

          if (derivs_h[i]) {
               /* Cn = Wn * Pn */
               dmat_copy(C, T_l[i], n, n);
               dmat_getrs(P, C, n, n, i1);

               /* Dn = An * (Un_l * Rn + Rn * Un_l) * Pn */
               dmat_mul(R_l[i], R, n, n, n, w1);
               dmat_mul(R, R_l[i], n, n, n, w2);
               dmat_add(w1, w2, w1, n, n);
               dmat_getrs(P, w1, n, n, i1);
               dmat_mul(A, w1, n, n, n, D);


               /*---------------------------------------------------------------
                *
                *-------------------------------------------------------------*/
if (S_p_l) {
               /* an = Rn * Snm + Snp * t */
               dm_v_mul(R, S_m, n, n, v1);
               dvec_scale(atran, S_p, v2, n);
               dvec_add(v1, v2, v1, n);

               /* v1 = Cn * a + Dn * a + An * (Un_l * Snm + Rn * Vnm + Vnp * t + Snp * L(t)) */
               dm_v_mul(C, v1, n, n, v2);
               dm_v_mul(D, v1, n, n, v3);
               dvec_add(v2, v3, v1, n);

               dm_v_mul(R_l[i], S_m, n, n, v2);
               dm_v_mul(R, S_m_l[i], n, n, v3);
               dvec_add(v2, v3, v2, n);
               dvec_scale(atran, S_p_l[i], v3, n);
               dvec_add(v2, v3, v2, n);
               dvec_scale(atran_l[i], S_p, v3, n);
               dvec_add(v2, v3, v2, n);

               dm_v_mul(A, v2, n, n, v3);
               dvec_add(v1, v3, v1, n);
}
if (S_m_l) {
               /* bn = Rn * Snm + Snp * t */
               dm_v_mul(R, S_p, n, n, v2);
               dvec_scale(atran, v2, v2, n);
               dvec_add(v2, S_m, v2, n);

               /* v2 = Cn * b + Dn * b + An * [Un_l * Snp * t + Rn * (Vnp * t + Snp * L(t)) + Vnm] */
               dm_v_mul(C, v2, n, n, v3);
               dm_v_mul(D, v2, n, n, v4);
               dvec_add(v3, v4, v2, n);

               dvec_scale(atran, S_p_l[i], v3, n);
               dvec_scale(atran_l[i], S_p, v4, n);
               dvec_add(v3, v4, v3, n);
               dm_v_mul(R, v3, n, n, v4);

               dm_v_mul(R_l[i], S_p, n, n, v3);
               dvec_scale(atran, v3, v3, n);
               dvec_add(v3, v4, v3, n);

               dvec_add(v3, S_m_l[i], v3, n);

               dm_v_mul(A, v3, n, n, v4);
               dvec_add(v2, v4, v2, n);
}
if (S_p_l) {
               /* V(n + 1)p = v1 + Vnp */
               dvec_add(v1, S_p_l[i], S_p_l[i], n);
}
if (S_m_l) {
               /* V(n + 1)m = v2 + Vnm * t + Snm * L(t) */
               dvec_scale(atran, S_m_l[i], v3, n);
               dvec_add(v2, v3, v2, n);
               dvec_scale(atran_l[i], S_m, v3, n);
               dvec_add(v2, v3, S_m_l[i], n);
}
               /* U(n + 1)_l = Un_l + (Cn + Dn) * Bn + An * (Un_l * Tn + Rn * Wn) */
               dmat_mul(R_l[i], T, n, n, n, w1);
               dmat_mul(R, T_l[i], n, n, n, w2);
               dmat_add(w1, w2, w1, n, n);
               dmat_mul(A, w1, n, n, n, w2);
               dmat_add(R_l[i], w2, R_l[i], n, n);

               dmat_add(C, D, w1, n, n);
               dmat_mul(w1, B, n, n, n, w2);
               dmat_add(R_l[i], w2, R_l[i], n, n);

               /* W(n + 1) = (Cn + Dn) * Tn + An * Wn */
               dmat_add(C, D, w1, n, n);
               dmat_mul(w1, T, n, n, n, w2);

               dmat_mul(A, T_l[i], n, n, n, w1);

               dmat_add(w2, w1, T_l[i], n, n);
}

          /*----------------------------------------------------------
           *
           *--------------------------------------------------------*/
          else
          if ((S_p_l || S_p_l) && derivs_p[i]) {

if (S_p_l) {
               /* v1 = Cn * a + Dn * a + An * (Un_l * Snm + Rn * Vnm + Vnp * t + Snp * L(t)) */
               /* v1 =                   An * (             Rn * Vnm + Vnp * t + Snp * L(t)) */
               dm_v_mul(R, S_m_l[i], n, n, v1);
               dvec_scale(atran, S_p_l[i], v2, n);
               dvec_add(v1, v2, v1, n);
               dvec_scale(atran_l[i], S_p, v2, n);
               dvec_add(v1, v2, v2, n);

               dm_v_mul(A, v2, n, n, v1);
}
if (S_m_l) {
               /* v2 = Cn * b + Dn * b + An * [Un_l * Snp * t + Rn * (Vnp * t + Snp * L(t)) + Vnm] */
               /* v2 =                   An * [                 Rn * (Vnp * t + Snp * L(t)) + Vnm] */
               dvec_scale(atran, S_p_l[i], v2, n);
               dvec_scale(atran_l[i], S_p, v3, n);
               dvec_add(v2, v3, v2, n);
               dm_v_mul(R, v2, n, n, v3);

               dvec_add(v3, S_m_l[i], v3, n);

               dm_v_mul(A, v3, n, n, v2);
}
if (S_p_l) {
               /* V(n + 1)p = v1 + Vnp */
               dvec_add(v1, S_p_l[i], S_p_l[i], n);
}
if (S_m_l) {
               /* V(n + 1)m = v2 + Vnm * t + Snm * L(t) */
               dvec_scale(atran, S_m_l[i], v3, n);
               dvec_add(v2, v3, v2, n);
               dvec_scale(atran_l[i], S_m, v3, n);
               dvec_add(v2, v3, S_m_l[i], n);
}

          }
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
if (S_p) {
     /* v1 = A * (Rn * Snm + Snp * t) */
     dm_v_mul(R, S_m, n, n, v1);
     dvec_scale(atran, S_p, v2, n);
     dvec_add(v1, v2, v2, n);
     dm_v_mul(A, v2, n, n, v1);
}
if (S_m) {
     /* v2 = A * (Rn * Snp * t + Snm) */
     dm_v_mul(R, S_p, n, n, v2);
     dvec_scale(atran, v2, v2, n);
     dvec_add(v2, S_m, v3, n);
     dm_v_mul(A, v3, n, n, v2);
}
if (S_p) {
     /* Snp = v1 + Snp */
     dvec_add(v1, S_p, S_p, n);
}
if (S_m) {
     /* Snm = v2 + Snm * t */
     dvec_scale(atran, S_m, S_m, n);
     dvec_add(v2, S_m, S_m, n);
}
     /* R(n+1) = Rn + An * Bn */
     dmat_mul(A, B, n, n, n, w1);
     dmat_add(R, w1, R, n, n);

     /* T(n+1) = An * Tn */
     dmat_mul(A, T, n, n, n, w1);
     dmat_copy(T, w1, n, n);
}



/*******************************************************************************
 *
 ******************************************************************************/
void layer_double_s_l(double **R, double **T, double *S_m, double *S_p,
                      double ***R_l, double ***T_l, double **S_m_l, double **S_p_l,
                      int n, int n_derivs, double atran, double *atran_l,
                      uchar *derivs_h, uchar *derivs_p, work_data work) {

     int i;

     int *i1;

     double *v1;
     double *v2;
     double *v3;
     double *v4;

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
     w2 = get_work1(&work, WORK_DXX);

     P  = get_work1(&work, WORK_DXX);
     A  = get_work1(&work, WORK_DXX);
     B  = get_work1(&work, WORK_DXX);

     if (n_derivs > 0) {
          v4 = get_work1(&work, WORK_DX);

          w3 = get_work1(&work, WORK_DXX);

          C  = get_work1(&work, WORK_DXX);
          D  = get_work1(&work, WORK_DXX);
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
if (1) {
     /* P = LU of (E - Rn * Rn)^T */
     dsym_mul(R, R, n, n, P, 1);
     dmat_i_sub(P, w1, n);
     dmat_trans(w1, P, n, n);
     dmat_getrf(P, n, n, i1);

     /* An = Tn * Pn */
     dmat_copy(A, T, n, n);
     dmat_getrs(P, A, n, n, i1);

     /* Bn = Rn * Tn */
     dmat_mul(R, T, n, n, n, B);
}
else {
     /* w1 = LU of (E - Rn * Rn)^T */
     dsym_mul(R, R, n, n, P, 1);
     dmat_i_sub(P, P, n);
     dmat_potrf(P, n);

     /* An = Tn * Pn */
     dmat_copy(A, T, n, n);
     dmat_potrs(P, A, n, n);

     /* Bn = Rn * Tn */
     dmat_mul(R, T, n, n, n, B);
}

     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     for (i = 0; i < n_derivs; ++i) {
          if (derivs_h[i]) {
 if (1) {
               /* Cn = Wn * Pn */
               dmat_copy(C, T_l[i], n, n);
               dmat_getrs(P, C, n, n, i1);

               /* Dn = An * (Un_l * Rn + Rn * Un_l) * Pn */
               dsym_mta(R_l[i], R, n, n, w2);
               dmat_getrs(P, w2, n, n, i1);
               dmat_mul(A, w2, n, n, n, D);
}
else {
               /* Cn = Wn * Pn */
               dmat_copy(C, T_l[i], n, n);
               dmat_potrs(P, C, n, n);

               /* Dn = An * (Un_l * Rn + Rn * Un_l) * Pn */
               dsym_mta(R_l[i], R, n, n, w2);
               dmat_potrs(P, w2, n, n);
               dmat_mul(A, w2, n, n, n, D);
}

               /*---------------------------------------------------------------
                *
                *-------------------------------------------------------------*/
if (S_p_l) {
               /* an = Rn * Snm + Snp * t */
               dm_v_mul(R, S_m, n, n, v1);
               dvec_scale(atran, S_p, v2, n);
               dvec_add(v1, v2, v1, n);

               /* v1 = Cn * a + Dn * a + An * (Un_l * Snm + Rn * Vnm + Vnp * t + Snp * L(t)) */
               dm_v_mul(C, v1, n, n, v2);
               dm_v_mul(D, v1, n, n, v3);
               dvec_add(v2, v3, v1, n);

               dm_v_mul(R_l[i], S_m, n, n, v2);
               dm_v_mul(R, S_m_l[i], n, n, v3);
               dvec_add(v2, v3, v2, n);
               dvec_scale(atran, S_p_l[i], v3, n);
               dvec_add(v2, v3, v2, n);
               dvec_scale(atran_l[i], S_p, v3, n);
               dvec_add(v2, v3, v2, n);

               dm_v_mul(A, v2, n, n, v3);
               dvec_add(v1, v3, v1, n);
}
if (S_m_l) {
               /* bn = Rn * Snm + Snp * t */
               dm_v_mul(R, S_p, n, n, v2);
               dvec_scale(atran, v2, v2, n);
               dvec_add(v2, S_m, v2, n);

               /* v2 = Cn * b + Dn * b + An * [Un_l * Snp * t + Rn * (Vnp * t + Snp * L(t)) + Vnm] */
               dm_v_mul(C, v2, n, n, v3);
               dm_v_mul(D, v2, n, n, v4);
               dvec_add(v3, v4, v2, n);

               dvec_scale(atran, S_p_l[i], v3, n);
               dvec_scale(atran_l[i], S_p, v4, n);
               dvec_add(v3, v4, v3, n);
               dm_v_mul(R, v3, n, n, v4);

               dm_v_mul(R_l[i], S_p, n, n, v3);
               dvec_scale(atran, v3, v3, n);
               dvec_add(v3, v4, v3, n);

               dvec_add(v3, S_m_l[i], v3, n);

               dm_v_mul(A, v3, n, n, v4);
               dvec_add(v2, v4, v2, n);
}
if (S_p_l) {
               /* V(n + 1)p = v1 + Vnp */
               dvec_add(v1, S_p_l[i], S_p_l[i], n);
}
if (S_m_l) {
               /* V(n + 1)m = v2 + Vnm * t + Snm * L(t) */
               dvec_scale(atran, S_m_l[i], v3, n);
               dvec_add(v2, v3, v2, n);
               dvec_scale(atran_l[i], S_m, v3, n);
               dvec_add(v2, v3, S_m_l[i], n);
}
               /* U(n + 1)_l = Un_l + (Cn + Dn) * Bn + An * (Un_l * Tn + Rn * Wn) */
               dmat_add(C, D, w1, n, n);

               dmat_mul(R_l[i], T, n, n, n, w2);
               dmat_mul(R, T_l[i], n, n, n, w3);
               dmat_add(w2, w3, w2, n, n);

               dsym_mma(w1, B, A, w2, n, n, w3, 0);

               dsym_add(R_l[i], w3, R_l[i], n, 1);

               /* W(n + 1) = (Cn + Dn) * Tn + An * Wn */
               dsym_mul(w1, T, n, n, w2, 0);
               dsym_mul(A , T_l[i], n, n, w3, 0);

               dsym_add(w2, w3, T_l[i], n, 1);
          }


          /*----------------------------------------------------------
           *
           *--------------------------------------------------------*/
          else
          if ((S_p_l || S_p_l) && derivs_p[i]) {
if (S_p_l) {
               /* v1 = Cn * a + Dn * a + An * (Un_l * Snm + Rn * Vnm + Vnp * t + Snp * L(t)) */
               /* v1 =                   An * (             Rn * Vnm + Vnp * t + Snp * L(t)) */
               dm_v_mul(R, S_m_l[i], n, n, v1);
               dvec_scale(atran, S_p_l[i], v2, n);
               dvec_add(v1, v2, v1, n);
               dvec_scale(atran_l[i], S_p, v2, n);
               dvec_add(v1, v2, v2, n);

               dm_v_mul(A, v2, n, n, v1);
}
if (S_m_l) {
               /* v2 = Cn * b + Dn * b + An * [Un_l * Snp * t + Rn * (Vnp * t + Snp * L(t)) + Vnm] */
               /* v2 =                   An * [                 Rn * (Vnp * t + Snp * L(t)) + Vnm] */
               dvec_scale(atran, S_p_l[i], v2, n);
               dvec_scale(atran_l[i], S_p, v3, n);
               dvec_add(v2, v3, v2, n);
               dm_v_mul(R, v2, n, n, v3);

               dvec_add(v3, S_m_l[i], v3, n);

               dm_v_mul(A, v3, n, n, v2);
}
if (S_p_l) {
               /* V(n + 1)p = v1 + Vnp */
               dvec_add(v1, S_p_l[i], S_p_l[i], n);
}
if (S_m_l) {
               /* V(n + 1)m = v2 + Vnm * t + Snm * L(t) */
               dvec_scale(atran, S_m_l[i], v3, n);
               dvec_add(v2, v3, v2, n);
               dvec_scale(atran_l[i], S_m, v3, n);
               dvec_add(v2, v3, S_m_l[i], n);
}
          }
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
if (S_p) {
     /* v1 = A * (Rn * Snm + Snp * t) */
     dm_v_mul(R, S_m, n, n, v1);
     dvec_scale(atran, S_p, v2, n);
     dvec_add(v1, v2, v2, n);
     dm_v_mul(A, v2, n, n, v1);
}
if (S_m) {
     /* v2 = A * (Rn * Snp * t + Snm) */
     dm_v_mul(R, S_p, n, n, v2);
     dvec_scale(atran, v2, v2, n);
     dvec_add(v2, S_m, v3, n);
     dm_v_mul(A, v3, n, n, v2);
}
if (S_p) {
     /* Snp = v1 + Snp */
     dvec_add(v1, S_p, S_p, n);
}
if (S_m) {
     /* Snm = v2 + Snm * t */
     dvec_scale(atran, S_m, S_m, n);
     dvec_add(v2, S_m, S_m, n);
}
     /* R(n+1) = Rn + An * Bn */
     dsym_mul(A, B, n, n, w1, 1);
     dmat_add(R, w1, R, n, n);

     /* T(n+1) = An * Tn */
     dsym_mul(A, T, n, n, w1, 1);
     dmat_copy(T, w1, n, n);
}



/*******************************************************************************
 *
 ******************************************************************************/
#ifdef INCLUDE_DEV_SOURCE
#include "../src2/xrtm_doubling2.c"
#endif
