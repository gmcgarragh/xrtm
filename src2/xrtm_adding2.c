/*******************************************************************************
**
**    Copyright (C) 2007-2020 Greg McGarragh <greg.mcgarragh@colostate.edu>
**
**    This source code is licensed under the GNU General Public License (GPL),
**    Version 3.  See the file COPYING for more details.
**
*******************************************************************************/


/*******************************************************************************
 * 13 may be 12 and/or 23, 31 may be 21 and/or 32 (can accumulate a stack or
   double without a copy) at the cost of 2n^2 of memory
 ******************************************************************************/
void layer_add2(double **R12, double **T12, double *Se12, double *Sl12, double *a12,
                double **R21, double **T21, double *Se21, double *Sl21, double *a21,
                double **R23, double **T23, double *Se23, double *Sl23, double *a23,
                double **R32, double **T32, double *Se32, double *Sl32, double *a32,
                double **R13, double **T13, double *Se13, double *Sl13, double *a13,
                double **R31, double **T31, double *Se31, double *Sl31, double *a31,
                int n, double atran, double lin_fac, int flag1, int flag2, int flag3,
                work_data work) {

     int *i1;

     double *v1;
     double *v2;
     double *v3;

     double **w1;
     double **w2;

     double **A31;
     double **A13;

     i1  = get_work1(&work, WORK_IX);

     v1  = get_work1(&work, WORK_DX);
     v2  = get_work1(&work, WORK_DX);
     v3  = get_work1(&work, WORK_DX);

     w1  = get_work1(&work, WORK_DXX);
     w2  = get_work1(&work, WORK_DXX);

     A31 = get_work1(&work, WORK_DXX);
     A13 = get_work1(&work, WORK_DXX);


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     /* LU of (E - R23 * R21)^T */
     dmat_mul(R23, R21, n, n, n, w1);
     dmat_i_sub(w1, w1, n);
     dmat_trans(w1, w2, n, n);
     dmat_getrf(w2, n, n, i1);

     /* A31 = T21 * P31 */
     dmat_copy(A31, T21, n, n);
     dmat_getrs(w2, A31, n, n, i1);


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     /* LU of (E - R21 * R23)^T */
     dmat_mul(R21, R23, n, n, n, w1);
     dmat_i_sub(w1, w1, n);
     dmat_trans(w1, w2, n, n);
     dmat_getrf(w2, n, n, i1);

     /* A13 = T23 * P13 */
     dmat_copy(A13, T23, n, n);
     dmat_getrs(w2, A13, n, n, i1);


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (flag1) {
if (Se31) {
          /* v1 = A31 * (Se32 * t12 + R23 * Se12) */
          dvec_scale(atran, Se32, v1, n);
          dm_v_mul(R23, Se12, n, n, v2);
          dvec_add(v1, v2, v2, n);
          dm_v_mul(A31, v2, n, n, v1);
}
if (Se13) {
          /* v2 = A13 * (Se12 + R21 * Se32 * t12) */
          dvec_scale(atran, Se32, v2, n);
          dm_v_mul(R21, v2, n, n, v3);
          dvec_add(Se12, v3, v2, n);
          dm_v_mul(A13, v2, n, n, v3);
}
if (Se31) {
          /* Se31 = Se21       + v1 */
          dvec_add(Se21, v1, Se31, n);
}
if (Se13) {
          /* Se13 = Se23 * t12 + v2 */
          dvec_scale(atran, Se23, v2, n);
          dvec_add(v2, v3, Se13, n);
}
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (flag2) {
          if (flag3) {
               dvec_copy(a31, Sl31, n);
               dvec_copy(a13, Sl13, n);
          }
if (Se31) {
          /* v1 = A31 * (R23 * S12 + S32 + a31 * f) */
          dm_v_mul(R23, Sl12, n, n, v1);
          dvec_add(v1, Sl32, v1, n);
          dvec_scale(lin_fac, a31, v2, n);
          dvec_add(v1, v2, v2, n);
          dm_v_mul(A31, v2, n, n, v1);
}
if (Sl13) {
          /* v2 = A13 * (R21 * (Sl32 + a32 * f) + Sl12) */
          dvec_scale(lin_fac, a32, v2, n);
          dvec_add(Sl32, v2, v2, n);
          dm_v_mul(R21, v2, n, n, v3);
          dvec_add(v3, Sl12, v3, n);
          dm_v_mul(A13, v3, n, n, v2);
}
if (Sl31) {
          /* S31 = v1 + S21 */
          dvec_add(v1, Sl21, Sl31, n);
}
if (Sl13) {
          /* S13 = v2 + S23 + a23 * f */
          dvec_add(v2, Sl23, Sl13, n);
          dvec_scale(lin_fac, a23, v1, n);
          dvec_add(Sl13, v1, Sl13, n);
}
if (a31) {
          /* v1 = A31 * (R23 * a12 + a32) */
          dm_v_mul(R23, a12, n, n, v1);
          dvec_add(v1, a32, v2, n);
          dm_v_mul(A31, v2, n, n, v1);
}
if (a13) {
          /* v2 = A13 * (R21 * a32 + a12) */
          dm_v_mul(R21, a32, n, n, v2);
          dvec_add(v2, a12, v3, n);
          dm_v_mul(A13, v3, n, n, v2);
}
if (a31) {
          /* a31 = v1 + a21 */
          dvec_add(v1, a21, a31, n);
}
if (a13) {
          /* a13 = v2 + a23 */
          dvec_add(v2, a23, a13, n);
}
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
if (R13) {
     /* R13 = R12 + A31 * B13 */
     dmat_mul(A31, R23, n, n, n, w1);
     dmat_mul(w1, T12, n, n, n, w2);
     dmat_add(R12, w2, R13, n, n);
}
if (R31) {
     /* R31 = R32 + A13 * B31 */
     dmat_mul(A13, R21, n, n, n, w1);
     dmat_mul(w1, T32, n, n, n, w2);
     dmat_add(R32, w2, R31, n, n);
}

if (T31) {
     /* T31 = A31 * T32 */
     dmat_mul(A31, T32, n, n, n, w1);
     dmat_copy(T31, w1, n, n);
}
if (T13) {
     /* T13 = A13 * T12 */
     dmat_mul(A13, T12, n, n, n, w1);
     dmat_copy(T13, w1, n, n);
}
}



/*******************************************************************************
 * 13 may be 12 and/or 23, 31 may be 21 and/or 32 (can accumulate a stack or
   double without a copy) at the cost of 2n^2 of memory
 ******************************************************************************/
void layer_add2_l(double **R12, double **T12, double *S12,
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
                  uchar *derivs, work_data work) {

     int i;

     int *i1;

     double *v1;
     double *v2;
     double *v3;

     double *C31;
     double *C13;

     double **w1;
     double **w2;

     double **P31;
     double **P13;
     double **A31;
     double **A13;
     double **B31;
     double **B13;
     double **D31;
     double **D13;

     i1  = get_work1(&work, WORK_IX);

     C31 = get_work1(&work, WORK_DX);
     C13 = get_work1(&work, WORK_DX);

     v1  = get_work1(&work, WORK_DX);
     v2  = get_work1(&work, WORK_DX);
     v3  = get_work1(&work, WORK_DX);

     w1  = get_work1(&work, WORK_DXX);
     w2  = get_work1(&work, WORK_DXX);

     P31 = get_work1(&work, WORK_DXX);
     P13 = get_work1(&work, WORK_DXX);
     A31 = get_work1(&work, WORK_DXX);
     A13 = get_work1(&work, WORK_DXX);
     B31 = get_work1(&work, WORK_DXX);
     B13 = get_work1(&work, WORK_DXX);

     if (n_derivs > 0) {
          D31 = get_work1(&work, WORK_DXX);
          D13 = get_work1(&work, WORK_DXX);
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     /* w1 = LU of (E - R23 * R21)^T */
     dmat_mul(R23, R21, n, n, n, w1);
     dmat_i_sub(w1, w2, n);
     dmat_trans(w2, P31, n, n);
     dmat_getrf(P31, n, n, i1);

     /* A31 = T21 * P31 */
     dmat_copy(A31, T21, n, n);
     dmat_getrs(P31, A31, n, n, i1);

     /* B13 = R23 * T12 */
     dmat_mul(R23, T12, n, n, n, B13);
if (S31 || S31_l) {
     /* C31 = S32 * t12 + R23 * S12 */
     dvec_scale(atran, S32, v1, n);
     dm_v_mul(R23, S12, n, n, v2);
     dvec_add(v1, v2, C31, n);
}

     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     /* w1 = LU of (E - R21 * R23)^T */
     dmat_mul(R21, R23, n, n, n, w1);
     dmat_i_sub(w1, w2, n);
     dmat_trans(w2, P13, n, n);
     dmat_getrf(P13, n, n, i1);

     /* A13 = T23 * P13 */
     dmat_copy(A13, T23, n, n);
     dmat_getrs(P13, A13, n, n, i1);

     /* B31 = R21 * T32 */
     dmat_mul(R21, T32, n, n, n, B31);
if (S13 || S13_l) {
     /* C13 = S12 + R21 * S32 * t12 */
     dm_v_mul(R21, S32, n, n, v1);
     dvec_scale(atran, v1, v1, n);
     dvec_add(S12, v1, C13, n);
}

     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     for (i = 0; i < n_derivs; ++i) {
          if (! derivs[i])
               continue;

          /* D31 = (T21_l + A31 * (R23_l * R21 + R23 * R21_l)) * P31 */
          dmat_mul(R23_l[i], R21, n, n, n, w1);
          dmat_mul(R23, R21_l[i], n, n, n, w2);
          dmat_add(w1, w2, w1, n, n);
          dmat_mul(A31, w1, n, n, n, w2);

          dmat_add(T21_l[i], w2, D31, n, n);

          dmat_getrs(P31, D31, n, n, i1);


          /*-------------------------------------------------------------------------
           *
           *-----------------------------------------------------------------------*/

          /* D13 = (T23_l + A13 * (R21_l * R23 + R21 * R23_l)) * P13 */
          dmat_mul(R21_l[i], R23, n, n, n, w1);
          dmat_mul(R21, R23_l[i], n, n, n, w2);
          dmat_add(w1, w2, w1, n, n);
          dmat_mul(A13, w1, n, n, n, w2);

          dmat_add(T23_l[i], w2, D13, n, n);

          dmat_getrs(P13, D13, n, n, i1);


          /*-------------------------------------------------------------------------
           *
           *-----------------------------------------------------------------------*/
if (S31_l) {
          /* v1 = D31 * C31 + A31 * (S32_l * t12 + S32 * L(t12) + R23_l * S12 + R23 * S12_l)  */
          dvec_scale(atran, S32_l[i], v1, n);
          dvec_scale(atran_l[i], S32, v2, n);
          dvec_add(v1, v2, v1, n);
          dm_v_mul(R23_l[i], S12, n, n, v2);
          dvec_add(v1, v2, v1, n);
          dm_v_mul(R23, S12_l[i], n, n, v2);
          dvec_add(v1, v2, v1, n);

          dm_v_mul(A31, v1, n, n, v2);
          dm_v_mul(D31, C31, n, n, v1);
          dvec_add(v1, v2, v1, n);
}
if (S13_l) {
          /* v2 = D13 * C13 + A13 * (S12_l + R21_l * S32 * t12 + R21 * (S32_l * t12 + S32 * L(t12))) */
          dvec_scale(atran, S32_l[i], v2, n);
          dvec_scale(atran_l[i], S32, v3, n);
          dvec_add(v2, v3, v2, n);
          dm_v_mul(R21, v2, n, n, v3);
          dvec_add(S12_l[i], v3, v2, n);
          dm_v_mul(R21_l[i], S32, n, n, v3);
          dvec_scale(atran, v3, v3, n);
          dvec_add(v2, v3, v2, n);

          dm_v_mul(A13, v2, n, n, v3);
          dm_v_mul(D13, C13, n, n, v2);
          dvec_add(v2, v3, v2, n);
}
if (S31_l) {
          /* S31_l = S21_l + v1  */
          dvec_add(S21_l[i], v1, S31_l[i], n);
}
if (S13_l) {
          /* S13_l = S23_l * t12 + S23 * L(t12) + v2 */
          dvec_scale(atran, S23_l[i], S13_l[i], n);
          dvec_scale(atran_l[i], S23, v1, n);
          dvec_add(S13_l[i], v1, S13_l[i], n);
          dvec_add(S13_l[i], v2, S13_l[i], n);
}
if (R13_l) {
          /* R13_l = R12_l + D31 * B13 + A31 * (R23_l * T12 + R23 * T12_l) */
          dmat_mul(R23_l[i], T12, n, n, n, w1);
          dmat_mul(R23, T12_l[i], n, n, n, w2);
          dmat_add(w1, w2, w1, n, n);
          dmat_mul(A31, w1, n, n, n, w2);
          dmat_add(R12_l[i], w2, w1, n, n);
          dmat_mul(D31, B13, n, n, n, w2);
          dmat_add(w1, w2, R13_l[i], n, n);
}
if (R31_l) {
          /* R31_l = R32_l + D13 * B31 + A13 * (R21_l * T32 + R21 * T32_l) */
          dmat_mul(R21_l[i], T32, n, n, n, w1);
          dmat_mul(R21, T32_l[i], n, n, n, w2);
          dmat_add(w1, w2, w1, n, n);
          dmat_mul(A13, w1, n, n, n, w2);
          dmat_add(R32_l[i], w2, w1, n, n);
          dmat_mul(D13, B31, n, n, n, w2);
          dmat_add(w1, w2, R31_l[i], n, n);
}

if (T31_l) {
          /* T31_l = D31 * T32 + A31 * T32_l */
          dmat_mul(D31, T32, n, n, n, w1);
          dmat_mul(A31, T32_l[i], n, n, n, w2);
          dmat_add(w1, w2, w1, n, n);
          dmat_copy(T31_l[i], w1, n, n);
}
if (T13_l) {
          /* T13_l = D13 * T12 + A13 * T12_l */
          dmat_mul(D13, T12, n, n, n, w1);
          dmat_mul(A13, T12_l[i], n, n, n, w2);
          dmat_add(w1, w2, w1, n, n);
          dmat_copy(T13_l[i], w1, n, n);
}
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
if (S31) {
     /* S31 = S21       + A31 * C31 */
     dm_v_mul(A31, C31, n, n, v1);
     dvec_add(S21, v1, S31, n);
}
if (R13) {
     /* R13 = R12 + A31 * B13 */
     dmat_mul(A31, B13, n, n, n, w1);
     dmat_add(R12, w1, R13, n, n);
}
if (T31) {
     /* T31 = A31 * T32 */
     dmat_mul(A31, T32, n, n, n, w1);
     dmat_copy(T31, w1, n, n);
}
if (S13) {
     /* S13 = S23 * t12 + A13 * C13 */
     dm_v_mul(A13, C13, n, n, v2);
     dvec_scale(atran, S23, v1, n);
     dvec_add(v1, v2, S13, n);
}
if (R31) {
     /* R31 = R32 + A13 * B31 */
     dmat_mul(A13, B31, n, n, n, w1);
     dmat_add(R32, w1, R31, n, n);
}
if (T13) {
     /* T13 = A13 * T12 */
     dmat_mul(A13, T12, n, n, n, w1);
     dmat_copy(T13, w1, n, n);
}
}
