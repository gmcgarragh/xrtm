/*******************************************************************************
**
**    Copyright (C) 2007-2020 Greg McGarragh <greg.mcgarragh@colostate.edu>
**
**    This source code is licensed under the GNU General Public License (GPL),
**    Version 3.  See the file COPYING for more details.
**
*******************************************************************************/

#include <cublas.h>


/*******************************************************************************
 *
 ******************************************************************************/
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


     int *d_i1;

     double *d_R;
     double *d_T;
     double *d_Se_m;
     double *d_Se_p;

     double *d_v1;
     double *d_v2;
     double *d_v3;

     double *d_w1;
     double *d_w2;

     double *d_A;


     cublasStatus status;


     cublasInit();


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
     status = cublasAlloc(n * n, sizeof(double), (void **) &d_R);
     if (status != CUBLAS_STATUS_SUCCESS) {
          printf ("cublasAlloc()");
          cublasShutdown();
          exit(1);
     }
     status = cublasAlloc(n * n, sizeof(double), (void **) &d_T);
     status = cublasAlloc(n,     sizeof(double), (void **) &d_Se_m);
     status = cublasAlloc(n,     sizeof(double), (void **) &d_Se_p);

     status = cublasAlloc(n,     sizeof(int),   (void **) &d_i1);

     status = cublasAlloc(n,     sizeof(double), (void **) &d_v1);
     status = cublasAlloc(n,     sizeof(double), (void **) &d_v2);
     status = cublasAlloc(n,     sizeof(double), (void **) &d_v3);

     status = cublasAlloc(n * n, sizeof(double), (void **) &d_w1);
     status = cublasAlloc(n * n, sizeof(double), (void **) &d_w2);

     status = cublasAlloc(n * n, sizeof(double), (void **) &d_A);


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
#ifdef CRAP
     /* R(n+1) = Rn + An * Bn */
     dmat_mul(A, R, n, n, n, w1);
     dmat_mul(w1, T, n, n, n, w2);
     dmat_add(R, w2, R, n, n);

     /* T(n+1) = An * Tn */
     dmat_mul(A, T, n, n, n, w1);
     dmat_copy(T, w1, n, n);
#endif

     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     cublasSetVector(n * n, sizeof(double), *R, 1, d_R, 1);
     cublasSetVector(n * n, sizeof(double), *T, 1, d_T, 1);

     cublasSetVector(n * n, sizeof(double), *A, 1, d_A, 1);


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/

     cublasDgemm('n', 'n', n, n, n, 1., d_R, n, d_A, n, 0., d_w1, n);
     cublasDgemm('n', 'n', n, n, n, 1., d_T, n, d_w1, n, 1., d_R, n);

     cublasDgemm('n', 'n', n, n, n, 1., d_T, n, d_A, n, 0., d_w1, n);


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     cublasGetVector(n * n, sizeof(double), d_R, 1, *R, 1);

     cublasGetVector(n * n, sizeof(double), d_w1, 1, *T, 1);


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     cublasFree(d_R);
     cublasFree(d_T);
     cublasFree(d_Se_m);
     cublasFree(d_Se_p);

     cublasFree(d_v1);
     cublasFree(d_v2);
     cublasFree(d_v3);

     cublasFree(d_w1);
     cublasFree(d_w2);

     cublasFree(d_A);


     cublasShutdown();
}
