/*******************************************************************************
**
**    Copyright (C) 2007-2020 Greg McGarragh <greg.mcgarragh@colostate.edu>
**
**    This source code is licensed under the GNU General Public License (GPL),
**    Version 3.  See the file COPYING for more details.
**
*******************************************************************************/


#ifdef __cplusplus
extern "C" {
#endif

void rs_(int *, int *, double *, double *, int *, double *, double *, double *, int *);

void dsyev_(const char *, const char *, int *, double *, int *, double *, double *, int *, int *);

void rs2_(int *, int *, double *, double *, int *, double *, double *, double *, int *, double *, double *, double *);

void dsyev2_(const char *, const char *, int *, double *, int *, double *, double *, int *, int *, double *, double *, double *);

void dchdc_(double *a, int *lda, int *p, double *work, int *jpvt, int *job, int *info);

void dtrsm_(const char *, const char *, const char *, const char *, int *, int *, double *, double *, int *, double *, int *);

#ifdef __cplusplus
}
#endif


/*******************************************************************************
 *
 ******************************************************************************/
static void eig_func_sym_real(int n_quad, double **A, double *evals,
                              double **evecs, int eigen_solver, work_data work) {
#ifdef HAVE_EISPACK_LIBRARY
if (eigen_solver == EIGEN_SOLVER_SYM_REAL_EISPACK) {
     int matz = 1;
     int ierr;

     double *fv1;
     double *fv2;

     double **w1;

     w1  = get_work1(&work, WORK_DXX);

     fv1 = get_work_d1(&work, 4 * n_quad);
     fv2 = get_work_d1(&work, 4 * n_quad);

     rs_(&n_quad, &n_quad, *A, evals, &matz, *w1, fv1, fv2, &ierr);
     if (ierr) {
          fprintf(stderr, "ERROR: rs() info = %d\n", ierr);
          exit(1);
     }

     dmat_trans(w1, evecs, n_quad, n_quad);
}
else
#endif
{
     int lwork = 4 * n_quad;
     int info;

     double *dwork;

     dwork = get_work_d1(&work, 4 * n_quad);

     dsyev_("V", "L", &n_quad, *A, &n_quad, evals, dwork, &lwork, &info);
     if (info) {
          fprintf(stderr, "ERROR: dsyev() info = %d\n", info);
          exit(1);
     }

     dmat_trans(A, evecs, n_quad, n_quad);
}

}



/*******************************************************************************
 *
 ******************************************************************************/
static void eig_func_sym_real2(int n_quad, double **A, double *evals,
                               double **evecs, int eigen_solver, work_data work,
                               double *d, double *e, double **q) {
#ifdef HAVE_EISPACK_LIBRARY
if (eigen_solver == EIGEN_SOLVER_SYM_REAL_EISPACK) {
     int i;

     int matz = 1;
     int ierr;

     double *fv1;
     double *fv2;

     double **w1;

     w1  = get_work1(&work, WORK_DXX);

     fv1 = get_work_d1(&work, 4 * n_quad);
     fv2 = get_work_d1(&work, 4 * n_quad);

     rs2_(&n_quad, &n_quad, *A, evals, &matz, *w1, fv1, fv2, &ierr, d, e, *q);
     if (ierr) {
          fprintf(stderr, "ERROR: rs() info = %d\n", ierr);
          exit(1);
     }

     dmat_trans(w1, evecs, n_quad, n_quad);

     for (i = 0; i < n_quad-1; ++i)
          e[i] = e[i+1];
}
else
#endif
{
     int lwork = 4 * n_quad;
     int info;

     double *dwork;

     dwork = get_work_d1(&work, 4 * n_quad);

     dsyev2_("V", "L", &n_quad, *A, &n_quad, evals, dwork, &lwork, &info, d, e, *q);
     if (info) {
          fprintf(stderr, "ERROR: dsyev() info = %d\n", info);
          exit(1);
     }

     dmat_trans(A, evecs, n_quad, n_quad);
}

}



/*******************************************************************************
 *
 ******************************************************************************/
void eig_1n_sym_real(int n_quad, int n_derivs,
                     double **tpr, double **tmr, double **gamma,
                     double *evals, double **evecs,
                     double **tpr_l, double **tmr_l, double **gamma_l,
                     double *evals_l, double **evecs_l,
                     double **aux, double **aux_l, int eigen_solver,
                     work_data work) {
if (0) {
     int i;
     int j;
     int k;

     int *i1;

     double a;
     double b;

     double *v1;

     double *B;

     double **w1;
     double **w2;
     double **w3;

     double **A;


     w1 = get_work1(&work, WORK_DXX);
     w2 = get_work1(&work, WORK_DXX);
     w3 = get_work1(&work, WORK_DXX);


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     dmat_scale(-1., tmr, w1, n_quad, n_quad);
     dmat_scale(-1., tpr, w2, n_quad, n_quad);


if (0) {
     dmat_chol(w2, n_quad);

     for (i = 0; i < n_quad; ++i) {
          for (j = i + 1; j < n_quad; ++j) {
               w2[i][j] = w2[j][i];
          }
     }
}
else
if (0) {
     int job = 0;
     int info;

     int *jpvt;

     double **work2;

     jpvt  = get_work1(&work, WORK_IX);

     work2 = get_work1(&work, WORK_DXX);

     dchdc_(*w2, &n_quad, &n_quad, *work2, jpvt, &job, &info);

     for (i = 0; i < n_quad; ++i) {
          for (j = i + 1; j < n_quad; ++j) {
               w2[i][j] = w2[j][i];
          }
     }
}
else
if (1) {
     dmat_potrf(w2, n_quad);

     for (i = 0; i < n_quad; ++i) {
          for (j = i + 1; j < n_quad; ++j) {
               w2[j][i] = w2[i][j];
          }
     }
}

/*
     dmat_dgxtlmx(w1, w2, w3, n_quad);
     dmat_dtugxmx(w2, w3, w1, n_quad);
*/
     dmat_dtugxtlms(w2, w1, w2, w1, w3, n_quad);


     eig_func_sym_real(n_quad, w1, evals, evecs, eigen_solver, work);

     if (aux) {
          dmat_dtlgxmx(w2, evecs, aux, n_quad);
          dmat_scale(-1., aux, aux, n_quad, n_quad);
     }

{
     double alpha = 1.;

     dtrsm_("r", "l", "n", "n",
            &n_quad, &n_quad, &alpha, *w2, &n_quad, *evecs, &n_quad);
}

     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (n_derivs > 0) {
          i1 = get_work1(&work, WORK_IX2);

          v1 = get_work1(&work, WORK_DX);

          B  = alloc_array1_d(n_quad + 1);
          A  = alloc_array2_d(n_quad + 1, n_quad + 1);


          for (i = 0; i < n_quad; ++i) {
               a = sqrt(evals[i]) * 2.;
               b = evals[i];

               for (j = 0; j < n_quad; ++j) {
                    A[j][0] = a * evecs[j][i];

                    for (k = 0; k < n_quad; ++k)
                         A[j][k+1] = -gamma[j][k];

                    A[j][j+1] += b;
               }
               A[n_quad][0] = 0.;

               for (j = 0; j < n_quad; ++j)
                    A[n_quad][j+1] = evecs[j][i];


               for (j = 0; j < n_quad; ++j)
                    v1[j] = evecs[j][i];
               dm_v_mul(gamma_l, v1, n_quad, n_quad, B);

               B[n_quad] = 0.;


               dmat_getrf(A, n_quad + 1, n_quad + 1, i1);

               dmat_getrs(A, &B, n_quad + 1, 1, i1);


               evals_l[i] = B[0];

               for (j = 0; j < n_quad; ++j)
                    evecs_l[j][i] = B[j+1];
          }

          free_array1_d(B);
          free_array2_d(A);
     }
}
else {
     int i;
     int j;

     int *i1;

     double a;

     double *v1;
     double *v2;
     double *v3;

     double **w1;
     double **w2;
     double **w3;
     double **w4;
     double **w5;
     double **w6;
     double **w7;
/*
     double **w8;
*/
     i1 = get_work1(&work, WORK_IX);

     v1 = get_work1(&work, WORK_DX);
     v2 = get_work1(&work, WORK_DX);
     v3 = get_work1(&work, WORK_DX);

     w1 = get_work1(&work, WORK_DXX);
     w2 = get_work1(&work, WORK_DXX);
     w3 = get_work1(&work, WORK_DXX);


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (n_derivs == 0) {
          dmat_scale(-1., tmr, w1, n_quad, n_quad);
          dmat_scale(-1., tpr, w2, n_quad, n_quad);

          dmat_chol(w2, n_quad);

          for (i = 0; i < n_quad; ++i) {
               for (j = i + 1; j < n_quad; ++j) {
                    w2[i][j] = w2[j][i];
               }
          }
/*
          dmat_dgxtlmx(w1, w2, w3, n_quad);
          dmat_dtugxmx(w2, w3, w1, n_quad);
*/
          dmat_dtugxtlms(w2, w1, w2, w1, w3, n_quad);

          for (i = 0; i < n_quad; ++i) {
               for (j = i + 1; j < n_quad; ++j) {
                    w1[j][i] = w1[i][j];
               }
          }

          eig_func_sym_real(n_quad, w1, evals, evecs, eigen_solver, work);

          if (aux) {
               dmat_dtlgxmx(w2, evecs, aux, n_quad);
               dmat_scale(-1., aux, aux, n_quad, n_quad);
          }

{
          double alpha = 1.;

          dtrsm_("r", "l", "n", "n",
                 &n_quad, &n_quad, &alpha, *w2, &n_quad, *evecs, &n_quad);
}
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     else {
          w4 = get_work1(&work, WORK_DXX);
          w5 = get_work1(&work, WORK_DXX);
          w6 = get_work1(&work, WORK_DXX);
          w7 = get_work1(&work, WORK_DXX);
/*
          w8 = get_work1(&work, WORK_DXX);
*/

          dmat_scale(-1., tmr, w1, n_quad, n_quad);
          dmat_scale(-1., tpr, w2, n_quad, n_quad);
          dmat_scale(-1., tmr_l, w3, n_quad, n_quad);
          dmat_scale(-1., tpr_l, w4, n_quad, n_quad);


          dmat_chol_l(w2, w4, n_quad);

          for (i = 0; i < n_quad; ++i) {
               for (j = i + 1; j < n_quad; ++j) {
                    w2[i][j] = w2[j][i];
                    w4[i][j] = w4[j][i];
               }
          }
/*
          dmat_dtugxtlms(w2, w1, w2, w5, w7, n_quad);

          for (i = 0; i < n_quad; ++i) {
               for (j = i + 1; j < n_quad; ++j) {
                    w5[j][i] = w5[i][j];
               }
          }
*/
/*
          dmat_dgxtlmx(w1, w2, w8, n_quad);
          dmat_dtugxmx(w4, w8, w6, n_quad);

          dmat_dgxtlmx(w3, w2, w8, n_quad);
          dmat_dtugxmx(w2, w8, w7, n_quad);

          dmat_add(w6, w7, w6, n_quad, n_quad);

          dmat_dgxtlmx(w1, w4, w8, n_quad);
          dmat_dtugxmx(w2, w8, w7, n_quad);

          dmat_add(w6, w7, w6, n_quad, n_quad);
*/
/*
          dmat_dtugxtlms(w4, w1, w2, w6, w8, n_quad);
          dmat_dtugxtlms(w2, w3, w2, w7, w8, n_quad);
          dmat_add(w6, w7, w6, n_quad, n_quad);
          dmat_dtugxtlms(w2, w1, w4, w7, w8, n_quad);
          dmat_add(w6, w7, w6, n_quad, n_quad);
*/
/*
          for (i = 0; i < n_quad; ++i) {
               for (j = i + 1; j < n_quad; ++j) {
                    w6[j][i] = w6[i][j];
               }
          }
*/
          dmat_dtugxtlms_l(w2, w1, w2, w5, w4, w3, w4, w6, w7, n_quad);

          for (i = 0; i < n_quad; ++i) {
               for (j = i + 1; j < n_quad; ++j) {
                    w5[j][i] = w5[i][j];
                    w6[j][i] = w6[i][j];
               }
          }

          eig_func_sym_real2(n_quad, w5, evals, evecs, eigen_solver, work, v1, v2, w1);

          if (aux) {
               dmat_dtlgxmx(w2, evecs, aux, n_quad);
               dmat_scale(-1., aux, aux, n_quad, n_quad);
          }
if (0) {
          for (i = 0; i < n_quad; ++i) {
               a = 2. * sqrt(evals[i]);

               for (j = 0; j < n_quad; ++j)
                    v1[j] = evecs[j][i];

               dm_v_mul(w6, v1, n_quad, n_quad, v2);

               evals_l[i] = dvec_dot(v1, v2, n_quad) / a;


               for (j = 0; j < n_quad; ++j)
                    v3[j] = evals[i];

               dmat_diag_sub(v3, w5, w7, n_quad);
               dmat_getrf(w7, n_quad, n_quad, i1);

               dvec_scale(a * evals_l[i], v1, v3, n_quad);
               dvec_sub(v2, v3, v2, n_quad);

               dmat_getrs(w7, &v2, n_quad, 1, i1);

               for (j = 0; j < n_quad; ++j)
                    evecs_l[j][i] = v2[j];
          }
}
else {
          double *v4;
          double *v5;
          double *v6;
          double *v7;
          double *v8;
          double *v9;

          v4 = get_work1(&work, WORK_DX);
          v5 = get_work1(&work, WORK_DX);
          v6 = get_work1(&work, WORK_DX);
          v7 = get_work1(&work, WORK_DX);
          v8 = get_work1(&work, WORK_DX);
          v9 = get_work1(&work, WORK_DX);

          for (i = 0; i < n_quad; ++i) {
               a = 2. * sqrt(evals[i]);

               for (j = 0; j < n_quad; ++j)
                    v3[j] = evecs[j][i];

               dm_v_mul(w6, v3, n_quad, n_quad, v4);

               evals_l[i] = dvec_dot(v3, v4, n_quad) / a;


               for (j = 0; j < n_quad  ; ++j)
                    v6[j] = evals[i] - v1[j];

               for (j = 0; j < n_quad-1; ++j) {
                    v7[j] = -v2[j];
                    v8[j] = -v2[j];
               }

               dmat_gttrf(v7, v6, v8, v9, n_quad, i1);

               dvec_scale(a * evals_l[i], v3, v5, n_quad);
               dvec_sub(v4, v5, v4, n_quad);

               dm_v_mul(w1, v4, n_quad, n_quad, v5);

               dmat_gttrs(v7, v6, v8, v9, &v5, n_quad, 1, i1);

               dmat_trans(w1, w3, n_quad, n_quad);
               dm_v_mul(w3, v5, n_quad, n_quad, v4);

               for (j = 0; j < n_quad; ++j)
                    evecs_l[j][i] = v4[j];
          }
}
          if (aux_l) {
               dmat_dtlgxmx(w2, evecs_l, w1, n_quad);
               dmat_dtlgxmx(w4, evecs, w3, n_quad);
               dmat_add(w1, w3, aux_l, n_quad, n_quad);
               dmat_scale(-1., aux_l, aux_l, n_quad, n_quad);
          }

{
          double alpha = 1.;

          dtrsm_("r", "l", "n", "n",
                 &n_quad, &n_quad, &alpha, *w2, &n_quad, *evecs, &n_quad);
          dmat_dtugxmx(w4, evecs, w7, n_quad);
          dmat_sub(evecs_l, w7, evecs_l, n_quad, n_quad);
          dtrsm_("r", "l", "n", "n",
                 &n_quad, &n_quad, &alpha, *w2, &n_quad, *evecs_l, &n_quad);
}
     }
}

}




/*******************************************************************************
 *
 ******************************************************************************/
static void eig_func_gen_complex2(int n_quad, double **A, double *evals_r, double *evals_i,
                                  double **l_evecs, double **r_evecs, work_data work) {

     int lwork = 4 * n_quad;
     int info;

     double *ework;

     double **w1;
     double **w2;
     double **w3;

     w1 = get_work1(&work, WORK_DXX);
     w2 = get_work1(&work, WORK_DXX);
     w3 = get_work1(&work, WORK_DXX);

     ework = get_work_d1(&work, 4 * n_quad);

     dmat_trans(A, w1, n_quad, n_quad);

     dgeev_("V", "V", &n_quad, *w1, &n_quad, evals_r, evals_i,
            *w2, &n_quad, *w3, &n_quad, ework, &lwork, &info);
     if (info) {
          fprintf(stderr, "ERROR: dgeev() info = %d\n", info);
          exit(1);
     }

     dmat_trans(w2, l_evecs, n_quad, n_quad);
     dmat_trans(w3, r_evecs, n_quad, n_quad);
}



/*******************************************************************************
 *
 ******************************************************************************/
void eig_1n_gen_complex2(int n_quad, int n_derivs,
                         double  **gamma,
                         double  *evals_r, double  *evals_i, double  **evecs,
                         double ***gamma_l,
                         double **evals_r_l, double **evals_i_l, double ***evecs_l,
                         uchar *derivs, work_data work) {
     int i;
     int i1;
     int i2;
     int j;
     int k;

     int *ip;

     double ar;
     double br;
     double cr;
     double dr;

     double *v1;
     double *v2;
     double *v3;
     double *v4;

     double **w1;
     double **w2;

     double **l_evecs;
     double **r_evecs;

     dcomplex ac;
     dcomplex bc;
     dcomplex cc;
     dcomplex dc;

     dcomplex *vc1;
     dcomplex *vc2;
     dcomplex *vc3;
     dcomplex *vc4;

     dcomplex **wc1;
     dcomplex **wc2;


     l_evecs = get_work1(&work, WORK_DXX);
     r_evecs = evecs;


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     eig_func_gen_complex2(n_quad, gamma, evals_r, evals_i, l_evecs, r_evecs, work);


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (flags_or(derivs, n_derivs)) {
          ip = get_work1(&work, WORK_IX);

          v1 = get_work1(&work, WORK_DX);
          v2 = get_work1(&work, WORK_DX);
          v3 = get_work1(&work, WORK_DX);
          v4 = get_work1(&work, WORK_DX);

          w1 = get_work1(&work, WORK_DXX);
          w2 = get_work1(&work, WORK_DXX);

          vc1 = get_work1(&work, WORK_ZX);
          vc2 = get_work1(&work, WORK_ZX);
          vc3 = get_work1(&work, WORK_ZX);
          vc4 = get_work1(&work, WORK_ZX);

          wc1 = get_work1(&work, WORK_ZXX);
          wc2 = get_work1(&work, WORK_ZXX);

          for (i = 0; i < n_quad; ++i) {
               if (evals_i[i] == 0.) {
                    if (is_degenerate_eigenvalue(i, n_quad, evals_r)) {
                         zero_eigen_derivs(i, n_quad, n_derivs, evals_r_l, evals_i_l, evecs_l, derivs);
                         continue;
                    }

                    ar = sqrt(evals_r[i]) * 2.;
                    br = evals_r[i];

                    for (j = 0; j < n_quad; ++j)
                         v1[j] = br;

                    dmat_diag_sub(v1, gamma, w1, n_quad);
/*
                    perturb_zero_elems(n_quad + 1, w1);
*/
                    dmat_getrf(w1, n_quad, n_quad, ip);

                    for (j = 0; j < n_quad; ++j) {
                         v1[j] = l_evecs[j][i];
                         v2[j] = r_evecs[j][i];
                    }

                    cr = dvec_dot(v1, v2, n_quad);

                    dv_m_mul(v2, &v1, n_quad, w2);
                    dmat_scale(1. / cr, w2, w2, n_quad, n_quad);
                    dmat_i_sub(w2, w2, n_quad);

                    for (j = 0; j < n_derivs; ++j) {
                         if (! derivs[j])
                              continue;

                         dm_v_mul(gamma_l[j], v2, n_quad, n_quad, v3);

                         dr = dvec_dot(v1, v3, n_quad);

                         evals_r_l[j][i] = dr / (ar * cr);
                         evals_i_l[j][i] = 0.;

                         dm_v_mul(w2, v3, n_quad, n_quad, v4);

                         dmat_getrs(w1, &v4, n_quad, 1, ip);

                         for (k = 0; k < n_quad; ++k)
                              evecs_l[j][k][i] = v4[k];
                    }
               }
               else
               if (evals_i[i] > 0.) {
                    i1 = i;
                    i2 = i + 1;

                    bc = evals_r[i] + _Complex_I * evals_i[i];
                    ac = csqrt(bc) * 2.;

                    for (j = 0; j < n_quad; ++j) {
                         for (k = 0; k < n_quad; ++k) {
                              wc1[j][k] = -gamma[j][k];
                         }
                         wc1[j][j] += bc;
                    }

                    zmat_getrf(wc1, n_quad, n_quad, ip);

                    for (j = 0; j < n_quad; ++j) {
                         v1[j]  = r_evecs[j][i1];
                         v2[j]  = r_evecs[j][i2];

                         vc1[j] = l_evecs[j][i1] - _Complex_I * l_evecs[j][i2];
                         vc2[j] = r_evecs[j][i1] + _Complex_I * r_evecs[j][i2];
                    }

                    cc = zvec_dot(vc1, vc2, n_quad);

                    zv_m_mul(vc2, &vc1, n_quad, wc2);
                    zmat_scale(1. / cc, wc2, wc2, n_quad, n_quad);
                    zmat_i_sub(wc2, wc2, n_quad);

                    for (j = 0; j < n_derivs; ++j) {
                         if (! derivs[j])
                              continue;

                         dm_v_mul(gamma_l[j], v1, n_quad, n_quad, v3);
                         dm_v_mul(gamma_l[j], v2, n_quad, n_quad, v4);

                         for (k = 0; k < n_quad; ++k)
                              vc3[k] = v3[k] + _Complex_I * v4[k];

                         dc = zvec_dot(vc1, vc3, n_quad);

                         evals_r_l[j][i1] =  creal(dc / (ac * cc));
                         evals_i_l[j][i1] =  cimag(dc / (ac * cc));

                         evals_r_l[j][i2] =  creal(dc / (ac * cc));
                         evals_i_l[j][i2] = -cimag(dc / (ac * cc));

                         zm_v_mul(wc2, vc3, n_quad, n_quad, vc4);

                         zmat_getrs(wc1, &vc4, n_quad, 1, ip);

                         for (k = 0; k < n_quad; ++k) {
                              evecs_l[j][k][i1] = creal(vc4[k]);
                              evecs_l[j][k][i2] = cimag(vc4[k]);
                         }
                    }
               }
          }
     }
}



/*******************************************************************************
 *
 ******************************************************************************/
void eig_2n_sym_real(int n_quad, int n_derivs,
                     double  **tpr, double  **tmr, double  **gamma,
                     double  *nu, double  **X_p, double  **X_m,
                     double ***tpr_l, double ***tmr_l, double ***gamma_l,
                     double **nu_l, double ***X_p_l, double ***X_m_l,
                     double **aux, double ***aux_l, int eigen_solver,
                     uchar *derivs, save_tree_data save_tree, work_data work) {

     double *evals;

     double **evals_l;

     double **evecs;

     double ***evecs_l;

     evals = get_work1(&work, WORK_DX);

     evecs = get_work1(&work, WORK_DXX);

     if (flags_or(derivs, n_derivs)) {
          evals_l = get_work2(&work, WORK_DX, WORK_DERIVS_V, derivs);

          evecs_l = get_work2(&work, WORK_DXX, WORK_DERIVS_V, derivs);
     }

     eig_1n_sym_real(n_quad, n_derivs, tpr, tmr, NULL, evals, evecs, tpr_l[0], tmr_l[0], gamma_l[0], evals_l[0], evecs_l[0], aux, aux_l[0], eigen_solver, work);

     eig_1n_to_2n_real(n_quad, n_derivs, tpr, tmr, evals, evecs, nu, X_p, X_m, tpr_l, tmr_l, evals_l, evecs_l, nu_l, X_p_l, X_m_l, aux, aux_l, derivs, save_tree, work);
}
