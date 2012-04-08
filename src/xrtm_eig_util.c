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
#include "xrtm_eig_util.h"
#include "xrtm_eig_util_a.h"
#include "xrtm_matrix.h"
#include "xrtm_save_tree.h"
#include "xrtm_utility.h"


#ifdef __cplusplus
extern "C" {
#endif

void dasymtx_(double *, double *, double *, int *, int *, int *, int *, double *);

void rg_(int *, int *, double *, double *, double *, int *, double *, int *, double *, int *);

void dgeev_(const char *, const char *, int *, double *, int *, double *, double *, double *, int *, double *, int *, double *, int *, int *);

#ifdef __cplusplus
}
#endif


/*******************************************************************************
 *
 ******************************************************************************/
int is_dummy_node_eigenvector(int i, int n_quad, double **evecs) {

     int j;

     int count1 = 0;
     int count2 = 0;

     for (j = 0; j < n_quad; ++j) {
          count1 += fabs(1. - evecs[j][i]) < DBL_EPSILON ? 1 : 0;
          count2 += fabs(     evecs[j][i]) < DBL_EPSILON ? 1 : 0;
     }

     return count1 == 1 && count2 == n_quad - 1;
}



void fill_dummy_node_eigen_derivs(int i, int n_quad, int n_derivs, double **evals_r_l, double **evals_i_l, double ***evecs_l, uchar *derivs) {

     int j;
     int k;

     for (j = 0; j < n_derivs; ++j) {
          if (! derivs[j])
               continue;

          evals_r_l[j][i] = 0.;
          if (evals_i_l)
               evals_i_l[j][i] = 0.;

          for (k = 0; k < n_quad; ++k)
               evecs_l[j][k][i] = 0.;
     }
}



void perturb_zero_elems(int n, double **A) {

     int i;
     int j;

     for (i = 0; i < n; ++i) {
          for (j = 0; j < n; ++j) {
               if (A[i][j] == 0.)
                    A[i][j] = DBL_EPSILON * (i + j + 1);
          }
     }
}



/*******************************************************************************
 *
 ******************************************************************************/
static void eig_func_gen_real(int n_quad, double **A, double *evals,
                              double **evecs, int eigen_solver, work_data work) {
#ifdef HAVE_FORTRAN_COMPILER
     if (eigen_solver == EIGEN_SOLVER_GEN_REAL_ASYMTX) {
          int i;
          int j;

          int ier;

          double a;

          double **w1;

          double *wkd;

          w1  = get_work1(&work, WORK_DXX);

          wkd = get_work1(&work, WORK_DX2);

          dmat_trans(A, evecs, n_quad, n_quad);

          dasymtx_(*evecs, *w1, evals, &n_quad, &n_quad, &n_quad, &ier, wkd);
          if (ier) {
               eprintf("ERROR: dasymtx() info = %d\n", ier);
               exit(1);
          }

          dmat_trans(w1, evecs, n_quad, n_quad);

          for (j = 0; j < n_quad; ++j) {
               a = 0.;
               for (i = 0; i < n_quad; ++i) {
                    a += evecs[i][j] * evecs[i][j];
               }
               a = sqrt(a);

               for (i = 0; i < n_quad; ++i) {
                    evecs[i][j] /= a;
               }
          }
     }
     else
#endif
#ifdef HAVE_EISPACK_LIBRARY
     if (eigen_solver == EIGEN_SOLVER_GEN_REAL_EISPACK) {
          int matz = 1;
          int ierr;

          int *iv1;

          double *v1;

          double *fv1;

          double **w1;

          v1  = get_work1(&work, WORK_DX);

          w1  = get_work1(&work, WORK_DXX);

          iv1 = get_work1(&work, WORK_IX);

          fv1 = get_work1(&work, WORK_DX);

          dmat_trans(A, evecs, n_quad, n_quad);

          rg_(&n_quad, &n_quad, *evecs, evals, v1, &matz, *w1, iv1, fv1, &ierr);
          if (ierr) {
               eprintf("ERROR: rg() info = %d\n", ierr);
               exit(1);
          }

          dmat_trans(w1, evecs, n_quad, n_quad);
     }
     else
#endif
     {
          int lwork = 4 * n_quad;
          int info;

          double *v1;

          double **w1;

          double *dwork;

          v1    = get_work1(&work, WORK_DX);

          w1    = get_work1(&work, WORK_DXX);

          dwork = get_work_d1(&work, 4 * n_quad);

          dmat_trans(A, evecs, n_quad, n_quad);

          dgeev_("N", "V", &n_quad, *evecs, &n_quad, evals, v1, NULL, &n_quad, *w1,
                 &n_quad, dwork, &lwork, &info);
          if (info) {
               eprintf("ERROR: dgeev() info = %d\n", info);
               exit(1);
          }

          dmat_trans(w1, evecs, n_quad, n_quad);
     }
}



/*******************************************************************************
 *
 ******************************************************************************/
static void eigen_derivatives_real(int i, int n_quad, int n_derivs, double **gamma, double *evals_r, double *evals_i, double **evecs, double ***gamma_l, double **evals_r_l, double **evals_i_l, double ***evecs_l, uchar *derivs, double **A, double *B, int *ip, double *v1) {

     int j;
     int k;

     double a;
     double b;

     if (is_dummy_node_eigenvector(i, n_quad, evecs)) {
          fill_dummy_node_eigen_derivs(i, n_quad, n_derivs, evals_r_l, evals_i_l, evecs_l, derivs);
          return;
     }

     a = sqrt(evals_r[i]) * 2.;
     b = evals_r[i];

     for (j = 0; j < n_quad; ++j) {
          A[j][0] = a * evecs[j][i];

          for (k = 0; k < n_quad; ++k)
               A[j][k+1] = -gamma[j][k];

          A[j][j+1] += b;
     }
     A[n_quad][0] = 0.;

     for (j = 0; j < n_quad; ++j)
          A[n_quad][j+1] = evecs[j][i];

     if (evals_i_l)
          perturb_zero_elems(n_quad + 1, A);

     dmat_getrf(A, n_quad + 1, n_quad + 1, ip);

     for (j = 0; j < n_quad; ++j)
          v1[j] = evecs[j][i];

     for (j = 0; j < n_derivs; ++j) {
          if (! derivs[j])
               continue;

          dm_v_mul(gamma_l[j], v1, n_quad, n_quad, B);

          B[n_quad] = 0.;

          dmat_getrs(A, &B, n_quad + 1, 1, ip);

          evals_r_l[j][i] = B[0];

          if (evals_i_l)
               evals_i_l[j][i] = 0.;

          for (k = 0; k < n_quad; ++k)
               evecs_l[j][k][i] = B[k+1];
     }
}



/*******************************************************************************
 *
 ******************************************************************************/
static void eig_1n_gen_real(int n_quad, int n_derivs,
                     double  **gamma, double  *evals, double  **evecs,
                     double ***gamma_l, double **evals_l, double ***evecs_l,
                     int eigen_solver, uchar *derivs, save_tree_data save_tree,
                     work_data work) {

     int i;

     int *i1;

     double *v1;

     double *B;

     double **A;

     forward_save_eig_1n_gen_real_data *save;

     if (save_tree.t) {
          save_tree_encode_s(&save_tree, "eig_1n_gen_real");

          if (save_tree_retrieve_data(&save_tree, forward_save_eig_1n_gen_real_data, &save))
               forward_save_eig_1n_gen_real_alloc(save, n_quad);

          save->flag = 0;
     }

     eig_func_gen_real(n_quad, gamma, evals, evecs, eigen_solver, work);

     if (flags_or(derivs, n_derivs)) {
          i1 = get_work_i1(&work, n_quad + 1);

          v1 = get_work1(&work, WORK_DX);

          B  = get_work_d1(&work, n_quad + 1);
          A  = get_work_d2(&work, n_quad + 1, n_quad + 1);

          for (i = 0; i < n_quad; ++i) {
               eigen_derivatives_real(i, n_quad, n_derivs, gamma, evals, NULL, evecs, gamma_l, evals_l, NULL, evecs_l, derivs, A, B, i1, v1);

               if (save_tree.t) {
                    copy_array1_i(save->ip[i], i1, n_quad + 1);

                    dmat_copy(save->A[i], A, n_quad + 1, n_quad + 1);
               }
          }
     }
#ifdef USE_AD_FOR_TL_CALC_EIG_1N_GEN_REAL
     eig_1n_gen_real_tl_with_ad(n_quad, n_derivs, gamma, evals, evecs, gamma_l, evals_l, evecs_l, eigen_solver, derivs, save_tree, work);
#endif
}



/*******************************************************************************
 *
 ******************************************************************************/
static void eig_func_gen_complex(int n_quad, double **A, double *evals_r, double *evals_i,
                                 double **evecs, int eigen_solver, work_data work) {
#ifdef HAVE_EISPACK_LIBRARY
     if (eigen_solver == EIGEN_SOLVER_GEN_COMPLEX_EISPACK) {
          int matz = 1;
          int ierr;

          int *iv1;

          double *fv1;

          double **w1;
          double **w2;

          w1  = get_work1(&work, WORK_DXX);
          w2  = get_work1(&work, WORK_DXX);

          iv1 = get_work1(&work, WORK_IX);

          fv1 = get_work1(&work, WORK_DX);

          dmat_trans(A, w1, n_quad, n_quad);

          rg_(&n_quad, &n_quad, *w1, evals_r, evals_i, &matz, *w2, iv1, fv1, &ierr);
          if (ierr) {
               eprintf("ERROR: rg() info = %d\n", ierr);
               exit(1);
          }

          dmat_trans(w2, evecs, n_quad, n_quad);
     }
     else
#endif
{
          int lwork = 4 * n_quad;
          int info;

          double *ework;

          double **w1;
          double **w2;

          w1 = get_work1(&work, WORK_DXX);
          w2 = get_work1(&work, WORK_DXX);

          ework = get_work_d1(&work, 4 * n_quad);

          dmat_trans(A, w1, n_quad, n_quad);

          dgeev_("N", "V", &n_quad, *w1, &n_quad, evals_r, evals_i, NULL, &n_quad,
                 *w2, &n_quad, ework, &lwork, &info);
          if (info) {
               eprintf("ERROR: dgeev() info = %d\n", info);
               exit(1);
          }

          dmat_trans(w2, evecs, n_quad, n_quad);
     }
}



/*******************************************************************************
 *
 ******************************************************************************/
void eig_1n_gen_complex(int n_quad, int n_derivs,
                        double  **gamma,
                        double  *evals_r, double  *evals_i, double  **evecs,
                        double ***gamma_l,
                        double **evals_r_l, double **evals_i_l, double ***evecs_l,
                        int eigen_solver, uchar *derivs, save_tree_data save_tree,
                        work_data work) {

     int i;
     int i1;
     int i2;
     int j;
     int k;

     int *ip;

     double *v1;
     double *v2;
     double *v3;
     double *v4;

     double *Br;

     double **Ar;

     dcomplex ac;
     dcomplex bc;

     dcomplex *Bc;

     dcomplex **Ac;


     forward_save_eig_1n_gen_complex_data *save;

     if (save_tree.t) {
          save_tree_encode_s(&save_tree, "eig_1n_gen_complex");

          if (save_tree_retrieve_data(&save_tree, forward_save_eig_1n_gen_complex_data, &save))
               forward_save_eig_1n_gen_complex_alloc(save, n_quad);

          save->flag = 0;
     }


     eig_func_gen_complex(n_quad, gamma, evals_r, evals_i, evecs, eigen_solver, work);


     if (flags_or(derivs, n_derivs)) {
          ip = get_work_i1(&work, n_quad + 1);

          v1 = get_work1(&work, WORK_DX);
          v2 = get_work1(&work, WORK_DX);
          v3 = get_work1(&work, WORK_DX);
          v4 = get_work1(&work, WORK_DX);

          Br = get_work_d1(&work, n_quad + 1);
          Ar = get_work_d2(&work, n_quad + 1, n_quad + 1);

          Bc = get_work_dc1(&work, n_quad + 1);
          Ac = get_work_dc2(&work, n_quad + 1, n_quad + 1);

          for (i = 0; i < n_quad; ++i) {
               if (evals_i[i] == 0.) {
                    eigen_derivatives_real(i, n_quad, n_derivs, gamma, evals_r, evals_i, evecs, gamma_l, evals_r_l, evals_i_l, evecs_l, derivs, Ar, Br, ip, v1);

                    if (save_tree.t) {
                         copy_array1_i(save->ip[i], ip, n_quad + 1);

                         dmat_copy(save->A[i], Ar, n_quad + 1, n_quad + 1);
                    }
               }
               else 
               if (evals_i[i] > 0.) {
                    i1 = i;
                    i2 = i + 1;

                    bc = evals_r[i] + _Complex_I * evals_i[i];
                    ac = csqrt(bc) * 2.;

                    for (j = 0; j < n_quad; ++j) {
                         Ac[j][0] = ac * (evecs[j][i1] + _Complex_I * evecs[j][i2]);

                         for (k = 0; k < n_quad; ++k)
                              Ac[j][k+1] = -gamma[j][k];

                         Ac[j][j+1] += bc;
                    }
                    Ac[n_quad][0] = 0.;

                    for (j = 0; j < n_quad; ++j)
                         Ac[n_quad][j+1] = (evecs[j][i1] + _Complex_I * evecs[j][i2]);

                    zmat_getrf(Ac, n_quad + 1, n_quad + 1, ip);

                    for (j = 0; j < n_quad; ++j) {
                         v1[j] = evecs[j][i1];
                         v2[j] = evecs[j][i2];
                    }

                    for (j = 0; j < n_derivs; ++j) {
                         if (! derivs[j])
                              continue;

                         dm_v_mul(gamma_l[j], v1, n_quad, n_quad, v3);
                         dm_v_mul(gamma_l[j], v2, n_quad, n_quad, v4);

                         for (k = 0; k < n_quad; ++k)
                              Bc[k] = v3[k] + _Complex_I * v4[k];

                         Bc[n_quad] = 0.;

                         zmat_getrs(Ac, &Bc, n_quad + 1, 1, ip);

                         evals_r_l[j][i1] =  creal(Bc[0]);
                         evals_i_l[j][i1] =  cimag(Bc[0]);

                         evals_r_l[j][i2] =  creal(Bc[0]);
                         evals_i_l[j][i2] = -cimag(Bc[0]);

                         for (k = 0; k < n_quad; ++k) {
                              evecs_l[j][k][i1] = creal(Bc[k+1]);
                              evecs_l[j][k][i2] = cimag(Bc[k+1]);
                         }
                    }
               }
          }
     }
#ifdef USE_AD_FOR_TL_CALC_EIG_1N_GEN_COMPLEX
     eig_1n_gen_complex_tl_with_ad(n_quad, n_derivs, gamma, evals_r, evals_i, evecs, gamma_l, evals_r_l, evals_i_l, evecs_l, eigen_solver, derivs, save_tree, work);
#endif
}



/*******************************************************************************
 *
 ******************************************************************************/
void eig_1n_to_2n_real(int n_quad, int n_derivs,
                       double  **tpr, double  **tmr,
                       double  *evals, double  **evecs,
                       double  *nu, double  **X_p, double  **X_m,
                       double ***tpr_l, double ***tmr_l,
                       double **evals_l, double ***evecs_l,
                       double **nu_l, double ***X_p_l, double ***X_m_l,
                       double **aux, double ***aul,
                       uchar *derivs, save_tree_data save_tree, work_data work) {

     int i;
     int j;
     int k;

     double a;
     double b;
     double c;

     double **w1;
     double **w2;
     double **w3;

     forward_save_eig_1n_to_2n_real_data *save;
if (0) {
     if (! aux)
          w1 = get_work1(&work, WORK_DXX);

     for (i = 0; i < n_quad; ++i)
          nu[i] = sqrt(fabs(evals[i]));

     if (! aux)
          dmat_mul(tpr, evecs, n_quad, n_quad, n_quad, w1);
     else
          w1 = aux;

     for (i = 0; i < n_quad; ++i) {
          a = 1. / nu[i];
          for (j = 0; j < n_quad; ++j) {
               b = w1[j][i] * a;
               X_p[j][i] = .5 * (evecs[j][i] + b);
               X_m[j][i] = .5 * (evecs[j][i] - b);
          }
     }

     if (flags_or(derivs, n_derivs)) {
          if (! aul) {
               w2 = get_work1(&work, WORK_DXX);
               w3 = get_work1(&work, WORK_DXX);
          }

          for (i = 0; i < n_derivs; ++i) {
               if (! derivs[i])
                    continue;

               for (j = 0; j < n_quad; ++j)
                    nu_l[i][j] = evals_l[i][j];

               if (! aul) {
                    dmat_mul(tpr_l[i], evecs, n_quad, n_quad, n_quad, w2);
                    dmat_mul(tpr, evecs_l[i], n_quad, n_quad, n_quad, w3);
                    dmat_add(w2, w3, w2, n_quad, n_quad);
               }
               else
                    w2 = aul[i];

               for (j = 0; j < n_quad; ++j) {
                    a = -nu_l[i][j] / (nu[j] * nu[j]);
                    b =  1.         /  nu[j];

                    for (k = 0; k < n_quad; ++k) {
                         c = a * w1[k][j] + b * w2[k][j];
                         X_p_l[i][k][j] = .5 * (evecs_l[i][k][j] + c);
                         X_m_l[i][k][j] = .5 * (evecs_l[i][k][j] - c);
                    }
               }
          }
     }
}
else {
     if (save_tree.t) {
          save_tree_encode_s(&save_tree, "eig_1n_to_2n_real");

          if (save_tree_retrieve_data(&save_tree, forward_save_eig_1n_to_2n_real_data, &save))
               forward_save_eig_1n_to_2n_real_alloc(save, n_quad);
     }

     w1 = get_work1(&work, WORK_DXX);

     for (i = 0; i < n_quad; ++i)
          nu[i] = sqrt(fabs(evals[i]));

     if (! aux)
          dmat_mul(tpr, evecs, n_quad, n_quad, n_quad, w1);
     else
          w1 = aux;

     dmat_mul_dinv(w1, nu, w1, n_quad, n_quad);

     if (save_tree.t)
          dmat_copy(save->b, w1, n_quad, n_quad);

     dmat_add(evecs, w1, X_p, n_quad, n_quad);
     dmat_scale(.5, X_p, X_p, n_quad, n_quad);

     dmat_sub(evecs, w1, X_m, n_quad, n_quad);
     dmat_scale(.5, X_m, X_m, n_quad, n_quad);

     if (flags_or(derivs, n_derivs)) {
          w2 = get_work1(&work, WORK_DXX);
          w3 = get_work1(&work, WORK_DXX);

          for (i = 0; i < n_derivs; ++i) {
               if (! derivs[i])
                    continue;

               for (j = 0; j < n_quad; ++j)
                    nu_l[i][j] = evals_l[i][j];

               if (! aul) {
                    dmat_mul(tpr_l[i], evecs, n_quad, n_quad, n_quad, w2);
                    dmat_mul(tpr, evecs_l[i], n_quad, n_quad, n_quad, w3);
                    dmat_add(w2, w3, w2, n_quad, n_quad);
               }
               else
                    w2 = aul[i];

               dmat_mul_diag(w1, nu_l[i], w3, n_quad, n_quad);
               dmat_sub(w2, w3, w2, n_quad, n_quad);

               dmat_mul_dinv(w2, nu,    w2, n_quad, n_quad);

               dmat_add(evecs_l[i], w2, X_p_l[i], n_quad, n_quad);
               dmat_scale(.5, X_p_l[i], X_p_l[i], n_quad, n_quad);

               dmat_sub(evecs_l[i], w2, X_m_l[i], n_quad, n_quad);
               dmat_scale(.5, X_m_l[i], X_m_l[i], n_quad, n_quad);
          }
     }
#ifdef USE_AD_FOR_TL_CALC_EIG_1N_TO_2N_REAL
     eig_1n_to_2n_real_tl_with_ad(n_quad, n_derivs, tpr, tmr, evals, evecs, nu, X_p, X_m, tpr_l, tmr_l, evals_l, evecs_l, nu_l, X_p_l, X_m_l, aux, aul, derivs, save_tree, work);
#endif
}
}



/*******************************************************************************
 *
 ******************************************************************************/
void eig_1n_to_2n_complex(int n_quad, int n_derivs,
                          double **tpr, double **tmr,
                          double  *evals_r, double  *evals_i, double  **evecs,
                          dcomplex  *nu, dcomplex  **X_p, dcomplex  **X_m,
                          double ***tpr_l, double ***tmr_l,
                          double **evals_r_l, double **evals_i_l, double ***evecs_l,
                          dcomplex **nu_l, dcomplex ***X_p_l, dcomplex ***X_m_l,
                          uchar *derivs, save_tree_data save_tree, work_data work) {

     int i;
     int j;
     int k;

     double a;
     double b;

     double *vr1;
     double *vr2;
     double *vr3;
     double *vr4;
     double *vr5;
     double *vr6;
     double *vr7;
     double *vr8;

     dcomplex ac1;
     dcomplex bc1;
     dcomplex ac2;
     dcomplex bc2;

     dcomplex *vc1;
     dcomplex *vc2;
     dcomplex *vc3;
     dcomplex *vc4;

     forward_save_eig_1n_to_2n_complex_data *save;


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (save_tree.t) {
          save_tree_encode_s(&save_tree, "eig_1n_to_2n_complex");

          if (save_tree_retrieve_data(&save_tree, forward_save_eig_1n_to_2n_complex_data, &save))
               forward_save_eig_1n_to_2n_complex_alloc(save, n_quad);
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     vr1 = get_work1(&work, WORK_DX);
     vr2 = get_work1(&work, WORK_DX);
     vr3 = get_work1(&work, WORK_DX);
     vr4 = get_work1(&work, WORK_DX);

     vc1 = get_work1(&work, WORK_ZX);
     vc2 = get_work1(&work, WORK_ZX);
     vc3 = get_work1(&work, WORK_ZX);
     vc4 = get_work1(&work, WORK_ZX);

     if (flags_or(derivs, n_derivs)) {
          vr5 = get_work1(&work, WORK_DX);
          vr6 = get_work1(&work, WORK_DX);
          vr7 = get_work1(&work, WORK_DX);
          vr8 = get_work1(&work, WORK_DX);
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     for (i = 0; i < n_quad; ++i)
          nu[i] = csqrt(evals_r[i] + _Complex_I * evals_i[i]);

     if (flags_or(derivs, n_derivs)) {
          for (i = 0; i < n_derivs; ++i) {
               if (! derivs[i])
                    continue;

               for (j = 0; j < n_quad; ++j) {
                    nu_l[i][j] = evals_r_l[i][j] + _Complex_I * evals_i_l[i][j];
               }
          }
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     for (i = 0; i < n_quad; ++i) {
          if (cimag(nu[i]) == 0.) {
               for (j = 0; j < n_quad; ++j)
                    vr1[j] = evecs[j][i];

               dm_v_mul(tpr, vr1, n_quad, n_quad, vr2);

               for (j = 0; j < n_quad; ++j)
                    vr3[j] = vr2[j] / creal(nu[i]);

               if (save_tree.t)
                    dvec_copy(save->b[i], vr3, n_quad);

               for (j = 0; j < n_quad; ++j) {
                    X_p[j][i] = .5 * (evecs[j][i] + vr3[j]);
                    X_m[j][i] = .5 * (evecs[j][i] - vr3[j]);
               }

               if (flags_or(derivs, n_derivs)) {
                    for (j = 0; j < n_derivs; ++j) {
                         if (! derivs[j])
                              continue;

                         for (k = 0; k < n_quad; ++k)
                              vr3[k] = evecs_l[j][k][i];

                         dm_v_mul(tpr_l[j], vr1, n_quad, n_quad, vr4);
                         dm_v_mul(tpr, vr3, n_quad, n_quad, vr5);
                         dvec_add(vr4, vr5, vr4, n_quad);
if (0) {
                         dvec_scale(creal(nu_l[j][i]) / creal(nu[i]), vr2, vr5, n_quad);

                         dvec_sub(vr4, vr5, vr4, n_quad);

                         dvec_scale(1. / creal(nu[i]), vr4, vr4, n_quad);
}
else {
                         a = creal(nu_l[j][i]);
                         b = creal(nu     [i]);
                         dvec_scale(-a  / (b * b), vr2, vr5, n_quad);
                         dvec_scale( 1. /  b,      vr4, vr6, n_quad);
                         dvec_add(vr5, vr6, vr4, n_quad);
}
                         for (k = 0; k < n_quad; ++k) {
                              X_p_l[j][k][i] = .5 * (evecs_l[j][k][i] + vr4[k]);
                              X_m_l[j][k][i] = .5 * (evecs_l[j][k][i] - vr4[k]);
                         }
                    }
               }
          }

          else
          if (cimag(nu[i]) > 0.) {
               for (j = 0; j < n_quad; ++j) {
                    vr1[j] = evecs[j][i    ];
                    vr2[j] = evecs[j][i + 1];
               }

               for (j = 0; j < n_quad; ++j) {
                    vc1[j] = evecs[j][i] + _Complex_I *  evecs[j][i + 1];
                    vc2[j] = evecs[j][i] + _Complex_I * -evecs[j][i + 1];
               }

               dm_v_mul(tpr, vr1, n_quad, n_quad, vr3);
               dm_v_mul(tpr, vr2, n_quad, n_quad, vr4);

               if (save_tree.t) {
                    dvec_copy(save->vr3[i], vr3, n_quad);
                    dvec_copy(save->vr4[i], vr4, n_quad);
               }

               for (j = 0; j < n_quad; ++j) {
                    vc3[j] = (vr3[j] + _Complex_I *  vr4[j]) / nu[i    ];
                    vc4[j] = (vr3[j] + _Complex_I * -vr4[j]) / nu[i + 1];
               }

               for (j = 0; j < n_quad; ++j) {
                    X_p[j][i    ] = .5 * (vc1[j] + vc3[j]);
                    X_m[j][i    ] = .5 * (vc1[j] - vc3[j]);
                    X_p[j][i + 1] = .5 * (vc2[j] + vc4[j]);
                    X_m[j][i + 1] = .5 * (vc2[j] - vc4[j]);
               }

               if (flags_or(derivs, n_derivs)) {
                    for (j = 0; j < n_derivs; ++j) {
                         if (! derivs[j])
                              continue;

                         for (k = 0; k < n_quad; ++k) {
                              vr5[k] = evecs_l[j][k][i];
                              vr6[k] = evecs_l[j][k][i + 1];
                         }

                         for (k = 0; k < n_quad; ++k) {
                              vc1[k] = evecs_l[j][k][i] + _Complex_I *  evecs_l[j][k][i + 1];
                              vc2[k] = evecs_l[j][k][i] + _Complex_I * -evecs_l[j][k][i + 1];
                         }

                         dm_v_mul(tpr_l[j], vr1, n_quad, n_quad, vr7);
                         dm_v_mul(tpr, vr5, n_quad, n_quad, vr8);
                         dvec_add(vr7, vr8, vr5, n_quad);

                         dm_v_mul(tpr_l[j], vr2, n_quad, n_quad, vr7);
                         dm_v_mul(tpr, vr6, n_quad, n_quad, vr8);
                         dvec_add(vr7, vr8, vr6, n_quad);
if (0) {
                         ac1 = nu_l[j][i    ] / nu[i    ];
                         bc1 = nu_l[j][i + 1] / nu[i + 1];
                         for (k = 0; k < n_quad; ++k) {
                              vc3[k] = ((vr5[k] + _Complex_I *  vr6[k]) -
                                        (vr3[k] + _Complex_I *  vr4[k]) * ac1) / nu[i    ];
                              vc4[k] = ((vr5[k] + _Complex_I * -vr6[k]) -
                                        (vr3[k] + _Complex_I * -vr4[k]) * bc1) / nu[i + 1];
                         }
}
else {
                         ac1 = -nu_l[j][i    ] / (nu[i    ] * nu[i    ]);
                         bc1 =  1.             /  nu[i    ];
                         ac2 = -nu_l[j][i + 1] / (nu[i + 1] * nu[i + 1]);
                         bc2 =  1.             /  nu[i + 1];
                         for (k = 0; k < n_quad; ++k) {
                              vc3[k] = ac1 * (vr3[k] + _Complex_I *  vr4[k]) +
                                       bc1 * (vr5[k] + _Complex_I *  vr6[k]);
                              vc4[k] = ac2 * (vr3[k] + _Complex_I * -vr4[k]) +
                                       bc2 * (vr5[k] + _Complex_I * -vr6[k]);
                         }
}
                         for (k = 0; k < n_quad; ++k) {
                              X_p_l[j][k][i    ] = .5 * (vc1[k] + vc3[k]);
                              X_m_l[j][k][i    ] = .5 * (vc1[k] - vc3[k]);
                              X_p_l[j][k][i + 1] = .5 * (vc2[k] + vc4[k]);
                              X_m_l[j][k][i + 1] = .5 * (vc2[k] - vc4[k]);
                         }
                    }
               }
          }
     }
#ifdef USE_AD_FOR_TL_CALC_EIG_1N_TO_2N_COMPLEX
     eig_1n_to_2n_complex_tl_with_ad(n_quad, n_derivs, tpr, tmr, evals_r, evals_i, evecs, nu, X_p, X_m, tpr_l, tmr_l, evals_r_l, evals_i_l, evecs_l, nu_l, X_p_l, X_m_l, derivs, save_tree, work);
#endif
}



/*******************************************************************************
 *
 ******************************************************************************/
void eig_2n_gen_real(int n_quad, int n_derivs,
                     double  **tpr, double  **tmr, double  **gamma,
                     double  *nu, double  **X_p, double  **X_m, 
                     double ***tpr_l, double ***tmr_l, double ***gamma_l,
                     double **nu_l, double ***X_p_l, double ***X_m_l, 
                     int eigen_solver, uchar *derivs, save_tree_data save_tree,
                     work_data work) {

     double *evals;
     double **evals_l;

     double **evecs;
     double ***evecs_l;

     forward_save_eig_2n_gen_real_data *save;


     if (save_tree.t) {
          save_tree_encode_s(&save_tree, "eig_2n_gen_real");

          if (save_tree_retrieve_data(&save_tree, forward_save_eig_2n_gen_real_data, &save))
               forward_save_eig_2n_gen_real_alloc(save, n_quad);
     }


     evals = get_work1(&work, WORK_DX);
     evecs = get_work1(&work, WORK_DXX);

     if (flags_or(derivs, n_derivs)) {
          evals_l = get_work2(&work, WORK_DX, WORK_DERIVS_V, derivs);
          evecs_l = get_work2(&work, WORK_DXX, WORK_DERIVS_V, derivs);
     }


     eig_1n_gen_real(n_quad, n_derivs, gamma, evals, evecs, gamma_l, evals_l, evecs_l, eigen_solver, derivs, save_tree, work);

     eig_1n_to_2n_real(n_quad, n_derivs, tpr, tmr, evals, evecs, nu, X_p, X_m, tpr_l, tmr_l, evals_l, evecs_l, nu_l, X_p_l, X_m_l, NULL, NULL, derivs, save_tree, work);


     if (save_tree.t) {
          dvec_copy(save->evals, evals, n_quad);
          dmat_copy(save->evecs, evecs, n_quad, n_quad);
     }
}



/*******************************************************************************
 *
 ******************************************************************************/
void eig_2n_gen_complex(int n_quad, int n_derivs,
                        double  **tpr, double  **tmr, double  **gamma,
                        dcomplex  *nu, dcomplex  **X_p, dcomplex  **X_m, 
                        double ***tpr_l, double ***tmr_l, double ***gamma_l,
                        dcomplex **nu_l, dcomplex ***X_p_l, dcomplex ***X_m_l, 
                        int eigen_solver, uchar *derivs, save_tree_data save_tree,
                        work_data work) {

     double *evals_r;
     double *evals_i;

     double **evals_r_l;
     double **evals_i_l;

     double **evecs;

     double ***evecs_l;

     forward_save_eig_2n_gen_complex_data *save;


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (save_tree.t) {
          save_tree_encode_s(&save_tree, "eig_2n_gen_complex");

          if (save_tree_retrieve_data(&save_tree, forward_save_eig_2n_gen_complex_data, &save))
               forward_save_eig_2n_gen_complex_alloc(save, n_quad);
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     evals_r = get_work1(&work, WORK_DX);
     evals_i = get_work1(&work, WORK_DX);

     evecs   = get_work1(&work, WORK_DXX);

     if (flags_or(derivs, n_derivs)) {
          evals_r_l = get_work2(&work, WORK_DX, WORK_DERIVS_V, derivs);
          evals_i_l = get_work2(&work, WORK_DX, WORK_DERIVS_V, derivs);

          evecs_l   = get_work2(&work, WORK_DXX, WORK_DERIVS_V, derivs);
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     eig_1n_gen_complex (n_quad, n_derivs, gamma, evals_r, evals_i, evecs, gamma_l, evals_r_l, evals_i_l, evecs_l, eigen_solver, derivs, save_tree, work);
/*
     eig_1n_gen_complex2(n_quad, n_derivs, gamma, evals_r, evals_i, evecs, gamma_l, evals_r_l, evals_i_l, evecs_l,               derivs, save_tree, work);
*/
     eig_1n_to_2n_complex(n_quad, n_derivs, tpr, tmr, evals_r, evals_i, evecs, nu, X_p, X_m, tpr_l, tmr_l, evals_r_l, evals_i_l, evecs_l, nu_l, X_p_l, X_m_l, derivs, save_tree, work);


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (save_tree.t) {
          dvec_copy(save->evals_r, evals_r, n_quad);
          dvec_copy(save->evals_i, evals_i, n_quad);
          dmat_copy(save->evecs, evecs, n_quad, n_quad);
     }
}



/*******************************************************************************
 *
 ******************************************************************************/
#ifdef INCLUDE_DEV_SOURCE
#include "../src2/xrtm_eig_util2.c"
#endif
