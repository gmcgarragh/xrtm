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
#include "xrtm_eig_util.h"
#include "xrtm_eig_util_a.h"
#include "xrtm_matrix.h"
#include "xrtm_save_tree.h"
#include "xrtm_utility.h"


/*******************************************************************************
 *
 ******************************************************************************/
static void eigen_derivatives_real_a(int i, int n_quad, double **gamma, double *evals_r, double *evals_i, double **evecs, double **gamma_a, double *evals_r_a, double *evals_i_a, double **evecs_a, double **A, double *B, int *ip, int flag) {

     int j;
     int k;

     double a;
     double b;

     if (is_trivial_eigenvector(i, n_quad, evecs))
/*
     if (is_dummy_node_eigenvector(i, n_quad, evecs))
*/
          return;

     if (is_degenerate_eigenvalue(i, n_quad, evals_r))
          return;
if (! flag) {
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
/*
     if (evals_i_a)
          perturb_zero_elems(n_quad + 1, A);
*/
     dmat_getrf(A, n_quad + 1, n_quad + 1, ip);
}
     B[0] = evals_r_a[i];

     for (j = 0; j < n_quad; ++j)
          B[j+1] = evecs_a[j][i];

     dmat_getrs2('t', A, &B, n_quad + 1, 1, ip);

     for (j = 0; j < n_quad; ++j) {
          for (k = 0; k < n_quad; ++k) {
               gamma_a[j][k] += B[j] * evecs[k][i];
          }
     }

     evals_r_a[i] = 0.;

     for (j = 0; j < n_quad; ++j)
          evecs_a[j][i] = 0.;
}



/*******************************************************************************
 *
 ******************************************************************************/
static void forward_save_eig_1n_gen_real_free(forward_save_eig_1n_gen_real_data *d);

int forward_save_eig_1n_gen_real_alloc(forward_save_eig_1n_gen_real_data *d, int n_quad) {

     d->free = (void (*)(void *)) forward_save_eig_1n_gen_real_free;

     d->ip = alloc_array2_i(n_quad, n_quad + 1);

     d->A  = alloc_array3_d(n_quad, n_quad + 1, n_quad + 1);

     return 0;
}



static void forward_save_eig_1n_gen_real_free(forward_save_eig_1n_gen_real_data *d) {

     free_array2_i(d->ip);

     free_array3_d(d->A);
}



/*******************************************************************************
 *
 ******************************************************************************/
void eig_1n_gen_real_a(int n_quad,
                       double **gamma, double *evals, double **evecs,
                       double **gamma_a, double *evals_a, double **evecs_a,
                       save_tree_data save_tree, work_data work) {

     int i;
/*
     int *i1;
*/
     double *B;
/*
     double **A;
*/
     forward_save_eig_1n_gen_real_data *save;

     save_tree_encode_s(&save_tree, "eig_1n_gen_real");

     save_tree_retrieve_data(&save_tree, forward_save_eig_1n_gen_real_data, &save);
/*
     i1 = get_work_i1(&work, n_quad + 1);
*/
     B  = get_work_d1(&work, n_quad + 1);
/*
     A  = get_work_d2(&work, n_quad + 1, n_quad + 1);
*/
     for (i = 0; i < n_quad; ++i)
          eigen_derivatives_real_a(i, n_quad, gamma, evals, NULL, evecs, gamma_a, evals_a, NULL, evecs_a, save->A[i], B, save->ip[i], save->flag);
save->flag = 1;
}



/*******************************************************************************
 *
 ******************************************************************************/
static void forward_save_eig_1n_gen_complex_free(forward_save_eig_1n_gen_complex_data *d);

int forward_save_eig_1n_gen_complex_alloc(forward_save_eig_1n_gen_complex_data *d, int n_quad) {

     d->free = (void (*)(void *)) forward_save_eig_1n_gen_complex_free;

     d->ip = alloc_array2_i(n_quad, n_quad + 1);

     d->A  = alloc_array3_d(n_quad, n_quad + 1, n_quad + 1);

     return 0;
}



static void forward_save_eig_1n_gen_complex_free(forward_save_eig_1n_gen_complex_data *d) {

     free_array2_i(d->ip);

     free_array3_d(d->A);
}



/*******************************************************************************
 *
 ******************************************************************************/
void eig_1n_gen_complex_a(int n_quad,
                          double **gamma,
                          double *evals_r, double *evals_i, double **evecs,
                          double **gamma_a,
                          double *evals_r_a, double *evals_i_a, double **evecs_a,
                          save_tree_data save_tree, work_data work) {

     int i;
     int i1;
     int i2;
     int j;
     int k;

     int *ip;

     double *Br;
/*
     double **Ar;
*/
     dcomplex ac;
     dcomplex bc;

     dcomplex *Bc;

     dcomplex **Ac;

     forward_save_eig_1n_gen_complex_data *save;

     save_tree_encode_s(&save_tree, "eig_1n_gen_complex");

     save_tree_retrieve_data(&save_tree, forward_save_eig_1n_gen_complex_data, &save);

     ip  = get_work_i1(&work, n_quad + 1);

     Br  = get_work_d1(&work, n_quad + 1);
/*
     Ar  = get_work_d2(&work, n_quad + 1, n_quad + 1);
*/
     Bc  = get_work_dc1(&work, n_quad + 1);
     Ac  = get_work_dc2(&work, n_quad + 1, n_quad + 1);

     for (i = 0; i < n_quad; ++i) {
          if (evals_i[i] == 0.)
               eigen_derivatives_real_a(i, n_quad, gamma, evals_r, evals_i, evecs, gamma_a, evals_r_a, evals_i_a, evecs_a,save->A[i], Br, save->ip[i], save->flag);
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

               Bc[0]  = evals_r_a[i1] + _Complex_I * evals_i_a[i1];
               Bc[0] += evals_r_a[i2] - _Complex_I * evals_i_a[i2];

               for (k = 0; k < n_quad; ++k)
                    Bc[k+1] = evecs_a[k][i1] + _Complex_I * evecs_a[k][i2];

               zmat_getrs2('t', Ac, &Bc, n_quad + 1, 1, ip);

               for (j = 0; j < n_quad; ++j) {
                    for (k = 0; k < n_quad; ++k) {
                         gamma_a[j][k] += creal(Bc[j] * (evecs[k][i1] + _Complex_I * evecs[k][i2]));

                    }
               }
          }
     }

     dvec_zero(evals_r_a, n_quad);
     dvec_zero(evals_i_a, n_quad);
     dmat_zero(evecs_a, n_quad, n_quad);
save->flag = 1;
}



/*******************************************************************************
 *
 ******************************************************************************/
static void forward_save_eig_1n_to_2n_real_free(forward_save_eig_1n_to_2n_real_data *d);

int forward_save_eig_1n_to_2n_real_alloc(forward_save_eig_1n_to_2n_real_data *d, int n_quad) {

     d->free = (void (*)(void *)) forward_save_eig_1n_to_2n_real_free;

     d->b = alloc_array2_d(n_quad, n_quad);

     return 0;
}



static void forward_save_eig_1n_to_2n_real_free(forward_save_eig_1n_to_2n_real_data *d) {

     free_array2_d(d->b);
}



/*******************************************************************************
 *
 ******************************************************************************/
void eig_1n_to_2n_real_a(int n_quad,
                         double **tpr, double **tmr,
                         double *evals, double **evecs,
                         double *nu, double **X_p, double **X_m,
                         double **tpr_a, double **tmr_a,
                         double *evals_a, double **evecs_a,
                         double *nu_a, double **X_p_a, double **X_m_a,
                         save_tree_data save_tree, work_data work) {

     int i;

     int j;

     double **w1;
/*
     double **w2;
     double **w3;
*/
     forward_save_eig_1n_to_2n_real_data *save;

     save_tree_encode_s(&save_tree, "eig_1n_to_2n_real");

     save_tree_retrieve_data(&save_tree, forward_save_eig_1n_to_2n_real_data, &save);

     w1 = get_work1(&work, WORK_DXX);
/*
     w2 = get_work1(&work, WORK_DXX);
     w3 = get_work1(&work, WORK_DXX);
*/
     dmat_add(X_p_a, X_m_a, w1, n_quad, n_quad);
     dmat_scale(.5, w1, w1, n_quad, n_quad);
     dmat_add(evecs_a, w1, evecs_a, n_quad, n_quad);

     dmat_sub(X_p_a, X_m_a, w1, n_quad, n_quad);
     dmat_scale(.5, w1, w1, n_quad, n_quad);

     dmat_mul_dinv(w1, nu, w1, n_quad, n_quad);
/*
     dmat_trans(evecs, w2, n_quad, n_quad);
     dmat_mul(w1, w2, n_quad, n_quad, n_quad, w3);
     dmat_add(tpr_a, w3, tpr_a, n_quad, n_quad);
*/
     dmat_gxgxmx(0, w1, 1, evecs, 1., tpr_a, 1., n_quad, n_quad, n_quad);
/*
     dmat_trans(tpr, w2, n_quad, n_quad);
     dmat_mul(w2, w1, n_quad, n_quad, n_quad, w3);
     dmat_add(evecs_a, w3, evecs_a, n_quad, n_quad);
*/
     dmat_gxgxmx(1, tpr, 0, w1, 1., evecs_a, 1., n_quad, n_quad, n_quad);

     for (i = 0; i < n_quad; ++i) {
          for (j = 0; j < n_quad; ++j) {
               nu_a[i] -= save->b[j][i] * w1[j][i];
          }
     }

     dvec_add(evals_a, nu_a, evals_a, n_quad);

     dvec_zero(nu_a,  n_quad);
     dmat_zero(X_p_a, n_quad, n_quad);
     dmat_zero(X_m_a, n_quad, n_quad);
}



/*******************************************************************************
 *
 ******************************************************************************/
static void forward_save_eig_1n_to_2n_complex_free(forward_save_eig_1n_to_2n_complex_data *d);

int forward_save_eig_1n_to_2n_complex_alloc(forward_save_eig_1n_to_2n_complex_data *d, int n_quad) {

     d->free = (void (*)(void *)) forward_save_eig_1n_to_2n_complex_free;

     d->b   = alloc_array2_d(n_quad, n_quad);
     d->vr3 = alloc_array2_d(n_quad, n_quad);
     d->vr4 = alloc_array2_d(n_quad, n_quad);

     return 0;
}



static void forward_save_eig_1n_to_2n_complex_free(forward_save_eig_1n_to_2n_complex_data *d) {

     free_array2_d(d->b);
     free_array2_d(d->vr3);
     free_array2_d(d->vr4);
}



/*******************************************************************************
 *
 ******************************************************************************/
void eig_1n_to_2n_complex_a(int n_quad,
                            double **tpr, double **tmr,
                            double *evals_r, double *evals_i, double **evecs,
                            dcomplex *nu, dcomplex **X_p, dcomplex **X_m,
                            double **tpr_a, double **tmr_a,
                            double *evals_r_a, double *evals_i_a, double **evecs_a,
                            dcomplex *nu_a, dcomplex **X_p_a, dcomplex **X_m_a,
                            save_tree_data save_tree, work_data work) {

     int i;
     int j;
     int k;

     double *v1;

     double **w1;

     double *t;
     double *d_a;

     double *vr5_a;
     double *vr6_a;
     double *vr7_a;
     double *vr8_a;

     dcomplex ac1_a;
     dcomplex bc1;
     dcomplex ac2_a;
     dcomplex bc2;

     dcomplex *vc1_a;
     dcomplex *vc2_a;
     dcomplex *vc3_a;
     dcomplex *vc4_a;

     forward_save_eig_1n_to_2n_complex_data *save;


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     save_tree_encode_s(&save_tree, "eig_1n_to_2n_complex");

     save_tree_retrieve_data(&save_tree, forward_save_eig_1n_to_2n_complex_data, &save);


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     v1    = get_work1(&work, WORK_DX);

     vr5_a = get_work1(&work, WORK_DX);
     vr6_a = get_work1(&work, WORK_DX);
     vr7_a = get_work1(&work, WORK_DX);
     vr8_a = get_work1(&work, WORK_DX);

     vc1_a = get_work1(&work, WORK_ZX);
     vc2_a = get_work1(&work, WORK_ZX);
     vc3_a = get_work1(&work, WORK_ZX);
     vc4_a = get_work1(&work, WORK_ZX);

     t     = get_work1(&work, WORK_DX);
     d_a   = get_work1(&work, WORK_DX);

     w1    = get_work1(&work, WORK_DXX);


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     for (i = 0; i < n_quad; ++i) {
          if (cimag(nu[i]) == 0.) {
               for (j = 0; j < n_quad; ++j) {
                   evecs_a[j][i] += .5 * creal(X_p_a[j][i] + X_m_a[j][i]);
                   d_a    [j]     = .5 * creal(X_p_a[j][i] - X_m_a[j][i]);
               }

               dvec_scale(1. / creal(nu[i]), d_a, t, n_quad);

               for (j = 0; j < n_quad; ++j) {
                    for (k = 0; k < n_quad; ++k) {
                         tpr_a[j][k] += t[j] * evecs[k][i];
                    }
               }

               dmat_trans(tpr, w1, n_quad, n_quad);
               dm_v_mul(w1, t, n_quad, n_quad, v1);
               for (j = 0; j < n_quad; ++j)
                    evecs_a[j][i] += v1[j];

                nu_a[i] -= dvec_dot(save->b[i], t, n_quad);

                evals_r_a[i] += creal(nu_a[i]);
          }

          else
          if (cimag(nu[i]) > 0.) {
               for (j = 0; j < n_quad; ++j) {
                   vc1_a[j] = .5 * (X_p_a[j][i    ] + X_m_a[j][i    ]);
                   vc3_a[j] = .5 * (X_p_a[j][i    ] - X_m_a[j][i    ]);
                   vc2_a[j] = .5 * (X_p_a[j][i + 1] + X_m_a[j][i + 1]);
                   vc4_a[j] = .5 * (X_p_a[j][i + 1] - X_m_a[j][i + 1]);
               }

               ac1_a = 0.;
               bc1   = 1. / nu[i    ];
               ac2_a = 0.;
               bc2   = 1. / nu[i + 1];
               for (j = 0; j < n_quad; ++j) {
                    ac1_a += vc3_a[j] * (save->vr3[i][j] + _Complex_I * save->vr4[i][j]);
                    vr5_a[j]  =  creal(bc1 * vc3_a[j]);
                    vr6_a[j]  =  cimag(bc1 * vc3_a[j]);

                    ac2_a += vc4_a[j] * (save->vr3[i][j] - _Complex_I * save->vr4[i][j]);
                    vr5_a[j] +=  creal(bc2 * vc4_a[j]);
                    vr6_a[j] += -cimag(bc2 * vc4_a[j]);
               }

               nu_a[i    ] -= ac1_a / (nu[i    ] * nu[i    ]);
               nu_a[i + 1] -= ac2_a / (nu[i + 1] * nu[i + 1]);

               dmat_trans(tpr, w1, n_quad, n_quad);

               for (j = 0; j < n_quad; ++j) {
                    for (k = 0; k < n_quad; ++k) {
                         tpr_a[j][k] += vr5_a[j] * evecs[k][i    ];
                    }
               }
               dm_v_mul(w1, vr5_a, n_quad, n_quad, vr7_a);

               for (j = 0; j < n_quad; ++j) {
                    for (k = 0; k < n_quad; ++k) {
                         tpr_a[j][k] += vr6_a[j] * -evecs[k][i + 1];
                    }
               }
               dm_v_mul(w1, vr6_a, n_quad, n_quad, vr8_a);

               for (j = 0; j < n_quad; ++j) {
                    evecs_a[j][i    ] += vr7_a[j];
                    evecs_a[j][i + 1] += vr8_a[j];
               }

               for (j = 0; j < n_quad; ++j) {
                    evecs_a[j][i    ] += creal(vc1_a[j]) + creal(vc2_a[j]);
                    evecs_a[j][i + 1] += cimag(vc1_a[j]) - cimag(vc2_a[j]);
               }

               evals_r_a[i    ] += creal(nu_a[i    ]);
               evals_i_a[i    ] += cimag(nu_a[i    ]);

               evals_r_a[i + 1] += creal(nu_a[i + 1]);
               evals_i_a[i + 1] += cimag(nu_a[i + 1]);
          }
     }

     zvec_zero(nu_a,  n_quad);
     zmat_zero(X_p_a, n_quad, n_quad);
     zmat_zero(X_m_a, n_quad, n_quad);
}



/*******************************************************************************
 *
 ******************************************************************************/
static void forward_save_eig_2n_gen_real_free(forward_save_eig_2n_gen_real_data *d);

int forward_save_eig_2n_gen_real_alloc(forward_save_eig_2n_gen_real_data *d, int n_quad) {

     d->free = (void (*)(void *)) forward_save_eig_2n_gen_real_free;

     d->evals = alloc_array1_d(n_quad);
     d->evecs = alloc_array2_d(n_quad, n_quad);

     return 0;
}



static void forward_save_eig_2n_gen_real_free(forward_save_eig_2n_gen_real_data *d) {

     free_array1_d(d->evals);
     free_array2_d(d->evecs);
}



/*******************************************************************************
 *
 ******************************************************************************/
void eig_2n_gen_real_a(int n_quad,
                       double **tpr, double **tmr, double **gamma,
                       double *nu, double **X_p, double **X_m,
                       double **tpr_a, double **tmr_a, double **gamma_a,
                       double *nu_a, double **X_p_a, double **X_m_a,
                       save_tree_data save_tree, work_data work) {

     double *evals_a;
     double **evecs_a;

     forward_save_eig_2n_gen_real_data *save;


     save_tree_encode_s(&save_tree, "eig_2n_gen_real");

     save_tree_retrieve_data(&save_tree, forward_save_eig_2n_gen_real_data, &save);


     evals_a = get_work1(&work, WORK_DX);
     evecs_a = get_work1(&work, WORK_DXX);


     dvec_zero(evals_a, n_quad);
     dmat_zero(evecs_a, n_quad, n_quad);

     eig_1n_to_2n_real_a(n_quad, tpr, tmr, save->evals, save->evecs, nu, X_p, X_m, tpr_a, tmr_a, evals_a, evecs_a, nu_a, X_p_a, X_m_a, save_tree, work);

     eig_1n_gen_real_a(n_quad, gamma, save->evals, save->evecs, gamma_a, evals_a, evecs_a, save_tree, work);
}



/*******************************************************************************
 *
 ******************************************************************************/
static void forward_save_eig_2n_gen_complex_free(forward_save_eig_2n_gen_complex_data *d);

int forward_save_eig_2n_gen_complex_alloc(forward_save_eig_2n_gen_complex_data *d, int n_quad) {

     d->free = (void (*)(void *)) forward_save_eig_2n_gen_complex_free;

     d->evals_r = alloc_array1_d(n_quad);
     d->evals_i = alloc_array1_d(n_quad);
     d->evecs   = alloc_array2_d(n_quad, n_quad);

     return 0;
}



static void forward_save_eig_2n_gen_complex_free(forward_save_eig_2n_gen_complex_data *d) {

     free_array1_d(d->evals_r);
     free_array1_d(d->evals_i);
     free_array2_d(d->evecs);
}



/*******************************************************************************
 *
 ******************************************************************************/
void eig_2n_gen_complex_a(int n_quad,
                          double  **tpr, double  **tmr, double  **gamma,
                          dcomplex  *nu, dcomplex  **X_p, dcomplex  **X_m,
                          double **tpr_a, double **tmr_a, double **gamma_a,
                          dcomplex *nu_a, dcomplex **X_p_a, dcomplex **X_m_a,
                          save_tree_data save_tree, work_data work) {

     double *evals_r_a;
     double *evals_i_a;
     double **evecs_a;

     forward_save_eig_2n_gen_complex_data *save;


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     save_tree_encode_s(&save_tree, "eig_2n_gen_complex");

     save_tree_retrieve_data(&save_tree, forward_save_eig_2n_gen_complex_data, &save);


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     evals_r_a = get_work1(&work, WORK_DX);
     evals_i_a = get_work1(&work, WORK_DX);

     evecs_a   = get_work1(&work, WORK_DXX);


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     dvec_zero(evals_r_a, n_quad);
     dvec_zero(evals_i_a, n_quad);
     dmat_zero(evecs_a, n_quad, n_quad);

     eig_1n_to_2n_complex_a(n_quad, tpr, tmr, save->evals_r, save->evals_i, save->evecs, nu, X_p, X_m, tpr_a, tmr_a, evals_r_a, evals_i_a, evecs_a, nu_a, X_p_a, X_m_a, save_tree, work);

     eig_1n_gen_complex_a (n_quad, gamma, save->evals_r, save->evals_i, save->evecs, gamma_a, evals_r_a, evals_i_a, evecs_a, save_tree, work);
/*
     eig_1n_gen_complex2_a(n_quad, gamma, save->evals_r, save->evals_i, save->evecs, gamma_a, evals_r_a, evals_i_a, evecs_a, save_tree, work);
*/
}



/*******************************************************************************
 *
 ******************************************************************************/
#ifdef INCLUDE_DEV_SOURCE
#include "../src2/xrtm_eig_util_a2.c"
#endif
