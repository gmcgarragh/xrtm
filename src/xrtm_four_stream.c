/*******************************************************************************
**
**    Copyright (C) 2007-2020 Greg McGarragh <greg.mcgarragh@colostate.edu>
**
**    This source code is licensed under the GNU General Public License (GPL),
**    Version 3.  See the file COPYING for more details.
**
*******************************************************************************/

#include <gutil.h>
#include <gmath.h>
#include <gmath_matrix.h>

#include "xrtm.h"
#include "xrtm_derivs.h"
#include "xrtm_four_stream.h"
#include "xrtm_matrix.h"
#include "xrtm_save_tree.h"
#include "xrtm_scatter.h"
#include "xrtm_sfi.h"


#define CODE_TYPE 2

#if CODE_TYPE == 0
#define USE_EXPANDED
#elif CODE_TYPE == 1
#define USE_MAT2D_FUNCTIONS
#elif CODE_TYPE == 2
#define USE_EXPANDED
#define USE_ARRAYS_OF_POINTERS
#elif CODE_TYPE == 3
#define USE_MATXD_FUNCTIONS
#define USE_ARRAYS_OF_POINTERS
#endif

#define USE_BANDED_SOLVER


#ifdef USE_BANDED_SOLVER
#ifdef __cplusplus
extern "C" {
#endif

void dgbtrf_(int *, int *, int *, int *, double *, int *, int *, int *);
void dgbtrs_(const char *, int *, int *, int *, int *, double *, int *, int *, double *, int *, int *);

#ifdef __cplusplus
}
#endif
#endif


/*******************************************************************************
 *
 ******************************************************************************/
static void calc_Y_1(int i_four, double *Y) {

     /* This code was automatically generated
        with function print_Y_x_constant_code() */

     if (i_four == 0) {
          Y[0] =  1.0000000000000000e+00;
          Y[1] =  2.1132486540518713e-01;
          Y[2] = -4.3301270189221930e-01;
          Y[3] = -2.9339382851364088e-01;
     }
     else if (i_four == 1) {
          Y[1] = -6.9113739634803417e-01;
          Y[2] = -2.5297384456881489e-01;
          Y[3] =  3.2872926407690672e-01;
     }
     else if (i_four == 2) {
          Y[2] =  5.8502498576049988e-01;
          Y[3] =  2.7644581385388761e-01;
     }
     else if (i_four == 3) {
          Y[3] = -5.2199120072772409e-01;
     }
}



static void calc_Y_2(int i_four, double *Y) {

     /* This code was automatically generated
        with function print_Y_x_constant_code() */

     if (i_four == 0) {
          Y[0] =  1.0000000000000000e+00;
          Y[1] =  7.8867513459481287e-01;
          Y[2] =  4.3301270189221941e-01;
          Y[3] =  4.3393828513641042e-02;
     }
     else if (i_four == 1) {
          Y[1] = -4.3473643283710062e-01;
          Y[2] = -5.9386101120610690e-01;
          Y[3] = -5.6173675511805921e-01;
     }
     else if (i_four == 2) {
          Y[2] =  2.3147159516722612e-01;
          Y[3] =  4.0820738302757004e-01;
     }
     else if (i_four == 3) {
          Y[3] = -1.2991165542275507e-01;
     }
}



/*******************************************************************************
 *
 ******************************************************************************/
void rtm_four_stream(int i_four,
                     int n_derivs, int n_layers,
                     double *qx, double *qw, double F_0, double mu_0,
                     int n_ulevels, int *ulevels, double *utaus,
                     int n_umus, double *umus,
                     double ***coef, double ****coef_l,
                     double *omega, double **omega_l, double *ltau, double **ltau_l,
                     double *btran, double **btran_l,
                     double *as_0, double **as_0_l, double *atran, double **atran_l,
                     double *Rs_q0, double **Rs_q0_l, double **Rs_qq, double ***Rs_qq_l,
                     double *Rs_u0, double **Rs_u0_l, double **Rs_uq, double ***Rs_uq_l,
                     double *I1_m, double **I1_m_l, double *In_p, double **In_p_l,
                     double **I_p, double **I_m, double ***I_p_l, double ***I_m_l,
                     int add_single_scattering, int solar, int thermal, int surface,
                     int upwelling, int downwelling, int utau_output,
                     derivs_data *derivs, save_tree_data save_tree, work_data work) {

     int i;
     int ii;
     int iii;
     int j;
     int jj;
     int k;
     int l;
#ifdef USE_BANDED_SOLVER
     int n_diags;
#endif
     int m_comp;
     int n_comp;
#ifdef USE_BANDED_SOLVER
     int info;
     int nrhs = 1;
#endif
     int *ipiv;

     double a1;
     double a2;
     double a3;
#ifdef USE_EXPANDED
     double a4;
#endif
     double uatran;

     double v1[2];
#if defined USE_MAT2D_FUNCTIONS || defined USE_MATXD_FUNCTIONS
     double v2[2];
     double v3[2];
#endif
     double Y_0[6];

     double evals[2];

     double P_q0_pp[2];
     double P_q0_mp[2];

     double d[2];
     double e[2];
     double g[2];
     double h[2];
     double p[2];

     double *b;

     double **pp1;
     double **pp2;

     double **w1;
     double **w2;

     double ***X_p;
     double ***X_m;

     double **A;
#ifndef USE_ARRAYS_OF_POINTERS
     double W1[2][2];

     double Y_q[6][6];

     double P_qq_pp[2][2];
     double P_qq_mp[2][2];

     double alpha[2][2];
     double beta [2][2];

     double apb[2][2];
     double amb[2][2];

     double gamma[2][2];

     double evecs[2][2];

     double nu[n_layers][2];

     double f[2][2];

     double F_p[n_layers][2];
     double F_m[n_layers][2];

     double Lambda[n_layers][2];
#else
     double **W1;

     double **Y_q;

     double **P_qq_pp;
     double **P_qq_mp;

     double **alpha;
     double **beta;

     double **apb;
     double **amb;

     double **gamma;

     double **evecs;

     double **nu;

     double **f;

     double **F_p;
     double **F_m;

     double **Lambda;
#endif

     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
#ifndef USE_BANDED_SOLVER
     m_comp    = 2 * 2 * n_layers;
#else
     n_diags   = 3 * 2 - 1;
     m_comp    = 3 * n_diags + 1;
#endif
     n_comp    = 2 * 2 * n_layers;


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     ipiv    = get_work_i1(&work, n_comp);

     w1      = get_work_d2(&work, 2, 2);
     w2      = get_work_d2(&work, 2, 2);

     X_p     = get_work_d3(&work, n_layers, 2, 2);
     X_m     = get_work_d3(&work, n_layers, 2, 2);

     b       = get_work_d1(&work, n_comp);

     A       = get_work_d2(&work, n_comp, m_comp);
#ifdef USE_ARRAYS_OF_POINTERS
     W1      = get_work_d2(&work, 2, 2);

     Y_q     = get_work_d2(&work, 6, 6);

     P_qq_pp = get_work_d2(&work, 2, 2);
     P_qq_mp = get_work_d2(&work, 2, 2);

     alpha   = get_work_d2(&work, 2, 2);
     beta    = get_work_d2(&work, 2, 2);

     apb     = get_work_d2(&work, 2, 2);
     amb     = get_work_d2(&work, 2, 2);

     gamma   = get_work_d2(&work, 2, 2);

     evecs   = get_work_d2(&work, 2, 2);

     nu      = get_work_d2(&work, n_layers, 2);

     f       = get_work_d2(&work, 2, 2);

     F_p     = get_work_d2(&work, n_layers, 2);
     F_m     = get_work_d2(&work, n_layers, 2);

     Lambda  = get_work_d2(&work, n_layers, 2);
#endif

     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
/*
     print_Y_x_constant_code(4, 4, qx[0]);
     printf("\n");

     print_Y_x_constant_code(4, 4, qx[1]);

     exit(1);
*/
     calc_Y_1(i_four, Y_q[0]);
     calc_Y_2(i_four, Y_q[1]);
/*
     calc_Y_x(i_four, qx[0], Y_q[0]);
     calc_Y_x(i_four, qx[1], Y_q[1]);
*/
     calc_Y_x(i_four, -mu_0, Y_0);

if (n_derivs > 0) {
     for (i = 0; i < n_ulevels; ++i) {
          for (j = 0; j < n_derivs; ++j) {
               for (k = 0; k < 2 + n_umus; ++k) {
                    I_p_l[i][j][k] = 0.;
                    I_m_l[i][j][k] = 0.;
               }
          }
     }
}

     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     for (i = 0; i < n_layers; ++i) {

          for (j = 0; j < 2; ++j) {
               for (k = 0; k < 2; ++k) {
                    P_qq_pp[j][k] = 0.;
                    P_qq_mp[j][k] = 0.;
               }
          }

          for (l = i_four    ; l < 4; l += 2) {
               for (j = 0; j < 2; ++j) {
                    a1 = coef[i][0][l] * Y_q[j][l];
                    for (k = 0; k < 2; ++k) {
                         P_qq_pp[j][k] += a1 * Y_q[k][l];
                    }
               }
          }

          for (l = i_four + 1; l < 4; l += 2) {
               for (j = 0; j < 2; ++j) {
                    a1 = coef[i][0][l] * Y_q[j][l];
                    for (k = 0; k < 2; ++k) {
                         P_qq_mp[j][k] += a1 * Y_q[k][l];
                    }
               }
          }

          a1 = 2. - (i_four == 0 ? 1. : 0.);
          for (j = 0; j < 2; ++j) {
               for (k = 0; k < 2; ++k) {
                    a2 = (P_qq_pp[j][k] + P_qq_mp[j][k]) * a1;
                    a3 = (P_qq_pp[j][k] - P_qq_mp[j][k]) * a1;
                    P_qq_pp[j][k] = a2;
                    P_qq_mp[j][k] = a3;
               }
          }
/*
          build_phase_mats_scalar(i_four, 4, 2, 2, Y_q, Y_q, 0, coef[i][0], P_qq_pp, P_qq_mp);
*/

          /*--------------------------------------------------------------------
           *
           *------------------------------------------------------------------*/
          a1 = (1. + (i_four == 0 ? 1. : 0.)) / 4.;

          a2 = a1 * qw[0] * omega[i];
          a3 = a1 * qw[1] * omega[i];

          alpha[0][0] = (-1. + a2 * P_qq_pp[0][0]) / qx[0];
          alpha[0][1] = (      a3 * P_qq_pp[0][1]) / qx[0];
          alpha[1][0] = (      a2 * P_qq_pp[1][0]) / qx[1];
          alpha[1][1] = (-1. + a3 * P_qq_pp[1][1]) / qx[1];

          beta [0][0] = (    - a2 * P_qq_mp[0][0]) / qx[0];
          beta [0][1] = (    - a3 * P_qq_mp[0][1]) / qx[0];
          beta [1][0] = (    - a2 * P_qq_mp[1][0]) / qx[1];
          beta [1][1] = (    - a3 * P_qq_mp[1][1]) / qx[1];


          /*--------------------------------------------------------------------
           *
           *------------------------------------------------------------------*/
#ifdef USE_EXPANDED
          apb[0][0] = alpha[0][0] + beta[0][0];
          apb[0][1] = alpha[0][1] + beta[0][1];
          apb[1][0] = alpha[1][0] + beta[1][0];
          apb[1][1] = alpha[1][1] + beta[1][1];

          amb[0][0] = alpha[0][0] - beta[0][0];
          amb[0][1] = alpha[0][1] - beta[0][1];
          amb[1][0] = alpha[1][0] - beta[1][0];
          amb[1][1] = alpha[1][1] - beta[1][1];
#endif
#ifdef USE_MAT2D_FUNCTIONS
          mat2d_add(alpha, beta, apb);
          mat2d_sub(alpha, beta, amb);
#endif
#ifdef USE_MATXD_FUNCTIONS
          dmat_add(alpha, beta, apb, 2, 2);
          dmat_sub(alpha, beta, amb, 2, 2);
#endif

          /*--------------------------------------------------------------------
           *
           *------------------------------------------------------------------*/
#ifdef USE_EXPANDED
          gamma[0][0] = amb[0][0] * apb[0][0] + amb[0][1] * apb[1][0];
          gamma[0][1] = amb[0][0] * apb[0][1] + amb[0][1] * apb[1][1];
          gamma[1][0] = amb[1][0] * apb[0][0] + amb[1][1] * apb[1][0];
          gamma[1][1] = amb[1][0] * apb[0][1] + amb[1][1] * apb[1][1];
#endif
#ifdef USE_MAT2D_FUNCTIONS
          mat2d_mul(amb, apb, gamma);
#endif
#ifdef USE_MATXD_FUNCTIONS
          dmat_mul(amb, apb, 2, 2, 2, gamma);
#endif

          /*--------------------------------------------------------------------
           *
           *------------------------------------------------------------------*/
          a1 = gamma[0][0] + gamma[1][1];
          a2 = gamma[0][0] - gamma[1][1];
          a3 = a2 * a2;

          evals[0] = (a1 + sqrt(a3 + 4. * gamma[0][1] * gamma[1][0])) / 2.,
          evals[1] = (a1 - sqrt(a3 + 4. * gamma[1][0] * gamma[0][1])) / 2.;

          evecs[0][0] = 1.;
          evecs[0][1] = gamma[0][1] / (evals[1] - gamma[0][0]);
          evecs[1][0] = gamma[1][0] / (evals[0] - gamma[1][1]);
          evecs[1][1] = 1.;


          /*--------------------------------------------------------------------
           *
           *------------------------------------------------------------------*/
          nu[i][0] = sqrt(evals[0]);
          nu[i][1] = sqrt(evals[1]);
#ifdef USE_EXPANDED
          W1[0][0] = apb[0][0] * evecs[0][0] + apb[0][1] * evecs[1][0];
          W1[0][1] = apb[0][0] * evecs[0][1] + apb[0][1] * evecs[1][1];
          W1[1][0] = apb[1][0] * evecs[0][0] + apb[1][1] * evecs[1][0];
          W1[1][1] = apb[1][0] * evecs[0][1] + apb[1][1] * evecs[1][1];
#endif
#ifdef USE_MAT2D_FUNCTIONS
          mat2d_mul(apb, evecs, W1);
#endif
#ifdef USE_MATXD_FUNCTIONS
          dmat_mul(apb, evecs, 2, 2, 2, W1);
#endif
          W1[0][0] /= nu[i][0];
          W1[0][1] /= nu[i][1];
          W1[1][0] /= nu[i][0];
          W1[1][1] /= nu[i][1];

          X_p[i][0][0] = (evecs[0][0] + W1[0][0]) / 2.;
          X_p[i][0][1] = (evecs[0][1] + W1[0][1]) / 2.;
          X_p[i][1][0] = (evecs[1][0] + W1[1][0]) / 2.;
          X_p[i][1][1] = (evecs[1][1] + W1[1][1]) / 2.;

          X_m[i][0][0] = (evecs[0][0] - W1[0][0]) / 2.;
          X_m[i][0][1] = (evecs[0][1] - W1[0][1]) / 2.;
          X_m[i][1][0] = (evecs[1][0] - W1[1][0]) / 2.;
          X_m[i][1][1] = (evecs[1][1] - W1[1][1]) / 2.;


          /*--------------------------------------------------------------------
           *
           *------------------------------------------------------------------*/
if (solar) {
          a1 = as_0[i] * as_0[i];

          gamma[0][0] -= a1;
          gamma[1][1] -= a1;

          a2 = gamma[0][0] * gamma[1][1] - gamma[0][1] * gamma[1][0];

          f[0][0] =  gamma[1][1] / a2;
          f[0][1] = -gamma[0][1] / a2;
          f[1][0] = -gamma[1][0] / a2;
          f[1][1] =  gamma[0][0] / a2;


          /*--------------------------------------------------------------------
           *
           *------------------------------------------------------------------*/
          for (j = 0; j < 2; ++j) {
               P_q0_pp[j] = 0.;
               P_q0_mp[j] = 0.;
          }

          for (l = i_four    ; l < 4; l += 2) {
               a1 = coef[i][0][l] * Y_0[l];
               for (j = 0; j < 2; ++j) {
                    P_q0_pp[j] += a1 * Y_q[j][l];
               }
          }

          for (l = i_four + 1; l < 4; l += 2) {
               a1 = coef[i][0][l] * Y_0[l];
               for (j = 0; j < 2; ++j) {
                    P_q0_mp[j] += a1 * Y_q[j][l];
               }
          }

          a1 = 2. - (i_four == 0 ? 1. : 0.);
          for (j = 0; j < 2; ++j) {
               a2 = (P_q0_pp[j] + P_q0_mp[j]) * a1;
               a3 = (P_q0_pp[j] - P_q0_mp[j]) * a1;
               P_q0_pp[j] = a2;
               P_q0_mp[j] = a3;
          }
/*
          build_phase_vecs_scalar(i_four, 4, 2, Y_q, Y_0, 0, coef[i][0], P_q0_pp, P_q0_mp);
*/

          /*--------------------------------------------------------------------
           *
           *------------------------------------------------------------------*/
          a1 = F_0 * omega[i] / (4. * PI);

          for (j = 0; j < 2; ++j) {
               a2 = a1 / qx[j];
               d[j] = a2 * (P_q0_pp[j] + P_q0_mp[j]);
               e[j] = a2 * (P_q0_pp[j] - P_q0_mp[j]);
          }
#ifdef USE_EXPANDED
          g[0] = -(amb[0][0] * e[0] + amb[0][1] * e[1]) - as_0[i] * d[0];
          g[1] = -(amb[1][0] * e[0] + amb[1][1] * e[1]) - as_0[i] * d[1];

          p[0] = f[0][0] * g[0] + f[0][1] * g[1];
          p[1] = f[1][0] * g[0] + f[1][1] * g[1];

          h[0] = apb[0][0] * p[0] + apb[0][1] * p[1] + e[0];
          h[1] = apb[1][0] * p[0] + apb[1][1] * p[1] + e[1];

          F_p[i][0] = (mu_0 * h[0] + p[0]) / 2.;
          F_p[i][1] = (mu_0 * h[1] + p[1]) / 2.;

          F_m[i][0] = F_p[i][0] - p[0];
          F_m[i][1] = F_p[i][1] - p[1];
#endif
#ifdef USE_MAT2D_FUNCTIONS
          mat2d_matvec(amb, e, v1);
          vec2d_scale(-as_0[i], d, v2);
          vec2d_sub(v2, v1, g);

          mat2d_matvec(f, g, p);

          mat2d_matvec(apb, p, v1);
          vec2d_add(v1, e, h);

          vec2d_scale(mu_0, h, v1);
          vec2d_add(v1, p, v1);
          vec2d_scale(.5, v1, F_p[i]);

          dvec_sub(F_p[i], p, F_m[i], 2);
#endif
#ifdef USE_MATXD_FUNCTIONS
          dm_v_mul(amb, e, 2, 2, v1);
          dvec_scale(-as_0[i], d, v2, 2);
          dvec_sub(v2, v1, g, 2);

          dm_v_mul(f, g, 2, 2, p);

          dm_v_mul(apb, p, 2, 2, v1);
          dvec_add(v1, e, h, 2);

          dvec_scale(mu_0, h, v1, 2);
          dvec_add(v1, p, v1, 2);
          dvec_scale(.5, v1, F_p[i], 2);

          dvec_sub(F_p[i], p, F_m[i], 2);
#endif

          /*--------------------------------------------------------------------
           *
           *------------------------------------------------------------------*/
          dvec_scale(btran[i], F_p[i], F_p[i], 2);
          dvec_scale(btran[i], F_m[i], F_m[i], 2);
}

          /*--------------------------------------------------------------------
           *
           *------------------------------------------------------------------*/
if (thermal) {
          w1[0][0] = alpha[0][0] * qx[0];
          w1[0][1] = alpha[0][1];
          w1[1][0] = alpha[1][0];
          w1[1][1] = alpha[1][1] * qx[1];

          w2[0][0] = beta [0][0] * qx[0];
          w2[0][1] = beta [0][1];
          w2[1][0] = beta [1][0];
          w2[1][1] = beta [1][1] * qx[1];

          a2 = gamma[0][0] * gamma[1][1] - gamma[0][1] * gamma[1][0];

          f[0][0] =  gamma[1][1] / a2;
          f[0][1] = -gamma[0][1] / a2;
          f[1][0] = -gamma[1][0] / a2;
          f[1][1] =  gamma[0][0] / a2;
}

          /*--------------------------------------------------------------------
           *
           *------------------------------------------------------------------*/
          Lambda[i][0] = exp(-nu[i][0] * ltau[i]);
          Lambda[i][1] = exp(-nu[i][1] * ltau[i]);
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     init_array2_d(A, n_comp, m_comp, 0.);


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
#ifndef USE_BANDED_SOLVER
     insert_matrix1_d(2, 2, 0, 0, A, -1., X_m[0]);
     insert_matrix2_d(2, 2, 0, 1, A, -1., Lambda[0], X_p[0]);
#else
     insert_matrix_band_storage1_d(2, 2, n_diags, n_diags, 0, 0, A, -1., X_m[0]);
     insert_matrix_band_storage2_d(2, 2, n_diags, n_diags, 0, 1, A, -1., Lambda[0], X_p[0]);
#endif

     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (n_layers > 1) {
#ifndef USE_BANDED_SOLVER
          insert_matrix2_d(2, 2, 1, 0, A,  1., Lambda[0], X_p[0]);
          insert_matrix1_d(2, 2, 1, 1, A,  1.,            X_m[0]);

          insert_matrix2_d(2, 2, 2, 0, A, -1., Lambda[0], X_m[0]);
          insert_matrix1_d(2, 2, 2, 1, A, -1.,            X_p[0]);
#else
          insert_matrix_band_storage2_d(2, 2, n_diags, n_diags, 1, 0, A,  1., Lambda[0], X_p[0]);
          insert_matrix_band_storage1_d(2, 2, n_diags, n_diags, 1, 1, A,  1.,            X_m[0]);

          insert_matrix_band_storage2_d(2, 2, n_diags, n_diags, 2, 0, A, -1., Lambda[0], X_m[0]);
          insert_matrix_band_storage1_d(2, 2, n_diags, n_diags, 2, 1, A, -1.,            X_p[0]);
#endif
          ii = 1;
          jj = 2;
          for (i = 1; i < n_layers - 1; ++i) {
#ifndef USE_BANDED_SOLVER
               insert_matrix1_d(2, 2, ii,     jj,     A, -1.,            X_p[i]);
               insert_matrix2_d(2, 2, ii,     jj + 1, A, -1., Lambda[i], X_m[i]);
               insert_matrix1_d(2, 2, ii + 1, jj,     A,  1.,            X_m[i]);
               insert_matrix2_d(2, 2, ii + 1, jj + 1, A,  1., Lambda[i], X_p[i]);

               insert_matrix2_d(2, 2, ii + 2, jj,     A,  1., Lambda[i], X_p[i]);
               insert_matrix1_d(2, 2, ii + 2, jj + 1, A,  1.,            X_m[i]);
               insert_matrix2_d(2, 2, ii + 3, jj,     A, -1., Lambda[i], X_m[i]);
               insert_matrix1_d(2, 2, ii + 3, jj + 1, A, -1.,            X_p[i]);
#else
               insert_matrix_band_storage1_d(2, 2, n_diags, n_diags, ii,     jj,     A, -1.,            X_p[i]);
               insert_matrix_band_storage2_d(2, 2, n_diags, n_diags, ii,     jj + 1, A, -1., Lambda[i], X_m[i]);
               insert_matrix_band_storage1_d(2, 2, n_diags, n_diags, ii + 1, jj,     A,  1.,            X_m[i]);
               insert_matrix_band_storage2_d(2, 2, n_diags, n_diags, ii + 1, jj + 1, A,  1., Lambda[i], X_p[i]);

               insert_matrix_band_storage2_d(2, 2, n_diags, n_diags, ii + 2, jj,     A,  1., Lambda[i], X_p[i]);
               insert_matrix_band_storage1_d(2, 2, n_diags, n_diags, ii + 2, jj + 1, A,  1.,            X_m[i]);
               insert_matrix_band_storage2_d(2, 2, n_diags, n_diags, ii + 3, jj,     A, -1., Lambda[i], X_m[i]);
               insert_matrix_band_storage1_d(2, 2, n_diags, n_diags, ii + 3, jj + 1, A, -1.,            X_p[i]);
#endif
               ii += 2;
               jj += 2;
          }
#ifndef USE_BANDED_SOLVER
          insert_matrix1_d(2, 2, ii,     jj,     A, -1.,            X_p[i]);
          insert_matrix2_d(2, 2, ii,     jj + 1, A, -1., Lambda[i], X_m[i]);

          insert_matrix1_d(2, 2, ii + 1, jj,     A,  1.,            X_m[i]);
          insert_matrix2_d(2, 2, ii + 1, jj + 1, A,  1., Lambda[i], X_p[i]);
#else
          insert_matrix_band_storage1_d(2, 2, n_diags, n_diags, ii,     jj,     A, -1.,            X_p[i]);
          insert_matrix_band_storage2_d(2, 2, n_diags, n_diags, ii,     jj + 1, A, -1., Lambda[i], X_m[i]);

          insert_matrix_band_storage1_d(2, 2, n_diags, n_diags, ii + 1, jj,     A,  1.,            X_m[i]);
          insert_matrix_band_storage2_d(2, 2, n_diags, n_diags, ii + 1, jj + 1, A,  1., Lambda[i], X_p[i]);
#endif
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (! surface) {
          pp1 = X_p[n_layers - 1];
          pp2 = X_m[n_layers - 1];
     }
     else {
          pp1 = w1;
          pp2 = w2;

          dmat_mul(Rs_qq, X_m[n_layers - 1], 2, 2, 2, w1);
          dmat_add(X_p[n_layers - 1], w1, w1, 2, 2);

          dmat_mul(Rs_qq, X_p[n_layers - 1], 2, 2, 2, w2);
          dmat_add(X_m[n_layers - 1], w2, w2, 2, 2);
     }
#ifndef USE_BANDED_SOLVER
     insert_matrix2_d(2, 2, 2 * n_layers - 1, 2 * n_layers - 2,     A, 1., Lambda[n_layers - 1], pp1);
     insert_matrix1_d(2, 2, 2 * n_layers - 1, 2 * n_layers - 2 + 1, A, 1.,                       pp2);
#else
     insert_matrix_band_storage2_d(2, 2, n_diags, n_diags, 2 * n_layers - 1, 2 * n_layers - 2,     A, 1., Lambda[n_layers - 1], pp1);
     insert_matrix_band_storage1_d(2, 2, n_diags, n_diags, 2 * n_layers - 1, 2 * n_layers - 2 + 1, A, 1.,                       pp2);
#endif

     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     insert_vector1_d(2, 0, b, -1., F_m[0]);

     ii = 1;
     for (i = 1; i < n_layers; ++i) {
          dvec_scale(atran[i - 1], F_p[i - 1], v1, 2);
          dvec_sub(F_p[i], v1, v1, 2);
          insert_vector1_d(2, ii,     b, 1., v1);

          dvec_scale(atran[i - 1], F_m[i - 1], v1, 2);
          dvec_sub(F_m[i], v1, v1, 2);
          insert_vector1_d(2, ii + 1, b, 1., v1);

          ii += 2;
     }

     if (! surface) {
          dvec_scale(-atran[n_layers - 1], F_p[n_layers - 1], v1, 2);

          insert_vector1_d(2, 2 * n_layers - 1, b, 1., v1);
     }
     else {
          dm_v_mul(Rs_qq, F_m[n_layers - 1], 2, 2, v1);
          dvec_sub(F_p[n_layers - 1], v1, v1, 2);
          dvec_scale(atran[n_layers - 1], v1, v1, 2);

          dvec_sub(In_p, v1, v1, 2);

          insert_vector1_d(2, 2 * n_layers - 1, b, 1., v1);
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
#ifndef USE_BANDED_SOLVER
     dmat_getrf(A, 4 * n_layers, 4 * n_layers, ipiv);

     dmat_getrs(A, &b, 4 * n_layers, 1, ipiv);
#else
     dgbtrf_(&n_comp, &n_comp, &n_diags, &n_diags, *A, &m_comp, ipiv, &info);
     if (info) {
          fprintf(stderr, "ERROR: dgbtrf() info = %d\n", info);
          exit(1);
     }

     dgbtrs_("N", &n_comp, &n_diags, &n_diags, &nrhs, *A, &m_comp, ipiv, b,
             &n_comp, &info);
     if (info) {
          fprintf(stderr, "ERROR: dgbtrs() info = %d\n", info);
          exit(1);
     }
#endif

     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
if (n_umus == 0) {
     if (! utau_output) {
          for (i = 0; i < n_ulevels; ++i) {
               ii = ulevels[i];

               if (ii != n_layers) {
                    iii = ii * 2 * 2;
#ifdef USE_EXPANDED
                    a1 = b[iii + 0    ];
                    a2 = b[iii + 0 + 2] * Lambda[ii][0];
                    a3 = b[iii + 1    ];
                    a4 = b[iii + 1 + 2] * Lambda[ii][1];

                    I_p[i][0] = a1 * X_p[ii][0][0] +
                                a2 * X_m[ii][0][0] +
                                a3 * X_p[ii][0][1] +
                                a4 * X_m[ii][0][1] + F_p[ii][0];
                    I_p[i][1] = a1 * X_p[ii][1][0] +
                                a2 * X_m[ii][1][0] +
                                a3 * X_p[ii][1][1] +
                                a4 * X_m[ii][1][1] + F_p[ii][1];
/*
                    for (j = 0; j < 2; ++j) {
                         a1 = 0.;
                         for (k = 0; k < 2; ++k) {
                              a1 += b[iii + k    ] *  X_p[ii][j][k];
                              a1 += b[iii + k + 2] *  X_m[ii][j][k] * Lambda[ii][k];
                         }

                         I_p[i][j] = a1 + F_p[ii][j];
                    }
*/
#endif
#ifdef USE_MAT2D_FUNCTIONS
#endif
#ifdef USE_MATXD_FUNCTIONS
                    dm_v_diag_mul(Lambda[ii], &b[iii + 2], v1, 2);

                    dm_v_mul(X_p[ii], &b[iii], 2, 2, v2);
                    dm_v_mul(X_m[ii], v1, 2, 2, v3);
                    dvec_add(v2, v3, v2, 2);
                    dvec_add(v2, F_p[ii], I_p[i], 2);
#endif
#ifdef USE_EXPANDED
                    I_m[i][0] = a1 * -X_m[ii][0][0] +
                                a2 * -X_p[ii][0][0] +
                                a3 * -X_m[ii][0][1] +
                                a4 * -X_p[ii][0][1] + F_m[ii][0];
                    I_m[i][1] = a1 * -X_m[ii][1][0] +
                                a2 * -X_p[ii][1][0] +
                                a3 * -X_m[ii][1][1] +
                                a4 * -X_p[ii][1][1] + F_m[ii][1];
/*
                    for (j = 0; j < 2; ++j) {
                         a1 = 0.;
                         for (k = 0; k < 2; ++k) {
                              a1 += b[iii + k    ] * -X_m[ii][j][k];
                              a1 += b[iii + k + 2] * -X_p[ii][j][k] * Lambda[ii][k];
                         }

                         I_m[i][j] = a1 + F_m[ii][j];
                    }
*/
#endif
#ifdef USE_MAT2D_FUNCTIONS
#endif
#ifdef USE_MATXD_FUNCTIONS
                    dm_v_mul(X_m[ii], &b[iii], 2, 2, v2);
                    dm_v_mul(X_p[ii], v1, 2, 2, v3);
                    dvec_sub(F_m[ii], v2, I_m[i], 2);
                    dvec_sub(I_m[i], v3, I_m[i], 2);
#endif
               }
               else {
                    ii = n_layers - 1;

                    iii = ii * 2 * 2;
#ifdef USE_EXPANDED
                    a1 = b[iii + 0    ] * Lambda[ii][0];
                    a2 = b[iii + 0 + 2];
                    a3 = b[iii + 1    ] * Lambda[ii][1];
                    a4 = b[iii + 1 + 2];

                    I_p[i][0] = a1 * X_p[ii][0][0] +
                                a2 * X_m[ii][0][0] +
                                a3 * X_p[ii][0][1] +
                                a4 * X_m[ii][0][1] + F_p[ii][0] * atran[ii];
                    I_p[i][1] = a1 * X_p[ii][1][0] +
                                a2 * X_m[ii][1][0] +
                                a3 * X_p[ii][1][1] +
                                a4 * X_m[ii][1][1] + F_p[ii][1] * atran[ii];
/*
                    for (j = 0; j < 2; ++j) {
                         a1 = 0.;
                         for (k = 0; k < 2; ++k) {
                              a1 += b[iii + k    ] *  X_p[ii][j][k] * Lambda[ii][k];
                              a1 += b[iii + k + 2] *  X_m[ii][j][k];
                         }

                         I_p[i][j] = a1 + F_p[ii][j] * atran[ii];
                    }
*/
#endif
#ifdef USE_MAT2D_FUNCTIONS
#endif
#ifdef USE_MATXD_FUNCTIONS
                    dm_v_diag_mul(Lambda[ii], &b[iii], v1, 2);

                    dm_v_mul(X_p[ii], v1, 2, 2, v2);
                    dm_v_mul(X_m[ii], &b[iii + 2], 2, 2, v3);
                    dvec_add(v2, v3, v2, 2);
                    dvec_scale(atran[ii], F_p[ii], v3, 2);
                    dvec_add(v2, v3, I_p[i], 2);
#endif
#ifdef USE_EXPANDED
                    I_m[i][0] = a1 * -X_m[ii][0][0] +
                                a2 * -X_p[ii][0][0] +
                                a3 * -X_m[ii][0][1] +
                                a4 * -X_p[ii][0][1] + F_m[ii][0] * atran[ii];
                    I_m[i][1] = a1 * -X_m[ii][1][0] +
                                a2 * -X_p[ii][1][0] +
                                a3 * -X_m[ii][1][1] +
                                a4 * -X_p[ii][1][1] + F_m[ii][1] * atran[ii];
/*
                    for (j = 0; j < 2; ++j) {
                         a1 = 0.;
                         for (k = 0; k < 2; ++k) {
                              a1 += b[iii + k    ] * -X_m[ii][j][k] * Lambda[ii][k];
                              a1 += b[iii + k + 2] * -X_p[ii][j][k];
                         }

                         I_m[i][j] = a1 + F_m[ii][j] * atran[ii];
                    }
*/
#endif
#ifdef USE_MAT2D_FUNCTIONS
#endif
#ifdef USE_MATXD_FUNCTIONS
                    dm_v_mul(X_m[ii], v1, 2, 2, v2);
                    dm_v_mul(X_p[ii], &b[iii + 2], 2, 2, v3);
                    dvec_scale(atran[ii], F_m[ii], v1, 2);
                    dvec_sub(v1, v2, I_m[i], 2);
                    dvec_sub(I_m[i], v3, I_m[i], 2);
#endif
               }
          }

          void calc_radiance_levels(int n_quad, int n_layers, int n_derivs,
                                    int n_ulevels, int *ulevels,
                                    double *ltau, double **ltau_l,
                                    double *atran, double **atran_l,
                                    double **nu, double ***X_p, double ***X_m,
                                    double **F_p, double **F_m,
                                    double **F0_p, double **F0_m, double **F1_p, double **F1_m,
                                    double ***nu_l, double ****X_p_l, double ****X_m_l,
                                    double ***F_p_l, double ***F_m_l,
                                    double ***F0_p_l, double ***F0_m_l, double ***F1_p_l, double ***F1_m_l,
                                    double *B, double **B_l,
                                    double **I_p, double **I_m, double ***I_p_l, double ***I_m_l,
                                    int thermal,
                                    uchar **derivs_layers, uchar **derivs_beam,
                                    save_tree_data save_tree, work_data work);

          calc_radiance_levels(2, n_layers, 0,
                               n_ulevels, ulevels,
                               ltau, ltau_l,
                               atran, atran_l,
                               nu, X_p, X_m,
                               F_p, F_m,
                               NULL, NULL, NULL, NULL,
                               NULL, NULL, NULL,
                               NULL, NULL,
                               NULL, NULL, NULL, NULL,
                               b, NULL,
                               I_p, I_m, I_p_l, I_m_l,
                               0,
                               NULL, NULL, save_tree, work);

     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     else {
          double lambda_p[2];
          double lambda_m[2];

          for (i = 0; i < n_ulevels; ++i) {
               ii = ulevels[i];

               for (j = 0; j < 2; ++j) {
                    lambda_p[j] = exp(-nu[ii][j] * (utaus[i] - 0.));
                    lambda_m[j] = exp(-nu[ii][j] * (ltau[ii] - utaus[i]));
               }

               uatran = exp(-utaus[i] * as_0[ii]);

               iii = ii * 2 * 2;
#ifdef USE_EXPANDED

               a1 = b[iii + 0    ] * lambda_p[0];
               a2 = b[iii + 0 + 2] * lambda_m[0];
               a3 = b[iii + 1    ] * lambda_p[1];
               a4 = b[iii + 1 + 2] * lambda_m[1];

               I_p[i][0] = a1 * X_p[ii][0][0] +
                           a2 * X_m[ii][0][0] +
                           a3 * X_p[ii][0][1] +
                           a4 * X_m[ii][0][1] + F_p[ii][0] * uatran;
               I_p[i][1] = a1 * X_p[ii][1][0] +
                           a2 * X_m[ii][1][0] +
                           a3 * X_p[ii][1][1] +
                           a4 * X_m[ii][1][1] + F_p[ii][1] * uatran;
/*
               for (j = 0; j < 2; ++j) {
                    a1 = 0.;
                    for (k = 0; k < 2; ++k) {
                         a1 += b[iii + k    ] *  X_p[ii][j][k] * lambda_p[k];
                         a1 += b[iii + k + 2] *  X_m[ii][j][k] * lambda_m[k];
                    }

                    I_p[i][j] = a1 + F_p[ii][j] * uatran;
               }
*/
#endif
#ifdef USE_MAT2D_FUNCTIONS
#endif
#ifdef USE_MATXD_FUNCTIONS
               dm_v_diag_mul(lambda_p, &b[iii], v1, 2);
               dm_v_mul(X_p[ii], v1, 2, 2, v2);
               dm_v_diag_mul(lambda_m, &b[iii + 2], v1, 2);
               dm_v_mul(X_m[ii], v1, 2, 2, v3);
               dvec_add(v2, v3, v1, 2);
               dvec_scale(uatran, F_p[ii], v2, 2);
               dvec_add(v1, v2, I_p[i], 2);
#endif
#ifdef USE_EXPANDED
               I_m[i][0] = a1 * -X_m[ii][0][0] +
                           a2 * -X_p[ii][0][0] +
                           a3 * -X_m[ii][0][1] +
                           a4 * -X_p[ii][0][1] + F_m[ii][0] * uatran;
               I_m[i][1] = a1 * -X_m[ii][1][0] +
                           a2 * -X_p[ii][1][0] +
                           a3 * -X_m[ii][1][1] +
                           a4 * -X_p[ii][1][1] + F_m[ii][1] * uatran;
/*
               for (j = 0; j < 2; ++j) {
                    a1 = 0.;
                    for (k = 0; k < 2; ++k) {
                         a1 += b[iii + k    ] * -X_m[ii][j][k] * lambda_p[k];
                         a1 += b[iii + k + 2] * -X_p[ii][j][k] * lambda_m[k];
                    }

                    I_m[i][j] = a1 + F_m[ii][j] * uatran;
               }
*/
#endif
#ifdef USE_MAT2D_FUNCTIONS
#endif
#ifdef USE_MATXD_FUNCTIONS
               dm_v_diag_mul(lambda_p, &b[iii], v1, 2);
               dm_v_mul(X_m[ii], v1, 2, 2, v2);
               dm_v_diag_mul(lambda_m, &b[iii + 2], v1, 2);
               dm_v_mul(X_p[ii], v1, 2, 2, v3);
               dvec_scale(exp(-utaus[i] * as_0[ii]), F_m[ii], v1, 2);
               dvec_sub(v1, v2, I_m[i], 2);
               dvec_sub(I_m[i], v3, I_m[i], 2);
#endif
          }
/*
          void calc_radiance_taus(int n_quad, int n_layers, int n_derivs,
                                  int n_ulevels, int *ulevels, double *utaus,
                                  double *ltau, double **ltau_l,
                                  double *as_0, double **as_0_l, double *atran, double **atran_l,
                                  double **nu, double ***X_p, double ***X_m,
                                  double **F_p, double **F_m,
                                  double ***nu_l, double ****X_p_l, double ****X_m_l,
                                  double ***F_p_l, double ***F_m_l,
                                  double *B, double **B_l,
                                  double **I_p, double **I_m, double ***I_p_l, double ***I_m_l,
                                  uchar **derivs_layers, uchar **derivs_beam,
                                  save_tree_data save_tree, work_data work);

          calc_radiance_taus(2, n_layers, 0,
                             n_ulevels, ulevels, utaus,
                             ltau, ltau_l,
                             as_0, as_0_l, atran, atran_l,
                             nu, X_p, X_m,
                             F_p, F_m,
                             NULL, NULL, NULL,
                             NULL, NULL,
                             b, NULL,
                             I_p, I_m, I_p_l, I_m_l,
                             NULL, NULL,
                             save_tree, work);
*/
     }
}

     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
if (n_umus != 0) {
#ifndef USE_ARRAYS_OF_POINTERS
     double Y_u[1][6];

     double P_uq_pp[n_layers][n_umus][2];
     double P_uq_mp[n_layers][n_umus][2];

     double P_u0_mm[n_layers][n_umus];
     double P_u0_pm[n_layers][n_umus];
#else
     double **Y_u;

     double ***P_uq_pp;
     double ***P_uq_mp;

     double **P_u0_mm;
     double **P_u0_pm;
#endif
     double **I_m2;


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     Y_u     = get_work_d2(&work, n_umus, 6);
#ifdef USE_ARRAYS_OF_POINTERS
     P_uq_pp = get_work_d3(&work, n_layers, n_umus, 2);
     P_uq_mp = get_work_d3(&work, n_layers, n_umus, 2);

     P_u0_mm = get_work_d2(&work, n_layers, n_umus);
     P_u0_pm = get_work_d2(&work, n_layers, n_umus);
#endif
     I_m2    = get_work_d2(&work, n_umus, 2);


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     for (i = 0; i < n_umus; ++i)
          calc_Y_x(i_four, umus[i], Y_u[i]);


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     for (i = 0; i < n_layers; ++i) {

          for (j = 0; j < n_umus; ++j) {
               for (k = 0; k < 2; ++k) {
                    P_uq_pp[i][j][k] = 0.;
                    P_uq_mp[i][j][k] = 0.;
               }
          }

          for (l = i_four    ; l < 4; l += 2) {
               for (j = 0; j < n_umus; ++j) {
                    a1 = coef[i][0][l] * Y_u[j][l];
                    for (k = 0; k < 2; ++k) {
                         P_uq_pp[i][j][k] += a1 * Y_q[k][l];
                    }
               }
          }

          for (l = i_four + 1; l < 4; l += 2) {
               for (j = 0; j < n_umus; ++j) {
                    a1 = coef[i][0][l] * Y_u[j][l];
                    for (k = 0; k < 2; ++k) {
                         P_uq_mp[i][j][k] += a1 * Y_q[k][l];
                    }
               }
          }

          a1 = 2. - (i_four == 0 ? 1. : 0.);
          for (j = 0; j < n_umus; ++j) {
               for (k = 0; k < 2; ++k) {
                    a2 = (P_uq_pp[i][j][k] + P_uq_mp[i][j][k]) * a1;
                    a3 = (P_uq_pp[i][j][k] - P_uq_mp[i][j][k]) * a1;
                    P_uq_pp[i][j][k] = a2;
                    P_uq_mp[i][j][k] = a3;
               }
          }
/*
          build_phase_mats_scalar(i_four, 4, n_umus, 2, Y_u, Y_q, 0, coef[i][0], P_uq_pp[i], P_uq_mp[i]);
*/

          /*--------------------------------------------------------------------
           *
           *------------------------------------------------------------------*/
          for (j = 0; j < n_umus; ++j) {
               P_u0_pm[i][j] = 0.;
               P_u0_mm[i][j] = 0.;
          }

          for (l = i_four    ; l < 4; l += 2) {
               a1 = coef[i][0][l] * Y_0[l];
               for (j = 0; j < n_umus; ++j) {
                    P_u0_pm[i][j] += a1 * Y_u[j][l];
               }
          }

          for (l = i_four + 1; l < 4; l += 2) {
               a1 = coef[i][0][l] * Y_0[l];
               for (j = 0; j < n_umus; ++j) {
                    P_u0_mm[i][j] += a1 * Y_u[j][l];
               }
          }

          a1 = 2. - (i_four == 0 ? 1. : 0.);
          for (j = 0; j < n_umus; ++j) {
               a2 = (P_u0_pm[i][j] + P_u0_mm[i][j]) * a1;
               a3 = (P_u0_pm[i][j] - P_u0_mm[i][j]) * a1;
               P_u0_pm[i][j] = a2;
               P_u0_mm[i][j] = a3;
          }
/*
          build_phase_vecs_scalar(i_four, 4, n_umus, Y_u, Y_0, 0, coef[i][0], P_u0_pm[i], P_u0_mm[i]);
*/
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     ii = n_layers - 1;

     iii = ii * 2 * 2;

     a1 = b[iii + 0    ] * Lambda[ii][0];
     a2 = b[iii + 0 + 2];
     a3 = b[iii + 1    ] * Lambda[ii][1];
     a4 = b[iii + 1 + 2];

     I_m2[0][0] = a1 * -X_m[ii][0][0] +
                  a2 * -X_p[ii][0][0] +
                  a3 * -X_m[ii][0][1] +
                  a4 * -X_p[ii][0][1] + F_m[ii][0] * atran[ii];
     I_m2[0][1] = a1 * -X_m[ii][1][0] +
                  a2 * -X_p[ii][1][0] +
                  a3 * -X_m[ii][1][1] +
                  a4 * -X_p[ii][1][1] + F_m[ii][1] * atran[ii];


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (upwelling)
          sfi_up(i_four,
                 2, 1, 0, n_layers,
                 1., qx, qw, F_0, mu_0,
                 n_ulevels, ulevels, utaus, n_umus, umus,
                 omega, omega_l, ltau, ltau_l,
                 0., NULL,
                 NULL, NULL,
                 btran, btran_l,
                 as_0, as_0_l, atran, atran_l,
                 NULL, NULL,
                 P_u0_pm,
                 P_uq_pp, P_uq_mp, P_uq_pp, P_uq_mp,
                 nu, X_p, X_m,
                 F_p, F_m,
                 NULL, NULL,
                 NULL, NULL, NULL, NULL,
                 NULL, NULL,
                 NULL,
                 NULL, NULL, NULL, NULL,
                 NULL, NULL, NULL,
                 NULL, NULL,
                 NULL, NULL,
                 Rs_u0, Rs_u0_l, Rs_uq, Rs_uq_l,
                 b, NULL,
                 0., NULL,
                 I_m2, NULL, I_p, NULL, 2,
                 1, 0, surface, 1, 0, utau_output,
                 NULL, work);

     if (downwelling)
          sfi_dn(i_four,
                 2, 1, 0, n_layers,
                 1., qx, qw, F_0,
                 n_ulevels, ulevels, utaus, n_umus, umus,
                 omega, omega_l, ltau, ltau_l,
                 NULL, NULL,
                 btran, btran_l,
                 as_0, as_0_l, atran, atran_l,
                 NULL, NULL,
                 P_u0_mm,
                 P_uq_pp, P_uq_mp, P_uq_pp, P_uq_mp,
                 nu, X_p, X_m,
                 F_p, F_m,
                 NULL, NULL,
                 NULL, NULL, NULL, NULL,
                 NULL, NULL,
                 NULL,
                 NULL, NULL, NULL, NULL,
                 NULL, NULL, NULL,
                 NULL, NULL,
                 NULL, NULL,
                 b, NULL,
                 0., NULL,
                 I_m2, NULL, I_m, NULL, 2,
                 1, 0, 1, 0, utau_output,
                 NULL, work);
}

}
