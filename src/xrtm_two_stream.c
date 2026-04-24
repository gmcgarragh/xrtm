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
#include "xrtm_derivs.h"
#include "xrtm_save_tree.h"
#include "xrtm_two_stream.h"
#include "xrtm_utility.h"


#define LAPACK_GENERAL_SOLVER	0
#define LAPACK_BANDED_SOLVER	1
#define PENTADIAGANOL_SOLVER	2

#define SOLVER LAPACK_GENERAL_SOLVER


#ifdef __cplusplus
extern "C" {
#endif

#if   SOLVER == LAPACK_GENERAL_SOLVER
void dgetrs_(const char *, int *, int *, double *, int *, int *, double *, int *, int *);
#elif   SOLVER == LAPACK_BANDED_SOLVER
void dgbtrf_(int *, int *, int *, int *, double *, int *, int *, int *);
void dgbtrs_(const char *, int *, int *, int *, int *, double *, int *, int *, double *, int *, int *);
#elif SOLVER == PENTADIAGANOL_SOLVER
void pentdag_(double *a, double *b, double *c, double *d, double *e, double *f, double *u, int *n);
void pentdag1_(double *a, double *b, double *c, double *d, double *e, double *p, double *q, int *n);
void pentdag2_(double *a, double *b, double *c, double *f, double *u, double *p, double *q, int *n);
#endif

#ifdef __cplusplus
}
#endif


/*******************************************************************************
 *
 ******************************************************************************/
void rtm_two_stream(int i_four,
                    int n_derivs, int n_layers,
                    double qx, double F_0, double mu_0,
                    int n_ulevels, int *ulevels, double *utaus,
                    int n_umus, double *umus,
                    double ***coef, double ****coef_l,
                    double *omega, double **omega_l, double *ltau, double **ltau_l,
                    double *btran, double **btran_l,
                    double *as_0, double **as_0_l, double *atran, double **atran_l,
                    double *Rs_q0, double **Rs_q0_l, double **Rs_qq, double ***Rs_qq_l,
                    double *Rs_u0, double **Rs_u0_l, double **Rs_uq, double ***Rs_uq_l,
                    double **I_p, double **I_m, double ***I_p_l, double ***I_m_l,
                    int add_single_scattering, int solar, int thermal, int surface,
                    int upwelling, int downwelling, int utau_output,
                    derivs_data *derivs, save_tree_data save_tree, work_data work) {

     int i;
     int ii;
     int iii;
     int j;
#if SOLVER == LAPACK_GENERAL_SOLVER
     int jj;
#endif
     int k;
     int kk;
#if SOLVER == LAPACK_GENERAL_SOLVER
     int m_comp;
     int n_comp;

     int info;
     int nrhs = 1;
#endif
#if SOLVER == LAPACK_BANDED_SOLVER
     int jj;

     int n_diags;
     int m_comp;
     int n_comp;

     int k2;

     int info;
     int nrhs = 1;
#endif
#if SOLVER == PENTADIAGANOL_SOLVER
     int n_comp;
#endif
     double a1;
     double a2;
     double a3;
     double a4;
     double a5;
     double a6;
     double a7;
     double a8;
     double a9;
     double a10;
     double a11;
     double a12;
     double a13;
     double a14;
     double a15;
     double a16;

     double a1_l;
     double a2_l;
     double a6_l;
     double a7_l;
     double a8_l;
     double a9_l;
     double a10_l;
     double a11_l;
     double a12_l;
     double a13_l;
     double a14_l;
     double a15_l;

     double solfac;

     double sqrt1;
     double sqrt2;

     double P_qq_pp;
     double P_qq_mp;

     double P_qq_pp_l;
     double P_qq_mp_l;

     double P_q0_mm;
     double P_q0_pm;

     double P_q0_mm_l;
     double P_q0_pm_l;

     double alpha;
     double beta;

     double alpha_l;
     double beta_l;

     double V_p;
     double V_m;

     double V_p_l;
     double V_m_l;

     double uatran;

     double *nu;

     double **nu_l;

     double *X_p;
     double *X_m;

     double **X_p_l;
     double **X_m_l;

     double *F_p;
     double *F_m;

     double **F_p_l;
     double **F_m_l;

     double *Lambda;

     double **Lambda_l;

     double *b;

     double **b_l;

     double *x;

     double **x_l;
#  if SOLVER == LAPACK_GENERAL_SOLVER || SOLVER == LAPACK_BANDED_SOLVER
     int *ipiv;

     double **A;
#elif SOLVER == PENTADIAGANOL_SOLVER
     int i_a;
     int i_b;
     int i_c;
     int i_d;
     int i_e;

     double *penta_a;
     double *penta_b;
     double *penta_c;
     double *penta_d;
     double *penta_e;
     double *penta_u;

     double *penta_p;
     double *penta_q;

     double **penta_u_l;
#endif
     double utau_l;
     double uatran_l;

     double *v1;
     double *v2;
     double *v3;
     double *v4;

     double **P_uq_pp;
     double **P_uq_mp;

     double ***P_uq_pp_l;
     double ***P_uq_mp_l;

     double **P_u0_mm;
     double **P_u0_pm;

     double ***P_u0_mm_l;
     double ***P_u0_pm_l;

     double *I_0;
     double **I_0_l;


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
#if SOLVER == LAPACK_GENERAL_SOLVER
     m_comp    = 2 * n_layers;
#elif SOLVER == LAPACK_BANDED_SOLVER
     n_diags   = 3 - 1;
     m_comp    = 3 * n_diags + 1;
#endif
     n_comp    = 2 * n_layers;


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     nu      = get_work_d1(&work, n_layers);

     X_p     = get_work_d1(&work, n_layers);
     X_m     = get_work_d1(&work, n_layers);

     F_p     = get_work_d1(&work, n_layers);
     F_m     = get_work_d1(&work, n_layers);

     Lambda  = get_work_d1(&work, n_layers);

     b       = get_work_d1(&work, n_comp);
#if SOLVER == LAPACK_GENERAL_SOLVER || SOLVER == LAPACK_BANDED_SOLVER
     A       = get_work_d2(&work, n_comp, m_comp);

     ipiv    = get_work_i1(&work, n_comp);
#elif SOLVER == PENTADIAGANOL_SOLVER
     penta_a = get_work_d1(&work, n_comp);
     penta_b = get_work_d1(&work, n_comp);
     penta_c = get_work_d1(&work, n_comp);
     penta_d = get_work_d1(&work, n_comp);
     penta_e = get_work_d1(&work, n_comp);
     penta_u = get_work_d1(&work, n_comp);

     penta_p = get_work_d1(&work, n_comp);
     penta_q = get_work_d1(&work, n_comp);
#endif
     if (n_derivs > 0) {
          b_l       = get_work_d2(&work, n_derivs, n_comp);
#if SOLVER == PENTADIAGANOL_SOLVER
          penta_u_l = get_work_d2(&work, n_derivs, n_comp);
#endif
     }

     if (flags_or2(derivs->layers, n_layers, n_derivs)) {
          nu_l      = get_work_d2(&work, n_layers, n_derivs);

          X_p_l     = get_work_d2(&work, n_layers, n_derivs);
          X_m_l     = get_work_d2(&work, n_layers, n_derivs);

          Lambda_l  = get_work_d2(&work, n_layers, n_derivs);
     }
/*
     if (flags_or2(derivs->beam, n_layers, n_derivs)) {
*/
          F_p_l     = get_work_d2(&work, n_layers, n_derivs);
          F_m_l     = get_work_d2(&work, n_layers, n_derivs);
/*
     }
*/

     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     solfac = F_0 * mu_0 / PI;


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     sqrt1 = sqrt(1. - qx * qx);
     sqrt2 = sqrt(1. - mu_0 * mu_0);

     a1 = qx * qx;
     a2 = 1. - a1;
     a3 = qx * mu_0;
     a4 = sqrt1 * sqrt2;
if (solar)
     a5 = F_0 / (4. * PI) / qx;

     for (i = 0; i < n_layers; ++i) {
          if (i_four == 0) {
               a6 = coef[i][0][1] * a1;
               P_qq_pp = (1. + a6);
               P_qq_mp = (1. - a6);

               a6 = coef[i][0][1] * a3;
               P_q0_mm = (1. + a6);
               P_q0_pm = (1. - a6);

               a6 = 1. / 2. / qx;
          }
          else {
               a6 = coef[i][0][1] * a2;
               P_qq_pp = (     a6); // / 2. * 2.;
               P_qq_mp = (     a6); // / 2. * 2.;

               a6 = coef[i][0][1] * a4;
               P_q0_mm = (     a6); //  / 2. * 2.;
               P_q0_pm = (     a6); //  / 2. * 2.;

               a6 = 1. / 4. / qx;
          }

          alpha = (-1. / qx + a6 * omega[i] * P_qq_pp);
          beta  =           - a6 * omega[i] * P_qq_mp ;

          a7 = alpha * alpha;
          a8 = beta  * beta;

          nu[i] = sqrt(a7 - a8);

          a9 = alpha + beta;

          X_p[i] = (nu[i] + a9) / (nu[i] - a9);
          X_m[i] = 1.;

          Lambda[i] = exp(-nu[i] * ltau[i]);

          F_p[i] = 0.;
          F_m[i] = 0.;
if (solar) {
          a10 = a5 * omega[i] * btran[i];

          a11 = a10 * P_q0_pm;
          a12 = a10 * P_q0_mm;

          a13 = as_0[i] + alpha;
          a14 = as_0[i] - alpha;

          a15 = a13 * a14 + a8;

          F_p[i] += (a11  * a13 + beta * a12) / a15;
          F_m[i] += (a12 * -a14 + beta * a11) / a15;
}
          a9  = nu[i] - a9;
          a16 = 2. * nu[i];

          for (j = 0; j < n_derivs; ++j) {
               if (derivs->layers[i][j]) {
                    if (i_four == 0) {
                         a6_l = coef_l[i][j][0][1] * a1;
                         P_qq_pp_l = (   + a6_l);
                         P_qq_mp_l = (   - a6_l);
                    }
                    else {
                         a6_l = coef_l[i][j][0][1] * a2;
                         P_qq_pp_l = (     a6_l); // / 2. * 2.;
                         P_qq_mp_l = (     a6_l); // / 2. * 2.;
                    }

                    alpha_l =   a6 * (omega_l[i][j] * P_qq_pp + omega[i] * P_qq_pp_l);
                    beta_l  = - a6 * (omega_l[i][j] * P_qq_mp + omega[i] * P_qq_mp_l);

                    a7_l = 2. * alpha_l * alpha;
                    a8_l = 2. * beta_l  * beta;

                    nu_l[i][j] = (a7_l - a8_l) / a16;

                    a9_l = alpha_l + beta_l;

                    X_p_l[i][j] = ((nu_l[i][j] + a9_l) - X_p[i] * (nu_l[i][j] - a9_l)) / a9;
                    X_m_l[i][j] = 0.;

                    Lambda_l[i][j] = (-nu_l[i][j] * ltau[i] - nu[i] * ltau_l[i][j]) * Lambda[i];

                    if (i_four == 0) {
                         a6_l = coef_l[i][j][0][1] * a3;
                         P_q0_mm_l = (   + a6_l);
                         P_q0_pm_l = (   - a6_l);
                    }
                    else {
                         a6_l = coef_l[i][j][0][1] * a4;
                         P_q0_mm_l = (     a6_l); //  / 2. * 2.;
                         P_q0_pm_l = (     a6_l); //  / 2. * 2.;
                    }
               }

               F_p_l[i][j] = 0.;
               F_m_l[i][j] = 0.;


               if (derivs->layers[i][j]) {
                    a10_l = a5 * (omega_l[i][j] * btran[i] + omega[i] * btran_l[i][j]);

                    a11_l = a10_l * P_q0_pm + a10 * P_q0_pm_l;
                    a12_l = a10_l * P_q0_mm + a10 * P_q0_mm_l;

                    a13_l = as_0_l[i][j] + alpha_l;
                    a14_l = as_0_l[i][j] - alpha_l;

                    a15_l = a13_l * a14 + a13 * a14_l + a8_l;

                    F_p_l[i][j] += ((a11_l *  a13 + a11  * a13_l + beta_l * a12 + beta * a12_l) - F_p[i] * a15_l) / a15;
                    F_m_l[i][j] += ((a12_l * -a14 + a12 * -a14_l + beta_l * a11 + beta * a11_l) - F_m[i] * a15_l) / a15;
               }
               else
               if (solar && derivs->beam[i][j]) {
                    a10_l = a5 *                             omega[i] * btran_l[i][j];

                    a11_l = a10_l * P_q0_pm;
                    a12_l = a10_l * P_q0_mm;

                    a13_l = as_0_l[i][j];
                    a14_l = as_0_l[i][j];

                    a15_l = a13_l * a14 + a13 * a14_l;

                    F_p_l[i][j] += ((a11_l *  a13 + a11  * a13_l +                beta * a12_l) - F_p[i] * a15_l) / a15;
                    F_m_l[i][j] += ((a12_l * -a14 + a12 * -a14_l +                beta * a11_l) - F_m[i] * a15_l) / a15;
               }
          }
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
#if SOLVER == LAPACK_GENERAL_SOLVER || SOLVER == LAPACK_BANDED_SOLVER
     init_array2_d(A, n_comp, m_comp, 0.);
#endif

     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
#if SOLVER == LAPACK_GENERAL_SOLVER
     A[0][0] = -X_m[0];
     A[0][1] = -X_p[0] * Lambda[0];


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (n_layers > 1) {
          A[1][0] =  X_p[0] * Lambda[0];
          A[1][1] =  X_m[0];
          A[2][0] = -X_m[0] * Lambda[0];
          A[2][1] = -X_p[0];

          ii = 1;
          jj = 2;
          for (i = 1; i < n_layers - 1; ++i) {
               A[ii    ][jj    ] = -X_p[i];
               A[ii    ][jj + 1] = -X_m[i] * Lambda[i];
               A[ii + 1][jj    ] =  X_m[i];
               A[ii + 1][jj + 1] =  X_p[i] * Lambda[i];

               A[ii + 2][jj    ] =  A[ii + 1][jj + 1];
               A[ii + 2][jj + 1] =  X_m[i];
               A[ii + 3][jj    ] =  A[ii    ][jj + 1];
               A[ii + 3][jj + 1] = -X_p[i];

               ii += 2;
               jj += 2;
          }

          A[ii    ][jj    ] = -X_p[i];
          A[ii    ][jj + 1] = -X_m[i] * Lambda[i];
          A[ii + 1][jj    ] =  X_m[i];
          A[ii + 1][jj + 1] =  X_p[i] * Lambda[i];
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (! surface) {
          V_p = X_p[n_layers - 1];
          V_m = X_m[n_layers - 1];
     }
     else {
          V_p = X_p[n_layers - 1] - Rs_qq[0][0] * -X_m[n_layers - 1];
          V_m = X_m[n_layers - 1] - Rs_qq[0][0] * -X_p[n_layers - 1];
     }

     A[2 * n_layers - 1][2 * n_layers - 2    ] = V_p * Lambda[n_layers - 1];
     A[2 * n_layers - 1][2 * n_layers - 2 + 1] = V_m;


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
#elif SOLVER == LAPACK_BANDED_SOLVER
     k2 = n_diags + n_diags + 1;

     ii = 0 + k2 - 1;
     jj = 0;
     A[jj    ][ii - (jj    )] = -X_m[0];
     A[jj + 1][ii - (jj + 1)] = -X_p[0] * Lambda[0];


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (n_layers > 1) {
          ii = 1 + k2 - 1;
          jj = 0;
          A[jj    ][ii       - (jj    )] =  X_p[0] * Lambda[0];
          A[jj + 1][ii       - (jj + 1)] =  X_m[0];
          A[jj    ][(ii + 1) - (jj    )] = -X_m[0] * Lambda[0];
          A[jj + 1][(ii + 1) - (jj + 1)] = -X_p[0];

          ii = 1 + k2 - 1;
          jj = 2;
          for (i = 1; i < n_layers - 1; ++i) {
               A[jj    ][(ii    ) - (jj    )] = -X_p[i];
               A[jj + 1][(ii    ) - (jj + 1)] = -X_m[i] * Lambda[i];
               A[jj    ][(ii + 1) - (jj    )] =  X_m[i];
               A[jj + 1][(ii + 1) - (jj + 1)] =  X_p[i] * Lambda[i];

               A[jj    ][(ii + 2) - (jj    )] =  A[jj + 1][(ii + 1) - (jj + 1)];
               A[jj + 1][(ii + 2) - (jj + 1)] =  X_m[i];
               A[jj    ][(ii + 3) - (jj    )] =  A[jj + 1][(ii    ) - (jj + 1)];
               A[jj + 1][(ii + 3) - (jj + 1)] = -X_p[i];

               ii += 2;
               jj += 2;
          }

          A[jj    ][(ii    ) - (jj    )] = -X_p[i];
          A[jj + 1][(ii    ) - (jj + 1)] = -X_m[i] * Lambda[i];
          A[jj    ][(ii + 1) - (jj    )] =  X_m[i];
          A[jj + 1][(ii + 1) - (jj + 1)] =  X_p[i] * Lambda[i];
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (! surface) {
          V_p = X_p[n_layers - 1];
          V_m = X_m[n_layers - 1];
     }
     else {
          V_p = X_p[n_layers - 1] - Rs_qq[0][0] * -X_m[n_layers - 1];
          V_m = X_m[n_layers - 1] - Rs_qq[0][0] * -X_p[n_layers - 1];
     }

     ii = 2 * n_layers - 1 + k2 - 1;
     jj = 2 * n_layers - 2;
     A[jj    ][ii - (jj    )] = V_p * Lambda[n_layers - 1];
     A[jj + 1][ii - (jj + 1)] = V_m;


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
#elif SOLVER == PENTADIAGANOL_SOLVER
     i_a = 2;
     i_b = 1;
     i_c = 0;
     i_d = 0;
     i_e = 0;

     penta_c[i_c++] = -X_m[0];
     penta_d[i_d++] = -X_p[0] * Lambda[0];
     penta_e[i_e++] = 0.;


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (n_layers > 1) {
          penta_b[i_b++] =  X_p[0] * Lambda[0];
          penta_c[i_c++] =  X_m[0];
          penta_a[i_a++] = -X_m[0] * Lambda[0];
          penta_b[i_b++] = -X_p[0];

          for (i = 1; i < n_layers - 1; ++i) {
               penta_d[i_d++] = -X_p[i];
               penta_e[i_e++] = -X_m[i] * Lambda[i];
               penta_c[i_c++] =  X_m[i];
               penta_d[i_d++] =  X_p[i] * Lambda[i];
               penta_e[i_e++] = 0.;

               penta_a[i_a++] = 0.;
               penta_b[i_b++] =  penta_d[i_d - 1];
               penta_c[i_c++] =  X_m[i];
               penta_a[i_a++] =  penta_e[i_e - 2];
               penta_b[i_b++] = -X_p[i];
          }

          penta_d[i_d++] = -X_p[i];
          penta_e[i_e++] = -X_m[i] * Lambda[i];
          penta_c[i_c++] =  X_m[i];
          penta_d[i_d++] =  X_p[i] * Lambda[i];
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (! surface) {
          V_p = X_p[n_layers - 1];
          V_m = X_m[n_layers - 1];
     }
     else {
          V_p = X_p[n_layers - 1] - Rs_qq[0][0] * -X_m[n_layers - 1];
          V_m = X_m[n_layers - 1] - Rs_qq[0][0] * -X_p[n_layers - 1];
     }

     penta_a[i_a++] = 0.;
     penta_b[i_b++] = V_p * Lambda[n_layers - 1];
     penta_c[i_c++] = V_m;
#endif

     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     dvec_zero(b, n_comp);

if (solar) {
     b[0] = -F_m[0];

     ii = 1;
     for (i = 1; i < n_layers; ++i) {
          b[ii    ] = F_p[i] - F_p[i - 1] * atran[i - 1];
          b[ii + 1] = F_m[i] - F_m[i - 1] * atran[i - 1];

          ii += 2;
     }

     if (! surface) {
          a1 = 0.;
          a2 = F_p[n_layers - 1];
     }
     else {
          a1 = solfac * btran[n_layers - 0] * Rs_q0[0];
          a2 = F_p[n_layers - 1] - Rs_qq[0][0] * F_m[n_layers - 1];
     }

     b[2 * n_layers - 1] = a1 - a2 * atran[n_layers - 1];
}

     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
#if SOLVER == LAPACK_GENERAL_SOLVER
     dmat_getrf(A, 2 * n_layers, 2 * n_layers, ipiv);
/*
     dmat_getrs(A, &b, 2 * n_layers, 1, ipiv);
*/
     dgetrs_("t", &n_comp, &nrhs, *A, &n_comp, ipiv, b, &n_comp, &info);

     if (info) {
          fprintf(stderr, "ERROR: dgetrs() info = %d\n", info);
          exit(1);
     }

     x = b;
#elif SOLVER == LAPACK_BANDED_SOLVER
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

     x = b;
#elif SOLVER == PENTADIAGANOL_SOLVER
if (0)
     pentdag_(penta_a, penta_b, penta_c, penta_d, penta_e, b, penta_u, &n_comp);
else {
     pentdag1_(penta_a, penta_b, penta_c, penta_d, penta_e, penta_p, penta_q, &n_comp);
     pentdag2_(penta_a, penta_b, penta_c, b, penta_u, penta_p, penta_q, &n_comp);
}
     x = penta_u;
#endif

     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     for (j = 0; j < n_derivs; ++j) {
          b_l[j][0] = 0.;
if (solar) {
          if (derivs->beam[0][j])
               b_l[j][0] += -F_m_l[0][j];
}
          if (derivs->layers[0][j])
               b_l[j][0] += X_m_l[0][j] * x[0] + (X_p_l[0][j] * Lambda[0] + X_p[0] * Lambda_l[0][j]) * x[1];

          ii = 1;
          kk = 0;
          for (i = 1; i < n_layers; ++i) {
if (0) {
               b_l[j][ii    ] = F_p_l[i][j] - F_p_l[i - 1][j] * atran[i - 1] - F_p[i - 1] * atran_l[i - 1][j]   -   (X_p_l[i - 1][j] * Lambda[i - 1] + X_p[i - 1] * Lambda_l[i - 1][j]) * x[kk]   -   X_m_l[i - 1][j] * x[kk + 1]   +   X_p_l[i][j] * x[kk + 2]   +   (X_m_l[i][j] * Lambda[i] + X_m[i] * Lambda_l[i][j]) * x[kk + 3];
               b_l[j][ii + 1] = F_m_l[i][j] - F_m_l[i - 1][j] * atran[i - 1] - F_m[i - 1] * atran_l[i - 1][j]   +   (X_m_l[i - 1][j] * Lambda[i - 1] + X_m[i - 1] * Lambda_l[i - 1][j]) * x[kk]   +   X_p_l[i - 1][j] * x[kk + 1]   -   X_m_l[i][j] * x[kk + 2]   -   (X_p_l[i][j] * Lambda[i] + X_p[i] * Lambda_l[i][j]) * x[kk + 3];
}
else {
               b_l[j][ii    ] = 0.;
               b_l[j][ii + 1] = 0.;
if (solar) {
               if (derivs->beam[i - 1][j]) {
                    b_l[j][ii    ] += - F_p_l[i - 1][j] * atran[i - 1] - F_p[i - 1] * atran_l[i - 1][j];
                    b_l[j][ii + 1] += - F_m_l[i - 1][j] * atran[i - 1] - F_m[i - 1] * atran_l[i - 1][j];
               }

               if (derivs->beam[i   ][j]) {
                    b_l[j][ii    ] +=   F_p_l[i    ][j];
                    b_l[j][ii + 1] +=   F_m_l[i    ][j];
               }
}
               if (derivs->layers[i - 1][j]) {
                    b_l[j][ii    ] += - (X_p_l[i - 1][j] * Lambda[i - 1] + X_p[i - 1] * Lambda_l[i - 1][j]) * x[kk] - X_m_l[i - 1][j] * x[kk + 1];
                    b_l[j][ii + 1] += + (X_m_l[i - 1][j] * Lambda[i - 1] + X_m[i - 1] * Lambda_l[i - 1][j]) * x[kk] + X_p_l[i - 1][j] * x[kk + 1];
               }

               if (derivs->layers[i    ][j]) {
                    b_l[j][ii    ] += +  X_p_l[i   ][j] * x[kk + 2] + (X_m_l[i    ][j] * Lambda[i] + X_m[i    ] * Lambda_l[i    ][j]) * x[kk + 3];
                    b_l[j][ii + 1] += -  X_m_l[i   ][j] * x[kk + 2] - (X_p_l[i    ][j] * Lambda[i] + X_p[i    ] * Lambda_l[i    ][j]) * x[kk + 3];
               }
}
               ii += 2;
               kk += 2;
          }

          if (! surface) {
if (0) {
               a1_l = 0.;
               a2_l = F_p_l[n_layers - 1][j];

               V_p_l = X_p_l[n_layers - 1][j];
               V_m_l = X_m_l[n_layers - 1][j];
}
else {
               b_l[j][2 * n_layers - 1] = 0.;
if (solar) {
               a1_l = 0.;
               a2_l = 0.;
               if (derivs->beam[n_layers - 1][j])
                    a2_l += F_p_l[n_layers - 1][j];

               b_l[j][2 * n_layers - 1] += a1_l - a2_l * atran[n_layers - 1] - a2 * atran_l[n_layers - 1][j];
}
               V_p_l = 0.;
               V_m_l = 0.;
               if (derivs->layers[n_layers - 1][j]) {
                    V_p_l = X_p_l[n_layers - 1][j];
                    V_m_l = X_m_l[n_layers - 1][j];

                    b_l[j][2 * n_layers - 1] += -(V_p_l * Lambda[n_layers - 1] + V_p * Lambda_l[n_layers - 1][j]) * x[2 * n_layers - 1 - 1] - V_m_l * x[2 * n_layers - 1];
               }
}
          }
          else {
if (0) {
               a1_l = solfac * (btran_l[n_layers - 0][j] * Rs_q0[0] + btran[n_layers - 0] * Rs_q0_l[j][0]);
               a2_l = F_p_l[n_layers - 1][j] - Rs_qq_l[j][0][0] * F_m[n_layers - 1] - Rs_qq[0][0] * F_m_l[n_layers - 1][j];

               V_p_l = X_p_l[n_layers - 1][j] - Rs_qq_l[j][0][0] * -X_m[n_layers - 1] - Rs_qq[0][0] * -X_m_l[n_layers - 1][j];
               V_m_l = X_m_l[n_layers - 1][j] - Rs_qq_l[j][0][0] * -X_p[n_layers - 1] - Rs_qq[0][0] * -X_p_l[n_layers - 1][j];
}
else {
               b_l[j][2 * n_layers - 1] = 0.;
if (solar) {
               a1_l = 0.;
               if (derivs->beam[n_layers - 0][j])
                    a1_l += btran_l[n_layers - 0][j] * Rs_q0[0];
               if (derivs->layers[n_layers - 0][j])
                    a1_l += btran[n_layers - 0] * Rs_q0_l[j][0];
               a1_l *= solfac;

               a2_l = 0.;
               if (derivs->beam[n_layers - 1][j])
                    a2_l += F_p_l[n_layers - 1][j] - Rs_qq     [0][0] * F_m_l[n_layers - 1][j];
               if (derivs->layers[n_layers - 0][j])
                    a2_l +=                        - Rs_qq_l[j][0][0] * F_m  [n_layers - 1];

               b_l[j][2 * n_layers - 1] += a1_l - a2_l * atran[n_layers - 1] - a2 * atran_l[n_layers - 1][j];
}
               a1_l  = 0.;
               V_p_l = 0.;
               V_m_l = 0.;
               if (derivs->layers[n_layers - 1][j]) {
                    a1_l   = V_p * Lambda_l[n_layers - 1][j];

                    V_p_l += X_p_l[n_layers - 1][j] - Rs_qq[0][0] * -X_m_l[n_layers - 1][j];
                    V_m_l += X_m_l[n_layers - 1][j] - Rs_qq[0][0] * -X_p_l[n_layers - 1][j];
               }
               if (derivs->layers[n_layers - 0][j]) {
                    V_p_l += - Rs_qq_l[j][0][0] * -X_m[n_layers - 1];
                    V_m_l += - Rs_qq_l[j][0][0] * -X_p[n_layers - 1];
               }

               b_l[j][2 * n_layers - 1] += -(V_p_l * Lambda[n_layers - 1] + a1_l) * x[2 * n_layers - 1 - 1] - V_m_l * x[2 * n_layers - 1];
}
          }
if (0)
          b_l[j][2 * n_layers - 1] = a1_l - a2_l * atran[n_layers - 1] - a2 * atran_l[n_layers - 1][j] - (V_p_l * Lambda[n_layers - 1] + V_p * Lambda_l[n_layers - 1][j]) * x[2 * n_layers - 1 - 1] - V_m_l * x[2 * n_layers - 1];


          /*--------------------------------------------------------------------
           *
           *------------------------------------------------------------------*/
#if SOLVER == LAPACK_GENERAL_SOLVER
/*
          dmat_getrs(A, &b_l[j], 2 * n_layers, 1, ipiv);
*/
          dgetrs_("t", &n_comp, &nrhs, *A, &n_comp, ipiv, b_l[j], &n_comp, &info);

          if (info) {
               fprintf(stderr, "ERROR: dgetrs() info = %d\n", info);
               exit(1);
          }

          x_l = b_l;
#elif SOLVER == LAPACK_BANDED_SOLVER
          dgbtrs_("N", &n_comp, &n_diags, &n_diags, &nrhs, *A, &m_comp, ipiv, b_l[j],
                  &n_comp, &info);
          if (info) {
               fprintf(stderr, "ERROR: dgbtrs() info = %d\n", info);
               exit(1);
          }

          x_l = b_l;
#elif SOLVER == PENTADIAGANOL_SOLVER
          pentdag2_(penta_a, penta_b, penta_c, b_l[j], penta_u_l[j], penta_p, penta_q, &n_comp);

          x_l = penta_u_l;
#endif
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (n_umus == 0) {
          if (! utau_output) {
               for (i = 0; i < n_ulevels; ++i) {
                    ii = ulevels[i];

                    if (ii != n_layers) {
                         iii = ii * 2 * 1;
if (1) {
                         a1 = x[iii    ];
                         a2 = x[iii + 1] * Lambda[ii];

                         I_p[i][0] = a1 *  X_p[ii] + a2 *  X_m[ii] + F_p[ii];
                         I_m[i][0] = a1 * -X_m[ii] + a2 * -X_p[ii] + F_m[ii];
     if (0) {
                         for (j = 0; j < n_derivs; ++j) {
                              a1_l = x_l[j][iii    ];
                              a2_l = x_l[j][iii + 1] * Lambda[ii] + x[iii + 1] * Lambda_l[ii][j];

                              I_p_l[i][j][0] = a1_l *  X_p[ii] + a1 *  X_p_l[ii][j]  +   a2_l *  X_m[ii] + a2 *  X_m_l[ii][j]   +   F_p_l[ii][j];
                              I_m_l[i][j][0] = a1_l * -X_m[ii] + a1 * -X_m_l[ii][j]  +   a2_l * -X_p[ii] + a2 * -X_p_l[ii][j]   +   F_m_l[ii][j];
                         }
     }
     else {
                         for (j = 0; j < n_derivs; ++j) {
                              a1_l  = x_l[j][iii    ];
                              a2_l  = x_l[j][iii + 1] * Lambda[ii];

                              if (derivs->layers[ii][j])
                                   a2_l += x[iii + 1] * Lambda_l[ii][j];

                              I_p_l[i][j][0]  = a1_l *  X_p[ii] + a2_l *  X_m[ii];
                              I_m_l[i][j][0]  = a1_l * -X_m[ii] + a2_l * -X_p[ii];

                              if (derivs->layers[ii][j]) {
                                   I_p_l[i][j][0] += a1 *  X_p_l[ii][j] + a2 *  X_m_l[ii][j];
                                   I_m_l[i][j][0] += a1 * -X_m_l[ii][j] + a2 * -X_p_l[ii][j];
                              }

                              if (derivs->beam[ii][j]) {
                                   I_p_l[i][j][0] += F_p_l[ii][j];
                                   I_m_l[i][j][0] += F_m_l[ii][j];
                              }
                         }
     }
}
else
                         calc_level_radiance(iii, n_derivs, 1., NULL, 1., Lambda[ii], X_p[ii], X_m[ii], F_p[ii], F_m[ii], NULL, Lambda_l[ii], X_p_l[ii], X_m_l[ii], F_p_l[ii], F_m_l[ii], x, x_l, I_p[i], I_m[i], I_p_l[i], I_m_l[i], derivs);
                    }
                    else {
                         ii = n_layers - 1;

                         iii = ii * 2 * 1;
if (1) {
                         a1 = x[iii    ] * Lambda[ii];
                         a2 = x[iii + 1];

                         I_p[i][0] = a1 *  X_p[ii] + a2 *  X_m[ii] + F_p[ii] * atran[ii];
                         I_m[i][0] = a1 * -X_m[ii] + a2 * -X_p[ii] + F_m[ii] * atran[ii];
     if (0) {
                         for (j = 0; j < n_derivs; ++j) {
                              a1_l = x_l[j][iii    ] * Lambda[ii] + x[iii    ] * Lambda_l[ii][j];
                              a2_l = x_l[j][iii + 1];

                              I_p_l[i][j][0] = a1_l *  X_p[ii] + a1 *  X_p_l[ii][j]   +   a2_l *  X_m[ii] + a2 *  X_m_l[ii][j]   +   F_p_l[ii][j] * atran[ii] + F_p[ii] * atran_l[ii][j];
                              I_m_l[i][j][0] = a1_l * -X_m[ii] + a1 * -X_m_l[ii][j]   +   a2_l * -X_p[ii] + a2 * -X_p_l[ii][j]   +   F_m_l[ii][j] * atran[ii] + F_m[ii] * atran_l[ii][j];
                         }
     }
     else {
                         for (j = 0; j < n_derivs; ++j) {
                              a1_l  = x_l[j][iii] * Lambda  [ii];
                              if (derivs->layers[ii][j])
                                   a1_l += x[iii] * Lambda_l[ii][j];
                              a2_l  = x_l[j][iii + 1];

                              I_p_l[i][j][0]  = a1_l *  X_p[ii] + a2_l *  X_m[ii];
                              I_m_l[i][j][0]  = a1_l * -X_m[ii] + a2_l * -X_p[ii];

                              if (derivs->layers[ii][j]) {
                                   I_p_l[i][j][0] += a1 *  X_p_l[ii][j] + a2 *  X_m_l[ii][j];
                                   I_m_l[i][j][0] += a1 * -X_m_l[ii][j] + a2 * -X_p_l[ii][j];
                              }

                              if (derivs->beam[ii][j]) {
                                   I_p_l[i][j][0] += F_p_l[ii][j] * atran[ii] + F_p[ii] * atran_l[ii][j];
                                   I_m_l[i][j][0] += F_m_l[ii][j] * atran[ii] + F_m[ii] * atran_l[ii][j];
                              }
                         }
     }
}
else
                         calc_level_radiance(iii, n_derivs, atran[ii], atran_l[ii], Lambda[ii], 1., X_p[ii], X_m[ii], F_p[ii], F_m[ii], Lambda_l[ii], NULL, X_p_l[ii], X_m_l[ii], F_p_l[ii], F_m_l[ii], x, x_l, I_p[i], I_m[i], I_p_l[i], I_m_l[i], derivs);
                    }
               }
          }


          /*--------------------------------------------------------------------
           *
           *------------------------------------------------------------------*/
          else {
               for (i = 0; i < n_ulevels; ++i) {
                    ii = ulevels[i];

                    iii = ii * 2 * 1;

                    uatran = exp(-utaus[i] * as_0[ii]);

                    a1 = x[iii]     * exp(-nu[ii] * (utaus[i] - 0.));
                    a2 = x[iii + 1] * exp(-nu[ii] * (ltau[ii] - utaus[i]));

                    I_p[i][0] = a1 *  X_p[ii] + a2 *  X_m[ii] + F_p[ii] * uatran;
                    I_m[i][0] = a1 * -X_m[ii] + a2 * -X_p[ii] + F_m[ii] * uatran;
if (0) {
                    for (j = 0; j < n_derivs; ++j) {
                         utau_l   =  utaus[i] * ltau_l[ii][j] / ltau[ii];

                         uatran_l = (-utau_l * as_0[ii] - utaus[i] * as_0_l[ii][j]) * uatran;

                         a1_l = x_l[j][iii]     * exp(-nu[ii] * (utaus[i] - 0.))       + x[iii]     * (-nu_l[ii][j] * (utaus[i] - 0.)       - nu[ii] * (utau_l        - 0.))     * exp(-nu[ii] * (utaus[i] - 0.));
                         a2_l = x_l[j][iii + 1] * exp(-nu[ii] * (ltau[ii] - utaus[i])) + x[iii + 1] * (-nu_l[ii][j] * (ltau[ii] - utaus[i]) - nu[ii] * (ltau_l[ii][j] - utau_l)) * exp(-nu[ii] * (ltau[ii] - utaus[i]));

                         I_p_l[i][j][0] = a1_l *  X_p[ii] + a1 *  X_p_l[ii][j]   +   a2_l *  X_m[ii] + a2 *  X_m_l[ii][j]   +   F_p_l[ii][j] * uatran + F_p[ii] * uatran_l;
                         I_m_l[i][j][0] = a1_l * -X_m[ii] + a1 * -X_m_l[ii][j]   +   a2_l * -X_p[ii] + a2 * -X_p_l[ii][j]   +   F_m_l[ii][j] * uatran + F_m[ii] * uatran_l;
                    }
}
else {
                    for (j = 0; j < n_derivs; ++j) {
                         utau_l   =  utaus[i] * ltau_l[ii][j] / ltau[ii];

                         uatran_l = (-utau_l * as_0[ii] - utaus[i] * as_0_l[ii][j]) * uatran;

                         a1_l = x_l[j][iii] * exp(-nu[ii] * (utaus[i] - 0.));

                         if (derivs->layers[ii][j])
                              a1_l += x[iii] * (-nu_l[ii][j] * (utaus[i] - 0.) - nu[ii] * (utau_l - 0.)) * exp(-nu[ii] * (utaus[i] - 0.));

                         a2_l = x_l[j][iii + 1] * exp(-nu[ii] * (ltau[ii] - utaus[i]));

                         if (derivs->layers[ii][j])
                              a2_l += x[iii + 1] * (-nu_l[ii][j] * (ltau[ii] - utaus[i]) - nu[ii] * (ltau_l[ii][j] - utau_l)) * exp(-nu[ii] * (ltau[ii] - utaus[i]));

                         I_p_l[i][j][0]  = a1_l *  X_p[ii] + a2_l *  X_m[ii];
                         I_m_l[i][j][0]  = a1_l * -X_m[ii] + a2_l * -X_p[ii];

                         if (derivs->layers[ii][j]) {
                              I_p_l[i][j][0] += a1 *  X_p_l[ii][j] + a2 *  X_m_l[ii][j];
                              I_m_l[i][j][0] += a1 * -X_m_l[ii][j] + a2 * -X_p_l[ii][j];
                         }

                         if (derivs->beam[ii][j]) {
                              I_p_l[i][j][0] += F_p_l[ii][j] * uatran + F_p[ii] * uatran_l;
                              I_m_l[i][j][0] += F_m_l[ii][j] * uatran + F_m[ii] * uatran_l;
                         }
                    }
}
               }
          }
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     else {
          v1      = get_work_d1(&work, n_umus);
          v2      = get_work_d1(&work, n_umus);
          v3      = get_work_d1(&work, n_umus);
          v4      = get_work_d1(&work, n_umus);

          P_uq_pp = get_work_d2(&work, n_layers, n_umus);
          P_uq_mp = get_work_d2(&work, n_layers, n_umus);

          P_u0_mm = get_work_d2(&work, n_layers, n_umus);
          P_u0_pm = get_work_d2(&work, n_layers, n_umus);

          I_0     = get_work_d1(&work, n_umus);

          if (n_derivs > 0) {
               P_uq_pp_l = get_work_d3(&work, n_layers, n_derivs, n_umus);
               P_uq_mp_l = get_work_d3(&work, n_layers, n_derivs, n_umus);

               P_u0_mm_l = get_work_d3(&work, n_layers, n_derivs, n_umus);
               P_u0_pm_l = get_work_d3(&work, n_layers, n_derivs, n_umus);

               I_0_l     = get_work_d2(&work, n_derivs, n_umus);
          }


          /*--------------------------------------------------------------------
           *
           *------------------------------------------------------------------*/
          for (i = 0; i < n_umus; ++i) {
               a1 = sqrt(1. - umus[i] * umus[i]);

               v1[i] = umus[i] * qx;
               v2[i] = a1 * sqrt1;
               v3[i] = umus[i] * mu_0;
               v4[i] = a1 * sqrt2;
          }

          for (i = 0; i < n_layers; ++i) {
               for (j = 0; j < n_umus; ++j) {
                    if (i_four == 0) {
                         a6 = coef[i][0][1] * v1[j];
                         P_uq_pp[i][j] = (1. + a6);
                         P_uq_mp[i][j] = (1. - a6);

                         a6 = coef[i][0][1] * v3[j];
                         P_u0_mm[i][j] = (1. + a6);
                         P_u0_pm[i][j] = (1. - a6);
                    }
                    else {
                         a6 = coef[i][0][1] * v2[j];
                         P_uq_pp[i][j] = (     a6); // / 2. * 2.;
                         P_uq_mp[i][j] = (     a6); // / 2. * 2.;

                         a6 = coef[i][0][1] * v4[j];
                         P_u0_mm[i][j] = (     a6); // / 2. * 2.;
                         P_u0_pm[i][j] = (     a6); // / 2. * 2.;
                    }
               }

               for (j = 0; j < n_derivs; ++j) {
                    for (k = 0; k < n_umus; ++k) {
                         if (i_four == 0) {
                              a6 = coef_l[i][j][0][1] * v1[k];
                              P_uq_pp_l[i][j][k] = (   + a6);
                              P_uq_mp_l[i][j][k] = (   - a6);

                              a6 = coef_l[i][j][0][1] * v3[k];
                              P_u0_mm_l[i][j][k] = (   + a6);
                              P_u0_pm_l[i][j][k] = (   - a6);
                         }
                         else {
                              a6 = coef_l[i][j][0][1] * v2[k];
                              P_uq_pp_l[i][j][k] = (     a6); // / 2. * 2.;
                              P_uq_mp_l[i][j][k] = (     a6); // / 2. * 2.;

                              a6 = coef_l[i][j][0][1] * v4[k];
                              P_u0_mm_l[i][j][k] = (     a6); // / 2. * 2.;
                              P_u0_pm_l[i][j][k] = (     a6); // / 2. * 2.;
                         }
                    }
               }
          }


          /*--------------------------------------------------------------------
           *
           *------------------------------------------------------------------*/
if (upwelling) {
          for (i = 0; i < n_umus; ++i)
               I_0[i] = 0.;

          for (i = 0; i < n_derivs; ++i) {
               for (j = 0; j < n_umus; ++j)
                    I_0_l[i][j] = 0.;
          }

          i  =     n_layers - 1;
          ii = 2 * n_layers - 2;

          if (surface) {
               a1 = x[ii] * -X_m[i] * Lambda[i] + x[ii + 1] * -X_p[i] + F_m[i] * atran[i];

               for (j = 0; j < n_umus; ++j) {
                    I_0[j] = a1 * Rs_uq[j][0];
if (solar) {
     if (add_single_scattering)
                    I_0[j] += solfac * btran[n_layers - 0] * Rs_u0[j];
}
               }

               for (j = 0; j < n_derivs; ++j) {
if (0)
                    a1_l = (x_l[j][ii] * -X_m[i] + x[ii] * -X_m_l[i][j]) * Lambda[i] + x[ii] * -X_m[i] * Lambda_l[i][j] + x_l[j][ii + 1] * -X_p[i] + x[ii + 1] * -X_p_l[i][j] + F_m_l[i][j] * atran[i] + F_m[i] * atran_l[i][j];
else {
                    a1_l = x_l[j][ii] * -X_m[i];
                    if (derivs->layers[i][j])
                         a1_l += x[ii] * -X_m_l[i][j];
                    a1_l *= Lambda[i];

                    if (derivs->layers[i][j])
                         a1_l += x[ii] * -X_m[i] * Lambda_l[i][j];

                    a1_l += x_l[j][ii + 1] * -X_p[i];

                    if (derivs->layers[i][j])
                         a1_l += x[ii + 1] * -X_p_l[i][j];

                    if (solar &&
                        derivs->beam[i][j])
                         a1_l += F_m_l[i][j] * atran[i] + F_m[i] * atran_l[i][j];
}
                    for (k = 0; k < n_umus; ++k) {
if (0)
                         I_0_l[j][k] = solfac * (btran_l[n_layers - 0][j] * Rs_u0[k] + btran[n_layers - 0] * Rs_u0_l[j][k]) + (a1_l * Rs_uq[k][0] + a1 * Rs_uq_l[j][k][0]);
else {
                         I_0_l[j][k] = 0.;
if (add_single_scattering) {
                         if (solar &&
                             derivs->beam[n_layers - 0][j])
                              I_0_l[j][k] += btran_l[n_layers - 0][j] * Rs_u0[k];
                         if (derivs->layers[n_layers - 0][j])
                              I_0_l[j][k] += btran[n_layers - 0] * Rs_u0_l[j][k];
                         I_0_l[j][k] *= solfac;
}
                         I_0_l[j][k] += a1_l * Rs_uq[k][0];
                         if (derivs->layers[n_layers - 0][j])
                              I_0_l[j][k] += a1 * Rs_uq_l[j][k][0];
}
                    }
               }
          }

          sfi(i_four, 1, 1, n_derivs, n_layers, F_0,
              n_ulevels, ulevels, utaus, n_umus, umus,
              omega, omega_l, ltau, ltau_l,
              btran, btran_l,
              as_0, as_0_l, atran, atran_l,
              P_u0_pm, P_uq_pp, P_uq_mp,
              nu, X_p, X_m,
              F_p, F_m,
              P_u0_pm_l, P_uq_pp_l, P_uq_mp_l,
              nu_l, X_p_l, X_m_l,
              F_p_l, F_m_l,
              x, x_l,
              I_0, I_0_l, 1, I_p, I_p_l,
              add_single_scattering, solar, thermal, utau_output,
              derivs, work, 0);
}

          /*--------------------------------------------------------------------
           *
           *------------------------------------------------------------------*/
if (downwelling) {
          for (i = 0; i < n_umus; ++i)
               I_0[i] = 0.;

          for (i = 0; i < n_derivs; ++i) {
               for (j = 0; j < n_umus; ++j)
                    I_0_l[i][j] = 0.;
          }

          sfi(i_four, 1, 1, n_derivs, n_layers, F_0,
              n_ulevels, ulevels, utaus, n_umus, umus,
              omega, omega_l, ltau, ltau_l,
              btran, btran_l,
              as_0, as_0_l, atran, atran_l,
              P_u0_mm, P_uq_mp, P_uq_pp,
              nu, X_p, X_m,
              F_p, F_m,
              P_u0_mm_l, P_uq_mp_l, P_uq_pp_l,
              nu_l, X_p_l, X_m_l,
              F_p_l, F_m_l,
              x, x_l,
              I_0, I_0_l, 1, I_m, I_m_l,
              add_single_scattering, solar, thermal, utau_output,
              derivs, work, 1);
}
     }
}



void calc_level_radiance(int i_const, int n_derivs, double atran, double *atran_l, double Lambda1, double Lambda2, double X_p, double X_m, double F_p, double F_m, double *Lambda1_l, double *Lambda2_l, double *X_p_l, double *X_m_l, double *F_p_l, double *F_m_l, double *x, double **x_l, double *I_p, double *I_m, double **I_p_l, double **I_m_l, derivs_data *derivs) {

     int j;

     double a1;
     double a1_l;
     double a2;
     double a2_l;
     double a3_l;

     a1 = x[i_const + 0] * Lambda1;
     a2 = x[i_const + 1] * Lambda2;

     I_p[0] = a1 *  X_p + a2 *  X_m + F_p * atran;
     I_m[0] = a1 * -X_m + a2 * -X_p + F_m * atran;
if (0) {
     for (j = 0; j < n_derivs; ++j) {
          a1_l = x_l[j][i_const + 0] * Lambda1 + x[i_const + 0] * (Lambda1_l ? Lambda1_l[j] : 0.);
          a2_l = x_l[j][i_const + 1] * Lambda2 + x[i_const + 1] * (Lambda2_l ? Lambda2_l[j] : 0.);

          a3_l = atran_l ? atran_l[j] : 0.;

          I_p_l[j][0] = a1_l *  X_p + a1 *  X_p_l[j]   +   a2_l *  X_m + a2 *  X_m_l[j]   +   F_p_l[j] * atran + F_p * a3_l;
          I_m_l[j][0] = a1_l * -X_m + a1 * -X_m_l[j]   +   a2_l * -X_p + a2 * -X_p_l[j]   +   F_m_l[j] * atran + F_m * a3_l;
     }
}
else {
     for (j = 0; j < n_derivs; ++j) {
          a1_l = x_l[j][i_const    ] * Lambda1;
          if (derivs->layers[j])
               a1_l += x[i_const    ] * (Lambda1_l ? Lambda1_l[j] : 0.);

          a2_l = x_l[j][i_const + 1] * Lambda2;
          if (derivs->layers[j])
               a2_l += x[i_const + 1] * (Lambda2_l ? Lambda2_l[j] : 0.);

          a3_l = atran_l ? atran_l[j] : 0.;

          I_p_l[j][0] = a1_l *  X_p + a2_l *  X_m;
          I_m_l[j][0] = a1_l * -X_m + a2_l * -X_p;

          if (derivs->layers[j]) {
               I_p_l[j][0] += a1 *  X_p_l[j] + a2 *  X_m_l[j];
               I_m_l[j][0] += a1 * -X_m_l[j] + a2 * -X_p_l[j];
          }
          if (derivs->beam[j]) {
               I_p_l[j][0] += F_p_l[j] * atran + F_p * a3_l;
               I_m_l[j][0] += F_m_l[j] * atran + F_m * a3_l;
          }
     }
}
}



/*******************************************************************************
 *
 ******************************************************************************/
static void sfi_layer(int i_layer,
                      double atau, double btau, double ptau,
                      int i_four, int n_derivs, double F_0,
                      int n_umus, double *umus,
                      double omega, double *omega_l, double ltau, double *ltau_l,
                      double btran, double *btran_l,
                      double as_0, double *as_0_l, double atran, double *atran_l,
                      double *P_u0_pm, double *P_uq_pp, double *P_uq_mp,
                      double nu, double X_p, double X_m,
                      double F_p, double F_m,
                      double **P_u0_pm_l, double **P_uq_pp_l, double **P_uq_mp_l,
                      double *nu_l, double *X_p_l, double *X_m_l,
                      double *F_p_l, double *F_m_l,
                      double *b, double **b_l,
                      double *I_0, double **I_0_l,
                      int offset, double *I_u, double **I_u_l,
                      int add_single_scattering, int solar, int thermal, int utau_output,
                      derivs_data *derivs, work_data work, int flag) {

     int ii;
     int j;
     int k;

     double a1;
     double a1_l;
     double a2;

     double solfac;

     double q1;
     double q2;
     double q3;

     double u1[n_umus];
     double u2[n_umus];
     double u3[n_umus];
     double u4[n_umus];

     double AA;

     double e1[n_umus];
     double e2;
     double e3;
     double e4;
     double e5;

     double e1_l[n_derivs][n_umus];
     double e2_l;
     double e3_l;
     double e4_l;
     double e5_l;

     double ptau_l[n_derivs];

     double atau_l[n_derivs];
     double btau_l[n_derivs];

     double Fu_0[n_umus];
     double Du_0[n_umus];
     double Gu_0[n_umus];
     double Eu_0;

     double qq1;
     double qq2;
     double qq3;
     double qq4;

     double uq1[n_umus];
     double uq2[n_umus];

     double Xu_p[n_umus];
     double Xu_m[n_umus];
     double Yu_p[n_umus];
     double Yu_m[n_umus];

     double Du_p[n_umus];
     double Du_m[n_umus];
     double Eu_p;
     double Eu_m;


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     solfac = F_0 / (4. * PI);


     a1 = (1. + (i_four == 0 ? 1. : 0.)) / 4.;

     AA = a1 * 1.;


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     for (j = 0; j < n_umus; ++j) {
/*
          if (! flag)
               e1[j] = exp(-(ltau - ptau) / umus[j]);
          else
               e1[j] = exp(-        ptau  / umus[j]);
*/
          e1[j] = exp(-atau / umus[j]);
     }

     for (j = 0; j < n_derivs; ++j) {
          if (derivs->layers[i_layer][j]) {
               ptau_l[j] = ptau * ltau_l[j] / ltau;

               atau_l[j] = atau * ltau_l[j] / ltau;
               btau_l[j] = btau * ltau_l[j] / ltau;

               for (k = 0; k < n_umus; ++k)
                    e1_l[j][k] = -atau_l[j] / umus[k] * e1[k];
          }
     }


     for (j = 0; j < n_derivs; ++j) {
          for (k = 0; k < n_umus; ++k) {
               I_u_l[j][offset + k] = I_0_l[j][k] * e1[k];

               if (derivs->layers[i_layer][j])
                    I_u_l[j][offset + k] += I_0[k] * e1_l[j][k];
          }
     }

     for (j = 0; j < n_umus; ++j)
          I_u[offset + j] = I_0[j] * e1[j];


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
if (solar) {
     q1 = AA * F_p;
     q2 = AA * F_m;

     for (j = 0; j < n_umus; ++j) {
          u1[j] = P_uq_pp[j] * q1;
          u2[j] = P_uq_mp[j] * q2;
          u1[j] = u1[j] + u2[j];
     }

     a1 = btran * solfac;

     for (j = 0; j < n_umus; ++j) {
if (add_single_scattering)
          u1[j] += a1 * P_u0_pm[j];
          Fu_0[j] = omega * u1[j];
     }

     e2 = exp(-ptau * as_0);

     for (j = 0; j < n_umus; ++j) {
          if (! flag)
               Du_0[j] = (e2 - atran * e1[j]) / (1. + umus[j] * as_0);
          else
               Du_0[j] = (e2 -         e1[j]) / (1. - umus[j] * as_0);
     }

     for (j = 0; j < n_derivs; ++j) {
          if (derivs->beam[i_layer][j]) {
               q3 = AA * F_p_l[j];
               for (k = 0; k < n_umus; ++k)
                    u3[k] = P_uq_pp[k] * q3;

               q3 = AA * F_m_l[j];
               for (k = 0; k < n_umus; ++k) {
                    u4[k] = P_uq_mp[k] * q3;

                    u3[k] = u3[k] + u4[k];
               }

               a2 = btran_l[j] * solfac;

               for (k = 0; k < n_umus; ++k) {
if (! add_single_scattering)
                    Gu_0[k]  =                      omega * (                       u3[k]);
else
                    Gu_0[k]  =                      omega * (a2 * P_u0_pm     [k] + u3[k]);
               }
          }

          if (derivs->layers[i_layer][j]) {
               for (k = 0; k < n_umus; ++k) {
                    u2[k] = P_uq_pp_l[j][k] * q1;
                    u3[k] = P_uq_mp_l[j][k] * q2;
                    u2[k] = u2[k] + u3[k];

if (! add_single_scattering)
                    Gu_0[k] += omega_l[j] * u1[k] + omega * (                       u2[k]);
else
                    Gu_0[k] += omega_l[j] * u1[k] + omega * (a1 * P_u0_pm_l[j][k] + u2[k]);
               }

               e2_l = (-ptau_l[j] * as_0 - ptau * as_0_l[j]) * e2;
          }

          for (k = 0; k < n_umus; ++k) {
               if (derivs->beam[i_layer][j])
                    I_u_l[j][offset + k] += Gu_0[k] * Du_0[k];

               if (derivs->layers[i_layer][j]) {
                    if (! flag)
                         Eu_0 = ((e2_l - atran_l[j] * e1[k] - atran * e1_l[j][k]) - Du_0[k] * ( umus[k] * as_0_l[j])) / (1. + umus[k] * as_0);
                    else
                         Eu_0 = ((e2_l -                              e1_l[j][k]) - Du_0[k] * (-umus[k] * as_0_l[j])) / (1. - umus[k] * as_0);

                    I_u_l[j][offset + k] += Fu_0[k] * Eu_0;
               }
          }
     }

     for (j = 0; j < n_umus; ++j)
          I_u[offset + j] += Fu_0[j] * Du_0[j];
}

     /*-------------------------------------------------------------------------
      *
      *------------------------------------------------------------------------*/
     qq1 = AA * X_p;
     qq2 = AA * X_m;

     for (j = 0; j < n_umus; ++j) {
          uq1[j] = P_uq_pp[j] * qq1;
          uq2[j] = P_uq_mp[j] * qq2;
          Xu_p[j] = uq1[j] - uq2[j];

          uq1[j] = P_uq_pp[j] * qq2;
          uq2[j] = P_uq_mp[j] * qq1;
          Xu_m[j] = uq1[j] - uq2[j];
     }

     e3 = exp(-nu * btau);
     e4 = exp(-nu * ltau);
     e5 = exp(-nu * atau);

     ii = 2 * i_layer;

     for (j = 0; j < n_umus; ++j) {
          a1 = umus[j] * nu;

          if (! flag) {
               Du_p[j] = (e3 - e4 * e1[j]) / (1. + a1);
               Du_m[j] = (e5 -      e1[j]) / (1. - a1);
          }
          else {
/*
               Du_p[j] = (e5 -      e1[j]) / (1. - a1);
               Du_m[j] = (e3 - e4 * e1[j]) / (1. + a1);
*/
               Du_m[j] = (e3 - e4 * e1[j]) / (1. + a1);
               Du_p[j] = (e5 -      e1[j]) / (1. - a1);
          }
     }

     for (j = 0; j < n_umus; ++j) {
          u1[j] = b[ii] * Xu_p[j] * Du_p[j] + b[ii + 1] * Xu_m[j] * Du_m[j];
          u2[j] = omega * Xu_p[j] * Du_p[j];
          u3[j] = omega * Xu_m[j] * Du_m[j];
     }

     for (j = 0; j < n_derivs; ++j) {
          if (derivs->layers[i_layer][j]) {
               qq3 = AA * X_p_l[j];
               qq4 = AA * X_m_l[j];

               for (k = 0; k < n_umus; ++k) {
                    uq1 [k] = P_uq_pp_l[j][k] * qq1;
                    uq2 [k] = P_uq_mp_l[j][k] * qq2;
                    Yu_p[k] = uq1[k] - uq2[k];

                    uq1 [k] = P_uq_pp[k] * qq3;
                    Yu_p[k] = Yu_p[k] + uq1[k];
                    uq1 [k] = P_uq_mp[k] * qq4;
                    Yu_p[k] = Yu_p[k] - uq1[k];

                    uq1 [k] = P_uq_pp_l[j][k] * qq2;
                    uq2 [k] = P_uq_mp_l[j][k] * qq1;
                    Yu_m[k] = uq1[k] - uq2[k];

                    uq1 [k] = P_uq_pp[k] * qq4;
                    Yu_m[k] = Yu_m[k] + uq1[k];
                    uq1 [k] = P_uq_mp[k] * qq3;
                    Yu_m[k] = Yu_m[k] - uq1[k];
               }
/*
               e3_l = (-nu_l[j] * ptau  - nu * ptau_l[j])  * e3;
               e4_l = (-nu_l[j] * ltau - nu * ltau_l[j]) * e4;
               e5_l = (-nu_l[j] * (ltau - ptau) - nu * (ltau_l[j] - ptau_l[j])) * e5;
*/
               e3_l = (-nu_l[j] * btau - nu * btau_l[j]) * e3;
               e4_l = (-nu_l[j] * ltau - nu * ltau_l[j]) * e4;
               e5_l = (-nu_l[j] * atau - nu * atau_l[j]) * e5;
          }

          for (k = 0; k < n_umus; ++k) {
               I_u_l[j][offset + k] += b_l[j][ii] * u2[k] + b_l[j][ii + 1] * u3[k];

               if (derivs->layers[i_layer][j]) {
                    I_u_l[j][offset + k] += omega_l[j] * u1[k];

                    a1   = umus[k] * nu;

                    a1_l = umus[k] * nu_l[j];

                    if (! flag) {
                         Eu_p = ((e3_l - (e4_l * e1[k] + e4 * e1_l[j][k])) + Du_p[k] * -a1_l) / (1. + a1);
                         Eu_m = ((e5_l -                      e1_l[j][k])  + Du_m[k] *  a1_l) / (1. - a1);
                    }
                    else {
/*
                         Eu_p = ((e3_l -                      e1_l[j][k])  + Du_p[k] *  a1_l) / (1. - a1);
                         Eu_m = ((e5_l - (e4_l * e1[k] + e4 * e1_l[j][k])) + Du_m[k] * -a1_l) / (1. + a1);
*/
                         Eu_m = ((e3_l - (e4_l * e1[k] + e4 * e1_l[j][k])) + Du_m[k] * -a1_l) / (1. + a1);
                         Eu_p = ((e5_l -                      e1_l[j][k])  + Du_p[k] *  a1_l) / (1. - a1);
                    }

                    I_u_l[j][offset + k] += omega * (b[ii    ] * Yu_p[k] * Du_p[k] + b[ii    ] * Xu_p[k] * Eu_p +
                                                     b[ii + 1] * Yu_m[k] * Du_m[k] + b[ii + 1] * Xu_m[k] * Eu_m);
               }
          }
     }

     for (j = 0; j < n_umus; ++j)
          I_u[offset + j] += omega * u1[j];
}



/*******************************************************************************
 *
 ******************************************************************************/
void sfi(int i_four, int n_quad, int n_stokes, int n_derivs, int n_layers, double F_0,
         int n_ulevels, int *ulevels, double *utaus, int n_umus, double *umus,
         double *omega, double **omega_l, double *ltau, double **ltau_l,
         double *btran, double **btran_l,
         double *as_0, double **as_0_l, double *atran, double **atran_l,
         double **P_u0_pm, double **P_uq_pp, double **P_uq_mp,
         double *nu, double *X_p, double *X_m,
         double *F_p, double *F_m,
         double ***P_u0_pm_l, double ***P_uq_pp_l, double ***P_uq_mp_l,
         double **nu_l, double **X_p_l, double **X_m_l,
         double **F_p_l, double **F_m_l,
         double *B, double **B_l,
         double *I_0, double **I_0_l, int offset, double **I_u, double ***I_u_l,
         int add_single_scattering, int solar, int thermal, int utau_output,
         derivs_data *derivs, work_data work, int flag) {

     int i;
     int i1;
     int i2;
     int j;
     int k;

     int inc;

     int i_out_level;

     double atau;
     double btau;
     double ptau;

     double **I_u_l2;


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (! flag) {
          i1          = n_layers - 1;
          i2          = ulevels[0] - 1;
          inc         = -1;
          i_out_level = n_ulevels - 1;
     }
     else {
          i1          = 0;
          if (! utau_output)
               i2 = ulevels[n_ulevels - 1];
          else
               i2 = ulevels[n_ulevels - 1] + 1;
          inc         =  1;
          i_out_level = 0;
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if ((! flag && ! utau_output && i1 + 1 == ulevels[i_out_level]) || (flag && ! utau_output && i1 == ulevels[i_out_level]) || (! flag && utau_output && i1 + 1 == ulevels[i_out_level] && utaus[i_out_level] == 0.) || (flag && utau_output && i1 == ulevels[i_out_level] && utaus[i_out_level] == 0.)) {
          for (i = 0; i < n_umus; ++i) {
               I_u[i_out_level][offset + i] = I_0[i];
          }
          for (i = 0; i < n_derivs; ++i) {
               for (j = 0; j < n_umus; ++j) {
                    I_u_l[i_out_level][i][offset + j] = I_0_l[i][j];
               }
          }

          i_out_level += inc;
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     for (i = i1; i != i2; i += inc) {

          if ((! flag && utau_output && i == ulevels[i_out_level]) || (flag && utau_output && i == ulevels[i_out_level])) {
               while (i_out_level >= 0 && i_out_level < n_ulevels && i == ulevels[i_out_level]) {
                    ptau = utaus[i_out_level];

                    if (! flag) {
                         atau = ltau[i] - utaus[i_out_level];
                         btau =           utaus[i_out_level];
                    }
                    else {
                         atau =           utaus[i_out_level];
                         btau = ltau[i] - utaus[i_out_level];
                    }

                    if (n_derivs > 0)
                         I_u_l2 = I_u_l[i_out_level];

                    if (n_derivs == 0)
                         sfi_layer(i, atau, btau, ptau, i_four, n_derivs, n_umus, F_0, umus, omega[i], NULL, ltau[i], NULL, btran[i], NULL, as_0[i], NULL, atran[i], NULL, P_u0_pm[i], P_uq_pp[i], P_uq_mp[i], nu[i], X_p[i], X_m[i], F_p[i], F_m[i], NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, B, B_l, I_0, I_0_l, offset, I_u[i_out_level], I_u_l2, add_single_scattering, solar, thermal, utau_output, NULL, work, flag);
                    else
                         sfi_layer(i, atau, btau, ptau, i_four, n_derivs, n_umus, F_0, umus, omega[i], omega_l[i], ltau[i], ltau_l[i], btran[i], btran_l[i], as_0[i], as_0_l[i], atran[i], atran_l[i], P_u0_pm[i], P_uq_pp[i], P_uq_mp[i], nu[i], X_p[i], X_m[i], F_p[i], F_m[i], P_u0_pm_l[i], P_uq_pp_l[i], P_uq_mp_l[i], nu_l[i], X_p_l[i], X_m_l[i], F_p_l[i], F_m_l[i], B, B_l, I_0, I_0_l, offset, I_u[i_out_level], I_u_l2, add_single_scattering, solar, thermal, utau_output, derivs, work, flag);

                    i_out_level += inc;
               }
          }

          if (! flag)
               ptau = 0.;
          else
               ptau = ltau[i];

          atau = ltau[i];
          btau = 0.;

          if (n_derivs == 0)
               sfi_layer(i, atau, btau, ptau, i_four, n_derivs, F_0, n_umus, umus, omega[i], NULL, ltau[i], NULL, btran[i], NULL, as_0[i], NULL, atran[i], NULL, P_u0_pm[i], P_uq_pp[i], P_uq_mp[i], nu[i], X_p[i], X_m[i], F_p[i], F_m[i], NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, B, B_l, I_0, I_0_l, 0, I_0, I_0_l, add_single_scattering, solar, thermal, utau_output, NULL, work, flag);
          else
               sfi_layer(i, atau, btau, ptau, i_four, n_derivs, F_0, n_umus, umus, omega[i], omega_l[i], ltau[i], ltau_l[i], btran[i], btran_l[i], as_0[i], as_0_l[i], atran[i], atran_l[i], P_u0_pm[i], P_uq_pp[i], P_uq_mp[i], nu[i], X_p[i], X_m[i], F_p[i], F_m[i], P_u0_pm_l[i], P_uq_pp_l[i], P_uq_mp_l[i], nu_l[i], X_p_l[i], X_m_l[i], F_p_l[i], F_m_l[i], B, B_l, I_0, I_0_l, 0, I_0, I_0_l, add_single_scattering, solar, thermal, utau_output, derivs, work, flag);

          if ((! flag && ! utau_output && i == ulevels[i_out_level]) || (flag && ! utau_output && i + 1 == ulevels[i_out_level])) {
               for (j = 0; j < n_umus; ++j) {
                    I_u[i_out_level][offset + j] = I_0[j];
               }
               for (j = 0; j < n_derivs; ++j) {
                    for (k = 0; k < n_umus; ++k) {
                         I_u_l[i_out_level][j][offset + k] = I_0_l[j][k];
                    }
               }

               i_out_level += inc;
          }
     }
}
