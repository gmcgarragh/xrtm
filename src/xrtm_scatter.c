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

#include <rtutil_math.h>

#include "xrtm.h"
#include "xrtm_scatter.h"
#include "xrtm_scatter_a.h"
#include "xrtm_utility.h"


/*******************************************************************************
 *
 ******************************************************************************/
double scat_angle(double mu_0, double phi_0, double mu, double phi) {

     return mu*mu_0 + sqrt(1. - mu*mu) * sqrt(1. - mu_0*mu_0) *
                      cos((phi - phi_0)*D2R);
}



/*******************************************************************************
 *
 ******************************************************************************/
void phase_func(int n_coef, double *p, double *chi, double *P) {

     int i;

     *P = 0.;
     for (i = 0; i < n_coef; ++i)
          *P += p[i] * chi[i];
}



/*******************************************************************************
 *
 ******************************************************************************/
void build_phase_vecs_scalar(int i_four, int n_coef, int n_mus, double **Y1, double *Y2, int lda, double *chi, double *P_pp, double *P_mp) {

     int i;
     int k;

     double a;
     double b;
     double c;
/*
     a = 2. - (i_four == 0 ? 1. : 0.);
     for (i = 0; i < n_mus; ++i) {
          b = 1.;
          P_pp[i] = 0.;
          P_mp[i] = 0.;
          for (k = i_four; k < n_coef; ++k) {
               c = chi[k] * Y1[i][k] * Y2[k];
               P_pp[i] += c;
               P_mp[i] += c * b;
               b *= -1.;
          }

          P_pp[i] *= a;
          P_mp[i] *= a;
     }
*/
/*
     for (i = 0; i < n_mus; ++i) {
          P_pp[i] = 0.;
          P_mp[i] = 0.;
     }

     a = 1.;
     for (k = i_four; k < n_coef; ++k) {
          b = chi[k] * Y2[k];
          for (i = 0; i < n_mus; ++i) {
               c = b * Y1[i][k];
               P_pp[i] += c;
               P_mp[i] += c * a;
          }
          a *= -1.;
     }

     a = 2. - (i_four == 0 ? 1. : 0.);
     for (i = 0; i < n_mus; ++i) {
          P_pp[i] *= a;
          P_mp[i] *= a;
     }
*/
     for (i = 0; i < n_mus; ++i) {
          P_pp[i] = 0.;
          P_mp[i] = 0.;
     }

     for (k = i_four    ; k < n_coef; k += 2) {
          a = chi[k] * Y2[k];
          for (i = 0; i < n_mus; ++i) {
               P_pp[i] += a * Y1[i][k];
          }
     }

     for (k = i_four + 1; k < n_coef; k += 2) {
          a = chi[k] * Y2[k];
          for (i = 0; i < n_mus; ++i) {
               P_mp[i] += a * Y1[i][k];
          }
     }

     a = 2. - (i_four == 0 ? 1. : 0.);
     for (i = 0; i < n_mus; ++i) {
          b = (P_pp[i] + P_mp[i]) * a;
          c = (P_pp[i] - P_mp[i]) * a;
          P_pp[i] = b;
          P_mp[i] = c;
     }
}



/*******************************************************************************
 *
 ******************************************************************************/
void build_phase_mats_scalar(int i_four, int n_coef, int n_mus1, int n_mus2, double **Y1, double **Y2, int lda, double *chi, double **P_pp, double **P_mp) {

     int i;
     int j;
     int k;

     double a;
     double b;
     double c;
/*
     a = 2. - (i_four == 0 ? 1. : 0.);
     for (i = 0; i < n_mus1; ++i) {
          for (j = 0; j < n_mus2; ++j) {
               b = 1.;
               P_pp[i][j] = 0.;
               P_mp[i][j] = 0.;
               for (k = i_four; k < n_coef; ++k) {
                    c = chi[k] * Y1[i][k] * Y2[j][k];
                    P_pp[i][j] += c;
                    P_mp[i][j] += c * b;
                    b *= -1.;
               }
               P_pp[i][j] *= a;
               P_mp[i][j] *= a;
          }
     }
*/
/*
     for (i = 0; i < n_mus1; ++i) {
          for (j = 0; j < n_mus2; ++j) {
               P_pp[i][j] = 0.;
               P_mp[i][j] = 0.;
          }
     }

     a = 1.;
     for (k = i_four; k < n_coef; ++k) {
          for (i = 0; i < n_mus1; ++i) {
               b = chi[k] * Y1[i][k];
               for (j = 0; j < n_mus2; ++j) {
                    c = b * Y2[j][k];
                    P_pp[i][j] += c;
                    P_mp[i][j] += c * a;
               }
          }
          a *= -1.;
     }

     a = 2. - (i_four == 0 ? 1. : 0.);
     for (i = 0; i < n_mus1; ++i) {
          for (j = 0; j < n_mus2; ++j) {
               P_pp[i][j] *= a;
               P_mp[i][j] *= a;
          }
     }
*/
     for (i = 0; i < n_mus1; ++i) {
          for (j = 0; j < n_mus2; ++j) {
               P_pp[i][j] = 0.;
               P_mp[i][j] = 0.;
          }
     }

     for (k = i_four    ; k < n_coef; k += 2) {
          for (i = 0; i < n_mus1; ++i) {
               a = chi[k] * Y1[i][k];
               for (j = 0; j < n_mus2; ++j) {
                    P_pp[i][j] += a * Y2[j][k];
               }
          }
     }

     for (k = i_four + 1; k < n_coef; k += 2) {
          for (i = 0; i < n_mus1; ++i) {
               a = chi[k] * Y1[i][k];
               for (j = 0; j < n_mus2; ++j) {
                    P_mp[i][j] += a * Y2[j][k];
               }
          }
     }

     a = 2. - (i_four == 0 ? 1. : 0.);
     for (i = 0; i < n_mus1; ++i) {
          for (j = 0; j < n_mus2; ++j) {
               b = (P_pp[i][j] + P_mp[i][j]) * a;
               c = (P_pp[i][j] - P_mp[i][j]) * a;
               P_pp[i][j] = b;
               P_mp[i][j] = c;
          }
     }
}



/*******************************************************************************
 *
 ******************************************************************************/
void build_scat_vector_gc(int n_coef, int n_stokes, double **gsf, double **chi, double *F) {

     int i;
     int ii;

     ii = 0;

     F[ii] = 0.;
     for (i = 0; i < n_coef; ++i)
          F[ii] += chi[0][i] * gsf[0][i];
     ii++;

     if (n_stokes > 1) {
          F[ii] = 0.;

          for (i = 2; i < n_coef; ++i)
               F[ii] += chi[4][i] * gsf[1][i];
          ii++;
     }

     if (n_stokes > 2)
          F[ii++] = 0.;

     if (n_stokes > 3)
          F[ii++] = 0.;
}



void build_scat_matrix_gc(int n_coef, double **gsf, double **chi, double **F, int flag) {

     int i;
     int j;

     double a[6];

     for (i = 0; i < 6; ++i)
          a[i] = 0.;

     for (i = 0; i < 2; ++i) {
          a[0] += chi[0][i] * gsf[0][i];
          if (flag)
               a[3] += chi[3][i] * gsf[0][i];
     }

     for ( ; i < n_coef; ++i) {
          a[0] += chi[0][i] * gsf[0][i];

          if (flag) {
               a[3] += chi[3][i] * gsf[0][i];

               a[4] += chi[4][i] * gsf[1][i];
               a[5] += chi[5][i] * gsf[1][i];

               a[1] += (chi[1][i] + chi[2][i]) * gsf[2][i];
               a[2] += (chi[1][i] - chi[2][i]) * gsf[3][i];
          }
     }

     for (i = 0; i < 4; ++i) {
          for (j = 0; j < 4; ++j) {
               F[i][j] = 0.;
          }
     }

     F[0][0] = a[0];

     if (flag) {
          F[3][3] = a[3];

          F[0][1] = a[4];
          F[1][0] = a[4];
          F[2][3] = a[5];
          F[3][2] = a[5];

          F[1][1] = F[0][0];
          F[2][2] = F[3][3];
/*
          F[1][1] = (a[1] + a[2]) / 2.;
          F[2][2] = (a[1] - a[2]) / 2.;
*/
     }
}


/*******************************************************************************
 *
 ******************************************************************************/
void build_scat_vector_lc(int n_coef, int n_stokes, double *p, double **chi, double *F) {

     int i;
     int ii;

     ii = 0;

     F[ii] = 0.;
     for (i = 0; i < n_coef; ++i)
          F[ii] += chi[0][i] * p[i];
     ii++;

     if (n_stokes > 1) {
          F[ii] = 0.;
          for (i = 0; i < n_coef; ++i)
               F[ii] += chi[4][i] * p[i];
          ii++;
     }

     if (n_stokes > 2)
          F[ii++] = 0.;

     if (n_stokes > 3)
          F[ii++] = 0.;
}



void build_scat_matrix_lc(int n_coef, double *p, double **chi, double **F, int flag) {

     int i;
     int j;

     int row[6] = {0, 1, 2, 3, 0, 2};
     int col[6] = {0, 1, 2, 3, 1, 3};

     for (i = 0; i < 4; ++i) {
          for (j = 0; j < 4; ++j) {
               F[i][j] = 0.;
          }
     }

     for (j = 0; j < n_coef; ++j) {
          F[row[0]][col[0]] += chi[0][j] * p[j];
     }

     if (flag) {
          for (i = 1; i < 6; ++i) {
               for (j = 0; j < n_coef; ++j) {
                    F[row[i]][col[i]] += chi[i][j] * p[j];
               }
          }

          F[1][0] =  F[0][1];
          F[3][2] = -F[2][3];
     }
}



/*******************************************************************************
 *
 ******************************************************************************/
static void iteration1(int i_four, int i, int n_stokes, double mu, double ***A) {

     double a;
     double b;
     double c;
     double f;

     a = i + 1.;

     if (n_stokes >= 3) {
          b = sqrt((i + 3.) * (i - 1.));
          c = sqrt(i * i - 4.);
     }

     f = (2. * i + 1.) * mu;

     A[i+1][0][0] = (f * A[i][0][0] - i * A[i-1][0][0]) / a;
     if (n_stokes >= 3) {
          A[i+1][1][1] = (f * A[i][1][1] - c * A[i-1][1][1]) / b;
          A[i+1][2][2] = (f * A[i][2][2] - c * A[i-1][2][2]) / b;
     }
     if (n_stokes >= 4)
          A[i+1][3][3] = (f * A[i][3][3] - i * A[i-1][3][3]) / a;
}



static void iteration2(int i_four, int i, int n_stokes, double mu, double ***A) {

     double a;
     double b;
     double e;
     double f;
     double g;

     a = i + 1. - i_four;

     g = 2. * i + 1.;

     if (n_stokes >= 3) {
          b = a / (i + 1.) * sqrt((i + 3.) * (i - 1.));
          e = 2. * i_four * g / (i * (i + 1.));
     }

     f = g * mu;

     A[i+1][0][0] = (f * A[i][0][0]) / a;
     if (n_stokes >= 3) {
          A[i+1][1][1] = (f * A[i][1][1] + e * A[i][2][1]) / b;
          A[i+1][2][2] = (f * A[i][2][2] + e * A[i][1][2]) / b;
          A[i+1][1][2] = (f * A[i][1][2] + e * A[i][2][2]) / b;
/*
          A[i+1][2][1] = (f * A[i][2][1] + e * A[i][1][1]) / b;
*/
          A[i+1][2][1] = A[i+1][1][2];
     }
     if (n_stokes >= 4)
          A[i+1][3][3] = (f * A[i][3][3]) / a;
}



static void iteration3(int i_four, int i, int n_stokes, double mu, double ***A) {

     double a;
     double b;
     double c;
     double d;
     double e;
     double f;
     double g;

     a = i + 1. - i_four;

     c = i + i_four;

     g = 2. * i + 1.;

     if (n_stokes >= 3) {
          b = a / (i + 1.) * sqrt((i + 3.) * (i - 1.));
          d = c / i * sqrt(i * i - 4.);
          e = 2. * i_four * g / (i * (i + 1.));
     }

     f = g * mu;

     A[i+1][0][0] = (f * A[i][0][0] - c * A[i-1][0][0]) / a;
     if (n_stokes >= 3) {
          A[i+1][1][1] = (f * A[i][1][1] - d * A[i-1][1][1] + e * A[i][2][1]) / b;
          A[i+1][2][2] = (f * A[i][2][2] - d * A[i-1][2][2] + e * A[i][1][2]) / b;
          A[i+1][1][2] = (f * A[i][1][2] - d * A[i-1][1][2] + e * A[i][2][2]) / b;
/*
          A[i+1][2][1] = (f * A[i][2][1] - d * A[i-1][2][1] + e * A[i][1][1]) / b;
*/
          A[i+1][2][1] = A[i+1][1][2];
     }
     if (n_stokes >= 4)
          A[i+1][3][3] = (f * A[i][3][3] - c * A[i-1][3][3]) / a;
}



void basic_matrix(int i_four, int n_coef, int n_stokes, double mu, double ***A) {

     int i;

     double a;
     double b;
     double c;
     double d;
     double e;
     double f;

     double P[3];

     if (i_four >= n_coef)
          return;

     init_array3_d(A + i_four, n_coef - i_four, n_stokes, n_stokes, 0.);

     if (i_four == 0) {
          A[0][0][0] = 1.;
          if (n_stokes >= 4)
               A[0][3][3] = 1.;

          A[1][0][0] = mu;
          if (n_stokes >= 4)
               A[1][3][3] = mu;

          if (n_coef >= 3) {
               leg_poly(3, mu, P);

               a = sqrt(6.) / 4. * (1. - mu*mu);

               A[2][0][0] = P[2];
               if (n_stokes >= 3) {
                    A[2][1][1] = a;
                    A[2][2][2] = a;
               }
               if (n_stokes >= 4)
                    A[2][3][3] = P[2];

               for (i = 2; i < n_coef - 1; ++i)
                    iteration1(i_four, i, n_stokes, mu, A);
          }
     }
     else
     if (i_four == 1) {
          a = sqrt(1. - mu*mu);

          A[1][0][0] = a;
          if (n_stokes >= 4)
               A[1][3][3] = a;

          if (n_coef >= 3) {
               b =  3. * mu * a;
               a *= sqrt(6.);
               c = -.5 * mu * a;
               d = -.5 *      a;
               A[2][0][0] =  b;
               if (n_stokes >= 3) {
                    A[2][1][1] =  c;
                    A[2][1][2] = -d;
                    A[2][2][1] = -d;
                    A[2][2][2] =  c;
               }
               if (n_stokes >= 4)
                    A[2][3][3] =  b;

               for (i = 2; i < n_coef - 1; ++i)
                    iteration3(i_four, i, n_stokes, mu, A);
          }
     }

     else {
          a = mu * mu;

          /* reformulate to deal with division by zero when mu = 1. */
/*
          f = pow((1. - a), i_four / 2.);

          b = 1.;
          for (i = 3; i <= 2 * i_four - 1; i += 2)
               b *= i;
          b *= f;

          c = f * factorial(2 * i_four) / pow(2, i_four) /
                  sqrt(factorial(i_four - 2) * factorial(i_four + 2));

          f = c / (1. - a);

          d = f * (1. + a);

          e = f *  2. * mu;
*/
          f = pow((1. - a), (i_four - 2.) / 2.);

          b = 1.;
          for (i = 3; i <= 2 * i_four - 1; i += 2)
               b *= i;
          b *= f * (1. - a);

          /* reformulate to deal with huge factorial(2 * i_four) */
/*
          c = f * factorial(2 * i_four) / pow(2, i_four) /
                  sqrt(factorial(i_four - 2) * factorial(i_four + 2));
*/
          c =                             pow(2., i_four) *
                  sqrt(factorial(i_four - 2) * factorial(i_four + 2));

          c = f * factorial_evaluate(1. / c, 2, 2 * i_four);

          d = c * (1. + a);

          e = c *  2. * mu;

          A[i_four][0][0] =  b;
          if (n_stokes >= 3) {
               A[i_four][1][1] =  d;
               A[i_four][1][2] = -e;
               A[i_four][2][1] = -e;
               A[i_four][2][2] =  d;
          }
          if (n_stokes >= 4)
               A[i_four][3][3] =  b;

          if (i_four < n_coef - 1) {
               iteration2(i_four, i_four, n_stokes, mu, A);
               for (i = i_four + 1; i < n_coef - 1; ++i)
                    iteration3(i_four, i, n_stokes, mu, A);
          }
     }
}



/*******************************************************************************
 *
 ******************************************************************************/
static void gc_to_B(int i, int n_stokes, double **gc, double **B, double c) {

     B[0][0] = c *  gc[0][i];
     if (n_stokes >= 3) {
          B[1][1] = c *  gc[1][i];
          B[0][1] = c * -gc[4][i];
          B[1][0] = c * -gc[4][i];
          B[2][2] = c *  gc[2][i];
     }
     if (n_stokes >= 4) {
          B[3][3] = c *  gc[3][i];
          B[2][3] = c *  gc[5][i];
          B[3][2] = c * -gc[5][i];
     }
}



/*******************************************************************************
 *
 ******************************************************************************/
void build_phase_vecs_vector_gc(int i_four, int n_coef, int n_mus, int n_stokes, double *mu, double mu_0, double **gc, double *P_pp, double *P_mp, work_data work) {

     int i;
     int ii;
     int j;

     int n_mus_v;

     int n_stokes1;
     int n_stokes2;

     double a;
     double b;
     double c;

     double D[] = {1., 1., -1., -1.};

     double *v1;

     double **w1;

     double **l0;

     double **B;

     double ***A;
     double ***ll;

     if (n_stokes == 2) {
          n_stokes1 = 2;
          n_stokes2 = 3;
     }
     else {
          n_stokes1 = n_stokes;
          n_stokes2 = n_stokes;
     }

     n_mus_v = n_mus * n_stokes1;

     if (i_four >= n_coef) {
          dvec_zero(P_pp, n_mus_v);
          dvec_zero(P_mp, n_mus_v);

          return;
     }

     v1 = get_work_d1(&work, n_mus_v);

     w1 = get_work_d2(&work, n_mus_v, n_stokes2);

     l0 = get_work_d2(&work, n_coef, n_stokes2);

     B  = get_work_d2(&work, n_stokes2, n_stokes2);

     A  = get_work_d3(&work, n_coef, n_stokes2, n_stokes2);

     ll = get_work_d3(&work, n_coef, n_mus_v,    n_stokes2);

     ii = 0;
     for (i = 0; i < n_mus; ++i) {
          basic_matrix(i_four, n_coef, n_stokes2, mu[i], A);

          for (j = i_four; j < n_coef; ++j) {
               dmat_copy(&ll[j][ii], A[j], n_stokes1, n_stokes2);
          }

          ii += n_stokes1;
     }

     basic_matrix(i_four, n_coef, 1, mu_0, A);
     for (i = i_four; i < n_coef; ++i) {
          l0[i][0] = A[i][0][0];
          for (j = 1; j < n_stokes2; ++j) {
               l0[i][j] = 0.;
          }
     }

     dmat_zero(B, n_stokes2, n_stokes2);

     dvec_zero(P_pp, n_mus_v);
     dvec_zero(P_mp, n_mus_v);

     a = 1.;
     b = 1. / factorial(2 * i_four);
     c = 2. - (i_four == 0 ? 1. : 0.);
     for (i = i_four; i < n_coef; ++i) {
          gc_to_B(i, n_stokes2, gc, B, b * c);

          dmat_mul(ll[i], B, n_mus_v, n_stokes2, n_stokes2, w1);

          dm_v_mul(w1, l0[i], n_mus_v, n_stokes2, v1);
          dvec_add(P_pp, v1, P_pp, n_mus_v);

          dmat_mul_diag(w1, D, w1, n_mus_v, n_stokes2);
          dm_v_mul(w1, l0[i], n_mus_v, n_stokes2, v1);
          if (a > 0.)
               dvec_add(P_mp, v1, P_mp, n_mus_v);
          else
               dvec_sub(P_mp, v1, P_mp, n_mus_v);

          a *= -1.;

          b *= (i + 1 - i_four) / (i + 1. + i_four);
     }

     dm_v_mul_D_A(n_mus, n_stokes1, P_mp, P_mp);
/*
     dm_v_mul_D_A_2(n_mus, n_stokes1, P_mp, P_mp);
*/
}



/*******************************************************************************
 *
 ******************************************************************************/
void build_phase_mats_vector_gc(int i_four, int n_coef, int n_mus1, int n_mus2, int n_stokes, double *mu1, double *mu2, double **gc, double **P_pp, double **P_mp, double **P_mm, double **P_pm, int vector, work_data work) {

     int i;
     int ii;
     int j;

     int n_mus_v1;
     int n_mus_v2;

     int n_stokes1;
     int n_stokes2;

     double a;
     double b;
     double c;

     double D[] = {1., 1., -1., -1.};

     double **w1;
     double **w2;
     double **w3;

     double **B;

     double ***A;
     double ***ll;
     double ***ll2;

     if (n_stokes == 2) {
          n_stokes1 = 2;
          n_stokes2 = 3;
     }
     else {
          n_stokes1 = n_stokes;
          n_stokes2 = n_stokes;
     }

     n_mus_v1 = n_mus1 * n_stokes1;
     n_mus_v2 = n_mus2 * n_stokes1;

     if (i_four >= n_coef) {
          dmat_zero(P_pp, n_mus_v1, n_mus_v2);
          dmat_zero(P_mp, n_mus_v1, n_mus_v2);
          dmat_zero(P_mm, n_mus_v1, n_mus_v2);
          dmat_zero(P_pm, n_mus_v1, n_mus_v2);

          return;
     }

     w1  = get_work_d2(&work, n_mus_v1,  n_stokes2);
     w2  = get_work_d2(&work, n_stokes2, n_mus_v2);
     w3  = get_work_d2(&work, n_mus_v1,  n_mus_v2);

     B   = get_work_d2(&work, n_stokes2, n_stokes2);

     A   = get_work_d3(&work, n_coef, n_stokes2, n_stokes2);

     ll  = get_work_d3(&work, n_coef, n_mus_v1, n_stokes2);
     ll2 = get_work_d3(&work, n_coef, n_mus_v2, n_stokes2);

     ii = 0;
     for (i = 0; i < n_mus1; ++i) {
          basic_matrix(i_four, n_coef, n_stokes2, mu1[i], A);

          for (j = i_four; j < n_coef; ++j) {
               dmat_copy(&ll[j][ii], A[j], n_stokes1, n_stokes2);
          }

          ii += n_stokes1;
     }

     ii = 0;
     for (i = 0; i < n_mus2; ++i) {
          basic_matrix(i_four, n_coef, n_stokes2, mu2[i], A);

          for (j = i_four; j < n_coef; ++j) {
               dmat_copy(&ll2[j][ii], A[j], n_stokes1, n_stokes2);
          }

          ii += n_stokes1;
     }


     dmat_zero(B, n_stokes2, n_stokes2);

     dmat_zero(P_pp, n_mus_v1, n_mus_v2);
     dmat_zero(P_mp, n_mus_v1, n_mus_v2);

     a =  1.;
     b =  1. / factorial(2 * i_four);
     c =  2. - (i_four == 0 ? 1. : 0.);
     for (i = i_four; i < n_coef; ++i) {
          gc_to_B(i, n_stokes2, gc, B, b * c);

          dmat_mul(ll[i], B, n_mus_v1, n_stokes2, n_stokes2, w1);

          dmat_trans(ll2[i], w2, n_mus_v2, n_stokes2);
          dmat_mul(w1, w2, n_mus_v1, n_mus_v2, n_stokes2, w3);
          dmat_add(P_pp, w3, P_pp, n_mus_v1, n_mus_v2);

          dmat_mul_diag(w1, D, w1, n_mus_v1, n_stokes2);
          dmat_mul(w1, w2, n_mus_v1, n_mus_v2, n_stokes2, w3);
          if (a > 0.)
               dmat_add(P_mp, w3, P_mp, n_mus_v1, n_mus_v2);
          else
               dmat_sub(P_mp, w3, P_mp, n_mus_v1, n_mus_v2);

          a *= -1.;

          b *= (i + 1 - i_four) / (i + 1. + i_four);
     }

     dmat_mul_D_A2(n_mus1, n_stokes, n_mus2, n_stokes, P_mp, P_mp);
/*
     dmat_mul_D_A2_2(n_mus1, n_stokes, n_mus2, n_stokes, P_mp, P_mp);
*/
     if (vector)
          phase_matrix_symmetry3(n_mus1, n_stokes, n_mus2, n_stokes, P_pp, P_mp, P_mm, P_pm, 1.);
}



/*******************************************************************************
 * mu1 = incident direction
 * mu2 = scattering direction
 ******************************************************************************/
void scat_vector_rotate(int n_stokes, double mu1, double mu2, double MU, double d_phi, double *P1, double *P2) {

     /* Hovenier 1983 */

     double sin_scat;

     double sin_theta2;

     double P1_01;

     double cos_sig1;

     double cos_sig1_2;

     double C1;
     double S1;

     P2[0] = P1[0];
     if (n_stokes == 1)
          return;

     P1_01 = P1[1];

     sin_scat   = sqrt(1. - MU *MU );

     if (sin_scat == 0.)
          cos_sig1 = 0.;
     else {
          sin_theta2 = sqrt(1. - mu2*mu2);

          if (sin_theta2 == 0.)
/*
               cos_sig1 = 1.;
*/
               cos_sig1 = cos(D2R*d_phi);
          else
               cos_sig1 = (mu1 - mu2 * MU) / (sin_theta2 * sin_scat);
     }

     cos_sig1_2 = cos_sig1 * cos_sig1;

     C1 = 2. * cos_sig1_2 - 1.;

     if (1. - cos_sig1_2 < 0.)
          S1 = 0.;
     else
          S1 = 2.* sqrt(1. - cos_sig1_2) * cos_sig1;

     if (n_stokes >= 2)
          P2[1] =  P1_01*C1;
     if (n_stokes >= 3)
          P2[2] = -P1_01*S1;
     if (n_stokes >= 4)
          P2[3] =  0.0;
}



void scat_vector_rotate_a(int n_stokes, double mu1, double mu2, double MU, double d_phi, double *P1_a, double *P2_a) {

     /* Hovenier 1983 */

     double sin_scat;

     double sin_theta2;

     double cos_sig1;

     double cos_sig1_2;

     double C1;
     double S1;

     P2_a[0] = P1_a[0];
     if (n_stokes == 1)
          return;

     sin_scat   = sqrt(1. - MU *MU );

     if (sin_scat == 0.)
          cos_sig1 = 0.;
     else {
          sin_theta2 = sqrt(1. - mu2*mu2);

          if (sin_theta2 == 0.)
/*
               cos_sig1 = 1.;
*/
               cos_sig1 = cos(D2R*d_phi);
          else
               cos_sig1 = (mu1 - mu2 * MU) / (sin_theta2 * sin_scat);
     }

     cos_sig1_2 = cos_sig1 * cos_sig1;

     C1 = 2. * cos_sig1_2 - 1.;

     if (1. - cos_sig1_2 < 0.)
          S1 = 0.;
     else
          S1 = 2.* sqrt(1. - cos_sig1_2) * cos_sig1;

     if (n_stokes >= 2)
          P1_a[1]  =  P2_a[1]*C1;
     if (n_stokes >= 3)
          P1_a[1] += -P2_a[2]*S1;
     if (n_stokes >= 4)
          P2_a[3]  =  0.0;
}



/*******************************************************************************
 *
 ******************************************************************************/
void phase_matrix_symmetry2(int n_quad, int n_stokes,
                            double **P_pp, double **P_mp,
                            double **P_mm, double **P_pm, double f) {

     phase_matrix_symmetry_ldx2(n_quad, n_stokes,
                                *P_pp, *P_mp, *P_mm, *P_pm, n_quad * n_stokes, f);

}



void phase_matrix_symmetry_ldx2(int n_quad, int n_stokes,
                               double *P_pp, double *P_mp,
                               double *P_mm, double *P_pm, int ldx, double f) {

     phase_matrix_symmetry_ldx3(n_quad, n_stokes, n_quad, n_stokes, P_pp, P_mp, P_mm, P_pm, ldx, f);
}



/*******************************************************************************
 *
 ******************************************************************************/
void phase_matrix_symmetry3(int n_quad1, int n_stokes1, int n_quad2, int n_stokes2,
                            double **P_pp, double **P_mp,
                            double **P_mm, double **P_pm, double f) {

     phase_matrix_symmetry_ldx3(n_quad1, n_stokes1, n_quad2, n_stokes2,
                              *P_pp, *P_mp, *P_mm, *P_pm, n_quad2 * n_stokes2, f);
}



void phase_matrix_symmetry_ldx3(int n_quad1, int n_stokes1, int n_quad2, int n_stokes2,
                               double *P_pp, double *P_mp,
                               double *P_mm, double *P_pm, int ldx, double f) {

     int i;
     int ii;
     int j;
     int jj;

     int i_offset;
     int j_offset;

     int n_quad_v1;
     int n_quad_v2;

     double a;
     double b;

     n_quad_v1 = n_quad1 * n_stokes1;
     n_quad_v2 = n_quad2 * n_stokes2;

     if (n_stokes1 != 4 || n_stokes2 != 4) {
          ii = 0;
          for (i = 0; i < n_quad_v1; ++i) {
               a = f * (i % n_stokes1 >= 2 ? -1. : 1.);
               jj = ii;
               for (j = 0; j < n_quad_v2; ++j) {
                    b = a * (j % n_stokes2 >= 2 ? -1. : 1.);
                    P_mm[jj] = b * P_pp[jj];
                    P_pm[jj] = b * P_mp[jj];
                    jj++;
               }
               ii += ldx;
          }
     }
     else {
          if (f > 0.) {
               i_offset = 0;
               for (i = 0; i < n_quad_v1; i += 4) {
                    for (ii = i; ii < i + 2; ++ii) {
                         j_offset = i_offset;
                         for (j = 0; j < n_quad_v2; j += 4) {
                              P_mm[j_offset    ] =  P_pp[j_offset    ];
                              P_mm[j_offset + 1] =  P_pp[j_offset + 1];
                              P_mm[j_offset + 2] = -P_pp[j_offset + 2];
                              P_mm[j_offset + 3] = -P_pp[j_offset + 3];
          
                              P_pm[j_offset    ] =  P_mp[j_offset    ];
                              P_pm[j_offset + 1] =  P_mp[j_offset + 1];
                              P_pm[j_offset + 2] = -P_mp[j_offset + 2];
                              P_pm[j_offset + 3] = -P_mp[j_offset + 3];

                              j_offset += 4;
                         }

                         i_offset += ldx;
                    }
          
                    for (      ; ii < i + 4; ++ii) {
                         j_offset = i_offset;
                         for (j = 0; j < n_quad_v2; j += 4) {
                              P_mm[j_offset    ] = -P_pp[j_offset    ];
                              P_mm[j_offset + 1] = -P_pp[j_offset + 1];
                              P_mm[j_offset + 2] =  P_pp[j_offset + 2];
                              P_mm[j_offset + 3] =  P_pp[j_offset + 3];
          
                              P_pm[j_offset    ] = -P_mp[j_offset    ];
                              P_pm[j_offset + 1] = -P_mp[j_offset + 1];
                              P_pm[j_offset + 2] =  P_mp[j_offset + 2];
                              P_pm[j_offset + 3] =  P_mp[j_offset + 3];

                              j_offset += 4;
                         }

                         i_offset += ldx;
                    }
               }
          }
          else {
               i_offset = 0;
               for (i = 0; i < n_quad_v1; i += 4) {
                    for (ii = i; ii < i + 2; ++ii) {
                         j_offset = i_offset;
                         for (j = 0; j < n_quad_v2; j += 4) {
                              P_mm[j_offset    ] = -P_pp[j_offset    ];
                              P_mm[j_offset + 1] = -P_pp[j_offset + 1];
                              P_mm[j_offset + 2] =  P_pp[j_offset + 2];
                              P_mm[j_offset + 3] =  P_pp[j_offset + 3];
          
                              P_pm[j_offset    ] = -P_mp[j_offset    ];
                              P_pm[j_offset + 1] = -P_mp[j_offset + 1];
                              P_pm[j_offset + 2] =  P_mp[j_offset + 2];
                              P_pm[j_offset + 3] =  P_mp[j_offset + 3];

                              j_offset += 4;
                         }

                         i_offset += ldx;
                    }
          
                    for (      ; ii < i + 4; ++ii) {
                         j_offset = i_offset;
                         for (j = 0; j < n_quad_v2; j += 4) {
                              P_mm[j_offset    ] =  P_pp[j_offset    ];
                              P_mm[j_offset + 1] =  P_pp[j_offset + 1];
                              P_mm[j_offset + 2] = -P_pp[j_offset + 2];
                              P_mm[j_offset + 3] = -P_pp[j_offset + 3];
          
                              P_pm[j_offset    ] =  P_mp[j_offset    ];
                              P_pm[j_offset + 1] =  P_mp[j_offset + 1];
                              P_pm[j_offset + 2] = -P_mp[j_offset + 2];
                              P_pm[j_offset + 3] = -P_mp[j_offset + 3];

                              j_offset += 4;
                         }

                         i_offset += ldx;
                    }
               }
          }
     }
}



/*******************************************************************************
 *
 ******************************************************************************/
#ifdef INCLUDE_DEV_SOURCE
#include "../src2/xrtm_scatter2.c"
#endif
