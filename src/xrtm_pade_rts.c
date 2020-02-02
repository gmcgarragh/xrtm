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

#include <rtutil_math.h>

#include "xrtm.h"
#include "xrtm_doubling.h"
#include "xrtm_matrix.h"
#include "xrtm_pade_rts.h"
#include "xrtm_save_tree.h"
#include "xrtm_scatter.h"
#include "xrtm_source.h"
#include "xrtm_utility.h"


matd1d2 NONE;


/*******************************************************************************
 *
 ******************************************************************************/
#include "pade_table.h"
#include "pade_table_l.h"

void pade_get_s_and_r(double **r_, double **t_, int n_quad,
                      double ltau, int n_umus, double *umus, int *s, int *r, int derivs) {

     int i;
     int j;

     double a;

     double angle;

     angle = 1.;
     for (i = 0; i < n_umus; ++i) {
          if (umus[i] < 0.)
/*
          if (umus[i] < angle)
*/
               angle = umus[i];
     }

     angle = acos(angle)*R2D;

     a = dmat_p_one_norm_A(r_, t_, n_quad) / pade_table_norm0 * ltau;

     i = pade_table_ltau_shift +
         (int) ceil(log(a) / pade_table_ltau_scale);
     if (i < 0)
          i = 0;
     else
     if (i > pade_table_n_ltaus) {
          fprintf(stderr, "ERROR: pade_get_s_and_r(): ltau index too high\n");
          exit(1);
     }

     j = pade_table_n_angles -
         ((int) ceil(sin((90. - angle)*D2R) * pade_table_n_angles) - 1);

     if (j < 0)
          j = 0;
     if (j > pade_table_n_angles - 1) {
          fprintf(stderr, "ERROR: pade_get_s_and_r(): angle index too high\n");
          exit(1);
     }

     if (! derivs) {
          *s = pade_table_s  [i][j];
          *r = pade_table_r  [i][j];
     }
     else {
          *s = pade_table_s_l[i][j];
          *r = pade_table_r_l[i][j];
     }
}



/*******************************************************************************
 *
 ******************************************************************************/
static void pade_coefs_x(int p, int q, double f, double *coefs) {

     int i;
     int ii;

     double a;
     double g;

     coefs[0] = 1.;

     a = 1.;
     g = f;
     for (i = 0, ii = 1; i < p; ++i, ++ii) {
          a *= (double) (p - i) / (double) ((2 * p - i) * ii);
          coefs[ii] = a * g;

          g *= f;
     }
}


/*
static void pade_coefs_n(int p, int q, double *coefs) {

     pade_coefs_x(p, q,  1., coefs);
}



static void pade_coefs_d(int p, int q, double *coefs) {

     pade_coefs_x(q, p, -1., coefs);
}



static void pade_coefs_symmetry(int n, double *cn, double *cd) {

     int i;

     double f;

     f = 1.;
     for (i = 0; i <= n; ++i) {
          cd[i] =  cn[i] * 1.;
          f *= -1.;
     }
}
*/


/*******************************************************************************
 *
 ******************************************************************************/
matd1d2 get_work_d2_d(work_data *work, size_t m, size_t n, size_t d) {

     matd1d2 a;

     a.a = get_work_d2(work, m, n);

     if (d > 0) {
          a.d1 = get_work_d2(work, d, n);
          a.d2 = get_work_d1(work, d);
     }

     return a;
}


matd1d2 *get_work_d2_d2(work_data *work, size_t v, size_t m, size_t n, size_t d) {

     size_t i;

     matd1d2 *a;

     a = (matd1d2 *) get_work_x1(work, v, sizeof(matd1d2));

     for (i = 0; i < v; ++i) {
          a[i].a = get_work_d2(work, m, n);

          if (d > 0) {
               a[i].d1 = get_work_d2(work, d, n);
               a[i].d2 = get_work_d1(work, d);
          }
     }

     return a;
}



/*******************************************************************************
 *
 ******************************************************************************/
/*
static void matrix_to_matd1d2(double **a_gen, matd1d2 a_d1d2, int n_quad, int n_umus) {

}
*/


static void matd1d2_to_matrix(matd1d2 a_d1d2, double **a_gen, int n_quad, int n_umus) {

     int i;
     int ii;
     int j;

     int n_quad_x = n_quad + n_umus;

     for (i = 0; i < n_quad; ++i) {
          for (j = 0; j < n_quad; ++j)
               a_gen[i][j] = a_d1d2.a[i][j];

          for (     ; j < n_quad_x; ++j)
               a_gen[i][j] = 0.;
     }

     if (n_umus > 0) {
          ii = 0;
          for (i = n_quad; i < n_quad_x; ++i) {
               for (j = 0; j < n_quad; ++j)
                    a_gen[i][j] = a_d1d2.d1[ii][j];

               for (     ; j < n_quad_x; ++j)
                    a_gen[i][j] = 0.;

               a_gen[i][i] = a_d1d2.d2[ii];

               ii++;
          }
     }
}



/*******************************************************************************
 *
 ******************************************************************************/
static void dmat_zero_d(double **a1, double *a2, long n_quad, long n_umus) {

     dmat_zero(a1, n_umus, n_quad);
     dvec_zero(a2, n_umus);
}



static void dmat_zero_d2(matd1d2 a, long n_quad, long n_umus) {

     dmat_zero(a.a, n_quad, n_quad);

     if (n_umus > 0)
          dmat_zero_d(a.d1, a.d2, n_quad, n_umus);
}



static void dmat_init_d(double **a1, double *a2, double alpha, long n_quad, long n_umus) {

     dmat_init(a1, 0.,    n_umus, n_quad);
     dvec_init(a2, alpha, n_umus);
}



static void dmat_init_d2(matd1d2 a, double alpha, long n_quad, long n_umus) {

     dmat_init(a.a, alpha, n_quad, n_quad);

     if (n_umus > 0)
          dmat_init_d(a.d1, a.d2, alpha, n_quad, n_umus);
}



static void dmat_copy_d(double **b1, double *b2, double **a1, double *a2, long n_quad, long n_umus) {

     dmat_copy(b1, a1, n_umus, n_quad);
     dvec_copy(b2, a2, n_umus);
}



static void dmat_copy_d2(matd1d2 b, matd1d2 a, long n_quad, long n_umus) {

     dmat_copy(b.a, a.a, n_quad, n_quad);

     if (n_umus > 0)
          dmat_copy_d(b.d1, b.d2, a.d1, a.d2, n_quad, n_umus);
}



static void dmat_add_d(double **a1, double *a2, double **b1, double *b2, double **c1, double *c2, long n_quad, long n_umus) {

     dmat_add(a1, b1, c1, n_umus, n_quad);
     dvec_add(a2, b2, c2, n_umus);
}



static void dmat_add_d2(matd1d2 a, matd1d2 b, matd1d2 c, long n_quad, long n_umus) {

     dmat_add(a.a, b.a, c.a, n_quad, n_quad);

     if (n_umus > 0)
          dmat_add_d(a.d1, a.d2, b.d1, b.d2, c.d1, c.d2, n_quad, n_umus);
}



static void dmat_sub_d(double **a1, double *a2, double **b1, double *b2, double **c1, double *c2, long n_quad, long n_umus) {

     dmat_sub(a1, b1, c1, n_umus, n_quad);
     dvec_sub(a2, b2, c2, n_umus);
}


/*
static void dmat_sub_d2(matd1d2 a, matd1d2 b, matd1d2 c, long n_quad, long n_umus) {

     dmat_sub(a.a, b.a, c.a, n_quad, n_quad);

     if (n_umus > 0)
          dmat_sub_d(a.d1, a.d2, b.d1, b.d2, c.d1, c.d2, n_quad, n_umus);
}
*/


static void dmat_scale_d(double alpha, double **a1, double *a2, double **b1, double *b2, long n_quad, long n_umus) {

     dmat_scale(alpha, a1, b1, n_umus, n_quad);
     dvec_scale(alpha, a2, b2, n_umus);
}



static void dmat_scale_d2(double alpha, matd1d2 a, matd1d2 b, long n_quad, long n_umus) {

     dmat_scale(alpha, a.a, b.a, n_quad, n_quad);

     if (n_umus > 0)
          dmat_scale_d(alpha, a.d1, a.d2, b.d1, b.d2, n_quad, n_umus);
}



static void dmat_scale_add_d(double **a1, double *a2, double **b1, double *b2, double **c1, double *c2, double alpha, long n_quad, long n_umus) {

     dmat_scale_add(a1, b1, c1, alpha, n_umus, n_quad);
     dvec_scale_add(a2, b2, c2, alpha, n_umus);
}



static void dmat_scale_add_d2(matd1d2 a, matd1d2 b, matd1d2 c, double alpha, long n_quad, long n_umus) {

     dmat_scale_add(a.a, b.a, c.a, alpha, n_quad, n_quad);

     if (n_umus > 0)
          dmat_scale_add_d(a.d1, a.d2, b.d1, b.d2, c.d1, c.d2, alpha, n_quad, n_umus);
}



static void dmat_mul_d(double **a1, double *a2, double **b, double **b1, double *b2, long n_quad, long n_umus, double **c1, double *c2, double **w1) {

     int i;

     dmat_mul(a1, b,  n_umus, n_quad, n_quad, c1);
/*
     dmat_mul_ldx(*a1, n_quad, *b, n_quad + n_umus, n_umus, n_quad, n_quad, *c1, n_quad);
*/
     dmat_diag_mul(a2, b1, w1, n_umus, n_quad);
     dmat_add(c1, w1, c1, n_umus, n_quad);

     for (i = 0; i < n_umus; ++i)
          c2[i] = a2[i] * b2[i];
}



static void dmat_mul_d2(matd1d2 a, matd1d2 b, long n_quad, long n_umus, matd1d2 c, double **w1) {

     dmat_mul(a.a, b.a, n_quad, n_quad, n_quad, c.a);

     if (n_umus > 0)
          dmat_mul_d(a.d1, a.d2, b.a, b.d1, b.d2, n_quad, n_umus, c.d1, c.d2, w1);
}



static void dmat_inv_mul_d(double **a_inv_b, double **a1, double *a2, double **b, double **b1, double *b2, long n_quad, long n_umus, double **c1, double *c2, double **w1) {

     int i;

     dmat_mul(a1, a_inv_b,  n_umus, n_quad, n_quad, w1);
/*
     dmat_mul_ldx(*a1, n_quad, *a_inv_b, n_quad + n_umus, n_umus, n_quad, n_quad, *w1, n_quad);
*/
     dmat_sub(b1, w1, w1, n_umus, n_quad);

     dmat_dinv_mul(a2, w1, c1, n_umus, n_quad);

     for (i = 0; i < n_umus; ++i)
          c2[i] = b2[i] / a2[i];
}



/*******************************************************************************
 *
 ******************************************************************************/
static void dmat_trans_stokes(int n_quad, int n_stokes, double **a, double **b) {

     int i;
     int j;

     int n_quad_v = n_quad * n_stokes;

     double x;
     double y;

     if (n_stokes <= 2)
          dmat_trans(a, b, n_quad_v, n_quad_v);
     else {
          for (i = 0; i < n_quad_v; ++i) {
               x = (i % n_stokes == 3) ? -1. : 1.;
               for (j = 0; j < n_quad_v; ++j) {
                    y = x * ((j % n_stokes == 3) ? -1. : 1.);
                    b[j][i] = a[i][j] * y;
               }
          }
     }
}



/*******************************************************************************
 *
 ******************************************************************************/
static void dmat_sym_stokes(int n_quad, int n_stokes, double **a) {

     int i;
     int j;

     int n_quad_v = n_quad * n_stokes;

     double x;
     double y;

     if (n_stokes <= 2)
          dmat_sym(a, n_quad_v, n_quad_v);
     else {
          for (i = 0; i < n_quad_v; ++i) {
               x = (i % n_stokes == 3) ? -1. : 1.;
               for (j = 0; j < i + 1; ++j) {
                    y = x * ((j % n_stokes == 3) ? -1. : 1.);
                    a[i][j] = a[j][i] * y;
               }
          }
     }
}



/*******************************************************************************
 *
 ******************************************************************************/
static void add_u_v_term_3_d2(int n_quad, int n_umus,
                              matd1d2 a1, matd1d2 a2,
                              matd1d2 b1, matd1d2 b2,
                              matd1d2 c1, matd1d2 c2,
                              matd1d2 u1, matd1d2 u2, matd1d2 v1, matd1d2 v2,
                              matd1d2 w1, matd1d2 w2, matd1d2 x1, matd1d2 x2,
                              double s0, double s1, double s2, double s3, double **w3u1) {

     dmat_mul_d2(a1, b1, n_quad, n_umus, c1, w3u1);
     dmat_mul_d2(a2, b2, n_quad, n_umus, c2, w3u1);

     if (s0 != 0.) {
          dmat_scale_add_d2(u1, c1, u1, s0, n_quad, n_umus);
          dmat_scale_add_d2(u2, c2, u2, s0, n_quad, n_umus);
     }
     if (s1 != 0.) {
          dmat_scale_add_d2(v1, c1, v1, s1, n_quad, n_umus);
          dmat_scale_add_d2(v2, c2, v2, s1, n_quad, n_umus);
     }
     if (s2 != 0.) {
          dmat_scale_add_d2(w1, c1, w1, s2, n_quad, n_umus);
          dmat_scale_add_d2(w2, c2, w2, s2, n_quad, n_umus);
     }
     if (s3 != 0.) {
          dmat_scale_add_d2(x1, c1, x1, s3, n_quad, n_umus);
          dmat_scale_add_d2(x2, c2, x2, s3, n_quad, n_umus);
     }
}



static void add_u_v_term_l_3_d2(int n_quad, int n_umus,
                                matd1d2 ax1, matd1d2 ax2,
                                matd1d2 bx1, matd1d2 bx2,
                                matd1d2 ay1, matd1d2 ay2,
                                matd1d2 by1, matd1d2 by2,
                                matd1d2 cy1, matd1d2 cy2,
                                matd1d2 uy1, matd1d2 uy2, matd1d2 vy1, matd1d2 vy2,
                                matd1d2 wy1, matd1d2 wy2, matd1d2 xy1, matd1d2 xy2,
                                matd1d2 w1, matd1d2 w2, double **w3u1,
                                double s0, double s1, double s2, double s3) {

     dmat_mul_d2(ay1, bx1, n_quad, n_umus, w1, w3u1);
     dmat_mul_d2(ax1, by1, n_quad, n_umus, w2, w3u1);
     dmat_add_d2(w1, w2, cy1, n_quad, n_umus);

     dmat_mul_d2(ay2, bx2, n_quad, n_umus, w1, w3u1);
     dmat_mul_d2(ax2, by2, n_quad, n_umus, w2, w3u1);
     dmat_add_d2(w1, w2, cy2, n_quad, n_umus);

     if (s0 != 0.) {
          dmat_scale_add_d2(uy1, cy1, uy1, s0, n_quad, n_umus);
          dmat_scale_add_d2(uy2, cy2, uy2, s0, n_quad, n_umus);
     }
     if (s1 != 0.) {
          dmat_scale_add_d2(vy1, cy1, vy1, s1, n_quad, n_umus);
          dmat_scale_add_d2(vy2, cy2, vy2, s1, n_quad, n_umus);
     }
     if (s2 != 0.) {
          dmat_scale_add_d2(wy1, cy1, wy1, s2, n_quad, n_umus);
          dmat_scale_add_d2(wy2, cy2, wy2, s2, n_quad, n_umus);
     }
     if (s3 != 0.) {
          dmat_scale_add_d2(xy1, cy1, xy1, s3, n_quad, n_umus);
          dmat_scale_add_d2(xy2, cy2, xy2, s3, n_quad, n_umus);
     }
}



/*******************************************************************************
 *
 ******************************************************************************/
static void reduced_recover(matd1d2 n11, matd1d2 n12, matd1d2 n21, matd1d2 n22,
                            matd1d2 d11, matd1d2 d12, matd1d2 d21, matd1d2 d22,
                            matd1d2 N11, matd1d2 N12, matd1d2 N21, matd1d2 N22,
                            matd1d2 D11, matd1d2 D12, matd1d2 D21, matd1d2 D22,
                            int n_quad, int n_umus) {

     int i;
     int j;

     double a;
     double b;
     double c;
     double d;
#ifdef CRAP
     double a1;
     double b1;
     double c1;
     double d1;

     double a2;
     double b2;
     double c2;
     double d2;
#endif
     for (i = 0; i < n_quad; ++i) {
          for (j = 0; j < n_quad; ++j) {

               /*---------------------------------------------------------------
                *
                *-------------------------------------------------------------*/
#ifdef CRAP
               a = .5 * (n11.a[i][j] + n21.a[i][j]);
               b = .5 * (n11.a[i][j] - n21.a[i][j]);
               c = .5 * (n12.a[i][j] + n22.a[i][j]);
               d = .5 * (n12.a[i][j] - n22.a[i][j]);

               N11.a[i][j] = a + c;
               N12.a[i][j] = a - c;
               N21.a[i][j] = b + d;
               N22.a[i][j] = b - d;

               a = .5 * (d11.a[i][j] + d21.a[i][j]);
               b = .5 * (d11.a[i][j] - d21.a[i][j]);
               c = .5 * (d12.a[i][j] + d22.a[i][j]);
               d = .5 * (d12.a[i][j] - d22.a[i][j]);

               D11.a[i][j] = a + c;
               D12.a[i][j] = a - c;
               D21.a[i][j] = b + d;
               D22.a[i][j] = b - d;
#endif

               /*---------------------------------------------------------------
                *
                *-------------------------------------------------------------*/
               a = .5 * (n11.a[i][j] + n21.a[i][j]);
               b = .5 * (n11.a[i][j] - n21.a[i][j]);
               c = .5 * (n12.a[i][j] + n22.a[i][j]);
               d = .5 * (n12.a[i][j] - n22.a[i][j]);

               N11.a[i][j] = a + c;
               N12.a[i][j] = a - c;
               N21.a[i][j] = b + d;
               N22.a[i][j] = b - d;
/*
               D11.a[i][j] = b - d;
               D12.a[i][j] = b + d;
               D21.a[i][j] = a - c;
               D22.a[i][j] = a + c;
*/

               /*---------------------------------------------------------------
                *
                *-------------------------------------------------------------*/
#ifdef CRAP
               a = .5 * (d11.a[i][j] + d21.a[i][j]);
               b = .5 * (d11.a[i][j] - d21.a[i][j]);
               c = .5 * (d12.a[i][j] + d22.a[i][j]);
               d = .5 * (d12.a[i][j] - d22.a[i][j]);

               D11.a[i][j] = a + c;
               D12.a[i][j] = a - c;
               D21.a[i][j] = b + d;
               D22.a[i][j] = b - d;
/*
               N11.a[i][j] = b - d;
               N12.a[i][j] = b + d;
               N21.a[i][j] = a - c;
               N22.a[i][j] = a + c;
*/
#endif
          }
     }

#ifdef CRAP
     for (i = 0; i < n_quad; ++i) {
          for (j = i; j < n_quad; ++j) {
#ifdef CRAP
               a1 = .5 * (n11.a[i][j] + n21.a[i][j]);
               b1 = .5 * (n11.a[i][j] - n21.a[i][j]);
               c1 = .5 * (n12.a[j][i] + n11.a[j][i]);
               d1 = .5 * (n12.a[j][i] - n11.a[j][i]);

               a2 = .5 * (n11.a[j][i] + n21.a[j][i]);
               b2 = .5 * (n11.a[j][i] - n21.a[j][i]);
               c2 = .5 * (n12.a[i][j] + n11.a[i][j]);
               d2 = .5 * (n12.a[i][j] - n11.a[i][j]);

               N11.a[i][j] = a1 + c1;
               N12.a[i][j] = a1 - c1;
               N21.a[i][j] = b1 + d1;
               N22.a[i][j] = b1 - d1;
/*
               D11.a[i][j] = b1 - d1;
               D12.a[i][j] = b1 + d1;
               D21.a[i][j] = a1 - c1;
               D22.a[i][j] = a1 + c1;
*/
               N11.a[j][i] = a2 + c2;
               N12.a[j][i] = a2 - c2;
               N21.a[j][i] = b2 + d2;
               N22.a[j][i] = b2 - d2;
/*
               D11.a[j][i] = b2 - d2;
               D12.a[j][i] = b2 + d2;
               D21.a[j][i] = a2 - c2;
               D22.a[j][i] = a2 + c2;
*/
#endif
               a1 = .5 * (d11.a[i][j] + d21.a[i][j]);
               b1 = .5 * (d11.a[i][j] - d21.a[i][j]);
               c1 = .5 * (d12.a[j][i] + d11.a[j][i]);
               d1 = .5 * (d12.a[j][i] - d11.a[j][i]);

               a2 = .5 * (d11.a[j][i] + d21.a[j][i]);
               b2 = .5 * (d11.a[j][i] - d21.a[j][i]);
               c2 = .5 * (d12.a[i][j] + d11.a[i][j]);
               d2 = .5 * (d12.a[i][j] - d11.a[i][j]);

               D11.a[i][j] = a1 + c1;
               D12.a[i][j] = a1 - c1;
               D21.a[i][j] = b1 + d1;
               D22.a[i][j] = b1 - d1;
/*
               N11.a[i][j] = b1 - d1;
               N12.a[i][j] = b1 + d1;
               N21.a[i][j] = a1 - c1;
               N22.a[i][j] = a1 + c1;
*/
               D11.a[j][i] = a2 + c2;
               D12.a[j][i] = a2 - c2;
               D21.a[j][i] = b2 + d2;
               D22.a[j][i] = b2 - d2;
/*
               N11.a[j][i] = b2 - d2;
               N12.a[j][i] = b2 + d2;
               N21.a[j][i] = a2 - c2;
               N22.a[j][i] = a2 + c2;
*/
          }
     }
#endif

     for (i = 0; i < n_umus; ++i) {
          for (j = 0; j < n_quad; ++j) {
               a = .5 * (n11.d1[i][j] + n21.d1[i][j]);
               b = .5 * (n11.d1[i][j] - n21.d1[i][j]);
               c = .5 * (n12.d1[i][j] + n22.d1[i][j]);
               d = .5 * (n12.d1[i][j] - n22.d1[i][j]);

               N11.d1[i][j] = a + c;
               N12.d1[i][j] = a - c;
               N21.d1[i][j] = b + d;
               N22.d1[i][j] = b - d;
          }
     }

     for (i = 0; i < n_umus; ++i) {
          a = .5 * (n11.d2[i] + n21.d2[i]);
          b = .5 * (n11.d2[i] - n21.d2[i]);
          c = .5 * (n12.d2[i] + n22.d2[i]);
          d = .5 * (n12.d2[i] - n22.d2[i]);

          N11.d2[i] = a + c;
          N12.d2[i] = a - c;
          N21.d2[i] = b + d;
          N22.d2[i] = b - d;
     }
}



/*******************************************************************************
 *
 ******************************************************************************/
static void reduced_update(double **a, double **b, double **c, double alpha, int m, int n) {

     int i;
     int j;

     if (alpha == 1.) {
          for (i = 0; i < m; ++i) {
               for (j = 0; j < n; ++j) {
                    c[i][j] = a[i][j] + b[i][j];
               }
          }
     }
     else {
          for (i = 0; i < m; ++i) {
               for (j = 0; j < n; ++j) {
                    c[i][j] = a[i][j] + b[i][j] * alpha;
               }
          }
     }
}



static void reduced_update_d(double **au1, double *au2, double **bu1, double *bu2,
                             double **cu1, double *cu2, double alpha, int n_quad, int n_umus) {

     reduced_update( au1,  bu1,  cu1, alpha, n_umus, n_quad);
     reduced_update(&au2, &bu2, &cu2, alpha, 1,      n_umus);
}


/*
static void reduced_update_d2(matd1d2 a, matd1d2 b,
                              matd1d2 c, double alpha, int n_quad, int n_umus) {

     reduced_update(a.a, b.a, c.a, alpha, n_quad, n_quad);

     if (n_umus) {
          reduced_update( a.d1,  b.d1,  c.d1, alpha, n_umus, n_quad);
          reduced_update(&a.d2, &b.d2, &c.d2, alpha, 1,      n_umus);
     }
}
*/


/*******************************************************************************
 *
 ******************************************************************************/
static void calc_n_d_gen_3(int n_quad, int n_umus, int n_stokes, int n_derivs, double *pade_c, int pade_r,
                           matd1d2 A12, matd1d2 A21, matd1d2 *B12, matd1d2 *B21,
                           matd1d2 N11, matd1d2 N12, matd1d2 N21, matd1d2 N22,
                           matd1d2 D11, matd1d2 D12, matd1d2 D21, matd1d2 D22,
                           matd1d2 *O11, matd1d2 *O12, matd1d2 *O21, matd1d2 *O22,
                           matd1d2 *E11, matd1d2 *E12, matd1d2 *E21, matd1d2 *E22,
                           uchar *derivs, work_data work) {

     int i;
     int j;
     int k;

     int n_quad_v;
     int n_umus_v;

     double **w1u1;

     matd1d2 w1;
     matd1d2 w2;
     matd1d2 w3;
     matd1d2 w4;

     matd1d2 *wd1;
     matd1d2 *wd2;

     matd1d2 w11[2];
     matd1d2 w22[2];

     matd1d2 *wd11[2];
     matd1d2 *wd22[2];


     n_quad_v = n_quad * n_stokes;
     n_umus_v = n_umus * n_stokes;


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (pade_r >= 2) {
          w1 = get_work_d2_d(&work, n_quad_v, n_quad_v, n_umus_v);
          w2 = get_work_d2_d(&work, n_quad_v, n_quad_v, n_umus_v);
     }

     if (pade_r >  2) {
          w11[0] = get_work_d2_d(&work, n_quad_v, n_quad_v, n_umus_v);
          w11[1] = get_work_d2_d(&work, n_quad_v, n_quad_v, n_umus_v);
          w22[0] = get_work_d2_d(&work, n_quad_v, n_quad_v, n_umus_v);
          w22[1] = get_work_d2_d(&work, n_quad_v, n_quad_v, n_umus_v);
     }

     if (n_umus > 0) {
          if (pade_r >= 2)
               w1u1 = get_work_d2(&work, n_umus_v, n_quad_v);
     }

     if (flags_or(derivs, n_derivs)) {
          if (pade_r >  2) {
               w3 = get_work_d2_d(&work, n_quad_v, n_quad_v, n_umus_v);
               w4 = get_work_d2_d(&work, n_quad_v, n_quad_v, n_umus_v);

               wd1 = get_work_d2_d2(&work, n_derivs, n_quad_v, n_quad_v, n_umus_v);
               wd2 = get_work_d2_d2(&work, n_derivs, n_quad_v, n_quad_v, n_umus_v);

               wd11[0] = get_work_d2_d2(&work, n_derivs, n_quad_v, n_quad_v, n_umus_v);
               wd11[1] = get_work_d2_d2(&work, n_derivs, n_quad_v, n_quad_v, n_umus_v);
               wd22[0] = get_work_d2_d2(&work, n_derivs, n_quad_v, n_quad_v, n_umus_v);
               wd22[1] = get_work_d2_d2(&work, n_derivs, n_quad_v, n_quad_v, n_umus_v);
          }
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (pade_r == 1 ||
         pade_r == 2) {
          dmat_init_d2(N11, pade_c[0], n_quad_v, n_umus_v);
          dmat_init_d2(N22, pade_c[0], n_quad_v, n_umus_v);

          dmat_scale_d2(pade_c[1], A12, N12, n_quad_v, n_umus_v);
          dmat_scale_d2(pade_c[1], A21, N21, n_quad_v, n_umus_v);

          for (i = 0; i < n_derivs; ++i) {
               if (! derivs[i])
                    continue;

               dmat_init_d2(O11[i], 0., n_quad_v, n_umus_v);
               dmat_init_d2(O22[i], 0., n_quad_v, n_umus_v);

               dmat_scale_d2(pade_c[1], B12[i], O12[i], n_quad_v, n_umus_v);
               dmat_scale_d2(pade_c[1], B21[i], O21[i], n_quad_v, n_umus_v);
          }
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (pade_r == 2) {
          dmat_mul_d2(A12, A21, n_quad_v, n_umus_v, w1, w1u1);
          dmat_scale_add_d2(N11, w1, N11, pade_c[2], n_quad_v, n_umus_v);

          dmat_mul_d2(A21, A12, n_quad_v, n_umus_v, w1, w1u1);
          dmat_scale_add_d2(N22, w1, N22, pade_c[2], n_quad_v, n_umus_v);

          for (i = 0; i < n_derivs; ++i) {
               if (! derivs[i])
                    continue;

               dmat_mul_d2(B12[i], A21, n_quad_v, n_umus_v, w1, w1u1);
               dmat_mul_d2(A12, B21[i], n_quad_v, n_umus_v, w2, w1u1);
               dmat_add_d2(w1, w2, w1, n_quad_v, n_umus_v);

               dmat_scale_add_d2(O11[i], w1, O11[i], pade_c[2], n_quad_v, n_umus_v);

               dmat_mul_d2(B21[i], A12, n_quad_v, n_umus_v, w1, w1u1);
               dmat_mul_d2(A21, B12[i], n_quad_v, n_umus_v, w2, w1u1);
               dmat_add_d2(w1, w2, w1, n_quad_v, n_umus_v);

               dmat_scale_add_d2(O22[i], w1, O22[i], pade_c[2], n_quad_v, n_umus_v);
          }
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (pade_r == 1 || pade_r == 2)
          return;


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (pade_r == 12 || pade_r == 13) {
          dmat_init_d2(N11, pade_c[0], n_quad_v, n_umus_v);
          dmat_init_d2(N22, pade_c[0], n_quad_v, n_umus_v);
          dmat_init_d2(D11, pade_c[1], n_quad_v, n_umus_v);
          dmat_init_d2(D22, pade_c[1], n_quad_v, n_umus_v);

          dmat_init_d2(N12, pade_c[6], n_quad_v, n_umus_v);
          dmat_init_d2(N21, pade_c[6], n_quad_v, n_umus_v);
          dmat_init_d2(D12, pade_c[7], n_quad_v, n_umus_v);
          dmat_init_d2(D21, pade_c[7], n_quad_v, n_umus_v);

          add_u_v_term_3_d2(n_quad_v, n_umus_v, A12, A21, A21, A12, w1, w2, N11, N22, D11, D22, N12, N21, D12, D21, pade_c[2], 2 < pade_r ? pade_c[3] : 0., pade_c[8], 8 < pade_r ? pade_c[9] : 0., w1u1);

          dmat_copy_d2(w11[0], w1, n_quad_v, n_umus_v);
          dmat_copy_d2(w22[0], w2, n_quad_v, n_umus_v);


          for (i = 0; i < n_derivs; ++i) {
               if (! derivs[i])
                    continue;

               dmat_zero_d2(O11[i], n_quad_v, n_umus_v);
               dmat_zero_d2(O22[i], n_quad_v, n_umus_v);
               dmat_zero_d2(E11[i], n_quad_v, n_umus_v);
               dmat_zero_d2(E22[i], n_quad_v, n_umus_v);

               dmat_zero_d2(O12[i], n_quad_v, n_umus_v);
               dmat_zero_d2(O21[i], n_quad_v, n_umus_v);
               dmat_zero_d2(E12[i], n_quad_v, n_umus_v);
               dmat_zero_d2(E21[i], n_quad_v, n_umus_v);

               add_u_v_term_l_3_d2(n_quad_v, n_umus_v, A12, A21, A21, A12, B12[i], B21[i], B21[i], B12[i], wd1[i], wd2[i], O11[i], O22[i], E11[i], E22[i], O12[i], O21[i], E12[i], E21[i], w3, w4, w1u1, pade_c[2], 2 < pade_r ? pade_c[3] : 0., pade_c[8], 8 < pade_r ? pade_c[9] : 0.);

               dmat_copy_d2(wd11[0][i], wd1[i], n_quad_v, n_umus_v);
               dmat_copy_d2(wd22[0][i], wd2[i], n_quad_v, n_umus_v);
          }


          i = 4;
          j = 1;
          add_u_v_term_3_d2(n_quad_v, n_umus_v, w11[1 - j], w22[1 - j], w1, w2, w11[j], w22[j], N11,  N22,  D11,  D22,  N12, N21, D12, D21, pade_c[i], i < pade_r ? pade_c[i+1] : 0., pade_c[i+6], i+6 < pade_r ? pade_c[i+7] : 0., w1u1);
          for (k = 0; k < n_derivs; ++k) {
               if (! derivs[k])
                    continue;

               add_u_v_term_l_3_d2(n_quad_v, n_umus_v, w11[1 - j], w22[1 - j], w1, w2, wd11[1 - j][k], wd22[1 - j][k], wd1[k], wd2[k], wd11[j][k], wd22[j][k], O11[k], O22[k], E11[k], E22[k], O12[k], O21[k], E12[k], E21[k], w3, w4, w1u1, pade_c[i], i < pade_r ? pade_c[i+1] : 0., pade_c[i+6], i+6 < pade_r ? pade_c[i+7] : 0.);
          }
          i = i + 2;
          j = 1 - j;

          add_u_v_term_3_d2(n_quad_v, n_umus_v, w11[1 - j], w22[1 - j], w1, w2, w11[j], w22[j], NONE, NONE, NONE, NONE, N12, N21, D12, D21, 0.,       0.,                             pade_c[i+6], i+6 < pade_r ? pade_c[i+7] : 0., w1u1);
          for (k = 0; k < n_derivs; ++k) {
               if (! derivs[k])
                    continue;

               add_u_v_term_l_3_d2(n_quad_v, n_umus_v, w11[1 - j], w22[1 - j], w1, w2, wd11[1 - j][k], wd22[1 - j][k], wd1[k], wd2[k], wd11[j][k], wd22[j][k], NONE,   NONE,   NONE,   NONE,   O12[k], O21[k], E12[k], E21[k], w3, w4, w1u1, 0.,        0.,                            pade_c[i+6], i+6 < pade_r ? pade_c[i+7] : 0.);
          }

          dmat_mul_d2(w11[j], N12, n_quad_v, n_umus_v, w1, w1u1);
          dmat_add_d2(w1, N11, N11, n_quad_v, n_umus_v);

          dmat_mul_d2(w22[j], N21, n_quad_v, n_umus_v, w1, w1u1);
          dmat_add_d2(w1, N22, N22, n_quad_v, n_umus_v);


          for (i = 0; i < n_derivs; ++i) {
               if (! derivs[i])
                    continue;

               dmat_mul_d2(wd11[j][i], N12, n_quad_v, n_umus_v, w1, w1u1);
               dmat_add_d2(w1, O11[i], O11[i], n_quad_v, n_umus_v);
               dmat_mul_d2(w11[j], O12[i], n_quad_v, n_umus_v, w1, w1u1);
               dmat_add_d2(w1, O11[i], O11[i], n_quad_v, n_umus_v);

               dmat_mul_d2(wd22[j][i], N21, n_quad_v, n_umus_v, w1, w1u1);
               dmat_add_d2(w1, O22[i], O22[i], n_quad_v, n_umus_v);
               dmat_mul_d2(w22[j], O21[i], n_quad_v, n_umus_v, w1, w1u1);
               dmat_add_d2(w1, O22[i], O22[i], n_quad_v, n_umus_v);
          }


          dmat_mul_d2(w11[j], D12, n_quad_v, n_umus_v, w1, w1u1);
          dmat_add_d2(w1, D11, D11, n_quad_v, n_umus_v);

          dmat_mul_d2(w22[j], D21, n_quad_v, n_umus_v, w1, w1u1);
          dmat_add_d2(w1, D22, D22, n_quad_v, n_umus_v);


          for (i = 0; i < n_derivs; ++i) {
               if (! derivs[i])
                    continue;

               dmat_mul_d2(wd11[j][i], D12, n_quad_v, n_umus_v, w1, w1u1);
               dmat_add_d2(w1, E11[i], E11[i], n_quad_v, n_umus_v);
               dmat_mul_d2(w11[j], E12[i], n_quad_v, n_umus_v, w1, w1u1);
               dmat_add_d2(w1, E11[i], E11[i], n_quad_v, n_umus_v);

               dmat_mul_d2(wd22[j][i], D21, n_quad_v, n_umus_v, w1, w1u1);
               dmat_add_d2(w1, E22[i], E22[i], n_quad_v, n_umus_v);
               dmat_mul_d2(w22[j], E21[i], n_quad_v, n_umus_v, w1, w1u1);
               dmat_add_d2(w1, E22[i], E22[i], n_quad_v, n_umus_v);
          }
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     else {
          dmat_init_d2(N11, pade_c[0], n_quad_v, n_umus_v);
          dmat_init_d2(N22, pade_c[0], n_quad_v, n_umus_v);

          dmat_init_d2(D11, pade_c[1], n_quad_v, n_umus_v);
          dmat_init_d2(D22, pade_c[1], n_quad_v, n_umus_v);

          add_u_v_term_3_d2(n_quad_v, n_umus_v, A12, A21, A21, A12, w1, w2, N11, N22, D11, D22, NONE, NONE, NONE, NONE, pade_c[2], 2 < pade_r ? pade_c[3] : 0., 0., 0., w1u1);

          dmat_copy_d2(w11[0], w1, n_quad_v, n_umus_v);
          dmat_copy_d2(w22[0], w2, n_quad_v, n_umus_v);

          for (i = 0; i < n_derivs; ++i) {
               if (! derivs[i])
                    continue;

               dmat_zero_d2(O11[i], n_quad_v, n_umus_v);
               dmat_zero_d2(O22[i], n_quad_v, n_umus_v);

               dmat_zero_d2(E11[i], n_quad_v, n_umus_v);
               dmat_zero_d2(E22[i], n_quad_v, n_umus_v);

               add_u_v_term_l_3_d2(n_quad_v, n_umus_v, A12, A21, A21, A12, B12[i], B21[i], B21[i], B12[i], wd1[i], wd2[i], O11[i], O22[i], E11[i], E22[i], NONE, NONE, NONE, NONE, w3, w4, w1u1, pade_c[2], 2 < pade_r ? pade_c[3] : 0., 0., 0.);

               dmat_copy_d2(wd11[0][i], wd1[i], n_quad_v, n_umus_v);
               dmat_copy_d2(wd22[0][i], wd2[i], n_quad_v, n_umus_v);
          }

          j = 1;
          for (i = 4; i <= pade_r; i += 2) {
               add_u_v_term_3_d2(n_quad_v, n_umus_v, w11[1 - j], w22[1 - j], w1, w2, w11[j], w22[j], N11, N22, D11, D22, NONE, NONE, NONE, NONE, pade_c[i], i < pade_r ? pade_c[i+1] : 0., 0., 0., w1u1);

               for (k = 0; k < n_derivs; ++k) {
                    if (! derivs[k])
                         continue;

                    add_u_v_term_l_3_d2(n_quad_v, n_umus_v, w11[1 - j], w22[1 - j], w1, w2, wd11[1 - j][k], wd22[1 - j][k], wd1[k], wd2[k], wd11[j][k], wd22[j][k], O11[k], O22[k], E11[k], E22[k], NONE, NONE, NONE, NONE, w3, w4, w1u1, pade_c[i], i < pade_r ? pade_c[i+1] : 0., 0., 0.);
               }

               j = 1 - j;
          }
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     dmat_mul_d2(A12, D22, n_quad_v, n_umus_v, N12, w1u1);
     dmat_mul_d2(A21, D11, n_quad_v, n_umus_v, N21, w1u1);

     for (i = 0; i < n_derivs; ++i) {
          if (! derivs[i])
               continue;

          dmat_mul_d2(B12[i], D22, n_quad_v, n_umus_v, w1, w1u1);
          dmat_mul_d2(A12, E22[i], n_quad_v, n_umus_v, w2, w1u1);
          dmat_add_d2(w1, w2, O12[i], n_quad_v, n_umus_v);

          dmat_mul_d2(B21[i], D11, n_quad_v, n_umus_v, w1, w1u1);
          dmat_mul_d2(A21, E11[i], n_quad_v, n_umus_v, w2, w1u1);
          dmat_add_d2(w1, w2, O21[i], n_quad_v, n_umus_v);
     }
}



/*******************************************************************************
 *
 ******************************************************************************/
static void calc_n_d_sym_3(int n_quad, int n_umus, int n_stokes, int n_derivs, double *pade_c, int pade_r,
                           matd1d2 A12, matd1d2 A21, matd1d2 *B12, matd1d2 *B21,
                           matd1d2 N11, matd1d2 N12, matd1d2 N21, matd1d2 N22,
                           matd1d2 D11, matd1d2 D12, matd1d2 D21, matd1d2 D22,
                           matd1d2 *O11, matd1d2 *O12, matd1d2 *O21, matd1d2 *O22,
                           matd1d2 *E11, matd1d2 *E12, matd1d2 *E21, matd1d2 *E22,
                           uchar *derivs, work_data work) {

     int i;
     int j;
     int k;

     int n_quad_v;
     int n_umus_v;

     double **w1u1;

     matd1d2 w1;
     matd1d2 w2;
     matd1d2 w3;
     matd1d2 w4;

     matd1d2 *wd1;
     matd1d2 *wd2;

     matd1d2 wx1[2];
     matd1d2 wx2[2];

     matd1d2 *wy1[2];
     matd1d2 *wy2[2];

     n_quad_v = n_quad * n_stokes;
     n_umus_v = n_umus * n_stokes;


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (pade_r >= 2) {
          w1 = get_work_d2_d(&work, n_quad_v, n_quad_v, n_umus_v);
          w2 = get_work_d2_d(&work, n_quad_v, n_quad_v, n_umus_v);
     }

     if (pade_r >  2) {
          wx1[0] = get_work_d2_d(&work, n_quad_v, n_quad_v, n_umus_v);
          wx1[1] = get_work_d2_d(&work, n_quad_v, n_quad_v, n_umus_v);
          wx2[0] = get_work_d2_d(&work, n_quad_v, n_quad_v, n_umus_v);
          wx2[1] = get_work_d2_d(&work, n_quad_v, n_quad_v, n_umus_v);
     }

     if (n_umus > 0) {
          if (pade_r >= 2 && ! flags_or(derivs, n_derivs))
               w3 = get_work_d2_d(&work, n_quad_v, n_quad_v, n_umus_v);

          if (pade_r >= 2)
               w1u1 = get_work_d2(&work, n_umus_v, n_quad_v);
     }

     if (flags_or(derivs, n_derivs)) {
          if (pade_r >= 2) {
               w3 = get_work_d2_d(&work, n_quad_v, n_quad_v, n_umus_v);
               w4 = get_work_d2_d(&work, n_quad_v, n_quad_v, n_umus_v);

               wd1 = get_work_d2_d2(&work, n_derivs, n_quad_v, n_quad_v, n_umus_v);
               wd2 = get_work_d2_d2(&work, n_derivs, n_quad_v, n_quad_v, n_umus_v);
          }

          if (pade_r > 2) {
               wy1[0] = get_work_d2_d2(&work, n_derivs, n_quad_v, n_quad_v, n_umus_v);
               wy1[1] = get_work_d2_d2(&work, n_derivs, n_quad_v, n_quad_v, n_umus_v);
               wy2[0] = get_work_d2_d2(&work, n_derivs, n_quad_v, n_quad_v, n_umus_v);
               wy2[1] = get_work_d2_d2(&work, n_derivs, n_quad_v, n_quad_v, n_umus_v);
          }
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (pade_r == 1 ||
         pade_r == 2) {

          dmat_init(N11.a, pade_c[0], n_quad_v, n_quad_v);
          dmat_init(N22.a, pade_c[0], n_quad_v, n_quad_v);

          dmat_scale(pade_c[1], A12.a, N12.a, n_quad_v, n_quad_v);
          dmat_scale(pade_c[1], A21.a, N21.a, n_quad_v, n_quad_v);

          if (n_umus > 0) {
               dmat_init_d(N11.d1, N11.d2, pade_c[0], n_quad_v, n_umus_v);
               dmat_init_d(N22.d1, N22.d2, pade_c[0], n_quad_v, n_umus_v);

               dmat_scale_d(pade_c[1], A12.d1, A12.d2, N12.d1, N12.d2, n_quad_v, n_umus_v);
               dmat_scale_d(pade_c[1], A21.d1, A21.d2, N21.d1, N21.d2, n_quad_v, n_umus_v);
          }

          for (i = 0; i < n_derivs; ++i) {
               if (! derivs[i])
                    continue;

               dmat_init(O11[i].a, 0., n_quad_v, n_quad_v);
               dmat_init(O22[i].a, 0., n_quad_v, n_quad_v);

               dmat_scale(pade_c[1], B12[i].a, O12[i].a, n_quad_v, n_quad_v);
               dmat_scale(pade_c[1], B21[i].a, O21[i].a, n_quad_v, n_quad_v);

               if (n_umus > 0) {
                    dmat_init_d(O11[i].d1, O11[i].d2, 0., n_quad_v, n_umus_v);
                    dmat_init_d(O22[i].d1, O22[i].d2, 0., n_quad_v, n_umus_v);

                    dmat_scale_d(pade_c[1], B12[i].d1, B12[i].d2, O12[i].d1, O12[i].d2, n_quad_v, n_umus_v);
                    dmat_scale_d(pade_c[1], B21[i].d1, B21[i].d2, O21[i].d1, O21[i].d2, n_quad_v, n_umus_v);
               }
          }
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (pade_r == 2) {
          dmat_mul(A12.a, A21.a, n_quad_v, n_quad_v, n_quad_v, w1.a);
          reduced_update(N11.a, w1.a, N11.a, pade_c[2], n_quad_v, n_quad_v);

          dmat_trans_stokes(n_quad, n_stokes, N11.a, N22.a);

          if (n_umus > 0) {
               dmat_mul_d(A12.d1, A12.d2, A21.a, A21.d1, A21.d2, n_quad_v, n_umus_v, w1.d1, w1.d2, w1u1);
               reduced_update_d(N11.d1, N11.d2, w1.d1, w1.d2, N11.d1, N11.d2, pade_c[2], n_quad_v, n_umus_v);

               dmat_mul_d(A21.d1, A21.d2, A12.a, A12.d1, A12.d2, n_quad_v, n_umus_v, w1.d1, w1.d2, w1u1);
               reduced_update_d(N22.d1, N22.d2, w1.d1, w1.d2, N22.d1, N22.d2, pade_c[2], n_quad_v, n_umus_v);
          }

          for (i = 0; i < n_derivs; ++i) {
               if (! derivs[i])
                    continue;

               dmat_mul(B12[i].a, A21.a, n_quad_v, n_quad_v, n_quad_v, w1.a);
               dmat_mul(A12.a, B21[i].a, n_quad_v, n_quad_v, n_quad_v, w2.a);
               dmat_add(w1.a, w2.a, w1.a, n_quad_v, n_quad_v);

               reduced_update(O11[i].a, w1.a, O11[i].a, pade_c[2], n_quad_v, n_quad_v);

               dmat_trans_stokes(n_quad, n_stokes, O11[i].a, O22[i].a);

               if (n_umus > 0) {
                    dmat_mul_d(B12[i].d1, B12[i].d2, A21.a, A21.d1, A21.d2,    n_quad_v, n_umus_v, w1.d1, w1.d2, w1u1);
                    dmat_mul_d(A12.d1, A12.d2, B21[i].a, B21[i].d1, B21[i].d2, n_quad_v, n_umus_v, w2.d1, w2.d2, w1u1);
                    dmat_add_d(w1.d1, w1.d2, w2.d1, w2.d2, w1.d1, w1.d2, n_quad_v, n_umus_v);

                    reduced_update_d(O11[i].d1, O11[i].d2, w1.d1, w1.d2, O11[i].d1, O11[i].d2, pade_c[2], n_quad_v, n_umus_v);

                    dmat_mul_d(B21[i].d1, B21[i].d2, A12.a, A12.d1, A12.d2,    n_quad_v, n_umus_v, w1.d1, w1.d2, w1u1);
                    dmat_mul_d(A21.d1, A21.d2, B12[i].a, B12[i].d1, B12[i].d2, n_quad_v, n_umus_v, w2.d1, w2.d2, w1u1);
                    dmat_add_d(w1.d1, w1.d2, w2.d1, w2.d2, w1.d1, w1.d2, n_quad_v, n_umus_v);

                    reduced_update_d(O22[i].d1, O22[i].d2, w1.d1, w1.d2, O22[i].d1, O22[i].d2, pade_c[2], n_quad_v, n_umus_v);
               }
          }
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (pade_r == 1 || pade_r == 2)
          return;


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (pade_r == 12 || pade_r == 13) {
          dmat_init(N11.a, pade_c[0], n_quad_v, n_quad_v);
          dmat_init(D11.a, pade_c[1], n_quad_v, n_quad_v);

          dmat_init(N12.a, pade_c[6], n_quad_v, n_quad_v);
          dmat_init(D12.a, pade_c[7], n_quad_v, n_quad_v);

          dmat_mul(A12.a, A21.a, n_quad_v, n_quad_v, n_quad_v, w1.a);

          reduced_update(N11.a, w1.a, N11.a, pade_c[2], n_quad_v, n_quad_v);
          if (2 < pade_r)
               reduced_update(D11.a, w1.a, D11.a, pade_c[3], n_quad_v, n_quad_v);
          reduced_update(N12.a, w1.a, N12.a, pade_c[8], n_quad_v, n_quad_v);
          if (8 < pade_r)
               reduced_update(D12.a, w1.a, D12.a, pade_c[9], n_quad_v, n_quad_v);

          dmat_copy(wx1[0].a, w1.a, n_quad_v, n_quad_v);

          if (n_umus > 0) {
               dmat_init_d(N11.d1, N11.d2, pade_c[0], n_quad_v, n_umus_v);
               dmat_init_d(N22.d1, N22.d2, pade_c[0], n_quad_v, n_umus_v);

               dmat_init_d(D11.d1, D11.d2, pade_c[1], n_quad_v, n_umus_v);
               dmat_init_d(D22.d1, D22.d2, pade_c[1], n_quad_v, n_umus_v);

               dmat_init_d(N12.d1, N12.d2, pade_c[6], n_quad_v, n_umus_v);
               dmat_init_d(N21.d1, N21.d2, pade_c[6], n_quad_v, n_umus_v);

               dmat_init_d(D12.d1, D12.d2, pade_c[7], n_quad_v, n_umus_v);
               dmat_init_d(D21.d1, D21.d2, pade_c[7], n_quad_v, n_umus_v);

               dmat_mul_d(A12.d1, A12.d2, A21.a, A21.d1, A21.d2, n_quad_v, n_umus_v, w1.d1, w1.d2, w1u1);
               dmat_mul_d(A21.d1, A21.d2, A12.a, A12.d1, A12.d2, n_quad_v, n_umus_v, w2.d1, w2.d2, w1u1);

               reduced_update_d(N11.d1, N11.d2, w1.d1, w1.d2, N11.d1, N11.d2, pade_c[2], n_quad_v, n_umus_v);
               reduced_update_d(N22.d1, N22.d2, w2.d1, w2.d2, N22.d1, N22.d2, pade_c[2], n_quad_v, n_umus_v);
               if (2 < pade_r) {
                    reduced_update_d(D11.d1, D11.d2, w1.d1, w1.d2, D11.d1, D11.d2, pade_c[3], n_quad_v, n_umus_v);
                    reduced_update_d(D22.d1, D22.d2, w2.d1, w2.d2, D22.d1, D22.d2, pade_c[3], n_quad_v, n_umus_v);
               }

               reduced_update_d(N12.d1, N12.d2, w1.d1, w1.d2, N12.d1, N12.d2, pade_c[8], n_quad_v, n_umus_v);
               reduced_update_d(N21.d1, N21.d2, w2.d1, w2.d2, N21.d1, N21.d2, pade_c[8], n_quad_v, n_umus_v);
               if (8 < pade_r) {
                    reduced_update_d(D12.d1, D12.d2, w1.d1, w1.d2, D12.d1, D12.d2, pade_c[9], n_quad_v, n_umus_v);
                    reduced_update_d(D21.d1, D21.d2, w2.d1, w2.d2, D21.d1, D21.d2, pade_c[9], n_quad_v, n_umus_v);
               }

               dmat_copy_d(wx1[0].d1, wx1[0].d2, w1.d1, w1.d2, n_quad_v, n_umus_v);
               dmat_copy_d(wx2[0].d1, wx2[0].d2, w2.d1, w2.d2, n_quad_v, n_umus_v);
          }

          for (i = 0; i < n_derivs; ++i) {
               if (! derivs[i])
                    continue;

               dmat_zero(O11[i].a, n_quad_v, n_quad_v);
               dmat_zero(E11[i].a, n_quad_v, n_quad_v);

               dmat_zero(O12[i].a, n_quad_v, n_quad_v);
               dmat_zero(E12[i].a, n_quad_v, n_quad_v);

               dmat_mul(B12[i].a, A21.a, n_quad_v, n_quad_v, n_quad_v, w2.a);
               dmat_mul(A12.a, B21[i].a, n_quad_v, n_quad_v, n_quad_v, w3.a);
               dmat_add(w2.a, w3.a, wd1[i].a, n_quad_v, n_quad_v);

               reduced_update(O11[i].a, wd1[i].a, O11[i].a, pade_c[2], n_quad_v, n_quad_v);
               if (2 < pade_r)
                    reduced_update(E11[i].a, wd1[i].a, E11[i].a, pade_c[3], n_quad_v, n_quad_v);
               reduced_update(O12[i].a, wd1[i].a, O12[i].a, pade_c[8], n_quad_v, n_quad_v);
               if (8 < pade_r)
                    reduced_update(E12[i].a, wd1[i].a, E12[i].a, pade_c[9], n_quad_v, n_quad_v);

               dmat_copy(wy1[0][i].a, wd1[i].a, n_quad_v, n_quad_v);

               if (n_umus > 0) {
                    dmat_zero_d(O11[i].d1, O11[i].d2, n_quad_v, n_umus_v);
                    dmat_zero_d(O22[i].d1, O22[i].d2, n_quad_v, n_umus_v);

                    dmat_zero_d(E11[i].d1, E11[i].d2, n_quad_v, n_umus_v);
                    dmat_zero_d(E22[i].d1, E22[i].d2, n_quad_v, n_umus_v);

                    dmat_zero_d(O12[i].d1, O12[i].d2, n_quad_v, n_umus_v);
                    dmat_zero_d(O21[i].d1, O21[i].d2, n_quad_v, n_umus_v);

                    dmat_zero_d(E12[i].d1, E12[i].d2, n_quad_v, n_umus_v);
                    dmat_zero_d(E21[i].d1, E21[i].d2, n_quad_v, n_umus_v);


                    dmat_mul_d(B12[i].d1, B12[i].d2, A21.a, A21.d1, A21.d2,    n_quad_v, n_umus_v, w3.d1, w3.d2, w1u1);
                    dmat_mul_d(A12.d1, A12.d2, B21[i].a, B21[i].d1, B21[i].d2, n_quad_v, n_umus_v, w4.d1, w4.d2, w1u1);
                    dmat_add_d(w3.d1, w3.d2, w4.d1, w4.d2, wd1[i].d1, wd1[i].d2, n_quad_v, n_umus_v);

                    dmat_mul_d(B21[i].d1, B21[i].d2, A12.a, A12.d1, A12.d2,    n_quad_v, n_umus_v, w3.d1, w3.d2, w1u1);
                    dmat_mul_d(A21.d1, A21.d2, B12[i].a, B12[i].d1, B12[i].d2, n_quad_v, n_umus_v, w4.d1, w4.d2, w1u1);
                    dmat_add_d(w3.d1, w3.d2, w4.d1, w4.d2, wd2[i].d1, wd2[i].d2, n_quad_v, n_umus_v);


                    reduced_update_d(O11[i].d1, O11[i].d2, wd1[i].d1, wd1[i].d2, O11[i].d1, O11[i].d2, pade_c[2], n_quad_v, n_umus_v);
                    reduced_update_d(O22[i].d1, O22[i].d2, wd2[i].d1, wd2[i].d2, O22[i].d1, O22[i].d2, pade_c[2], n_quad_v, n_umus_v);
                    if (2 < pade_r) {
                         reduced_update_d(E11[i].d1, E11[i].d2, wd1[i].d1, wd1[i].d2, E11[i].d1, E11[i].d2, pade_c[3], n_quad_v, n_umus_v);
                         reduced_update_d(E22[i].d1, E22[i].d2, wd2[i].d1, wd2[i].d2, E22[i].d1, E22[i].d2, pade_c[3], n_quad_v, n_umus_v);
                    }

                    dmat_copy_d(wy1[0][i].d1, wy1[0][i].d2, wd1[i].d1, wd1[i].d2, n_quad_v, n_umus_v);


                    reduced_update_d(O12[i].d1, O12[i].d2, wd1[i].d1, wd1[i].d2, O12[i].d1, O12[i].d2, pade_c[8], n_quad_v, n_umus_v);
                    reduced_update_d(O21[i].d1, O21[i].d2, wd2[i].d1, wd2[i].d2, O21[i].d1, O21[i].d2, pade_c[8], n_quad_v, n_umus_v);
                    if (8 < pade_r) {
                         reduced_update_d(E12[i].d1, E12[i].d2, wd1[i].d1, wd1[i].d2, E12[i].d1, E12[i].d2, pade_c[9], n_quad_v, n_umus_v);
                         reduced_update_d(E21[i].d1, E21[i].d2, wd2[i].d1, wd2[i].d2, E21[i].d1, E21[i].d2, pade_c[9], n_quad_v, n_umus_v);
                    }

                    dmat_copy_d(wy2[0][i].d1, wy2[0][i].d2, wd2[i].d1, wd2[i].d2, n_quad_v, n_umus_v);
               }
          }

          i = 4;
          j = 1;

          dmat_mul(wx1[1 - j].a, w1.a, n_quad_v, n_quad_v, n_quad_v, wx1[j].a);

          reduced_update(N11.a, wx1[j].a, N11.a, pade_c[i], n_quad_v, n_quad_v);
          if (i < pade_r)
               reduced_update(D11.a, wx1[j].a, D11.a, pade_c[i+1], n_quad_v, n_quad_v);
          reduced_update(N12.a, wx1[j].a, N12.a, pade_c[i+6], n_quad_v, n_quad_v);
          if (i+6 < pade_r)
               reduced_update(D12.a, wx1[j].a, D12.a, pade_c[i+7], n_quad_v, n_quad_v);

          if (n_umus > 0) {
               dmat_mul_d(wx1[1 - j].d1, wx1[1 - j].d2, w1.a, w1.d1, w1.d2, n_quad_v, n_umus_v, wx1[j].d1, wx1[j].d2, w1u1);
               dmat_trans_stokes(n_quad, n_stokes, w1.a, w2.a);
               dmat_mul_d(wx2[1 - j].d1, wx2[1 - j].d2, w2.a, w2.d1, w2.d2, n_quad_v, n_umus_v, wx2[j].d1, wx2[j].d2, w1u1);

               reduced_update_d(N11.d1, N11.d2, wx1[j].d1, wx1[j].d2, N11.d1, N11.d2, pade_c[i], n_quad_v, n_umus_v);
               reduced_update_d(N22.d1, N22.d2, wx2[j].d1, wx2[j].d2, N22.d1, N22.d2, pade_c[i], n_quad_v, n_umus_v);
               if (i < pade_r) {
                    reduced_update_d(D11.d1, D11.d2, wx1[j].d1, wx1[j].d2, D11.d1, D11.d2, pade_c[i+1], n_quad_v, n_umus_v);
                    reduced_update_d(D22.d1, D22.d2, wx2[j].d1, wx2[j].d2, D22.d1, D22.d2, pade_c[i+1], n_quad_v, n_umus_v);
               }

               reduced_update_d(N12.d1, N12.d2, wx1[j].d1, wx1[j].d2, N12.d1, N12.d2, pade_c[i+6], n_quad_v, n_umus_v);
               reduced_update_d(N21.d1, N21.d2, wx2[j].d1, wx2[j].d2, N21.d1, N21.d2, pade_c[i+6], n_quad_v, n_umus_v);
               if (i+6 < pade_r) {
                    reduced_update_d(D12.d1, D12.d2, wx1[j].d1, wx1[j].d2, D12.d1, D12.d2, pade_c[i+7], n_quad_v, n_umus_v);
                    reduced_update_d(D21.d1, D21.d2, wx2[j].d1, wx2[j].d2, D21.d1, D21.d2, pade_c[i+7], n_quad_v, n_umus_v);
               }
          }

          for (k = 0; k < n_derivs; ++k) {
               if (! derivs[k])
                    continue;

               dmat_mul(wy1[1 - j][k].a, w1.a, n_quad_v, n_quad_v, n_quad_v, w2.a);
               dmat_mul(wx1[1 - j].a, wd1[k].a, n_quad_v, n_quad_v, n_quad_v, w3.a);
               dmat_add(w2.a, w3.a, wy1[j][k].a, n_quad_v, n_quad_v);

               reduced_update(O11[k].a, wy1[j][k].a, O11[k].a, pade_c[i], n_quad_v, n_quad_v);
               if (i < pade_r)
                    reduced_update(E11[k].a, wy1[j][k].a, E11[k].a, pade_c[i+1], n_quad_v, n_quad_v);
               reduced_update(O12[k].a, wy1[j][k].a, O12[k].a, pade_c[i+6], n_quad_v, n_quad_v);
               if (i+6 < pade_r)
                    reduced_update(E12[k].a, wy1[j][k].a, E12[k].a, pade_c[i+7], n_quad_v, n_quad_v);

               if (n_umus > 0) {
                    dmat_mul_d(wy1[1 - j][k].d1, wy1[1 - j][k].d2, w1.a, w1.d1, w1.d2, n_quad_v, n_umus_v, w3.d1, w3.d2, w1u1);
                    dmat_mul_d(wx1[1 - j].d1, wx1[1 - j].d2, wd1[k].a, wd1[k].d1, wd1[k].d2, n_quad_v, n_umus_v, w4.d1, w4.d2, w1u1);
                    dmat_add_d(w3.d1, w3.d2, w4.d1, w4.d2, wy1[j][k].d1, wy1[j][k].d2, n_quad_v, n_umus_v);

                    dmat_trans_stokes(n_quad, n_stokes, w1.a, w2.a);
                    dmat_mul_d(wy2[1 - j][k].d1, wy2[1 - j][k].d2, w2.a, w2.d1, w2.d2, n_quad_v, n_umus_v, w3.d1, w3.d2, w1u1);
                    dmat_trans_stokes(n_quad, n_stokes, wd1[k].a, wd2[k].a);
                    dmat_mul_d(wx2[1 - j].d1, wx2[1 - j].d2, wd2[k].a, wd2[k].d1, wd2[k].d2, n_quad_v, n_umus_v, w4.d1, w4.d2, w1u1);
                    dmat_add_d(w3.d1, w3.d2, w4.d1, w4.d2, wy2[j][k].d1, wy2[j][k].d2, n_quad_v, n_umus_v);


                    reduced_update_d(O11[k].d1, O11[k].d2, wy1[j][k].d1, wy1[j][k].d2, O11[k].d1, O11[k].d2, pade_c[i], n_quad_v, n_umus_v);
                    reduced_update_d(O22[k].d1, O22[k].d2, wy2[j][k].d1, wy2[j][k].d2, O22[k].d1, O22[k].d2, pade_c[i], n_quad_v, n_umus_v);
                    if (i < pade_r) {
                         reduced_update_d(E11[k].d1, E11[k].d2, wy1[j][k].d1, wy1[j][k].d2, E11[k].d1, E11[k].d2, pade_c[i+1], n_quad_v, n_umus_v);
                         reduced_update_d(E22[k].d1, E22[k].d2, wy2[j][k].d1, wy2[j][k].d2, E22[k].d1, E22[k].d2, pade_c[i+1], n_quad_v, n_umus_v);
                    }

                    reduced_update_d(O12[k].d1, O12[k].d2, wy1[j][k].d1, wy1[j][k].d2, O12[k].d1, O12[k].d2, pade_c[i+6], n_quad_v, n_umus_v);
                    reduced_update_d(O21[k].d1, O21[k].d2, wy2[j][k].d1, wy2[j][k].d2, O21[k].d1, O21[k].d2, pade_c[i+6], n_quad_v, n_umus_v);
                    if (i+6 < pade_r) {
                         reduced_update_d(E12[k].d1, E12[k].d2, wy1[j][k].d1, wy1[j][k].d2, E12[k].d1, E12[k].d2, pade_c[i+7], n_quad_v, n_umus_v);
                         reduced_update_d(E21[k].d1, E21[k].d2, wy2[j][k].d1, wy2[j][k].d2, E21[k].d1, E21[k].d2, pade_c[i+7], n_quad_v, n_umus_v);
                    }
               }
          }

          i = i + 2;
          j = 1 - j;


          dmat_mul(wx1[1 - j].a, w1.a, n_quad_v, n_quad_v, n_quad_v, wx1[j].a);

          reduced_update(N12.a, wx1[j].a, N12.a, pade_c[i+6], n_quad_v, n_quad_v);
          if (i+6 < pade_r)
               reduced_update(D12.a, wx1[j].a, D12.a, pade_c[i+7], n_quad_v, n_quad_v);

          if (n_umus > 0) {
               dmat_mul_d(wx1[1 - j].d1, wx1[1 - j].d2, w1.a, w1.d1, w1.d2, n_quad_v, n_umus_v, wx1[j].d1, wx1[j].d2, w1u1);
               dmat_trans_stokes(n_quad, n_stokes, w1.a, w2.a);
               dmat_mul_d(wx2[1 - j].d1, wx2[1 - j].d2, w2.a, w2.d1, w2.d2, n_quad_v, n_umus_v, wx2[j].d1, wx2[j].d2, w1u1);

               reduced_update_d(N12.d1, N12.d2, wx1[j].d1, wx1[j].d2, N12.d1, N12.d2, pade_c[i+6], n_quad_v, n_umus_v);
               reduced_update_d(N21.d1, N21.d2, wx1[j].d1, wx1[j].d2, N21.d1, N21.d2, pade_c[i+6], n_quad_v, n_umus_v);

               if (i+6 < pade_r)  {
                    reduced_update_d(D12.d1, D12.d2, wx1[j].d1, wx1[j].d2, D12.d1, D12.d2, pade_c[i+7], n_quad_v, n_umus_v);
                    reduced_update_d(D21.d1, D21.d2, wx1[j].d1, wx1[j].d2, D21.d1, D21.d2, pade_c[i+7], n_quad_v, n_umus_v);
               }
          }

          for (k = 0; k < n_derivs; ++k) {
               if (! derivs[k])
                    continue;

               dmat_mul(wy1[1 - j][k].a, w1.a, n_quad_v, n_quad_v, n_quad_v, w2.a);
               dmat_mul(wx1[1 - j].a, wd1[k].a, n_quad_v, n_quad_v, n_quad_v, w3.a);
               dmat_add(w2.a, w3.a, wy1[j][k].a, n_quad_v, n_quad_v);

               reduced_update(O12[k].a, wy1[j][k].a, O12[k].a, pade_c[i+6], n_quad_v, n_quad_v);
               if (i+6 < pade_r)
                    reduced_update(E12[k].a, wy1[j][k].a, E12[k].a, pade_c[i+7], n_quad_v, n_quad_v);

               if (n_umus > 0) {
                    dmat_mul_d(wy1[1 - j][k].d1, wy1[1 - j][k].d2, w1.a, w1.d1, w1.d2, n_quad_v, n_umus_v, w3.d1, w3.d2, w1u1);
                    dmat_mul_d(wx1[1 - j].d1, wx1[1 - j].d2, wd1[k].a, wd1[k].d1, wd1[k].d2,  n_quad_v, n_umus_v, w4.d1, w4.d2, w1u1);
                    dmat_add_d(w3.d1, w3.d2, w4.d1, w4.d2, wy1[j][k].d1, wy1[j][k].d2, n_quad_v, n_umus_v);

                    dmat_trans_stokes(n_quad, n_stokes, w1.a, w2.a);
                    dmat_mul_d(wy2[1 - j][k].d1, wy2[1 - j][k].d2, w2.a, w2.d1, w2.d2, n_quad_v, n_umus_v, w3.d1, w3.d2, w1u1);
                    dmat_trans_stokes(n_quad, n_stokes, wd1[k].a, wd2[k].a);
                    dmat_mul_d(wx2[1 - j].d1, wx2[1 - j].d2, wd2[k].a, wd2[k].d1, wd2[k].d2,  n_quad_v, n_umus_v, w4.d1, w4.d2, w1u1);
                    dmat_add_d(w3.d1, w3.d2, w4.d1, w4.d2, wy2[j][k].d1, wy2[j][k].d2, n_quad_v, n_umus_v);

                    reduced_update_d(O12[k].d1, O12[k].d2, wy1[j][k].d1, wy1[j][k].d2, O12[k].d1, O12[k].d2, pade_c[i+6], n_quad_v, n_umus_v);
                    reduced_update_d(O21[k].d1, O21[k].d2, wy2[j][k].d1, wy2[j][k].d2, O21[k].d1, O21[k].d2, pade_c[i+6], n_quad_v, n_umus_v);
                    if (i+6 < pade_r) {
                         reduced_update_d(E12[k].d1, E12[k].d2, wy1[j][k].d1, wy1[j][k].d2, E12[k].d1, E12[k].d2, pade_c[i+7], n_quad_v, n_umus_v);
                         reduced_update_d(E21[k].d1, E21[k].d2, wy2[j][k].d1, wy2[j][k].d2, E21[k].d1, E21[k].d2, pade_c[i+7], n_quad_v, n_umus_v);
                    }
               }

          }


          dmat_mul(wx1[j].a, N12.a, n_quad_v, n_quad_v, n_quad_v, w1.a);
          dmat_add(w1.a, N11.a, N11.a, n_quad_v, n_quad_v);

          if (n_umus > 0) {
               dmat_mul_d(wx1[j].d1, wx1[j].d2, N12.a, N12.d1, N12.d2, n_quad_v, n_umus_v, w1.d1, w1.d2, w1u1);
               dmat_trans_stokes(n_quad, n_stokes, N12.a, N21.a);
               dmat_mul_d(wx2[j].d1, wx2[j].d2, N21.a, N21.d1, N21.d2, n_quad_v, n_umus_v, w2.d1, w2.d2, w1u1);

               dmat_add_d(w1.d1, w1.d2, N11.d1, N11.d2, N11.d1, N11.d2, n_quad_v, n_umus_v);
               dmat_add_d(w2.d1, w2.d2, N22.d1, N22.d2, N22.d1, N22.d2, n_quad_v, n_umus_v);
          }

          for (i = 0; i < n_derivs; ++i) {
               if (! derivs[i])
                    continue;

               dmat_mul(wy1[j][i].a, N12.a, n_quad_v, n_quad_v, n_quad_v, w1.a);
               dmat_add(w1.a, O11[i].a, O11[i].a, n_quad_v, n_quad_v);

               dmat_mul(wx1[j].a, O12[i].a, n_quad_v, n_quad_v, n_quad_v, w1.a);
               dmat_add(w1.a, O11[i].a, O11[i].a, n_quad_v, n_quad_v);

               if (n_umus > 0) {
                    dmat_mul_d(wy1[j][i].d1, wy1[j][i].d2, N12.a, N12.d1, N12.d2,  n_quad_v, n_umus_v, w3.d1, w3.d2, w1u1);
                    dmat_add_d(w3.d1, w3.d2, O11[i].d1, O11[i].d2, O11[i].d1, O11[i].d2, n_quad_v, n_umus_v);

                    dmat_mul_d(wx1[j].d1, wx1[j].d2, O12[i].a, O12[i].d1, O12[i].d2, n_quad_v, n_umus_v, w4.d1, w4.d2, w1u1);
                    dmat_add_d(w4.d1, w4.d2, O11[i].d1, O11[i].d2, O11[i].d1, O11[i].d2, n_quad_v, n_umus_v);

                    dmat_mul_d(wy2[j][i].d1, wy2[j][i].d2, N21.a, N21.d1, N21.d2,  n_quad_v, n_umus_v, w3.d1, w3.d2, w1u1);
                    dmat_add_d(w3.d1, w3.d2, O22[i].d1, O22[i].d2, O22[i].d1, O22[i].d2, n_quad_v, n_umus_v);

                    dmat_trans_stokes(n_quad, n_stokes, O12[i].a, O21[i].a);
                    dmat_mul_d(wx2[j].d1, wx2[j].d2, O21[i].a, O21[i].d1, O21[i].d2, n_quad_v, n_umus_v, w4.d1, w4.d2, w1u1);
                    dmat_add_d(w4.d1, w4.d2, O22[i].d1, O22[i].d2, O22[i].d1, O22[i].d2, n_quad_v, n_umus_v);
               }
          }


          dmat_mul(wx1[j].a, D12.a, n_quad_v, n_quad_v, n_quad_v, w1.a);
          dmat_add(w1.a, D11.a, D11.a, n_quad_v, n_quad_v);

          if (n_umus > 0) {
               dmat_mul_d(wx1[j].d1, wx1[j].d2, D12.a, D12.d1, D12.d2, n_quad_v, n_umus_v, w1.d1, w1.d2, w1u1);
               dmat_trans_stokes(n_quad, n_stokes, D12.a, D21.a);
               dmat_mul_d(wx2[j].d1, wx2[j].d2, D21.a, D21.d1, D21.d2, n_quad_v, n_umus_v, w2.d1, w2.d2, w1u1);

               dmat_add_d(w1.d1, w1.d2, D11.d1, D11.d2, D11.d1, D11.d2, n_quad_v, n_umus_v);
               dmat_add_d(w2.d1, w2.d2, D22.d1, D22.d2, D22.d1, D22.d2, n_quad_v, n_umus_v);
          }

          for (i = 0; i < n_derivs; ++i) {
               if (! derivs[i])
                    continue;

               dmat_mul(wy1[j][i].a, D12.a, n_quad_v, n_quad_v, n_quad_v, w1.a);
               dmat_add(w1.a, E11[i].a, E11[i].a, n_quad_v, n_quad_v);

               dmat_mul(wx1[j].a, E12[i].a, n_quad_v, n_quad_v, n_quad_v, w1.a);
               dmat_add(w1.a, E11[i].a, E11[i].a, n_quad_v, n_quad_v);

               if (n_umus > 0) {
                    dmat_mul_d(wy1[j][i].d1, wy1[j][i].d2, D12.a, D12.d1, D12.d2,  n_quad_v, n_umus_v, w3.d1, w3.d2, w1u1);
                    dmat_add_d(w3.d1, w3.d2, E11[i].d1, E11[i].d2, E11[i].d1, E11[i].d2, n_quad_v, n_umus_v);

                    dmat_mul_d(wx1[j].d1, wx1[j].d2, E12[i].a, E12[i].d1, E12[i].d2, n_quad_v, n_umus_v, w4.d1, w4.d2, w1u1);
                    dmat_add_d(w4.d1, w4.d2, E11[i].d1, E11[i].d2, E11[i].d1, E11[i].d2, n_quad_v, n_umus_v);

                    dmat_mul_d(wy2[j][i].d1, wy2[j][i].d2, D21.a, D21.d1, D21.d2,  n_quad_v, n_umus_v, w3.d1, w3.d2, w1u1);
                    dmat_add_d(w3.d1, w3.d2, E22[i].d1, E22[i].d2, E22[i].d1, E22[i].d2, n_quad_v, n_umus_v);

                    dmat_trans_stokes(n_quad, n_stokes, E12[i].a, E21[i].a);
                    dmat_mul_d(wx2[j].d1, wx2[j].d2, E21[i].a, E21[i].d1, E21[i].d2, n_quad_v, n_umus_v, w4.d1, w4.d2, w1u1);
                    dmat_add_d(w4.d1, w4.d2, E22[i].d1, E22[i].d2, E22[i].d1, E22[i].d2, n_quad_v, n_umus_v);
               }
          }
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     else {
          dmat_init(N11.a, pade_c[0], n_quad_v, n_quad_v);
          dmat_init(D11.a, pade_c[1], n_quad_v, n_quad_v);

          dmat_mul(A12.a, A21.a, n_quad_v, n_quad_v, n_quad_v, w1.a);

          reduced_update(N11.a, w1.a, N11.a, pade_c[2], n_quad_v, n_quad_v);
          if (2 < pade_r)
               reduced_update(D11.a, w1.a, D11.a, pade_c[3], n_quad_v, n_quad_v);

          dmat_copy(wx1[0].a, w1.a, n_quad_v, n_quad_v);

          if (n_umus > 0) {
               dmat_init_d(N11.d1, N11.d2, pade_c[0], n_quad_v, n_umus_v);
               dmat_init_d(N22.d1, N22.d2, pade_c[0], n_quad_v, n_umus_v);

               dmat_init_d(D11.d1, D11.d2, pade_c[1], n_quad_v, n_umus_v);
               dmat_init_d(D22.d1, D22.d2, pade_c[1], n_quad_v, n_umus_v);

               dmat_mul_d(A12.d1, A12.d2, A21.a, A21.d1, A21.d2, n_quad_v, n_umus_v, w1.d1, w1.d2, w1u1);
               dmat_mul_d(A21.d1, A21.d2, A12.a, A12.d1, A12.d2, n_quad_v, n_umus_v, w2.d1, w2.d2, w1u1);

               reduced_update_d(N11.d1, N11.d2, w1.d1, w1.d2, N11.d1, N11.d2, pade_c[2], n_quad_v, n_umus_v);
               reduced_update_d(N22.d1, N22.d2, w2.d1, w2.d2, N22.d1, N22.d2, pade_c[2], n_quad_v, n_umus_v);

               if (2 < pade_r) {
                    reduced_update_d(D11.d1, D11.d2, w1.d1, w1.d2, D11.d1, D11.d2, pade_c[3], n_quad_v, n_umus_v);
                    reduced_update_d(D22.d1, D22.d2, w2.d1, w2.d2, D22.d1, D22.d2, pade_c[3], n_quad_v, n_umus_v);
               }

               dmat_copy_d(wx1[0].d1, wx1[0].d2, w1.d1, w1.d2, n_quad_v, n_umus_v);
               dmat_copy_d(wx2[0].d1, wx2[0].d2, w2.d1, w2.d2, n_quad_v, n_umus_v);
          }

          for (i = 0; i < n_derivs; ++i) {
               if (! derivs[i])
                    continue;

               dmat_zero(O11[i].a, n_quad_v, n_quad_v);
               dmat_zero(E11[i].a, n_quad_v, n_quad_v);

               dmat_mul(B12[i].a, A21.a, n_quad_v, n_quad_v, n_quad_v, w3.a);
               dmat_mul(A12.a, B21[i].a, n_quad_v, n_quad_v, n_quad_v, w4.a);
               dmat_add(w3.a, w4.a, wd1[i].a, n_quad_v, n_quad_v);

               reduced_update(O11[i].a, wd1[i].a, O11[i].a, pade_c[2], n_quad_v, n_quad_v);
               if (2 < pade_r)
                    reduced_update(E11[i].a, wd1[i].a, E11[i].a, pade_c[3], n_quad_v, n_quad_v);

               dmat_copy(wy1[0][i].a, wd1[i].a, n_quad_v, n_quad_v);

               if (n_umus > 0) {
                    dmat_zero_d(O11[i].d1, O11[i].d2, n_quad_v, n_umus_v);
                    dmat_zero_d(O22[i].d1, O22[i].d2, n_quad_v, n_umus_v);

                    dmat_zero_d(E11[i].d1, E11[i].d2, n_quad_v, n_umus_v);
                    dmat_zero_d(E22[i].d1, E22[i].d2, n_quad_v, n_umus_v);


                    dmat_mul_d(B12[i].d1, B12[i].d2, A21.a, A21.d1, A21.d2,    n_quad_v, n_umus_v, w3.d1, w3.d2, w1u1);
                    dmat_mul_d(A12.d1, A12.d2, B21[i].a, B21[i].d1, B21[i].d2, n_quad_v, n_umus_v, w4.d1, w4.d2, w1u1);
                    dmat_add_d(w3.d1, w3.d2, w4.d1, w4.d2, wd1[i].d1, wd1[i].d2, n_quad_v, n_umus_v);

                    reduced_update_d(O11[i].d1, O11[i].d2, wd1[i].d1, wd1[i].d2, O11[i].d1, O11[i].d2, pade_c[2], n_quad_v, n_umus_v);
                    if (2 < pade_r)
                         reduced_update_d(E11[i].d1, E11[i].d2, wd1[i].d1, wd1[i].d2, E11[i].d1, E11[i].d2, pade_c[3], n_quad_v, n_umus_v);

                    dmat_copy_d(wy1[0][i].d1, wy1[0][i].d2, wd1[i].d1, wd1[i].d2, n_quad_v, n_umus_v);


                    dmat_mul_d(B21[i].d1, B21[i].d2, A12.a, A12.d1, A12.d2,    n_quad_v, n_umus_v, w3.d1, w3.d2, w1u1);
                    dmat_mul_d(A21.d1, A21.d2, B12[i].a, B12[i].d1, B12[i].d2, n_quad_v, n_umus_v, w4.d1, w4.d2, w1u1);
                    dmat_add_d(w3.d1, w3.d2, w4.d1, w4.d2, wd2[i].d1, wd2[i].d2, n_quad_v, n_umus_v);

                    reduced_update_d(O22[i].d1, O22[i].d2, wd2[i].d1, wd2[i].d2, O22[i].d1, O22[i].d2, pade_c[2], n_quad_v, n_umus_v);
                    if (2 < pade_r)
                         reduced_update_d(E22[i].d1, E22[i].d2, wd2[i].d1, wd2[i].d2, E22[i].d1, E22[i].d2, pade_c[3], n_quad_v, n_umus_v);

                    dmat_copy_d(wy2[0][i].d1, wy2[0][i].d2, wd2[i].d1, wd2[i].d2, n_quad_v, n_umus_v);
               }
          }

          j = 1;
          for (i = 4; i <= pade_r; i += 2) {
               dmat_mul(wx1[1 - j].a, w1.a, n_quad_v, n_quad_v, n_quad_v, wx1[j].a);

               reduced_update(N11.a, wx1[j].a, N11.a, pade_c[i], n_quad_v, n_quad_v);
               if (i < pade_r)
                    reduced_update(D11.a, wx1[j].a, D11.a, pade_c[i+1], n_quad_v, n_quad_v);

               if (n_umus > 0) {
                    dmat_mul_d(wx1[1 - j].d1, wx1[1 - j].d2, w1.a, w1.d1, w1.d2, n_quad_v, n_umus_v, wx1[j].d1, wx1[j].d2, w1u1);
                    dmat_trans_stokes(n_quad, n_stokes, w1.a, w2.a);
                    dmat_mul_d(wx2[1 - j].d1, wx2[1 - j].d2, w2.a, w2.d1, w2.d2, n_quad_v, n_umus_v, wx2[j].d1, wx2[j].d2, w1u1);

                    reduced_update_d(N11.d1, N11.d2, wx1[j].d1, wx1[j].d2, N11.d1, N11.d2, pade_c[i], n_quad_v, n_umus_v);
                    reduced_update_d(N22.d1, N22.d2, wx2[j].d1, wx2[j].d2, N22.d1, N22.d2, pade_c[i], n_quad_v, n_umus_v);

                    if (i < pade_r) {
                         reduced_update_d(D11.d1, D11.d2, wx1[j].d1, wx1[j].d2, D11.d1, D11.d2, pade_c[i+1], n_quad_v, n_umus_v);
                         reduced_update_d(D22.d1, D22.d2, wx2[j].d1, wx2[j].d2, D22.d1, D22.d2, pade_c[i+1], n_quad_v, n_umus_v);
                    }
               }

               for (k = 0; k < n_derivs; ++k) {
                    if (! derivs[k])
                         continue;

                    dmat_mul(wy1[1 - j][k].a, w1.a, n_quad_v, n_quad_v, n_quad_v, w3.a);
                    dmat_mul(wx1[1 - j].a, wd1[k].a, n_quad_v, n_quad_v, n_quad_v, w4.a);
                    dmat_add(w3.a, w4.a, wy1[j][k].a, n_quad_v, n_quad_v);

                    reduced_update(O11[k].a, wy1[j][k].a, O11[k].a, pade_c[i], n_quad_v, n_quad_v);
                    if (i < pade_r)
                         reduced_update(E11[k].a, wy1[j][k].a, E11[k].a, pade_c[i+1], n_quad_v, n_quad_v);

                    if (n_umus > 0) {
                         dmat_mul_d(wy1[1 - j][k].d1, wy1[1 - j][k].d2, w1.a, w1.d1, w1.d2, n_quad_v, n_umus_v, w3.d1, w3.d2, w1u1);
                         dmat_mul_d(wx1[1 - j].d1, wx1[1 - j].d2, wd1[k].a, wd1[k].d1, wd1[k].d2, n_quad_v, n_umus_v, w4.d1, w4.d2, w1u1);
                         dmat_add_d(w3.d1, w3.d2, w4.d1, w4.d2, wy1[j][k].d1, wy1[j][k].d2, n_quad_v, n_umus_v);

                         reduced_update_d(O11[k].d1, O11[k].d2, wy1[j][k].d1, wy1[j][k].d2, O11[k].d1, O11[k].d2, pade_c[i], n_quad_v, n_umus_v);
                         if (i < pade_r)
                              reduced_update_d(E11[k].d1, E11[k].d2, wy1[j][k].d1, wy1[j][k].d2, E11[k].d1, E11[k].d2, pade_c[i+1], n_quad_v, n_umus_v);

                         dmat_trans_stokes(n_quad, n_stokes, w1.a, w2.a);
                         dmat_mul_d(wy2[1 - j][k].d1, wy2[1 - j][k].d2, w2.a, w2.d1, w2.d2, n_quad_v, n_umus_v, w3.d1, w3.d2, w1u1);
                         dmat_trans_stokes(n_quad, n_stokes, wd1[k].a, wd2[k].a);
                         dmat_mul_d(wx2[1 - j].d1, wx2[1 - j].d2, wd2[k].a, wd2[k].d1, wd2[k].d2, n_quad_v, n_umus_v, w4.d1, w4.d2, w1u1);
                         dmat_add_d(w3.d1, w3.d2, w4.d1, w4.d2, wy2[j][k].d1, wy2[j][k].d2, n_quad_v, n_umus_v);

                         reduced_update_d(O22[k].d1, O22[k].d2, wy2[j][k].d1, wy2[j][k].d2, O22[k].d1, O22[k].d2, pade_c[i], n_quad_v, n_umus_v);
                         if (i < pade_r)
                              reduced_update_d(E22[k].d1, E22[k].d2, wy2[j][k].d1, wy2[j][k].d2, E22[k].d1, E22[k].d2, pade_c[i+1], n_quad_v, n_umus_v);
                    }
               }

               j = 1 - j;
          }
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     dmat_trans_stokes(n_quad, n_stokes, N11.a, N22.a);
     dmat_trans_stokes(n_quad, n_stokes, D11.a, D22.a);
if (n_stokes <= 2) {
     dsym_mul(A12.a, D22.a, n_quad_v, n_quad_v, N12.a, 1);
     dsym_mul(A21.a, D11.a, n_quad_v, n_quad_v, N21.a, 1);
}
else {
     dsym_mul(A12.a, D22.a, n_quad_v, n_quad_v, N12.a, 0);
     dmat_sym_stokes(n_quad, n_stokes, N12.a);
     dsym_mul(A21.a, D11.a, n_quad_v, n_quad_v, N21.a, 0);
     dmat_sym_stokes(n_quad, n_stokes, N21.a);
}
     if (n_umus > 0) {
          dmat_mul_d(A12.d1, A12.d2, D22.a, D22.d1, D22.d2, n_quad_v, n_umus_v, N12.d1, N12.d2, w1u1);
          dmat_mul_d(A21.d1, A21.d2, D11.a, D11.d1, D11.d2, n_quad_v, n_umus_v, N21.d1, N21.d2, w1u1);
     }

     for (i = 0; i < n_derivs; ++i) {
          if (! derivs[i])
               continue;

          dmat_trans_stokes(n_quad, n_stokes, O11[i].a, O22[i].a);
          dmat_trans_stokes(n_quad, n_stokes, E11[i].a, E22[i].a);
if (n_stokes <= 2) {
          dsym_mul(B12[i].a, D22.a, n_quad_v, n_quad_v, w1.a, 0);
          dsym_mul(A12.a, E22[i].a, n_quad_v, n_quad_v, w2.a, 0);
          dsym_add(w1.a, w2.a, O12[i].a, n_quad_v, 1);

          dsym_mul(B21[i].a, D11.a, n_quad_v, n_quad_v, w1.a, 0);
          dsym_mul(A21.a, E11[i].a, n_quad_v, n_quad_v, w2.a, 0);
          dsym_add(w1.a, w2.a, O21[i].a, n_quad_v, 1);
}
else {
          dsym_mul(B12[i].a, D22.a, n_quad_v, n_quad_v, w1.a, 0);
          dsym_mul(A12.a, E22[i].a, n_quad_v, n_quad_v, w2.a, 0);
          dsym_add(w1.a, w2.a, O12[i].a, n_quad_v, 0);
          dmat_sym_stokes(n_quad, n_stokes, O12[i].a);

          dsym_mul(B21[i].a, D11.a, n_quad_v, n_quad_v, w1.a, 0);
          dsym_mul(A21.a, E11[i].a, n_quad_v, n_quad_v, w2.a, 0);
          dsym_add(w1.a, w2.a, O21[i].a, n_quad_v, 0);
          dmat_sym_stokes(n_quad, n_stokes, O21[i].a);
}
          if (n_umus > 0) {
               dmat_mul_d(B12[i].d1, B12[i].d2, D22.a, D22.d1, D22.d2,    n_quad_v, n_umus_v, w1.d1, w1.d2, w1u1);
               dmat_mul_d(A12.d1, A12.d2, E22[i].a, E22[i].d1, E22[i].d2, n_quad_v, n_umus_v, w2.d1, w2.d2, w1u1);
               dmat_add_d(w1.d1, w1.d2, w2.d1, w2.d2, O12[i].d1, O12[i].d2, n_quad_v, n_umus_v);

               dmat_mul_d(B21[i].d1, B21[i].d2, D11.a, D11.d1, D11.d2,    n_quad_v, n_umus_v, w1.d1, w1.d2, w1u1);
               dmat_mul_d(A21.d1, A21.d2, E11[i].a, E11[i].d1, E11[i].d2, n_quad_v, n_umus_v, w2.d1, w2.d2, w1u1);
               dmat_add_d(w1.d1, w1.d2, w2.d1, w2.d2, O21[i].d1, O21[i].d2, n_quad_v, n_umus_v);
          }
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
/*
     for (i = 0; i < n_quad_v; ++i) {
          for (j = 0; j < n_quad_v; ++j) {
               D11[i][j] =  N11[i][j];
               D12[i][j] = -N12[i][j];
               D21[i][j] = -N21[i][j];
               D22[i][j] =  N22[i][j];
          }
     }

     for (i = 0; i < n_derivs; ++i) {
          if (! derivs[i])
               continue;

          for (j = 0; j < n_quad_v; ++j) {
               for (k = 0; k < n_quad_v; ++k) {
                    E11[i].a[j][k] =  O11[i].a[j][k];
                    E12[i].a[j][k] = -O12[i].a[j][k];
                    E21[i].a[j][k] = -O21[i].a[j][k];
                    E22[i].a[j][k] =  O22[i].a[j][k];
               }
          }
     }
*/
}



/*******************************************************************************
 *
 ******************************************************************************/
void pade_calc_R_and_T_2(matd1d2 N11, matd1d2 N12, matd1d2 N21, matd1d2 N22,
                         matd1d2 D11, matd1d2 D12, matd1d2 D21, matd1d2 D22,
                         matd1d2 R_p, matd1d2 T_p, matd1d2 R_m, matd1d2 T_m,
                         matd1d2 *O11, matd1d2 *O12, matd1d2 *O21, matd1d2 *O22,
                         matd1d2 *E11, matd1d2 *E12, matd1d2 *E21, matd1d2 *E22,
                         matd1d2 *R_p_l, matd1d2 *T_p_l, matd1d2 *R_m_l, matd1d2 *T_m_l,
                         int n_quad, int n_stokes, int n_derivs,
                         int check_condition, int reduced, int symmetric, int vector,
                         double *condition, uchar *derivs, work_data work, int n_umus) {

     int i;

     int n_quad_v;
     int n_umus_v;
int symmetric2;
     int *i1;
     int *i2;
     int *i3;

     double aa;

     double **w1;
     double **w2;

     double **lu1;
     double **lu2;
     double **b;
     double **c;
     double **d;

     double **a_l;
     double **b_l;
     double **c_l;
     double **d_l;

     double **lu3;
     double **e;
     double **f;

     double **e_l;
     double **f_l;

     double **w1u1;
     double  *w1u2;
     double **w2u1;
     double  *w2u2;
     double **w3u1;

     double **au1;
     double  *au2;
     double **cu1;
     double  *cu2;
     double **du1;
     double  *du2;

     double **au1_l;
     double  *au2_l;
     double **bu1_l;
     double  *bu2_l;
     double **cu1_l;
     double  *cu2_l;
     double **du1_l;
     double  *du2_l;
symmetric2 = symmetric;
symmetric  = 0;
     n_quad_v = n_quad * n_stokes;
     n_umus_v = n_umus * n_stokes;


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     i2  = get_work_i1(&work, n_quad_v);

     w1  = get_work_d2(&work, n_quad_v, n_quad_v);

     lu1 = get_work_d2(&work, n_quad_v, n_quad_v);
     lu2 = get_work_d2(&work, n_quad_v, n_quad_v);
     c   = get_work_d2(&work, n_quad_v, n_quad_v);
     d   = get_work_d2(&work, n_quad_v, n_quad_v);

     if (! symmetric) {
          b  = get_work_d2(&work, n_quad_v, n_quad_v);

          i1 = get_work_i1(&work, n_quad_v);
     }

     if (n_umus > 0) {
          w1u1 = get_work_d2(&work, n_umus_v, n_quad_v);
          w1u2 = get_work_d1(&work, n_umus_v);
          w2u1 = get_work_d2(&work, n_umus_v, n_quad_v);
          w2u2 = get_work_d1(&work, n_umus_v);
          w3u1 = get_work_d2(&work, n_umus_v, n_quad_v);

          au1  = get_work_d2(&work, n_umus_v, n_quad_v);
          au2  = get_work_d1(&work, n_umus_v);

          cu1  = get_work_d2(&work, n_umus_v, n_quad_v);
          cu2  = get_work_d1(&work, n_umus_v);
          du1  = get_work_d2(&work, n_umus_v, n_quad_v);
          du2  = get_work_d1(&work, n_umus_v);
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (reduced) {
          D22 = N11;
          D21 = N12;
          D12 = N21;
          D11 = N22;
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     /* lu1 = LU D11 */
     dmat_copy(lu1, D11.a, n_quad_v, n_quad_v);
if (! symmetric) {
     dmat_getrf(lu1, n_quad_v, n_quad_v, i1);
     if (check_condition) {
          aa = dmat_p_one_norm(D11.a, n_quad_v, n_quad_v);
          aa = dmat_gecon('1', lu1, n_quad_v, aa, i2, *w1);
     }
}
else {
     dmat_potrf(lu1, n_quad_v);
     if (check_condition) {
          aa = dmat_p_one_norm(D11.a, n_quad_v, n_quad_v);
          aa = dmat_pocon(lu1, n_quad_v, aa, i2, *w1);
     }
}
     if (check_condition) {
          *condition = aa;
/*
          if (aa < DBL_EPSILON) {
               fprintf(stderr, "ERROR: pade condition < machine precision\n");
               return -1;
          }
*/
     }

     /* c = D11-1 * N12 */
     dmat_trans(N12.a, w1, n_quad_v, n_quad_v);
if (! symmetric)
     dmat_getrs(lu1, w1, n_quad_v, n_quad_v, i1);
else
     dmat_potrs(lu1, w1, n_quad_v, n_quad_v);
     dmat_trans(w1, c, n_quad_v, n_quad_v);

     /* d = D11-1 * D12 */
     dmat_trans(D12.a, w1, n_quad_v, n_quad_v);
if (! symmetric)
     dmat_getrs(lu1, w1, n_quad_v, n_quad_v, i1);
else
     dmat_potrs(lu1, w1, n_quad_v, n_quad_v);
     dmat_trans(w1, d, n_quad_v, n_quad_v);

     /* lu2 = LU a = (D21 * c - N22) */
     dmat_mul(D21.a, c, n_quad_v, n_quad_v, n_quad_v, lu2);
     dmat_sub(lu2, N22.a, lu2, n_quad_v, n_quad_v);
     dmat_getrf(lu2, n_quad_v, n_quad_v, i2);
if (! symmetric2) {
     /*          b = (D21 * d - D22) */
     dmat_mul(D21.a, d, n_quad_v, n_quad_v, n_quad_v, b);
     dmat_sub(b, D22.a, b, n_quad_v, n_quad_v);

     /* T = a-1 * b */
     dmat_trans(b, w1, n_quad_v, n_quad_v);
     dmat_getrs(lu2, w1, n_quad_v, n_quad_v, i2);
     dmat_trans(w1, T_m.a, n_quad_v, n_quad_v);
}
else {
     /*          b = (D21 * d - D22) */
     dsym_m_s(D21.a, d, D22.a, n_quad_v, n_quad_v, T_m.a, 1);

     /* T = a-1 * b */
     dmat_getrs(lu2, T_m.a, n_quad_v, n_quad_v, i2);
}

     /* R = c * T - d */
if (! symmetric2) {
     dmat_mul(c, T_m.a, n_quad_v, n_quad_v, n_quad_v, R_m.a);
     dmat_sub(R_m.a, d, R_m.a, n_quad_v, n_quad_v);
}
else
     dsym_m_s(c, T_m.a, d, n_quad_v, n_quad_v, R_m.a, 1);

     if (n_umus > 0) {
          dmat_inv_mul_d(c, D11.d1, D11.d2, N12.a, N12.d1, N12.d2, n_quad_v, n_umus_v, cu1, cu2, w1u1);

          dmat_inv_mul_d(d, D11.d1, D11.d2, D12.a, D12.d1, D12.d2, n_quad_v, n_umus_v, du1, du2, w1u1);

          dmat_mul_d(D21.d1, D21.d2, c, cu1, cu2, n_quad_v, n_umus_v, w1u1, w1u2, w3u1);
          dmat_sub_d(w1u1, w1u2, N22.d1, N22.d2, w1u1, w1u2, n_quad_v, n_umus_v);
dmat_copy_d(au1, au2, w1u1, w1u2, n_quad_v, n_umus_v);
          dmat_mul_d(D21.d1, D21.d2, d, du1, du2, n_quad_v, n_umus_v, w2u1, w2u2, w3u1);
          dmat_sub_d(w2u1, w2u2, D22.d1, D22.d2, w2u1, w2u2, n_quad_v, n_umus_v);

          dmat_inv_mul_d(T_m.a, w1u1, w1u2, b, w2u1, w2u2, n_quad_v, n_umus_v, T_m.d1, T_m.d2, w3u1);

          dmat_mul_d(cu1, cu2, T_m.a, T_m.d1, T_m.d2, n_quad_v, n_umus_v, R_m.d1, R_m.d2, w3u1);
          dmat_sub_d(R_m.d1, R_m.d2, du1, du2, R_m.d1, R_m.d2, n_quad_v, n_umus_v);
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (flags_or(derivs, n_derivs)) {
          w2  = get_work_d2(&work, n_quad_v, n_quad_v);

          a_l = get_work_d2(&work, n_quad_v, n_quad_v);
          b_l = get_work_d2(&work, n_quad_v, n_quad_v);
          c_l = get_work_d2(&work, n_quad_v, n_quad_v);
          d_l = get_work_d2(&work, n_quad_v, n_quad_v);

          if (n_umus > 0) {
               au1_l  = get_work_d2(&work, n_umus_v, n_quad_v);
               au2_l  = get_work_d1(&work, n_umus_v);
               bu1_l  = get_work_d2(&work, n_umus_v, n_quad_v);
               bu2_l  = get_work_d1(&work, n_umus_v);
               cu1_l  = get_work_d2(&work, n_umus_v, n_quad_v);
               cu2_l  = get_work_d1(&work, n_umus_v);
               du1_l  = get_work_d2(&work, n_umus_v, n_quad_v);
               du2_l  = get_work_d1(&work, n_umus_v);
          }

          if (reduced) {
               E22 = O11;
               E21 = O12;
               E12 = O21;
               E11 = O22;
          }

          for (i = 0; i < n_derivs; ++i) {

               if (! derivs[i])
                    continue;

               /* L(c) = D11-1 * (-L(D11) * c + L(N12)) */
               dmat_mul(E11[i].a, c, n_quad_v, n_quad_v, n_quad_v, w1);
               dmat_sub(O12[i].a, w1, w1, n_quad_v, n_quad_v);
               dmat_trans(w1, w2, n_quad_v, n_quad_v);
if (! symmetric)
               dmat_getrs(lu1, w2, n_quad_v, n_quad_v, i1);
else
               dmat_potrs(lu1, w2, n_quad_v, n_quad_v);
               dmat_trans(w2, c_l, n_quad_v, n_quad_v);

               if (n_umus > 0) {
                    /* L(c) = D11-1 * (-L(D11) * c + L(N12)) */
                    dmat_mul_d(E11[i].d1, E11[i].d2, c, cu1, cu2, n_quad_v, n_umus_v, w1u1, w1u2, w3u1);
                    dmat_sub_d(O12[i].d1, O12[i].d2, w1u1, w1u2, w1u1, w1u2, n_quad_v, n_umus_v);

                    dmat_inv_mul_d(c_l, D11.d1, D11.d2, w1, w1u1, w1u2, n_quad_v, n_umus_v, cu1_l, cu2_l, w3u1);
               }

               /* L(d) = D11-1 * (-L(D11) * d + L(D12)) */
               dmat_mul(E11[i].a, d, n_quad_v, n_quad_v, n_quad_v, w1);
               dmat_sub(E12[i].a, w1, w1, n_quad_v, n_quad_v);
               dmat_trans(w1, w2, n_quad_v, n_quad_v);
if (! symmetric)
               dmat_getrs(lu1, w2, n_quad_v, n_quad_v, i1);
else
               dmat_potrs(lu1, w2, n_quad_v, n_quad_v);
               dmat_trans(w2, d_l, n_quad_v, n_quad_v);

               if (n_umus > 0) {
                    /* L(d) = D11-1 * (-L(D11) * d + L(D12)) */
                    dmat_mul_d(E11[i].d1, E11[i].d2, d, du1, du2, n_quad_v, n_umus_v, w1u1, w1u2, w3u1);
                    dmat_sub_d(E12[i].d1, E12[i].d2, w1u1, w1u2, w1u1, w1u2, n_quad_v, n_umus_v);

                    dmat_inv_mul_d(d_l, D11.d1, D11.d2, w1, w1u1, w1u2, n_quad_v, n_umus_v, du1_l, du2_l, w3u1);
               }

               /* L(a) = L(D21) * c + D21 * L(c) - L(N22) */
               dmat_mul(E21[i].a, c, n_quad_v, n_quad_v, n_quad_v, w1);
               dmat_mul(D21.a, c_l, n_quad_v, n_quad_v, n_quad_v, w2);
               dmat_add(w1, w2, w1, n_quad_v, n_quad_v);
               dmat_sub(w1, O22[i].a, a_l, n_quad_v, n_quad_v);

               if (n_umus > 0) {
                    /* L(a) = L(D21) * c + D21 * L(c) - L(N22) */
                    dmat_mul_d(E21[i].d1, E21[i].d2, c, cu1, cu2, n_quad_v, n_umus_v, w1u1, w1u2, w3u1);
                    dmat_mul_d(D21.d1, D21.d2, c_l, cu1_l, cu2_l, n_quad_v, n_umus_v, w2u1, w2u2, w3u1);
                    dmat_add_d(w1u1, w1u2, w2u1, w2u2, w1u1, w1u2, n_quad_v, n_umus_v);
                    dmat_sub_d(w1u1, w1u2, O22[i].d1, O22[i].d2, au1_l, au2_l, n_quad_v, n_umus_v);
               }

               /* L(b) = L(D21) * d + D21 * L(d) - L(D22) */
if (! symmetric2) {
               dmat_mul(E21[i].a, d, n_quad_v, n_quad_v, n_quad_v, w1);
               dmat_mul(D21.a, d_l, n_quad_v, n_quad_v, n_quad_v, w2);
               dmat_add(w1, w2, b_l, n_quad_v, n_quad_v);
               dmat_sub(b_l, E22[i].a, b_l, n_quad_v, n_quad_v);
}
else {
               dsym_mma(E21[i].a, d, D21.a, d_l, n_quad_v, n_quad_v, b_l, 0);
               dsym_sub(b_l, E22[i].a, b_l, n_quad_v, 1);
}
               if (n_umus > 0) {
                    /* L(b) = L(D21) * d + D21 * L(d) - L(D22) */
                    dmat_mul_d(E21[i].d1, E21[i].d2, d, du1, du2, n_quad_v, n_umus_v, w1u1, w1u2, w3u1);
                    dmat_mul_d(D21.d1, D21.d2, d_l, du1_l, du2_l, n_quad_v, n_umus_v, w2u1, w2u2, w3u1);
                    dmat_add_d(w1u1, w1u2, w2u1, w2u2, w1u1, w1u2, n_quad_v, n_umus_v);
                    dmat_sub_d(w1u1, w1u2, E22[i].d1, E22[i].d2, bu1_l, bu2_l, n_quad_v, n_umus_v);
               }

               /* L(T) = a-1 * (-L(a) * T + L(b)) */
               dmat_mul(a_l, T_m.a, n_quad_v, n_quad_v, n_quad_v, w1);
               dmat_sub(b_l, w1, w1, n_quad_v, n_quad_v);
               dmat_trans(w1, w2, n_quad_v, n_quad_v);
               dmat_getrs(lu2, w2, n_quad_v, n_quad_v, i2);
               dmat_trans(w2, T_m_l[i].a, n_quad_v, n_quad_v);

               if (n_umus > 0) {
                    /* L(T) = a-1 * (-L(a) * T + L(b)) */
                    dmat_mul_d(au1_l, au2_l, T_m.a, T_m.d1, T_m.d2, n_quad_v, n_umus_v, w1u1, w1u2, w3u1);
                    dmat_sub_d(bu1_l, bu2_l, w1u1, w1u2, w1u1, w1u2, n_quad_v, n_umus_v);

                    dmat_inv_mul_d(T_m_l[i].a, au1, au2, w1, w1u1, w1u2, n_quad_v, n_umus_v, T_m_l[i].d1, T_m_l[i].d2, w3u1);
               }

               /* L(R) = L(c) * T + c * L(T) - L(d) */
if (! symmetric2) {
               dmat_mul(c_l, T_m.a, n_quad_v, n_quad_v, n_quad_v, w1);
               dmat_mul(c, T_m_l[i].a, n_quad_v, n_quad_v, n_quad_v, w2);
               dmat_add(w1, w2, R_m_l[i].a, n_quad_v, n_quad_v);
               dmat_sub(R_m_l[i].a, d_l, R_m_l[i].a, n_quad_v, n_quad_v);
}
else {
               dsym_mma(c_l, T_m.a, c, T_m_l[i].a, n_quad_v, n_quad_v, R_m_l[i].a, 0);
               dsym_sub(R_m_l[i].a, d_l, R_m_l[i].a, n_quad_v, 1);
}
               if (n_umus > 0) {
                    /* L(R) = L(c) * T + c * L(T) - L(d) */
                    dmat_mul_d(cu1_l, cu2_l, T_m.a, T_m.d1, T_m.d2, n_quad_v, n_umus_v, w1u1, w1u2, w3u1);
                    dmat_mul_d(cu1, cu2, T_m_l[i].a, T_m_l[i].d1, T_m_l[i].d2, n_quad_v, n_umus_v, w2u1, w2u2, w3u1);
                    dmat_add_d(w1u1, w1u2, w2u1, w2u2, R_m_l[i].d1, R_m_l[i].d2, n_quad_v, n_umus_v);
                    dmat_sub_d(R_m_l[i].d1, R_m_l[i].d2, du1_l, du2_l, R_m_l[i].d1, R_m_l[i].d2, n_quad_v, n_umus_v);
               }
          }
     }

     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (vector) {
if (0) {
          i3  = get_work1(&work, WORK_IX);

          lu3 = get_work1(&work, WORK_DXX);
          e   = get_work1(&work, WORK_DXX);
          f   = get_work1(&work, WORK_DXX);

          /* f = D11-1 * N11 */
          dmat_trans(N11.a, w1, n_quad_v, n_quad_v);
if (! symmetric)
          dmat_getrs(lu1, w1, n_quad_v, n_quad_v, i1);
else
          dmat_potrs(lu1, w1, n_quad_v, n_quad_v);
          dmat_trans(w1, f, n_quad_v, n_quad_v);

          /* lu3 = LU b */
          dmat_copy(lu3, b, n_quad_v, n_quad_v);
          dmat_getrf(lu3, n_quad_v, n_quad_v, i3);

          /* e = b-1 * (N12 - D21 * f) */
          dmat_mul(D21.a, f, n_quad_v, n_quad_v, n_quad_v, w1);
          dmat_sub(N21.a, w1, e, n_quad_v, n_quad_v);
          dmat_trans(e, w1, n_quad_v, n_quad_v);
          dmat_getrs(lu3, w1, n_quad_v, n_quad_v, i3);
          dmat_trans(w1, e, n_quad_v, n_quad_v);

          /* T_p = f + (R_m + d) * e */
          dmat_add(R_m.a, d, w1, n_quad_v, n_quad_v);
          dmat_mul(w1, e, n_quad_v, n_quad_v, n_quad_v, T_p.a);
          dmat_add(f, T_p.a, T_p.a, n_quad_v, n_quad_v);

          /* R_p = T_m * e */
          dmat_mul(T_m.a, e, n_quad_v, n_quad_v, n_quad_v, R_p.a);
}
else {
          phase_matrix_symmetry2(n_quad, n_stokes, T_m.a, R_m.a, T_p.a, R_p.a, 1.);
}
          /*--------------------------------------------------------------------
           *
           *------------------------------------------------------------------*/
          if (flags_or(derivs, n_derivs)) {
if (0) {
               e_l = get_work1(&work, WORK_DXX);
               f_l = get_work1(&work, WORK_DXX);

               for (i = 0; i < n_derivs; ++i) {

                    if (! derivs[i])
                         continue;

                    /* L(f) = D11-1 * (-L(D11) * f + L(N11)) */
                    dmat_mul(E11[i].a, f, n_quad_v, n_quad_v, n_quad_v, w1);
                    dmat_sub(O11[i].a, w1, w1, n_quad_v, n_quad_v);
                    dmat_trans(w1, w2, n_quad_v, n_quad_v);
if (! symmetric)
                    dmat_getrs(lu1, w2, n_quad_v, n_quad_v, i1);
else
                    dmat_potrs(lu1, w2, n_quad_v, n_quad_v);
                    dmat_trans(w2, f_l, n_quad_v, n_quad_v);

                    /* L(e) = b-1 * (-L(b) * e + L(N21) - L(D21) * f - D21 * L(f)) */
                    dmat_mul(b_l, e, n_quad_v, n_quad_v, n_quad_v, w1);
                    dmat_sub(O21[i].a, w1, w1, n_quad_v, n_quad_v);
                    dmat_mul(E21[i].a, f, n_quad_v, n_quad_v, n_quad_v, w2);
                    dmat_sub(w1, w2, w1, n_quad_v, n_quad_v);
                    dmat_mul(D21.a, f_l, n_quad_v, n_quad_v, n_quad_v, w2);
                    dmat_sub(w1, w2, w1, n_quad_v, n_quad_v);
                    dmat_trans(w1, w2, n_quad_v, n_quad_v);
                    dmat_getrs(lu3, w2, n_quad_v, n_quad_v, i3);
                    dmat_trans(w2, e_l, n_quad_v, n_quad_v);

                    /* L(T_p) = L(f) + (L(R_m) + L(d)) * e + (R_m + d) * L(e) */
                    dmat_add(R_m_l[i].a, d_l, w1, n_quad_v, n_quad_v);
                    dmat_mul(w1, e, n_quad_v, n_quad_v, n_quad_v, w2);
                    dmat_add(f_l, w2, T_p_l[i].a, n_quad_v, n_quad_v);
                    dmat_add(R_m.a, d, w1, n_quad_v, n_quad_v);
                    dmat_mul(w1, e_l, n_quad_v, n_quad_v, n_quad_v, w2);
                    dmat_add(T_p_l[i].a, w2, T_p_l[i].a, n_quad_v, n_quad_v);

                    /* L(R_p) = L(T_m) * e + T_m * L(e)) */
                    dmat_mul(T_m_l[i].a, e, n_quad_v, n_quad_v, n_quad_v, w1);
                    dmat_mul(T_m.a, e_l, n_quad_v, n_quad_v, n_quad_v, w2);
                    dmat_add(w1, w2, R_p_l[i].a, n_quad_v, n_quad_v);
               }
}
else {
               for (i = 0; i < n_derivs; ++i) {
                    if (! derivs[i])
                         continue;

                    phase_matrix_symmetry2(n_quad, n_stokes, T_m_l[i].a, R_m_l[i].a, T_p_l[i].a, R_p_l[i].a, 1.);
               }
}
          }
     }
}



/*******************************************************************************
 *
 ******************************************************************************/
void rtm_pade_rts2(int n_quad, int n_stokes, int n_derivs, double F_0,
                  double *qxv, double *qwv, double *umus, int n_umus,
                  double planck0, double planck1, double *planck0_l, double *planck1_l,
                  double omega, double *omega_l, double ltau, double *ltau_l,
                  double as_0, double *as_0_l, double atran, double *atran_l,
                  double  *P_x0_p, double  *P_x0_m,
                  double  **r_p, double  **t_p, double  **r_m, double  **t_m,
                  double  **R_p, double  **T_p, double  **R_m, double  **T_m,
                  double  *S_p, double  *S_m, double  *Sl_p, double  *Sl_m,
                  double **P_x0_p_l, double **P_x0_m_l,
                  double ***r_p_l, double ***t_p_l, double ***r_m_l, double ***t_m_l,
                  double ***R_p_l, double ***T_p_l, double ***R_m_l, double ***T_m_l,
                  double **S_p_l, double **S_m_l, double **Sl_p_l, double **Sl_m_l,
                  int pade_s, int pade_r, int check_condition, int symmetric,
                  int solar, int thermal, int vector, int gamma_init, double *condition,
                  uchar *derivs_layers, uchar *derivs_beam, uchar *derivs_thermal,
                  uchar *derivs_sources, save_tree_data save_tree, work_data work) {

     int i;
     int ii;
     int j;
     int jj;
     int k;
     int kk;

     int n_quad_v;

     int n_quad2;
     int n_umus2;
     int n_umus_v2;
     int n_quad_v2;

     int m_pade;

     int symmetric2;

     double a;
     double b;

     double pade_c[MAX_PADE_R];

     double *F_p;
     double *F_m;

     double *Ft0_p;
     double *Ft0_m;
     double *Ft1_p;
     double *Ft1_m;

     double **F_p_l;
     double **F_m_l;

     double **Ft0_p_l;
     double **Ft0_m_l;
     double **Ft1_p_l;
     double **Ft1_m_l;

     double **tpr;
     double **tmr;
     double **gamma;

     double ***tpr_l;
     double ***tmr_l;
     double ***gamma_l;

     matd1d2 A12;
     matd1d2 A21;

     matd1d2 *B12;
     matd1d2 *B21;

     matd1d2 N11;
     matd1d2 N12;
     matd1d2 N21;
     matd1d2 N22;

     matd1d2 D11;
     matd1d2 D12;
     matd1d2 D21;
     matd1d2 D22;

     matd1d2 *O11;
     matd1d2 *O12;
     matd1d2 *O21;
     matd1d2 *O22;

     matd1d2 *E11;
     matd1d2 *E12;
     matd1d2 *E21;
     matd1d2 *E22;
/*
     matd1d2 R_p2;
     matd1d2 T_p2;
*/
     matd1d2 R_m2;
     matd1d2 T_m2;
/*
     matd1d2 *R_p_l2;
     matd1d2 *T_p_l2;
*/
     matd1d2 *R_m_l2;
     matd1d2 *T_m_l2;


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     n_quad_v = n_quad * n_stokes;


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (pade_s < 0 || pade_r <= 0)
          pade_get_s_and_r(r_p, t_p, n_quad_v, ltau, n_umus, umus, &pade_s, &pade_r, flags_or(derivs_layers, n_derivs));

     pade_coefs_x(pade_r, pade_r, 1., pade_c);


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     n_umus2   = n_umus;
     n_quad2   = n_quad  - n_umus;

     n_umus_v2 = n_umus2 * n_stokes;
     n_quad_v2 = n_quad2 * n_stokes;


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (solar) {
          F_p = get_work1(&work, WORK_DX);
          F_m = get_work1(&work, WORK_DX);
     }

     if (thermal) {
          Ft0_p = get_work1(&work, WORK_DX);
          Ft0_m = get_work1(&work, WORK_DX);
          Ft1_p = get_work1(&work, WORK_DX);
          Ft1_m = get_work1(&work, WORK_DX);
     }

     tpr = get_work1(&work, WORK_DXX);
     tmr = get_work1(&work, WORK_DXX);

     A12 = get_work_d2_d(&work, n_quad_v2, n_quad_v2, n_umus_v2);
     A21 = get_work_d2_d(&work, n_quad_v2, n_quad_v2, n_umus_v2);

     N11 = get_work_d2_d(&work, n_quad_v2, n_quad_v2, n_umus_v2);
     N12 = get_work_d2_d(&work, n_quad_v2, n_quad_v2, n_umus_v2);
     N21 = get_work_d2_d(&work, n_quad_v2, n_quad_v2, n_umus_v2);
     N22 = get_work_d2_d(&work, n_quad_v2, n_quad_v2, n_umus_v2);

     D11 = get_work_d2_d(&work, n_quad_v2, n_quad_v2, n_umus_v2);
     D12 = get_work_d2_d(&work, n_quad_v2, n_quad_v2, n_umus_v2);
     D21 = get_work_d2_d(&work, n_quad_v2, n_quad_v2, n_umus_v2);
     D22 = get_work_d2_d(&work, n_quad_v2, n_quad_v2, n_umus_v2);

     R_m2 = get_work_d2_d(&work, n_quad_v2, n_quad_v2, n_umus_v2);
     T_m2 = get_work_d2_d(&work, n_quad_v2, n_quad_v2, n_umus_v2);

     if (solar) {
          if (flags_or(derivs_beam, n_derivs)) {
               F_p_l = get_work2(&work, WORK_DX, WORK_DERIVS_V, derivs_beam);
               F_m_l = get_work2(&work, WORK_DX, WORK_DERIVS_V, derivs_beam);
          }
     }

     if (thermal) {
          if (flags_or(derivs_thermal, n_derivs)) {
               Ft0_p_l = get_work2(&work, WORK_DX, WORK_DERIVS_V, derivs_thermal);
               Ft0_m_l = get_work2(&work, WORK_DX, WORK_DERIVS_V, derivs_thermal);
               Ft1_p_l = get_work2(&work, WORK_DX, WORK_DERIVS_V, derivs_thermal);
               Ft1_m_l = get_work2(&work, WORK_DX, WORK_DERIVS_V, derivs_thermal);
          }
     }

     if (flags_or(derivs_layers, n_derivs)) {
          tpr_l = get_work2(&work, WORK_DXX, WORK_DERIVS_V, derivs_layers);
          tmr_l = get_work2(&work, WORK_DXX, WORK_DERIVS_V, derivs_layers);

          B12 = get_work_d2_d2(&work, n_derivs, n_quad_v2, n_quad_v2, n_umus_v2);
          B21 = get_work_d2_d2(&work, n_derivs, n_quad_v2, n_quad_v2, n_umus_v2);

          O11 = get_work_d2_d2(&work, n_derivs, n_quad_v2, n_quad_v2, n_umus_v2);
          O12 = get_work_d2_d2(&work, n_derivs, n_quad_v2, n_quad_v2, n_umus_v2);
          O21 = get_work_d2_d2(&work, n_derivs, n_quad_v2, n_quad_v2, n_umus_v2);
          O22 = get_work_d2_d2(&work, n_derivs, n_quad_v2, n_quad_v2, n_umus_v2);

          E11 = get_work_d2_d2(&work, n_derivs, n_quad_v2, n_quad_v2, n_umus_v2);
          E12 = get_work_d2_d2(&work, n_derivs, n_quad_v2, n_quad_v2, n_umus_v2);
          E21 = get_work_d2_d2(&work, n_derivs, n_quad_v2, n_quad_v2, n_umus_v2);
          E22 = get_work_d2_d2(&work, n_derivs, n_quad_v2, n_quad_v2, n_umus_v2);

          R_m_l2 = get_work_d2_d2(&work, n_derivs, n_quad_v2, n_quad_v2, n_umus_v2);
          T_m_l2 = get_work_d2_d2(&work, n_derivs, n_quad_v2, n_quad_v2, n_umus_v2);
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     build_txr(n_quad, n_stokes, n_derivs, r_p, t_p, tpr, tmr, r_p_l, t_p_l, tpr_l, tmr_l, derivs_layers, work);


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     m_pade = 1 << pade_s;

     a = ltau / m_pade;

     for (i = 0; i < n_quad_v2; ++i) {
          for (j = 0; j < n_quad_v2; ++j) {
               A12.a[i][j] = tpr[i][j] * a;
               A21.a[i][j] = tmr[i][j] * a;
          }
     }

     if (n_umus2 > 0) {
          ii = 0;
          for (i = n_quad_v2; i < n_quad_v; ++i) {
               jj = 0;
               for (j = 0; j < n_quad_v2; ++j) {
                    A12.d1[ii][jj] = tpr[i][j] * a;
                    A21.d1[ii][jj] = tmr[i][j] * a;
                    jj++;
               }

               A12.d2[ii] = tpr[i][i] * a;
               A21.d2[ii] = tmr[i][i] * a;

               ii++;
          }
     }

     if (flags_or(derivs_layers, n_derivs)) {
          for (i = 0; i < n_derivs; ++i) {
               if (! derivs_layers[i])
                    continue;

               b = ltau_l[i] / m_pade;

               for (j = 0; j < n_quad_v2; ++j) {
                    for (k = 0; k < n_quad_v2; ++k) {
                         B12[i].a[j][k] = tpr_l[i][j][k] * a + tpr[j][k] * b;
                         B21[i].a[j][k] = tmr_l[i][j][k] * a + tmr[j][k] * b;
                    }
               }

               if (n_umus2 > 0) {
                    jj = 0;
                    for (j = n_quad_v2; j < n_quad_v; ++j) {
                         kk = 0;
                         for (k = 0; k < n_quad_v2; ++k) {
                              B12[i].d1[jj][kk] = tpr_l[i][j][k] * a + tpr[j][k] * b;
                              B21[i].d1[jj][kk] = tmr_l[i][j][k] * a + tmr[j][k] * b;
                              kk++;
                         }

                         B12[i].d2[jj] = tpr_l[i][j][j] * a + tpr[j][j] * b;
                         B21[i].d2[jj] = tmr_l[i][j][j] * a + tmr[j][j] * b;

                         jj++;
                    }
               }
          }
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (! symmetric)
          calc_n_d_gen_3(n_quad2, n_umus2, n_stokes, n_derivs,
                          pade_c, pade_r, A12, A21, B12, B21,
                          N11, N12, N21, N22, D11, D12, D21, D22,
                          O11, O12, O21, O22, E11, E12, E21, E22,
                          derivs_layers, work);
     else
          calc_n_d_sym_3(n_quad2, n_umus2, n_stokes, n_derivs,
                          pade_c, pade_r, A12, A21, B12, B21,
                          N11, N12, N21, N22, D11, D12, D21, D22,
                          O11, O12, O21, O22, E11, E12, E21, E22,
                          derivs_layers, work);


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     reduced_recover(N11, N12, N21, N22, D11, D12, D21, D22,
                     N11, N12, N21, N22, D11, D12, D21, D22, n_quad_v2, n_umus_v2);

     if (flags_or(derivs_layers, n_derivs)) {
          for (i = 0; i < n_derivs; ++i) {
               if (! derivs_layers[i])
                    continue;

               reduced_recover(O11[i], O12[i], O21[i], O22[i], E11[i], E12[i], E21[i], E22[i],
                               O11[i], O12[i], O21[i], O22[i], E11[i], E12[i], E21[i], E22[i], n_quad_v2, n_umus_v2);
          }
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (! vector && n_umus == 0)
          symmetric2 = symmetric;
     else
          symmetric2 = 0;


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
/*
     pade_calc_R_and_T_1(N11, N12, N21, N22, D11, D12, D21, D22,
                         R_p2, T_p2, R_m2, T_m2,
                         O11, O12, O21, O22, E11, E12, E21, E22,
                         R_p_l2, T_p_l2, R_m_l2, T_m_l2,
                         n_quad, n_stokes, n_derivs,
                         check_condition, 1, symmetric2, 0, condition, derivs_layers, work);
*/
     pade_calc_R_and_T_2(N11, N12, N21, N22, D11, D12, D21, D22,
                         NONE, NONE, R_m2, T_m2,
                         O11, O12, O21, O22, E11, E12, E21, E22,
                         NULL, NULL, R_m_l2, T_m_l2,
                         n_quad2, n_stokes, n_derivs,
                         check_condition, 1, symmetric2, 0, condition, derivs_layers, work, n_umus2);

     matd1d2_to_matrix(R_m2, R_p, n_quad_v2, n_umus_v2);
     matd1d2_to_matrix(T_m2, T_p, n_quad_v2, n_umus_v2);

     if (flags_or(derivs_layers, n_derivs)) {
          for (i = 0; i < n_derivs; ++i) {
               if (! derivs_layers[i])
                    continue;

               matd1d2_to_matrix(R_m_l2[i], R_p_l[i], n_quad_v2, n_umus_v2);
               matd1d2_to_matrix(T_m_l2[i], T_p_l[i], n_quad_v2, n_umus_v2);
          }
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     for (i = 0; i < pade_s; ++i) {
          if (! symmetric2) {
               if (! flags_or(derivs_layers, n_derivs))
                    layer_double  (R_p, T_p, NULL, NULL, NULL, NULL, NULL, NULL, n_quad_v, 0., 0., 0, 0, i == 0, work);
               else
                    layer_double_l(R_p, T_p, NULL, NULL, NULL, NULL, NULL, NULL, R_p_l, T_p_l, NULL, NULL, NULL, NULL, NULL, NULL, n_quad_v, n_derivs, 0., NULL, 0., NULL, 0, 0, i == 0, derivs_layers, NULL, NULL, work);
          }
          else {
               if (! flags_or(derivs_layers, n_derivs))
                    layer_double_s  (R_p, T_p, NULL, NULL, NULL, NULL, NULL, NULL, n_quad_v, 0., 0., 0, 0, i == 0, work);
               else
                    layer_double_s_l(R_p, T_p, NULL, NULL, NULL, NULL, NULL, NULL, R_p_l, T_p_l, NULL, NULL, NULL, NULL, NULL, NULL, n_quad_v, n_derivs, 0., NULL, 0., NULL, 0, 0, i == 0, derivs_layers, NULL, NULL, work);
          }
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (vector) {
          dmat_mul_D_A(n_quad, n_stokes, R_p, R_p);
          phase_matrix_symmetry2(n_quad, n_stokes, T_p, R_p, T_m, R_m, 1.);

          if (flags_or(derivs_layers, n_derivs)) {
               for (i = 0; i < n_derivs; ++i) {
                    if (! derivs_layers[i])
                    continue;

                    dmat_mul_D_A(n_quad, n_stokes, R_p_l[i], R_p_l[i]);
                    phase_matrix_symmetry2(n_quad, n_stokes, T_p_l[i], R_p_l[i], T_m_l[i], R_m_l[i], 1.);
               }
          }
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (solar) {
          if (! CLASSICAL_PARTICULAR_SOLUTION_USE_2D) {
               gamma = get_work1(&work, WORK_DXX);

               if (flags_or(derivs_layers, n_derivs))
                    gamma_l = get_work2(&work, WORK_DXX, WORK_DERIVS_V, derivs_layers);

               if (gamma_init)
                    dmat_init(gamma, n_quad_v, n_quad_v, 1.);
               else
                    build_gamma(n_quad_v, n_derivs, tpr, tmr, gamma, tpr_l, tmr_l, gamma_l, derivs_layers, work);

               build_source_vectors_solar_classic_1n(n_quad, n_stokes, n_derivs, qxv, F_0, omega, omega_l, as_0, as_0_l, P_x0_p, P_x0_m, tpr, tmr, gamma, F_p, F_m, P_x0_p_l, P_x0_m_l, tpr_l, tmr_l, gamma_l, F_p_l, F_m_l, derivs_layers, derivs_beam, save_tree, work);
          }
          else {
               if (! vector)
                    build_source_vectors_solar_classic_2n(n_quad_v, n_derivs, qxv, F_0, omega, omega_l, as_0, as_0_l, P_x0_p, P_x0_m, r_p, t_p, F_p, F_m, P_x0_p_l, P_x0_m_l, r_p_l, t_p_l, F_p_l, F_m_l, derivs_layers, derivs_beam, work);
               else
                    build_source_vectors_solar_classic_2n2(n_quad_v, n_derivs, qxv, F_0, omega, omega_l, as_0, as_0_l, P_x0_p, P_x0_m, r_p, t_p, r_m, t_m, F_p, F_m, P_x0_p_l, P_x0_m_l, r_p_l, t_p_l, r_m_l, t_m_l, F_p_l, F_m_l, derivs_layers, derivs_beam, work);
          }

          if (! vector)
               build_global_source_solar (n_quad_v, n_derivs, atran, atran_l, R_p, T_p, F_p, F_m, S_p, S_m, R_p_l, T_p_l, F_p_l, F_m_l, S_p_l, S_m_l, derivs_layers, derivs_beam, work, NULL, NULL, NULL, NULL);
          else
               build_global_source_solar2(n_quad_v, n_derivs, atran, atran_l, R_p, T_p, R_m, T_m, F_p, F_m, S_p, S_m, R_p_l, T_p_l, R_m_l, T_m_l, F_p_l, F_m_l, S_p_l, S_m_l, derivs_layers, derivs_beam, work, NULL, NULL, NULL, NULL);
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (thermal) {
          build_source_vectors_thermal2(n_quad, n_stokes, n_derivs, qxv, planck0, planck1, planck0_l, planck1_l, omega, omega_l, ltau, ltau_l, r_p, t_p, r_m, t_m, Ft0_p, Ft0_m, Ft1_p, Ft1_m, r_p_l, t_p_l, r_m_l, t_m_l, Ft0_p_l, Ft0_m_l, Ft1_p_l, Ft1_m_l, derivs_layers, derivs_thermal, work);

          if (! vector)
               build_global_source_thermal(n_quad_v, n_derivs, R_p, T_p, R_p, T_p, Ft0_p, Ft0_m, Ft1_p, Ft1_m, Sl_p, Sl_m, R_p_l, T_p_l, R_p_l, T_p_l, Ft0_p_l, Ft0_m_l, Ft1_p_l, Ft1_m_l, Sl_p_l, Sl_m_l, derivs_layers, derivs_thermal, work);
          else
               build_global_source_thermal(n_quad_v, n_derivs, R_p, T_p, R_m, T_m, Ft0_p, Ft0_m, Ft1_p, Ft1_m, Sl_p, Sl_m, R_p_l, T_p_l, R_m_l, T_m_l, Ft0_p_l, Ft0_m_l, Ft1_p_l, Ft1_m_l, Sl_p_l, Sl_m_l, derivs_layers, derivs_thermal, work);
     }
}



/*******************************************************************************
 *
 ******************************************************************************/
#ifdef INCLUDE_DEV_SOURCE
#include "../src2/xrtm_pade_rts2.c"
#endif
