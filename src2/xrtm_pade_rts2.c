/*******************************************************************************
**
**    Copyright (C) 2007-2020 Greg McGarragh <greg.mcgarragh@colostate.edu>
**
**    This source code is licensed under the GNU General Public License (GPL),
**    Version 3.  See the file COPYING for more details.
**
*******************************************************************************/
/*
#include <mpfr.h>
*/

/*******************************************************************************
 *
 ******************************************************************************/
/*
static double taylor_exp2(int i, int n, double x);

static double taylor_exp(int n, double x) {

     return taylor_exp2(0, n, x);
}



static double taylor_exp2(int i, int n, double x) {

     double a;
     double b;
     double c;

     if (i == 0) {
          a = 1.;
          b = 1.;
          c = 1.;
     }
     else {
          a = pow(x, i);
          b = factorial(i);
          c = a / b;
     }

     for (i = i + 1; i <= n; ++i) {
          a *= x;
          b *= (double) i;
          c += a / b;
     }

     return c;
}
*/


/*******************************************************************************
 *
 ******************************************************************************/
#ifdef CRAP
static double pade_exp(int r, double *coefs, double x) {

     int i;

     double a;
     double f;
     double n;
     double d;

     n = coefs[0];
     d = coefs[0];

     a =  1;
     f = -1;
     for (i = 1; i <= r; ++i) {
          a *= x;
          n += coefs[i] * a;
          d += coefs[i] * a * f;
          f *= -1;
     }

     return n / d;
/*
     int i;

     double a;
     double x2;
     double u;
     double v;

     u = coefs[0];
     v = coefs[1];

     x2 = x * x;

     a = 1;
     for (i = 2; i <= r; i += 2) {
          a *= x2;
          u += coefs[i] * a;
          if (i < r)
               v += coefs[i+1] * a;
     }

     v *= x;

     return (u + v) / (u - v);
*/
}



static double najfeld_exp(int r, double *coefs, double x) {

     int i;

     double a;
     double b;
     double n;
     double d;

     double coefs1[] = {1., 8. / 17., 7. / 255., 4. /  9945., 1. /   765765.};
     double coefs2[] = {1., 7. / 51., 1. / 255., 2. / 69615., 1. / 34459425.};

     x /= 2.;

     n = coefs1[0];
     d = coefs2[0];

     a = b = x * x;

     n += coefs1[1] * a;
     d += coefs2[1] * a;

     for (i = 2; i <= 4; ++i) {
          a *= b;
          n += coefs1[i] * a;
          d += coefs2[i] * a;
     }

     a = x * d;

     return (n + a) / (n - a);
}
#endif

/*
static void pade_coefs_x_mprf(int p, int q, double f, mpfr_t *coefs) {

     int i;
     int ii;

     long ia;
     long ib;

     mpfr_t ma;
     mpfr_t mb;
     mpfr_t mc;
     mpfr_t mf;
     mpfr_t mg;

     mpfr_init(ma);
     mpfr_init(mb);
     mpfr_init(mc);
     mpfr_init(mf);
     mpfr_init(mg);

     mpfr_set_d(coefs[0], 1., GMP_RNDN);

     mpfr_set_d(ma, 1., GMP_RNDN);
     mpfr_set_d(mf, f, GMP_RNDN);
     mpfr_set_d(mg, f, GMP_RNDN);
     for (i = 0, ii = 1; i < p; ++i, ++ii) {
          ia = p - i;
          ib = (2 * p - i) * ii;

          mpfr_set_si(mb, ia, GMP_RNDN);
          mpfr_set_si(mc, ib, GMP_RNDN);
          mpfr_div(mb, mb, mc, GMP_RNDN);
          mpfr_mul(ma, ma, mb, GMP_RNDN);

          mpfr_mul(coefs[ii], ma, mg, GMP_RNDN);

          mpfr_mul(mg, mg, mf, GMP_RNDN);
     }
}



static void pade_exp_mpfr(int r, mpfr_t *coefs, mpfr_t x, mpfr_t e) {

     int i;

     mpfr_t a;
     mpfr_t b;
     mpfr_t x2;
     mpfr_t u;
     mpfr_t v;

     mpfr_init(a);
     mpfr_init(b);
     mpfr_init(x2);
     mpfr_init(u);
     mpfr_init(v);

     mpfr_set(u, coefs[0], GMP_RNDN);
     mpfr_set(v, coefs[1], GMP_RNDN);

     mpfr_mul(x2, x, x, GMP_RNDN);

     mpfr_set_d(a, 1., GMP_RNDN);
     for (i = 2; i <= r; i += 2) {
          mpfr_mul(a, a, x2, GMP_RNDN);

          mpfr_mul(b, coefs[i], a, GMP_RNDN);
          mpfr_add(u, u, b, GMP_RNDN);
          if (i < r) {
               mpfr_mul(b, coefs[i+1], a, GMP_RNDN);
               mpfr_add(v, v, b, GMP_RNDN);
          }
     }

     mpfr_mul(v, v, x, GMP_RNDN);

     mpfr_add(a, u, v, GMP_RNDN);
     mpfr_sub(b, u, v, GMP_RNDN);
     mpfr_div(e, a, b, GMP_RNDN);
}
*/


/*******************************************************************************
 *
 ******************************************************************************/
static void dmat_mul_a2(double **ax1, double **ax2,
                        double **b11, double **b12, double **b21, double **b22,
                        long n, double **cx1, double **cx2) {

     dmat_mul_brow(ax1, ax2, b11, b12, b21, b22, n, cx1, cx2);
}



/*******************************************************************************
 *
 ******************************************************************************/
static void dmat_mul_a3(double **ax1, double **ax2,
                        double **b11, double **b12, double **b21, double **b22,
                        int n, double **cx1, double **cx2, double f, int flag) {

     dsym_mul_kernel(ax1, b11, 1., n, n, cx1, 0., 0,    0);
     dsym_mul_kernel(ax2, b21, 1., n, n, cx1, 1., flag, 0);

     dsym_mul_kernel(ax1, b12, 1., n, n, cx2, 0., 0,    0);
     dsym_mul_kernel(ax2, b22, 1., n, n, cx2, 1., 0,    0);

     if (flag)
          dmat_sym(cx2, n, f);
}



/*******************************************************************************
 *
 ******************************************************************************/
static void r_t_matrix_symmetry(int n_quad, int n_stokes,
          double **A11, double **A12, double **A21, double **A22, double f) {

     int i;
     int j;

     int n;

     n = n_quad * n_stokes;

     if (n_stokes <= 2) {
          for (i = 0; i < n; ++i) {
               for (j = 0; j < n; ++j) {
                    A21[i][j] = f * A12[i][j];
                    A22[i][j] = f * A11[i][j];
               }
          }
     }
     else
          phase_matrix_symmetry2(n_quad, n_stokes, A11, A12, A22, A21, f);
}



/*******************************************************************************
 *
 ******************************************************************************/
static void add_u_v_term_2(int n_quad, int n_stokes,
                           double **a11, double **a12,
                           double **b11, double **b12, double **b21, double **b22,
                           double **c11, double **c12,
                           double **u11, double **u12, double **v11, double **v12,
                           double **w11, double **w12, double **x11, double **x12,
                           double c0, double c1, double c2, double c3, int symmetric) {

     int n_quad_v;

     n_quad_v = n_quad * n_stokes;

     if (! symmetric)
          dmat_mul_a2(a11, a12, b11, b12, b21, b22, n_quad_v, c11, c12);
     else
          dmat_mul_a3(a11, a12, b11, b12, b21, b22, n_quad_v, c11, c12, -1., 1);

     if (c0 != 0.)
          dmat_scale_add_1x2(u11, u12, c11, c12, u11, u12, c0, n_quad_v);
     if (c1 != 0.)
          dmat_scale_add_1x2(v11, v12, c11, c12, v11, v12, c1, n_quad_v);
     if (c2 != 0.)
          dmat_scale_add_1x2(w11, w12, c11, c12, w11, w12, c2, n_quad_v);
     if (c3 != 0.)
          dmat_scale_add_1x2(x11, x12, c11, c12, x11, x12, c3, n_quad_v);
}



static void add_u_v_term_l_2(int n_quad, int n_stokes,
                             double **ax11, double **ax12,
                             double **bx11, double **bx12, double **bx21, double **bx22,
                             double **ay11, double **ay12,
                             double **by11, double **by12, double **by21, double **by22,
                             double **cy11, double **cy12,
                             double **uy11, double **uy12, double **vy11, double **vy12,
                             double **wy11, double **wy12, double **xy11, double **xy12,
                             double **w111, double **w112, double **w121, double **w122,
                             double **w211, double **w212, double **w221, double **w222,
                             double c0, double c1, double c2, double c3, int symmetric) {

     int n_quad_v;

     n_quad_v = n_quad * n_stokes;

     if (! symmetric) {
          dmat_mul_a2(ay11, ay12, bx11, bx12, bx21, bx22, n_quad_v, w111, w112);
          dmat_mul_a2(ax11, ax12, by11, by12, by21, by22, n_quad_v, w211, w212);
          dmat_scale_add_1x2(w111, w112, w211, w212, cy11, cy12, 1., n_quad_v);
     }
     else {
          dmat_mul_a3(ay11, ay12, bx11, bx12, bx21, bx22, n_quad_v, w111, w112, -1., 0);
          dmat_mul_a3(ax11, ax12, by11, by12, by21, by22, n_quad_v, w211, w212, -1., 0);
          dsym_add(w111, w211, cy11, n_quad_v, 1);
          dsym_add(w112, w212, cy12, n_quad_v, 0);
          dmat_sym(cy12, n_quad_v, -1.);
     }

     if (c0 != 0.)
          dmat_scale_add_1x2(uy11, uy12, cy11, cy12, uy11, uy12, c0, n_quad_v);
     if (c1 != 0.)
          dmat_scale_add_1x2(vy11, vy12, cy11, cy12, vy11, vy12, c1, n_quad_v);
     if (c2 != 0.)
          dmat_scale_add_1x2(wy11, wy12, cy11, cy12, wy11, wy12, c2, n_quad_v);
     if (c3 != 0.)
          dmat_scale_add_1x2(xy11, xy12, cy11, cy12, xy11, xy12, c3, n_quad_v);
}



/*******************************************************************************
 *
 ******************************************************************************/
static void calc_n_d_2(int n_quad, int n_stokes, int n_derivs, double *pade_c, int pade_r,
                       double  **A11, double  **A12, double  **A21, double  **A22,
                       double ***B11, double ***B12, double ***B21, double ***B22,
                       double  **N11, double  **N12, double  **N21, double  **N22,
                       double  **D11, double  **D12, double  **D21, double  **D22,
                       double ***O11, double ***O12, double ***O21, double ***O22,
                       double ***E11, double ***E12, double ***E21, double ***E22,
                       int symmetric, uchar *derivs, work_data work) {

     int i;
     int j;
     int k;

     int n_quad_v;

     double **w111;
     double **w112;
     double **w121;
     double **w122;
     double **w311;
     double **w312;
     double **w321;
     double **w322;
     double **w411;
     double **w412;
     double **w421;
     double **w422;

     double **wx11[2];
     double **wx12[2];

     double ***w211;
     double ***w212;
     double ***w221;
     double ***w222;

     double ***wy11[2];
     double ***wy12[2];


     n_quad_v = n_quad * n_stokes;


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (pade_r > 1) {
          w111 = get_work1(&work, WORK_DXX);
          w112 = get_work1(&work, WORK_DXX);
          w121 = get_work1(&work, WORK_DXX);
          w122 = get_work1(&work, WORK_DXX);
     }

     if (pade_r > 2) {
          wx11[0] = get_work1(&work, WORK_DXX);
          wx12[0] = get_work1(&work, WORK_DXX);
          wx11[1] = get_work1(&work, WORK_DXX);
          wx12[1] = get_work1(&work, WORK_DXX);
     }

     if (flags_or(derivs, n_derivs)) {
          if (pade_r > 1) {
               w311 = get_work1(&work, WORK_DXX);
               w312 = get_work1(&work, WORK_DXX);
               w321 = get_work1(&work, WORK_DXX);
               w322 = get_work1(&work, WORK_DXX);
               w411 = get_work1(&work, WORK_DXX);
               w412 = get_work1(&work, WORK_DXX);
               w421 = get_work1(&work, WORK_DXX);
               w422 = get_work1(&work, WORK_DXX);
          }

          if (pade_r > 2) {
               w211    = get_work2(&work, WORK_DXX, WORK_DERIVS_V, derivs);
               w212    = get_work2(&work, WORK_DXX, WORK_DERIVS_V, derivs);
               w221    = get_work2(&work, WORK_DXX, WORK_DERIVS_V, derivs);
               w222    = get_work2(&work, WORK_DXX, WORK_DERIVS_V, derivs);

               wy11[0] = get_work2(&work, WORK_DXX, WORK_DERIVS_V, derivs);
               wy12[0] = get_work2(&work, WORK_DXX, WORK_DERIVS_V, derivs);
               wy11[1] = get_work2(&work, WORK_DXX, WORK_DERIVS_V, derivs);
               wy12[1] = get_work2(&work, WORK_DXX, WORK_DERIVS_V, derivs);
          }
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (pade_r == 1 ||
         pade_r == 2) {
          dmat_init_2x2(N11, N12, N21, N22, pade_c[0], n_quad_v);

          dmat_scale_add_2x2(N11, N12, N21, N22,
                            A11, A12, A21, A22,
                            D11, D12, D21, D22, -pade_c[1], n_quad_v);
          dmat_scale_add_2x2(N11, N12, N21, N22,
                            A11, A12, A21, A22,
                            N11, N12, N21, N22,  pade_c[1], n_quad_v);

          for (i = 0; i < n_derivs; ++i) {
               if (! derivs[i])
                    continue;

               dmat_init_2x2(O11[i], O12[i], O21[i], O22[i], 0., n_quad_v);

               dmat_scale_add_2x2(O11[i], O12[i], O21[i], O22[i],
                                 B11[i], B12[i], B21[i], B22[i],
                                 E11[i], E12[i], E21[i], E22[i], -pade_c[1], n_quad_v);
               dmat_scale_add_2x2(O11[i], O12[i], O21[i], O22[i],
                                 B11[i], B12[i], B21[i], B22[i],
                                 O11[i], O12[i], O21[i], O22[i],  pade_c[1], n_quad_v);
          }
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (pade_r == 2) {
          if (! symmetric)
               dmat_mul_a2(A11, A12, A11, A12,
                          A21, A22, n_quad_v, w111, w112);
          else
               dmat_mul_a3(A11, A12, A11, A12,
                          A21, A22, n_quad_v, w111, w112, -1., 1);
          r_t_matrix_symmetry(n_quad, n_stokes,
                              w111, w112, w121, w122, 1.);

          dmat_scale_add_2x2(D11, D12, D21, D22,
                            w111, w112, w121, w122,
                            D11, D12, D21, D22, pade_c[2], n_quad_v);
          dmat_scale_add_2x2(N11, N12, N21, N22,
                            w111, w112, w121, w122,
                            N11, N12, N21, N22, pade_c[2], n_quad_v);

          for (i = 0; i < n_derivs; ++i) {
               if (! derivs[i])
                    continue;

               dmat_mul_a2(B11[i], B12[i], A11, A12,
                          A21, A22, n_quad_v, w311, w312);

               if (! symmetric)
                    dmat_mul_a2(A11, A12, B11[i], B12[i],
                               B21[i], B22[i], n_quad_v, w411, w412);
               else {
                    dmat_scale_trans(w311, w411,  1., n_quad_v, n_quad_v);
                    dmat_scale_trans(w312, w412, -1., n_quad_v, n_quad_v);
               }

               dmat_scale_add_1x2(w311, w312, w411, w412,
                                 w311, w312, 1., n_quad_v);

               r_t_matrix_symmetry(n_quad, n_stokes,
                                   w311, w312, w321, w322, 1.);

               dmat_scale_add_2x2(E11[i], E12[i], E21[i], E22[i],
                                 w311, w312, w321, w322,
                                 E11[i], E12[i], E21[i], E22[i], pade_c[2], n_quad_v);
               dmat_scale_add_2x2(O11[i], O12[i], O21[i], O22[i],
                                 w311, w312, w321, w322,
                                 O11[i], O12[i], O21[i], O22[i], pade_c[2], n_quad_v);
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
          dmat_init_1x2(N11, N12, pade_c[0], n_quad_v);
          dmat_init_1x2(D11, D12, pade_c[1], n_quad_v);
          dmat_init_1x2(N21, N22, pade_c[6], n_quad_v);
          dmat_init_1x2(D21, D22, pade_c[7], n_quad_v);

          add_u_v_term_2(n_quad, n_stokes, A11, A12, A11, A12, A21, A22, w111, w112, N11, N12, D11, D12, N21, N22, D21, D22, pade_c[2], 2 < pade_r ? pade_c[3] : 0., pade_c[8], 8 < pade_r ? pade_c[9] : 0., symmetric);

          r_t_matrix_symmetry(n_quad, n_stokes,
                              w111, w112, w121, w122, 1.);

          dmat_copy_1x2(wx11[0], wx12[0], w111, w112, n_quad_v);


          for (i = 0; i < n_derivs; ++i) {
               if (! derivs[i])
                    continue;

               dmat_init_1x2(O11[i], O12[i], 0., n_quad_v);
               dmat_init_1x2(E11[i], E12[i], 0., n_quad_v);
               dmat_init_1x2(O21[i], O22[i], 0., n_quad_v);
               dmat_init_1x2(E21[i], E22[i], 0., n_quad_v);

               add_u_v_term_l_2(n_quad, n_stokes, A11, A12, A11, A12, A21, A22, B11[i], B12[i], B11[i], B12[i], B21[i], B22[i], w211[i], w212[i], O11[i], O12[i], E11[i], E12[i], O21[i], O22[i], E21[i], E22[i], w311, w312, w321, w322, w411, w412, w421, w422, pade_c[2], 2 < pade_r ? pade_c[3] : 0., pade_c[8], 8 < pade_r ? pade_c[9] : 0., symmetric);

               r_t_matrix_symmetry(n_quad, n_stokes,
                                   w211[i], w212[i], w221[i], w222[i], 1.);

               dmat_copy_1x2(wy11[0][i], wy12[0][i], w211[i], w212[i], n_quad_v);
          }


          i = 4;
          j = 1;
          add_u_v_term_2(n_quad, n_stokes, wx11[1 - j], wx12[1 - j], w111, w112, w121, w122, wx11[j], wx12[j], N11,  N12,  D11,  D12,  N21, N22, D21, D22, pade_c[i], i < pade_r ? pade_c[i+1] : 0., pade_c[i+6], i+6 < pade_r ? pade_c[i+7] : 0., symmetric);
          for (k = 0; k < n_derivs; ++k) {
               if (! derivs[k])
                    continue;

               add_u_v_term_l_2(n_quad, n_stokes, wx11[1 - j], wx12[1 - j], w111, w112, w121, w122, wy11[1 - j][k], wy12[1 - j][k], w211[k], w212[k], w221[k], w222[k], wy11[j][k], wy12[j][k], O11[k], O12[k], E11[k], E12[k], O21[k], O22[k], E21[k], E22[k], w311, w312, w321, w322, w411, w412, w421, w422, pade_c[i], i < pade_r ? pade_c[i+1] : 0., pade_c[i+6], i+6 < pade_r ? pade_c[i+7] : 0., symmetric);
          }
          i = i + 2;
          j = 1 - j;

          add_u_v_term_2(n_quad, n_stokes, wx11[1 - j], wx12[1 - j], w111, w112, w121, w122, wx11[j], wx12[j], NULL, NULL, NULL, NULL, N21, N22, D21, D22, 0.,       0.,                             pade_c[i+6], i+6 < pade_r ? pade_c[i+7] : 0., symmetric);
          for (k = 0; k < n_derivs; ++k) {
               if (! derivs[k])
                    continue;

               add_u_v_term_l_2(n_quad, n_stokes, wx11[1 - j], wx12[1 - j], w111, w112, w121, w122, wy11[1 - j][k], wy12[1 - j][k], w211[k], w212[k], w221[k], w222[k], wy11[j][k], wy12[j][k], NULL, NULL, NULL, NULL, O21[k], O22[k], E21[k], E22[k], w311, w312, w321, w322, w411, w412, w421, w422, 0.,       0.,                           pade_c[i+6], i+6 < pade_r ? pade_c[i+7] : 0., symmetric);
          }

          r_t_matrix_symmetry(n_quad, n_stokes,
                              N21, N22, w121, w122, 1.);
          if (! symmetric)
               dmat_mul_a2(wx11[j], wx12[j], N21, N22, w121, w122, n_quad_v, w111, w112);
          else
               dmat_mul_a3(wx11[j], wx12[j], N21, N22, w121, w122, n_quad_v, w111, w112, -1., 1);
          dmat_scale_add_1x2(w111, w112, N11, N12, N11, N12, 1., n_quad_v);


          for (k = 0; k < n_derivs; ++k) {
               if (! derivs[k])
                    continue;

               r_t_matrix_symmetry(n_quad, n_stokes, O21[k], O22[k], w221[k], w222[k], 1.);

               if (! symmetric) {
                    dmat_mul_a2(wy11[j][k], wy12[j][k], N21, N22, w121, w122, n_quad_v, w111, w112);
                    dmat_mul_a2(wx11[j], wx12[j], O21[k], O22[k], w221[k], w222[k], n_quad_v, w211[k], w212[k]);
                    dmat_scale_add_1x2(w111, w112, w211[k], w212[k], w111, w112, 1., n_quad_v);
               }

               else {
                    dmat_mul_a3(wy11[j][k], wy12[j][k], N21, N22, w121, w122, n_quad_v, w111, w112, -1., 0);
                    dmat_mul_a3(wx11[j], wx12[j], O21[k], O22[k], w221[k], w222[k], n_quad_v, w211[k], w212[k], -1., 0);
                    dsym_add(w111, w211[k], w111, n_quad_v, 1);
                    dsym_add(w112, w212[k], w112, n_quad_v, 0);
                    dmat_sym(w112, n_quad_v, -1.);
               }

               dmat_scale_add_1x2(w111, w112, O11[k], O12[k], O11[k], O12[k], 1., n_quad_v);
          }


          r_t_matrix_symmetry(n_quad, n_stokes,
                              D21, D22, w121, w122, 1.);
          if (! symmetric)
               dmat_mul_a2(wx11[j], wx12[j], D21, D22, w121, w122, n_quad_v, w111, w112);
          else
               dmat_mul_a3(wx11[j], wx12[j], D21, D22, w121, w122, n_quad_v, w111, w112, -1., 1);
          dmat_scale_add_1x2(w111, w112, D11, D12, D11, D12, 1., n_quad_v);


          for (k = 0; k < n_derivs; ++k) {
               if (! derivs[k])
                    continue;

               r_t_matrix_symmetry(n_quad, n_stokes, E21[k], E22[k], w221[k], w222[k], 1.);

               if (! symmetric) {
                    dmat_mul_a2(wy11[j][k], wy12[j][k], D21, D22, w121, w122, n_quad_v, w111, w112);
                    dmat_mul_a2(wx11[j], wx12[j], E21[k], E22[k], w221[k], w222[k], n_quad_v, w211[k], w212[k]);
                    dmat_scale_add_1x2(w111, w112, w211[k], w212[k], w111, w112, 1., n_quad_v);
               }

               else {
                    dmat_mul_a3(wy11[j][k], wy12[j][k], D21, D22, w121, w122, n_quad_v, w111, w112, -1., 0);
                    dmat_mul_a3(wx11[j], wx12[j], E21[k], E22[k], w221[k], w222[k], n_quad_v, w211[k], w212[k], -1., 0);
                    dsym_add(w111, w211[k], w111, n_quad_v, 1);
                    dsym_add(w112, w212[k], w112, n_quad_v, 0);
                    dmat_sym(w112, n_quad_v, -1.);
               }

               dmat_scale_add_1x2(w111, w112, E11[k], E12[k], E11[k], E12[k], 1., n_quad_v);
          }
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     else {
          dmat_init_1x2(N11, N12, pade_c[0], n_quad_v);
          dmat_init_1x2(D11, D12, pade_c[1], n_quad_v);

          add_u_v_term_2(n_quad, n_stokes, A11, A12, A11, A12, A21, A22, w111, w112, N11, N12, D11, D12, NULL, NULL, NULL, NULL, pade_c[2], 2 < pade_r ? pade_c[3] : 0., 0., 0., symmetric);

          r_t_matrix_symmetry(n_quad, n_stokes,
                              w111, w112, w121, w122, 1.);

          dmat_copy_1x2(wx11[0], wx12[0], w111, w112, n_quad_v);

          for (i = 0; i < n_derivs; ++i) {
               if (! derivs[i])
                    continue;

               dmat_init_1x2(O11[i], O12[i], 0., n_quad_v);
               dmat_init_1x2(E11[i], E12[i], 0., n_quad_v);

               add_u_v_term_l_2(n_quad, n_stokes, A11, A12, A11, A12, A21, A22, B11[i], B12[i], B11[i], B12[i], B21[i], B22[i], w211[i], w212[i], O11[i], O12[i], E11[i], E12[i], NULL, NULL, NULL, NULL, w311, w312, w321, w322, w411, w412, w421, w422, pade_c[2], 2 < pade_r ? pade_c[3] : 0., 0., 0., symmetric);

               r_t_matrix_symmetry(n_quad, n_stokes,
                                   w211[i], w212[i], w221[i], w222[i], 1.);

               dmat_copy_1x2(wy11[0][i], wy12[0][i], w211[i], w212[i], n_quad_v);
          }

          j = 1;
          for (i = 4; i <= pade_r; i += 2) {
               add_u_v_term_2(n_quad, n_stokes, wx11[1 - j], wx12[1 - j], w111, w112, w121, w122, wx11[j], wx12[j], N11, N12, D11, D12, NULL, NULL, NULL, NULL, pade_c[i], i < pade_r ? pade_c[i+1] : 0., 0., 0., symmetric);

               for (k = 0; k < n_derivs; ++k) {
                    if (! derivs[k])
                         continue;

                    add_u_v_term_l_2(n_quad, n_stokes, wx11[1 - j], wx12[1 - j], w111, w112, w121, w122, wy11[1 - j][k], wy12[1 - j][k], w211[k], w212[k], w221[k], w222[k], wy11[j][k], wy12[j][k], O11[k], O12[k], E11[k], E12[k], NULL, NULL, NULL, NULL, w311, w312, w321, w322, w411, w412, w421, w422, pade_c[i], i < pade_r ? pade_c[i+1] : 0., 0., 0., symmetric);
               }

               j = 1 - j;
          }
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     r_t_matrix_symmetry(n_quad, n_stokes,
                         N11, N12, N21, N22, 1.);
     r_t_matrix_symmetry(n_quad, n_stokes,
                         D11, D12, D21, D22, 1.);

     if (! symmetric)
          dmat_mul_a2(A11, A12, D11, D12,
                     D21, D22, n_quad_v, w111, w112);
     else
          dmat_mul_a3(A11, A12, D11, D12,
                     D21, D22, n_quad_v, w111, w112, 1., 1);
     r_t_matrix_symmetry(n_quad, n_stokes,
                         w111, w112, w121, w122, -1.);

     for (i = 0; i < n_derivs; ++i) {
          if (! derivs[i])
               continue;

          r_t_matrix_symmetry(n_quad, n_stokes,
                              O11[i],  O12[i],  O21[i],  O22[i], 1.);
          r_t_matrix_symmetry(n_quad, n_stokes,
                              E11[i],  E12[i],  E21[i],  E22[i], 1.);

          if (! symmetric) {

               dmat_mul_a2(B11[i], B12[i], D11, D12,
                          D21, D22, n_quad_v, w311, w312);
               dmat_mul_a2(A11, A12, E11[i], E12[i],
                          E21[i], E22[i], n_quad_v, w411, w412);
               dmat_scale_add_1x2(w311, w312, w411, w412,
                                 w211[i], w212[i], 1., n_quad_v);

          }
          else {
               dmat_mul_a3(B11[i], B12[i], D11, D12,
                          D21, D22, n_quad_v, w311, w312, 1., 0);
               dmat_mul_a3(A11, A12, E11[i], E12[i],
                          E21[i], E22[i], n_quad_v, w411, w412, 1., 0);
               dsym_add(w311, w411, w211[i], n_quad_v, 1);
               dsym_add(w312, w412, w212[i], n_quad_v, 0);
               dmat_sym(w212[i], n_quad_v, 1.);
          }

          r_t_matrix_symmetry(n_quad, n_stokes,
                              w211[i], w212[i], w221[i], w222[i], -1.);
     }

     dmat_scale_add_2x2(N11, N12, N21, N22,
                       w111, w112, w121, w122,
                       D11, D12, D21, D22, -1., n_quad_v);
     dmat_scale_add_2x2(N11, N12, N21, N22,
                       w111, w112, w121, w122,
                       N11, N12, N21, N22,  1., n_quad_v);

     for (i = 0; i < n_derivs; ++i) {
          if (! derivs[i])
               continue;

          dmat_scale_add_2x2(O11[i], O12[i], O21[i], O22[i],
                            w211[i], w212[i], w221[i], w222[i],
                            E11[i], E12[i], E21[i], E22[i], -1., n_quad_v);
          dmat_scale_add_2x2(O11[i], O12[i], O21[i], O22[i],
                            w211[i], w212[i], w221[i], w222[i],
                            O11[i], O12[i], O21[i], O22[i],  1., n_quad_v);
     }
}



/*******************************************************************************
 *
 ******************************************************************************/
void pade_calc_R_and_T_1(double  **N11, double  **N12, double  **N21, double  **N22,
                         double  **D11, double  **D12, double  **D21, double  **D22,
                         double  **R_p, double  **T_p, double  **R_m, double  **T_m,
                         double ***O11, double ***O12, double ***O21, double ***O22,
                         double ***E11, double ***E12, double ***E21, double ***E22,
                         double ***R_p_l, double ***T_p_l, double ***R_m_l, double ***T_m_l,
                         int n_quad, int n_stokes, int n_derivs,
                         int symmetric, int reduced, int vector,
                         int check_condition, double *condition,
                         uchar *derivs, work_data work) {

     int i;
     int ii;
     int j;
     int jj;
     int k;
     int kk;

     int n_quad_v;
     int n_quad_v2;

     int *i1;
     int *i2;

     double a;

     double **w11;

     double **w21;
     double **w22;
     double **w23;
     double **w24;

     double **lu2;

     double **P11;
     double **P12;
     double **P21;
     double **P22;

     double **Q11;
     double **Q12;
     double **Q21;
     double **Q22;


     n_quad_v = n_quad * n_stokes;

     n_quad_v2 = n_quad_v * 2;


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     i1  = get_work1(&work, WORK_IX2);

     w21 = get_work1(&work, WORK_DXX2);
     w22 = get_work1(&work, WORK_DXX2);

     lu2 = get_work1(&work, WORK_DXX2);

     P11 = get_work1(&work, WORK_DXX);
     P12 = get_work1(&work, WORK_DXX);
     P21 = get_work1(&work, WORK_DXX);
     P22 = get_work1(&work, WORK_DXX);

     if (derivs) {
          w23 = get_work1(&work, WORK_DXX2);
          w24 = get_work1(&work, WORK_DXX2);

          Q11 = get_work1(&work, WORK_DXX);
          Q12 = get_work1(&work, WORK_DXX);
          Q21 = get_work1(&work, WORK_DXX);
          Q22 = get_work1(&work, WORK_DXX);
     }

     if (derivs || vector)
          w11 = get_work1(&work, WORK_DXX);

     if (check_condition)
          i2 = get_work1(&work, WORK_IX);


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
     ii = n_quad_v;
     for (i = 0; i < n_quad_v; ++i) {
          jj = n_quad_v;
          for (j = 0; j < n_quad_v; ++j) {
/*
               w21[i ][j ] = N11[i][j];
               w21[i ][jj] = N12[i][j];
               w21[ii][j ] = N21[i][j];
               w21[ii][jj] = N22[i][j];
*/
               w21[i ][ j] = N12[i][j];
               w21[ii][ j] = N22[i][j];

               lu2[i ][j ] = D11[i][j];
               lu2[i ][jj] = D12[i][j];
               lu2[ii][j ] = D21[i][j];
               lu2[ii][jj] = D22[i][j];

               jj++;
          }
          ii++;
     }

     /* w24 = P = D-1 * N */
     if (check_condition)
          a = dmat_p_one_norm(lu2, n_quad_v, n_quad_v);

     dmat_getrf(lu2, n_quad_v2, n_quad_v2, i1);

     if (check_condition) {
          a = dmat_gecon('1', lu2, n_quad_v, a, i2, *w22);
          *condition = a;
/*
          if (a < DBL_EPSILON) {
               fprintf(stderr, "ERROR: pade condition < machine precision\n");
               return -1;
          }
*/
     }

     dmat_trans(w21, w22, n_quad_v2, n_quad_v2);
/*
     dmat_getrs(lu2, w22, n_quad_v2, n_quad_v2, i1);
*/
     dmat_getrs(lu2, w22, n_quad_v2, n_quad_v , i1);
     dmat_trans(w22, w21, n_quad_v2, n_quad_v2);

     ii = n_quad_v;
     for (i = 0; i < n_quad_v; ++i) {
          jj = n_quad_v;
          for (j = 0; j < n_quad_v; ++j) {
/*
               P11[i][j] = w21[i ][j ];
               P12[i][j] = w21[i ][jj];
               P21[i][j] = w21[ii][j ];
               P22[i][j] = w21[ii][jj];
*/
               P12[i][j] = w21[i ][ j];
               P22[i][j] = w21[ii][ j];
               jj++;
          }
          ii++;
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/

if (1) {
/*
if (! symmetric) {
*/
     /* T = P22-1 */
     dmat_copy(T_m, P22, n_quad_v, n_quad_v);
     dmat_getrf(T_m, n_quad_v, n_quad_v, i1);
     dmat_getri(T_m, n_quad_v, i1);
}
else {
     dmat_copy(T_m, P22, n_quad_v, n_quad_v);
     dmat_potrf(T_m, n_quad_v);
     dmat_potri(T_m, n_quad_v);
     for (i = 0; i < n_quad_v; ++i) {
          for (j = i; j < n_quad_v; ++j) {
               T_m[j][i] = T_m[i][j];
          }
     }
}
          /* R = P12 * T */
if (! symmetric)
     dmat_mul(P12, T_m, n_quad_v, n_quad_v, n_quad_v, R_m);
else
     dsym_mul(P12, T_m, n_quad_v, n_quad_v, R_m, 1);


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (vector) {
if (0) {
          w11 = get_work1(&work, WORK_DXX);

          /* R_p = -T_m * P21 */
          dmat_mul(T_m, P21, n_quad_v, n_quad_v, n_quad_v, R_p);
          dmat_scale(-1., R_p, R_p, n_quad_v, n_quad_v);

          /* T_p = P11 + P12 * R_p */
          dmat_mul(P12, R_p, n_quad_v, n_quad_v, n_quad_v, w11);
          dmat_add(P11, w11, T_p, n_quad_v, n_quad_v);
}
else {
          phase_matrix_symmetry2(n_quad, n_stokes, T_m, R_m, T_p, R_p, 1.);
}
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (flags_or(derivs, n_derivs)) {
          if (reduced) {
               E22 = O11;
               E21 = O12;
               E12 = O21;
               E11 = O22;
          }

          for (i = 0; i < n_derivs; ++i) {

               if (! derivs[i])
                    continue;

               jj = n_quad_v;
               for (j = 0; j < n_quad_v; ++j) {
                    kk = n_quad_v;
                    for (k = 0; k < n_quad_v; ++k) {
                         w22[j ][k ] = O11[i][j][k];
                         w22[j ][kk] = O12[i][j][k];
                         w22[jj][k ] = O21[i][j][k];
                         w22[jj][kk] = O22[i][j][k];

                         w23[j ][k ] = E11[i][j][k];
                         w23[j ][kk] = E12[i][j][k];
                         w23[jj][k ] = E21[i][j][k];
                         w23[jj][kk] = E22[i][j][k];

                         kk++;
                    }
                    jj++;
               }

               /* w28 = L(P) = L(D-1 * N) = D-1 * (L(N) - L(D) * P) */
               dmat_mul(w23, w21, n_quad_v2, n_quad_v2, n_quad_v2, w24);
               dmat_sub(w22, w24, w23, n_quad_v2, n_quad_v2);
               dmat_trans(w23, w24, n_quad_v2, n_quad_v2);
               dmat_getrs(lu2, w24, n_quad_v2, n_quad_v2, i1);
               dmat_trans(w24, w22, n_quad_v2, n_quad_v2);

               jj = n_quad_v;
               for (j = 0; j < n_quad_v; ++j) {
                    kk = n_quad_v;
                    for (k = 0; k < n_quad_v; ++k) {
                         Q11[j][k] = w22[j ][k ];
                         Q12[j][k] = w22[j ][kk];
                         Q21[j][k] = w22[jj][k ];
                         Q22[j][k] = w22[jj][kk];
                         kk++;
                    }
                    jj++;
               }


               /*---------------------------------------------------------------
                *
                *-------------------------------------------------------------*/
               /* W = T * Q22 * T */
               dmat_mul(T_m, Q22, n_quad_v, n_quad_v, n_quad_v, w11);
if (! symmetric)
               dmat_mul(w11, T_m, n_quad_v, n_quad_v, n_quad_v, T_m_l[i]);
else
               dsym_mul(w11, T_m, n_quad_v, n_quad_v, T_m_l[i], 1);
               dmat_scale(-1., T_m_l[i], T_m_l[i], n_quad_v, n_quad_v);

               /* U = Q12 * T + P22 * W */
if (! symmetric) {
               dmat_mul(Q12, T_m, n_quad_v, n_quad_v, n_quad_v, R_m_l[i]);
               dmat_mul(P12, T_m_l[i], n_quad_v, n_quad_v, n_quad_v, w11);
               dmat_add(R_m_l[i], w11, R_m_l[i], n_quad_v, n_quad_v);
}
else
               dsym_mma(Q12, T_m, P12, T_m_l[i], n_quad_v, n_quad_v, R_m_l[i], 1);


              /*---------------------------------------------------------------
               *
               *-------------------------------------------------------------*/
               if (vector) {
if (0) {
                    /* R_p_l = -T_m_l * P21 - T_m * Q21 */
                    dmat_mul(T_m_l[i], P21, n_quad_v, n_quad_v, n_quad_v, R_p_l[i]);
                    dmat_scale(-1., R_p_l[i], R_p_l[i], n_quad_v, n_quad_v);
                    dmat_mul(T_m, Q21, n_quad_v, n_quad_v, n_quad_v, w11);
                    dmat_sub(R_p_l[i], w11, R_p_l[i], n_quad_v, n_quad_v);

                    /* T_p_l = Q11 + Q12 * R_p + P12 * R_p_l */
                    dmat_mul(Q12, R_p, n_quad_v, n_quad_v, n_quad_v, w11);
                    dmat_add(Q11, w11, T_p_l[i], n_quad_v, n_quad_v);
                    dmat_mul(P12, R_p_l[i], n_quad_v, n_quad_v, n_quad_v, w11);
                    dmat_add(T_p_l[i], w11, T_p_l[i], n_quad_v, n_quad_v);
}
else {
                    phase_matrix_symmetry2(n_quad, n_stokes, T_m_l[i], R_m_l[i], T_p_l[i], R_p_l[i], 1.);
}
               }
          }
     }
}



/*******************************************************************************
 *
 ******************************************************************************/
void rtm_pade_rts(int n_quad, int n_stokes, int n_derivs, double F_0,
                  double *qx_v, double *qwv, double *umus, int n_umus,
                  double planck0, double planck1, double *planck0_l, double *planck1_l,
                  double omega, double *omega_l, double ltau, double *ltau_l,
                  double btran, double *btran_l,
                  double as_0, double *as_0_l, double atran, double *atran_l,
                  double  *P_0p, double  *P_0m,
                  double  **r_p, double  **t_p, double  **r_m, double  **t_m,
                  double  **R_p, double  **T_p, double  **R_m, double  **T_m,
                  double  *S_p, double  *S_m,
                  double **P_0p_l, double **P_0m_l,
                  double ***r_p_l, double ***t_p_l, double ***r_m_l, double ***t_m_l,
                  double ***R_p_l, double ***T_p_l, double ***R_m_l, double ***T_m_l,
                  double **S_p_l, double **S_m_l,
                  int pade_s, int pade_r, int symmetric, int vector,
                  int thermal, int check_condition, double *condition,
                  uchar *derivs_layers, uchar *derivs_beam_down, save_tree_data save_tree, work_data work) {

     int i;
     int j;
     int k;
/*
     int flag;
*/
     int n_quad_v;

     int m_pade;

     int symmetric2;

     double a;
     double b;

     double pade_c[MAX_PADE_R];

     double *F_p;
     double *F_m;

     double *F0_p;
     double *F0_m;
     double *F1_p;
     double *F1_m;

     double *St_p;
     double *St_m;

     double **F_p_l;
     double **F_m_l;

     double **A11;
     double **A12;
     double **A21;
     double **A22;

     double ***B11;
     double ***B12;
     double ***B21;
     double ***B22;

     double **N11;
     double **N12;
     double **N21;
     double **N22;

     double **D11;
     double **D12;
     double **D21;
     double **D22;

     double ***O11;
     double ***O12;
     double ***O21;
     double ***O22;

     double ***E11;
     double ***E12;
     double ***E21;
     double ***E22;
/*
     double **R_p2;
     double **T_p2;
     double **R_m2;
     double **T_m2;

     double ***R_p_l2;
     double ***T_p_l2;
     double ***R_m_l2;
     double ***T_m_l2;
*/
     double **tpr;
     double **tmr;
     double **gamma;

     double ***tpr_l;
     double ***tmr_l;
     double ***gamma_l;


     n_quad_v = n_quad * n_stokes;


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     F_p = get_work1(&work, WORK_DX);
     F_m = get_work1(&work, WORK_DX);

     if (thermal) {
          F0_p = get_work1(&work, WORK_DX);
          F0_m = get_work1(&work, WORK_DX);
          F1_p = get_work1(&work, WORK_DX);
          F1_m = get_work1(&work, WORK_DX);

          St_p = get_work1(&work, WORK_DX);
          St_m = get_work1(&work, WORK_DX);
     }

     A11 = get_work1(&work, WORK_DXX);
     A12 = get_work1(&work, WORK_DXX);
     A21 = get_work1(&work, WORK_DXX);
     A22 = get_work1(&work, WORK_DXX);

     N11 = get_work1(&work, WORK_DXX);
     N12 = get_work1(&work, WORK_DXX);
     N21 = get_work1(&work, WORK_DXX);
     N22 = get_work1(&work, WORK_DXX);

     D11 = get_work1(&work, WORK_DXX);
     D12 = get_work1(&work, WORK_DXX);
     D21 = get_work1(&work, WORK_DXX);
     D22 = get_work1(&work, WORK_DXX);

     if (flags_or(derivs_layers, n_derivs)) {
          B11 = get_work2(&work, WORK_DXX, WORK_DERIVS_V, derivs_layers);
          B12 = get_work2(&work, WORK_DXX, WORK_DERIVS_V, derivs_layers);
          B21 = get_work2(&work, WORK_DXX, WORK_DERIVS_V, derivs_layers);
          B22 = get_work2(&work, WORK_DXX, WORK_DERIVS_V, derivs_layers);

          O11 = get_work2(&work, WORK_DXX, WORK_DERIVS_V, derivs_layers);
          O12 = get_work2(&work, WORK_DXX, WORK_DERIVS_V, derivs_layers);
          O21 = get_work2(&work, WORK_DXX, WORK_DERIVS_V, derivs_layers);
          O22 = get_work2(&work, WORK_DXX, WORK_DERIVS_V, derivs_layers);

          E11 = get_work2(&work, WORK_DXX, WORK_DERIVS_V, derivs_layers);
          E12 = get_work2(&work, WORK_DXX, WORK_DERIVS_V, derivs_layers);
          E21 = get_work2(&work, WORK_DXX, WORK_DERIVS_V, derivs_layers);
          E22 = get_work2(&work, WORK_DXX, WORK_DERIVS_V, derivs_layers);

     }

     if (flags_or(derivs_beam_down, n_derivs)) {
          F_p_l = get_work2(&work, WORK_DX, WORK_DERIVS_V, derivs_beam_down);
          F_m_l = get_work2(&work, WORK_DX, WORK_DERIVS_V, derivs_beam_down);
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (pade_s < 0 || pade_r <= 0)
          pade_get_s_and_r(r_p, t_p, n_quad_v, ltau, n_umus, umus, &pade_s, &pade_r, flags_or(derivs_beam_down, n_derivs));


     pade_coefs_x(pade_r, pade_r, 1., pade_c);


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     m_pade = 1 << pade_s;

     a = ltau / m_pade;

     for (i = 0; i < n_quad_v; ++i) {
          for (j = 0; j < n_quad_v; ++j) {
               A11[i][j] =  t_p[i][j] * a;
               A12[i][j] = -r_m[i][j] * a;
               A21[i][j] =  r_p[i][j] * a;
               A22[i][j] = -t_m[i][j] * a;
          }
     }

     if (flags_or(derivs_layers, n_derivs)) {
          for (i = 0; i < n_derivs; ++i) {
               if (! derivs_layers[i])
                    continue;

               b = ltau_l[i] / m_pade;

               for (j = 0; j < n_quad_v; ++j) {
                    for (k = 0; k < n_quad_v; ++k) {
                         B11[i][j][k] =  t_p_l[i][j][k] * a + t_p[j][k] * b;
                         B12[i][j][k] = -r_m_l[i][j][k] * a - r_m[j][k] * b;
                         B21[i][j][k] =  r_p_l[i][j][k] * a + r_p[j][k] * b;
                         B22[i][j][k] = -t_m_l[i][j][k] * a - t_m[j][k] * b;
                    }
               }
          }
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     calc_n_d_2(n_quad, n_stokes, n_derivs, pade_c, pade_r,
                A11, A12, A21, A22, B11, B12, B21, B22,
                N11, N12, N21, N22, D11, D12, D21, D22,
                O11, O12, O21, O22, E11, E12, E21, E22,
                symmetric, derivs_layers, work);


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
if (! vector)
     symmetric2 = symmetric;
else
     symmetric2 = 0;
/*
     if (! vector) {
          flag = 0;

          R_p2 = NULL;
          T_p2 = NULL;
          R_m2 = R_p;
          T_m2 = T_p;

          if (flags_or(derivs_layers, n_derivs)) {
               R_p_l2 = NULL;
               T_p_l2 = NULL;
               R_m_l2 = R_p_l;
               T_m_l2 = T_p_l;
          }
     }
     else {
          flag = 1;

          R_p2 = R_p;
          T_p2 = T_p;
          R_m2 = R_m;
          T_m2 = T_m;

          if (flags_or(derivs_layers, n_derivs)) {
               R_p_l2 = R_p_l;
               T_p_l2 = T_p_l;
               R_m_l2 = R_m_l;
               T_m_l2 = T_m_l;
          }
     }
*/
/*
     pade_calc_R_and_T_1(N11, N12, N21, N22, D11, D12, D21, D22,
                         R_p2, T_p2, R_m2, T_m2,
                         O11, O12, O21, O22, E11, E12, E21, E22,
                         R_p_l2, T_p_l2, R_m_l2, T_m_l2,
                         n_quad, n_stokes, n_derivs,
                         symmetric2, 0, flag, check_condition, condition, derivs_layers, work);

     pade_calc_R_and_T_2(N11, N12, N21, N22, D11, D12, D21, D22,
                         R_p2, T_p2, R_m2, T_m2,
                         O11, O12, O21, O22, E11, E12, E21, E22,
                         R_p_l2, T_p_l2, R_m_l2, T_m_l2,
                         n_quad, n_stokes, n_derivs,
                         symmetric2, 0, flag, check_condition, condition, derivs_layers, work);
*/

     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     for (i = 0; i < pade_s; ++i) {
          if (! vector) {
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
          else {
/*
               if (! symmetric2) {
*/
                    if (! flags_or(derivs_layers, n_derivs))
                         variant_double  (R_m, T_m, NULL, NULL, NULL, R_p, T_p, NULL, NULL, NULL, n_quad, n_stokes, 0., 0., 0, 0, i == 0, work);
                    else
                         variant_double_l(R_m, T_m, NULL, R_p, T_p, NULL, R_m_l, T_m_l, NULL, R_p_l, T_p_l, NULL, n_quad, n_stokes, n_derivs, 0., NULL, derivs_layers, work);
/*
               }
               else {
                    if (! flags_or(derivs_layers, n_derivs))
                         variant_double_s  (R_m, T_m, NULL, NULL, NULL, R_p, T_p, NULL, NULL, NULL, n_quad, n_stokes, 0., 0., 0, 0, i == 0, work);
                    else
                         variant_double_s_l(R_m, T_m, NULL, R_p, T_p, NULL, R_m_l, T_m_l, NULL, R_p_l, T_p_l, NULL, n_quad, n_stokes, n_derivs, 0., NULL, derivs_layers, work);
               }
*/
          }
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
if (1) {
     tpr   = get_work1(&work, WORK_DXX);
     tmr   = get_work1(&work, WORK_DXX);
     gamma = get_work1(&work, WORK_DXX);

     if (flags_or(derivs_layers, n_derivs)) {
          tpr_l   = get_work2(&work, WORK_DXX, WORK_DERIVS_V, derivs_layers);
          tmr_l   = get_work2(&work, WORK_DXX, WORK_DERIVS_V, derivs_layers);
          gamma_l = get_work2(&work, WORK_DXX, WORK_DERIVS_V, derivs_layers);
     }

     build_txr(n_quad, n_stokes, n_derivs, r_p, t_p, tpr, tmr, r_p_l, t_p_l, tpr_l, tmr_l, derivs_layers, work);
     build_gamma(n_quad_v, n_derivs, tpr, tmr, gamma, tpr_l, tmr_l, gamma_l, derivs_layers, work);

     build_source_vectors_solar_classic_1n(n_quad, n_stokes, n_derivs, qx_v, F_0, omega, omega_l, as_0, as_0_l, P_0p, P_0m, tpr, tmr, gamma, F_p, F_m, P_0p_l, P_0m_l, tpr_l, tmr_l, gamma_l, F_p_l, F_m_l, derivs_layers, derivs_beam_down, save_tree, work);
}
else {
     if (! vector)
          build_source_vectors_solar_classic_2n (n_quad_v, n_derivs, qx_v, F_0, omega, omega_l, as_0, as_0_l, P_0p, P_0m, r_p, t_p, F_p, F_m, P_0p_l, P_0m_l, r_p_l, t_p_l, F_p_l, F_m_l, derivs_layers, derivs_beam_down, work);
     else
          build_source_vectors_solar_classic_2n2(n_quad_v, n_derivs, qx_v, F_0, omega, omega_l, as_0, as_0_l, P_0p, P_0m, r_p, t_p, r_m, t_m, F_p, F_m, P_0p_l, P_0m_l, r_p_l, t_p_l, r_m_l, t_m_l, F_p_l, F_m_l, derivs_layers, derivs_beam_down, work);
}
     if (! vector)
          build_global_source_solar (n_quad_v, n_derivs, atran, atran_l, R_p, T_p, F_p, F_m, S_p, S_m, R_p_l, T_p_l, F_p_l, F_m_l, S_p_l, S_m_l, derivs_layers, derivs_beam_down, work, NULL, NULL, NULL, NULL);
     else
          build_global_source_solar2(n_quad_v, n_derivs, atran, atran_l, R_p, T_p, R_m, T_m, F_p, F_m, S_p, S_m, R_p_l, T_p_l, R_m_l, T_m_l, F_p_l, F_m_l, S_p_l, S_m_l, derivs_layers, derivs_beam_down, work, NULL, NULL, NULL, NULL);


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
if (thermal) {
     build_source_vectors_thermal2(n_quad, n_stokes, n_derivs, qx_v, planck0, planck1, planck0_l, planck1_l, omega, omega_l, ltau, ltau_l, r_p, t_p, r_m, t_m, F0_p, F0_m, F1_p, F1_m, r_p_l, t_p_l, r_m_l, t_m_l, NULL, NULL, NULL, NULL, derivs_layers, derivs_beam_down, work);

     if (! vector)
          build_global_source_thermal(n_quad_v, n_derivs, R_p, T_p, R_p, T_p, F0_p, F0_m, F1_p, F1_m, St_p, St_m, R_p_l, T_p_l, R_p_l, T_p_l, NULL, NULL, NULL, NULL, S_p_l, S_m_l, derivs_layers, derivs_beam_down, work);
     else
          build_global_source_thermal(n_quad_v, n_derivs, R_p, T_p, R_m, T_m, F0_p, F0_m, F1_p, F1_m, St_p, St_m, R_p_l, T_p_l, R_m_l, T_m_l, NULL, NULL, NULL, NULL, S_p_l, S_m_l, derivs_layers, derivs_beam_down, work);
}

     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (F_0 == 0. && ! thermal) {
          dvec_zero(S_p, n_quad_v);
          dvec_zero(S_m, n_quad_v);
     }
     else
     if (F_0 == 0. && thermal) {
          dvec_copy(S_p, St_p, n_quad_v);
          dvec_copy(S_m, St_m, n_quad_v);
     }
     else
     if (F_0 != 0. && thermal) {
          dvec_add (S_p, St_p, S_p, n_quad_v);
          dvec_add (S_m, St_m, S_m, n_quad_v);
     }
}
