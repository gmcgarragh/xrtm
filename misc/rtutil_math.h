/*******************************************************************************
**
**    Copyright (C) 2007-2020 Greg McGarragh <greg.mcgarragh@colostate.edu>
**
**    This source code is licensed under the GNU General Public License (GPL),
**    Version 3.  See the file COPYING for more details.
**
*******************************************************************************/

#ifndef RTUTIL_MATH_H
#define RTUTIL_MATH_H

#ifdef __cplusplus
extern "C" {
#endif


double factorial(int x);
double factorial_evaluate(double y, int x1, int x2);

void extended_trapezoidal_quad(int n, double a, double b, double *x, double *w);

void gauss_leg_quad(int n, double *x, double *w);
void gauss_leg_quad2(int n, double *x, double *w);
void gauss_leg_quadx(int n, double x1, double x2, double *x, double *w);
void gauss_leg_quadx_l(int n, int n_derivs,
                       double x1, double x2, double *x1_l, double *x2_l,
                       double *x, double *w, double **x_l, double **w_l);
void gauss_leg_quad_fit(int n, double x1, double x2,
                       double *xa, double *wa, double *xb, double *wb);
void gauss_leg_quad_fit_l(int n, int n_derivs,
                          double x1, double x2, double *x1_l, double *x2_l,
                          double *xa, double *wa, double *xb, double *wb,
                          double **xa_l, double **wa_l, double **xb_l, double **wb_l);

int simpsons_rule_quad(int n, double a, double b, double *x, double *w);

void trapezoidal_quad(int n, double a, double b, double *x, double *w);

void leg_poly(int n_p, double mu, double *p);
void leg_poly2(int n_mu, int n_l, double *mu, double **p);
void leg_poly_assoc(int m, int n_y, double mu,
                    double *c, double *d, double *y, int flag);
void leg_poly_assoc2(int m, int n_mu, int n_y, double *mu,
                     double *c, double *d, double **y, int flag);
void gen_spher_funcs(int p, int q, int n_l, double mu, dcomplex *z);
void gen_spher_funcs2(int p, int q, int n_mu, int n_l, double *mu, dcomplex **z);
void gen_spher_funcs_0q(int q, int n_l, double mu, dcomplex *gsf);
void gen_spher_funcs_0q_2(int q, int n_mu, int n_l, double *mu, dcomplex **gsf);
void gen_spher_funcs_p0(int p, int n_l, double mu, dcomplex *gsf);
void gen_spher_funcs_p0_2(int p, int n_mu, int n_l, double *mu, dcomplex **gsf);
void gen_spher_funcs_prt(int n_l, double mu,
                         double *p00, double *p0p2, double *p2p2, double *p2m2);
void gen_spher_funcs_prt2(int n_mu, int n_l, double *mu,
                          double **p00, double **p0p2, double **p2p2, double **p2m2);


#ifdef __cplusplus
}
#endif

#endif /* RTUTIL_MATH_H */
