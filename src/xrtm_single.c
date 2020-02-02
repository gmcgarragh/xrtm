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
#include "xrtm_save_tree.h"
#include "xrtm_single.h"
#include "xrtm_single_a.h"
#include "xrtm_utility.h"


/*******************************************************************************
 *
 ******************************************************************************/
static void ssr_up_layer(int i_layer, int n_stokes, int n_derivs, double utau, double *umus, int n_umus, double *omega, double **omega_l, double *ltau, double **ltau_l, double *btran, double **btran_l, double *as_0, double **as_0_l, double **P, double ***P_l, double *I, double **I_l, double *I_ss, double **I_ss_l, uchar **derivs_layers, work_data work) {

     int i;
     int ii;
     int j;
     int jj;
     int k;
     int kk;

     double a;
     double b;
     double *b_l;
     double c;
     double *c_l;
     double d;
     double *d_l;
     double e;
     double e_l;
     double f;
     double g;
     double h;
     double p;
     double p_l;

     double utau_l;
#ifdef USE_AD_FOR_TL_SSR_UP_LAYER
     ssr_up_layer_tl_with_ad(i_layer, n_stokes, n_derivs, utau, umus, n_umus, omega, omega_l, ltau, ltau_l, btran, btran_l, as_0, as_0_l, P, P_l, I, I_l, I_ss, I_ss_l, derivs_layers, work);
#endif
     b_l = get_work1(&work, WORK_DDERIVS);
     c_l = get_work1(&work, WORK_DDERIVS);
     d_l = get_work1(&work, WORK_DDERIVS);

     a = exp(-utau * as_0[i_layer]);
     b = btran[i_layer] * a * omega[i_layer];
     c = ltau[i_layer] - utau;
     d = exp(-c * as_0[i_layer]);

     f = a * omega[i_layer];
     g = btran[i_layer] * a;
     for (i = 0; i < n_derivs; ++i) {
          utau_l = utau * ltau_l[i_layer][i] / ltau[i_layer];

          b_l[i] = (btran_l[i_layer][i] + btran[i_layer] * (-utau_l * as_0[i_layer] -utau * as_0_l[i_layer][i])) * f + g * omega_l[i_layer][i];
          c_l[i] = ltau_l[i_layer][i] - utau_l;
          d_l[i] = (-c_l[i] * as_0[i_layer] - c * as_0_l[i_layer][i]) * d;
     }

     for (i = 0; i < n_umus; ++i) {
          ii = i * n_stokes;

          e = exp(-c / umus[i]);
          f = 1. / umus[i] + as_0[i_layer];
          g = (1. - e * d) / f;
          h = b * g;

          p = e / umus[i];
          for (j = 0; j < n_derivs; ++j) {
               e_l = -c_l[j] * p;
               p_l = g * b_l[j] + b * (-e_l * d - e * d_l[j] - g * as_0_l[i_layer][j]) / f;

               for (k = 0; k < n_stokes; ++k) {
                    kk = ii + k;
#ifndef USE_AD_FOR_TL_SSR_UP_LAYER
                    I_ss_l[j][kk] = e_l * I[kk] + e * I_l[j][kk] + P[i_layer][kk] * p_l;
                    if (derivs_layers[i_layer][j])
                         I_ss_l[j][kk] += b * P_l[i_layer][j][kk] * g;
#endif
               }
          }

          for (j = 0; j < n_stokes; ++j) {
               jj = ii + j;
               I_ss[jj] = e * I[jj] + h * P[i_layer][jj];
          }
     }
}



/*******************************************************************************
 *
 ******************************************************************************/
void single_scattered_radiance_up(int n_stokes, int n_derivs, int n_layers, double F_0, int n_ulevels, int *ulevels, double *utaus, double *umus, int n_umus, double *omega, double **omega_l, double *ltau, double **ltau_l, double *btran, double **btran_l, double *as_0, double **as_0_l, double *atran, double **atran_l, double **P, double ***P_l, double *I_in, double **I_in_l, double **I_ss, double ***I_ss_l, int utau_output, uchar **derivs_layers, uchar **derivs_beam, save_tree_data save_tree, work_data work) {

     int i;
     int ii;
     int i1;
     int i2;
     int j;
     int jj;
     int k;
     int kk;
     int l;

     int n_umus_v;

     int i_ulevel;

     double a;
     double b;

     double *I;
     double **I_l;

     double **I_ss_l2;

     forward_save_single_scattered_radiance_data *save;


     n_umus_v = n_umus * n_stokes;


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (F_0 == 0.) {
          init_array2_d(I_ss, n_ulevels, n_umus_v, 0.);
          if (n_derivs > 0)
               init_array3_d(I_ss_l, n_ulevels, n_derivs, n_umus_v, 0.);
          return;
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     I = get_work_d1(&work, n_umus * n_stokes);
     if (n_derivs > 0)
          I_l = get_work_d2(&work, n_derivs, n_umus * n_stokes);


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (save_tree.t) {
          save_tree_encode_s(&save_tree, "single_scattered_radiance_up");

          if (save_tree_retrieve_data(&save_tree, forward_save_single_scattered_radiance_data, &save))
               forward_save_single_scattered_radiance_alloc(save, n_layers, n_umus * n_stokes);
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     a = F_0 / (4. * PI);


     for (i = 0; i < n_umus; ++i) {
          ii = i * n_stokes;

          b = 1. / umus[i] * a;

          for (j = 0; j < n_stokes; ++j) {
               jj = ii + j;
               I[jj] = I_in[jj] / b;
          }

          for (j = 0; j < n_derivs; ++j)
               for (k = 0; k < n_stokes; ++k) {
                    kk = ii + k;
                    I_l[j][kk] = I_in_l[j][kk] / b;
          }
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     i1       = n_layers - 1;
     i2       = ulevels[0];
     i_ulevel = n_ulevels - 1;


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if ((! utau_output && i1 + 1 == ulevels[i_ulevel]) || (utau_output && i1 + 1 == ulevels[i_ulevel] && utaus[i_ulevel] == 0.)) {
          for (i = 0; i < n_umus_v; ++i)
               I_ss[i_ulevel][i] = I[i];

          for (i = 0; i < n_derivs; ++i) {
               for (j = 0; j < n_umus_v; ++j)
                    I_ss_l[i_ulevel][i][j] = I_l[i][j];
          }

          i_ulevel--;
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     for (i = i1; i >= i2; --i) {
          if (  utau_output && i_ulevel >= 0 && i == ulevels[i_ulevel]) {
               while (i_ulevel >= 0 && i == ulevels[i_ulevel]) {
                    if (n_derivs > 0)
                         I_ss_l2 = I_ss_l[i_ulevel];

                    ssr_up_layer(i, n_stokes, n_derivs, utaus[i_ulevel], umus, n_umus, omega, omega_l, ltau, ltau_l, btran, btran_l, as_0, as_0_l, P, P_l, I, I_l, I_ss[i_ulevel], I_ss_l2, derivs_layers, work);

                    i_ulevel--;
               }
          }

          if (save_tree.t)
               dvec_copy(save->I[i], I, n_umus * n_stokes);

          ssr_up_layer(i, n_stokes, n_derivs, 0., umus, n_umus, omega, omega_l, ltau, ltau_l, btran, btran_l, as_0, as_0_l, P, P_l, I, I_l, I, I_l, derivs_layers, work);

          if (! utau_output && i == ulevels[i_ulevel]) {
               for (j = 0; j < n_umus_v; ++j)
                    I_ss[i_ulevel][j] = I[j];

              for (j = 0; j < n_derivs; ++j) {
                    for (k = 0; k < n_umus_v; ++k)
                         I_ss_l[i_ulevel][j][k] = I_l[j][k];
               }

               i_ulevel--;
          }
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     for (i = 0; i < n_ulevels; ++i) {
          for (j = 0; j < n_umus; ++j) {
               ii = j * n_stokes;

               b = 1. / umus[j] * a;

               for (k = 0; k < n_stokes; ++k)
                    I_ss[i][ii + k] *= b;

               for (k = 0; k < n_derivs; ++k) {
                    for (l = 0; l < n_stokes; ++l) {
                         I_ss_l[i][k][ii + l] *= b;
                    }
               }
          }
     }
#ifdef USE_AD_FOR_TL_SINGLE_SCATTERED_RADIANCE_UP
     single_scattered_radiance_up_tl_with_ad(n_stokes, n_derivs, n_layers, F_0, n_ulevels, ulevels, utaus, umus, n_umus, omega, omega_l, ltau, ltau_l, btran, btran_l, as_0, as_0_l, atran, atran_l, P, P_l, I_in, I_in_l, I_ss, I_ss_l, utau_output, derivs_layers, derivs_beam, save_tree, work);
#endif
}



/*******************************************************************************
 *
 ******************************************************************************/
static void ssr_dn_layer(int i_layer, int n_stokes, int n_derivs, double utau, double *umus, int n_umus, double *omega, double **omega_l, double *ltau, double **ltau_l, double *btran, double **btran_l, double *as_0, double **as_0_l, double **P, double ***P_l, double *I, double **I_l, double *I_ss, double **I_ss_l, uchar **derivs_layers, work_data work) {

     int i;
     int ii;
     int j;
     int jj;
     int k;
     int kk;

     double a;
     double *a_l;
     double b;
     double *b_l;
     double c;
     double c_l;
     double d;
     double e;
     double f;
     double g;
     double g_l;

     double *utau_l;
#ifdef USE_AD_FOR_TL_SSR_DN_LAYER
     ssr_dn_layer_tl_with_ad(i_layer, n_stokes, n_derivs, utau, umus, n_umus, omega, omega_l, ltau, ltau_l, btran, btran_l, as_0, as_0_l, P, P_l, I, I_l, I_ss, I_ss_l, derivs_layers, work);
#endif
     a_l    = get_work1(&work, WORK_DDERIVS);
     b_l    = get_work1(&work, WORK_DDERIVS);
     utau_l = get_work1(&work, WORK_DDERIVS);

     a = btran[i_layer] * omega[i_layer];

     b = exp(-utau * as_0[i_layer]);

     for (i = 0; i < n_derivs; ++i) {
          utau_l[i] = utau * ltau_l[i_layer][i] / ltau[i_layer];

          a_l[i] = btran_l[i_layer][i] * omega[i_layer] + btran[i_layer] * omega_l[i_layer][i];
          b_l[i] = (-utau_l[i] * as_0[i_layer] - utau * as_0_l[i_layer][i]) * b;
     }

     for (i = 0; i < n_umus; ++i) {
          ii = i * n_stokes;

          c = exp(-utau    / umus[i]);
          d = as_0[i_layer] - 1. / umus[i];
          e = (c - b) / d;
          f = a * e;

          g = c / umus[i];
          for (j = 0; j < n_derivs; ++j) {
               c_l = -utau_l[j] * g;

               g_l = e * a_l[j] + a * (c_l - b_l[j] - e * as_0_l[i_layer][j]) / d;

               for (k = 0; k < n_stokes; ++k) {
                    kk = ii + k;
#ifndef USE_AD_FOR_TL_SSR_DN_LAYER
                    I_ss_l[j][kk] = c_l * I[kk] + c * I_l[j][kk] + P[i_layer][kk] * g_l;
                    if (derivs_layers[i_layer][j])
                         I_ss_l[j][kk] += a * P_l[i_layer][j][kk] * e;
#endif
               }
          }

          for (j = 0; j < n_stokes; ++j) {
               jj = ii + j;
               I_ss[jj] = c * I[jj] + f * P[i_layer][jj];
          }
     }
}



/*******************************************************************************
 *
 ******************************************************************************/
void single_scattered_radiance_dn(int n_stokes, int n_derivs, int n_layers, double F_0, int n_ulevels, int *ulevels, double *utaus, double *umus, int n_umus, double *omega, double **omega_l, double *ltau, double **ltau_l, double *btran, double **btran_l, double *as_0, double **as_0_l, double *atran, double **atran_l, double **P, double ***P_l, double *I_in, double **I_in_l, double **I_ss, double ***I_ss_l, int utau_output, uchar **derivs_layers, uchar **derivs_beam, save_tree_data save_tree, work_data work) {

     int i;
     int ii;
     int i1;
     int i2;
     int j;
     int jj;
     int k;
     int kk;
     int l;

     int n_umus_v;

     int i_ulevel;

     double a;
     double b;

     double *I;
     double **I_l;

     double **I_ss_l2;

     forward_save_single_scattered_radiance_data *save;


     n_umus_v = n_umus * n_stokes;


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (F_0 == 0.) {
          init_array2_d(I_ss, n_ulevels, n_umus_v, 0.);
          if (n_derivs > 0)
               init_array3_d(I_ss_l, n_ulevels, n_derivs, n_umus_v, 0.);
          return;
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     I = get_work_d1(&work, n_umus * n_stokes);
     if (n_derivs > 0)
          I_l = get_work_d2(&work, n_derivs, n_umus * n_stokes);


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (save_tree.t) {
          save_tree_encode_s(&save_tree, "single_scattered_radiance_dn");

          if (save_tree_retrieve_data(&save_tree, forward_save_single_scattered_radiance_data, &save))
               forward_save_single_scattered_radiance_alloc(save, n_layers, n_umus * n_stokes);
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     a = F_0 / (4. * PI);


     for (i = 0; i < n_umus; ++i) {
          ii = i * n_stokes;

          b = 1. / umus[i] * a;

          for (j = 0; j < n_stokes; ++j) {
               jj = ii + j;
               I[jj] = I_in[jj] / b;
          }

          for (j = 0; j < n_derivs; ++j)
               for (k = 0; k < n_stokes; ++k) {
                    kk = ii + k;
                    I_l[j][kk] = I_in_l[j][kk] / b;
          }
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     i1       = 0;
     if (! utau_output)
          i2 = ulevels[n_ulevels - 1];
     else
          i2 = ulevels[n_ulevels - 1] + 1;
     i_ulevel = 0;


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if ((! utau_output && i1 == ulevels[i_ulevel]) || (utau_output && i1 == ulevels[i_ulevel] && utaus[i_ulevel] == 0.)) {
          for (i = 0; i < n_umus_v; ++i)
               I_ss[i_ulevel][i] = I[i];

          for (i = 0; i < n_derivs; ++i) {
               for (j = 0; j < n_umus_v; ++j)
                    I_ss_l[i_ulevel][i][j] = I_l[i][j];
          }

          i_ulevel++;
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     for (i = i1; i < i2; ++i) {
          if (  utau_output && i_ulevel < n_ulevels && i == ulevels[i_ulevel]) {
               while (i_ulevel < n_ulevels && i == ulevels[i_ulevel]) {
                    if (n_derivs > 0)
                         I_ss_l2 = I_ss_l[i_ulevel];

                    ssr_dn_layer(i, n_stokes, n_derivs, utaus[i_ulevel], umus, n_umus, omega, omega_l, ltau, ltau_l, btran, btran_l, as_0, as_0_l, P, P_l, I, I_l, I_ss[i_ulevel], I_ss_l2, derivs_layers, work);

                    i_ulevel++;
               }
          }

          if (save_tree.t)
               dvec_copy(save->I[i], I, n_umus * n_stokes);

          ssr_dn_layer(i, n_stokes, n_derivs, ltau[i], umus, n_umus, omega, omega_l, ltau, ltau_l, btran, btran_l, as_0, as_0_l, P, P_l, I, I_l, I, I_l, derivs_layers, work);

          if (! utau_output && i + 1 == ulevels[i_ulevel]) {
               for (j = 0; j < n_umus_v; ++j)
                    I_ss[i_ulevel][j] = I[j];

              for (j = 0; j < n_derivs; ++j) {
                    for (k = 0; k < n_umus_v; ++k)
                         I_ss_l[i_ulevel][j][k] = I_l[j][k];
               }

               i_ulevel++;
          }
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     for (i = 0; i < n_ulevels; ++i) {
          for (j = 0; j < n_umus; ++j) {
               ii = j * n_stokes;

               b = 1. / umus[j] * a;

               for (k = 0; k < n_stokes; ++k)
                    I_ss[i][ii + k] *= b;

               for (k = 0; k < n_derivs; ++k) {
                    for (l = 0; l < n_stokes; ++l) {
                         I_ss_l[i][k][ii + l] *= b;
                    }
               }
          }
     }
#ifdef USE_AD_FOR_TL_SINGLE_SCATTERED_RADIANCE_DN
     single_scattered_radiance_dn_tl_with_ad(n_stokes, n_derivs, n_layers, F_0, n_ulevels, ulevels, utaus, umus, n_umus, omega, omega_l, ltau, ltau_l, btran, btran_l, as_0, as_0_l, atran, atran_l, P, P_l, I_in, I_in_l, I_ss, I_ss_l, utau_output, derivs_layers, derivs_beam, save_tree, work);
#endif
}



/*******************************************************************************
 *
 ******************************************************************************/
void n_t_tms_correction_up(int n_stokes, int n_derivs, int n_layers, double F_0, int n_ulevels, int *ulevels, double *utaus, double *umus, int n_umus, double *omega, double **omega_l, double *omega2, double **omega2_l, double *ltau, double **ltau_l, double *btran, double **btran_l, double *as_0, double **as_0_l, double *atran, double **atran_l, double **P_trun, double **P_full, double ***P_trun_l, double ***P_full_l, double *I_in, double **I_in_l, double **I_c, double ***I_c_l, int utau_output, uchar **derivs_layers, uchar **derivs_beam, save_tree_data save_tree, work_data work) {

     int i;
     int j;

     double **a;

     double ***a_l;

     if (save_tree.t)
          save_tree_encode_s(&save_tree, "n_t_tms_correction_up");

     a = get_work_d2(&work, n_ulevels, n_umus * n_stokes);

     if (n_derivs > 0)
          a_l = get_work_d3(&work, n_ulevels, n_derivs, n_umus * n_stokes);

     if (save_tree.t) save_tree_encode_s(&save_tree, "trun");
     single_scattered_radiance_up(n_stokes, n_derivs, n_layers, F_0, n_ulevels, ulevels, utaus, umus, n_umus, omega,  omega_l,  ltau, ltau_l, btran, btran_l, as_0, as_0_l, atran, atran_l, P_trun, P_trun_l, I_in, I_in_l, a,   a_l, utau_output, derivs_layers, derivs_beam, save_tree, work);

     if (save_tree.t) save_tree_recode_s(&save_tree, "_full");
     single_scattered_radiance_up(n_stokes, n_derivs, n_layers, F_0, n_ulevels, ulevels, utaus, umus, n_umus, omega2, omega2_l, ltau, ltau_l, btran, btran_l, as_0, as_0_l, atran, atran_l, P_full, P_full_l, I_in, I_in_l, I_c, I_c_l, utau_output, derivs_layers, derivs_beam, save_tree, work);

     for (i = 0; i < n_ulevels; ++i) {
          dvec_sub(I_c[i], a[i], I_c[i], n_umus * n_stokes);
          for (j = 0; j < n_derivs; ++j)
               dvec_sub(I_c_l[i][j], a_l[i][j], I_c_l[i][j], n_umus * n_stokes);
     }
}



/*******************************************************************************
 *
 ******************************************************************************/
void n_t_tms_correction_dn(int n_stokes, int n_derivs, int n_layers, double F_0, int n_ulevels, int *ulevels, double *utaus, double *umus, int n_umus, double *omega, double **omega_l, double *omega2, double **omega2_l, double *ltau, double **ltau_l, double *btran, double **btran_l, double *as_0, double **as_0_l, double *atran, double **atran_l, double **P_trun, double **P_full, double ***P_trun_l, double ***P_full_l, double *I_in, double **I_in_l, double **I_c, double ***I_c_l, int utau_output, uchar **derivs_layers, uchar **derivs_beam, save_tree_data save_tree, work_data work) {

     int i;
     int j;

     double **a;

     double ***a_l;

     if (save_tree.t)
          save_tree_encode_s(&save_tree, "n_t_tms_correction_dn");

     a = get_work_d2(&work, n_ulevels, n_umus * n_stokes);

     if (n_derivs > 0)
          a_l = get_work_d3(&work, n_ulevels, n_derivs, n_umus * n_stokes);

     if (save_tree.t) save_tree_encode_s(&save_tree, "trun");
     single_scattered_radiance_dn(n_stokes, n_derivs, n_layers, F_0, n_ulevels, ulevels, utaus, umus, n_umus, omega,  omega_l,  ltau, ltau_l, btran, btran_l, as_0, as_0_l, atran, atran_l, P_trun, P_trun_l, I_in, I_in_l, a,   a_l, utau_output, derivs_layers, derivs_beam, save_tree, work);

     if (save_tree.t) save_tree_recode_s(&save_tree, "_full");
     single_scattered_radiance_dn(n_stokes, n_derivs, n_layers, F_0, n_ulevels, ulevels, utaus, umus, n_umus, omega2, omega2_l, ltau, ltau_l, btran, btran_l, as_0, as_0_l, atran, atran_l, P_full, P_full_l, I_in, I_in_l, I_c, I_c_l, utau_output, derivs_layers, derivs_beam, save_tree, work);

     for (i = 0; i < n_ulevels; ++i) {
          dvec_sub(I_c[i], a[i], I_c[i], n_umus * n_stokes);
          for (j = 0; j < n_derivs; ++j)
               dvec_sub(I_c_l[i][j], a_l[i][j], I_c_l[i][j], n_umus * n_stokes);
     }
}
