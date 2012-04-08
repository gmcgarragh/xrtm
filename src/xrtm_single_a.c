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
#include "xrtm_save_tree.h"
#include "xrtm_single_a.h"
#include "xrtm_utility.h"


/*******************************************************************************
 *
 ******************************************************************************/
static void forward_save_single_scattered_radiance_free(forward_save_single_scattered_radiance_data *d);

int forward_save_single_scattered_radiance_alloc(forward_save_single_scattered_radiance_data *d, int n_layers, int n_stokes) {

     d->free = (void (*)(void *)) forward_save_single_scattered_radiance_free;

     d->I = alloc_array2_d(n_layers, n_stokes);

     return 0;
}



static void forward_save_single_scattered_radiance_free(forward_save_single_scattered_radiance_data *d) {

     free_array2_d(d->I);
}



/*******************************************************************************
 *
 ******************************************************************************/
static void ssr_up_layer_a(int i_layer, int n_stokes, double utau, double *umus, int n_umus, double *omega, double *omega_a, double *ltau, double *ltau_a, double *btran, double *btran_a, double *as_0, double *as_0_a, double **P, double **P_a, double *I, double *I_a, double *I_ss, double *I_ss_a, work_data work) {

     int i;
     int ii;
     int j;
     int jj;

     double a;
/*
     double b;
*/
     double c;
     double d;
     double e;

     double t;

     double utau_a;

     double a_save;
     double a_save_a;

     double b_save_a;

     a = exp(-utau * as_0[i_layer]);
/*
     b = btran[i_layer] * a * omega[i_layer];
*/
     c = ltau[i_layer] - utau;
     d = exp(-c * as_0[i_layer]);

     utau_a = 0.;

     for (i = 0; i < n_umus; ++i) {
          ii = i * n_stokes;

          e = exp(-c / umus[i]);

          a_save = (1. - e * d) / (1. / umus[i] + as_0[i_layer]);

          for (j = 0; j < n_stokes; ++j) {
               jj = ii + j;

               I_a[jj] += I_ss_a[jj] * e;

               ltau_a[i_layer] -= I[jj] * I_ss_a[jj] / umus[i] * e;

               utau_a += I[jj] * I_ss_a[jj] / umus[i] * e;

               b_save_a = I_ss_a[jj];

               btran_a[i_layer] += b_save_a * a * omega[i_layer] * P[i_layer][jj] * a_save;

               utau_a -= b_save_a * btran[i_layer] * as_0[i_layer] * a * omega[i_layer] * P[i_layer][jj] * a_save;

               as_0_a[i_layer] -= b_save_a * btran[i_layer] * utau * a * omega[i_layer] * P[i_layer][jj] * a_save;

               omega_a[i_layer] += b_save_a * btran[i_layer] * a * P[i_layer][jj] * a_save;

               P_a[i_layer][jj] += b_save_a * btran[i_layer] * a * omega[i_layer] * a_save;

               a_save_a = b_save_a * btran[i_layer] * a * omega[i_layer] * P[i_layer][jj];

               t = a_save_a / (1. / umus[i] + as_0[i_layer]);

               ltau_a[i_layer] += t / umus[i] * e * d;

               utau_a -= t / umus[i] * e * d;

               ltau_a[i_layer] += t * e * as_0[i_layer] * d;

               utau_a -= t * e * as_0[i_layer] * d;

               as_0_a[i_layer] += t * e * c * d;

               as_0_a[i_layer] -= t * a_save;

               I_ss_a[jj] = 0.;
          }
     }

     ltau_a[i_layer] += utau_a * utau / ltau[i_layer];
}



/*******************************************************************************
 *
 ******************************************************************************/
void single_scattered_radiance_up_a(int n_stokes, int n_layers, double F_0, int n_ulevels, int *ulevels, double *utaus, double *umus, int n_umus, double *omega, double *omega_a, double *ltau, double *ltau_a, double *btran, double *btran_a, double *as_0, double *as_0_a, double *atran, double *atran_a, double **P, double **P_a, double *I_in, double *I_in_a, double **I_ss, double **I_ss_a, int utau_output, save_tree_data save_tree, work_data work) {

     int i;
     int ii;
     int i1;
     int i2;
     int j;
     int jj;
     int k;

     int n_umus_v;

     int i_ulevel;

     double a;
     double b;

     double *I_a;
     double *I_a2;

     forward_save_single_scattered_radiance_data *save;


     n_umus_v = n_umus * n_stokes;


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     I_a  = get_work_d1(&work, n_umus_v);
     I_a2 = get_work_d1(&work, n_umus_v);

     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     save_tree_encode_s(&save_tree, "single_scattered_radiance_up");

     save_tree_retrieve_data(&save_tree, forward_save_single_scattered_radiance_data, &save);


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     a = F_0 / (4. * PI);


     for (i = 0; i < n_ulevels; ++i) {
          for (j = 0; j < n_umus; ++j) {
               ii = j * n_stokes;

               b = 1. / umus[j] * a;

               for (k = 0; k < n_stokes; ++k)
                    I_ss_a[i][ii + k] *= b;
          }
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     i1       = n_layers - 1;
     i2       = ulevels[0];
     i_ulevel = 0;


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     init_array1_d(I_a, n_umus_v, 0.);

     for (i = i2; i <= i1; ++i) {
          if (! utau_output && i_ulevel < n_ulevels && i == ulevels[i_ulevel]) {
               for (j = 0; j < n_umus_v; ++j)
                    I_a[j] += I_ss_a[i_ulevel][j];

               i_ulevel++;
          }

          init_array1_d(I_a2, n_umus_v, 0.);

          ssr_up_layer_a(i, n_stokes, 0., umus, n_umus, omega, omega_a, ltau, ltau_a, btran, btran_a, as_0, as_0_a, P, P_a, save->I[i], I_a2, NULL, I_a, work);

          dvec_copy(I_a, I_a2, n_umus_v);

          if (  utau_output && i_ulevel < n_ulevels && i == ulevels[i_ulevel]) {
               while (i_ulevel < n_ulevels && i == ulevels[i_ulevel]) {
                    ssr_up_layer_a(i, n_stokes, utaus[i_ulevel], umus, n_umus, omega, omega_a, ltau, ltau_a, btran, btran_a, as_0, as_0_a, P, P_a, save->I[i], I_a, NULL, I_ss_a[i_ulevel], work);

                    i_ulevel++;
               }
          }
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if ((! utau_output && i_ulevel < n_ulevels && i1 + 1 == ulevels[i_ulevel]) || (utau_output && i_ulevel < n_ulevels && i1 + 1 == ulevels[i_ulevel] && utaus[i_ulevel] == 0.)) {
          for (i = 0; i < n_umus_v; ++i)
               I_a[i] += I_ss_a[i_ulevel][i];

          i_ulevel++;
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     for (i = 0; i < n_umus; ++i) {
          ii = i * n_stokes;

          b = 1. / umus[i] * a;

          for (j = 0; j < n_stokes; ++j) {
               jj = ii + j;
               I_in_a[jj] += I_a[jj] / b;
          }
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     init_array2_d(I_ss_a, n_ulevels, n_umus_v, 0.);
}



/*******************************************************************************
 *
 ******************************************************************************/
static void ssr_dn_layer_a(int i_layer, int n_stokes, double utau, double *umus, int n_umus, double *omega, double *omega_a, double *ltau, double *ltau_a, double *btran, double *btran_a, double *as_0, double *as_0_a, double **P, double **P_a, double *I, double *I_a, double *I_ss, double *I_ss_a, work_data work) {

     int i;
     int ii;
     int j;
     int jj;
/*
     double a;
*/
     double b;
     double c;

     double t;

     double utau_a;

     double c_save;
     double c_save_a;

     double d_save_a;
/*
     a = btran[i_layer] * omega[i_layer];
*/
     b = exp(-utau * as_0[i_layer]);

     utau_a = 0.;

     for (i = 0; i < n_umus; ++i) {
          ii = i * n_stokes;

          c = exp(-utau / umus[i]);

          c_save = (c - b) / (as_0[i_layer] - 1. / umus[i]);

          for (j = 0; j < n_stokes; ++j) {
               jj = ii + j;

               I_a[jj] += I_ss_a[jj] * c;

               utau_a -= I[jj] * I_ss_a[jj] / umus[i] * c;
/*
               ltau_a[i_layer] += I[jj] * I_ss_a[jj] / umus[i] * c;
*/
               d_save_a = I_ss_a[jj];

               btran_a[i_layer] += d_save_a * omega[i_layer] * P[i_layer][jj] * c_save;

               omega_a[i_layer] += d_save_a * btran[i_layer] * P[i_layer][jj] * c_save;

               P_a[i_layer][jj] += d_save_a * btran[i_layer] * omega[i_layer] * c_save;

               c_save_a = d_save_a * btran[i_layer] * omega[i_layer] * P[i_layer][jj];

               t = c_save_a / (as_0[i_layer] - 1. / umus[i]);

               utau_a -= t / umus[i] * c;
/*
               ltau_a[i_layer] += t / umus[i] * c;
*/
               utau_a += t * as_0[i_layer] * b;
/*
               ltau_a[i_layer] -= t * as_0[i_layer] * b;
*/
               as_0_a[i_layer] += t * utau * b;

               as_0_a[i_layer] -= t * c_save;

               I_ss_a[jj] = 0.;
          }

     }

     ltau_a[i_layer] += utau_a * utau / ltau[i_layer];
}



/*******************************************************************************
 *
 ******************************************************************************/
void single_scattered_radiance_dn_a(int n_stokes, int n_layers, double F_0, int n_ulevels, int *ulevels, double *utaus, double *umus, int n_umus, double *omega, double *omega_a, double *ltau, double *ltau_a, double *btran, double *btran_a, double *as_0, double *as_0_a, double *atran, double *atran_a, double **P, double **P_a, double *I_in, double *I_in_a, double **I_ss, double **I_ss_a, int utau_output, save_tree_data save_tree, work_data work) {

     int i;
     int ii;
     int i1;
     int i2;
     int j;
     int jj;
     int k;

     int n_umus_v;

     int i_ulevel;

     double a;
     double b;

     double *I_a;
     double *I_a2;

     forward_save_single_scattered_radiance_data *save;


     n_umus_v = n_umus * n_stokes;


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     I_a  = get_work_d1(&work, n_umus * n_stokes);
     I_a2 = get_work_d1(&work, n_umus_v);


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     save_tree_encode_s(&save_tree, "single_scattered_radiance_dn");

     save_tree_retrieve_data(&save_tree, forward_save_single_scattered_radiance_data, &save);


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     a = F_0 / (4. * PI);


     for (i = 0; i < n_ulevels; ++i) {
          for (j = 0; j < n_umus; ++j) {
               ii = j * n_stokes;

               b = 1. / umus[j] * a;

               for (k = 0; k < n_stokes; ++k)
                    I_ss_a[i][ii + k] *= b;
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
     i_ulevel = n_ulevels - 1;


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     init_array1_d(I_a, n_umus_v, 0.);

     for (i = i2 - 1; i >= i1; --i) {
          if (! utau_output && i_ulevel >= 0 && i + 1 == ulevels[i_ulevel]) {
               for (j = 0; j < n_umus_v; ++j)
                    I_a[j] += I_ss_a[i_ulevel][j];

               i_ulevel--;
          }

          init_array1_d(I_a2, n_umus_v, 0.);

          ssr_dn_layer_a(i, n_stokes, ltau[i], umus, n_umus, omega, omega_a, ltau, ltau_a, btran, btran_a, as_0, as_0_a, P, P_a, save->I[i], I_a2, NULL, I_a, work);

          dvec_copy(I_a, I_a2, n_umus_v);

          if (  utau_output && i_ulevel < n_ulevels && i == ulevels[i_ulevel]) {
               while (i_ulevel >= 0 && i == ulevels[i_ulevel]) {
                    ssr_dn_layer_a(i, n_stokes, utaus[i_ulevel], umus, n_umus, omega, omega_a, ltau, ltau_a, btran, btran_a, as_0, as_0_a, P, P_a, save->I[i], I_a, NULL, I_ss_a[i_ulevel], work);

                    i_ulevel--;
               }
          }
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if ((! utau_output && i_ulevel >= 0 && i1 == ulevels[i_ulevel]) || (utau_output && i_ulevel >= 0 && i1 == ulevels[i_ulevel] && utaus[i_ulevel] == 0.)) {
          for (i = 0; i < n_umus_v; ++i)
               I_a[i] += I_ss_a[i_ulevel][i];

          i_ulevel--;
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     for (i = 0; i < n_umus; ++i) {
          ii = i * n_stokes;

          b = 1. / umus[i] * a;

          for (j = 0; j < n_stokes; ++j) {
               jj = ii + j;
               I_in_a[jj] = I_a[jj] / b;
          }
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     init_array2_d(I_ss_a, n_ulevels, n_umus_v, 0.);
}



/*******************************************************************************
 *
 ******************************************************************************/
void n_t_tms_correction_up_a(int n_stokes, int n_layers, double F_0, int n_ulevels, int *ulevels, double *utaus, double *umus, int n_umus, double *omega, double *omega_a, double *omega2, double *omega2_a, double *ltau, double *ltau_a, double *btran, double *btran_a, double *as_0, double *as_0_a, double *atran, double *atran_a, double **P_trun, double **P_full, double **P_trun_a, double **P_full_a, double *I_in, double *I_in_a, double **I_c, double **I_c_a, int utau_output, save_tree_data save_tree, work_data work) {

     int i;

     double **a_a;

     save_tree_encode_s(&save_tree, "n_t_tms_correction_up");

     a_a = get_work_d2(&work, n_ulevels, n_umus * n_stokes);

     for (i = 0; i < n_ulevels; ++i)
          dvec_scale(-1., I_c_a[i], a_a[i], n_umus * n_stokes);

     if (save_tree.t) save_tree_encode_s(&save_tree, "trun");
     single_scattered_radiance_up_a(n_stokes, n_layers, F_0, n_ulevels, ulevels, utaus, umus, n_umus, omega,  omega_a,  ltau, ltau_a, btran, btran_a, as_0, as_0_a, atran, atran_a, P_trun, P_trun_a, I_in, I_in_a, NULL, a_a,   utau_output, save_tree, work);

     if (save_tree.t) save_tree_recode_s(&save_tree, "_full");
     single_scattered_radiance_up_a(n_stokes, n_layers, F_0, n_ulevels, ulevels, utaus, umus, n_umus, omega2, omega2_a, ltau, ltau_a, btran, btran_a, as_0, as_0_a, atran, atran_a, P_full, P_full_a, I_in, I_in_a, NULL, I_c_a, utau_output, save_tree, work);
}



/*******************************************************************************
 *
 ******************************************************************************/
void n_t_tms_correction_dn_a(int n_stokes, int n_layers, double F_0, int n_ulevels, int *ulevels, double *utaus, double *umus, int n_umus, double *omega, double *omega_a, double *omega2, double *omega2_a, double *ltau, double *ltau_a, double *btran, double *btran_a, double *as_0, double *as_0_a, double *atran, double *atran_a, double **P_trun, double **P_full, double **P_trun_a, double **P_full_a, double *I_in, double *I_in_a, double **I_c, double **I_c_a, int utau_output, save_tree_data save_tree, work_data work) {

     int i;

     double **a_a;

     save_tree_encode_s(&save_tree, "n_t_tms_correction_dn");

     a_a = get_work_d2(&work, n_ulevels, n_umus * n_stokes);

     for (i = 0; i < n_ulevels; ++i)
          dvec_scale(-1., I_c_a[i], a_a[i], n_umus * n_stokes);

     if (save_tree.t) save_tree_encode_s(&save_tree, "trun");
     single_scattered_radiance_dn_a(n_stokes, n_layers, F_0, n_ulevels, ulevels, utaus, umus, n_umus, omega,  omega_a,  ltau, ltau_a, btran, btran_a, as_0, as_0_a, atran, atran_a, P_trun, P_trun_a, I_in, I_in_a, NULL, a_a,   utau_output, save_tree, work);

     if (save_tree.t) save_tree_recode_s(&save_tree, "_full");
     single_scattered_radiance_dn_a(n_stokes, n_layers, F_0, n_ulevels, ulevels, utaus, umus, n_umus, omega2, omega2_a, ltau, ltau_a, btran, btran_a, as_0, as_0_a, atran, atran_a, P_full, P_full_a, I_in, I_in_a, NULL, I_c_a, utau_output, save_tree, work);
}



/*******************************************************************************
 *
 ******************************************************************************/
#ifdef INCLUDE_DEV_SOURCE
#include "../src2/xrtm_single_a2.c"
#endif
