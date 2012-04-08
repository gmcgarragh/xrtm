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
#include "xrtm_brdf.h"
#include "xrtm_brdf_a.h"
#include "xrtm_utility.h"


/*******************************************************************************
 *
 ******************************************************************************/
void brdf_build_kernel_vecs_a(int i_offset, int n_quad, int j_offset, int n_stokes, int n_kernels, int n_kernel_quad, enum xrtm_kernel_type *kernels, double *ampfac, double *ampfac_a, double **params, double **params_a, double *kernel_qx, double *kernel_qw, brdf_aux_data *aux, double **kernel_f, double **kernel_f_a, work_data work) {

     int i;
     int j;

     int flag;
/*
     double f_a;
*/
     for (i = 0; i < n_kernels; ++i) {
          flag = brdf_needs_fourier(kernels + i, 1);

          if (! flag) {
               for (j = 0; j < n_quad; ++j) {
                    ampfac_a[i] += kernel_f_a[0][j];
               }
          }
          else {

          }
     }
}



/*******************************************************************************
 *
 ******************************************************************************/
void brdf_build_ref_vec_a(int i_four, int i_offset, int n_quad, int j_offset, int n_stokes, double qf, double *qx_v, double *qw_v, int n_kernels, int n_kernel_quad, enum xrtm_kernel_type *kernels, double *kernel_qx, double *kernel_qw, double **kernel_f_a, double *R_s_a, work_data work) {

     int i;
     int ii;

     int flag;

     double a;

     flag = brdf_needs_fourier(kernels, n_kernels);

     a = (2. - (i_four == 0 ? 1. : 0.));

     for (i = 0, ii = 0; i < n_quad; ++i, ii += n_stokes)
          R_s_a[ii] *= a;

     if (! flag) {
          for (i = 0, ii = 0; i < n_quad; ++i, ii += n_stokes) {
               kernel_f_a[0][i] += R_s_a[ii];
          }
     }
     else {

     }
}



/*******************************************************************************
 *
 ******************************************************************************/
void brdf_build_kernel_mats_a(int i_offset, int n_quad1, int j_offset, int n_quad2, int n_stokes, int n_kernels, int n_kernel_quad, enum xrtm_kernel_type *kernels, double *ampfac, double *ampfac_a, double **params, double **params_a, double *kernel_qx, double *kernel_qw, brdf_aux_data *aux, double ***kernel_f, double ***kernel_f_a, work_data work) {

     int i;
     int j;
     int k;

     int flag;
/*
     double f_a;
*/
     for (i = 0; i < n_kernels; ++i) {
          flag = brdf_needs_fourier(kernels + i, 1);

          if (! flag) {
               for (j = 0; j < n_quad1; ++j) {
                    for (k = 0; k < n_quad2; ++k) {
                         ampfac_a[i] += kernel_f_a[0][j][k];
                    }
               }
          }
          else {

          }
     }
}



/*******************************************************************************
 *
 ******************************************************************************/
void brdf_build_ref_mat_a(int i_four, int i_offset, int n_quad1, int j_offset, int n_quad2, int n_stokes, double qf, double *qx_v, double *qw_v, int n_kernels, int n_kernel_quad, enum xrtm_kernel_type *kernels, double *kernel_qx, double *kernel_qw, double ***kernel_f_a, double **R_s_a, work_data work) {

     int ii;
     int i;
     int j;
     int jj;

     int flag;

     int jj_offset;

     double a;
     double b;

     flag = brdf_needs_fourier(kernels, n_kernels);

     a = (2. - (i_four == 0 ? 1. : 0.)) * (1. + (i_four == 0 ? 1. : 0.)) * qf;

     jj_offset = j_offset * n_stokes;
     for (j = 0, jj = 0; j < n_quad2; ++j, jj += n_stokes) {
          b = a * qx_v[jj_offset + jj] * qw_v[jj_offset + jj];
          for (i = 0, ii = 0; i < n_quad1; ++i, ii += n_stokes) {
               R_s_a[ii][jj] *= b;
          }
     }

     if (! flag) {
          for (i = 0, ii = 0; i < n_quad1; ++i, ii += n_stokes) {
               for (j = 0, jj = 0; j < n_quad2; ++j, jj += n_stokes) {
                    kernel_f_a[0][i][j] += R_s_a[ii][jj];
               }
          }
     }
     else {

     }
}
