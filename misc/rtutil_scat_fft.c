/*******************************************************************************
**
**    Copyright (C) 2007-2020 Greg McGarragh <greg.mcgarragh@colostate.edu>
**
**    This source code is licensed under the GNU General Public License (GPL),
**    Version 3.  See the file COPYING for more details.
**
*******************************************************************************/

int rtutil_scat_fft_dummy(void);

int rtutil_scat_fft_dummy(void) { return 0; }


#ifdef FFTW_FLAG

#include <gutil.h>

#include <fftw3.h>

#include "rtutil_scat_fft.h"


/*******************************************************************************
 *
 ******************************************************************************/
int scat_matrix_n_phis(int n_four, int n_coef) {

     if (n_four - 1 == 0)
          return 2 * (int) ((n_coef - 1 + 1) / 2) + 4;
     else
          return 2 * (int) pow(2., (int) (log(n_coef - 1 + 4.) / log(2.) + 1.));
}



/*******************************************************************************
 *
 ******************************************************************************/
void fft_data_alloc(fft_data *fft, int n_phis) {

     fft->in   = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * n_phis);
     fft->out  = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * n_phis);

     fft->plan = fftw_plan_dft_1d(n_phis, fft->in, fft->out, FFTW_BACKWARD, FFTW_ESTIMATE);
}



void fft_data_free(fft_data *fft) {

     fftw_free(fft->in);
     fftw_free(fft->out);

     fftw_destroy_plan(fft->plan);
}



void azimuth_to_fourier(int n_four, int n_stokes1, int n_stokes2, int n_phis, double ***F, double *P, int accel1, int accel2, fft_data *fft) {

     int i;
     int j;
     int k;

     int i_offset;
     int j_offset;

     double f;

     int r_or_i[4][4] = { {0, 0, 1, 1},
                          {0, 0, 1, 1},
                          {1, 1, 0, 0},
                          {1, 1, 0, 0} };

     double p_or_n[4][4] = { { 1.,  1., -1., -1.},
                             { 1.,  1., -1., -1.},
                             { 1.,  1.,  1.,  1.},
                             { 1.,  1.,  1.,  1.} };

     f = 2. / n_phis;

     for (i = 0; i < n_stokes1; ++i) {

          i_offset = i * accel2;

          for (j = 0; j < n_stokes2; ++j) {

               j_offset = i_offset + j;

               for (k = 0; k < n_phis; ++k) {
                    fft->in[k][0] = F[k][i][j];
                    fft->in[k][1] = 0.;
               }

               fftw_execute(fft->plan);

               P[j_offset] = p_or_n[i][j] * fft->out[0][r_or_i[i][j]] / n_phis;

               for (k = 1; k < n_four; ++k)
                    P[k * accel1 + j_offset] = p_or_n[i][j] * fft->out[k][r_or_i[i][j]] * f;
          }
     }
}

#endif
