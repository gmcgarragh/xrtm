/*******************************************************************************
**
**    Copyright (C) 2007-2020 Greg McGarragh <greg.mcgarragh@colostate.edu>
**
**    This source code is licensed under the GNU General Public License (GPL),
**    Version 3.  See the file COPYING for more details.
**
*******************************************************************************/

#ifndef RTUTIL_SCAT_FFT_H
#define RTUTIL_SCAT_FFT_H

#include <fftw3.h>

#ifdef __cplusplus
extern "C" {
#endif


typedef struct {
     fftw_complex *in;
     fftw_complex *out;
     fftw_plan plan;
} fft_data;


int scat_matrix_n_phis(int n_four, int n_coef);

void fft_data_alloc(fft_data *fft, int n_phis);
void fft_data_free(fft_data *fft);

void azimuth_to_fourier(int n_four, int n_stokes1, int n_stokes2, int n_phis, double ***F, double *P, int accel1, int accel2, fft_data *fft);


#ifdef __cplusplus
}
#endif

#endif /* RTUTIL_SCAT_FFT_H */
