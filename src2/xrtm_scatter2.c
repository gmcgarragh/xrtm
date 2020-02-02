/*******************************************************************************
**
**    Copyright (C) 2007-2020 Greg McGarragh <greg.mcgarragh@colostate.edu>
**
**    This source code is licensed under the GNU General Public License (GPL),
**    Version 3.  See the file COPYING for more details.
**
*******************************************************************************/

#ifndef CPROTO
#include <rtutil_scat_fft.h>
#endif

/*******************************************************************************
 *
 ******************************************************************************/
void scat_vector_rotate_evans(int n_stokes, double mu1, double mu2, double MU, double d_phi, double *P1, double *P2) {

     double sin_theta1;
     double sin_theta2;

     double sin_d_phi;
     double cos_d_phi;

     double sin_scat;

     double P1_01;

     double sin1;
     double cos1;

     double sin21;
     double cos21;

     P2[0] = P1[0];
     if (n_stokes == 1)
          return;

     P1_01 = P1[1];

     if (fabs(MU) == 1.) {
          if (n_stokes >= 2)
               P2[1] =  P1_01;
          if (n_stokes >= 3)
               P2[2] =  0.;
          if (n_stokes >= 4)
               P2[3] =  0.;
     }
     else {
          sin_theta1 = sqrt(1. - mu1*mu1);
          sin_theta2 = sqrt(1. - mu2*mu2);

          sin_scat   = sqrt(1. - MU *MU );

          sin_d_phi  = sin(d_phi);
          cos_d_phi  = cos(d_phi);

          sin1 =                   sin_theta2    *sin_d_phi  / sin_scat;
          cos1 = (sin_theta1*mu2 - sin_theta2*mu1*cos_d_phi) / sin_scat;

          sin21 =      2.*sin1*cos1;
          cos21 = 1. - 2.*sin1*sin1;

          if (n_stokes >= 2)
               P2[1] =  P1_01*cos21;
          if (n_stokes >= 3)
               P2[2] = -P1_01*sin21;
          if (n_stokes >= 4)
               P2[3] =  0.0;
     }
}



void scat_matrix_rotate_evans(int n_stokes, double mu1, double mu2, double MU, double d_phi, double **P1, double **P2) {

     double sin_theta1;
     double sin_theta2;

     double sin_d_phi;
     double cos_d_phi;

     double sin_scat;

     double P1_01;
     double P1_11;
     double P1_22;

     double P1_23;
     double P1_33;

     double sin1;
     double sin2;
     double cos1;
     double cos2;

     double sin21;
     double cos21;
     double sin22;
     double cos22;

     double sin21_sin22;
     double cos21_cos22;
     double sin21_cos22;
     double cos21_sin22;

     P2[0][0] = P1[0][0];
     if (n_stokes == 1)
          return;

     P1_01 = P1[0][1];
     P1_11 = P1[1][1];
     P1_22 = P1[2][2];

     if (n_stokes >= 4) {
        P1_23 = P1[2][3];
        P1_33 = P1[3][3];
     }

     if (fabs(MU) == 1.) {
          if (n_stokes >= 2) {
               P2[0][1] =  P1_01;
               P2[1][0] =  P1_01;
               P2[1][1] =  P1_11;
          }
          if (n_stokes >= 3) {
               P2[0][2] =  0.;
               P2[1][2] =  0.;
               P2[2][0] =  0.;
               P2[2][1] =  0.;
               P2[2][2] =  P1_22;
          }
          if (n_stokes >= 4) {
               P2[0][3] =  0.;
               P2[1][3] =  0.;
               P2[2][3] =  P1_23;
               P2[3][0] =  0.;
               P2[3][1] =  0.;
               P2[3][2] = -P1_23;
               P2[3][3] =  P1_33;
          }
     }
     else {
          sin_theta1 = sqrt(1. - mu1*mu1);
          sin_theta2 = sqrt(1. - mu2*mu2);

          sin_scat   = sqrt(1. - MU *MU );

          sin_d_phi  = sin(d_phi);
          cos_d_phi  = cos(d_phi);
/*
if (sin_scat == 0.) {
          sin1 =  0.;
          sin2 =  0.;
          cos1 =  1.;
          cos2 = -1.;
}
else {
*/
          sin1 =                   sin_theta2    *sin_d_phi  / sin_scat;
          sin2 =                   sin_theta1    *sin_d_phi  / sin_scat;
          cos1 = (sin_theta1*mu2 - sin_theta2*mu1*cos_d_phi) / sin_scat;
          cos2 = (sin_theta2*mu1 - sin_theta1*mu2*cos_d_phi) / sin_scat;
/*
}
*/
          sin21 =      2.*sin1*cos1;
          cos21 = 1. - 2.*sin1*sin1;
          sin22 =      2.*sin2*cos2;
          cos22 = 1. - 2.*sin2*sin2;

          sin21_sin22 = sin21*sin22;
          cos21_cos22 = cos21*cos22;

          if (n_stokes >= 2) {
               P2[0][1] =  P1_01*cos21;
               P2[1][0] =  P1_01*cos22;
               P2[1][1] =  P1_11*cos21_cos22 - P1_22*sin21_sin22;
          }
          if (n_stokes >= 3) {
               sin21_cos22 = sin21*cos22;
               cos21_sin22 = cos21*sin22;

               P2[0][2] = -P1_01*sin21;
               P2[1][2] = -P1_11*sin21_cos22 - P1_22*cos21_sin22;
               P2[2][0] =  P1_01*sin22;
               P2[2][1] =  P1_11*cos21_sin22 + P1_22*sin21_cos22;
               P2[2][2] = -P1_11*sin21_sin22 + P1_22*cos21_cos22;
          }
          if (n_stokes >= 4) {
               P2[0][3] =  0.0;
               P2[1][3] = -P1_23*sin22;
               P2[2][3] =  P1_23*cos22;
               P2[3][0] =  0.0;
               P2[3][1] = -P1_23*sin21;
               P2[3][2] = -P1_23*cos21;
               P2[3][3] =  P1_33;
          }
     }
}



void phase_matrix_symmetry(int n_stokes, double **P1, double **P2) {

     P2[0][0] = P1[0][0];
     if (n_stokes >= 2) {
          P2[0][1] =  P1[0][1];
          P2[1][0] =  P1[1][0];
          P2[1][1] =  P1[1][1];
     }
     if (n_stokes >= 3) {
          P2[0][2] = -P1[0][2];
          P2[1][2] = -P1[1][2];
          P2[2][0] = -P1[2][0];
          P2[2][1] = -P1[2][1];
          P2[2][2] =  P1[2][2];
     }
     if (n_stokes >= 4) {
          P2[0][3] = -P1[0][3];
          P2[1][3] = -P1[1][3];
          P2[2][3] =  P1[2][3];
          P2[3][0] = -P1[3][0];
          P2[3][1] = -P1[3][1];
          P2[3][2] =  P1[3][2];
          P2[3][3] =  P1[3][3];
     }
}



static void build_phase_pair_vector_lc(int n_four, int n_coef, int n_mus1, double *mu1_array, double mu1_units, int n_stokes1, int n_mus2, double *mu2_array, double mu2_units, int n_stokes2, int azi_sym, double **coefs, double *P, work_data work, int accel) {

     int i;
     int j;
     int l;

     int flag;

     int n_phis;
     int h_phis;

     int accel1;
     int accel2;

     int i_offset;
     int j_offset;

     double a;
     double b;
     double c;

     double twopi = 2.*PI;

     double mu1;
     double mu2;
     double phi;

     double d_phi;

     double MU;

     double *p;

     double ***F;

     fft_data fft;


     n_phis = scat_matrix_n_phis(n_four, n_coef);

     if (azi_sym)
          h_phis = n_phis / 2 + 1;
     else
          h_phis = n_phis;

     d_phi = twopi / n_phis;


     p = alloc_array1_d(n_coef);
     F = alloc_array3_d(n_phis + 1, 4, 4);

     fft_data_alloc(&fft, n_phis);


     init_array3_d(F, n_phis + 1, 4, 4, 0.);


     accel1 = accel;
     accel2 = n_stokes2 * n_mus2;

     flag = n_stokes1 > 1 || n_stokes2 > 1;

     for (i = 0; i < n_mus1; ++i) {
          i_offset = i * n_stokes1;

          mu1 = mu1_units * mu1_array[i];

          a = 1. - mu1 * mu1;

          for (j = 0; j < n_mus2; ++j) {
               j_offset = j * n_stokes2;

               mu2 = mu2_units * mu2_array[j];

               b = mu1 * mu2;
               c = sqrt(a * (1. - mu2 * mu2));

               for (l = 0; l < h_phis; ++l) {

                    phi = l * d_phi;

                    MU = b + c * cos(phi);

                    leg_poly(n_coef, MU, p);

                    build_scat_matrix_lc(n_coef, p, coefs, F[l], flag);

                    scat_matrix_rotate_evans(n_stokes1, mu2, mu1, MU, phi, F[l], F[l]);

                    if (azi_sym)
                         phase_matrix_symmetry(n_stokes1, F[l], F[n_phis-l-1 + 1]);
               }

               azimuth_to_fourier(n_four, n_stokes1, n_stokes2, n_phis, F, &P[i_offset * accel2 + j_offset], accel1, accel2, &fft);
          }
     }

     free_array1_d(p);
     free_array3_d(F);

     fft_data_free(&fft);
}



void build_phase_vecs_vector_lc(int n_four, int n_coef, int n_mus, int n_stokes, double *mu, double mu_0, double **coefs,  double **P_pp,  double **P_mp,  work_data work, int accel) {

     build_phase_pair_vector_lc(n_four, n_coef, n_mus, mu, 1., n_stokes, 1, &mu_0,  1., 1, 0, coefs, *P_pp, work, accel);

     build_phase_pair_vector_lc(n_four, n_coef, n_mus, mu, -1., n_stokes, 1, &mu_0,  1., 1, 0, coefs, *P_mp, work, accel);
}



void build_phase_mats_vector_lc(int n_four, int n_coef, int n_mus1, int n_mus2, int n_stokes, double *mu1, double *mu2, double **coefs, double ***P_pp, double ***P_mp, double ***P_mm, double ***P_pm, int vector, work_data work, int accel) {

     int i;
     int ii;

     build_phase_pair_vector_lc(n_four, n_coef, n_mus1, mu1,  1., n_stokes, n_mus2, mu2, 1., n_stokes, 1, coefs, **P_pp, work, accel);

     build_phase_pair_vector_lc(n_four, n_coef, n_mus1, mu1, -1., n_stokes, n_mus2, mu2, 1., n_stokes, 1, coefs, **P_mp, work, accel);

     if (vector) {
          for (i = 0; i < n_four; ++i) {
               ii = i * accel;

               phase_matrix_symmetry_ldx3(n_mus1, n_stokes, n_mus2, n_stokes, (**P_pp)+ii, (**P_mp)+ii, (**P_mm)+ii, (**P_pm)+ii, n_mus2*n_stokes, 1.);
          }
     }
}
