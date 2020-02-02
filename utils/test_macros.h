/*******************************************************************************
**
**    Copyright (C) 2007-2020 Greg McGarragh <greg.mcgarragh@colostate.edu>
**
**    This source code is licensed under the GNU General Public License (GPL),
**    Version 3.  See the file COPYING for more details.
**
*******************************************************************************/

#ifndef TEST_MACROS_H
#define TEST_MACROS_H

#include <xrtm_interface.h>

#ifdef __cplusplus
extern "C" {
#endif


/*******************************************************************************
 *
 ******************************************************************************/
#define HANDLE_RETURN(CALL, S) do {							\
     int r;										\
											\
     r = CALL;										\
											\
     if (r < 0) {									\
          fprintf(stderr, "ERROR: %s\n", S);						\
          return -1;									\
     }											\
     else										\
     if (r > 0) {									\
          fprintf(stderr, "FAIL: %s\n", S);						\
          return  1;									\
     }											\
} while (0)


#define SET_SELF_CHECK(A, B, S) do {							\
     if (A != B) {									\
          fprintf (stderr, "FAIL: xrtm_set self check failed: %s\n", S);		\
          return 1;									\
     }											\
} while (0)


#define SET_SELF_CLOSE(A, B, E, S) do {							\
     if (fabs(A - B) > E) {								\
          fprintf (stderr, "FAIL: xrtm_set self check failed: %s\n", S);		\
          return 1;									\
     }											\
} while (0)


#define SET_SELF_CHECK2(A, B, N, S) do {						\
     int i_SET_SELF_CHECK2;								\
											\
     for (i_SET_SELF_CHECK2 = 0; i_SET_SELF_CHECK2 < N; ++i_SET_SELF_CHECK2 ) {		\
          if (A[i_SET_SELF_CHECK2] != B[i_SET_SELF_CHECK2]) {				\
               fprintf (stderr, "FAIL: xrtm_set self check failed: %s\n", S);		\
               return 1;								\
          }										\
     }											\
} while (0)


#define SET_SELF_CLOSE2(A, B, N, E, S) do {						\
     int i_SET_SELF_CLOSE2;								\
											\
     for (i_SET_SELF_CLOSE2 = 0; i_SET_SELF_CLOSE2 < N; ++i_SET_SELF_CLOSE2) {		\
          if (fabs(A[i_SET_SELF_CLOSE2] - B[i_SET_SELF_CLOSE2]) > E) {			\
               fprintf (stderr, "FAIL: xrtm_set self check failed: %s\n", S);		\
               return 1;								\
          }										\
     }											\
} while (0)



/*******************************************************************************
 *
 ******************************************************************************/
#define XRTM_CREATE(D, OPTIONS, SOLVERS, N_COEF, N_QUAD, N_STOKES, N_DERIVS, N_LAYERS, N_THETA_0S, N_KERNELS, N_KERNEL_QUAD, KERNELS, N_OUT_LEVELS, N_MUS) do {	\
     if (xrtm_create(D, OPTIONS, SOLVERS, N_COEF, N_QUAD, N_STOKES, N_DERIVS, N_LAYERS, N_THETA_0S, N_KERNELS, N_KERNEL_QUAD, KERNELS, N_OUT_LEVELS, N_MUS)) {	\
          fprintf(stderr, "ERROR: xrtm_create()\n");														\
          return -1;																	\
     }																			\
} while (0)


#define XRTM_DESTROY(D) do {						\
     xrtm_destroy(D);							\
} while (0)



/*******************************************************************************
 *
 ******************************************************************************/
#define XRTM_SET_DOUB_D_TAU(D, D_TAU) do {				\
     double a;								\
									\
     if (xrtm_set_doub_d_tau(D, D_TAU)) {				\
          fprintf(stderr, "ERROR: xrtm_set_doub_d_tau()\n");		\
          return -1;							\
     }									\
									\
     if ((a = xrtm_get_doub_d_tau(D)) < -1) {				\
          fprintf(stderr, "ERROR: xrtm_get_doub_d_tau()\n");		\
          return -1;							\
     }									\
									\
     SET_SELF_CHECK(D_TAU, a, "doub_d_tau");				\
} while (0)


#define XRTM_SET_PADE_PARAMS(D, S_PADE, R_PADE) do {			\
     int a;								\
     int b;								\
									\
     if (xrtm_set_pade_params(D, S_PADE, R_PADE)) {			\
          fprintf(stderr, "ERROR: xrtm_set_pade_params()\n");		\
          return -1;							\
     }									\
									\
     if ((xrtm_get_pade_params(D, &a, &b)) < 0) {			\
          fprintf(stderr, "ERROR: xrtm_get_pade_params()\n");		\
          return -1;							\
     }									\
									\
     SET_SELF_CHECK(S_PADE, a, "s_pade");				\
     SET_SELF_CHECK(R_PADE, b, "r_pade");				\
} while (0)


#define XRTM_SET_SOS_PARAMS(D, MAX_OS, MAX_TAU, SOS_TOL) do {		\
     int a;								\
     double b;								\
     double c;								\
									\
     if (xrtm_set_sos_params(D, MAX_OS, MAX_TAU, SOS_TOL)) {		\
          fprintf(stderr, "ERROR: xrtm_set_sos_params()\n");		\
          return -1;							\
     }									\
									\
     if ((xrtm_get_sos_params(D, &a, &b, &c)) < 0) {			\
          fprintf(stderr, "ERROR: xrtm_get_sos_params()\n");		\
          return -1;							\
     }									\
									\
     SET_SELF_CHECK(MAX_OS,  a, "sos_max_os");				\
     SET_SELF_CHECK(MAX_TAU, b, "max_tau");				\
     SET_SELF_CHECK(SOS_TOL, c, "sos_tol");				\
} while (0)


#define XRTM_SET_FOURIER_TOL(D, FOURIER_TOL) do {			\
     double a;								\
									\
     if (xrtm_set_fourier_tol(D, FOURIER_TOL)) {			\
          fprintf(stderr, "ERROR: xrtm_set_fourier_tol()\n");		\
          return -1;							\
     }									\
									\
     if ((a = xrtm_get_fourier_tol(D)) < -1) {				\
          fprintf(stderr, "ERROR: xrtm_get_fourier_tol()\n");		\
          return -1;							\
     }									\
									\
     SET_SELF_CHECK(FOURIER_TOL, a, "fourier_tol");			\
} while (0)


#define XRTM_SET_PLANET_R(D, PLANET_R) do {				\
     double a;								\
									\
     if (xrtm_set_planet_r(D, PLANET_R)) {				\
          fprintf(stderr, "ERROR: xrtm_set_planet_r()\n");		\
          return -1;							\
     }									\
									\
     if ((a = xrtm_get_planet_r(D)) < -1) {				\
          fprintf(stderr, "ERROR: xrtm_get_PLANET_R()\n");		\
          return -1;							\
     }									\
									\
     SET_SELF_CHECK(PLANET_R, a, "phi_0");				\
} while (0)


#define XRTM_SET_LEVELS_Z(D, LEVELS_Z) do {				\
     int n;								\
     double *a;								\
									\
     if (xrtm_set_levels_z(D, LEVELS_Z)) {				\
          fprintf(stderr, "ERROR: xrtm_set_levels_z()\n");		\
          return -1;							\
     }									\
									\
     n = xrtm_get_n_layers(D) + 1;					\
									\
     a = alloc_array1_d(n);						\
									\
     if (xrtm_get_levels_z(D, a)) {					\
          fprintf(stderr, "ERROR: xrtm_get_levels_z()\n");		\
          return -1;							\
     }									\
									\
     SET_SELF_CLOSE2((LEVELS_Z), a, n, 1.e-4, "levels_z");		\
									\
     free_array1_d(a);							\
} while (0)


#define XRTM_SET_OUT_LEVELS(D, OUT_LEVELS) do {				\
     int n;								\
     int *a;								\
									\
     if (xrtm_set_out_levels(D, OUT_LEVELS)) {				\
          fprintf(stderr, "ERROR: xrtm_set_out_levels()\n");		\
          return -1;							\
     }									\
									\
     n = xrtm_get_n_out_levels(D);					\
									\
     a = alloc_array1_i(n);						\
									\
     if (xrtm_get_out_levels(D, a)) {					\
          fprintf(stderr, "ERROR: xrtm_get_out_levels()\n");		\
          return -1;							\
     }									\
									\
     SET_SELF_CLOSE2((OUT_LEVELS), a, n, 1.e-4, "out_levels");		\
									\
     free_array1_i(a);							\
} while (0)

#define XRTM_SET_OUT_TAUS(D, OUT_TAUS) do {				\
     int n;								\
     double *a;								\
									\
     if (xrtm_set_out_taus(D, OUT_TAUS)) {				\
          fprintf(stderr, "ERROR: xrtm_set_out_taus()\n");		\
          return -1;							\
     }									\
									\
     n = xrtm_get_n_out_levels(D);					\
									\
     a = alloc_array1_d(n);						\
									\
     if (xrtm_get_out_taus(D, a)) {					\
          fprintf(stderr, "ERROR: xrtm_get_out_taus()\n");		\
          return -1;							\
     }									\
									\
     SET_SELF_CLOSE2((OUT_TAUS), a, n, 1.e-4, "out_taus");		\
									\
     free_array1_d(a);							\
} while (0)


#define XRTM_SET_OUT_THETAS(D, OUT_THETAS) do {				\
     int n;								\
     double *a;								\
									\
     if (xrtm_set_out_thetas(D, OUT_THETAS)) {				\
          fprintf(stderr, "ERROR: xrtm_set_out_thetas()\n");		\
          return -1;							\
     }									\
									\
     n = xrtm_get_n_out_thetas(D);					\
									\
     a = alloc_array1_d(n);						\
									\
     if (xrtm_get_out_thetas(D, a)) {					\
          fprintf(stderr, "ERROR: xrtm_get_out_thetas()\n");		\
          return -1;							\
     }									\
									\
     SET_SELF_CLOSE2((OUT_THETAS), a, n, 1.e-4, "out_thetas");		\
									\
     free_array1_d(a);							\
} while (0)


#define XRTM_SET_F_ISO_TOP(D, F_ISO_TOP) do {				\
     double a;								\
									\
     if (xrtm_set_F_iso_top(D, F_ISO_TOP)) {				\
          fprintf(stderr, "ERROR: xrtm_set_F_iso_top()\n");		\
          return -1;							\
     }									\
									\
     if ((a = xrtm_get_F_iso_top(D)) < -1) {				\
          fprintf(stderr, "ERROR: xrtm_get_F_iso_top()\n");		\
          return -1;							\
     }									\
									\
     SET_SELF_CHECK(F_ISO_TOP, a, "F_iso_top");				\
} while (0)


#define XRTM_SET_F_ISO_BOT(D, F_ISO_BOT) do {				\
     double a;								\
									\
     if (xrtm_set_F_iso_bot(D, F_ISO_BOT)) {				\
          fprintf(stderr, "ERROR: xrtm_set_F_iso_bot()\n");		\
          return -1;							\
     }									\
									\
     if ((a = xrtm_get_F_iso_bot(D)) < -1) {				\
          fprintf(stderr, "ERROR: xrtm_get_F_iso_bot()\n");		\
          return -1;							\
     }									\
									\
     SET_SELF_CHECK(F_ISO_BOT, a, "F_iso_bot");				\
} while (0)


#define XRTM_SET_F_0(D, F_0) do {					\
     double a;								\
									\
     if (xrtm_set_F_0(D, F_0)) {					\
          fprintf(stderr, "ERROR: xrtm_set_F_0()\n");			\
          return -1;							\
     }									\
									\
     if ((a = xrtm_get_F_0(D)) < -1) {					\
          fprintf(stderr, "ERROR: xrtm_get_F_0()\n");			\
          return -1;							\
     }									\
									\
     SET_SELF_CHECK(F_0, a, "F_0");					\
} while (0)


#define XRTM_SET_THETA_0(D, THETA_0) do {				\
     double a;								\
									\
     if (xrtm_set_theta_0(D, THETA_0)) {				\
          fprintf(stderr, "ERROR: xrtm_set_theta_0()\n");		\
          return -1;							\
     }									\
									\
     if ((a = xrtm_get_theta_0(D)) < -1) {				\
          fprintf(stderr, "ERROR: xrtm_get_theta_0()\n");		\
          return -1;							\
     }									\
									\
     SET_SELF_CLOSE(THETA_0, a, (D)->misc_input.threshold_mu_0_singlarity*1.e2, "theta_0");				\
} while (0)


#define XRTM_SET_PHI_0(D, PHI_0) do {					\
     double a;								\
									\
     if (xrtm_set_phi_0(D, PHI_0)) {					\
          fprintf(stderr, "ERROR: xrtm_set_phi_0()\n");			\
          return -1;							\
     }									\
									\
     if ((a = xrtm_get_phi_0(D)) < -1) {				\
          fprintf(stderr, "ERROR: xrtm_get_phi_0()\n");			\
          return -1;							\
     }									\
									\
     SET_SELF_CHECK(PHI_0, a, "phi_0");					\
} while (0)


#define XRTM_SET_LAMBDA(D, LAMBDA) do {					\
     double a;								\
									\
     if (xrtm_set_lambda(D, LAMBDA)) {					\
          fprintf(stderr, "ERROR: xrtm_set_lambda()\n");		\
          return -1;							\
     }									\
									\
     if ((a = xrtm_get_lambda(D)) < -1) {				\
          fprintf(stderr, "ERROR: xrtm_get_lambda()\n");		\
          return -1;							\
     }									\
									\
     SET_SELF_CHECK(LAMBDA, a, "LAMBDA");				\
} while (0)


#define XRTM_SET_LEVELS_B(D, LEVELS_B) do {				\
     int n;								\
     double *a;								\
									\
     if (xrtm_set_levels_b(D, LEVELS_B)) {				\
          fprintf(stderr, "ERROR: xrtm_set_levels_b()\n");		\
          return -1;							\
     }									\
									\
     n = xrtm_get_n_layers(D) + 1;					\
									\
     a = alloc_array1_d(n);						\
									\
     if (xrtm_get_levels_b(D, a)) {					\
          fprintf(stderr, "ERROR: xrtm_get_levels_b()\n");		\
          return -1;							\
     }									\
									\
     SET_SELF_CLOSE2((LEVELS_B), a, n, 1.e-4, "levels_b");		\
									\
     free_array1_d(a);							\
} while (0)

#define XRTM_SET_SURFACE_B(D, SURFACE_B) do {				\
     double a;								\
									\
     if (xrtm_set_surface_b(D, SURFACE_B)) {				\
          fprintf(stderr, "ERROR: xrtm_set_surface_b()\n");		\
          return -1;							\
     }									\
									\
     if ((a = xrtm_get_surface_b(D)) < -1) {				\
          fprintf(stderr, "ERROR: xrtm_get_surface_b()\n");		\
          return -1;							\
     }									\
									\
     SET_SELF_CHECK(SURFACE_B, a, "SURFACE_B");				\
} while (0)


/*******************************************************************************
 *
 ******************************************************************************/
#define XRTM_SET_COEF_1(D, I_LAYER, N_COEF_LAYER, COEF) do {		\
     if (xrtm_set_coef_1(D, I_LAYER, N_COEF_LAYER, COEF)) {		\
          fprintf(stderr, "ERROR: xrtm_set_coef_1()\n");		\
          return -1;							\
     }									\
} while (0)


#define XRTM_SET_COEF_N(D, N_COEF_LAYER, COEF) do {			\
     if (xrtm_set_coef_n(D, N_COEF_LAYER, COEF)) {			\
          fprintf(stderr, "ERROR: xrtm_set_coef_n()\n");		\
          return -1;							\
     }									\
} while (0)


#define XRTM_SET_COEF_L_11(D, I_LAYER, I_DERIV, COEF_L) do {		\
     if (xrtm_set_coef_l_11(D, I_LAYER, I_DERIV, COEF_L)) {		\
          fprintf(stderr, "ERROR: xrtm_set_coef_l_11()\n");		\
          return -1;							\
     }									\
} while (0)


#define XRTM_SET_COEF_L_1N(D, I_LAYER, COEF_L) do {			\
     if (xrtm_set_coef_l_1n(D, I_LAYER, COEF_L)) {			\
          fprintf(stderr, "ERROR: xrtm_set_coef_l_1n()\n");		\
          return -1;							\
     }									\
} while (0)



/*******************************************************************************
 *
 ******************************************************************************/
#define XRTM_SET_OMEGA_1(D, I_LAYER, OMEGA) do {			\
     double a;								\
									\
     if (xrtm_set_omega_1(D, I_LAYER, OMEGA)) {				\
          fprintf(stderr, "ERROR: xrtm_set_omega_1()\n");		\
          return -1;							\
     }									\
									\
     if ((a = xrtm_get_omega(D, I_LAYER)) < -1) {			\
          fprintf(stderr, "ERROR: xrtm_get_omega()\n");			\
          return -1;							\
     }									\
									\
     SET_SELF_CLOSE(OMEGA, a, D->misc_input.threshold_omega_singlarity*1.e0, "omega");					\
} while (0)


#define XRTM_SET_OMEGA_N(D, OMEGA) do {					\
     if (xrtm_set_omega_n(D, OMEGA)) {					\
          fprintf(stderr, "ERROR: xrtm_set_omega_n()\n");		\
          return -1;							\
     }									\
} while (0)


#define XRTM_SET_OMEGA_L_1N(D, I_LAYER, OMEGA_L) do {			\
     if (xrtm_set_omega_l_1n(D, I_LAYER, OMEGA_L)) {			\
          fprintf(stderr, "ERROR: xrtm_set_omega_l_1n()\n");		\
          return -1;							\
     }									\
} while (0)



/*******************************************************************************
 *
 ******************************************************************************/
#define XRTM_SET_LTAU_1(D, I_LAYER, LTAU) do {				\
     double a;								\
									\
     if (xrtm_set_ltau_1(D, I_LAYER, LTAU)) {				\
          fprintf(stderr, "ERROR: xrtm_set_ltau_1()\n");		\
          return -1;							\
     }									\
									\
     if ((a = xrtm_get_ltau(D, I_LAYER)) < -1) {			\
          fprintf(stderr, "ERROR: xrtm_get_ltau()\n");			\
          return -1;							\
     }									\
									\
     SET_SELF_CHECK(LTAU, a, "ltau");					\
} while (0)


#define XRTM_SET_LTAU_N(D, LTAU) do {					\
     if (xrtm_set_ltau_n(D, LTAU)) {					\
          fprintf(stderr, "ERROR: xrtm_set_ltau_n()\n");		\
          return -1;							\
     }									\
} while (0)


#define XRTM_SET_LTAU_L_1N(D, I_LAYER, LTAU_L) do {			\
     if (xrtm_set_ltau_l_1n(D, I_LAYER, LTAU_L)) {			\
          fprintf(stderr, "ERROR: xrtm_set_ltau_l_1n()\n");		\
          return -1;							\
     }									\
} while (0)



/*******************************************************************************
 *
 ******************************************************************************/
#define XRTM_SET_KERNEL_AMPFAC(D, I_KERNEL, AMPFAC) do {		\
     double a;								\
									\
     if (xrtm_set_kernel_ampfac(D, I_KERNEL, AMPFAC)) {			\
          fprintf(stderr, "ERROR: xrtm_set_kernel_ampfac()\n");		\
          return -1;							\
     }									\
									\
     if ((a = xrtm_get_kernel_ampfac(D, I_KERNEL)) < -1) {		\
          fprintf(stderr, "ERROR: xrtm_get_kernel_ampfac()\n");		\
          return -1;							\
     }									\
									\
     SET_SELF_CHECK(AMPFAC, a, "kernel_ampfac");			\
} while (0)


#define XRTM_SET_KERNEL_AMPFAC_L_N(D, I_KERNEL, AMPFAC) do {		\
     if (xrtm_set_kernel_ampfac_l_n(D, I_KERNEL, AMPFAC)) {		\
          fprintf(stderr, "ERROR: xrtm_set_kernel_ampfac_l_n()\n");	\
          return -1;							\
     }									\
} while (0)


#ifdef __cplusplus
}
#endif

#endif /* TEST_MACROS_H */
