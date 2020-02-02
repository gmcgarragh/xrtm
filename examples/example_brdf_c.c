/*******************************************************************************
 * This is an example of creating an XRTM instance, setting inputs, running the
 * actual model and destroying the instance.  The example includes a 4 layer
 * atmosphere with Rayleigh scattering in layers 1, 2, and 4, aerosols in layer
 * 3, and a lambertian surface.  The aerosol is the accumulation mode mineral
 * dust described by d'Almeida, 1991 except with an effective variance of 0.2.
 * Both the Rayleigh and aerosol scattering are for 2.25 um.  Along with
 * intensities, derivatives for intensity with respect to the optical depth and
 * single scattering albedo in the aerosol layer and with intensity with respect
 * to the surface albedo are also determined.
 ******************************************************************************/
#include <stdio.h>
#include <stdlib.h>

#include <gutil.h>

#include <xrtm_interface.h>

/*
#define XRTM_BRDF
*/

typedef struct {
     double a;
     double b;
     double c;
     double d;
     double e;
     double g;
     double h;
     double R_ir;
     double T_i;
     double T_r;
} hapke_aux_data;


void hapke_aux(brdf_aux_data *aux, hapke_aux_data *aux2, int i, int j, int n_derivs, uchar *derivs, double *p, double **p_l);
double hapke_kernel(brdf_aux_data *aux, hapke_aux_data *aux2, enum xrtm_brdf_geometry brdf_geom, int i, int j, int k, double phi, int n_derivs, uchar *derivs, double *p, double **p_l, double *f_l);


int main(int argc, char *argv[]) {

     int i;
     int j;
     int k;

     double ****I_p;
     double ****I_m;
     double *****K_p;
     double *****K_m;

     /* Structure the represents the XRTM instance the members of which should
        never be modfied by the user. */
     xrtm_data xrtm;


     /*-------------------------------------------------------------------------
      * Define inputs.
      *-----------------------------------------------------------------------*/
     /* Bit mask of XRTM options.  Use Delta-m scaling and the Nakajima and
        Tanka TMS correction.  See xrtm_model.h for other supported options. */
     int options        = XRTM_OPTION_CALC_DERIVS | XRTM_OPTION_DELTA_M |
                          XRTM_OPTION_N_T_TMS | XRTM_OPTION_OUTPUT_AT_LEVELS |
                          XRTM_OPTION_SOURCE_SOLAR;

     /* Bit mask of solvers which will be used with this XRTM instance.  Use the
        eigen/adding solver (same method as radiant).  See xrtm_model.h for
        other supported solvers. */
     int solvers         = XRTM_SOLVER_EIG_ADD;

     int max_coef        = 34;
     int n_quad          = 8;	/* Quadrature points per hemisphere */
     int n_stokes        = 1;
     int n_derivs        = 3;
     int n_layers        = 4;
     int n_theta_0s      = 1;
     int n_kernels       = 1;
     int n_kernel_quad   = 16;	/* Not used for the Lambertian kernel */
#ifdef XRTM_BRDF
     int kernels[]       = {XRTM_KERNEL_HAPKE};
#else
     int kernels[]       = {XRTM_KERNEL_USER_DEFINED};
#endif
     int n_out_levels    = 2;
     int n_out_thetas    = 3;
     int n_out_phis      = 1;

     double F_0          = 1.;
     double theta_0      = 35.;

     int out_levels[]    = {0, 4};

     double out_thetas[] = {0., 30., 60.};
     double phi[]        = {45.};
     double *out_phis[]  = {phi, phi, phi};

     double ltau []      = {.02, .05, 1., .1};

     double omega[]      = {1., 1., .9, 1.};

     int n_ray_coef      = 3;

     double ray_coef1[]  = {1.000000e+00,
                            0.000000e+00,
                            4.798741e-01};
     double *ray_coef[]  = {ray_coef1};

     int n_aer_coef      = 34;

     double aer_coef1[]  = {1.000000e+00,
                            1.865569e+00,
                            1.789985e+00,
                            1.220838e+00,
                            7.472409e-01,
                            4.017337e-01,
                            2.173326e-01,
                            1.054020e-01,
                            5.737447e-02,
                            2.570752e-02,
                            1.527185e-02,
                            6.202491e-03,
                            4.278587e-03,
                            1.529611e-03,
                            1.276447e-03,
                            3.964385e-04,
                            4.036524e-04,
                            1.112891e-04,
                            1.338887e-04,
                            3.468251e-05,
                            4.611093e-05,
                            1.204792e-05,
                            1.637357e-05,
                            4.577401e-06,
                            5.975423e-06,
                            1.849954e-06,
                            2.241820e-06,
                            7.774087e-07,
                            8.673507e-07,
                            3.351400e-07,
                            3.476180e-07,
                            1.472730e-07,
                            1.448956e-07,
                            6.591328e-08};
     double *aer_coef[]  = {aer_coef1};

     int n_coef[]        = {n_ray_coef, n_ray_coef, n_aer_coef, n_ray_coef};

     double **coef[4]    = {ray_coef, ray_coef, aer_coef, ray_coef};

     double brdf_ampfac  = {.5};

     double brdf_params[]= {.25, .5, .75};


     /*-------------------------------------------------------------------------
      * Create an XRTM instance.
      *-----------------------------------------------------------------------*/
     if (xrtm_create(&xrtm, options, solvers, max_coef, n_quad, n_stokes, n_derivs,
                     n_layers, n_theta_0s, n_kernels, n_kernel_quad,
                     (enum xrtm_kernel_type *) kernels, n_out_levels, n_out_thetas)) {
          fprintf(stderr, "error calling xrtm_create()\n");
          exit(1);
     }


     /*-------------------------------------------------------------------------
      * Set inputs.
      *
      * Inputs must be set before the first model run.  For subsequent runs only
      * the inputs that change need to be set.  For example calculating the
      * radiance across the O2-A band spectrum, assuming constant scattering
      * properites, would require only updating ltau and omega for each point.
      *-----------------------------------------------------------------------*/
     if (xrtm_set_fourier_tol(&xrtm, .0001)) {
          fprintf(stderr, "ERROR: xrtm_set_fourier_tol()\n");
          exit(1);
     }
     if (xrtm_set_out_levels(&xrtm, out_levels)) {
          fprintf(stderr, "ERROR: xrtm_set_out_levels()\n");
          exit(1);
     }
     if (xrtm_set_out_thetas(&xrtm, out_thetas)) {
          fprintf(stderr, "ERROR: xrtm_set_out_thetas()\n");
          exit(1);
     }
     if (xrtm_set_F_iso_top(&xrtm, 0.)) {
          fprintf(stderr, "ERROR: xrtm_set_F_iso_top()\n");
          exit(1);
     }
     if (xrtm_set_F_iso_bot(&xrtm, 0.)) {
          fprintf(stderr, "ERROR: xrtm_set_F_iso_bot()\n");
          exit(1);
     }
     if (xrtm_set_F_0(&xrtm, F_0)) {
          fprintf(stderr, "ERROR: xrtm_set_F_0()\n");
          exit(1);
     }
     if (xrtm_set_theta_0(&xrtm, theta_0)) {
          fprintf(stderr, "ERROR: xrtm_set_theta_0()\n");
          exit(1);
     }
     if (xrtm_set_phi_0(&xrtm, 0.)) {
          fprintf(stderr, "ERROR: xrtm_set_phi_0()\n");
          exit(1);
     }

     /* Set optical property inputs */
     if (xrtm_set_ltau_n(&xrtm, ltau)) {
          fprintf(stderr, "ERROR: xrtm_set_ltau_n()\n");
          exit(1);
     }
     if (xrtm_set_omega_n(&xrtm, omega)) {
          fprintf(stderr, "ERROR: xrtm_set_omega_n()\n");
          exit(1);
     }
     if (xrtm_set_coef_n(&xrtm, n_coef, coef)) {
          fprintf(stderr, "ERROR: xrtm_set_coef_n()\n");
          exit(1);
     }

     /* Alternatively the inputs can be set one layer at a time. */
     for (i = 0; i < n_layers; ++i) {
          if (xrtm_set_ltau_1(&xrtm, i, ltau[i])) {
               fprintf(stderr, "ERROR: xrtm_set_ltau_1()\n");
               exit(1);
          }
          if (xrtm_set_omega_1(&xrtm, i, omega[i])) {
               fprintf(stderr, "ERROR: xrtm_set_omega_1()\n");
               exit(1);
          }
          if (xrtm_set_coef_1(&xrtm, i, n_coef[i], coef[i])) {
               fprintf(stderr, "ERROR: xrtm_set_coef_1()\n");
               exit(1);
          }
     }

     /* Set surface albedo */
     if (xrtm_set_kernel_ampfac(&xrtm, 0, brdf_ampfac)) {
          fprintf(stderr, "ERROR: xrtm_set_kernel_ampfac()\n");
          exit(1);
     }

     if (xrtm_set_kernel_params_n(&xrtm, 0, brdf_params)) {
          fprintf(stderr, "ERROR: xrtm_set_kernel_params()\n");
          exit(1);
     }
#ifndef XRTM_BRDF
     if (xrtm_set_kernel_func(&xrtm, 0, (void (*)(brdf_aux_data *aux, void *aux2, int i, int j, int n_derivs, uchar *derivs, double *p, double **p_l)) hapke_aux,
                              (double (*)(brdf_aux_data *aux, void *aux2, enum xrtm_brdf_geometry brdf_geom, int i, int j, int k, double phi, int n_derivs, uchar *derivs, double *p, double **p_l, double *f_l)) hapke_kernel)) {
          fprintf(stderr, "ERROR: xrtm_set_kernel_func()\n");
          exit(1);
     }
#endif

     /*-------------------------------------------------------------------------
      * Set linearized inputs.
      *-----------------------------------------------------------------------*/
     if (xrtm_set_ltau_l_11(&xrtm, 2, 0, 1.)) {
          fprintf(stderr, "ERROR: xrtm_set_ltau_l_11()\n");
          exit(1);
     }
     if (xrtm_set_omega_l_11(&xrtm, 2, 1, 1.)) {
          fprintf(stderr, "ERROR: xrtm_set_omega_l_11()\n");
          exit(1);
     }

     if (xrtm_set_kernel_ampfac_l_1(&xrtm, 0, 2, 1.)) {
          fprintf(stderr, "ERROR: xrtm_set_kernel_ampfac_l_1()\n");
          exit(1);
     }

     /* Must be called everytime linearized values are changed */
     if (xrtm_update_varied_layers(&xrtm)) {
          fprintf(stderr, "ERROR: xrtm_update_varied_layers()\n");
          exit(1);
     }


     /*-------------------------------------------------------------------------
      * Currently the XRTM C interface uses arrays of pointers to arrays to
      * pointers and so on ... for multidimensioanl array io.  This may change
      * if true multidimensioanl arrays (row major) prove to be more convient.
      * alloc_array*_d() allocates pointer based multidimensional arrays.
      *-----------------------------------------------------------------------*/
     I_p = alloc_array4_d(n_out_levels,           n_out_thetas, n_out_phis, n_stokes);
     I_m = alloc_array4_d(n_out_levels,           n_out_thetas, n_out_phis, n_stokes);
     K_p = alloc_array5_d(n_out_levels, n_derivs, n_out_thetas, n_out_phis, n_stokes);
     K_m = alloc_array5_d(n_out_levels, n_derivs, n_out_thetas, n_out_phis, n_stokes);


     /*-------------------------------------------------------------------------
      * Run the model for radiances and associated derivatives.  If this is the
      * initial run and all the required inputs have not been initialized then
      * XRTM will print a appropriate message and return < 0.
      *-----------------------------------------------------------------------*/
     if (xrtm_radiance(&xrtm, XRTM_SOLVER_EIG_ADD, n_out_phis, out_phis, I_p, I_m, K_p, K_m)) {
          fprintf(stderr, "ERROR: xrtm_radiance()\n");
          exit(1);
     }


     /*-------------------------------------------------------------------------
      * Output results.
      *-----------------------------------------------------------------------*/
#if PLATFORM == WIN32_MSVC
     _set_output_format(_TWO_DIGIT_EXPONENT);
#endif
     for (i = 0; i < n_out_levels; ++i) {
          printf("level: %d\n", i);
          printf("     intensity:\n");
          for (j = 0; j < n_out_thetas; ++j) {
               printf("          theta = %9.2E, I_p = %13.6E, I_m = %13.6E\n",
                      out_thetas[j], I_p[i][j][0][0], I_m[i][j][0][0]);
          }
          for (j = 0; j < n_derivs; ++j) {
               printf("     derivative: %d\n", j);
               for (k = 0; k < n_out_thetas; ++k) {
                    printf("          theta = %9.2E, K_p = %13.6E, K_m = %13.6E\n",
                           out_thetas[k], K_p[i][j][k][0][0], K_m[i][j][k][0][0]);
               }
          }
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     free_array4_d(I_p);
     free_array4_d(I_m);
     free_array5_d(K_p);
     free_array5_d(K_m);


     /*-------------------------------------------------------------------------
      * Destroy the model.
      *-----------------------------------------------------------------------*/
     xrtm_destroy(&xrtm);


     exit(0);
}



/*******************************************************************************
 * p[0] = Sigma (single scattering albedo)
 * p[1] = B_0   (amplitude of hot-spot)
 * p[2] = Delta (angular width of hot-spot)
 ******************************************************************************/
void hapke_aux(brdf_aux_data *aux, hapke_aux_data *aux2, int i, int j, int n_derivs, uchar *derivs, double *p, double **p_l) {

     double a;
     double b;

     double gamma;

     a = (aux->cos_theta[i] + aux->cos_theta[j]) * 4.;
     aux2->R_ir = p[0] / a;

     aux2->a    =  aux->cos_theta[i] * aux->cos_theta[j];
     aux2->b    =  aux->sin_theta[i] * aux->sin_theta[j];

     aux2->c    = p[1] * p[2];

     gamma      = sqrt(1. - p[0]);

     a          = 2. * aux->cos_theta[i];
     b          = 2. * aux->cos_theta[j];

     aux2->d    = ((1. + a) / (1. + a * gamma)) * ((1. + b) / (1. + b * gamma));
}



double hapke_kernel(brdf_aux_data *aux, hapke_aux_data *aux2, enum xrtm_brdf_geometry brdf_geom, int i, int j, int k, double phi, int n_derivs, uchar *derivs, double *p, double **p_l, double *f_l) {

     /* Spurr 2004 (contains a typo) */

     double cos_ksi;
     double ksi;
     double P;
     double B;

     cos_ksi = aux2->a + aux2->b * aux->cos_phi[k];

     P       = 1. + .5 * cos_ksi;

     ksi     = acos(cos_ksi);

     B       = aux2->c / (p[2] + tan(ksi / 2.));

     return aux2->R_ir * ((1. + B) * P + aux2->d - 1.);
}
