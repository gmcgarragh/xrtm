#include <iostream>
using namespace std;

#include <xrtm_int_cpp.h>


int main(int argc, char *argv[]) {

     int i;
     int j;
     int k;

     double ****I_p;
     double ****I_m;
     double *****K_p;
     double *****K_m;

     /* Class that represents the XRTM instance. */
     xrtm_int_cpp *xrtm;


     /*-------------------------------------------------------------------------
      * Define inputs.
      *-----------------------------------------------------------------------*/
     /* Bit mask of XRTM options.  Use Delta-m scaling and the Nakajima and
        Tanka TMS correction.  See xrtm_model.h for other supported options. */
     int options        = xrtm_int_cpp::OPTION_CALC_DERIVS |
                          xrtm_int_cpp::OPTION_DELTA_M |
                          xrtm_int_cpp::OPTION_N_T_TMS |
                          xrtm_int_cpp::OPTION_OUTPUT_AT_LEVELS |
                          xrtm_int_cpp::OPTION_SOURCE_SOLAR;

     /* Bit mask of solvers which will be used with this XRTM instance.  Use the
        eigen/adding solver (same method as radiant).  See xrtm_model.h for
        other supported solvers. */
     int solvers         = xrtm_int_cpp::SOLVER_EIG_ADD;

     int max_coef        = 34;
     int n_quad          = 8;	/* Quadrature points per hemisphere */
     int n_stokes        = 1;
     int n_derivs        = 3;
     int n_layers        = 4;
     int n_theta_0s      = 1;
     int n_kernels       = 1;
     int n_kernel_quad   = 16;	/* Not used for the Lambertian kernel */
     xrtm_int_cpp::kernel_type
          kernels[]      = {xrtm_int_cpp::KERNEL_LAMBERTIAN};
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

     double albedo       = .2;


     /*-------------------------------------------------------------------------
      * Create an XRTM instance.
      *-----------------------------------------------------------------------*/
     xrtm = new xrtm_int_cpp(options, solvers, max_coef, n_quad, n_stokes,
                             n_derivs, n_layers, n_theta_0s, n_kernels,
                             n_kernel_quad, kernels, n_out_levels, n_out_thetas);


     /*-------------------------------------------------------------------------
      * Set inputs.
      *
      * Inputs must be set before the first model run.  For subsequent runs only
      * the inputs that change need to be set.  For example calculating the
      * radiance across the O2-A band spectrum, assuming constant scattering
      * properites, would require only updating ltau and omega for each point.
      *-----------------------------------------------------------------------*/
     try {xrtm->set_fourier_tol(.0001); }
     catch (xrtm_int_cpp::xrtm_errors &) {
          cerr << "ERROR: xrtm_set_fourier_tol()\n";
          exit(1);
     }
     try {xrtm->set_out_levels(out_levels); }
     catch (xrtm_int_cpp::xrtm_errors &) {
          cerr << "ERROR: xrtm_set_out_levels()\n";
          exit(1);
     }
     try {xrtm->set_out_thetas(out_thetas); }
     catch (xrtm_int_cpp::xrtm_errors &) {
          cerr << "ERROR: xrtm_set_out_thetas()\n";
          exit(1);
     }
     try {xrtm->set_F_iso_top(0.); }
     catch (xrtm_int_cpp::xrtm_errors &) {
          cerr << "ERROR: xrtm_set_F_iso_top()\n";
          exit(1);
     }
     try {xrtm->set_F_iso_bot(0.); }
     catch (xrtm_int_cpp::xrtm_errors &) {
          cerr << "ERROR: xrtm_set_F_iso_bot()\n";
          exit(1);
     }
     try {xrtm->set_F_0(F_0); }
     catch (xrtm_int_cpp::xrtm_errors &) {
          cerr << "ERROR: xrtm_set_F_0()\n";
          exit(1);
     }
     try {xrtm->set_theta_0(theta_0); }
     catch (xrtm_int_cpp::xrtm_errors &) {
          cerr << "ERROR: xrtm_set_theta_0()\n";
          exit(1);
     }
     try {xrtm->set_phi_0(0.); }
     catch (xrtm_int_cpp::xrtm_errors &) {
          cerr << "ERROR: xrtm_set_phi_0()\n";
          exit(1);
     }


     /* Set optical property inputs */
     try {xrtm->set_ltau_n(ltau); }
     catch (xrtm_int_cpp::xrtm_errors &) {
          cerr << "ERROR: xrtm_set_ltau_n()\n";
          exit(1);
     }
     try {xrtm->set_omega_n(omega); }
     catch (xrtm_int_cpp::xrtm_errors &) {
          cerr << "ERROR: xrtm_set_omega_n()\n";
          exit(1);
     }
     try {xrtm->set_coef_n(n_coef, coef); }
     catch (xrtm_int_cpp::xrtm_errors &) {
          cerr << "ERROR: xrtm_set_coef_n()\n";
          exit(1);
     }

     /* Alternatively the inputs can be set one layer at a time. */
     for (i = 0; i < n_layers; ++i) {
          try {xrtm->set_ltau_1(i, ltau[i]); }
     catch (xrtm_int_cpp::xrtm_errors &) {
               cerr << "ERROR: xrtm_set_ltau_1()\n";
               exit(1);
          }
          try {xrtm->set_omega_1(i, omega[i]); }
     catch (xrtm_int_cpp::xrtm_errors &) {
               cerr << "ERROR: xrtm_set_omega_1()\n";
               exit(1);
          }
          try {xrtm->set_coef_1(i, n_coef[i], coef[i]); }
     catch (xrtm_int_cpp::xrtm_errors &) {
               cerr << "ERROR: xrtm_set_coef_1()\n";
               exit(1);
          }
     }

     /* Set surface albedo */
     try {xrtm->set_kernel_ampfac(0, albedo); }
     catch (xrtm_int_cpp::xrtm_errors &) {
          cerr << "ERROR: xrtm_set_kernel_ampfac()\n";
          exit(1);
     }


     /*-------------------------------------------------------------------------
      * Set linearized inputs.
      *-----------------------------------------------------------------------*/
     try {xrtm->set_ltau_l_11(2, 0, 1.); }
     catch (xrtm_int_cpp::xrtm_errors &) {
               cerr << "ERROR: xrtm_set_ltau_l_11()\n";
          exit(1);
     }
     try {xrtm->set_omega_l_11(2, 1, 1.); }
     catch (xrtm_int_cpp::xrtm_errors &) {
               cerr << "ERROR: xrtm_set_omega_l_11()\n";
          exit(1);
     }

     try {xrtm->set_kernel_ampfac_l_1(0, 2, 1.); }
     catch (xrtm_int_cpp::xrtm_errors &) {
               cerr << "ERROR: xrtm_set_kernel_ampfac_l_1()\n";
          exit(1);
     }

     /* Must be called everytime linearized values are changed */
     try {xrtm->update_varied_layers(); }
     catch (xrtm_int_cpp::xrtm_errors &) {
               cerr << "ERROR: xrtm_update_varied_layers()\n";
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
     try {xrtm->radiance(xrtm_int_cpp::SOLVER_EIG_ADD, n_out_phis, out_phis,
                         I_p, I_m, K_p, K_m); }
     catch (xrtm_int_cpp::xrtm_errors &) {
          cerr << "ERROR: xrtm_radiance()\n";
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
      * Delete the model.
      *-----------------------------------------------------------------------*/
     delete xrtm;


     exit(0);
}
