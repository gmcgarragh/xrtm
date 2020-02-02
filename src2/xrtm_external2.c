/*******************************************************************************
**
**    Copyright (C) 2007-2020 Greg McGarragh <greg.mcgarragh@colostate.edu>
**
**    This source code is licensed under the GNU General Public License (GPL),
**    Version 3.  See the file COPYING for more details.
**
*******************************************************************************/

#define RADIANCE_SCALE 1.e28


double plkavg_(double *wvnmlo, double *wvnmhi, double *t);

void radtran3_planck_function_(double *temp, char *units, double *wavelength, double *planck);


/*******************************************************************************
 *
 ******************************************************************************/
#ifdef __cplusplus
extern "C" {
#endif
void disort_(
       int *nlyr, double *dtauc, double *ssalb,
       int *nmom, double *pmom, double *temper,
       double *wvnmlo,  double *wvnmhi,
       int *usrtau, int *ntau, double *utau, int *nstr,
       int *usrang, int *numu, double *umu, int *nphi, double *phi,
       int *ibcnd, double *fbeam, double *umu0, double *phi0,
       double *fisot, int *lamber, double *albedo,
       double *btemp, double *ttemp, double *temis,
       int *plank, int *onlyfl, double *accur, int *prnt, char *header,
       int *maxcly, int *maxulv, int *maxumu, int *maxphi, int *maxmom,
       double *rfldir, double *rfldn, double *flup, double *dfdt,
       double *uavg, double *uu, double *albmed, double *trnmed,
       int *deltam, int *corint);
#ifdef __cplusplus
}
#endif

int call_disort(int n_coef, int n_quad, int n_layers, double *qx, double lambda, double F_0, double mu_0, double phi_0, int *ulevels, double *utaus, int n_ulevels, double *umus, int n_umus, double *phis, int n_phis, double F_iso_top, double *levels_b, double surface_b, double albedo_, int *n_coef_layer, double ***coefs, double *omega, double *ltau, double ****I_p, double ****I_m, double *mean_p, double *mean_m, double *flux_p, double *flux_m, double fourier_tol, int delta_m, int n_t_tms, int radiance, int mean, int flux, int solar, int thermal, int utau_output) {

     int i;
     int ii;
     int j;
     int k;

     char header[127];

     int nlyr   = n_layers;
     int nmom   = n_coef - 1;
     int usrtau = 1;
     int ntau;
     int nstr   = 2 * n_quad;
     int usrang;
     int numu;
     int nphi   = n_phis;
     int ibcnd  = 0;
     int lamber = 1;
     int plank  = thermal;
     int onlyfl = ! radiance;
     int maxcly = n_layers;
     int maxulv = n_ulevels;
     int maxumu = MAX(nstr, 2 * n_umus);
     int maxphi = n_phis;
     int maxmom = n_coef - 1;
     int deltam = delta_m;
     int corint = n_t_tms;

     int prnt[5] = {0, 0, 0, 0, 0};

     int *index;

     double ttau;

     double wvnmhi = 1. / ((lambda - .5) * 1.e-6 * 1.e2);
     double wvnmlo = 1. / ((lambda + .5) * 1.e-6 * 1.e2);
     double fbeam  = 0.;
     double umu0   = mu_0;
     double phi0   = phi_0;
     double fisot  = F_iso_top;
     double albedo = albedo_;
     double btemp;
     double ttemp;
     double temis  = 1.;
     double accur  = fourier_tol;
/*
     double *dtauc;
     double *ssalb;
*/
     double **pmom;
     double *temper;
     double *utau;
     double *umu;
/*
     double *phi;
*/
     double *rfldir;
     double *rfldn;
     double *flup;
     double *dfdt;
     double *uavg;
     double *albmed;
     double *trnmed;

     double ***uu;
/*
     dtauc  = alloc_array1_d(maxcly);
     ssalb  = alloc_array1_d(maxcly);
*/
     pmom   = alloc_array2_d(maxcly, maxmom + 1);
     temper = alloc_array1_d(maxcly + 1);
     utau   = alloc_array1_d(maxulv);
     umu    = alloc_array1_d(maxumu);
/*
     phi    = alloc_array1_d(maxphi);
*/
     rfldir = alloc_array1_d(maxulv);
     rfldn  = alloc_array1_d(maxulv);
     flup   = alloc_array1_d(maxulv);
     dfdt   = alloc_array1_d(maxulv);
     uavg   = alloc_array1_d(maxulv);
     albmed = alloc_array1_d(maxumu);
     trnmed = alloc_array1_d(maxumu);

     uu = alloc_array3_d(maxphi, maxulv, maxumu);

     memset(header, ' ', 127);

/*
     for (i = 0; i < n_layers; ++i) {
          dtauc[i] = ltau [i];
          ssalb[i] = omega[i];
     }
*/
     if (solar)
          fbeam  = F_0;
     if (thermal) {
          ttemp = 0.;

          for (i = 0; i < n_layers + 1; ++i)
               temper[i] = levels_b[i];

          btemp = surface_b;
     }

     for (i = 0; i < n_layers; ++i) {
          for (j = 0; j < n_coef_layer[i]; ++j) {
               pmom[i][j] = coefs[i][0][j] / (2. * j + 1.);
          }

          memset(&pmom[i][j], 0, (n_coef - j) * sizeof(double));
     }

     if (n_umus == 0)
          usrang = 0;
     else {
          usrang = 1;

          numu = n_umus * 2;

          index = alloc_array1_i(n_umus);

          for (i = 0; i < n_umus; ++i)
               index[i] = i;

          quick_sort_index(umus, index, n_umus, sizeof(double), (int (*)(const void *, const void *)) cmp_p_d);

          ii = 0;
          for (i = n_umus-1; i >= 0; --i)
               umu[ii++] = -umus[index[i]];
          for (i = 0; i < n_umus; ++i)
               umu[ii++] =  umus[index[i]];
     }


     nphi = n_phis;
/*
     for (i = 0; i < n_phis; ++i)
          phi[i] = phis[i];
*/
     ntau = n_ulevels;
     if (! utau_output) {
          ii = 0;
          ttau = 0.;
          for (i = 0; i < n_layers + 1; ++i) {
               if (i == ulevels[ii])
                    utau[ii++] = ttau;
               if (ii >= n_ulevels)
                    break;
               ttau += ltau[i];
          }
     }
     else {
          for (i = 0; i < n_ulevels; ++i)
               utau[i] = utaus[i];
     }

     if (radiance) {
          if (! uu) {
               fprintf(stderr, "ERROR: memory allocation failed\n");
               return -1;
          }
     }

     disort_(&nlyr, ltau, omega,
             &nmom, *pmom, temper,
             &wvnmlo,  &wvnmhi,
             &usrtau, &ntau, utau, &nstr,
             &usrang, &numu, umu, &nphi, phis,
             &ibcnd, &fbeam, &umu0, &phi0,
             &fisot, &lamber, &albedo,
             &btemp, &ttemp, &temis,
             &plank, &onlyfl, &accur, prnt, header,
             &maxcly, &maxulv, &maxumu, &maxphi, &maxmom,
             rfldir, rfldn, flup, dfdt,
             uavg, **uu, albmed, trnmed, &deltam, &corint);

     if (radiance) {
          for (i = 0; i < n_ulevels; ++i) {
               for (j = 0; j < n_phis; ++j) {
                    if (n_umus == 0) {
                         for (k = 0; k < n_quad; ++k) {
                              I_p[i][k][j][0] = uu[j][i][n_quad + k    ];
                              I_m[i][k][j][0] = uu[j][i][n_quad - k - 1];
                         }
                    }
                    else {
                         for (k = 0; k < n_umus ; ++k) {
                              I_p[i][k][j][0] = uu[j][i][n_umus  + index[k]    ];
                              I_m[i][k][j][0] = uu[j][i][n_umus  - index[k] - 1];
                         }
                    }
               }
          }
     }

     if (mean) {
          for (i = 0; i < n_ulevels; ++i) {
               mean_p[i] = uavg[i];
               mean_m[i] = uavg[i];
          }
     }

     if (flux) {
          for (i = 0; i < n_ulevels; ++i) {
               flux_p[i] = flup[i];
               flux_m[i] = rfldir[i] + rfldn[i];
          }
     }
/*
     free_array1_d(dtauc);
     free_array1_d(ssalb);
*/
     free_array2_d(pmom);
     free_array1_d(temper);
     free_array1_d(utau);
     free_array1_d(umu);
/*
     free_array1_d(phi);
*/
     free_array1_d(rfldir);
     free_array1_d(rfldn);
     free_array1_d(flup);
     free_array1_d(dfdt);
     free_array1_d(uavg);
     free_array1_d(albmed);
     free_array1_d(trnmed);

     free_array3_d(uu);

     if (n_umus != 0)
          free_array1_i(index);

     return 0;
}



/*******************************************************************************
 *
 ******************************************************************************/
#ifdef __cplusplus
extern "C" {
#endif
void call_lidort_f_(int *n_four,
                    int *n_elem,
                    int *n_coef,
                    int *n_quad,
                    int *n_derivs,
                    int *n_layers,
                    double *qx,
                    double *F_0,
                    double *sza,
                    double *phi_0,
                    int *ulevels,
                    double *utaus,
                    int *n_ulevels,
                    double *umus,
                    int *n_umus,
                    double *phis,
                    int *n_phis,
                    double *planet_r,
                    double *levels_z,
                    int *n_kernels,
                    int *n_kernel_quad,
                    int *kernels,
                    double *ampfac,
                    double *ampfac_l,
                    double *params,
                    double *params_l,
                    int *n_coef_layer,
                    double *coefs,
                    double *coefs_l,
                    double *omega,
                    double *omega_l,
                    double *ltau,
                    double *ltau_l,
                    double *I_p,
                    double *I_m,
                    double *I_p_l,
                    double *I_m_l,
                    double *mean_p,
                    double *mean_m,
                    double *mean_p_l,
                    double *mean_m_l,
                    double *flux_p,
                    double *flux_m,
                    double *flux_p_l,
                    double *flux_m_l,
                    int *delta_m,
                    int *n_t_tms,
                    int *psa,
                    int *quad_output,
                    int *radiance,
                    int *mean,
                    int *flux,
                    int *utau_output,
                    uchar *derivs,
                    double *epsilon,
                    int *info,
                    int *n_mus2);
#ifdef __cplusplus
}
#endif

int call_lidort(int n_four, int n_elem, int n_coef, int n_quad, int n_derivs, int n_layers, double *qx, double F_0, double mu_0, double phi_0, int *ulevels, double *utaus, int n_ulevels, double *umus, int n_umus, double *phis, int n_phis, double planet_r, double *levels_z, int n_kernels, int n_kernel_quad, enum xrtm_kernel_type *kernels, double *ampfac, double **ampfac_l, double **params, double ***params_l, int *n_coef_layer, double ***coefs, double ****coefs_l, double *omega, double **omega_l, double *ltau, double **ltau_l, double ****I_p, double ****I_m, double *****I_p_l, double *****I_m_l, double *mean_p, double *mean_m, double **mean_p_l, double **mean_m_l, double *flux_p, double *flux_m, double **flux_p_l, double **flux_m_l, double fourier_tol, int delta_m, int n_t_tms, int psa, int quad_output, int radiance, int mean, int flux, int utau_output, uchar **derivs) {

     uchar *derivs2;

     int info;

     int n_mus2;

     double sza;

     double *coefs_l2;
     double *omega_l2;
     double *ltau_l2;

     double *params2;

     double *ampfac_l2;
     double *params_l2;

     double *I_p2;
     double *I_m2;

     double *I_p_l2;
     double *I_m_l2;

     double *mean_p_l2;
     double *mean_m_l2;
     double *flux_p_l2;
     double *flux_m_l2;

     sza = acos(mu_0)*R2D;

     if (mu_0 == 1.)
          psa = 0;

     if (n_umus == 0)
          n_mus2 = n_quad;
     else
          n_mus2 = n_umus;

     if (radiance) {
          I_p2 = ***I_p;
          I_m2 = ***I_m;
     }

     if (n_kernels > 0)
          params2 = *params;

     if (n_derivs > 0) {
          coefs_l2 = ***coefs_l;
          omega_l2 = *omega_l;
          ltau_l2  = *ltau_l;
          derivs2  = *derivs;

          if (radiance) {
               I_p_l2    = ****I_p_l;
               I_m_l2    = ****I_m_l;
          }
          if (mean) {
               mean_p_l2 = *mean_p_l;
               mean_m_l2 = *mean_m_l;
          }
          if (flux) {
               flux_p_l2 = *flux_p_l;
               flux_m_l2 = *flux_m_l;
          }

          if (n_kernels > 0) {
               ampfac_l2 = *ampfac_l;
               params_l2 = **params_l;
          }
     }

     call_lidort_f_(&n_four, &n_elem, &n_coef, &n_quad, &n_derivs, &n_layers, qx, &F_0, &sza, &phi_0, ulevels, utaus, &n_ulevels, umus, &n_umus, phis, &n_phis, &planet_r, levels_z, &n_kernels, &n_kernel_quad, (int *) kernels, ampfac, ampfac_l2, params2, params_l2, n_coef_layer, **coefs, coefs_l2, omega, omega_l2, ltau, ltau_l2, I_p2, I_m2, I_p_l2, I_m_l2, mean_p, mean_m, mean_p_l2, mean_m_l2, flux_p, flux_m, flux_p_l2, flux_m_l2, &delta_m, &n_t_tms, &psa, &quad_output, &radiance, &mean, &flux, &utau_output, derivs2, &fourier_tol, &info, &n_mus2);

     if (info) {
          fprintf(stderr, "ERROR: call_lidort_f()\n");
          return -1;
     }

     return 0;
}



/*******************************************************************************
 *
 ******************************************************************************/
#ifdef __cplusplus
extern "C" {
#endif
void call_lrad_f_(int *n_four,
                  int *n_elem,
                  int *n_coef,
                  int *n_quad,
                  int *n_stokes,
                  int *n_derivs,
                  int *n_layers,
                  double *qx,
                  double *F_0,
                  double *mu_0,
                  double *phi_0,
                  int *n_ulevels,
                  double *umus,
                  int *n_umus,
                  double *phis,
                  int *n_phis,
                  double *levels_z,
                  int *n_kernels,
                  int *n_kernel_quad,
                  int *kernels,
                  double *ampfac,
                  double *ampfac_l,
                  double *params,
                  double *params_l,
                  int *n_coef_layer,
                  double *coefs0,
                  double *coefs0_l,
                  double *coefs,
                  double *coefs_l,
                  double *omega0,
                  double *omega0_l,
                  double *omega,
                  double *omega_l,
                  double *ltau,
                  double *ltau_l,
                  double *I_p,
                  double *I_m,
                  double *I_p_l,
                  double *I_m_l,
                  int *n_t_tms,
                  int *psa,
                  uchar *derivs,
                  double *epsilon,
                  int *info,
                  int *n_mus2);
#ifdef __cplusplus
}
#endif

int call_lrad(int n_four, int n_elem, int n_coef, int n_quad, int n_stokes, int n_derivs, int n_layers, double *qx, double lambda, double F_0, double mu_0, double phi_0, int n_ulevels, double *umus, int n_umus, double *phis, int n_phis, double planet_r, double *levels_z, double *levels_t, double surface_t, int n_kernels, int n_kernel_quad, enum xrtm_kernel_type *kernels, double *ampfac, double **ampfac_l, double **params, double ***params_l, int *n_coef_layer, double ***coefs0, double ****coefs0_l, double ***coefs, double ****coefs_l, double *omega0, double **omega0_l, double *omega, double **omega_l, double *ltau, double **ltau_l, double ****I_p, double ****I_m, double *****I_p_l, double *****I_m_l, double fourier_tol, int delta_m, int n_t_tms, int psa, uchar **derivs) {

     uchar *derivs2;

     int info;

     int n_mus2;

     double *coefs0_l2;
     double *coefs_l2;
     double *omega0_l2;
     double *omega_l2;
     double *ltau_l2;

     double *params2;

     double *ampfac_l2;
     double *params_l2;

     double *I_p_l2;
     double *I_m_l2;

     if (mu_0 == 1.)
          psa = 0;

     if (n_umus == 0)
          n_mus2 = n_quad;
     else
          n_mus2 = n_umus;

     if (n_kernels > 0)
          params2 = *params;

     if (n_derivs > 0) {
          coefs0_l2 = ***coefs0_l;
          coefs_l2  = ***coefs_l;
          omega0_l2 = *omega0_l;
          omega_l2  = *omega_l;
          ltau_l2   = *ltau_l;
          derivs2   = *derivs;

          I_p_l2    = ****I_p_l;
          I_m_l2    = ****I_m_l;

          if (n_kernels > 0) {
               ampfac_l2 = *ampfac_l;
               params_l2 = **params_l;
          }
     }

     call_lrad_f_(&n_four, &n_elem, &n_coef, &n_quad, &n_stokes, &n_derivs, &n_layers, qx, &F_0, &mu_0, &phi_0, &n_ulevels, umus, &n_umus, phis, &n_phis, levels_z, &n_kernels, &n_kernel_quad, (int *) kernels, ampfac, ampfac_l2, params2, params_l2, n_coef_layer, **coefs0, coefs0_l2, **coefs, coefs_l2, omega0, omega0_l2, omega, omega_l2, ltau, ltau_l2, ***I_p, ***I_m, I_p_l2, I_m_l2, &n_t_tms, &psa, derivs2, &fourier_tol, &info, &n_mus2);

     if (info) {
          fprintf(stderr, "ERROR: call_lrad_f()\n");
          return -1;
     }

     return 0;
}



/*******************************************************************************
 *
 ******************************************************************************/
static void xrtm_to_polrad_coefs(int n_coef, int n_stokes, double **coefs, double ***coefs2) {

     int j;
     int k;
     int l;

     for (j = 0; j < n_coef; ++j) {
          for (k = 0; k < 4; ++k) {
               for (l = 0; l < 4; ++l) {
                    coefs2[j][k][l] = 0.;
               }
          }

          coefs2[j][0][0] =  coefs[0][j];        /* 0,0  */
if (n_stokes > 1) {
          coefs2[j][1][1] =  coefs[1][j];        /* 1,1  */
          coefs2[j][1][0] = -coefs[4][j];        /* 0,1  */
          coefs2[j][0][1] = -coefs[4][j];        /* 1,0  */
          coefs2[j][2][2] =  coefs[2][j];        /* 2,2  */
          coefs2[j][3][3] =  coefs[3][j];        /* 3,3  */
          coefs2[j][3][2] = -coefs[5][j];        /* 2,3  */
          coefs2[j][2][3] = -coefs[5][j];        /* 3,2  */
}
     }
}


#ifdef __cplusplus
extern "C" {
#endif
void call_polrad_f_(int *n_four,
                    int *n_quad,
                    int *n_stokes,
                    int *n_layers,
                    double *qx,
                    double *F_0,
                    double *theta_0,
                    double *phi_0,
                    int *n_ulevels,
                    double *umus,
                    int *n_umus,
                    double *phis,
                    int *n_phis,
                    double *albedo,
                    int *n_coef_layer,
                    double *coefs,
                    double *omega,
                    double *ltau,
                    double *I_p,
                    double *I_m,
                    int *delta_m,
                    double *epsilon,
                    int *info,
                    int *n_mus2);
#ifdef __cplusplus
}
#endif

int call_polrad(int n_four, int n_coef, int n_quad, int n_stokes, int n_layers, double *qx, double F_0, double mu_0, double phi_0, int n_ulevels, double *umus, int n_umus, double *phis, int n_phis, double albedo, int *n_coef_layer, double ***coefs, double *omega, double *ltau, double ****I_p, double ****I_m, double fourier_tol, int delta_m) {

     int i;

     int info;

     int n_mus2;

     int *n_coef_layer2;

     double theta_0;

     double *omega2;
     double *ltau2;

     double ****coefs2;
fourier_tol = DBL_EPSILON;
     if (n_umus == 0)
          n_mus2 = n_quad;
     else
          n_mus2 = n_umus;

     theta_0 = acos(mu_0)*R2D;

     n_coef_layer2 = alloc_array1_i(n_layers);
     coefs2        = alloc_array4_d(40, 5100, 4, 4);
     omega2        = alloc_array1_d(n_layers);
     ltau2         = alloc_array1_d(n_layers);

     init_array4_d(coefs2, 40, 5100, 4, 4, 0.);

     for (i = 0; i < n_layers; ++i) {
          n_coef_layer2[n_layers - i - 1] = n_coef_layer[i];

          xrtm_to_polrad_coefs(n_coef_layer[i], n_stokes, coefs[i], coefs2[n_layers - i - 1]);

          omega2[n_layers - i - 1] = omega[i];
          ltau2 [n_layers - i - 1] = ltau [i];
     }

     call_polrad_f_(&n_four, &n_quad, &n_stokes, &n_layers, qx, &F_0, &theta_0, &phi_0, &n_ulevels, umus, &n_umus, phis, &n_phis, &albedo, n_coef_layer2, ***coefs2, omega2, ltau2, ***I_p,***I_m, &delta_m, &fourier_tol, &info, &n_mus2);

     if (info) {
          fprintf(stderr, "ERROR: call_polrad_f()\n");
          return -1;
     }

     free_array1_i(n_coef_layer2);
     free_array4_d(coefs2);
     free_array1_d(omega);
     free_array1_d(ltau);

     return 0;
}



/*******************************************************************************
 *
 ******************************************************************************/
#ifdef __cplusplus
extern "C" {
#endif
void call_radiant_f_(int *n_four,
                     int *n_elem,
                     int *n_coef,
                     int *n_quad,
                     int *n_derivs,
                     int *n_layers,
                     double *qx,
                     double *lambda,
                     double *F_0,
                     double *theta_0,
                     double *phi_0,
                     double *utaus,
                     int *n_ulevels,
                     double *mus,
                     int *n_mus,
                     double *phis,
                     int *n_phis,
                     double *top_t,
                     double *planet_r,
                     double *levels_z,
                     double *levels_t,
                     double *surface_t,
                     int *n_kernels,
                     int *n_kernel_quad,
                     int *kernels,
                     double *ampfac,
                     double *ampfac_l,
                     double *params,
                     double *params_l,
                     int *n_coef_layer,
                     double *coefs,
                     double *coefs_l,
                     double *omega,
                     double *omega_l,
                     double *ltau,
                     double *ltau_l,
                     double *I_p,
                     double *I_m,
                     double *I_p_l,
                     double *I_m_l,
                     double *flux_p,
                     double *flux_m,
                     int *quadrature,
                     int *delta_m,
                     int *n_t_tms,
                     int *psa,
                     int *radiance,
                     int *flux,
                     int *thermal,
                     int *utau_output,
                     uchar *derivs,
                     double *epsilon,
                     int *info,
                     int *n_mus2);
#ifdef __cplusplus
}
#endif

int call_radiant(int n_four, int n_elem, int n_coef, int n_quad, int n_derivs, int n_layers, double *qx, double lambda, double F_0, double mu_0, double phi_0, double *utaus, int n_ulevels, double *umus, int n_umus, double *phis, int n_phis, double top_b, double planet_r, double *levels_z, double *levels_b, double surface_b, int n_kernels, int n_kernel_quad, enum xrtm_kernel_type *kernels, double *ampfac, double **ampfac_l, double **params, double ***params_l, int *n_coef_layer, double ***coefs, double ****coefs_l, double *omega, double **omega_l, double *ltau, double **ltau_l, double ****I_p, double ****I_m, double *****I_p_l, double *****I_m_l, double *flux_p, double *flux_m, double fourier_tol, int quad_type_, int delta_m, int n_t_tms, int psa, int radiance, int flux, int thermal, int utau_output, uchar **derivs) {

     uchar *derivs2;

     int i;

     int info;

     int quadrature;

     int n_mus2;

     double theta_0;

     double top_t;
     double *levels_t;
     double surface_t;

     double *coefs_l2;
     double *omega_l2;
     double *ltau_l2;

     double *params2;

     double *ampfac_l2;
     double *params_l2;

     double *I_p2;
     double *I_m2;

     double *I_p_l2;
     double *I_m_l2;


     if (quad_type_ == QUAD_NORM_GAUS_LEG)
          quadrature = 1;
     else
     if (quad_type_ == QUAD_DOUB_GAUS_LEG)
          quadrature = 2;
     else
     if (quad_type_ == QUAD_LOBATTO)
          quadrature = 3;

     theta_0 = acos(mu_0)*R2D;

     if (mu_0 == 1.)
          psa = 0;

     if (thermal) {
          top_t = top_b;

          levels_t = alloc_array1_d(n_layers + 1);
          for (i = 0; i < n_layers + 1; ++i)
               levels_t[i] = levels_b[i];

          surface_t = surface_b;
     }

     if (n_umus == 0)
          n_mus2 = n_quad;
     else
          n_mus2 = n_umus;

     if (radiance) {
          I_p2 = ***I_p;
          I_m2 = ***I_m;
     }

     if (n_kernels > 0) {
          params2   = *params;
     }

     if (n_derivs > 0) {
          coefs_l2  = ***coefs_l;
          omega_l2  = *omega_l;
          ltau_l2   = *ltau_l;
          derivs2   = *derivs;

          if (radiance) {
               I_p_l2    = ****I_p_l;
               I_m_l2    = ****I_m_l;
          }

          if (n_kernels > 0) {
               ampfac_l2 = *ampfac_l;
               params_l2 = **params_l;
          }
     }

     call_radiant_f_(&n_four, &n_elem, &n_coef, &n_quad, &n_derivs, &n_layers, qx, &lambda, &F_0, &theta_0, &phi_0, utaus, &n_ulevels, umus, &n_umus, phis, &n_phis, &top_t, &planet_r, levels_z, levels_t, &surface_t, &n_kernels, &n_kernel_quad, (int *) kernels, ampfac, ampfac_l2, params2, params_l2, n_coef_layer, **coefs, coefs_l2, omega, omega_l2, ltau, ltau_l2, I_p2, I_m2, I_p_l2, I_m_l2, flux_p, flux_m, &quadrature, &delta_m, &n_t_tms, &psa, &radiance, &flux, &thermal, &utau_output, derivs2, &fourier_tol, &info, &n_mus2);

     if (info) {
          fprintf(stderr, "ERROR: call_radiant_f()\n");
          return -1;
     }

     if (thermal)
          free_array1_d(levels_t);

     return 0;
}



/*******************************************************************************
 *
 ******************************************************************************/
#ifdef __cplusplus
extern "C" {
#endif
void radtran_(int    *nstokes,
              int    *nummu,
              int    *aziorder,
              double *max_delta_tau,
              int    *src_code,
              char   *quad_type,
              char   *extra_angles,
              char   *deltam,
              double *direct_flux,
              double *direct_mu,
              double *ground_temp,
              char   *ground_type,
              double *ground_albedo,
              double *ground_index,
              double *sky_temp,
              double *wavelength,
              int    *num_layers,
              double *height,
              double *temperatures,
              double *gas_extinct,
              int    *numlegen,
              double *legendre_coef,
              double *part_extinct,
              double *scat_extinct,
              int    *noutlevels,
              int    *outlevels,
              double *mu_values,
              double *up_flux,
              double *down_flux,
              double *up_rad,
              double *down_rad);
#ifdef __cplusplus
}
#endif

int call_radtran3(int n_four, int n_elem, int n_coef, int n_quad, int n_stokes, int n_layers, double *qx, double lambda, double F_0, double mu_0, double phi_0, int *ulevels, int n_ulevels, double *umus, int n_umus, double *phis, int n_phis, double top_b, double *levels_b, double surface_b, double albedo, int *n_coef_layer, double ***coefs, double *omega, double *ltau, double ****I_p, double ****I_m, double *flux_p, double *flux_m, double d_tau, int quad_type_, int delta_m, int radiance, int flux, int solar, int thermal) {

     int i;
     int j;
     int k;
     int kk;
     int l;
     int m;
     int n;
     int nn;

     int i_mus2;
     int n_mus2;

     double a;

     int    nstokes         = n_stokes;
     int    nummu           = n_quad + n_umus;
     int    aziorder;
     double max_delta_tau   = d_tau;
     int    src_code;
     char   quad_type;
     char   extra_angles;
     char   deltam;
     double direct_flux;
     double direct_mu;
     double ground_temp;
     char   ground_type     = 'L';
     double ground_albedo   = albedo;
     double ground_index[2] = {0., 0.};
     double sky_temp;
     double wavelength      = lambda;
     int    num_layers      = n_layers;
     double *height;
     double *temperatures;
     double *gas_extinct;
     int    numlegen        = n_coef - 1;
     double ***legendre_coef;
     double *part_extinct;
     double *scat_extinct;
     int    noutlevels      = n_ulevels;
     int    *outlevels;
     double *mu_values;
     double **up_flux;
     double **down_flux;
     double ****up_rad;
     double ****down_rad;

     if (F_0 == 0.)
          n_four = 1;
     aziorder = n_four - 1;

     height        = alloc_array1_d(num_layers + 1);
     temperatures  = alloc_array1_d(num_layers + 1);
     gas_extinct   = alloc_array1_d(num_layers);
     legendre_coef = alloc_array3_d(num_layers, n_coef + 1, 6);
     part_extinct  = alloc_array1_d(num_layers);
     scat_extinct  = alloc_array1_d(num_layers);
     outlevels     = alloc_array1_i(noutlevels);
     mu_values     = alloc_array1_d(nummu);
     up_flux       = alloc_array2_d(noutlevels, nstokes);
     down_flux     = alloc_array2_d(noutlevels, nstokes);
     up_rad        = alloc_array4_d(noutlevels, aziorder + 1, nummu, nstokes);
     down_rad      = alloc_array4_d(noutlevels, aziorder + 1, nummu, nstokes);

     if (quad_type_ == QUAD_NORM_GAUS_LEG)
          quad_type = 'G';
     else
     if (quad_type_ == QUAD_DOUB_GAUS_LEG)
          quad_type = 'D';
     else
     if (quad_type_ == QUAD_LOBATTO)
          quad_type = 'L';
#ifdef DEBUG
     else {
         fprintf(stderr, "ERROR: call_radtran3(): end of if / else if\n");
         exit(1);
     }
#endif
     if (! delta_m) {
          deltam = 'N';
          n = n_coef;
     }
     else {
          deltam = 'Y';
          n = n_coef + 1;
     }

     if (! solar && ! thermal) {
          src_code    = 0;
          ground_temp = 0.;
          sky_temp    = 0.;
     }
     else
     if (solar && ! thermal) {
          src_code    = 1;
          direct_flux = F_0 * mu_0;
          direct_mu   = mu_0;
          ground_temp = 0.;
          sky_temp    = 0.;
     }
     else
     if (! solar && thermal) {
          src_code    = 2;
          ground_temp = surface_b;
          sky_temp    = top_b;
     }
     else {
          src_code    = 3;
          direct_flux = F_0 * mu_0;
          direct_mu   = mu_0;
          ground_temp = surface_b;
          sky_temp    = top_b;
     }

     height[0] = 0.;
     if (thermal)
          temperatures[0] = levels_b[0];
     for (i = 0; i < n_layers; ++i) {
          height[i+1] = height[i] + ltau[i];
          if (! thermal)
               temperatures[i+1] = 0.;
          else
               temperatures[i+1] = levels_b[i+1];

          gas_extinct [i] = 0.;
          part_extinct[i] = 1.;
          scat_extinct[i] = omega[i];

          nn = n <= n_coef_layer[i] ? n : n_coef_layer[i];

          for (j = 0; j < nn; ++j) {
               legendre_coef[i][j][0] = coefs[i][0][j];
               if (n_elem > 1) {
                    legendre_coef[i][j][4] = coefs[i][1][j];
                    legendre_coef[i][j][2] = coefs[i][2][j];
                    legendre_coef[i][j][5] = coefs[i][3][j];
                    legendre_coef[i][j][1] = coefs[i][4][j];
                    legendre_coef[i][j][3] = coefs[i][5][j];
               }
               else {
                    legendre_coef[i][j][4] = 0.;
                    legendre_coef[i][j][2] = 0.;
                    legendre_coef[i][j][5] = 0.;
                    legendre_coef[i][j][1] = 0.;
                    legendre_coef[i][j][3] = 0.;
               }
          }

          memset(&legendre_coef[i][j][0], 0, (n_coef - j + 1) * 6 * sizeof(double));
     }

     if (n_umus == 0)
          extra_angles = 'N';
     else {
          extra_angles = 'Y';

          for (i = 0; i < n_quad; ++i)
               mu_values[i] = 0.;

          for ( ; i < n_quad + n_umus; ++i)
               mu_values[i] = umus[i - n_quad];
     }

     for (i = 0; i < n_ulevels; ++i)
          outlevels[i] = ulevels[i] + 1;

     radtran_(&nstokes,
              &nummu,
              &aziorder,
              &max_delta_tau,
              &src_code,
              &quad_type,
              &extra_angles,
              &deltam,
              &direct_flux,
              &direct_mu,
              &ground_temp,
              &ground_type,
              &ground_albedo,
              ground_index,
              &sky_temp,
              &wavelength,
              &num_layers,
              height,
              temperatures,
              gas_extinct,
              &numlegen,
              **legendre_coef,
              part_extinct,
              scat_extinct,
              &noutlevels,
              outlevels,
              mu_values,
              *up_flux,
              *down_flux,
              ***up_rad,
              ***down_rad);

     if (n_umus == 0) {
          i_mus2 = 0;
          n_mus2 = n_quad;
     }
     else {
          i_mus2 = n_quad;
          n_mus2 = n_umus;
     }

     if (radiance) {
          for (i = 0; i < n_ulevels; ++i) {
               for (j = 0; j < n_phis; ++j) {
                    for (k = i_mus2, kk = 0; k < i_mus2 + n_mus2; ++k, ++kk) {
                         for (l = 0; l < 2 && l < n_stokes; ++l) {
                              I_p[i][kk][j][l] = 0.;
                              I_m[i][kk][j][l] = 0.;
                              for (m = 0; m < n_four; ++m) {
                                   a = cos(m * (phis[j] - phi_0)*D2R);
                                   I_p[i][kk][j][l] += a * up_rad  [i][m][k][l];
                                   I_m[i][kk][j][l] += a * down_rad[i][m][k][l];
                              }
                         }
                         for (     ; l < 4 && l < n_stokes; ++l) {
                              I_p[i][kk][j][l] = 0.;
                              I_m[i][kk][j][l] = 0.;
                              for (m = 0; m < n_four; ++m) {
                                   a = sin(m * (phis[j] - phi_0)*D2R);
                                   I_p[i][kk][j][l] += a * up_rad  [i][m][k][l];
                                   I_m[i][kk][j][l] += a * down_rad[i][m][k][l];
                              }
                         }
                    }
               }
          }
     }

     if (flux) {
          for (i = 0; i < n_ulevels; ++i) {
               for (j = 0; j < n_stokes; ++j) {
                    flux_p[i] = up_flux[i][0];
                    flux_m[i] = down_flux[i][0];
               }
          }

     }

     free_array1_d(height);
     free_array1_d(temperatures);
     free_array1_d(gas_extinct);
     free_array3_d(legendre_coef);
     free_array1_d(part_extinct);
     free_array1_d(scat_extinct);
     free_array1_i(outlevels);
     free_array1_d(mu_values);
     free_array2_d(up_flux);
     free_array2_d(down_flux);
     free_array4_d(up_rad);
     free_array4_d(down_rad);

     return 0;
}



/*******************************************************************************
 *
 ******************************************************************************/
#ifdef __cplusplus
extern "C" {
#endif
void call_soi_f_(int *n_four,
                 int *n_elem,
                 int *n_coef,
                 int *n_quad,
                 int *n_layers,
                 double *qx,
                 double *F_0,
                 double *theta_0,
                 double *phi_0,
                 double *umus,
                 int *n_umus,
                 double *phis,
                 int *n_phis,
                 double *planet_r,
                 double *levels_z,
                 int *n_kernels,
                 int *n_kernel_quad,
                 int *kernels,
                 double *ampfac,
                 double *params,
                 int *n_coef_layer,
                 double *coefs,
                 double *omega,
                 double *ltau,
                 double *I_p,
                 double *I_m,
                 double *d_tau,
                 int *delta_m,
                 int *n_t_tms,
                 int *psa,
                 double *epsilon,
                 int *info,
                 int *n_mus2);
#ifdef __cplusplus
}
#endif

int call_soi(int n_four, int n_elem, int n_coef, int n_quad, int n_layers, double *qx, double F_0, double mu_0, double phi_0, double *umus, int n_umus, double *phis, int n_phis, double planet_r, double *levels_z, int n_kernels, int n_kernel_quad, enum xrtm_kernel_type *kernels, double *ampfac, double **params, int *n_coef_layer, double ***coefs, double *omega, double *ltau, double ****I_p, double ****I_m, double d_tau, double fourier_tol, int delta_m, int n_t_tms, int psa) {

     int info;

     int n_mus2;

     double theta_0;

     double *params2;

     theta_0 = acos(mu_0)*R2D;

     if (n_umus == 0)
          n_mus2 = n_quad;
     else
          n_mus2 = n_umus;

     if (n_kernels > 0)
          params2 = *params;

     call_soi_f_(&n_four, &n_elem, &n_coef, &n_quad, &n_layers, qx, &F_0, &theta_0, &phi_0, umus, &n_umus, phis, &n_phis, &planet_r, levels_z, &n_kernels, &n_kernel_quad, (int *) kernels, ampfac, params2, n_coef_layer, **coefs, omega, ltau, ***I_p, ***I_m, &d_tau, &delta_m, &n_t_tms, &psa, &fourier_tol, &info, &n_mus2);

     if (info) {
          fprintf(stderr, "ERROR: call_soi_f()\n");
          return -1;
     }

     return 0;
}



/*******************************************************************************
 *
 ******************************************************************************/
#ifdef __cplusplus
extern "C" {
#endif
void call_2stream_f_(int *n_four,
                     int *n_elem,
                     int *n_coef,
                     int *n_quad,
                     int *n_derivs,
                     int *n_layers,
                     double *qx,
                     double *lambda,
                     double *F_0,
                     double *theta_0,
                     double *phi_0,
                     double *mus,
                     int *n_mus,
                     double *phis,
                     int *n_phis,
                     double *top_t,
                     double *planet_r,
                     double *levels_z,
                     double *levels_t,
                     double *surface_t,
                     int *n_kernels,
                     int *n_kernel_quad,
                     int *kernels,
                     double *ampfac,
                     double *ampfac_l,
                     double *params,
                     double *params_l,
                     double *g,
                     double *g_l,
                     double *omega,
                     double *omega_l,
                     double *ltau,
                     double *ltau_l,
                     double *I_p,
                     double *I_m,
                     double *I_p_l,
                     double *I_m_l,
                     double *flux_p,
                     double *flux_m,
                     int *quadrature,
                     int *delta_m,
                     int *n_t_tms,
                     int *psa,
                     int *radiance,
                     int *flux,
                     int *thermal,
                     uchar *derivs,
                     double *epsilon,
                     int *info,
                     int *n_mus2);
#ifdef __cplusplus
}
#endif



/*******************************************************************************
 *
 ******************************************************************************/
#ifdef __cplusplus
extern "C" {
#endif
void twostr_(double *albedo,
             double *btemp,
             int *deltam,
             double *dtauc,
             double *fbeam,
             double *fisot,
             double *gg,
             char *header,
             int *ierror,
             int *maxcly,
             int *maxulv,
             int *nlyr,
             int *planck,
             int *ntau,
             int *prnt,
             int *quiet,
             double *radius,
             int *spher,
             double *ssalb,
             double *temis,
             double *temper,
             double *ttemp,
             double *umu0,
             int *usrtau,
             double *utau,
             double *wvnmlo,
             double *wvnmhi,
             double *zd,
             double *dfdt,
             double *flup,
             double *rfldir,
             double *rfldn,
             double *uavg);
#ifdef __cplusplus
}
#endif

int call_twostr(int n_layers, double lambda, double F_0, double mu_0, double phi_0, int *ulevels, double *utaus, int n_ulevels, double top_t, double planet_r, double *levels_z, double *levels_t, double surface_t, double albedo_, double *g, double *omega, double *ltau, double *mean_p, double *mean_m, double *flux_p, double *flux_m, double *flux_div, int delta_m, int mean, int flux, int div, int psa, int thermal, int utau_output) {

     int i;
     int ii;

     char   header[127];

     int    deltam  = delta_m;
     int    ierror[22];
     int    maxcly  = n_layers;
     int    maxulv  = n_layers + 1;
     int    nlyr    = n_layers;
     int    planck  = thermal;
     int    ntau    = 2;
     int    prnt[2] = {0, 0};
     int    quiet   = 0;
     int    spher   = psa;
     int    usrtau  = 1;

     double ttau;

     double albedo  = albedo_;
     double btemp;
     double *dtauc;
     double fbeam   = F_0;
     double fisot;
     double *gg;
     double radius  = planet_r;
     double *ssalb;
     double temis   = 0.;
     double *temper;
     double ttemp   = 0.;
     double umu0    = mu_0;
     double *utau;
     double wvnmhi = 1. / ((lambda - .5) / 1.e6 * 1.e2);
     double wvnmlo = 1. / ((lambda + .5) / 1.e6 * 1.e2);
     double *zd;
     double *dfdt;
     double *flup;
     double *rfldir;
     double *rfldn;
     double *uavg;

     dtauc  = alloc_array1_d(maxcly);
     gg     = alloc_array1_d(maxcly);
     ssalb  = alloc_array1_d(maxcly);
     temper = alloc_array1_d(maxcly + 1);
     utau   = alloc_array1_d(maxulv);
     zd     = alloc_array1_d(maxcly + 1);
     dfdt   = alloc_array1_d(maxulv);
     flup   = alloc_array1_d(maxulv);
     rfldir = alloc_array1_d(maxulv);
     rfldn  = alloc_array1_d(maxulv);
     uavg   = alloc_array1_d(maxulv);


     memset(header, ' ', 127);


     for (i = 0; i < n_layers; ++i) {
          dtauc[i] = ltau [i];
          gg   [i] = g    [i];
          ssalb[i] = omega[i];
     }

     if (! thermal)
          fisot = 0.;
     else {

          fisot = top_t;
/*
          fisot = planck(lambda * 1.e-6, top_t) * 1.e-6;
*/
          for (i = 0; i < n_layers + 1; ++i) {
               temper[i] = levels_t[i];
/*
               temper[i] = planck(lambda * 1.e-6, levels_t[i]) * 1.e-6;
*/
          }

          btemp = surface_t;
/*
          btemp = planck(lambda * 1.e-6, surface_t) * 1.e-6;
*/
     }

     ntau = n_ulevels;
     if (! utau_output) {
          ii = 0;
          ttau = 0.;
          for (i = 0; i < n_layers + 1; ++i) {
               if (i == ulevels[ii])
                    utau[ii++] = ttau;
               if (ii >= n_ulevels)
                    break;
               ttau += ltau[i];
          }
     }
     else {
          for (i = 0; i < n_ulevels; ++i)
               utau[i] = utaus[i];
     }
if (0) {
     twostr_(&albedo,
             &btemp,
             &deltam,
             dtauc,
             &fbeam,
             &fisot,
             gg,
             header,
             ierror,
             &maxcly,
             &maxulv,
             &nlyr,
             &planck,
             &ntau,
             prnt,
             &quiet,
             &radius,
             &spher,
             ssalb,
             &temis,
             temper,
             &ttemp,
             &umu0,
             &usrtau,
             utau,
             &wvnmlo,
             &wvnmhi,
             zd,
             dfdt,
             flup,
             rfldir,
             rfldn,
             uavg);
}
     for (i = 0; i < 22; ++i) {
          if (ierror[i] != 0) {
               fprintf(stderr, "ERROR: twostr ierror[%d] = %d\n", i, ierror[i]);
               return -1;
          }
     }

     for (i = 0; i < n_ulevels; ++i) {
          if (mean) {
               mean_p[i] = uavg[i];
               mean_m[i] = uavg[i];
          }

          if (flux) {
               flux_p[i] = flup [i];
               flux_m[i] = rfldn[i] + rfldir[i];
          }

          if (div)
               flux_div[i] = dfdt[i];
     }

     free_array1_d(dtauc);
     free_array1_d(gg);
     free_array1_d(ssalb);
     free_array1_d(temper);
     free_array1_d(utau);
     free_array1_d(zd);
     free_array1_d(dfdt);
     free_array1_d(flup);
     free_array1_d(rfldir);
     free_array1_d(rfldn);
     free_array1_d(uavg);

     return 0;
}



/*******************************************************************************
 *
 ******************************************************************************/
#ifdef __cplusplus
extern "C" {
#endif
void call_vlidort_f_(int *n_four,
                     int *n_elem,
                     int *n_coef,
                     int *n_quad,
                     int *n_stokes,
                     int *n_derivs,
                     int *n_layers,
                     double *qx,
                     double *F_0,
                     double *sza,
                     double *phi_0,
                     int *ulevels,
                     double *utaus,
                     int *n_ulevels,
                     double *umus,
                     int *n_umus,
                     double *phis,
                     int *n_phis,
                     double *planet_r,
                     double *levels_z,
                     double *levels_b,
                     double *surface_b,
                     int *n_kernels,
                     int *n_kernel_quad,
                     int *kernels,
                     double *ampfac,
                     double *ampfac_l,
                     double *params,
                     double *params_l,
                     int *n_coef_layer,
                     double *coefs,
                     double *coefs_l,
                     double *omega,
                     double *omega_l,
                     double *ltau,
                     double *ltau_l,
                     double *I_p,
                     double *I_m,
                     double *I_p_l,
                     double *I_m_l,
                     double *mean_p,
                     double *mean_m,
                     double *mean_p_l,
                     double *mean_m_l,
                     double *flux_p,
                     double *flux_m,
                     double *flux_p_l,
                     double *flux_m_l,
                     int *delta_m,
                     int *n_t_tms,
                     int *psa,
                     int *quad_output,
                     int *radiance,
                     int *mean,
                     int *flux,
                     int *thermal,
                     int *utau_output,
                     uchar *derivs,
                     double *epsilon,
                     int *info,
                     int *n_mus2);
#ifdef __cplusplus
}
#endif

int call_vlidort(int n_four, int n_elem, int n_coef, int n_quad, int n_stokes, int n_derivs, int n_layers, double *qx, double lambda, double F_0, double mu_0, double phi_0, int *ulevels, double *utaus, int n_ulevels, double *umus, int n_umus, double *phis, int n_phis, double top_b, double planet_r, double *levels_z, double *levels_b, double surface_b, int n_kernels, int n_kernel_quad, enum xrtm_kernel_type *kernels, double *ampfac, double **ampfac_l, double **params, double ***params_l, int *n_coef_layer, double ***coefs, double ****coefs_l, double *omega, double **omega_l, double *ltau, double **ltau_l, double ****I_p, double ****I_m, double *****I_p_l, double *****I_m_l, double *mean_p, double *mean_m, double **mean_p_l, double **mean_m_l, double *flux_p, double *flux_m, double **flux_p_l, double **flux_m_l, double fourier_tol, int delta_m, int n_t_tms, int psa, int quad_output, int radiance, int mean, int flux, int thermal, int utau_output, uchar **derivs) {

     uchar *derivs2;

     int info;

     int n_mus2;

     double sza;

     double *coefs_l2;
     double *omega_l2;
     double *ltau_l2;

     double *params2;

     double *ampfac_l2;
     double *params_l2;

     double *I_p2;
     double *I_m2;

     double *I_p_l2;
     double *I_m_l2;

     double *mean_p_l2;
     double *mean_m_l2;
     double *flux_p_l2;
     double *flux_m_l2;

     sza = acos(mu_0)*R2D;

     if (mu_0 == 1.)
          psa = 0;

     if (n_umus == 0)
          n_mus2 = n_quad;
     else
          n_mus2 = n_umus;

     if (radiance) {
          I_p2 = ***I_p;
          I_m2 = ***I_m;
     }

     if (n_kernels > 0)
          params2 = *params;

     if (n_derivs > 0) {
          coefs_l2 = ***coefs_l;
          omega_l2 = *omega_l;
          ltau_l2  = *ltau_l;
          derivs2  = *derivs;

          if (radiance) {
               I_p_l2    = ****I_p_l;
               I_m_l2    = ****I_m_l;
          }
          if (mean) {
               mean_p_l2 = *mean_p_l;
               mean_m_l2 = *mean_m_l;
          }
          if (flux) {
               flux_p_l2 = *flux_p_l;
               flux_m_l2 = *flux_m_l;
          }

          if (n_kernels > 0) {
               ampfac_l2 = *ampfac_l;
               params_l2 = **params_l;
          }
     }

     call_vlidort_f_(&n_four, &n_elem, &n_coef, &n_quad, &n_stokes, &n_derivs, &n_layers, qx, &F_0, &sza, &phi_0, ulevels, utaus, &n_ulevels, umus, &n_umus, phis, &n_phis, &planet_r, levels_z, levels_b, &surface_b, &n_kernels, &n_kernel_quad, (int *) kernels, ampfac, ampfac_l2, params2, params_l2, n_coef_layer, **coefs, coefs_l2, omega, omega_l2, ltau, ltau_l2, I_p2, I_m2, I_p_l2, I_m_l2, mean_p, mean_m, mean_p_l2, mean_m_l2, flux_p, flux_m, flux_p_l2, flux_m_l2, &delta_m, &n_t_tms, &psa, &quad_output, &radiance, &mean, &flux, &thermal, &utau_output, derivs2, &fourier_tol, &info, &n_mus2);

     if (info) {
          fprintf(stderr, "ERROR: call_vlidort_f()\n");
          return -1;
     }

     return 0;
}
