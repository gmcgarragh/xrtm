/*******************************************************************************
**
**    Copyright (C) 2007-2020 Greg McGarragh <greg.mcgarragh@colostate.edu>
**
**    This source code is licensed under the GNU General Public License (GPL),
**    Version 3.  See the file COPYING for more details.
**
*******************************************************************************/

#include <gutil.h>
#include <gindex_name_value.h>
#include <gmath_matrix.h>

#include <rtutil_math.h>

#include "xrtm.h"
#include "xrtm_brdf.h"
#include "xrtm_brdf_a.h"
#include "xrtm_utility.h"


/*******************************************************************************
 *
 ******************************************************************************/
static const char *xrtm_kernel_names[] = {
     "lambertian",
     "direct_and_diffuse",
     "roujean",
     "li_sparse",
     "li_dense",
     "ross_thin",
     "ross_thick",
     "hapke",
     "rahman",
     "cox_munk",
     "user_defined"
};

static long xrtm_kernel_types[] = {
     XRTM_KERNEL_LAMBERTIAN,
     XRTM_KERNEL_DIRECT_AND_DIFFUSE,
     XRTM_KERNEL_ROUJEAN,
     XRTM_KERNEL_LI_SPARSE,
     XRTM_KERNEL_LI_DENSE,
     XRTM_KERNEL_ROSS_THIN,
     XRTM_KERNEL_ROSS_THICK,
     XRTM_KERNEL_HAPKE,
     XRTM_KERNEL_RAHMAN,
     XRTM_KERNEL_COX_MUNK,
     XRTM_KERNEL_USER_DEFINED
};


GINDEX_NAME_VALUE_TEMPLATE(xrtm_kernel, "xrtm kernel", N_XRTM_KERNEL_TYPES)



/*******************************************************************************
 *
 ******************************************************************************/
int kernel_is_valid(enum xrtm_kernel_type type) {

    if (type < 0 || type >= N_XRTM_KERNEL_TYPES) {
         fprintf(stderr, "ERROR: invalid kernel type: %d\n", type);
         return 0;
    }

    return 1;
}



/*******************************************************************************
 *
 ******************************************************************************/
static int n_params_kernel[] = {0, 4, 0, 2, 2, 0, 0, 3, 3, 2, MAX_KERNEL_PARAMS};

int kernel_n_params(enum xrtm_kernel_type type) {

    if (! kernel_is_valid(type)) {
         fprintf(stderr, "ERROR: kernel_is_valid()\n");
         return -1;
    }

    return n_params_kernel[type];
}



/*******************************************************************************
 *
 ******************************************************************************/
static int kernel_needs_fourier[] = {0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1};

int brdf_needs_fourier(enum xrtm_kernel_type *kernels, int n_kernels) {

    int i;

    for (i = 0; i < n_kernels; ++i) {
         if (kernel_needs_fourier[kernels[i]])
              return 1;
    }

    return 0;
}



/*******************************************************************************
 *
 ******************************************************************************/
int brdf_aux_alloc(brdf_aux_data *aux, int n_quad, int n_kernel_quad) {

     aux->cos_theta  = alloc_array1_d(n_quad + 1);
     aux->sin_theta  = alloc_array1_d(n_quad + 1);
     aux->tan_theta  = alloc_array1_d(n_quad + 1);
     aux->tan_theta2 = alloc_array1_d(n_quad + 1);

     if (! (aux->cos_theta && aux->sin_theta && aux->tan_theta && aux->tan_theta2)) {
          fprintf(stderr, "ERROR: alloc_array1_d(%d)\n", n_quad + 1);
          return -1;
     }

     aux->cos_phi    = alloc_array1_d(n_kernel_quad);
     aux->sin_phi    = alloc_array1_d(n_kernel_quad);
     aux->sin_phi2   = alloc_array1_d(n_kernel_quad);

     if (! (aux->cos_phi && aux->sin_phi && aux->sin_phi2)) {
          fprintf(stderr, "ERROR: alloc_array1_d(%d)\n", n_kernel_quad);
          return -1;
     }

     return 0;
}



void brdf_aux_calc_base(brdf_aux_data *aux, int n_quad, double *qx, int n_kernel_quad, double *kernel_qx) {

     int i;

     for (i = 0; i < n_quad; ++i)
          aux->cos_theta[i] = qx[i];

     for (i = 0; i < n_quad; ++i) {
          aux->sin_theta [i] = sqrt(1. - aux->cos_theta[i]*aux->cos_theta[i]);
          aux->tan_theta [i] = aux->sin_theta[i] / aux->cos_theta[i];
          aux->tan_theta2[i] = aux->tan_theta[i] * aux->tan_theta[i];
     }

     for (i = 0; i < n_kernel_quad; ++i) {
          aux->cos_phi [i] = cos(kernel_qx[i]);
          aux->sin_phi [i] = sin(kernel_qx[i]);
          aux->sin_phi2[i] = aux->sin_phi[i]*aux->sin_phi[i];
     }
}



void brdf_aux_calc_part(brdf_aux_data *aux, int i_quad, int n_quad, double *qx, int n_kernel_quad, double *kernel_qx) {

     int i;

     for (i = i_quad; i < n_quad; ++i)
          aux->cos_theta[i] = qx[i];

     for (i = i_quad; i < n_quad; ++i) {
          aux->sin_theta [i] = sqrt(1. - aux->cos_theta[i]*aux->cos_theta[i]);
          aux->tan_theta [i] = aux->sin_theta[i] / aux->cos_theta[i];
          aux->tan_theta2[i] = aux->tan_theta[i] * aux->tan_theta[i];
     }
}



void brdf_aux_calc_mu_0(brdf_aux_data *aux, int n_quad, double mu_0, int n_kernel_quad, double *kernel_qx) {

     int i;

     i = n_quad;

     aux->cos_theta[i] = mu_0;

     aux->sin_theta [i] = sqrt(1. - aux->cos_theta[i]*aux->cos_theta[i]);
     aux->tan_theta [i] = aux->sin_theta[i] / aux->cos_theta[i];
     aux->tan_theta2[i] = aux->tan_theta[i] * aux->tan_theta[i];
}



void brdf_aux_free(brdf_aux_data *aux) {

     free_array1_d(aux->cos_theta);
     free_array1_d(aux->sin_theta);
     free_array1_d(aux->tan_theta);
     free_array1_d(aux->tan_theta2);

     free_array1_d(aux->cos_phi);
     free_array1_d(aux->sin_phi);
     free_array1_d(aux->sin_phi2);
}



/*******************************************************************************
 *
 ******************************************************************************/
typedef union {
     struct /* direct_and_diffuse_data */ {
          double rho_0v;
          double rho_0d;
          double rho_vd;
          double rho_dd;
     } direct_and_diffuse;

     struct /* roujean_aux_data */ {
          double a;
          double b;
          double c;
          double d;
     } roujean;

     struct /* li_aux_data */ {
          double a;
          double b;
          double c;
          double d;
          double e;
          double g;

          double tan_theta_i_p;
          double tan_theta_r_p;

          double tan_theta_i_p2;
          double tan_theta_r_p2;

          double theta_i_p;
          double theta_r_p;

          double cos_theta_i_p;
          double sin_theta_i_p;

          double cos_theta_r_p;
          double sin_theta_r_p;

          double r;
     } li;

     struct /* ross_aux_data */ {
          double a;
          double b;
          double c;
     } ross;

     struct /* hapke_aux_data */ {
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
     } hapke;

     struct /* rahman_aux_data  */ {
          double a;
          double b;
          double c;
          double d;
          double e;
          double g;
          double gamma1_2;
     } rahman;

     struct /* cox_munk_aux_data */ {
          double a;
          double b;
          double c;
          double d;
     } cox_munk;

} kernel_aux_data;



/*******************************************************************************
 *
 ******************************************************************************/
static void lambertian_aux(brdf_aux_data *aux, void *aux2, int i, int j, int n_derivs, uchar *derivs, double *p, double **p_l) {

}



static double lambertian_kernel(brdf_aux_data *aux, void *aux2, enum xrtm_brdf_geometry brdf_geom, int i, int j, int k, double phi, int n_derivs, uchar *derivs, double *p, double **p_l, double *f_l) {

     int ii;

     if (flags_or(derivs, n_derivs)) {
          for (ii = 0; ii < n_derivs; ++ii) {
               if (derivs[ii]) {
                    f_l[ii] = 0.;
               }
          }
     }

     return 1.;
}


/*
static void lambertian_kernel_a(brdf_aux_data *aux, void *aux2, enum xrtm_brdf_geometry brdf_geom, int i, int j, int k, double phi, double *p, double *p_a, double f_a) {

}
*/


/*******************************************************************************
 *
 ******************************************************************************/
static void direct_and_diffuse_aux(brdf_aux_data *aux, void *aux2, int i, int j, int n_derivs, uchar *derivs, double *p, double **p_l) {

}



static double direct_diffuse_param(enum xrtm_brdf_geometry brdf_geom, double *p) {

     switch(brdf_geom) {
          case XRTM_BRDF_GEOMETRY_QQ:
               return p[3];
          case XRTM_BRDF_GEOMETRY_Q0:
               return p[1];
          case XRTM_BRDF_GEOMETRY_UQ:
               return p[2];
          case XRTM_BRDF_GEOMETRY_U0:
               return p[0];
          default:
               fprintf(stderr, "ERROR: invalid brdf geometry: %d\n", brdf_geom);
               exit(1);
     }
}

static double direct_and_diffuse_kernel(brdf_aux_data *aux, void *aux2, enum xrtm_brdf_geometry brdf_geom, int i, int j, int k, double phi, int n_derivs, uchar *derivs, double *p, double **p_l, double *f_l) {

     int ii;

     if (flags_or(derivs, n_derivs)) {
          for (ii = 0; ii < n_derivs; ++ii) {
               if (derivs[ii]) {
                    f_l[ii] = direct_diffuse_param(brdf_geom, p_l[ii]);
               }
          }
     }

     return direct_diffuse_param(brdf_geom, p);
}


/*
static void direct_and_diffuse_kernel_a(brdf_aux_data *aux, void *aux2, enum xrtm_brdf_geometry brdf_geom, int i, int j, int k, double phi, double *p, double *p_a, double f_a) {

}
*/


/*******************************************************************************
 *
 ******************************************************************************/
static void roujean_aux(brdf_aux_data *aux, kernel_aux_data *aux2, int i, int j, int n_derivs, uchar *derivs, double *p, double **p_l) {

     double a;

     aux2->roujean.a = aux->tan_theta [i] + aux->tan_theta [j];
     aux2->roujean.b = aux->tan_theta2[i] + aux->tan_theta2[j];

     a = aux->tan_theta[i] * aux->tan_theta[j];
     aux2->roujean.c = 2. * a;
     aux2->roujean.d = a / (2. * PI);
}



static double roujean_kernel(brdf_aux_data *aux, kernel_aux_data *aux2, enum xrtm_brdf_geometry brdf_geom, int i, int j, int k, double phi, int n_derivs, uchar *derivs, double *p, double **p_l, double *f_l) {

     /* Roujean 1992, Wanner 1995, Spurr 2004 (contains a typo in eq for f) */

     double a = 1.;

     double delta;

     if (phi < 0.) a = -1.;

     delta = (aux2->roujean.a + sqrt(aux2->roujean.b - aux2->roujean.c * aux->cos_phi[k])) / PI;

     return ((PI - a * phi) * aux->cos_phi[k] + aux->sin_phi[k]) * aux2->roujean.d - delta;
}



/*******************************************************************************
 * p[0] = x = b / r (crown shape)
 * p[1] = y = h / b (relative height)
 ******************************************************************************/
static void li_aux(brdf_aux_data *aux, kernel_aux_data *aux2, int i, int j, int n_derivs, uchar *derivs, double *p, double **p_l) {

     aux2->li.tan_theta_i_p  = p[0] * aux->tan_theta[i];
     aux2->li.tan_theta_r_p  = p[0] * aux->tan_theta[j];

     aux2->li.tan_theta_i_p2 = aux2->li.tan_theta_i_p * aux2->li.tan_theta_i_p;
     aux2->li.tan_theta_r_p2 = aux2->li.tan_theta_r_p * aux2->li.tan_theta_r_p;

     aux2->li.theta_i_p      = atan(aux2->li.tan_theta_i_p);
     aux2->li.theta_r_p      = atan(aux2->li.tan_theta_r_p);

     aux2->li.cos_theta_i_p  = cos(aux2->li.theta_i_p);
     aux2->li.sin_theta_i_p  = sin(aux2->li.theta_i_p);

     aux2->li.cos_theta_r_p  = cos(aux2->li.theta_r_p);
     aux2->li.sin_theta_r_p  = sin(aux2->li.theta_r_p);

     aux2->li.a              = aux2->li.cos_theta_i_p * aux2->li.cos_theta_r_p;
     aux2->li.b              = aux2->li.sin_theta_i_p * aux2->li.sin_theta_r_p;

     aux2->li.c              = aux2->li.tan_theta_i_p2 + aux2->li.tan_theta_r_p2;
     aux2->li.d              = aux2->li.tan_theta_i_p  * aux2->li.tan_theta_r_p * 2.;
     aux2->li.e              = aux2->li.tan_theta_i_p2 * aux2->li.tan_theta_r_p2;

     aux2->li.r              = 1. / aux2->li.cos_theta_i_p + 1. / aux2->li.cos_theta_r_p;

     aux2->li.g              = p[1] / aux2->li.r;
}



static void li_common(brdf_aux_data *aux, kernel_aux_data *aux2, enum xrtm_brdf_geometry brdf_geom, int i, int j, int k, double phi, int n_derivs, uchar *derivs, double *p, double **p_l, double *p_, double *q) {

     /* Wanner 1995, Spurr 2004 (contains a typo in eq for d(x)) */

     double cos_ksi_p;
     double d;
     double h;
     double cos_t;

     cos_ksi_p = aux2->li.a + aux2->li.b * aux->cos_phi[k];

    *p_        = (1. + cos_ksi_p) / aux2->li.cos_theta_r_p;

     d         = sqrt(aux2->li.c - aux2->li.d * aux->cos_phi[k]);

     h         = sqrt(d * d + aux2->li.e * aux->sin_phi2[k]);

     cos_t     = aux2->li.g * h;

     if (cos_t > 1.) /* overlap area = 0 */
          *q = 1.;
     else
          *q = 1. - (acos(cos_t) - sqrt(1. - cos_t * cos_t) * cos_t) / PI ;
}




static double li_sparse_kernel(brdf_aux_data *aux, kernel_aux_data *aux2, enum xrtm_brdf_geometry brdf_geom, int i, int j, int k, double phi, int n_derivs, uchar *derivs, double *p, double **p_l, double *f_l) {

     /* Wanner 1995, Spurr 2004 (contains a typo) */

     double p_;
     double q;

     li_common(aux, aux2, brdf_geom, i, j, k, phi, n_derivs, derivs, p, p_l, &p_, &q);

     return .5 * p_ - q * aux2->li.r;
}



static double li_dense_kernel(brdf_aux_data *aux, kernel_aux_data *aux2, enum xrtm_brdf_geometry brdf_geom, int i, int j, int k, double phi, int n_derivs, uchar *derivs, double *p, double **p_l, double *f_l) {

     /* Wanner 1995, Spurr 2004 (contains a typo) */

     double p_;
     double q;

     li_common(aux, aux2, brdf_geom, i, j, k, phi, n_derivs, derivs, p, p_l, &p_, &q);

     return p_ / (q * aux2->li.r) - 2.;
}



/*******************************************************************************
 *
 ******************************************************************************/
static void ross_aux(brdf_aux_data *aux, kernel_aux_data *aux2, int i, int j, int n_derivs, uchar *derivs, double *p, double **p_l) {

     aux2->ross.a = aux->cos_theta[i] * aux->cos_theta[j];
     aux2->ross.b = aux->sin_theta[i] * aux->sin_theta[j];
     aux2->ross.c = aux->cos_theta[i] + aux->cos_theta[j];
}



static double ross_thin_kernel(brdf_aux_data *aux, kernel_aux_data *aux2, enum xrtm_brdf_geometry brdf_geom, int i, int j, int k, double phi, int n_derivs, uchar *derivs, double *p, double **p_l, double *f_l) {

     /* Wanner 1995, Spurr 2004 */

     double cos_ksi;
     double ksi;

     cos_ksi = aux2->ross.a + aux2->ross.b * aux->cos_phi[k];

     ksi     = acos(cos_ksi);

     return ((PI / 2. - ksi) * cos_ksi + sin(ksi)) / aux2->ross.a - PI / 2.;
}



static double ross_thick_kernel(brdf_aux_data *aux, kernel_aux_data *aux2, enum xrtm_brdf_geometry brdf_geom, int i, int j, int k, double phi, int n_derivs, uchar *derivs, double *p, double **p_l, double *f_l) {

     /* Wanner 1995, Spurr 2004 */

     double cos_ksi;
     double ksi;

     cos_ksi     = aux2->ross.a + aux2->ross.b * aux->cos_phi[k];

     ksi         = acos(cos_ksi);

     return ((PI / 2. - ksi) * cos_ksi + sin(ksi)) / aux2->ross.c - PI / 4.;
}



/*******************************************************************************
 * p[0] = Sigma (single scattering albedo)
 * p[1] = B_0   (amplitude of hot-spot)
 * p[2] = Delta (angular width of hot-spot)
 ******************************************************************************/
static void hapke_aux(brdf_aux_data *aux, kernel_aux_data *aux2, int i, int j, int n_derivs, uchar *derivs, double *p, double **p_l) {

     double a;
     double b;

     double gamma;

     a = (aux->cos_theta[i] + aux->cos_theta[j]) * 4.;
     aux2->hapke.R_ir = p[0] / a;

     aux2->hapke.a    =  aux->cos_theta[i] * aux->cos_theta[j];
     aux2->hapke.b    =  aux->sin_theta[i] * aux->sin_theta[j];

     aux2->hapke.c    = p[1] * p[2];

     gamma      = sqrt(1. - p[0]);

     a          = 2. * aux->cos_theta[i];
     b          = 2. * aux->cos_theta[j];

     aux2->hapke.d    = ((1. + a) / (1. + a * gamma)) * ((1. + b) / (1. + b * gamma));
}



static double hapke_kernel(brdf_aux_data *aux, kernel_aux_data *aux2, enum xrtm_brdf_geometry brdf_geom, int i, int j, int k, double phi, int n_derivs, uchar *derivs, double *p, double **p_l, double *f_l) {

     /* Spurr 2004 (contains a typo) */

     double cos_ksi;
     double ksi;
     double P;
     double B;

     cos_ksi = aux2->hapke.a + aux2->hapke.b * aux->cos_phi[k];

     P       = 1. + .5 * cos_ksi;

     ksi     = acos(cos_ksi);

     B       = aux2->hapke.c / (p[2] + tan(ksi / 2.));

     return aux2->hapke.R_ir * ((1. + B) * P + aux2->hapke.d - 1.);
}



/*******************************************************************************
 * p[0] = gamma0 (amplitude)
 * p[1] = gamma1 (asymmetry parameter)
 * p[2] = gamma2 (angular spread)
 ******************************************************************************/
static void rahman_aux(brdf_aux_data *aux, kernel_aux_data *aux2, int i, int j, int n_derivs, uchar *derivs, double *p, double **p_l) {

     aux2->rahman.a        = aux->cos_theta[i] * aux->cos_theta[j];
     aux2->rahman.b        = aux->sin_theta[i] * aux->sin_theta[j];

     aux2->rahman.gamma1_2 = p[1]*p[1];

     aux2->rahman.c        = 2. * p[1];

     aux2->rahman.d        = aux->tan_theta2[i] + aux->tan_theta2[j];
     aux2->rahman.e        = aux->tan_theta [i] * aux->tan_theta [j] * 2.;

     aux2->rahman.g        = aux2->rahman.a * (aux->cos_theta[i] + aux->cos_theta[j]);
}



static double rahman_kernel(brdf_aux_data *aux, kernel_aux_data *aux2, enum xrtm_brdf_geometry brdf_geom, int i, int j, int k, double phi, int n_derivs, uchar *derivs, double *p, double **p_l, double *f_l) {

     /* Rahman 1993b, Spurr 2004 (contains 2 errors) */

     double cos_ksi;
     double P;
     double delta;
     double R;

     cos_ksi     = aux2->rahman.a + aux2->rahman.b * aux->cos_phi[k];

     P           = (1. - aux2->rahman.gamma1_2) / pow((1. + aux2->rahman.gamma1_2 + aux2->rahman.c * cos_ksi), 1.5);

     delta       = sqrt(aux2->rahman.d - aux2->rahman.e * aux->cos_phi[k]);

     R           = (1. - p[0]) / (1. + delta);

     return p[0] * P * (1. + R) * pow(aux2->rahman.g, p[2] - 1.);
}



/*******************************************************************************
 * p[0] = sigma^2 = alpha + beta * W (W = wind speed in m / s)
 * p[1] = n = mr (real part of the complex index of refraction)
 ******************************************************************************/
static void cox_munk_aux(brdf_aux_data *aux, kernel_aux_data *aux2, int i, int j, int n_derivs, uchar *derivs, double *p, double **p_l) {

     aux2->cox_munk.a =  aux->cos_theta[i] * aux->cos_theta[j];
     aux2->cox_munk.b =  aux->sin_theta[i] * aux->sin_theta[j];
     aux2->cox_munk.c = (aux->cos_theta[i] + aux->cos_theta[j]) / 2.;
     aux2->cox_munk.d =  4. * aux2->cox_munk.a * PI * p[0];
/*
     aux2->cox_munk.d =  4. * aux2->cox_munk.a *      p[0];
*/
}



static double cox_munk_kernel(brdf_aux_data *aux, kernel_aux_data *aux2, enum xrtm_brdf_geometry brdf_geom, int i, int j, int k, double phi, int n_derivs, uchar *derivs, double *p, double **p_l, double *f_l) {

     double a;
     double cos_ksi;
     double lambda;
     double r_p;
     double r_m;
     double c;
     double R;
     double gamma;
     double gamma2;
     double tau;
     double P;

     cos_ksi     = aux2->cox_munk.a + aux2->cox_munk.b *  aux->cos_phi[k];
/*
     cos_ksi     = aux2->cox_munk.a + aux2->cox_munk.b * -aux->cos_phi[k];
*/
     lambda      = cos(acos(cos_ksi) / 2.);

     a           = p[1] * lambda;

     c           = sqrt(p[1] + (lambda*lambda - 1.));

     r_p         = (a - c) / (a + c);
     r_m         = (lambda - c) / (lambda + c);

     R           = (r_p*r_p + r_m*r_m) / 2.;

     gamma       = aux2->cox_munk.c / lambda;
     gamma2      = gamma * gamma;

     tau         = tan(PI / 2. - asin(gamma));

     P           = 1. / (aux2->cox_munk.d * gamma2*gamma2) * exp(-tau*tau / p[0]);
/*
     P           = 1. / (aux2->cox_munk.d * gamma2*gamma2) * exp(-tau*tau / p[0]);
*/
     return R * P;
}



/*******************************************************************************
 *
 ******************************************************************************/
static double call_aux_func(brdf_aux_data *aux, kernel_aux_data *aux2, int kernel, int i, int j, int n_derivs, uchar *derivs, double *p, double **p_l, brdf_kernel_func_data *func) {

     switch(kernel) {
          case XRTM_KERNEL_LAMBERTIAN:
               lambertian_aux(aux, aux2, i, j, n_derivs, derivs, p, p_l);
               break;
          case XRTM_KERNEL_DIRECT_AND_DIFFUSE:
               direct_and_diffuse_aux(aux, aux2, i, j, n_derivs, derivs, p, p_l);
               break;
          case XRTM_KERNEL_ROUJEAN:
               roujean_aux   (aux, aux2, i, j, n_derivs, derivs, p, p_l);
               break;
          case XRTM_KERNEL_LI_SPARSE:
          case XRTM_KERNEL_LI_DENSE:
               li_aux        (aux, aux2, i, j, n_derivs, derivs, p, p_l);
               break;
          case XRTM_KERNEL_ROSS_THIN:
          case XRTM_KERNEL_ROSS_THICK:
               ross_aux      (aux, aux2, i, j, n_derivs, derivs, p, p_l);
               break;
          case XRTM_KERNEL_HAPKE:
               hapke_aux     (aux, aux2, i, j, n_derivs, derivs, p, p_l);
               break;
          case XRTM_KERNEL_RAHMAN:
               rahman_aux    (aux, aux2, i, j, n_derivs, derivs, p, p_l);
               break;
          case XRTM_KERNEL_COX_MUNK:
               cox_munk_aux  (aux, aux2, i, j, n_derivs, derivs, p, p_l);
               break;
          case XRTM_KERNEL_USER_DEFINED:
               func->aux     (aux, aux2, i, j, n_derivs, derivs, p, p_l);
               break;
          default:
               fprintf(stderr, "ERROR: invalid brdf kernel: %d\n", kernel);
               exit(1);
     }

     return 0.;
}



/*******************************************************************************
 *
 ******************************************************************************/
static double call_kernel_func(brdf_aux_data *aux, kernel_aux_data *aux2, enum xrtm_brdf_geometry brdf_geom, int kernel, int i, int j, int k, double phi, int n_derivs, uchar *derivs, double *p, double **p_l, brdf_kernel_func_data *func, double *f_l) {

     switch(kernel) {
          case XRTM_KERNEL_LAMBERTIAN:
               return lambertian_kernel(aux, aux2, brdf_geom, i, j, k, phi, n_derivs, derivs, p, p_l, f_l);
          case XRTM_KERNEL_DIRECT_AND_DIFFUSE:
               return direct_and_diffuse_kernel(aux, aux2, brdf_geom, i, j, k, phi, n_derivs, derivs, p, p_l, f_l);
          case XRTM_KERNEL_ROUJEAN:
               return roujean_kernel   (aux, aux2, brdf_geom, i, j, k, phi, n_derivs, derivs, p, p_l, f_l);
          case XRTM_KERNEL_LI_SPARSE:
               return li_sparse_kernel (aux, aux2, brdf_geom, i, j, k, phi, n_derivs, derivs, p, p_l, f_l);
          case XRTM_KERNEL_LI_DENSE:
               return li_dense_kernel  (aux, aux2, brdf_geom, i, j, k, phi, n_derivs, derivs, p, p_l, f_l);
          case XRTM_KERNEL_ROSS_THIN:
               return ross_thin_kernel (aux, aux2, brdf_geom, i, j, k, phi, n_derivs, derivs, p, p_l, f_l);
          case XRTM_KERNEL_ROSS_THICK:
               return ross_thick_kernel(aux, aux2, brdf_geom, i, j, k, phi, n_derivs, derivs, p, p_l, f_l);
          case XRTM_KERNEL_HAPKE:
               return hapke_kernel     (aux, aux2, brdf_geom, i, j, k, phi, n_derivs, derivs, p, p_l, f_l);
          case XRTM_KERNEL_RAHMAN:
               return rahman_kernel    (aux, aux2, brdf_geom, i, j, k, phi, n_derivs, derivs, p, p_l, f_l);
          case XRTM_KERNEL_COX_MUNK:
               return cox_munk_kernel  (aux, aux2, brdf_geom, i, j, k, phi, n_derivs, derivs, p, p_l, f_l);
          case XRTM_KERNEL_USER_DEFINED:
               return func->kernel     (aux, aux2, brdf_geom, i, j, k, phi, n_derivs, derivs, p, p_l, f_l);
          default:
               fprintf(stderr, "ERROR: invalid brdf kernel: %d\n", kernel);
               exit(1);
     }

     return 0.;
}



/*******************************************************************************
 *
 ******************************************************************************/
void brdf_build_kernel_vecs(int i_offset, int n_quad, int j_offset, int n_stokes, int n_derivs, int n_kernels, int n_kernel_quad, enum xrtm_kernel_type *kernels, double *ampfac, double **ampfac_l, double **params, double ***params_l, brdf_kernel_func_data *func, double *kernel_qx, double *kernel_qw, brdf_aux_data *aux, enum xrtm_brdf_geometry brdf_geom, double **kernel_f, double ***kernel_f_l, uchar *derivs, work_data work) {

     int i;
     int j;
     int k;
     int kk;
     int l;

     int flag;

     double f;
     double *f_l;

     double **params_l2;

     kernel_aux_data aux2;

     init_array2_d(kernel_f, n_kernel_quad * 2 + 1, n_quad, 0.);

     if (flags_or(derivs, n_derivs)) {
          f_l = get_work1(&work, WORK_DDERIVS);

          for (i = 0; i < n_derivs; ++i) {
               if (derivs[i]) {
                    init_array2_d(kernel_f_l[i], n_kernel_quad * 2 + 1, n_quad, 0.);
               }
          }
     }

     for (i = 0; i < n_kernels; ++i) {
          if (flags_or(derivs, n_derivs))
               params_l2 = params_l[i];

          flag = brdf_needs_fourier(kernels + i, 1);

          if (! flag) {
               for (j = 0; j < n_quad; ++j) {
                    f = call_kernel_func(aux, &aux2, brdf_geom, kernels[i], j_offset, i_offset+j, 0, 0., n_derivs, derivs, params[i], params_l2, &func[i], f_l);
                    kernel_f[0][j] += f * ampfac[i];
                    if (flags_or(derivs, n_derivs)) {
                         for (k = 0; k < n_derivs; ++k) {
                              if (derivs[k]) {
                                   kernel_f_l[k][0][j] += f_l[k] * ampfac[i] + f * ampfac_l[i][k];
                              }
                         }
                    }
               }
          }
          else {
               for (j = 0; j < n_quad; ++j) {
                    call_aux_func(aux, &aux2, kernels[i], j_offset, i_offset+j, n_derivs, derivs, params[i], params_l2, &func[i]);

                    for (k = 0, kk = n_kernel_quad; k < n_kernel_quad; ++k, ++kk) {
                         f = call_kernel_func(aux, &aux2, brdf_geom, kernels[i], j_offset, i_offset+j, k,  kernel_qx[k], n_derivs, derivs, params[i], params_l2, &func[i], f_l);
                         kernel_f[k +1][j] += f * ampfac[i];
                         if (flags_or(derivs, n_derivs)) {
                              for (l = 0; l < n_derivs; ++l) {
                                   if (derivs[l]) {
                                        kernel_f_l[l][k +1][j] += f_l[l] * ampfac[i] + f * ampfac_l[i][l];
                                   }
                              }
                         }

                         f = call_kernel_func(aux, &aux2, brdf_geom, kernels[i], j_offset, i_offset+j, k, -kernel_qx[k], n_derivs, derivs, params[i], params_l2, &func[i], f_l);
                         kernel_f[kk+1][j] += f * ampfac[i];
                         if (flags_or(derivs, n_derivs)) {
                              for (l = 0; l < n_derivs; ++l) {
                                   if (derivs[l]) {
                                        kernel_f_l[l][kk+1][j] += f_l[l] * ampfac[i] + f * ampfac_l[i][l];
                                   }
                              }
                         }
                    }
               }
          }
     }
}



/*******************************************************************************
 *
 ******************************************************************************/
void brdf_build_ref_vec(int i_four, int i_offset, int n_quad, int j_offset, int n_stokes, double qf, double *qx_v, double *qw_v, int n_kernels, int n_kernel_quad, enum xrtm_kernel_type *kernels, double *kernel_qx, double *kernel_qw, double **kernel_f, double *R_s, work_data work) {

     int i;
     int ii;
     int j;
     int jj;

     int flag;

     int n_quad_v;

     double a;

     double *azmfac;

     n_quad_v = n_quad * n_stokes;

     flag = brdf_needs_fourier(kernels, n_kernels);

     a = (2. - (i_four == 0 ? 1. : 0.));

     if (i_four > 0 || n_stokes > 1)
          dvec_zero(R_s, n_quad_v);

     if (! flag) {
          for (i = 0, ii = 0; i < n_quad; ++i, ii += n_stokes) {
               R_s[ii] = kernel_f[0][i];
          }
     }
     else {
          azmfac = get_work_d1(&work, n_kernel_quad * 2);

          for (i = 0, ii = n_kernel_quad; i < n_kernel_quad; ++i, ++ii) {
               azmfac[i ] = kernel_qw[i] * cos(i_four *  kernel_qx[i]) / 2.;
               azmfac[ii] = kernel_qw[i] * cos(i_four * -kernel_qx[i]) / 2.;
          }

          for (i = 0, ii = 0; i < n_quad; ++i, ii += n_stokes) {
               if (i_four == 0)
                    R_s[ii] = kernel_f[0][i];

               for (j = 0, jj = n_kernel_quad; j < n_kernel_quad; ++j, ++jj) {
                    R_s[ii] += kernel_f[j +1][i] * azmfac[j ];
                    R_s[ii] += kernel_f[jj+1][i] * azmfac[jj];
               }
          }
     }

     for (i = 0, ii = 0; i < n_quad; ++i, ii += n_stokes)
          R_s[ii] *= a;
}



/*******************************************************************************
 *
 ******************************************************************************/
void brdf_build_kernel_mats(int i_offset, int n_quad1, int j_offset, int n_quad2, int n_stokes, int n_derivs, int n_kernels, int n_kernel_quad, enum xrtm_kernel_type *kernels, double *ampfac, double **ampfac_l, double **params, double ***params_l, brdf_kernel_func_data *func, double *kernel_qx, double *kernel_qw, brdf_aux_data *aux, enum xrtm_brdf_geometry brdf_geom, double ***kernel_f, double ****kernel_f_l, uchar *derivs, work_data work) {

     int i;
     int j;
     int k;
     int l;
     int ll;
     int m;

     int flag;

     double f;
     double *f_l;

     double **params_l2;

     kernel_aux_data aux2;

     init_array3_d(kernel_f, n_kernel_quad * 2 + 1, n_quad1, n_quad2, 0.);

     if (flags_or(derivs, n_derivs)) {
          f_l  = get_work1(&work, WORK_DDERIVS);

          for (i = 0; i < n_derivs; ++i) {
               if (derivs[i]) {
                    init_array3_d(kernel_f_l[i], n_kernel_quad * 2 + 1, n_quad1, n_quad2, 0.);
               }
          }
     }

     for (i = 0; i < n_kernels; ++i) {
          if (flags_or(derivs, n_derivs))
               params_l2 = params_l[i];

          flag = brdf_needs_fourier(kernels + i, 1);

          if (! flag) {
               for (j = 0; j < n_quad1; ++j) {
                    for (k = 0; k < n_quad2; ++k) {
                         f = call_kernel_func(aux, &aux2, brdf_geom, kernels[i], j_offset+k, i_offset+j, 0, 0., n_derivs, derivs, params[i], params_l2, &func[i], f_l);
                         kernel_f[0][j][k] += f * ampfac[i];
                         if (flags_or(derivs, n_derivs)) {
                              for (l = 0; l < n_derivs; ++l) {
                                   if (derivs[l]) {
                                        kernel_f_l[l][0][j][k] += f_l[l] * ampfac[i] + f * ampfac_l[i][l];
                                   }
                              }
                         }
                    }
               }
          }
          else {
               for (j = 0; j < n_quad1; ++j) {
                    for (k = 0; k < n_quad2; ++k) {
                         call_aux_func(aux, &aux2, kernels[i], j_offset+k, i_offset+j, n_derivs, derivs, params[i], params_l2, &func[i]);

                         for (l = 0, ll = n_kernel_quad; l < n_kernel_quad; ++l, ++ll) {
                              f = call_kernel_func(aux, &aux2, brdf_geom, kernels[i], j_offset+k, i_offset+j, l,  kernel_qx[l], n_derivs, derivs, params[i], params_l2, &func[i], f_l);
                              kernel_f[l +1][j][k] += f * ampfac[i];
                              if (flags_or(derivs, n_derivs)) {
                                   for (m = 0; m < n_derivs; ++m) {
                                        if (derivs[m]) {
                                             kernel_f_l[m][l +1][j][k] += f_l[m] * ampfac[i] + f * ampfac_l[i][m];
                                        }
                                   }
                              }

                              f = call_kernel_func(aux, &aux2, brdf_geom, kernels[i], j_offset+k, i_offset+j, l, -kernel_qx[l], n_derivs, derivs, params[i], params_l2, &func[i], f_l);
                              kernel_f[ll+1][j][k] += f * ampfac[i];
                              if (flags_or(derivs, n_derivs)) {
                                   for (m = 0; m < n_derivs; ++m) {
                                        if (derivs[m]) {
                                             kernel_f_l[m][ll+1][j][k] += f_l[m] * ampfac[i] + f * ampfac_l[i][m];
                                        }
                                   }
                              }
                         }
                    }
               }
          }
     }
}



/*******************************************************************************
 *
 ******************************************************************************/
void brdf_build_ref_mat(int i_four, int i_offset, int n_quad1, int j_offset, int n_quad2, int n_stokes, double qf, double *qx_v, double *qw_v, int n_kernels, int n_kernel_quad, enum xrtm_kernel_type *kernels, double *kernel_qx, double *kernel_qw, double ***kernel_f, double **R_s, work_data work) {

     int ii;
     int i;
     int j;
     int jj;
     int k;
     int kk;

     int flag;

     int jj_offset;

     double a;
     double b;

     double *azmfac;

     flag = brdf_needs_fourier(kernels, n_kernels);

     a = (2. - (i_four == 0 ? 1. : 0.)) * (1. + (i_four == 0 ? 1. : 0.)) * qf;

     if (i_four > 0 || n_stokes > 1)
          dmat_zero(R_s, n_quad1 * n_stokes, n_quad2 * n_stokes);

     if (! flag) {
          for (i = 0, ii = 0; i < n_quad1; ++i, ii += n_stokes) {
               for (j = 0, jj = 0; j < n_quad2; ++j, jj += n_stokes) {
                    R_s[ii][jj] = kernel_f[0][i][j];
               }
          }
     }
     else {
          azmfac = get_work_d1(&work, n_kernel_quad * 2);

          for (i = 0, ii = n_kernel_quad; i < n_kernel_quad; ++i, ++ii) {
               azmfac[i ] = kernel_qw[i] * cos(i_four *  kernel_qx[i]) / 2.;
               azmfac[ii] = kernel_qw[i] * cos(i_four * -kernel_qx[i]) / 2.;
          }

          for (i = 0, ii = 0; i < n_quad1; ++i, ii += n_stokes) {
               for (j = 0, jj = 0; j < n_quad2; ++j, jj += n_stokes) {
                    if (i_four == 0)
                         R_s[ii][jj] = kernel_f[0][i][j];

                    for (k = 0, kk = n_kernel_quad; k < n_kernel_quad; ++k, ++kk) {
                         R_s[ii][jj] += kernel_f[k +1][i][j] * azmfac[k ];
                         R_s[ii][jj] += kernel_f[kk+1][i][j] * azmfac[kk];
                    }
               }
          }
     }

     jj_offset = j_offset * n_stokes;
     for (j = 0, jj = 0; j < n_quad2; ++j, jj += n_stokes) {
          b = a * qx_v[jj_offset + jj] * qw_v[jj_offset + jj];
          for (i = 0, ii = 0; i < n_quad1; ++i, ii += n_stokes) {
               R_s[ii][jj] *= b;
          }
     }
}
