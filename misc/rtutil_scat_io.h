/*******************************************************************************
**
**    Copyright (C) 2007-2020 Greg McGarragh <greg.mcgarragh@colostate.edu>
**
**    This source code is licensed under the GNU General Public License (GPL),
**    Version 3.  See the file COPYING for more details.
**
*******************************************************************************/

#ifndef RTUTIL_SCAT_IO_H
#define RTUTIL_SCAT_IO_H

#include <gutil.h>

#include "rtutil_size_dist.h"

#ifdef __cplusplus
extern "C" {
#endif


typedef struct {
     enum size_dist_type dist_type;
     int n_int1;
     int n_int2;
     int n_quad;
     double lambda;
     double mr;
     double mi;
     double a1;
     double a2;
     double a3;
     double a4;
     double a5;
     double r1;
     double r2;
     double accuracy;
     double norm;
     double reff;
     double veff;
     double gavg;
     double vavg;
     double ravg;
     double rvw;
     double cext;
     double csca;
     double cbak;
     double g;
} scat_header_data;


void scat_header_import(scat_header_data *h, double lambda,
                        double mr, double mi, enum size_dist_type dist_type,
                        double a1, double a2, double a3, double a4, double a5,
                        double r1, double r2,
                        int n_int1, int n_int2, int n_quad, double accuracy,
                        double norm, double reff, double veff,
                        double gavg, double vavg, double ravg, double rvw,
                        double cext, double csca, double cbak, double g);
void scat_header_export(scat_header_data *h, double *lambda,
                        double *mr, double *mi, enum size_dist_type *dist_type,
                        double *a1, double *a2, double *a3, double *a4, double *a5,
                        double *r1, double *r2,
                        int *n_int1, int *n_int2, int *n_quad, double *accuracy,
                        double *norm, double *reff, double *veff,
                        double *gavg, double *vavg, double *ravg, double *rvw,
                        double *cext, double *csca, double *cbak, double *g);
void scat_header_zero(scat_header_data *h);

int scat_header_size_bin(void);
int scat_header_read_fn(const char *filename, scat_header_data *h);
int scat_header_read_fp(FILE *fp, scat_header_data *h);
int scat_header_skip(FILE *fp);
int scat_header_write_fn(const char *filename, scat_header_data *h);
int scat_header_write_fp(FILE *fp, scat_header_data *h);
int scat_header_blank(FILE *fp);

int scat_header_read_bin_fn(const char *filename, scat_header_data *h);
int scat_header_read_bin_fp(FILE *fp, scat_header_data *h);
int scat_header_skip_bin(FILE *fp);
int scat_header_write_bin_fn(const char *filename, scat_header_data *h);
int scat_header_write_bin_fp(FILE *fp, scat_header_data *h);
int scat_header_blank_bin(FILE *fp);

int scat_header_interp(double a, scat_header_data *h1, scat_header_data *h2,
                       scat_header_data *h);

int phase_func_read_fn(const char *filename, scat_header_data *h,
                       double **theta, double ***phase, int allsix);
int phase_func_read_fp(FILE *fp, scat_header_data *h,
                       double **theta_, double ***phase_, int allsix);
int phase_func_write_fn(const char *filename, scat_header_data *h,
                        double *theta, double **phase, int n_angles, int allsix);
int phase_func_write_fp(FILE *fp, scat_header_data *h,
                        double *theta, double **phase, int n_angles, int allsix);

int scat_coefs_read_fn(const char *filename, int n_coef,
                        scat_header_data *h, double ***coefs, int allsix);
int scat_coefs_read_fp(FILE *fp, int n_coef,
                        scat_header_data *h, double ***coefs_, int allsix);
int scat_coefs_write_fn(const char *filename, int n_coef2,
                         scat_header_data *h, double **coefs, int allsix);
int scat_coefs_write_fp(FILE *fp, int n_coef2,
                         scat_header_data *h, double **coefs, int allsix);

int scat_coefs_size_bin(int n_coef, int size, int allsix);
int scat_coefs_read_bin_fn(const char *filename, int n_coef,
                           scat_header_data *h, void ***coefs, int size, int n_elem);
int scat_coefs_read_bin_fp(FILE *fp, int n_coef,
                           scat_header_data *h, void ***coefs_, int size, int n_elem);
int scat_coefs_write_bin_fn(const char *filename, int n_coef2,
                            scat_header_data *h, double **coefs, int size, int n_elem);
int scat_coefs_write_bin_fp(FILE *fp, int n_coef2,
                            scat_header_data *h, double **coefs, int size, int n_elem);

int scat_coefs_interp(double a,
                      int n_coef1, scat_header_data *h1, double **coefs1,
                      int n_coef2, scat_header_data *h2, double **coefs2,
                      scat_header_data *h , double **coefs, int allsix);

void scat_coefs_hg(int n_coef, double g, double *coefs);
void scat_coefs_hg_l(int n_coef, double g, double g_l, double *coefs_l);
void scat_coefs_dhg(int n_coef, double g1, double g2, double f, double *coefs);
void scat_coefs_dhg_l(int n_coef, double g1, double g2, double f,
                      double g1_l, double g2_l, double f_l, double *coefs_l);
void scat_coefs_ray(double rho, double **coefs);
void scat_coefs_ray_l(double rho, double rho_l, double **coefs_l);

int scat_total_read_fn(const char *filename, scat_header_data *h, int *n_quad,
                       int *n_four, int *quad_type, double *cextx, double *cexty,
                       double *cscatx, double *cscaty, double *cbackx, double *cbacky,
                       double *gx, double *gy);
int scat_total_read_fp(FILE *fp, scat_header_data *h, int *n_quad,
                       int *n_four, int *quad_type, double *cextx, double *cexty,
                       double *cscatx, double *cscaty, double *cbackx, double *cbacky,
                       double *gx, double *gy);
int scat_total_write_fn(const char *filename, scat_header_data *h, int n_quad,
                        int n_four, int quad_type, double *cextx, double *cexty,
                        double *cscatx, double *cscaty, double *cbackx, double *cbacky,
                        double *gx, double *gy);
int scat_total_write_fp(FILE *fp, scat_header_data *h, int n_quad,
                        int n_four, int quad_type, double *cextx, double *cexty,
                        double *cscatx, double *cscaty, double *cbackx, double *cbacky,
                        double *gx, double *gy);

int phase_mats_read_fn(const char *filename, scat_header_data *h, int *n_quad,
                       int *n_four, int *quad_type, double ******P,
                       double ****Kext, double ***Kemi);
int phase_mats_read_fp(FILE *fp, scat_header_data *h, int *n_quad,
                       int *n_four, int *quad_type, double ******P_,
                       double ****Kext_, double ***Kemi_);
int phase_mats_write_fn(const char *filename, scat_header_data *h, int n_quad,
                        int n_four, int quad_type, double *qx, double *****P,
                        double ***Kext, double **Kemi);
int phase_mats_write_fp(FILE *fp, scat_header_data *h, int n_quad,
                        int n_four, int quad_type, double *qx, double *****P,
                        double ***Kext, double **Kemi);

int user_defined_n_coef(const char *filename);

int load_scat_coefs(const char *filename, int n_coef, double ***coefs_, int *flag);

int load_scat_coefs2(const char *filename, int n_coef, int n_coef2, double ***coefs_, int *flag);


#ifdef __cplusplus
}
#endif

#endif /* RTUTIL_SCAT_IO_H */
