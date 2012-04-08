/******************************************************************************%
**
**    Copyright (C) 2007-2012 Greg McGarragh <gregm@atmos.colostate.edu>
**
**    This source code is licensed under the GNU General Public License (GPL),
**    Version 3.  See the file COPYING for more details.
**
*******************************************************************************/

#include <gutil.h>

#include "rtutil_scat_io.h"
#include "rtutil_support.h"


/*******************************************************************************
 * Phase function (pf1)
 * 	F11
 * 	a1
 *
 * Phase function (pf6)
 * 	F11, F22, F33, F44, F12=F21, F34=F43
 * 	a1,  a2,  a3,  a4,  b1,      b2
 *
 * Legendre coefficients (pm1):
 * 	F11
 *
 * Legendre coefficients (pm6):
 * 	F11, F22, F33, F44, F12=F21, F34=F43
 *
 * Greek (spherical) constants (gc6):
 * 	alpha1, alpha2, alpha3, alpha4, -beta1, -beta2
 * 	beta,   alpha,  zeta,   delta,   gamma,  epsilon
 *
 ******************************************************************************/


/*******************************************************************************
 *
 ******************************************************************************/
void scat_header_import(scat_header_data *h, double lambda,
                        double mr, double mi, enum size_dist_type dist_type,
                        double a1, double b1, double a2, double b2,
                        double gamma, double r1, double r2,
                        int n_int1, int n_int2, int n_quad, double accuracy,
                        double norm, double reff, double veff,
                        double gavg, double vavg, double ravg, double rvw,
                        double cext, double csca, double cbak, double g) {

     h->lambda    = lambda;
     h->mr        = mr;
     h->mi        = mi;
     h->dist_type = dist_type;
     h->a1        = a1;
     h->b1        = b1;
     h->a2        = a2;
     h->b2        = b2;
     h->gamma     = gamma;
     h->r1        = r1;
     h->r2        = r2;
     h->n_int1    = n_int1;
     h->n_int2    = n_int2;
     h->n_quad    = n_quad;
     h->accuracy  = accuracy;
     h->norm      = norm;
     h->reff      = reff;
     h->veff      = veff;
     h->gavg      = gavg;
     h->vavg      = vavg;
     h->ravg      = ravg;
     h->rvw       = rvw;
     h->cext      = cext;
     h->csca      = csca;
     h->cbak      = cbak;
     h->g         = g;
}



void scat_header_export(scat_header_data *h, double *lambda,
                        double *mr, double *mi, enum size_dist_type *dist_type,
                        double *a1, double *b1, double *a2, double *b2,
                        double *gamma, double *r1, double *r2,
                        int *n_int1, int *n_int2, int *n_quad, double *accuracy,
                        double *norm, double *reff, double *veff,
                        double *gavg, double *vavg, double *ravg, double *rvw,
                        double *cext, double *csca, double *cbak, double *g) {

     *lambda    = h->lambda;
     *mr        = h->mr;
     *mi        = h->mi;
     *dist_type = h->dist_type;
     *a1        = h->a1;
     *b1        = h->b1;
     *a2        = h->a2;
     *b2        = h->b2;
     *gamma     = h->gamma;
     *r1        = h->r1;
     *r2        = h->r2;
     *n_int1    = h->n_int1;
     *n_int2    = h->n_int2;
     *n_quad    = h->n_quad;
     *accuracy  = h->accuracy;
     *norm      = h->norm;
     *reff      = h->reff;
     *veff      = h->veff;
     *gavg      = h->gavg;
     *vavg      = h->vavg;
     *ravg      = h->ravg;
     *rvw       = h->rvw;
     *cext      = h->cext;
     *csca      = h->csca;
     *cbak      = h->cbak;
     *g         = h->g;
}



void scat_header_zero(scat_header_data *h) {

     scat_header_import(h, 0., 0., 0., SIZE_DIST_NONE, 0., 0., 0., 0., 0., 0.,
                        0., 0, 0, 0, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.);
}



/*******************************************************************************
 *
 ******************************************************************************/
int scat_header_read_fn(const char *filename, scat_header_data *h) {

     FILE *fp;

     if ((fp = fopen(filename, "r")) == NULL) {
          eprintf("ERROR: Error opening file for reading: %s ... %s\n",
                  filename, strerror(errno));
          return -1;
     }

     if (scat_header_read_fp(fp, h)) {
          eprintf("ERROR: scat_header_read_fp()\n");
          return -1;
     }

     fclose(fp);

     return 0;
}



int scat_header_read_fp(FILE *fp, scat_header_data *h) {

     char *name;

     char format[15 + 20 + 1];

     int max_length;

     if (fscanf(fp, "lambda    = %lf\n", &h->lambda)    < 1) goto L1;
     if (fscanf(fp, "mr        = %lf\n", &h->mr)        < 1) goto L1;
     if (fscanf(fp, "mi        = %lf\n", &h->mi)        < 1) goto L1;

     max_length = size_dist_name_max_length();
     sprintf(format, "dist_type = %%%ds\n", max_length);
     name = (char *) malloc(max_length + 1);
     if (fscanf(fp, format,  name)                      < 1) goto L1;
     h->dist_type = size_dist_code(name);
     free(name);

     if (fscanf(fp, "a1        = %lf\n", &h->a1)        < 1) goto L1;
     if (fscanf(fp, "b1        = %lf\n", &h->b1)        < 1) goto L1;
     if (fscanf(fp, "a2        = %lf\n", &h->a2)        < 1) goto L1;
     if (fscanf(fp, "b2        = %lf\n", &h->b2)        < 1) goto L1;
     if (fscanf(fp, "gamma     = %lf\n", &h->gamma)     < 1) goto L1;
     if (fscanf(fp, "r1        = %lf\n", &h->r1)        < 1) goto L1;
     if (fscanf(fp, "r2        = %lf\n", &h->r2)        < 1) goto L1;
     if (fscanf(fp, "n_int1    = %d\n",  &h->n_int1)    < 1) goto L1;
     if (fscanf(fp, "n_int2    = %d\n",  &h->n_int2)    < 1) goto L1;
     if (fscanf(fp, "n_quad    = %d\n",  &h->n_quad)    < 1) goto L1;
     if (fscanf(fp, "norm      = %lf\n", &h->norm)      < 1) goto L1;
     if (fscanf(fp, "reff      = %lf\n", &h->reff)      < 1) goto L1;
     if (fscanf(fp, "veff      = %lf\n", &h->veff)      < 1) goto L1;
     if (fscanf(fp, "gavg      = %lf\n", &h->gavg)      < 1) goto L1;
     if (fscanf(fp, "vavg      = %lf\n", &h->vavg)      < 1) goto L1;
     if (fscanf(fp, "ravg      = %lf\n", &h->ravg)      < 1) goto L1;
     if (fscanf(fp, "rvw       = %lf\n", &h->rvw)       < 1) goto L1;
     if (fscanf(fp, "cext      = %lf\n", &h->cext)      < 1) goto L1;
     if (fscanf(fp, "csca      = %lf\n", &h->csca)      < 1) goto L1;
     if (fscanf(fp, "cbak      = %lf\n", &h->cbak)      < 1) goto L1;
     if (fscanf(fp, "g         = %lf\n", &h->g)         < 1) goto L1;

     return 0;

L1:  eprintf("ERROR: problem reading header from file: %s\n", strerror(errno));
     return -1;
}



int scat_header_skip(FILE *fp) {

     scat_header_data h;

     if (scat_header_read_fp(fp, &h)) {
          eprintf("ERROR: scat_header_read_fp()\n");
          return -1;
     }

     return 0;
}



int scat_header_write_fn(const char *filename, scat_header_data *h) {

     FILE *fp;

     if ((fp = fopen(filename, "r")) == NULL) {
          eprintf("ERROR: Error opening file for writing: %s ... %s\n",
                  filename, strerror(errno));
          return -1;
     }

     if (scat_header_write_fp(fp, h)) {
          eprintf("ERROR: scat_header_write_fp()\n");
          return -1;
     }

     fclose(fp);

     return 0;
}



int scat_header_write_fp(FILE *fp, scat_header_data *h) {

     const char *name;

     if (fprintf(fp, "lambda    = %e\n", h->lambda)    < 1) goto L1;
     if (fprintf(fp, "mr        = %e\n", h->mr)        < 1) goto L1;
     if (fprintf(fp, "mi        = %e\n", h->mi)        < 1) goto L1;
     name = size_dist_name(h->dist_type);
     if (fprintf(fp, "dist_type = %s\n", name) < 1) goto L1;
     if (fprintf(fp, "a1        = %e\n", h->a1)        < 1) goto L1;
     if (fprintf(fp, "b1        = %e\n", h->b1)        < 1) goto L1;
     if (fprintf(fp, "a2        = %e\n", h->a2)        < 1) goto L1;
     if (fprintf(fp, "b2        = %e\n", h->b2)        < 1) goto L1;
     if (fprintf(fp, "gamma     = %e\n", h->gamma)     < 1) goto L1;
     if (fprintf(fp, "r1        = %e\n", h->r1)        < 1) goto L1;
     if (fprintf(fp, "r2        = %e\n", h->r2)        < 1) goto L1;
     if (fprintf(fp, "n_int1    = %d\n", h->n_int1)    < 1) goto L1;
     if (fprintf(fp, "n_int2    = %d\n", h->n_int2)    < 1) goto L1;
     if (fprintf(fp, "n_quad    = %d\n", h->n_quad)    < 1) goto L1;
     if (fprintf(fp, "norm      = %e\n", h->norm)      < 1) goto L1;
     if (fprintf(fp, "reff      = %e\n", h->reff)      < 1) goto L1;
     if (fprintf(fp, "veff      = %e\n", h->veff)      < 1) goto L1;
     if (fprintf(fp, "gavg      = %e\n", h->gavg)      < 1) goto L1;
     if (fprintf(fp, "vavg      = %e\n", h->vavg)      < 1) goto L1;
     if (fprintf(fp, "ravg      = %e\n", h->ravg)      < 1) goto L1;
     if (fprintf(fp, "rvw       = %e\n", h->rvw)       < 1) goto L1;
     if (fprintf(fp, "cext      = %e\n", h->cext)      < 1) goto L1;
     if (fprintf(fp, "csca      = %e\n", h->csca)      < 1) goto L1;
     if (fprintf(fp, "cbak      = %e\n", h->cbak)      < 1) goto L1;
     if (fprintf(fp, "g         = %e\n", h->g)         < 1) goto L1;

     return 0;

L1:  eprintf("ERROR: problem writing header from file: %s\n", strerror(errno));
     return -1;
}



int scat_header_blank(FILE *fp) {

     scat_header_data h;

     scat_header_zero(&h);

     if (scat_header_write_fp(fp, &h)) {
          eprintf("ERROR: scat_header_write_fp()\n");
          return -1;
     }

     return 0;
}



/*******************************************************************************
 *
 ******************************************************************************/
int scat_header_size_bin() {

     return 1 * sizeof(int) + 20 * sizeof(double);
}



int scat_header_read_bin_fn(const char *filename, scat_header_data *h) {

     FILE *fp;

     if ((fp = fopen(filename, "r")) == NULL) {
          eprintf("ERROR: Error opening file for reading: %s ... %s\n",
                  filename, strerror(errno));
          return -1;
     }

     if (scat_header_read_bin_fp(fp, h)) {
          eprintf("ERROR: scat_header_read_bin_fp()\n");
          return -1;
     }

     fclose(fp);

     return 0;
}



int scat_header_read_bin_fp(FILE *fp, scat_header_data *h) {

     if (fread(&h->lambda,    sizeof(double), 1, fp) < 1) goto L1;
     if (fread(&h->mr,        sizeof(double), 1, fp) < 1) goto L1;
     if (fread(&h->mi,        sizeof(double), 1, fp) < 1) goto L1;
     if (fread(&h->dist_type, sizeof(int),    1, fp) < 1) goto L1;
     if (fread(&h->a1,        sizeof(double), 1, fp) < 1) goto L1;
     if (fread(&h->b1,        sizeof(double), 1, fp) < 1) goto L1;
     if (fread(&h->a2,        sizeof(double), 1, fp) < 1) goto L1;
     if (fread(&h->b2,        sizeof(double), 1, fp) < 1) goto L1;
     if (fread(&h->gamma,     sizeof(double), 1, fp) < 1) goto L1;
     if (fread(&h->r1,        sizeof(double), 1, fp) < 1) goto L1;
     if (fread(&h->r2,        sizeof(double), 1, fp) < 1) goto L1;
     if (fread(&h->n_int1,    sizeof(int),    1, fp) < 1) goto L1;
     if (fread(&h->n_int2,    sizeof(int),    1, fp) < 1) goto L1;
     if (fread(&h->n_quad,    sizeof(int),    1, fp) < 1) goto L1;
     if (fread(&h->accuracy,  sizeof(double), 1, fp) < 1) goto L1;
     if (fread(&h->norm,      sizeof(double), 1, fp) < 1) goto L1;
     if (fread(&h->reff,      sizeof(double), 1, fp) < 1) goto L1;
     if (fread(&h->veff,      sizeof(double), 1, fp) < 1) goto L1;
     if (fread(&h->gavg,      sizeof(double), 1, fp) < 1) goto L1;
     if (fread(&h->vavg,      sizeof(double), 1, fp) < 1) goto L1;
     if (fread(&h->ravg,      sizeof(double), 1, fp) < 1) goto L1;
     if (fread(&h->rvw,       sizeof(double), 1, fp) < 1) goto L1;
     if (fread(&h->cext,      sizeof(double), 1, fp) < 1) goto L1;
     if (fread(&h->csca,      sizeof(double), 1, fp) < 1) goto L1;
     if (fread(&h->cbak,      sizeof(double), 1, fp) < 1) goto L1;
     if (fread(&h->g,         sizeof(double), 1, fp) < 1) goto L1;

     return 0;

L1:  eprintf("ERROR: problem reading header from file: %s\n", strerror(errno));
     return -1;
}



int scat_header_skip_bin(FILE *fp) {

     scat_header_data h;

     if (scat_header_read_bin_fp(fp, &h)) {
          eprintf("ERROR: scat_header_read_bin_fp()\n");
          return -1;
     }

     return 0;
}



int scat_header_write_bin_fn(const char *filename, scat_header_data *h) {

     FILE *fp;

     if ((fp = fopen(filename, "r")) == NULL) {
          eprintf("ERROR: Error opening file for writing: %s ... %s\n",
                  filename, strerror(errno));
          return -1;
     }

     if (scat_header_write_bin_fp(fp, h)) {
          eprintf("ERROR: scat_header_write_bin_fp()\n");
          return -1;
     }

     fclose(fp);

     return 0;
}



int scat_header_write_bin_fp(FILE *fp, scat_header_data *h) {

     if (fwrite(&h->lambda,    sizeof(double), 1, fp) < 1) goto L1;
     if (fwrite(&h->mr,        sizeof(double), 1, fp) < 1) goto L1;
     if (fwrite(&h->mi,        sizeof(double), 1, fp) < 1) goto L1;
     if (fwrite(&h->dist_type, sizeof(int),    1, fp) < 1) goto L1;
     if (fwrite(&h->a1,        sizeof(double), 1, fp) < 1) goto L1;
     if (fwrite(&h->b1,        sizeof(double), 1, fp) < 1) goto L1;
     if (fwrite(&h->a2,        sizeof(double), 1, fp) < 1) goto L1;
     if (fwrite(&h->b2,        sizeof(double), 1, fp) < 1) goto L1;
     if (fwrite(&h->gamma,     sizeof(double), 1, fp) < 1) goto L1;
     if (fwrite(&h->r1,        sizeof(double), 1, fp) < 1) goto L1;
     if (fwrite(&h->r2,        sizeof(double), 1, fp) < 1) goto L1;
     if (fwrite(&h->n_int1,    sizeof(int),    1, fp) < 1) goto L1;
     if (fwrite(&h->n_int2,    sizeof(int),    1, fp) < 1) goto L1;
     if (fwrite(&h->n_quad,    sizeof(int),    1, fp) < 1) goto L1;
     if (fwrite(&h->accuracy,  sizeof(double), 1, fp) < 1) goto L1;
     if (fwrite(&h->norm,      sizeof(double), 1, fp) < 1) goto L1;
     if (fwrite(&h->reff,      sizeof(double), 1, fp) < 1) goto L1;
     if (fwrite(&h->veff,      sizeof(double), 1, fp) < 1) goto L1;
     if (fwrite(&h->gavg,      sizeof(double), 1, fp) < 1) goto L1;
     if (fwrite(&h->vavg,      sizeof(double), 1, fp) < 1) goto L1;
     if (fwrite(&h->ravg,      sizeof(double), 1, fp) < 1) goto L1;
     if (fwrite(&h->rvw,       sizeof(double), 1, fp) < 1) goto L1;
     if (fwrite(&h->cext,      sizeof(double), 1, fp) < 1) goto L1;
     if (fwrite(&h->csca,      sizeof(double), 1, fp) < 1) goto L1;
     if (fwrite(&h->cbak,      sizeof(double), 1, fp) < 1) goto L1;
     if (fwrite(&h->g,         sizeof(double), 1, fp) < 1) goto L1;

     return 0;

L1:  eprintf("ERROR: problem writing header to file: %s\n", strerror(errno));
     return -1;
}



int scat_header_blank_bin(FILE *fp) {

     scat_header_data h;

     scat_header_zero(&h);

     if (scat_header_write_bin_fp(fp, &h)) {
          eprintf("ERROR: scat_header_write_bin_fp()\n");
          return -1;
     }

     return 0;
}


/*******************************************************************************
 *
 ******************************************************************************/
int scat_header_interp(double a, scat_header_data *h1, scat_header_data *h2,
                       scat_header_data *h) {

     double b = 1. - a;

     if (h1->dist_type != h1->dist_type) {
          eprintf("ERROR: cannot interpolate different size distribution types\n");
          return -1;
     }

     h->lambda    = b * h1->lambda   + a * h2->lambda;
     h->mr        = b * h1->mr       + a * h2->mr;
     h->mi        = b * h1->mi       + a * h2->mi;
     h->dist_type = h1->dist_type;
     h->a1        = b * h1->a1       + a * h2->a1;
     h->b1        = b * h1->b1       + a * h2->b1;
     h->a2        = b * h1->a2       + a * h2->a2;
     h->b2        = b * h1->b2       + a * h2->b2;
     h->gamma     = b * h1->gamma    + a * h2->gamma;
     h->r1        = b * h1->r1       + a * h2->r1;
     h->r2        = b * h1->r2       + a * h2->r2;
     h->n_int1    = MAX(h1->n_int1, h2->n_int1);
     h->n_int2    = MAX(h1->n_int2, h2->n_int2);
     h->n_quad    = MAX(h1->n_quad, h2->n_quad);
     h->accuracy  = b * h1->accuracy + a * h2->accuracy;
     h->norm      = b * h1->norm     + a * h2->norm;
     h->reff      = b * h1->reff     + a * h2->reff;
     h->veff      = b * h1->veff     + a * h2->veff;
     h->gavg      = b * h1->gavg     + a * h2->gavg;
     h->vavg      = b * h1->vavg     + a * h2->vavg;
     h->ravg      = b * h1->ravg     + a * h2->ravg;
     h->rvw       = b * h1->rvw      + a * h2->rvw;
     h->cext      = b * h1->cext     + a * h2->cext;
     h->csca      = b * h1->csca     + a * h2->csca;
     h->cbak      = b * h1->cbak     + a * h2->cbak;
     h->g         = b * h1->g        + a * h2->g;

     return 0;
}



/*******************************************************************************
 *
 ******************************************************************************/
int phase_func_read_fn(const char *filename, scat_header_data *h,
                       double **theta, double ***phase, int allsix) {

     int n_angles;

     FILE *fp;

     if ((fp = fopen(filename, "r")) == NULL) {
          eprintf("ERROR: Error opening file for reading: %s ... %s\n",
                  filename, strerror(errno));
          return -1;
     }

     if ((n_angles = phase_func_read_fp(fp, h, theta, phase, allsix)) < 0) {
          eprintf("ERROR: phase_func_read_fp()\n");
          return -1;
     }

     fclose(fp);

     return n_angles;
}



int phase_func_read_fp(FILE *fp, scat_header_data *h,
                       double **theta_, double ***phase_, int allsix) {

     int i;

     int n_angles;

     double *theta;

     double **phase;

     if (! h) {
          if (scat_header_skip(fp)) {
               eprintf("ERROR: scat_header_skip()\n");
               return -1;
          }
     }
     else {
          if (scat_header_read_fp(fp, h)) {
               eprintf("ERROR: scat_header_read_fp()\n");
               return -1;
          }
     }

     if (fscanf(fp, "%d", &n_angles) < 1) goto L1;

     theta = alloc_array1_d(n_angles);
     phase = alloc_array2_d(6, n_angles);

     for (i = 0; i < n_angles; ++i) {
          if (fscanf(fp, "%lf", &theta[i]) < 1) goto L1;

          if (! allsix) {
               if (fscanf(fp, "%lf", &phase[0][i]) < 1) goto L1;
               phase[1][i] = 0.;
               phase[2][i] = 0.;
               phase[3][i] = 0.;
               phase[4][i] = 0.;
               phase[5][i] = 0.;
          }
          else {
               if (fscanf(fp, "%lf%lf%lf%lf%lf%lf",
                          &phase[0][i], &phase[1][i], &phase[2][i],
                          &phase[3][i], &phase[4][i], &phase[5][i]) < 6) goto L1;
          }
     }

     *theta_ = theta;

     *phase_ = phase;

     return n_angles;

L1:  eprintf("ERROR: problem reading phase function from file: %s\n",
             strerror(errno));
     return -1;
}



int phase_func_write_fn(const char *filename, scat_header_data *h,
                        double *theta, double **phase, int n_angles, int allsix) {

     FILE *fp;

     if ((fp = fopen(filename, "w")) == NULL) {
          eprintf("ERROR: Error opening file for writing: %s ... %s\n",
                  filename, strerror(errno));
          return -1;
     }

     if (phase_func_write_fp(fp, h, theta, phase, n_angles, allsix) < 0) {
          eprintf("ERROR: phase_func_write_fp()\n");
          return -1;
     }

     fclose(fp);

     return 0;
}



int phase_func_write_fp(FILE *fp, scat_header_data *h,
                        double *theta, double **phase, int n_angles, int allsix) {

     int i;

     if (! h) {
          if (scat_header_blank(fp)) {
               eprintf("ERROR: scat_header_blank()\n");
               return -1;
          }
     }
     else {
          if (scat_header_write_fp(fp, h)) {
               eprintf("ERROR: scat_header_write_fp()\n");
               return -1;
          }
     }

     if (fprintf(fp, "%d\n", n_angles) < 1) goto L1;

     for (i = 0; i < n_angles; ++i) {
          fprintf(fp, "%e ", theta[i]);
          if (! allsix) {
               if (fprintf(fp, "%e\n", phase[0][i]) < 1) goto L1;
          }
          else {
               if (fprintf(fp, "%e %e %e %e %e %e\n",
                           phase[0][i], phase[1][i], phase[2][i],
                           phase[3][i], phase[4][i], phase[5][i]) < 6) goto L1;
          }
     }

     return 0;

L1:  eprintf("ERROR: problem writing phase function to file: %s\n",
             strerror(errno));
     return -1;
}



/*******************************************************************************
 *
 ******************************************************************************/
int scat_coefs_read_fn(const char *filename, int n_coef,
                       scat_header_data *h, double ***chi, int allsix) {

     int n_total;

     FILE *fp;

     if ((fp = fopen(filename, "r")) == NULL) {
          eprintf("ERROR: Error opening file for reading: %s ... %s\n",
                  filename, strerror(errno));
          return -1;
     }

     n_total = scat_coefs_read_fp(fp, n_coef, h, chi, allsix);

     fclose(fp);

     return n_total;
}



int scat_coefs_read_fp(FILE *fp, int n_coef,
                       scat_header_data *h, double ***chi_, int allsix) {

     int i;
     int j;

     int n_elem;
     int n_elem2;

     int n_coef2;
     int n_total;

     double **chi;

     n_elem = allsix ? 6 : 1;

     if (! h) {
          if (scat_header_skip(fp)) {
               eprintf("ERROR: scat_header_skip()\n");
               return -1;
          }
     }
     else {
          if (scat_header_read_fp(fp, h)) {
               eprintf("ERROR: scat_header_read_fp()\n");
               return -1;
          }
     }

     n_elem2 = 6;
/*
     if (fscanf(fp, "%d", &n_elem2) < 1) goto L1;
*/
     if (n_elem2 < n_elem)  {
          eprintf("ERROR: n_elem in file < n_elem requested\n");
          return -1;
     }

     if (fscanf(fp, "%d", &n_coef2) < 1) goto L1;

     if (n_coef <= 0)
          n_total = n_coef2;
     else
          n_total = MIN(n_coef2, n_coef);

     chi = alloc_array2_d(n_elem, n_total);
     init_array2_d(chi, n_elem, n_total, 0.);

     for (i = 0; i < n_total; ++i) {
          for (j = 0; j < n_elem; ++j)
               if (fscanf(fp, "%lf", &chi[j][i]) < 1) goto L1;

          for (     ; j < n_elem2; ++j)
               if (fscanf(fp, "%*s"            ) < 0) goto L1;
     }

     *chi_ = chi;

     return n_total;

L1:  eprintf("ERROR: problem reading phase coefficients from file: %s\n",
             strerror(errno));
     return -1;
}



int scat_coefs_write_fn(const char *filename, int n_coef2,
                        scat_header_data *h, double **chi, int allsix) {

     int ret_val;

     FILE *fp;

     if ((fp = fopen(filename, "w")) == NULL) {
          eprintf("ERROR: Error opening file for writing: %s ... %s\n",
                  filename, strerror(errno));
          return -1;
     }

     ret_val = scat_coefs_write_fp(fp, n_coef2, h, chi, allsix);

     fclose(fp);

     return ret_val;
}



int scat_coefs_write_fp(FILE *fp, int n_coef2,
                        scat_header_data *h, double **chi, int allsix) {

     int i;

     if (! h) {
          if (scat_header_blank(fp)) {
               eprintf("ERROR: scat_header_blank()\n");
               return -1;
          }
     }
     else {
          if (scat_header_write_fp(fp, h)) {
               eprintf("ERROR: scat_header_write_fp()\n");
               return -1;
          }
     }

     if (fprintf(fp, "%d\n", n_coef2) < 1) goto L1;

     for (i = 0; i < n_coef2; ++i) {
          if (! allsix) {
               if (fprintf(fp, "%e\n", chi[0][i]) < 1) goto L1;
          }
          else {
               if (fprintf(fp, "%e %e %e %e %e %e\n",
                           chi[0][i], chi[1][i], chi[2][i],
                           chi[3][i], chi[4][i], chi[5][i]) < 6) goto L1;
          }
     }

     return 0;


L1:  eprintf("ERROR: problem writing phase coefficients to file: %s\n",
             strerror(errno));
     return -1;
}



/*******************************************************************************
 *
 ******************************************************************************/
int scat_coefs_read_bin_fn(const char *filename, int n_coef,
                           scat_header_data *h, double ***chi, int allsix) {

     int n_total;

     FILE *fp;

     if ((fp = fopen(filename, "r")) == NULL) {
          eprintf("ERROR: Error opening file for reading: %s ... %s\n",
                  filename, strerror(errno));
          return -1;
     }

     n_total = scat_coefs_read_bin_fp(fp, n_coef, h, chi, allsix);

     fclose(fp);

     return n_total;
}



int scat_coefs_read_bin_fp(FILE *fp, int n_coef,
                           scat_header_data *h, double ***chi_, int allsix) {

     int n_coef2;
     int n_total;

     double **chi;

     if (! h) {
          if (scat_header_skip_bin(fp)) {
               eprintf("ERROR: scat_header_skip_bin()\n");
               return -1;
          }
     }
     else {
          if (scat_header_read_bin_fp(fp, h)) {
               eprintf("ERROR: scat_header_read_bin_fp()\n");
               return -1;
          }
     }

     if (fread(&n_coef2, sizeof(int), 1, fp) < 1) goto L1;

     if (n_coef <= 0)
          n_total = n_coef2;
     else
          n_total = MIN(n_coef2, n_coef);

     chi = alloc_array2_d(6, n_total);
     init_array2_d(chi, 6, n_total, 0.);

     if (! allsix) {
          if (fread(chi[0], sizeof(double), n_total, fp) < (size_t) n_total) goto L1;

          init_array1_d(chi[1], n_total, 0.);
          init_array1_d(chi[2], n_total, 0.);
          init_array1_d(chi[3], n_total, 0.);
          init_array1_d(chi[4], n_total, 0.);
          init_array1_d(chi[5], n_total, 0.);
     }
     else {
          if (fread(chi[0], sizeof(double), n_total, fp) < (size_t) n_total) goto L1;
          if (fread(chi[1], sizeof(double), n_total, fp) < (size_t) n_total) goto L1;
          if (fread(chi[2], sizeof(double), n_total, fp) < (size_t) n_total) goto L1;
          if (fread(chi[3], sizeof(double), n_total, fp) < (size_t) n_total) goto L1;
          if (fread(chi[4], sizeof(double), n_total, fp) < (size_t) n_total) goto L1;
          if (fread(chi[5], sizeof(double), n_total, fp) < (size_t) n_total) goto L1;
     }

     *chi_ = chi;

     return n_total;

L1:  eprintf("ERROR: problem reading phase coefficients from file: %s\n",
             strerror(errno));
     return -1;
}



int scat_coefs_write_bin_fn(const char *filename, int n_coef2,
                            scat_header_data *h, double **chi, int allsix) {

     int ret_val;

     FILE *fp;

     if ((fp = fopen(filename, "w")) == NULL) {
          eprintf("ERROR: Error opening file for writing: %s ... %s\n",
                  filename, strerror(errno));
          return -1;
     }

     ret_val = scat_coefs_write_bin_fp(fp, n_coef2, h, chi, allsix);

     fclose(fp);

     return ret_val;
}



int scat_coefs_write_bin_fp(FILE *fp, int n_coef2,
                            scat_header_data *h, double **chi, int allsix) {

     if (! h) {
          if (scat_header_blank_bin(fp)) {
               eprintf("ERROR: scat_header_blank_bin()\n");
               return -1;
          }
     }
     else {
          if (scat_header_write_bin_fp(fp, h)) {
               eprintf("ERROR: scat_header_write_bin_fp()\n");
               return -1;
          }
     }

     if (fwrite(&n_coef2, sizeof(int), 1, fp) < 1) goto L1;

     if (! allsix) {
          if (fwrite(chi[0], sizeof(double), n_coef2, fp) < (size_t) n_coef2) goto L1;
     }
     else {
          if (fwrite(chi[0], sizeof(double), n_coef2, fp) < (size_t) n_coef2) goto L1;
          if (fwrite(chi[1], sizeof(double), n_coef2, fp) < (size_t) n_coef2) goto L1;
          if (fwrite(chi[2], sizeof(double), n_coef2, fp) < (size_t) n_coef2) goto L1;
          if (fwrite(chi[3], sizeof(double), n_coef2, fp) < (size_t) n_coef2) goto L1;
          if (fwrite(chi[4], sizeof(double), n_coef2, fp) < (size_t) n_coef2) goto L1;
          if (fwrite(chi[5], sizeof(double), n_coef2, fp) < (size_t) n_coef2) goto L1;
     }

     return 0;


L1:  eprintf("ERROR: problem writing phase coefficients to file: %s\n",
             strerror(errno));
     return -1;
}



/*******************************************************************************
 *
 ******************************************************************************/
int scat_coefs_interp(double a,
                      int n_coef1, scat_header_data *h1, double **chi1,
                      int n_coef2, scat_header_data *h2, double **chi2,
                      scat_header_data *h , double **chi, int allsix) {

     int i;
     int j;

     int n_elem;
     int n_coef;

     double b = 1. - a;

     if (scat_header_interp(a, h1, h2, h)) {
          eprintf("ERROR: scat_header_interp()\n");
          return -1;
     }

     n_coef = MAX(n_coef1, n_coef2);

     n_elem = allsix ? 6 : 1;

     for (i = 0; i < n_elem; ++i) {
          for (j = 0; j < n_coef; ++j) {
               chi[i][j] = 0.;

               if (j < n_coef1)
                    chi[i][j] += b * chi1[i][j];
               if (j < n_coef2 && a != 0.)
                    chi[i][j] += a * chi2[i][j];
          }
     }

     return n_coef;
}



/*******************************************************************************
 *
 ******************************************************************************/
void scat_coefs_hg(int n_coef, double g, double *chi) {

     int i;

     double a;

     a = 1.;
     for (i = 0; i < n_coef; ++i) {
          chi[i] = (2. * i + 1.) * a;
          a *= g;
     }
}



void scat_coefs_hg_l(int n_coef, double g, double g_l, double *chi_l) {

     int i;

     double a;
     double a_l;

     a   = 1.;
     a_l = 0.;
     for (i = 0; i < n_coef; ++i) {
          chi_l[i] = (2. * i + 1.) * a_l;
          a_l = a_l * g + a * g_l;
          a   = a * g;
     }
}



void scat_coefs_dhg(int n_coef, double g1, double g2, double f, double *chi) {

     int i;

     double a1;
     double a2;

     a1 = 1.;
     a2 = 1.;
     for (i = 0; i < n_coef; ++i) {
          chi[i] = (2. * i + 1.) * (f * a1 + (1. - f) * a2);
          a1 *= g1;
          a2 *= g2;
     }
}




void scat_coefs_dhg_l(int n_coef, double g1, double g2, double f,
                      double g1_l, double g2_l, double f_l, double *chi_l) {

     int i;

     double a1;
     double a1_l;
     double a2;
     double a2_l;

     a1   = 1.;
     a2   = 1.;
     a1_l = 0.;
     a2_l = 0.;
     for (i = 0; i < n_coef; ++i) {
          chi_l[i] = (2. * i + 1.) * (f_l * a1 + f * a1_l +
                                     -f_l * a2 + (1. - f) * a2_l);
          a1_l = a1_l * g1 + a1 * g1_l;
          a1   = a1 * g1;
          a2_l = a2_l * g2 + a2 * g2_l;
          a2   = a2 * g2;
     }
}



void scat_coefs_ray(double rho, double **chi) {

     double a;

     a = (1. - rho) / (2. + rho);

     chi[0][0] = 1.;
     chi[1][0] = 0.;
     chi[2][0] = 0.;
     chi[3][0] = 0.;
     chi[4][0] = 0.;
     chi[5][0] = 0.;

     chi[0][1] = 0.;
     chi[1][1] = 0.;
     chi[2][1] = 0.;
     chi[3][1] = 3. * (1. - 2. * rho) / (2. + rho);
     chi[4][1] = 0.;
     chi[5][1] = 0.;

     chi[0][2] = a;
     chi[1][2] =       6. * a;
     chi[2][2] = 0.;
     chi[3][2] = 0.;
     chi[4][2] = sqrt(6.) * a;
     chi[5][2] = 0.;
}



void scat_coefs_ray_l(double rho, double rho_l, double **chi) {

     double a_l;

     a_l = -rho_l / (2. + rho) * (1. + (1. - rho) / (2. + rho));

     chi[0][0] = 0.;
     chi[1][0] = 0.;
     chi[2][0] = 0.;
     chi[3][0] = 0.;
     chi[4][0] = 0.;
     chi[5][0] = 0.;

     chi[0][1] = 0.;
     chi[1][1] = 0.;
     chi[2][1] = 0.;
     chi[3][1] = 3. * -rho_l / (2. + rho) * (2. + (1. - 2. * rho) / (2. + rho));
     chi[4][1] = 0.;
     chi[5][1] = 0.;

     chi[0][2] = a_l;
     chi[1][2] =       6. * a_l;
     chi[2][2] = 0.;
     chi[3][2] = 0.;
     chi[4][2] = sqrt(6.) * a_l;
     chi[5][2] = 0.;
}



/*******************************************************************************
 * 
 ******************************************************************************/
int scat_total_read_fn(const char *filename, scat_header_data *h, int *n_quad,
                       int *n_four, int *quad_type, double *cextx, double *cexty,
                       double *cscatx, double *cscaty, double *cbackx, double *cbacky,
                       double *gx, double *gy) {

     FILE *fp;

     if ((fp = fopen(filename, "r")) == NULL) {
          eprintf("ERROR: Error opening file for reading: %s ... %s\n",
                  filename, strerror(errno));
          return -1;
     }

     if (scat_total_read_fp(fp, h, n_quad, n_four, quad_type, cextx, cexty,
                            cscatx, cscaty, cbackx, cbacky, gx, gy)) {
          eprintf("ERROR: scat_total_read_fp(), %s\n", filename);
          return -1;
     }

     fclose(fp);

     return 0;
}



int scat_total_read_fp(FILE *fp, scat_header_data *h, int *n_quad,
                       int *n_four, int *quad_type, double *cextx, double *cexty,
                       double *cscatx, double *cscaty, double *cbackx, double *cbacky,
                       double *gx, double *gy) {

     int i;

     char *quad_name;

     int line;

     int i_four;

     line = 1;

     if (! h) {
          if (scat_header_skip(fp)) {
               eprintf("ERROR: scat_header_skip()\n");
               return -1;
          }
     }
     else {
          if (scat_header_read_fp(fp, h)) {
               eprintf("ERROR: scat_header_read_fp()\n");
               return -1;
          }
     }

     quad_name = (char *) malloc(size_dist_name_max_length());

     if (fscanf(fp, "%d%d%s%*[\n]", n_quad, &i_four, quad_name) < 3) {
          eprintf("ERROR: problem reading from file, line %d\n", line);
          return -1;
     }

     *n_four = i_four + 1;

     if ((*quad_type = quadrature_code(quad_name)) < 0) {
          eprintf("ERROR: quadrature_code(), line %d\n", line);
          return -1;
     }

     free(quad_name);

     line++;

     for (i = 0; i < 2 * *n_quad; ++i) {
           if (fscanf(fp, "%lf%lf%lf%lf%lf%lf%lf%lf",
                      &cextx[i], &cexty[i], &cscatx[i], &cscaty[i],
                      &cbackx[i], &cbacky[i], &gx[i], &gy[i]) < 6) {
                eprintf("ERROR: problem reading from file, line %d\n", line);
                return -1;
           }

           line++;
     }

     return 0;
}



int scat_total_write_fn(const char *filename, scat_header_data *h, int n_quad,
                        int n_four, int quad_type, double *cextx, double *cexty,
                        double *cscatx, double *cscaty, double *cbackx, double *cbacky,
                        double *gx, double *gy) {

     FILE *fp;

     if ((fp = fopen(filename, "w")) == NULL) {
          eprintf("ERROR: Error opening file for writing: %s ... %s\n",
                  filename, strerror(errno));
          return -1;
     }

     if (scat_total_write_fp(fp, h, n_quad, n_four, quad_type, cextx, cexty,
                             cscatx, cscaty, cbackx, cbacky, gx, gy)) {
          eprintf("ERROR: scat_total_write_fp(), %s\n", filename);
          return -1;
     }

     fclose(fp);

     return 0;
}



int scat_total_write_fp(FILE *fp, scat_header_data *h, int n_quad,
                        int n_four, int quad_type, double *cextx, double *cexty,
                        double *cscatx, double *cscaty, double *cbackx, double *cbacky,
                        double *gx, double *gy) {

     int i;

     if (! h) {
          if (scat_header_blank(fp)) {
               eprintf("ERROR: scat_header_blank()\n");
               return -1;
          }
     }
     else {
          if (scat_header_write_fp(fp, h)) {
               eprintf("ERROR: scat_header_write_fp()\n");
               return -1;
          }
     }

     if (fprintf(fp, "%d %d %s\n", n_quad, n_four - 1,
                 quadrature_name((enum quadrature_type) quad_type)) < 0) {
          eprintf("ERROR: problem writing to file\n");
          return -1;
     }

     for (i = 0; i < 2 * n_quad; ++i) {
          if (fprintf(fp, "%e %e %e %e %e %e %e %e\n",
                      cextx[i], cexty[i], cscatx[i], cscaty[i],
                      cbackx[i], cbacky[i], gx[i], gy[i]) < 0) {
               eprintf("ERROR: problem writing to file\n");
               return -1;
          }
     }

     return 0;
}



/*******************************************************************************
 *
 ******************************************************************************/
int phase_mats_read_fn(const char *filename, scat_header_data *h, int *n_quad,
                       int *n_four, int *quad_type, double ******P,
                       double ****Kext, double ***Kemi) {

     FILE *fp;

     if ((fp = fopen(filename, "r")) == NULL) {
          eprintf("ERROR: Error opening file for reading: %s, %s\n",
                  filename, strerror(errno));
          return -1;
     }

     if (phase_mats_read_fp(fp, h, n_quad, n_four, quad_type, P, Kext, Kemi)) {
          eprintf("ERROR: phase_mats_read_fp(), %s\n", filename);
          return -1;
     }

     fclose(fp);

     return 0;
}



int phase_mats_read_fp(FILE *fp, scat_header_data *h, int *n_quad,
                       int *n_four, int *quad_type, double ******P_,
                       double ****Kext_, double ***Kemi_) {

     char *quad_name;

     int i;
     int j;
     int k;
     int l;
/*
     int c;
*/
     int line;

     int i_four;

     double mu1;
     double mu2;

     double *****P;
     double ***Kext;
     double **Kemi;

     line = 1;
/*
     while (1) {
          if ((c = fgetc(fp)) == EOF) {
               eprintf("ERROR: problem reading from file, line %d\n", line);
               return -1;
          }

          if (c != 'C')
               break;

          read_newline(fp);

          line++;
     }

     ungetc(c, fp);
*/
     if (! h) {
          if (scat_header_skip(fp)) {
               eprintf("ERROR: scat_header_skip()\n");
               return -1;
          }
     }
     else {
          if (scat_header_read_fp(fp, h)) {
               eprintf("ERROR: scat_header_read_fp()\n");
               return -1;
          }
     }

     quad_name = (char *) malloc(size_dist_name_max_length());

     if (fscanf(fp, "%d%d%s%*[\n]", n_quad, &i_four, quad_name) < 3) {
          eprintf("ERROR: problem reading from file, line %d\n", line);
          return -1;
     }

     *n_four = i_four + 1;

     if ((*quad_type = quadrature_code(quad_name)) < 0) {
          eprintf("ERROR: quadrature_code(), line %d\n", line);
          return -1;
     }

     free(quad_name);

     line++;

     P    = alloc_array5_d(*n_four, 2 * *n_quad, 2 * *n_quad, 4, 4);
     Kext = alloc_array3_d(         2 * *n_quad, 4, 4);
     Kemi = alloc_array2_d(         2 * *n_quad, 4);

     read_newline(fp);
     line++;

     for (k = 0; k < 2 * *n_quad; ++k) {
          for (j = 0; j < 2 * *n_quad; ++j) {
               for (i = 0; i < *n_four; ++i) {
                    if (fscanf(fp, "%lf%lf%d", &mu1, &mu2, &i_four) < 3) {
                         eprintf("ERROR: problem reading from file, line %d\n", line);
                         return -1;
                    }

                    line++;

                    if (i_four != i) {
                         eprintf("ERROR: invalid fourier order encountered\n");
                         return -1;
                    }

                    for (l = 0; l < 4; ++l) {
                         if (fscanf(fp, "%lf%lf%lf%lf\n",
                                    &P[i][j][k][l][0], &P[i][j][k][l][1],
                                    &P[i][j][k][l][2], &P[i][j][k][l][3]) < 4) {
                              eprintf("ERROR: problem reading from file, line %d\n", line);
                              return -1;
                         }

                         line++;
                    }
               }
          }
     }

     read_newline(fp);

     for (i = 0; i < 2 * *n_quad; ++i) {
          if (fscanf(fp, "%lf", &mu1) < 1) {
               eprintf("ERROR: problem reading from file, line %d\n", line);
               return -1;
          }

          line++;

          for (j = 0; j < 4; ++j) {
               if (fscanf(fp, "%lf%lf%lf%lf\n",
                          &Kext[i][j][0], &Kext[i][j][1],
                          &Kext[i][j][2], &Kext[i][j][3]) < 4) {
                    eprintf("ERROR: problem reading from file, line %d\n", line);
                    return -1;
               }

               line++;
          }
     }

     read_newline(fp);

     for (i = 0; i < 2 * *n_quad; ++i) {
          if (fscanf(fp, "%lf%lf%lf%lf%lf\n", &mu1,
                     &Kemi[i][0], &Kemi[i][1],
                     &Kemi[i][2], &Kemi[i][3]) < 5) {
               eprintf("ERROR: problem reading from file, line %d\n", line);
               return -1;
          }

          line++;
     }

     *P_    = P;
     *Kext_ = Kext;
     *Kemi_ = Kemi;

     return 0;
}



int phase_mats_write_fn(const char *filename, scat_header_data *h, int n_quad,
                        int n_four, int quad_type, double *qx, double *****P,
                        double ***Kext, double **Kemi) {

     FILE *fp;

     if ((fp = fopen(filename, "w")) == NULL) {
          eprintf("ERROR: Error opening file for writing: %s, %s\n",
                  filename, strerror(errno));
          return -1;
     }

     if (phase_mats_write_fp(fp, h, n_quad, n_four, quad_type, qx, P, Kext, Kemi)) {
          eprintf("ERROR: phase_mats_write_fp(), %s\n", filename);
          return -1;
     }

     fclose(fp);

     return 0;
}



int phase_mats_write_fp(FILE *fp, scat_header_data *h, int n_quad,
                        int n_four, int quad_type, double *qx, double *****P,
                        double ***Kext, double **Kemi) {

     int i;
     int j;
     int k;
     int l;

     double mu1 = 0.;
     double mu2 = 0.;

     if (! h) {
          if (scat_header_blank(fp)) {
               eprintf("ERROR: scat_header_blank()\n");
               return -1;
          }
     }
     else {
          if (scat_header_write_fp(fp, h)) {
               eprintf("ERROR: scat_header_write_fp()\n");
               return -1;
          }
     }

     if (fprintf(fp, "%d %d %s\n", n_quad, n_four - 1,
                 quadrature_name((enum quadrature_type) quad_type)) < 0) {
          eprintf("ERROR: problem writing to file\n");
          return -1;
     }

     if (fprintf(fp, "C SCATTERING MATRIX\n") < 0) {
          eprintf("ERROR: problem writing to file\n");
          return -1;
     }

     for (k = 0; k < 2 * n_quad; ++k) {
          for (j = 0; j < 2 * n_quad; ++j) {
               for (i = 0; i < n_four; ++i) {
                    if (qx) {
                         mu1 = k < n_quad ? qx[k] : -qx[k - n_quad];
                         mu2 = j < n_quad ? qx[j] : -qx[j - n_quad];
                    }
                    if (fprintf(fp, "%16e%16e %d\n", mu1, mu2, i) < 0) {
                         eprintf("ERROR: problem writing to file\n");
                         return -1;
                    }

                    for (l = 0; l < 4; ++l) {
                         if (fprintf(fp, "%18e%16e%16e%16e\n",
                                     P[i][j][k][l][0], P[i][j][k][l][1],
                                     P[i][j][k][l][2], P[i][j][k][l][3]) < 0) {
                              eprintf("ERROR: problem writing to file\n");
                              return -1;
                         }
                    }
               }
          }
     }

     if (fprintf(fp, "C EXTINCTION MATRIX\n") < 0) {
          eprintf("ERROR: problem writing to file\n");
          return -1;
     }

     for (i = 0; i < 2 * n_quad; ++i) {
          if (qx)
               mu1 = i < n_quad ? qx[i] : -qx[i - n_quad];

          if (fprintf(fp, "%16e\n", mu1) < 0) {
               eprintf("ERROR: problem writing to file\n");
               return -1;
          }

          for (j = 0; j < 4; ++j) {
               if (fprintf(fp, "%18e%16e%16e%16e\n",
                           Kext[i][j][0], Kext[i][j][1],
                           Kext[i][j][2], Kext[i][j][3]) < 0) {
                    eprintf("ERROR: problem writing to file\n");
                    return -1;
               }
          }
     }

     if (fprintf(fp, "C EMISION VECTOR\n") < 0) {
          eprintf("ERROR: problem writing to file\n");
          return -1;
     }

     for (i = 0; i < 2 * n_quad; ++i) {
          if (qx)
               mu1 = i < n_quad ? qx[i] : -qx[i - n_quad];

          if (fprintf(fp, "%16e%16e%16e%16e%16e\n", mu1,
                      Kemi[i][0], Kemi[i][1],
                      Kemi[i][2], Kemi[i][3]) < 0) {
               eprintf("ERROR: problem writing to file\n");
               return -1;
          }
     }

     return 0;
}



/*******************************************************************************
 *
 ******************************************************************************/
int is_user_defined_n_coef(char *token) {

     if (strcmp(token, "hg")    == 0 || strcmp(token, "hg_l")  == 0 ||
         strcmp(token, "dhg")   == 0 || strcmp(token, "dhg_l") == 0 ||
         strcmp(token, "none")  == 0) {
          return 1;
     }

     return 0;
}



int load_scat_coefs(const char *filename, int n_coef, double ***chi_, int *flag) {

     char *temp;

     char *token;
     char *token2;
     char *lasts = NULL;

     int i;
     int j;

     int n_coef2;

     double g;
     double g1;
     double g2;
     double f;
     double rho;

     double g_l;
     double g1_l;
     double g2_l;
     double f_l;
     double rho_l;

     double **chi;

     temp = strdup(filename);

     token = strtok_r(temp, ":", &lasts);

     if (is_user_defined_n_coef(token)) {
          if (n_coef <= 0) {
               eprintf("ERROR: must specifiy the number of phase function "
                      "coefficients for the %s phase function()\n", token);
               return -1;
          }

          chi = alloc_array2_d(6, n_coef);
          init_array2_d(chi, 6, n_coef, 0.);
     }

     if (strcmp(token, "hg") == 0) {
          token2 = strtok_r(NULL, "", &lasts);

          if (strtod_errmsg(token2, "asymmetry parameter g", &g)) {
               eprintf("ERROR: strtod_errmsg()\n");
               goto L1;
          }
/*
          if (sscanf(token2, "%lf", &g) != 1) {
               eprintf("ERROR: invalid value for asymmetry parameter g: %s\n", token2);
               goto L1;
          }
*/
          scat_coefs_hg(n_coef, g, chi[0]);

          n_coef2 = n_coef;

          *flag = SCAT_COEF_GC;
     }
     else if (strcmp(token, "hg_l") == 0) {
          token2 = strtok_r(NULL, ":", &lasts);

          if (strtod_errmsg(token2, "asymmetry parameter g", &g)) {
               eprintf("ERROR: strtod_errmsg()\n");
               goto L1;
          }
/*
          if (sscanf(token2, "%lf", &g) != 1) {
               eprintf("ERROR: invalid value for asymmetry parameter g: %s\n", token2);
               goto L1;
          }
*/
          token2 = strtok_r(NULL, "", &lasts);

          if (strtod_errmsg(token2, "linearized asymmetry parameter g_l", &g_l)) {
               eprintf("ERROR: strtod_errmsg()\n");
               goto L1;
          }
/*
          if (sscanf(token2, "%lf", &g_l) != 1) {
               eprintf("ERROR: invalid value for asymmetry parameter g_l: %s\n", token2);
               goto L1;
          }
*/
          scat_coefs_hg_l(n_coef, g, g_l, chi[0]);

          n_coef2 = n_coef;

          *flag = SCAT_COEF_GC;
     }
     else if (strcmp(token, "dhg") == 0) {
          token2 = strtok_r(NULL, ":", &lasts);

          if (sscanf(token2, "%lf", &g1) != 1) {
               eprintf("ERROR: invalid value for forward asymmetry parameter g1: %s\n", token2);
               goto L1;
          }
/*
          if (sscanf(token2, "%lf", &g1) != 1) {
               eprintf("ERROR: invalid value for forward asymmetry parameter g1: %s\n", token2);
               goto L1;
          }
*/
          token2 = strtok_r(NULL, ":", &lasts);

          if (sscanf(token2, "%lf", &g2) != 1) {
               eprintf("ERROR: invalid value for backward asymmetry parameter g2: %s\n", token2);
               goto L1;
          }
/*
          if (sscanf(token2, "%lf", &g2) != 1) {
               eprintf("ERROR: invalid value for backward asymmetry parameter g2: %s\n", token2);
               goto L1;
          }
*/
          token2 = strtok_r(NULL, "", &lasts);

          if (sscanf(token2, "%lf", &f) != 1) {
               eprintf("ERROR: invalid value for forward/backward weight f: %s\n", token2);
               goto L1;
          }
/*
          if (sscanf(token2, "%lf", &f) != 1) {
               eprintf("ERROR: invalid value for forward/backward weight f: %s\n", token2);
               goto L1;
          }
*/
          scat_coefs_dhg(n_coef, g1, g2, f, chi[0]);

          n_coef2 = n_coef;

          *flag = SCAT_COEF_GC;
     }
     else if (strcmp(token, "dhg_l") == 0) {
          token2 = strtok_r(NULL, ":", &lasts);

          if (sscanf(token2, "%lf", &g1) != 1) {
               eprintf("ERROR: invalid value for forward asymmetry parameter g1: %s\n", token2);
               goto L1;
          }
/*
          if (sscanf(token2, "%lf", &g1) != 1) {
               eprintf("ERROR: invalid value for forward asymmetry parameter g1: %s\n", token2);
               goto L1;
          }
*/
          token2 = strtok_r(NULL, ":", &lasts);

          if (sscanf(token2, "%lf", &g2) != 1) {
               eprintf("ERROR: invalid value for backward asymmetry parameter g2: %s\n", token2);
               goto L1;
          }
/*
          if (sscanf(token2, "%lf", &g2) != 1) {
               eprintf("ERROR: invalid value for backward asymmetry parameter g2: %s\n", token2);
               goto L1;
          }
*/
          token2 = strtok_r(NULL, ":", &lasts);

          if (sscanf(token2, "%lf", &f) != 1) {
               eprintf("ERROR: invalid value for forward/backward weight f: %s\n", token2);
               goto L1;
          }
/*
          if (sscanf(token2, "%lf", &f) != 1) {
               eprintf("ERROR: invalid value for forward/backward weight f: %s\n", token2);
               goto L1;
          }
*/
          token2 = strtok_r(NULL, ":", &lasts);

          if (sscanf(token2, "%lf", &g1_l) != 1) {
               eprintf("ERROR: invalid value for linearized forward asymmetry parameter g1_l: %s\n", token2);
               goto L1;
          }
/*
          if (sscanf(token2, "%lf", &g1_l) != 1) {
               eprintf("ERROR: invalid value for linearized forward asymmetry parameter g1_l: %s\n", token2);
               goto L1;
          }
*/
          token2 = strtok_r(NULL, ":", &lasts);

          if (sscanf(token2, "%lf", &g2_l) != 1) {
               eprintf("ERROR: invalid value for linearized backward asymmetry parameter g2_l: %s\n", token2);
               goto L1;
          }
/*
          if (sscanf(token2, "%lf", &g2_l) != 1) {
               eprintf("ERROR: invalid value for linearized backward asymmetry parameter g2_l: %s\n", token2);
               goto L1;
          }
*/
          token2 = strtok_r(NULL, "", &lasts);

          if (sscanf(token2, "%lf", &f_l) != 1) {
               eprintf("ERROR: invalid value for linearized forward/backward weight f_l: %s\n", token2);
               goto L1;
          }
/*
          if (sscanf(token2, "%lf", &f_l) != 1) {
               eprintf("ERROR: invalid value for linearized forward/backward weight f_l: %s\n", token2);
               goto L1;
          }
*/
          scat_coefs_dhg_l(n_coef, g1, g2, f, g1_l, g2_l, f_l, chi[0]);

          n_coef2 = n_coef;

          *flag = SCAT_COEF_GC;
     }
     else if (strcmp(token, "gc") == 0 || strcmp(token, "gc6") == 0) {
          token2 = strtok_r(NULL, ":", &lasts);

          if ((n_coef2 = scat_coefs_read_fn(token2, n_coef, NULL, &chi, token[2]))  < 0) {
               eprintf("ERROR: scat_coefs_read_fn()\n");
               return -1;
          }

          *flag = SCAT_COEF_GC;
     }
     else if (strcmp(token, "lc") == 0 || strcmp(token, "lc6") == 0) {
          token2 = strtok_r(NULL, ":", &lasts);

          if ((n_coef2 = scat_coefs_read_fn(token2, n_coef, NULL, &chi, token[2]))  < 0) {
               eprintf("ERROR: scat_coefs_read_fn()\n");
               return -1;
          }

          *flag = SCAT_COEF_LC;
     }
     else if (strcmp(token, "none") == 0) {
          init_array2_d(chi, 6, n_coef, 0.);

          n_coef2 = n_coef;

          *flag = SCAT_COEF_NONE;
     }
     else if (strcmp(token, "ray") == 0) {
          token2 = strtok_r(NULL, ":", &lasts);

          if (sscanf(token2, "%lf", &rho) != 1) {
               eprintf("ERROR: invalid value for ray depolarization factor rho: %s\n", token2);
               goto L1;
          }
/*
          if (sscanf(token2, "%lf", &rho) != 1) {
               eprintf("ERROR: invalid value for ray depolarization factor rho: %s\n", token2);
               return -1;
          }
*/
          chi = alloc_array2_d(6, 3);
          init_array2_d(chi, 6, 3, 0.);

          scat_coefs_ray(rho, chi);

          n_coef2 = 3;

          *flag = SCAT_COEF_GC;
     }
     else if (strcmp(token, "ray_l") == 0) {
          token2 = strtok_r(NULL, ":", &lasts);

          if (sscanf(token2, "%lf", &rho) != 1) {
               eprintf("ERROR: invalid value for ray depolarization factor rho: %s\n", token2);
               goto L1;
          }
/*
          if (sscanf(token2, "%lf", &rho) != 1) {
               eprintf("ERROR: invalid value for ray depolarization factor rho: %s\n", token2);
               return -1;
          }
*/
          token2 = strtok_r(NULL, "", &lasts);

          if (sscanf(token2, "%lf", &rho_l) != 1) {
               eprintf("ERROR: invalid value for linearized ray depolarization factor rho_l: %s\n", token2);
               goto L1;
          }
/*
          if (sscanf(token2, "%lf", &rho_l) != 1) {
               eprintf("ERROR: invalid value for linearized ray depolarization factor rho_l: %s\n", token2);
               return -1;
          }
*/
          chi = alloc_array2_d(6, 3);
          init_array2_d(chi, 6, 3, 0.);

          scat_coefs_ray_l(rho, rho_l, chi);

          n_coef2 = 3;

          *flag = SCAT_COEF_GC;
     }
     else {
          eprintf("ERROR: invalid phase input type: %s\n", token);
          return -1;
     }

     if (chi[0][0] != 0. && chi[0][0] != 1.) {
          for (i = 0; i < 6; ++i) {
               for (j = 0; j < n_coef2; ++j) {
                    chi[i][j] /= chi[0][0];
               }
          }
     }

     *chi_ = chi;

     free(temp);

     return n_coef2;

L1:  free_array2_d(chi);

     return -1;
}



int load_scat_coefs2(const char *filename, int n_coef, int n_coef2, double ***chi_, int *flag) {

     char *temp;

     char *token;
     char *lasts = NULL;

     temp = strdup(filename);

     token = strtok_r(temp, ":", &lasts);

     if (! is_user_defined_n_coef(token))
          n_coef = n_coef2;

     free(temp);

     if ((n_coef2 = load_scat_coefs(filename, n_coef, chi_, flag)) < 0) {
          eprintf("ERROR: load_scat_coefs()\n");
          return -1;
     }

     return n_coef2;
}
