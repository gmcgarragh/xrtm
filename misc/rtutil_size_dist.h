/*******************************************************************************
**
**    Copyright (C) 2007-2020 Greg McGarragh <greg.mcgarragh@colostate.edu>
**
**    This source code is licensed under the GNU General Public License (GPL),
**    Version 3.  See the file COPYING for more details.
**
*******************************************************************************/

#ifndef RTUTIL_SIZE_DIST_H
#define RTUTIL_SIZE_DIST_H

#include <gutil.h>

#ifdef __cplusplus
extern "C" {
#endif


#define SIZE_DIST_MAX_PARAMS 5

enum size_dist_type {
     SIZE_DIST_NONE,
     SIZE_DIST_MODIFIED_GAMMA,
     SIZE_DIST_LOG_NORMAL,
     SIZE_DIST_POWER_LAW,
     SIZE_DIST_GAMMA,
     SIZE_DIST_MODIFIED_POWER_LAW,
     SIZE_DIST_MODIFIED_BIMODAL_LOG_NORMAL,
     SIZE_DIST_EXPONENTIAL,
     SIZE_DIST_MONO,

     N_SIZE_DIST_TYPES
};


int get_particle_dist_n(enum size_dist_type dist_type, int n_int1, int n_int2, int n_quad);

void get_power_law_range(double a1, double a2, double *r1, double *r2);

int get_particle_dist(enum size_dist_type dist_type,
                      int n_int1, int n_int2, int n_quad, int n_derivs,
                      double a1, double a2, double a3, double a4,
                      double a5, double r1, double r2,
                      double *a1_l, double *a2_l, double *a3_l, double *a4_l,
                      double *a5_l, double *r1_l, double *r2_l,
                      double *qx, double *qw, double *nr,
                      double *r1_, double *r2_,
                      double **qx_l, double **qw_l,double **nr_l,
                      double *r1_l_, double *r2_l_,
                      double *norm, double *norm_l);

int get_dist_parameters(int n_quad, int n_derivs,
                        double *qx,   double *qw,   double *nr,
                        double **qx_l, double **qw_l, double **nr_l,
                        double *reff,   double *veff,   double *gavg,
                        double *vavg,   double *ravg,   double *rvw,
                        double *reff_l, double *veff_l, double *gavg_l,
                        double *vavg_l, double *ravg_l, double *rvw_l);


#ifdef __cplusplus
}
#endif

#endif /* RTUTIL_SIZE_DIST_H */
