/*******************************************************************************
**
**    Copyright (C) 2007-2020 Greg McGarragh <greg.mcgarragh@colostate.edu>
**
**    This source code is licensed under the GNU General Public License (GPL),
**    Version 3.  See the file COPYING for more details.
**
*******************************************************************************/

#ifndef XRTM_BRDF_H
#define XRTM_BRDF_H

#include "xrtm_work.h"

#ifdef __cplusplus
extern "C" {
#endif


#define MAX_KERNEL_PARAMS 4

enum xrtm_kernel_type {
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
     XRTM_KERNEL_USER_DEFINED,

     N_XRTM_KERNEL_TYPES
};


enum xrtm_brdf_geometry {
     XRTM_BRDF_GEOMETRY_QQ,
     XRTM_BRDF_GEOMETRY_Q0,
     XRTM_BRDF_GEOMETRY_UQ,
     XRTM_BRDF_GEOMETRY_U0,

     N_XRTM_BRDF_GEOMETRIES
};


typedef struct {
     double *cos_theta;
     double *sin_theta;
     double *tan_theta;
     double *tan_theta2;
     double *cos_phi;
     double *sin_phi;
     double *sin_phi2;
} brdf_aux_data;


typedef struct {
     void (*aux)(brdf_aux_data *aux, void *aux2, int i, int j, int n_derivs, uchar *derivs, double *p, double **p_l);
     double (*kernel)(brdf_aux_data *aux, void *aux2, enum xrtm_brdf_geometry brdf_geom, int i, int j, int k, double phi, int n_derivs, uchar *derivs, double *p, double **p_l, double *f_l);
} brdf_kernel_func_data;


#include "prototypes/xrtm_brdf_p.h"


#ifdef __cplusplus
}
#endif

#endif /* XRTM_BRDF_H */
