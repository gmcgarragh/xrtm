/******************************************************************************%
**
**    Copyright (C) 2007-2012 Greg McGarragh <gregm@atmos.colostate.edu>
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


#define MAX_KERNEL_PARAMS 3

enum xrtm_kernel_type {
     XRTM_KERNEL_LAMBERTIAN,
     XRTM_KERNEL_ROUJEAN,
     XRTM_KERNEL_LI_SPARSE,
     XRTM_KERNEL_LI_DENSE,
     XRTM_KERNEL_ROSS_THIN,
     XRTM_KERNEL_ROSS_THICK,
     XRTM_KERNEL_HAPKE,
     XRTM_KERNEL_RAHMAN,
     XRTM_KERNEL_COX_MUNK,

     N_XRTM_KERNEL_TYPES
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


#include "prototypes/xrtm_brdf_p.h"


#ifdef __cplusplus
}
#endif

#endif /* XRTM_BRDF_H */
