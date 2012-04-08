/******************************************************************************%
**
**    Copyright (C) 2007-2012 Greg McGarragh <gregm@atmos.colostate.edu>
**
**    This source code is licensed under the GNU General Public License (GPL),
**    Version 3.  See the file COPYING for more details.
**
*******************************************************************************/

#ifndef XRTM_FD_H
#define XRTM_FD_H

#include "xrtm_model.h"

#ifdef __cplusplus
extern "C" {
#endif


enum fd_method_type {
     FD_METHOD_FORWARD,
     FD_METHOD_BACKWARD,
     FD_METHOD_CENTRAL,

     N_FD_METHOD_TYPES
};


typedef struct {
     int n_derivs;

     double *delta;

     double **g_p;
     double ****coef_p;
     double **ltau_p;
     double **omega_p;
     double **kernel_ampfac_p;
     double ***kernel_params_p;

     xrtm_data *d;
} xrtm_fd_data;


#include "prototypes/xrtm_fd_p.h"


#ifdef __cplusplus
}
#endif

#endif /* XRTM_FD_H */
