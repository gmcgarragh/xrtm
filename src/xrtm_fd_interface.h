/*******************************************************************************
**
**    Copyright (C) 2007-2020 Greg McGarragh <greg.mcgarragh@colostate.edu>
**
**    This source code is licensed under the GNU General Public License (GPL),
**    Version 3.  See the file COPYING for more details.
**
*******************************************************************************/

#ifndef XRTM_FD_INTERFACE_H
#define XRTM_FD_INTERFACE_H

#include "xrtm_interface.h"
#include "xrtm_model.h"

#ifdef __cplusplus
extern "C" {
#endif


enum xrtm_fd_method_type {
     XRTM_FD_METHOD_FORWARD,
     XRTM_FD_METHOD_BACKWARD,
     XRTM_FD_METHOD_CENTRAL,

     N_XRTM_FD_METHOD_TYPES
};


typedef struct {
     int n_derivs;

     double *delta;

     double *F_iso_top_p;
     double *F_iso_bot_p;

     double **levels_b_p;
     double *surface_b_p;

     double **g_p;
     double ****coef_p;
     double **ltau_p;
     double **omega_p;

     double **kernel_ampfac_p;
     double ***kernel_params_p;

     xrtm_data *d;
} xrtm_fd_data;


#include "prototypes/xrtm_fd_interface_p.h"


#ifdef __cplusplus
}
#endif

#endif /* XRTM_FD_INTERFACE_H */
