/******************************************************************************%
**
**    Copyright (C) 2007-2012 Greg McGarragh <gregm@atmos.colostate.edu>
**
**    This source code is licensed under the GNU General Public License (GPL),
**    Version 3.  See the file COPYING for more details.
**
*******************************************************************************/

#ifndef XRTM_EIG_RTS_A_H
#define XRTM_EIG_RTS_A_H

#include "xrtm_save_tree.h"
#include "xrtm_work.h"

#ifdef __cplusplus
extern "C" {
#endif


typedef struct {
     void (*free)(void *);
     int *X_m_ip;
     int *c_ip;
     double **X_m_lu;
     double  *lambda;
     double **a;
     double **b;
     double **c_lu;
     double **d;
     double **e;
} forward_save_calc_global_r_and_t_data;


typedef struct {
     void (*free)(void *);
     int *X_m_ip;
     int *c_ip;
     dcomplex **X_m_lu;
     dcomplex  *lambda;
     dcomplex **a;
     dcomplex **b;
     dcomplex **c_lu;
     dcomplex **d;
     dcomplex **e;
} forward_save_calc_global_r_and_t_data2;


typedef struct {
     void (*free)(void *);

     int flag;

     double   *F_p;
     double   *F_m;

     double  **tpr;
     double  **tmr;
     double  **gamma;

     double   *nu;
     double  **X_p;
     double  **X_m;

     dcomplex  *nu_c;
     dcomplex **X_p_c;
     dcomplex **X_m_c;
} forward_save_rtm_eig_rts_data;


#include "prototypes/xrtm_eig_rts_a_p.h"


#ifdef __cplusplus
}
#endif

#endif /* XRTM_EIG_RTS_A_H */
