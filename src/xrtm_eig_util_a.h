/******************************************************************************%
**
**    Copyright (C) 2007-2012 Greg McGarragh <gregm@atmos.colostate.edu>
**
**    This source code is licensed under the GNU General Public License (GPL),
**    Version 3.  See the file COPYING for more details.
**
*******************************************************************************/

#ifndef XRTM_EIG_UTIL_A_H
#define XRTM_EIG_UTIL_A_H

#include "xrtm_save_tree.h"
#include "xrtm_work.h"

#ifdef __cplusplus
extern "C" {
#endif


typedef struct {
     void (*free)(void *);
     int flag;
     int **ip;
     double ***A;
} forward_save_eig_1n_gen_real_data;


typedef struct {
     void (*free)(void *);
     double **b;
} forward_save_eig_1n_to_2n_real_data;


typedef struct {
     void (*free)(void *);
     double  *evals;
     double **evecs;
} forward_save_eig_2n_gen_real_data;


typedef struct {
     void (*free)(void *);
     int flag;
     int **ip;
     double ***A;
} forward_save_eig_1n_gen_complex_data;

typedef struct {
     void (*free)(void *);
     double **b;
     double **vr3;
     double **vr4;
} forward_save_eig_1n_to_2n_complex_data;


typedef struct {
     void (*free)(void *);
     double  *evals_r;
     double  *evals_i;
     double **evecs;
} forward_save_eig_2n_gen_complex_data;


#include "prototypes/xrtm_eig_util_a_p.h"


#ifdef __cplusplus
}
#endif

#endif /* XRTM_EIG_UTIL_A_H */
