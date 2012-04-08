/******************************************************************************%
**
**    Copyright (C) 2007-2012 Greg McGarragh <gregm@atmos.colostate.edu>
**
**    This source code is licensed under the GNU General Public License (GPL),
**    Version 3.  See the file COPYING for more details.
**
*******************************************************************************/

#ifndef XRTM_RADIANCE_A_H
#define XRTM_RADIANCE_A_H

#include "xrtm_save_tree.h"
#include "xrtm_work.h"

#ifdef __cplusplus
extern "C" {
#endif


typedef struct {
     void (*free)(void *);
     int *ip;
     double **P;
} forward_save_radiance_boa_all_data;


#include "prototypes/xrtm_radiance_a_p.h"


#ifdef __cplusplus
}
#endif

#endif /* XRTM_RADIANCE_A_H */
