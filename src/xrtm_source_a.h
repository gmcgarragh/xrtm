/******************************************************************************%
**
**    Copyright (C) 2007-2012 Greg McGarragh <gregm@atmos.colostate.edu>
**
**    This source code is licensed under the GNU General Public License (GPL),
**    Version 3.  See the file COPYING for more details.
**
*******************************************************************************/

#ifndef XRTM_SOURCE_A_H
#define XRTM_SOURCE_A_H

#include "xrtm_save_tree.h"
#include "xrtm_work.h"

#ifdef __cplusplus
extern "C" {
#endif


typedef struct {
     void (*free)(void *);
     int *ip;
     double  *d;
     double  *e;
     double  *p;
     double  *h;
     double **f;
} forward_save_build_source_vectors_1n_data;


#include "prototypes/xrtm_source_a_p.h"


#ifdef __cplusplus
}
#endif

#endif /* XRTM_SOURCE_A_H */
