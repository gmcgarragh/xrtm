/*******************************************************************************
**
**    Copyright (C) 2007-2020 Greg McGarragh <greg.mcgarragh@colostate.edu>
**
**    This source code is licensed under the GNU General Public License (GPL),
**    Version 3.  See the file COPYING for more details.
**
*******************************************************************************/

#ifndef XRTM_SINGLE_A_H
#define XRTM_SINGLE_A_H

#include "xrtm_save_tree.h"
#include "xrtm_work.h"

#ifdef __cplusplus
extern "C" {
#endif


typedef struct {
     void (*free)(void *);
     double **I;
} forward_save_single_scattered_radiance_data;


#include "prototypes/xrtm_single_a_p.h"


#ifdef __cplusplus
}
#endif

#endif /* XRTM_SINGLE_A_H */
