/******************************************************************************%
**
**    Copyright (C) 2007-2012 Greg McGarragh <gregm@atmos.colostate.edu>
**
**    This source code is licensed under the GNU General Public License (GPL),
**    Version 3.  See the file COPYING for more details.
**
*******************************************************************************/

#ifndef XRTM_ADDING_A_H
#define XRTM_ADDING_A_H

#include "xrtm_save_tree.h"
#include "xrtm_work.h"

#ifdef __cplusplus
extern "C" {
#endif


typedef struct {
     void (*free)(void *);
     int     *i13;
     double  *C13;
     double **P13;
     double **A13;
     double **B31;
} forward_save_layer_add_ref_data;


typedef struct {
     void (*free)(void *);

     int     *i31;
     double  *C31;
     double **P31;
     double **A31;
     double **B13;

     int     *i13;
     double  *C13;
     double **P13;
     double **A13;
     double **B31;
} forward_save_layer_add_all_data;


#include "prototypes/xrtm_adding_a_p.h"


#ifdef __cplusplus
}
#endif

#endif /* XRTM_ADDING_A_H */
