/******************************************************************************%
**
**    Copyright (C) 2007-2012 Greg McGarragh <gregm@atmos.colostate.edu>
**
**    This source code is licensed under the GNU General Public License (GPL),
**    Version 3.  See the file COPYING for more details.
**
*******************************************************************************/

#ifndef XRTM_SAVE_TREE_H
#define XRTM_SAVE_TREE_H

#include <gstruct.h>

#ifdef __cplusplus
extern "C" {
#endif


#define save_tree_retrieve_data(D, TYPE, V)     __save_tree_retrieve_data(D, sizeof(TYPE), (void **) V) 

#define save_tree_retrieve_proxy(D, S, TYPE, V) __save_tree_retrieve_proxy(D, S, sizeof(TYPE), (void **) V) 


typedef struct {
     void (*free)(void *);
} forward_save_data;


typedef struct {
     int is_proxy;
     ulong hash;
     void *d;
} save_node_data;


typedef struct {
     char *temp;
     int i;
     ulong *hash;
     gaatree *t;
} save_tree_data;


#include "prototypes/xrtm_save_tree_p.h"


#ifdef __cplusplus
}
#endif

#endif /* XRTM_SAVE_TREE_H */
