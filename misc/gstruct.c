/******************************************************************************%
**
**    Copyright (C) 1998-2012 Greg McGarragh <gregm@atmos.colostate.edu>
**
**    This source code is licensed under the GNU General Public License (GPL),
**    Version 3.  See the file COPYING for more details.
**
*******************************************************************************/

#include "gutil.h"
#include "gstruct.h"


#ifdef GSTRUCT_BASIC_TYPE
#undef GSTRUCT_BASIC_TYPE
#endif
#define GSTRUCT_BASIC_TYPE

#ifdef GSTRUCT_NILL
#undef GSTRUCT_NILL
#endif
#define GSTRUCT_NILL NULL

#ifdef GSTRUCT_POINTER_TYPE
#undef GSTRUCT_POINTER_TYPE
#endif
#define GSTRUCT_POINTER_TYPE

#ifdef GSTRUCT_PREFIX
#undef GSTRUCT_PREFIX
#endif
#define GSTRUCT_PREFIX g

#ifdef GSTRUCT_TYPE
#undef GSTRUCT_TYPE
#endif
#define GSTRUCT_TYPE void *

/*
#ifndef HTREE_USE_SUB_LIST
#define HTREE_USE_SUB_LIST
#endif
*/
#ifndef HTREE_USE_SUB_TREE
#define HTREE_USE_SUB_TREE
#endif

#define IMPLEMENTATION
#include "gstruct_aatree.h"

/*
#include "gstruct_btree.h"
#include "gstruct_dlist.h"
#include "gstruct_dclist.h"
#include "gstruct_htree.h"
#include "gstruct_queue.h"
#include "gstruct_slist.h"
#include "gstruct_sclist.h"
#include "gstruct_stack.h"

#include "gstruct_extra.h"
*/
#undef IMPLEMENTATION
/*
#undef HTREE_USE_SUB_LIST
*/
#undef HTREE_USE_SUB_TREE


#undef GSTRUCT_BASIC
#undef GSTRUCT_NILL
#undef GSTRUCT_POINTER_TYPE
#undef GSTRUCT_PREFIX
#undef GSTRUCT_TYPE


void stdfree(void *v) {

     free(v);
}
