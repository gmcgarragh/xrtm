/******************************************************************************%
**
**    Copyright (C) 1998-2012 Greg McGarragh <gregm@atmos.colostate.edu>
**
**    This source code is licensed under the GNU General Public License (GPL),
**    Version 3.  See the file COPYING for more details.
**
*******************************************************************************/

#include "gutil.h"

#include "gmath_matrix.h"


/*******************************************************************************
 * 
 ******************************************************************************/
#define prefix_ d
#define cprefix_ f
#define type_ double
#include "gvector.h"
#undef prefix_
#undef cprefix_
#undef type_

#define prefix_ z
#define cprefix_ c
#define type_ dcomplex
#include "gvector.h"
#undef prefix_
#undef cprefix_
#undef type_
