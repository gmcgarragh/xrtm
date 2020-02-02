/*******************************************************************************
**
**    Copyright (C) 1998-2020 Greg McGarragh <greg.mcgarragh@colostate.edu>
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
#define is_complex_ 0
#define prefix_ d
#define type_ double
#define xsqrt_ sqrt
#include "gvector.h"
#undef is_complex_
#undef prefix_
#undef cprefix_
#undef type_
#undef xsqrt_

#define is_complex_ 1
#define prefix_ z
#define type_ dcomplex
#define xsqrt_ csqrt
#include "gvector.h"
#undef is_complex_
#undef prefix_
#undef cprefix_
#undef type_
#undef xsqrt_
