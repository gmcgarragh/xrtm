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


static uchar max_pow_of_two[] = {
     0, 0, 1, 1, 2, 2, 2, 2, 3, 3,	/*   0 -   9 */
     3, 3, 3, 3, 3, 3, 4, 4, 4, 4,	/*  10 -  19 */
     4, 4, 4, 4, 4, 4, 4, 4, 4, 4,	/*  20 -  29 */
     4, 4, 5, 5, 5, 5, 5, 5, 5, 5,	/*  30 -  39 */
     5, 5, 5, 5, 5, 5, 5, 5, 5, 5,	/*  40 -  49 */
     5, 5, 5, 5, 5, 5, 5, 5, 5, 5,	/*  50 -  59 */
     5, 5, 5, 5, 6, 6, 6, 6, 6, 6,	/*  60 -  69 */
     6, 6, 6, 6, 6, 6, 6, 6, 6, 6,	/*  70 -  79 */
     6, 6, 6, 6, 6, 6, 6, 6, 6, 6,	/*  80 -  89 */
     6, 6, 6, 6, 6, 6, 6, 6, 6, 6,	/*  90 -  99 */
     6, 6, 6, 6, 6, 6, 6, 6, 6, 6,	/* 100 - 109 */
     6, 6, 6, 6, 6, 6, 6, 6, 6, 6,	/* 110 - 119 */
     6, 6, 6, 6, 6, 6, 6, 6, 7, 7};	/* 120 - 129 */

static uchar powers_of_two[] = {1, 2, 4, 8, 16, 32, 64, 128};



/*******************************************************************************
 *
 ******************************************************************************/
#define is_complex_ 0
#define prefix_ d
#define type_ double
#define xreal_(x) (x)
#define xabs_ fabs
#include "gmatrix.h"
#undef is_complex_
#undef prefix_
#undef type_
#undef xreal_
#undef xabs_

#define is_complex_ 1
#define prefix_ z
#define type_ dcomplex
#define xreal_(x) (creal(x))
#define xabs_ cabs
#include "gmatrix.h"
#undef is_complex_
#undef prefix_
#undef type_
#undef xreal_
#undef xabs_
