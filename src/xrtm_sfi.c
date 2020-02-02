/*******************************************************************************
**
**    Copyright (C) 2007-2020 Greg McGarragh <greg.mcgarragh@colostate.edu>
**
**    This source code is licensed under the GNU General Public License (GPL),
**    Version 3.  See the file COPYING for more details.
**
*******************************************************************************/

#include <gutil.h>
#include <gmath_matrix.h>

#include "xrtm.h"
#include "xrtm_derivs.h"
#include "xrtm_matrix.h"
#include "xrtm_utility.h"


#define LEGACY_MODE 1


/*******************************************************************************
 *
 ******************************************************************************/
#define REAL

#define TYPE			double
#define TYPE_PREFIX		d
#define TYPE_TYPE_PREFIX	d
#define TYPE_POSTFIX		d
#define WORK_XX			WORK_DX
#define WORK_XXX		WORK_DXX
#define WORK_XUX		WORK_DUX
#define XEXP			exp
#define XREAL(x)		(x)

#define SFI_LAYER		sfi_layer
#define SFI			sfi
#define SFI_UP			sfi_up
#define SFI_DN			sfi_dn

#include "type_set.h"

#include "xrtm_sfi_x.c"

#include "type_unset.h"

#undef REAL

#undef TYPE
#undef TYPE_PREFIX
#undef TYPE_TYPE_PREFIX
#undef TYPE_POSTFIX
#undef WORK_XX
#undef WORK_XXX
#undef WORK_XUX
#undef XEXP
#undef XREAL

#undef SFI_LAYER
#undef SFI
#undef SFI_UP
#undef SFI_DN



/*******************************************************************************
 *
 ******************************************************************************/
#define TYPE			dcomplex
#define TYPE_PREFIX		z
#define TYPE_TYPE_PREFIX	dz
#define TYPE_POSTFIX		dc
#define WORK_XX			WORK_ZX
#define WORK_XXX		WORK_ZXX
#define WORK_XUX		WORK_ZUX
#define XEXP			cexp
#define XREAL(x)		(creal(x))

#define SFI_LAYER		sfi_layer2
#define SFI			sfi2
#define SFI_UP			sfi_up2
#define SFI_DN			sfi_dn2

#include "type_set.h"

#include "xrtm_sfi_x.c"

#include "type_unset.h"

#undef TYPE
#undef TYPE_PREFIX
#undef TYPE_TYPE_PREFIX
#undef TYPE_POSTFIX
#undef WORK_XX
#undef WORK_XXX
#undef WORK_XUX
#undef XEXP
#undef XREAL

#undef SFI_LAYER
#undef SFI
#undef SFI_UP
#undef SFI_DN
