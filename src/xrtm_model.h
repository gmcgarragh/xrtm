/******************************************************************************%
**
**    Copyright (C) 2007-2012 Greg McGarragh <gregm@atmos.colostate.edu>
**
**    This source code is licensed under the GNU General Public License (GPL),
**    Version 3.  See the file COPYING for more details.
**
*******************************************************************************/

#ifndef XRTM_SOLVERS_MODEL_H
#define XRTM_SOLVERS_MODEL_H

#include "xrtm_interface.h"
#include "xrtm_save_tree.h"
#include "xrtm_work.h"

#ifdef __cplusplus
extern "C" {
#endif


enum xrtm_dep_mask {
     DEP_MASK_INIT        = (1<<0),
     DEP_MASK_DOUB_D_TAU  = (1<<1),
     DEP_MASK_PADE_PARAMS = (1<<2),
     DEP_MASK_SOS_PARAMS  = (1<<3),
     DEP_MASK_FOURIER_TOL = (1<<4),
     DEP_MASK_LAMBDA      = (1<<5),
     DEP_MASK_F_0         = (1<<6),
     DEP_MASK_MU_0        = (1<<7),
     DEP_MASK_PHI_0       = (1<<8),
     DEP_MASK_ULEVELS     = (1<<9),
     DEP_MASK_UTAUS       = (1<<10),
     DEP_MASK_UMUS        = (1<<11),
     DEP_MASK_TOP_B       = (1<<12),
     DEP_MASK_TOP_B_L     = (1<<13),
     DEP_MASK_PLANET_R    = (1<<14),
     DEP_MASK_LEVELS_Z    = (1<<15),
     DEP_MASK_LEVELS_B    = (1<<16),
     DEP_MASK_LEVELS_B_L  = (1<<17),
     DEP_MASK_G           = (1<<18),
     DEP_MASK_G_L         = (1<<19),
     DEP_MASK_COEF        = (1<<20),
     DEP_MASK_COEF_L      = (1<<21),
     DEP_MASK_OMEGA       = (1<<22),
     DEP_MASK_OMEGA_L     = (1<<23),
     DEP_MASK_LTAU        = (1<<24),
     DEP_MASK_LTAU_L      = (1<<25),
     DEP_MASK_SURFACE_B   = (1<<26),
     DEP_MASK_SURFACE_B_L = (1<<27),
     DEP_MASK_BRDF        = (1<<28),
     DEP_MASK_BRDF_L      = (1<<29),
/*
     DEP_MASK_MASK        = 0xffffffff
*/
     DEP_MASK_MASK        = -1
};


#include "prototypes/xrtm_model_p.h"


#ifdef __cplusplus
}
#endif

#endif /* XRTM_SOLVERS_MODEL_H */
