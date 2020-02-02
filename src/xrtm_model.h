/*******************************************************************************
**
**    Copyright (C) 2007-2020 Greg McGarragh <greg.mcgarragh@colostate.edu>
**
**    This source code is licensed under the GNU General Public License (GPL),
**    Version 3.  See the file COPYING for more details.
**
*******************************************************************************/

#ifndef XRTM_MODEL_H
#define XRTM_MODEL_H

#include "xrtm_interface.h"
#include "xrtm_save_tree.h"
#include "xrtm_work.h"

#ifdef __cplusplus
extern "C" {
#endif


enum xrtm_dep_mask {
     DEP_MASK_INIT        = (1L<<0),
     DEP_MASK_DOUB_D_TAU  = (1L<<1),
     DEP_MASK_PADE_PARAMS = (1L<<2),
     DEP_MASK_SOS_PARAMS  = (1L<<3),
     DEP_MASK_FOURIER_TOL = (1L<<4),
     DEP_MASK_LAMBDA      = (1L<<5),
     DEP_MASK_PLANET_R    = (1L<<6),
     DEP_MASK_LEVELS_Z    = (1L<<7),
     DEP_MASK_ULEVELS     = (1L<<8),
     DEP_MASK_UTAUS       = (1L<<9),
     DEP_MASK_UMUS        = (1L<<10),
     DEP_MASK_F_0         = (1L<<11),
     DEP_MASK_MU_0        = (1L<<12),
     DEP_MASK_PHI_0       = (1L<<13),
     DEP_MASK_F_ISO_TOP   = (1L<<14),
     DEP_MASK_F_ISO_TOP_L = (1L<<15),
     DEP_MASK_F_ISO_BOT   = (1L<<16),
     DEP_MASK_F_ISO_BOT_L = (1L<<17),
     DEP_MASK_LEVELS_B    = (1L<<18),
     DEP_MASK_LEVELS_B_L  = (1L<<19),
     DEP_MASK_SURFACE_B   = (1L<<20),
     DEP_MASK_SURFACE_B_L = (1L<<21),
     DEP_MASK_G           = (1L<<22),
     DEP_MASK_G_L         = (1L<<23),
     DEP_MASK_COEF        = (1L<<24),
     DEP_MASK_COEF_L      = (1L<<25),
     DEP_MASK_OMEGA       = (1L<<26),
     DEP_MASK_OMEGA_L     = (1L<<27),
     DEP_MASK_LTAU        = (1L<<28),
     DEP_MASK_LTAU_L      = (1L<<29),
     DEP_MASK_BRDF        = (1L<<30),
     DEP_MASK_BRDF_L      = (1L<<31),
/*
     DEP_MASK_MASK        = 0xffffffff
*/
     DEP_MASK_MASK        = -1L
};


#include "prototypes/xrtm_model_p.h"


#ifdef __cplusplus
}
#endif

#endif /* XRTM_MODEL_H */
