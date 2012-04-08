/******************************************************************************%
**
**    Copyright (C) 2007-2012 Greg McGarragh <gregm@atmos.colostate.edu>
**
**    This source code is licensed under the GNU General Public License (GPL),
**    Version 3.  See the file COPYING for more details.
**
*******************************************************************************/

#ifndef XRTM_PADE_RTS_H
#define XRTM_PADE_RTS_H

#include "xrtm_save_tree.h"
#include "xrtm_work.h"

#ifdef __cplusplus
extern "C" {
#endif


#define MAX_PADE_R 32


typedef struct {
     double **a;
     double **d1;
     double *d2;
} matd1d2;


#include "prototypes/xrtm_pade_rts_p.h"


#ifdef __cplusplus
}
#endif

#endif /* XRTM_PADE_RTS_H */
