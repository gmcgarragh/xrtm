/******************************************************************************%
**
**    Copyright (C) 2007-2012 Greg McGarragh <gregm@atmos.colostate.edu>
**
**    This source code is licensed under the GNU General Public License (GPL),
**    Version 3.  See the file COPYING for more details.
**
*******************************************************************************/

#ifndef XRTM_ADDING_H
#define XRTM_ADDING_H

#include "xrtm_save_tree.h"
#include "xrtm_work.h"

#ifdef __cplusplus
extern "C" {
#endif


enum xrtm_deriv_type {
     DERIV_TOP_P    = 1 << 2,
     DERIV_TOP_L    = 1 << 3,

     DERIV_BOTTOM_P = 1 << 0,
     DERIV_BOTTOM_L = 1 << 1
};


enum xrtm_adding_type {
     ADDING_U_U = 0,

     ADDING_U_P = 0                         | DERIV_BOTTOM_P,			/* 1  */
     ADDING_P_U = DERIV_TOP_P               | 0,				/* 4  */
     ADDING_P_P = DERIV_TOP_P               | DERIV_BOTTOM_P,			/* 5  */

     ADDING_U_L = 0                         | DERIV_BOTTOM_P | DERIV_BOTTOM_L,	/* 3  */
     ADDING_L_U = DERIV_TOP_P | DERIV_TOP_L | 0,				/* 12 */

     ADDING_P_L = DERIV_TOP_P               | DERIV_BOTTOM_P | DERIV_BOTTOM_L,	/* 7  */
     ADDING_L_P = DERIV_TOP_P | DERIV_TOP_L | DERIV_BOTTOM_P,			/* 13 */

     ADDING_L_L = DERIV_TOP_P | DERIV_TOP_L | DERIV_BOTTOM_P | DERIV_BOTTOM_L	/* 15 */
};


typedef struct {
     double **P13;
     double **P31;
     double **A13;
     double **A31;
     double **B13;
     double **B31;
     double  *C13;
     double  *C31;
} layer_add_aux_data;


#include "prototypes/xrtm_adding_p.h"


#ifdef __cplusplus
}
#endif

#endif /* XRTM_ADDING_H */
