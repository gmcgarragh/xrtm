/*******************************************************************************
**
**    Copyright (C) 2007-2020 Greg McGarragh <greg.mcgarragh@colostate.edu>
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
     DERIV_TOP_S    = 1 << 2,
     DERIV_TOP_L    = 1 << 3,

     DERIV_BOTTOM_S = 1 << 0,
     DERIV_BOTTOM_L = 1 << 1
};


enum xrtm_adding_type {
     ADDING_U_U = 0,								/* 0 */
     ADDING_U_S =                                              DERIV_BOTTOM_S,	/* 1 */
     ADDING_U_L =                             DERIV_BOTTOM_L,			/* 2 */
     ADDING_U_B =                             DERIV_BOTTOM_L | DERIV_BOTTOM_S,	/* 3 */
     ADDING_S_U =               DERIV_TOP_S,					/* 4 */
     ADDING_S_S =               DERIV_TOP_S |                  DERIV_BOTTOM_S,	/* 5 */
     ADDING_S_L =               DERIV_TOP_S | DERIV_BOTTOM_L,			/* 6 */
     ADDING_S_B =               DERIV_TOP_S | DERIV_BOTTOM_L | DERIV_BOTTOM_S,	/* 7 */
     ADDING_L_U = DERIV_TOP_L,							/* 8 */
     ADDING_L_S = DERIV_TOP_L |                                DERIV_BOTTOM_S,	/* 9 */
     ADDING_L_L = DERIV_TOP_L |               DERIV_BOTTOM_L,			/* 10 */
     ADDING_L_B = DERIV_TOP_L |               DERIV_BOTTOM_L | DERIV_BOTTOM_S,	/* 11 */
     ADDING_B_U = DERIV_TOP_L | DERIV_TOP_S,					/* 12 */
     ADDING_B_S = DERIV_TOP_L | DERIV_TOP_S |                  DERIV_BOTTOM_S,	/* 13 */
     ADDING_B_L = DERIV_TOP_L | DERIV_TOP_S | DERIV_BOTTOM_L,			/* 14 */
     ADDING_B_B = DERIV_TOP_L | DERIV_TOP_S | DERIV_BOTTOM_L | DERIV_BOTTOM_S	/* 15 */
};



typedef struct {
     int *i13;
     int *i31;
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
