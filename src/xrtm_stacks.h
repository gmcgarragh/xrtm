/*******************************************************************************
**
**    Copyright (C) 2007-2020 Greg McGarragh <greg.mcgarragh@colostate.edu>
**
**    This source code is licensed under the GNU General Public License (GPL),
**    Version 3.  See the file COPYING for more details.
**
*******************************************************************************/

#ifndef XRTM_STACKS_H
#define XRTM_STACKS_H

#include "xrtm_adding.h"

#ifdef __cplusplus
extern "C" {
#endif


typedef struct stack_data_ {
     int i1;
     int i2;

     int stack;
     int count;

     int new_[1024];
     int touched;
     int allocated;

     double **S_p;
     double **S_m;

     double **Sl_p;
     double **Sl_m;

     double ***R_p;
     double ***T_p;
     double ***R_m;
     double ***T_m;

     double ***S_p_l;
     double ***S_m_l;

     double ***Sl_p_l;
     double ***Sl_m_l;

     double ****R_p_l;
     double ****T_p_l;
     double ****R_m_l;
     double ****T_m_l;

     layer_add_aux_data *d;

     struct stack_data_ *p1;
     struct stack_data_ *p2;
} stack_data;


#include "prototypes/xrtm_stacks_p.h"


#ifdef __cplusplus
}
#endif

#endif /* XRTM_STACKS_H */
