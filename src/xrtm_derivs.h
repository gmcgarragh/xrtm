/*******************************************************************************
**
**    Copyright (C) 2007-2020 Greg McGarragh <greg.mcgarragh@colostate.edu>
**
**    This source code is licensed under the GNU General Public License (GPL),
**    Version 3.  See the file COPYING for more details.
**
*******************************************************************************/

#ifndef XRTM_DERIVS_H
#define XRTM_DERIVS_H

#ifdef __cplusplus
extern "C" {
#endif


typedef struct {
     uchar **layers;
     uchar **total;
     uchar **stacks;
     uchar **beam;
     uchar **thermal;
     uchar **sources;

     uchar **adding_up;
     uchar **adding_down;

     uchar **layers_m;
     uchar **total_m;
     uchar **stacks_m;
     uchar **beam_m;
     uchar **thermal_m;
     uchar **sources_m;

     uchar **adding_up_m;
     uchar **adding_down_m;

     uchar *layers_union;
     uchar *total_union;
     uchar *stacks_union;
     uchar *beam_union;
     uchar *thermal_union;
     uchar *sources_union;

     uchar *adding_up_union;
     uchar *adding_down_union;

     uchar *layers_m_union;
     uchar *total_m_union;
     uchar *stacks_m_union;
     uchar *beam_m_union;
     uchar *thermal_m_union;
     uchar *sources_m_union;

     uchar *adding_up_m_union;
     uchar *adding_down_m_union;
} derivs_data;


#include "prototypes/xrtm_derivs_p.h"


#ifdef __cplusplus
}
#endif

#endif /* XRTM_DERIVS_H */
