/*******************************************************************************
**
**    Copyright (C) 2007-2020 Greg McGarragh <greg.mcgarragh@colostate.edu>
**
**    This source code is licensed under the GNU General Public License (GPL),
**    Version 3.  See the file COPYING for more details.
**
*******************************************************************************/

#ifndef XRTM_WORK_H
#define XRTM_WORK_H

#ifdef __cplusplus
extern "C" {
#endif


#define WORK_ALIGNMENT	16

#define WORK_MAX_PAGES	1024

#define WORK_I_TYPE	int
#define WORK_LIST_INC	1024

#ifndef __cplusplus
#define get_work1(WORK, TYPE)                 work_get1(WORK, TYPE)
#define get_work2(WORK, TYPE, V_TYPE, DERIVS) work_get2(WORK, TYPE, V_TYPE, DERIVS)
#define get_work3(WORK, TYPE, V_TYPE, DERIVS) work_get3(WORK, TYPE, V_TYPE, DERIVS)
#else
#ifdef NO_WORK_CAST
#define get_work1 work_get1
#define get_work2 work_get2
#define get_work3 work_get3
#else
#define GET_WORK1_WORK_IX(WORK)      (int *)       work_get1(WORK, WORK_IX)
#define GET_WORK1_WORK_IX2(WORK)     (int *)       work_get1(WORK, WORK_IX2)
#define GET_WORK1_WORK_DQQ(WORK)     (double **)   work_get1(WORK, WORK_DQQ)
#define GET_WORK1_WORK_DX(WORK)      (double *)    work_get1(WORK, WORK_DX)
#define GET_WORK1_WORK_DX2(WORK)     (double *)    work_get1(WORK, WORK_DX2)
#define GET_WORK1_WORK_DXX(WORK)     (double **)   work_get1(WORK, WORK_DXX)
#define GET_WORK1_WORK_DXX2(WORK)    (double **)   work_get1(WORK, WORK_DXX2)
#define GET_WORK1_WORK_DD(WORK)      (double *)    work_get1(WORK, WORK_DD)
#define GET_WORK1_WORK_DS(WORK)      (double *)    work_get1(WORK, WORK_DS)
#define GET_WORK1_WORK_DU(WORK)      (double *)    work_get1(WORK, WORK_DU)
#define GET_WORK1_WORK_DUX(WORK)     (double **)   work_get1(WORK, WORK_DUX)
#define GET_WORK1_WORK_DDERIVS(WORK) (double *)    work_get1(WORK, WORK_DDERIVS)
#define GET_WORK1_WORK_DLAYERS(WORK) (double *)    work_get1(WORK, WORK_DLAYERS)
#define GET_WORK1_WORK_DBOTH(WORK)   (double **)   work_get1(WORK, WORK_DBOTH)
#define GET_WORK1_WORK_ZX(WORK)      (dcomplex *)  work_get1(WORK, WORK_ZX)
#define GET_WORK1_WORK_ZXX(WORK)     (dcomplex **) work_get1(WORK, WORK_ZXX)
#define GET_WORK1_WORK_ZU(WORK)      (dcomplex *)  work_get1(WORK, WORK_ZU)
#define GET_WORK1_WORK_ZUX(WORK)     (dcomplex **) work_get1(WORK, WORK_ZUX)

#define xet_work1(WORK, TYPE) GET_WORK1_##TYPE(WORK)
#define get_work1(WORK, TYPE) xet_work1(WORK, TYPE)


#define GET_WORK2_WORK_IX(WORK,      V_TYPE, DERIVS) (int **)       work_get2(WORK, WORK_IX,      V_TYPE, DERIV)
#define GET_WORK2_WORK_IX2(WORK,     V_TYPE, DERIVS) (int **)       work_get2(WORK, WORK_IX2,     V_TYPE, DERIVS)
#define GET_WORK2_WORK_DQQ(WORK,     V_TYPE, DERIVS) (double ***)   work_get2(WORK, WORK_DQQ,     V_TYPE, DERIVS)
#define GET_WORK2_WORK_DX(WORK,      V_TYPE, DERIVS) (double **)    work_get2(WORK, WORK_DX,      V_TYPE, DERIVS)
#define GET_WORK2_WORK_DX2(WORK,     V_TYPE, DERIVS) (double **)    work_get2(WORK, WORK_DX2,     V_TYPE, DERIVS)
#define GET_WORK2_WORK_DXX(WORK,     V_TYPE, DERIVS) (double ***)   work_get2(WORK, WORK_DXX,     V_TYPE, DERIVS)
#define GET_WORK2_WORK_DXX2(WORK,    V_TYPE, DERIVS) (double ***)   work_get2(WORK, WORK_DXX2,    V_TYPE, DERIVS)
#define GET_WORK2_WORK_DD(WORK,      V_TYPE, DERIVS) (double **)    work_get2(WORK, WORK_DD,      V_TYPE, DERIVS)
#define GET_WORK2_WORK_DS(WORK,      V_TYPE, DERIVS) (double **)    work_get2(WORK, WORK_DS,      V_TYPE, DERIVS)
#define GET_WORK2_WORK_DU(WORK,      V_TYPE, DERIVS) (double **)    work_get2(WORK, WORK_DU,      V_TYPE, DERIVS)
#define GET_WORK2_WORK_DUX(WORK,     V_TYPE, DERIVS) (double ***)   work_get2(WORK, WORK_DUX,     V_TYPE, DERIVS)
#define GET_WORK2_WORK_DDERIVS(WORK, V_TYPE, DERIVS) (double **)    work_get2(WORK, WORK_DDERIVS, V_TYPE, DERIVS)
#define GET_WORK2_WORK_DLAYERS(WORK, V_TYPE, DERIVS) (double **)    work_get2(WORK, WORK_DLAYERS, V_TYPE, DERIVS)
#define GET_WORK2_WORK_DBOTH(WORK,   V_TYPE, DERIVS) (double **)    work_get2(WORK, WORK_DBOTH,   V_TYPE, DERIVS)
#define GET_WORK2_WORK_ZX(WORK,      V_TYPE, DERIVS) (dcomplex **)  work_get2(WORK, WORK_ZX,      V_TYPE, DERIVS)
#define GET_WORK2_WORK_ZXX(WORK,     V_TYPE, DERIVS) (dcomplex ***) work_get2(WORK, WORK_ZXX,     V_TYPE, DERIVS)
#define GET_WORK2_WORK_ZU(WORK,      V_TYPE, DERIVS) (dcomplex **)  work_get2(WORK, WORK_ZU,      V_TYPE, DERIVS)
#define GET_WORK2_WORK_ZUX(WORK,     V_TYPE, DERIVS) (dcomplex ***) work_get2(WORK, WORK_ZUX,     V_TYPE, DERIVS)

#define xet_work2(WORK, TYPE, V_TYPE, DERIVS) GET_WORK2_##TYPE(WORK, V_TYPE, DERIVS)
#define get_work2(WORK, TYPE, V_TYPE, DERIVS) xet_work2(WORK, TYPE, V_TYPE, DERIVS)


#define GET_WORK3_WORK_IX(WORK,      V_TYPE, DERIVS) (int ***)       work_get3(WORK, WORK_IX,      V_TYPE, DERIV)
#define GET_WORK3_WORK_IX2(WORK,     V_TYPE, DERIVS) (int ***)       work_get3(WORK, WORK_IX2,     V_TYPE, DERIVS)
#define GET_WORK3_WORK_DQQ(WORK,     V_TYPE, DERIVS) (double ****)   work_get3(WORK, WORK_DQQ,     V_TYPE, DERIVS)
#define GET_WORK3_WORK_DX(WORK,      V_TYPE, DERIVS) (double ***)    work_get3(WORK, WORK_DX,      V_TYPE, DERIVS)
#define GET_WORK3_WORK_DX2(WORK,     V_TYPE, DERIVS) (double ***)    work_get3(WORK, WORK_DX2,     V_TYPE, DERIVS)
#define GET_WORK3_WORK_DXX(WORK,     V_TYPE, DERIVS) (double ****)   work_get3(WORK, WORK_DXX,     V_TYPE, DERIVS)
#define GET_WORK3_WORK_DXX2(WORK,    V_TYPE, DERIVS) (double ****)   work_get3(WORK, WORK_DXX2,    V_TYPE, DERIVS)
#define GET_WORK3_WORK_DD(WORK,      V_TYPE, DERIVS) (double ***)    work_get3(WORK, WORK_DD,      V_TYPE, DERIVS)
#define GET_WORK3_WORK_DS(WORK,      V_TYPE, DERIVS) (double ***)    work_get3(WORK, WORK_DS,      V_TYPE, DERIVS)
#define GET_WORK3_WORK_DU(WORK,      V_TYPE, DERIVS) (double ***)    work_get3(WORK, WORK_DU,      V_TYPE, DERIVS)
#define GET_WORK3_WORK_DUX(WORK,     V_TYPE, DERIVS) (double ****)   work_get3(WORK, WORK_DUX,     V_TYPE, DERIVS)
#define GET_WORK3_WORK_DDERIVS(WORK, V_TYPE, DERIVS) (double ***)    work_get3(WORK, WORK_DDERIVS, V_TYPE, DERIVS)
#define GET_WORK3_WORK_DLAYERS(WORK, V_TYPE, DERIVS) (double ***)    work_get3(WORK, WORK_DLAYERS, V_TYPE, DERIVS)
#define GET_WORK3_WORK_DBOTH(WORK,   V_TYPE, DERIVS) (double ***)    work_get3(WORK, WORK_DBOTH,   V_TYPE, DERIVS)
#define GET_WORK3_WORK_ZX(WORK,      V_TYPE, DERIVS) (dcomplex ***)  work_get3(WORK, WORK_ZX,      V_TYPE, DERIVS)
#define GET_WORK3_WORK_ZXX(WORK,     V_TYPE, DERIVS) (dcomplex ****) work_get3(WORK, WORK_ZXX,     V_TYPE, DERIVS)
#define GET_WORK3_WORK_ZU(WORK,      V_TYPE, DERIVS) (dcomplex ***)  work_get3(WORK, WORK_ZU,      V_TYPE, DERIVS)
#define GET_WORK3_WORK_ZUX(WORK,     V_TYPE, DERIVS) (dcomplex ****) work_get3(WORK, WORK_ZUX,     V_TYPE, DERIVS)

#define xet_work3(WORK, TYPE, V_TYPE, DERIVS) GET_WORK3_##TYPE(WORK, V_TYPE, DERIVS)
#define get_work3(WORK, TYPE, V_TYPE, DERIVS) xet_work3(WORK, TYPE, V_TYPE, DERIVS)
#endif
#endif

enum work_type {
     WORK_IX = 0,
     WORK_IX2,

     WORK_DQQ,	/* n_quad_v */

     WORK_DX,	/* n_quad_v_x */
     WORK_DX2,
     WORK_DXX,
     WORK_DXX2,

     WORK_DD,	/* n_quad_v_d */

     WORK_DS,

     WORK_DU,
     WORK_DUX,

     WORK_DDERIVS,
     WORK_DLAYERS,
     WORK_DBOTH,

     WORK_ZX,
     WORK_ZXX,
     WORK_ZU,
     WORK_ZUX,

     N_WORK_TYPES
};

enum work_v_type {
     WORK_DERIVS_V,
     WORK_LAYERS_V,
     WORK_BOTH_V,

     N_WORK_TYPES_V
};


typedef struct {
     int n_quad_v;
     int n_quad_v_x;
     int n_stokes;
     int n_derivs;
     int n_layers;
     int n_umus_v;

     uchar *pages[WORK_MAX_PAGES];

     size_t min_page_size;

     size_t n_bytes[WORK_MAX_PAGES];

     int n_list  [N_WORK_TYPES];
     int n_list_v[N_WORK_TYPES_V];

     int max_list  [N_WORK_TYPES];
     int max_list_v[N_WORK_TYPES_V];

     void **lists   [N_WORK_TYPES];
     void ***lists_v[N_WORK_TYPES_V];
} work_shared;


typedef struct {
     size_t i_page;
     size_t i_byte;

     WORK_I_TYPE i_list  [N_WORK_TYPES];
     WORK_I_TYPE i_list_v[N_WORK_TYPES_V];

     work_shared *w;
} work_data;


#include "prototypes/xrtm_work_p.h"


#ifdef __cplusplus
}
#endif

#endif /* XRTM_WORK_H */
