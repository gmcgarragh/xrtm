/******************************************************************************%
**
**    Copyright (C) 2007-2012 Greg McGarragh <gregm@atmos.colostate.edu>
**
**    This source code is licensed under the GNU General Public License (GPL),
**    Version 3.  See the file COPYING for more details.
**
*******************************************************************************/

#ifndef INPUT_LOCAL_H
#define INPUT_LOCAL_H

#include <gutil.h>

#include <xrtm_fd.h>
#include <xrtm_model.h>

#include "input_util.h"

#ifdef __cplusplus
extern "C" {
#endif


#define MAX_DIMENS 8


enum lex_type {
     LEX_TYPE_STRING = 1,
     LEX_TYPE_INT,
     LEX_TYPE_DOUBLE,
     LEX_TYPE_IDENTIFIER
};


enum input_type {
     INPUT_TYPE_VOID = 128,
     INPUT_TYPE_STRING,
     INPUT_TYPE_INT,
     INPUT_TYPE_DOUBLE
};


enum input_option {
     INPUT_OPTION_FD = 256,

     N_INPUT_OPTIONS
};


enum xrtm_input {
     XRTM_INPUT_OPTIONS = 384,
     XRTM_INPUT_SOLVERS,
     XRTM_INPUT_MAX_COEF,
     XRTM_INPUT_N_QUAD,
     XRTM_INPUT_N_STOKES,
     XRTM_INPUT_N_DERIVS,
     XRTM_INPUT_N_LAYERS,
     XRTM_INPUT_N_KERNELS,
     XRTM_INPUT_N_KERNEL_QUAD,
     XRTM_INPUT_KERNELS,
     XRTM_INPUT_N_OUT_LEVELS,
     XRTM_INPUT_N_OUT_THETAS,

     XRTM_INPUT_DOUB_D_TAU,
     XRTM_INPUT_PADE_PARAMS,
     XRTM_INPUT_SOS_PARAMS,
     XRTM_INPUT_FOURIER_TOL,
     XRTM_INPUT_LAMBDA,
     XRTM_INPUT_F_0,
     XRTM_INPUT_THETA_0,
     XRTM_INPUT_PHI_0,
     XRTM_INPUT_OUT_LEVELS,
     XRTM_INPUT_OUT_TAUS,
     XRTM_INPUT_OUT_THETAS,
     XRTM_INPUT_TOP_B,
     XRTM_INPUT_TOP_B_L,
     XRTM_INPUT_PLANET_R,
     XRTM_INPUT_LEVELS_Z,
     XRTM_INPUT_LEVELS_B,
     XRTM_INPUT_LEVELS_B_L,
     XRTM_INPUT_G,
     XRTM_INPUT_G_L,
     XRTM_INPUT_CHI,
     XRTM_INPUT_COEF_L,
     XRTM_INPUT_COEF_FILES,
     XRTM_INPUT_COEF_FILES_L,
     XRTM_INPUT_OMEGA,
     XRTM_INPUT_OMEGA_L,
     XRTM_INPUT_LTAU,
     XRTM_INPUT_LTAU_L,
     XRTM_INPUT_SURFACE_B,
     XRTM_INPUT_SURFACE_B_L,
     XRTM_INPUT_KERNEL_AMPFAC,
     XRTM_INPUT_KERNEL_AMPFAC_L,
     XRTM_INPUT_KERNEL_PARAMS,
     XRTM_INPUT_KERNEL_PARAMS_L,

     XRTM_INPUT_DELTA,
     XRTM_INPUT_TOP_B_P,
     XRTM_INPUT_LEVELS_B_P,
     XRTM_INPUT_G_P,
     XRTM_INPUT_COEF_P,
     XRTM_INPUT_COEF_FILES_P,
     XRTM_INPUT_OMEGA_P,
     XRTM_INPUT_LTAU_P,
     XRTM_INPUT_SURFACE_B_P,
     XRTM_INPUT_KERNEL_AMPFAC_P,
     XRTM_INPUT_KERNEL_PARAMS_P,

     N_XRTM_INPUTS
};


typedef union {
     char *s;
     int i;
     double d;
} lex_type_data;


typedef struct {
     const char *file;
     int  line;
     const char *statement;
} locus_data;


typedef struct {
     char nl;
     int use_aligned;
     int use_array2;
     int use_dash;
     int use_equal;
} format_data;


typedef struct {
     int use_aligned;
     int use_array2;
     int use_dash;
     int use_equal;
     xrtm_data *d;
     xrtm_fd_data *fd;
     misc_data *md;
     locus_data locus;
} input_data;


#undef  GSTRUCT_BASIC_TYPE
#define GSTRUCT_NILL list_null
#undef  GSTRUCT_POINTER_TYPE
#define GSTRUCT_PREFIX f
#define GSTRUCT_TYPE input_type_data

#define DEC_BEFORE
#include <gstruct_dclist.h>
#undef DEC_BEFORE


typedef struct {
     char *name;

     int flag;
     int type;

     int order;
     int dimens[MAX_DIMENS];

     union {
          char *s;
          int i;
          double d;
          void *v;
          fdclist *list;
     } d;
} input_type_data;


#define DEC_AFTER
#include <gstruct_dclist.h>
#undef DEC_AFTER

#undef GSTRUCT_BASIC_TYPE
#undef GSTRUCT_NILL
#undef GSTRUCT_POINTER_TYPE
#undef GSTRUCT_PREFIX
#undef GSTRUCT_TYPE


#define YY_DECL int yylex(locus_data *locus, lex_type_data *lex_type)


YY_DECL;
char *get_yytext();
void yypreinclude(FILE *fp);
void yypostinclude();
void yyrewind(int r, lex_type_data *lex_type);

void input_error(locus_data *locus, const char *s, ...);


#ifdef __cplusplus
}
#endif

#endif /* INPUT_LOCAL_H */
