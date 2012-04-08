/******************************************************************************%
**
**    Copyright (C) 2007-2012 Greg McGarragh <gregm@atmos.colostate.edu>
**
**    This source code is licensed under the GNU General Public License (GPL),
**    Version 3.  See the file COPYING for more details.
**
*******************************************************************************/

#ifndef XRTM_SUPPORT_H
#define XRTM_SUPPORT_H

#include <setjmp.h>

#include "xrtm.h"
#include "xrtm_brdf.h"
#include "xrtm_model.h"
#include "xrtm_eig_util.h"

#ifdef __cplusplus
extern "C" {
#endif


#ifdef INCLUDE_DEV_SOURCE
#include "../src2/xrtm_support2.h"
#endif


/*******************************************************************************
 *
 ******************************************************************************/
#ifndef INCLUDE_DEV_SOURCE
#define XRTM_SOLVERS_EXT_EXTERNAL   (0)
#define XRTM_SOLVERS_EXT_ALL        (0)
#define XRTM_SOLVERS_EXT_VECTOR     (0)
#define XRTM_SOLVERS_EXT_USE_G      (0)
#define XRTM_SOLVERS_EXT_USE_COEF   (0)
#define XRTM_SOLVERS_EXT_LINEARIZED (0)
#define XRTM_SOLVERS_EXT_DOUBLING   (0)
#define XRTM_SOLVERS_EXT_ADDING     (0)
#define XRTM_SOLVERS_EXT_BVP        (0)
#define XRTM_SOLVERS_EXT_SFI        (0)
#endif


#define XRTM_SOLVERS_INTERNAL   (XRTM_SOLVER_DOUB_ADD | XRTM_SOLVER_EIG_ADD | \
                                 XRTM_SOLVER_EIG_BVP | XRTM_SOLVER_MEM_BVP | \
                                 XRTM_SOLVER_PADE_ADD | XRTM_SOLVER_SINGLE | \
                                 XRTM_SOLVER_SOS | XRTM_SOLVER_TWO_OS)

#define XRTM_SOLVERS_ALL        (XRTM_SOLVERS_INTERNAL | XRTM_SOLVERS_EXT_EXTERNAL)

#define XRTM_SOLVERS_VECTOR     (XRTM_SOLVER_DOUB_ADD | XRTM_SOLVER_EIG_ADD | \
                                 XRTM_SOLVER_EIG_BVP | XRTM_SOLVER_MEM_BVP | \
                                 XRTM_SOLVER_PADE_ADD | XRTM_SOLVER_SINGLE | \
                                 XRTM_SOLVER_SOS | XRTM_SOLVER_TWO_OS | \
                                 XRTM_SOLVERS_EXT_VECTOR)

#define XRTM_SOLVERS_USE_G      (XRTM_SOLVERS_EXT_USE_G)

#define XRTM_SOLVERS_USE_COEF   (XRTM_SOLVER_DOUB_ADD | XRTM_SOLVER_EIG_ADD | \
                                 XRTM_SOLVER_EIG_BVP | XRTM_SOLVER_MEM_BVP | \
                                 XRTM_SOLVER_PADE_ADD | XRTM_SOLVER_SINGLE | \
                                 XRTM_SOLVER_SOS | XRTM_SOLVER_TWO_OS | \
                                 XRTM_SOLVERS_EXT_USE_COEF)

#define XRTM_SOLVERS_LINEARIZED (XRTM_SOLVER_DOUB_ADD | XRTM_SOLVER_EIG_ADD | \
                                 XRTM_SOLVER_EIG_BVP | XRTM_SOLVER_MEM_BVP | \
                                 XRTM_SOLVER_PADE_ADD | XRTM_SOLVER_SINGLE | \
                                 XRTM_SOLVER_SOS | XRTM_SOLVER_TWO_OS | \
                                 XRTM_SOLVERS_EXT_LINEARIZED)

#define XRTM_SOLVERS_DOUBLING   (XRTM_SOLVER_DOUB_ADD | XRTM_SOLVERS_EXT_DOUBLING)

#define XRTM_SOLVERS_ADDING     (XRTM_SOLVER_DOUB_ADD | XRTM_SOLVER_EIG_ADD | \
                                 XRTM_SOLVER_PADE_ADD | XRTM_SOLVERS_EXT_ADDING)

#define XRTM_SOLVERS_BVP        (XRTM_SOLVER_EIG_BVP | XRTM_SOLVER_MEM_BVP | \
                                 XRTM_SOLVER_TWO_OS | XRTM_SOLVER_SOS | \
                                 XRTM_SOLVERS_EXT_BVP)

#define XRTM_SOLVERS_HAS_SFI    (XRTM_SOLVER_EIG_BVP | XRTM_SOLVER_MEM_BVP | \
                                 XRTM_SOLVERS_EXT_HAS_SFI)



/*******************************************************************************
 *
 ******************************************************************************/
#define XRTM_SOLUTION_ALL (XRTM_OUTPUT_RADIANCE | XRTM_OUTPUT_RADIANCE_MEAN | \
                           XRTM_OUTPUT_FLUX | XRTM_OUTPUT_FLUX_DIVERGENCE)



/*******************************************************************************
 *
 ******************************************************************************/
#define ALLOC_ARRAY1(m, type)             ((type *)       jalloc_array1(env, m, sizeof(type)))
#define ALLOC_ARRAY1_UC(m)                ((uchar *)      jalloc_array1(env, m, sizeof(uchar)))
#define ALLOC_ARRAY1_I(m)                 ((int *)        jalloc_array1(env, m, sizeof(int)))
#define ALLOC_ARRAY1_D(m)                 ((double *)     jalloc_array1(env, m, sizeof(double)))

#define ALLOC_ARRAY2(m, n, type)          ((type **)      jalloc_array2(env, m, n, sizeof(type)))
#define ALLOC_ARRAY2_UC(m, n)             ((uchar **)     jalloc_array2(env, m, n, sizeof(uchar)))
#define ALLOC_ARRAY2_I(m, n)              ((int **)       jalloc_array2(env, m, n, sizeof(int)))
#define ALLOC_ARRAY2_D(m, n)              ((double **)    jalloc_array2(env, m, n, sizeof(double)))

#define ALLOC_ARRAY3(m, n, o, type)       ((type ***)     jalloc_array3(env, m, n, o, sizeof(type)))
#define ALLOC_ARRAY3_UC(m, n, o)          ((uchar ***)    jalloc_array3(env, m, n, o, sizeof(uchar)))
#define ALLOC_ARRAY3_I(m, n, o)           ((int ***)      jalloc_array3(env, m, n, o, sizeof(int)))
#define ALLOC_ARRAY3_D(m, n, o)           ((double ***)   jalloc_array3(env, m, n, o, sizeof(double)))

#define ALLOC_ARRAY4(m, n, o, p, type)    ((type ****)    jalloc_array4(env, m, n, o, p, sizeof(type)))
#define ALLOC_ARRAY4_UC(m, n, o, p)       ((uchar ****)   jalloc_array4(env, m, n, o, p, sizeof(uchar)))
#define ALLOC_ARRAY4_I(m, n, o, p)        ((int ****)     jalloc_array4(env, m, n, o, p, sizeof(int)))
#define ALLOC_ARRAY4_D(m, n, o, p)        ((double ****)  jalloc_array4(env, m, n, o, p, sizeof(double)))

#define ALLOC_ARRAY5(m, n, o, p, q, type) ((type *****)   jalloc_array5(env, m, n, o, p, q, sizeof(type)))
#define ALLOC_ARRAY5_UC(m, n, o, p, q)    ((uchar *)      jalloc_array5(env, m, n, o, p, q, sizeof(uchar)))
#define ALLOC_ARRAY5_I(m, n, o, p, q)     ((int *****)    jalloc_array5(env, m, n, o, p, q, sizeof(int)))
#define ALLOC_ARRAY5_D(m, n, o, p, q)     ((double *****) jalloc_array5(env, m, n, o, p, q, sizeof(double)))



/*******************************************************************************
 *
 ******************************************************************************/
#define CHECK_N_KERNESL_NOT_ZERO(D, R) do {					\
     if (D->n_kernels == 0) {							\
          eprintf("ERROR: can't set/get brdf inputs if n_kernels = 0\n");	\
          return R;								\
     }										\
} while (0)

#define CHECK_INIT_FOR_DERIVS() do {						\
     if (! (d->options & XRTM_OPTION_CALC_DERIVS)) {				\
          eprintf("ERROR: can't set linearized inputs when model is "		\
                  "not initialized to evaluate derivatives\n");			\
          return -1;								\
     }										\
} while (0)

#define CHECK_INDEX_RANGE(X, N, I_NAME, N_NAME, R) do {				\
     if (X < 0 || X >= N) {							\
          eprintf("ERROR: invalid value for %s: %d, must "			\
                  "be >= 0 and < %s\n", I_NAME, X, N_NAME);			\
          return R;								\
     }										\
} while (0)

#define GET_INPUTS_CHECK_SET_FLAGS(FLAGS, NAME, R) do {				\
     if (d->initial_inputs && ! FLAGS) {					\
          eprintf("ERROR: can't get %s until it is set\n", NAME);		\
          return R;								\
     }										\
} while (0)

#define SOLUTION_PREPROCESS()  do {						\
     if (solution_preprocess(d)) {						\
          eprintf("ERROR: solution_preprocess()\n");				\
          return -1;								\
     }										\
} while (0)


#include "prototypes/xrtm_support_p.h"


#ifdef __cplusplus
}
#endif

#endif /* XRTM_SUPPORT_H */
