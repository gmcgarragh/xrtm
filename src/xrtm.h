/******************************************************************************%
**
**    Copyright (C) 2007-2012 Greg McGarragh <gregm@atmos.colostate.edu>
**
**    This source code is licensed under the GNU General Public License (GPL),
**    Version 3.  See the file COPYING for more details.
**
*******************************************************************************/

#ifndef XRTM_H
#define XRTM_H

#ifdef __cplusplus
extern "C" {
#endif


#include "../config.h"


#ifdef INCLUDE_DEV_SOURCE
#include "../src2/xrtm2.h"
#endif


#define USE_PADE_CHECK_CONDITION	0

#define USE_REBUILD_STACKS		0

#define USE_SYMMETRIC_FORM		0

#ifdef HAVE_FORTRAN_COMPILER
#define EIGEN_SOLVER_GEN_REAL 		EIGEN_SOLVER_GEN_REAL_ASYMTX
#else
#define EIGEN_SOLVER_GEN_REAL 		EIGEN_SOLVER_GEN_REAL_LAPACK
#endif
#define EIGEN_SOLVER_GEN_COMPLEX 	EIGEN_SOLVER_GEN_COMPLEX_LAPACK

#define THRESHOLD_SYM_MUL_BLOCK		8 + 4

#define THRESHOLD_MU_0_SINGLARITY	1.e-5
#define THRESHOLD_OMEGA_SINGLARITY	1.e-5


#ifdef __cplusplus
}
#endif

#endif /* XRTM_H */
