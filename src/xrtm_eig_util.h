/*******************************************************************************
**
**    Copyright (C) 2007-2020 Greg McGarragh <greg.mcgarragh@colostate.edu>
**
**    This source code is licensed under the GNU General Public License (GPL),
**    Version 3.  See the file COPYING for more details.
**
*******************************************************************************/

#ifndef XRTM_EIG_UTIL_H
#define XRTM_EIG_UTIL_H

#include "xrtm_save_tree.h"
#include "xrtm_work.h"

#ifdef __cplusplus
extern "C" {
#endif


#ifdef INCLUDE_DEV_SOURCE
#include "../src2/xrtm_eig_util2.h"
#endif


enum eigen_solver_gen_real_type {
#ifdef HAVE_FORTRAN_COMPILER
     EIGEN_SOLVER_GEN_REAL_ASYMTX,
#endif
#ifdef HAVE_EISPACK_LIBRARY
     EIGEN_SOLVER_GEN_REAL_EISPACK,
#endif
     EIGEN_SOLVER_GEN_REAL_LAPACK,

     N_EIGEN_SOLVER_GEN_REAL_TYPES
};


enum eigen_solver_gen_complex_type {
#ifdef HAVE_EISPACK_LIBRARY
     EIGEN_SOLVER_GEN_COMPLEX_EISPACK,
#endif
     EIGEN_SOLVER_GEN_COMPLEX_LAPACK,

     N_EIGEN_SOLVER_GEN_COMPLEX_TYPES
};


#include "prototypes/xrtm_eig_util_p.h"


#ifdef __cplusplus
}
#endif

#endif /* XRTM_EIG_UTIL_H */
