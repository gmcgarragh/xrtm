/*******************************************************************************
**
**    Copyright (C) 2007-2020 Greg McGarragh <greg.mcgarragh@colostate.edu>
**
**    This source code is licensed under the GNU General Public License (GPL),
**    Version 3.  See the file COPYING for more details.
**
*******************************************************************************/


/*******************************************************************************
 *
 ******************************************************************************/
#define USE_LIDORT_QUAD_OUTPUT		0

#define USE_PADE_GAMMA_INIT		0

#define EIGEN_SOLVER_SYM_REAL 		EIGEN_SOLVER_SYM_REAL_EISPACK


/*******************************************************************************
 *
 ******************************************************************************/
#define XRTM_SOLVERS_EXT_EXTERNAL      (XRTM_SOLVER_DISORT | \
                                        XRTM_SOLVER_LIDORT | \
                                        XRTM_SOLVER_LRAD | \
                                        XRTM_SOLVER_POLRAD | \
                                        XRTM_SOLVER_RADIANT | \
                                        XRTM_SOLVER_RADTRAN3 | \
                                        XRTM_SOLVER_SOI | \
                                        XRTM_SOLVER_TWOSTR | \
                                        XRTM_SOLVER_VLIDORT)

#define XRTM_SOLVERS_EXT_VECTOR        (XRTM_SOLVER_LRAD | \
                                        XRTM_SOLVER_POLRAD | \
                                        XRTM_SOLVER_RADTRAN3 | \
                                        XRTM_SOLVER_VLIDORT)

#define XRTM_SOLVERS_EXT_USE_G         (XRTM_SOLVER_TWOSTR)

#define XRTM_SOLVERS_EXT_USE_COEF      (XRTM_SOLVER_DISORT | \
                                        XRTM_SOLVER_LIDORT | \
                                        XRTM_SOLVER_LRAD | \
                                        XRTM_SOLVER_POLRAD | \
                                        XRTM_SOLVER_RADIANT | \
                                        XRTM_SOLVER_RADTRAN3 | \
                                        XRTM_SOLVER_SOI | \
                                        XRTM_SOLVER_VLIDORT)

#define XRTM_SOLVERS_EXT_LINEARIZED    (XRTM_SOLVER_LIDORT | \
                                        XRTM_SOLVER_LRAD | \
                                        XRTM_SOLVER_RADIANT | \
                                        XRTM_SOLVER_VLIDORT)

#define XRTM_SOLVERS_EXT_DOUBLING      (XRTM_SOLVER_POLRAD | \
                                        XRTM_SOLVER_RADTRAN3)

#define XRTM_SOLVERS_EXT_ADDING        (XRTM_SOLVER_POLRAD | \
                                        XRTM_SOLVER_RADIANT | \
                                        XRTM_SOLVER_RADTRAN3)

#define XRTM_SOLVERS_EXT_BVP           (XRTM_SOLVER_DISORT | \
                                        XRTM_SOLVER_LIDORT | \
                                        XRTM_SOLVER_TWOSTR | \
                                        XRTM_SOLVER_VLIDORT)

#define XRTM_SOLVERS_EXT_HAS_SFI       (XRTM_SOLVER_DISORT | \
                                        XRTM_SOLVER_LIDORT | \
                                        XRTM_SOLVER_TWOSTR | \
                                        XRTM_SOLVER_VLIDORT)

#define XRTM_SOLVERS_EXT_REQUIRES_WORK (XRTM_SOLVER_LRAD)
