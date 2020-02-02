/*******************************************************************************
**
**    Copyright (C) 2007-2020 Greg McGarragh <greg.mcgarragh@colostate.edu>
**
**    This source code is licensed under the GNU General Public License (GPL),
**    Version 3.  See the file COPYING for more details.
**
*******************************************************************************/


#define GROUP_SOLVERS_DEV			\
     XRTM_SOLVER_DISORT,   DEF_TOL,		\
     XRTM_SOLVER_LIDORT,   DEF_TOL,		\
     XRTM_SOLVER_VLIDORT,  DEF_TOL,

#define GROUP_SOLVERS_DEV_DERIVS1		\
     XRTM_SOLVER_LIDORT,   DEF_TOL, DEF_TOL,	\
     XRTM_SOLVER_VLIDORT,  DEF_TOL, DEF_TOL,

#define GROUP_SOLVERS_DEV_DERIVS12		\
     XRTM_SOLVER_LIDORT,   DEF_TOL, DEF_TOL,


#define DEV_RANGE_OVER_N_QUAD			(XRTM_SOLVER_LRAD | XRTM_SOLVER_RADIANT | XRTM_SOLVER_RADTRAN3 | XRTM_SOLVER_SOI)

#define DEV_RANGE_OVER_N_STOKES			(XRTM_SOLVER_RADTRAN3)

#define DEV_RANGE_OVER_N_DERIVS			(XRTM_SOLVER_LRAD | XRTM_SOLVER_RADIANT)

#define DEV_RANGE_OVER_N_LAYERS			(XRTM_SOLVER_LRAD | XRTM_SOLVER_RADIANT)

#define DEV_RANGE_OVER_N_OUT_THETAS		(XRTM_SOLVER_DISORT | XRTM_SOLVER_LIDORT | XRTM_SOLVER_RADTRAN3 | XRTM_SOLVER_VLIDORT)

#define DEV_RANGE_OVER_OPTICAL_PROPERTY_INPUTS	(XRTM_SOLVER_DISORT | XRTM_SOLVER_LIDORT | XRTM_SOLVER_LRAD | XRTM_SOLVER_RADIANT | XRTM_SOLVER_RADTRAN3 | XRTM_SOLVER_SOI | XRTM_SOLVER_VLIDORT)

#define DEV_RANGE_OVER_EIGEN_SOVERS		(0)

#define DEV_OPT_OUT_RANGE_OVER_OPTIONS(MASK)															\
																				\
     if (i & XRTM_OPTION_NO_AZIMUTHAL)																\
          MASK &= 0xffffffff ^ (XRTM_SOLVER_DISORT | XRTM_SOLVER_LRAD | XRTM_SOLVER_TWOSTR | XRTM_SOLVER_VLIDORT);						\
     if ((i & XRTM_OPTION_N_T_TMS) && (i & XRTM_OPTION_NO_AZIMUTHAL))												\
          MASK &= 0xffffffff ^ (XRTM_SOLVER_LIDORT);														\
     if (i & XRTM_OPTION_QUAD_NORM_GAUS_LEG)															\
          MASK &= 0xffffffff ^ (XRTM_SOLVER_DISORT | XRTM_SOLVER_LIDORT | XRTM_SOLVER_LRAD | XRTM_SOLVER_SOI | XRTM_SOLVER_TWOSTR | XRTM_SOLVER_VLIDORT);	\
     if (i & XRTM_OPTION_QUAD_LOBATTO)																\
          MASK &= 0xffffffff ^ (XRTM_SOLVER_DISORT | XRTM_SOLVER_LIDORT | XRTM_SOLVER_LRAD | XRTM_SOLVER_SOI | XRTM_SOLVER_TWOSTR | XRTM_SOLVER_VLIDORT);	\
     if (i & XRTM_OPTION_PSA)																	\
          MASK &= 0xffffffff ^ (XRTM_SOLVER_DISORT | XRTM_SOLVER_RADTRAN3);											\
     if (i & XRTM_OPTION_SFI)																	\
          MASK &= 0xffffffff ^ (XRTM_SOLVER_LRAD | XRTM_SOLVER_RADIANT | XRTM_SOLVER_RADTRAN3 | XRTM_SOLVER_SOI);

#define DEV_RANGE_OVER_SAVING_OPTIONS		(0)


int test_errors_dev_solvers(test_data *t);
