/*******************************************************************************
**
**    Copyright (C) 2007-2020 Greg McGarragh <greg.mcgarragh@colostate.edu>
**
**    This source code is licensed under the GNU General Public License (GPL),
**    Version 3.  See the file COPYING for more details.
**
*******************************************************************************/


#define N_XRTM_DEV_SOLVERS 9

#define XRTM_DEV_SOLVER_ENUMS		\
     XRTM_SOLVER_DISORT   = (1L<<16),	\
     XRTM_SOLVER_LIDORT   = (1L<<17),	\
     XRTM_SOLVER_LRAD     = (1L<<18),	\
     XRTM_SOLVER_POLRAD   = (1L<<19),	\
     XRTM_SOLVER_RADIANT  = (1L<<20),	\
     XRTM_SOLVER_RADTRAN3 = (1L<<21),	\
     XRTM_SOLVER_SOI      = (1L<<22),	\
     XRTM_SOLVER_TWOSTR   = (1L<<23),	\
     XRTM_SOLVER_VLIDORT  = (1L<<24)

#define XRTM_DEV_SOLVER_NAMES	\
     "disort",			\
     "lidort",			\
     "lrad",			\
     "polrad",			\
     "radiant",			\
     "radtran3",		\
     "soi",			\
     "twostr",			\
     "vlidort"

#define XRTM_DEV_SOLVER_MASKS	\
     XRTM_SOLVER_DISORT,	\
     XRTM_SOLVER_LIDORT,	\
     XRTM_SOLVER_LRAD,		\
     XRTM_SOLVER_POLRAD,	\
     XRTM_SOLVER_RADIANT,	\
     XRTM_SOLVER_RADTRAN3,	\
     XRTM_SOLVER_SOI,		\
     XRTM_SOLVER_TWOSTR,	\
     XRTM_SOLVER_VLIDORT

#define XRTM_DEV_MISC_INPUT		\
     int use_lidort_quad_output;	\
					\
     int use_pade_gamma_init;		\
					\
     int eigen_solver_sym_real;
