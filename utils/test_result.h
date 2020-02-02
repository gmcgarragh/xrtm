/*******************************************************************************
**
**    Copyright (C) 2007-2020 Greg McGarragh <greg.mcgarragh@colostate.edu>
**
**    This source code is licensed under the GNU General Public License (GPL),
**    Version 3.  See the file COPYING for more details.
**
*******************************************************************************/

#ifndef TEST_RESULT_H
#define TEST_RESULT_H

#include <glist.h>

#include <xrtm_interface.h>

#ifdef __cplusplus
extern "C" {
#endif


typedef struct {
      char *name;
      struct list_data *prev;
      struct list_data *next;

      int index;
      int index2;

      char *cmd;
} test_cmd_data;


typedef struct {
      enum xrtm_solver_mask mask;
      double x_p;
      double x_m;
      double x_p_delta;
      double x_m_delta;
} test_result_solver_data;


typedef struct {
      char *name;
      struct list_data *prev;
      struct list_data *next;

      int index;
      int index2;
      int i_out_level;
      int i_deriv;
      double out_theta;
      double out_phi;
      int i_stokes;

      int n_solvers;

      test_result_solver_data *solvers;
} test_result_data;


#include "prototypes/test_result_p.h"


#ifdef __cplusplus
}
#endif

#endif /* TEST_RESULT_H */
