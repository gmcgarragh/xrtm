/*******************************************************************************
**
**    Copyright (C) 2007-2020 Greg McGarragh <greg.mcgarragh@colostate.edu>
**
**    This source code is licensed under the GNU General Public License (GPL),
**    Version 3.  See the file COPYING for more details.
**
*******************************************************************************/

#include <gutil.h>
#include <gindex_name_value.h>

#include "xrtm_fd_interface.h"
#include "xrtm_interface.h"
#include "xrtm_support.h"


/*******************************************************************************
 *
 ******************************************************************************/
static const char *xrtm_fd_method_names[] = {
     "forward",
     "backward",
     "central"
};

static long xrtm_fd_method_types[] = {
     XRTM_FD_METHOD_FORWARD,
     XRTM_FD_METHOD_BACKWARD,
     XRTM_FD_METHOD_CENTRAL
};


GINDEX_NAME_VALUE_TEMPLATE(xrtm_fd_method, "xrtm fd_method", N_XRTM_FD_METHOD_TYPES)


/*******************************************************************************
 *
 ******************************************************************************/
int xrtm_fd_create(xrtm_fd_data *fd, xrtm_data *xrtm, int n_derivs) {

     int r;

     jmp_buf env;


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     fd->d        = xrtm;

     fd->n_derivs = n_derivs;


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (fd->d->options & XRTM_OPTION_CALC_DERIVS) {
          fprintf(stderr, "ERROR: can't use xrtm option \"calc_derivs\" with xrtm_fd\n");
          return -1;

     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if ((r = setjmp(env)) == 0) {


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     fd->delta = ALLOC_ARRAY1_D(fd->n_derivs);

     fd->F_iso_top_p = ALLOC_ARRAY1_D(fd->n_derivs);
     fd->F_iso_bot_p = ALLOC_ARRAY1_D(fd->n_derivs);

     if (fd->d->options & XRTM_OPTION_SOURCE_THERMAL)
          fd->levels_b_p = ALLOC_ARRAY2_D(fd->d->n_layers + 1, fd->n_derivs);

     fd->ltau_p  = ALLOC_ARRAY2_D(fd->d->n_layers, fd->n_derivs);
     fd->omega_p = ALLOC_ARRAY2_D(fd->d->n_layers, fd->n_derivs);

     if (fd->d->solvers & XRTM_SOLVERS_USE_G)
          fd->g_p    = ALLOC_ARRAY2_D(fd->d->n_layers, fd->n_derivs);
     if (fd->d->solvers & XRTM_SOLVERS_USE_COEF)
          fd->coef_p = ALLOC_ARRAY4_D(fd->d->n_layers, fd->n_derivs, fd->d->n_elem, fd->d->n_coef);

     if (fd->d->options & XRTM_OPTION_SOURCE_THERMAL)
          fd->surface_b_p = ALLOC_ARRAY1_D(fd->n_derivs);

     if (fd->d->n_kernels > 0) {
          fd->kernel_ampfac_p = ALLOC_ARRAY2_D(fd->d->n_kernels, fd->n_derivs);

          fd->kernel_params_p = ALLOC_ARRAY3_D(fd->d->n_kernels, fd->n_derivs, MAX_KERNEL_PARAMS);
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     }
     else if (r < 0)
          return -1;


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     init_array1_d(fd->F_iso_top_p, fd->n_derivs, 0.);
     init_array1_d(fd->F_iso_bot_p, fd->n_derivs, 0.);

     if (fd->d->options & XRTM_OPTION_SOURCE_THERMAL)
          init_array2_d(fd->levels_b_p, fd->d->n_layers + 1, fd->n_derivs, 0.);

     init_array2_d(fd->ltau_p, fd->d->n_layers, fd->n_derivs, 0.);
     init_array2_d(fd->omega_p, fd->d->n_layers, fd->n_derivs, 0.);

     if (fd->d->solvers & XRTM_SOLVERS_USE_G)
          init_array2_d(fd->g_p, fd->d->n_layers, fd->n_derivs, 0.);
     if (fd->d->solvers & XRTM_SOLVERS_USE_COEF)
          init_array4_d(fd->coef_p, fd->d->n_layers, fd->n_derivs, fd->d->n_elem, fd->d->n_coef, 0.);

     if (fd->d->options & XRTM_OPTION_SOURCE_THERMAL)
          init_array1_d(fd->surface_b_p, fd->n_derivs, 0.);

     if (fd->d->n_kernels > 0) {
          init_array2_d(fd->kernel_ampfac_p, fd->d->n_kernels, fd->n_derivs, 0.);

          init_array3_d(fd->kernel_params_p, fd->d->n_kernels, fd->n_derivs, MAX_KERNEL_PARAMS, 0.);
     }


     return 0;
}



/*******************************************************************************
 *
 ******************************************************************************/
int xrtm_fd_destroy(xrtm_fd_data *d) {

     free_array1_d(d->delta);

     free_array1_d(d->F_iso_top_p);
     free_array1_d(d->F_iso_bot_p);

     if (d->d->options & XRTM_OPTION_SOURCE_THERMAL)
          free_array2_d(d->levels_b_p);

     free_array2_d(d->ltau_p);
     free_array2_d(d->omega_p);

     if (d->d->solvers & XRTM_SOLVERS_USE_G)
          free_array2_d(d->g_p);
     if (d->d->solvers & XRTM_SOLVERS_USE_COEF)
          free_array4_d(d->coef_p);

     if (d->d->options & XRTM_OPTION_SOURCE_THERMAL)
          free_array1_d(d->surface_b_p);

     if (d->d->n_kernels > 0) {
          free_array2_d(d->kernel_ampfac_p);

          free_array3_d(d->kernel_params_p);
     }

     return 0;
}



/*******************************************************************************
 *
 ******************************************************************************/
/*
static int check_inputs_all_set(xrtm_fd_data *d) {

     return 0;
}



static int solution_preprocess(xrtm_fd_data *d) {

     return 0;
}
*/


/*******************************************************************************
 *
 ******************************************************************************/
int xrtm_fd_get_n_derivs(xrtm_fd_data *d) {

     return d->n_derivs;
}



/*******************************************************************************
 *
 ******************************************************************************/
int xrtm_fd_set_delta(xrtm_fd_data *d, double *delta) {

     int i;

     for (i = 0; i < d->n_derivs; ++i)
          d->delta[i] = delta[i];

     return 0;
}



int xrtm_fd_get_delta(xrtm_fd_data *d, double *delta) {

     int i;

     for (i = 0; i < d->n_derivs; ++i)
         delta[i] = d->delta[i];

     return 0;
}


/*******************************************************************************
 *
 ******************************************************************************/
int xrtm_fd_set_F_iso_top_p_1(xrtm_fd_data *d, int i_deriv, double F_iso_top_p) {

     CHECK_INDEX_RANGE(i_deriv, d->n_derivs, "i_deriv", "n_derivs", -1);

     d->F_iso_top_p[i_deriv] = F_iso_top_p;
/*
     set_deps(d, DEP_MASK_F_iso_top_L, 0, d->n_layers);
*/
     return 0;
}



int xrtm_fd_set_F_iso_top_p_n(xrtm_fd_data *d, double *F_iso_top_p) {

     int i;

     for (i = 0; i < d->n_derivs; ++i)
          d->F_iso_top_p[i] = F_iso_top_p[i];
/*
     set_deps(d, DEP_MASK_F_iso_top_L, 0, d->n_layers);
*/
     return 0;
}



double xrtm_fd_get_F_iso_top_p(xrtm_fd_data *d, int i_deriv) {

     CHECK_INDEX_RANGE(i_deriv, d->n_derivs, "i_deriv", "n_derivs", -1);

     return d->F_iso_top_p[i_deriv];
}



/*******************************************************************************
 *
 ******************************************************************************/
int xrtm_fd_set_F_iso_bot_p_1(xrtm_fd_data *d, int i_deriv, double F_iso_bot_p) {

     CHECK_INDEX_RANGE(i_deriv, d->n_derivs, "i_deriv", "n_derivs", -1);

     d->F_iso_bot_p[i_deriv] = F_iso_bot_p;
/*
     set_deps(d, DEP_MASK_F_iso_top_L, 0, d->n_layers);
*/
     return 0;
}



int xrtm_fd_set_F_iso_bot_p_n(xrtm_fd_data *d, double *F_iso_bot_p) {

     int i;

     for (i = 0; i < d->n_derivs; ++i)
          d->F_iso_bot_p[i] = F_iso_bot_p[i];
/*
     set_deps(d, DEP_MASK_F_iso_top_L, 0, d->n_layers);
*/
     return 0;
}



double xrtm_fd_get_F_iso_bot_p(xrtm_fd_data *d, int i_deriv) {

     CHECK_INDEX_RANGE(i_deriv, d->n_derivs, "i_deriv", "n_derivs", -1);

     return d->F_iso_bot_p[i_deriv];
}



/*******************************************************************************
 *
 ******************************************************************************/
int xrtm_fd_set_levels_b_p_1(xrtm_fd_data *d, int i_deriv, double *levels_b_p) {

     int i;
     int j;

     if (! (d->d->options & XRTM_OPTION_SOURCE_THERMAL)) {
          fprintf(stderr, "ERROR: model not initialized for thermal sources\n");
          return XRTM_INT_ERROR;
     }
/*
     CHECK_INIT_FOR_DERIVS();
*/
     for (i = 0; i < d->d->n_layers + 1; ++i)
          d->levels_b_p[i][j] = levels_b_p[i_deriv];
/*
     set_deps(d, DEP_MASK_LEVELS_B_L, 0, d->n_layers);
*/
     return 0;
}



int xrtm_fd_set_levels_b_p_n(xrtm_fd_data *d, double **levels_b_p) {

     int i;
     int j;

     if (! (d->d->options & XRTM_OPTION_SOURCE_THERMAL)) {
          fprintf(stderr, "ERROR: model not initialized for thermal sources\n");
          return XRTM_INT_ERROR;
     }
/*
     CHECK_INIT_FOR_DERIVS();
*/
     for (i = 0; i < d->d->n_layers + 1; ++i) {
          for (j = 0; j < d->n_derivs; ++j)
               d->levels_b_p[i][j] = levels_b_p[i][j];
     }
/*
     set_deps(d, DEP_MASK_LEVELS_B_L, 0, d->n_layers);
*/
     return 0;
}



double xrtm_fd_get_levels_b_p(xrtm_fd_data *d, int i_deriv, double *levels_b_p) {

     int i;

     if (! (d->d->options & XRTM_OPTION_SOURCE_THERMAL)) {
          fprintf(stderr, "ERROR: model not initialized for thermal sources\n");
          return XRTM_INT_ERROR;
     }

     CHECK_INDEX_RANGE(i_deriv, d->n_derivs, "i_deriv", "n_derivs", -1);

     for (i = 0; i < d->d->n_layers + 1; ++i)
          levels_b_p[i] = d->levels_b_p[i][i_deriv];

     return 0;
}



/*******************************************************************************
 *
 ******************************************************************************/
int xrtm_fd_set_surface_b_p_1(xrtm_fd_data *d, int i_deriv, double surface_b_p) {

     if (! (d->d->options & XRTM_OPTION_SOURCE_THERMAL)) {
          fprintf(stderr, "ERROR: model not initialized for thermal sources\n");
          return XRTM_INT_ERROR;
     }

     CHECK_INDEX_RANGE(i_deriv, d->n_derivs, "i_deriv", "n_derivs", -1);

     d->surface_b_p[i_deriv] = surface_b_p;
/*
     set_deps(d, DEP_MASK_SURFACE_B_L, 0, d->n_layers);
*/
     return 0;
}



int xrtm_fd_set_surface_b_p_n(xrtm_fd_data *d, double *surface_b_p) {

     int i;

     if (! (d->d->options & XRTM_OPTION_SOURCE_THERMAL)) {
          fprintf(stderr, "ERROR: model not initialized for thermal sources\n");
          return XRTM_INT_ERROR;
     }

     for (i = 0; i < d->n_derivs; ++i)
          d->surface_b_p[i] = surface_b_p[i];
/*
     set_deps(d, DEP_MASK_SURFACE_B_L, 0, d->n_layers);
*/
     return 0;
}



double xrtm_fd_get_surface_b_p(xrtm_fd_data *d, int i_deriv) {

     if (! (d->d->options & XRTM_OPTION_SOURCE_THERMAL)) {
          fprintf(stderr, "ERROR: model not initialized for thermal sources\n");
          return XRTM_DBL_ERROR;
     }

     CHECK_INDEX_RANGE(i_deriv, d->n_derivs, "i_deriv", "n_derivs", -1);

     return d->surface_b_p[i_deriv];

     return 0;
}



/*******************************************************************************
 *
 ******************************************************************************/
static int xrtm_fd_set_g_x_p(xrtm_fd_data *d, void *g_p, int i_layer,
                             int n_layers, int i_deriv, int n_derivs, int type) {
/*
     int deps;
*/
     if (! (d->d->solvers & XRTM_SOLVERS_USE_G)) {
          fprintf(stderr, "ERROR: no two-stream solvers have been specified\n");
          return -1;
     }

     CHECK_INDEX_RANGE(i_layer, d->d->n_layers, "i_layer", "n_layers", -1);
     CHECK_INDEX_RANGE(i_deriv,    d->n_derivs, "i_deriv", "n_derivs", -1);

     if (xrtm_set_layer_x_l(d->d, d->g_p, g_p,
                            i_layer, n_layers, i_deriv, n_derivs, type)) {
          fprintf(stderr, "ERROR: xrtm_set_layer_x_l()\n");
          return 1;
     }
/*
     deps = DEP_MASK_G_L;

     set_deps(d, deps, i_layer, n_layers);
*/
     return 0;
}



int xrtm_fd_set_g_p_11(xrtm_fd_data *d, int i_layer, int i_deriv, double g_p) {

     if (xrtm_fd_set_g_x_p(d, &g_p, i_layer, i_layer + 1, i_deriv, i_deriv + 1, 0)) {
          fprintf(stderr, "ERROR: xrtm_fd_set_g_x_p()\n");
          return -1;
     }

     return 0;
}



int xrtm_fd_set_g_p_n1(xrtm_fd_data *d, int i_deriv, double *g_p) {

     if (xrtm_fd_set_g_x_p(d, g_p, 0, d->d->n_layers, i_deriv, i_deriv + 1, 1)) {
          fprintf(stderr, "ERROR: xrtm_fd_set_g_x_p()\n");
          return -1;
     }

     return 0;
}



int xrtm_fd_set_g_p_1n(xrtm_fd_data *d, int i_layer, double *g_p) {

     if (xrtm_fd_set_g_x_p(d, g_p, i_layer, i_layer + 1, 0, d->n_derivs, 2)) {
          fprintf(stderr, "ERROR: xrtm_fd_set_g_x_p()\n");
          return -1;
     }

     return 0;
}



int xrtm_fd_set_g_p_nn(xrtm_fd_data *d, double **g_p) {

     if (xrtm_fd_set_g_x_p(d, g_p, 0, d->d->n_layers, 0, d->n_derivs, 3)) {
          fprintf(stderr, "ERROR: xrtm_fd_set_g_x_p()\n");
          return -1;
     }

     return 0;
}



double xrtm_fd_get_g_p(xrtm_fd_data *d, int i_layer, int i_deriv) {

     if (! (d->d->solvers & XRTM_SOLVERS_USE_G)) {
          fprintf(stderr, "ERROR: no two-stream solvers have been specified\n");
          return XRTM_DBL_ERROR;
     }

     CHECK_INDEX_RANGE(i_layer, d->d->n_layers, "i_layer", "n_layers", XRTM_DBL_ERROR);
     CHECK_INDEX_RANGE(i_deriv,    d->n_derivs, "i_deriv", "n_derivs", XRTM_DBL_ERROR);

     return d->g_p[i_layer][i_deriv];
}



/*******************************************************************************
 *
 ******************************************************************************/
static int xrtm_fd_set_coef_x_p(xrtm_fd_data *d, void *coef_p, int i_layer, int n_layers,
                                int i_deriv, int n_derivs, int type) {

     int i;
     int ii;
     int j;
     int jj;
     int k;
     int l;
/*
     int deps;
*/
     double a;

     if (! (d->d->solvers & XRTM_SOLVERS_USE_COEF)) {
          fprintf(stderr, "ERROR: no multi-stream solvers have been specified\n");
          return -1;
     }

     CHECK_INDEX_RANGE(i_layer, d->d->n_layers, "i_layer", "n_layers", -1);
     CHECK_INDEX_RANGE(i_deriv,    d->n_derivs, "i_deriv", "n_derivs", -1);

     for (i = i_layer, ii = 0; i < n_layers; ++i, ++ii) {
          for (j = i_deriv, jj = 0; j < n_derivs; ++j, ++jj) {
               for (k = 0; k < d->d->n_elem; ++k) {
                    for (l = 0; l < d->d->n_coef_layer[i]; ++l) {
                         switch(type) {
                              case 0:
                                   a = ((double **)   coef_p)[k][l];
                                   break;
                              case 1:
                                   a = ((double ***)  coef_p)[ii][k][l];
                                   break;
                              case 2:
                                   a = ((double ***)  coef_p)[jj][k][l];
                                   break;
                              case 3:
                                   a = ((double ****) coef_p)[ii][jj][k][l];
                                   break;
                              default:
#ifdef DEBUG
                                   fprintf(stderr, "ERROR: xrtm_fd_set_coef_x_p(): end of switch\n");
                                   exit(1);
#endif
                                   break;
                         }

                         d->coef_p[i][j][k][l] = a;
                    }
               }
          }
     }
/*
     deps = DEP_MASK_COEF_L;

     if (d->solvers & xrtm_fd_INTERNAL && d->options & xrtm_fd_DELTA_M) {
          deps |= DEP_MASK_OMEGA_L;
          deps |= DEP_MASK_LTAU_L;
     }

     set_deps(d, deps, i_layer, n_layers);
*/
     return 0;
}



int xrtm_fd_set_coef_p_11(xrtm_fd_data *d, int i_layer, int i_deriv, double **coef_p) {

     if (xrtm_fd_set_coef_x_p(d, coef_p, i_layer, i_layer + 1, i_deriv, i_deriv + 1, 0)) {
          fprintf(stderr, "ERROR: xrtm_fd_set_coef_x_p()\n");
          return -1;
     }

     return 0;
}



int xrtm_fd_set_coef_p_n1(xrtm_fd_data *d, int i_deriv, double ***coef_p) {

     if (xrtm_fd_set_coef_x_p(d, coef_p, 0, d->d->n_layers, i_deriv, i_deriv + 1, 1)) {
          fprintf(stderr, "ERROR: xrtm_fd_set_coef_x_p()\n");
          return -1;
     }

     return 0;
}



int xrtm_fd_set_coef_p_1n(xrtm_fd_data *d, int i_layer, double ***coef_p) {

     if (xrtm_fd_set_coef_x_p(d, coef_p, i_layer, i_layer + 1, 0, d->n_derivs, 2)) {
          fprintf(stderr, "ERROR: xrtm_fd_set_coef_x_p()\n");
          return -1;
     }

     return 0;
}



int xrtm_fd_set_coef_p_nn(xrtm_fd_data *d, double ****coef_p) {

     if (xrtm_fd_set_coef_x_p(d, coef_p, 0, d->d->n_layers, 0, d->n_derivs, 3)) {
          fprintf(stderr, "ERROR: xrtm_fd_set_coef_x_p()\n");
          return -1;
     }

     return 0;
}



double xrtm_fd_get_coef_p(xrtm_fd_data *d, int i_layer, int i_deriv, int i_elem, int i_coef) {

     if (! (d->d->solvers & XRTM_SOLVERS_USE_COEF)) {
          fprintf(stderr, "ERROR: no multi-stream solvers have been specified\n");
          return XRTM_DBL_ERROR;
     }

     CHECK_INDEX_RANGE(i_layer, d->d->n_layers,              "i_layer", "n_layers",              XRTM_DBL_ERROR);
     CHECK_INDEX_RANGE(i_deriv, d->n_derivs,                 "i_deriv", "n_derivs",              XRTM_DBL_ERROR);
     CHECK_INDEX_RANGE(i_elem,  6,                           "i_elem",  "6",                     XRTM_DBL_ERROR);
     CHECK_INDEX_RANGE(i_layer, d->d->n_coef_layer[i_layer], "i_coef",  "n_coef_layer[i_layer]", XRTM_DBL_ERROR);

     return d->coef_p[i_layer][i_deriv][i_elem][i_coef];
}



/*******************************************************************************
 *
 ******************************************************************************/
static int xrtm_fd_set_omega_x_p(xrtm_fd_data *d, void *omega_p, int i_layer,
                                 int n_layers, int i_deriv, int n_derivs, int type) {
/*
     int deps;
*/
     CHECK_INDEX_RANGE(i_layer, d->d->n_layers, "i_layer", "n_layers", -1);
     CHECK_INDEX_RANGE(i_deriv,    d->n_derivs, "i_deriv", "n_derivs", -1);

     if (xrtm_set_layer_x_l(d->d, d->omega_p, omega_p,
                            i_layer, n_layers, i_deriv, n_derivs, type)) {
          fprintf(stderr, "ERROR: xrtm_set_layer_x_l()\n");
          return 1;
     }
/*
     deps = DEP_MASK_OMEGA_L;

     if (d->options & xrtm_fd_DELTA_M)
          deps |= DEP_MASK_LTAU_L;

     set_deps(d, deps, i_layer, n_layers);
*/
     return 0;
}



int xrtm_fd_set_omega_p_11(xrtm_fd_data *d, int i_layer, int i_deriv, double omega_p) {

     if (xrtm_fd_set_omega_x_p(d, &omega_p, i_layer, i_layer + 1, i_deriv, i_deriv + 1, 0)) {
          fprintf(stderr, "ERROR: xrtm_fd_set_omega_x_p()\n");
          return -1;
     }

     return 0;
}



int xrtm_fd_set_omega_p_n1(xrtm_fd_data *d, int i_deriv, double *omega_p) {

     if (xrtm_fd_set_omega_x_p(d, omega_p, 0, d->d->n_layers, i_deriv, i_deriv + 1, 1)) {
          fprintf(stderr, "ERROR: xrtm_fd_set_omega_x_p()\n");
          return -1;
     }

     return 0;
}



int xrtm_fd_set_omega_p_1n(xrtm_fd_data *d, int i_layer, double *omega_p) {

     if (xrtm_fd_set_omega_x_p(d, omega_p, i_layer, i_layer + 1, 0, d->n_derivs, 2)) {
          fprintf(stderr, "ERROR: xrtm_fd_set_omega_x_p()\n");
          return -1;
     }

     return 0;
}



int xrtm_fd_set_omega_p_nn(xrtm_fd_data *d, double **omega_p) {

     if (xrtm_fd_set_omega_x_p(d, omega_p, 0, d->d->n_layers, 0, d->n_derivs, 3)) {
          fprintf(stderr, "ERROR: xrtm_fd_set_omega_x_p()\n");
          return -1;
     }

     return 0;
}



double xrtm_fd_get_omega_p(xrtm_fd_data *d, int i_layer, int i_deriv) {

     CHECK_INDEX_RANGE(i_layer, d->d->n_layers, "i_layer", "n_layers", XRTM_DBL_ERROR);
     CHECK_INDEX_RANGE(i_deriv,    d->n_derivs, "i_deriv", "n_derivs", XRTM_DBL_ERROR);

     return d->omega_p[i_layer][i_deriv];
}



/*******************************************************************************
 *
 ******************************************************************************/
static int xrtm_fd_set_ltau_x_p(xrtm_fd_data *d, void *ltau_p, int i_layer,
                                int n_layers, int i_deriv, int n_derivs, int type) {

     CHECK_INDEX_RANGE(i_layer, d->d->n_layers, "i_layer", "n_layers", -1);
     CHECK_INDEX_RANGE(i_deriv,    d->n_derivs, "i_deriv", "n_derivs", -1);

     if (xrtm_set_layer_x_l(d->d, d->ltau_p, ltau_p,
                            i_layer, n_layers, i_deriv, n_derivs, type)) {
          fprintf(stderr, "ERROR: xrtm_set_layer_x_l()\n");
          return 1;
     }
/*
     set_deps(d, DEP_MASK_LTAU_L, i_layer, n_layers);
*/
     return 0;
}



int xrtm_fd_set_ltau_p_11(xrtm_fd_data *d, int i_layer, int i_deriv, double ltau_p) {

     if (xrtm_fd_set_ltau_x_p(d, &ltau_p, i_layer, i_layer + 1, i_deriv, i_deriv + 1, 0)) {
          fprintf(stderr, "ERROR: xrtm_fd_set_ltau_x_p()\n");
          return -1;
     }

     return 0;
}



int xrtm_fd_set_ltau_p_n1(xrtm_fd_data *d, int i_deriv, double *ltau_p) {

     if (xrtm_fd_set_ltau_x_p(d, ltau_p, 0, d->d->n_layers, i_deriv, i_deriv + 1, 1)) {
          fprintf(stderr, "ERROR: xrtm_fd_set_ltau_x_p()\n");
          return -1;
     }

     return 0;
}



int xrtm_fd_set_ltau_p_1n(xrtm_fd_data *d, int i_layer, double *ltau_p) {

     if (xrtm_fd_set_ltau_x_p(d, ltau_p, i_layer, i_layer + 1, 0, d->n_derivs, 2)) {
          fprintf(stderr, "ERROR: xrtm_fd_set_ltau_x_p()\n");
          return -1;
     }

     return 0;
}



int xrtm_fd_set_ltau_p_nn(xrtm_fd_data *d, double **ltau_p) {

     if (xrtm_fd_set_ltau_x_p(d, ltau_p, 0, d->d->n_layers, 0, d->n_derivs, 3)) {
          fprintf(stderr, "ERROR: xrtm_fd_set_ltau_x_p()\n");
          return -1;
     }

     return 0;
}



double xrtm_fd_get_ltau_p(xrtm_fd_data *d, int i_layer, int i_deriv) {

     CHECK_INDEX_RANGE(i_layer, d->d->n_layers, "i_layer", "n_layers", XRTM_DBL_ERROR);
     CHECK_INDEX_RANGE(i_deriv,    d->n_derivs, "i_deriv", "n_derivs", XRTM_DBL_ERROR);

     return d->ltau_p[i_layer][i_deriv];
}



/*******************************************************************************
 *
 ******************************************************************************/
int xrtm_fd_set_kernel_ampfac_p_1(xrtm_fd_data *d, int i_kernel, int i_deriv, double ampfac_p) {

     CHECK_N_KERNELS_NOT_ZERO(d->d, -1);

     CHECK_INDEX_RANGE(i_kernel, d->d->n_kernels, "i_kernel", "n_kernel", -1);
     CHECK_INDEX_RANGE(i_deriv,     d->n_derivs,  "i_deriv",  "n_derivs", -1);

     d->kernel_ampfac_p[i_kernel][i_deriv] = ampfac_p;
/*
     set_deps(d, DEP_MASK_BRDF_L, 0, d->n_layers);
*/
     return 0;
}



int xrtm_fd_set_kernel_ampfac_p_n(xrtm_fd_data *d, int i_kernel, double *ampfac_p) {

     int i;

     CHECK_N_KERNELS_NOT_ZERO(d->d, -1);

     CHECK_INDEX_RANGE(i_kernel, d->d->n_kernels, "i_kernel", "n_kernel", -1);

     for (i = 0; i < d->n_derivs; ++i)
          d->kernel_ampfac_p[i_kernel][i] = ampfac_p[i];
/*
     set_deps(d, DEP_MASK_BRDF_L, 0, d->n_layers);
*/
     return 0;
}



double xrtm_fd_get_kernel_ampfac_p(xrtm_fd_data *d, int i_kernel, int i_deriv) {

     CHECK_N_KERNELS_NOT_ZERO(d->d, XRTM_DBL_ERROR);

     CHECK_INDEX_RANGE(i_kernel, d->d->n_kernels, "i_kernel", "n_kernel", XRTM_DBL_ERROR);
     CHECK_INDEX_RANGE(i_deriv,     d->n_derivs,  "i_deriv",  "n_derivs", XRTM_DBL_ERROR);

     return d->kernel_ampfac_p[i_kernel][i_deriv];
}



int xrtm_fd_set_kernel_params_p_11(xrtm_fd_data *d, int i_kernel, int i_deriv, int i_param, double param_p) {

     int n_params;

     CHECK_N_KERNELS_NOT_ZERO(d->d, -1);

     CHECK_INDEX_RANGE(i_kernel, d->d->n_kernels, "i_kernel", "n_kernel", -1);

     n_params = kernel_n_params((enum xrtm_kernel_type) d->d->kernels[i_kernel]);

     CHECK_INDEX_RANGE(i_deriv,  d->n_derivs,  "i_deriv",  "n_derivs", -1);
     CHECK_INDEX_RANGE(i_param,     n_params,  "i_param",  "n_params", -1);

     d->kernel_params_p[i_kernel][i_deriv][i_param] = param_p;
/*
     set_deps(d, DEP_MASK_BRDF_L, 0, d->n_layers);
*/
     return 0;
}



int xrtm_fd_set_kernel_params_p_1n(xrtm_fd_data *d, int i_kernel, int i_deriv, double *params_p) {

     int i;

     int n_params;

     CHECK_N_KERNELS_NOT_ZERO(d->d, -1);

     CHECK_INDEX_RANGE(i_kernel, d->d->n_kernels, "i_kernel", "n_kernel", -1);
     CHECK_INDEX_RANGE(i_deriv,     d->n_derivs,  "i_deriv",  "n_derivs", -1);

     n_params = kernel_n_params((enum xrtm_kernel_type) d->d->kernels[i_kernel]);

     for (i = 0; i < n_params; ++i)
          d->kernel_params_p[i_kernel][i_deriv][i] = params_p[i];
/*
     set_deps(d, DEP_MASK_BRDF_L, 0, d->n_layers);
*/
     return 0;
}


int xrtm_fd_set_kernel_params_p_n1(xrtm_fd_data *d, int i_kernel, int i_param, double *params_p) {

     int i;

     int n_params;

     CHECK_N_KERNELS_NOT_ZERO(d->d, -1);

     CHECK_INDEX_RANGE(i_kernel, d->d->n_kernels, "i_kernel", "n_kernel", -1);

     n_params = kernel_n_params((enum xrtm_kernel_type) d->d->kernels[i_kernel]);

     CHECK_INDEX_RANGE(i_param,     n_params,  "i_param",  "n_params", -1);

     for (i = 0; i < d->n_derivs; ++i)
          d->kernel_params_p[i_kernel][i][i_param] = params_p[i];
/*
     set_deps(d, DEP_MASK_BRDF_L, 0, d->n_layers);
*/
     return 0;
}



int xrtm_fd_set_kernel_params_p_nn(xrtm_fd_data *d, int i_kernel, double **params_p) {

     int i;
     int j;

     int n_params;

     CHECK_N_KERNELS_NOT_ZERO(d->d, -1);

     CHECK_INDEX_RANGE(i_kernel, d->d->n_kernels, "i_kernel", "n_kernel", -1);

     n_params = kernel_n_params((enum xrtm_kernel_type) d->d->kernels[i_kernel]);

     for (i = 0; i < d->n_derivs; ++i) {
          for (j = 0; j < n_params; ++j) {
               d->kernel_params_p[i_kernel][i][j] = params_p[i][j];
          }
     }
/*
     set_deps(d, DEP_MASK_BRDF_L, 0, d->n_layers);
*/
     return 0;
}



double xrtm_fd_get_kernel_params_p(xrtm_fd_data *d, int i_kernel, int i_deriv, int i_param) {

     int n_params;

     CHECK_N_KERNELS_NOT_ZERO(d->d, XRTM_DBL_ERROR);

     CHECK_INDEX_RANGE(i_kernel, d->d->n_kernels, "i_kernel", "n_kernel", XRTM_DBL_ERROR);

     n_params = kernel_n_params((enum xrtm_kernel_type) d->d->kernels[i_kernel]);

     CHECK_INDEX_RANGE(i_deriv,  d->n_derivs,  "i_deriv",  "n_derivs", XRTM_DBL_ERROR);
     CHECK_INDEX_RANGE(i_param,     n_params,  "i_param",  "n_params", XRTM_DBL_ERROR);

     return d->kernel_params_p[i_kernel][i_deriv][i_param];
}



/*******************************************************************************
 *
 ******************************************************************************/
static int perturb_values(xrtm_fd_data *d, int j, enum xrtm_solver_mask solver, double f, double **coef) {

     int i;
     int k;
     int l;

     int n_coef_layer;

     int n_params;

     double *levels_b;

     if (d->d->options & XRTM_OPTION_SOURCE_THERMAL) {
           if (xrtm_set_F_iso_top(d->d, xrtm_get_F_iso_top(d->d) + f * d->F_iso_top_p[j])) {
                fprintf(stderr, "ERROR: xrtm_set_F_iso_top()\n");
                return -1;
           }

           if (xrtm_set_F_iso_bot(d->d, xrtm_get_F_iso_bot(d->d) + f * d->F_iso_bot_p[j])) {
                fprintf(stderr, "ERROR: xrtm_set_F_iso_bot()\n");
                return -1;
           }

           levels_b = alloc_array1_d(d->d->n_layers + 1);

           if (xrtm_get_levels_b(d->d, levels_b)) {
                fprintf(stderr, "ERROR: xrtm_get_levels_b()\n");
                return -1;
           }

           for (i = 0; i < d->d->n_layers + 1; ++i) {
                levels_b[i] += f * d->levels_b_p[i][j];
           }

           if (xrtm_set_levels_b(d->d, levels_b)) {
                fprintf(stderr, "ERROR: xrtm_set_levels_b()\n");
                return -1;
           }

           free_array1_d(levels_b);


           if (xrtm_set_surface_b(d->d, xrtm_get_surface_b(d->d) + f * d->surface_b_p[j])) {
                fprintf(stderr, "ERROR: xrtm_set_surface_b()\n");
                return -1;
           }
     }

     if (solver & XRTM_SOLVERS_USE_G) {
           for (i = 0; i < d->d->n_layers; ++i) {
                if (xrtm_set_g_1(d->d, i, xrtm_get_g(d->d, i) + f * d->g_p[i][j])) {
                     fprintf(stderr, "ERROR: xrtm_set_g_1()\n");
                     return -1;
                }
           }
      }
      else {
           for (i = 0; i < d->d->n_layers; ++i) {
                n_coef_layer = xrtm_get_n_coef(d->d, i);

                for (k = 0; k < d->d->n_elem; ++k) {
                     for (l = 0; l < n_coef_layer; ++l) {
                          coef[k][l] = xrtm_get_coef(d->d, i, k, l) + f * d->coef_p[i][j][k][l];
                     }
                }

                if (xrtm_set_coef_1(d->d, i, n_coef_layer, coef)) {
                     fprintf(stderr, "ERROR: xrtm_set_coef_1()\n");
                     return -1;
                }
           }
      }

      for (i = 0; i < d->d->n_layers; ++i) {
           if (xrtm_set_omega_1(d->d, i, xrtm_get_omega(d->d, i) + f * d->omega_p[i][j])) {
                fprintf(stderr, "ERROR: xrtm_set_omega_1()\n");
                return -1;
           }
      }

      for (i = 0; i < d->d->n_layers; ++i) {
           if (xrtm_set_ltau_1(d->d,  i, xrtm_get_ltau(d->d, i)  + f * d->ltau_p [i][j])) {
                fprintf(stderr, "ERROR: xrtm_set_ltau_1()\n");
                return -1;
           }
     }

     for (i = 0; i < d->d->n_kernels; ++i) {
          if (xrtm_set_kernel_ampfac(d->d, i, xrtm_get_kernel_ampfac(d->d, i) + f * d->kernel_ampfac_p[i][j])) {
                fprintf(stderr, "ERROR: xrtm_set_ampfac()\n");
                return -1;
          }

          n_params = kernel_n_params((enum xrtm_kernel_type) d->d->kernels[i]);

          for (k = 0; k < n_params; ++k) {
                if (xrtm_set_kernel_params_1(d->d, i, k, xrtm_get_kernel_params(d->d, i, k) + f * d->kernel_params_p[i][j][k])) {
                      fprintf(stderr, "ERROR: xrtm_set_kernel_param_1()\n");
                      return -1;
                }
          }
     }

     return 0;
}



/*******************************************************************************
 *
 ******************************************************************************/
int xrtm_fd_solution(xrtm_fd_data *d, enum xrtm_solver_mask solver, int solutions, enum xrtm_fd_method_type fd_method, int n_out_phis, double **out_phis, double ****I_p, double ****I_m, double *****I_p_l, double *****I_m_l, double *mean_p, double *mean_m, double **mean_p_l, double **mean_m_l, double *flux_p, double *flux_m, double **flux_p_l, double **flux_m_l, double *flux_div, double **flux_div_l) {

     int i;
     int j;
     int k;
     int l;
     int m;
/*
     int n_quad_x;
*/
     int n_mus2;

     double f;

     double **coef;

     double ****I_p2;
     double ****I_m2;


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
/*
     n_quad_x = d->d->n_quad + d->d->n_umus;
*/
     if (d->d->n_umus == 0)
          n_mus2  = d->d->n_quad;
     else
          n_mus2  = d->d->n_umus;
/*
     else
          n_mus2  = n_quad_x;
*/


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (xrtm_radiance(d->d, solver, n_out_phis, out_phis, I_p, I_m, NULL, NULL)) {
          fprintf(stderr, "ERROR: xrtm_radiance()\n");
          return -1;
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     coef = alloc_array2_d(d->d->n_elem, d->d->n_coef);

     I_p2 = alloc_array4_d(d->d->n_ulevels, n_mus2, n_out_phis, d->d->n_stokes);
     I_m2 = alloc_array4_d(d->d->n_ulevels, n_mus2, n_out_phis, d->d->n_stokes);
/*
     coef = get_work2_d(&work, d->d->n_elem, d->d->n_coef);

     I_p2 = get_work4_d(&work, d->d->n_ulevels, n_out_phis, n_mus2, d->d->n_stokes);
     I_m2 = get_work4_d(&work, d->d->n_ulevels, n_out_phis, n_mus2, d->d->n_stokes);
*/

     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     for (j = 0; j < d->n_derivs; ++j) {
          if (fd_method == XRTM_FD_METHOD_FORWARD || fd_method == XRTM_FD_METHOD_BACKWARD) {
               if (fd_method == XRTM_FD_METHOD_FORWARD)
                    f =  1.;
               else
               if (fd_method == XRTM_FD_METHOD_BACKWARD)
                    f = -1.;

               if (perturb_values(d, j, solver, f, coef)) {
                    fprintf(stderr, "ERROR: perturb_values()\n");
                    return -1;
               }

               if (xrtm_solution(d->d, solver, solutions, n_out_phis, out_phis, I_p2, I_m2, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL)) {
                    fprintf(stderr, "ERROR: xrtm_solution()\n");
                    return -1;
               }

               for (i = 0; i < d->d->n_ulevels; ++i) {
                    for (k = 0; k < n_mus2; ++k) {
                         for (l = 0; l < n_out_phis; ++l) {
                              for (m = 0; m < d->d->n_stokes; ++m) {
                                   I_p_l[i][j][k][l][m] = (I_p2[i][k][l][m] - I_p[i][k][l][m]) / d->delta[j];
                                   I_m_l[i][j][k][l][m] = (I_m2[i][k][l][m] - I_m[i][k][l][m]) / d->delta[j];
                              }
                         }
                    }
               }
          }


          /*--------------------------------------------------------------------
           *
           *------------------------------------------------------------------*/
          else {
               if (perturb_values(d, j, solver,  1., coef)) {
                    fprintf(stderr, "ERROR: perturb_values()\n");
                    return -1;
               }

               if (xrtm_solution(d->d, solver, solutions, n_out_phis, out_phis, I_p2, I_m2, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL)) {
                    fprintf(stderr, "ERROR: xrtm_solution()\n");
                    return -1;
               }

               for (i = 0; i < d->d->n_ulevels; ++i) {
                    for (k = 0; k < n_mus2; ++k) {
                         for (l = 0; l < n_out_phis; ++l) {
                              for (m = 0; m < d->d->n_stokes; ++m) {
                                   I_p_l[i][j][k][l][m]  = I_p2[i][k][l][m] / d->delta[j];
                                   I_m_l[i][j][k][l][m]  = I_m2[i][k][l][m] / d->delta[j];
                              }
                         }
                    }
               }

               if (perturb_values(d, j, solver, -2., coef)) {
                    fprintf(stderr, "ERROR: perturb_values()\n");
                    return -1;
               }

               if (xrtm_solution(d->d, solver, solutions, n_out_phis, out_phis, I_p2, I_m2, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL)) {
                    fprintf(stderr, "ERROR: xrtm_solution()\n");
                    return -1;
               }

               for (i = 0; i < d->d->n_ulevels; ++i) {
                    for (k = 0; k < n_mus2; ++k) {
                         for (l = 0; l < n_out_phis; ++l) {
                              for (m = 0; m < d->d->n_stokes; ++m) {
                                   I_p_l[i][j][k][l][m] -= I_p2[i][k][l][m] / d->delta[j];
                                   I_m_l[i][j][k][l][m] -= I_m2[i][k][l][m] / d->delta[j];
                              }
                         }
                    }
               }

               if (perturb_values(d, j, solver,  1., coef)) {
                    fprintf(stderr, "ERROR: perturb_values()\n");
                    return -1;
               }
          }
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/

     free_array2_d(coef);

     free_array4_d(I_p2);
     free_array4_d(I_m2);


     return 0;
}



/*******************************************************************************
 *
 ******************************************************************************/
int xrtm_fd_radiance(xrtm_fd_data *d, enum xrtm_solver_mask solver, enum xrtm_fd_method_type fd_method, int n_out_phis, double **out_phis, double ****I_p, double ****I_m, double *****I_p_l, double *****I_m_l) {

     if (xrtm_fd_solution(d, solver, XRTM_OUTPUT_RADIANCE, fd_method, n_out_phis, out_phis, I_p, I_m, I_p_l, I_m_l, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL)) {
          fprintf(stderr, "ERROR: xrtm_fd_solution()\n");
          return -1;
     }

     return 0;
}
