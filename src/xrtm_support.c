/******************************************************************************%
**
**    Copyright (C) 2007-2012 Greg McGarragh <gregm@atmos.colostate.edu>
**
**    This source code is licensed under the GNU General Public License (GPL),
**    Version 3.  See the file COPYING for more details.
**
*******************************************************************************/

#include <gutil.h>

#include <rtutil_support.h>

#include "xrtm_support.h"


/*******************************************************************************
 *
 ******************************************************************************/
static const char *kernel_names[] = {
     "lambertian",
     "roujean",
     "li_sparse",
     "li_dense",
     "ross_thin",
     "ross_thick",
     "hapke",
     "rahman",
     "cox_munk",
};

int kernel_code(char *name) {
     return name_to_code(name, kernel_names,
                         N_XRTM_KERNEL_TYPES, "kernel type");
}

const char *kernel_name(enum xrtm_kernel_type code) {
     return code_to_name(code, kernel_names,
                         N_XRTM_KERNEL_TYPES, "kernel type");
}



static const char *eigen_solver_gen_real_names[] = {
     "asymtx",
     "eispack",
     "lapack"
};

int eigen_solver_gen_real_code(char *name) {
     return name_to_code(name, eigen_solver_gen_real_names,
          N_EIGEN_SOLVER_GEN_REAL_TYPES, "general real eigen solver types");
}

const char *eigen_solver_gen_real_name(enum eigen_solver_gen_real_type code) {
     return code_to_name(code, eigen_solver_gen_real_names,
          N_EIGEN_SOLVER_GEN_REAL_TYPES, "general real eigen solver types");
}



static const char *eigen_solver_gen_complex_names[] = {
     "eispack",
     "lapack"
};

int eigen_solver_gen_complex_code(char *name) {
     return name_to_code(name, eigen_solver_gen_complex_names,
          N_EIGEN_SOLVER_GEN_COMPLEX_TYPES, "general complex eigen solver types");
}

const char *eigen_solver_gen_complex_name(enum eigen_solver_gen_complex_type code) {
     return code_to_name(code, eigen_solver_gen_complex_names,
          N_EIGEN_SOLVER_GEN_COMPLEX_TYPES, "general complex eigen solver types");
}



/*******************************************************************************
 *
 ******************************************************************************/
void *jalloc_array1(jmp_buf env, long m, int size) {

     void *x;

     if (! (x = alloc_array1(m, size))) {
          eprintf("ERROR: alloc_array1_d(%ld)\n", m);
          longjmp(env, -1);
     }

     return x;
}



void *jalloc_array2(jmp_buf env, long m, long n, int size) {

     void *x;

     if (! (x = alloc_array2(m, n, size))) {
          eprintf("ERROR: alloc_array2_d(%ld, %ld)\n", m, n);
          longjmp(env, -1);
     }

     return x;
}



void *jalloc_array3(jmp_buf env, long m, long n, long o, int size) {

     void *x;

     if (! (x = alloc_array3(m, n, o, size))) {
          eprintf("ERROR: alloc_array3_d(%ld, %ld, %ld)\n", m, n, o);
          longjmp(env, -1);
     }

     return x;
}



void *jalloc_array4(jmp_buf env, long m, long n, long o, long p, int size) {

     void *x;

     if (! (x = alloc_array4(m, n, o, p, size))) {
          eprintf("ERROR: alloc_array4_d(%ld, %ld, %ld, %ld)\n", m, n, o, p);
          longjmp(env, -1);
     }

     return x;
}



void *jalloc_array5(jmp_buf env, long m, long n, long o, long p, long q, int size) {

     void *x;

     if (! (x = alloc_array5(m, n, o, p, q, size))) {
          eprintf("ERROR: alloc_array5_d(%ld, %ld, %ld, %ld, %d)\n", m, n, o, p, q);
          longjmp(env, -1);
     }

     return x;
}



/*******************************************************************************
 *
 ******************************************************************************/
int xrtm_set_layer_x(xrtm_data *d, double *x, double *x_, int i_layer, int n_layers, uchar *set_flags) {

     int i;

     for (i = i_layer; i < n_layers; ++i) {
          x[i] = x_[i - i_layer];

          if (d->initial_inputs)
               set_flags[i] = 1;
     }

     return 0;
}



int xrtm_set_layer_x_l(xrtm_data *d, double **x, void *x_, int i_layer, int n_layers, int i_deriv, int n_derivs, int type) {

     int i;
     int ii;
     int j;
     int jj;

     double a;

     for (i = i_layer, ii = 0; i < n_layers; ++i, ++ii) {
          for (j = i_deriv, jj = 0; j < n_derivs; ++j, ++jj) {
               switch(type) {
                    case 0:
                         a = ((double *)  x_)[ 0];
                         break;
                    case 1:
                         a = ((double *)  x_)[ii];
                         break;
                    case 2:
                         a = ((double *)  x_)[jj];
                         break;
                    case 3:
                         a = ((double **) x_)[ii][jj];
                         break;
                    default:
#ifdef DEBUG
                         eprintf("ERROR: xrtm_set_layer_x_l(): end of switch\n");
                         exit(1);
#endif
                         break;
               }

               x[i][j] = a;
          }
     }

     return 0;
}



/*******************************************************************************
 *
 ******************************************************************************/
#ifdef INCLUDE_DEV_SOURCE
#include "../src2/xrtm_support2.c"
#endif
