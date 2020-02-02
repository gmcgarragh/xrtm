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

#include "xrtm_support.h"


/*******************************************************************************
 *
 ******************************************************************************/
static const char *eigen_solver_gen_real_names[] = {
     "asymtx",
     "eispack",
     "lapack"
};

static long eigen_solver_gen_real_types[] = {
#ifdef HAVE_FORTRAN_COMPILER
     EIGEN_SOLVER_GEN_REAL_ASYMTX,
#endif
#ifdef HAVE_EISPACK_LIBRARY
     EIGEN_SOLVER_GEN_REAL_EISPACK,
#endif
     EIGEN_SOLVER_GEN_REAL_LAPACK,

     N_EIGEN_SOLVER_GEN_REAL_TYPES
};


GINDEX_NAME_VALUE_TEMPLATE(eigen_solver_gen_real, "eigen solver general real", N_EIGEN_SOLVER_GEN_REAL_TYPES)



static const char *eigen_solver_gen_complex_names[] = {
     "eispack",
     "lapack"
};

static long eigen_solver_gen_complex_types[] = {
#ifdef HAVE_EISPACK_LIBRARY
     EIGEN_SOLVER_GEN_COMPLEX_EISPACK,
#endif
     EIGEN_SOLVER_GEN_COMPLEX_LAPACK,

     N_EIGEN_SOLVER_GEN_COMPLEX_TYPES
};


GINDEX_NAME_VALUE_TEMPLATE(eigen_solver_gen_complex, "eigen solver general complex", N_EIGEN_SOLVER_GEN_COMPLEX_TYPES)




/*******************************************************************************
 *
 ******************************************************************************/
void *jalloc_array1(jmp_buf env, long m, int size) {

     void *x;

     if (! (x = alloc_array1(m, size))) {
          fprintf(stderr, "ERROR: alloc_array1_d(%ld)\n", m);
          longjmp(env, -1);
     }

     return x;
}



void *jalloc_array2(jmp_buf env, long m, long n, int size) {

     void *x;

     if (! (x = alloc_array2(m, n, size))) {
          fprintf(stderr, "ERROR: alloc_array2_d(%ld, %ld)\n", m, n);
          longjmp(env, -1);
     }

     return x;
}



void *jalloc_array3(jmp_buf env, long m, long n, long o, int size) {

     void *x;

     if (! (x = alloc_array3(m, n, o, size))) {
          fprintf(stderr, "ERROR: alloc_array3_d(%ld, %ld, %ld)\n", m, n, o);
          longjmp(env, -1);
     }

     return x;
}



void *jalloc_array4(jmp_buf env, long m, long n, long o, long p, int size) {

     void *x;

     if (! (x = alloc_array4(m, n, o, p, size))) {
          fprintf(stderr, "ERROR: alloc_array4_d(%ld, %ld, %ld, %ld)\n", m, n, o, p);
          longjmp(env, -1);
     }

     return x;
}



void *jalloc_array5(jmp_buf env, long m, long n, long o, long p, long q, int size) {

     void *x;

     if (! (x = alloc_array5(m, n, o, p, q, size))) {
          fprintf(stderr, "ERROR: alloc_array5_d(%ld, %ld, %ld, %ld, %ld)\n", m, n, o, p, q);
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

          if (d->inputs_initialized)
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
                         fprintf(stderr, "ERROR: xrtm_set_layer_x_l(): end of switch\n");
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
