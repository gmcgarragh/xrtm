/******************************************************************************%
**
**    Copyright (C) 2007-2012 Greg McGarragh <gregm@atmos.colostate.edu>
**
**    This source code is licensed under the GNU General Public License (GPL),
**    Version 3.  See the file COPYING for more details.
**
*******************************************************************************/

#include <gutil.h>

#include "xrtm.h"
#include "xrtm_work.h"
#include "xrtm_utility.h"


/*******************************************************************************
 *
 ******************************************************************************/
int work_init(work_data *work, int n_quad, int n_quad_x, int n_stokes, int n_derivs, int n_layers, int n_umus) {

     int i;

     work_shared *w;

     work->i_page = 0;
     work->i_byte = 0;

     for (i = 0; i < N_WORK_TYPES  ; ++i)
          work->i_list  [i] = 0;

     for (i = 0; i < N_WORK_TYPES_V; ++i)
          work->i_list_v[i] = 0;

     work->w = (work_shared *) malloc(sizeof(work_shared));

     if (! work->w) {
          eprintf("ERROR: memory allocation failed\n");
          return -1;
     }

     w = work->w;

     w->n_quad_v   =  n_quad   * n_stokes;
     w->n_quad_v_x =  n_quad_x * n_stokes;
     w->n_stokes   =  n_stokes;
     w->n_derivs   =  n_derivs;
     w->n_layers   =  n_layers + 1;
     w->n_umus_v   =  n_umus * n_stokes;

     for (i = 0; i < WORK_MAX_PAGES; ++i) {
          w->n_bytes[i] = 0;
          w->pages  [i] = NULL;
     }

     w->min_page_size = (n_derivs > 0 ? n_derivs : 1) * w->n_quad_v * w->n_quad_v * sizeof(double);

     for (i = 0; i < N_WORK_TYPES  ; ++i) {
          w->n_list    [i] = 0;
          w->max_list  [i] = 0;
          w->lists     [i] = NULL;
     }

     for (i = 0; i < N_WORK_TYPES_V; ++i) {
          w->n_list_v  [i] = 0;
          w->max_list_v[i] = 0;
          w->lists_v   [i] = NULL;
     }

     return 0;
}



/*******************************************************************************
 *
 ******************************************************************************/
static void free_work(work_shared *w, int i_type, void (*free_array)(void *)) {

     int i;

     for (i = 0; i < w->n_list[i_type]; ++i)
          free_array((void *) w->lists[i_type][i]);

     free((void *) w->lists[i_type]);
}



static void free_work_v(work_shared *w, int i_type, void (*free_array)(void *)) {

     int i;

     for (i = 0; i < w->n_list_v[i_type]; ++i)
          free_array((void *) w->lists_v[i_type][i]);

     free((void *) w->lists_v[i_type]);
}



void work_free(work_data *work) {

     int i;

     work_shared *w;

     w = work->w;

     for (i = 0; i < WORK_MAX_PAGES; ++i) {
          if (w->n_bytes[i] > 0)
               free(w->pages[i]);
     }

     free_work(w, WORK_IX,      (void (*)(void *)) free_array1_i);
     free_work(w, WORK_IX2,     (void (*)(void *)) free_array1_i);

     free_work(w, WORK_DX,      (void (*)(void *)) free_array1_d);
     free_work(w, WORK_DX2,     (void (*)(void *)) free_array1_d);

     free_work(w, WORK_DXX,     (void (*)(void *)) free_array2_d);
     free_work(w, WORK_DXX2,    (void (*)(void *)) free_array2_d);

     free_work(w, WORK_DU,      (void (*)(void *)) free_array1_d);
     free_work(w, WORK_DUX,     (void (*)(void *)) free_array2_d);

     free_work(w, WORK_DD,      (void (*)(void *)) free_array1_d);
     free_work(w, WORK_DS,      (void (*)(void *)) free_array1_d);

     free_work(w, WORK_DDERIVS, (void (*)(void *)) free_array1_d);
     free_work(w, WORK_DLAYERS, (void (*)(void *)) free_array1_d);
     free_work(w, WORK_DBOTH,   (void (*)(void *)) free_array2_d);

     free_work(w, WORK_ZX,      (void (*)(void *)) free_array1_d);
     free_work(w, WORK_ZXX,     (void (*)(void *)) free_array2_d);

     free_work(w, WORK_ZU,      (void (*)(void *)) free_array1_d);
     free_work(w, WORK_ZUX,     (void (*)(void *)) free_array2_d);

     free_work_v(w, WORK_DERIVS_V, (void (*)(void *)) free_array1);
     free_work_v(w, WORK_LAYERS_V, (void (*)(void *)) free_array1);
     free_work_v(w, WORK_BOTH_V,   (void (*)(void *)) free_array2);

     free(w);
}



/*******************************************************************************
 *
 ******************************************************************************/
static size_t get_work_common(work_data *work, size_t length) {

     size_t i;

     length = length % 16 ? length / WORK_ALIGNMENT * WORK_ALIGNMENT + WORK_ALIGNMENT : length;

     if (work->i_byte == 0 && work->w->n_bytes[work->i_page] == 0) {
          work->w->n_bytes[work->i_page] = MAX(work->w->min_page_size, length);

          work->w->pages  [work->i_page] = (uchar *) malloc(work->w->n_bytes[work->i_page]);

          i = 0;

          work->i_byte = length;
     }
     else
     if (work->i_byte == 0 && length > work->w->n_bytes[work->i_page]) {
          work->w->n_bytes[work->i_page] = MAX(work->w->min_page_size, length);

          free(work->w->pages[work->i_page]);
          work->w->pages  [work->i_page] = (uchar *) malloc(work->w->n_bytes[work->i_page]);

          i = 0;

          work->i_byte = length;
     }
     else
     if (work->i_byte + length > work->w->n_bytes[work->i_page]) {
          work->i_page++;

          if (work->i_page >= WORK_MAX_PAGES) {
               eprintf("ERROR: number of required pages greater than WORK_MAX_PAGES\n");
               exit(1);
          }

          if (length > work->w->n_bytes[work->i_page]) {
               work->w->n_bytes[work->i_page] = MAX(work->w->min_page_size, length);

               free(work->w->pages[work->i_page]);
               work->w->pages  [work->i_page] = (uchar *) malloc(work->w->n_bytes[work->i_page]);
          }

          i = 0;

          work->i_byte = length;
     }
     else {
          i = work->i_byte;

          work->i_byte += length;
     }

     return i;
}


/*******************************************************************************
 *
 ******************************************************************************/
void *get_work_x1(work_data *work, size_t m, size_t size) {

     size_t i_byte;
     size_t length;

     length = array_mem_length1(m, size);

     i_byte = get_work_common(work, length);

     return (void *) array_from_mem1(work->w->pages[work->i_page] + i_byte, m, size, 0);
}



uchar *get_work_uc1(work_data *work, size_t m) {

     return (uchar *) get_work_x1(work, m, sizeof(uchar));
}



int *get_work_i1(work_data *work, size_t m) {

     return (int *) get_work_x1(work, m, sizeof(int));
}



double *get_work_d1(work_data *work, size_t m) {

     return (double *) get_work_x1(work, m, sizeof(double));
}



dcomplex *get_work_dc1(work_data *work, size_t m) {

     return (dcomplex *) get_work_x1(work, m, sizeof(dcomplex));
}



/*******************************************************************************
 *
 ******************************************************************************/
void **get_work_x2(work_data *work, size_t m, size_t n, size_t size) {

     size_t i_byte;
     size_t length;

     length = array_mem_length2(m, n, size);

     i_byte = get_work_common(work, length);

     return (void **) array_from_mem2(work->w->pages[work->i_page] + i_byte, m, n, size, 0);
}



int **get_work_i2(work_data *work, size_t m, size_t n) {

     return (int **) get_work_x2(work, m, n, sizeof(int));
}



double **get_work_d2(work_data *work, size_t m, size_t n) {

     return (double **) get_work_x2(work, m, n, sizeof(double));
}



dcomplex **get_work_dc2(work_data *work, size_t m, size_t n) {

     return (dcomplex **) get_work_x2(work, m, n, sizeof(dcomplex));
}



/*******************************************************************************
 *
 ******************************************************************************/
void ***get_work_x3(work_data *work, size_t m, size_t n, size_t o, size_t size) {

     size_t i_byte;
     size_t length;

     length = array_mem_length3(m, n, o, size);

     i_byte = get_work_common(work, length);

     return (void ***) array_from_mem3(work->w->pages[work->i_page] + i_byte, m, n, o, size, 0);
}



int ***get_work_i3(work_data *work, size_t m, size_t n, size_t o) {

     return (int ***) get_work_x3(work, m, n, o, sizeof(int));
}



double ***get_work_d3(work_data *work, size_t m, size_t n, size_t o) {

     return (double ***) get_work_x3(work, m, n, o, sizeof(double));
}



dcomplex ***get_work_dc3(work_data *work, size_t m, size_t n, size_t o) {

     return (dcomplex ***) get_work_x3(work, m, n, o, sizeof(dcomplex));
}



/*******************************************************************************
 *
 ******************************************************************************/
void ****get_work_x4(work_data *work, size_t m, size_t n, size_t o, size_t p, size_t size) {

     size_t i_byte;
     size_t length;

     length = array_mem_length4(m, n, o, p, size);

     i_byte = get_work_common(work, length);

     return (void ****) array_from_mem4(work->w->pages[work->i_page] + i_byte, m, n, o, p, size, 0);
}



int ****get_work_i4(work_data *work, size_t m, size_t n, size_t o, size_t p) {

     return (int ****) get_work_x4(work, m, n, o, p, sizeof(int));
}



double ****get_work_d4(work_data *work, size_t m, size_t n, size_t o, size_t p) {

     return (double ****) get_work_x4(work, m, n, o, p, sizeof(double));
}



dcomplex ****get_work_dc4(work_data *work, size_t m, size_t n, size_t o, size_t p) {

     return (dcomplex ****) get_work_x4(work, m, n, o, p, sizeof(dcomplex));
}



/*******************************************************************************
 *
 ******************************************************************************/
static void check_work_bounds(int i, int *n, void **data) {

     if (i >= *n) {
          *n += ((i - *n + 1) / WORK_LIST_INC + 1) * WORK_LIST_INC;
          *data = realloc(*data, *n * sizeof(void *));
     }
}



/*******************************************************************************
 *
 ******************************************************************************/
#define GET_WORK1_VECTOR(i, n, max, alloc_vector, mm, data) do {		\
										\
     if (i >= n)								\
          get_work_vector(&i, &n, &max, alloc_vector, mm, (void ***) &data);	\
										\
     return data[i++];								\
} while (0)



static void get_work_vector(WORK_I_TYPE *i, int *n, int *max, void **(*alloc_vector_x)(long m), int mm, void ***data) {

     check_work_bounds(*n, max, (void **) data);

     (*data)[*i] = alloc_vector_x(mm);
     if (! (*data)[*i]) {
          eprintf("ERROR: memory allocation failed\n");
          exit(1);
     }

     (*n)++;
}



static void **get_work_vectors(WORK_I_TYPE *iq, int *n1, int *max1, WORK_I_TYPE *iq2, int *n2, int *max2, void **(*alloc_vector_x)(long m), int ll, int mm, void ****data1, void ***data2, uchar *flags) {

     int i;

     int count;

     if (*iq >= *n1) {

          check_work_bounds(*n1, max1, (void **) data1);

          (*data1)[*iq] = (void **) alloc_array1(ll, sizeof(void **));
          if (! (*data1)[*iq]) {
               eprintf("ERROR: memory allocation failed\n");
               exit(1);
          }

          (*n1)++;
     }

     count = flags ? flags_count(flags, ll) : ll;

     if (*iq2 + count - 1 >= *n2) {
          check_work_bounds(*iq2 + count - 1, max2, (void **) data2);

          for (i = *n2; i < *iq2 + count; ++i) {
               (*data2)[i] = alloc_vector_x(mm);
               if (! (*data2)[i]) {
                    eprintf("ERROR: memory allocation failed\n");
                    exit(1);
               }
          }

          *n2 = *iq2 + count;
     }

     if (! flags) {
          for (i = 0; i < ll; ++i) {
               (*data1)[*iq][i] = (*data2)[(*iq2)++];
          }
     }
     else {
          for (i = 0; i < ll; ++i) {
               if (flags[i])
                    (*data1)[*iq][i] = (*data2)[(*iq2)++];
#ifdef DEBUG
               else
                    (*data1)[*iq][i] = NULL;
#endif
          }
     }

     return (*data1)[(*iq)++];
}



static void ***get_work_vectors2(WORK_I_TYPE *iq, int *n1, int *max1, WORK_I_TYPE *iq2, int *n2, int *max2, void **(*alloc_vector_x)(long m), int ll1, int ll2, int mm, void *****data1, void ***data2, uchar **flags) {

     int i;
     int j;

     int count;

     if (*iq >= *n1) {
          check_work_bounds(*n1, max1, (void **) data1);

          (*data1)[*iq] = (void ***) alloc_array2(ll1, ll2, sizeof(void ***));
          if (! (*data1)[*iq]) {
               eprintf("ERROR: memory allocation failed\n");
               exit(1);
          }

          (*n1)++;
     }

     if (! flags)
          count = ll1 * ll2;
     else {
          count = 0;
          for (i = 0; i < ll1; ++i)
               count += flags_count(flags[i], ll2);
     }

     if (*iq2 + count - 1 >= *n2) {
          check_work_bounds(*iq2 + count - 1, max2, (void **) data2);

          for (i = *n2; i < *iq2 + count; ++i) {
               (*data2)[i] = alloc_vector_x(mm);
               if (! (*data2)[i]) {
                    eprintf("ERROR: memory allocation failed\n");
                    exit(1);
               }
          }

          *n2 = *iq2 + count;
     }

     if (! flags) {
          for (i = 0; i < ll1; ++i) {
               for (j = 0; j < ll2; ++j) {
                    (*data1)[*iq][i][j] = (*data2)[(*iq2)++];
               }
          }
     }
     else {
          for (i = 0; i < ll1; ++i) {
               for (j = 0; j < ll2; ++j) {
                    if (flags[i][j])
                         (*data1)[*iq][i][j] = (*data2)[(*iq2)++];
#ifdef DEBUG
                    else
                         (*data1)[*iq][i][j] = NULL;
#endif
               }
          }
     }

     return (*data1)[(*iq)++];
}



/*******************************************************************************
 *
 ******************************************************************************/
#define GET_WORK1_MATRIX(i, n, max, alloc_matrix, mm, nn, data) do {			\
											\
     if (i >= n)									\
           get_work_matrix(&i, &n, &max, alloc_matrix, mm, nn, (void ****) &data);	\
											\
     return data[i++];									\
} while (0)



static void get_work_matrix(WORK_I_TYPE *i, int *n, int *max, void ** (*alloc_matrix_x)(long m, long n), int mm, int nn, void ****data) {

     check_work_bounds(*n, max, (void **) data);

     (*data)[*i] = alloc_matrix_x(mm, nn);
     if (! (*data)[*i]) {
          eprintf("ERROR: memory allocation failed\n");
          exit(1);
     }

     (*n)++;
}



static void **get_work_matrices(WORK_I_TYPE *iq, int *n1, int *max1, WORK_I_TYPE *iq2, int *n2, int *max2, void **(*alloc_matrix_x)(long m, long n), int ll, int mm, int nn, void ****data1, void ****data2, uchar *flags) {

     int i;

     int count;

     if (*iq >= *n1) {
          check_work_bounds(*n1, max1, (void **) data1);

          (*data1)[*iq] = (void **) alloc_array1(ll, sizeof(void **));
          if (! (*data1)[*iq]) {
               eprintf("ERROR: memory allocation failed\n");
               exit(1);
          }

          (*n1)++;
     }

     count = flags ? flags_count(flags, ll) : ll;

     if (*iq2 + count - 1 >= *n2) {
          check_work_bounds(*iq2 + count - 1, max2, (void **) data2);

          for (i = *n2; i < *iq2 + count; ++i) {
               (*data2)[i] = alloc_matrix_x(mm, nn);
               if (! (*data2)[i]) {
                    eprintf("ERROR: memory allocation failed\n");
                    exit(1);
               }
          }

          *n2 = *iq2 + count;
     }

     if (! flags) {
          for (i = 0; i < ll; ++i) {
               (*data1)[*iq][i] = (*data2)[(*iq2)++];
          }
     }
     else {
          for (i = 0; i < ll; ++i) {
               if (flags[i])
                    (*data1)[*iq][i] = (*data2)[(*iq2)++];
#ifdef DEBUG
               else
                    (*data1)[*iq][i] = NULL;
#endif
          }
     }

     return (*data1)[(*iq)++];
}



static void ***get_work_matrices2(WORK_I_TYPE *iq, int *n1, int *max1, WORK_I_TYPE *iq2, int *n2, int *max2, void **(*alloc_matrix_x)(long m, long n), int ll1, int ll2, int mm, int nn, void *****data1, void ****data2, uchar **flags) {

     int i;
     int j;

     int count;

     if (*iq >= *n1) {
          check_work_bounds(*n1, max1, (void **) data1);

          (*data1)[*iq] = (void ***) alloc_array2(ll1, ll2, sizeof(void **));
          if (! (*data1)[*iq]) {
               eprintf("ERROR: memory allocation failed\n");
               exit(1);
          }

          (*n1)++;
     }

     if (! flags)
          count = ll1 * ll2;
     else {
          count = 0;
          for (i = 0; i < ll1; ++i)
               count += flags_count(flags[i], ll2);
     }

     if (*iq2 + count - 1 >= *n2) {
          check_work_bounds(*iq2 + count - 1, max2, (void **) data2);

          for (i = *n2; i < *iq2 + count; ++i) {
               (*data2)[i] = alloc_matrix_x(mm, nn);
               if (! (*data2)[i]) {
                    eprintf("ERROR: memory allocation failed\n");
                    exit(1);
               }
          }

          *n2 = *iq2 + count;
     }

     if (! flags) {
          for (i = 0; i < ll1; ++i) {
               for (j = 0; j < ll2; ++j) {
                    (*data1)[*iq][i][j] = (*data2)[(*iq2)++];
               }
          }
     }
     else {
          for (i = 0; i < ll1; ++i) {
               for (j = 0; j < ll2; ++j) {
                    if (! flags || flags[i][j])
                         (*data1)[*iq][i][j] = (*data2)[(*iq2)++];
#ifdef DEBUG
                    else
                         (*data1)[*iq][i][j] = NULL;
#endif
               }
          }
     }

     return (*data1)[(*iq)++];
}



/*******************************************************************************
 *
 ******************************************************************************/
void *work_get1(work_data *work, enum work_type type) {

     work_shared *w;

     w = work->w;

     switch(type) {
          case WORK_IX:
               GET_WORK1_VECTOR(work->i_list[WORK_IX], w->n_list[WORK_IX], w->max_list[WORK_IX], (void **(*)(long)) alloc_array1_i, w->n_quad_v_x, w->lists[WORK_IX]);
          case WORK_IX2:
               GET_WORK1_VECTOR(work->i_list[WORK_IX2], w->n_list[WORK_IX2], w->max_list[WORK_IX2], (void **(*)(long)) alloc_array1_i, 2 * w->n_quad_v_x, w->lists[WORK_IX2]);
          case WORK_DQQ:
               GET_WORK1_MATRIX(work->i_list[WORK_DQQ], w->n_list[WORK_DQQ], w->max_list[WORK_DQQ], (void **(*)(long, long)) alloc_array2_d, w->n_quad_v, w->n_quad_v, w->lists[WORK_DQQ]);
          case WORK_DX:
               GET_WORK1_VECTOR(work->i_list[WORK_DX], w->n_list[WORK_DX], w->max_list[WORK_DX], (void **(*)(long)) alloc_array1_d, w->n_quad_v_x, w->lists[WORK_DX]);
          case WORK_DX2:
               GET_WORK1_VECTOR(work->i_list[WORK_DX2], w->n_list[WORK_DX2], w->max_list[WORK_DX2], (void **(*)(long)) alloc_array1_d, 2 * w->n_quad_v_x, w->lists[WORK_DX2]);
          case WORK_DXX:
               GET_WORK1_MATRIX(work->i_list[WORK_DXX], w->n_list[WORK_DXX], w->max_list[WORK_DXX], (void **(*)(long, long)) alloc_array2_d, w->n_quad_v_x, w->n_quad_v_x, w->lists[WORK_DXX]);
          case WORK_DXX2:
               GET_WORK1_MATRIX(work->i_list[WORK_DXX2], w->n_list[WORK_DXX2], w->max_list[WORK_DXX2], (void **(*)(long, long)) alloc_array2_d, 2 * w->n_quad_v_x, 2 * w->n_quad_v_x, w->lists[WORK_DXX2]);
          case WORK_DD:
               GET_WORK1_VECTOR(work->i_list[WORK_DD], w->n_list[WORK_DD], w->max_list[WORK_DD], (void **(*)(long)) alloc_array1_d, w->n_quad_v + w->n_umus_v, w->lists[WORK_DD]);
          case WORK_DS:
               GET_WORK1_VECTOR(work->i_list[WORK_DS], w->n_list[WORK_DS], w->max_list[WORK_DS], (void **(*)(long)) alloc_array1_d, w->n_stokes, w->lists[WORK_DS]);
          case WORK_DU:
               GET_WORK1_VECTOR(work->i_list[WORK_DU], w->n_list[WORK_DU], w->max_list[WORK_DU], (void **(*)(long)) alloc_array1_d, w->n_umus_v, w->lists[WORK_DU]);
          case WORK_DUX:
               GET_WORK1_MATRIX(work->i_list[WORK_DUX], w->n_list[WORK_DUX], w->max_list[WORK_DUX], (void **(*)(long, long)) alloc_array2_d, w->n_umus_v, w->n_quad_v_x, w->lists[WORK_DUX]);
          case WORK_DDERIVS:
               GET_WORK1_VECTOR(work->i_list[WORK_DDERIVS], w->n_list[WORK_DDERIVS], w->max_list[WORK_DDERIVS], (void **(*)(long)) alloc_array1_d, w->n_derivs, w->lists[WORK_DDERIVS]);
          case WORK_DLAYERS:
               GET_WORK1_VECTOR(work->i_list[WORK_DLAYERS], w->n_list[WORK_DLAYERS], w->max_list[WORK_DLAYERS], (void **(*)(long)) alloc_array1_d, w->n_layers, w->lists[WORK_DLAYERS]);
          case WORK_DBOTH:
               GET_WORK1_MATRIX(work->i_list[WORK_DBOTH], w->n_list[WORK_DBOTH], w->max_list[WORK_DBOTH], (void **(*)(long, long)) alloc_array2_d, w->n_layers, w->n_derivs, w->lists[WORK_DBOTH]);
          case WORK_ZX:
               GET_WORK1_VECTOR(work->i_list[WORK_ZX], w->n_list[WORK_ZX], w->max_list[WORK_ZX], (void **(*)(long)) alloc_array1_dc, w->n_quad_v_x, w->lists[WORK_ZX]);
          case WORK_ZXX:
               GET_WORK1_MATRIX(work->i_list[WORK_ZXX], w->n_list[WORK_ZXX],  w->max_list[WORK_ZXX],  (void **(*)(long, long)) alloc_array2_dc, w->n_quad_v_x, w->n_quad_v_x, w->lists[WORK_ZXX]);
          case WORK_ZU:
               GET_WORK1_VECTOR(work->i_list[WORK_ZU], w->n_list[WORK_ZU], w->max_list[WORK_ZU], (void **(*)(long)) alloc_array1_dc, w->n_umus_v, w->lists[WORK_ZU]);
          case WORK_ZUX:
               GET_WORK1_MATRIX(work->i_list[WORK_ZUX], w->n_list[WORK_ZUX], w->max_list[WORK_ZUX], (void **(*)(long, long)) alloc_array2_dc, w->n_umus_v, w->n_quad_v_x, w->lists[WORK_ZUX]);
          default:
#ifdef DEBUG
               eprintf("ERROR: invalid work type\n");
               exit(1);
#endif
               break;
     }

     return NULL;
}



/*******************************************************************************
 *
 ******************************************************************************/
void *work_get2(work_data *work, enum work_type type,
                enum work_v_type v_type, uchar *flags) {

     WORK_I_TYPE *i;

     int ll;

     int *max;

     int *n;

     void ****data;

     work_shared *w;

     w = work->w;

     switch(v_type) {
          case WORK_DERIVS_V:
               i    = &(work->i_list_v[WORK_DERIVS_V]);
               n    = &(w->n_list_v[WORK_DERIVS_V]);
               max  = &(w->max_list_v[WORK_DERIVS_V]);
               ll   = w->n_derivs;
               data = (void ****) &(w->lists_v[WORK_DERIVS_V]);
               break;
          case WORK_LAYERS_V:
               i    = &(work->i_list_v[WORK_LAYERS_V]);
               n    = &(w->n_list_v[WORK_LAYERS_V]);
               max  = &(w->max_list_v[WORK_LAYERS_V]);
               ll   = w->n_layers + 1;
               data = (void ****) &(w->lists_v[WORK_LAYERS_V]);
               break;
          default:
#ifdef DEBUG
               eprintf("ERROR: invalid work multi type\n");
               exit(1);
#endif
               break;
     }

     switch(type) {
          case WORK_IX:
               return get_work_vectors(i, n, max, &(work->i_list[WORK_IX]), &(w->n_list[WORK_IX]), &(w->max_list[WORK_IX]), (void **(*)(long)) alloc_array1_i, ll, w->n_quad_v_x, data, (void ***) &w->lists[WORK_IX], flags);
          case WORK_IX2:
               return get_work_vectors(i, n, max, &(work->i_list[WORK_IX2]), &(w->n_list[WORK_IX2]), &(w->max_list[WORK_IX2]), (void **(*)(long)) alloc_array1_i, ll, 2 * w->n_quad_v_x, data, (void ***) &w->lists[WORK_IX2], flags);
          case WORK_DQQ:
               return get_work_matrices(i, n, max, &(work->i_list[WORK_DQQ]), &(w->n_list[WORK_DQQ]), &(w->max_list[WORK_DQQ]), (void **(*)(long, long)) alloc_array2_d, ll, w->n_quad_v, w->n_quad_v, data, (void ****) &w->lists[WORK_DQQ], flags);
          case WORK_DX:
               return get_work_vectors(i, n, max, &(work->i_list[WORK_DX]), &(w->n_list[WORK_DX]), &(w->max_list[WORK_DX]), (void **(*)(long)) alloc_array1_d, ll, w->n_quad_v_x, data, (void ***) &w->lists[WORK_DX], flags);
          case WORK_DX2:
               return get_work_vectors(i, n, max, &(work->i_list[WORK_DX2]), &(w->n_list[WORK_DX2]), &(w->max_list[WORK_DX2]), (void **(*)(long)) alloc_array1_d, ll, 2 * w->n_quad_v_x, data, (void ***) &w->lists[WORK_DX2], flags);
          case WORK_DXX:
               return get_work_matrices(i, n, max, &(work->i_list[WORK_DXX]), &(w->n_list[WORK_DXX]), &(w->max_list[WORK_DXX]), (void **(*)(long, long)) alloc_array2_d, ll, w->n_quad_v_x, w->n_quad_v_x, data, (void ****) &w->lists[WORK_DXX], flags);
          case WORK_DXX2:
               return get_work_matrices(i, n, max, &(work->i_list[WORK_DXX2]), &(w->n_list[WORK_DXX2]), &(w->max_list[WORK_DXX2]), (void **(*)(long, long)) alloc_array2_d, ll, 2 * w->n_quad_v_x, 2 * w->n_quad_v_x, data, (void ****) &w->lists[WORK_DXX2], flags);
          case WORK_DD:
               return get_work_vectors(i, n, max, &(work->i_list[WORK_DD]), &(w->n_list[WORK_DD]), &(w->max_list[WORK_DD]), (void **(*)(long)) alloc_array1_d, ll, w->n_quad_v + w->n_umus_v, data, (void ***) &w->lists[WORK_DD], flags);
          case WORK_DS:
               return get_work_vectors(i, n, max, &(work->i_list[WORK_DS]), &(w->n_list[WORK_DS]), &(w->max_list[WORK_DS]), (void **(*)(long)) alloc_array1_d, ll, w->n_stokes, data, (void ***) &w->lists[WORK_DS], flags);
          case WORK_DU:
               return get_work_vectors(i, n, max, &(work->i_list[WORK_DU]), &(w->n_list[WORK_DU]), &(w->max_list[WORK_DU]), (void **(*)(long)) alloc_array1_d, ll, w->n_umus_v, data, (void ***) &w->lists[WORK_DU], flags);
          case WORK_DUX:
               return get_work_matrices(i, n, max, &(work->i_list[WORK_DUX]), &(w->n_list[WORK_DUX]), &(w->max_list[WORK_DUX]), (void **(*)(long, long)) alloc_array2_d, ll, w->n_umus_v, w->n_quad_v_x, data, (void ****) &w->lists[WORK_DUX], flags);
          case WORK_DDERIVS:
               return get_work_vectors(i, n, max, &(work->i_list[WORK_DDERIVS]), &(w->n_list[WORK_DDERIVS]), &(w->max_list[WORK_DDERIVS]), (void **(*)(long)) alloc_array1_d, ll, w->n_derivs, data, (void ***) &w->lists[WORK_DDERIVS], flags);
          case WORK_DLAYERS:
               return get_work_vectors(i, n, max, &(work->i_list[WORK_DLAYERS]), &(w->n_list[WORK_DLAYERS]), &(w->max_list[WORK_DLAYERS]), (void **(*)(long)) alloc_array1_d, ll, w->n_layers, data, (void ***) &w->lists[WORK_DLAYERS], flags);
          case WORK_ZX:
               return get_work_vectors(i, n, max, &(work->i_list[WORK_ZX]), &(w->n_list[WORK_ZX]), &(w->max_list[WORK_ZX]), (void **(*)(long)) alloc_array1_dc, ll, w->n_quad_v_x, data, (void ***) &w->lists[WORK_ZX], flags);
          case WORK_ZXX:
               return get_work_matrices(i, n, max, &(work->i_list[WORK_ZXX]), &(w->n_list[WORK_ZXX]), &(w->max_list[WORK_ZXX]), (void **(*)(long, long)) alloc_array2_dc, ll, w->n_quad_v_x, w->n_quad_v_x, data, (void ****) &w->lists[WORK_ZXX], flags);
          case WORK_ZU:
               return get_work_vectors(i, n, max, &(work->i_list[WORK_ZU]), &(w->n_list[WORK_ZU]), &(w->max_list[WORK_ZU]), (void **(*)(long)) alloc_array1_dc, ll, w->n_umus_v, data, (void ***) &w->lists[WORK_ZU], flags);
          case WORK_ZUX:
               return get_work_matrices(i, n, max, &(work->i_list[WORK_ZUX]), &(w->n_list[WORK_ZUX]), &(w->max_list[WORK_ZUX]), (void **(*)(long, long)) alloc_array2_dc, ll, w->n_umus_v, w->n_quad_v_x, data, (void ****) &w->lists[WORK_ZUX], flags);
          default:
#ifdef DEBUG
               eprintf("ERROR: invalid work single type\n");
               exit(1);
#endif
               break;
     }

     return NULL;
}



/*******************************************************************************
 *
 ******************************************************************************/
void *work_get3(work_data *work, enum work_type type,
                enum work_v_type v_type, uchar **flags) {

     WORK_I_TYPE *i;

     int ll1;
     int ll2;

     int *max;

     int *n;

     void *****data2;

     work_shared *w;

     w = work->w;

     if (v_type != WORK_BOTH_V) {
          eprintf("ERROR: invalid work v type\n");
          exit(1);
     }

     i     = &(work->i_list_v[WORK_BOTH_V]);
     n     = &(w->n_list_v[WORK_BOTH_V]);
     max   = &(w->max_list_v[WORK_BOTH_V]);
     ll1   = w->n_layers;
     ll2   = w->n_derivs;
     data2 = (void *****) &(w->lists_v[WORK_BOTH_V]);

     switch(type) {
          case WORK_IX:
               return get_work_vectors2(i, n, max, &(work->i_list[WORK_IX]), &(w->n_list[WORK_IX]), &(w->max_list[WORK_IX]), (void **(*)(long)) alloc_array1_i, ll1, ll2, w->n_quad_v_x, data2, (void ***) &w->lists[WORK_IX], flags);
          case WORK_IX2:
               return get_work_vectors2(i, n, max, &(work->i_list[WORK_IX2]), &(w->n_list[WORK_IX2]), &(w->max_list[WORK_IX2]), (void **(*)(long)) alloc_array1_i, ll1, ll2, 2 * w->n_quad_v_x, data2, (void ***) &w->lists[WORK_IX2], flags);
          case WORK_DQQ:
               return get_work_matrices2(i, n, max, &(work->i_list[WORK_DQQ]), &(w->n_list[WORK_DQQ]), &(w->max_list[WORK_DQQ]), (void **(*)(long, long)) alloc_array2_d, ll1, ll2, w->n_quad_v, w->n_quad_v, data2, (void ****) &w->lists[WORK_DQQ], flags);
          case WORK_DX:
               return get_work_vectors2(i, n, max, &(work->i_list[WORK_DX]), &(w->n_list[WORK_DX]), &(w->max_list[WORK_DX]), (void **(*)(long)) alloc_array1_d, ll1, ll2, w->n_quad_v_x, data2, (void ***) &w->lists[WORK_DX], flags);
          case WORK_DX2:
               return get_work_vectors2(i, n, max, &(work->i_list[WORK_DX2]), &(w->n_list[WORK_DX2]), &(w->max_list[WORK_DX2]), (void **(*)(long)) alloc_array1_d, ll1, ll2, 2 * w->n_quad_v_x, data2, (void ***) &w->lists[WORK_DX2], flags);
          case WORK_DXX:
               return get_work_matrices2(i, n, max, &(work->i_list[WORK_DXX]), &(w->n_list[WORK_DXX]), &(w->max_list[WORK_DXX]), (void **(*)(long, long)) alloc_array2_d, ll1, ll2, w->n_quad_v_x, w->n_quad_v_x, data2, (void ****) &w->lists[WORK_DXX], flags);
          case WORK_DXX2:
               return get_work_matrices2(i, n, max, &(work->i_list[WORK_DXX2]), &(w->n_list[WORK_DXX2]), &(w->max_list[WORK_DXX2]), (void **(*)(long, long)) alloc_array2_d, ll1, ll2, 2 * w->n_quad_v_x, 2 * w->n_quad_v_x, data2, (void ****) &w->lists[WORK_DXX2], flags);
          case WORK_DD:
               return get_work_vectors2(i, n, max, &(work->i_list[WORK_DD]), &(w->n_list[WORK_DD]), &(w->max_list[WORK_DD]), (void **(*)(long)) alloc_array1_d, ll1, ll2, w->n_quad_v + w->n_umus_v, data2, (void ***) &w->lists[WORK_DD], flags);
          case WORK_DS:
               return get_work_vectors2(i, n, max, &(work->i_list[WORK_DS]), &(w->n_list[WORK_DS]), &(w->max_list[WORK_DS]), (void **(*)(long)) alloc_array1_d, ll1, ll2, w->n_stokes, data2, (void ***) &w->lists[WORK_DS], flags);
          case WORK_DU:
               return get_work_vectors2(i, n, max, &(work->i_list[WORK_DU]), &(w->n_list[WORK_DU]), &(w->max_list[WORK_DU]), (void **(*)(long)) alloc_array1_d, ll1, ll2, w->n_umus_v, data2, (void ***) &w->lists[WORK_DU], flags);
          case WORK_DUX:
               return get_work_matrices2(i, n, max, &(work->i_list[WORK_DUX]), &(w->n_list[WORK_DUX]), &(w->max_list[WORK_DUX]), (void **(*)(long, long)) alloc_array2_d, ll1, ll2, w->n_umus_v, w->n_quad_v_x, data2, (void ****) &w->lists[WORK_DUX], flags);
          case WORK_DDERIVS:
               return get_work_vectors2(i, n, max, &(work->i_list[WORK_DDERIVS]), &(w->n_list[WORK_DDERIVS]), &(w->max_list[WORK_DDERIVS]), (void **(*)(long)) alloc_array1_d, ll1, ll2, w->n_derivs, data2, (void ***) &w->lists[WORK_DDERIVS], flags);
          case WORK_DLAYERS:
               return get_work_vectors2(i, n, max, &(work->i_list[WORK_DLAYERS]), &(w->n_list[WORK_DLAYERS]), &(w->max_list[WORK_DLAYERS]), (void **(*)(long)) alloc_array1_d, ll1, ll2, w->n_layers, data2, (void ***) &w->lists[WORK_DLAYERS], flags);
          case WORK_ZX:
               return get_work_vectors2(i, n, max, &(work->i_list[WORK_ZX]), &(w->n_list[WORK_ZX]), &(w->max_list[WORK_ZX]), (void **(*)(long)) alloc_array1_dc, ll1, ll2, w->n_quad_v_x, data2, (void ***) &w->lists[WORK_ZX], flags);
          case WORK_ZXX:
               return get_work_matrices2(i, n, max, &(work->i_list[WORK_ZXX]), &(w->n_list[WORK_ZXX]), &(w->max_list[WORK_ZXX]), (void **(*)(long, long)) alloc_array2_dc, ll1, ll2, w->n_quad_v_x, w->n_quad_v_x, data2, (void ****) &w->lists[WORK_ZXX], flags);
          case WORK_ZU:
               return get_work_vectors2(i, n, max, &(work->i_list[WORK_ZU]), &(w->n_list[WORK_ZU]), &(w->max_list[WORK_ZU]), (void **(*)(long)) alloc_array1_dc, ll1, ll2, w->n_umus_v, data2, (void ***) &w->lists[WORK_ZU], flags);
          case WORK_ZUX:
               return get_work_matrices2(i, n, max, &(work->i_list[WORK_ZUX]), &(w->n_list[WORK_ZUX]), &(w->max_list[WORK_ZUX]), (void **(*)(long, long)) alloc_array2_dc, ll1, ll2, w->n_umus_v, w->n_quad_v_x, data2, (void ****) &w->lists[WORK_ZUX], flags);
          default:
#ifdef DEBUG
               eprintf("ERROR: invalid work single type\n");
               exit(1);
#endif
               break;
     }

     return NULL;
}
