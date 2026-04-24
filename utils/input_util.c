/*******************************************************************************
**
**    Copyright (C) 2007-2020 Greg McGarragh <greg.mcgarragh@colostate.edu>
**
**    This source code is licensed under the GNU General Public License (GPL),
**    Version 3.  See the file COPYING for more details.
**
*******************************************************************************/

#include <gutil.h>

#include <rtutil_scat_io.h>

#include <xrtm_fd_interface.h>
#include <xrtm_interface.h>
#include <xrtm_support.h>

#include "input_local.h"
#include "input_util.h"
#define YY_NO_UNISTD_H
#include "input_yylex.h"


/*******************************************************************************
 *
 ******************************************************************************/
void input_error(locus_data *locus, const char *s, ...) {

     char *format;

     size_t size;

     va_list ap;

     size = 7 + strlen(locus->file) + 1 + 16 + 1 + strlen(locus->statement) + 2 + strlen(s) + 2 + 1;

     format = (char *) malloc(size);

     snprintf(format, size, "ERROR: %s:%d:%s, %s\n", locus->file, locus->line, locus->statement, s);

     va_start(ap, s);

     vfprintf(stderr, format, ap);

     va_end(ap);
}


/*******************************************************************************
 *
 ******************************************************************************/
input_type_data list_null = {NULL, 0, INPUT_TYPE_VOID, 0, {0, 0, 0, 0}, {NULL}};


#undef  GSTRUCT_BASIC_TYPE
#define GSTRUCT_NILL list_null
#undef  GSTRUCT_POINTER_TYPE
#define GSTRUCT_PREFIX f
#define GSTRUCT_TYPE input_type_data

#define IMPLEMENTATION
#include <gstruct_dclist.h>
#undef IMPLEMENTATION

#undef GSTRUCT_BASIC_TYPE
#undef GSTRUCT_NILL
#undef GSTRUCT_POINTER_TYPE
#undef GSTRUCT_PREFIX
#undef GSTRUCT_TYPE



/*******************************************************************************
 *
 ******************************************************************************/
static input_type_data input_type_strdup(input_data *input, const char *s) {

     input_type_data input_type;

     input_type.name  = NULL;
     input_type.flag  = 0;
     input_type.type  = INPUT_TYPE_STRING;
     input_type.order = 0;

     input_type.d.s = strdup(s);

     return input_type;
}



/*******************************************************************************
 *
 ******************************************************************************/
static void parse_char(locus_data *locus, char c) {

     lex_type_data lex_type;

     if (yylex(locus, &lex_type) != c) {
          input_error(locus, "expected '%c' at %s", c, get_yytext());
          exit(1);
     }
}



static void parse_equal(input_data *input) {

     if (input->use_equal)
          parse_char(&input->locus, '=');
}



static void parse_comma(input_data *input) {

     parse_char(&input->locus, ',');
}



/*******************************************************************************
 *
 ******************************************************************************/
static void copy_dimens(int *dimens2, int *dimens1, int n) {

     int i;

     for (i = 0; i < n; ++i)
          dimens2[i] = dimens1[i];
}



static int parse_scaler_int(input_data *input);

static int dimens_parse(input_data *input, int *dimens) {

     int r;

     lex_type_data data;

     r = yylex(&input->locus, &data);

     if (r != '[') {
          yyrewind(r, &data);
          return  0;
     }
     else {
          r = yylex(&input->locus, &data);

          if (r == ':')
               dimens[0] = -1;
          else {
               yyrewind(r, &data);
               dimens[0] = parse_scaler_int(input);
          }

          r = yylex(&input->locus, &data);

          if (r != ']') {
               input_error(&input->locus, "expected ']' at %s", get_yytext());
               exit(1);
          }

          return 1 + dimens_parse(input, dimens + 1);
     }
}



/*******************************************************************************
 *
 ******************************************************************************/
#undef TYPE
#undef NAME
#undef X_VAL
#undef FORMAT

#define C_TYPE char *
#define LEX_TYPE LEX_TYPE_STRING
#define LEX_MEMBER s
#define FILE_NAME string
#define FILE_STRING "string"
#define FILE_TYPE INPUT_TYPE_STRING
#define FILE_MEMBER s
#define FORMAT "\'%s\'"
#define FREE_SCALER
#include "parse_data.h"
#undef C_TYPE
#undef LEX_TYPE
#undef LEX_MEMBER
#undef FILE_NAME
#undef FILE_STRING
#undef FILE_TYPE
#undef FILE_MEMBER
#undef FORMAT
#undef FREE_SCALER

#define C_TYPE int
#define LEX_TYPE LEX_TYPE_INT
#define LEX_MEMBER i
#define FILE_NAME int
#define FILE_STRING "int"
#define FILE_TYPE INPUT_TYPE_INT
#define FILE_MEMBER i
#define FORMAT "%d"
#include "parse_data.h"
#undef C_TYPE
#undef LEX_TYPE
#undef LEX_MEMBER
#undef FILE_NAME
#undef FILE_STRING
#undef FILE_TYPE
#undef FILE_MEMBER
#undef FORMAT

#define C_TYPE double
#define LEX_TYPE  LEX_TYPE_INT
#define LEX_TYPE2 LEX_TYPE_DOUBLE
#define LEX_MEMBER  i
#define LEX_MEMBER2 d
#define FILE_NAME double
#define FILE_STRING "double"
#define FILE_TYPE INPUT_TYPE_DOUBLE
#define FILE_MEMBER d
#define FORMAT "%e"
#include "parse_data.h"
#undef C_TYPE
#undef LEX_TYPE
#undef LEX_TYPE2
#undef LEX_MEMBER
#undef LEX_MEMBER2
#undef FILE_NAME
#undef FILE_STRING
#undef FILE_TYPE
#undef FILE_MEMBER
#undef FORMAT


/*******************************************************************************
 *
 ******************************************************************************/
static input_type_data parse_type(input_data *input, int type, int n, int *dimens) {

     switch(type) {
          case INPUT_TYPE_STRING:
               return parse_string(input, n, dimens);
               break;
          case INPUT_TYPE_INT:
               return parse_int(input, n, dimens);
               break;
          case INPUT_TYPE_DOUBLE:
               return parse_double(input, n, dimens);
               break;
          default:
               fprintf(stderr, "ERROR: end of switch block error\n");
               exit(1);
     }
}



static void free_type(input_type_data *input_type) {

     switch(input_type->type) {
          case INPUT_TYPE_STRING:
               free_string(input_type);
               break;
          case INPUT_TYPE_INT:
               free_int(input_type);
               break;
          case INPUT_TYPE_DOUBLE:
               free_double(input_type);
               break;
          default:
               fprintf(stderr, "ERROR: end of switch block error\n");
               exit(1);
     }
}




static void print_type(FILE *fp, format_data *d, input_type_data *input_type) {

     switch(input_type->type) {
          case INPUT_TYPE_STRING:
               print_string(fp, d, input_type);
               break;
          case INPUT_TYPE_INT:
               print_int(fp, d, input_type);
               break;
          case INPUT_TYPE_DOUBLE:
               print_double(fp, d, input_type);
               break;
          default:
               fprintf(stderr, "ERROR: end of switch block error\n");
               exit(1);
     }
}



/*******************************************************************************
 *
 ******************************************************************************/
static fdclist *list_parse1(input_data *input, int type, int n, int *dimens) {

     int r;

     lex_type_data lex_type;

     fdclist *list;

     list = fdclist_create();
     fdclist_set_free(list, (void (*)(void *)) free_type);
     fdclist_set_print(list, (int (*)(FILE *, const void *)) print_type);

     r = yylex(&input->locus, &lex_type);

     if (r != '(') {
          input_error(&input->locus, "expected '(' at %s", get_yytext());
          exit(1);
     }

     do {
          fdclist_insert_end(list, parse_type(input, type, n, dimens));
     } while ((r = yylex(&input->locus, &lex_type)) == ',');

     if (r != ')') {
          input_error(&input->locus, "expected ')' at %s", get_yytext());
          exit(1);
     }

     return list;
}



static fdclist *list_prase2(input_data *input, int type, int n, int *dimens) {

     int r;

     lex_type_data lex_type;

     fdclist *list;

     list = fdclist_create();
     fdclist_set_free(list, (void (*)(void *)) free_type);
     fdclist_set_print(list, (int (*)(FILE *, const void *)) print_type);

     do {
          fdclist_insert_end(list, parse_type(input, type, n, dimens));
     } while ((r = yylex(&input->locus, &lex_type)) == ',');

     yyrewind(r, &lex_type);

     return list;
}



static fdclist *parse_list(input_data *input, int type, int n, int *dimens) {

     if (! input->use_array2)
          return list_parse1(input, type, n, dimens);
     else
          return list_prase2(input, type, n, dimens);

}



/*******************************************************************************
 *
 ******************************************************************************/
static void print_data(FILE *fp, format_data *d, input_type_data *data) {

     switch(data->type) {
          case INPUT_TYPE_STRING:
               print_string(fp, d, data);
               break;
          case INPUT_TYPE_INT:
               print_int(fp, d, data);
               break;
          case INPUT_TYPE_DOUBLE:
               print_double(fp, d, data);
               break;
          default:
               break;
     }
}


static void print_list(FILE *fp, format_data *d, fdclist *list) {

     fdcelem *elem;
     fdcelem *first;

     elem = first = fdclist_first(list);

     if (! d->use_array2)
          fprintf(fp, "(");

     if (first) {
          print_data(fp, d, fdcelem_pointer(elem));

          while ((elem = fdcelem_advance(elem)) != first) {
              fprintf(fp, ", ");
              print_data(fp, d, fdcelem_pointer(elem));
          }
     }

     if (! d->use_array2)
          fprintf(fp, ")");
}



/*******************************************************************************
 *
 ******************************************************************************/
static long xrtm_string_list_to_mask(input_data *input, fdclist *list, long (*name_to_mask)(const char *)) {

     long mask;

     long options;

     fdcelem *elem;
     fdcelem *first;

     options = 0;

     if ((elem = first = fdclist_first(list))) {
          do {
               if ((mask = name_to_mask(fdcelem_value(elem).d.s)) == 0)
                    return -1;
               options |= mask;
          } while ((elem = fdcelem_advance(elem)) != first);
     }

     return options;
}



static long *xrtm_string_list_to_vector(input_data *input, fdclist *list, long (*name_to_mask)(const char *)) {

     int i;

     long mask;

     long *options;

     fdcelem *elem;
     fdcelem *first;

     options = 0;

     options = alloc_array1_l(fdclist_count(list));

     if ((elem = first = fdclist_first(list))) {
          i = 0;
          do {
               if ((mask = name_to_mask(fdcelem_value(elem).d.s)) == 0)
                    return NULL;
               options[i++] = mask;
          } while ((elem = fdcelem_advance(elem)) != first);
     }

     return options;
}



static fdclist *xrtm_mask_to_string_list(input_data *input, long options, int (*count)(), long (*index_to_mask)(int), const char *(*index_to_name)(int)) {

     int i;
     int n;

     fdclist *list;

     list = fdclist_create();
     fdclist_set_free(list, (void (*)(void *)) free_type);
     fdclist_set_print(list, (int (*)(FILE *, const void *)) print_type);

     n = count();

     for (i = 0; i < n; ++i) {
          if (index_to_mask(i) & options) {
               fdclist_insert_end(list, input_type_strdup(input, index_to_name(i)));
          }
     }

     return list;
}



static long xrtm_string_list_to_options_mask(input_data *input, fdclist *list) {

     return xrtm_string_list_to_mask(input, list, (long (*)(const char*)) xrtm_option_name_to_mask);
}



static long *xrtm_string_list_to_options_vector(input_data *input, fdclist *list) {

     return xrtm_string_list_to_vector(input, list, (long (*)(const char*)) xrtm_option_name_to_mask);
}



static fdclist *xrtm_options_mask_to_string_list(input_data *input, long options) {

     return xrtm_mask_to_string_list(input, options, xrtm_option_n, (long (*)(int)) xrtm_option_index_to_mask, xrtm_option_index_to_name);
}



static long xrtm_string_list_to_solvers_mask(input_data *input, fdclist *list) {

     return xrtm_string_list_to_mask(input, list, (long (*)(const char*)) xrtm_solver_name_to_mask);
}



static long *xrtm_string_list_to_solvers_vector(input_data *input, fdclist *list) {

     return xrtm_string_list_to_vector(input, list, (long (*)(const char*)) xrtm_solver_name_to_mask);
}



static fdclist *xrtm_solvers_mask_to_string_list(input_data *input, long solvers) {

     return xrtm_mask_to_string_list(input, solvers, xrtm_solver_n, (long (*)(int)) xrtm_solver_index_to_mask, xrtm_solver_index_to_name);
}



static int *xrtm_string_list_to_kernels_vector(input_data *input, fdclist *list) {

     int i;

     int code;

     int *kernels;

     fdcelem *elem;
     fdcelem *first;

     kernels = alloc_array1_i(fdclist_count(list));

     if ((elem = first = fdclist_first(list))) {
          i = 0;
          do {
               if ((code = xrtm_kernel_name_to_value(fdcelem_value(elem).d.s)) < 0)
                    return NULL;
               kernels[i++] = code;
          } while ((elem = fdcelem_advance(elem)) != first);
     }

     return kernels;
}



static fdclist *xrtm_kernels_vector_to_string_list(input_data *input, enum xrtm_kernel_type *kernels, int n) {

     int i;

     fdclist *list;

     list = fdclist_create();
     fdclist_set_free(list, (void (*)(void *)) free_type);
     fdclist_set_print(list, (int (*)(FILE *, const void *)) print_type);

     for (i = 0; i < n; ++i)
          fdclist_insert_end(list, input_type_strdup(input, xrtm_kernel_value_to_name(kernels[i])));

     return list;
}



/*******************************************************************************
 *
 ******************************************************************************/
static void parse_scaler_param_int(input_data *input, const char *name, int (*xrtm_set_x)(xrtm_data *, int)) {

     int i;

     parse_equal(input);

     i = parse_scaler_int(input);

     if (xrtm_set_x(input->d, i)) {
          input_error(&input->locus, "xrtm_set_%s()\n", name);
          exit(1);
     }
}



static void parse_scaler_param_int2(input_data *input, const char *name, int (*xrtm_set_x)(xrtm_data *, int, int)) {

     int i1;
     int i2;

     parse_equal(input);

     i1 = parse_scaler_int(input);
     parse_comma(input);
     i2 = parse_scaler_int(input);

     if (xrtm_set_x(input->d, i1, i2)) {
          input_error(&input->locus, "xrtm_set_%s()\n", name);
          exit(1);
     }
}



static void parse_sos_params(input_data *input) {

     int i;
     double d1;
     double d2;

     parse_equal(input);

     i  = parse_scaler_int(input);
     parse_comma(input);
     d1 = parse_scaler_double(input);
     parse_comma(input);
     d2 = parse_scaler_double(input);

     if (xrtm_set_sos_params(input->d, i, d1, d2)) {
          input_error(&input->locus, "xrtm_set_sos_params()\n");
          exit(1);
     }
}



static void parse_scaler_param_double(input_data *input, const char *name, int (*xrtm_set_x)(xrtm_data *, double)) {

     double x;

     parse_equal(input);

     x = parse_scaler_double(input);

     if (xrtm_set_x(input->d, x)) {
          input_error(&input->locus, "xrtm_set_%s()\n", name);
          exit(1);
     }
}



static void parse_array_param_int(input_data *input, const char *name, int n, int (*xrtm_set_x)(xrtm_data *, void *)) {

     int *x1;

     parse_equal(input);

     x1 = parse_array1_int(input, n);

     if (xrtm_set_x(input->d, x1)) {
          input_error(&input->locus, "xrtm_set_%s()\n", name);
          exit(1);
     }

     free_array1_int(x1, n);
}



static void parse_array_param_double(input_data *input, void *d, const char *name, int n, int (*xrtm_set_x)(void *, void *)) {

     double *x1;

     parse_equal(input);

     x1 = parse_array1_double(input, n);

     if (xrtm_set_x(d, x1)) {
          input_error(&input->locus, "xrtm_set_%s()\n", name);
          exit(1);
     }

     free_array1_double(x1, n);
}



static void parse_array_param_double2(input_data *input, void *d, const char *name, int m, int n, int (*xrtm_set_x)(void *, void *)) {

     double **x1;

     parse_equal(input);

     x1 = parse_array2_double(input, m, n);

     if (xrtm_set_x(d, x1)) {
          input_error(&input->locus, "xrtm_set_%s()\n", name);
          exit(1);
     }

     free_array2_double(x1, m, n);
}



static void parse_layer_param_double(input_data *input, const char *name, int (*xrtm_set_x_1)(xrtm_data *, int, double), int (*xrtm_set_x_n)(xrtm_data *, double *)) {

     int n;

     int n_layers;

     int dimens[MAX_DIMENS];

     double x;
     double *x1;

     n = dimens_parse(input, dimens);

     if (n == 1 && dimens[0] > -1) {
          parse_equal(input);
          x = parse_scaler_double(input);
          if (xrtm_set_x_1(input->d, dimens[0], x)) {
               input_error(&input->locus, "xrtm_set_%s_1()\n", name);
               exit(1);
          }
     }
     else
     if (n == 0 || (n == 1 && dimens[0] == -1)) {
          parse_equal(input);
          n_layers = xrtm_get_n_layers(input->d);
          x1 = parse_array1_double(input, n_layers);
          if (xrtm_set_x_n(input->d, x1)) {
               input_error(&input->locus, "xrtm_set_%s_n()\n", name);
               exit(1);
          }
          free_array1_double(x1, n_layers);
     }
     else {
          input_error(&input->locus, "invalid dimensions");
          exit(1);
     }
}



static void parse_layer_param_l_double(input_data *input, void *d, int n_derivs, const char *name, const char *type, int (*xrtm_set_x_x_11)(void *, int, int, double), int (*xrtm_set_x_x_n1)(void *, int, double *), int (*xrtm_set_x_x_1n)(void *, int, double *), int (*xrtm_set_x_x_nn)(void *, double **)) {

     int n;

     int n_layers;

     int dimens[MAX_DIMENS];

     double x;
     double *x1;
     double **x2;

     n = dimens_parse(input, dimens);

     if (n == 2 && dimens[0] > -1 && dimens[1] > -1) {
          parse_equal(input);
          x = parse_scaler_double(input);
          if (xrtm_set_x_x_11(d, dimens[0], dimens[1], x)) {
               input_error(&input->locus, "xrtm_set_%s_%s_11()\n", name, type);
               exit(1);
          }
     }
     else
     if ( n == 2 && dimens[0] == -1 && dimens[1] > -1) {
          parse_equal(input);
          n_layers = xrtm_get_n_layers(input->d);
          x1 = parse_array1_double(input, n_layers);
          if (xrtm_set_x_x_n1(d, dimens[1], x1)) {
               input_error(&input->locus, "xrtm_set_%s_%s_n1()\n", name, type);
               exit(1);
          }
          free_array1_double(x1, n_layers);
     }
     else
     if ((n == 1 && dimens[0] > -1) || (n == 2 && dimens[0] > -1 && dimens[1] == -1)) {
          parse_equal(input);
          x1 = parse_array1_double(input, n_derivs);
          if (xrtm_set_x_x_1n(d, dimens[0], x1)) {
               input_error(&input->locus, "xrtm_set_%s_%s_1n()\n", name, type);
               exit(1);
          }
          free_array1_double(x1, n_derivs);
     }
     else
     if ((n == 0) || (n == 1 && dimens[0] == -1 && dimens[1] == -1) || (n == 1 && dimens[0] == -1) || (n == 2 && dimens[0] == -1 && dimens[1] == -1)) {
          parse_equal(input);
          n_layers = xrtm_get_n_layers(input->d);
          x2 = parse_array2_double(input, n_layers, n_derivs);
          if (xrtm_set_x_x_nn(d, x2)) {
               input_error(&input->locus, "xrtm_set_%s_%s_nn()\n", name, type);
               exit(1);
          }
          free_array2_double(x2, n_layers, n_derivs);
     }
     else {
          input_error(&input->locus, "invalid dimensions");
          exit(1);
     }
}



static void trans_array_double(double **a, double **b, int m, int n) {

     int i;
     int j;

     for (i = 0; i < m; ++i) {
          for (j = 0; j < n; ++j) {
               b[j][i] = a[i][j];
          }
     }
}



static double **parse_coef_1(input_data *input, int n_elem, int n_coef) {

     double **a1;
     double **a2;

     a1 = parse_array2_double(input, n_coef, n_elem);

     a2 = alloc_array2_d(n_elem, n_coef);

     trans_array_double(a1, a2, n_coef, n_elem);

     free_array2_double(a1, n_coef, n_elem);

     return a2;
}



static double ***parse_coef_2(input_data *input, int n_x, int n_elem, int n_coef) {

     int i;

     double ***a1;
     double ***a2;

     a1 = parse_array3_double(input, n_x, n_coef, n_elem);

     a2 = alloc_array3_d(n_x, n_elem, n_coef);

     for (i = 0; i < n_x; ++i)
          trans_array_double(a1[i], a2[i], n_coef, n_elem);

     free_array3_double(a1, n_x, n_coef, n_elem);

     return a2;
}



static double ****parse_coef_3(input_data *input, int n_x, int n_y, int n_elem, int n_coef) {

     int i;
     int j;

     double ****a1;
     double ****a2;

     a1 = parse_array4_double(input, n_x, n_y, n_coef, n_elem);

     a2 = alloc_array4_d(n_x, n_y, n_elem, n_coef);

     for (i = 0; i < n_x; ++i) {
          for (j = 0; j < n_y; ++j) {
               trans_array_double(a1[i][j], a2[i][j], n_coef, n_elem);
          }
     }

     free_array4_double(a1, n_x, n_y, n_coef, n_elem);

     return a2;
}



static void parse_coef(input_data *input) {

     int i;
     int n;

     int n_elem;
     int n_coef;
     int n_layers;

     int dimens[MAX_DIMENS];

     int *n_coef_layer;

     double **x2;
     double ***x3;

     n = dimens_parse(input, dimens);

     n_elem = xrtm_get_options(input->d) & XRTM_OPTION_VECTOR ? 6 : 1;

     if (n == 1 && dimens[0] > -1) {
          parse_equal(input);
          n_coef = parse_scaler_int(input);
          parse_comma(input);
          x2 = parse_coef_1(input, n_elem, n_coef);
          if (xrtm_set_coef_1(input->d, dimens[0], n_coef, x2)) {
               input_error(&input->locus, "xrtm_set_coef_1()");
               exit(1);
          }
          free_array2_d(x2);

          free(input->md->coef_files[dimens[0]]);
          input->md->coef_files[dimens[0]] = NULL;
     }
     else
     if (n == 0 || (n == 1 && dimens[0] == -1)) {
          parse_equal(input);
          n_coef = xrtm_get_max_coef(input->d);
          n_layers = xrtm_get_n_layers(input->d);
          n_coef_layer = alloc_array1_i(n_layers);
          init_array1_i(n_coef_layer, n_layers, n_coef);
          x3 = parse_coef_2(input, n_layers, n_elem, n_coef);
          if (xrtm_set_coef_n(input->d, n_coef_layer, x3)) {
               input_error(&input->locus, "xrtm_set_coef_n()");
               exit(1);
          }
          free_array1_i(n_coef_layer);
          free_array3_d(x3);

          for (i = 0; i < n_layers; ++i) {
               free(input->md->coef_files[i]);
               input->md->coef_files[i] = NULL;
          }
     }
     else {
          input_error(&input->locus, "invalid dimensions");
          exit(1);
     }
}



static void parse_coef_l(input_data *input, void *d, int n_derivs, const char *type, int (*xrtm_set_coef_x_11)(void *, int, int, double **), int (*xrtm_set_coef_x_n1)(void *, int, double ***), int (*xrtm_set_coef_x_1n)(void *, int, double ***), int (*xrtm_set_coef_x_nn)(void *, double ****)) {

     int i;
     int j;
     int n;

     int n_elem;
     int n_coef;
     int n_layers;

     int dimens[MAX_DIMENS];

     double **x2;
     double ***x3;
     double ****x4;

     n = dimens_parse(input, dimens);

     n_elem = xrtm_get_options(input->d) & XRTM_OPTION_VECTOR ? 6 : 1;
     n_coef = xrtm_get_max_coef(input->d);

     if (n == 2 && dimens[0] > -1 && dimens[1] > -1) {
          parse_equal(input);
          n_coef = xrtm_get_n_coef(input->d, dimens[0]);
          x2 = parse_coef_1(input, n_elem, n_coef);
          if (xrtm_set_coef_l_11(input->d, dimens[0], dimens[1], x2)) {
               input_error(&input->locus, "xrtm_set_coef_%s_11()", type);
               exit(1);
          }
          free_array2_d(x2);

          free(input->md->coef_files_l[dimens[0]][dimens[1]]);
          input->md->coef_files_l[dimens[0]][dimens[1]] = NULL;
     }
     else
     if ( n == 2 && dimens[0] == -1 && dimens[1] > -1) {
          parse_equal(input);
          n_layers = xrtm_get_n_layers(input->d);
          x3 = parse_coef_2(input, n_layers, n_elem, n_coef);
          if (xrtm_set_coef_l_n1(input->d, dimens[1], x3)) {
               input_error(&input->locus, "xrtm_set_coef_%s_n1()", type);
               exit(1);
          }
          free_array3_d(x3);

          for (i = 0; i < n_layers; ++i) {
               free(input->md->coef_files_l[i][dimens[1]]);
               input->md->coef_files_l[i][dimens[1]] = NULL;
          }
     }
     else
     if ((n == 1 && dimens[0] > -1) || (n == 2 && dimens[0] > -1 && dimens[1] == -1)) {
          parse_equal(input);
          n_coef = xrtm_get_n_coef(input->d, dimens[0]);
          x3 = parse_coef_2(input, n_derivs, n_elem, n_coef);
          if (xrtm_set_coef_l_1n(input->d, dimens[0], x3)) {
               input_error(&input->locus, "xrtm_set_coef_%s_1n()", type);
               exit(1);
          }
          free_array3_d(x3);

          for (i = 0; i < n_derivs; ++i) {
               free(input->md->coef_files_l[dimens[0]][i]);
               input->md->coef_files_l[dimens[0]][i] = NULL;
          }
     }
     else
     if ((n == 0) || (n == 1 && dimens[0] == -1 && dimens[1] == -1) || (n == 1 && dimens[0] == -1) || (n == 2 && dimens[0] == -1 && dimens[1] == -1)) {
          parse_equal(input);
          n_layers = xrtm_get_n_layers(input->d);
          x4 = parse_coef_3(input, n_layers, n_derivs, n_elem, n_coef);
          if (xrtm_set_coef_l_nn(input->d, x4)) {
               input_error(&input->locus, "xrtm_set_coef_%s_nn()", type);
               exit(1);
          }
          free_array4_d(x4);

          for (i = 0; i < n_layers; ++i) {
               for (j = 0; j < n_derivs; ++j) {
                    free(input->md->coef_files_l[i][j]);
                    input->md->coef_files_l[i][j] = NULL;
               }
          }
     }
     else {
          input_error(&input->locus, "invalid dimensions");
          exit(1);
     }
}



static void parse_coef_files(input_data *input) {

     char *x;
     char **x1;

     int i;
     int n;

     int flag;
/*
     int n_elem;
*/
     int n_coef;
     int n_layers;

     int dimens[MAX_DIMENS];

     int n_coef_layer;
     int *n_coef_layer2;

     double **coef1;
     double ***coef2;

     n = dimens_parse(input, dimens);
/*
     n_elem = xrtm_get_options(input->d) & XRTM_OPTION_VECTOR ? 6 : 1;
*/
     n_coef = xrtm_get_max_coef(input->d);

     if (n == 1 && dimens[0] > -1) {
          parse_equal(input);

          x = parse_scaler_string(input);

          if ((n_coef_layer = load_scat_coefs(x, n_coef, &coef1, &flag)) < 0) {
               fprintf(stderr, "ERROR: load_scat_coefs(), layer %d\n", dimens[0]);
               exit(1);
          }

          if (xrtm_set_coef_1(input->d, dimens[0], n_coef_layer, coef1)) {
               input_error(&input->locus, "xrtm_set_coef_1()");
               exit(1);
          }

          free(input->md->coef_files[dimens[0]]);
          input->md->coef_files[dimens[0]] = strdup(x);

          free(x);
          free_array2_d(coef1);
     }
     else
     if (n == 0 || (n == 1 && dimens[0] == -1)) {
          parse_equal(input);

          n_layers = xrtm_get_n_layers(input->d);

          x1 = parse_array1_string(input, n_layers);

          n_coef_layer2 = alloc_array1_i(n_layers);

          coef2 = (double ***) alloc_array1(n_layers, sizeof(double **));

          for (i = 0; i < n_layers; ++i) {
               if ((n_coef_layer2[i] = load_scat_coefs(x1[i], n_coef, &coef2[i], &flag)) < 0) {
                    fprintf(stderr, "ERROR: load_scat_coefs(), layer %d\n", i);
                    exit(1);
               }
          }

          if (xrtm_set_coef_n(input->d, n_coef_layer2, coef2)) {
               input_error(&input->locus, "xrtm_set_coef_n()");
               exit(1);
          }

          for (i = 0; i < n_layers; ++i) {
               free(input->md->coef_files[i]);
               input->md->coef_files[i] = strdup(x1[i]);
          }

          free_array1_string(x1, n_layers);

          free_array1_i(n_coef_layer2);
          for (i = 0; i < n_layers; ++i) {
               free_array2_d(coef2[i]);
          }
          free_array1(coef2);
     }
     else {
          input_error(&input->locus, "invalid dimensions");
          exit(1);
     }
}



static void parse_coef_files_l(input_data *input, void *d, int n_derivs, const char *type, int (*xrtm_set_coef_x_11)(void *, int, int, double **), int (*xrtm_set_coef_x_n1)(void *, int, double ***), int (*xrtm_set_coef_x_1n)(void *, int, double ***), int (*xrtm_set_coef_x_nn)(void *, double ****)) {

     char *x;
     char **x1;
     char ***x2;

     int i;
     int j;
     int n;

     int flag;
/*
     int n_elem;
*/
     int n_coef;
     int n_layers;

     int dimens[MAX_DIMENS];

     int n_coef2;

     int n_coef_layer;

     double **coef;
     double ***coef1;
     double ****coef2;

     n = dimens_parse(input, dimens);
/*
     n_elem = xrtm_get_options(input->d) & XRTM_OPTION_VECTOR ? 6 : 1;
*/
     n_coef = xrtm_get_max_coef(input->d);

     if (n == 2 && dimens[0] > -1 && dimens[1] > -1) {
          parse_equal(input);

          x = parse_scaler_string(input);

          n_coef_layer = xrtm_get_n_coef(input->d, dimens[0]);

          if ((n_coef2 = load_scat_coefs2(x, n_coef_layer, n_coef, &coef, &flag)) < 0) {
               fprintf(stderr, "ERROR: load_scat_coefs2(), layer %d, deriv %d\n", dimens[0], dimens[1]);
               exit(1);
          }

          if (n_coef2 != n_coef_layer) {
/*
               input_error(&input->locus, "# of phase function expansion coefficient derivatives != to that for corresponding base values");
               exit(1);
*/
          }

          if (xrtm_set_coef_x_11(d, dimens[0], dimens[1], coef)) {
               input_error(&input->locus, "xrtm_set_coef_%s_11()", type);
               exit(1);
          }

          free(input->md->coef_files_l[dimens[0]][dimens[1]]);
          input->md->coef_files_l[dimens[0]][dimens[1]] = strdup(x);

          free(x);
          free_array2_d(coef);
     }
     else
     if ( n == 2 && dimens[0] == -1 && dimens[1] > -1) {
          parse_equal(input);

          n_layers = xrtm_get_n_derivs(input->d);

          x1 = parse_array1_string(input, n_layers);

          coef1 = (double ***) alloc_array1(n_layers, sizeof(double **));

          for (i = 0; i < n_layers; ++i) {
               n_coef_layer = xrtm_get_n_coef(input->d, i);

               if ((n_coef2 = load_scat_coefs2(x1[i], n_coef_layer, n_coef, &coef1[i], &flag)) < 0) {
                    fprintf(stderr, "ERROR: load_scat_coefs2(), layer %d, deriv %d\n", i, dimens[1]);
                    exit(1);
               }

               if (n_coef2 != n_coef_layer) {
/*
                    input_error(&input->locus, "# of phase function expansion coefficient derivatives != to that for corresponding base values");
                    exit(1);
*/
               }
          }

          if (xrtm_set_coef_x_n1(d, dimens[1], coef1)) {
               input_error(&input->locus, "xrtm_set_coef_%s_n1()", type);
               exit(1);
          }

          for (i = 0; i < n_layers; ++i) {
               free(input->md->coef_files_l[i][dimens[1]]);
               input->md->coef_files_l[i][dimens[1]] = strdup(x1[i]);
          }

          free_array1_string(x1, n_layers);

          for (i = 0; i < n_layers; ++i)
               free_array2_d(coef1[i]);
          free_array1(coef1);
     }
     else
     if ((n == 1 && dimens[0] > -1) || (n == 2 && dimens[0] > -1 && dimens[1] == -1)) {
          parse_equal(input);

          x1 = parse_array1_string(input, n_derivs);

          coef1 = (double ***) alloc_array1(n_derivs, sizeof(double **));

          n_coef_layer = xrtm_get_n_coef(input->d, dimens[0]);

          for (i = 0; i < n_derivs; ++i) {
               if ((n_coef2 = load_scat_coefs2(x1[i], n_coef_layer, n_coef, &coef1[i], &flag)) < 0) {
                    fprintf(stderr, "ERROR: load_scat_coefs2(), layer %d, deriv %d\n", dimens[0], i);
                    exit(1);
               }

               if (n_coef2 != n_coef_layer) {
/*
                    input_error(&input->locus, "# of phase function expansion coefficient derivatives != to that for corresponding base values");
                    exit(1);
*/
               }
          }

          if (xrtm_set_coef_x_1n(d, dimens[0], coef1)) {
               input_error(&input->locus, "xrtm_set_coef_%s_1n()", type);
               exit(1);
          }

          for (i = 0; i < n_derivs; ++i) {
               free(input->md->coef_files_l[dimens[0]][i]);
               input->md->coef_files_l[dimens[0]][i] = strdup(x1[i]);
          }

          free_array1_string(x1, n_derivs);

          for (i = 0; i < n_derivs; ++i)
               free_array2_d(coef1[i]);
          free_array1(coef1);
     }
     else
     if ((n == 0) || (n == 1 && dimens[0] == -1 && dimens[1] == -1) || (n == 1 && dimens[0] == -1) || (n == 2 && dimens[0] == -1 && dimens[1] == -1)) {
          parse_equal(input);

          n_layers = xrtm_get_n_layers(input->d);

          x2 = parse_array2_string(input, n_layers, n_derivs);

          coef2 = (double ****) alloc_array2(n_layers, n_derivs, sizeof(double **));

          for (i = 0; i < n_layers; ++i) {
               n_coef_layer = xrtm_get_n_coef(input->d, i);

               for (j = 0; j < n_derivs; ++j) {
                    if ((n_coef2 = load_scat_coefs2(x2[i][j], n_coef_layer, n_coef, &coef2[i][j], &flag)) < 0) {
                         fprintf(stderr, "ERROR: load_scat_coefs2(), layer %d, deriv %d\n", i, j);
                         exit(1);
                    }

                    if (n_coef2 != n_coef_layer) {
/*
                         input_error(&input->locus, "# of phase function expansion coefficient derivatives != to that for corresponding base values");
                         exit(1);
*/
                    }
               }
          }

          if (xrtm_set_coef_x_nn(d, coef2)) {
               input_error(&input->locus, "xrtm_set_coef_%s_nn()", type);
               exit(1);
          }

          for (i = 0; i < n_layers; ++i) {
               for (j = 0; j < n_derivs; ++j) {
                    free(input->md->coef_files_l[i][j]);
                    input->md->coef_files_l[i][j] = strdup(x2[i][j]);
               }
          }

          free_array2_string(x2, n_layers, n_derivs);

          for (i = 0; i < n_layers; ++i) {
               for (j = 0; j < n_derivs; ++j) {
                    free_array2_d(coef2[i][j]);
               }
          }
          free_array2((void **) coef2);
     }
     else {
          input_error(&input->locus, "invalid dimensions");
          exit(1);
     }
}



static void parse_kernel_ampfac(input_data *input) {

     int n;

     int dimens[MAX_DIMENS];

     double x;

     n = dimens_parse(input, dimens);

     if (n == 1 && dimens[0] > -1) {
          parse_equal(input);
          x = parse_scaler_double(input);
          if (xrtm_set_kernel_ampfac(input->d, dimens[0], x)) {
               input_error(&input->locus, "xrtm_set_kernel_ampfac()");
               exit(1);
          }
     }
     else {
          input_error(&input->locus, "invalid dimensions");
          exit(1);
     }
}



static void parse_kernel_ampfac_l(input_data *input, void *d, int n_derivs, const char *type, int (*xrtm_set_kernel_ampfac_x_1)(void *, int, int, double), int (*xrtm_set_kernel_ampfac_x_n)(void *, int, double *)) {

     int n;

     int dimens[MAX_DIMENS];

     double x;
     double *x1;

     n = dimens_parse(input, dimens);

     if (n == 2 && dimens[0] > -1 && dimens[1] > -1) {
          parse_equal(input);
          x = parse_scaler_double(input);
          if (xrtm_set_kernel_ampfac_x_1(d, dimens[0], dimens[1], x)) {
               input_error(&input->locus, "xrtm_set_kernel_ampfac_%s_1()", type);
               exit(1);
          }
     }
     else
     if ((n == 1 && dimens[0] > -1) || (n == 2 && dimens[0] > -1 && dimens[1] == -1)) {
          parse_equal(input);
          x1 = parse_array1_double(input, n_derivs);
          if (xrtm_set_kernel_ampfac_x_n(d, dimens[0], x1)) {
               input_error(&input->locus, "xrtm_set_kernel_ampfac_%s_n()", type);
               exit(1);
          }
          free_array1_double(x1, n_derivs);
     }
     else {
          input_error(&input->locus, "invalid dimensions");
          exit(1);
     }
}



static void parse_kernel_params(input_data *input) {

     int n;

     int n_params;

     int dimens[MAX_DIMENS];

     double x;
     double *x1;

     n = dimens_parse(input, dimens);

     if (n <= 0 || dimens[0] <= -1) {
          input_error(&input->locus, "invalid dimensions");
          exit(1);
     }

     if (n == 2 && dimens[1] > -1) {
          parse_equal(input);
          x = parse_scaler_double(input);
          if (xrtm_set_kernel_params_1(input->d, dimens[0], dimens[1], x)) {
               input_error(&input->locus, "xrtm_set_kernel_params_1()");
               exit(1);
          }
     }
     else
     if (n == 1 || (n == 2 && dimens[1] == -1)) {
          parse_equal(input);
          n_params = kernel_n_params((enum xrtm_kernel_type) xrtm_get_kernel(input->d, dimens[0]));
          if (n_params == 0) {
               input_error(&input->locus, "kernel %d requires no parameters\n", dimens[0]);
               exit(1);
          }
          x1 = parse_array1_double(input, n_params);
          if (xrtm_set_kernel_params_n(input->d, dimens[0], x1)) {
               input_error(&input->locus, "xrtm_set_kernel_params_n()");
               exit(1);
          }
          free_array1_double(x1, n_params);
     }
     else {
          input_error(&input->locus, "invalid dimensions");
          exit(1);
     }
}



static void parse_kernel_params_l(input_data *input, void *d, int n_derivs, const char *type, int (*xrtm_set_kernel_params_x_11)(void *, int, int, int, double), int (*xrtm_set_kernel_params_x_n1)(void *, int, int, double *), int (*xrtm_set_kernel_params_x_1n)(void *, int, int, double *), int (*xrtm_set_kernel_params_x_nn)(void *, int, double **)) {

     int n;

     int n_params;

     int dimens[MAX_DIMENS];

     double x;
     double *x1;
     double **x2;

     n = dimens_parse(input, dimens);

     if (n == 3 && dimens[1] > -1 && dimens[2] > -1) {
          parse_equal(input);
          x = parse_scaler_double(input);
          if (xrtm_set_kernel_params_x_11(d, dimens[0], dimens[1], dimens[2], x)) {
               input_error(&input->locus, "xrtm_set_kernel_params_%s_11()", type);
               exit(1);
          }
     }
     else
     if ( n == 3 && dimens[1] == -1 && dimens[2] > -1) {
          parse_equal(input);
          x1 = parse_array1_double(input, n_derivs);
          if (xrtm_set_kernel_params_x_n1(d, dimens[0], dimens[2], x1)) {
               input_error(&input->locus, "xrtm_set_kernel_params_%s_n1()", type);
               exit(1);
          }
          free_array1_double(x1, n_derivs);
     }
     else
     if ((n == 2 && dimens[1] > -1) || (n == 3 && dimens[1] > -1 && dimens[2] == -1)) {
          parse_equal(input);
          n_params = kernel_n_params((enum xrtm_kernel_type) xrtm_get_kernel(input->d, dimens[0]));
          x1 = parse_array1_double(input, n_params);
          if (xrtm_set_kernel_params_x_1n(d, dimens[0], dimens[1], x1)) {
               input_error(&input->locus, "xrtm_set_kernel_params_%s_1n()", type);
               exit(1);
          }
          free_array1_double(x1, n_params);
     }
     else
     if ((n == 1) || (n == 2 && dimens[1] == -1 && dimens[2] == -1) || (n == 2 && dimens[1]) || (n == 3 && dimens[1] == -1 && dimens[2] == -1)) {
          parse_equal(input);
          n_params = kernel_n_params((enum xrtm_kernel_type) xrtm_get_kernel(input->d, dimens[0]));
          x2 = parse_array2_double(input, n_derivs, n_params);
          if (xrtm_set_kernel_params_x_nn(d, dimens[0], x2)) {
               input_error(&input->locus, "xrtm_set_kernel_params_%s_nn()", type);
               exit(1);
          }
          free_array2_double(x2, n_derivs, n_params);
     }
     else {
          input_error(&input->locus, "invalid dimensions");
          exit(1);
     }
}



/*******************************************************************************
 *
 ******************************************************************************/
static int print_lhs(FILE *fp, format_data *d, const char *name) {

     int i = 0;

     char format[sizeof("-%-18s = ")];

     if (d->use_dash)
          format[i++] = '-';

     format[i++] = '%';

     if (d->use_aligned) {
          format[i++] = '-';
          format[i++] = '1';
          format[i++] = '8';
     }

     format[i++] = 's';
     format[i++] = ' ';

     if (d->use_equal) {
          format[i++] = '=';
          format[i++] = ' ';
     }

     format[i++] = '\0';

     return fprintf(fp, format, name);
}



static int print_lhs_dimens(FILE *fp, format_data *d, const char *name, int n_dimens, int *dimens) {

     int i;
     int n;
     int r;

     char *temp;

     int length;

     length = strlen(name) + n_dimens * (1 + 20 + 1) + 1;

     temp = (char *) malloc(length);

     n = 0;

     n = snprintf(temp + n, length - n, "%s", name);

     for (i = 0; i < n_dimens; ++i)
          n += snprintf(temp + n, length - n, "[%d]", dimens[i]);

     r = print_lhs(fp, d, temp);

     free(temp);

     return r;
}



static int print_array_param_int(FILE *fp, format_data *f, xrtm_data *d, const char *name, int n, int (*xrtm_get_x)(xrtm_data *, int *)) {

     int *i1;

     i1 = alloc_array1_i(n);

     if (xrtm_get_x(d, i1) == XRTM_INT_ERROR) {
          fprintf(stderr, "ERROR: xrtm_get_x()\n");
          return -1;
     }

     print_lhs(fp, f, name);
     print_array1_int(fp, f, i1, n);
     fprintf(fp, "%c", f->nl);

     free_array1_i(i1);

     return 0;
}



static int print_array_param_double(FILE *fp, format_data *f, xrtm_data *d, const char *name, int n, int (*xrtm_get_x)(xrtm_data *, double *)) {

     double *x1;

     x1 = alloc_array1_d(n);

     if (xrtm_get_x(d, x1) == XRTM_DBL_ERROR) {
          fprintf(stderr, "ERROR: xrtm_get_x()\n");
          return -1;
     }

     print_lhs(fp, f, name);
     print_array1_double(fp, f, x1, n);
     fprintf(fp, "%c", f->nl);

     free_array1_d(x1);

     return 0;
}



static int print_layer_param_double(FILE *fp, format_data *f, xrtm_data *d, const char *name, double (*xrtm_get_x)(xrtm_data *, int)) {

     int i;

     int n_layers;

     double *x1;

     n_layers = xrtm_get_n_layers(d);

     x1 = alloc_array1_d(n_layers);
     for (i = 0; i < n_layers; ++i) {
          if ((x1[i] = xrtm_get_x(d, i)) == XRTM_DBL_ERROR) {
               fprintf(stderr, "ERROR: xrtm_get_%s(), i_layer = %d\n", name, i);
               return -1;
          }
     }

     print_lhs(fp, f, name);
     print_array1_double(fp, f, x1, n_layers);
     free_array1_d(x1);
     fprintf(fp, "%c", f->nl);

     return 0;
}



static int print_layer_param_l_double(FILE *fp, format_data *f, xrtm_data *d, const char *name, double (*xrtm_get_x_l)(xrtm_data *, int, int)) {

     int i;
     int j;

     int n_layers;
     int n_derivs;

     double **x2;

     n_layers = xrtm_get_n_layers(d);
     n_derivs = xrtm_get_n_derivs(d);

     x2 = alloc_array2_d(n_layers, n_derivs);
     for (i = 0; i < n_layers; ++i) {
          for (j = 0; j < n_derivs; ++j) {
               if ((x2[i][j] = xrtm_get_x_l(d, i, j)) == XRTM_DBL_ERROR) {
                    fprintf(stderr, "ERROR: xrtm_get_%s(), i_layer = %d, i_deriv = %d\n", name, i, j);
                    return -1;
               }
          }
     }

     print_lhs(fp, f, name);
     print_array2_double(fp, f, x2, n_layers, n_derivs);
     free_array2_d(x2);
     fprintf(fp, "%c", f->nl);

     return 0;
}



static int print_coef(FILE *fp, format_data *f, xrtm_data *d, misc_data *md) {

     int i;
     int j;
     int k;

     int n_elem;
     int n_coef;
     int n_coef2;
     int n_layers;

     double **x2;

     n_elem   = xrtm_get_options(d) & XRTM_OPTION_VECTOR ? 6 : 1;
     n_coef   = xrtm_get_max_coef(d);
     n_layers = xrtm_get_n_layers(d);

     x2 = alloc_array2_d(n_coef, n_elem);

     for (i = 0; i < n_layers; ++i) {
          if (! md || ! md->coef_files[i]) {
               n_coef2 = xrtm_get_n_coef(d, i);

               for (j = 0; j < n_elem; ++j) {
                    for (k = 0; k < n_coef2; ++k) {
                         if ((x2[k][j] = xrtm_get_coef(d, i, j, k)) == XRTM_DBL_ERROR) {
                              fprintf(stderr, "ERROR: xrtm_get_coef(), i_layer = %d, i_elem = %d, i_coef = %d\n", i, j, k);
                              return -1;
                         }
                    }
               }

               print_lhs_dimens(fp, f, "coef", 1, &i);
               fprintf(fp, "%d, ", n_coef2);
               print_array2_double(fp, f, x2, n_coef2, n_elem);
               fprintf(fp, "%c", f->nl);
          }
          else {
               print_lhs_dimens(fp, f, "coef_files", 1, &i);
               fprintf(fp, "\'%s\'", md->coef_files[i]);
               fprintf(fp, "%c", f->nl);
          }
     }

     free_array2_d(x2);

     return 0;
}



static int print_coef_l(FILE *fp, format_data *f, xrtm_data *d, misc_data *md) {

     int i;
     int j;
     int k;
     int l;

     int n_elem;
     int n_coef;
     int n_coef2;
     int n_derivs;
     int n_layers;

     double ***x3;

     n_elem   = xrtm_get_options(d) & XRTM_OPTION_VECTOR ? 6 : 1;
     n_coef   = xrtm_get_max_coef(d);
     n_derivs = xrtm_get_n_derivs(d);
     n_layers = xrtm_get_n_layers(d);

     x3 = alloc_array3_d(n_derivs, n_coef, n_elem);

     for (i = 0; i < n_layers; ++i) {
          if (! md || ! md->coef_files[i]) {
               n_coef2 = xrtm_get_n_coef(d, i);

               for (j = 0; j < n_derivs; ++j) {
                    for (k = 0; k < n_elem; ++k) {
                         for (l = 0; l < n_coef2; ++l) {
                              if ((x3[j][l][k] = xrtm_get_coef_l(d, i, j, k, l)) == XRTM_DBL_ERROR) {
                                   fprintf(stderr, "ERROR: xrtm_get_coef_l(), i_layer = %d, i_deriv = %d, i_elem = %d, i_coef = %d\n", i, j, k, l);
                                   return -1;
                              }
                         }
                    }
               }

               print_lhs_dimens(fp, f, "coef_l", 1, &i);
               print_array3_double(fp, f, x3, n_derivs, n_coef2, n_elem);
               fprintf(fp, "%c", f->nl);
          }
          else {
               print_lhs_dimens(fp, f, "coef_files_l", 1, &i);
               print_array1_string(fp, f, md->coef_files_l[i], n_derivs);
               fprintf(fp, "%c", f->nl);
          }
     }

     free_array3_d(x3);

     return 0;
}



static int print_kernel_ampfac_l(FILE *fp, format_data *f, xrtm_data *d, int i_kernel) {

     int i;

     int n_derivs;

     double *x1;

     n_derivs = xrtm_get_n_derivs(d);

     x1 = alloc_array1_d(n_derivs);
     for (i = 0; i < n_derivs; ++i) {
          if ((x1[i] = xrtm_get_kernel_ampfac_l(d, i_kernel, i)) == XRTM_DBL_ERROR) {
               fprintf(stderr, "ERROR: xrtm_get_kernel_ampfac_l(), i_kernel = %d, i_deriv = %d\n", i_kernel, i);
               return -1;
          }
     }

     print_lhs_dimens(fp, f, "kernel_ampfac_l", 1, &i_kernel);
     print_array1_double(fp, f, x1, n_derivs);
     free_array1_d(x1);
     fprintf(fp, "%c", f->nl);

     return 0;
}



static int print_kernel_params(FILE *fp, format_data *f, xrtm_data *d, int i_kernel) {

     int i;

     int n_params;

     double *x1;

     n_params = kernel_n_params((enum xrtm_kernel_type) xrtm_get_kernel(d, i_kernel));

     x1 = alloc_array1_d(n_params);
     for (i = 0; i < n_params; ++i) {
          if ((x1[i] = xrtm_get_kernel_params(d, i_kernel, i)) == XRTM_DBL_ERROR) {
               fprintf(stderr, "ERROR: xrtm_get_kernel_params(), i_kernel = %d, i_param = %d\n", i_kernel, i);
               return -1;
          }
     }

     print_lhs_dimens(fp, f, "kernel_params", 1, &i_kernel);
     print_array1_double(fp, f, x1, n_params);
     free_array1_d(x1);
     fprintf(fp, "%c", f->nl);

     return 0;
}



static int print_kernel_params_l(FILE *fp, format_data *f, xrtm_data *d, int i_kernel) {

     int i;
     int j;

     int n_params;
     int n_derivs;

     double **x2;

     n_params = kernel_n_params((enum xrtm_kernel_type) xrtm_get_kernel(d, i_kernel));
     n_derivs = xrtm_get_n_derivs(d);

     x2 = alloc_array2_d(n_derivs, n_params);
     for (i = 0; i < n_derivs; ++i) {
          for (j = 0; j < n_params; ++j) {
               if ((x2[i][j] = xrtm_get_kernel_params_l(d, i_kernel, i, j)) == XRTM_DBL_ERROR) {
                    fprintf(stderr, "ERROR: xrtm_get_kernel_params(), i_kernel = %d, i_param = %d, i_deriv = %d\n", i_kernel, i, j);
                    return -1;
               }
          }
     }

     print_lhs_dimens(fp, f, "kernel_params_l", 1, &i_kernel);
     print_array2_double(fp, f, x2, n_derivs, n_params);
     free_array2_d(x2);
     fprintf(fp, "%c", f->nl);

     return 0;
}



/*******************************************************************************
 *
 ******************************************************************************/
static void input_parse1(input_data *input) {

     int r;

     long options      = 0;
     long solvers      = 0;
     int n_coef        = -999;
     int n_quad        = -999;
     int n_stokes      = -999;
     int n_derivs      = -999;
     int n_layers      = -999;
     int n_kernels     = -999;
     int n_kernel_quad = -999;
     int *kernels      = NULL;
     int n_out_levels  = -999;
     int n_out_thetas  = -999;

     int n_derivs2;
     int n_derivs3;

     int flag = 0;

     int fd_method = -1;

     fdclist *list;

     lex_type_data lex_type;

     while ((r = yylex(&input->locus, &lex_type))) {
          switch(r) {
               case INPUT_OPTION_FD:
                    parse_equal(input);
                    fd_method = parse_scaler_int(input);
                    flag = 1;
                    break;
               case XRTM_INPUT_OPTIONS:
                    parse_equal(input);
                    list = parse_list(input, INPUT_TYPE_STRING, 0, NULL);
                    options = xrtm_string_list_to_options_mask(input, list);
                    if (options < 0) {
                         input_error(&input->locus, "invalid option string");
                         exit(1);
                    }
                    fdclist_free(list);
                    break;
               case XRTM_INPUT_SOLVERS:
                    parse_equal(input);
                    list = parse_list(input, INPUT_TYPE_STRING, 0, NULL);
                    solvers = xrtm_string_list_to_solvers_mask(input, list);
                    if (solvers < 0) {
                         input_error(&input->locus, "invalid solver string");
                         exit(1);
                    }
                    fdclist_free(list);
                    break;
               case XRTM_INPUT_MAX_COEF:
                    parse_equal(input);
                    n_coef = parse_scaler_int(input);
                    break;
               case XRTM_INPUT_N_QUAD:
                    parse_equal(input);
                    n_quad = parse_scaler_int(input);
                    break;
               case XRTM_INPUT_N_STOKES:
                    parse_equal(input);
                    n_stokes = parse_scaler_int(input);
                    break;
               case XRTM_INPUT_N_DERIVS:
                    parse_equal(input);
                    n_derivs = parse_scaler_int(input);
                    break;
               case XRTM_INPUT_N_LAYERS:
                    parse_equal(input);
                    n_layers = parse_scaler_int(input);
                    break;
               case XRTM_INPUT_N_KERNELS:
                    parse_equal(input);
                    n_kernels = parse_scaler_int(input);
                    break;
               case XRTM_INPUT_N_KERNEL_QUAD:
                    parse_equal(input);
                    n_kernel_quad = parse_scaler_int(input);
                    break;
               case XRTM_INPUT_KERNELS:
                    parse_equal(input);
                    list = parse_list(input, INPUT_TYPE_STRING, 0, NULL);
                    kernels = xrtm_string_list_to_kernels_vector(input, list);
                    if (kernels == NULL) {
                         input_error(&input->locus, "invalid kernel string");
                         exit(1);
                    }
                    fdclist_free(list);
                    break;
               case XRTM_INPUT_N_OUT_LEVELS:
                    parse_equal(input);
                    n_out_levels = parse_scaler_int(input);
                    break;
               case XRTM_INPUT_N_OUT_THETAS:
                    parse_equal(input);
                    n_out_thetas = parse_scaler_int(input);
                    break;
               default:
                    if (r >= XRTM_INPUT_DOUB_D_TAU && r <= XRTM_INPUT_KERNEL_PARAMS_L) {
                         yyrewind(r, &lex_type);
                         goto L1;
                    }
                    else {
                         input_error(&input->locus, "invalid input parameter: %s", get_yytext());
                         exit(1);
                    }
          }
     }

L1:  if (flag && options & XRTM_OPTION_CALC_DERIVS) {
          fprintf(stderr, "ERROR: cannot use \"fd\" and xrtm option \"calc_derivs\" together\n");
          exit(1);
     }

     if (! flag)
          n_derivs2 = n_derivs;
     else {
          n_derivs2 = 0;
          n_derivs3 = n_derivs;
     }

     if (xrtm_create(input->d, (enum xrtm_option_mask) options, (enum xrtm_solver_mask) solvers, n_coef, n_quad, n_stokes, n_derivs2, n_layers, 1, n_kernels, n_kernel_quad, (enum xrtm_kernel_type *) kernels, n_out_levels, n_out_thetas)) {
          fprintf(stderr, "ERROR: invalid model creation parameters\n");
          exit(1);
     }

     if (flag) {
          if (xrtm_fd_create(input->fd, input->d, n_derivs3)) {
               fprintf(stderr, "ERROR: invalid fd model creation parameters\n");
               exit(1);
          }
     }

     misc_init(input->d, input->md, fd_method);

     free_array1_i(kernels);

     return;
}



static void input_parse2(input_data *input) {

     int m;
     int n;
     int r;

     int n_derivs;

     lex_type_data lex_type;

     if (input->md->fd_method < 0)
          n_derivs = xrtm_get_n_derivs(input->d);
     else
          n_derivs = xrtm_fd_get_n_derivs(input->fd);

     while ((r = yylex(&input->locus, &lex_type))) {
          switch(r) {
               case XRTM_INPUT_DOUB_D_TAU:
                    parse_scaler_param_double(input, "doub_d_tau", xrtm_set_doub_d_tau);
                    break;
               case XRTM_INPUT_PADE_PARAMS:
                    parse_scaler_param_int2(input, "pade_params", xrtm_set_pade_params);
                    break;
               case XRTM_INPUT_SOS_PARAMS:
                    parse_sos_params(input);
                    break;
               case XRTM_INPUT_FOURIER_TOL:
                    parse_scaler_param_double(input, "fourier_tol", xrtm_set_fourier_tol);
                    break;
               case XRTM_INPUT_LAMBDA:
                    parse_scaler_param_double(input, "lambda", xrtm_set_lambda);
                    break;
               case XRTM_INPUT_PLANET_R:
                    parse_scaler_param_double(input, "planet_r", xrtm_set_planet_r);
                    break;
               case XRTM_INPUT_LEVELS_Z:
                    n = xrtm_get_n_layers(input->d) + 1;
                    parse_array_param_double(input, input->d, "levels_z", n, (int (*)(void *, void *)) xrtm_set_levels_z);
                    break;
               case XRTM_INPUT_OUT_LEVELS:
                    n = xrtm_get_n_out_levels(input->d);
                    parse_array_param_int(input, "out_levels", n, (int (*)(xrtm_data *, void *)) xrtm_set_out_levels);
                    break;
               case XRTM_INPUT_OUT_TAUS:
                    n = xrtm_get_n_out_levels(input->d);
                    parse_array_param_double(input, input->d, "out_taus", n, (int (*)(void *, void *)) xrtm_set_out_taus);
                    break;
               case XRTM_INPUT_OUT_THETAS:
                    n = xrtm_get_n_out_thetas(input->d);
                    parse_array_param_double(input, input->d, "out_thetas", n, (int (*)(void *, void *)) xrtm_set_out_thetas);
                    break;
               case XRTM_INPUT_F_0:
                    parse_scaler_param_double(input, "F_0", xrtm_set_F_0);
                    break;
               case XRTM_INPUT_THETA_0:
                    parse_scaler_param_double(input, "theta_0", xrtm_set_theta_0);
                    break;
               case XRTM_INPUT_PHI_0:
                    parse_scaler_param_double(input, "phi_0", xrtm_set_phi_0);
                    break;
               case XRTM_INPUT_F_ISO_TOP:
                    parse_scaler_param_double(input, "F_iso_top", xrtm_set_F_iso_top);
                    break;
               case XRTM_INPUT_F_ISO_TOP_L:
                    parse_array_param_double(input, input->d, "F_iso_top_l", n_derivs, (int (*)(void *, void *)) xrtm_set_F_iso_top_l_n);
                    break;
               case XRTM_INPUT_F_ISO_TOP_P:
                    parse_array_param_double(input, input->fd, "F_iso_top_p", n_derivs, (int (*)(void *, void *)) xrtm_fd_set_F_iso_top_p_n);
                    break;
               case XRTM_INPUT_F_ISO_BOT:
                    parse_scaler_param_double(input, "F_iso_bot", xrtm_set_F_iso_bot);
                    break;
               case XRTM_INPUT_F_ISO_BOT_L:
                    parse_array_param_double(input, input->d, "F_iso_bot_l", n_derivs, (int (*)(void *, void *)) xrtm_set_F_iso_bot_l_n);
                    break;
               case XRTM_INPUT_F_ISO_BOT_P:
                    parse_array_param_double(input, input->fd, "F_iso_bot_p", n_derivs, (int (*)(void *, void *)) xrtm_fd_set_F_iso_bot_p_n);
                    break;
               case XRTM_INPUT_LEVELS_B:
                    n = xrtm_get_n_layers(input->d) + 1;
                    parse_array_param_double(input, input->d, "levels_b", n, (int (*)(void *, void *)) xrtm_set_levels_b);
                    break;
               case XRTM_INPUT_LEVELS_B_L:
                    m = xrtm_get_n_layers(input->d) + 1;
                    parse_array_param_double2(input, input->d, "levels_b_l", m, n_derivs, (int (*)(void *, void *)) xrtm_set_levels_b_l_n);
                    break;
               case XRTM_INPUT_LEVELS_B_P:
                    m = xrtm_get_n_layers(input->d) + 1;
                    parse_array_param_double2(input, input->fd, "levels_b_p", m, n_derivs, (int (*)(void *, void *)) xrtm_fd_set_levels_b_p_n);
                    break;
               case XRTM_INPUT_SURFACE_B:
                    parse_scaler_param_double(input, "surface_b", xrtm_set_surface_b);
                    break;
               case XRTM_INPUT_SURFACE_B_L:
                    parse_array_param_double(input, input->d, "surface_b_l_n", n_derivs, (int (*)(void *, void *)) xrtm_set_surface_b_l_n);
                    break;
               case XRTM_INPUT_SURFACE_B_P:
                    parse_array_param_double(input, input->fd, "surface_b_p_n", n_derivs, (int (*)(void *, void *)) xrtm_fd_set_surface_b_p_n);
                    break;
               case XRTM_INPUT_G:
                    parse_layer_param_double(input, "g", xrtm_set_g_1, xrtm_set_g_n);
                    break;
               case XRTM_INPUT_G_L:
                    parse_layer_param_l_double(input, input->d, n_derivs, "g", "l", (int (*)(void *, int, int, double)) xrtm_set_g_l_11, (int (*)(void *, int, double *)) xrtm_set_g_l_n1, (int (*)(void *, int, double *)) xrtm_set_g_l_1n, (int (*)(void *, double **)) xrtm_set_g_l_nn);
                    break;
               case XRTM_INPUT_G_P:
                    parse_layer_param_l_double(input, input->fd, n_derivs, "g", "p", (int (*)(void *, int, int, double)) xrtm_fd_set_g_p_11, (int (*)(void *, int, double *)) xrtm_fd_set_g_p_n1, (int (*)(void *, int, double *)) xrtm_fd_set_g_p_1n, (int (*)(void *, double **)) xrtm_fd_set_g_p_nn);
                    break;
               case XRTM_INPUT_CHI:
                    parse_coef(input);
                    break;
               case XRTM_INPUT_COEF_L:
                    parse_coef_l(input, input->d, n_derivs, "l", (int (*)(void *, int, int, double **)) xrtm_set_coef_l_11, (int (*)(void *, int, double ***)) xrtm_set_coef_l_n1, (int (*)(void *, int, double ***)) xrtm_set_coef_l_1n, (int (*)(void *, double ****)) xrtm_set_coef_l_nn);
                    break;
               case XRTM_INPUT_COEF_P:
                    parse_coef_l(input, input->fd, n_derivs, "p", (int (*)(void *, int, int, double **)) xrtm_fd_set_coef_p_11, (int (*)(void *, int, double ***)) xrtm_fd_set_coef_p_n1, (int (*)(void *, int, double ***)) xrtm_fd_set_coef_p_1n, (int (*)(void *, double ****)) xrtm_fd_set_coef_p_nn);
                    break;
               case XRTM_INPUT_COEF_FILES:
                    parse_coef_files(input);
                    break;
               case XRTM_INPUT_COEF_FILES_L:
                    parse_coef_files_l(input, input->d, n_derivs, "l", (int (*)(void *, int, int, double **)) xrtm_set_coef_l_11, (int (*)(void *, int, double ***)) xrtm_set_coef_l_n1, (int (*)(void *, int, double ***)) xrtm_set_coef_l_1n, (int (*)(void *, double ****)) xrtm_set_coef_l_nn);
                    break;
               case XRTM_INPUT_COEF_FILES_P:
                    parse_coef_files_l(input, input->fd, n_derivs, "p", (int (*)(void *, int, int, double **)) xrtm_fd_set_coef_p_11, (int (*)(void *, int, double ***)) xrtm_fd_set_coef_p_n1, (int (*)(void *, int, double ***)) xrtm_fd_set_coef_p_1n, (int (*)(void *, double ****)) xrtm_fd_set_coef_p_nn);
                    break;
               case XRTM_INPUT_OMEGA:
                    parse_layer_param_double(input, "omega", xrtm_set_omega_1, xrtm_set_omega_n);
                    break;
               case XRTM_INPUT_OMEGA_L:
                    parse_layer_param_l_double(input, input->d, n_derivs, "omega", "l", (int (*)(void *, int, int, double)) xrtm_set_omega_l_11, (int (*)(void *, int, double *)) xrtm_set_omega_l_n1, (int (*)(void *, int, double *)) xrtm_set_omega_l_1n, (int (*)(void *, double **)) xrtm_set_omega_l_nn);
                    break;
               case XRTM_INPUT_OMEGA_P:
                    parse_layer_param_l_double(input, input->fd, n_derivs, "omega", "p", (int (*)(void *, int, int, double)) xrtm_fd_set_omega_p_11, (int (*)(void *, int, double *)) xrtm_fd_set_omega_p_n1, (int (*)(void *, int, double *)) xrtm_fd_set_omega_p_1n, (int (*)(void *, double **)) xrtm_fd_set_omega_p_nn);
                    break;
               case XRTM_INPUT_LTAU:
                    parse_layer_param_double(input, "ltau", xrtm_set_ltau_1, xrtm_set_ltau_n);
                    break;
               case XRTM_INPUT_LTAU_L:
                    parse_layer_param_l_double(input, input->d, n_derivs, "ltau", "l", (int (*)(void *, int, int, double)) xrtm_set_ltau_l_11, (int (*)(void *, int, double *)) xrtm_set_ltau_l_n1, (int (*)(void *, int, double *)) xrtm_set_ltau_l_1n, (int (*)(void *, double **)) xrtm_set_ltau_l_nn);
                    break;
               case XRTM_INPUT_LTAU_P:
                    parse_layer_param_l_double(input, input->fd, n_derivs, "ltau", "p", (int (*)(void *, int, int, double)) xrtm_fd_set_ltau_p_11, (int (*)(void *, int, double *)) xrtm_fd_set_ltau_p_n1, (int (*)(void *, int, double *)) xrtm_fd_set_ltau_p_1n, (int (*)(void *, double **)) xrtm_fd_set_ltau_p_nn);
                    break;
               case XRTM_INPUT_KERNEL_AMPFAC:
                    parse_kernel_ampfac(input);
                    break;
               case XRTM_INPUT_KERNEL_AMPFAC_L:
                    parse_kernel_ampfac_l(input, input->d, n_derivs, "l", (int (*)(void *, int, int, double)) xrtm_set_kernel_ampfac_l_1, (int (*)(void *, int, double *)) xrtm_set_kernel_ampfac_l_n);
                    break;
               case XRTM_INPUT_KERNEL_AMPFAC_P:
                    parse_kernel_ampfac_l(input, input->fd, n_derivs, "p", (int (*)(void *, int, int, double)) xrtm_fd_set_kernel_ampfac_p_1, (int (*)(void *, int, double *)) xrtm_fd_set_kernel_ampfac_p_n);
                    break;
               case XRTM_INPUT_KERNEL_PARAMS:
                    parse_kernel_params(input);
                    break;
               case XRTM_INPUT_KERNEL_PARAMS_L:
                    parse_kernel_params_l(input, input->d, n_derivs, "l", (int (*)(void *, int, int, int, double)) xrtm_set_kernel_params_l_11, (int (*)(void *, int, int, double *)) xrtm_set_kernel_params_l_n1, (int (*)(void *, int, int, double *)) xrtm_set_kernel_params_l_1n, (int (*)(void *, int, double **)) xrtm_set_kernel_params_l_nn);
                    break;
               case XRTM_INPUT_KERNEL_PARAMS_P:
                    parse_kernel_params_l(input, input->fd, n_derivs, "p", (int (*)(void *, int, int, int, double)) xrtm_fd_set_kernel_params_p_11, (int (*)(void *, int, int, double *)) xrtm_fd_set_kernel_params_p_n1, (int (*)(void *, int, int, double *)) xrtm_fd_set_kernel_params_p_1n, (int (*)(void *, int, double **)) xrtm_fd_set_kernel_params_p_nn);
                    break;
               case XRTM_INPUT_DELTA:
                    n = xrtm_fd_get_n_derivs(input->fd);
                    parse_array_param_double(input, input->fd, "delta", n, (int (*)(void *, void *)) xrtm_fd_set_delta);
                    break;
               default:
                    input_error(&input->locus, "invalid input parameter: %s", get_yytext());
                    exit(1);
          }
     }
}



static void input_parse(input_data *input) {

     input_parse1(input);
     input_parse2(input);
}



/*******************************************************************************
 *
 ******************************************************************************/
static void input_init(input_data *input, xrtm_data *d, xrtm_fd_data *fd, misc_data *md, int use_aligned, int use_array2, int use_dash, int use_equal, const char *file) {

     input->use_aligned  = use_aligned;
     input->use_array2   = use_array2;
     input->use_dash     = use_dash;
     input->use_equal    = use_equal;

     input->d            = d;
     input->fd           = fd;
     input->md           = md;

     input->locus.file   = file;
     input->locus.line   = 1;
}



static void input_free(input_data *input) {

}



/*******************************************************************************
 *
 ******************************************************************************/
void misc_init(xrtm_data *d, misc_data *md, int fd_method) {

     int i;
     int j;

     md->n_layers = xrtm_get_n_layers(d);
     md->n_derivs = xrtm_get_n_derivs(d);

     md->fd_method = fd_method;

     md->coef_files = (char **) alloc_array1(md->n_layers, sizeof(char *));
     for (i = 0; i < md->n_layers; ++i)
          md->coef_files[i] = NULL;

     md->coef_files_l = NULL;
     if (md->n_derivs > 0 || md->fd_method >= 0) {
          md->coef_files_l = (char ***) alloc_array2(md->n_layers, md->n_derivs, sizeof(char *));
          for (i = 0; i < md->n_layers; ++i) {
               for (j = 0; j < md->n_derivs; ++j) {
                    md->coef_files_l[i][j] = NULL;
               }
          }
     }
}



void misc_free(misc_data *md) {

     int i;
     int j;

     for (i = 0; i <  md->n_layers; ++i)
          free(md->coef_files[i]);
     free_array1(md->coef_files);

     if (md->n_derivs > 0 || md->fd_method >= 0) {
          for (i = 0; i < md->n_layers; ++i) {
               for (j = 0; j < md->n_derivs; ++j) {
                    free(md->coef_files_l[i][j]);
               }
          }
          free_array2((void **) md->coef_files_l);
     }
}



/*******************************************************************************
 *
 ******************************************************************************/
int xrtm_sread_input(xrtm_data *d, xrtm_fd_data *fd, misc_data *md, char *s, int use_aligned, int use_array2, int use_dash, int use_equal) {

     input_data input;

     YY_BUFFER_STATE buffer;

     buffer = yy_scan_string(s);

     input_init(&input, d, fd, md, use_aligned, use_array2, use_dash, use_equal, "string");

     input_parse(&input);

     yy_delete_buffer(buffer);

     input_free(&input);

     return 0;
}



int xrtm_fread_input(xrtm_data *d, xrtm_fd_data *fd, misc_data *md, const char *filename, int use_aligned, int use_array2, int use_dash, int use_equal) {

     const char *file;

     input_data input;

     if (strcmp(filename, "-") == 0) {
          yyin = stdin;

          file = "stdin";
     }
     else {
          if ((yyin = fopen(filename, "r")) == NULL) {
               fprintf(stderr, "ERROR: Error opening file for reading: %s ... %s\n",
                      filename, strerror(errno));
               return -1;
          }

          file = filename;
     }

     input_init(&input, d, fd, md, use_aligned, use_array2, use_dash, use_equal, file);

     input_parse(&input);

     if (yyin != stdin)
          fclose(yyin);

     input_free(&input);

     return 0;
}



int xrtm_swrite_input(xrtm_data *d, xrtm_fd_data *fd, misc_data *md, char *s, char nl, int use_aligned, int use_array2, int use_dash, int use_equal) {

     return 0;
}



int xrtm_fwrite_input_fn(xrtm_data *d, xrtm_fd_data *fd, misc_data *md, const char *filename, char nl, int use_aligned, int use_array2, int use_dash, int use_equal) {

     FILE *fp;

     if (strcmp(filename, "stdout") == 0)
          fp = stdout;
     else
     if (strcmp(filename, "stderr") == 0)
          fp = stderr;
     else {
          if ((fp = fopen(filename, "w")) == NULL) {
               fprintf(stderr, "ERROR: Error opening file for writing: %s ... %s\n",
                      filename, strerror(errno));
               return -1;
          }
     }

     if (xrtm_fwrite_input_fp(d, fd, md, fp, nl, use_aligned, use_array2, use_dash, use_equal)) {
          fprintf(stderr, "ERROR: xrtm_fwrite_input_fp()\n");
          return -1;
     }

     if (fp != stdin && fp != stderr)
          fclose(fp);

     return 0;
}



int xrtm_fwrite_input_fp(xrtm_data *d, xrtm_fd_data *fd, misc_data *md, FILE *fp, char nl, int use_aligned, int use_array2, int use_dash, int use_equal) {

     int i;
     int j;

     int options;
     int solvers;
     int n_coef;
     int n_quad;
     int n_stokes;
     int n_derivs;
     int n_layers;
     int n_kernels;
     int n_kernel_quad;
     int n_out_levels;
     int n_out_thetas;

     int *i1;

     double x;
     double y;

     fdclist *list;

     format_data f;

     f.nl          = nl;
     f.use_aligned = use_aligned;
     f.use_array2  = use_array2;
     f.use_dash    = use_dash;
     f.use_equal   = use_equal;

     options       = xrtm_get_options(d);
     solvers       = xrtm_get_solvers(d);
     n_coef        = xrtm_get_max_coef(d);
     n_quad        = xrtm_get_n_quad(d);
     n_stokes      = xrtm_get_n_stokes(d);
     n_derivs      = xrtm_get_n_derivs(d);
     n_layers      = xrtm_get_n_layers(d);
     n_kernels     = xrtm_get_n_kernels(d);
     n_kernel_quad = xrtm_get_n_kernel_quad(d);
     n_out_levels  = xrtm_get_n_out_levels(d);
     n_out_thetas  = xrtm_get_n_out_thetas(d);

     list = xrtm_options_mask_to_string_list(NULL, options);
     print_lhs(fp, &f, "options");
     print_list(fp, &f, list);
     fprintf(fp, "%c", f.nl);
     fdclist_free(list);

     list = xrtm_solvers_mask_to_string_list(NULL, solvers);
     print_lhs(fp, &f, "solvers");
     print_list(fp, &f, list);
     fprintf(fp, "%c", f.nl);
     fdclist_free(list);

     print_lhs(fp, &f, "max_coef");      fprintf(fp, "%d%c", n_coef,        f.nl);
     print_lhs(fp, &f, "n_quad");        fprintf(fp, "%d%c", n_quad,        f.nl);
     print_lhs(fp, &f, "n_stokes");      fprintf(fp, "%d%c", n_stokes,      f.nl);
     print_lhs(fp, &f, "n_derivs");      fprintf(fp, "%d%c", n_derivs,      f.nl);
     print_lhs(fp, &f, "n_layers");      fprintf(fp, "%d%c", n_layers,      f.nl);
     print_lhs(fp, &f, "n_kernels");     fprintf(fp, "%d%c", n_kernels,     f.nl);
     print_lhs(fp, &f, "n_kernel_quad"); fprintf(fp, "%d%c", n_kernel_quad, f.nl);

     i1 = alloc_array1_i(n_kernels);
     for (i = 0; i < n_kernels; ++i)
          i1[i] = xrtm_get_kernel(d, i);
     list = xrtm_kernels_vector_to_string_list(NULL, (enum xrtm_kernel_type *) i1, n_kernels);
     free_array1_i(i1);
     print_lhs(fp, &f, "kernels");
     print_list(fp, &f, list);
     fprintf(fp, "%c", f.nl);
     fdclist_free(list);

     print_lhs(fp, &f, "n_out_levels"); fprintf(fp, "%d%c", n_out_levels, f.nl);
     print_lhs(fp, &f, "n_out_thetas"); fprintf(fp, "%d%c", n_out_thetas, f.nl);

     if (f.use_dash)
          fprintf(fp, ":: ");

     if (solvers & XRTM_SOLVERS_DOUBLING) {
          if ((x = xrtm_get_doub_d_tau(d)) < 0.) {
               fprintf(stderr, "ERROR: xrtm_get_doub_d_tau()\n");
               return -1;
          }
          print_lhs(fp, &f, "doub_d_tau"); fprintf(fp, "%e%c", x, f.nl);
     }

     if (solvers & XRTM_SOLVER_PADE_ADD) {
          if (xrtm_get_pade_params(d, &i, &j) < 0.) {
               fprintf(stderr, "ERROR: xrtm_get_pade_params()\n");
               return -1;
          }
          print_lhs(fp, &f, "pade_params"); fprintf(fp, "%d, %d%c", i, j, f.nl);
     }

     if (solvers & XRTM_SOLVER_SOS) {
          if (xrtm_get_sos_params(d, &i, &x, &y) < 0.) {
               fprintf(stderr, "ERROR: xrtm_get_sos_params()\n");
               return -1;
          }
          print_lhs(fp, &f, "sos_params"); fprintf(fp, "%d, %e, %e%c", i, x, y, f.nl);
     }

     if ((x = xrtm_get_fourier_tol(d)) < 0.) {
          fprintf(stderr, "ERROR: xrtm_get_fourier_tol()\n");
          return -1;
     }
     print_lhs(fp, &f, "fourier_tol"); fprintf(fp, "%e%c", x, f.nl);

     if (options & XRTM_OPTION_SOURCE_THERMAL) {
          if ((x = xrtm_get_lambda(d)) < 0.) {
               fprintf(stderr, "ERROR: xrtm_get_lambda()\n");
               return -1;
          }
          print_lhs(fp, &f, "lambda"); fprintf(fp, "%e%c", x, f.nl);
     }

     if (options & XRTM_OPTION_PSA) {
          if ((x = xrtm_get_planet_r(d)) < 0.) {
               fprintf(stderr, "ERROR: xrtm_get_planet_r()\n");
               return -1;
          }
          print_lhs(fp, &f, "planet_r"); fprintf(fp, "%e%c", x, f.nl);

          if (print_array_param_double(fp, &f, d, "levels_z", n_layers + 1, xrtm_get_levels_z)) {
               fprintf(stderr, "ERROR: print_array_param_double()\n");
               return -1;
          }
     }

     if (options & XRTM_OPTION_OUTPUT_AT_LEVELS) {
          if (print_array_param_int(fp, &f, d, "out_levels", n_out_levels, xrtm_get_out_levels)) {
               fprintf(stderr, "ERROR: print_array_param_int()\n");
               return -1;
          }
     }

     if (options & XRTM_OPTION_OUTPUT_AT_TAUS) {
          if (print_array_param_double(fp, &f, d, "out_taus", n_out_levels, xrtm_get_out_taus)) {
               fprintf(stderr, "ERROR: print_array_param_double()\n");
               return -1;
          }
     }

     if (n_out_thetas > 0) {
          if (print_array_param_double(fp, &f, d, "out_thetas", n_out_thetas, xrtm_get_out_thetas)) {
               fprintf(stderr, "ERROR: print_array_param_double()\n");
               return -1;
          }
     }

     if ((x = xrtm_get_F_iso_top(d)) < 0.) {
          fprintf(stderr, "ERROR: xrtm_get_F_iso_top()\n");
          return -1;
     }
     print_lhs(fp, &f, "F_iso_top"); fprintf(fp, "%e%c", x, f.nl);

     if ((x = xrtm_get_F_iso_bot(d)) < 0.) {
          fprintf(stderr, "ERROR: xrtm_get_F_iso_bot()\n");
          return -1;
     }
     print_lhs(fp, &f, "F_iso_bot"); fprintf(fp, "%e%c", x, f.nl);


     if (options & XRTM_OPTION_SOURCE_SOLAR) {
          if ((x = xrtm_get_F_0(d)) < 0.) {
               fprintf(stderr, "ERROR: xrtm_get_F_0()\n");
               return -1;
          }
          print_lhs(fp, &f, "F_0"); fprintf(fp, "%e%c", x, f.nl);

          if ((x = xrtm_get_theta_0(d)) < 0.) {
               fprintf(stderr, "ERROR: xrtm_get_theta_0()\n");
               return -1;
          }
          print_lhs(fp, &f, "theta_0"); fprintf(fp, "%e%c", x, f.nl);

          if ((x = xrtm_get_phi_0(d)) < 0.) {
               fprintf(stderr, "ERROR: xrtm_get_phi_0()\n");
               return -1;
          }
          print_lhs(fp, &f, "phi_0"); fprintf(fp, "%e%c", x, f.nl);
     }

     if (options & XRTM_OPTION_SOURCE_THERMAL) {
          if (print_array_param_double(fp, &f, d, "levels_b", n_layers + 1, xrtm_get_levels_b)) {
               fprintf(stderr, "ERROR: print_array_param_double()\n");
               return -1;
          }

          if ((x = xrtm_get_surface_b(d)) < 0.) {
               fprintf(stderr, "ERROR: xrtm_get_surface_b()\n");
               return -1;
          }
          print_lhs(fp, &f, "surface_b"); fprintf(fp, "%e%c", x, f.nl);
     }

     if (solvers & XRTM_SOLVERS_USE_G) {
          if (print_layer_param_double(fp, &f, d, "g", xrtm_get_g)) {
               fprintf(stderr, "ERROR: print_layer_param_double()\n");
               return -1;
          }

          if (n_derivs > 0) {
               if (print_layer_param_l_double(fp, &f, d, "g_l", xrtm_get_g_l)) {
                    fprintf(stderr, "ERROR: print_layer_param_l_double()\n");
                    return -1;
               }
          }
     }

     if (print_layer_param_double(fp, &f, d, "omega", xrtm_get_omega)) {
          fprintf(stderr, "ERROR: print_layer_param_double()\n");
          return -1;
     }

     if (n_derivs > 0) {
          if (print_layer_param_l_double(fp, &f, d, "omega_l", xrtm_get_omega_l)) {
               fprintf(stderr, "ERROR: print_layer_param_l_double()\n");
               return -1;
          }
     }

     if (print_layer_param_double(fp, &f, d, "ltau", xrtm_get_ltau)) {
          fprintf(stderr, "ERROR: print_layer_param_double()\n");
          return -1;
     }

     if (n_derivs > 0) {
          if (print_layer_param_l_double(fp, &f, d, "ltau_l", xrtm_get_ltau_l)) {
               fprintf(stderr, "ERROR: print_layer_param_l_double()\n");
               return -1;
          }
     }

     if (print_coef(fp, &f, d, md)) {
          fprintf(stderr, "ERROR: print_coef()\n");
          return -1;
     }

     if (n_derivs > 0) {
          if (print_coef_l(fp, &f, d, md)) {
               fprintf(stderr, "ERROR: print_coef_l()\n");
               return -1;
          }
     }

     for (i = 0; i < n_kernels; ++i) {
          if ((x = xrtm_get_kernel_ampfac(d, i)) < 0.) {
               fprintf(stderr, "ERROR: xrtm_get_kernel_ampfac(), i_kernel = %d", i);
               return -1;
          }
          print_lhs_dimens(fp, &f, "kernel_ampfac", 1, &i); fprintf(fp, "%e%c", x, f.nl);
     }

     if (n_derivs > 0) {
          for (i = 0; i < n_kernels; ++i) {
               if (print_kernel_ampfac_l(fp, &f, d, i)) {
                    fprintf(stderr, "ERROR: print_kernel_ampfac_l()\n");
                    return -1;
               }
          }
     }

     for (i = 0; i < n_kernels; ++i) {
          if (kernel_n_params((enum xrtm_kernel_type) xrtm_get_kernel(d, i))) {
               if (print_kernel_params(fp, &f, d, i)) {
                    fprintf(stderr, "ERROR: print_kernel_params()\n");
                    return -1;
               }
          }
     }

     if (n_derivs > 0) {
          for (i = 0; i < n_kernels; ++i) {
               if (kernel_n_params((enum xrtm_kernel_type) xrtm_get_kernel(d, i))) {
                    if (print_kernel_params_l(fp, &f, d, i)) {
                         fprintf(stderr, "ERROR: print_kernel_params_l()\n");
                         return -1;
                    }
               }
          }
     }

     return 0;
}
