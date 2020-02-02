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

#include "rtutil_size_dist.h"
#include "rtutil_support.h"


/*******************************************************************************
 *
 ******************************************************************************/
static int count_list_entries(char *list) {

     char *lasts;

     int i;

     lasts = NULL;

     i = 0;
     strtok_r(list, " \n,", &lasts);
     do {
          i++;
     } while ((strtok_r(NULL, " \n,", &lasts)));

     return i;
}



static int count_list_entries2(char *list, int *m_, int *n_) {

     char *token;
     char *lasts;

     int i;
     int j;

     int n;

     lasts = NULL;

     i =  0;
     j = -1;
     token = strtok_r(list, " \n;", &lasts);
     do {
          i++;

          n = count_list_entries(token);

          if (j != -1 && n != j) {
               fprintf(stderr, "ERROR: leading dimensions do not match\n");
               return -1;
          }
          else
               j = n;
     } while ((token = strtok_r(NULL, " \n;", &lasts)));

     *m_ = i;
     *n_ = j;

     return 0;
}



static void parse_list_s(const char *name, char *list, char **array) {

     char *token;
     char *lasts;

     int i;

     lasts = NULL;

     i = 0;
     token = strtok_r(list, " \n,", &lasts);
     do {
          array[i] = strdup(token);
          i++;
     } while ((token = strtok_r(NULL, " \n,", &lasts)));
}



int get_value_list_s(const char *name, char *list, char ***array) {

     char *temp;

     int n;

     temp = strdup(list);

     n = count_list_entries(temp);

     free(temp);

     if (n == 0)
          return 0;

     *array = (char **) malloc(n * sizeof(char *));

     parse_list_s(name, list, *array);

     return n;
}



int put_value_list_s(char **array, int m, char *list, int length) {

     int i;

     int nn = 0;

     list[0] = '\0';

     for (i = 0; i < m; ++i)
          nn += snprintf(list + nn, length - nn, "%s,", array[i]);

     if (nn > 0)
          list[nn-1] = '\0';

     return nn - 1;
}



int get_value_list2_s(const char *name, char *list, char ****array, int *m_, int *n_) {

     char *temp;

     char *token;
     char *lasts;

     int i;
     int m;
     int n;

     temp = strdup(list);

     if (count_list_entries2(temp, &m, &n) < 0) {
          fprintf(stderr, "ERROR: count_list_entries2()\n");
          return -1;
     }

     free(temp);

     if (m == 0 || n == 0)
          return -1;

     *array = (char ***) alloc_array2(m, n, sizeof(char *));

     lasts = NULL;

     i = 0;
     token = strtok_r(list, " \n;", &lasts);
     do {
          parse_list_s(name, token, (*array)[i]);
          i++;
     } while ((token = strtok_r(NULL, " \n;", &lasts)));

     *m_ = m;
     *n_ = n;

     return 0;
}



int put_value_list2_s(char ***array, int m, int n, char *list, int length) {

     int i;

     int nn = 0;

     for (i = 0; i < m; ++i) {
          nn += put_value_list_s(array[i], n, list + nn, length - nn);
          nn += snprintf(list + nn, length - nn, ";");
     }

     if (nn > 0)
          list[nn-1] = '\0';

     return nn - 1;
}



static void parse_list_i(const char *name, char *list, int *array) {

     char *token;
     char *lasts;

     int i;

     lasts = NULL;

     i = 0;
     token = strtok_r(list, " \n,", &lasts);
     do {
          array[i] = strtoi_errmsg_exit(token, name);
          i++;
     } while ((token = strtok_r(NULL, " \n,", &lasts)));
}



int get_value_list_i(const char *name, char *list, int **array) {

     char *temp;

     int n;

     temp = strdup(list);

     n = count_list_entries(temp);

     free(temp);

     if (n == 0)
          return 0;

     *array = (int *) malloc(n * sizeof(int));

     parse_list_i(name, list, *array);

     return n;
}



int put_value_list_i(int *array, int m, char *list, int length) {

     int i;

     int nn = 0;

     list[0] = '\0';

     for (i = 0; i < m; ++i)
          nn += snprintf(list + nn, length - nn, "%d,", array[i]);

     if (nn > 0)
          list[nn-1] = '\0';

     return nn - 1;
}



int get_value_list2_i(const char *name, char *list, int ***array, int *m_, int *n_) {

     char *temp;

     char *token;
     char *lasts;

     int i;
     int m;
     int n;

     temp = strdup(list);

     if (count_list_entries2(temp, &m, &n) < 0) {
          fprintf(stderr, "ERROR: count_list_entries2()\n");
          return -1;
     }

     free(temp);

     if (m == 0 || n == 0)
          return -1;

     *array = alloc_array2_i(m, n);

     lasts = NULL;

     i = 0;
     token = strtok_r(list, " \n;", &lasts);
     do {
          parse_list_i(name, token, (*array)[i]);
          i++;
     } while ((token = strtok_r(NULL, " \n;", &lasts)));

     *m_ = m;
     *n_ = n;

     return 0;
}



int put_value_list2_i(int **array, int m, int n, char *list, int length) {

     int i;

     int nn = 0;

     for (i = 0; i < m; ++i) {
          nn += put_value_list_i(array[i], n, list + nn, length - nn);
          nn += snprintf(list + nn, length - nn, ";");
     }

     if (nn > 0)
          list[nn-1] = '\0';

     return nn - 1;
}



static void parse_list_d(const char *name, char *list, double *array) {

     char *token;
     char *lasts;

     int i;

     lasts = NULL;

     i = 0;
     token = strtok_r(list, " \n,", &lasts);
     do {
          array[i] = strtod_errmsg_exit(token, name);
          i++;
     } while ((token = strtok_r(NULL, " \n,", &lasts)));
}



int get_value_list_d(const char *name, char *list, double **array) {

     char *temp;

     int n;

     temp = strdup(list);

     n = count_list_entries(temp);

     free(temp);

     if (n == 0)
          return 0;

     *array = (double *) malloc(n * sizeof(double));

     parse_list_d(name, list, *array);

     return n;
}



int put_value_list_d(double *array, int m, char *list, int length, const char *format) {

     int i;

     int nn = 0;

     list[0] = '\0';

     for (i = 0; i < m; ++i) {
          nn += snprintf(list + nn, length - nn, format, array[i]);
          nn += snprintf(list + nn, length - nn, ",");
     }

     if (nn > 0)
          list[nn-1] = '\0';

     return nn - 1;
}



int get_value_list2_d(const char *name, char *list, double ***array, int *m_, int *n_) {

     char *temp;

     char *token;
     char *lasts;

     int i;
     int m;
     int n;

     temp = strdup(list);

     if (count_list_entries2(temp, &m, &n) < 0) {
          fprintf(stderr, "ERROR: count_list_entries2()\n");
          return -1;
     }

     free(temp);

     if (m == 0 || n == 0)
          return -1;

     *array = alloc_array2_d(m, n);

     lasts = NULL;

     i = 0;
     token = strtok_r(list, " \n;", &lasts);
     do {
          parse_list_d(name, token, (*array)[i]);
          i++;
     } while ((token = strtok_r(NULL, " \n;", &lasts)));

     *m_ = m;
     *n_ = n;

     return 0;
}



int put_value_list2_d(double **array, int m, int n, char *list, int length, const char *format) {

     int i;

     int nn = 0;

     for (i = 0; i < m; ++i) {
          nn += put_value_list_d(array[i], n, list + nn, length - nn, format);
          nn += snprintf(list + nn, length - nn, ";");
     }

     if (nn > 0)
          list[nn-1] = '\0';

     return nn - 1;
}



/*******************************************************************************
 *
 ******************************************************************************/
char **get_derivs_list_s(int n_derivs, const char *name, char *list) {

     int n;

     char **x;

     if (n_derivs <= 0) {
          fprintf(stderr, "ERROR: Must specifiy -derivs <n_derivs> option before %s\n",
                  name);
          exit(1);
     }

     if ((n = get_value_list_s(name, list, &x)) < 0) {
          fprintf(stderr, "ERROR: get_value_list_d()\n");
          exit(1);
     }

     if (n != n_derivs) {
          fprintf(stderr, "ERROR: Number of values in %s list is not consistent "
                  "with the specified number of derivs: %d\n", name, n_derivs);
          exit(1);
     }

     return x;
}



double *get_derivs_list_d(int n_derivs, const char *name, char *list) {

     int n;

     double *x;

     if (n_derivs <= 0) {
          fprintf(stderr, "ERROR: Must specifiy -derivs <n_derivs> option before %s\n",
                  name);
          exit(1);
     }

     if ((n = get_value_list_d(name, list, &x)) < 0) {
          fprintf(stderr, "ERROR: get_value_list_d()\n");
          exit(1);
     }

     if (n != n_derivs) {
          fprintf(stderr, "ERROR: Number of values in %s list is not consistent "
                  "with the specified number of derivs: %d\n", name, n_derivs);
          exit(1);
     }

     return x;
}



/*******************************************************************************
 *
 ******************************************************************************/
static const size_t __size_dist_name_max_length = sizeof("modified_bimodal_log_normal");

static const char *size_dist_names[] = {
     "none",
     "modified_gamma",
     "log_normal",
     "power_law",
     "gamma",
     "modified_power_law",
     "modified_bimodal_log_normal",
     "exponential",
     "mono"
};

size_t size_dist_name_max_length(void) {
     return __size_dist_name_max_length;
}

enum size_dist_type size_dist_code(const char *name) {
     return (enum size_dist_type)
          name_to_index(name, size_dist_names, N_SIZE_DIST_TYPES,
                        "size distribution type");
}

const char *size_dist_name(enum size_dist_type code) {
     return index_to_name(code, size_dist_names, N_SIZE_DIST_TYPES,
                         "size distribution type");
}



static const char *quadrature_names[] = {
     "norm_gaus_leg",
     "doub_gaus_leg",
     "lobatto"
};

enum quadrature_type quadrature_code(const char *name) {
     return (enum quadrature_type)
          name_to_index(name, quadrature_names, N_QUADRATURE_TYPES,
                        "quadrature type");
}

const char *quadrature_name(enum quadrature_type code) {
     return index_to_name(code, quadrature_names, N_QUADRATURE_TYPES,
                          "quadrature type");
}



static const char *scat_coef_names[] = {
     "gc",
     "lc",
     "none"
};

enum scat_coef_type scat_coef_code(const char *name) {
     return (enum scat_coef_type)
          name_to_index(name, scat_coef_names, N_SCAT_COEF_TYPES,
                        "phase coefficent type");
}

const char *scat_coef_name(enum scat_coef_type code) {
     return index_to_name(code, scat_coef_names, N_SCAT_COEF_TYPES,
                          "phase coefficent type");
}



/*******************************************************************************
 *
 ******************************************************************************/
int parse_dist_values(int argc, char *argv[], int *i_,
                      char *dist_name, enum size_dist_type *dist_type,
                      double *a1, double *a2, double *a3, double *a4,
                      double *a5, double *r1, double *r2, int flag) {

     char *temp;

     int length;

     const char *option = "-dist";

     int i = *i_;

     if ((*dist_type = size_dist_code(dist_name)) < 0) {
          fprintf(stderr, "ERROR: (-dist) invalid distribution name: %s\n", dist_name);
          return -1;
     }

     length = strlen(option) + 1 + strlen(dist_name) + 1;

     temp = (char *) malloc(length);

     snprintf(temp, length, "%s %s", option, dist_name);

     if (*dist_type == SIZE_DIST_MONO) {
          check_arg_count(i, argc, 1, temp);
          *a1 = strtod_errmsg_exit(argv[++i], temp);
          *r1 = *a1;
          *r2 = *a1;
     }
     else
     if (*dist_type == SIZE_DIST_MODIFIED_GAMMA) {
          check_arg_count(i, argc, 5, temp);
          *a1 = strtod_errmsg_exit(argv[++i], temp);
          *a2 = strtod_errmsg_exit(argv[++i], temp);
          *a3 = strtod_errmsg_exit(argv[++i], temp);
          *r1 = strtod_errmsg_exit(argv[++i], temp);
          *r2 = strtod_errmsg_exit(argv[++i], temp);
     }
     else
     if (*dist_type == SIZE_DIST_LOG_NORMAL) {
          check_arg_count(i, argc, 4, temp);
          *a1 = strtod_errmsg_exit(argv[++i], temp);
          *a2 = strtod_errmsg_exit(argv[++i], temp);
          *r1 = strtod_errmsg_exit(argv[++i], temp);
          *r2 = strtod_errmsg_exit(argv[++i], temp);
     }
     else
     if (*dist_type == SIZE_DIST_POWER_LAW) {
          check_arg_count(i, argc, 2, temp);
          *a1 = strtod_errmsg_exit(argv[++i], temp);
          *a2 = strtod_errmsg_exit(argv[++i], temp);
     }
     else
     if (*dist_type == SIZE_DIST_GAMMA) {
          check_arg_count(i, argc, 4, temp);
          *a1 = strtod_errmsg_exit(argv[++i], temp);
          *a2 = strtod_errmsg_exit(argv[++i], temp);
          *r1 = strtod_errmsg_exit(argv[++i], temp);
          *r2 = strtod_errmsg_exit(argv[++i], temp);
     }
     else
     if (*dist_type == SIZE_DIST_MODIFIED_POWER_LAW) {
          check_arg_count(i, argc, 3, temp);
          *a1 = strtod_errmsg_exit(argv[++i], temp);
          *r1 = strtod_errmsg_exit(argv[++i], temp);
          *r2 = strtod_errmsg_exit(argv[++i], temp);
     }
     else
     if (*dist_type == SIZE_DIST_MODIFIED_BIMODAL_LOG_NORMAL) {
          check_arg_count(i, argc, 7, temp);
          *a1 = strtod_errmsg_exit(argv[++i], temp);
          *a2 = strtod_errmsg_exit(argv[++i], temp);
          *a3 = strtod_errmsg_exit(argv[++i], temp);
          *a4 = strtod_errmsg_exit(argv[++i], temp);
          *a5 = strtod_errmsg_exit(argv[++i], temp);
          *r1 = strtod_errmsg_exit(argv[++i], temp);
          *r2 = strtod_errmsg_exit(argv[++i], temp);
     }
     else
     if (flag && *dist_type == SIZE_DIST_EXPONENTIAL) {
          check_arg_count(i, argc, 3, temp);
          *a1 = strtod_errmsg_exit(argv[++i], temp); /* lambda */
          *r1 = strtod_errmsg_exit(argv[++i], temp);
          *r2 = strtod_errmsg_exit(argv[++i], temp);
     }

     *i_ = i;

     free(temp);

     return 0;
}



int parse_dist_values_l(int argc, char *argv[], int *i_,
                        enum size_dist_type dist_type,
                        double **a1, double **a2, double **a3, double **a4,
                        double **a5, double **r1, double **r2, int n_derivs,
                        int flag) {

     int i = *i_;

     if (dist_type == SIZE_DIST_MONO) {
          check_arg_count(i, argc, 1, "-dist_l");
          *a1 = get_derivs_list_d(n_derivs, "-dist_l", argv[++i]);
     }
     else
     if (dist_type == SIZE_DIST_MODIFIED_GAMMA) {
          check_arg_count(i, argc, 5, "-dist_l");
          *a1 = get_derivs_list_d(n_derivs, "-dist_l", argv[++i]);
          *a2 = get_derivs_list_d(n_derivs, "-dist_l", argv[++i]);
          *a3 = get_derivs_list_d(n_derivs, "-dist_l", argv[++i]);
          *r1 = get_derivs_list_d(n_derivs, "-dist_l", argv[++i]);
          *r2 = get_derivs_list_d(n_derivs, "-dist_l", argv[++i]);
     }
     else
     if (dist_type == SIZE_DIST_LOG_NORMAL) {
          check_arg_count(i, argc, 4, "-dist_l");
          *a1 = get_derivs_list_d(n_derivs, "-dist_l", argv[++i]);
          *a2 = get_derivs_list_d(n_derivs, "-dist_l", argv[++i]);
          *r1 = get_derivs_list_d(n_derivs, "-dist_l", argv[++i]);
          *r2 = get_derivs_list_d(n_derivs, "-dist_l", argv[++i]);
     }
     else
     if (dist_type == SIZE_DIST_POWER_LAW) {
          check_arg_count(i, argc, 2, "-dist_l");
          *a1 = get_derivs_list_d(n_derivs, "-dist_l", argv[++i]);
          *a2 = get_derivs_list_d(n_derivs, "-dist_l", argv[++i]);
     }
     else
     if (dist_type == SIZE_DIST_GAMMA) {
          check_arg_count(i, argc, 4, "-dist_l");
          *a1 = get_derivs_list_d(n_derivs, "-dist_l", argv[++i]);
          *a2 = get_derivs_list_d(n_derivs, "-dist_l", argv[++i]);
          *r1 = get_derivs_list_d(n_derivs, "-dist_l", argv[++i]);
          *r2 = get_derivs_list_d(n_derivs, "-dist_l", argv[++i]);
     }
     else
     if (dist_type == SIZE_DIST_MODIFIED_POWER_LAW) {
          check_arg_count(i, argc, 3, "-dist_l");
          *a1 = get_derivs_list_d(n_derivs, "-dist_l", argv[++i]);
          *r1 = get_derivs_list_d(n_derivs, "-dist_l", argv[++i]);
          *r2 = get_derivs_list_d(n_derivs, "-dist_l", argv[++i]);
     }
     else
     if (dist_type == SIZE_DIST_MODIFIED_BIMODAL_LOG_NORMAL) {
          check_arg_count(i, argc, 7, "-dist_l");
          *a1 = get_derivs_list_d(n_derivs, "-dist_l", argv[++i]);
          *a2 = get_derivs_list_d(n_derivs, "-dist_l", argv[++i]);
          *a3 = get_derivs_list_d(n_derivs, "-dist_l", argv[++i]);
          *a4 = get_derivs_list_d(n_derivs, "-dist_l", argv[++i]);
          *a5 = get_derivs_list_d(n_derivs, "-dist_l", argv[++i]);
          *r1 = get_derivs_list_d(n_derivs, "-dist_l", argv[++i]);
          *r2 = get_derivs_list_d(n_derivs, "-dist_l", argv[++i]);
     }
     else
     if (flag && dist_type == SIZE_DIST_EXPONENTIAL) {
          check_arg_count(i, argc, 3, "-dist_l");
          *a1 = get_derivs_list_d(n_derivs, "-dist_l", argv[++i]);
          *r1 = get_derivs_list_d(n_derivs, "-dist_l", argv[++i]);
          *r2 = get_derivs_list_d(n_derivs, "-dist_l", argv[++i]);
     }
     else {
          fprintf(stderr, "ERROR: (-dist) invalid distribution type: %d\n", dist_type);
          return -1;
     }

     *i_ = i;

     return 0;
}
