/*******************************************************************************
**
**    Copyright (C) 2007-2020 Greg McGarragh <greg.mcgarragh@colostate.edu>
**
**    This source code is licensed under the GNU General Public License (GPL),
**    Version 3.  See the file COPYING for more details.
**
*******************************************************************************/

#ifndef RTUTIL_SUPPORT_H
#define RTUTIL_SUPPORT_H

#include "rtutil_size_dist.h"

#ifdef __cplusplus
extern "C" {
#endif


enum quadrature_type {
     QUAD_NORM_GAUS_LEG,
     QUAD_DOUB_GAUS_LEG,
     QUAD_LOBATTO,

     N_QUADRATURE_TYPES
};


enum scat_coef_type {
     SCAT_COEF_GC,
     SCAT_COEF_LC,
     SCAT_COEF_NONE,

     N_SCAT_COEF_TYPES
};


int get_value_list_s(const char *name, char *list, char ***array);
int put_value_list_s(char **array, int m, char *list, int length);
int get_value_list2_s(const char *name, char *list, char ****array, int *m_, int *n_);
int put_value_list2_s(char ***array, int m, int n, char *list, int length);
int get_value_list_i(const char *name, char *list, int **array);
int put_value_list_i(int *array, int m, char *list, int length);
int get_value_list2_i(const char *name, char *list, int ***array, int *m_, int *n_);
int put_value_list2_i(int **array, int m, int n, char *list, int length);
int get_value_list_d(const char *name, char *list, double **array);
int put_value_list_d(double *array, int m, char *list, int length, const char *format);
int get_value_list2_d(const char *name, char *list, double ***array, int *m_, int *n_);
int put_value_list2_d(double **array, int m, int n, char *list, int length, const char *format);

char **get_derivs_list_s(int n_derivs, const char *name, char *list);
double *get_derivs_list_d(int n_derivs, const char *name, char *list);

size_t size_dist_name_max_length(void);
enum size_dist_type size_dist_code(const char *name);
const char *size_dist_name(enum size_dist_type code);

enum quadrature_type quadrature_code(const char *name);
const char *quadrature_name(enum quadrature_type code);

enum scat_coef_type scat_coef_code(const char *name);
const char *scat_coef_name(enum scat_coef_type code);

int parse_dist_values(int argc, char *argv[], int *i_,
                      char *dist_name, enum size_dist_type *dist_type,
                      double *a1, double *a2, double *a3, double *a4,
                      double *a5, double *r1, double *r2, int flag);
int parse_dist_values_l(int argc, char *argv[], int *i_,
                        enum size_dist_type dist_type,
                        double **a1, double **a2, double **a3, double **a4,
                        double **a5, double **r1, double **r2, int n_derivs,
                        int flag);


#ifdef __cplusplus
}
#endif

#endif /* RTUTIL_SUPPORT_H */
