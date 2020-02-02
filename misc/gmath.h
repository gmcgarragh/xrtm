/*******************************************************************************
**
**    Copyright (C) 1998-2020 Greg McGarragh <greg.mcgarragh@colostate.edu>
**
**    This source code is licensed under the GNU General Public License (GPL),
**    Version 3.  See the file COPYING for more details.
**
*******************************************************************************/

#ifndef GMATH_H
#define GMATH_H

#include "gutil.h"

#ifdef __cplusplus
extern "C" {
#endif


typedef struct {
     double m11; double m12;
     double m21; double m22;

     double dx;  double dy;
} affine2d;


typedef struct {
     double m11; double m12; double m13;
     double m21; double m22; double m23;
     double m31; double m32; double m33;

     double dx;  double dy;  double dz;
} affine3d;


/* **** gmath_affine.c **** */

void affine2d_init(affine2d *);
void affine3d_init(affine3d *);
void affine2d_copy(affine2d *, affine2d *);
void affine3d_copy(affine3d *, affine3d *);
int affine2d_equality(affine2d *, affine2d *);
int affine3d_equality(affine3d *, affine3d *);
void affine2d_print(affine2d *);
void affine3d_print(affine3d *);
void affine2d_multiply(affine2d *, affine2d *, affine2d *);
void affine3d_multiply(affine3d *, affine3d *, affine3d *);
int affine2d_invert(affine2d *, affine2d *);
int affine3d_invert(affine3d *, affine3d *);
/*
void affine2d_translate(affine2d *, double, double);
void affine3d_translate(affine3d *, double, double, double);
void affine2d_scale(affine2d *, double, double);
void affine3d_scale(affine3d *, double, double, double);
void affine2d_shear(affine2d *, double, double);
void affine3d_shear(affine3d *, double, double, double, double, double, double);
void affine2d_rotate(affine2d *, double a);
void affine3d_rotate(affine3d *, double, double, double);
*/
void affine2d_translate_x(affine2d *, double);
void affine2d_translate_y(affine2d *, double);
void affine2d_translate(affine2d *, double, double);
void affine3d_translate_x(affine3d *, double);
void affine3d_translate_y(affine3d *, double);
void affine3d_translate_z(affine3d *, double);
void affine3d_translate(affine3d *, double, double, double);
void affine2d_scale_x(affine2d *, double);
void affine2d_scale_y(affine2d *, double);
void affine2d_scale(affine2d *, double, double);
void affine3d_scale_x(affine3d *, double);
void affine3d_scale_y(affine3d *, double);
void affine3d_scale_z(affine3d *, double);
void affine3d_scale(affine3d *, double, double, double);
void affine2d_shear_g(affine2d *, double);
void affine2d_shear_h(affine2d *, double);
void affine2d_shear(affine2d *, double, double);
void affine3d_shear_xy(affine3d *, double);
void affine3d_shear_xz(affine3d *, double);
void affine3d_shear_yx(affine3d *, double);
void affine3d_shear_yz(affine3d *, double);
void affine3d_shear_zx(affine3d *, double);
void affine3d_shear_zy(affine3d *, double);
void affine3d_shear(affine3d *, double, double, double, double, double, double);
void affine2d_rotate(affine2d *, double a);
void affine3d_rotate_x(affine3d *, double);
void affine3d_rotate_y(affine3d *, double);
void affine3d_rotate_z(affine3d *, double);
void affine3d_rotate(affine3d *, double, double, double);
void map2d_1(affine2d *, double *, double *);
void map2d_2(affine2d *, double, double, double *, double *);
void map3d_1(affine3d *, double *, double *, double *);
void map3d_2(affine3d *, double, double, double, double *, double *, double *);
void affine3d_read(FILE *, affine3d *);
void affine3d_write(FILE *, affine3d *);


/* **** gmath_geom.c **** */

double prnval(double);
double point_line_min_dist(double, double, double, double, double, double);
void xyztorpz(double, double, double, double *, double *, double *);
void rpztoxyz(double *, double *, double *, double, double, double);
void xyztortp(double, double, double, double *, double *, double *);
void rtptoxyz(double *, double *, double *, double, double, double);
void xyztorll(double, double, double, double *, double *, double *);
void rlltoxyz(double *, double *, double *, double, double, double);
void rpztortp(double, double, double, double *, double *, double *);
void rtptorpz(double *, double *, double *, double, double, double);
void rpztorll(double, double, double, double *, double *, double *);
void rlltorpz(double *, double *, double *, double, double, double);
void rtptorll(double, double, double, double *, double *, double *);
void rlltortp(double *, double *, double *, double, double, double);
double great_circle_dist(double, double, double, double, double);


/* **** gmath_interp.c **** */

int lin_interp_first(double *, double, double *, int, int *, double *, double *);
int lin_interp_sorted(double *, double, double *, int, int *, double *, double *);
int lin_interp_sorted2(double *x, int n, double *xx, int nn, double *y, int *i0, double *f0, double *yy);
int lin_interp_bsearch(double *, double, double *, int, int *, double *, double *);
int lin_interp_bsearch2(double *x, int n, double *xx, int nn, double *y, int *i0, double *f0, double *yy);
int log_interp_sorted(double *, double, double *, int, int *, double *, double *);
int log_interp_bsearch(double *, double, double *, int, int *, double *, double *);
double lin_int(double, double, double, double, double);
double log_int(double, double, double, double, double, double);
int neville(double *, double *, int, double, double *, double *, double *, double *);


/* **** gmath_mat2d3d.c **** */

void vec2d_zero(double [2]);
void vec2d_copy(double [2], double [2]);
void vec2d_print(double [2]);
void vec2d_scale(double, double [2], double [2]);
void vec2d_add(double [2], double [2], double [2]);
void vec2d_sub(double [2], double [2], double [2]);
void vec2d_lincmb(double [2], double, double [2], double, double [2]);
double vec2d_dot(double [2], double [2]);
double vec2d_mag(double [2]);
void vec2d_unit(double [2], double [2]);
void vec2d_cross(double [2], double [2], double [2]);
void vec2d_dyadic(double [2], double [2], double [2][2]);

void vec3d_zero(double [3]);
void vec3d_copy(double [3], double [3]);
void vec3d_print(double [3]);
void vec3d_scale(double, double [3], double [3]);
void vec3d_add(double [3], double [3], double [3]);
void vec3d_sub(double [3], double [3], double [3]);
void vec3d_lincmb(double [3], double, double [3], double, double [3]);
double vec3d_dot(double [3], double [3]);
double vec3d_mag(double [3]);
void vec3d_unit(double [3], double [3]);
void vec3d_cross(double [3], double [3], double [3]);
void vec3d_dyadic(double [3], double [3], double [3][3]);

void mat2d_zero(double [2][2]);
void mat2d_init(double [2][2]);
void mat2d_copy(double [2][2], double [2][2]);
void mat2d_print(double [2][2]);
void mat2d_trans(double [2][2], double [2][2]);
void mat2d_add(double [2][2], double [2][2], double [2][2]);
void mat2d_sub(double [2][2], double [2][2], double [2][2]);
void mat2d_scale(double, double [2][2], double [2][2]);
void mat2d_matvec(double [2][2], double [2], double [2]);
void mat2d_mul(double [2][2], double [2][2], double [2][2]);
int mat2d_inv(double [2][2], double [2][2]);
double mat2d_det(double [2][2]);
void mat2d_rot(double [2][2], double);

void mat3d_zero(double [3][3]);
void mat3d_init(double [3][3]);
void mat3d_copy(double [3][3], double [3][3]);
void mat3d_print(double [3][3]);
void mat3d_trans(double [3][3], double [3][3]);
void mat3d_add(double [3][3], double [3][3], double [3][3]);
void mat3d_sub(double [3][3], double [3][3], double [3][3]);
void mat3d_scale(double, double [3][3], double [3][3]);
void mat3d_matvec(double [3][3], double [3], double [3]);
void mat3d_mul(double [3][3], double [3][3], double [3][3]);
int mat3d_inv(double [3][3], double [3][3]);
double mat3d_det(double [3][3]);
void mat3d_rot_x(double [3][3], double);
void mat3d_rot_y(double [3][3], double);
void mat3d_rot_z(double [3][3], double);
void mat3d_rot(double [3][3], double, double, double);


/* **** gutil_stats.c **** */

#define all_stats_proto(postfix_, type_)					\
type_ XCAT(min_, postfix_)(type_, type_);					\
type_ XCAT(min_1, postfix_)(type_ *, long);					\
type_ XCAT(min_2, postfix_)(type_ *, long, long, long);				\
type_ XCAT(min_index_1, postfix_)(type_ *, long *, long);			\
long  XCAT(min_pos_1, postfix_)(type_ *, long);					\
type_ XCAT(min_pos_2, postfix_)(type_ *, long, long, long, long *, long *);	\
long  XCAT(min_pos_index_1, postfix_)(type_ *, long *, long);			\
type_ XCAT(max_, postfix_)(type_, type_);					\
type_ XCAT(max_1, postfix_)(type_ *, long);					\
type_ XCAT(max_2, postfix_)(type_ *, long, long, long);				\
type_ XCAT(max_index_1, postfix_)(type_ *, long *, long);			\
long  XCAT(max_pos_1, postfix_)(type_ *, long);					\
type_ XCAT(max_pos_2, postfix_)(type_ *, long, long, long, long *, long *);	\
long XCAT(max_pos_index_1, postfix_)(type_ *, long *, long);			\
double XCAT(mean_1, postfix_)(type_ *, long);					\
double XCAT(mean_2, postfix_)(type_ *, long, long, long);			\
double XCAT(mean_index_1, postfix_)(type_ *, long *, long);			\
double XCAT(variance_1, postfix_)(type_ *, long, int);				\
double XCAT(variance_2, postfix_)(type_ *, long, long, long, int);		\
double XCAT(variance_index_1, postfix_)(type_ *, long *, long, int);		\
double XCAT(stddev_1, postfix_)(type_ *, long, int);				\
double XCAT(stddev_2, postfix_)(type_ *, long, long, long, int);		\
double XCAT(stddev_index_1, postfix_)(type_ *, long *, long, int);

all_stats_proto(uc,  uchar)
all_stats_proto(c,   char)
all_stats_proto(us,  ushort)
all_stats_proto(s,   short)
all_stats_proto(ui,  uint)
all_stats_proto(i,   int)
all_stats_proto(ul,  ulong)
all_stats_proto(l,   long)
all_stats_proto(ull, uint64_t)
all_stats_proto(ll,  int64_t)
all_stats_proto(f,   float)
all_stats_proto(d,   double)


/* **** gmath_util.c **** */

float gnanf(void);
double gnand(void);
long double gnanl(void);
double grint(double x);
int prec_zero(double);
int prec_equal(double, double);


#ifdef __cplusplus
}
#endif

#endif /* GMATH_H */
